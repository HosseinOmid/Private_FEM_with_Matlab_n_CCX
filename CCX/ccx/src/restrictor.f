!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!     
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!     
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!     
      subroutine restrictor(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &     nodef,idirf,df,cp,r,physcon,dvi,numf,set
     &     ,shcon,nshcon,rhcon,nrhcon,ntmat_,mi,iaxial)
!     
!     pressure loss element with partial total head loss 
!
!     author: Yannick Muller
!      
      implicit none
!     
      logical identity,crit,isothermal
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(5),idirf(5),index,iflag,iaxial,
     &     inv,ipkon(*),kon(*),kgas,icase,k_oil,nshcon(*),
     &     nrhcon(*),ntmat_,mi(*)
!     
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(5),kappa,R,d,
     &     Tt1,Tt2,pt1,pt2,cp,physcon(3),km1,dvi,
     &     kp1,kdkm1,reynolds,kdkp1,
     &     pt2pt1,pt1pt2,pt1pt2_crit,qred_crit,qred1,qred2,zeta,
     &     A1,A2,root, expon1,expon2,expon3,fact1,fact2,sqrt,pi,
     &     pt2_lim,M2,M1,xflow_oil,T1,T2,phi,
     &     shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),zeta_phi,Aeff,
     &     C2,tdkp1
!     
      phi=0.d0
      if (iflag.eq.0) then
         identity=.true.
!     
         if(nactdog(2,node1).ne.0)then
            identity=.false.
         elseif(nactdog(2,node2).ne.0)then
            identity=.false.
         elseif(nactdog(1,nodem).ne.0)then
            identity=.false.
         endif
!     
      elseif (iflag.eq.1)then
!     
         isothermal=.false.
         index=ielprop(nelem)
         kappa=(cp/(cp-R))
         kp1=kappa+1d0
         km1=kappa-1d0
!     
!     defining surfaces for branch elements
!     
         if(lakon(nelem)(2:6).eq.'REBRJ') then
            if(nelem.eq.int(prop(index+2))) then
               A1=prop(index+5)
               A2=A1
            elseif(nelem.eq.int(prop(index+3)))then
               A1=prop(index+6)
               A2=A1
            endif
            zeta=1.d0
         elseif(lakon(nelem)(2:6).eq.'REBRS') then
            if(nelem.eq.int(prop(index+2))) then
               A1=prop(index+5)
               A2=A1
            elseif(nelem.eq.int(prop(index+3)))then
               A1=prop(index+6)
               A2=A1
            endif
            zeta=1.d0
!     
!     for other Restrictor elements         
!     
         else if (lakon(nelem)(2:5).eq.'REUS' ) then
            A1=prop(index+1)
            A2=prop(index+2)
            zeta=prop(index+4)
            if(A1.gt.A2) then
               A1=A2
            endif
         else
            A1=prop(index+1)
            A2=prop(index+2) 
            zeta=1.d0
         endif
!     
         pt1=v(2,node1)
         pt2=v(2,node2)
!     
         if(pt1.ge.pt2) then
            inv=1
            Tt1=v(0,node1)+physcon(1)
            Tt2=v(0,node2)+physcon(1)
         else
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
            Tt1=v(0,node2)+physcon(1)
            Tt2=v(0,node1)+physcon(1)
         endif
!     
         pt1pt2=pt1/pt2
         pt2pt1=1/pt1pt2
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
!     
         if(.not.isothermal) then
            pt1pt2_crit=(0.5d0*kp1)**(zeta*kdkm1)
         else
            pt1pt2_crit=0.5d0*(3*kappa-1)**(zeta*kdkm1)
         endif
!     
         if(pt1pt2.gt.pt1pt2_crit) then
            crit=.true.
            pt1pt2=pt1pt2_crit
         endif
!     
         if(A1.le.A2) then
!     

            Qred1=dsqrt(kappa/R)*pt1pt2**(-0.5d0*kp1/(kappa*zeta))
     &           *dsqrt(2.d0/km1*(pt1pt2**(km1/(kappa*zeta))-1d0))
!     
            Qred2=pt1pt2*A1/A2*Qred1
!     
            if(.not.isothermal) then
               Qred_crit=dsqrt(kappa/R)*(1.d0+0.5d0*km1)
     &              **(-0.5d0*kp1/km1)
            else
               Qred_crit=dsqrt(1/R)*(1+0.5*km1/kappa)
     &              **(-0.5d0*kp1/km1)
            endif
!     
            if (Qred2.lt.Qred_crit) then
               if((Qred1.gt.Qred_crit).or.(pt1pt2.gt.pt1pt2_crit)) then
                  xflow=inv*A1*pt1*Qred_crit/dsqrt(Tt1)
               else
                  xflow=inv*A1*pt1*Qred1/dsqrt(Tt1)
               endif
            else
               call pt2_lim_calc(pt1,a2,a1,kappa,zeta,pt2_lim)
!     
               xflow=inv*A2*pt2_lim*Qred_crit/dsqrt(Tt2)
!     
            endif
!     
         else     
            Qred_crit=dsqrt(kappa/R)*(1.d0+0.5d0*km1)
     &              **(-0.5d0*kp1/km1)
            call pt2_lim_calc(pt1,a2,a1,kappa,zeta,pt2_lim)
!     
            xflow=inv*A2*pt2_lim*Qred_crit/dsqrt(Tt2)
         endif

         pt2pt1=pt2/pt1
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         tdkp1=2.d0/kp1
         C2=tdkp1**kdkm1
         if(A1.gt.A2) then
            Aeff=A2
         else
            Aeff=A1
         endif
         if(pt2pt1.gt.C2) then
            xflow=inv*pt1*Aeff*dsqrt(2.d0*kdkm1*pt2pt1**(2.d0/kappa)
     &           *(1.d0-pt2pt1**(1.d0/kdkm1))/r)/dsqrt(Tt1)
         else
            xflow=inv*pt1*Aeff*dsqrt(kappa/r)*tdkp1**(kp1/(2.d0*km1))/
     &           dsqrt(Tt1)
         endif
         if(lakon(nelem)(2:5).ne.'RECO') then
            xflow=0.75*xflow
         else
            xflow=xflow
         endif
!     
      elseif (iflag.eq.2)then
!     
         numf=4
         isothermal=.false.
         pi=4.d0*datan(1.d0)
         kappa=(cp/(cp-R))
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         kdkp1=kappa/kp1
         index=ielprop(nelem)
!
         pt1=v(2,node1)
         pt2=v(2,node2)
!
         if(pt1.ge.pt2) then
            inv=1
         else
            inv=-1
         endif     
!     
!     defining surfaces and oil properties for branches elements
!     
         if(lakon(nelem)(2:6).eq.'REBRJ') then
            if(nelem.eq.int(prop(index+2))) then
               A1=prop(index+5)
               A2=A1
               xflow_oil=prop(index+9)
               k_oil=int(prop(index+11))
            elseif(nelem.eq.int(prop(index+3)))then
               A1=prop(index+6)
               A2=A1
               xflow_oil=prop(index+10)
               k_oil=int(prop(index+11))
            endif
         elseif(lakon(nelem)(2:6).eq.'REBRS') then
            if(nelem.eq.int(prop(index+2))) then
               A1=prop(index+5)
               A2=A1
               if(lakon(nelem)(2:8).eq.'REBRSI1') then
                  xflow_oil=prop(index+11)
                  k_oil=int(prop(index+13))
               else
                  xflow_oil=prop(index+9)
                  k_oil=int(prop(index+11))
               endif
            elseif(nelem.eq.int(prop(index+3)))then
               A1=prop(index+6)
               A2=A1
               if(lakon(nelem)(2:8).eq.'REBRSI1') then
                  xflow_oil=prop(index+12)
                  k_oil=int(prop(index+13))
               else 
                  xflow_oil=prop(index+10)
                  k_oil=int(prop(index+11))
               endif
            endif
!     
!     for other Restrictor elements         
!     
       
         else                 
            if(inv.gt.0.d0) then
               A1=prop(index+1)
               A2=prop(index+2)
            else              
               A1=prop(index+2)
               A2=prop(index+1)              
            endif
!           
            if(lakon(nelem)(2:5).eq.'REEL') then
               xflow_oil=prop(index+4)
               k_oil=int(prop(index+5))
            elseif((lakon(nelem)(2:7).eq.'RELOID').or.
     &              (lakon(nelem)(2:5).eq.'REUS').or. 
     &              (lakon(nelem)(2:5).eq.'REEN').or.
     &              (lakon(nelem)(2:5).eq.'REEX').or.
     &              (lakon(nelem)(2:7).eq.'REWAOR').or.
     &              (lakon(nelem)(2:7).eq.'RELOLI')) then
               xflow_oil=prop(index+5)
               k_oil=int(prop(index+6))
            elseif((lakon(nelem)(2:5).eq.'RECO').or.
     &              (lakon(nelem)(2:7).eq.'REBEMA').or.
     &              (lakon(nelem)(2:7).eq.'REBEMI').or.
     &              (lakon(nelem)(2:8).eq.'REBEIDC')) then
               xflow_oil=prop(index+6)
               k_oil=int(prop(index+7))
            elseif(lakon(nelem)(2:8).eq.'REBEIDR') then
               xflow_oil=prop(index+8)
               k_oil=int(prop(index+9))
            endif
         endif
!
         if(pt1.gt.pt2) then
            inv=1
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)+physcon(1)
            Tt2=v(0,node2)+physcon(1)
!     
            icase=0
            call ts_calc(xflow,Tt1,Pt1,kappa,r,a1,T1,icase)
            call ts_calc(xflow,Tt2,Pt2,kappa,r,a2,T2,icase)
!
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
            
         elseif(pt1.eq.pt2) then
            inv=1
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)+physcon(1)
            Tt2=v(0,node2)+physcon(1)
!     
            pt2=pt2-0.01*pt2
            icase=0
            call ts_calc(xflow,Tt1,Pt1,kappa,r,a1,T1,icase)
            call ts_calc(xflow,Tt2,Pt2,kappa,r,a2,T2,icase)
!
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
!     
         else
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
            xflow=-v(1,nodem)*iaxial
            Tt1=v(0,node2)+physcon(1)
            Tt2=v(0,node1)+physcon(1)
            icase=0
            call ts_calc(xflow,Tt1,Pt1,kappa,r,a1,T1,icase)
            call ts_calc(xflow,Tt2,Pt2,kappa,r,a2,T2,icase)
            nodef(1)=node2
            nodef(2)=node2
            nodef(3)=nodem
            nodef(4)=node1 
         endif

!     
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!     
!     calculation of the dynamic viscosity 
!     
         if( lakon(nelem)(2:3).eq.'RE') then
            icase=0
         endif
!     
         if (A1.le.A2) then
            if(dabs(dvi).lt.1E-30) then
               kgas=0
               call dynamic_viscosity(kgas,T1,dvi)
            endif
         else
            if(dabs(dvi).lt.1E-30) then
               kgas=0
               call dynamic_viscosity(kgas,T2,dvi)
            endif
         endif
!     
!     Reynolds number calculation
!     
         if (lakon(nelem)(2:5).eq.'REBR') then
            d=dsqrt(4d0*A1/Pi)
            reynolds=dabs(xflow)*d/(dvi*A1)     
         else
            d=prop(index+3)
            if(A1.le.A2) then
               reynolds=dabs(xflow)*d/(dvi*A1)
            else
               reynolds=dabs(xflow)*d/(dvi*A2)
            endif
         endif

         if(xflow_oil.lt.1E-10) then
            xflow_oil=0d0
         endif 
!
!     BEND MILLER with oil
!     
         if(lakon(nelem)(2:7).eq.'REBEMI') then
            if(xflow_oil.ne.0d0) then
!
               call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &              xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &              v,dvi,cp,r,k_oil,phi,zeta,nshcon,nrhcon,
     &              shcon,rhcon,ntmat_,mi)
!
               zeta=phi*zeta
            else

               call zeta_calc(nelem,prop,ielprop,lakon,reynolds,zeta,
     &              isothermal,kon,ipkon,R,Kappa,v,mi)
               phi=1.d0
            endif
!
!     long orifice idelchick with oil
!
         elseif(lakon(nelem)(2:7).eq.'RELOID') then
            if(xflow_oil.ne.0d0) then
!               
               call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &              xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &              v,dvi,cp,r,k_oil,phi,zeta,nshcon,nrhcon,
     &              shcon,rhcon,ntmat_,mi)
               zeta=phi*zeta

            else

               call zeta_calc(nelem,prop,ielprop,lakon,reynolds,zeta,
     &              isothermal,kon,ipkon,R,Kappa,v,mi)
               phi=1.d0
            endif
!
!     every other zeta elements with/without oil
!
         else               
!    
            if(xflow_oil.ne.0d0) then
               call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &              xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &              v,dvi,cp,r,k_oil,phi,zeta,nshcon,nrhcon,
     &              shcon,rhcon,ntmat_,mi)
               call zeta_calc(nelem,prop,ielprop,lakon,reynolds,zeta,
     &              isothermal,kon,ipkon,R,Kappa,v,mi)
               zeta=phi*zeta
            else
               phi=1.d0
               call zeta_calc(nelem,prop,ielprop,lakon,reynolds,zeta,
     &              isothermal,kon,ipkon,R,Kappa,v,mi)
               zeta=phi*zeta
            endif
         endif
!
         if(zeta.lt.0) then
            pt1=v(2,node1)
            pt2=v(2,node2)
            xflow=v(1,nodem)*iaxial
            Tt2=v(0,node2)
            Tt1=v(0,node1)
           call ts_calc(xflow,Tt1,Pt1,kappa,r,A1,T1,icase)
           call ts_calc(xflow,Tt2,Pt2,kappa,r,A2,T2,icase)
!     
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
!     
         endif
!     
         if(.not.isothermal) then
            pt1pt2_crit=(0.5d0*kp1)**(zeta*kdkm1)
         else
            pt1pt2_crit=0.5d0*(3*kappa-1)**(zeta*kdkm1)
         endif
         pt1pt2=pt1/pt2
!     
!     Mach number caclulation
!    
         M1=dsqrt(2d0/km1*(Tt1/T1-1d0))
         if((1.d0-M1).le.1E-6) then
            if(zeta.gt.0d0) then
               call limit_case_calc(a2,pt1,Tt2,xflow,zeta,r,kappa,
     &              pt2_lim,M2)
!     
            endif
         else
            M2=dsqrt(2d0/km1*(Tt2/T2-1d0))
         endif
!     
!     Section A1 smaller than or equal to section A2
!     or for all BRANCHES ELEMENTS
!     
         if (A1.le.A2) then
!     
!     definition of the reduced mass flows
!     
            if(zeta.gt.0) then
!     
               Qred1=dsqrt(kappa/R)*pt1pt2**(-0.5d0*kp1/(kappa*zeta))
     &              *dsqrt(2.d0/km1*(pt1pt2**(km1/(kappa*zeta))-1d0))
!     
            elseif(zeta.lt.0d0) then
!     
               Qred1=dabs(xflow)*dsqrt(Tt1)/(pt1*A1)
!     
            endif
!     
            Qred2=pt1pt2*A1/A2*Qred1
!     
            if(.not.isothermal) then
               Qred_crit=dsqrt(kappa/R)*(1.d0+0.5d0*km1)
     &              **(-0.5d0*kp1/km1)
            else
               Qred_crit=dsqrt(1/R)*(1+0.5*km1/kappa)
     &              **(-0.5d0*kp1/km1)
            endif
!     
!     icase zeta greater than zero
!     
            if(zeta.gt.0) then
!     
!     definition of the coefficients 
!     
               sqrt=dsqrt(R*Tt1/kappa)
               expon1=-0.5d0*kp1/(zeta*kappa)
               fact1=pt1pt2**expon1
               expon2=km1/(zeta*kappa)
               fact2=pt1pt2**expon2
               expon3=1d0/(zeta*kappa)
               root=2d0/km1*(fact2-1d0)
!     
               if(Qred2.lt.Qred_crit) then
!     
                  if((Qred1.lt.Qred_crit)
     &                 .and.(pt1pt2.lt.pt1pt2_crit))then
!     
!     section 1 is not critical
!     
!     residual
!     
                     f=xflow*sqrt/(A1*Pt1)-fact1*dsqrt(root)
!     
!     pressure node1
!     
                     df(1)=-xflow*sqrt/(A1*Pt1**2)+
     &                    fact1/pt1*dsqrt(root)
     &                    *(-expon1-expon3*fact2/root)
!     
!     temperature node1
!     
                     df(2)=0.5d0*xflow*dsqrt(R/(kappa*Tt1))/(A1*Pt1)
!     
!     mass flow
!     
                     df(3)=inv*sqrt/(A1*Pt1)
!     
!     pressure node2
!     
                     df(4)=fact1/pt2*dsqrt(root)*
     &                    (expon1+expon3*fact2/root)
!     
                  else
!     
!     section1 is critical
!     
                     f=xflow*sqrt/(pt1*A1)-dsqrt(R/kappa)*qred_crit
!     
!     pressure node1
!     
                     df(1)=-xflow*sqrt/(A1*pt1**2)
!     
!     temperature node1
!     
                     df(2)=0.5d0*xflow*dsqrt(R/kappa)
     &                    /(pt1*A1*dsqrt(Tt1))
!     
!     mass flow
!     
                     df(3)=inv*sqrt/(A1*pt1)
!     
!     pressure node2
!     
                     df(4)=0.d0
!     
                  endif
!     
               else
!     
!     section A2 critical
!
                  call pt2_lim_calc(pt1,a2,a1,kappa,zeta,pt2_lim)
                  pt1pt2=pt1/pt2_lim
!     
                  fact1=pt1pt2**expon1
!     
                  fact2=pt1pt2**expon2
!     
                  root=2d0/km1*(fact2-1d0)
!     
                  f=xflow*sqrt/(A1*Pt1)-fact1*dsqrt(root)
!     
!     pressure node1
!     
                  df(1)=-xflow*sqrt/(A1*Pt1**2)+
     &                 fact1/pt1*dsqrt(root)
     &                 *(-expon1-expon3*fact2/root)
!     
!     temperature node1
!     
                  df(2)=0.5d0*xflow*dsqrt(R/(kappa*Tt1))/(A1*Pt1)
!     
!     mass flow
!     
                  df(3)=inv*sqrt/(A1*Pt1)
!     
!     pressure node2
!     
                  df(4)=0
!     
               endif
!     
!     icase zeta less than zero
!     
            elseif(zeta.lt.0) then
!     
               expon1=-kp1/(zeta*kappa)
               fact1=pt1pt2**expon1
               expon2=km1/(zeta*kappa)
               fact2=pt1pt2**expon2
               expon3=1d0/(zeta*kappa)
               root=2d0/km1*(fact2-1d0) 
!     
               if(Qred1.lt.Qred_crit) then
!     
!     section 1 is not critical
!     
!     residual
!     
                  f=xflow**2*R*Tt1/(A1**2*Pt1**2*Kappa)
     &                 -fact1*root
!     
!     pressure node1
!     
                  df(1)=-2*xflow**2*R*Tt1/(A1**2*Pt1**3*Kappa)
     &                 -1/pt1*fact1*(expon1*root
     &                 +2/(zeta*kappa)*fact2)
!     
!     temperature node1
!     
                  df(2)=xflow**2*R/(A1**2*Pt1**2*Kappa)
!     
!     mass flow
!     
                  df(3)=2*xflow*R*Tt1/(A1**2*Pt1**2*Kappa)
!     
!     pressure node2
!     
                  df(4)=-(1/Pt2*fact1)
     &                 *(-expon1*root-2/(zeta*kappa)*fact2)
!     
!     section1 is critical
!     
               else
!     
                  f=xflow**2*R*Tt1/(A1**2*Pt1**2*Kappa)
     &                 -R/kappa*qred_crit**2
!     
!     pressure node1
!     
                  df(1)=-2*xflow**2*R*Tt1/(A1**2*pt1**3*kappa)
!     
!     temperature node1
!     
                  df(2)=xflow**2*R/(A1**2*Pt1**2*Kappa)
!     
!     mass flow
!     
                  df(3)=2*xflow*R*Tt1/(A1**2*Pt1**2*Kappa)
!     
!     pressure node2
!     
                  df(4)=0.d0
!     
               endif
!
!     zeta = 0
!
            elseif(zeta.eq.0d0) then
!     
               f=pt1-pt2
!     
               df(1)=1
!     
               df(2)=0
!     
               df(3)=0
!     
               df(4)=-1
!     
            endif
!     
         else
!     
!     A1 greater than A2 
!     
            Qred2=dabs(xflow)*dsqrt(Tt2)/(A2*Pt2)
!     
            Qred1=1/pt1pt2*A2/A1*Qred2
!
            Qred_crit=dsqrt(kappa/R)*(1.d0+0.5d0*km1)
     &           **(-0.5d0*kp1/km1)

!     definition of the coefficients 
!     
            if(zeta.gt.0d0) then
!     
               sqrt=dsqrt(R*Tt1/kappa)
!
               expon1=-0.5d0*kp1/(zeta*kappa)
               fact1=pt1pt2**expon1
               expon2=km1/(zeta*kappa)
               fact2=pt1pt2**expon2
               expon3=1d0/(zeta*kappa)
               root=2d0/km1*(fact2-1d0)
!     
               if(pt1pt2.ge.pt1pt2_crit) then
                  pt1pt2=pt1pt2_crit
                  pt2=pt1/pt1pt2_crit
               endif
!     
               if((Qred2.lt.Qred_crit)
     &              .and.(pt1/pt2.lt.pt1pt2_crit)) then
!     
!     section 2 is not critical
!     
!     residual
!     
                  f=xflow*sqrt/(A2*Pt2)-fact1*dsqrt(root)
!     
!     pressure node1
!     
                  df(1)=-fact1/pt1*dsqrt(root)
     &                 *(expon1+0.5*dsqrt(2/km1)*expon2*fact2/root)
!     
!     temperature node1
!     
                  df(2)=0.5d0*xflow*sqrt/(A2*Pt2*Tt1)
!     
!     mass flow
!     
                  df(3)=inv*sqrt/(A2*Pt2)
!     
!     pressure node2
!     
                  df(4)=-xflow*sqrt/(A2*Pt2**2)
     &                 -fact1/pt2*dsqrt(root)*
     &                 (-expon1-0.5*dsqrt(2/km1)*expon2*fact2/root)
!     
               else
                  write(*,*) 
     &               '*WARNING in restrictor: A1 greater A2 critical'
!     
!     section2 is critical
!     
                  pt2=pt1/pt1pt2_crit
!     
                  f=xflow*dsqrt(Tt1)/(pt2*A2)-qred_crit
!     
!     pressure node1
!     
                  df(1)=0
!     
!     temperature node1
!     
                  df(2)=0.5d0*xflow/(A2*pt2*dsqrt(Tt2))
!     
!     mass flow
!     
                  df(3)=inv*dsqrt(Tt1)/(A2*pt2)
!     
!     pressure node2
!     
                  df(4)=-xflow*dsqrt(Tt1)/(A2*pt2**2)
!     
               endif
!     
            elseif(zeta.eq.0d0) then
!     
               Qred1=dabs(xflow)*dsqrt(Tt1*kappa/R)/(A1*Pt1)
               Qred2=dabs(xflow)*dsqrt(Tt2*kappa/R)/(A2*Pt2)
               Qred_crit=dsqrt(kappa/R)*(1.d0+0.5d0*km1)
     &              **(-0.5d0*kp1/km1)
!     
               f=pt1/pt2-1.d0
!     
               df(1)=1/pt2
!     
               df(2)=0
!     
               df(3)=0
!     
               df(4)=-pt1/pt2**2
!     
            endif
         endif
!     
      elseif((iflag.eq.3).or.(iflag.eq.4)) then
!
         isothermal=.false.
         pi=4.d0*datan(1.d0)
         kappa=(cp/(cp-R))
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         kdkp1=kappa/kp1
         index=ielprop(nelem)
!     
         pt1=v(2,node1)
         pt2=v(2,node2)
         if(pt1.ge.pt2) then
            inv=1
         else
            inv=-1
         endif
!     
!     defining surfaces for branches elements
!     
         if(lakon(nelem)(2:6).eq.'REBRJ') then
            if(nelem.eq.int(prop(index+2))) then
               A1=prop(index+5)
               A2=A1
               xflow_oil=prop(index+9)
               k_oil=int(prop(index+11))
            elseif(nelem.eq.int(prop(index+3)))then
               A1=prop(index+6)
               A2=A1
               xflow_oil=prop(index+10)
               k_oil=int(prop(index+11))
            endif
         elseif(lakon(nelem)(2:6).eq.'REBRS') then
            if(nelem.eq.int(prop(index+2))) then
               A1=prop(index+5)
               A2=A1
               if(lakon(nelem)(2:8).eq.'REBRSI1') then
                  xflow_oil=prop(index+11)
                  k_oil=int(prop(index+13))
               else
                  xflow_oil=prop(index+9)
                  k_oil=int(prop(index+11))
               endif
            elseif(nelem.eq.int(prop(index+3)))then
               A1=prop(index+6)
               A2=A1
               if(lakon(nelem)(2:8).eq.'REBRSI1') then
                  xflow_oil=prop(index+12)
                  k_oil=int(prop(index+13))
               else 
                  xflow_oil=prop(index+10)
                  k_oil=int(prop(index+11))
               endif
            endif
!     
!     for other Restrictor elements         
!     
         else
            A1=prop(index+1)
            A2=prop(index+2)
            if(lakon(nelem)(2:5).eq.'REEL') then
               xflow_oil=prop(index+4)
               k_oil=int(prop(index+5))
            elseif((lakon(nelem)(2:7).eq.'RELOID').or.
     &              (lakon(nelem)(2:5).eq.'REUS').or.  
     &              (lakon(nelem)(2:5).eq.'REEN').or.  
     &              (lakon(nelem)(2:5).eq.'REEX').or.   
     &              (lakon(nelem)(2:7).eq.'REWAOR').or.
     &              (lakon(nelem)(2:7).eq.'RELOLI')) then
               xflow_oil=prop(index+5)
               k_oil=int(prop(index+6))
            elseif((lakon(nelem)(2:5).eq.'RECO').or.
     &              (lakon(nelem)(2:7).eq.'REBEMA').or.
     &              (lakon(nelem)(2:7).eq.'REBEMI').or.
     &              (lakon(nelem)(2:8).eq.'REBEIDC')) then
               xflow_oil=prop(index+6)
               k_oil=int(prop(index+7))
            elseif(lakon(nelem)(2:7).eq.'REBEIDR') then
               xflow_oil=prop(index+8)
               k_oil=int(prop(index+9))
            endif
         endif
!     
         if(pt1.ge.pt2) then
            inv=1
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)+physcon(1)
            Tt2=v(0,node2)+physcon(1)     
            icase=0
            call ts_calc(xflow,Tt1,Pt1,kappa,r,a1,T1,icase)
            call ts_calc(xflow,Tt2,Pt2,kappa,r,a2,T2,icase)
!     
         else
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
            xflow=-v(1,nodem)*iaxial
            Tt1=v(0,node2)+physcon(1)
            Tt2=v(0,node1)+physcon(1)
            icase=0
            call ts_calc(xflow,Tt1,Pt1,kappa,r,a1,T1,icase)
            call ts_calc(xflow,Tt2,Pt2,kappa,r,a2,T2,icase)
!            
         endif
!     
!     calculation of the dynamic viscosity 
!     
         if( lakon(nelem)(2:3).eq.'RE') then
            icase=0
         elseif(lakon(nelem)(2:5).eq.'REEX') then
            if(lakon(int(prop(index+4)))(2:6).eq.'GAPFA') then
               icase=0
            elseif(lakon(int(prop(index+4)))(2:6).eq.'GAPFI') then
               icase=1
            endif
         endif
!     
         if (A1.le.A2) then
            if(dabs(dvi).lt.1E-30) then
               kgas=0
               call dynamic_viscosity(kgas,T1,dvi)
            endif
         else
            if(dabs(dvi).lt.1E-30) then
               kgas=0
               call dynamic_viscosity(kgas,T2,dvi)
            endif
         endif
!     
!     Reynolds number calculation
!     
         if (lakon(nelem)(2:5).eq.'REBR') then
            d=dsqrt(4d0*A1/Pi)
            reynolds=dabs(xflow)*d/(dvi*A1)
         else
            d=prop(index+3)    
            if(A1.le.A2) then
               reynolds=dabs(xflow)*d/(dvi*A1)
            else
               reynolds=dabs(xflow)*d/(dvi*A2)
            endif
         endif

         if(xflow_oil.lt.1E-10) then
            xflow_oil=0d0
         endif
!
!     BEND MILLER with oil
!  
         if(lakon(nelem)(2:7).eq.'REBEMI') then
            if(xflow_oil.ne.0d0) then
               call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &              xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &              v,dvi,cp,r,k_oil,phi,zeta,nshcon,nrhcon,
     &              shcon,rhcon,ntmat_,mi)

               call zeta_calc(nelem,prop,ielprop,lakon,reynolds,zeta,
     &              isothermal,kon,ipkon,R,Kappa,v,mi)
! 
               zeta_phi=phi*zeta
            else
!
               call zeta_calc(nelem,prop,ielprop,lakon,reynolds,zeta,
     &              isothermal,kon,ipkon,R,Kappa,v,mi)
               phi=1.d0
               zeta_phi=phi*zeta
!               
            endif
!
!   long  orifice in a wall with oil after Idelchik
!
         elseif(lakon(nelem)(2:7).eq.'RELOID') then
            if(xflow_oil.ne.0d0) then
               call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &              xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &              v,dvi,cp,r,k_oil,phi,zeta,nshcon,nrhcon,
     &              shcon,rhcon,ntmat_,mi)
!
               zeta_phi=phi*zeta
            else
!
               call zeta_calc(nelem,prop,ielprop,lakon,reynolds,zeta,
     &              isothermal,kon,ipkon,R,Kappa,v,mi)
               phi=1.d0
               zeta_phi=phi*zeta
            endif
!
!     every other zeta elements with/without oil
!
         else
!
            if(xflow_oil.ne.0) then
               call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &              xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &              v,dvi,cp,r,k_oil,phi,zeta,nshcon,nrhcon,
     &              shcon,rhcon,ntmat_,mi)
!
               call zeta_calc(nelem,prop,ielprop,lakon,reynolds,zeta,
     &              isothermal,kon,ipkon,R,Kappa,v,mi)
! 
               zeta_phi=phi*zeta
            else
               phi=1.d0
               call zeta_calc(nelem,prop,ielprop,lakon,reynolds,zeta,
     &              isothermal,kon,ipkon,R,Kappa,v,mi)
               zeta_phi=phi*zeta
            endif
         endif
!     
         if(zeta.le.0) then
            pt1=v(2,node1)
            pt2=v(2,node2)
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)
            Tt2=v(0,node2)
!     
         endif
!     
         if(.not.isothermal) then
            pt1pt2_crit=(0.5d0*kp1)**(zeta*kdkm1)
         else
            pt1pt2_crit=0.5d0*(3*kappa-1)**(zeta*kdkm1)
         endif
         pt1pt2=pt1/pt2
!     
!     Mach number calculation
!     
         M1=dsqrt(2d0/km1*(Tt1/T1-1d0))
         if((1.d0-M1).le.1E-3) then
  
            if(zeta.gt.0d0)then
               if(xflow_oil.eq.0) then
                  call limit_case_calc(a2,pt1,Tt2,xflow,zeta,r,kappa,
     &                 pt2_lim,M2)
               else
                  call limit_case_calc(a2,pt1,Tt2,xflow,zeta_phi,r,kappa
     &                 ,pt2_lim,M2)
               endif
            endif
         else
            M2=dsqrt(2d0/km1*(Tt2/T2-1d0))
         endif
!     
         if(iflag.eq.3) then
            write(1,*) ''
            write(1,55) 'In line ',int(nodem/1000),' from node ',node1,
     &           ' to node ', node2,' :   air massflow rate = '
     &           ,xflow,' kg/s'
     &           , ' , oil massflow rate = ',xflow_oil,' kg/s'
!            
            if(lakon(nelem)(4:5).ne.'BR') then
!     
!     for restrictors
!     
               if(inv.eq.1) then
                write(1,56)'      Inlet node ',node1,' :    Tt1=  ',Tt1,
     &              ' K , Ts1 = ',T1,' K , Pt1 = ',Pt1/1E5,
     &              ' Bar , M1 = ',M1
                write(1,*)'             Element F    ',set(numf)
     &              (1:30)
                write(1,57)'             eta = ',dvi,' kg/(m*s) , Re = '
     &           ,reynolds
                write(1,58)'             PHI = ',phi,' , ZETA = ',zeta,
     &        ' , ZETA_PHI = ',phi*zeta
                write(1,56)'       Outlet node ',node2,' :   Tt2 = ',
     &              Tt2,
     &              ' K, Ts2 = ',T2,' K , Pt2 = ',Pt2/1e5,
     &              ' Bar , M2= ',M2
!     
             else if(inv.eq.-1) then
                write(1,56)'       Inlet node ',node2,' :    Tt1= ',Tt1,
     &              ' K, Ts1= ',T1,' K, Pt1= ',Pt1/1E5,
     &              ' Bar, M1= ',M1
                write(1,*)'            Element F    ',set(numf)
     &              (1:30)
                write(1,57)'             eta =',dvi,' kg/(m*s) , Re = '
     &           ,reynolds
                write(1,58)'             PHI = ',phi,' , ZETA = ',zeta,
     &         ' , ZETA_PHI = ',phi*zeta
                write(1,56)'       Outlet node ',node1,' :   Tt2 = ',
     &              Tt2,
     &              ' K , Ts2 = ',T2,' K , Pt2 = ',Pt2/1e5,
     &              ' Bar , M2 = ',M2
             endif
          else
!     
!     for branches
!     
             if(inv.eq.1) then
                write(1,56)'       Inlet node ',node1,' :    Tt1 = ',
     &              Tt1,
     &              ' K , Ts1 = ',T1,' K , Pt1 = ',Pt1/1E5,
     &              ' Bar , M1 = ',M1
                write(1,*)'             Element B    ',set(numf)
     &              (1:20)
                write(1,57)'             eta = ',dvi,' kg/(m*s) , Re = '
     &           ,reynolds
                write(1,58)'             PHI = ',phi,' , ZETA = ',zeta,
     &        ' , ZETA_PHI = ',phi*zeta              
                write(1,56)'       Outlet node ',node2,' :   Tt2 = ',
     &              Tt2,
     &              ' K , Ts2 = ',T2,' K , Pt2 = ',Pt2/1E5,
     &              ' Bar , M2 = ',M2
!     
             else if(inv.eq.-1) then
                write(1,56)'       Inlet node ',node2,' :    Tt1 = ',
     &              Tt1,
     &              'K, Ts1 = ',T1,' K, Pt1 = ',Pt1/1E5,
     &              ' Bar , M1 = ',M1
                write(1,*)'             Element B    ',set(numf)
     &              (1:20)
                write(1,57)'             eta =',dvi,' kg/(m*s) , Re = '
     &           ,reynolds
                write(1,58)'             PHI = ',phi,' , ZETA = ',zeta,
     &' , ZETA_PHI = ',phi*zeta
                write(1,56)'       Outlet node ',node1,' :   Tt2 = ',
     &              Tt2,
     &              ' K , Ts2 = ',T2,' K , Pt2 = ',Pt2/1E5,
     &              ' Bar , M2 = ',M2
             endif
          endif
!     
       elseif (iflag.eq.4) then
!
!        Write the main information about the element
          write(1,*) ''
!         
          write(1,78)'Element nr.= ',nelem,', type=',lakon(nelem),
     &                 ', name= ',set(numf)(1:30)
!         
          write(1,79)'Nodes: ',node1,',',nodem,',',node2
!
          if(lakon(nelem)(4:5).ne.'BR') then
!     
!           for restrictors
!
             write(1,80)'Inlet: Tt1= ',Tt1,
     &              ', pt1= ',pt1,', M1= ',M1
!
             write(1,77)'mass flow = ',xflow,
     &              ', oil mass flow =',xflow_oil,
     &              ', kappa = ',kappa,
     &              ', phi= ',phi,
     &              ', zeta= ',zeta,
     &              ', eta= ',dvi,
     &              ', Re= ',reynolds,
     &              ', phi = ',phi,
     &              ', zeta_phi = ',zeta_phi

             write(1,80)'Outlet: Tt2= ',Tt2,
     &              ', pt2= ',pt2,', M2= ',M2
!
          else
!     
!     for branches
!     
             if(inv.eq.1) then
                write(1,56)'       Inlet node ',node1,':    Tt1= ',Tt1,
     &              ' K, Ts1= ',T1,' K, Pt1= ',Pt1/1E5,
     &              ' Bar, M1= ',M1
                write(1,*)'             Element B    ',set(numf)
     &              (1:20)
                write(1,57)'             eta= ',dvi,' kg/(m*s), Re= '
     &,reynolds,', PHI= ',phi,', ZETA= ',zeta
                write(1,56)'       Outlet node ',node2,':   Tt2= ',Tt2,
     &              ' K, Ts2= ',T2,' K, Pt2= ',Pt2/1E5,
     &              ' Bar, M2= ',M2
!     
             else if(inv.eq.-1) then
                write(1,56)'       Inlet node ',node2,':    Tt1= ',Tt1,
     &              'K, Ts1= ',T1,'K, Pt1= ',Pt1/1E5,
     &              'Bar, M1= ',M1
                write(1,*)'             Element B    ',set(numf)
     &              (1:20)
                write(1,57)'                 eta=',dvi,' kg/(m*s), Re= '
     &              ,reynolds,', PHI= ',phi,', ZETA= ',zeta
                write(1,56)'       Outlet node ',node1,':   Tt2= ',Tt2,
     &              ' K, Ts2= ',T2,' K, Pt2= ',Pt2/1E5,
     &              ' Bar, M2= ',M2
             endif
          endif
       endif
      endif
!     
 55   format(1x,a,i6.3,a,i6.3,a,i6.3,a,f9.6,a,a,f9.6,a) 
 56   format(1x,a,i6.3,a,f6.1,a,f6.1,a,f8.5,a,f8.6)
 57   format(1x,a,g9.4,a,g11.4)
 58   format(1x,a,f10.5,a,f10.5,a,f10.5)
 77   format(3x,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,
     &     e10.4,a,f10.2,a,f10.6,a,f10.6)
 78   format(a,i4,a,a,a,a)
 79   format(3x,a,i4,a,i4,a,i4)
 80   format(3x,a,f10.6,a,f10.2,a,f10.6)
!     
      xflow=xflow/iaxial
      df(3)=df(3)*iaxial
!     
      return
      end
