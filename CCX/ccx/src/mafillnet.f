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
!     This subroutine creates the matrix ac for gas problems
!     
      subroutine mafillnet(itg,ieg,ntg,ac,nload,sideload,nelemload,
     &     xloadact,lakon,ntmat_,v,shcon,nshcon,ipkon,kon,co,nflow,iinc,
     &     istep,dtime,ttime,time,
     &     ielmat,nteq,prop,ielprop,nactdog,nacteq,physcon,
     &     rhcon,nrhcon,ipobody,ibody,xbodyact,nbody,vold,xloadold,
     &     reltime,nmethod,set,mi,nmpc,nodempc,ipompc,coefmpc,labmpc,
     &     iaxial)
!     
      implicit none
!     
      logical identity
      character*8 lakonl,lakon(*)
      character*20 sideload(*),labmpc(*)
      character*81 set(*)
!     
      integer mi(*),itg(*),ieg(*),ntg,nteq,nflow,nload,
     &     ielmat(mi(3),*),iaxial,
     &     nelemload(2,*),nope,nopes,mint2d,i,j,k,l,iflag,
     &     node,imat,ntmat_,id,ifaceq(9,6),ifacet(7,4),
     &     ifacew(8,5),node1,node2,nshcon(*),nelem,ig,index,konl(20),
     &     ipkon(*),kon(*),idof,iinc,ibody(3,*),istep,jltyp,nfield,
     &     ipobody(2,*),nodem,ieq,kflag,nrhcon(*),numf,
     &     idofp1,idofp2,idofm,idoft1,idoft2,idoft,nactdog(0:3,*),
     &     nacteq(0:3,*),ielprop(*),nodef(8),idirf(8),nbody,
     &     nmethod,icase,nmpc,nodempc(3,*),ipompc(*),idir,ider
!     
      real*8 ac(nteq,*),xloadact(2,*),cp,h(2),physcon(*),dvi,
     &     xl2(3,8),coords(3),dxsj2,temp,xi,et,weight,xsj2(3),
     &     gastemp,v(0:mi(2),*),shcon(0:3,ntmat_,*),co(3,*),shp2(7,8),
     &     ftot,field,prop(*),f,df(8),tg1,tg2,r,rho,tl2(8),
     &     dtime,ttime,time,areaj,xflow,tvar(2),g(3),coefmpc(*),
     &     rhcon(0:1,ntmat_,*),xbodyact(7,*),sinktemp,ts1,ts2,xs2(3,7),
     &     xdenom1,xdenom2,xcst,xk1,xk2,expon,a,dt1,dt2,kappa,
     &     pt1,pt2,inv,vold(0:mi(2),*),xloadold(2,*),reltime,pi
!     
      include "gauss.f"
!     
      data ifaceq /4,3,2,1,11,10,9,12,21,
     &            5,6,7,8,13,14,15,16,22,
     &            1,2,6,5,9,18,13,17,23,
     &            2,3,7,6,10,19,14,18,24,
     &            3,4,8,7,11,20,15,19,25,
     &            4,1,5,8,12,17,16,20,26/
      data ifacet /1,3,2,7,6,5,11,
     &     1,2,4,5,9,8,12,
     &     2,3,4,6,10,9,13,
     &     1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/
!     
      data iflag /2/
!     
      kflag=2
      ider=1
!
      Pi=4.d0*datan(1.d0)
      tvar(1)=time
      tvar(2)=ttime+time
!     
!     reinitialisation of the Ac matrix
!
      do i=1,nteq
         do j=1,nteq
            ac(i,j)=0.d0
         enddo
      enddo
!
!     solving for the gas temperatures in forced convection
!     
      ftot=0.d0
!     
!     element contribution.
!              
 
      do i=1,nflow
         nelem=ieg(i)
         index=ipkon(nelem)
         node1=kon(index+1)
         nodem=kon(index+2)
         node2=kon(index+3)
!
         xflow=v(1,nodem)
!
        if((lakon(nelem)(2:3).ne.'LP').and.
     &     (lakon(nelem)(2:3).ne.'LI')) then
            if(node1.eq.0) then
               tg1=v(0,node2)
               tg2=tg1
               ts1=v(3,node2)
               ts2=ts1
            elseif(node2.eq.0) then
               tg1=v(0,node1)
               tg2=tg1
               ts1=v(3,node1)
               ts2=ts1
            else
               tg1=v(0,node1)
               tg2=v(0,node2)
               ts1=v(3,node1)
               ts2=v(3,node2)
            endif
!
            gastemp=(ts1+ts2)/2.d0
!
!     for liquid pipe element only the upstream temperature is used to 
!     determine thematerial properties
!
         else
            
            if(xflow.gt.0) then
               if(node1.eq.0) then
                  gastemp=v(0,node2)
               else
                  gastemp=v(0,node1)
               endif
            else
               if(node2.eq.0) then
                  gastemp=v(0,node1)
               else
                  gastemp=v(0,node2)
               endif
            endif
!
            if(node1.eq.0) then
               tg2=v(0,node2)
               tg1=tg2
            elseif(node2.eq.0) then
               tg1=v(0,node1)
               tg2=tg1
            else
               tg1=v(0,node1)
               tg2=v(0,node2)
            endif
         endif
!
         imat=ielmat(1,nelem)
!
         call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,cp,r,dvi,
     &        rhcon,nrhcon,rho)
!     
         kappa=(cp/(cp-R))
!
!     Definitions of the constant for isothermal flow elements
!
         if((lakon(nelem)(2:6).eq.'GAPFI')
     &        .or.(lakon(nelem)(2:6).eq.'GAPII'))then
            if((node1.ne.0).and.(node2.ne.0)) then
!
               icase=1
               A=prop(ielprop(nelem)+1)
               pt1=v(2,node1)
               pt2=v(2,node2)
               if(pt1.ge.pt2)then
                  inv=1.d0
                  pt1=v(2,node1)
                  pt2=v(2,node2)
                  if(dabs(tg2/ts2-(1+0.5*(kappa-1)/kappa)).lt.1E-5) then
                    
                     pt2=dabs(xflow)*dsqrt(Tg2*R)/A
     &                    *(1+0.5*(kappa-1)/kappa)
     &                    **(0.5*(kappa+1)/(kappa-1)) 
                  endif
                  tg1=v(0,node1)
                  call ts_calc(xflow,Tg1,Pt1,kappa,r,a,Ts1,icase)
!
                  tg2=v(0,node2)
                  call ts_calc(xflow,Tg2,Pt2,kappa,r,a,Ts2,icase)
               else
!                 
                  inv=-1.d0
                  pt1=v(2,node2)
                  pt2=v(2,node1)
                  if(dabs(tg2/ts2-(1+0.5*(kappa-1)/kappa)).lt.1E-5) then
                    
                     pt2=dabs(xflow)*dsqrt(Tg2*R)/A
     &                    *(1+0.5*(kappa-1)/kappa)
     &                    **(0.5*(kappa+1)/(kappa-1)) 
                  endif
                  tg1=v(0,node2)
                  call ts_calc(xflow,Tg1,Pt1,kappa,r,a,Ts1,icase)
                  tg2=v(0,node1)
                  call ts_calc(xflow,Tg2,Pt2,kappa,r,a,Ts2,icase)
               endif
               dt1=tg1/ts1-1d0
               dt2=tg2/ts2-1d0
               expon=2.d0*kappa/(kappa-1.d0)
               xcst=2.d0*Cp*A**2/r**2
               xk1=pt1**2*(ts1/tg1)**expon
               xdenom1=xcst*xk1*(1.d0-expon*(tg1/ts1-1.d0))
     &              /ts1+2.d0*xflow**2
               xk2=pt2**2*(ts2/tg2)**expon
               xdenom2=xcst*xk2*(1.d0-expon*(tg2/ts2-1.d0))
     &              /ts2+2.d0*xflow**2
            endif
         endif
!     
         if(node1.ne.0) then
            idoft1=nactdog(0,node1)
            idofp1=nactdog(2,node1)
         else
            idoft1=0
            idofp1=0
         endif
         if(node2.ne.0) then
            idoft2=nactdog(0,node2)
            idofp2=nactdog(2,node2)
         else
            idoft2=0
            idofp2=0
         endif
         idofm=nactdog(1,nodem)
!     
         if(node1.ne.0) then
!     
!     energy equation contribution node1
!     
            if (nacteq(0,node1).ne.0) then
               ieq=nacteq(0,node1)
               if ((xflow.le.0d0).and.(nacteq(3,node1).eq.0))then
!     
!     adiabatic element
!     
                  if(idoft1.ne.0) then
                     ac(ieq,idoft1)=ac(ieq,idoft1)-cp*xflow
                  endif
!     
                  if(idoft2.ne.0)then
                     ac(ieq,idoft2)=ac(ieq,idoft2)+cp*xflow
                  endif
!     
                  if(idofm.ne.0) then
                     ac(ieq,idofm)=ac(ieq,idofm)-cp*(tg1-tg2)
                  endif
!
               elseif(nacteq(3,node1).ne.0)then
!     
!     isothermal element
!     
                  if(nacteq(3,node1).eq.node2) then
!     
                     if(inv.eq.-1d0) then
                        if(idoft1.ne.0) then
                           ac(ieq,idoft1)=-xcst*xk1*(1.d0-expon
     &                          *(1.d0-ts1/tg1))/(xdenom1*ts1) 
                        endif
!     
                        if(idoft2.ne.0)then
                           ac(ieq,idoft2)=xcst*xk2*(1.d0-expon
     &                          *(1.d0-ts2/tg2))/(xdenom2*ts2) 
                        endif
!     
                        if(idofm.ne.0) then
                           ac(ieq,idofm)=(-2.d0*xflow*ts2
     &                          /xdenom2+2.d0*xflow*ts1/xdenom1)
                        endif
!     
                        if(idofp1.ne.0) then
                           ieq=nacteq(2,idofp1)
                           ac(ieq,idofp1)=2.d0*xcst*dt1*xk1
     &                          /(pt1*xdenom1)
                        endif
!     
                        if(idofp2.ne.0) then
                           ac(ieq,idofp2)=-2.d0*xcst*dt2*xk2
     &                          /(pt2*xdenom2)
                        endif
!     
                     elseif(inv.eq.1d0)then
                        if(idoft1.ne.0) then
                           ac(ieq,idoft1)=xcst*xk1*(1.d0-expon
     &                          *(1.d0-ts1/tg1))/(xdenom1*ts1) 
                        endif
!     
                        if(idoft2.ne.0)then
                           ac(ieq,idoft2)=-xcst*xk2*(1.d0-expon
     &                          *(1.d0-ts2/tg2))/(xdenom2*ts2) 
                        endif
!     
                        if(idofm.ne.0) then
                           ac(ieq,idofm)=-(-2.d0*xflow*ts2
     &                          /xdenom2+2.d0*xflow*ts1/xdenom1)
                        endif
!     
                        if(idofp1.ne.0) then
                           ieq=nacteq(2,idofp1)
                           ac(ieq,idofp1)=-2.d0*xcst*dt1*xk1
     &                          /(pt1*xdenom1)
                        endif
!     
                        if(idofp2.ne.0) then
                           ac(ieq,idofp2)=2.d0*xcst*dt2*xk2
     &                          /(pt2*xdenom2)
                        endif
!     
                     endif
                  endif
               endif
            endif
!     
!     mass equation contribution node1
!     
               if (nacteq(1,node1).ne.0) then
                  ieq=nacteq(1,node1)
                  if (idofm.ne.0) then
                  ac(ieq,idofm)=1.d0
               endif
            endif
         endif
!     
         if(node2.ne.0) then
!     
!     energy equation contribution node2
!     
            if (nacteq(0,node2).ne.0) then
               ieq=nacteq(0,node2)
               if ((xflow.ge.0d0).and.(nacteq(3,node2).eq.0))then
!
!     adiabatic element
!     
                  if(idoft1.ne.0)then
                     ac(ieq,idoft1)=ac(ieq,idoft1)-cp*xflow
                  endif
!     
                  if(idoft2.ne.0) then
                     ac(ieq,idoft2)=ac(ieq,idoft2)+cp*xflow
                  endif
!     
                  if(idofm.ne.0) then
                     ac(ieq,idofm)=ac(ieq,idofm)+cp*(tg2-tg1)
                  endif
!     
               elseif((nacteq(3,node2).eq.node1))then
!     
!     isothermal element
!     
                  if(inv.eq.-1d0) then
                     if(idoft1.ne.0)then
                        ac(ieq,idoft1)=-xcst*xk1*(1.d0-expon
     &                       *(1.d0-ts1/tg1))/(xdenom1*ts1) 
                     endif
!     
                     if(idoft2.ne.0) then
                        ac(ieq,idoft2)=(xcst*xk2*(1.d0-expon
     &                       *(1.d0-ts2/tg2))/(xdenom2*ts2)) 
                     endif
!     
                     if(idofm.ne.0) then
                        ac(ieq,idofm)=(-2.d0*xflow*ts2
     &                       /xdenom2+2.d0*xflow*ts1/xdenom1)
                     endif
!     
                     if(idofp1.ne.0) then
                        ac(ieq,idofp1)=+2.d0*xcst*dt1*xk1
     &                       /(pt1*xdenom1)
                     endif
!     
                     if(idofp2.ne.0) then
                        ac(ieq,idofp2)=-2.d0*xcst*dt2*xk2
     &                       /(pt2*xdenom2)
                     endif
!
                  elseif(inv.eq.1d0) then
     
                     if(idoft1.ne.0)then
                        ac(ieq,idoft1)=xcst*xk1*(1.d0-expon
     &                       *(1.d0-ts1/tg1))/(xdenom1*ts1) 
                     endif
!     
                     if(idoft2.ne.0) then
                        ac(ieq,idoft2)=-(xcst*xk2*(1.d0-expon
     &                       *(1.d0-ts2/tg2))/(xdenom2*ts2)) 
                     endif
!     
                     if(idofm.ne.0) then
                        ac(ieq,idofm)=-(-2.d0*xflow*ts2
     &                       /xdenom2+2.d0*xflow*ts1/xdenom1)
                     endif
!     
                     if(idofp1.ne.0) then
                        ac(ieq,idofp1)=+2.d0*xcst*dt1*xk1
     &                       /(pt1*xdenom1)
                     endif
!     
                     if(idofp2.ne.0) then
                        ac(ieq,idofp2)=-2.d0*xcst*dt2*xk2
     &                       /(pt2*xdenom2)
                    endif
                 endif
!     
              endif 
           endif
!     
!          mass equation contribution node2
!     
           if (nacteq(1,node2).ne.0) then
              ieq=nacteq(1,node2)
              if(idofm.ne.0)then
                 ac(ieq,idofm)=-1.d0
              endif
           endif
        endif
!     
!     element equation
!     
        if (nacteq(2,nodem).ne.0) then
           ieq=nacteq(2,nodem)
!     
!     for liquids: determine the gravity vector
!     
           if(lakon(nelem)(2:3).eq.'LI') then
              do j=1,3
                 g(j)=0.d0
               enddo
               if(nbody.gt.0) then
                  index=nelem
                  do
                     j=ipobody(1,index)
                     if(j.eq.0) exit
                     if(ibody(1,j).eq.2) then
                        g(1)=g(1)+xbodyact(1,j)*xbodyact(2,j)
                        g(2)=g(2)+xbodyact(1,j)*xbodyact(3,j)
                        g(3)=g(3)+xbodyact(1,j)*xbodyact(4,j)
                     endif
                     index=ipobody(2,index)
                     if(index.eq.0) exit
                  enddo
               endif
            endif
!     
            call flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &           nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &           nodef,idirf,df,cp,R,rho,physcon,g,co,dvi,numf,
     &           vold,set,shcon,nshcon,rhcon,nrhcon,ntmat_,mi,ider,
     &           iaxial)
!
            do k=1,numf
               idof=nactdog(idirf(k),nodef(k))
               if(idof.ne.0)then
                  ac(ieq,idof)=df(k)
               endif
            enddo
         endif
      enddo
!     
!     convection with the walls
!     
      do i=1,nload
         if(sideload(i)(3:4).eq.'FC') then
            nelem=nelemload(1,i)
            index=ipkon(nelem)
            if(index.lt.0) cycle
            lakonl=lakon(nelem)
            node=nelemload(2,i)
            ieq=nacteq(0,node)
            if(ieq.eq.0) then 
               cycle
            endif
!     
            call nident(itg,node,ntg,id)
!     
!     calculate the area
!     
            read(sideload(i)(2:2),'(i1)') ig
!     
!     number of nodes and integration points in the face
!     
            if(lakonl(4:4).eq.'2') then
               nope=20
               nopes=8
            elseif(lakonl(4:4).eq.'8') then
               nope=8
               nopes=4
            elseif(lakonl(4:5).eq.'10') then
               nope=10
               nopes=6
            elseif(lakonl(4:4).eq.'4') then
               nope=4
               nopes=3
            elseif(lakonl(4:5).eq.'15') then
               nope=15
            elseif(lakonl(4:4).eq.'6') then
               nope=6
            else
               cycle
            endif
!     
            if(lakonl(4:5).eq.'8R') then
               mint2d=1
            elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R'))
     &              then
               if(lakonl(7:7).eq.'A') then
                  mint2d=2
               else
                  mint2d=4
               endif
            elseif(lakonl(4:4).eq.'2') then
               mint2d=9
            elseif(lakonl(4:5).eq.'10') then
               mint2d=3
            elseif(lakonl(4:4).eq.'4') then
               mint2d=1
            endif
!     
            if(lakonl(4:4).eq.'6') then
               mint2d=1
               if(ig.le.2) then
                  nopes=3
               else
                  nopes=4
               endif
            endif
            if(lakonl(4:5).eq.'15') then
               if(ig.le.2) then
                  mint2d=3
                  nopes=6
               else
                  mint2d=4
                  nopes=8
               endif
            endif
!     
!     connectivity of the element
!     
            index=ipkon(nelem)
            if(index.lt.0) then
               write(*,*) '*ERROR in mafillnet: element ',nelem
               write(*,*) '       is not defined'
               call exit(201)
            endif
            do k=1,nope
               konl(k)=kon(index+k)
            enddo
!     
!     coordinates of the nodes belonging to the face
!     
            if((nope.eq.20).or.(nope.eq.8)) then
               do k=1,nopes
                  tl2(k)=v(0,konl(ifaceq(k,ig)))
                  do j=1,3
                     xl2(j,k)=co(j,konl(ifaceq(k,ig)))+
     &                    v(j,konl(ifaceq(k,ig)))
                  enddo
               enddo
            elseif((nope.eq.10).or.(nope.eq.4)) then
               do k=1,nopes
                  tl2(k)=v(0,konl(ifacet(k,ig)))
                  do j=1,3
                     xl2(j,k)=co(j,konl(ifacet(k,ig)))+
     &                    v(j,konl(ifacet(k,ig)))
                  enddo
               enddo
            else
               do k=1,nopes
                  tl2(k)=v(0,konl(ifacew(k,ig)))
                  do j=1,3
                     xl2(j,k)=co(j,konl(ifacew(k,ig)))+
     &                    v(j,konl(ifacew(k,ig)))
                  enddo
               enddo
            endif
!     
!     integration to obtain the area and the mean
!     temperature
!     
            do l=1,mint2d
               if((lakonl(4:5).eq.'8R').or.
     &              ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
                  xi=gauss2d1(1,l)
                  et=gauss2d1(2,l)
                  weight=weight2d1(l)
               elseif((lakonl(4:4).eq.'8').or.
     &                 (lakonl(4:6).eq.'20R').or.
     &                 ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
                  xi=gauss2d2(1,l)
                  et=gauss2d2(2,l)
                  weight=weight2d2(l)
               elseif(lakonl(4:4).eq.'2') then
                  xi=gauss2d3(1,l)
                  et=gauss2d3(2,l)
                  weight=weight2d3(l)
               elseif((lakonl(4:5).eq.'10').or.
     &                 ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
                  xi=gauss2d5(1,l)
                  et=gauss2d5(2,l)
                  weight=weight2d5(l)
               elseif((lakonl(4:4).eq.'4').or.
     &                 ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
                  xi=gauss2d4(1,l)
                  et=gauss2d4(2,l)
                  weight=weight2d4(l)
               endif
!     
               if(nopes.eq.8) then
                  call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.4) then
                  call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
               elseif(nopes.eq.6) then
                  call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
               else
                  call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
               endif
!     
               dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &              xsj2(3)*xsj2(3))
               areaj=dxsj2*weight
!     
               temp=0.d0
               do k=1,3
                  coords(k)=0.d0
               enddo
               do j=1,nopes
                  temp=temp+tl2(j)*shp2(4,j)
                  do k=1,3
                     coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
                  enddo
               enddo
!     
               if(sideload(i)(5:6).ne.'NU') then
                  h(1)=xloadact(1,i)
               else
                  read(sideload(i)(2:2),'(i1)') jltyp
                  jltyp=jltyp+10
                  sinktemp=v(0,node)
                  call film(h,sinktemp,temp,istep,
     &                 iinc,tvar,nelem,l,coords,jltyp,field,nfield,
     &                 sideload(i),node,areaj,v,mi)
                  if(nmethod.eq.1) h(1)=xloadold(1,i)+
     &                 (h(1)-xloadold(1,i))*reltime
               endif
!     
               idoft=nactdog(0,node)
               if(idoft.gt.0) then
                  if(lakonl(5:7).eq.'0RA') then
                     ac(ieq,idoft)=ac(ieq,idoft)+2.d0*h(1)*dxsj2*weight
                  else
                     ac(ieq,idoft)=ac(ieq,idoft)+h(1)*dxsj2*weight
                  endif
               endif
            enddo
         endif
      enddo
!
!     additional multiple point constraints
!
      j=nteq+1
      do i=nmpc,1,-1
         if(labmpc(i)(1:7).ne.'NETWORK') cycle
         j=j-1
         index=ipompc(i)
!
         do
            node=nodempc(1,index)
            idir=nodempc(2,index)
            if(nactdog(idir,node).ne.0) then
               ac(j,nactdog(idir,node))=coefmpc(index)
            endif
            index=nodempc(3,index)
            if(index.eq.0) exit
         enddo
      enddo
!
!      write(30,*) nteq
!      do i=1,nteq
!         write(30,'(17(1x,e11.4))') (ac(i,j),j=1,nteq)
!      enddo
!
      return
      end
      
      
