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
      subroutine massflow_percent(node1,node2,nodem,nelem,lakon,kon,
     &        ipkon,nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,iaxial)
!     
!     partial massflow  element
!
!     author: Yannick Muller
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(4),idirf(4),index,iflag,
     &     inv,ipkon(*),kon(*),number,kgas,iaxial,
     &     nodea,nodeb,mi(*),i,itype,nodemup,
     &     nrhcon(*),ntmat_,nshcon(*)
!
      real*4 ofvidg
!     
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(4),kappa,R,a,d,xl,
     &     p1,p2,T1,physcon(3),pi,xflow_oil,T2,co(3,*),vold(0:mi(2),*),
     &     xflow_sum,percent_xflow,cp,dvi,pt1,pt2,Tt1,Tt2,
     &     shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*)
!
      external ofvidg
!
      pi=4.d0*datan(1.d0) 
  
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
         index=ielprop(nelem)
         percent_xflow=prop(index+1)
         xflow_sum=0
!         
         do i=2,10 
            if(int(prop(index+i)).ne.0) then
               nodemup=kon(ipkon(int(prop(index+i)))+2)
               if(v(1,nodemup).gt.0)then         
                  xflow_sum=xflow_sum+v(1,nodemup)*iaxial
               endif
            endif
         enddo
!
         if(xflow_sum.eq.0d0) then
            xflow_sum=0.001d0
         endif
!
         xflow=xflow_sum*percent_xflow
!     
      elseif ((iflag.eq.2).or.(iflag.eq.3))then
!     
        
         index=ielprop(nelem)
         percent_xflow=prop(index+1)
         xflow_sum=0
         do i=2,10
            nodemup=kon(ipkon(int(prop(index+i)))+2)
            if(int(prop(index+i)).ne.0) then
               nodemup=kon(ipkon(int(prop(index+i)))+2)
               if(v(1,nodemup).gt.0)then        
                  xflow_sum=xflow_sum+v(1,nodemup)*iaxial
               endif
            endif
         enddo
         
         if(xflow_sum.eq.0.d0) then
            xflow_sum=1E-5
         endif
!
         inv=1
!
         pt1=v(2,node1)
         pt2=v(2,node2)
         xflow=v(1,nodem)*iaxial
         Tt1=v(0,node1)+physcon(1)
         Tt2=v(0,node2)+physcon(1)
!
         nodef(1)=node1
         nodef(2)=node1
         nodef(3)=nodem
         nodef(4)=node2
!     
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!
         if(iflag.eq.2) then
            numf=4
!
            f=xflow/xflow_sum-percent_xflow
!     
            df(1)=0
            df(2)=0
            df(3)=1/xflow_sum
            df(4)=0
!     
!     output
!     
         elseif(iflag.eq.3)then
!         
            xflow_oil=0
!
            write(1,*) ''
            write(1,55) 'In line ',int(nodem/1000),' from node ',node1,
     &           ' to node ', node2,' :   air massflow rate = '
     &           ,inv*xflow,' kg/s',
     &           ', oil massflow rate = ',xflow_oil,' kg/s'
 55         format(1X,A,I6.3,A,I6.3,A,I6.3,A,F9.6,A,A,F9.6,A)
!            
            write(1,56)'       Inlet node ',node1,' :   Tt1 = ',Tt1,
     &           ' K , Ts1 = ',Tt1,' K , Pt1 = ',Pt1/1E5, ' Bar'
!            
            write(1,*)'             Element PMASS ',set(numf)(1:30)
            write(1,57)'        Massflow upstream = ',xflow_sum,
     &        ' [kg/s]'
            write(1,58)'        Massflow fraction = ', percent_xflow
            write(1,56)'       Outlet node ',node2,':    Tt2=',Tt2,
     &           'K, Ts2=',Tt2,'K, Pt2=',Pt2/1e5, 'Bar'
!            
         endif
      endif
!     
 56   format(1X,A,I6.3,A,f6.1,A,f6.1,A,f9.5,A)
 57   format(1X,A,G9.4,A)
 58   format(1X,A,G9.4)
!     
      xflow=xflow/iaxial
      df(3)=df(3)*iaxial
!     
      return
      end
