!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine calcgamma(nface,ielfa,vel,gradvfa,gamma,xlet,
     &  xxn,xxj,ipnei,betam,nef)
!
!     determine gamma:
!        upwind difference: gamma=0
!        central difference: gamma=1
!
      implicit none
!
      integer nface,ielfa(4,*),i,j,indexf,ipnei(*),iel1,iel2,nef
!
      real*8 vel(nef,0:5),gradvfa(3,3,*),xxn(3,*),xxj(3,*),vud,vcd,
     &  gamma(*),phic,xlet(*),betam
!
      do i=1,nface
         iel2=ielfa(2,i)
!
!        faces with only one neighbor need not be treated
!
         if(iel2.le.0) cycle
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
!
         vcd=(vel(iel2,1)-vel(iel1,1))*xxn(1,indexf)+
     &       (vel(iel2,2)-vel(iel1,2))*xxn(2,indexf)+
     &       (vel(iel2,3)-vel(iel1,3))*xxn(3,indexf)
!
         vud=2.d0*xlet(indexf)*
     &       (xxn(1,indexf)*(gradvfa(1,1,i)*xxj(1,indexf)+
     &                       gradvfa(1,2,i)*xxj(2,indexf)+
     &                       gradvfa(1,3,i)*xxj(3,indexf))+
     &        xxn(2,indexf)*(gradvfa(2,1,i)*xxj(1,indexf)+
     &                       gradvfa(2,2,i)*xxj(2,indexf)+
     &                       gradvfa(2,3,i)*xxj(3,indexf))+
     &        xxn(3,indexf)*(gradvfa(3,1,i)*xxj(1,indexf)+
     &                       gradvfa(3,2,i)*xxj(2,indexf)+
     &                       gradvfa(3,3,i)*xxj(3,indexf)))
c         write(*,*) xlet(indexf)
c         write(*,*) xxn(1,indexf),xxn(2,indexf),xxn(3,indexf)
c         write(*,*) xxj(1,indexf),xxj(2,indexf),xxj(3,indexf)
c         write(*,*) 'calcgamma ',vcd,vud
!
         if(dabs(vud).lt.1.d-20) then
            gamma(i)=1.d0
            cycle
         endif
!            
         phic=1.d0-vcd/vud
!
         if(phic.ge.1.d0) then
            gamma(i)=0.d0
         elseif(phic.le.0.d0) then
            gamma(i)=0.d0
         elseif(betam.le.phic) then
            gamma(i)=1.d0
         else
            gamma(i)=phic/betam
         endif
c         write(*,*) 'calcgamma ',i,gamma(i)
      enddo
!            
      return
      end
