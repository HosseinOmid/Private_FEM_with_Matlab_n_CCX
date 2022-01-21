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
      subroutine correctvel(hel,adv,vfa,ipnei,area,bv,xxn,neifa,
     &  lakon,ne,neq)
!
!     correction of the velocity at the element centers due to the
!     pressure change (balance of mass)
!
!     the solution is stored in field bv.
!
      implicit none
!
      character*8 lakon(*)
!
      integer i,j,k,ne,jdof1,indexf,ipnei(*),neifa(*),ifa,
     &  neq,numfaces
!
      real*8 bv(neq,3),hel(3,*),adv(*),xxn(3,*),area(*),vfa(0:5,*)
!
      do i=1,ne
!
         jdof1=i
         do k=1,3
            bv(jdof1,k)=0.d0
         enddo
         indexf=ipnei(i)
!
         if(lakon(i)(4:4).eq.'8') then
            numfaces=6
         elseif(lakon(i)(4:4).eq.'6') then
            numfaces=5
         else
            numfaces=4
         endif
         do j=1,numfaces
            indexf=indexf+1
            ifa=neifa(indexf)
            do k=1,3
               bv(jdof1,k)=bv(jdof1,k)
     &              +vfa(4,ifa)*area(ifa)*xxn(k,indexf)
            enddo
         enddo
c         write(*,*) 'correctvel1 ',i,bv(jdof1,1),bv(jdof1,2),bv(jdof1,3)
!
         do k=1,3
            bv(jdof1,k)=(hel(k,jdof1)-bv(jdof1,k))/adv(jdof1)
         enddo
c         write(*,*) 'correctvel2 ',i,bv(jdof1,1),bv(jdof1,2),bv(jdof1,3)
      enddo
!  
      return
      end
