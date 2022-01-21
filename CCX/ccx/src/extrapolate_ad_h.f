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
      subroutine extrapolate_ad_h(nface,ielfa,xrlfa,adv,advfa,
     &                  hel,hfa)
!
!     inter/extrapolation of adv at the center of the elements
!     to the center of the faces
!
!     inter/extrapolation of h at the center of the elements 
!     to the center of the faces; division through advfa to obtain
!     the face velocity
!
      implicit none
!
      integer nface,ielfa(4,*),ipo1,iel2,ipo3,i,j
!
      real*8 xrlfa(3,*),xl1,advfa(*),adv(*),hel(3,*),hfa(3,*)
!     
      do i=1,nface
         ipo1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!
!           internal face
!
            advfa(i)=xl1*adv(ipo1)+xrlfa(2,i)*adv(iel2)
            do j=1,3
               hfa(j,i)=(xl1*hel(j,ipo1)
     &                   +xrlfa(2,i)*hel(j,iel2))/advfa(i)
            enddo
         elseif(ielfa(3,i).gt.0) then
!
!           external face; linear extrapolation
!
            ipo3=ielfa(3,i)
            advfa(i)=xl1*adv(ipo1)+xrlfa(3,i)*adv(ipo3)
            do j=1,3
               hfa(j,i)=(xl1*hel(j,ipo1)+xrlfa(3,i)*hel(j,ipo3))
     &                  /advfa(i)
            enddo
         else
!
!           external face: constant extrapolation (only one adjacent
!           element layer)
!
            advfa(i)=adv(ipo1)
            do j=1,3
               hfa(j,i)=hel(j,ipo1)/advfa(i)
            enddo
         endif
c         write(*,*) 'extrapolate_ad ',ielfa(1,i),ielfa(4,i),advfa(i)
      enddo
!            
      return
      end
