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
      subroutine extrapolate_gradtel(nface,ielfa,xrlfa,gradtel,gradtfa)
!
!     inter/extrapolation of p element values to the faces
!
      implicit none
!
      integer nface,ielfa(4,*),i,iel1,iel2,l
!
      real*8 xrlfa(3,*),gradtel(3,*),gradtfa(3,*),xl1
!     
      do i=1,nface
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!
!           face in between two elements
!
            do l=1,3
               gradtfa(l,i)=xl1*gradtel(l,iel1)+
     &              xrlfa(2,i)*gradtel(l,iel2)
            enddo
         elseif(ielfa(3,i).gt.0) then
!     
!     boundary face; linear extrapolation
!     
            do l=1,3
               gradtfa(l,i)=xl1*gradtel(l,iel1)+
     &              xrlfa(3,i)*gradtel(l,ielfa(3,i))
            enddo
         else
!     
!     boundary face; constant extrapolation (one element layer)
!     
            do l=1,3
               gradtfa(l,i)=gradtel(l,iel1)
            enddo
         endif
c         write(*,*) 'extrapolate_gradtel ',(gradtfa(l,i),l=1,3)
      enddo
!            
      return
      end
