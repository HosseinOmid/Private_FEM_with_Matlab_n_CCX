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
      subroutine forcadd(node,i,val,nodeforc,ndirforc,xforc,
     &  nforc,nforc_,iamforc,iamplitude,nam,ntrans,trab,inotr,co,
     &  ikforc,ilforc,isector,add,user,idefforc)
!
!     adds a cload condition to the data base
!
      implicit none
!
      logical add,user
!
      integer nodeforc(2,*),ndirforc(*),node,i,nforc,nforc_,j,
     &  iamforc(*),iamplitude,nam,ntrans,inotr(2,*),itr,idf(3),
     &  ikforc(*),ilforc(*),idof,id,k,isector,idefforc(*)
!
      real*8 xforc(*),val,trab(7,*),a(3,3),co(3,*)
!
      if(ntrans.eq.0) then
         itr=0
      else
         itr=inotr(1,node)
      endif
!
      if(itr.eq.0) then
!
!        no transformation applies to the node
!
         idof=8*(node-1)+i
         call nident(ikforc,idof,nforc,id)
         if(id.gt.0) then
            do
               if(ikforc(id).eq.idof) then
                  k=ilforc(id)
                  if(nodeforc(2,k).eq.isector) then
                     if(add.or.(idefforc(k).eq.1)) then
                        xforc(k)=xforc(k)+val
                     else
                        xforc(k)=val
                        if(.not.user) idefforc(k)=1
                     endif
                     if(nam.gt.0) iamforc(k)=iamplitude
                     return
                  endif
                  id=id-1
                  if(id.eq.0) exit
               else
                  exit
               endif
            enddo
         endif
c
         nforc=nforc+1
         if(nforc.gt.nforc_) then
            write(*,*) '*ERROR in forcadd: increase nforc_'
            call exit(201)
         endif
         nodeforc(1,nforc)=node
         nodeforc(2,nforc)=isector
         ndirforc(nforc)=i
         xforc(nforc)=val
         if(.not.user) idefforc(nforc)=1
         if(nam.gt.0) iamforc(nforc)=iamplitude
!
!        updating ikforc and ilforc
!            
         do j=nforc,id+2,-1
            ikforc(j)=ikforc(j-1)
            ilforc(j)=ilforc(j-1)
         enddo
         ikforc(id+1)=idof
         ilforc(id+1)=nforc
      else
!
!        a transformation applies
!
         call transformatrix(trab(1,itr),co(1,node),a)
!
         do j=1,3
            idf(j)=0
            idof=8*(node-1)+j
            call nident(ikforc,idof,nforc,id)
            if(id.gt.0) then
               do
                  if(ikforc(id).eq.idof) then
                     k=ilforc(id)
                     if(nodeforc(2,k).eq.isector) then
                        idf(j)=ilforc(id)
                        exit
                     endif
                     id=id-1
                     if(id.eq.0) exit
                  else
                     exit
                  endif
               enddo
            endif
         enddo
!
         if((idf(1).ne.0).and.(.not.user)) then
!
!        a force was previously applied to this node. The component
!        in direction i is filtered out and replaced by the new
!        value
!
!        if an amplitude is selected, it applies to all components
!        of the force in the node. No separate amplitudes are allowed.
!
            if((.not.add).and.(idefforc(idf(i)).ne.1))
     &        val=val-xforc(idf(1))*a(1,i)-xforc(idf(2))*a(2,i)
     &             -xforc(idf(3))*a(3,i)
!
            xforc(idf(1))=xforc(idf(1))+val*a(1,i)
            xforc(idf(2))=xforc(idf(2))+val*a(2,i)
            xforc(idf(3))=xforc(idf(3))+val*a(3,i)
!
!           only first entry is tagged
!
            idefforc(idf(i))=1
!
            if(nam.gt.0) then
               do j=1,3
                  iamforc(idf(j))=iamplitude
               enddo
            endif
         else
            do j=1,3
               nforc=nforc+1
               if(nforc.gt.nforc_) then
                  write(*,*) '*ERROR in forcadd: increase nforc_'
                  call exit(201)
               endif
               nodeforc(1,nforc)=node
               nodeforc(2,nforc)=isector
               ndirforc(nforc)=j
               if(user) then
                  xforc(nforc)=val
               else
                  xforc(nforc)=val*a(j,i)
               endif
               idefforc(nforc)=1
               if(nam.gt.0) iamforc(nforc)=iamplitude
!
!              updating ikforc and ilforc
! 
               idof=8*(node-1)+j
               call nident(ikforc,idof,nforc-1,id)
               do k=nforc,id+2,-1
                  ikforc(k)=ikforc(k-1)
                  ilforc(k)=ilforc(k-1)
               enddo
               ikforc(id+1)=idof
               ilforc(id+1)=nforc
            enddo
         endif
      endif
!
      return
      end

