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
      subroutine physicalconstantss(inpc,textpart,physcon,
     &  istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc)
!
!     reading the input deck: *PHYSICAL CONSTANTS
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer i,istep,istat,n,key,iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  ipoinpc(0:*)
!
      real*8 physcon(*)
!
      if(istep.gt.0) then
         write(*,*) '*ERROR in physicalconstants: *PHYSICAL CONSTANTS'
         write(*,*) '        should only be used before the first STEP'
         call exit(201)
      endif
!
      do i=2,n
         if(textpart(i)(1:13).eq.'ABSOLUTEZERO=') then
            read(textpart(i)(14:33),'(f20.0)',iostat=istat) physcon(1)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*PHYSICAL CONSTANTS%")
         elseif(textpart(i)(1:16).eq.'STEFANBOLTZMANN=') then
            read(textpart(i)(17:36),'(f20.0)',iostat=istat) physcon(2)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*PHYSICAL CONSTANTS%")
         elseif(textpart(i)(1:14).eq.'NEWTONGRAVITY=') then
            read(textpart(i)(15:24),'(f20.0)',iostat=istat) physcon(3)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*PHYSICAL CONSTANTS%")
         else
            write(*,*) 
     &        '*WARNING in physicalconstants: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*PHYSICAL CONSTANTS%")
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end







