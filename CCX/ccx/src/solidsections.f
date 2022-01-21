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
      subroutine solidsections(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,ielmat,matname,nmat,ielorien,orname,norien,
     &  lakon,thicke,kon,ipkon,irstrt,istep,istat,n,iline,ipol,inl,
     &  ipoinp,inp,cs,mcs,iaxial,ipoinpc,mi)
!
!     reading the input deck: *SOLID SECTION
!
      implicit none
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*80 matname(*),orname(*),material,orientation
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer mi(*),istartset(*),iendset(*),ialset(*),ielmat(mi(3),*),
     &  ielorien(mi(3),*),kon(*),ipkon(*),indexe,irstrt,nset,nmat,
     &  norien,
     &  istep,istat,n,key,i,j,k,l,imaterial,iorientation,ipos,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),mcs,iaxial,ipoinpc(0:*)
!
      real*8 thicke(mi(3),*),thickness,pi,cs(17,*)
!
      if((istep.gt.0).and.(irstrt.ge.0)) then
         write(*,*) '*ERROR in solidsections: *SOLID SECTION should'
         write(*,*) '  be placed before all step definitions'
         call exit(201)
      endif
!
      pi=4.d0*datan(1.d0)
!
      orientation='
     &                           '
      elset='
     &                      '
      ipos=0
!
      do i=2,n
         if(textpart(i)(1:9).eq.'MATERIAL=') then
            material=textpart(i)(10:89)
         elseif(textpart(i)(1:12).eq.'ORIENTATION=') then
            orientation=textpart(i)(13:92)
         elseif(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         else
            write(*,*) 
     &        '*WARNING in solidsections: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*SOLID SECTION%")
         endif
      enddo
!
!     check for the existence of the set,the material and orientation
!
      do i=1,nmat
         if(matname(i).eq.material) exit
      enddo
      if(i.gt.nmat) then
         do i=1,nmat
            if(matname(i)(1:11).eq.'ANISO_CREEP') then
               if(matname(i)(12:20).eq.material(1:9)) exit
            elseif(matname(i)(1:10).eq.'ANISO_PLAS') then
               if(matname(i)(11:20).eq.material(1:10)) exit
            endif
         enddo
      endif
      if(i.gt.nmat) then
         write(*,*) '*ERROR in solidsections: nonexistent material'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &"*SOLID SECTION%")
         call exit(201)
      endif
      imaterial=i
!
      if(orientation.eq.'                    ') then
         iorientation=0
      else
         do i=1,norien
            if(orname(i).eq.orientation) exit
         enddo
         if(i.gt.norien) then
            write(*,*)'*ERROR in solidsections: nonexistent orientation'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline,
     &"*SOLID SECTION%")
            call exit(201)
         endif
         iorientation=i
      endif
!
      if(ipos.eq.0) then
         write(*,*) '*ERROR in solidsections: no element set ',elset
         write(*,*) '       was been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &"*SOLID SECTION%")
         call exit(201)
      endif
      do i=1,nset
         if(set(i).eq.elset) exit
      enddo
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*) '*ERROR in solidsections: element set ',elset
         write(*,*) '  has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &"*SOLID SECTION%")
         call exit(201)
      endif
!
!     assigning the elements of the set the appropriate material
!     and orientation number
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            if((lakon(ialset(j))(1:1).eq.'B').or.
     &         (lakon(ialset(j))(1:1).eq.'S')) then
               write(*,*) '*ERROR in solidsections: *SOLID SECTION can'
               write(*,*) '       not be used for beam or shell elements
     &'
               write(*,*) '       Faulty element: ',ialset(j)
               call exit(201)
            endif
            ielmat(1,ialset(j))=imaterial
            ielorien(1,ialset(j))=iorientation
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               if((lakon(k)(1:1).eq.'B').or.
     &              (lakon(k)(1:1).eq.'S')) then
                  write(*,*) '*ERROR in solidsections: *SOLID SECTION ca
     &n'
                  write(*,*) '       not be used for beam or shell eleme
     &nts'
                  write(*,*) '       Faulty element: ',k
                  call exit(201)
               endif
               ielmat(1,k)=imaterial
               ielorien(1,k)=iorientation
            enddo
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
!     assigning a thickness to plane stress elements and an angle to
!     axisymmetric elements
!
      if((key.eq.0).or.(lakon(ialset(istartset(i)))(1:2).eq.'CA')) then
         if(key.eq.0) then
            read(textpart(1)(1:20),'(f20.0)',iostat=istat) thickness
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*SOLID SECTION%")
!
!        for axial symmetric structures:
!           thickness for axial symmetric elements: 2 degrees
!           thickness for plane stress elements: reduced by 180
!           thickness for plane strain elements: reduced by 180
!
            if(iaxial.eq.180) then
               if(lakon(ialset(istartset(i)))(1:2).eq.'CA') then
                  thickness=datan(1.d0)*8.d0/iaxial
               elseif(lakon(ialset(istartset(i)))(1:3).eq.'CPS') then
                  thickness=thickness/iaxial
               elseif(lakon(ialset(istartset(i)))(1:3).eq.'CPE') then
                  thickness=thickness/iaxial
               endif
            endif
         else
            thickness=datan(1.d0)*8.d0/iaxial
         endif
!
         do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
               if((lakon(ialset(j))(1:2).eq.'CP').or.
     &            (lakon(ialset(j))(1:2).eq.'CA')) then
                  indexe=ipkon(ialset(j))
                  do l=1,8
                     thicke(1,indexe+l)=thickness
                  enddo
               endif
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  if((lakon(k)(1:2).eq.'CP').or.
     &               (lakon(k)(1:2).eq.'CA')) then
                     indexe=ipkon(k)
                     do l=1,8
                        thicke(1,indexe+l)=thickness
                     enddo
                  endif
               enddo
            endif
         enddo
!
!        defining cyclic symmetric conditions for axisymmetric
!        elements (needed for cavity radiation)
!
         do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
               if(lakon(ialset(j))(1:2).eq.'CA') then
                  if(mcs.gt.1) then
                     write(*,*) '*ERROR in solidsections: '
                     write(*,*) '       axisymmetric elements cannot be
     &combined with cyclic symmetry'
                     call exit(201)
                  elseif(mcs.eq.1) then
                     if(int(cs(1,1)).ne.int(2.d0*pi/thickness+0.5d0)) 
     &                 then
                        write(*,*) '*ERROR in solidsections: '
                        write(*,*) '       it is not allowed to define t
     &wo different'
                        write(*,*) '       angles for an axisymmetric st
     &ructure'
                        call exit(201)
                     else
                        exit
                     endif
                  endif
                  mcs=1
                  cs(1,1)=2.d0*pi/thickness+0.5d0
                  cs(2,1)=-0.5d0
                  cs(3,1)=-0.5d0
                  cs(5,1)=1.5d0
                  do k=6,9
                     cs(k,1)=0.d0
                  enddo
                  cs(10,1)=1.d0
                  cs(11,1)=0.d0
                  cs(12,1)=-1.d0
                  cs(14,1)=0.5
                  cs(15,1)=dcos(thickness)
                  cs(16,1)=dsin(thickness)
                  exit
               endif
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
c                  if(lakon(ialset(j))(1:2).eq.'CA') then
                  if(lakon(k)(1:2).eq.'CA') then
                     if(mcs.gt.1) then
                        write(*,*) '*ERROR in solidsections: '
                        write(*,*) '       axisymmetric elements cannot 
     &be combined with cyclic symmetry'
                        call exit(201)
                     elseif(mcs.eq.1) then
                        if(int(cs(1,1)).ne.int(2.d0*pi/thickness+0.5d0)) 
     &                       then
                           write(*,*) '*ERROR in solidsections: '
                           write(*,*) '       it is not allowed to defin
     &e two different'
                           write(*,*) '       angles for an axisymmetric
     & structure'
                           call exit(201)
                        else
                           exit
                        endif
                     endif
                     mcs=1
                     cs(1,1)=2.d0*pi/thickness+0.5d0
                     cs(2,1)=-0.5d0
                     cs(3,1)=-0.5d0
                     cs(5,1)=1.5d0
                     do k=6,9
                        cs(k,1)=0.d0
                     enddo
                     cs(10,1)=1.d0
                     cs(11,1)=0.d0
                     cs(12,1)=-1.d0
                     cs(14,1)=0.5
                     cs(15,1)=dcos(thickness)
                     cs(16,1)=dsin(thickness)
                     exit
                  endif
              enddo
            endif
         enddo
!
         if(key.eq.0) then
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
         endif
      endif
!
      return
      end

