
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
      subroutine fluidsections(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,ielmat,matname,nmat,
     &  irstrt,istep,istat,n,iline,ipol,
     &  inl,ipoinp,inp,lakon,ielprop,nprop,nprop_,prop,iaxial,
     &  ipoinpc,mi)
!
!     reading the input deck: *FLUID SECTION
!
      implicit none
!
      logical liquid,manning
!
      character*1 inpc(*)
      character*7 elname
      character*8 lakon(*)
      character*80 matname(*),material,typename,typename_oil
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer mi(*),istartset(*),iendset(*),ialset(*),
     &  ielmat(mi(3),*),iaxial,
     &  irstrt,nset,nmat,ndprop,npropstart,ndpropread,
     &  istep,istat,n,key,i,j,k,imaterial,ipos,lprop,ipoinpc(0:*),
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),ielprop(*),nprop,nprop_,
     &  nodea,nodeb,noil _mat,npu,nfix,nstart
!
      real*8 prop(*)

      noil_mat=0
      liquid=.false.
      manning=.false.
!
      if((istep.gt.0).and.(irstrt.ge.0)) then
         write(*,*) 
     &      '*ERROR reading *FLUID SECTION: *FLUID SECTION should'
         write(*,*) '  be placed before all step definitions'
         call exit(201)
      endif
!
      typename='
     &                        '
!
      do i=2,n
         if(textpart(i)(1:9).eq.'MATERIAL=') then
            material=textpart(i)(10:89)
         elseif(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
!           
         elseif(textpart(i)(1:5).eq.'TYPE=') then
            read(textpart(i)(6:85),'(a80)',iostat=istat) typename
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
!           
         elseif(textpart(i)(1:4).eq.'OIL=') then
            read(textpart(i)(5:85),'(a80)',iostat=istat) typename_oil
            do j=1, nmat
               if(matname(j)(1:80).eq. typename_oil(1:80)) then
                  noil_mat=j
                  exit
               endif
            enddo
            if(j.gt.nmat) then
               write(*,*) 
     &            '*ERROR reading *FLUID SECTION: no oil with the'
               write(*,*) '       name',typename_oil,'has been defined'
               call exit(201)
            endif
!
         elseif(textpart(i)(1:6).eq.'LIQUID') then
            liquid=.true.
         elseif(textpart(i)(1:7).eq.'MANNING') then
            manning=.true.
         else
            write(*,*) 
     &      '*WARNING reading *FLUID SECTION: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
         endif
      enddo
!
!     types of sections
!
      if(typename(1:5).eq.'INOUT')then
         elname='IO     '
         ndprop=0
      elseif(typename(1:9).eq.'PIPEINOUT')then
         elname='LIPIO  '
         ndprop=0
      elseif(typename(1:12).eq.'CHANNELINOUT')then
         elname='LICHIO '
         ndprop=0
      elseif(typename(1:18).eq.'ABSOLUTETORELATIVE') then
         elname='ATR    '
         ndprop=3
!
!     version of Yavor Dobrev
!
      elseif(typename(1:10).eq.'ACCTUBEONE') then
         elname ='ACCTUBO'
!
         ndprop=19 
!
!     version of Stefanie Asenbauer
!
      elseif(typename(1:7).eq.'ACCTUBE') then
         elname ='ACCTUBE'
!
         ndprop=24
!
!     inlets and outlets for the general case (flows, pressures and
!     temperatures unknown)
!
      elseif(typename(1:5).eq.'INLET') then
         elname='INLT   '
         ndprop=0
      elseif(typename(1:6).eq.'OUTLET') then
         elname='OUTLT  '
         ndprop=0
!
!     elements "BLINDLINK" with boundary massflows
!     (pressures and temperatures unknown)
!
      elseif(typename(1:9).eq.'BLINDLINK') then
         elname='BCMF  '
         ndprop=0
!
      elseif(typename(1:12).eq.'BLEEDTAPPING') then
         elname='ORBT   '
         ndprop=39
!
      elseif(typename(1:10).eq.'CARBONSEAL') then
         ndprop=3
         if(typename(11:12).eq.'GE') then
            elname='CARBSGE'
         else
            elname='CARBS  '
         endif
!
      elseif(typename(1:7).eq.'CHANNEL') then
         if(typename(8:17).eq.'SLUICEGATE') then
            elname='LICHSG '
            ndprop=7
         elseif(typename(8:20).eq.'SLUICEOPENING') then
            elname='LICHSO '
            ndprop=6
         elseif(typename(8:16).eq.'WEIRCREST') then
            elname='LICHWE '
            ndprop=7
         elseif(typename(8:16).eq.'WEIRSLOPE') then
            elname='LICHWO '
            ndprop=6
         elseif(typename(8:17).eq.'STRAIGHT') then
            elname='LICH   '
            ndprop=7
         elseif(typename(8:16).eq.'RESERVOIR') then
            elname='LICHRE '
            ndprop=6
         elseif(typename(8:25).eq.'DISCONTINUOUSSLOPE') then
            elname='LICHDS '
            ndprop=9
         elseif(typename(8:27).eq.'DISCONTINUOUSOPENING') then
            elname='LICHDO '
            ndprop=6
         elseif(typename(8:18).eq.'CONTRACTION') then
            elname='LICHCO '
            ndprop=6
         elseif(typename(8:18).eq.'ENLARGEMENT') then
            elname='LICHEL '
            ndprop=6
         elseif(typename(8:11).eq.'STEP') then
            elname='LICHST '
            ndprop=6
         elseif(typename(8:11).eq.'DROP') then
            elname='LICHDR '
            ndprop=6
         else
            write(*,*) '*ERROR reading *FLUID SECTION:'
            write(*,*) '       unknown channel section'
            call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
         endif
!
      elseif(typename(1:14).eq.'CHARACTERISTIC') then
         ndprop=22
         elname='CHAR   '
!      
      elseif(typename(1:10).eq.'CROSSSPLIT')then
         elname='CROSPL'
         ndprop=12
!
      elseif(typename(1:12).eq.'GASPIPEFANNO') then
!     
!     version Fanno(friction and oil)
!
         if(typename(13:21).eq.'ADIABATIC') then
            if(typename(22:27).eq.'ALBERS') then
               elname='GAPFAA '
               ndprop=7
            elseif(typename(22:28).eq.'FRIEDEL') then
               elname='GAPFAF '
               ndprop=7
            else
               elname='GAPFA  '
               ndprop=7
            endif
            
         elseif(typename(13:22).eq.'ISOTHERMAL') then
            if(typename(23:28).eq.'ALBERS') then
               elname='GAPFIA '
               ndprop=7
            elseif(typename(23:29).eq.'FRIEDEL') then
               elname='GAPFIF '
               ndprop=7
            else
               elname='GAPFI  '
               ndprop=7
            endif
         endif
!
!     generalised version(friction,rotation,oil and section variation)
!
      elseif(typename(1:7).eq.'GASPIPE') then
!     
         if(typename(8:16).eq.'ADIABATIC') then
            if(typename(17:22).eq.'ALBERS') then
               elname='GAPIAA '
               ndprop=10
            elseif(typename(17:23).eq.'FRIEDEL') then
               elname='GAPIAF '
               ndprop=10
            else
               elname='GAPIA  '
               ndprop=10
            endif
         elseif(typename(8:17).eq.'ISOTHERMAL') then
            if(typename(18:23).eq.'ALBERS') then
               elname='GAPIIA '
               ndprop=10
            elseif(typename(18:24).eq.'FRIEDEL') then
               elname='GAPIIF '
               ndprop=10
            else
               elname='GAPII  '
               ndprop=10
            endif
         endif
!       
      elseif (typename(1:16).eq.'FREEDISCPUMPFLOW') then
         elname='FDPF '
         ndprop=6
!     
      elseif (typename(1:18).eq.'FREECONVECTIONFLOW') then
         elname='FCVF '
         ndprop=5
!     
      elseif (typename(1:15).eq.'MASSFLOWPERCENT') then
         elname='MFPC '
         ndprop=10
!     
      elseif(typename(1:9).eq.'LABYRINTH') then
         if(typename(10:17).eq.'FLEXIBLE') then
            ndprop=14 
            if(typename(18:23).eq.'SINGLE') then
               elname='LABFSN '
            elseif(typename(18:25).eq.'STRAIGHT') then
               elname='LABFSR '
            elseif(typename(18:24).eq.'STEPPED') then
               elname='LABFSP '
            endif
         elseif(typename(10:14).eq.'DUMMY') then
            ndprop=1
            elname='LABD   '
         else
            ndprop=11 
            if(typename(10:15).eq.'SINGLE') then
               elname='LABSN  '
            elseif(typename(10:17).eq.'STRAIGHT') then
               elname='LABSR  '
            elseif(typename(10:16).eq.'STEPPED') then
               elname='LABSP  '
            endif 
         endif
!     
      elseif(typename(1:10).eq.'LIQUIDPUMP') then
         ndprop=20
         elname='LIPU   '
!     
      elseif(typename(1:8).eq.'MOEHRING') then
         if(typename(9:19).eq.'CENTRIFUGAL') then
            elname='MRGF   '
         elseif(typename(9:19).eq.'CENTRIPETAL') then
            elname='MRGP   '
         endif
         ndprop=10
!
      elseif(typename(1:7).eq.'ORIFICE') then
         ndprop=8
         if(typename(8:13).eq.'OWN_AL') then
            elname='ORMA   '
         elseif(typename(8:12).eq.'PK_AL') then
            elname='ORPA   '
         elseif(typename(8:12).eq.'MS_MS') then
            elname='ORMM   '
         elseif(typename(8:12).eq.'PK_MS') then
            elname='ORPM   '
         elseif(typename(8:12).eq.'BRAGG') then
             elname='ORBG   '
         else 
            elname='ORC1   '
         endif
!
      elseif(typename(1:4).eq.'PIPE') then
         if(typename(5:8).eq.'BEND') then
            elname='LIPIBE '
            ndprop=4
         elseif(typename(5:10).eq.'BRANCH') then
            elname='LIPIBR '
            ndprop=6
         elseif(typename(5:15).eq.'CONTRACTION') then
            elname='LIPICO '
            ndprop=2
         elseif(typename(5:13).eq.'DIAPHRAGM') then
            elname='LIPIDI '
            ndprop=2
         elseif(typename(5:15).eq.'ENLARGEMENT') then
            elname='LIPIEL '
            ndprop=2
         elseif(typename(5:12).eq.'ENTRANCE') then
            elname='LIPIEN '
            ndprop=2
         elseif(typename(5:13).eq.'GATEVALVE') then
            elname='LIPIGV '
            ndprop=2
         elseif(typename(5:11).eq.'MANNING') then
            if(typename(12:19).eq.'FLEXIBLE') then
               elname='LIPIMAF'
               ndprop=4
            else
               elname='LIPIMA '
               ndprop=3
            endif
         elseif(typename(5:19).eq.'WHITE-COLEBROOK') then
            if(typename(20:27).eq.'FLEXIBLE') then
               elname='LIPIWCF'
               ndprop=6
            else
               elname='LIPIWC '
               ndprop=5
            endif
         endif
!
      elseif(typename(1:14).eq.'PRESWIRLNOZZLE') then
         elname='ORPN   '
         ndprop=39
!
      elseif(typename(1:18).eq.'RELATIVETOABSOLUTE') then
         elname='RTA    '
         ndprop=3
!
      elseif(typename(1:10).eq.'RESTRICTOR') then
         if(typename(11:15).eq.'USER') then
            elname='REUS   '
              ndprop=6
!          
         elseif(typename(11:29).eq.'LONGORIFICEIDELCHIK') then
            elname='RELOID '
            ndprop=6
!            
         elseif(typename(11:21).eq.'WALLORIFICE') then 
            elname='REWAOR '
            ndprop=6
         elseif(typename(11:18).eq.'ENTRANCE') then
            elname='REEN   '
            ndprop=6
         elseif(typename(11:21).eq.'ENLARGEMENT') then
            elname='REEL   '
             ndprop=5
         elseif(typename(11:21).eq.'CONTRACTION') then
            elname='RECO   '
            ndprop=7
         elseif(typename(11:23).eq.'BENDIDELCIRC') then
            elname='REBEIDC'
            ndprop=7
         elseif(typename(11:23).eq.'BENDIDELRECT') then
            elname='REBEIDR'
            ndprop=9
         elseif(typename(11:20).eq.'BENDMILLER') then
            elname='REBEMI '
            ndprop=7
         elseif(typename(11:17).eq.'BENDOWN') then
            elname='REBEMA '
            ndprop=7
         elseif(typename(11:14).eq.'EXIT') then
            elname='REEX   '
            ndprop=6
         elseif(typename(11:33).eq.'LONGORIFICELICHTAROWICZ') then
            elname='RELOLI '
            ndprop=6
         endif
!     
      elseif(typename(1:6).eq.'BRANCH') then
         if(typename(7:13).eq.'JOINTGE')then
            elname='REBRJG '
            ndprop=11
         elseif(typename(7:20).eq.'JOINTIDELCHIK1')then
            elname='REBRJI1'
            ndprop=11
         elseif(typename(7:20).eq.'JOINTIDELCHIK2')then
            elname='REBRJI2'
            ndprop=11 
         elseif(typename(7:13).eq.'SPLITGE')then
            elname='REBRSG '
            ndprop=11  
         elseif(typename(7:20).eq.'SPLITIDELCHIK1')then
            elname='REBRSI1'
            ndprop=15
         elseif(typename(7:20).eq.'SPLITIDELCHIK2')then
            elname='REBRSI2'
            ndprop=13
         endif
!
      elseif(typename(1:7).eq.'RIMSEAL') then
         elname='RIMS   '
         ndprop=5
! 
      elseif(typename(1:14).eq.'ROTATINGCAVITY') then
         if(typename(15:20).eq.'LINEAR')then
            elname='RCVL   '
         else if(typename(15:23).eq.'NONLINEAR')then
            elname='RCVN   '
         endif
         ndprop=5
!
      elseif(typename(1:13).eq.'RADIALOUTFLOW') then
         if(typename(14:24).eq.'RADIALINLET')then
            elname='RORI   '

         else if(typename(14:23).eq.'AXIALINLET')then
            elname='ROAI   '
         endif
         ndprop=6
!
      else if (typename(1:6).eq.'S-PUMP') then
         elname='SPUMP  '
         ndprop=5
!
      else if (typename(1:6).eq.'O-PUMP') then
         elname='LDOP   '
         ndprop=3
!         
      elseif(typename(1:6).eq.'VORTEX') then
         if(typename(7:10).eq.'FREE') then
            elname='VOFR   '
            ndprop=9
         elseif(typename(7:12).eq.'FORCED') then
            elname='VOFO   '
            ndprop=7
         endif
!     
      elseif(typename(1:1).eq.' ') then
         elname='       '
         ndprop=0
      else
         write(*,*) '*ERROR reading *FLUID SECTION: ',typename
         write(*,*) '       is an unknown fluid section type'
         call exit(201)
      endif
!
!     check for the existence of the set and the material
!
      do i=1,nmat
         if(matname(i).eq.material) exit
      enddo
      if(i.gt.nmat) then
         write(*,*) 
     &     '*ERROR reading *FLUID SECTION: nonexistent material'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
         call exit(201)
      endif
      imaterial=i
!
      do i=1,nset
         if(set(i).eq.elset) exit
      enddo
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*) '*ERROR reading *FLUID SECTION: element set ',elset
         write(*,*) '  has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
         call exit(201)
      endif
!
      npropstart=nprop
!
!     reading the element properties
!
!     liquid pump or gas characteristic
!
      if((elname(1:4).eq.'LIPU').or.(elname(1:4).eq.'CHAR').or.
     &   (elname(1:4).eq.'ORPN').or.(elname(1:4).eq.'ORBT')) then
!
!        determine the number of fixed inputs (before entering
!        a curve with arbitary many data points)
!
         if(elname(1:4).eq.'LIPU') then
            nfix=0
         elseif(elname(1:4).eq.'CHAR') then
            nfix=1
         elseif(elname(1:4).eq.'ORPN') then
            nfix=5
         elseif(elname(1:4).eq.'ORBT') then
            nfix=3
         endif
         ndprop=0
         npu=0
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            if(ndprop.eq.0) then
               do j=1,nfix
                  ndprop=ndprop+1
                  read(textpart(j),'(f40.0)',iostat=istat) 
     &                 prop(nprop+ndprop)
                  if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
               enddo
               nstart=nfix+1
!
!              adding the point (0.,0.) for characteristics
!
               if(elname(1:4).eq.'CHAR') then
                  prop(nprop+ndprop+2)=0.d0
                  prop(nprop+ndprop+3)=0.d0
                  ndprop=ndprop+2
                  npu=npu+2
               endif
            else
               nstart=1
            endif
            do j=nstart,n
               npu=npu+1
               ndprop=ndprop+1
               read(textpart(j),'(f40.0)',iostat=istat) 
     &              prop(nprop+ndprop+1)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
            enddo
         enddo
!
!        check whether data points are paired
!
         if(2*(npu/2).ne.npu) then
            write(*,*) '*ERROR reading *FLUID SECTION: less y data'
            write(*,*) '       points x data points for fluid'
            write(*,*) '       section type',elname(1:4)
            call exit(201)
         endif
!
         prop(nprop+nfix+1)=npu/2
!
!        check whether data points are in appropriate order
!        (only for liquid pump so far)
!
         if(elname(1:4).eq.'LIPU') then
            do j=2,npu/2,2
               if(prop(nprop+nfix+2*j).le.prop(nprop+nfix+2*j-2)) then
                  write(*,*) '*ERROR reading *FLUID SECTION: volumetric'
                  write(*,*) '       flow data points must be in'
                  write(*,*) '       strict increasing order'
                  call exit(201)
               elseif(prop(nprop+nfix+2*j+1).gt.prop(nprop+nfix+2*j-1))
     &              then
                  write(*,*) '*ERROR reading *FLUID SECTION: total'
                  write(*,*) '       head data points must be '
                  write(*,*) '       decreasing'
                  call exit(201)
               endif
            enddo
         endif
!
c         if(iaxial.ne.0) then
c            do j=1,npu/2,2
c               prop(nprop+nfix+j+1)=prop(nprop+nfix+j+1)/iaxial
c            enddo
c         endif
         nprop=nprop+ndprop+1
! 
      elseif(elname.eq.'ACCTUBE') then
!         
!        Read the first 20 elements
! 
         lprop=0
         ndpropread=ndprop
!         
         do j=1,(ndpropread-1)/8+1
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            do k=1,8
               lprop=lprop+1
               if(lprop.gt.ndpropread) exit
               read(textpart(k),'(f40.0)',iostat=istat)
     &                 prop(nprop+lprop)
               
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
               
!              If 20 elements have been read, check how many more
!              are to be read
               if (lprop.eq.20) then
                  ndpropread = lprop +
     &                        (prop(nprop+19)+prop(nprop+20))*2
                  ndprop = ndpropread
               endif
            enddo
         enddo
!         
!        Until now 3 lines have been read
!        If necessary read the rest
!
         if(ndpropread.gt.lprop) then
            do j=3,(ndpropread-1)/8+1
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               if((istat.lt.0).or.(key.eq.1)) exit
               do k=1,8
                  lprop=lprop+1
                  if(lprop.gt.ndpropread) exit
                  read(textpart(k),'(f40.0)',iostat=istat)
     &                 prop(nprop+lprop)
               
                  if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
               enddo
            enddo
         endif
!
         nprop=nprop+ndprop
!
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
!         
      elseif(ndprop.gt.0) then
!
!        flexible liquid pipe
!
         if((elname.eq.'LIPIMAF').or.(elname.eq.'LIPIWCF')) then
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
            read(textpart(1),'(i10)',iostat=istat) nodea
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
            prop(nprop+1)=nodea+0.5d0
            read(textpart(2),'(i10)',iostat=istat) nodeb
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
            prop(nprop+2)=nodeb+0.5d0
            if(elname.eq.'LIPIMAF') then
               read(textpart(3),'(f40.0)',iostat=istat)
     &              prop(nprop+3)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
               prop(nprop+4)=iaxial+0.5d0
            else
               do j=3,5
                  read(textpart(j),'(f40.0)',iostat=istat)
     &                 prop(nprop+j)
                  if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
               enddo            
               prop(nprop+6)=iaxial+0.5d0
            endif
            nprop=nprop+ndprop
         else
!
!           reading the properties of all other elements
!
            lprop=0
            ndpropread=ndprop
            if(elname.eq.'REBEIDR') then
               ndpropread=ndpropread-1
            endif
            do j=1,(ndpropread-1)/8+1
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               do k=1,8
                  lprop=lprop+1
                  if(lprop.gt.ndpropread) exit
                  read(textpart(k),'(f40.0)',iostat=istat)
     &                 prop(nprop+lprop)
                  if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*FLUID SECTION%")
               enddo
            enddo
!     
!           reducing the area in case of an axisymmetric model
!
c            if((elname(1:3).ne.'LAB').and.(elname(1:5).ne.'GAPFF')
c     &           .and.(elname(1:3).ne.'ORF')
c     &           .and.(elname(1:5).ne.'GAPIF')) then
c               if(iaxial.ne.0) prop(nprop+1)=prop(nprop+1)/iaxial
c            endif
            nprop=nprop+ndprop
         endif
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
      else
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
      endif
!
      if(nprop.gt.nprop_) then
         write(*,*) '*ERROR reading *FLUID SECTION: increase nprop_'
         call exit(201)
      endif
!
      if(elname(1:4).eq.'RIMS') then
         prop(npropstart+1)= prop(npropstart+1)+0.5d0
         prop(npropstart+2)= prop(npropstart+2)+0.5d0
      endif
!
!     complementing the properties of restrictor elements
!
      if((elname(1:6).eq.'REBEMI').or.
     &     (elname(1:6).eq.'REBEMA')) then
         prop(npropstart+7)=noil_mat+0.5d0
!
      elseif(elname(1:7).eq.'REBEIDC') then
         prop(npropstart+7)=noil_mat+0.5d0
!     
      elseif(elname(1:7).eq.'REBEIDR') then
         prop(npropstart+9)=noil_mat+0.5d0
!     
      elseif(elname(1:4).eq.'REEX') then
!       
            prop(npropstart+6)=noil_mat+0.5d0
            prop(npropstart+5)=prop(npropstart+4)
            prop(npropstart+4)=prop(npropstart+3)
            prop(npropstart+3)=prop(npropstart+2)
            prop(npropstart+2)=100000*prop(npropstart+1)
!     
      elseif(elname(1:6).eq.'REWAOR') then
!        
            prop(npropstart+6)=noil_mat+0.5d0
            prop(npropstart+5)=prop(npropstart+4)
            prop(npropstart+4)=prop(npropstart+3)
            prop(npropstart+3)=prop(npropstart+2)
            prop(npropstart+2)=prop(npropstart+1)
            prop(npropstart+1)=100000*prop(npropstart+2)
!
      elseif(elname(1:4).eq.'REEN') then
!        
            prop(npropstart+6)=noil_mat+0.5d0
            prop(npropstart+5)=prop(npropstart+3)
            prop(npropstart+4)=0.5d0
            prop(npropstart+3)=prop(npropstart+2)
            prop(npropstart+2)=prop(npropstart+1)
            prop(npropstart+1)=100000*prop(npropstart+2)
!
      elseif(elname(1:7).eq.'REBRJI1') then
         prop(npropstart+11)=noil_mat+0.5d0
!
      elseif(elname(1:7).eq.'REBRJI2') then
         prop(npropstart+11)=noil_mat+0.5d0
         if(1.d0-(prop(npropstart+5)+prop(npropstart+6))/
     &        prop(npropstart+4).gt.0.01d0)then
            write(*,*) '*ERROR: reading *FLUID SECTION:'
            write(*,*) '        in element type RESTRICTOR 
     &                           BRANCH JOINT IDELCHIK2'
            write(*,*) '        A0 ist not equal to A1+A2'
            call exit(201)
         endif
      elseif(elname(1:6).eq.'REBRJG') then
         prop(npropstart+11)=noil_mat+0.5d0
      elseif(elname(1:6).eq.'REBRSG') then
         prop(npropstart+11)=noil_mat+0.5d0
      elseif(elname(1:7).eq.'REBRSI1') then
         prop(npropstart+15)=noil_mat+0.5d0
         if(dabs(prop(npropstart+13)).le.1E-5) then
            prop(npropstart+13)=1.d0
         endif
         if(dabs(prop(npropstart+14)).le.1E-5) then
            prop(npropstart+14)=1.d0
         endif
!     
      elseif(elname(1:7).eq.'REBRSI2') then
         prop(npropstart+13)=noil_mat+0.5d0
         if(dabs(prop(npropstart+11)).le.1E-5) then
            prop(npropstart+11)=1.d0
         endif
         if(dabs(prop(npropstart+12)).le.1E-5) then
            prop(npropstart+12)=1.d0
         endif
!     
      endif
!
!     "transforming" integer inputs into real inputs
!
      if(elname(1:4).eq.'REBR') then
!
!     for branches element the three first inputs should be integers 
!     since they are respectively the label number of the element
!     modeling branch 0 , 1 and 2.
!     Since all inputs are real the "integer input" are augmented by 0.5
!     and will be converted in integer using fucntion INT
         prop(npropstart+3)=prop(npropstart+3)+0.5d0
         prop(npropstart+2)=prop(npropstart+2)+0.5d0
         prop(npropstart+1)=prop(npropstart+1)+0.5d0
!
!     number of the table to be used for Preswirlnozzles
!
      elseif(elname(1:4).eq.'ORPN') then
         prop(npropstart+4)=prop(npropstart+4)+0.5d0
!
!     number of the table to be used for bleed tappings
!
      elseif(elname(1:4).eq.'ORBT') then
         prop(npropstart+3)=prop(npropstart+3)+0.5d0
!
!     label number of the swirl generating element for free vortices
!
      elseif(elname(1:4).eq.'VOFR') then
         if(prop(npropstart+6).ne.0d0) then
            prop(npropstart+6)=prop(npropstart+6)+0.5d0
         endif
!
!     k_oil for restrictors type USER, ENTRENCE, LONG ORIFICE IDELCHIK,
!     EXIT, ORIFICE IN A WALL, LONG ORIFICE LICHTAROWICZ
!
      elseif((elname(1:4).eq.'REUS').or.
     &         (elname(1:6).eq.'RELOID').or.  
c     &         (elname(1:4).eq.'REEN').or.  
c     &         (elname(1:4).eq.'REEX').or.   
c     &         (elname(1:6).eq.'REWAOR').or. 
     &         (elname(1:6).eq.'RELOLI')) then

         prop(npropstart+6)= noil_mat+0.5d0

!
!     k_oil for restrictors type ENLARGEMENT
!
      elseif (elname(1:4).eq.'REEL') then
         prop(npropstart+5)= noil_mat+0.5d0
!
!     k_oil for restrictors type CONTRACTION, BEN MILLER, BEND MA
!     (BEND IDELCHICK has already been treated)
!
      elseif(elname(1:4).eq.'RECO') then
        prop(npropstart+7)= noil_mat+0.5d0
!
!     label number of the upstream and downstream nodes for Moehrings
!
      elseif(elname(1:2).eq.'MR') then
         prop(npropstart+5)= prop(npropstart+5)+0.1d0
         prop(npropstart+6)= prop(npropstart+6)+0.1d0
!
!     k_oil for pipe elements
!
      elseif(elname(1:3).eq.'GAP') then
         prop(npropstart+7)=noil_mat+0.5
!
!     Node1 and Node2 for LABYRINTH FLEXIBLE
!
      elseif(elname(1:3).eq.'LAB') then
         if((elname(1:4).eq.'LABF')) then
            prop(npropstart+14)= prop(npropstart+13)+0.5d0
            prop(npropstart+13)= prop(npropstart+12)
            prop(npropstart+12)= prop(npropstart+11)
            prop(npropstart+11)= prop(npropstart+10)
            prop(npropstart+10)= prop(npropstart+9)
            prop(npropstart+9)= prop(npropstart+8)
            prop(npropstart+8)= prop(npropstart+7)
            prop(npropstart+7)= prop(npropstart+6)
            prop(npropstart+6)= prop(npropstart+5)
            prop(npropstart+5)= prop(npropstart+4)
            prop(npropstart+4)= prop(npropstart+3)
            prop(npropstart+3)=iaxial+0.5d0
         elseif(elname(1:4).eq.'LABD') then
            prop(npropstart+1)=prop(npropstart+1)+0.5d0
         else
            prop(npropstart+11)= prop(npropstart+10)
            prop(npropstart+10)= prop(npropstart+9)
            prop(npropstart+9)= prop(npropstart+8)
            prop(npropstart+8)= prop(npropstart+7)
            prop(npropstart+7)= prop(npropstart+6)
            prop(npropstart+6)= prop(npropstart+5)
            prop(npropstart+5)= prop(npropstart+4)+0.5d0
            prop(npropstart+4)= prop(npropstart+3)
            prop(npropstart+3)= iaxial+0.5d0
         endif
      elseif(elname(1:4).eq.'LICH') then
!
!        upstream and downstream reference elements for channels
!
         if((elname(5:6).eq.'SO').or.
     &      (elname(5:6).eq.'WO').or.
     &      (elname(5:6).eq.'  ').or.
     &      (elname(5:6).eq.'RE').or.
     &      (elname(5:6).eq.'DS').or.
     &      (elname(5:6).eq.'DO')) then
            prop(npropstart+6)=prop(npropstart+6)+0.5d0
         elseif(elname(5:6).eq.'SG') then
            prop(npropstart+5)=prop(npropstart+5)+0.5d0
            prop(npropstart+7)=prop(npropstart+7)+0.5d0
         elseif(elname(5:6).eq.'WE') then
            prop(npropstart+5)=prop(npropstart+5)+0.5d0
         endif
         if(elname(5:6).eq.'DS') then
            prop(npropstart+8)=prop(npropstart+8)+0.5d0
         endif
!
!        if Manning: change sign of 5th entry as marker
!
         if(manning) then
            if((elname(5:6).eq.'SO').or.
     &         (elname(5:6).eq.'WO').or.
     &         (elname(5:6).eq.'  ').or.
     &         (elname(5:6).eq.'RE').or.
     &         (elname(5:6).eq.'DS').or.
     &         (elname(5:6).eq.'DO')) then
               prop(npropstart+5)=-prop(npropstart+5)
            endif
         endif
      endif         
!     
!     check the range of the parameters
!
      if((elname(1:4).ne.'LIPU').and.(elname.ne.'       ').and.
     &   (elname(1:3).ne.'LAB').and.(elname(1:4).ne.'CHAR').and.
     &   (elname(1:5).ne.'CARBS')) then
         if(prop(npropstart+1).lt.0.d0) then
            write(*,*) '*ERROR reading *FLUID SECTION: section area'
            write(*,*) '       is not positive'
            call exit(201)
         endif
      endif
!
      if((elname(1:2).eq.'OR').and.(elname(1:4).ne.'ORC1').and.
     &     (elname(1:4).ne.'ORBT') .and.(elname(1:4).ne.'ORPN').and.
     &     (elname(1:4).ne.'ORFL')) then
         if(prop(npropstart+2).lt.0.d0) then
            write(*,*) '*ERROR reading *FLUID SECTION: diameter of the'
            write(*,*) '       orifice is not positive'
            call exit(201)
         endif
         if(prop(npropstart+3).lt.0.d0) then
            write(*,*) '*ERROR reading *FLUID SECTION: length of the'
            write(*,*) '       orifice is not positive'
            call exit(201)
         endif
         if((prop(npropstart+4).gt.1.d-20).and.
     &      (prop(npropstart+5).gt.1.d-20)) then
            write(*,*) 
     &         '*ERROR reading *FLUID SECTION: either the radius'
            write(*,*) '       of the orifice must be zero or the'
            write(*,*) '       chamfer angle'
            call exit(201)
         endif
         if(prop(npropstart+4)/prop(npropstart+2).lt.0.d0) then
            write(*,*) 
     &          '*ERROR reading *FLUID SECTION: r/d of an orifice'
            write(*,*) '       must not be negative'
            call exit(201)
         endif
         if(prop(npropstart+4)/prop(npropstart+2).gt.0.82) then
            write(*,*) 
     &         '*ERROR reading *FLUID SECTION: r/d of an orifice'
            write(*,*) '       must not not exceed 0.82'
c            call exit(201)
         endif
         if(prop(npropstart+5).lt.0.d0) then
            write(*,*) 
     &         '*ERROR reading *FLUID SECTION: the chamfer angle'
            write(*,*) '       of an orifice must not be negative'
            call exit(201)
         endif
         if(prop(npropstart+5).gt.90.d0) then
            write(*,*) 
     &         '*ERROR reading *FLUID SECTION: the chamfer angle'
            write(*,*) '       of an orifice must not exceed 90�'
c            call exit(201)
         endif
         if(prop(npropstart+6).lt.0.d0) then
            write(*,*) '*ERROR reading *FLUID SECTION: d/D (orifice)'
            write(*,*) '       must not be negative'
            call exit(201)
         endif
         if(prop(npropstart+6).gt.0.5d0) then
            write(*,*) '*ERROR reading *FLUID SECTION: d/D (orifice)'
            write(*,*) '       must not exceed 0.5'
c            call exit(201)
         endif
         if(prop(npropstart+3)/prop(npropstart+2).lt.0.d0) then
            write(*,*) 
     &         '*ERROR reading *FLUID SECTION: L/d of an orifice'
            write(*,*) '       must not be negative'
            call exit(201)
         endif
         if(elname(3:3).eq.'P') then
            if(prop(npropstart+3)/prop(npropstart+2).gt.2.d0) then
               write(*,*) 
     &           '*ERROR reading *FLUID SECTION: L/d of an orifice'
               write(*,*) '       with Parker must not exceed 2'
               call exit(201)
            endif
!         else
!            if(prop(npropstart+3)/prop(npropstart+2).gt.10.d0) then
!               write(*,*) '*ERROR reading *FLUID SECTION: L/d of an orifice'
!               write(*,*) '        must not exceed 10'
!               call exit(201)
!            endif
         endif
      elseif(elname(1:4).eq.'ORBT') then
         if(prop(npropstart+2).lt.0.d0) then
            write(*,*) '*ERROR reading *FLUID SECTION: ps1/pt1 (bleed'
            write(*,*) '       tapping) must not be negative'
            call exit(201)
         endif
         if(prop(npropstart+2).gt.1.d0) then
            write(*,*) '*ERROR reading *FLUID SECTION: ps1/pt1 (bleed'
            write(*,*) '       tapping) must not exceed 1'
c            call exit(201)
         endif
      elseif(elname(1:4).eq.'ORPN') then
         if(prop(npropstart+2).lt.0.d0) then
            write(*,*) '*ERROR reading *FLUID SECTION: theta (preswirl'
            write(*,*) '       nozzle) must not be negative'
            call exit(201)
         endif
         if(prop(npropstart+2).gt.90.d0) then
            write(*,*) '*ERROR reading *FLUID SECTION: theta (preswirl'
            write(*,*) '       nozzle) must not exceed 90�'
c            call exit(201)
         endif
         if(prop(npropstart+3).lt.0.d0) then
            write(*,*) '*ERROR reading *FLUID SECTION: k_phi (preswirl'
            write(*,*) '       nozzle) must not be negative'
            call exit(201)
         endif
         if(prop(npropstart+3).gt.1.05d0) then
            write(*,*) '*ERROR reading *FLUID SECTION: k_phi (preswirl'
            write(*,*) '       nozzle) must not exceed 1.05'
c            call exit(201)
         endif
      elseif(elname(1:4).eq.'ORBG') then
         if(prop(npropstart+1).lt.0.d0) then
            write(*,*) '*ERROR reading *FLUID SECTION: section area'
            write(*,*) '       is not positive'
            call exit(201)
         endif
         if(prop(npropstart+2).lt.0.d0 .or.
     &      prop(npropstart+2).ge.1.d0) then
            write(*,*) 
     &         '*ERROR reading *FLUID SECTION: using Bragg Method'
            write(*,*) 'Cd by crtitical pressure ratio '
            write(*,*) '*FLUID SECTIONS position 2'
            write(*,*) '0 < Cd _crit < 1'
c            call exit(201)
         endif
      endif
!
      if((elname(1:3).eq.'LAB').and.(elname(1:4).ne.'LABF').and.
     &    (elname(1:4).ne.'LABD') ) then
         if((prop(npropstart+1).gt.1000.d0)
     &        .or.(prop(npropstart+1).lt.0.d0)) then
            write(*,*) 
     &         '*ERROR reading *FLUID SECTION: the selected pitch t'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=t<=1000 mm'
c            call exit(201)
         endif
         if((prop(npropstart+2).gt.100d0)
     &        .or.(prop(npropstart+2).lt.0.d0)) then
            write(*,*) 
     &          '*ERROR reading *FLUID SECTION: the selected gap s'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=s<=100 mm'
c            call exit(201)
         endif
         if((prop(npropstart+3).gt.5000.d0)
     &        .or.(prop(npropstart+3).lt.0.d0)) then
            write(*,*) '*ERROR reading *FLUID SECTION:'
            write(*,*) 'the selected diameter d'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=d<=5000'
c            call exit(201)
         endif
         if((prop(npropstart+5).gt.9.d0)
     &        .or.(prop(npropstart+5).lt.0.d0)) then
            write(*,*) '*ERROR reading *FLUID SECTION:'
            write(*,*)'the selected spike number n'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=n<=9'
            write(*,*) prop(npropstart+4)
c           call exit(201)
        endif
         if((prop(npropstart+5).gt.100.d0)
     &       .or.(prop(npropstart+5).lt.0.d0)) then
            write(*,*) '*ERROR reading *FLUID SECTION:'
            write(*,*) 'the selected spike breadth'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=b<=100 mm'
c            call exit(201)
         endif
         if((prop(npropstart+6).gt.9.d0)
     &        .or.(prop(npropstart+6).lt.0.d0)) then
            write(*,*) '*ERROR reading *FLUID SECTION:'
            write(*,*) 'the selected spike height'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=b<=20 mm'
c            call exit(201)
         endif
         if((prop(npropstart+7).gt.4.d0)
     &        .or.(prop(npropstart+7).lt.0.d0)) then
            write(*,*) '*ERROR reading *FLUID SECTION:'
            write(*,*) 'the selected Honeycomb cell width'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=L<=4 mm'
c            call exit(201)
         endif
         if((prop(npropstart+8).gt.5.d0)
     &        .or.(prop(npropstart+8).lt.0.d0)) then
            write(*,*) '*ERROR reading *FLUID SECTION:'
            write(*,*) 'the selected edge radius'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=r<=5 mm'
c            call exit(201)
         endif
         if((prop(npropstart+9).gt.100.d0)
     &        .or.(prop(npropstart+9).lt.0.d0)) then
            write(*,*) '*ERROR reading *FLUID SECTION:'
            write(*,*) 'the selected position of the spike'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=X<=0.1 mm'
c            call exit(201)
         endif
         if((prop(npropstart+10).gt.100.d0)
     &        .or.(prop(npropstart+10).lt.0.d0)) then
            write(*,*) '*ERROR reading *FLUID SECTION:'
            write(*,*) 'the selected height of the step'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=Hst<=0.1 mm'
c            call exit(201)
         endif
      endif
!
      if(elname(1:5).eq.'CARBS') then
         if(elname(1:7).eq.'CARBSGE') then
            if(prop(npropstart+1).le.0.d0) then
               write(*,*) '*ERROR reading *FLUID SECTION:'
               write(*,*) 'the selected diameter of the carbon seal'
               write(*,*) 'has been defined as less or equal to 0'
               call exit(201)
            endif
         else
            if((prop(npropstart+1).le.0.d0)
     &           .or.(prop(npropstart+2).le.0.d0)
     &           .or.(prop(npropstart+3).le.0.d0)) then
               write(*,*) '*ERROR reading *FLUID SECTION:'
               write(*,*) 'the selected diameter'
               write(*,*) 'or the selected length' 
               write(*,*) 'or the selected gap of the carbon seal'
               write(*,*) 'has been defined as less or equal to 0'
               call exit(201)
            endif
         endif
      endif
!
      if(elname(1:3).eq.'RCVL') then
         if(prop(npropstart+2).lt.(prop(npropstart+3))) then
             write(*,*) '*ERROR reading *FLUID SECTION: '
             write(*,*) 'element TYPE=ROTATING CAVITY(Radial inflow)'
             write(*,*) 'the specified upstream radius is smaller than'
             write(*,*) 'the specified downstream radius!'
             write(*,*) 'Please check the element definition.'
             call exit(201)
         endif
      endif 
      if(elname(1:3).eq.'RCVN') then
         if(prop(npropstart+1).lt.(prop(npropstart+2))) then
             write(*,*) '*ERROR reading *FLUID SECTION: '
             write(*,*) 'element TYPE=ROTATING CAVITY(Radial inflow)'
             write(*,*) 'the specified upstream radius is smaller than'
             write(*,*) 'the specified downstream radius!'
             write(*,*) 'Please check the element definition.'
             call exit(201)
         endif
      endif
!
      if(elname(1:2).eq.'RE') then
         if((prop(npropstart+1).le.0)
     &        .or.(prop(npropstart+2).le.0) 
     &        .or.(prop(npropstart+3).le.0))then
            write(*,*) '*ERROR reading *FLUID SECTION:'
            write(*,*) 'A1,A2 or Dh less or equal 0'
            call exit(201)
!     
         elseif(elname(1:4).eq.'REEL') then
            if(prop(npropstart+1).ge.(prop(npropstart+2))) then
               write(*,*) '*ERROR reading *FLUID SECTION:'
               write(*,*) 'Section A1 is greater than section A2'
               call exit(201)
            endif
!
         elseif(elname(1:4).eq.'RECO') then
            if(prop(npropstart+1).lt.(prop(npropstart+2))) then
               write(*,*) '*ERROR reading *FLUID SECTION:'
               write(*,*) 'Section A2 is greater than section A1'
               call exit(201)
            endif
         endif
!     
         if(elname(1:4).eq.'REBR') then
            if((prop(npropstart+1).le.0)
     &           .or.(prop(npropstart+2).le.0) 
     &           .or.(prop(npropstart+3).le.0))then
               write(*,*) '*ERROR reading *FLUID SECTION:'
               write(*,*) 'trying to define a branch '
               write(*,*) 'all three elements must be different from 0'
               call exit(201)
!     
            elseif((prop(npropstart+4).le.0)
     &              .or.(prop(npropstart+5).le.0) 
     &              .or.(prop(npropstart+6).le.0))then 
               write(*,*) '*ERROR reading *FLUID SECTION:'
               write(*,*) 'trying to define a branch '
               write(*,*) 'all sections must be positive'
               call exit(201)
!
            elseif((prop(npropstart+7).lt.0)
     &              .or.(prop(npropstart+8).lt.0))then
               write(*,*) '*ERROR reading *FLUID SECTION:'
               write(*,*) 'trying to define a branch '
               write(*,*) 'alpha1 & 2 cannot be negative'
               call exit(201)
!     
            elseif((prop(npropstart+7).gt.90)
     &              .or.(prop(npropstart+8).gt.90)) then
               write(*,*) '*ERROR reading *FLUID SECTION:'
               write(*,*) 'trying to define a branch '
               write(*,*) 'alpha1 & 2 cannot greater than 90 gegrees'
               call exit(201)
!     
            elseif((elname(5:7).eq.'SI1')
     &              .or.(elname(5:7).eq.'JI1')
     &              .or.(elname(5:7).eq.'JI2')) then
               if(prop(npropstart+7).gt.0) then
                  write(*,*) '*ERROR reading *FLUID SECTION:'
                  write(*,*) 'trying to define a branch '
                  write(*,*) 'Type IDELCHIK JOINT1 or SPLIT1&2'
                  write(*,*) 'alpha1 must be 0 degrees'
                  call exit(201)
               endif
!
            elseif((elname(5:7).eq.'SI1').or.
     &              (elname(5:7).eq.'JI1'))then
               if(prop(npropstart+4).ne.(prop(npropstart+5)))then
                  write(*,*) '*ERROR reading *FLUID SECTION:'
                  write(*,*) 'trying to define a branch '
                  write(*,*) 'Type IDELCHIK SPLIT1 or JOINT1 '
                  write(*,*) 'A1=A0'
                  call exit(201)
               endif
! 
            elseif(elname(5:7).eq.'JI2') then
               if((prop(npropstart+5)+(prop(npropstart+6)))
     &              .ne.prop(npropstart+4))then
                  write(*,*) '*ERROR reading *FLUID SECTION:'
                  write(*,*) 'trying to define a branch '
                  write(*,*) 'Type IDELCHIK JOINT2 '
                  write(*,*) 'A1+A2 must be equal to A0'
                  call exit(201)
               endif
            endif
!
!     General Vortex
!     
         elseif(elname(1:2).eq.'VO') then
!
!     inner and outer radius less or equal to 0
            if((prop(npropstart+1).le.0) .or.
     &           (prop(npropstart+2).le.0)) then
               write(*,*)'*ERROR reading *FLUID SECTION:'
               write(*,*)'trying to define a VORTEX'
               write(*,*)'R1 and R2 must be positive'
               call exit(201)
!     
            elseif(prop(npropstart+3).le.0) then
               write(*,*)'*ERROR reading *FLUID SECTION:'
               write(*,*)'trying to define a VORTEX'
               write(*,*)'eta must be different positive'
               call exit(201)
            endif
!
!     FREE VORTEX
!
         elseif(elname(1:4).eq.'VOFR') then
!     
!     the swirl comes from another upstream  element
            if(prop(npropstart+6).ne.0d0) then
!     
!     the rotation speed must be 0
               if(prop(npropstart+7).ne.0) then
                  write(*,*)'*ERROR reading *FLUID SECTION:'
                  write(*,*)'trying to define a FREE VORTEX'
                  write(*,*)'rotation speed and upstream element'
                  write(*,*)'cannot be simultaneously used '
                  call exit(201)
               endif
            endif
!
!     FORCED VORTEX
!
         elseif(elname(1:4).eq.'VOFO') then   
!
!     Core swirl ratio must be defined and positive
            if(prop(npropstart+4).le.0) then
               write(*,*)'*ERROR reading *FLUID SECTION:'
               write(*,*)'trying to define a FORCED VORTEX'
               write(*,*)'Core swirl ratio Kr is strictly positive'
               call exit(201)
            endif
            if(prop(npropstart+4).gt.1) then
               write(*,*)'*ERROR reading *FLUID SECTION:'
               write(*,*)'trying to define a FORCED VORTEX'
               write(*,*)'Core swirl ratio Kr cannot be greater than 1'
               call exit(201)
            endif
!
!     Rotation speed must be defined and positive
            if(prop(npropstart+5).le.0) then
               write(*,*)'*ERROR reading *FLUID SECTION:'
               write(*,*)'trying to define a FORCED VORTEX'
               write(*,*)'Rotation speed n is strictly positive'
            endif
         endif
!
!     Absolute/relative system
!
         if((elname(1:3).eq.'ATR').or.
     &        (elname(1:3).eq.'RTA')) then
            if(prop(npropstart+1).eq.0) then
               write(*,*)'*ERROR reading *FLUID SECTION:'
               write(*,*)'trying to define an element' 
               write(*,*)'TYPE= ABSOLUTE TO RELATIVE or'
               write(*,*)'TYPE= RELATIVE TO ABSOLUTE'
               write(*,*)'Rotation speed is strictly positive'
               call exit(201)
            elseif(prop(npropstart+2).eq.0) then
               if(prop(npropstart+3).eq.0) then
                  write(*,*)'*ERROR reading *FLUID SECTION:'
                  write(*,*)'trying to define an element' 
                  write(*,*)'TYPE= ABSOLUTE TO RELATIVE or'
                  write(*,*)'TYPE= RELATIVE TO ABSOLUTE'
                  write(*,*)'Tengential velocity is 0 but'
                  write(*,*)'no reference element has been provided'
                  call exit(201)
               endif
            elseif(prop(npropstart+3).ne.0) then
               if(prop(npropstart+2).ne.0) then
                  write(*,*)'*ERROR reading *FLUID SECTION:'
                  write(*,*)'trying to define an element' 
                  write(*,*)'TYPE= ABSOLUTE TO RELATIVE or'
                  write(*,*)'TYPE= RELATIVE TO ABSOLUTE'
                  write(*,*)'reference element has been provided'
                  write(*,*)'but tengential velocity is already defined'
                  call exit(201)
               endif 
            endif
         endif
      endif
!
!     assigning the elements of the set the appropriate material
!     and property pointer
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            if(lakon(ialset(j))(1:1).ne.'D') then
               write(*,*) '*ERROR reading *FLUID SECTION: element ',
     &            ialset(j),' is no fluid element.'
               call exit(201)
            endif
            lakon(ialset(j))(2:8)=elname(1:7)
!
!           gas type elements used for liquids are labeled with LP
!
            if((liquid.and.(lakon(ialset(j))(2:3).eq.'RE')).or.
     &           (liquid.and.(lakon(ialset(j))(2:3).eq.'OR'))) 
     &             lakon(ialset(j))(2:3)='LP'
!
!
             if(liquid.and.(lakon(ialset(j))(2:3).eq.'VO')) then
               lakon(ialset(j))(2:3)='LP'
               if(lakon(ialset(j))(4:5).eq.'FR') then 
                  lakon(ialset(j))(4:5)='VF'
               else if(lakon(ialset(j))(4:5).eq.'FO') then 
                  lakon(ialset(j))(4:5)='VS'
               endif
            endif
!     
            if(liquid.and.((lakon(ialset(j))(4:5).eq.'BG').or.
     &           (lakon(ialset(j))(4:5).eq.'BT').or.
     &           (lakon(ialset(j))(4:5).eq.'MA').or.
     &           (lakon(ialset(j))(4:5).eq.'MM').or.
     &           (lakon(ialset(j))(4:5).eq.'PA').or.
     &           (lakon(ialset(j))(4:5).eq.'PM').or.
     &           (lakon(ialset(j))(4:5).eq.'PN'))) then
               write(*,*) ''
               write(*,*) '*ERROR reading *FLUID SECTION: element ',k,
     &' is no valid incompressible orifice element.'
               write(*,*) ' Please change element type', 
     &' and/or definition.'
               write(*,*) ' Calculation stopped.'
               call exit(201)
            endif
!
            ielmat(1,ialset(j))=imaterial
            if(ndprop.gt.0) ielprop(ialset(j))=npropstart
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               if(lakon(k)(1:1).ne.'D') then
                  write(*,*) '*ERROR reading *FLUID SECTION: element ',
     &                 k,' is no fluid element.'
                  call exit(201)
               endif
               lakon(k)(2:8)=elname(1:7)
!
!     gas type elements used for liquids are labeled with LP
!     
               if((liquid.and.(lakon(ialset(j))(2:3).eq.'RE')).or.
     &              (liquid.and.(lakon(ialset(j))(2:3).eq.'OR'))) 
     &              lakon(ialset(j))(2:3)='LP'
!     
!     
               if(liquid.and.(lakon(ialset(j))(2:3).eq.'VO')) then
                  lakon(ialset(j))(2:3)='LP'
                  if(lakon(ialset(j))(4:5).eq.'FR') then 
                     lakon(ialset(j))(4:5)='VF'
                  else if(lakon(ialset(j))(4:5).eq.'FO') then 
                     lakon(ialset(j))(4:5)='VS'
                  endif
               endif
!     
               if(liquid.and.((lakon(ialset(j))(4:5).eq.'BG').or.
     &              (lakon(ialset(j))(4:5).eq.'BT').or.
     &              (lakon(ialset(j))(4:5).eq.'MA').or.
     &              (lakon(ialset(j))(4:5).eq.'MM').or.
     &              (lakon(ialset(j))(4:5).eq.'PA').or.
     &              (lakon(ialset(j))(4:5).eq.'PM').or.
     &              (lakon(ialset(j))(4:5).eq.'PN'))) then
                  write(*,*) ''
                  write(*,*)'*ERROR reading *FLUID SECTION: element ',k,
     &                 ' is no valid incompressible orifice element.'
                  write(*,*) ' Please change element type ',
     &'and/or definition.'
                  write(*,*) ' Calculation stopped.'
                  call exit(201)
               endif
               ielmat(1,k)=imaterial
               if(ndprop.gt.0) ielprop(k)=npropstart
            enddo
         endif
      enddo
!
      return
      end


