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
      subroutine printoutelem(prlab,ipkon,lakon,kon,co,
     &     ener,mi,ii,nelem,energytot,volumetot,enerkintot,nkin,ne,
     &     stx,nodes,thicke,ielmat,ielem,iface,mortar)
!
!     stores whole element results for element "nelem" in the .dat file
!
      implicit none
!
      character*6 prlab(*)
      character*8 lakon(*)
!
      integer ipkon(*),nelem,ii,kon(*),mi(*),nope,indexe,i,j,k,
     &  konl(20),iface,mortar,ielem,
     &  mint3d,jj,nener,iflag,nkin,ne,nodes,ki,kl,ilayer,nlayer,kk,
     &  nopes,ielmat(mi(3),*),mint2d
!
      real*8 ener(mi(1),*),energytot,volumetot,energy,volume,co(3,*),
     &  xl(3,20),xi,et,ze,xsj,shp(4,20),weight,enerkintot,enerkin,
     &  stx(6,mi(1),*),a,gs(8,4),dlayer(4),tlayer(4),thickness,
     &  thicke(mi(3),*),xlayer(mi(3),4),shp2(7,8),xs2(3,7),xsj2(3),
     &  xl2(3,8)
!
      include "gauss.f"
!
      data iflag /2/
!
      if(ipkon(nelem).lt.0) return
      indexe=ipkon(nelem)
c      if(ipkon(nelem).eq.-1) then
c!
c!        linear element corresponding to a remeshing of a quadratic
c!        element adjacent to a contact surface
c!
c         return
c      elseif(ipkon(nelem).lt.-1) then
c!
c!        element is quadratic and adjacent to a contact surface
c!        -> it has been remeshed; the first node of the topology has
c!           been replaced by a pointer to the first linear element
c!           of the remeshing, the first node of which is identical to
c!           the first node of the original quadratic element
c!
c         indexe=-ipkon(nelem)-2
c         ielemremesh=kon(indexe+1)
c         kon(indexe+1)=kon(ipkon(ielemremesh)+1)
c      else
c         indexe=ipkon(nelem)
c      endif
!
      if((prlab(ii)(1:4).eq.'ELSE').or.(prlab(ii)(1:4).eq.'CELS')) then
         nener=1
      else
         nener=0
      endif
!
      if(lakon(nelem)(1:5).eq.'C3D8I') then
         nope=11
      elseif(lakon(nelem)(4:4).eq.'2') then
         nope=20
      elseif(lakon(nelem)(4:4).eq.'8') then
         nope=8
      elseif(lakon(nelem)(4:5).eq.'10') then
         nope=10
      elseif(lakon(nelem)(4:4).eq.'4') then
         nope=4
      elseif(lakon(nelem)(4:5).eq.'15') then
         nope=15
      elseif(lakon(nelem)(4:5).eq.'6') then
         nope=6
      else
         nope=0
      endif
!
!        composite materials
!
      if(lakon(nelem)(7:8).eq.'LC') then
!
!        determining the number of layers
!
         nlayer=0
         do k=1,mi(3)
            if(ielmat(k,nelem).ne.0) then
               nlayer=nlayer+1
            endif
         enddo
         mint2d=4
         nopes=8
!
!        determining the layer thickness and global thickness
!        at the shell integration points
!
         iflag=1
         indexe=ipkon(nelem)
         do kk=1,mint2d
            xi=gauss3d2(1,kk)
            et=gauss3d2(2,kk)
            call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            tlayer(kk)=0.d0
            do i=1,nlayer
               thickness=0.d0
               do j=1,nopes
                  thickness=thickness+thicke(i,indexe+j)*shp2(4,j)
               enddo
               tlayer(kk)=tlayer(kk)+thickness
               xlayer(i,kk)=thickness
            enddo
         enddo
         iflag=2
!
         ilayer=0
         do i=1,4
            dlayer(i)=0.d0
         enddo
!     
      endif
!
      do j=1,nope
         konl(j)=kon(indexe+j)
         do k=1,3
            xl(k,j)=co(k,konl(j))
         enddo
      enddo
!
      energy=0.d0
      volume=0.d0
      enerkin=0.d0
!
      if(lakon(nelem)(4:5).eq.'8R') then
         mint3d=1
      elseif((lakon(nelem)(4:4).eq.'8').or.
     &        (lakon(nelem)(4:6).eq.'20R')) then
         if(lakon(nelem)(7:8).eq.'LC') then
            mint3d=8*nlayer
         else
            mint3d=8
         endif
      elseif(lakon(nelem)(4:4).eq.'2') then
         mint3d=27
      elseif(lakon(nelem)(4:5).eq.'10') then
         mint3d=4
      elseif(lakon(nelem)(4:4).eq.'4') then
         mint3d=1
      elseif(lakon(nelem)(4:5).eq.'15') then
         mint3d=9
      elseif(lakon(nelem)(4:5).eq.'6') then
         mint3d=2
      else
         if(nener.eq.1)then
            energy=ener(1,nelem)
         endif
         mint3d=0
      endif
!
      do jj=1,mint3d
         if(lakon(nelem)(4:5).eq.'8R') then
            xi=gauss3d1(1,jj)
            et=gauss3d1(2,jj)
            ze=gauss3d1(3,jj)
            weight=weight3d1(jj)
         elseif((lakon(nelem)(4:4).eq.'8').or.
     &           (lakon(nelem)(4:6).eq.'20R'))
     &           then
            if(lakon(nelem)(7:8).ne.'LC') then
               xi=gauss3d2(1,jj)
               et=gauss3d2(2,jj)
               ze=gauss3d2(3,jj)
               weight=weight3d2(jj)
            else
               kl=mod(jj,8)
               if(kl.eq.0) kl=8
!     
               xi=gauss3d2(1,kl)
               et=gauss3d2(2,kl)
               ze=gauss3d2(3,kl)
               weight=weight3d2(kl)
!     
               ki=mod(jj,4)
               if(ki.eq.0) ki=4
!     
               if(kl.eq.1) then
                  ilayer=ilayer+1
                  if(ilayer.gt.1) then
                     do i=1,4
                        dlayer(i)=dlayer(i)+xlayer(ilayer-1,i)
                     enddo
                  endif
               endif
               ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &              tlayer(ki)-1.d0
               weight=weight*xlayer(ilayer,ki)/tlayer(ki)
            endif
         elseif(lakon(nelem)(4:4).eq.'2') then
            xi=gauss3d3(1,jj)
            et=gauss3d3(2,jj)
            ze=gauss3d3(3,jj)
            weight=weight3d3(jj)
         elseif(lakon(nelem)(4:5).eq.'10') then
            xi=gauss3d5(1,jj)
            et=gauss3d5(2,jj)
            ze=gauss3d5(3,jj)
            weight=weight3d5(jj)
         elseif(lakon(nelem)(4:4).eq.'4') then
            xi=gauss3d4(1,jj)
            et=gauss3d4(2,jj)
            ze=gauss3d4(3,jj)
            weight=weight3d4(jj)
         elseif(lakon(nelem)(4:5).eq.'15') then
            xi=gauss3d8(1,jj)
            et=gauss3d8(2,jj)
            ze=gauss3d8(3,jj)
            weight=weight3d8(jj)
         else
            xi=gauss3d7(1,jj)
            et=gauss3d7(2,jj)
            ze=gauss3d7(3,jj)
            weight=weight3d7(jj)
         endif
!
         if(lakon(nelem)(1:5).eq.'C3D8R') then
            call shape8hr(xl,xsj,shp,gs,a)
         elseif(lakon(nelem)(1:5).eq.'C3D8I') then
            call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.20) then
            call shape20h(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.8) then
            call shape8h(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.10) then
            call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.4) then
            call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
         elseif(nope.eq.15) then
            call shape15w(xi,et,ze,xl,xsj,shp,iflag)
         else
            call shape6w(xi,et,ze,xl,xsj,shp,iflag)
         endif
!
         if(nener.eq.1) energy=energy+weight*xsj*ener(jj,nelem)
         if(nkin.eq.1) enerkin=enerkin+weight*xsj*ener(jj,nelem+ne)
         volume=volume+weight*xsj
      enddo
!
      volumetot=volumetot+volume
      if(nener.eq.1) energytot=energytot+energy
      if(nkin.eq.1) enerkintot=enerkintot+enerkin
!     
!     writing to file
!     
      if((prlab(ii)(1:5).eq.'ELSE ').or.
     &     (prlab(ii)(1:5).eq.'ELSET')) then
         write(5,'(i10,1p,1x,e13.6)') nelem,energy
      elseif((prlab(ii)(1:5).eq.'CELS ').or.
     &        (prlab(ii)(1:5).eq.'CELST')) then
         if(mortar.eq.0) then
            write(5,'(i10,1p,1x,e13.6)') nodes,energy
         elseif(mortar.eq.1) then
            write(5,'(i10,1x,i10,1p,1x,e13.6)') ielem,iface,energy
         endif
      elseif((prlab(ii)(1:5).eq.'CDIS ').or.
     &        (prlab(ii)(1:5).eq.'CDIST')) then
         if(mortar.eq.0) then
            write(5,'(i10,1p,1x,e13.6,1p,1x,e13.6,1p,1x,e13.6)') nodes,
     &           stx(1,1,nelem),stx(2,1,nelem),stx(3,1,nelem)
         elseif(mortar.eq.1) then
            write(5,'(i10,1x,i10,1p,1x,e13.6,1p,1x,e13.6,1p,1x,e13.6)') 
     &           ielem,iface,
     &           stx(1,1,nelem),stx(2,1,nelem),stx(3,1,nelem)
         endif
      elseif((prlab(ii)(1:5).eq.'CSTR ').or.
     &        (prlab(ii)(1:5).eq.'CSTRT')) then
         if(mortar.eq.0) then
            write(5,'(i10,1p,1x,e13.6,1p,1x,e13.6,1p,1x,e13.6)') nodes,
     &           stx(4,1,nelem),stx(5,1,nelem),stx(6,1,nelem)
         elseif(mortar.eq.1) then
            write(5,'(i10,1x,i10,1p,1x,e13.6,1p,1x,e13.6,1p,1x,e13.6)') 
     &           ielem,iface,
     &           stx(4,1,nelem),stx(5,1,nelem),stx(6,1,nelem)
         endif
      elseif((prlab(ii)(1:5).eq.'EVOL ').or.
     &        (prlab(ii)(1:5).eq.'EVOLT')) then
         write(5,'(i10,1p,1x,e13.6)') nelem,volume
      elseif((prlab(ii)(1:5).eq.'ELKE ').or.
     &        (prlab(ii)(1:5).eq.'ELKET')) then
         write(5,'(i10,1p,1x,e13.6)') nelem,enerkin
      endif
c!
c!     restoring the topology of a quadratic element which has been
c!     remeshed because of its adjacency to a contact surface
c!    
c      if(ipkon(nelem).lt.-1) then
c         kon(indexe+1)=ielemremesh
c      endif
!
      return
      end
      
