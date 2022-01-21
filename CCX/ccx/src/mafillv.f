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
      subroutine mafillv(nef,ipnei,neifa,neiel,vfa,xxn,area,
     &  auv,adv,jq,irow,nzs,bv,vel,cosa,umfa,xlet,xle,gradvfa,xxi,
     &  body,volume,compressible,ielfa,lakon,ifabou,nbody,neq,
     &  dtimef,velo,veloo,sel,xrlfa,gamma,xxj,nactdohinv,a1,
     &  a2,a3)
!
      implicit none
!
      character*2 one,two,three
      character*8 lakon(*)
!
      integer i,nef,jdof1,indexf,ipnei(*),j,ifa,iel,neifa(*),
     &  neiel(*),jdof2,jq(*),irow(*),nzs,iwall,compressible,ielfa(4,*),
     &  ipointer,ifabou(*),nbody,neq,k,indexb,numfaces,nactdohinv(*)
!
      real*8 flux,vfa(0:5,*),xxn(3,*),area(*),auv(*),adv(*),bv(neq,3),
     &  vel(nef,0:5),cosa(*),umfa(*),xlet(*),xle(*),coef,gradvfa(3,3,*),
     &  xxi(3,*),body(0:3,*),volume(*),coef2,dtimef,velo(nef,0:5),
     &  veloo(nef,0:5),rhovel,constant,sel(3,*),xrlfa(3,*),gamma(*),
     &  xxj(3,*),a1,a2,a3
!
cd      do i=1,nzs
cd         write(*,*) 'mafillv auv,irow',i,auv(i),irow(i)
cd      enddo
cd      do i=1,ne
cd         write(*,*) 'mafillv jq nactdoh',i,jq(i),nactdoh(i)
cd      enddo
      one='b1'
      two='b2'
      three='b3'
!
      do i=1,nef
         jdof1=i
         indexf=ipnei(i)
         if(lakon(i)(4:4).eq.'8') then
            numfaces=6
         elseif(lakon(i)(4:4).eq.'6') then
            numfaces=5
         else
            numfaces=4
         endif
         do j=1,numfaces
c            write(*,*) 'mafillv1 elem face ',i,j
!
!              convection
!
            indexf=indexf+1
            ifa=neifa(indexf)
            iel=neiel(indexf)
            if(iel.ne.0) jdof2=iel
            flux=(vfa(1,ifa)*xxn(1,indexf)+
     &           vfa(2,ifa)*xxn(2,indexf)+
     &           vfa(3,ifa)*xxn(3,indexf))
     &           *vfa(5,ifa)*area(ifa)
!
!              vfa sometimes undefined?
!
            if(flux.ge.0.d0) then
!
!                 outflowing flux
!
               call add_sm_fl_as(auv,adv,jq,irow,jdof1,jdof1,
     &              flux,nzs)
c               write(*,*) 'mafillv2 adv ',jdof1,flux
!centdiff
               bv(jdof1,1)=bv(jdof1,1)-(vfa(1,ifa)-vel(i,1))*flux
               bv(jdof1,2)=bv(jdof1,2)-(vfa(2,ifa)-vel(i,2))*flux
               bv(jdof1,3)=bv(jdof1,3)-(vfa(3,ifa)-vel(i,3))*flux
c                  write(*,*) 'mafillv3 ',-(vfa(1,ifa)-vel(i,1))*flux
c                  write(*,*) 'mafillv4 ',-(vfa(2,ifa)-vel(i,2))*flux
c                  write(*,*) 'mafillv5 ',-(vfa(3,ifa)-vel(i,3))*flux
!end centdiff
            else
               if(iel.gt.0) then
!
!                    incoming flux from neighboring element
!
                  call add_sm_fl_as(auv,adv,jq,irow,jdof1,jdof2,flux,
     &                 nzs)
c                  write(*,*) 'mafillv6 auv ',jdof1,jdof2,flux
!centdiff
                  bv(jdof1,1)=bv(jdof1,1)-(vfa(1,ifa)-vel(iel,1))*flux
                  bv(jdof1,2)=bv(jdof1,2)-(vfa(2,ifa)-vel(iel,2))*flux
                  bv(jdof1,3)=bv(jdof1,3)-(vfa(3,ifa)-vel(iel,3))*flux
c                  write(*,*) 'mafillv7 ',-(vfa(1,ifa)-vel(iel,1))*flux
c                  write(*,*) 'mafillv8 ',-(vfa(2,ifa)-vel(iel,2))*flux
c                  write(*,*) 'mafillv9 ',-(vfa(3,ifa)-vel(iel,3))*flux
!end centdiff
               else
!
!                    incoming flux through boundary
!
                  if(ielfa(2,ifa).lt.0) then
                     indexb=-ielfa(2,ifa)
                     if(((ifabou(indexb+1).ne.0).and.
     &                    (ifabou(indexb+2).ne.0).and.
     &                    (ifabou(indexb+3).ne.0)).or.
     &                    (dabs(flux).lt.1.d-10)) then
                        bv(jdof1,1)=bv(jdof1,1)-vfa(1,ifa)*flux
                        bv(jdof1,2)=bv(jdof1,2)-vfa(2,ifa)*flux
                        bv(jdof1,3)=bv(jdof1,3)-vfa(3,ifa)*flux
c                  write(*,*) 'mafillv10 ',-vfa(1,ifa)*flux
c                  write(*,*) 'mafillv11 ',-vfa(2,ifa)*flux
c                  write(*,*) 'mafillv12 ',-vfa(3,ifa)*flux
                     else
                        write(*,*) '*ERROR in mafillv: not all'
                        write(*,*) '       components of an incoming'
                        write(*,*) '       flux through face ',j
                        write(*,*)'       of element ',nactdohinv(i),
     &                        ' are given'
c                           write(*,*) indexb,flux
c                           write(*,*) ifabou(indexb+1)
c                           write(*,*) ifabou(indexb+2)
c                           write(*,*) ifabou(indexb+3)
c                           write(*,*) xxn(1,indexf),xxn(2,indexf),
c     &                              xxn(3,indexf)
c                           write(*,*) vfa(1,ifa),vfa(2,ifa),vfa(3,ifa)
cccc                        call exit(201)
                     endif
                  else
                     write(*,*) '*ERROR in mafillv: not all'
                     write(*,*) '       components of an incoming'
                     write(*,*) '       flux through face ',j
                     write(*,*)'       of element ',nactdohinv(i),
     &                     ' are given'
cccc                     call exit(201)
                  endif
               endif
            endif
!
!              diffusion
!
            if(iel.ne.0) then
!
!                 neighboring element
!
               coef=umfa(ifa)*area(ifa)/xlet(indexf)
               call add_sm_fl_as(auv,adv,jq,irow,jdof1,jdof1,
     &              coef,nzs)
c               write(*,*) 'mafillv13 adv ',jdof1,coef
               call add_sm_fl_as(auv,adv,jq,irow,jdof1,jdof2,
     &              -coef,nzs)
c               write(*,*) 'mafillv14 auv ',jdof1,jdof2,-coef
!
!                 correction for non-orthogonal grid
!
               bv(jdof1,1)=bv(jdof1,1)+umfa(ifa)*area(ifa)*
     &              (gradvfa(1,1,ifa)*(xxn(1,indexf)-xxj(1,indexf))+
     &              gradvfa(1,2,ifa)*(xxn(2,indexf)-xxj(2,indexf))+
     &              gradvfa(1,3,ifa)*(xxn(3,indexf)-xxj(3,indexf)))
               bv(jdof1,2)=bv(jdof1,2)+umfa(ifa)*area(ifa)*
     &              (gradvfa(2,1,ifa)*(xxn(1,indexf)-xxj(1,indexf))+
     &              gradvfa(2,2,ifa)*(xxn(2,indexf)-xxj(2,indexf))+
     &              gradvfa(2,3,ifa)*(xxn(3,indexf)-xxj(3,indexf)))
               bv(jdof1,3)=bv(jdof1,3)+umfa(ifa)*area(ifa)*
     &              (gradvfa(3,1,ifa)*(xxn(1,indexf)-xxj(1,indexf))+
     &              gradvfa(3,2,ifa)*(xxn(2,indexf)-xxj(2,indexf))+
     &              gradvfa(3,3,ifa)*(xxn(3,indexf)-xxj(3,indexf)))
c               write(*,*) 'mafillv15 ',jdof1,
c     &         umfa(ifa)*area(ifa)*
c     &              (gradvfa(1,1,ifa)*(xxn(1,indexf)-xxj(1,indexf))+
c     &              gradvfa(1,2,ifa)*(xxn(2,indexf)-xxj(2,indexf))+
c     &              gradvfa(1,3,ifa)*(xxn(3,indexf)-xxj(3,indexf))),
c     &         umfa(ifa)*area(ifa)*
c     &              (gradvfa(2,1,ifa)*(xxn(1,indexf)-xxj(1,indexf))+
c     &              gradvfa(2,2,ifa)*(xxn(2,indexf)-xxj(2,indexf))+
c     &              gradvfa(2,3,ifa)*(xxn(3,indexf)-xxj(3,indexf))),
c     &         umfa(ifa)*area(ifa)*
c     &              (gradvfa(3,1,ifa)*(xxn(1,indexf)-xxj(1,indexf))+
c     &              gradvfa(3,2,ifa)*(xxn(2,indexf)-xxj(2,indexf))+
c     &              gradvfa(3,3,ifa)*(xxn(3,indexf)-xxj(3,indexf)))
c                  write(*,*) 'mafillv4 ',i,j,jdof1,bv(jdof1,1)
            else
!
!                 boundary; check whether wall (specified by user),
!                           outlet (no velocity boundary conditions or
!                           none of those
!
               iwall=0
               ipointer=abs(ielfa(2,ifa))
               if(ipointer.gt.0) then
                  iwall=ifabou(ipointer+5)
               endif
               if(iwall.eq.0) then
!
!                    external face, but no wall
!
                  if((ifabou(ipointer+1).ne.0).or.
     &                 (ifabou(ipointer+2).ne.0).or.
     &                 (ifabou(ipointer+3).ne.0)) then
!
!                       no outlet: face velocity fixed
!
                     coef=umfa(ifa)*area(ifa)/xle(indexf)
                     call add_sm_fl_as(auv,adv,jq,irow,jdof1,jdof1,
     &                    coef,nzs)
c                     write(*,*) 'mafillv16 adv ',jdof1,coef
                     bv(jdof1,1)=bv(jdof1,1)+coef*vfa(1,ifa)
                     bv(jdof1,2)=bv(jdof1,2)+coef*vfa(2,ifa)
                     bv(jdof1,3)=bv(jdof1,3)+coef*vfa(3,ifa)
c                     write(*,*) 'mafillv17 ',coef*vfa(1,ifa)
c                     write(*,*) 'mafillv18 ',coef*vfa(2,ifa)
c                     write(*,*) 'mafillv19 ',coef*vfa(3,ifa)
c                  write(*,*) 'mafillv5 ',i,j,jdof1,bv(jdof1,1)
                  else
!
!                       outlet: no diffusion
!
                  endif
!
!                    correction for non-orthogonal grid
!
                  bv(jdof1,1)=bv(jdof1,1)+umfa(ifa)*area(ifa)*
     &                 (gradvfa(1,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &                 gradvfa(1,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &                 gradvfa(1,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
                  bv(jdof1,2)=bv(jdof1,2)+umfa(ifa)*area(ifa)*
     &                 (gradvfa(2,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &                 gradvfa(2,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &                 gradvfa(2,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
                  bv(jdof1,3)=bv(jdof1,3)+umfa(ifa)*area(ifa)*
     &                 (gradvfa(3,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &                 gradvfa(3,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &                 gradvfa(3,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
c                  write(*,*) 'mafillv20',umfa(ifa)*area(ifa)*
c     &                 (gradvfa(1,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
c     &                 gradvfa(1,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
c     &                 gradvfa(1,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
c                  write(*,*) 'mafillv21',umfa(ifa)*area(ifa)*
c     &                 (gradvfa(2,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
c     &                 gradvfa(2,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
c     &                 gradvfa(2,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
c                  write(*,*) 'mafillv22',umfa(ifa)*area(ifa)*
c     &                 (gradvfa(3,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
c     &                 gradvfa(3,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
c     &                 gradvfa(3,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
c                  write(*,*) 'mafillv6 ',i,j,jdof1,bv(jdof1,1)
c                  write(*,*) umfa(ifa),area(ifa)
c                  write(*,*)gradvfa(1,1,ifa),xxn(1,indexf),xxi(1,indexf)
c                  write(*,*)gradvfa(1,2,ifa),xxn(2,indexf),xxi(2,indexf)
c                  write(*,*)gradvfa(1,3,ifa),xxn(3,indexf),xxi(3,indexf)
               else
!     
!                    wall
!     
                  coef=umfa(ifa)*area(ifa)/(xle(indexf)*cosa(indexf))
                  call add_sm_fl_as(auv,adv,jq,irow,jdof1,jdof1,coef,
     &                 nzs)
c                  write(*,*) 'mafillv23 adv ',jdof1,coef
!
!                    correction for non-orthogonal grid and nonzero
!                    wall velocity
!
c                     coef2=((vel(1,i))*xxn(1,indexf)+
c     &                     (vel(2,i))*xxn(2,indexf)+
c     &                     (vel(3,i))*xxn(3,indexf))*coef
c                     bv(jdof1,1)=bv(jdof1,1)+coef*vfa(1,ifa)+
c     &                                     coef2*xxn(1,indexf)
c                     bv(jdof1,2)=bv(jdof1,2)+coef*vfa(2,ifa)+
c     &                                     coef2*xxn(2,indexf)
c                     bv(jdof1,3)=bv(jdof1,3)+coef*vfa(3,ifa)+
c     &                                     coef2*xxn(3,indexf)
                  coef2=((vel(i,1)-vfa(1,ifa))*xxn(1,indexf)+
     &                 (vel(i,2)-vfa(2,ifa))*xxn(2,indexf)+
     &                 (vel(i,3)-vfa(3,ifa))*xxn(3,indexf))*coef
                  bv(jdof1,1)=bv(jdof1,1)+coef*vfa(1,ifa)+
     &                 coef2*xxn(1,indexf)
                  bv(jdof1,2)=bv(jdof1,2)+coef*vfa(2,ifa)+
     &                 coef2*xxn(2,indexf)
                  bv(jdof1,3)=bv(jdof1,3)+coef*vfa(3,ifa)+
     &                 coef2*xxn(3,indexf)
c                  write(*,*) 'mafillv24 ',coef*vfa(1,ifa)+
c     &                 coef2*xxn(1,indexf)
c                  write(*,*) 'mafillv25 ',coef*vfa(2,ifa)+
c     &                 coef2*xxn(2,indexf)
c                  write(*,*) 'mafillv26 ',coef*vfa(3,ifa)+
c     &                 coef2*xxn(3,indexf)
c                  write(*,*) 'mafillv7 ',i,j,jdof1,bv(jdof1,1)
               endif
            endif
!     
!     compressible
!     
c               if(compressible.eq.1) then
c                  bv(jdof1,1)=bv(jdof1,1)+umfa(ifa)*area(ifa)*
c     &                 (gradvfa(1,1,ifa)*xxn(1,indexf)+
c     &                 gradvfa(2,1,ifa)*xxn(2,indexf)+
c     &                 gradvfa(3,1,ifa)*xxn(3,indexf)-
c     &                 2.d0*(gradvfa(1,1,ifa)+gradvfa(2,2,ifa)+
c     &                 gradvfa(3,3,ifa))*xxn(1,indexf)/3.d0)
c                  bv(jdof1,2)=bv(jdof1,2)+umfa(ifa)*area(ifa)*
c     &                 (gradvfa(1,2,ifa)*xxn(1,indexf)+
c     &                 gradvfa(2,2,ifa)*xxn(2,indexf)+
c     &                 gradvfa(3,2,ifa)*xxn(3,indexf)-
c     &                 2.d0*(gradvfa(1,1,ifa)+gradvfa(2,2,ifa)+
c     &                 gradvfa(3,3,ifa))*xxn(2,indexf)/3.d0)
c                  bv(jdof1,3)=bv(jdof1,3)+umfa(ifa)*area(ifa)*
c     &                 (gradvfa(1,3,ifa)*xxn(1,indexf)+
c     &                 gradvfa(2,3,ifa)*xxn(2,indexf)+
c     &                 gradvfa(3,3,ifa)*xxn(3,indexf)-
c     &                 2.d0*(gradvfa(1,1,ifa)+gradvfa(2,2,ifa)+
c     &                 gradvfa(3,3,ifa))*xxn(3,indexf)/3.d0)
c                  write(*,*) 'mafillv8 ',jdof1,j,bv(jdof1,3)
c               endif
!     
!     pressure
!     
c               bv(jdof1,1)=bv(jdof1,1)-vfa(4,ifa)*xxn(1,indexf)*area(ifa)
ccd         write(*,*) 'press ',i,j,one,-vfa(4,ifa)*xxn(1,indexf)*area(ifa)
c               bv(jdof1,2)=bv(jdof1,2)-vfa(4,ifa)*xxn(2,indexf)*area(ifa)
cc        write(*,*) 'press ',i,j,two,-vfa(4,ifa)*xxn(2,indexf)*area(ifa)
c               bv(jdof1,3)=bv(jdof1,3)-vfa(4,ifa)*xxn(3,indexf)*area(ifa)
cc       write(*,*) 'press ',i,j,three,-vfa(4,ifa)*xxn(3,indexf)*area(ifa)
         enddo
!     
!           body force
!     
         rhovel=vel(i,5)*volume(i)
!
         if(nbody.gt.0) then
            bv(jdof1,1)=bv(jdof1,1)+rhovel*body(1,i)
            bv(jdof1,2)=bv(jdof1,2)+rhovel*body(2,i)
            bv(jdof1,3)=bv(jdof1,3)+rhovel*body(3,i)
c                  write(*,*) 'mafillv8 ',i,j,jdof1,bv(jdof1,1)
         endif
!
!           transient term
!
         constant=rhovel/dtimef
         bv(jdof1,1)=bv(jdof1,1)-(a2*velo(i,1)+a3*veloo(i,1))*constant
         bv(jdof1,2)=bv(jdof1,2)-(a2*velo(i,2)+a3*veloo(i,2))*constant
         bv(jdof1,3)=bv(jdof1,3)-(a2*velo(i,3)+a3*veloo(i,3))*constant
c         write(*,*) 'mafillv27',(4.d0*velo(i,1)-veloo(i,1))*constant
c         write(*,*) 'mafillv28',(4.d0*velo(i,2)-veloo(i,2))*constant
c         write(*,*) 'mafillv29',(4.d0*velo(i,3)-veloo(i,3))*constant
c                  write(*,*) 'mafillv9 ',i,j,jdof1,bv(jdof1,1)
         constant=a1*constant
         call add_sm_fl(auv,adv,jq,irow,jdof1,jdof1,constant,nzs)
c         write(*,*) 'mafillv30 adv ',jdof1,constant
!
!           copying b into sel (rhs without pressure)
!
         do j=1,3
            sel(j,jdof1)=bv(jdof1,j)
         enddo
!
!           pressure contribution to b
!
         indexf=ipnei(i)
         do j=1,numfaces
            indexf=indexf+1
            ifa=neifa(indexf)
            bv(jdof1,1)=bv(jdof1,1)
     &           -vfa(4,ifa)*xxn(1,indexf)*area(ifa)
            bv(jdof1,2)=bv(jdof1,2)
     &           -vfa(4,ifa)*xxn(2,indexf)*area(ifa)
            bv(jdof1,3)=bv(jdof1,3)
     &           -vfa(4,ifa)*xxn(3,indexf)*area(ifa)
c           write(*,*) 'mafillv31 ',j,-vfa(4,ifa)*xxn(1,indexf)*area(ifa)
c           write(*,*) 'mafillv32 ',j,-vfa(4,ifa)*xxn(2,indexf)*area(ifa)
c           write(*,*) 'mafillv33 ',j,-vfa(4,ifa)*xxn(3,indexf)*area(ifa)
c                  write(*,*) 'mafillv11 ',i,j,jdof1,bv(jdof1,1)
         enddo
!            
      enddo
!
c      do i=1,nzs
c         write(*,*) 'mafillv auv,irow',i,auv(i),auv(i+nzs)
c      enddo
c      do i=1,nef
c         write(*,*) 'mafillv b adv',i,bv(i,1),adv(i)
c      enddo
!     
      return
      end
