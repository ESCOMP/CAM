!-----------------------------------------------------------------------
!BOP
! !ROUTINE: pkez --- Calculate solution to hydrostatic equation
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine pkez(nx, im, km, jfirst, jlast, kfirst, klast,    &
                      ifirst, ilast, pe, pk, cap3v, ks, peln, pkz, eta, high_alt)
!****6***0*********0*********0*********0*********0*********0**********72
!
! !USES:
      use shr_kind_mod, only: r8 => shr_kind_r8

      implicit none

!
! This routine may be called assuming either yz or xy decompositions.
! For xy decomposition, the effective "nx" is 1.
!

! !INPUT PARAMETERS:
      integer, intent(in) :: nx                          ! SMP decomposition in x
      integer, intent(in) :: im, km                      ! Dimensions
      integer, intent(in) :: jfirst, jlast               ! Latitude strip
      integer, intent(in) :: kfirst, klast               ! Vertical strip
      integer, intent(in) :: ifirst, ilast               ! Longitude strip
      real (r8), intent(in) ::  pe(ifirst:ilast, kfirst:klast+1, jfirst:jlast)    ! Edge pressure
      integer, intent(in) :: ks
      logical, intent(in) :: eta     ! Is on ETA coordinate?
                      ! True:  input pe    ; output pk, pkz, peln
                      ! False: input pe, pk; output     pkz, peln
      real (r8), intent(in) :: cap3v(ifirst:ilast,jfirst:jlast,km)
      logical, intent(in) :: high_alt

! !INPUT/OUTPUT PARAMETERS:
      real (r8), intent(inout) :: pk(ifirst:ilast,jfirst:jlast,kfirst:klast+1)

! !OUTPUT PARAMETERS
      real (r8), intent(out) :: pkz(ifirst:ilast,jfirst:jlast,kfirst:klast)
      real (r8), intent(out) :: peln(ifirst:ilast, kfirst:klast+1, jfirst:jlast)   ! log pressure (pe) at layer edges

! !DESCRIPTION:
!
!
! !CALLED FROM:
!     te_map and fvccm3
!
! !REVISION HISTORY:
!
!     WS  99.05.19 : Removed fvcore.h
!     WS  99.07.27 : Limited region to jfirst:jlast
!     WS  99.10.22 : Deleted cp as argument (was not used)
!     WS  99.11.05 : Documentation; pruning of arguments
!     SJL 00.01.02 : SMP decomposition in i
!     AAM 00.08.10 : Add kfirst:klast
!     AAM 01.06.27 : Add ifirst:ilast
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local
      real (r8) pk2(ifirst:ilast, kfirst:klast+1)
      real (r8) pek
      real (r8) lnp
      real (r8) lnpk
      real (r8) cap3vi(ifirst:ilast,jfirst:jlast,km+1)
      real (r8) pkln(ifirst:ilast,km+1,jfirst:jlast)  ! log pk at layer edges
      integer i, j, k, itot, nxu
      integer ixj, jp, it, i1, i2

      itot = ilast - ifirst + 1
! Use smaller block sizes only if operating on full i domain
      nxu = 1
      if (itot .eq. im) nxu = nx

      it = itot / nxu
      jp = nxu * ( jlast - jfirst + 1 )

      if ( eta ) then
         if (high_alt) then
            !$omp parallel do private(i,j,k)
            do k=2,km
               do j=jfirst,jlast
                  do i=ifirst,ilast
                     cap3vi(i,j,k) = 0.5_r8*(cap3v(i,j,k-1)+cap3v(i,j,k))
                  enddo
               enddo
            enddo
            cap3vi(:,:,1) = 1.5_r8 * cap3v(:,:,1) - 0.5_r8 * cap3v(:,:,2)
            cap3vi(:,:,km+1) = 1.5_r8 * cap3v(:,:,km) - 0.5_r8 * cap3v(:,:,km-1)
         else
            cap3vi(:,:,:) =  cap3v(ifirst,jfirst,1)
         endif
      endif

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(ixj, i1, i2, i, j, k, pek, lnp, pk2)

! WS 99.07.27 : Limited region to jfirst:jlast

      do 1000 ixj=1,jp

         j  = jfirst + (ixj-1) / nxu
         i1 = ifirst + it * mod(ixj-1, nxu)
         i2 = i1 + it - 1

        if ( eta ) then

! <<<<<<<<<<< Eta cordinate Coordinate  >>>>>>>>>>>>>>>>>>>
          if (kfirst .eq. 1) then
            pek =     pe(i1,1,j)**cap3vi(i1,j,1)
            lnp = log(pe(i1,1,j))
            lnpk = log(pek)
            do i=i1,i2
               pk2(i,1)   = pek
              peln(i,1,j) = lnp
              pkln(i,1,j) = lnpk
            enddo
          endif

          if(ks .ne. 0) then
            do k=max(2,kfirst), min(ks+1,klast+1)
              pek = pe(i1,k,j)**cap3vi(i1,j,k)
              lnp = log(pe(i1,k,j))
              lnpk = log(pek)
              do i=i1,i2
                pk2(i,k)   = pek
                peln(i,k,j) =  lnp
                pkln(i,k,j) = lnpk
              enddo
            enddo

            do k=kfirst, min(ks,klast)
              pek = (       pk2(i1,k+1)   - pk2(i1,k))   /     &
                    (pkln(i1,k+1,j) - pkln(i1,k,j))
              do i=i1,i2
                 pkz(i,j,k) = pek
              enddo
            enddo
          endif

          do k=max(ks+2,kfirst), klast+1
#if !defined( VECTOR_MATH )
            do i=i1,i2
               pk2(i,k) = pe(i,k,j)**cap3vi(i,j,k)
            enddo
#else
            call vlog(pk2(i1,k), pe(i1,k,j), it)
            do i=i1,i2
               pk2(i,k) = cap3vi(i,j,k) * pk2(i,k)
            enddo
            call vexp(pk2(i1,k), pk2(i1,k), it)
#endif
          enddo

          do k=max(ks+2,kfirst), klast+1
            do i=i1,i2
               peln(i,k,j) =  log(pe(i,k,j))
               pkln(i,k,j) =  log(pk2(i,k))
            enddo
          enddo

          do k=max(ks+1,kfirst), klast
            do i=i1,i2
               pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k)) /         &
                            (pkln(i,k+1,j) - pkln(i,k,j))
            enddo
          enddo

          do k=kfirst, klast+1
            do i=i1,i2
               pk(i,j,k) = pk2(i,k)
            enddo
          enddo

        else

! <<<<<<<<<<< General Coordinate  >>>>>>>>>>>>>>>>>>>

          if (kfirst .eq. 1) then
            lnp = log(pe(i1,1,j)) ! do log only one time at top -- assumes pe is constant at top

            do i=i1,i2
               peln(i,1,j) = lnp
            enddo
          endif

          do k=max(2,kfirst), klast+1
             do i=i1,i2
                peln(i,k,j) = log(pe(i,k,j))
             enddo
          enddo
          do k=kfirst, klast+1 ! variable pk at the top interface --> 
             do i=i1,i2
                pk2(i,k) = pk(i,j,k)
             enddo
          enddo
          if (high_alt) then
             do k=kfirst, klast+1 ! variable pk at the top interface --> 
                do i=i1,i2
                   pkln(i,k,j) = log(pk(i,j,k))
                enddo
             enddo
          endif

          if (high_alt) then
             do k=kfirst, klast
                do i=i1,i2
                   pkz(i,j,k) = ( pk2(i,k+1) - pk2(i,k) )  /    &
                        (pkln(i,k+1,j) - pkln(i,k,j))
                enddo
             enddo
          else
             do k=kfirst, klast
                do i=i1,i2
                   pkz(i,j,k) = ( pk2(i,k+1) - pk2(i,k) )  /    &
                        (cap3v(i,j,k)*(peln(i,k+1,j) - peln(i,k,j)))
                enddo
             enddo
          endif

       endif

1000  continue

      return
!EOC
      end
!-----------------------------------------------------------------------
