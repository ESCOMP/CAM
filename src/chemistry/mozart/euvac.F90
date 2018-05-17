!-----------------------------------------------------------------------
! An emperical model which uses F10.7 index to provide EUV solar spectrum
!  Richards, P.G., et al, EUVAC, A solar EUV flux model for aeronomic 
!  calculations, J.Geophys. Res., 99, 8981 - 8992, 1994
!-----------------------------------------------------------------------
      module euvac

      use shr_kind_mod,     only : r8 => shr_kind_r8
      use cam_logfile,      only : iulog
      implicit none

      private
      public :: euvac_init
      public :: euvac_set_etf
      public :: euvac_etf

      save

      integer               :: nbins
      real(r8), allocatable :: refmin(:)
      real(r8), allocatable :: afac(:)
      real(r8), target, allocatable :: euvac_etf(:)

      logical :: euvac_on 

      contains

      subroutine euvac_init (euvac_file)
!---------------------------------------------------------------
!       ... initialize euvac etf module
!---------------------------------------------------------------

      use cam_pio_utils,  only : cam_pio_openfile
      use pio,            only : pio_nowrite, pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, &
                                 pio_get_var, file_desc_t, pio_closefile
      use error_messages, only : alloc_err
      use ioFileMod,      only : getfil
      implicit none

      character(len=*), intent(in) :: euvac_file

!---------------------------------------------------------------
!       ... local variables
!---------------------------------------------------------------
      type(file_desc_t)  :: ncid
      integer  :: ierr
      integer  :: dimid
      integer  :: varid
      integer  :: astat
      character(len=256) :: locfn

      euvac_on = len_trim(euvac_file)>0 .and. euvac_file.ne.'NONE'
      if (.not.euvac_on) return

!-----------------------------------------------------------------------
!       ... readin the etf data
!-----------------------------------------------------------------------
      call getfil( euvac_file, locfn, 0 )
      call cam_pio_openfile (ncid, trim(locfn), PIO_NOWRITE)
!-----------------------------------------------------------------------
!       ... get number of bins
!-----------------------------------------------------------------------
      ierr = pio_inq_dimid( ncid, 'bin', dimid )
      ierr = pio_inq_dimlen( ncid, dimid, nbins )

!-----------------------------------------------------------------------
!       ... allocate primary arrays
!-----------------------------------------------------------------------
      allocate( refmin(nbins), afac(nbins), euvac_etf(nbins), stat=astat )
      if( astat /= 0 ) then
         call alloc_err( astat, 'euvac_init', 'wc ... euvac_etf', nbins )
      end if
!-----------------------------------------------------------------------
!       ... read primary arrays
!-----------------------------------------------------------------------
      ierr = pio_inq_varid( ncid, 'REFMIN', varid )
      ierr = pio_get_var( ncid, varid, refmin )
      ierr = pio_inq_varid( ncid, 'AFAC', varid )
      ierr = pio_get_var( ncid, varid, afac )
      
      call pio_closefile( ncid )

      end subroutine euvac_init

      subroutine euvac_set_etf( f107, f107a )
!---------------------------------------------------------------
!       ... set euvac etf
!---------------------------------------------------------------

      use spmd_utils,     only : masterproc

      implicit none

!---------------------------------------------------------------
!       ... dummy arguments
!---------------------------------------------------------------
      real(r8), intent(in) :: f107
      real(r8), intent(in) :: f107a

!---------------------------------------------------------------
!       ... local variables
!---------------------------------------------------------------
      real(r8), parameter :: factor = 80._r8
      real(r8) :: pindex

      if (.not.euvac_on) return

      pindex = .5_r8*(f107 + f107a) - factor
      euvac_etf(:) = refmin(:) * max( .8_r8,(1._r8 + afac(:)*pindex) )

      if( masterproc ) then
         write(iulog,*) ' '
         write(iulog,*) '--------------------------------------------------------'
         write(iulog,*) 'euvac_set_etf: f107,f107a = ',f107,f107a
#ifdef EUVAC_DIAGS
         write(iulog,*) 'euvac_set_etf: etf'
         do w = 1,nbins
            write(iulog,'(1p,2g15.7)') euvac_etf(w)
         end do
#endif
         write(iulog,*) '--------------------------------------------------------'
         write(iulog,*) ' '
      end if

      end subroutine euvac_set_etf

      end module euvac
