!-----------------------------------------------------------------------
! Solar irradiance / photon flux data
!-----------------------------------------------------------------------
module solar_irrad_data
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun
  use pio,              only: file_desc_t,pio_inq_dimid,pio_inq_varid,pio_inq_dimlen, pio_noerr,pio_internal_error,pio_bcast_error
  use pio,              only: pio_get_var, pio_seterrorhandling
  use cam_pio_utils,    only: cam_pio_openfile
  use cam_logfile,      only: iulog
  use infnan,           only: nan, assignment(=)
  use input_data_utils, only: time_coordinate

  implicit none

  save
  private
  public :: solar_irrad_init
  public :: solar_irrad_advance

  integer,  public, protected :: nbins  ! number of wavelength samples of spectrum, wavelength endpoints
  real(r8), public, protected, allocatable :: we(:)
  real(r8), public, protected, allocatable :: sol_etf(:)
  real(r8), public, protected, allocatable :: ssi_ref(:)  ! a reference spectrum constructed from 3 solar cycles of data
  real(r8), public, protected, allocatable :: sol_irrad(:)
  real(r8), public, protected              :: sol_tsi = -1.0_r8
  real(r8), public, protected              :: ref_tsi
  logical,  public, protected :: do_spctrl_scaling = .false.
  logical,  public, protected :: has_spectrum = .false.
  logical,  public, protected :: has_ref_spectrum = .false.

  type(file_desc_t) :: file_id
  integer :: ssi_vid
  integer :: tsi_vid
  integer :: ref_vid
  integer :: tsi_ref_vid

  logical  :: initialized = .false.
  logical  :: has_tsi = .false.
  real(r8) :: itsi(2)
  real(r8), allocatable :: irradi(:,:)
  real(r8), allocatable :: irrad_fac(:)
  real(r8), allocatable :: etf_fac(:)
  real(r8), allocatable :: dellam(:)

  logical  :: fixed_scon

  type(time_coordinate) :: time_coord

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine solar_irrad_init(filepath, fixed, fixed_ymd, fixed_tod, const_tsi, heatng_spctrl_scl )

    use ioFileMod, only : getfil
    use physconst, only : c0, planck

    !---------------------------------------------------------------
    ! arguments
    !---------------------------------------------------------------
    character(len=*), intent(in) :: filepath
    logical,  intent(in) :: fixed
    integer,  intent(in) :: fixed_ymd
    integer,  intent(in) :: fixed_tod
    real(r8), intent(in) :: const_tsi
    logical,  intent(in) :: heatng_spctrl_scl

    !---------------------------------------------------------------
    ! local vars
    !---------------------------------------------------------------
    integer :: astat, dimid, vid
    character(len=256) :: filen   
    real(r8), allocatable :: lambda(:)
    integer :: i, wvl_vid
    real(r8), parameter :: c = c0     ! speed of light (m/s)
    real(r8), parameter :: h = planck ! Planck's constant (Joule sec)
    real(r8), parameter :: fac = 1._r8/(h*c)
    integer :: ierr

    has_spectrum = .false.

    if ( filepath.ne.'NONE' ) then
       fixed_scon = .false.
    else
       fixed_scon = .true.
    endif

    if ( const_tsi>0._r8 ) then
       sol_tsi = const_tsi
    endif
    ref_tsi = nan

    if ( fixed_scon ) return

    call time_coord%initialize( filepath, fixed=fixed, fixed_ymd=fixed_ymd, fixed_tod=fixed_tod, &
                                force_time_interp=.true., try_dates=.true. )

    call getfil( filepath, filen, 0 )
    call cam_pio_openfile( file_id, filen, 0 )
    if(masterproc)   write(iulog,*)'solar_data_init: data file = ',trim(filen)
    call pio_seterrorhandling(file_id, pio_bcast_error)
    ierr = pio_inq_varid( file_id, 'ssi', ssi_vid )
    has_spectrum = ierr==PIO_NOERR

    ierr = pio_inq_varid( file_id, 'tsi', tsi_vid )
    has_tsi = ierr==PIO_NOERR .and. const_tsi<0._r8

    ierr = pio_inq_varid( file_id, 'ssi_ref', ref_vid )
    has_ref_spectrum = ierr==PIO_NOERR
    call pio_seterrorhandling(file_id, pio_internal_error)

    if ( has_spectrum ) then
       call pio_seterrorhandling(file_id, pio_bcast_error)
       ierr = pio_inq_varid( file_id, 'wavelength', wvl_vid )
       call pio_seterrorhandling(file_id, pio_internal_error)
       
       if ( ierr==PIO_NOERR ) then
          ierr = pio_inq_dimid( file_id, 'wavelength', dimid )
       else ! for backwards compatibility
          ierr = pio_inq_varid( file_id, 'wvl', wvl_vid  )
          ierr = pio_inq_dimid( file_id, 'wvl', dimid )
       endif
       ierr = pio_inq_dimlen( file_id, dimid, nbins )
       if ( has_ref_spectrum ) then
          ierr = pio_inq_varid( file_id, 'tsi_ref', tsi_ref_vid )
       endif
    endif

    do_spctrl_scaling = has_spectrum .and. heatng_spctrl_scl

    allocate(lambda(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate lambda; error = ',astat
       call endrun('solar_data_init')
    end if
    allocate(dellam(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate dellam; error = ',astat
       call endrun('solar_data_init')
    end if
    allocate(irrad_fac(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate irrad_fac; error = ',astat
       call endrun('solar_data_init')
    end if
    allocate(etf_fac(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate etf_fac; error = ',astat
       call endrun('solar_data_init')
    end if
    allocate(sol_irrad(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate sol_irrad; error = ',astat
       call endrun('solar_data_init')
    end if
    allocate(ssi_ref(nbins), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'solar_data_init: failed to allocate ssi_ref; error = ',astat
       call endrun('solar_data_init')
    end if
    ssi_ref(:) = nan

    if (has_spectrum) then
       ierr = pio_get_var( file_id, wvl_vid, lambda )
       ierr = pio_inq_varid( file_id, 'band_width', vid  )
       ierr = pio_get_var( file_id, vid, dellam )
    endif

    if(masterproc)   write(iulog,*)'solar_data_init: has_ref_spectrum',has_ref_spectrum
    if ( has_ref_spectrum ) then
       ierr = pio_inq_varid( file_id, 'ssi_ref', vid  )
       ierr = pio_get_var( file_id, vid, ssi_ref )
       ierr = pio_get_var( file_id, tsi_ref_vid, ref_tsi )
    endif

    if ( has_spectrum ) then
       allocate(sol_etf(nbins), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'solar_data_init: failed to allocate sol_etf; error = ',astat
          call endrun('solar_data_init')
       end if
       allocate(irradi(nbins,2), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'solar_data_init: failed to allocate irradi; error = ',astat
          call endrun('solar_data_init')
       end if

       allocate(we(nbins+1), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'solar_data_init: failed to allocate we; error = ',astat
          call endrun('solar_data_init')
       end if

       we(:nbins)  = lambda(:nbins) - 0.5_r8*dellam(:nbins)
       we(nbins+1) = lambda(nbins)  + 0.5_r8*dellam(nbins)
       do i = 1,nbins
          irrad_fac(i) = 1.e-3_r8                ! mW/m2/nm --> W/m2/nm
          etf_fac(i)   = 1.e-16_r8*lambda(i)*fac ! mW/m2/nm --> photons/cm2/sec/nm
       enddo
       if(has_ref_spectrum) then
          ssi_ref = ssi_ref * 1.e-3_r8        ! mW/m2/nm --> W/m2/nm
       endif
    endif

    deallocate(lambda)
    deallocate(dellam)

    ! need to force data loading when the model starts at a time =/ 00:00:00.000
    ! -- may occur in restarts also
    call solar_irrad_advance()
    initialized = .true.

  end subroutine solar_irrad_init

!-----------------------------------------------------------------------
! Reads in the ETF data for the current date.  
!-----------------------------------------------------------------------
  subroutine solar_irrad_advance( )

    integer  :: i, index, nt
    integer  :: offset(2), count(2)
    logical  :: read_data
    real(r8) :: data(nbins)
    integer  :: ierr
    real(r8) :: delt

    if ( fixed_scon ) return
    if ( time_coord%fixed .and. initialized ) return

    index = -1

    read_data = time_coord%read_more() .or. .not.initialized
    call time_coord%advance()

    if ( read_data ) then
       nt = 2
       index = time_coord%indxs(1)

       ! get the surrounding time slices
       offset = (/ 1, index /)
       count =  (/ nbins, nt /)

       if (has_spectrum) then
          ierr = pio_get_var( file_id, ssi_vid, offset, count, irradi )
       endif
       if (has_tsi .and. (.not.do_spctrl_scaling)) then
          ierr = pio_get_var( file_id, tsi_vid, (/index/), (/nt/), itsi )
          if ( any(itsi(:nt) < 0._r8) ) then
             call endrun( 'solar_data_advance: invalid or missing tsi data  ' )
          endif
       endif
    endif

    delt = time_coord%wghts(2)

    if (has_spectrum) then
       data(:) = irradi(:,1) + delt*( irradi(:,2) - irradi(:,1) )
       
       do i = 1,nbins
          sol_irrad(i) = data(i)*irrad_fac(i) ! W/m2/nm
          sol_etf(i)   = data(i)*etf_fac(i)   ! photons/cm2/sec/nm 
       enddo
    endif
    if (has_tsi .and. (.not.do_spctrl_scaling)) then
       sol_tsi = itsi(1) + delt*( itsi(2) - itsi(1) )
    endif

  end subroutine solar_irrad_advance

end module solar_irrad_data
