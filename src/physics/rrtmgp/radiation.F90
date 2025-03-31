module radiation

!---------------------------------------------------------------------------------
!
! CAM interface to RRTMGP radiation parameterization.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,        only: r8=>shr_kind_r8, cl=>shr_kind_cl
use spmd_utils,          only: masterproc
use ppgrid,              only: pcols, pver, pverp, begchunk, endchunk
use ref_pres,            only: pref_edge
use physics_types,       only: physics_state, physics_ptend
use phys_control,        only: phys_getopts
use physics_buffer,      only: physics_buffer_desc, pbuf_add_field, dtype_r8, pbuf_get_index, &
                               pbuf_set_field, pbuf_get_field, pbuf_old_tim_idx
use camsrfexch,          only: cam_out_t, cam_in_t
use physconst,           only: cappa, cpair, gravit, stebol

use time_manager,        only: get_nstep, is_first_step, is_first_restart_step, &
                               get_curr_calday, get_step_size

use rad_constituents,    only: N_DIAG, rad_cnst_get_call_list, rad_cnst_get_gas, rad_cnst_out

use radconstants,        only: nradgas, gasnamelength, nswbands, nlwbands, &
                               gaslist, radconstants_init
use rad_solar_var,       only: rad_solar_var_init, get_variability

use cloud_rad_props,     only: cloud_rad_props_init

use cospsimulator_intr,  only: docosp, cospsimulator_intr_init, &
                               cospsimulator_intr_run, cosp_nradsteps

use scamMod,             only: scm_crm_mode, single_column, have_cld, cldobs

use cam_history,         only: addfld, add_default, horiz_only, outfld, hist_fld_active

use radiation_data,      only: rad_data_register, rad_data_init

use ioFileMod,           only: getfil
use cam_pio_utils,       only: cam_pio_openfile
use pio,                 only: file_desc_t, var_desc_t,                       &
                               pio_int, pio_double, PIO_NOERR,                &
                               pio_seterrorhandling, PIO_BCAST_ERROR,         &
                               pio_inq_dimlen, pio_inq_dimid, pio_inq_varid,  &
                               pio_def_var, pio_put_var, pio_get_var,         &
                               pio_put_att, PIO_NOWRITE, pio_closefile

use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
use ccpp_optical_props,      only: ty_optical_props_1scl_ccpp, ty_optical_props_2str_ccpp
use ccpp_source_functions,   only: ty_source_func_lw_ccpp
use ccpp_fluxes,             only: ty_fluxes_broadband_ccpp
use ccpp_fluxes_byband,    only: ty_fluxes_byband_ccpp

use string_utils,        only: to_lower
use cam_abortutils,      only: endrun, handle_allocate_error
use cam_logfile,         only: iulog
use rrtmgp_pre,          only: radiation_do_ccpp


implicit none
private
save

public :: &
   radiation_readnl,         &! read namelist variables
   radiation_do,             &! query which radiation calcs are done this timestep
   radiation_register,       &! registers radiation physics buffer fields
   radiation_init,           &! initialization
   radiation_define_restart, &! define variables for restart
   radiation_write_restart,  &! write variables to restart
   radiation_read_restart,   &! read variables from restart
   radiation_tend,           &! compute heating rates and fluxes
   rad_out_t                  ! type for diagnostic outputs

integer,public, allocatable :: cosp_cnt(:)       ! counter for cosp
integer,public              :: cosp_cnt_init = 0 !initial value for cosp counter

real(r8), public, protected :: nextsw_cday       ! future radiation calday for surface models

type rad_out_t
   real(r8) :: solin(pcols)         ! Solar incident flux

   real(r8) :: qrsc(pcols,pver)

   real(r8) :: fsnsc(pcols)         ! Clear sky surface abs solar flux
   real(r8) :: fsntc(pcols)         ! Clear sky total column abs solar flux
   real(r8) :: fsdsc(pcols)         ! Clear sky surface downwelling solar flux
   
   real(r8) :: fsntoa(pcols)        ! Net solar flux at TOA
   real(r8) :: fsntoac(pcols)       ! Clear sky net solar flux at TOA
   real(r8) :: fsutoa(pcols)        ! upwelling solar flux at TOA

   real(r8) :: fsnirt(pcols)        ! Near-IR flux absorbed at toa
   real(r8) :: fsnrtc(pcols)        ! Clear sky near-IR flux absorbed at toa
   real(r8) :: fsnirtsq(pcols)      ! Near-IR flux absorbed at toa >= 0.7 microns

   real(r8) :: fsn200(pcols)        ! Net SW flux interpolated to 200 mb
   real(r8) :: fsn200c(pcols)       ! Net clear-sky SW flux interpolated to 200 mb
   real(r8) :: fsnr(pcols)          ! Net SW flux interpolated to tropopause
   
   real(r8) :: flux_sw_up(pcols,pverp)     ! upward shortwave flux on interfaces
   real(r8) :: flux_sw_clr_up(pcols,pverp) ! upward shortwave clearsky flux
   real(r8) :: flux_sw_dn(pcols,pverp)     ! downward flux
   real(r8) :: flux_sw_clr_dn(pcols,pverp) ! downward clearsky flux

   real(r8) :: flux_lw_up(pcols,pverp)     ! upward longwave flux on interfaces
   real(r8) :: flux_lw_clr_up(pcols,pverp) ! upward longwave clearsky flux
   real(r8) :: flux_lw_dn(pcols,pverp)     ! downward flux
   real(r8) :: flux_lw_clr_dn(pcols,pverp) ! downward clearsky flux

   real(r8) :: qrlc(pcols,pver)

   real(r8) :: flntc(pcols)         ! Clear sky lw flux at model top
   real(r8) :: flut(pcols)          ! Upward flux at top of model
   real(r8) :: flutc(pcols)         ! Upward Clear Sky flux at top of model
   real(r8) :: lwcf(pcols)          ! longwave cloud forcing

   real(r8) :: fln200(pcols)        ! net longwave flux interpolated to 200 mb
   real(r8) :: fln200c(pcols)       ! net clearsky longwave flux interpolated to 200 mb
   real(r8) :: flnr(pcols)          ! net longwave flux interpolated to tropopause

   real(r8) :: flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
   real(r8) :: fldsc(pcols)         ! Clear sky lw flux at srf (down)

   real(r8) :: tot_cld_vistau(pcols,pver)   ! gbx water+ice cloud optical depth (only during day, night = fillvalue)
   real(r8) :: tot_icld_vistau(pcols,pver)  ! in-cld water+ice cloud optical depth (only during day, night = fillvalue)
   real(r8) :: liq_icld_vistau(pcols,pver)  ! in-cld liq cloud optical depth (only during day, night = fillvalue)
   real(r8) :: ice_icld_vistau(pcols,pver)  ! in-cld ice cloud optical depth (only during day, night = fillvalue)
   real(r8) :: snow_icld_vistau(pcols,pver) ! snow in-cloud visible sw optical depth for output on history files
   real(r8) :: grau_icld_vistau(pcols,pver) ! Graupel in-cloud visible sw optical depth for output on history files
end type rad_out_t

! Control variables set via namelist
character(len=cl) :: coefs_lw_file ! filepath for lw coefficients
character(len=cl) :: coefs_sw_file ! filepath for sw coefficients

integer :: iradsw = -1     ! freq. of shortwave radiation calc in time steps (positive)
                           ! or hours (negative).
integer :: iradlw = -1     ! frequency of longwave rad. calc. in time steps (positive)
                           ! or hours (negative).

integer :: irad_always = 0 ! Specifies length of time in timesteps (positive)
                           ! or hours (negative) SW/LW radiation will be
                           ! run continuously from the start of an
                           ! initial or restart run
logical :: use_rad_dt_cosz  = .false. ! if true, use radiation dt for all cosz calculations
logical :: spectralflux     = .false. ! calculate fluxes (up and down) per band.
logical :: graupel_in_rad   = .false. ! graupel in radiation code
logical :: use_rad_uniform_angle = .false. ! if true, use the namelist rad_uniform_angle for the coszrs calculation

! Gathered indices of day and night columns 
!  chunk_column_index = IdxDay(daylight_column_index)
integer :: nday           ! Number of daylight columns
integer :: nnite          ! Number of night columns
integer :: idxday(pcols)   ! chunk indices of daylight columns
integer :: idxnite(pcols) ! chunk indices of night columns
real(r8) :: coszrs(pcols)   ! Cosine solar zenith angle
real(r8) :: eccf            ! Earth orbit eccentricity factor

integer :: band2gpt_sw(2,nswbands)

! active_calls is set by a rad_constituents method after parsing namelist input
! for the rad_climate and rad_diag_N entries.
logical :: active_calls(0:N_DIAG)

! Physics buffer indices
integer :: qrs_idx      = 0 
integer :: qrl_idx      = 0 
integer :: su_idx       = 0 
integer :: sd_idx       = 0 
integer :: lu_idx       = 0 
integer :: ld_idx       = 0 
integer :: fsds_idx     = 0
integer :: fsns_idx     = 0
integer :: fsnt_idx     = 0
integer :: flns_idx     = 0
integer :: flnt_idx     = 0
integer :: cld_idx      = 0 
integer :: cldfsnow_idx = 0 
integer :: cldfgrau_idx = 0    
integer :: dei_idx
integer :: mu_idx
integer :: lambda_idx
integer :: iciwp_idx
integer :: iclwp_idx
integer :: des_idx
integer :: icswp_idx
integer :: icgrauwp_idx
integer :: degrau_idx

character(len=4) :: diag(0:N_DIAG) =(/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ',&
                                      '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

! averaging time interval for zenith angle
real(r8) :: dt_avg = 0._r8
real(r8) :: rad_uniform_angle = -99._r8

! Number of layers in radiation calculations.
integer :: nlay
! Number of interfaces in radiation calculations.
integer :: nlayp

! Number of CAM layers in radiation calculations.  Is either equal to nlay, or is
! 1 less than nlay if "extra layer" is used in the radiation calculations.
integer :: nlaycam

! Indices for copying data between CAM/WACCM and RRTMGP arrays.  Since RRTMGP is
! vertical order agnostic we can send data using the top to bottom order used
! in CAM/WACCM.  But the number of layers that RRTMGP does computations for
! may not match the number of layers in CAM/WACCM for two reasons:
! 1. If the CAM model top is below 1 Pa, then RRTMGP does calculations for an
!    extra layer that is added between 1 Pa and the model top.
! 2. If the WACCM model top is above 1 Pa, then RRMTGP only does calculations
!    for those model layers that are below 1 Pa.
integer :: ktopcam ! Index in CAM arrays of top level (layer or interface) at which
                   ! RRTMGP is active.
integer :: ktoprad ! Index in RRTMGP arrays of the layer or interface corresponding
                   ! to CAM's top layer or interface.
                   ! Note: for CAM's top to bottom indexing, the index of a given layer
                   ! (midpoint) and the upper interface of that layer, are the same.

integer :: nlwgpts
integer :: nswgpts

! Band indices for bands containing specific wavelengths
integer :: idx_sw_diag
integer :: idx_nir_diag
integer :: idx_uv_diag
integer :: idx_sw_cloudsim
integer :: idx_lw_diag
integer :: idx_lw_cloudsim

real(r8) :: sw_low_bounds(nswbands)
real(r8) :: sw_high_bounds(nswbands)

! Flag to perform shortwave or longwave on current timestep
logical :: dosw
logical :: dolw

! Gas optics objects contain the data read from the coefficients files.
type(ty_gas_optics_rrtmgp_ccpp) :: kdist_sw
type(ty_gas_optics_rrtmgp_ccpp) :: kdist_lw

! lower case version of gaslist for RRTMGP
character(len=gasnamelength) :: gaslist_lc(nradgas)

type(var_desc_t) :: cospcnt_desc  ! cosp
type(var_desc_t) :: nextsw_cday_desc

!=========================================================================================
contains
!=========================================================================================

subroutine radiation_readnl(nlfile)

   ! Read radiation_nl namelist group.

   use namelist_utils,  only: find_group_name
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_logical, &
                              mpi_character, mpi_real8

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   integer :: dtime      ! timestep size
   character(len=32) :: errmsg
   character(len=*), parameter :: sub = 'radiation_readnl'

   character(len=cl) :: rrtmgp_coefs_lw_file, rrtmgp_coefs_sw_file

   namelist /radiation_nl/ &
      rrtmgp_coefs_lw_file, rrtmgp_coefs_sw_file, iradsw, iradlw,        &
      irad_always, use_rad_dt_cosz, spectralflux, use_rad_uniform_angle, &
      rad_uniform_angle, graupel_in_rad
   !-----------------------------------------------------------------------------

   if (masterproc) then
      open( newunit=unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'radiation_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, radiation_nl, iostat=ierr)
         if (ierr /= 0) then
            write(errmsg,'(a,i5)') 'iostat =', ierr
            call endrun(sub//': ERROR reading namelist: '//trim(errmsg))
         end if
      end if
      close(unitn)
   end if

   ! Broadcast namelist variables
   call mpi_bcast(rrtmgp_coefs_lw_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rrtmgp_coefs_lw_file")
   call mpi_bcast(rrtmgp_coefs_sw_file, cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rrtmgp_coefs_sw_file")
   call mpi_bcast(iradsw, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: iradsw")
   call mpi_bcast(iradlw, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: iradlw")
   call mpi_bcast(irad_always, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: irad_always")
   call mpi_bcast(use_rad_dt_cosz, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: use_rad_dt_cosz")
   call mpi_bcast(spectralflux, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: spectralflux")
   call mpi_bcast(use_rad_uniform_angle, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: use_rad_uniform_angle")
   call mpi_bcast(rad_uniform_angle, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rad_uniform_angle")   
   call mpi_bcast(graupel_in_rad, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: graupel_in_rad")

   if (use_rad_uniform_angle .and. rad_uniform_angle == -99._r8) then
      call endrun(sub//': ERROR - use_rad_uniform_angle is set to .true,' &
                     //' but rad_uniform_angle is not set ')
   end if

   ! Set module data
   coefs_lw_file   = rrtmgp_coefs_lw_file
   coefs_sw_file   = rrtmgp_coefs_sw_file

   ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
   dtime  = get_step_size()
   if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
   if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
   if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

   !----------------------------------------------------------------------- 
   ! Print runtime options to log.
   !-----------------------------------------------------------------------

   if (masterproc) then
      write(iulog,*) 'RRTMGP radiation scheme parameters:'
      write(iulog,10) trim(coefs_lw_file), trim(coefs_sw_file), nlwbands, nswbands, &
         iradsw, iradlw, irad_always, use_rad_dt_cosz, spectralflux, graupel_in_rad
   end if

10 format('  LW coefficents file: ',                                a/, &
          '  SW coefficents file: ',                                a/, &
          '  Number of LW bands:                                 ',i5/, &
          '  Number of SW bands:                                 ',i5/, &
          '  Frequency (timesteps) of Shortwave Radiation calc:  ',i5/, &
          '  Frequency (timesteps) of Longwave Radiation calc:   ',i5/, &
          '  SW/LW calc done every timestep for first N steps. N=',i5/, &
          '  Use average zenith angle:                           ',l5/, &
          '  Output spectrally resolved fluxes:                  ',l5/, &
          '  Graupel in Radiation Code:                          ',l5/)

end subroutine radiation_readnl

!================================================================================================

subroutine radiation_register

   ! Register radiation fields in the physics buffer

   call pbuf_add_field('QRS' , 'global',dtype_r8,(/pcols,pver/), qrs_idx) ! shortwave radiative heating rate 
   call pbuf_add_field('QRL' , 'global',dtype_r8,(/pcols,pver/), qrl_idx) ! longwave  radiative heating rate 

   call pbuf_add_field('FSDS' , 'global',dtype_r8,(/pcols/), fsds_idx) ! Surface solar downward flux
   call pbuf_add_field('FSNS' , 'global',dtype_r8,(/pcols/), fsns_idx) ! Surface net shortwave flux
   call pbuf_add_field('FSNT' , 'global',dtype_r8,(/pcols/), fsnt_idx) ! Top-of-model net shortwave flux

   call pbuf_add_field('FLNS' , 'global',dtype_r8,(/pcols/), flns_idx) ! Surface net longwave flux
   call pbuf_add_field('FLNT' , 'global',dtype_r8,(/pcols/), flnt_idx) ! Top-of-model net longwave flux

   ! If the namelist has been configured for preserving the spectral fluxes, then create
   ! physics buffer variables to store the results.  This data is accessed by CARMA.
   if (spectralflux) then
      call pbuf_add_field('SU'  , 'global',dtype_r8,(/pcols,pverp,nswbands/), su_idx) ! shortwave upward flux (per band)
      call pbuf_add_field('SD'  , 'global',dtype_r8,(/pcols,pverp,nswbands/), sd_idx) ! shortwave downward flux (per band)
      call pbuf_add_field('LU'  , 'global',dtype_r8,(/pcols,pverp,nlwbands/), lu_idx) ! longwave upward flux (per band)
      call pbuf_add_field('LD'  , 'global',dtype_r8,(/pcols,pverp,nlwbands/), ld_idx) ! longwave downward flux (per band)
   end if

   ! Register fields for offline radiation driver.
   call rad_data_register()

end subroutine radiation_register

!================================================================================================

function radiation_do(op)

   ! Return true if the specified operation is done this timestep.

   character(len=*), intent(in) :: op             ! name of operation
   logical                      :: radiation_do   ! return value

   ! Local variables
   integer :: nstep             ! current timestep number
   integer :: errcode
   character(len=512) :: errmsg
   !-----------------------------------------------------------------------

   nstep = get_nstep()

   select case (op)
      case ('sw') ! do a shortwave heating calc this timestep?
         call radiation_do_ccpp(op, nstep, iradsw, irad_always, radiation_do, errmsg, errcode)
      case ('lw') ! do a longwave heating calc this timestep?
         call radiation_do_ccpp(op, nstep, iradlw, irad_always, radiation_do, errmsg, errcode)
      case default
         call endrun('radiation_do: unknown operation:'//op)
   end select

end function radiation_do

!================================================================================================

subroutine radiation_init(pbuf2d)
   use rrtmgp_inputs, only: rrtmgp_inputs_init
   use rrtmgp_inputs_cam, only: rrtmgp_inputs_cam_init

   ! Initialize the radiation and cloud optics.
   ! Add fields to the history buffer.

   ! arguments
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local variables
   character(len=512) :: errmsg

   ! names of gases that are available in the model 
   ! -- needed for the kdist initialization routines 
   type(ty_gas_concs_ccpp) :: available_gases

   real(r8) :: qrl_unused(1,1)

   integer :: i, icall
   integer :: nstep                       ! current timestep number
   logical :: history_amwg                ! output the variables used by the AMWG diag package
   logical :: history_vdiag               ! output the variables used by the AMWG variability diag package
   logical :: history_budget              ! output tendencies and state variables for CAM4
                                          ! temperature, water vapor, cloud ice and cloud
                                          ! liquid budgets.
   integer :: history_budget_histfile_num ! history file number for budget fields
   integer :: ierr, istat, errflg

   integer :: dtime

   character(len=*), parameter :: sub = 'radiation_init'
   !-----------------------------------------------------------------------
   
   ! Create lowercase version of the gaslist for RRTMGP.  The ty_gas_concs_ccpp objects
   ! work with CAM's uppercase names, but other objects that get input from the gas
   ! concs objects don't work.
   do i = 1, nradgas
      gaslist_lc(i) = to_lower(gaslist(i))
   end do

   ! PEVERWHEE - add this to new rrtmgp_pre_iinit routine (possible also above code?)
   errmsg = available_gases%gas_concs%init(gaslist_lc)
   call stop_on_err(errmsg, sub, 'available_gases%init')

   ! Read RRTMGP coefficients files and initialize kdist objects.
   call coefs_init(coefs_sw_file, available_gases, kdist_sw)
   call coefs_init(coefs_lw_file, available_gases, kdist_lw)

   ! Set up inputs to RRTMGP
   call rrtmgp_inputs_init(ktopcam, ktoprad, nlaycam, sw_low_bounds, sw_high_bounds, nswbands,               &
                   pref_edge, nlay, pver, pverp, kdist_sw, kdist_lw, qrl_unused, is_first_step(), use_rad_dt_cosz, &
                   get_step_size(), get_nstep(), iradsw, dt_avg, irad_always, is_first_restart_step(),       &
                   nlwbands, nradgas, gasnamelength, iulog, idx_sw_diag, idx_nir_diag, idx_uv_diag,          &
                   idx_sw_cloudsim, idx_lw_diag, idx_lw_cloudsim, gaslist, nswgpts, nlwgpts, nlayp,          &
                   nextsw_cday, get_curr_calday(), band2gpt_sw, errmsg, errflg)


   call rrtmgp_inputs_cam_init(ktopcam, ktoprad, idx_sw_diag, idx_nir_diag, idx_uv_diag, idx_sw_cloudsim, idx_lw_diag, &
           idx_lw_cloudsim)

   call radconstants_init(idx_sw_diag, idx_nir_diag, idx_uv_diag, idx_lw_diag)

   call rad_solar_var_init(nswbands)

   ! initialize output fields for offline driver
   call rad_data_init(pbuf2d)

   call cloud_rad_props_init()
  
   cld_idx      = pbuf_get_index('CLD')
   cldfsnow_idx = pbuf_get_index('CLDFSNOW', errcode=ierr)
   cldfgrau_idx = pbuf_get_index('CLDFGRAU', errcode=ierr)

   if (is_first_step()) then
      call pbuf_set_field(pbuf2d, qrl_idx, 0._r8)
   end if

   ! Set the radiation timestep for cosz calculations if requested using
   ! the adjusted iradsw value from radiation
   !if (use_rad_dt_cosz)  then
   !   dtime  = get_step_size()
   !   dt_avg = iradsw*dtime
   !end if

   ! Surface components to get radiation computed today
   !if (.not. is_first_restart_step()) then
   !   nextsw_cday = get_curr_calday()
   !end if

   call phys_getopts(history_amwg_out   = history_amwg,    &
                     history_vdiag_out  = history_vdiag,   &
                     history_budget_out = history_budget,  &
                     history_budget_histfile_num_out = history_budget_histfile_num)

   ! "irad_always" is number of time steps to execute radiation continuously from
   ! start of initial OR restart run
   !nstep = get_nstep()
   !if (irad_always > 0) then
   !   irad_always = irad_always + nstep
   !end if

   if (docosp) call cospsimulator_intr_init()

   allocate(cosp_cnt(begchunk:endchunk), stat=istat)
   call handle_allocate_error(istat, sub, 'cosp_cnt')
   if (is_first_restart_step()) then
      cosp_cnt(begchunk:endchunk) = cosp_cnt_init
   else
      cosp_cnt(begchunk:endchunk) = 0     
   end if

   ! Add fields to history buffer

   call addfld('TOT_CLD_VISTAU',  (/ 'lev' /), 'A',   '1',             &
               'Total gbx cloud extinction visible sw optical depth',  &
               sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('TOT_ICLD_VISTAU', (/ 'lev' /), 'A',  '1',              &
               'Total in-cloud extinction visible sw optical depth',   &
               sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('LIQ_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', &
               'Liquid in-cloud extinction visible sw optical depth',  &
               sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('ICE_ICLD_VISTAU', (/ 'lev' /), 'A',  '1',              &
               'Ice in-cloud extinction visible sw optical depth',     &
               sampling_seq='rad_lwsw', flag_xyfill=.true.)

   if (cldfsnow_idx > 0) then
      call addfld('SNOW_ICLD_VISTAU', (/ 'lev' /), 'A', '1',           &
                  'Snow in-cloud extinction visible sw optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
   end if
   if (cldfgrau_idx > 0 .and. graupel_in_rad) then
      call addfld('GRAU_ICLD_VISTAU', (/ 'lev' /), 'A', '1',           &
                  'Graupel in-cloud extinction visible sw optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
   endif

   ! get list of active radiation calls
   call rad_cnst_get_call_list(active_calls)

   ! Add shortwave radiation fields to history master field list.

   do icall = 0, N_DIAG

      if (active_calls(icall)) then

         call addfld('SOLIN'//diag(icall),    horiz_only,   'A', 'W/m2', &
                     'Solar insolation', sampling_seq='rad_lwsw')
         call addfld('QRS'//diag(icall),      (/ 'lev' /),  'A', 'K/s',  &
                     'Solar heating rate', sampling_seq='rad_lwsw')
         call addfld('QRSC'//diag(icall),     (/ 'lev' /),  'A', 'K/s',  &
                     'Clearsky solar heating rate', sampling_seq='rad_lwsw')
         call addfld('FSNT'//diag(icall),     horiz_only,   'A', 'W/m2', &
                     'Net solar flux at top of model', sampling_seq='rad_lwsw')
         call addfld('FSNTC'//diag(icall),    horiz_only,   'A', 'W/m2', &
                     'Clearsky net solar flux at top of model', sampling_seq='rad_lwsw')
         call addfld('FSNTOA'//diag(icall),   horiz_only,   'A', 'W/m2', &
                     'Net solar flux at top of atmosphere', sampling_seq='rad_lwsw')
         call addfld('FSNTOAC'//diag(icall),  horiz_only,   'A', 'W/m2', &
                     'Clearsky net solar flux at top of atmosphere', sampling_seq='rad_lwsw')
         call addfld('SWCF'//diag(icall),     horiz_only,   'A', 'W/m2', &
                     'Shortwave cloud forcing', sampling_seq='rad_lwsw')
         call addfld('FSUTOA'//diag(icall),   horiz_only,   'A', 'W/m2', &
                     'Upwelling solar flux at top of atmosphere', sampling_seq='rad_lwsw')
         call addfld('FSNIRTOA'//diag(icall), horiz_only,   'A', 'W/m2', &
                     'Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', sampling_seq='rad_lwsw')
         call addfld('FSNRTOAC'//diag(icall), horiz_only,   'A', 'W/m2', &
                      'Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', sampling_seq='rad_lwsw')
         call addfld('FSNRTOAS'//diag(icall), horiz_only,   'A', 'W/m2', &
                     'Net near-infrared flux (>= 0.7 microns) at top of atmosphere', sampling_seq='rad_lwsw')
         call addfld('FSN200'//diag(icall),   horiz_only,   'A', 'W/m2', &
                     'Net shortwave flux at 200 mb', sampling_seq='rad_lwsw')
         call addfld('FSN200C'//diag(icall),  horiz_only,   'A', 'W/m2', &
                     'Clearsky net shortwave flux at 200 mb', sampling_seq='rad_lwsw')
         call addfld('FSNR'//diag(icall),     horiz_only,   'A', 'W/m2', &
                     'Net solar flux at tropopause', sampling_seq='rad_lwsw')
         call addfld('SOLL'//diag(icall),     horiz_only,   'A', 'W/m2', &
                     'Solar downward near infrared direct  to surface', sampling_seq='rad_lwsw')
         call addfld('SOLS'//diag(icall),     horiz_only,   'A', 'W/m2', &
                     'Solar downward visible direct  to surface', sampling_seq='rad_lwsw')
         call addfld('SOLLD'//diag(icall),    horiz_only,   'A', 'W/m2', &
                     'Solar downward near infrared diffuse to surface', sampling_seq='rad_lwsw')
         call addfld('SOLSD'//diag(icall),    horiz_only,   'A', 'W/m2', &
                     'Solar downward visible diffuse to surface', sampling_seq='rad_lwsw')
         call addfld('FSNS'//diag(icall),     horiz_only,   'A', 'W/m2', &
                     'Net solar flux at surface', sampling_seq='rad_lwsw')
         call addfld('FSNSC'//diag(icall),    horiz_only,   'A', 'W/m2', &
                     'Clearsky net solar flux at surface', sampling_seq='rad_lwsw')
         call addfld('FSDS'//diag(icall),     horiz_only,   'A', 'W/m2', &
                     'Downwelling solar flux at surface', sampling_seq='rad_lwsw')
         call addfld('FSDSC'//diag(icall),    horiz_only,   'A', 'W/m2', &
                     'Clearsky downwelling solar flux at surface', sampling_seq='rad_lwsw')

         ! Fluxes on CAM grid
         call addfld('FUS'//diag(icall),      (/ 'ilev' /), 'I', 'W/m2', &
                     'Shortwave upward flux', sampling_seq='rad_lwsw')
         call addfld('FDS'//diag(icall),      (/ 'ilev' /), 'I', 'W/m2', &
                     'Shortwave downward flux', sampling_seq='rad_lwsw')
         call addfld('FUSC'//diag(icall),     (/ 'ilev' /), 'I', 'W/m2', &
                     'Shortwave clear-sky upward flux', sampling_seq='rad_lwsw')
         call addfld('FDSC'//diag(icall),     (/ 'ilev' /), 'I', 'W/m2', &
                     'Shortwave clear-sky downward flux', sampling_seq='rad_lwsw')

         if (history_amwg) then
            call add_default('SOLIN'//diag(icall),   1, ' ')
            call add_default('QRS'//diag(icall),     1, ' ')
            call add_default('FSNT'//diag(icall),    1, ' ')
            call add_default('FSNTC'//diag(icall),   1, ' ')
            call add_default('FSNTOA'//diag(icall),  1, ' ')
            call add_default('FSNTOAC'//diag(icall), 1, ' ')
            call add_default('SWCF'//diag(icall),    1, ' ')
            call add_default('FSNS'//diag(icall),    1, ' ')
            call add_default('FSNSC'//diag(icall),   1, ' ')
            call add_default('FSUTOA'//diag(icall),  1, ' ')
            call add_default('FSDSC'//diag(icall),   1, ' ')
            call add_default('FSDS'//diag(icall),    1, ' ')
         endif

      end if
   end do

   if (scm_crm_mode) then
      call add_default('FUS     ', 1, ' ')
      call add_default('FUSC    ', 1, ' ')
      call add_default('FDS     ', 1, ' ')
      call add_default('FDSC    ', 1, ' ')
   endif

   ! Add longwave radiation fields to history master field list.

   do icall = 0, N_DIAG

      if (active_calls(icall)) then
         call addfld('QRL'//diag(icall),     (/ 'lev' /), 'A', 'K/s',  &
                     'Longwave heating rate', sampling_seq='rad_lwsw')
         call addfld('QRLC'//diag(icall),    (/ 'lev' /), 'A', 'K/s',  &
                     'Clearsky longwave heating rate', sampling_seq='rad_lwsw')
         call addfld('FLNT'//diag(icall),    horiz_only,  'A', 'W/m2', &
                     'Net longwave flux at top of model', sampling_seq='rad_lwsw')
         call addfld('FLNTC'//diag(icall),   horiz_only,  'A', 'W/m2', &
                     'Clearsky net longwave flux at top of model', sampling_seq='rad_lwsw')
         call addfld('FLUT'//diag(icall),    horiz_only,  'A', 'W/m2', &
                     'Upwelling longwave flux at top of model', sampling_seq='rad_lwsw')
         call addfld('FLUTC'//diag(icall),   horiz_only,  'A', 'W/m2', &
                     'Clearsky upwelling longwave flux at top of model', sampling_seq='rad_lwsw')
         call addfld('LWCF'//diag(icall),    horiz_only,  'A', 'W/m2', &
                     'Longwave cloud forcing', sampling_seq='rad_lwsw')
         call addfld('FLN200'//diag(icall),  horiz_only,  'A', 'W/m2', &
                     'Net longwave flux at 200 mb', sampling_seq='rad_lwsw')
         call addfld('FLN200C'//diag(icall), horiz_only,  'A', 'W/m2', &
                     'Clearsky net longwave flux at 200 mb', sampling_seq='rad_lwsw')
         call addfld('FLNR'//diag(icall),    horiz_only,  'A', 'W/m2', &
                     'Net longwave flux at tropopause', sampling_seq='rad_lwsw')
         call addfld('FLNS'//diag(icall),    horiz_only,  'A', 'W/m2', &
                     'Net longwave flux at surface', sampling_seq='rad_lwsw')
         call addfld('FLNSC'//diag(icall),   horiz_only,  'A', 'W/m2', &
                     'Clearsky net longwave flux at surface', sampling_seq='rad_lwsw')
         call addfld('FLDS'//diag(icall),    horiz_only,  'A', 'W/m2', &
                     'Downwelling longwave flux at surface', sampling_seq='rad_lwsw')
         call addfld('FLDSC'//diag(icall),   horiz_only,  'A', 'W/m2', &
                     'Clearsky Downwelling longwave flux at surface', sampling_seq='rad_lwsw')

         ! Fluxes on CAM grid
         call addfld('FUL'//diag(icall),    (/ 'ilev' /), 'I', 'W/m2', &
                     'Longwave upward flux', sampling_seq='rad_lwsw')
         call addfld('FDL'//diag(icall),    (/ 'ilev' /), 'I', 'W/m2', &
                     'Longwave downward flux', sampling_seq='rad_lwsw')
         call addfld('FULC'//diag(icall),   (/ 'ilev' /), 'I', 'W/m2', &
                     'Longwave clear-sky upward flux', sampling_seq='rad_lwsw')
         call addfld('FDLC'//diag(icall),   (/ 'ilev' /), 'I', 'W/m2', &
                     'Longwave clear-sky downward flux', sampling_seq='rad_lwsw')

         if (history_amwg) then
            call add_default('QRL'//diag(icall),   1, ' ')
            call add_default('FLNT'//diag(icall),  1, ' ')
            call add_default('FLNTC'//diag(icall), 1, ' ')
            call add_default('FLUT'//diag(icall),  1, ' ')
            call add_default('FLUTC'//diag(icall), 1, ' ')
            call add_default('LWCF'//diag(icall),  1, ' ')
            call add_default('FLNS'//diag(icall),  1, ' ')
            call add_default('FLNSC'//diag(icall), 1, ' ')
            call add_default('FLDS'//diag(icall),  1, ' ')
         end if

      end if
   end do

   call addfld('EMIS', (/ 'lev' /), 'A', '1', 'Cloud longwave emissivity')

   if (scm_crm_mode) then
      call add_default ('FUL     ', 1, ' ')
      call add_default ('FULC    ', 1, ' ')
      call add_default ('FDL     ', 1, ' ')
      call add_default ('FDLC    ', 1, ' ')
   endif

   ! Heating rate needed for d(theta)/dt computation
   call addfld ('HR',(/ 'lev' /), 'A','K/s','Heating rate needed for d(theta)/dt computation')

   if ( history_budget .and. history_budget_histfile_num > 1 ) then
      call add_default ('QRL     ', history_budget_histfile_num, ' ')
      call add_default ('QRS     ', history_budget_histfile_num, ' ')
   end if

   if (history_vdiag) then
      call add_default('FLUT', 2, ' ')
      call add_default('FLUT', 3, ' ')
   end if
   
end subroutine radiation_init

!===============================================================================

subroutine radiation_define_restart(file)

   ! define variables to be written to restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables
   integer :: ierr
   !----------------------------------------------------------------------------

   call pio_seterrorhandling(file, PIO_BCAST_ERROR)

   ierr = pio_def_var(file, 'nextsw_cday', pio_double, nextsw_cday_desc)
   ierr = pio_put_att(file, nextsw_cday_desc, 'long_name', 'future radiation calday for surface models')
   if (docosp) then
      ierr = pio_def_var(File, 'cosp_cnt_init', pio_int, cospcnt_desc)
   end if

end subroutine radiation_define_restart
  
!===============================================================================

subroutine radiation_write_restart(file)

   ! write variables to restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables
   integer :: ierr
   !----------------------------------------------------------------------------
   ierr = pio_put_var(File, nextsw_cday_desc, (/ nextsw_cday /))
   if (docosp) then
      ierr = pio_put_var(File, cospcnt_desc, (/cosp_cnt(begchunk)/))
   end if

end subroutine radiation_write_restart
  
!===============================================================================

subroutine radiation_read_restart(file)

   ! read variables from restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables
   integer :: ierr
   type(var_desc_t) :: vardesc
   integer :: err_handling

   !----------------------------------------------------------------------------
   if (docosp) then
      call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)
      ierr = pio_inq_varid(File, 'cosp_cnt_init', vardesc)
      call pio_seterrorhandling(File, err_handling)
      if (ierr /= PIO_NOERR) then
         cosp_cnt_init = 0
      else
         ierr = pio_get_var(File, vardesc, cosp_cnt_init)
      end if
   end if

   ierr = pio_inq_varid(file, 'nextsw_cday', vardesc)
   ierr = pio_get_var(file, vardesc, nextsw_cday)


end subroutine radiation_read_restart
  
!===============================================================================

subroutine radiation_tend( &
   state, ptend, pbuf, cam_out, cam_in, net_flx, rd_out)

   !----------------------------------------------------------------------- 
   ! 
   ! CAM driver for radiation computation.
   ! 
   !-----------------------------------------------------------------------

   ! Location/Orbital Parameters for cosine zenith angle
   use phys_grid,          only: get_rlat_all_p, get_rlon_all_p
   use cam_control_mod,    only: eccen, mvelpp, lambm0, obliqr
   use shr_orb_mod,        only: shr_orb_decl, shr_orb_cosz

   use rrtmgp_inputs,      only: rrtmgp_inputs_run
   use rrtmgp_pre,         only: rrtmgp_pre_run
   use rrtmgp_lw_cloud_optics,      only: rrtmgp_lw_cloud_optics_run
   use rrtmgp_lw_mcica_subcol_gen,  only: rrtmgp_lw_mcica_subcol_gen_run
   use rrtmgp_lw_gas_optics_pre,    only: rrtmgp_lw_gas_optics_pre_run
   use rrtmgp_lw_gas_optics,        only: rrtmgp_lw_gas_optics_run
   use rrtmgp_lw_main,              only: rrtmgp_lw_main_run
   use rrtmgp_dry_static_energy_tendency, only: rrtmgp_dry_static_energy_tendency_run
   use rrtmgp_post,                 only: rrtmgp_post_run

   use rrtmgp_inputs_cam,  only: rrtmgp_get_gas_mmrs, &
                                 rrtmgp_set_aer_lw, rrtmgp_set_gases_sw, rrtmgp_set_cloud_sw, &
                                 rrtmgp_set_aer_sw

   ! RRTMGP drivers for flux calculations.
   use mo_rte_lw,          only: rte_lw
   use mo_rte_sw,          only: rte_sw

   use radheat,            only: radheat_tend

   use radiation_data,     only: rad_data_write

   use interpolate_data,   only: vertinterp
   use tropopause,         only: tropopause_find_cam, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE
   use cospsimulator_intr, only: docosp, cospsimulator_intr_run, cosp_nradsteps


   ! Arguments
   type(physics_state), intent(in), target :: state
   type(physics_ptend), intent(out)        :: ptend
   type(physics_buffer_desc), pointer      :: pbuf(:)
   type(cam_out_t),     intent(inout)      :: cam_out
   type(cam_in_t),      intent(in)         :: cam_in
   real(r8),            intent(out)        :: net_flx(pcols)

   type(rad_out_t), target, optional, intent(out) :: rd_out


   ! Local variables
   type(rad_out_t), pointer :: rd  ! allow rd_out to be optional by allocating a local object
                                   ! if the argument is not present
   logical  :: write_output
  
   integer  :: i, k, gas_idx, istat
   integer  :: lchnk, ncol
   logical  :: dosw, dolw
   integer  :: icall           ! loop index for climate/diagnostic radiation calls

   real(r8) :: calday          ! current calendar day
   real(r8) :: delta           ! Solar declination angle  in radians
   real(r8) :: eccf            ! Earth orbit eccentricity factor
   real(r8) :: clat(pcols)     ! current latitudes(radians)
   real(r8) :: clon(pcols)     ! current longitudes(radians)
   real(r8) :: coszrs(pcols)   ! Cosine solar zenith angle

   integer :: itim_old

   real(r8), pointer :: cld(:,:)      ! cloud fraction
   real(r8), pointer :: cldfsnow(:,:) ! cloud fraction of just "snow clouds"
   real(r8), pointer :: cldfgrau(:,:) ! cloud fraction of just "graupel clouds"
   real(r8)          :: cldfprime(pcols,pver)   ! combined cloud fraction
   real(r8), pointer :: qrs(:,:) ! shortwave radiative heating rate 
   real(r8), pointer :: qrl(:,:) ! longwave  radiative heating rate 
   real(r8), pointer :: fsds(:)  ! Surface solar down flux
   real(r8), pointer :: fsns(:)  ! Surface solar absorbed flux
   real(r8), pointer :: fsnt(:)  ! Net column abs solar flux at model top
   real(r8), pointer :: flns(:)  ! Srf longwave cooling (up-down) flux
   real(r8), pointer :: flnt(:)  ! Net outgoing lw flux at model top

   real(r8), pointer :: dei(:,:)
   real(r8), pointer :: mu(:,:)
   real(r8), pointer :: lambda(:,:)
   real(r8), pointer :: iciwp(:,:)
   real(r8), pointer :: iclwp(:,:)
   real(r8), pointer :: des(:,:)
   real(r8), pointer :: icswp(:,:)
   real(r8), pointer :: icgrauwp(:,:)
   real(r8), pointer :: degrau(:,:)

   real(r8), pointer, dimension(:,:,:) :: su => NULL()  ! shortwave spectral flux up
   real(r8), pointer, dimension(:,:,:) :: sd => NULL()  ! shortwave spectral flux down
   real(r8), pointer, dimension(:,:,:) :: lu => NULL()  ! longwave  spectral flux up
   real(r8), pointer, dimension(:,:,:) :: ld => NULL()  ! longwave  spectral flux down

   ! tropopause diagnostic
   integer  :: troplev(pcols)
   real(r8) :: p_trop(pcols)

   ! state data passed to radiation calc
   real(r8), allocatable :: t_sfc(:)
   real(r8), allocatable :: emis_sfc(:,:)
   real(r8), allocatable :: t_rad(:,:)
   real(r8), allocatable :: pmid_rad(:,:)
   real(r8), allocatable :: pint_rad(:,:)
   real(r8), allocatable :: t_day(:,:)
   real(r8), allocatable :: pmid_day(:,:)
   real(r8), allocatable :: pint_day(:,:)
   real(r8), allocatable :: coszrs_day(:)
   real(r8), allocatable :: alb_dir(:,:)
   real(r8), allocatable :: alb_dif(:,:)
   real(r8), allocatable :: tauc(:,:,:)
   real(r8), allocatable :: cldf(:,:)

   real(r8), allocatable :: gas_mmrs(:,:,:)

   ! in-cloud optical depths for COSP
   real(r8) :: cld_tau_cloudsim(pcols,pver)    ! liq + ice
   real(r8) :: snow_tau_cloudsim(pcols,pver)   ! snow
   real(r8) :: grau_tau_cloudsim(pcols,pver)   ! graupel
   real(r8) :: cld_lw_abs_cloudsim(pcols,pver) ! liq + ice
   real(r8) :: snow_lw_abs_cloudsim(pcols,pver)! snow
   real(r8) :: grau_lw_abs_cloudsim(pcols,pver)! graupel

   ! Set vertical indexing in RRTMGP to be the same as CAM (top to bottom).
   logical, parameter :: top_at_1 = .true.

   logical :: do_graupel, do_snow

   ! TOA solar flux on RRTMGP g-points
   real(r8), allocatable :: toa_flux(:,:)
   ! Scale factors based on spectral distribution from input irradiance dataset
   real(r8), allocatable :: sfac(:,:)
   
   ! Planck sources for LW.
   type(ty_source_func_lw_ccpp) :: sources_lw

   ! Gas volume mixing ratios.  Use separate objects for LW and SW because SW only does
   ! calculations for daylight columns.
   ! These objects have a final method which deallocates the internal memory when they
   ! go out of scope (i.e., when radiation_tend returns), so no need for explicit deallocation.
   type(ty_gas_concs_ccpp) :: gas_concs_lw
   type(ty_gas_concs_ccpp) :: gas_concs_sw

   ! Atmosphere optics.  This object is initialized with gas optics, then is incremented
   ! by the aerosol optics for the clear-sky radiative flux calculations, and then
   ! incremented again by the cloud optics for the all-sky radiative flux calculations.
   type(ty_optical_props_1scl_ccpp) :: atm_optics_lw
   type(ty_optical_props_2str_ccpp) :: atm_optics_sw

   ! Cloud optical properties objects (McICA sampling of cloud optical properties).
   type(ty_optical_props_1scl_ccpp) :: cloud_lw
   type(ty_optical_props_2str_ccpp) :: cloud_sw

   ! Aerosol optical properties objects.
   type(ty_optical_props_1scl_ccpp) :: aer_lw
   type(ty_optical_props_2str_ccpp) :: aer_sw

   ! Flux objects contain all fluxes computed by RRTMGP.
   ! SW allsky fluxes always include spectrally resolved fluxes needed for surface models.
   type(ty_fluxes_byband_ccpp) :: fsw
   ! LW allsky fluxes only need spectrally resolved fluxes when spectralflux=.true.
   type(ty_fluxes_byband_ccpp) :: flw
   ! Only broadband fluxes needed for clear sky (diagnostics).
   type(ty_fluxes_broadband_ccpp) :: fswc, flwc

   ! Arrays for output diagnostics on CAM grid.
   real(r8) :: fns(pcols,pverp)     ! net shortwave flux
   real(r8) :: fcns(pcols,pverp)    ! net clear-sky shortwave flux
   real(r8) :: fnl(pcols,pverp)     ! net longwave flux
   real(r8) :: fcnl(pcols,pverp)    ! net clear-sky longwave flux

   ! Unused variables for rte_lw
   real(r8) :: fluxlwup_jac(1,1)
   real(r8) :: lw_ds(1,1)

   ! for COSP
   real(r8) :: emis(pcols,pver)        ! Cloud longwave emissivity
   real(r8) :: gb_snow_tau(pcols,pver) ! grid-box mean snow_tau
   real(r8) :: gb_snow_lw(pcols,pver)  ! grid-box mean LW snow optical depth

   real(r8) :: ftem(pcols,pver)        ! Temporary workspace for outfld variables

   character(len=128) :: errmsg
   integer            :: errflg, err
   character(len=*), parameter :: sub = 'radiation_tend'
   !--------------------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol

   if (present(rd_out)) then
      rd => rd_out
      write_output = .false.
   else
      allocate(rd, stat=istat)
      call handle_allocate_error(istat, sub, 'rd')
      write_output = .true.
   end if

   ! Cosine solar zenith angle for current time step
   calday = get_curr_calday()
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)

   call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, &
                     delta, eccf)

   if (use_rad_uniform_angle) then
      do i = 1, ncol
         coszrs(i) = shr_orb_cosz(calday, clat(i), clon(i), delta, dt_avg, &
                                  uniform_angle=rad_uniform_angle)
      end do
   else
      do i = 1, ncol
         ! if dt_avg /= 0, it triggers using avg coszrs
         coszrs(i) = shr_orb_cosz(calday, clat(i), clon(i), delta, dt_avg)
      end do
   end if

   call rrtmgp_pre_run(coszrs, get_nstep(), get_step_size(), iradsw, iradlw, irad_always, &
           ncol, nextsw_cday, idxday, nday, idxnite, nnite, dosw, dolw, nlay, nlwbands,   &
           nswbands, spectralflux, fsw, fswc, flw, flwc, errmsg, errflg)
   if (errflg /= 0) then
      call endrun(sub//': '//errmsg)
   end if

   ! Associate pointers to physics buffer fields
   itim_old = pbuf_old_tim_idx()
   nullify(cldfsnow)
   if (cldfsnow_idx > 0) then
      call pbuf_get_field(pbuf, cldfsnow_idx, cldfsnow, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   end if
   nullify(cldfgrau)
   if (cldfgrau_idx > 0 .and. graupel_in_rad) then
      call pbuf_get_field(pbuf, cldfgrau_idx, cldfgrau, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   endif

   call pbuf_get_field(pbuf, cld_idx, cld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, qrs_idx, qrs)
   call pbuf_get_field(pbuf, qrl_idx, qrl)

   call pbuf_get_field(pbuf, fsns_idx, fsns)
   call pbuf_get_field(pbuf, fsnt_idx, fsnt)
   call pbuf_get_field(pbuf, fsds_idx, fsds)
   call pbuf_get_field(pbuf, flns_idx, flns)
   call pbuf_get_field(pbuf, flnt_idx, flnt)

   if (spectralflux) then
      call pbuf_get_field(pbuf, su_idx, su)
      call pbuf_get_field(pbuf, sd_idx, sd)
      call pbuf_get_field(pbuf, lu_idx, lu)
      call pbuf_get_field(pbuf, ld_idx, ld)
   end if

   !  For CRM, make cloud equal to input observations:
   if (scm_crm_mode .and. have_cld) then
      do k = 1, pver
         cld(:ncol,k)= cldobs(k)
      end do
   end if

   ! Find tropopause height if needed for diagnostic output
   if (hist_fld_active('FSNR') .or. hist_fld_active('FLNR')) then
      !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
      troplev(:) = 0
      p_trop(:) = 0._r8
      !REMOVECAM_END
      call tropopause_find_cam(state, troplev, tropP=p_trop, primary=TROP_ALG_HYBSTOB, &
                           backup=TROP_ALG_CLIMATE)
   end if

   if (dosw .or. dolw) then

      allocate( &
         t_sfc(ncol), emis_sfc(nlwbands,ncol), toa_flux(nday,nswgpts),     &
         sfac(nday,nswgpts),                                               &
         t_rad(ncol,nlay), pmid_rad(ncol,nlay), pint_rad(ncol,nlay+1),     &
         t_day(nday,nlay), pmid_day(nday,nlay), pint_day(nday,nlay+1),     &
         coszrs_day(nday), alb_dir(nswbands,nday), alb_dif(nswbands,nday), &
         cldf(ncol,nlaycam), tauc(nlwbands,ncol,nlaycam), stat=istat)
      call handle_allocate_error(istat, sub, 't_sfc,..,alb_dif')
      allocate(gas_mmrs(ncol, pver, nradgas), stat=istat, errmsg=errmsg)
      if (errflg /= 0) then
         call handle_allocate_error(istat, sub, 'gas_mmrs, message: '//errmsg)
      end if

      ! Prepares state variables, daylit columns, albedos for RRTMGP
      ! Also calculates modified cloud fraction
      call rrtmgp_inputs_run(dosw, dolw, associated(cldfsnow), associated(cldfgrau), &
                  state%pmid, state%pint, state%t, &
                  nday, idxday, cldfprime, coszrs, kdist_sw, t_sfc,       &
                  emis_sfc, t_rad, pmid_rad, pint_rad, t_day, pmid_day,   &
                  pint_day, coszrs_day, alb_dir, alb_dif, cam_in%lwup, stebol,  &
                  ncol, ktopcam, ktoprad, nswbands, cam_in%asdir, cam_in%asdif, &
                  sw_low_bounds, sw_high_bounds, cam_in%aldir, cam_in%aldif, nlay,     &
                  pverp, pver, cld, cldfsnow, cldfgrau, graupel_in_rad,   &
                  gasnamelength, gaslist, gas_concs_lw, aer_lw, atm_optics_lw, &
                  kdist_lw, sources_lw, aer_sw, atm_optics_sw, gas_concs_sw,   &
                  errmsg, errflg)
      if (errflg /= 0) then
         call endrun(sub//': '//errmsg)
      end if

      ! Output the mass per layer, and total column burdens for gas and aerosol
      ! constituents in the climate list.
      call rad_cnst_out(0, state, pbuf)

      !========================!
      ! SHORTWAVE calculations !
      !========================!

      if (dosw) then

         ! Set cloud optical properties in cloud_sw object.
         call rrtmgp_set_cloud_sw( &
            state, pbuf, nlay, nday, idxday, nswgpts,            &
            nnite, idxnite, pmid_day, cld, cldfsnow,                      &
            cldfgrau, cldfprime, graupel_in_rad, kdist_sw, cloud_sw,      &
            rd%tot_cld_vistau, rd%tot_icld_vistau, rd%liq_icld_vistau,    &
            rd%ice_icld_vistau, rd%snow_icld_vistau, rd%grau_icld_vistau, &
            cld_tau_cloudsim, snow_tau_cloudsim, grau_tau_cloudsim )

         if (write_output) then
            call radiation_output_cld(lchnk, rd)
         end if

         ! The climate (icall==0) calculation must occur last.
         do icall = N_DIAG, 0, -1
            if (active_calls(icall)) then

               if (nday > 0) then

                  ! Set gas volume mixing ratios for this call in gas_concs_sw.
                  call rrtmgp_set_gases_sw( &
                     icall, state, pbuf, nlay, nday, &
                     idxday, gas_concs_sw)

                  ! Compute the gas optics (stored in atm_optics_sw).
                  ! toa_flux is the reference solar source from RRTMGP data.
                  errmsg = kdist_sw%gas_props%gas_optics( &
                     pmid_day, pint_day, t_day, gas_concs_sw%gas_concs, atm_optics_sw%optical_props, &
                     toa_flux)
                  call stop_on_err(errmsg, sub, 'kdist_sw%gas_props%gas_optics')

                  ! Scale the solar source
                  call get_variability(toa_flux, sfac, band2gpt_sw, nswbands)
                  toa_flux = toa_flux * sfac * eccf

               end if

               ! Set SW aerosol optical properties in the aer_sw object.
               ! This call made even when no daylight columns because it does some
               ! diagnostic aerosol output.
               call rrtmgp_set_aer_sw( &
                  icall, state, pbuf, nday, idxday, nnite, idxnite, aer_sw)
                  
               if (nday > 0) then

                  ! Increment the gas optics (in atm_optics_sw) by the aerosol optics in aer_sw.
                  errmsg = aer_sw%optical_props%increment(atm_optics_sw%optical_props)
                  call stop_on_err(errmsg, sub, 'aer_sw%optical_props%increment')

                  ! Compute clear-sky fluxes.
                  errmsg = rte_sw(&
                     atm_optics_sw%optical_props, top_at_1, coszrs_day, toa_flux, &
                     alb_dir, alb_dif, fswc%fluxes)
                  call stop_on_err(errmsg, sub, 'clear-sky rte_sw')

                  ! Increment the aerosol+gas optics (in atm_optics_sw) by the cloud optics in cloud_sw.
                  errmsg = cloud_sw%optical_props%increment(atm_optics_sw%optical_props)
                  call stop_on_err(errmsg, sub, 'cloud_sw%optical_props%increment')

                  ! Compute all-sky fluxes.
                  errmsg = rte_sw(&
                     atm_optics_sw%optical_props, top_at_1, coszrs_day, toa_flux, &
                     alb_dir, alb_dif, fsw%fluxes)
                  call stop_on_err(errmsg, sub, 'all-sky rte_sw')

               end if

               ! Transform RRTMGP outputs to CAM outputs and compute heating rates.
               call set_sw_diags()

               if (write_output) then
                  call radiation_output_sw(lchnk, ncol, icall, rd, pbuf, cam_out)
               end if

            end if ! (active_calls(icall))
         end do    ! loop over diagnostic calcs (icall)
      end if  ! if (dosw)

      !=======================!
      ! LONGWAVE calculations !
      !=======================!

      if (dolw) then

         ! Grab additional pbuf fields for LW cloud optics
         dei_idx = pbuf_get_index('DEI',errcode=err)
         mu_idx  = pbuf_get_index('MU',errcode=err)
         lambda_idx = pbuf_get_index('LAMBDAC',errcode=err)
         iciwp_idx  = pbuf_get_index('ICIWP',errcode=err)
         iclwp_idx  = pbuf_get_index('ICLWP',errcode=err)
         des_idx    = pbuf_get_index('DES',errcode=err)
         icswp_idx  = pbuf_get_index('ICSWP',errcode=err)
         icgrauwp_idx  = pbuf_get_index('ICGRAUWP',errcode=err) ! Available when using MG3
         degrau_idx    = pbuf_get_index('DEGRAU',errcode=err)   ! Available when using MG3
         call pbuf_get_field(pbuf, lambda_idx,  lambda)
         call pbuf_get_field(pbuf, mu_idx,      mu)
         call pbuf_get_field(pbuf, iclwp_idx,   iclwp)
         call pbuf_get_field(pbuf, iciwp_idx, iciwp)
         call pbuf_get_field(pbuf, dei_idx,   dei)
         call pbuf_get_field(pbuf, icswp_idx, icswp)
         call pbuf_get_field(pbuf, des_idx,   des)
         if (icgrauwp_idx > 0) then
            call pbuf_get_field(pbuf, icgrauwp_idx, icgrauwp)
         end if
         if (degrau_idx > 0) then
            call pbuf_get_field(pbuf, degrau_idx,   degrau)
         end if
         do_graupel = ((icgrauwp_idx > 0) .and. (degrau_idx > 0) .and. associated(cldfgrau))
         do_snow = associated(cldfsnow)
         ! Set cloud optical properties in cloud_lw object.
         call rrtmgp_lw_cloud_optics_run(ncol, nlay, nlaycam, cld, cldfsnow, cldfgrau,       &
             cldfprime, graupel_in_rad, kdist_lw, cloud_lw, lambda, mu, iclwp, iciwp,  &
             dei, icswp, des, icgrauwp, degrau, nlwbands, do_snow, &
             do_graupel, cld_lw_abs_cloudsim, snow_lw_abs_cloudsim, pver, ktopcam, &
             grau_lw_abs_cloudsim, idx_lw_cloudsim, tauc, cldf, errmsg, errflg)

         if (errflg /= 0) then
            call endrun(sub//': '//errmsg)
         end if
         call rrtmgp_lw_mcica_subcol_gen_run(ktoprad, &
                 kdist_lw, nlwbands, nlwgpts, ncol, pver, nlaycam, nlwgpts, &
                 state%pmid, cldf, tauc, cloud_lw, errmsg, errflg )
         if (errflg /= 0) then
            call endrun(sub//': '//errmsg)
         end if

         ! The climate (icall==0) calculation must occur last.
         do icall = N_DIAG, 0, -1

            if (active_calls(icall)) then

               ! Grab the gas mass mixing ratios from rad_constituents
               call rrtmgp_get_gas_mmrs(icall, state, pbuf, nlay, gas_mmrs)

               ! Set gas volume mixing ratios for this call in gas_concs_lw
               call rrtmgp_lw_gas_optics_pre_run(icall, gas_mmrs, state%pmid, state%pint, nlay, ncol, gaslist, idxday, &
                  pverp, ktoprad, ktopcam, dolw, nradgas, gas_concs_lw, errmsg, errflg)
               if (errflg /= 0) then
                  call endrun(sub//': '//errmsg)
               end if

               ! Compute the gas optics and Planck sources.
               call rrtmgp_lw_gas_optics_run(dolw, 1, ncol, ncol, pmid_rad, pint_rad, t_rad,  &
                  t_sfc, gas_concs_lw, atm_optics_lw, sources_lw, t_rad, .false., kdist_lw, errmsg, &
                  errflg)
               if (errflg /= 0) then
                  call endrun(sub//': '//errmsg)
               end if

               ! Set LW aerosol optical properties in the aer_lw object.
               call rrtmgp_set_aer_lw(icall, state, pbuf, aer_lw)

               call rrtmgp_lw_main_run(dolw, dolw, .true., .false., .false., &
                                 0, ncol, 1, ncol, atm_optics_lw, &
                                 cloud_lw, top_at_1, sources_lw, emis_sfc, kdist_lw, &
                                 aer_lw, fluxlwup_jac, lw_ds, flwc, flw, errmsg, errflg)
               if (errflg /= 0) then
                  call endrun(sub//': '//errmsg)
               end if
               
               ! Transform RRTMGP outputs to CAM outputs and compute heating rates.
               call set_lw_diags()

               if (write_output) then
                  call radiation_output_lw(lchnk, ncol, icall, rd, pbuf, cam_out)
               end if

            end if ! (active_calls(icall))
         end do    ! loop over diagnostic calcs (icall)
      end if  ! if (dolw)

      deallocate( &
         t_sfc, emis_sfc, toa_flux, sfac, t_rad, pmid_rad, pint_rad,  &
         t_day, pmid_day, pint_day, coszrs_day, alb_dir, alb_dif)

      !================!
      ! COSP simulator !
      !================!

      if (docosp) then

         emis(:,:) = 0._r8
         emis(:ncol,:) = 1._r8 - exp(-cld_lw_abs_cloudsim(:ncol,:))
         call outfld('EMIS', emis, pcols, lchnk)

         ! compute grid-box mean SW and LW snow optical depth for use by COSP
         gb_snow_tau(:,:) = 0._r8
         gb_snow_lw(:,:)  = 0._r8
         if (cldfsnow_idx > 0) then
            do i = 1, ncol
               do k = 1, pver
                  if (cldfsnow(i,k) > 0._r8) then

                     ! Add graupel to snow tau for cosp
                     if (cldfgrau_idx > 0 .and. graupel_in_rad) then
                        gb_snow_tau(i,k) = snow_tau_cloudsim(i,k)*cldfsnow(i,k) + &
                                           grau_tau_cloudsim(i,k)*cldfgrau(i,k)
                        gb_snow_lw(i,k)  = snow_lw_abs_cloudsim(i,k)*cldfsnow(i,k) + &
                                           grau_lw_abs_cloudsim(i,k)*cldfgrau(i,k)
                     else
                        gb_snow_tau(i,k) = snow_tau_cloudsim(i,k)*cldfsnow(i,k)
                        gb_snow_lw(i,k)  = snow_lw_abs_cloudsim(i,k)*cldfsnow(i,k)
                     end if
                  end if
               end do
            end do
         end if

         ! advance counter for this timestep (chunk dimension required for thread safety)
         cosp_cnt(lchnk) = cosp_cnt(lchnk) + 1

         ! if counter is the same as cosp_nradsteps, run cosp and reset counter
         if (cosp_nradsteps .eq. cosp_cnt(lchnk)) then

            ! N.B.: For snow optical properties, the GRID-BOX MEAN shortwave and longwave
            !       optical depths are passed.
            call cospsimulator_intr_run( &
               state,  pbuf, cam_in, emis, coszrs,                     &
               cld_swtau_in=cld_tau_cloudsim, snow_tau_in=gb_snow_tau, &
               snow_emis_in=gb_snow_lw)
            cosp_cnt(lchnk) = 0
         end if
      end if   ! docosp
   end if   ! if (dosw .or. dolw) then

   ! Calculate dry static energy if LW calc wasn't done; needed before calling radheat_run
   call rrtmgp_dry_static_energy_tendency_run(ncol, state%pdel, (.not. dosw), (.not. dolw), &
             qrs, qrl, errmsg, errflg)
   if (errflg /= 0) then
      call endrun(sub//': '//errmsg)
   end if

   ! Output for PORT: Parallel Offline Radiative Transport
   call rad_data_write(pbuf, state, cam_in, coszrs)

   ! Compute net radiative heating tendency.  Note that the WACCM version
   ! of radheat_tend merges upper atmosphere heating rates with those calculated
   ! by RRTMGP.
   call radheat_tend(state, pbuf,  ptend, qrl, qrs, fsns, &
                     fsnt, flns, flnt, cam_in%asdir, net_flx)

   if (write_output) then
      ! Compute heating rate for dtheta/dt
      do k = 1, pver
         do i = 1, ncol
            ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5_r8/state%pmid(i,k))**cappa
         end do
      end do
      call outfld('HR', ftem, pcols, lchnk)
   end if

   if (.not. present(rd_out)) then
      deallocate(rd)
   end if

   call rrtmgp_post_run(ncol, qrs, qrl, fsns, state%pdel, atm_optics_sw, cloud_sw, aer_sw, &
                  fsw, fswc, sources_lw, cloud_lw, aer_lw, flw, flwc, cam_out%netsw, errmsg, errflg)
   if (errflg /= 0) then
     call endrun(sub//': '//errmsg)
   end if

   !-------------------------------------------------------------------------------
   contains
   !-------------------------------------------------------------------------------

   subroutine set_sw_diags()

      ! Transform RRTMGP output for CAM and compute heating rates.
      ! SW fluxes from RRTMGP are on daylight columns only, so expand to
      ! full chunks for output to CAM history.

      integer :: i
      real(r8), dimension(size(fsw%fluxes%bnd_flux_dn,1), &
                          size(fsw%fluxes%bnd_flux_dn,2), &
                          size(fsw%fluxes%bnd_flux_dn,3)) :: flux_dn_diffuse
      !-------------------------------------------------------------------------

      ! Initialize to provide 0.0 values for night columns.
      fns               = 0._r8 ! net sw flux
      fcns              = 0._r8 ! net sw clearsky flux
      fsds              = 0._r8 ! downward sw flux at surface
      rd%fsdsc          = 0._r8 ! downward sw clearsky flux at surface
      rd%fsutoa         = 0._r8 ! upward sw flux at TOA
      rd%fsntoa         = 0._r8 ! net sw at TOA
      rd%fsntoac        = 0._r8 ! net sw clearsky flux at TOA
      rd%solin          = 0._r8 ! solar irradiance at TOA
      rd%flux_sw_up     = 0._r8
      rd%flux_sw_dn     = 0._r8
      rd%flux_sw_clr_up = 0._r8
      rd%flux_sw_clr_dn = 0._r8

      qrs      = 0._r8
      fsns     = 0._r8
      fsnt     = 0._r8
      rd%qrsc  = 0._r8
      rd%fsnsc = 0._r8
      rd%fsntc = 0._r8

      do i = 1, nday
         fns(idxday(i),ktopcam:)  = fsw%fluxes%flux_net(i, ktoprad:)
         fcns(idxday(i),ktopcam:) = fswc%fluxes%flux_net(i,ktoprad:)
         fsds(idxday(i))          = fsw%fluxes%flux_dn(i, nlay+1)
         rd%fsdsc(idxday(i))      = fswc%fluxes%flux_dn(i, nlay+1)
         rd%fsutoa(idxday(i))     = fsw%fluxes%flux_up(i, 1)
         rd%fsntoa(idxday(i))     = fsw%fluxes%flux_net(i, 1)
         rd%fsntoac(idxday(i))    = fswc%fluxes%flux_net(i, 1)
         rd%solin(idxday(i))      = fswc%fluxes%flux_dn(i, 1)
         rd%flux_sw_up(idxday(i),ktopcam:)     = fsw%fluxes%flux_up(i,ktoprad:)
         rd%flux_sw_dn(idxday(i),ktopcam:)     = fsw%fluxes%flux_dn(i,ktoprad:)
         rd%flux_sw_clr_up(idxday(i),ktopcam:) = fswc%fluxes%flux_up(i,ktoprad:)
         rd%flux_sw_clr_dn(idxday(i),ktopcam:) = fswc%fluxes%flux_dn(i,ktoprad:)
      end do

      ! Compute heating rate as a dry static energy tendency.
      call heating_rate('SW', ncol, fns, qrs)
      call heating_rate('SW', ncol, fcns, rd%qrsc)

      fsns(:ncol)     = fns(:ncol,pverp)    ! net sw flux at surface
      fsnt(:ncol)     = fns(:ncol,ktopcam)  ! net sw flux at top-of-model (w/o extra layer)
      rd%fsnsc(:ncol) = fcns(:ncol,pverp)   ! net sw clearsky flux at surface
      rd%fsntc(:ncol) = fcns(:ncol,ktopcam) ! net sw clearsky flux at top

      ! Output fluxes at 200 mb
      call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fns,  rd%fsn200)
      call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcns, rd%fsn200c)
      if (hist_fld_active('FSNR')) then
         do i = 1,ncol
            call vertinterp(1, 1, pverp, state%pint(i,:), p_trop(i), fns(i,:), rd%fsnr(i))
         end do
      end if

      if (spectralflux) then
         su  = 0._r8
         sd  = 0._r8
         do i = 1, nday
            su(idxday(i),ktopcam:,:) = fsw%fluxes%bnd_flux_up(i,ktoprad:,:)
            sd(idxday(i),ktopcam:,:) = fsw%fluxes%bnd_flux_dn(i,ktoprad:,:)
         end do
      end if

      ! Export surface fluxes
      ! sols(pcols)      Direct solar rad on surface (< 0.7)
      ! soll(pcols)      Direct solar rad on surface (>= 0.7)
      ! RRTMGP: Near-IR bands (1-10), 820-16000 cm-1, 0.625-12.195 microns
      ! Put half of band 10 in each of the UV/visible and near-IR values,
      ! since this band straddles 0.7 microns:
      ! UV/visible bands 10-13, 16000-50000 cm-1, 0.200-0.625 micron

      ! reset fluxes
      cam_out%sols  = 0.0_r8
      cam_out%soll  = 0.0_r8
      cam_out%solsd = 0.0_r8
      cam_out%solld = 0.0_r8

      ! Calculate diffuse flux from total and direct
      flux_dn_diffuse = fsw%fluxes%bnd_flux_dn - fsw%fluxes%bnd_flux_dn_dir

      do i = 1, nday
         cam_out%soll(idxday(i)) = sum(fsw%fluxes%bnd_flux_dn_dir(i,nlay+1,1:9))      &
                                   + 0.5_r8 * fsw%fluxes%bnd_flux_dn_dir(i,nlay+1,10)

         cam_out%sols(idxday(i)) = 0.5_r8 * fsw%fluxes%bnd_flux_dn_dir(i,nlay+1,10)   &
                                   + sum(fsw%fluxes%bnd_flux_dn_dir(i,nlay+1,11:14))

         cam_out%solld(idxday(i)) = sum(flux_dn_diffuse(i,nlay+1,1:9))         &
                                    + 0.5_r8 * flux_dn_diffuse(i,nlay+1,10)
         
         cam_out%solsd(idxday(i)) = 0.5_r8 * flux_dn_diffuse(i, nlay+1, 10)    &
                                    + sum(flux_dn_diffuse(i,nlay+1,11:14))
      end do

   end subroutine set_sw_diags

   !-------------------------------------------------------------------------------

   subroutine set_lw_diags()

      ! Set CAM LW diagnostics
      !----------------------------------------------------------------------------
 
      fnl = 0._r8
      fcnl = 0._r8

      ! RTE-RRTMGP convention for net is (down - up) **CAM assumes (up - down) !!
      fnl(:ncol,ktopcam:)  = -1._r8 * flw%fluxes%flux_net(    :, ktoprad:)
      fcnl(:ncol,ktopcam:) = -1._r8 * flwc%fluxes%flux_net(   :, ktoprad:)

      rd%flux_lw_up(:ncol,ktopcam:)     = flw%fluxes%flux_up( :, ktoprad:)
      rd%flux_lw_clr_up(:ncol,ktopcam:) = flwc%fluxes%flux_up(:, ktoprad:)
      rd%flux_lw_dn(:ncol,ktopcam:)     = flw%fluxes%flux_dn( :, ktoprad:)
      rd%flux_lw_clr_dn(:ncol,ktopcam:) = flwc%fluxes%flux_dn(:, ktoprad:)

      call heating_rate('LW', ncol, fnl, qrl)
      call heating_rate('LW', ncol, fcnl, rd%qrlc)

      flns(:ncol) = fnl(:ncol, pverp)
      flnt(:ncol) = fnl(:ncol, ktopcam)

      rd%flnsc(:ncol) = fcnl(:ncol, pverp)
      rd%flntc(:ncol) = fcnl(:ncol, ktopcam)    ! net lw flux at top-of-model

      cam_out%flwds(:ncol) = flw%fluxes%flux_dn(:, nlay+1)
      rd%fldsc(:ncol)      = flwc%fluxes%flux_dn(:, nlay+1)

      rd%flut(:ncol)  = flw%fluxes%flux_up(:, ktoprad)
      rd%flutc(:ncol) = flwc%fluxes%flux_up(:, ktoprad)

      ! Output fluxes at 200 mb
      call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fnl,  rd%fln200)
      call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcnl, rd%fln200c)
      if (hist_fld_active('FLNR')) then
         do i = 1,ncol
            call vertinterp(1, 1, pverp, state%pint(i,:), p_trop(i), fnl(i,:), rd%flnr(i))
         end do
      end if

      if (spectralflux) then
         lu  = 0._r8
         ld  = 0._r8
         lu(:ncol, ktopcam:, :)  = flw%fluxes%bnd_flux_up(:, ktoprad:, :)
         ld(:ncol, ktopcam:, :)  = flw%fluxes%bnd_flux_dn(:, ktoprad:, :)
      end if

   end subroutine set_lw_diags

   !-------------------------------------------------------------------------------

   subroutine heating_rate(type, ncol, flux_net, hrate)

      ! Compute heating rate as a dry static energy tendency
      
      ! arguments
      character(2), intent(in)  :: type ! either LW or SW
      integer,      intent(in)  :: ncol
      real(r8),     intent(in)  :: flux_net(pcols,pverp)  ! W/m^2
      real(r8),     intent(out) :: hrate(pcols,pver)      ! J/kg/s

      ! local vars
      integer :: k

      ! Initialize for layers where RRTMGP is not providing fluxes.
      hrate = 0.0_r8

      select case (type)
      case ('LW')

         do k = ktopcam, pver
            ! (flux divergence as bottom-MINUS-top) * g/dp
            hrate(:ncol,k) = (flux_net(:ncol,k+1) - flux_net(:ncol,k)) * &
                          gravit * state%rpdel(:ncol,k)
         end do

      case ('SW')

         do k = ktopcam, pver
            ! top - bottom
            hrate(:ncol,k) = (flux_net(:ncol,k) - flux_net(:ncol,k+1)) * &
                          gravit * state%rpdel(:ncol,k)
         end do

      end select

   end subroutine heating_rate

   !----------------------------------------------------------------------------
   !            -- end contains statement of radiation_tend --
   !----------------------------------------------------------------------------
end subroutine radiation_tend

!===============================================================================

subroutine radiation_output_sw(lchnk, ncol, icall, rd, pbuf, cam_out)

   ! Dump shortwave radiation information to history buffer.

   integer ,               intent(in) :: lchnk
   integer,                intent(in) :: ncol
   integer,                intent(in) :: icall
   type(rad_out_t),        intent(in) :: rd
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cam_out_t),        intent(in) :: cam_out

   ! local variables
   real(r8), pointer :: qrs(:,:)
   real(r8), pointer :: fsnt(:)
   real(r8), pointer :: fsns(:)
   real(r8), pointer :: fsds(:)

   real(r8) :: ftem(pcols)
   !----------------------------------------------------------------------------

   call pbuf_get_field(pbuf, qrs_idx,  qrs)
   call pbuf_get_field(pbuf, fsnt_idx, fsnt)
   call pbuf_get_field(pbuf, fsns_idx, fsns)
   call pbuf_get_field(pbuf, fsds_idx, fsds)

   call outfld('SOLIN'//diag(icall),    rd%solin,      pcols, lchnk)

   ! QRS is output as temperature tendency.
   call outfld('QRS'//diag(icall),      qrs(:ncol,:)/cpair,     ncol, lchnk)
   call outfld('QRSC'//diag(icall),     rd%qrsc(:ncol,:)/cpair, ncol, lchnk)

   call outfld('FSNT'//diag(icall),    fsnt,               pcols, lchnk)
   call outfld('FSNTC'//diag(icall),   rd%fsntc,           pcols, lchnk)
   call outfld('FSNTOA'//diag(icall),  rd%fsntoa,          pcols, lchnk)
   call outfld('FSNTOAC'//diag(icall), rd%fsntoac,         pcols, lchnk)

   ftem(:ncol) = rd%fsntoa(:ncol) - rd%fsntoac(:ncol)
   call outfld('SWCF'//diag(icall),     ftem,          pcols, lchnk)

   call outfld('FSUTOA'//diag(icall),   rd%fsutoa,     pcols, lchnk)

   call outfld('FSNIRTOA'//diag(icall), rd%fsnirt,     pcols, lchnk)
   call outfld('FSNRTOAC'//diag(icall), rd%fsnrtc,     pcols, lchnk)
   call outfld('FSNRTOAS'//diag(icall), rd%fsnirtsq,   pcols, lchnk)

   call outfld('FSN200'//diag(icall),   rd%fsn200,     pcols, lchnk)
   call outfld('FSN200C'//diag(icall),  rd%fsn200c,    pcols, lchnk)

   call outfld('FSNR'//diag(icall),     rd%fsnr,       pcols, lchnk)

   call outfld('SOLS'//diag(icall),     cam_out%sols,  pcols, lchnk)
   call outfld('SOLL'//diag(icall),     cam_out%soll,  pcols, lchnk)
   call outfld('SOLSD'//diag(icall),    cam_out%solsd, pcols, lchnk)
   call outfld('SOLLD'//diag(icall),    cam_out%solld, pcols, lchnk)

   call outfld('FSNS'//diag(icall),     fsns,          pcols, lchnk)
   call outfld('FSNSC'//diag(icall),    rd%fsnsc,      pcols, lchnk)

   call outfld('FSDS'//diag(icall),     fsds,          pcols, lchnk)
   call outfld('FSDSC'//diag(icall),    rd%fsdsc,      pcols, lchnk)

   call outfld('FUS'//diag(icall),  rd%flux_sw_up,     pcols, lchnk)
   call outfld('FUSC'//diag(icall), rd%flux_sw_clr_up, pcols, lchnk)
   call outfld('FDS'//diag(icall),  rd%flux_sw_dn,     pcols, lchnk)
   call outfld('FDSC'//diag(icall), rd%flux_sw_clr_dn, pcols, lchnk)

end subroutine radiation_output_sw

!===============================================================================

subroutine radiation_output_cld(lchnk, rd)

   ! Dump shortwave cloud optics information to history buffer.

   integer ,        intent(in) :: lchnk
   type(rad_out_t), intent(in) :: rd
   !----------------------------------------------------------------------------

   call outfld('TOT_CLD_VISTAU',  rd%tot_cld_vistau,  pcols, lchnk)
   call outfld('TOT_ICLD_VISTAU', rd%tot_icld_vistau, pcols, lchnk)
   call outfld('LIQ_ICLD_VISTAU', rd%liq_icld_vistau, pcols, lchnk)
   call outfld('ICE_ICLD_VISTAU', rd%ice_icld_vistau, pcols, lchnk)
   if (cldfsnow_idx > 0) then
      call outfld('SNOW_ICLD_VISTAU', rd%snow_icld_vistau, pcols, lchnk)
   endif
   if (cldfgrau_idx > 0 .and. graupel_in_rad) then
      call outfld('GRAU_ICLD_VISTAU', rd%grau_icld_vistau , pcols, lchnk)
   endif

end subroutine radiation_output_cld

!===============================================================================

subroutine radiation_output_lw(lchnk, ncol, icall, rd, pbuf, cam_out)

   ! Dump longwave radiation information to history buffer

   integer,                intent(in) :: lchnk
   integer,                intent(in) :: ncol
   integer,                intent(in) :: icall  ! icall=0 for climate diagnostics
   type(rad_out_t),        intent(in) :: rd
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cam_out_t),        intent(in) :: cam_out

   ! local variables
   real(r8), pointer :: qrl(:,:)
   real(r8), pointer :: flnt(:)
   real(r8), pointer :: flns(:)

   real(r8) :: ftem(pcols)
   !----------------------------------------------------------------------------

   call pbuf_get_field(pbuf, qrl_idx,  qrl)
   call pbuf_get_field(pbuf, flnt_idx, flnt)
   call pbuf_get_field(pbuf, flns_idx, flns)

   call outfld('QRL'//diag(icall),     qrl(:ncol,:)/cpair,     ncol, lchnk)
   call outfld('QRLC'//diag(icall),    rd%qrlc(:ncol,:)/cpair, ncol, lchnk)

   call outfld('FLNT'//diag(icall),    flnt,          pcols, lchnk)
   call outfld('FLNTC'//diag(icall),   rd%flntc,      pcols, lchnk)

   call outfld('FLUT'//diag(icall),    rd%flut,       pcols, lchnk)
   call outfld('FLUTC'//diag(icall),   rd%flutc,      pcols, lchnk)
   
   ftem(:ncol) = rd%flutc(:ncol) - rd%flut(:ncol)
   call outfld('LWCF'//diag(icall),    ftem,          pcols, lchnk)

   call outfld('FLN200'//diag(icall),  rd%fln200,     pcols, lchnk)
   call outfld('FLN200C'//diag(icall), rd%fln200c,    pcols, lchnk)

   call outfld('FLNR'//diag(icall),    rd%flnr,       pcols, lchnk)

   call outfld('FLNS'//diag(icall),    flns,          pcols, lchnk)
   call outfld('FLNSC'//diag(icall),   rd%flnsc,      pcols, lchnk)

   call outfld('FLDS'//diag(icall),    cam_out%flwds, pcols, lchnk)
   call outfld('FLDSC'//diag(icall),   rd%fldsc,      pcols, lchnk)

   call outfld('FDL'//diag(icall),  rd%flux_lw_dn,     pcols, lchnk)
   call outfld('FDLC'//diag(icall), rd%flux_lw_clr_dn, pcols, lchnk)
   call outfld('FUL'//diag(icall),  rd%flux_lw_up,     pcols, lchnk)
   call outfld('FULC'//diag(icall), rd%flux_lw_clr_up, pcols, lchnk)

end subroutine radiation_output_lw

!===============================================================================

subroutine coefs_init(coefs_file, available_gases, kdist)
   use rrtmgp_lw_gas_optics_data, only: rrtmgp_lw_gas_optics_data_init

   ! Read data from coefficients file.  Initialize the kdist object.
   ! available_gases object provides the gas names that CAM provides.

   ! arguments
   character(len=*),                  intent(in)  :: coefs_file
   class(ty_gas_concs_ccpp),          intent(in)  :: available_gases
   class(ty_gas_optics_rrtmgp_ccpp),  intent(inout) :: kdist

   ! local variables
   type(file_desc_t) :: fh    ! pio file handle
   character(len=cl) :: locfn ! path to file on local storage

   ! File dimensions
   integer ::            &
      absorber,          &
      atmos_layer,       &
      bnd,               &
      pressure,          &
      temperature,       &
      absorber_ext,      &
      pressure_interp,   &
      mixing_fraction,   &
      gpt,               &
      temperature_Planck
   
   integer :: i
   integer :: did, vid
   integer :: ierr, istat

   character(32), dimension(:),  allocatable :: gas_names
   integer,  dimension(:,:,:),   allocatable :: key_species
   integer,  dimension(:,:),     allocatable :: band2gpt
   real(r8), dimension(:,:),     allocatable :: band_lims_wavenum
   real(r8), dimension(:),       allocatable :: press_ref, temp_ref
   real(r8)                                  :: press_ref_trop, temp_ref_t, temp_ref_p
   real(r8), dimension(:,:,:),   allocatable :: vmr_ref
   real(r8), dimension(:,:,:,:), allocatable :: kmajor
   real(r8), dimension(:,:,:),   allocatable :: kminor_lower, kminor_upper
   real(r8), dimension(:,:),     allocatable :: totplnk
   real(r8), dimension(:,:,:,:), allocatable :: planck_frac
   real(r8), dimension(:),       allocatable :: solar_src_quiet, solar_src_facular, solar_src_sunspot
   real(r8)                                  :: tsi_default
   real(r8), dimension(:,:,:),   allocatable :: rayl_lower, rayl_upper
   character(len=32), dimension(:),  allocatable :: gas_minor,         &
                                                    identifier_minor,  &
                                                    minor_gases_lower, &
                                                    minor_gases_upper, &
                                                    scaling_gas_lower, &
                                                    scaling_gas_upper
   integer, dimension(:,:),          allocatable :: minor_limits_gpt_lower, &
                                                    minor_limits_gpt_upper
   ! Send these to RRTMGP as logicals,
   ! but they have to be read from the netCDF as integers
   logical, dimension(:),            allocatable :: minor_scales_with_density_lower, &
                                                    minor_scales_with_density_upper
   logical, dimension(:),            allocatable :: scale_by_complement_lower, &
                                                    scale_by_complement_upper
   integer, dimension(:), allocatable :: int2log   ! use this to convert integer-to-logical.
   integer, dimension(:),            allocatable :: kminor_start_lower, kminor_start_upper
   real(r8), dimension(:,:),         allocatable :: optimal_angle_fit
   real(r8)                                      :: mg_default, sb_default

   integer :: pairs, &
              minorabsorbers, &
              minor_absorber_intervals_lower, &
              minor_absorber_intervals_upper, &
              contributors_lower, &
              contributors_upper, &
              fit_coeffs

   character(len=128) :: error_msg
   character(len=512) :: errmsg
   character(len=*), parameter :: sub = 'coefs_init'
   !----------------------------------------------------------------------------

   ! Open file
   call getfil(coefs_file, locfn, 0)
   call cam_pio_openfile(fh, locfn, PIO_NOWRITE)

   call pio_seterrorhandling(fh, PIO_BCAST_ERROR)

   ! Get dimensions

   ierr = pio_inq_dimid(fh, 'absorber', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': absorber not found')
   ierr = pio_inq_dimlen(fh, did, absorber)

   ierr = pio_inq_dimid(fh, 'atmos_layer', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': atmos_layer not found')
   ierr = pio_inq_dimlen(fh, did, atmos_layer)

   ierr = pio_inq_dimid(fh, 'bnd', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': bnd not found')
   ierr = pio_inq_dimlen(fh, did, bnd)

   ierr = pio_inq_dimid(fh, 'pressure', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': pressure not found')
   ierr = pio_inq_dimlen(fh, did, pressure)

   ierr = pio_inq_dimid(fh, 'temperature', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': temperature not found')
   ierr = pio_inq_dimlen(fh, did, temperature)

   ierr = pio_inq_dimid(fh, 'absorber_ext', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': absorber_ext not found')
   ierr = pio_inq_dimlen(fh, did, absorber_ext)

   ierr = pio_inq_dimid(fh, 'pressure_interp', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': pressure_interp not found')
   ierr = pio_inq_dimlen(fh, did, pressure_interp)

   ierr = pio_inq_dimid(fh, 'mixing_fraction', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': mixing_fraction not found')
   ierr = pio_inq_dimlen(fh, did, mixing_fraction)
   
   ierr = pio_inq_dimid(fh, 'gpt', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': gpt not found')
   ierr = pio_inq_dimlen(fh, did, gpt)

   temperature_Planck = 0
   ierr = pio_inq_dimid(fh, 'temperature_Planck', did)
   if (ierr == PIO_NOERR) then
      ierr = pio_inq_dimlen(fh, did, temperature_Planck)
   end if
   ierr = pio_inq_dimid(fh, 'pair', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': pair not found')
   ierr = pio_inq_dimlen(fh, did, pairs)
   ierr = pio_inq_dimid(fh, 'minor_absorber', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_absorber not found')
   ierr = pio_inq_dimlen(fh, did, minorabsorbers)
   ierr = pio_inq_dimid(fh, 'minor_absorber_intervals_lower', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_absorber_intervals_lower not found')
   ierr = pio_inq_dimlen(fh, did, minor_absorber_intervals_lower)
   ierr = pio_inq_dimid(fh, 'minor_absorber_intervals_upper', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_absorber_intervals_upper not found')
   ierr = pio_inq_dimlen(fh, did, minor_absorber_intervals_upper)
   ierr = pio_inq_dimid(fh, 'contributors_lower', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': contributors_lower not found')
   ierr = pio_inq_dimlen(fh, did, contributors_lower)
   ierr = pio_inq_dimid(fh, 'contributors_upper', did)
   if (ierr /= PIO_NOERR) call endrun(sub//': contributors_upper not found')
   ierr = pio_inq_dimlen(fh, did, contributors_upper)

   ierr = pio_inq_dimid(fh, 'fit_coeffs', did)
   if (ierr == PIO_NOERR) then
      ierr = pio_inq_dimlen(fh, did, fit_coeffs)
   end if

   ! Get variables
   
   ! names of absorbing gases
   allocate(gas_names(absorber), stat=istat)
   call handle_allocate_error(istat, sub, 'gas_names')
   ierr = pio_inq_varid(fh, 'gas_names', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': gas_names not found')
   ierr = pio_get_var(fh, vid, gas_names)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading gas_names')

   ! key species pair for each band
   allocate(key_species(2,atmos_layer,bnd), stat=istat)
   call handle_allocate_error(istat, sub, 'key_species')
   ierr = pio_inq_varid(fh, 'key_species', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': key_species not found')
   ierr = pio_get_var(fh, vid, key_species)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading key_species')

   ! beginning and ending gpoint for each band
   allocate(band2gpt(2,bnd), stat=istat)
   call handle_allocate_error(istat, sub, 'band2gpt')
   ierr = pio_inq_varid(fh, 'bnd_limits_gpt', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': bnd_limits_gpt not found')
   ierr = pio_get_var(fh, vid, band2gpt)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading bnd_limits_gpt')

   ! beginning and ending wavenumber for each band
   allocate(band_lims_wavenum(2,bnd), stat=istat)
   call handle_allocate_error(istat, sub, 'band_lims_wavenum')
   ierr = pio_inq_varid(fh, 'bnd_limits_wavenumber', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': bnd_limits_wavenumber not found')
   ierr = pio_get_var(fh, vid, band_lims_wavenum)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading bnd_limits_wavenumber')

   ! pressures [hPa] for reference atmosphere; press_ref(# reference layers)
   allocate(press_ref(pressure), stat=istat)
   call handle_allocate_error(istat, sub, 'press_ref')
   ierr = pio_inq_varid(fh, 'press_ref', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': press_ref not found')
   ierr = pio_get_var(fh, vid, press_ref)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading press_ref')

   ! reference pressure separating the lower and upper atmosphere
   ierr = pio_inq_varid(fh, 'press_ref_trop', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': press_ref_trop not found')
   ierr = pio_get_var(fh, vid, press_ref_trop)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading press_ref_trop')

   ! temperatures [K] for reference atmosphere; temp_ref(# reference layers)
   allocate(temp_ref(temperature), stat=istat)
   call handle_allocate_error(istat, sub, 'temp_ref')
   ierr = pio_inq_varid(fh, 'temp_ref', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': temp_ref not found')
   ierr = pio_get_var(fh, vid, temp_ref)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading temp_ref')

   ! standard spectroscopic reference temperature [K]
   ierr = pio_inq_varid(fh, 'absorption_coefficient_ref_T', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': absorption_coefficient_ref_T not found')
   ierr = pio_get_var(fh, vid, temp_ref_t)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading absorption_coefficient_ref_T')

   ! standard spectroscopic reference pressure [hPa]
   ierr = pio_inq_varid(fh, 'absorption_coefficient_ref_P', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': absorption_coefficient_ref_P not found')
   ierr = pio_get_var(fh, vid, temp_ref_p)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading absorption_coefficient_ref_P')

   ! volume mixing ratios for reference atmosphere
   allocate(vmr_ref(atmos_layer, absorber_ext, temperature), stat=istat)
   call handle_allocate_error(istat, sub, 'vmr_ref')
   ierr = pio_inq_varid(fh, 'vmr_ref', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': vmr_ref not found')
   ierr = pio_get_var(fh, vid, vmr_ref)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading vmr_ref')

   ! absorption coefficients due to major absorbing gases
   allocate(kmajor(gpt,mixing_fraction,pressure_interp,temperature), stat=istat)
   call handle_allocate_error(istat, sub, 'kmajor')
   ierr = pio_inq_varid(fh, 'kmajor', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kmajor not found')
   ierr = pio_get_var(fh, vid, kmajor)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kmajor')

   ! absorption coefficients due to minor absorbing gases in lower part of atmosphere
   allocate(kminor_lower(contributors_lower, mixing_fraction, temperature), stat=istat)
   call handle_allocate_error(istat, sub, 'kminor_lower')
   ierr = pio_inq_varid(fh, 'kminor_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kminor_lower not found')
   ierr = pio_get_var(fh, vid, kminor_lower)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kminor_lower')

   ! absorption coefficients due to minor absorbing gases in upper part of atmosphere
   allocate(kminor_upper(contributors_upper, mixing_fraction, temperature), stat=istat)
   call handle_allocate_error(istat, sub, 'kminor_upper')
   ierr = pio_inq_varid(fh, 'kminor_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kminor_upper not found')
   ierr = pio_get_var(fh, vid, kminor_upper)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kminor_upper')

   ! integrated Planck function by band
   ierr = pio_inq_varid(fh, 'totplnk', vid)
   if (ierr == PIO_NOERR) then
      allocate(totplnk(temperature_Planck,bnd), stat=istat)
      call handle_allocate_error(istat, sub, 'totplnk')
      ierr = pio_get_var(fh, vid, totplnk)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading totplnk')
   end if

   ! Planck fractions
   ierr = pio_inq_varid(fh, 'plank_fraction', vid)
   if (ierr == PIO_NOERR) then
      allocate(planck_frac(gpt,mixing_fraction,pressure_interp,temperature), stat=istat)
      call handle_allocate_error(istat, sub, 'planck_frac')
      ierr = pio_get_var(fh, vid, planck_frac)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading plank_fraction')
   end if

   ierr = pio_inq_varid(fh, 'optimal_angle_fit', vid)
   if (ierr == PIO_NOERR) then
      allocate(optimal_angle_fit(fit_coeffs, bnd), stat=istat)
      call handle_allocate_error(istat, sub, 'optiman_angle_fit')
      ierr = pio_get_var(fh, vid, optimal_angle_fit)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading optimal_angle_fit')
   end if

   ierr = pio_inq_varid(fh, 'solar_source_quiet', vid)
   if (ierr == PIO_NOERR) then
      allocate(solar_src_quiet(gpt), stat=istat)
      call handle_allocate_error(istat, sub, 'solar_src_quiet')
      ierr = pio_get_var(fh, vid, solar_src_quiet)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading solar_source_quiet')
   end if

   ierr = pio_inq_varid(fh, 'solar_source_facular', vid)
   if (ierr == PIO_NOERR) then
      allocate(solar_src_facular(gpt), stat=istat)
      call handle_allocate_error(istat, sub, 'solar_src_facular')
      ierr = pio_get_var(fh, vid, solar_src_facular)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading solar_source_facular')
   end if

   ierr = pio_inq_varid(fh, 'solar_source_sunspot', vid)
   if (ierr == PIO_NOERR) then
      allocate(solar_src_sunspot(gpt), stat=istat)
      call handle_allocate_error(istat, sub, 'solar_src_sunspot')
      ierr = pio_get_var(fh, vid, solar_src_sunspot)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading solar_source_sunspot')
   end if

   ierr = pio_inq_varid(fh, 'tsi_default', vid)
   if (ierr == PIO_NOERR) then
      ierr = pio_get_var(fh, vid, tsi_default)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading tsi_default')
   end if

   ierr = pio_inq_varid(fh, 'mg_default', vid)
   if (ierr == PIO_NOERR) then
      ierr = pio_get_var(fh, vid, mg_default)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading mg_default')
   end if

   ierr = pio_inq_varid(fh, 'sb_default', vid)
   if (ierr == PIO_NOERR) then
      ierr = pio_get_var(fh, vid, sb_default)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading sb_default')
   end if

   ! rayleigh scattering contribution in lower part of atmosphere
   ierr = pio_inq_varid(fh, 'rayl_lower', vid)
   if (ierr == PIO_NOERR) then
      allocate(rayl_lower(gpt,mixing_fraction,temperature), stat=istat)
      call handle_allocate_error(istat, sub, 'rayl_lower')
      ierr = pio_get_var(fh, vid, rayl_lower)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading rayl_lower')
   end if

   ! rayleigh scattering contribution in upper part of atmosphere
   ierr = pio_inq_varid(fh, 'rayl_upper', vid)
   if (ierr == PIO_NOERR) then
      allocate(rayl_upper(gpt,mixing_fraction,temperature), stat=istat)
      call handle_allocate_error(istat, sub, 'rayl_upper')
      ierr = pio_get_var(fh, vid, rayl_upper)
      if (ierr /= PIO_NOERR) call endrun(sub//': error reading rayl_upper')
   end if

   allocate(gas_minor(minorabsorbers), stat=istat)
   call handle_allocate_error(istat, sub, 'gas_minor')
   ierr = pio_inq_varid(fh, 'gas_minor', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': gas_minor not found')
   ierr = pio_get_var(fh, vid, gas_minor)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading gas_minor')

   allocate(identifier_minor(minorabsorbers), stat=istat)
   call handle_allocate_error(istat, sub, 'identifier_minor')
   ierr = pio_inq_varid(fh, 'identifier_minor', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': identifier_minor not found')
   ierr = pio_get_var(fh, vid, identifier_minor)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading identifier_minor')
   
   allocate(minor_gases_lower(minor_absorber_intervals_lower), stat=istat)
   call handle_allocate_error(istat, sub, 'minor_gases_lower')
   ierr = pio_inq_varid(fh, 'minor_gases_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_gases_lower not found')
   ierr = pio_get_var(fh, vid, minor_gases_lower)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading minor_gases_lower')

   allocate(minor_gases_upper(minor_absorber_intervals_upper), stat=istat)
   call handle_allocate_error(istat, sub, 'minor_gases_upper')
   ierr = pio_inq_varid(fh, 'minor_gases_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_gases_upper not found')
   ierr = pio_get_var(fh, vid, minor_gases_upper)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading minor_gases_upper')

   allocate(minor_limits_gpt_lower(pairs,minor_absorber_intervals_lower), stat=istat)
   call handle_allocate_error(istat, sub, 'minor_limits_gpt_lower')
   ierr = pio_inq_varid(fh, 'minor_limits_gpt_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_limits_gpt_lower not found')
   ierr = pio_get_var(fh, vid, minor_limits_gpt_lower)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading minor_limits_gpt_lower')

   allocate(minor_limits_gpt_upper(pairs,minor_absorber_intervals_upper), stat=istat)
   call handle_allocate_error(istat, sub, 'minor_limits_gpt_upper')
   ierr = pio_inq_varid(fh, 'minor_limits_gpt_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_limits_gpt_upper not found')
   ierr = pio_get_var(fh, vid, minor_limits_gpt_upper)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading minor_limits_gpt_upper')

   ! Read as integer and convert to logical
   allocate(int2log(minor_absorber_intervals_lower), stat=istat)
   call handle_allocate_error(istat, sub, 'int2log for lower')

   allocate(minor_scales_with_density_lower(minor_absorber_intervals_lower), stat=istat)
   call handle_allocate_error(istat, sub, 'minor_scales_with_density_lower')
   ierr = pio_inq_varid(fh, 'minor_scales_with_density_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_scales_with_density_lower not found')
   ierr = pio_get_var(fh, vid, int2log)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading minor_scales_with_density_lower')
   do i = 1,minor_absorber_intervals_lower
      if (int2log(i) .eq. 0) then
         minor_scales_with_density_lower(i) = .false.
      else
         minor_scales_with_density_lower(i) = .true.
      end if
   end do

   allocate(scale_by_complement_lower(minor_absorber_intervals_lower), stat=istat)
   call handle_allocate_error(istat, sub, 'scale_by_complement_lower')
   ierr = pio_inq_varid(fh, 'scale_by_complement_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': scale_by_complement_lower not found')
   ierr = pio_get_var(fh, vid, int2log)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading scale_by_complement_lower')
   do i = 1,minor_absorber_intervals_lower
      if (int2log(i) .eq. 0) then
         scale_by_complement_lower(i) = .false.
      else
         scale_by_complement_lower(i) = .true.
      end if
   end do

   deallocate(int2log)

   ! Read as integer and convert to logical
   allocate(int2log(minor_absorber_intervals_upper), stat=istat)
   call handle_allocate_error(istat, sub, 'int2log for upper')

   allocate(minor_scales_with_density_upper(minor_absorber_intervals_upper), stat=istat)
   call handle_allocate_error(istat, sub, 'minor_scales_with_density_upper')
   ierr = pio_inq_varid(fh, 'minor_scales_with_density_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': minor_scales_with_density_upper not found')
   ierr = pio_get_var(fh, vid, int2log)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading minor_scales_with_density_upper')
   do i = 1,minor_absorber_intervals_upper
      if (int2log(i) .eq. 0) then
         minor_scales_with_density_upper(i) = .false.
      else
         minor_scales_with_density_upper(i) = .true.
      end if
   end do

   allocate(scale_by_complement_upper(minor_absorber_intervals_upper), stat=istat)
   call handle_allocate_error(istat, sub, 'scale_by_complement_upper')
   ierr = pio_inq_varid(fh, 'scale_by_complement_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': scale_by_complement_upper not found')
   ierr = pio_get_var(fh, vid, int2log)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading scale_by_complement_upper')
   do i = 1,minor_absorber_intervals_upper
      if (int2log(i) .eq. 0) then
         scale_by_complement_upper(i) = .false.
      else
         scale_by_complement_upper(i) = .true.
      end if
   end do

   deallocate(int2log)

   allocate(scaling_gas_lower(minor_absorber_intervals_lower), stat=istat)
   call handle_allocate_error(istat, sub, 'scaling_gas_lower')
   ierr = pio_inq_varid(fh, 'scaling_gas_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': scaling_gas_lower not found')
   ierr = pio_get_var(fh, vid, scaling_gas_lower)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading scaling_gas_lower')

   allocate(scaling_gas_upper(minor_absorber_intervals_upper), stat=istat)
   call handle_allocate_error(istat, sub, 'scaling_gas_upper')
   ierr = pio_inq_varid(fh, 'scaling_gas_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': scaling_gas_upper not found')
   ierr = pio_get_var(fh, vid, scaling_gas_upper)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading scaling_gas_upper')

   allocate(kminor_start_lower(minor_absorber_intervals_lower), stat=istat)
   call handle_allocate_error(istat, sub, 'kminor_start_lower')
   ierr = pio_inq_varid(fh, 'kminor_start_lower', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kminor_start_lower not found')
   ierr = pio_get_var(fh, vid, kminor_start_lower)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kminor_start_lower')

   allocate(kminor_start_upper(minor_absorber_intervals_upper), stat=istat)
   call handle_allocate_error(istat, sub, 'kminor_start_upper')
   ierr = pio_inq_varid(fh, 'kminor_start_upper', vid)
   if (ierr /= PIO_NOERR) call endrun(sub//': kminor_start_upper not found')
   ierr = pio_get_var(fh, vid, kminor_start_upper)
   if (ierr /= PIO_NOERR) call endrun(sub//': error reading kminor_start_upper')

   ! Close file
   call pio_closefile(fh)

   ! Initialize the gas optics object with data. The calls are slightly different depending
   ! on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)

   if (allocated(totplnk) .and. allocated(planck_frac)) then
      call rrtmgp_lw_gas_optics_data_init(kdist, available_gases, gas_names,           &
                  key_species, band2gpt, band_lims_wavenum, press_ref, press_ref_trop, &
                  temp_ref, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower,     &
                  kminor_upper, gas_minor, identifier_minor, minor_gases_lower,        &
                  minor_gases_upper, minor_limits_gpt_lower, minor_limits_gpt_upper,   &
                  minor_scales_with_density_lower, minor_scales_with_density_upper,    &
                  scaling_gas_lower, scaling_gas_upper, scale_by_complement_lower,     &
                  scale_by_complement_upper, kminor_start_lower, kminor_start_upper,   &
                  totplnk, planck_frac, rayl_lower, rayl_upper, optimal_angle_fit,     &
                  errmsg, ierr)
      if (ierr /= 0) then
         call endrun(sub//': ERROR message: '//errmsg)
      end if
   else if (allocated(solar_src_quiet)) then
      error_msg = kdist%gas_props%load( &
         available_gases%gas_concs, gas_names, key_species,     &
         band2gpt, band_lims_wavenum,                           &
         press_ref, press_ref_trop, temp_ref,                   &
         temp_ref_p, temp_ref_t, vmr_ref,                       & 
         kmajor, kminor_lower, kminor_upper,                    &
         gas_minor, identifier_minor,                           &
         minor_gases_lower, minor_gases_upper,                  &
         minor_limits_gpt_lower, minor_limits_gpt_upper,        &
         minor_scales_with_density_lower,                       &
         minor_scales_with_density_upper,                       &
         scaling_gas_lower, scaling_gas_upper,                  &
         scale_by_complement_lower,                             &
         scale_by_complement_upper,                             &
         kminor_start_lower,                                    &
         kminor_start_upper,                                    &
         solar_src_quiet, solar_src_facular, solar_src_sunspot, &
         tsi_default, mg_default, sb_default,                   &
         rayl_lower, rayl_upper)
   else
      error_msg = 'must supply either totplnk and planck_frac, or solar_src_[*]'
   end if

   call stop_on_err(error_msg, sub, 'kdist%gas_props%load')

   deallocate( &
      gas_names, key_species,               &
      band2gpt, band_lims_wavenum,          &
      press_ref, temp_ref, vmr_ref,         &
      kmajor, kminor_lower, kminor_upper,   &
      gas_minor, identifier_minor,          &
      minor_gases_lower, minor_gases_upper, &
      minor_limits_gpt_lower,               & 
      minor_limits_gpt_upper,               &
      minor_scales_with_density_lower,      &
      minor_scales_with_density_upper,      &
      scale_by_complement_lower,            & 
      scale_by_complement_upper,            &
      scaling_gas_lower, scaling_gas_upper, & 
      kminor_start_lower, kminor_start_upper)

   if (allocated(totplnk))           deallocate(totplnk)
   if (allocated(planck_frac))       deallocate(planck_frac)
   if (allocated(optimal_angle_fit)) deallocate(optimal_angle_fit)
   if (allocated(solar_src_quiet))   deallocate(solar_src_quiet)
   if (allocated(solar_src_facular)) deallocate(solar_src_facular)
   if (allocated(solar_src_sunspot)) deallocate(solar_src_sunspot)
   if (allocated(rayl_lower))        deallocate(rayl_lower)
   if (allocated(rayl_upper))        deallocate(rayl_upper)

end subroutine coefs_init

!=========================================================================================

subroutine stop_on_err(errmsg, sub, info)

! call endrun if RRTMGP function returns non-empty error message.

   character(len=*), intent(in) :: errmsg    ! return message from RRTMGP function
   character(len=*), intent(in) :: sub       ! name of calling subroutine
   character(len=*), intent(in) :: info      ! name of called function

   if (len_trim(errmsg) > 0) then
      call endrun(trim(sub)//': ERROR: '//trim(info)//': '//trim(errmsg))
   end if

end subroutine stop_on_err

!=========================================================================================

end module radiation

