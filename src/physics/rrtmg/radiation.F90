module radiation

!---------------------------------------------------------------------------------
!
! CAM interface to RRTMG radiation parameterization
!
!---------------------------------------------------------------------------------

use shr_kind_mod,        only: r8=>shr_kind_r8
use spmd_utils,          only: masterproc
use ppgrid,              only: pcols, pver, pverp, begchunk, endchunk
use physics_types,       only: physics_state, physics_ptend
use physics_buffer,      only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
use camsrfexch,          only: cam_out_t, cam_in_t
use physconst,           only: cappa, cpair

use time_manager,        only: get_nstep, is_first_restart_step, &
                               get_curr_calday, get_step_size

use rad_constituents,    only: N_DIAG, rad_cnst_get_call_list, rad_cnst_get_info, &
                               rad_cnst_get_gas, rad_cnst_out, oldcldoptics, &
                               liqcldoptics, icecldoptics

use radconstants,        only: nswbands, nlwbands, rrtmg_sw_cloudsim_band, rrtmg_lw_cloudsim_band, &
                               idx_sw_diag

use cospsimulator_intr,  only: docosp, cospsimulator_intr_init, &
                               cospsimulator_intr_run, cosp_nradsteps

use scamMod,             only: scm_crm_mode, single_column, have_cld, cldobs

use cam_history,         only: addfld, add_default, horiz_only, outfld, hist_fld_active
use cam_history_support, only: fillvalue

use pio,                 only: file_desc_t, var_desc_t,               &
                               pio_int, pio_noerr,                    &
                               pio_seterrorhandling, pio_bcast_error, &
                               pio_inq_varid, pio_def_var,            &
                               pio_put_var, pio_get_var

use cam_abortutils,      only: endrun
use error_messages,      only: handle_err
use perf_mod,            only: t_startf, t_stopf
use cam_logfile,         only: iulog

implicit none
private
save

public :: &
   radiation_readnl,         &! read namelist variables
   radiation_register,       &! registers radiation physics buffer fields
   radiation_nextsw_cday,    &! calendar day of next radiation calculation
   radiation_do,             &! query which radiation calcs are done this timestep
   radiation_init,           &! initialization
   radiation_define_restart, &! define variables for restart
   radiation_write_restart,  &! write variables to restart
   radiation_read_restart,   &! read variables from restart
   radiation_tend,           &! compute heating rates and fluxes
   rad_out_t                  ! type for diagnostic outputs

integer,public, allocatable :: cosp_cnt(:)       ! counter for cosp
integer,public              :: cosp_cnt_init = 0 !initial value for cosp counter

type rad_out_t
   real(r8) :: solin(pcols)         ! Solar incident flux

   real(r8) :: qrsc(pcols,pver)

   real(r8) :: fsntc(pcols)         ! Clear sky total column abs solar flux
   real(r8) :: fsntoa(pcols)        ! Net solar flux at TOA
   real(r8) :: fsntoac(pcols)       ! Clear sky net solar flux at TOA
   real(r8) :: fsutoa(pcols)        ! upwelling solar flux at TOA

   real(r8) :: fsnirt(pcols)        ! Near-IR flux absorbed at toa
   real(r8) :: fsnrtc(pcols)        ! Clear sky near-IR flux absorbed at toa
   real(r8) :: fsnirtsq(pcols)      ! Near-IR flux absorbed at toa >= 0.7 microns

   real(r8) :: fsn200(pcols)        ! fns interpolated to 200 mb
   real(r8) :: fsn200c(pcols)       ! fcns interpolated to 200 mb
   real(r8) :: fsnr(pcols)          ! fns interpolated to tropopause

   real(r8) :: fsnsc(pcols)         ! Clear sky surface abs solar flux
   real(r8) :: fsdsc(pcols)         ! Clear sky surface downwelling solar flux

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

   real(r8) :: cld_tau_cloudsim(pcols,pver)
   real(r8) :: aer_tau400(pcols,0:pver)
   real(r8) :: aer_tau550(pcols,0:pver)
   real(r8) :: aer_tau700(pcols,0:pver)
   
end type rad_out_t

! Namelist variables

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
logical :: use_rad_uniform_angle = .false. ! if true, use the namelist rad_uniform_angle for the coszrs calculation

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
integer :: cldfsnow_idx = 0 
integer :: cld_idx      = 0 

character(len=4) :: diag(0:N_DIAG) =(/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ','_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

! averaging time interval for zenith angle
real(r8) :: dt_avg = 0._r8

real(r8) :: rad_uniform_angle = -99._r8

! PIO descriptors (for restarts)
type(var_desc_t) :: cospcnt_desc

!===============================================================================
contains
!===============================================================================

subroutine radiation_readnl(nlfile)

   ! Read radiation_nl namelist group.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_logical, mpi_real8

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   integer :: dtime      ! timestep size
   character(len=*), parameter :: sub = 'radiation_readnl'

   namelist /radiation_nl/ iradsw, iradlw, irad_always, &
                           use_rad_dt_cosz, spectralflux, use_rad_uniform_angle, rad_uniform_angle
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'radiation_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, radiation_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   ! Broadcast namelist variables
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
   call mpi_bcast(rad_uniform_angle, 1, mpi_real8,  mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: rad_uniform_angle")

   if (use_rad_uniform_angle .and. rad_uniform_angle == -99._r8) then
      call endrun(sub // ' ERROR - use_rad_uniform_angle is set to .true, but rad_uniform_angle is not set ')
   end if

   ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
   dtime  = get_step_size()
   if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
   if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
   if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

   !----------------------------------------------------------------------- 
   ! Print runtime options to log.
   !-----------------------------------------------------------------------

   if (masterproc) then
      write(iulog,*) 'RRTMG radiation scheme parameters:'
      write(iulog,10) iradsw, iradlw, irad_always, use_rad_dt_cosz, spectralflux
   end if

10 format('  Frequency (timesteps) of Shortwave Radiation calc:  ',i5/, &
          '  Frequency (timesteps) of Longwave Radiation calc:   ',i5/, &
          '  SW/LW calc done every timestep for first N steps. N=',i5/, &
          '  Use average zenith angle:                           ',l5/, &
          '  Output spectrally resolved fluxes:                  ',l5/)

end subroutine radiation_readnl

!================================================================================================

subroutine radiation_register

   ! Register radiation fields in the physics buffer

   use physics_buffer, only: pbuf_add_field, dtype_r8
   use radiation_data, only: rad_data_register

   call pbuf_add_field('QRS' , 'global',dtype_r8,(/pcols,pver/), qrs_idx) ! shortwave radiative heating rate 
   call pbuf_add_field('QRL' , 'global',dtype_r8,(/pcols,pver/), qrl_idx) ! longwave  radiative heating rate 

   call pbuf_add_field('FSDS' , 'global',dtype_r8,(/pcols/), fsds_idx) ! Surface solar downward flux
   call pbuf_add_field('FSNS' , 'global',dtype_r8,(/pcols/), fsns_idx) ! Surface net shortwave flux
   call pbuf_add_field('FSNT' , 'global',dtype_r8,(/pcols/), fsnt_idx) ! Top-of-model net shortwave flux

   call pbuf_add_field('FLNS' , 'global',dtype_r8,(/pcols/), flns_idx) ! Surface net longwave flux
   call pbuf_add_field('FLNT' , 'global',dtype_r8,(/pcols/), flnt_idx) ! Top-of-model net longwave flux

   ! If the namelist has been configured for preserving the spectral fluxes, then create
   ! physics buffer variables to store the results.
   if (spectralflux) then
      call pbuf_add_field('SU'  , 'global',dtype_r8,(/pcols,pverp,nswbands/), su_idx) ! shortwave upward flux (per band)
      call pbuf_add_field('SD'  , 'global',dtype_r8,(/pcols,pverp,nswbands/), sd_idx) ! shortwave downward flux (per band)
      call pbuf_add_field('LU'  , 'global',dtype_r8,(/pcols,pverp,nlwbands/), lu_idx) ! longwave upward flux (per band)
      call pbuf_add_field('LD'  , 'global',dtype_r8,(/pcols,pverp,nlwbands/), ld_idx) ! longwave downward flux (per band)
   end if

   call rad_data_register()

end subroutine radiation_register

!================================================================================================

function radiation_do(op, timestep)

   ! Return true if the specified operation is done this timestep.

   character(len=*), intent(in) :: op             ! name of operation
   integer, intent(in), optional:: timestep
   logical                      :: radiation_do   ! return value

   ! Local variables
   integer :: nstep             ! current timestep number
   !-----------------------------------------------------------------------

   if (present(timestep)) then
      nstep = timestep
   else
      nstep = get_nstep()
   end if

   select case (op)

   case ('sw') ! do a shortwave heating calc this timestep?
      radiation_do = nstep == 0  .or.  iradsw == 1                     &
                    .or. (mod(nstep-1,iradsw) == 0  .and.  nstep /= 1) &
                    .or. nstep <= irad_always

   case ('lw') ! do a longwave heating calc this timestep?
      radiation_do = nstep == 0  .or.  iradlw == 1                     &
                    .or. (mod(nstep-1,iradlw) == 0  .and.  nstep /= 1) &
                    .or. nstep <= irad_always

   case default
      call endrun('radiation_do: unknown operation:'//op)

   end select
end function radiation_do

!================================================================================================

real(r8) function radiation_nextsw_cday()
  
   ! Return calendar day of next sw radiation calculation

   ! Local variables
   integer :: nstep      ! timestep counter
   logical :: dosw       ! true => do shosrtwave calc   
   integer :: offset     ! offset for calendar day calculation
   integer :: dTime      ! integer timestep size
   real(r8):: calday     ! calendar day of 
   !-----------------------------------------------------------------------

   radiation_nextsw_cday = -1._r8
   dosw   = .false.
   nstep  = get_nstep()
   dtime  = get_step_size()
   offset = 0
   do while (.not. dosw)
      nstep = nstep + 1
      offset = offset + dtime
      if (radiation_do('sw', nstep)) then
         radiation_nextsw_cday = get_curr_calday(offset=offset) 
         dosw = .true.
      end if
   end do
   if(radiation_nextsw_cday == -1._r8) then
      call endrun('error in radiation_nextsw_cday')
   end if
        
end function radiation_nextsw_cday

!================================================================================================

subroutine radiation_init(pbuf2d)

   ! Initialize the radiation parameterization, add fields to the history buffer

   use physics_buffer,  only: pbuf_get_index, pbuf_set_field
   use phys_control,    only: phys_getopts
   use radsw,           only: radsw_init
   use radlw,           only: radlw_init
   use rad_solar_var,   only: rad_solar_var_init
   use radiation_data,  only: rad_data_init
   use cloud_rad_props, only: cloud_rad_props_init
   use modal_aer_opt,   only: modal_aer_opt_init
   use rrtmg_state,     only: rrtmg_state_init
   use time_manager,    only: is_first_step


   ! arguments
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local variables
   integer :: icall, nmodes
   logical :: active_calls(0:N_DIAG)
   integer :: nstep                       ! current timestep number
   logical :: history_amwg                ! output the variables used by the AMWG diag package
   logical :: history_vdiag               ! output the variables used by the AMWG variability diag package
   logical :: history_budget              ! output tendencies and state variables for CAM4
                                          ! temperature, water vapor, cloud ice and cloud
                                          ! liquid budgets.
   integer :: history_budget_histfile_num ! output history file number for budget fields
   integer :: err

   integer :: dtime
   !-----------------------------------------------------------------------
    
   call rad_solar_var_init()
   call rrtmg_state_init()
   call rad_data_init(pbuf2d) ! initialize output fields for offline driver
   call radsw_init()
   call radlw_init()
   call cloud_rad_props_init()

   cld_idx      = pbuf_get_index('CLD')
   cldfsnow_idx = pbuf_get_index('CLDFSNOW',errcode=err)

   if (is_first_step()) then
      call pbuf_set_field(pbuf2d, qrl_idx, 0._r8)
   end if

   ! Set the radiation timestep for cosz calculations if requested using the adjusted iradsw value from radiation
   if (use_rad_dt_cosz)  then
      dtime  = get_step_size()
      dt_avg = iradsw*dtime
   end if

   call phys_getopts(history_amwg_out   = history_amwg,    &
                     history_vdiag_out  = history_vdiag,   &
                     history_budget_out = history_budget,  &
                     history_budget_histfile_num_out = history_budget_histfile_num)

   ! Determine whether modal aerosols are affecting the climate, and if so
   ! then initialize the modal aerosol optics module
   call rad_cnst_get_info(0, nmodes=nmodes)
   if (nmodes > 0) call modal_aer_opt_init()

   ! "irad_always" is number of time steps to execute radiation continuously from start of
   ! initial OR restart run
   nstep = get_nstep()
   if (irad_always > 0) then
      nstep       = get_nstep()
      irad_always = irad_always + nstep
   end if

   if (docosp) call cospsimulator_intr_init
    
   allocate(cosp_cnt(begchunk:endchunk))
   if (is_first_restart_step()) then
      cosp_cnt(begchunk:endchunk) = cosp_cnt_init
   else
      cosp_cnt(begchunk:endchunk) = 0     
   end if

   call addfld('O3colAbove',    horiz_only,   'A', 'DU', 'Column O3 above model top', sampling_seq='rad_lwsw')

   call addfld('TOT_CLD_VISTAU',  (/ 'lev' /), 'A',   '1', 'Total gbx cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('TOT_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Total in-cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('LIQ_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Liquid in-cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
   call addfld('ICE_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Ice in-cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)

   if (cldfsnow_idx > 0) then
      call addfld('SNOW_ICLD_VISTAU', (/ 'lev' /), 'A', '1', 'Snow in-cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
   endif

   ! get list of active radiation calls
   call rad_cnst_get_call_list(active_calls)

   ! Add shortwave radiation fields to history master field list.

   do icall = 0, N_DIAG

      if (active_calls(icall)) then

         call addfld('SOLIN'//diag(icall),    horiz_only,   'A', 'W/m2', 'Solar insolation', sampling_seq='rad_lwsw')

         call addfld('QRS'//diag(icall),      (/ 'lev' /),  'A', 'K/s',  'Solar heating rate', sampling_seq='rad_lwsw')
         call addfld('QRSC'//diag(icall),     (/ 'lev' /),  'A', 'K/s',  'Clearsky solar heating rate',                     &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNT'//diag(icall),     horiz_only,   'A', 'W/m2', 'Net solar flux at top of model',                  &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNTC'//diag(icall),    horiz_only,   'A', 'W/m2', 'Clearsky net solar flux at top of model',         &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNTOA'//diag(icall),   horiz_only,   'A', 'W/m2', 'Net solar flux at top of atmosphere',             &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNTOAC'//diag(icall),  horiz_only,   'A', 'W/m2', 'Clearsky net solar flux at top of atmosphere',    &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('SWCF'//diag(icall),     horiz_only,   'A', 'W/m2', 'Shortwave cloud forcing',                         &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSUTOA'//diag(icall),   horiz_only,   'A', 'W/m2', 'Upwelling solar flux at top of atmosphere',       &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNIRTOA'//diag(icall), horiz_only,   'A', 'W/m2',                                                    &
                               'Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', sampling_seq='rad_lwsw')
         call addfld('FSNRTOAC'//diag(icall), horiz_only,   'A', 'W/m2',                                                    &
                      'Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', sampling_seq='rad_lwsw')
         call addfld('FSNRTOAS'//diag(icall), horiz_only,   'A', 'W/m2',                                                    &
                              'Net near-infrared flux (>= 0.7 microns) at top of atmosphere', sampling_seq='rad_lwsw')

         call addfld('FSN200'//diag(icall),   horiz_only,   'A', 'W/m2', 'Net shortwave flux at 200 mb',                    &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSN200C'//diag(icall),  horiz_only,   'A', 'W/m2', 'Clearsky net shortwave flux at 200 mb',           &
                                                                                 sampling_seq='rad_lwsw')

         call addfld('FSNR'//diag(icall),     horiz_only,   'A', 'W/m2', 'Net solar flux at tropopause',                    &
                                                                                 sampling_seq='rad_lwsw')

         call addfld('SOLL'//diag(icall),     horiz_only,   'A', 'W/m2', 'Solar downward near infrared direct  to surface', &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('SOLS'//diag(icall),     horiz_only,   'A', 'W/m2', 'Solar downward visible direct  to surface',       &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('SOLLD'//diag(icall),    horiz_only,   'A', 'W/m2', 'Solar downward near infrared diffuse to surface', &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('SOLSD'//diag(icall),    horiz_only,   'A', 'W/m2', 'Solar downward visible diffuse to surface',       &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNS'//diag(icall),     horiz_only,   'A', 'W/m2', 'Net solar flux at surface',                       &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSNSC'//diag(icall),    horiz_only,   'A', 'W/m2', 'Clearsky net solar flux at surface',              &
                                                                                 sampling_seq='rad_lwsw')

         call addfld('FSDS'//diag(icall),     horiz_only,   'A', 'W/m2', 'Downwelling solar flux at surface',               &
                                                                                 sampling_seq='rad_lwsw')
         call addfld('FSDSC'//diag(icall),    horiz_only,   'A', 'W/m2', 'Clearsky downwelling solar flux at surface',      &
                                                                                 sampling_seq='rad_lwsw')

         call addfld('FUS'//diag(icall),      (/ 'ilev' /), 'I', 'W/m2', 'Shortwave upward flux')
         call addfld('FDS'//diag(icall),      (/ 'ilev' /), 'I', 'W/m2', 'Shortwave downward flux')
         call addfld('FUSC'//diag(icall),     (/ 'ilev' /), 'I', 'W/m2', 'Shortwave clear-sky upward flux')
         call addfld('FDSC'//diag(icall),     (/ 'ilev' /), 'I', 'W/m2', 'Shortwave clear-sky downward flux')

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

         call addfld('QRL'//diag(icall),     (/ 'lev' /), 'A', 'K/s',  'Longwave heating rate', sampling_seq='rad_lwsw')
         call addfld('QRLC'//diag(icall),    (/ 'lev' /), 'A', 'K/s',  'Clearsky longwave heating rate',                   &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNT'//diag(icall),    horiz_only,  'A', 'W/m2', 'Net longwave flux at top of model',                &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNTC'//diag(icall),   horiz_only,  'A', 'W/m2', 'Clearsky net longwave flux at top of model',       &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNTCLR'//diag(icall), horiz_only,  'A', 'W/m2', 'Clearsky ONLY points net longwave flux at top of model',&
                                                                           sampling_seq='rad_lwsw')
         call addfld('FREQCLR'//diag(icall), horiz_only,  'A', 'Frac', 'Frequency of Occurrence of Clearsky',       &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLUT'//diag(icall),    horiz_only,  'A', 'W/m2', 'Upwelling longwave flux at top of model',          &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLUTC'//diag(icall),   horiz_only,  'A', 'W/m2', 'Clearsky upwelling longwave flux at top of model', &
                                                                           sampling_seq='rad_lwsw')
         call addfld('LWCF'//diag(icall),    horiz_only,  'A', 'W/m2', 'Longwave cloud forcing', sampling_seq='rad_lwsw')

         call addfld('FLN200'//diag(icall),  horiz_only,  'A', 'W/m2', 'Net longwave flux at 200 mb',                      &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLN200C'//diag(icall), horiz_only,  'A', 'W/m2', 'Clearsky net longwave flux at 200 mb',             &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNR'//diag(icall),    horiz_only,  'A', 'W/m2', 'Net longwave flux at tropopause',                  &
                                                                           sampling_seq='rad_lwsw')

         call addfld('FLNS'//diag(icall),    horiz_only,  'A', 'W/m2', 'Net longwave flux at surface',                     &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLNSC'//diag(icall),   horiz_only,  'A', 'W/m2', 'Clearsky net longwave flux at surface',            &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLDS'//diag(icall),    horiz_only,  'A', 'W/m2', 'Downwelling longwave flux at surface',             &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FLDSC'//diag(icall),   horiz_only,  'A', 'W/m2', 'Clearsky Downwelling longwave flux at surface',    &
                                                                           sampling_seq='rad_lwsw')
         call addfld('FUL'//diag(icall),     (/ 'ilev' /),'I', 'W/m2', 'Longwave upward flux')
         call addfld('FDL'//diag(icall),     (/ 'ilev' /),'I', 'W/m2', 'Longwave downward flux')
         call addfld('FULC'//diag(icall),    (/ 'ilev' /),'I', 'W/m2', 'Longwave clear-sky upward flux')
         call addfld('FDLC'//diag(icall),    (/ 'ilev' /),'I', 'W/m2', 'Longwave clear-sky downward flux')

         if (history_amwg) then
            call add_default('QRL'//diag(icall),   1, ' ')
            call add_default('FLNT'//diag(icall),  1, ' ')
            call add_default('FLNTC'//diag(icall), 1, ' ')
            call add_default('FLNTCLR'//diag(icall), 1, ' ')
            call add_default('FREQCLR'//diag(icall), 1, ' ')
            call add_default('FLUT'//diag(icall),  1, ' ')
            call add_default('FLUTC'//diag(icall), 1, ' ')
            call add_default('LWCF'//diag(icall),  1, ' ')
            call add_default('FLNS'//diag(icall),  1, ' ')
            call add_default('FLNSC'//diag(icall), 1, ' ')
            call add_default('FLDS'//diag(icall),  1, ' ')
         endif

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

   call pio_seterrorhandling(File, PIO_BCAST_ERROR)

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

   integer :: err_handling
   integer :: ierr

   type(var_desc_t) :: vardesc
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

end subroutine radiation_read_restart
  
!===============================================================================

subroutine radiation_tend( &
   state, ptend, pbuf, cam_out, cam_in, net_flx, rd_out)

   !----------------------------------------------------------------------- 
   ! 
   ! Driver for radiation computation.
   ! 
   ! Revision history:
   ! 2007-11-05  M. Iacono        Install rrtmg_lw and sw as radiation model.
   ! 2007-12-27  M. Iacono        Modify to use CAM cloud optical properties with rrtmg.
   !-----------------------------------------------------------------------
    
   use phys_grid,          only: get_rlat_all_p, get_rlon_all_p
   use cam_control_mod,    only: eccen, mvelpp, lambm0, obliqr
   use shr_orb_mod,        only: shr_orb_decl, shr_orb_cosz

   use aer_rad_props,      only: aer_rad_props_sw, aer_rad_props_lw

   use cloud_rad_props,    only: get_ice_optics_sw, get_liquid_optics_sw, liquid_cloud_get_rad_props_lw, &
                                 ice_cloud_get_rad_props_lw, cloud_rad_props_get_lw, &
                                 snow_cloud_get_rad_props_lw, get_snow_optics_sw
   use slingo,             only: slingo_liq_get_rad_props_lw, slingo_liq_optics_sw
   use ebert_curry,        only: ec_ice_optics_sw, ec_ice_get_rad_props_lw

   use rad_solar_var,      only: get_variability
   use radsw,              only: rad_rrtmg_sw
   use radlw,              only: rad_rrtmg_lw
   use radheat,            only: radheat_tend

   use radiation_data,     only: rad_data_write
   use rrtmg_state,        only: rrtmg_state_create, rrtmg_state_update, rrtmg_state_destroy, rrtmg_state_t, &
                                 num_rrtmg_levs

   use interpolate_data,   only: vertinterp
   use tropopause,         only: tropopause_find, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE

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
  
   integer  :: i, k
   integer  :: lchnk, ncol
   logical  :: dosw, dolw

   real(r8) :: calday          ! current calendar day
   real(r8) :: delta           ! Solar declination angle  in radians
   real(r8) :: eccf            ! Earth orbit eccentricity factor
   real(r8) :: clat(pcols)     ! current latitudes(radians)
   real(r8) :: clon(pcols)     ! current longitudes(radians)
   real(r8) :: coszrs(pcols)   ! Cosine solar zenith angle

   ! Gathered indices of day and night columns 
   !  chunk_column_index = IdxDay(daylight_column_index)
   integer :: Nday           ! Number of daylight columns
   integer :: Nnite          ! Number of night columns
   integer :: IdxDay(pcols)  ! Indices of daylight columns
   integer :: IdxNite(pcols) ! Indices of night columns

   integer :: itim_old

   real(r8), pointer :: cld(:,:)      ! cloud fraction
   real(r8), pointer :: cldfsnow(:,:) ! cloud fraction of just "snow clouds- whatever they are"
   real(r8), pointer :: qrs(:,:)      ! shortwave radiative heating rate 
   real(r8), pointer :: qrl(:,:)      ! longwave  radiative heating rate 
   real(r8), pointer :: fsds(:)  ! Surface solar down flux
   real(r8), pointer :: fsns(:)  ! Surface solar absorbed flux
   real(r8), pointer :: fsnt(:)  ! Net column abs solar flux at model top
   real(r8), pointer :: flns(:)  ! Srf longwave cooling (up-down) flux
   real(r8), pointer :: flnt(:)  ! Net outgoing lw flux at model top

   real(r8), pointer, dimension(:,:,:) :: su => NULL()  ! shortwave spectral flux up
   real(r8), pointer, dimension(:,:,:) :: sd => NULL()  ! shortwave spectral flux down
   real(r8), pointer, dimension(:,:,:) :: lu => NULL()  ! longwave  spectral flux up
   real(r8), pointer, dimension(:,:,:) :: ld => NULL()  ! longwave  spectral flux down

   ! tropopause diagnostic
   integer  :: troplev(pcols)
   real(r8) :: p_trop(pcols)

   type(rrtmg_state_t) :: r_state ! contains the atm concentrations in layers needed for RRTMG

   ! cloud radiative parameters are "in cloud" not "in cell"
   real(r8) :: ice_tau    (nswbands,pcols,pver) ! ice extinction optical depth
   real(r8) :: ice_tau_w  (nswbands,pcols,pver) ! ice single scattering albedo * tau
   real(r8) :: ice_tau_w_g(nswbands,pcols,pver) ! ice assymetry parameter * tau * w
   real(r8) :: ice_tau_w_f(nswbands,pcols,pver) ! ice forward scattered fraction * tau * w
   real(r8) :: ice_lw_abs (nlwbands,pcols,pver)   ! ice absorption optics depth (LW)

   ! cloud radiative parameters are "in cloud" not "in cell"
   real(r8) :: liq_tau    (nswbands,pcols,pver) ! liquid extinction optical depth
   real(r8) :: liq_tau_w  (nswbands,pcols,pver) ! liquid single scattering albedo * tau
   real(r8) :: liq_tau_w_g(nswbands,pcols,pver) ! liquid assymetry parameter * tau * w
   real(r8) :: liq_tau_w_f(nswbands,pcols,pver) ! liquid forward scattered fraction * tau * w
   real(r8) :: liq_lw_abs (nlwbands,pcols,pver) ! liquid absorption optics depth (LW)

   ! cloud radiative parameters are "in cloud" not "in cell"
   real(r8) :: cld_tau    (nswbands,pcols,pver) ! cloud extinction optical depth
   real(r8) :: cld_tau_w  (nswbands,pcols,pver) ! cloud single scattering albedo * tau
   real(r8) :: cld_tau_w_g(nswbands,pcols,pver) ! cloud assymetry parameter * w * tau
   real(r8) :: cld_tau_w_f(nswbands,pcols,pver) ! cloud forward scattered fraction * w * tau
   real(r8) :: cld_lw_abs (nlwbands,pcols,pver) ! cloud absorption optics depth (LW)

   ! cloud radiative parameters are "in cloud" not "in cell"
   real(r8) :: snow_tau    (nswbands,pcols,pver) ! snow extinction optical depth
   real(r8) :: snow_tau_w  (nswbands,pcols,pver) ! snow single scattering albedo * tau
   real(r8) :: snow_tau_w_g(nswbands,pcols,pver) ! snow assymetry parameter * tau * w
   real(r8) :: snow_tau_w_f(nswbands,pcols,pver) ! snow forward scattered fraction * tau * w
   real(r8) :: snow_lw_abs (nlwbands,pcols,pver)! snow absorption optics depth (LW)

   ! combined cloud radiative parameters are "in cloud" not "in cell"
   real(r8) :: cldfprime(pcols,pver)              ! combined cloud fraction (snow plus regular)
   real(r8) :: c_cld_tau    (nswbands,pcols,pver) ! combined cloud extinction optical depth
   real(r8) :: c_cld_tau_w  (nswbands,pcols,pver) ! combined cloud single scattering albedo * tau
   real(r8) :: c_cld_tau_w_g(nswbands,pcols,pver) ! combined cloud assymetry parameter * w * tau
   real(r8) :: c_cld_tau_w_f(nswbands,pcols,pver) ! combined cloud forward scattered fraction * w * tau
   real(r8) :: c_cld_lw_abs (nlwbands,pcols,pver) ! combined cloud absorption optics depth (LW)

   real(r8) :: sfac(1:nswbands)     ! time varying scaling factors due to Solar Spectral Irrad at 1 A.U. per band

   integer :: icall                 ! index through climate/diagnostic radiation calls
   logical :: active_calls(0:N_DIAG)

   ! Aerosol radiative properties
   real(r8) :: aer_tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
   real(r8) :: aer_tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
   real(r8) :: aer_tau_w_g(pcols,0:pver,nswbands) ! aerosol assymetry parameter * w * tau
   real(r8) :: aer_tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * w * tau
   real(r8) :: aer_lw_abs (pcols,pver,nlwbands)   ! aerosol absorption optics depth (LW)

   real(r8) :: fns(pcols,pverp)     ! net shortwave flux
   real(r8) :: fcns(pcols,pverp)    ! net clear-sky shortwave flux
   real(r8) :: fnl(pcols,pverp)     ! net longwave flux
   real(r8) :: fcnl(pcols,pverp)    ! net clear-sky longwave flux
  
   ! for COSP
   real(r8) :: emis(pcols,pver)        ! Cloud longwave emissivity
   real(r8) :: gb_snow_tau(pcols,pver) ! grid-box mean snow_tau
   real(r8) :: gb_snow_lw(pcols,pver)  ! grid-box mean LW snow optical depth

   real(r8) :: ftem(pcols,pver)        ! Temporary workspace for outfld variables

   real(r8) :: freqclr(pcols)          ! Frequency of occurrence of clear sky columns
   real(r8) :: flntclr(pcols)          ! Clearsky only columns (zero if cloudy)

   character(*), parameter :: name = 'radiation_tend'
   !--------------------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol

   if (present(rd_out)) then
      rd => rd_out
      write_output = .false.
   else
      allocate(rd)
      write_output=.true.
   end if

   dosw = radiation_do('sw')      ! do shortwave heating calc this timestep?
   dolw = radiation_do('lw')      ! do longwave heating calc this timestep?

   ! Cosine solar zenith angle for current time step
   calday = get_curr_calday()
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)

   call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, &
                     delta, eccf)
   if (use_rad_uniform_angle) then
      do i = 1, ncol
         coszrs(i) = shr_orb_cosz(calday, clat(i), clon(i), delta, dt_avg, uniform_angle=rad_uniform_angle)
      end do
   else
      do i = 1, ncol
         coszrs(i) = shr_orb_cosz(calday, clat(i), clon(i), delta, dt_avg)
      end do
   end if

   ! Gather night/day column indices.
   Nday = 0
   Nnite = 0
   do i = 1, ncol
      if ( coszrs(i) > 0.0_r8 ) then
         Nday = Nday + 1
         IdxDay(Nday) = i
      else
         Nnite = Nnite + 1
         IdxNite(Nnite) = i
      end if
   end do

   ! Associate pointers to physics buffer fields
   itim_old = pbuf_old_tim_idx()
   if (cldfsnow_idx > 0) then
      call pbuf_get_field(pbuf, cldfsnow_idx, cldfsnow, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   endif
   call pbuf_get_field(pbuf, cld_idx, cld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, qrs_idx, qrs)
   call pbuf_get_field(pbuf, qrl_idx, qrl)

   call pbuf_get_field(pbuf, fsnt_idx, fsnt)
   call pbuf_get_field(pbuf, fsds_idx, fsds)
   call pbuf_get_field(pbuf, fsns_idx, fsns)
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
      call tropopause_find(state, troplev, tropP=p_trop, primary=TROP_ALG_HYBSTOB, backup=TROP_ALG_CLIMATE)
   endif

   if (dosw .or. dolw) then

      ! construct an RRTMG state object
      r_state = rrtmg_state_create( state, cam_in )

      call t_startf('cldoptics')

      if (cldfsnow_idx > 0) then
         do k = 1, pver
            do i = 1, ncol
               cldfprime(i,k) = max(cld(i,k), cldfsnow(i,k))
            end do
         end do
      else
         cldfprime(:ncol,:) = cld(:ncol,:)
      end if
      
      
      if (dosw) then

         if (oldcldoptics) then
            call ec_ice_optics_sw(state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp=.false.)
            call slingo_liq_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp=.false.)
         else
            select case (icecldoptics)
            case ('ebertcurry')
               call  ec_ice_optics_sw(state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp=.true.)
            case ('mitchell')
               call get_ice_optics_sw(state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f)
            case default
               call endrun('iccldoptics must be one either ebertcurry or mitchell')
            end select

            select case (liqcldoptics)
            case ('slingo')
               call slingo_liq_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp=.true.)
            case ('gammadist')
               call get_liquid_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f)
            case default
               call endrun('liqcldoptics must be either slingo or gammadist')
            end select
         end if

         cld_tau(:,:ncol,:)     =  liq_tau(:,:ncol,:)     + ice_tau(:,:ncol,:)
         cld_tau_w(:,:ncol,:)   =  liq_tau_w(:,:ncol,:)   + ice_tau_w(:,:ncol,:)
         cld_tau_w_g(:,:ncol,:) =  liq_tau_w_g(:,:ncol,:) + ice_tau_w_g(:,:ncol,:)
         cld_tau_w_f(:,:ncol,:) =  liq_tau_w_f(:,:ncol,:) + ice_tau_w_f(:,:ncol,:)
 
         if (cldfsnow_idx > 0) then
            ! add in snow
            call get_snow_optics_sw(state, pbuf, snow_tau, snow_tau_w, snow_tau_w_g, snow_tau_w_f)
            do i = 1, ncol
               do k = 1, pver

                  if (cldfprime(i,k) > 0._r8) then

                     c_cld_tau(:,i,k)     = ( cldfsnow(i,k)*snow_tau(:,i,k) &
                                             + cld(i,k)*cld_tau(:,i,k) )/cldfprime(i,k)

                     c_cld_tau_w(:,i,k)   = ( cldfsnow(i,k)*snow_tau_w(:,i,k)  &
                                             + cld(i,k)*cld_tau_w(:,i,k) )/cldfprime(i,k)

                     c_cld_tau_w_g(:,i,k) = ( cldfsnow(i,k)*snow_tau_w_g(:,i,k) &
                                             + cld(i,k)*cld_tau_w_g(:,i,k) )/cldfprime(i,k)

                     c_cld_tau_w_f(:,i,k) = ( cldfsnow(i,k)*snow_tau_w_f(:,i,k) &
                                             + cld(i,k)*cld_tau_w_f(:,i,k) )/cldfprime(i,k)
                  else
                     c_cld_tau(:,i,k)     = 0._r8
                     c_cld_tau_w(:,i,k)   = 0._r8
                     c_cld_tau_w_g(:,i,k) = 0._r8
                     c_cld_tau_w_f(:,i,k) = 0._r8
                  end if
               end do
            end do
         else
            c_cld_tau(:,:ncol,:)     = cld_tau(:,:ncol,:)
            c_cld_tau_w(:,:ncol,:)   = cld_tau_w(:,:ncol,:)
            c_cld_tau_w_g(:,:ncol,:) = cld_tau_w_g(:,:ncol,:)
            c_cld_tau_w_f(:,:ncol,:) = cld_tau_w_f(:,:ncol,:)
         end if

         ! Output cloud optical depth fields for the visible band
         rd%tot_icld_vistau(:ncol,:) = c_cld_tau(idx_sw_diag,:ncol,:)
         rd%liq_icld_vistau(:ncol,:) = liq_tau(idx_sw_diag,:ncol,:)
         rd%ice_icld_vistau(:ncol,:) = ice_tau(idx_sw_diag,:ncol,:)

         if (cldfsnow_idx > 0) then
            rd%snow_icld_vistau(:ncol,:) = snow_tau(idx_sw_diag,:ncol,:)
         endif

         ! multiply by total cloud fraction to get gridbox value
         rd%tot_cld_vistau(:ncol,:) = c_cld_tau(idx_sw_diag,:ncol,:)*cldfprime(:ncol,:)

         ! add fillvalue for night columns
         do i = 1, Nnite
            rd%tot_cld_vistau(IdxNite(i),:)   = fillvalue
            rd%tot_icld_vistau(IdxNite(i),:)  = fillvalue
            rd%liq_icld_vistau(IdxNite(i),:)  = fillvalue
            rd%ice_icld_vistau(IdxNite(i),:)  = fillvalue
            if (cldfsnow_idx > 0) then
               rd%snow_icld_vistau(IdxNite(i),:) = fillvalue
            end if
         end do

         if (write_output) call radiation_output_cld(lchnk, ncol, rd)

      end if   ! if (dosw)

      if (dolw) then

         if (oldcldoptics) then
            call cloud_rad_props_get_lw(state, pbuf, cld_lw_abs, oldcloud=.true.)
         else
            select case (icecldoptics)
            case ('ebertcurry')
               call ec_ice_get_rad_props_lw(state, pbuf, ice_lw_abs, oldicewp=.true.)
            case ('mitchell')
               call ice_cloud_get_rad_props_lw(state, pbuf, ice_lw_abs)
            case default
               call endrun('iccldoptics must be one either ebertcurry or mitchell')
            end select

            select case (liqcldoptics)
            case ('slingo')
               call slingo_liq_get_rad_props_lw(state, pbuf, liq_lw_abs, oldliqwp=.true.)
            case ('gammadist')
               call liquid_cloud_get_rad_props_lw(state, pbuf, liq_lw_abs)
            case default
               call endrun('liqcldoptics must be either slingo or gammadist')
            end select

            cld_lw_abs(:,:ncol,:) = liq_lw_abs(:,:ncol,:) + ice_lw_abs(:,:ncol,:)

         end if

         if (cldfsnow_idx > 0) then

            ! add in snow
            call snow_cloud_get_rad_props_lw(state, pbuf, snow_lw_abs)

            do i = 1, ncol
               do k = 1, pver
                  if (cldfprime(i,k) > 0._r8) then
                     c_cld_lw_abs(:,i,k) = ( cldfsnow(i,k)*snow_lw_abs(:,i,k) &
                                            + cld(i,k)*cld_lw_abs(:,i,k) )/cldfprime(i,k)
                  else
                     c_cld_lw_abs(:,i,k) = 0._r8
                  end if
               end do
            end do
         else
            c_cld_lw_abs(:,:ncol,:) = cld_lw_abs(:,:ncol,:)
         end if

      end if   ! if (dolw)

      call t_stopf('cldoptics')

      ! Solar radiation computation

      if (dosw) then

         call get_variability(sfac)

         ! Get the active climate/diagnostic shortwave calculations
         call rad_cnst_get_call_list(active_calls)

         ! The climate (icall==0) calculation must occur last.
         do icall = N_DIAG, 0, -1

            if (active_calls(icall)) then

               ! update the concentrations in the RRTMG state object
               call rrtmg_state_update(state, pbuf, icall, r_state)

               call aer_rad_props_sw(icall, state, pbuf, nnite, idxnite, &
                                     aer_tau, aer_tau_w, aer_tau_w_g, aer_tau_w_f)
               
               rd%cld_tau_cloudsim(:ncol,:) = cld_tau(rrtmg_sw_cloudsim_band,:ncol,:)
               rd%aer_tau550(:ncol,:)       = aer_tau(:ncol,:,idx_sw_diag)
               rd%aer_tau400(:ncol,:)       = aer_tau(:ncol,:,idx_sw_diag+1)
               rd%aer_tau700(:ncol,:)       = aer_tau(:ncol,:,idx_sw_diag-1)

               call rad_rrtmg_sw( &
                  lchnk, ncol, num_rrtmg_levs, r_state, state%pmid,          &
                  cldfprime, aer_tau, aer_tau_w, aer_tau_w_g,  aer_tau_w_f,  &
                  eccf, coszrs, rd%solin, sfac, cam_in%asdir,                &
                  cam_in%asdif, cam_in%aldir, cam_in%aldif, qrs, rd%qrsc,    &
                  fsnt, rd%fsntc, rd%fsntoa, rd%fsutoa, rd%fsntoac,          &
                  rd%fsnirt, rd%fsnrtc, rd%fsnirtsq, fsns, rd%fsnsc,         &
                  rd%fsdsc, fsds, cam_out%sols, cam_out%soll, cam_out%solsd, &
                  cam_out%solld, fns, fcns, Nday, Nnite,                     &
                  IdxDay, IdxNite, su, sd, E_cld_tau=c_cld_tau,              &
                  E_cld_tau_w=c_cld_tau_w, E_cld_tau_w_g=c_cld_tau_w_g,      &
                  E_cld_tau_w_f=c_cld_tau_w_f, old_convert=.false.)

               ! Output net fluxes at 200 mb
               call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcns, rd%fsn200c)
               call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fns,  rd%fsn200)
               if (hist_fld_active('FSNR')) then
                  do i = 1,ncol
                     call vertinterp(1, 1, pverp, state%pint(i,:), p_trop(i), fns(i,:), rd%fsnr(i))
                  end do
               end if

               if (write_output) call radiation_output_sw(lchnk, ncol, icall, rd, pbuf, cam_out)

            end if
         end do

      end if

      ! Output aerosol mmr
      call rad_cnst_out(0, state, pbuf)
                 
      ! Longwave radiation computation

      if (dolw) then

         call rad_cnst_get_call_list(active_calls)

         ! The climate (icall==0) calculation must occur last.
         do icall = N_DIAG, 0, -1

            if (active_calls(icall)) then

               ! update the conctrations in the RRTMG state object
               call  rrtmg_state_update( state, pbuf, icall, r_state)

               call aer_rad_props_lw(icall, state, pbuf,  aer_lw_abs)
                  
               call rad_rrtmg_lw( &
                  lchnk, ncol, num_rrtmg_levs, r_state, state%pmid,  &
                  aer_lw_abs, cldfprime, c_cld_lw_abs, qrl, rd%qrlc, &
                  flns, flnt, rd%flnsc, rd%flntc, cam_out%flwds,     &
                  rd%flut, rd%flutc, fnl, fcnl, rd%fldsc,            &
                  lu, ld)

               !  Output fluxes at 200 mb
               call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fnl,  rd%fln200)
               call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcnl, rd%fln200c)
               if (hist_fld_active('FLNR')) then
                  do i = 1,ncol
                     call vertinterp(1, 1, pverp, state%pint(i,:), p_trop(i), fnl(i,:), rd%flnr(i))
                  end do
               end if

               flntclr(:) = 0._r8
               freqclr(:) = 0._r8
               do i = 1, ncol
                  if (maxval(cldfprime(i,:)) <= 0.1_r8) then
                     freqclr(i) = 1._r8
                     flntclr(i) = rd%flntc(i)
                  end if
               end do

               if (write_output) call radiation_output_lw(lchnk, ncol, icall, rd, pbuf, cam_out, freqclr, flntclr)

            end if
         end do

      end if

      ! deconstruct the RRTMG state object
      call rrtmg_state_destroy(r_state)

      if (docosp) then

         ! initialize and calculate emis
         emis(:,:) = 0._r8
         emis(:ncol,:) = 1._r8 - exp(-cld_lw_abs(rrtmg_lw_cloudsim_band,:ncol,:))
         call outfld('EMIS', emis, pcols, lchnk)

         ! compute grid-box mean SW and LW snow optical depth for use by COSP
         gb_snow_tau(:,:) = 0._r8
         gb_snow_lw(:,:)  = 0._r8
         if (cldfsnow_idx > 0) then
            do i = 1, ncol
               do k = 1, pver
                  if (cldfsnow(i,k) > 0._r8) then
                     gb_snow_tau(i,k) = snow_tau(rrtmg_sw_cloudsim_band,i,k)*cldfsnow(i,k)
                     gb_snow_lw(i,k)  = snow_lw_abs(rrtmg_lw_cloudsim_band,i,k)*cldfsnow(i,k)
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
            call cospsimulator_intr_run(state,  pbuf, cam_in, emis, coszrs, &
               cld_swtau_in=cld_tau(rrtmg_sw_cloudsim_band,:,:),&
               snow_tau_in=gb_snow_tau, snow_emis_in=gb_snow_lw)
            cosp_cnt(lchnk) = 0
         end if
      end if

   else   !  if (dosw .or. dolw) then

      ! convert radiative heating rates from Q*dp to Q for energy conservation
      do k =1 , pver
         do i = 1, ncol
            qrs(i,k) = qrs(i,k)/state%pdel(i,k)
            qrl(i,k) = qrl(i,k)/state%pdel(i,k)
         end do
      end do

   end if   ! if (dosw .or. dolw) then

    ! output rad inputs and resulting heating rates
    call rad_data_write(  pbuf, state, cam_in, coszrs )

   ! Compute net radiative heating tendency
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

   ! convert radiative heating rates to Q*dp for energy conservation
   do k = 1, pver
      do i = 1, ncol
         qrs(i,k) = qrs(i,k)*state%pdel(i,k)
         qrl(i,k) = qrl(i,k)*state%pdel(i,k)
      end do
   end do

   cam_out%netsw(:ncol) = fsns(:ncol)

   if (.not. present(rd_out)) then
      deallocate(rd)
   end if

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

   call outfld('QRS'//diag(icall),      qrs(:ncol,:)/cpair,     ncol, lchnk)
   call outfld('QRSC'//diag(icall),     rd%qrsc(:ncol,:)/cpair, ncol, lchnk)

   call outfld('FSNT'//diag(icall),     fsnt,          pcols, lchnk)
   call outfld('FSNTC'//diag(icall),    rd%fsntc,      pcols, lchnk)
   call outfld('FSNTOA'//diag(icall),   rd%fsntoa,     pcols, lchnk)
   call outfld('FSNTOAC'//diag(icall),  rd%fsntoac,    pcols, lchnk)

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

end subroutine radiation_output_sw


!===============================================================================

subroutine radiation_output_cld(lchnk, ncol, rd)

   ! Dump shortwave cloud optics information to history buffer.

   integer ,        intent(in) :: lchnk
   integer,         intent(in) :: ncol
   type(rad_out_t), intent(in) :: rd
   !----------------------------------------------------------------------------

   call outfld('TOT_CLD_VISTAU',  rd%tot_cld_vistau,  pcols, lchnk)
   call outfld('TOT_ICLD_VISTAU', rd%tot_icld_vistau, pcols, lchnk)
   call outfld('LIQ_ICLD_VISTAU', rd%liq_icld_vistau, pcols, lchnk)
   call outfld('ICE_ICLD_VISTAU', rd%ice_icld_vistau, pcols, lchnk)
   if (cldfsnow_idx > 0) then
      call outfld('SNOW_ICLD_VISTAU', rd%snow_icld_vistau, pcols, lchnk)
   endif

end subroutine radiation_output_cld

!===============================================================================

subroutine radiation_output_lw(lchnk, ncol, icall, rd, pbuf, cam_out, freqclr, flntclr)

   ! Dump longwave radiation information to history buffer

   integer,                intent(in) :: lchnk
   integer,                intent(in) :: ncol
   integer,                intent(in) :: icall  ! icall=0 for climate diagnostics
   type(rad_out_t),        intent(in) :: rd
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cam_out_t),        intent(in) :: cam_out
   real(r8),               intent(in) :: freqclr(pcols)
   real(r8),               intent(in) :: flntclr(pcols)

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

   call outfld('FREQCLR'//diag(icall), freqclr,       pcols, lchnk)
   call outfld('FLNTCLR'//diag(icall), flntclr,       pcols, lchnk)

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

end subroutine radiation_output_lw

!===============================================================================

subroutine calc_col_mean(state, mmr_pointer, mean_value)

   ! Compute the column mean mass mixing ratio.  

   type(physics_state),        intent(in)  :: state
   real(r8), dimension(:,:),   pointer     :: mmr_pointer  ! mass mixing ratio (lev)
   real(r8), dimension(pcols), intent(out) :: mean_value   ! column mean mmr

   integer  :: i, k, ncol
   real(r8) :: ptot(pcols)
   !-----------------------------------------------------------------------

   ncol         = state%ncol
   mean_value   = 0.0_r8
   ptot         = 0.0_r8

   do k=1,pver
      do i=1,ncol
         mean_value(i) = mean_value(i) + mmr_pointer(i,k)*state%pdeldry(i,k)
         ptot(i)         = ptot(i) + state%pdeldry(i,k)
      end do
   end do
   do i=1,ncol
      mean_value(i) = mean_value(i) / ptot(i)
   end do

end subroutine calc_col_mean

!===============================================================================

end module radiation

