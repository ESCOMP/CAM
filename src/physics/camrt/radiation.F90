module radiation

!---------------------------------------------------------------------------------
!
! CAM interface to the legacy 'camrt' radiation code
!
!---------------------------------------------------------------------------------

use shr_kind_mod,        only: r8=>shr_kind_r8, cl=>shr_kind_cl
use spmd_utils,          only: masterproc
use ppgrid,              only: pcols, pver, pverp, begchunk, endchunk
use physics_types,       only: physics_state, physics_ptend
use phys_grid,           only: get_ncols_p
use camsrfexch,          only: cam_out_t, cam_in_t    
use physconst,           only: cpair, cappa
use time_manager,        only: get_nstep, is_first_restart_step, &
                               get_curr_calday, get_step_size
use cam_control_mod,     only: lambm0, obliqr, mvelpp, eccen

use radae,               only: abstot_3d, absnxt_3d, emstot_3d, initialize_radbuffer, ntoplw

use scamMod,             only: scm_crm_mode, single_column,have_cld,cldobs,&
                               have_clwp,clwpobs,have_tg,tground

use cam_grid_support,    only: cam_grid_write_attr, cam_grid_id,            &
                               cam_grid_header_info_t, cam_grid_dimensions, &
                               cam_grid_write_dist_array, cam_grid_read_dist_array

use cam_history,         only: outfld, hist_fld_active
use cam_history_support, only: fillvalue

use cam_pio_utils,       only: cam_pio_def_dim

use pio,                 only: file_desc_t, var_desc_t, &
                               pio_double, pio_int, pio_noerr, &
                               pio_seterrorhandling, pio_bcast_error, &
                               pio_inq_varid, &
                               pio_def_var, pio_def_dim, &
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
   radiation_define_restart, &!
   radiation_write_restart,  &!
   radiation_read_restart,   &!
   radiation_tend,           &! compute heating rates and fluxes
   rad_out_t                  ! type for diagnostic outputs

type rad_out_t
   real(r8) :: solin(pcols)         ! Solar incident flux
   real(r8) :: fsntoa(pcols)        ! Net solar flux at TOA
   real(r8) :: fsutoa(pcols)        ! upwelling solar flux at TOA
   real(r8) :: fsntoac(pcols)       ! Clear sky net solar flux at TOA
   real(r8) :: fsnirt(pcols)        ! Near-IR flux absorbed at toa
   real(r8) :: fsnrtc(pcols)        ! Clear sky near-IR flux absorbed at toa
   real(r8) :: fsnirtsq(pcols)      ! Near-IR flux absorbed at toa >= 0.7 microns
   real(r8) :: fsntc(pcols)         ! Clear sky total column abs solar flux
   real(r8) :: fsnsc(pcols)         ! Clear sky surface abs solar flux
   real(r8) :: fsdsc(pcols)         ! Clear sky surface downwelling solar flux
   real(r8) :: flut(pcols)          ! Upward flux at top of model
   real(r8) :: flutc(pcols)         ! Upward Clear Sky flux at top of model
   real(r8) :: flntc(pcols)         ! Clear sky lw flux at model top
   real(r8) :: flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
   real(r8) :: fldsc(pcols)         ! Clear sky lw flux at srf (down)
   real(r8) :: flwds(pcols)         ! Down longwave flux at surface
   real(r8) :: fsnr(pcols)
   real(r8) :: flnr(pcols)
   real(r8) :: fsds(pcols)          ! Surface solar down flux
   real(r8) :: fln200(pcols)        ! net longwave flux interpolated to 200 mb
   real(r8) :: fln200c(pcols)       ! net clearsky longwave flux interpolated to 200 mb
   real(r8) :: fsn200(pcols)        ! fns interpolated to 200 mb
   real(r8) :: fsn200c(pcols)       ! fcns interpolated to 200 mb
   real(r8) :: sols(pcols)          ! Solar downward visible direct  to surface
   real(r8) :: soll(pcols)          ! Solar downward near infrared direct  to surface
   real(r8) :: solsd(pcols)         ! Solar downward visible diffuse to surface
   real(r8) :: solld(pcols)         ! Solar downward near infrared diffuse to surface
   real(r8) :: qrsc(pcols,pver)     ! clearsky shortwave radiative heating rate
   real(r8) :: qrlc(pcols,pver)     ! clearsky longwave  radiative heating rate
   real(r8) :: fsdtoa(pcols)        ! Solar input = Flux Solar Downward Top of Atmosphere
   real(r8) :: swcf(pcols)          ! shortwave cloud forcing
   real(r8) :: lwcf(pcols)          ! longwave cloud forcing

   real(r8) :: tot_cld_vistau(pcols,pver)  
   real(r8) :: tot_icld_vistau(pcols,pver) 
   real(r8) :: liq_icld_vistau(pcols,pver)  ! in-cld liq cloud optical depth (only during day, night = fillvalue)
   real(r8) :: ice_icld_vistau(pcols,pver)  ! in-cld ice cloud optical depth (only during day, night = fillvalue)
end type rad_out_t

! Namelist variables

character(len=cl) :: absems_data
integer :: iradsw = -1     ! freq. of shortwave radiation calc in time steps (positive)
                           ! or hours (negative).
integer :: iradlw = -1     ! frequency of longwave rad. calc. in time steps (positive)
                           ! or hours (negative).
integer :: iradae = -12    ! frequency of absorp/emis calc in time steps (positive)
                           ! or hours (negative).
integer :: irad_always = 0 ! Specifies length of time in timesteps (positive)
                           ! or hours (negative) SW/LW radiation will be
                           ! run continuously from the start of an
                           ! initial or restart run
logical :: use_rad_dt_cosz = .false. ! if true use zenith angle averaged over 
                                     ! interval between radiation calculations

! Physics buffer indices
integer :: qrs_idx      = 0
integer :: qrl_idx      = 0 
integer :: fsds_idx     = 0
integer :: fsns_idx     = 0
integer :: fsnt_idx     = 0
integer :: flns_idx     = 0
integer :: flnt_idx     = 0
integer :: cld_idx      = 0
integer :: rel_idx      = 0
integer :: rei_idx      = 0
integer :: cicewp_idx   = -1
integer :: cliqwp_idx   = -1
integer :: cldemis_idx  = -1
integer :: cldtau_idx   = -1
integer :: nmxrgn_idx   = -1
integer :: pmxrgn_idx   = -1

! averaging time interval for zenith angle
real(r8) :: dt_avg = 0._r8

real(r8), parameter :: cgs2mks = 1.e-3_r8

! PIO descriptors (for restarts)

type(var_desc_t), allocatable :: abstot_desc(:)
type(var_desc_t) :: emstot_desc, absnxt_desc(4)

logical  :: use_rad_uniform_angle = .false. ! if true, use the namelist rad_uniform_angle for the zenith calculation
real(r8) :: rad_uniform_angle = -99._r8

!===============================================================================
contains
!===============================================================================

subroutine radiation_readnl(nlfile)

   ! Read radiation_nl namelist group.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_logical, &
                              mpi_character, mpi_real8

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   integer :: dtime      ! timestep size
   character(len=*), parameter :: sub = 'radiation_readnl'

   namelist /radiation_nl/ absems_data, iradsw, iradlw, iradae, irad_always, &
                           use_rad_dt_cosz, use_rad_uniform_angle, rad_uniform_angle
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
   call mpi_bcast(absems_data, len(absems_data), mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: absems_data")
   call mpi_bcast(iradsw, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: iradsw")
   call mpi_bcast(iradlw, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: iradlw")
   call mpi_bcast(iradae, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: iradae")
   call mpi_bcast(irad_always, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: irad_always")
   call mpi_bcast(use_rad_dt_cosz, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: use_rad_dt_cosz")
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

   ! Convert iradae from hours to timesteps if necessary and check that
   ! iradae must be an even multiple of iradlw
   if (iradae < 0) iradae = nint((-iradae*3600._r8)/dtime)
   if (mod(iradae,iradlw)/=0) then
      write(iulog,*) sub//': iradae must be an even multiple of iradlw.'
      write(iulog,*)'     iradae = ',iradae,', iradlw = ',iradlw
      call endrun(sub//': iradae must be an even multiple of iradlw.')
   end if

   !----------------------------------------------------------------------- 
   ! Print runtime options to log.
   !-----------------------------------------------------------------------

   if (masterproc) then
      write(iulog,*) 'CAMRT radiation scheme parameters:'
      write(iulog,10) iradsw, iradlw, iradae, irad_always, use_rad_dt_cosz
      write(iulog,*) '  Abs/Emis dataset: ', trim(absems_data)
   end if

10 format('  Frequency (timesteps) of Shortwave Radiation calc:    ',i5/, &
          '  Frequency (timesteps) of Longwave Radiation calc:     ',i5/, &
          '  Frequency (timesteps) of Absorptivity/Emissivity calc:',i5/, &
          '  SW/LW calc done every timestep for first N steps. N=  ',i5/, &
          '  Use average zenith angle:                             ',l5)


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

   call rad_data_register()

end subroutine radiation_register

!================================================================================================

function radiation_do(op, timestep)

   ! Returns true if the specified operation is done this timestep.

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

   case ('absems') ! do an absorptivity/emissivity calculation this timestep?
      radiation_do = nstep == 0  .or.  iradae == 1                     &
                    .or. (mod(nstep-1,iradae) == 0  .and.  nstep /= 1)

   case ('aeres') ! write absorptivity/emissivity to restart file this timestep?
      radiation_do = mod(nstep,iradae) /= 0
         
   case default
      call endrun('radiation_do: unknown operation:'//op)

   end select
end function radiation_do

!================================================================================================

real(r8) function radiation_nextsw_cday()
  
   ! Returns calendar day of next sw radiation calculation

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

   use cam_history,    only: addfld, add_default, horiz_only
   use physconst,      only: gravit, cpair, epsilo, stebol, &
                             pstd, mwdry, mwco2, mwo3
    
   use physics_buffer, only: physics_buffer_desc, pbuf_get_index
   use radsw,          only: radsw_init
   use radlw,          only: radlw_init
   use radae,          only: radae_init
   use radconstants,   only: radconstants_init
   use rad_solar_var,  only: rad_solar_var_init
   use radiation_data, only: rad_data_init
   use phys_control,   only: phys_getopts
   use time_manager,   only: get_step_size

   ! args
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)


   ! Local variables
   integer :: nstep                       ! current timestep number
   logical :: history_amwg                ! output the variables used by the AMWG diag package
   logical :: history_vdiag               ! output the variables used by the AMWG variability diag package
   logical :: history_budget              ! output tendencies and state variables for CAM4
                                          ! temperature, water vapor, cloud ice and cloud
                                          ! liquid budgets.
   integer :: history_budget_histfile_num ! output history file number for budget fields

   integer :: dtime
   !-----------------------------------------------------------------------

   call radconstants_init()
   call rad_solar_var_init()

   call radsw_init(gravit)
   call radlw_init(gravit, stebol)
   call radae_init( &
      gravit, epsilo, stebol, pstd, mwdry, &
      mwco2, mwo3, absems_data)

   call rad_data_init(pbuf2d)

   ! Set the radiation timestep for cosz calculations if requested using the adjusted iradsw value from radiation
   if (use_rad_dt_cosz)  then
      dtime  = get_step_size()
      dt_avg = real(iradsw*dtime, r8)
   end if

   ! Get physics buffer indices
   cld_idx    = pbuf_get_index('CLD')
   rel_idx    = pbuf_get_index('REL')
   rei_idx    = pbuf_get_index('REI')

   ! "irad_always" is number of time steps to execute radiation continuously from start of
   ! initial OR restart run
   nstep = get_nstep()
   if ( irad_always > 0) then
      nstep       = get_nstep()
      irad_always = irad_always + nstep
   end if

   ! Shortwave radiation
   call addfld ('SOLIN',           horiz_only,   'A','W/m2','Solar insolation',                            sampling_seq='rad_lwsw')
   call addfld ('SOLL',            horiz_only,   'A','W/m2','Solar downward near infrared direct  to surface',              &
                                                                                                           sampling_seq='rad_lwsw')
   call addfld ('SOLS',            horiz_only,   'A','W/m2','Solar downward visible direct  to surface',   sampling_seq='rad_lwsw')
   call addfld ('SOLLD',           horiz_only,   'A','W/m2','Solar downward near infrared diffuse to surface',              &
                                                                                                           sampling_seq='rad_lwsw')
   call addfld ('SOLSD',           horiz_only,   'A','W/m2','Solar downward visible diffuse to surface',   sampling_seq='rad_lwsw')
   call addfld ('QRS',             (/ 'lev' /),  'A','K/s', 'Solar heating rate',                          sampling_seq='rad_lwsw')
   call addfld ('QRSC',            (/ 'lev' /),  'A','K/s', 'Clearsky solar heating rate',                 sampling_seq='rad_lwsw')
   call addfld ('FSNS',            horiz_only,   'A','W/m2','Net solar flux at surface',                   sampling_seq='rad_lwsw')
   call addfld ('FSNT',            horiz_only,   'A','W/m2','Net solar flux at top of model',              sampling_seq='rad_lwsw')
   call addfld ('FSNTOA',          horiz_only,   'A','W/m2','Net solar flux at top of atmosphere',         sampling_seq='rad_lwsw')
   call addfld ('FSUTOA',          horiz_only,   'A','W/m2','Upwelling solar flux at top of atmosphere',   sampling_seq='rad_lwsw')
   call addfld ('FSNTOAC',         horiz_only,   'A','W/m2','Clearsky net solar flux at top of atmosphere',                 &
                                                                                                           sampling_seq='rad_lwsw')
   call addfld ('FSDTOA',          horiz_only,   'A','W/m2','Downwelling solar flux at top of atmosphere', sampling_seq='rad_lwsw')
   call addfld ('FSN200',          horiz_only,   'A','W/m2','Net shortwave flux at 200 mb',                sampling_seq='rad_lwsw')
   call addfld ('FSN200C',         horiz_only,   'A','W/m2','Clearsky net shortwave flux at 200 mb',       sampling_seq='rad_lwsw')
   call addfld ('FSNTC',           horiz_only,   'A','W/m2','Clearsky net solar flux at top of model',     sampling_seq='rad_lwsw')
   call addfld ('FSNSC',           horiz_only,   'A','W/m2','Clearsky net solar flux at surface',          sampling_seq='rad_lwsw')
   call addfld ('FSDSC',           horiz_only,   'A','W/m2','Clearsky downwelling solar flux at surface',                   &
                                                                                                           sampling_seq='rad_lwsw')
   call addfld ('FSDS',            horiz_only,   'A','W/m2','Downwelling solar flux at surface',           sampling_seq='rad_lwsw')
   call addfld ('FUS',             (/ 'ilev' /), 'I','W/m2','Shortwave upward flux')
   call addfld ('FDS',             (/ 'ilev' /), 'I','W/m2','Shortwave downward flux')
   call addfld ('FUSC',            (/ 'ilev' /), 'I','W/m2','Shortwave clear-sky upward flux')
   call addfld ('FDSC',            (/ 'ilev' /), 'I','W/m2','Shortwave clear-sky downward flux')
   call addfld ('FSNIRTOA',        horiz_only,   'A','W/m2','Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere',  &
                                                                                                           sampling_seq='rad_lwsw')
   call addfld ('FSNRTOAC',        horiz_only,   'A','W/m2',                                                                &
         'Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere',                           sampling_seq='rad_lwsw')
   call addfld ('FSNRTOAS',        horiz_only,   'A','W/m2','Net near-infrared flux (>= 0.7 microns) at top of atmosphere', &
                                                                                                           sampling_seq='rad_lwsw')
   call addfld ('FSNR',            horiz_only,   'A','W/m2','Net solar flux at tropopause',                sampling_seq='rad_lwsw')
   call addfld ('SWCF',            horiz_only,   'A','W/m2','Shortwave cloud forcing',                     sampling_seq='rad_lwsw')

   call addfld ('TOT_CLD_VISTAU',  (/ 'lev' /),  'A','1',   'Total gbx cloud visible sw optical depth',                     &
                                                                                         sampling_seq='rad_lwsw',flag_xyfill=.true.)
   call addfld ('TOT_ICLD_VISTAU', (/ 'lev' /),  'A','1',   'Total in-cloud visible sw optical depth',                      &
                                                                                         sampling_seq='rad_lwsw',flag_xyfill=.true.)
   call addfld ('LIQ_ICLD_VISTAU', (/ 'lev' /),  'A','1',   'Liquid in-cloud visible sw optical depth',                     &
                                                                                         sampling_seq='rad_lwsw',flag_xyfill=.true.)
   call addfld ('ICE_ICLD_VISTAU', (/ 'lev' /),  'A','1',   'Ice in-cloud visible sw optical depth',                        &
                                                                                         sampling_seq='rad_lwsw',flag_xyfill=.true.)

   ! Longwave radiation
   call addfld ('QRL',             (/ 'lev' /),  'A','K/s', 'Longwave heating rate',                       sampling_seq='rad_lwsw')
   call addfld ('QRLC',            (/ 'lev' /),  'A','K/s', 'Clearsky longwave heating rate',              sampling_seq='rad_lwsw')
   call addfld ('FLNS',            horiz_only,   'A','W/m2','Net longwave flux at surface',                sampling_seq='rad_lwsw')
   call addfld ('FLDS',            horiz_only,   'A','W/m2','Downwelling longwave flux at surface',        sampling_seq='rad_lwsw')
   call addfld ('FLNT',            horiz_only,   'A','W/m2','Net longwave flux at top of model',           sampling_seq='rad_lwsw')
   call addfld ('FLUT',            horiz_only,   'A','W/m2','Upwelling longwave flux at top of model',     sampling_seq='rad_lwsw')
   call addfld ('FLUTC',           horiz_only,   'A','W/m2','Clearsky upwelling longwave flux at top of model',             &
                                                                                                           sampling_seq='rad_lwsw')
   call addfld ('FLNTC',           horiz_only,   'A','W/m2','Clearsky net longwave flux at top of model',  sampling_seq='rad_lwsw')
   call addfld ('FLN200',          horiz_only,   'A','W/m2','Net longwave flux at 200 mb',                 sampling_seq='rad_lwsw')
   call addfld ('FLN200C',         horiz_only,   'A','W/m2','Clearsky net longwave flux at 200 mb',        sampling_seq='rad_lwsw')
   call addfld ('FLNR',            horiz_only,   'A','W/m2','Net longwave flux at tropopause',             sampling_seq='rad_lwsw')
   call addfld ('FLNSC',           horiz_only,   'A','W/m2','Clearsky net longwave flux at surface',       sampling_seq='rad_lwsw')
   call addfld ('FLDSC',           horiz_only,   'A','W/m2','Clearsky downwelling longwave flux at surface',                &
                                                                                                           sampling_seq='rad_lwsw')
   call addfld ('LWCF',            horiz_only,   'A','W/m2','Longwave cloud forcing',                      sampling_seq='rad_lwsw')
   call addfld ('FUL',             (/ 'ilev' /), 'I','W/m2','Longwave upward flux')
   call addfld ('FDL',             (/ 'ilev' /), 'I','W/m2','Longwave downward flux')
   call addfld ('FULC',            (/ 'ilev' /), 'I','W/m2','Longwave clear-sky upward flux')
   call addfld ('FDLC',            (/ 'ilev' /), 'I','W/m2','Longwave clear-sky downward flux')

   ! Heating rate needed for d(theta)/dt computation
   call addfld ('HR',              (/ 'lev' /),  'A','K/s', 'Heating rate needed for d(theta)/dt computation')
   
   ! determine default variables
   call phys_getopts(history_amwg_out   = history_amwg,   &
                      history_vdiag_out  = history_vdiag,  &
                      history_budget_out = history_budget, &
                      history_budget_histfile_num_out = history_budget_histfile_num)

   if (history_amwg) then
      ! Shortwave variables
      call add_default ('SOLIN   ', 1, ' ')
      call add_default ('QRS     ', 1, ' ')
      call add_default ('FSNS    ', 1, ' ')
      call add_default ('FSNT    ', 1, ' ')
      call add_default ('FSDTOA  ', 1, ' ')
      call add_default ('FSNTOA  ', 1, ' ')
      call add_default ('FSUTOA  ', 1, ' ')
      call add_default ('FSNTOAC ', 1, ' ')
      call add_default ('FSNTC   ', 1, ' ')
      call add_default ('FSNSC   ', 1, ' ')
      call add_default ('FSDSC   ', 1, ' ')
      call add_default ('FSDS    ', 1, ' ')
      call add_default ('SWCF    ', 1, ' ')
      ! Longwave variables
      call add_default ('QRL     ', 1, ' ')
      call add_default ('FLNS    ', 1, ' ')
      call add_default ('FLDS    ', 1, ' ')
      call add_default ('FLNT    ', 1, ' ')
      call add_default ('FLUT    ', 1, ' ')
      call add_default ('FLUTC   ', 1, ' ')
      call add_default ('FLNTC   ', 1, ' ')
      call add_default ('FLNSC   ', 1, ' ')
      call add_default ('FLDSC   ', 1, ' ')
      call add_default ('LWCF    ', 1, ' ')
   endif
   if (single_column.and.scm_crm_mode) then
      ! Shortwave variables
      call add_default ('FUS     ', 1, ' ')
      call add_default ('FUSC    ', 1, ' ')
      call add_default ('FDS     ', 1, ' ')
      call add_default ('FDSC    ', 1, ' ')
      ! Longwave variables
      call add_default ('FUL     ', 1, ' ')
      call add_default ('FULC    ', 1, ' ')
      call add_default ('FDL     ', 1, ' ')
      call add_default ('FDLC    ', 1, ' ')
   endif
    
   if ( history_budget .and. history_budget_histfile_num > 1 ) then
      call add_default ('QRL     ', history_budget_histfile_num, ' ')
      call add_default ('QRS     ', history_budget_histfile_num, ' ')
   end if
 
   if (history_vdiag) then
      call add_default('FLUT',2,' ')
      call add_default('FLUT',3,' ')
   end if
   
   cicewp_idx = pbuf_get_index('CICEWP')
   cliqwp_idx = pbuf_get_index('CLIQWP')
   cldemis_idx= pbuf_get_index('CLDEMIS')
   cldtau_idx = pbuf_get_index('CLDTAU')
   nmxrgn_idx = pbuf_get_index('NMXRGN')
   pmxrgn_idx = pbuf_get_index('PMXRGN')

end subroutine radiation_init

!===============================================================================

subroutine radiation_define_restart(file)

   ! define variables to be written to restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables
   integer :: i, ierr
   integer :: grid_id
   integer :: hdimcnt
   integer :: pver_id, pverp_id
   integer :: vsize
   integer :: dimids(4)

   type(cam_grid_header_info_t) :: info

   character(len=16) :: pname
   !----------------------------------------------------------------------------

   call pio_seterrorhandling(File, PIO_BCAST_ERROR)

   if (radiation_do('aeres')) then

      grid_id = cam_grid_id('physgrid')
      call cam_grid_write_attr(File, grid_id, info)
      hdimcnt = info%num_hdims()
      do i = 1, hdimcnt
         dimids(i) = info%get_hdimid(i)
      end do

      call cam_pio_def_dim(File, 'lev',  pver,  pver_id,  existOK=.true.)
      call cam_pio_def_dim(File, 'ilev', pverp, pverp_id, existOK=.true.)

      vsize = pverp - ntoplw + 1
      if (vsize /= pverp) then
         ierr = pio_def_dim(File, 'lwcols', vsize, dimids(hdimcnt+1))
      else
         dimids(hdimcnt+1) = pverp_id
      end if

      ! split into vsize variables to avoid excessive memory usage in IO

      allocate(abstot_desc(ntoplw:pverp))

      do i = ntoplw, pverp
         write(pname,'(a,i3.3)') 'NAL_absorp', i
         ierr = pio_def_var(File, trim(pname), pio_double, dimids(1:hdimcnt+1), abstot_desc(i))
      end do
	
      dimids(hdimcnt+1) = pverp_id
      ierr = pio_def_var(File, 'Emissivity', pio_double, dimids(1:hdimcnt+1), emstot_desc)

      dimids(hdimcnt+1) = pver_id
      do i=1,4
         write(pname,'(a,i3.3)') 'NN_absorp',i
         ierr = pio_def_var(File, pname, pio_double, dimids(1:hdimcnt+1), absnxt_desc(i))
      end do

   end if

end subroutine radiation_define_restart
  
!===============================================================================

subroutine radiation_write_restart(file)

   ! write variables to restart file

   ! arguments
   type(file_desc_t), intent(inout) :: file

   ! local variables
   integer :: i, ierr
   integer :: physgrid
   integer :: dims(3), gdims(3)
   integer :: nhdims
   integer :: ncol
   !----------------------------------------------------------------------------

   if ( radiation_do('aeres')  ) then
         
      physgrid = cam_grid_id('physgrid')
      call cam_grid_dimensions(physgrid, gdims(1:2), nhdims)

      do i = begchunk, endchunk
         ncol = get_ncols_p(i)
         if (ncol < pcols) then
            abstot_3d(ncol+1:pcols,:,:,i) = fillvalue
            absnxt_3d(ncol+1:pcols,:,:,i) = fillvalue
            emstot_3d(ncol+1:pcols,:,i)   = fillvalue
         end if
      end do
      
      ! abstot_3d is written as a series of 3D variables

      dims(1) = size(abstot_3d, 1) ! Should be pcols
      dims(2) = size(abstot_3d, 2) ! Should be (pverp-ntoplw+1)
      dims(3) = size(abstot_3d, 4) ! Should be endchunk - begchunk + 1
      gdims(nhdims+1) = dims(2)
      do i = ntoplw, pverp
         call cam_grid_write_dist_array(File, physgrid, dims(1:3),             &
              gdims(1:nhdims+1), abstot_3d(:,:,i,:), abstot_desc(i))
      end do

      dims(1) = size(emstot_3d, 1) ! Should be pcols
      dims(2) = size(emstot_3d, 2) ! Should be pverp
      dims(3) = size(emstot_3d, 3) ! Should be endchunk - begchunk + 1
      gdims(nhdims+1) = dims(2)
      call cam_grid_write_dist_array(File, physgrid, dims(1:3),               &
           gdims(1:nhdims+1), emstot_3d, emstot_desc)

      dims(1) = size(absnxt_3d, 1) ! Should be pcols
      dims(2) = size(absnxt_3d, 2) ! Should be pver
      dims(3) = size(absnxt_3d, 4) ! Should be endchunk - begchunk + 1
      gdims(nhdims+1) = dims(2)
      do i = 1, 4
         call cam_grid_write_dist_array(File, physgrid, dims(1:3),             &
              gdims(1:nhdims+1), absnxt_3d(:,:,i,:), absnxt_desc(i))
      end do

      ! module data was allocated in radiation_define_restart
      deallocate(abstot_desc)
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
   integer :: physgrid
   integer :: dims(3), gdims(3), nhdims
   integer :: vsize
   integer :: i

   type(var_desc_t) :: vardesc
   character(len=16) :: pname
   !----------------------------------------------------------------------------

   ! Put this call here for now.  It should move to an init method when the
   ! initialization and restart sequencing is unified.
   call initialize_radbuffer()

   if ( radiation_do('aeres')  ) then

      call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)
      ierr = pio_inq_varid(File, 'Emissivity', vardesc)
      call pio_seterrorhandling(File, err_handling)
      if (ierr /= PIO_NOERR) then
         if (masterproc) write(iulog,*) 'Warning: Emissivity variable not found on restart file.'
         return
      end if

      physgrid = cam_grid_id('physgrid')
      call cam_grid_dimensions(physgrid, gdims(1:2), nhdims)

      dims(1) = pcols
      dims(2) = pverp
      dims(3) = endchunk - begchunk + 1
      gdims(nhdims+1) = dims(2)

      call cam_grid_read_dist_array(File, physgrid, dims(1:3), &
           gdims(1:nhdims+1), emstot_3d, vardesc)
        
      vsize = pverp - ntoplw + 1
      dims(2) = vsize
      gdims(nhdims+1) = dims(2)
        
      do i = ntoplw, pverp
         write(pname,'(a,i3.3)') 'NAL_absorp', i
         ierr = pio_inq_varid(File, trim(pname), vardesc)
         call cam_grid_read_dist_array(File, physgrid, dims(1:3), &
              gdims(1:nhdims+1), abstot_3d(:,:,i,:), vardesc)
      end do

      dims(2) = pver
      gdims(nhdims+1) = dims(2)
      do i = 1, 4
         write(pname,'(a,i3.3)') 'NN_absorp', i
         ierr = pio_inq_varid(File, trim(pname), vardesc)
         call cam_grid_read_dist_array(File, physgrid, dims(1:3), &
              gdims(1:nhdims+1), absnxt_3d(:,:,i,:), vardesc)
      end do
   end if

end subroutine radiation_read_restart
  
!===============================================================================

subroutine radiation_tend( &
   state, ptend, pbuf, cam_out, cam_in, net_flx, rd_out)

   !-----------------------------------------------------------------------
   ! Driver for radiation computation.
   !
   ! NOTE: Radiation uses cgs units, so conversions must be done from
   !       model fields to radiation fields.
   !-----------------------------------------------------------------------

   use physics_buffer,      only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
   use phys_grid,           only: get_rlat_all_p, get_rlon_all_p
   use physics_types,       only: physics_state, physics_ptend
   use time_manager,        only: get_curr_calday
   use radheat,             only: radheat_tend
   use physconst,           only: cpair, stebol
   use radconstants,        only: nlwbands, nswbands
   use radsw,               only: radcswmx
   use radlw,               only: radclwmx
   use rad_constituents,    only: rad_cnst_get_gas, rad_cnst_out
   use aer_rad_props,       only: aer_rad_props_sw, aer_rad_props_lw
   use interpolate_data,    only: vertinterp
   use radiation_data,      only: rad_data_write
   use cloud_cover_diags,   only: cloud_cover_diags_out
   use tropopause,          only: tropopause_find, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE
   use orbit,               only: zenith

   ! Arguments
   type(physics_state),       target, intent(in)    :: state
   type(physics_ptend),               intent(out)   :: ptend
   type(physics_buffer_desc),         pointer       :: pbuf(:)
   type(cam_out_t),                   intent(inout) :: cam_out
   type(cam_in_t),                    intent(in)    :: cam_in
   real(r8),                          intent(out)   :: net_flx(pcols)
   type(rad_out_t), target, optional, intent(out)   :: rd_out

   ! Local variables
   type(rad_out_t), pointer :: rd  ! allow rd_out to be optional by allocating a local object
                                   ! if the argument is not present

   integer :: i, k
   integer :: lchnk, ncol

   logical :: dosw, dolw, doabsems
   integer, pointer :: nmxrgn(:)              ! pbuf pointer to Number of maximally overlapped regions
   real(r8),pointer :: pmxrgn(:,:)            ! Maximum values of pressure for each
                                              !    maximally overlapped region.
                                              !    0->pmxrgn(i,1) is range of pressure for
                                              !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
                                              !    2nd region, etc

   real(r8),pointer :: emis(:,:)              ! Cloud longwave emissivity
   real(r8),pointer :: cldtau(:,:)            ! Cloud longwave optical depth
   real(r8),pointer :: cicewp(:,:)            ! in-cloud cloud ice water path
   real(r8),pointer :: cliqwp(:,:)            ! in-cloud cloud liquid water path

   real(r8) :: cltot(pcols)                   ! Diagnostic total cloud cover
   real(r8) :: cllow(pcols)                   !       "     low  cloud cover
   real(r8) :: clmed(pcols)                   !       "     mid  cloud cover
   real(r8) :: clhgh(pcols)                   !       "     hgh  cloud cover

   real(r8) :: ftem(pcols,pver)               ! Temporary workspace for outfld variables

   integer :: itim_old
   real(r8), pointer, dimension(:,:) :: rel     ! liquid effective drop radius (microns)
   real(r8), pointer, dimension(:,:) :: rei     ! ice effective drop size (microns)
   real(r8), pointer, dimension(:,:) :: cld     ! cloud fraction
   real(r8), pointer, dimension(:,:) :: qrs     ! shortwave radiative heating rate
   real(r8), pointer, dimension(:,:) :: qrl     ! longwave  radiative heating rate

   real(r8) :: calday                        ! current calendar day
   real(r8) :: clat(pcols)                   ! current latitudes(radians)
   real(r8) :: clon(pcols)                   ! current longitudes(radians)
   real(r8) :: coszrs(pcols)                 ! Cosine solar zenith angle

   real(r8) :: fns(pcols,pverp)     ! net shortwave flux
   real(r8) :: fcns(pcols,pverp)    ! net clear-sky shortwave flux
   real(r8) :: fnl(pcols,pverp)     ! net longwave flux
   real(r8) :: fcnl(pcols,pverp)    ! net clear-sky longwave flux

   ! This is used by the chemistry.
   real(r8), pointer :: fsds(:)  ! Surface solar down flux

   ! This is used for the energy checker and the Eulerian dycore.
   real(r8), pointer :: fsns(:)  ! Surface solar absorbed flux
   real(r8), pointer :: fsnt(:)  ! Net column abs solar flux at model top
   real(r8), pointer :: flns(:)  ! Srf longwave cooling (up-down) flux
   real(r8), pointer :: flnt(:)  ! Net outgoing lw flux at model top

   real(r8) :: pbr(pcols,pver)      ! Model mid-level pressures (dynes/cm2)
   real(r8) :: pnm(pcols,pverp)     ! Model interface pressures (dynes/cm2)
   real(r8) :: eccf                 ! Earth/sun distance factor
   real(r8) :: lwupcgs(pcols)       ! Upward longwave flux in cgs units
 
   real(r8), pointer, dimension(:,:) :: n2o    ! nitrous oxide mass mixing ratio
   real(r8), pointer, dimension(:,:) :: ch4    ! methane mass mixing ratio
   real(r8), pointer, dimension(:,:) :: cfc11  ! cfc11 mass mixing ratio
   real(r8), pointer, dimension(:,:) :: cfc12  ! cfc12 mass mixing ratio
   real(r8), pointer, dimension(:,:) :: o3     ! Ozone mass mixing ratio
   real(r8), pointer, dimension(:,:) :: o2     ! Oxygen mass mixing ratio
   real(r8), dimension(pcols) :: o2_col        ! column oxygen mmr
   real(r8), pointer, dimension(:,:) :: co2    ! co2   mass mixing ratio
   real(r8), dimension(pcols) :: co2_col_mean  ! co2 column mean mmr
   real(r8), pointer, dimension(:,:) :: sp_hum ! specific humidity

   ! Aerosol shortwave radiative properties
   real(r8) :: aer_tau    (pcols,0:pver,nswbands) ! aerosol extinction optical depth
   real(r8) :: aer_tau_w  (pcols,0:pver,nswbands) ! aerosol single scattering albedo * tau
   real(r8) :: aer_tau_w_g(pcols,0:pver,nswbands) ! aerosol assymetry parameter * w * tau
   real(r8) :: aer_tau_w_f(pcols,0:pver,nswbands) ! aerosol forward scattered fraction * w * tau

   ! Aerosol longwave absorption optical depth
   real(r8) :: odap_aer(pcols,pver,nlwbands)

   ! Gathered indicies of day and night columns
   !  chunk_column_index = IdxDay(daylight_column_index)
   integer :: Nday                      ! Number of daylight columns
   integer :: Nnite                     ! Number of night columns
   integer, dimension(pcols) :: IdxDay  ! Indicies of daylight coumns
   integer, dimension(pcols) :: IdxNite ! Indicies of night coumns

   character(*), parameter :: name = 'radiation_tend'

   ! tropopause diagnostic
   integer :: troplev(pcols)
   real(r8):: p_trop(pcols)

   logical :: write_output ! switch for outfld calls
   !----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol

   calday = get_curr_calday()

   if (present(rd_out)) then
      rd => rd_out
      write_output = .false.
   else
      allocate(rd)
      write_output=.true.
   end if

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,    cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, qrs_idx,qrs)
   call pbuf_get_field(pbuf, qrl_idx,qrl)

   call pbuf_get_field(pbuf, fsds_idx, fsds)

   call pbuf_get_field(pbuf, fsns_idx, fsns)
   call pbuf_get_field(pbuf, fsnt_idx, fsnt)
   call pbuf_get_field(pbuf, flns_idx, flns)
   call pbuf_get_field(pbuf, flnt_idx, flnt)

   call pbuf_get_field(pbuf, rel_idx, rel)
   call pbuf_get_field(pbuf, rei_idx, rei)
   
   !  For CRM, make cloud equal to input observations:
   if (single_column.and.scm_crm_mode.and.have_cld) then
      do k = 1,pver
         cld(:ncol,k)= cldobs(k)
      enddo
   endif

   ! Cosine solar zenith angle for current time step
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)
   if (use_rad_uniform_angle) then
     call zenith (calday, clat, clon, coszrs, ncol, dt_avg, uniform_angle=rad_uniform_angle)
   else
     call zenith (calday, clat, clon, coszrs, ncol, dt_avg)
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

   dosw     = radiation_do('sw')      ! do shortwave heating calc this timestep?
   dolw     = radiation_do('lw')      ! do longwave heating calc this timestep?

   doabsems = radiation_do('absems')  ! do absorptivity/emissivity calc this timestep?

   if (dosw .or. dolw) then

      ! pbuf cloud properties set in cloud_diagnostics
      call pbuf_get_field(pbuf, cicewp_idx, cicewp)
      call pbuf_get_field(pbuf, cliqwp_idx, cliqwp)
      call pbuf_get_field(pbuf, cldemis_idx, emis)

      call pbuf_get_field(pbuf, cldtau_idx, cldtau)

      call pbuf_get_field(pbuf, pmxrgn_idx, pmxrgn)
      call pbuf_get_field(pbuf, nmxrgn_idx, nmxrgn)

      ! For CRM, make cloud liquid water path equal to input observations
      if(single_column.and.scm_crm_mode.and.have_clwp)then
         do k=1,pver
            cliqwp(:ncol,k) = clwpobs(k)
         end do
      endif

      ! Get specific humidity
      call rad_cnst_get_gas(0,'H2O', state, pbuf,  sp_hum)

      ! Get ozone mass mixing ratio.
      call rad_cnst_get_gas(0,'O3',  state, pbuf,  o3)

      ! Get CO2 mass mixing ratio and compute column mean values
      call rad_cnst_get_gas(0,'CO2', state, pbuf,  co2)
      call calc_col_mean(state, co2, co2_col_mean)

      ! construct cgs unit reps of pmid and pint and get "eccf" - earthsundistancefactor
      call radinp(ncol, state%pmid, state%pint, pbr, pnm, eccf)

      ! Solar radiation computation

      if (hist_fld_active('FSNR') .or. hist_fld_active('FLNR')) then
         call tropopause_find(state, troplev, tropP=p_trop, primary=TROP_ALG_HYBSTOB, backup=TROP_ALG_CLIMATE)
      endif

      if (dosw) then

         call t_startf('rad_sw')

         ! Get Oxygen mass mixing ratio.
         call rad_cnst_get_gas(0,'O2', state, pbuf,  o2)
         call calc_col_mean(state, o2, o2_col)
   
         ! Get aerosol radiative properties.
         call t_startf('aero_optics_sw')
         call aer_rad_props_sw(0, state, pbuf,  nnite, idxnite, &
            aer_tau, aer_tau_w, aer_tau_w_g, aer_tau_w_f)
         call t_stopf('aero_optics_sw')

         call radcswmx(lchnk, &
               ncol,       pnm,        pbr,        sp_hum,     o3,         &
               o2_col,     cld,        cicewp,     cliqwp,     rel,        &
               rei,        eccf,       coszrs,     rd%solin,      &
               cam_in%asdir, cam_in%asdif, cam_in%aldir, cam_in%aldif, nmxrgn, &
               pmxrgn,     qrs,        rd%qrsc,       fsnt,       rd%fsntc,      rd%fsdtoa, &
               rd%fsntoa,     rd%fsutoa,     rd%fsntoac,    rd%fsnirt,     rd%fsnrtc,     rd%fsnirtsq,   &
               fsns,       rd%fsnsc,      rd%fsdsc,      fsds,       cam_out%sols, &
               cam_out%soll, cam_out%solsd, cam_out%solld, fns, fcns,      &
               Nday,       Nnite,      IdxDay,     IdxNite,    co2_col_mean, &
               aer_tau,    aer_tau_w,  aer_tau_w_g, aer_tau_w_f , rd%liq_icld_vistau, rd%ice_icld_vistau  ) 

         call t_stopf('rad_sw')

         !  Output net fluxes at 200 mb
         call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcns, rd%fsn200c)
         call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fns,  rd%fsn200)
         if (hist_fld_active('FSNR')) then
            do i = 1,ncol
               call vertinterp(1, 1, pverp, state%pint(i,:), p_trop(i), fns(i,:), rd%fsnr(i))
               rd%fsnr(i) = rd%fsnr(i)*cgs2mks
            enddo
         else
            rd%fsnr(:) = 0._r8
         endif

         ! Convert units of shortwave fields needed by rest of model from CGS to MKS

         do i=1,ncol
            rd%solin(i)   = rd%solin(i)   *cgs2mks
            fsds(i)       = fsds(i)       *cgs2mks
            rd%fsnirt(i)  = rd%fsnirt(i)  *cgs2mks
            rd%fsnrtc(i)  = rd%fsnrtc(i)  *cgs2mks
            rd%fsnirtsq(i)= rd%fsnirtsq(i)*cgs2mks
            fsnt(i)       = fsnt(i)       *cgs2mks
            rd%fsdtoa(i)  = rd%fsdtoa(i)  *cgs2mks
            fsns(i)       = fsns(i)       *cgs2mks
            rd%fsntc(i)   = rd%fsntc(i)   *cgs2mks
            rd%fsnsc(i)   = rd%fsnsc(i)   *cgs2mks
            rd%fsdsc(i)   = rd%fsdsc(i)   *cgs2mks
            rd%fsntoa(i)  = rd%fsntoa(i)  *cgs2mks
            rd%fsutoa(i)  = rd%fsutoa(i)  *cgs2mks
            rd%fsntoac(i) = rd%fsntoac(i) *cgs2mks
            rd%fsn200(i)  = rd%fsn200(i)  *cgs2mks
            rd%fsn200c(i) = rd%fsn200c(i) *cgs2mks
            rd%swcf(i)    = rd%fsntoa(i) - rd%fsntoac(i)
         end do

         ! initialize tau_cld_vistau and tau_icld_vistau as fillvalue, they will stay fillvalue for night columns
         rd%tot_icld_vistau(1:pcols,1:pver) = fillvalue
         rd%tot_cld_vistau(1:pcols,1:pver)  = fillvalue

         ! only do calcs for tot_cld_vistau and tot_icld_vistau on daytime columns
         do i=1,Nday
            ! sum the water and ice optical depths to get total in-cloud cloud optical depth
            rd%tot_icld_vistau(IdxDay(i),1:pver) = rd%liq_icld_vistau(IdxDay(i),1:pver) + &
                                                   rd%ice_icld_vistau(IdxDay(i),1:pver)

            ! sum wat and ice, multiply by cloud fraction to get grid-box value
            rd%tot_cld_vistau(IdxDay(i),1:pver) = (rd%liq_icld_vistau(IdxDay(i),1:pver) + &
                                                   rd%ice_icld_vistau(IdxDay(i),1:pver))*cld(IdxDay(i),1:pver)
         end do

         ! add fillvalue for night columns
         do i = 1, Nnite
            rd%liq_icld_vistau(IdxNite(i),:)  = fillvalue
            rd%ice_icld_vistau(IdxNite(i),:)  = fillvalue
         end do

         if (write_output) call radiation_output_sw(state, rd, cam_out, fsns, fsnt, fsds, qrs)

      end if   ! dosw

      ! Longwave radiation computation

      if (dolw) then

         call t_startf("rad_lw")

         ! Convert upward longwave flux units to CGS

         do i=1,ncol
            lwupcgs(i) = cam_in%lwup(i)*1000._r8
            if (single_column .and. scm_crm_mode .and. have_tg) &
               lwupcgs(i) = 1000*stebol*tground(1)**4
         end do

         ! Get gas phase constituents.
         call rad_cnst_get_gas(0,'N2O',   state, pbuf,  n2o)
         call rad_cnst_get_gas(0,'CH4',   state, pbuf,  ch4)
         call rad_cnst_get_gas(0,'CFC11', state, pbuf,  cfc11)
         call rad_cnst_get_gas(0,'CFC12', state, pbuf,  cfc12)

         ! absems requires lw absorption optical depth and transmission through aerosols
         call t_startf('aero_optics_lw')
         if (doabsems) call aer_rad_props_lw(0, state, pbuf, odap_aer)
         call t_stopf('aero_optics_lw')

         call radclwmx(lchnk, ncol, doabsems, &
                 lwupcgs, state%t, sp_hum, o3, pbr, &
                 pnm, state%lnpmid, state%lnpint, n2o, ch4, &
                 cfc11, cfc12, cld, emis, pmxrgn, &
                 nmxrgn, qrl, rd%qrlc, flns, flnt, rd%flnsc, &
                 rd%flntc, cam_out%flwds, rd%fldsc, rd%flut, rd%flutc, &
                 fnl, fcnl, co2_col_mean, odap_aer)

         call t_stopf("rad_lw")

         !  Output fluxes at 200 mb
         call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fnl,  rd%fln200)
         call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcnl, rd%fln200c)
         if (hist_fld_active('FLNR')) then
            do i = 1,ncol
               call vertinterp(1, 1, pverp, state%pint(i,:), p_trop(i), fnl(i,:), rd%flnr(i))
            enddo
         else
            rd%flnr(:) = 0._r8
         endif

         ! Convert units of longwave fields needed by rest of model from CGS to MKS

         do i = 1, ncol
            flnt(i)          = flnt(i)            *cgs2mks
            rd%flut(i)       = rd%flut(i)         *cgs2mks
            rd%flutc(i)      = rd%flutc(i)        *cgs2mks
            rd%lwcf(i)       = rd%flutc(i) - rd%flut(i)
            flns(i)          = flns(i)            *cgs2mks
            rd%fldsc(i)      = rd%fldsc(i)        *cgs2mks
            rd%flntc(i)      = rd%flntc(i)        *cgs2mks
            rd%fln200(i)     = rd%fln200(i)       *cgs2mks
            rd%fln200c(i)    = rd%fln200c(i)      *cgs2mks
            rd%flnsc(i)      = rd%flnsc(i)        *cgs2mks
            cam_out%flwds(i) = cam_out%flwds(i)   *cgs2mks
            rd%flnr(i)       = rd%flnr(i)         *cgs2mks
         end do

         if (write_output) call radiation_output_lw(state,  rd, cam_out, flns, flnt, qrl)

      end if  ! dolw

      ! Output aerosol mmr
      if (write_output) call rad_cnst_out(0, state, pbuf)

      ! Cloud cover diagnostics
      ! radsw can change pmxrgn and nmxrgn so cldsav needs to follow radsw
      if (write_output) call cloud_cover_diags_out(lchnk, ncol, cld, state%pmid, nmxrgn, pmxrgn )

   else   !  if (dosw .or. dolw) then

      ! convert radiative heating rates from Q*dp to Q for energy conservation
      do k =1 , pver
         do i = 1, ncol
            qrs(i,k) = qrs(i,k)/state%pdel(i,k)
            qrl(i,k) = qrl(i,k)/state%pdel(i,k)
         end do
      end do

   end if   !  if (dosw .or. dolw) then

   ! output rad inputs and resulting heating rates
   call rad_data_write( pbuf, state, cam_in, coszrs )

   ! Compute net radiative heating tendency
   call radheat_tend(state, pbuf,  ptend, qrl, qrs, fsns, &
                     fsnt, flns, flnt, cam_in%asdir, net_flx)

   if (write_output) then
      ! Compute heating rate for dtheta/dt
      do k=1,pver
         do i=1,ncol
            ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5_r8/state%pmid(i,k))**cappa
         end do
      end do
      call outfld('HR      ',ftem    ,pcols   ,lchnk   )
   end if
       
   ! convert radiative heating rates to Q*dp for energy conservation
   do k =1 , pver
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
  
subroutine radiation_output_sw(state, rd, cam_out, fsns, fsnt, fsds, qrs)

   ! Dump shortwave radiation information to history buffer (diagnostics)

   type(physics_state), intent(in) :: state
   type(rad_out_t),     intent(in) :: rd
   type(cam_out_t),     intent(in) :: cam_out
   real(r8),            intent(in) :: fsns(pcols) ! Surface solar absorbed flux
   real(r8),            intent(in) :: fsnt(pcols) ! Net column abs solar flux at model top
   real(r8),            intent(in) :: fsds(pcols) ! Surface solar down flux
   real(r8),            pointer    :: qrs(:,:)    ! shortwave radiative heating rate

   ! Local variables
   integer :: lchnk, ncol
   real(r8) :: ftem(pcols,pver)               ! Temporary workspace for outfld variables
   !----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol

   ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair
   call outfld('QRS     ',ftem  ,pcols,lchnk)
   
   ftem(:ncol,:pver) = rd%qrsc(:ncol,:pver)/cpair
   call outfld('QRSC    ',ftem  ,pcols,lchnk)

   call outfld('SOLIN   ',rd%solin      ,pcols,lchnk)
   call outfld('FSDS    ',fsds       ,pcols,lchnk)
   call outfld('FSNIRTOA',rd%fsnirt     ,pcols,lchnk)
   call outfld('FSNRTOAC',rd%fsnrtc     ,pcols,lchnk)
   call outfld('FSNRTOAS',rd%fsnirtsq   ,pcols,lchnk)
   call outfld('FSNT    ',fsnt       ,pcols,lchnk)
   call outfld('FSDTOA  ',rd%fsdtoa     ,pcols,lchnk)
   call outfld('FSNS    ',fsns       ,pcols,lchnk)
   call outfld('FSNTC   ',rd%fsntc      ,pcols,lchnk)
   call outfld('FSNSC   ',rd%fsnsc      ,pcols,lchnk)
   call outfld('FSDSC   ',rd%fsdsc      ,pcols,lchnk)
   call outfld('FSNTOA  ',rd%fsntoa     ,pcols,lchnk)
   call outfld('FSUTOA  ',rd%fsutoa     ,pcols,lchnk)
   call outfld('FSNTOAC ',rd%fsntoac    ,pcols,lchnk)
   call outfld('SOLS    ',cam_out%sols  ,pcols,lchnk)
   call outfld('SOLL    ',cam_out%soll  ,pcols,lchnk)
   call outfld('SOLSD   ',cam_out%solsd ,pcols,lchnk)
   call outfld('SOLLD   ',cam_out%solld ,pcols,lchnk)
   call outfld('FSN200  ',rd%fsn200     ,pcols,lchnk)
   call outfld('FSN200C ',rd%fsn200c    ,pcols,lchnk)
   call outfld('FSNR'    ,rd%fsnr       ,pcols,lchnk)
   call outfld('SWCF    ',rd%swcf       ,pcols,lchnk)

   call outfld('TOT_CLD_VISTAU    ',rd%tot_cld_vistau     ,pcols,lchnk)
   call outfld('TOT_ICLD_VISTAU   ',rd%tot_icld_vistau    ,pcols,lchnk)
   call outfld('LIQ_ICLD_VISTAU   ',rd%liq_icld_vistau ,pcols,lchnk)
   call outfld('ICE_ICLD_VISTAU   ',rd%ice_icld_vistau ,pcols,lchnk)
   
end subroutine radiation_output_sw

!===============================================================================
  
subroutine radiation_output_lw(state,  rd, cam_out, flns, flnt, qrl)

   ! Dump longwave radiation information to history tape buffer (diagnostics)

   type(physics_state), intent(in) :: state
   type(rad_out_t),     intent(in) :: rd
   type(cam_out_t),     intent(in) :: cam_out
   real(r8),            intent(in) :: flns(pcols) ! Srf longwave cooling (up-down) flux
   real(r8),            intent(in) :: flnt(pcols) ! Net outgoing lw flux at model top
   real(r8),            pointer    :: qrl(:,:)    ! longwave  radiative heating rate

   ! Local variables
   integer :: lchnk, ncol
   real(r8) :: ftem(pcols,pver)
   !----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol

   call outfld('QRL     ',qrl(:ncol,:)/cpair,ncol,lchnk)
   call outfld('QRLC    ',rd%qrlc(:ncol,:)/cpair,ncol,lchnk)
   call outfld('FLNT    ',flnt  ,pcols,lchnk)
   call outfld('FLUT    ',rd%flut  ,pcols,lchnk)
   call outfld('FLUTC   ',rd%flutc ,pcols,lchnk)
   call outfld('FLNTC   ',rd%flntc ,pcols,lchnk)
   call outfld('FLNS    ',flns  ,pcols,lchnk)
   call outfld('FLDS    ',cam_out%flwds ,pcols,lchnk)
   call outfld('FLNSC   ',rd%flnsc ,pcols,lchnk)
   call outfld('FLDSC   ',rd%fldsc ,pcols,lchnk)
   call outfld('LWCF    ',rd%lwcf  ,pcols,lchnk)
   call outfld('FLN200  ',rd%fln200,pcols,lchnk)
   call outfld('FLN200C ',rd%fln200c,pcols,lchnk)
   call outfld('FLNR '   ,rd%flnr,pcols,lchnk)

end subroutine radiation_output_lw

!===============================================================================

subroutine radinp(ncol, pmid, pint, pmidrd, pintrd, eccf)

   use shr_orb_mod
   use time_manager, only: get_curr_calday

   !------------------------------Arguments--------------------------------
   integer, intent(in)   :: ncol                ! number of atmospheric columns

   real(r8), intent(in)  :: pmid(pcols,pver)    ! Pressure at model mid-levels (pascals)
   real(r8), intent(in)  :: pint(pcols,pverp)   ! Pressure at model interfaces (pascals)

   real(r8), intent(out) :: pmidrd(pcols,pver)  ! Pressure at mid-levels (dynes/cm*2)
   real(r8), intent(out) :: pintrd(pcols,pverp) ! Pressure at interfaces (dynes/cm*2)
   real(r8), intent(out) :: eccf                ! Earth-sun distance factor

   !---------------------------Local variables-----------------------------
   integer :: i, k
   real(r8) :: calday       ! current calendar day
   real(r8) :: delta        ! Solar declination angle
   !-----------------------------------------------------------------------

   calday = get_curr_calday()
   call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                      delta   ,eccf)

   ! Convert pressure from pascals to dynes/cm2
   do k=1,pver
      do i=1,ncol
         pmidrd(i,k) = pmid(i,k)*10.0_r8
         pintrd(i,k) = pint(i,k)*10.0_r8
      end do
   end do
   do i=1,ncol
      pintrd(i,pverp) = pint(i,pverp)*10.0_r8
   end do
   
end subroutine radinp

!===============================================================================

subroutine calc_col_mean(state, mmr_pointer, mean_value)

   ! Compute the column mean.  

   use cam_logfile,  only: iulog

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

