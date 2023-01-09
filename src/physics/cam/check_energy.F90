
module check_energy

!---------------------------------------------------------------------------------
! Purpose:
!
! Module to check
!   1. vertically integrated total energy and water conservation for each
!      column within the physical parameterizations
!
!   2. global mean total energy conservation between the physics output state
!      and the input state on the next time step.
!
!   3. add a globally uniform heating term to account for any change of total energy in 2.
!
! Author: Byron Boville  Oct 31, 2002
!
! Modifications:
!   03.03.29  Boville  Add global energy check and fixer.
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver, begchunk, endchunk
  use spmd_utils,      only: masterproc

  use gmean_mod,       only: gmean
  use physconst,       only: gravit, latvap, latice, cpair, rair
  use air_composition, only: cpairv, rairv
  use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_ptend_init
  use constituents,    only: cnst_get_ind, pcnst, cnst_name, cnst_get_type_byind
  use time_manager,    only: is_first_step
  use cam_logfile,     only: iulog

  implicit none
  private

! Public types:
  public check_tracers_data

! Public methods
  public :: check_energy_readnl    ! read namelist values
  public :: check_energy_register  ! register fields in physics buffer
  public :: check_energy_get_integrals ! get energy integrals computed in check_energy_gmean
  public :: check_energy_init      ! initialization of module
  public :: check_energy_timestep_init  ! timestep initialization of energy integrals and cumulative boundary fluxes
  public :: check_energy_chng      ! check changes in integrals against cumulative boundary fluxes
  public :: check_energy_gmean     ! global means of physics input and output total energy
  public :: check_energy_budgets_init ! initialization of energy budgets (addflds and budget_adds)
  public :: check_energy_budget_state_init ! initialization of energy budget integrals
  public :: check_energy_phys_budget_update    ! global budgets of physics energies
  public :: check_energy_phys_cnt_update    ! global budgets of physics energies
  public :: check_energy_fix       ! add global mean energy difference as a heating
  public :: check_tracers_init      ! initialize tracer integrals and cumulative boundary fluxes
  public :: check_tracers_chng      ! check changes in integrals against cumulative boundary fluxes

  public :: calc_te_and_aam_budgets ! calculate and output total energy and axial angular momentum diagnostics

  integer, public, parameter                         :: num_stages = 8
  character (len = 4), dimension(num_stages) :: stage = (/"phBF","phBP","phAP","phAM","dyBF","dyBP","dyAP","dyAM"/)
  character (len = 45),dimension(num_stages) :: stage_txt = (/&
       " before energy fixer                     ",& !phBF - physics energy
       " before parameterizations                ",& !phBF - physics energy
       " after parameterizations                 ",& !phAP - physics energy
       " after dry mass correction               ",& !phAM - physics energy
       " before energy fixer (dycore)            ",& !dyBF - dynamics energy
       " before parameterizations (dycore)       ",& !dyBF - dynamics energy
       " after parameterizations (dycore)        ",& !dyAP - dynamics energy
       " after dry mass correction (dycore)      " & !dyAM - dynamics energy
       /)

! Private module data

  logical  :: print_energy_errors = .false.

  real(r8) :: teout_glob           ! global mean energy of output state
  real(r8) :: teinp_glob           ! global mean energy of input state
  real(r8) :: tedif_glob           ! global mean energy difference
  real(r8) :: psurf_glob           ! global mean surface pressure
  real(r8) :: ptopb_glob           ! global mean top boundary pressure
  real(r8) :: heat_glob            ! global mean heating rate

! Physics buffer indices

  integer  :: teout_idx  = 0       ! teout index in physics buffer
  integer  :: dtcore_idx = 0       ! dtcore index in physics buffer
  integer  :: dqcore_idx = 0       ! dqcore index in physics buffer
  integer  :: ducore_idx = 0       ! ducore index in physics buffer
  integer  :: dvcore_idx = 0       ! dvcore index in physics buffer

  type check_tracers_data
     real(r8) :: tracer(pcols,pcnst)       ! initial vertically integrated total (kinetic + static) energy
     real(r8) :: tracer_tnd(pcols,pcnst)   ! cumulative boundary flux of total energy
     integer :: count(pcnst)               ! count of values with significant imbalances
  end type check_tracers_data


!===============================================================================
contains
!===============================================================================

subroutine check_energy_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_logical
   use cam_abortutils,  only: endrun

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: sub = 'check_energy_readnl'

   namelist /check_energy_nl/ print_energy_errors
   !-----------------------------------------------------------------------------

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'check_energy_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, check_energy_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub//': FATAL: reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(print_energy_errors, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: print_energy_errors")

   if (masterproc) then
      write(iulog,*) 'check_energy options:'
      write(iulog,*) '  print_energy_errors =', print_energy_errors
   end if

end subroutine check_energy_readnl

!===============================================================================

  subroutine check_energy_register()
!
! Register fields in the physics buffer.
!
!-----------------------------------------------------------------------

    use physics_buffer, only : pbuf_add_field, dtype_r8, dyn_time_lvls
    use physics_buffer, only : pbuf_register_subcol
    use subcol_utils,   only : is_subcol_on

!-----------------------------------------------------------------------

! Request physics buffer space for fields that persist across timesteps.

    call pbuf_add_field('TEOUT', 'global',dtype_r8 , (/pcols,dyn_time_lvls/),      teout_idx)
    call pbuf_add_field('DTCORE','global',dtype_r8,  (/pcols,pver,dyn_time_lvls/),dtcore_idx)
    ! DQCORE refers to dycore tendency of water vapor
    call pbuf_add_field('DQCORE','global',dtype_r8,  (/pcols,pver,dyn_time_lvls/),dqcore_idx)
    call pbuf_add_field('DUCORE','global',dtype_r8,  (/pcols,pver,dyn_time_lvls/),ducore_idx)
    call pbuf_add_field('DVCORE','global',dtype_r8,  (/pcols,pver,dyn_time_lvls/),dvcore_idx)
    if(is_subcol_on()) then
      call pbuf_register_subcol('TEOUT', 'phys_register', teout_idx)
      call pbuf_register_subcol('DTCORE', 'phys_register', dtcore_idx)
      call pbuf_register_subcol('DQCORE', 'phys_register', dqcore_idx)
      call pbuf_register_subcol('DUCORE', 'phys_register', ducore_idx)
      call pbuf_register_subcol('DVCORE', 'phys_register', dvcore_idx)
    end if

  end subroutine check_energy_register

!===============================================================================

subroutine check_energy_get_integrals( tedif_glob_out, heat_glob_out )

!-----------------------------------------------------------------------
! Purpose: Return energy integrals
!-----------------------------------------------------------------------

     real(r8), intent(out), optional :: tedif_glob_out
     real(r8), intent(out), optional :: heat_glob_out

!-----------------------------------------------------------------------

   if ( present(tedif_glob_out) ) then
      tedif_glob_out = tedif_glob
   endif
   if ( present(heat_glob_out) ) then
      heat_glob_out = heat_glob
   endif

end subroutine check_energy_get_integrals

!================================================================================================

  subroutine check_energy_init()
!
! Initialize the energy conservation module
!
!-----------------------------------------------------------------------
    use cam_history,       only: addfld, add_default, horiz_only
    use phys_control,      only: phys_getopts

    implicit none

    logical          :: history_budget, history_waccm
    integer          :: history_budget_histfile_num ! output history file number for budget fields
    integer          :: m                           ! budget array index into te_budgets
    character(len=32):: budget_name     ! budget names
    character(len=3) :: budget_pkgtype     ! budget type phy or dyn
    character(len=128):: budget_longname ! long name of budgets
!-----------------------------------------------------------------------

    call phys_getopts( history_budget_out = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num, &
                       history_waccm_out = history_waccm )

! register history variables
    call addfld('TEINP',  horiz_only,  'A', 'J/m2', 'Total energy of physics input')
    call addfld('TEOUT',  horiz_only,  'A', 'J/m2', 'Total energy of physics output')
    call addfld('TEFIX',  horiz_only,  'A', 'J/m2', 'Total energy after fixer')
    call addfld('EFIX',   horiz_only,  'A', 'W/m2', 'Effective sensible heat flux due to energy fixer')
    call addfld('DTCORE', (/ 'lev' /), 'A', 'K/s' , 'T tendency due to dynamical core')
    call addfld('DQCORE', (/ 'lev' /), 'A', 'kg/kg/s' , 'Water vapor tendency due to dynamical core')

    if ( history_budget ) then
       call add_default ('DTCORE', history_budget_histfile_num, ' ')
    end if
    if ( history_waccm ) then
       call add_default ('DTCORE', 1, ' ')
    end if

  end subroutine check_energy_init

!===============================================================================

  subroutine check_energy_timestep_init(state, tend, pbuf, col_type)
    use cam_thermo,      only: get_hydrostatic_energy
    use physics_buffer,  only: physics_buffer_desc, pbuf_set_field
    use cam_abortutils,  only: endrun
    use dyn_tests_utils, only: vc_physics, vc_dycore, vc_height
    use physics_types,   only: phys_te_idx, dyn_te_idx
!-----------------------------------------------------------------------
! Compute initial values of energy and water integrals,
! zero cumulative tendencies
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state),   intent(inout)    :: state
    type(physics_tend ),   intent(inout)    :: tend
    type(physics_buffer_desc), pointer      :: pbuf(:)
    integer, optional                       :: col_type  ! Flag inidicating whether using grid or subcolumns
!---------------------------Local storage-------------------------------
    real(r8)              :: cp_or_cv(state%psetcols,pver)
    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol

    ! cp_or_cv needs to be allocated to a size which matches state and ptend
    ! If psetcols == pcols, cpairv is the correct size and just copy into cp_or_cv
    ! If psetcols > pcols and all cpairv match cpair, then assign the constant cpair

    if (state%psetcols == pcols) then
       cp_or_cv(:,:) = cpairv(:,:,lchnk)
    else if (state%psetcols > pcols .and. all(cpairv(:,:,lchnk) == cpair)) then
       cp_or_cv(:,:) = cpair
    else
       call endrun('check_energy_timestep_init: cpairv is not allowed to vary when subcolumns are turned on')
    end if
    !
    ! CAM physics total energy
    !
    call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),&
         state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
         state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), state%T(1:ncol,1:pver),     &
         vc_physics, ps = state%ps(1:ncol), phis = state%phis(1:ncol),               &
         te = state%te_ini(1:ncol,phys_te_idx), H2O = state%tw_ini(1:ncol,phys_te_idx))
    !
    ! Dynamical core total energy
    !
    state%temp_ini(:ncol,:) = state%T(:ncol,:)
    state%z_ini(:ncol,:)    = state%zm(:ncol,:)
    if (vc_dycore == vc_height) then
      !
      ! compute cv if vertical coordinate is height: cv = cp - R
      !
      if (state%psetcols == pcols) then
        cp_or_cv(:,:) = cpairv(:,:,lchnk)-rairv(:,:,lchnk)
      else
        cp_or_cv(:,:) = cpair-rair
      endif

      call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),&
           state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
           state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), state%T(1:ncol,1:pver),     &
           vc_dycore, ps = state%ps(1:ncol), phis = state%phis(1:ncol),                &
           z_mid = state%z_ini(1:ncol,:),                                                  &
           te = state%te_ini(1:ncol,dyn_te_idx), H2O = state%tw_ini(1:ncol,dyn_te_idx))
    else
      state%te_ini(1:ncol,dyn_te_idx) = state%te_ini(1:ncol,phys_te_idx)
      state%tw_ini(1:ncol,dyn_te_idx) = state%tw_ini(1:ncol,phys_te_idx)
    end if

    state%te_cur(:ncol,:) = state%te_ini(:ncol,:)
    state%tw_cur(:ncol,:) = state%tw_ini(:ncol,:)

! zero cummulative boundary fluxes
    tend%te_tnd(:ncol) = 0._r8
    tend%tw_tnd(:ncol) = 0._r8

    state%count = 0

! initialize physics buffer
    if (is_first_step()) then
       call pbuf_set_field(pbuf, teout_idx, state%te_ini(:,dyn_te_idx), col_type=col_type)
    end if

  end subroutine check_energy_timestep_init

!===============================================================================

  subroutine check_energy_budget_state_init(state)
!-----------------------------------------------------------------------
! Compute initial values of energy and water integrals,
! zero cumulative tendencies
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------

    type(physics_state),   intent(inout)    :: state
!-----------------------------------------------------------------------

! zero cummulative boundary fluxes
    state%te_budgets(:,:,:) = 0._r8

  end subroutine check_energy_budget_state_init

!===============================================================================

  subroutine check_energy_budgets_init()
!-----------------------------------------------------------------------
! Compute initial values of energy and water integrals,
! zero cumulative tendencies
!-----------------------------------------------------------------------
    use budgets,            only: budget_add, budget_info, budget_num
    use cam_history,        only: addfld, horiz_only
    use cam_thermo,         only: thermo_budget_num_vars,thermo_budget_vars, &
                                  thermo_budget_vars_descriptor,thermo_budget_vars_unit
!---------------------------Local storage-------------------------------
    !
    ! variables for energy diagnostics
    !
    integer                                    :: istage, ivars, i
    character (len=256)                        :: str1, str2, str3
    character(len=32)                          :: budget_name     ! budget names
    character(len=3)                           :: budget_pkgtype,budget_optype
    character(len=128)                         :: budget_longname ! long name of budgets
!-----------------------------------------------------------------------

!
! energy diagnostics addflds for vars_stage combinations plus budget_adds for
! just the stages as the vars portion is accounted for via an extra array
! dimension in the state%te_budgets array.
!
    do istage = 1, num_stages
      do ivars=1, thermo_budget_num_vars
        write(str1,*) TRIM(ADJUSTL(thermo_budget_vars(ivars))),"_",TRIM(ADJUSTL(stage(istage)))
        write(str2,*) TRIM(ADJUSTL(thermo_budget_vars_descriptor(ivars)))," ", &
                           TRIM(ADJUSTL(stage_txt(istage)))
        write(str3,*) TRIM(ADJUSTL(thermo_budget_vars_unit(ivars)))
        call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'A', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)))
      end do
      call budget_add(TRIM(ADJUSTL(stage(istage))),'phy',longname=TRIM(ADJUSTL(stage_txt(istage))),outfld=.true.)
      write(iulog,*)'Calling addfld for ',TRIM(ADJUSTL(stage(istage)))
      call addfld (TRIM(ADJUSTL(stage(istage))),   horiz_only, 'A', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(stage_txt(istage))))
    end do

    ! Create budgets that are a sum/dif of 2 stages

    call budget_add('BP_param_and_efix','phAP','phBF','phy','dif',longname='dE/dt CAM physics parameterizations + efix dycore E (phAP-phBF)',outfld=.true.)
    call budget_add('BD_param_and_efix','dyAP','dyBF','phy','dif',longname='dE/dt CAM physics parameterizations + efix dycore E (dyAP-dyBF)',outfld=.true.)
    call budget_add('BP_phy_params','phAP','phBP','phy','dif',longname='dE/dt CAM physics parameterizations (phAP-phBP)',outfld=.true.)
    call budget_add('BD_phy_params','dyAP','dyBP','phy','dif',longname='dE/dt CAM physics parameterizations using dycore E (dyAP-dyBP)',outfld=.true.)
    call budget_add('BP_pwork','phAM','phAP','phy','dif',longname='dE/dt dry mass adjustment (phAM-phAP)',outfld=.true.)
    call budget_add('BD_pwork','dyAM','dyAP','phy','dif',longname='dE/dt dry mass adjustment using dycore E (dyAM-dyAP)',outfld=.true.)
    call budget_add('BP_efix','phBP','phBF','phy','dif',longname='dE/dt energy fixer (phBP-phBF)',outfld=.true.)
    call budget_add('BD_efix','dyBP','dyBF','phy','dif',longname='dE/dt energy fixer using dycore E (dyBP-dyBF)',outfld=.true.)
    call budget_add('BP_phys_tot','phAM','phBF','phy','dif',longname='dE/dt physics total (phAM-phBF)',outfld=.true.)
    call budget_add('BD_phys_tot','dyAM','dyBF','phy','dif',longname='dE/dt physics total using dycore E (dyAM-dyBF)',outfld=.true.)

    ! create addfld calls for all two stage budgets
    do i=1,budget_num
       call budget_info(i,name=budget_name,longname=budget_longname,pkgtype=budget_pkgtype,optype=budget_optype)
       if (budget_pkgtype=='phy'.and.(budget_optype=='dif'.or.budget_optype=='sum')) then
          do ivars=1, thermo_budget_num_vars
             write(str1,*) TRIM(ADJUSTL(thermo_budget_vars(ivars))),"_",TRIM(ADJUSTL(budget_name))
             write(str2,*) TRIM(ADJUSTL(thermo_budget_vars_descriptor(ivars)))," ", &
                  TRIM(ADJUSTL(budget_longname))
             write(str3,*) TRIM(ADJUSTL(thermo_budget_vars_unit(ivars)))
             write(iulog,*)'Calling addfld for ',TRIM(ADJUSTL(str1))
             call addfld (TRIM(ADJUSTL(str1)),   horiz_only, 'A', TRIM(ADJUSTL(str3)),TRIM(ADJUSTL(str2)))
             if (TRIM(ADJUSTL(thermo_budget_vars(ivars)))=='TE') then
                write(iulog,*)'Calling addfld for ',TRIM(ADJUSTL(budget_name))
                call addfld (TRIM(ADJUSTL(budget_name)),   horiz_only, 'A', 'J/m2',TRIM(ADJUSTL(budget_longname)))
             end if
          end do
       end if
    end do

  end subroutine check_energy_budgets_init

!===============================================================================

  subroutine check_energy_chng(state, tend, name, nstep, ztodt,        &
       flx_vap, flx_cnd, flx_ice, flx_sen)
    use cam_thermo,      only: get_hydrostatic_energy
    use dyn_tests_utils, only: vc_physics, vc_dycore, vc_height
    use cam_abortutils,  only: endrun
    use physics_types,   only: phys_te_idx, dyn_te_idx
!-----------------------------------------------------------------------
! Check that the energy and water change matches the boundary fluxes
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state)    , intent(inout) :: state
    type(physics_tend )    , intent(inout) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: flx_vap(:)          ! (pcols) - boundary flux of vapor         (kg/m2/s)
    real(r8), intent(in   ) :: flx_cnd(:)          ! (pcols) -boundary flux of liquid+ice    (m/s) (precip?)
    real(r8), intent(in   ) :: flx_ice(:)          ! (pcols) -boundary flux of ice           (m/s) (snow?)
    real(r8), intent(in   ) :: flx_sen(:)          ! (pcols) -boundary flux of sensible heat (w/m2)

!******************** BAB ******************************************************
!******* Note that the precip and ice fluxes are in precip units (m/s). ********
!******* I would prefer to have kg/m2/s.                                ********
!******* I would also prefer liquid (not total) and ice fluxes          ********
!*******************************************************************************

!---------------------------Local storage-------------------------------

    real(r8) :: te_xpd(state%ncol)                 ! expected value (f0 + dt*boundary_flux)
    real(r8) :: te_dif(state%ncol)                 ! energy of input state - original energy
    real(r8) :: te_tnd(state%ncol)                 ! tendency from last process
    real(r8) :: te_rer(state%ncol)                 ! relative error in energy column

    real(r8) :: tw_xpd(state%ncol)                 ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tw_dif(state%ncol)                 ! tw_inp - original water
    real(r8) :: tw_tnd(state%ncol)                 ! tendency from last process
    real(r8) :: tw_rer(state%ncol)                 ! relative error in water column

    real(r8) :: te(state%ncol)                     ! vertical integral of total energy
    real(r8) :: tw(state%ncol)                     ! vertical integral of total water
    real(r8) :: cp_or_cv(state%psetcols,pver)      ! cp or cv depending on vcoord
    real(r8) :: scaling(state%psetcols,pver)       ! scaling for conversion of temperature increment
    real(r8) :: temp(state%ncol,pver)              ! temperature

    real(r8) :: se(pcols)                          ! enthalpy or internal energy (J/m2)
    real(r8) :: po(pcols)                          ! surface potential or potential energy (J/m2)
    real(r8) :: ke(pcols)                          ! kinetic energy    (J/m2)
    real(r8) :: wv(pcols)                          ! column integrated vapor       (kg/m2)
    real(r8) :: liq(pcols)                         ! column integrated liquid      (kg/m2)
    real(r8) :: ice(pcols)                         ! column integrated ice         (kg/m2)

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i                                   ! column
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol

    ! If psetcols == pcols, cpairv is the correct size and just copy into cp_or_cv
    ! If psetcols > pcols and all cpairv match cpair, then assign the constant cpair

    if (state%psetcols == pcols) then
       cp_or_cv(:,:) = cpairv(:,:,lchnk)
    else if (state%psetcols > pcols .and. all(cpairv(:,:,:) == cpair)) then
       cp_or_cv(:,:) = cpair
    else
       call endrun('check_energy_chng: cpairv is not allowed to vary when subcolumns are turned on')
    end if

    call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),&
         state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
         state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), state%T(1:ncol,1:pver),     &
         vc_physics, ps = state%ps(1:ncol), phis = state%phis(1:ncol),               &
         te = te, H2O = tw, se=se,po=po,ke=ke,wv=wv,liq=liq,ice=ice)
    ! compute expected values and tendencies
    do i = 1, ncol
       ! change in static energy and total water
       te_dif(i) = te(i) - state%te_cur(i,phys_te_idx)
       tw_dif(i) = tw(i) - state%tw_cur(i,phys_te_idx)

       ! expected tendencies from boundary fluxes for last process
       te_tnd(i) = flx_vap(i)*(latvap+latice) - (flx_cnd(i) - flx_ice(i))*1000._r8*latice + flx_sen(i)
       tw_tnd(i) = flx_vap(i) - flx_cnd(i) *1000._r8

       ! cummulative tendencies from boundary fluxes
       tend%te_tnd(i) = tend%te_tnd(i) + te_tnd(i)
       tend%tw_tnd(i) = tend%tw_tnd(i) + tw_tnd(i)

       ! expected new values from previous state plus boundary fluxes
       te_xpd(i) = state%te_cur(i,phys_te_idx) + te_tnd(i)*ztodt
       tw_xpd(i) = state%tw_cur(i,phys_te_idx) + tw_tnd(i)*ztodt

       ! relative error, expected value - input state / previous state
       te_rer(i) = (te_xpd(i) - te(i)) / state%te_cur(i,phys_te_idx)
    end do

    ! relative error for total water (allow for dry atmosphere)
    tw_rer = 0._r8
    where (state%tw_cur(:ncol,phys_te_idx) > 0._r8)
       tw_rer(:ncol) = (tw_xpd(:ncol) - tw(:ncol)) / state%tw_cur(:ncol,1)
    end where

    ! error checking
    if (print_energy_errors) then
       if (any(abs(te_rer(1:ncol)) > 1.E-14_r8 .or. abs(tw_rer(1:ncol)) > 1.E-10_r8)) then
          do i = 1, ncol
             ! the relative error threshold for the water budget has been reduced to 1.e-10
             ! to avoid messages generated by QNEG3 calls
             ! PJR- change to identify if error in energy or water
             if (abs(te_rer(i)) > 1.E-14_r8 ) then
                state%count = state%count + 1
                write(iulog,*) "significant energy conservation error after ", name,        &
                      " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
                write(iulog,*) te(i),te_xpd(i),te_dif(i),tend%te_tnd(i)*ztodt,  &
                      te_tnd(i)*ztodt,te_rer(i)
             endif
             if ( abs(tw_rer(i)) > 1.E-10_r8) then
                state%count = state%count + 1
                write(iulog,*) "significant water conservation error after ", name,        &
                      " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
                write(iulog,*) tw(i),tw_xpd(i),tw_dif(i),tend%tw_tnd(i)*ztodt,  &
                      tw_tnd(i)*ztodt,tw_rer(i)
             end if
          end do
       end if
    end if

    ! copy new value to state

    do i = 1, ncol
      state%te_cur(i,phys_te_idx) = te(i)
      state%tw_cur(i,phys_te_idx) = tw(i)
    end do

    !
    ! Dynamical core total energy
    !
    if (vc_dycore == vc_height) then
      !
      ! compute cv if vertical coordinate is height: cv = cp - R
      !
      ! Note: cp_or_cv set above for pressure coordinate
      !
      if (state%psetcols == pcols) then
        cp_or_cv(:,:) = cpairv(:,:,lchnk)-rairv(:,:,lchnk)
      else
        cp_or_cv(:,:) = cpair-rair
      endif
      scaling(:,:) = cpairv(:,:,lchnk)/cp_or_cv(:,:) !cp/cv scaling

      temp(1:ncol,:) = state%temp_ini(1:ncol,:)+scaling(1:ncol,:)*(state%T(1:ncol,:)-state%temp_ini(1:ncol,:))
      call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),&
           state%pdel(1:ncol,1:pver), cp_or_cv(1:ncol,1:pver),                         &
           state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), temp(1:ncol,1:pver),        &
           vc_dycore, ps = state%ps(1:ncol), phis = state%phis(1:ncol),                &
           z_mid = state%z_ini(1:ncol,:),                                                  &
           te = state%te_cur(1:ncol,dyn_te_idx), H2O = state%tw_cur(1:ncol,dyn_te_idx))
    else
      state%te_cur(1:ncol,dyn_te_idx) = te(1:ncol)
      state%tw_cur(1:ncol,dyn_te_idx) = tw(1:ncol)
    end if
  end subroutine check_energy_chng


  subroutine check_energy_gmean(state, pbuf2d, dtime, nstep)

    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    use physics_types,   only: dyn_te_idx
!-----------------------------------------------------------------------
! Compute global mean total energy of physics input and output states
! computed consistently with dynamical core vertical coordinate
! (under hydrostatic assumption)
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
    type(physics_buffer_desc),    pointer    :: pbuf2d(:,:)

    real(r8), intent(in) :: dtime        ! physics time step
    integer , intent(in) :: nstep        ! current timestep number

!---------------------------Local storage-------------------------------
    integer :: ncol                      ! number of active columns
    integer :: lchnk                     ! chunk index

    real(r8) :: te(pcols,begchunk:endchunk,4)
                                         ! total energy of input/output states (copy)
    real(r8) :: te_glob(4)               ! global means of total energy
    real(r8), pointer :: teout(:)
!-----------------------------------------------------------------------

    ! Copy total energy out of input and output states
    do lchnk = begchunk, endchunk
       ncol = state(lchnk)%ncol
       ! input energy using dynamical core energy formula
       te(:ncol,lchnk,1) = state(lchnk)%te_ini(:ncol,dyn_te_idx)
       ! output energy
       call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk),teout_idx, teout)

       te(:ncol,lchnk,2) = teout(1:ncol)
       ! surface pressure for heating rate
       te(:ncol,lchnk,3) = state(lchnk)%pint(:ncol,pver+1)
       ! model top pressure for heating rate (not constant for z-based vertical coordinate!)
       te(:ncol,lchnk,4) = state(lchnk)%pint(:ncol,1)
    end do

    ! Compute global means of input and output energies and of
    ! surface pressure for heating rate (assume uniform ptop)
    call gmean(te, te_glob, 4)

    if (begchunk .le. endchunk) then
       teinp_glob = te_glob(1)
       teout_glob = te_glob(2)
       psurf_glob = te_glob(3)
       ptopb_glob = te_glob(4)

       ! Global mean total energy difference
       tedif_glob =  teinp_glob - teout_glob
       heat_glob  = -tedif_glob/dtime * gravit / (psurf_glob - ptopb_glob)
       if (masterproc) then
          write(iulog,'(1x,a9,1x,i8,5(1x,e25.17))') "nstep, te", nstep, teinp_glob, teout_glob, &
               heat_glob, psurf_glob, ptopb_glob
       end if
    else
       heat_glob = 0._r8
    end if  !  (begchunk .le. endchunk)

  end subroutine check_energy_gmean

  subroutine check_energy_phys_budget_update(state, dtime, nstep)

    use cam_history,     only: outfld
    use budgets,         only: budget_num, budget_info, &
                               budget_outfld, budget_num_phy, &
                               budget_put_global
    use cam_thermo,      only: thermo_budget_num_vars, thermo_budget_vars_massv
    use cam_abortutils,  only: endrun
    use dycore_budget,   only: print_budget
!-----------------------------------------------------------------------
! Compute global mean total energy of physics input and output states
! computed consistently with dynamical core vertical coordinate
! (under hydrostatic assumption)
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), intent(inout   ), dimension(begchunk:endchunk) :: state

    real(r8), intent(in) :: dtime        ! physics time step
    integer , intent(in) :: nstep        ! current timestep number

!---------------------------Local storage-------------------------------
    integer :: ncol                      ! number of active columns
    integer :: lchnk                     ! chunk index

    real(r8),allocatable :: te(:,:,:,:) ! total energy of input/output states (copy)
    real(r8),allocatable :: te_glob(:,:) ! global means of total energy
    integer  :: i,ii,s_ind,is1,is2,is1b,is2b
    character*32 :: budget_name               ! parameterization name for fluxes
    character*3 :: budget_pkgtype             ! parameterization type phy or dyn
    character*3 :: budget_optype              ! dif or stg
!-----------------------------------------------------------------------
    if (.not.allocated (te)) then
       allocate( te(pcols,begchunk:endchunk,budget_num_phy,thermo_budget_num_vars))
    end if
    if (.not.allocated (te_glob)) then
       allocate( te_glob(budget_num_phy,thermo_budget_num_vars))
    else
       write(iulog,*)'no alloc call shape te_glob=',shape(te_glob)
    end if
    te=0.0_r8
    te_glob=0.0_r8

    ! calculate energy budget differences
    do lchnk = begchunk, endchunk
       ncol = state(lchnk)%ncol
       do ii=1,budget_num
          call budget_info(ii,name=budget_name,pkgtype=budget_pkgtype,optype=budget_optype,state_ind=i)
          if (budget_pkgtype=='phy') then
             if (budget_optype=='dif') then
                call budget_info(ii,stg1stateidx=is1, stg2stateidx=is2,stg1index=is1b,stg2index=is2b)
                if (state(lchnk)%budget_cnt(is1b).ne.state(lchnk)%budget_cnt(is2b)) then
                   if (lchnk==begchunk.and.masterproc) write(iulog,*)'budget_cnt mismatch stage1=',state(lchnk)%budget_cnt(is1b),'stage2=',state(lchnk)%budget_cnt(is2b)
                   call endrun()
                end if
                if (state(lchnk)%budget_cnt(is1b)==0.or.state(lchnk)%budget_cnt(is2b)==0) then
                   te(:,lchnk,i,:)=0._r8
                   if (lchnk==begchunk.and.masterproc) write(iulog,*)'zeroing:',budget_name,' cnt1:',state(lchnk)%budget_cnt(is1b),' cnt2 ',state(lchnk)%budget_cnt(is2b)
                else
                   if (lchnk==begchunk.and.masterproc) write(iulog,*)'dif and norm into te:',budget_name,' cnt:',state(lchnk)%budget_cnt(ii),' cnt1:',state(lchnk)%budget_cnt(is1b),'budget index ii/is1b=',ii,'/',is1b
                   te(:,lchnk,i,:) = (state(lchnk)%te_budgets(:,:,is1)-state(lchnk)%te_budgets(:,:,is2))/state(lchnk)%budget_cnt(is1b)
                end if
             else if (budget_optype=='sum') then
                call budget_info(ii,stg1stateidx=is1, stg2stateidx=is2,stg1index=is1b,stg2index=is2b)
                if (state(lchnk)%budget_cnt(is1b).ne.state(lchnk)%budget_cnt(is2b)) then
                   if (lchnk==begchunk.and.masterproc) write(iulog,*)'budget_cnt mismatch stage1=',state(lchnk)%budget_cnt(is1b),'stage2=',state(lchnk)%budget_cnt(is2b)
                   call endrun()
                end if
                if (state(lchnk)%budget_cnt(is1b)==0.or.state(lchnk)%budget_cnt(is2b)==0) then
                   te(:,lchnk,i,:)=0._r8
                else          
                   if (lchnk==begchunk.and.masterproc) write(iulog,*)'sum and norm into te:',budget_name,' cnt:',state(lchnk)%budget_cnt(ii),' cnt1:',state(lchnk)%budget_cnt(is1b),'budget index ii/is1b=',ii,'/',is1b
                   te(:,lchnk,i,:) = (state(lchnk)%te_budgets(:,:,is1)+state(lchnk)%te_budgets(:,:,is2))/state(lchnk)%budget_cnt(is1b)
                end if
             else
                if (state(lchnk)%budget_cnt(ii)==0) then
                   te(:,lchnk,i,:)=0._r8
                   if (lchnk==begchunk.and.masterproc) write(iulog,*)'zeroing:',budget_name,' cnt:',state(lchnk)%budget_cnt(ii),' ii=',ii,"current vals=",state(lchnk)%te_budgets(:,:,ii)
                else
                   if (lchnk==begchunk.and.masterproc) write(iulog,*)'norm and read into te:',budget_name,' cnt:',state(lchnk)%budget_cnt(ii),'budget index=',ii
                   te(:,lchnk,i,:)=state(lchnk)%te_budgets(:,:,i)/state(lchnk)%budget_cnt(ii)
                end if
             end if
             if (budget_outfld(i).and.budget_pkgtype=='phy') call outfld(trim(budget_name), te(:ncol,lchnk,i,1), pcols, lchnk)
          end if
       end do
    end do
    ! Compute global means of budgets
    do i=1,thermo_budget_num_vars
       call gmean(te(:,:,:,i), te_glob(:,i), budget_num_phy)
       !divide by time to get flux if not a mass budget
       if (.not.thermo_budget_vars_massv(i)) te_glob(:,i)=te_glob(:,i)/dtime
    end do

    if (masterproc) then
       do ii=1,budget_num
          call budget_info(ii,name=budget_name,pkgtype=budget_pkgtype,optype=budget_optype,state_ind=s_ind)
          if (budget_pkgtype=='phy') then
             do i=1,thermo_budget_num_vars
                call budget_put_global(trim(budget_name),i,te_glob(s_ind,i))
                if (budget_optype=='dif'.or.budget_optype=='sum') then
                   call budget_info(ii,stg1index=is1b)
                   write(iulog,*)"putting global ",trim(budget_name)," m_cnst=",i," ",te_glob(s_ind,i)," cnt=",state(begchunk)%budget_cnt(is1b),"is1b=",is1b
                else
                   write(iulog,*)"putting global ",trim(budget_name)," m_cnst=",i," ",te_glob(s_ind,i)," cnt=",state(begchunk)%budget_cnt(ii),"ii=",ii
                end if
             end do
          end if
       end do
    end if
    call print_budget()
    do lchnk = begchunk, endchunk
       state(lchnk)%budget_cnt(:)=0._r8
       state(lchnk)%te_budgets(:,:,:)=0._r8
    end do

  end subroutine check_energy_phys_budget_update
  subroutine check_energy_phys_cnt_update(state)

    use budgets,         only: budget_num, budget_info, &
                               budget_outfld, budget_num_phy, &
                               budget_put_global
    use cam_abortutils,  only: endrun
    use dycore_budget,   only: print_budget
!-----------------------------------------------------------------------
! Compute global mean total energy of physics input and output states
! computed consistently with dynamical core vertical coordinate
! (under hydrostatic assumption)
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), intent(inout   ), dimension(begchunk:endchunk) :: state

!---------------------------Local storage-------------------------------
    integer :: ncol                      ! number of active columns
    integer :: lchnk                     ! chunk index

    integer  :: i,ii,s_ind,is1,is2,is1b,is2b
    character*32 :: budget_name               ! parameterization name for fluxes
    character*3 :: budget_pkgtype             ! parameterization type phy or dyn
    character*3 :: budget_optype              ! dif or stg
!-----------------------------------------------------------------------
    do lchnk = begchunk, endchunk
       ncol = state(lchnk)%ncol
       do ii=1,budget_num
          call budget_info(ii,name=budget_name,pkgtype=budget_pkgtype,optype=budget_optype,state_ind=i)
          if (budget_pkgtype=='phy'.and.(budget_optype=='dif'.or.budget_optype=='sub')) then
             state(lchnk)%budget_cnt(ii)=state(lchnk)%budget_cnt(ii)+1
             call budget_info(ii,stg1index=is1b,stg2index=is2b)
             if (lchnk==begchunk .and. masterproc) write(iulog,*)trim(budget_name)," cnt(",ii,") updated to ",state(lchnk)%budget_cnt(ii),'stage1=',state(lchnk)%budget_cnt(is1b),'stage2=',state(lchnk)%budget_cnt(is2b),' is1b/is2b=',is1b,'/',is2b
          end if
       end do
    end do
  end subroutine check_energy_phys_cnt_update

!===============================================================================
  subroutine check_energy_fix(state, ptend, nstep, eshflx)

!-----------------------------------------------------------------------
! Add heating rate required for global mean total energy conservation
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), intent(in   ) :: state
    type(physics_ptend), intent(out)   :: ptend

    integer , intent(in   ) :: nstep          ! time step number
    real(r8), intent(out  ) :: eshflx(pcols)  ! effective sensible heat flux

!---------------------------Local storage-------------------------------
    integer  :: i                        ! column
    integer  :: ncol                     ! number of atmospheric columns in chunk
!-----------------------------------------------------------------------
    ncol = state%ncol

    call physics_ptend_init(ptend, state%psetcols, 'chkenergyfix', ls=.true.)

#if ( defined OFFLINE_DYN )
    ! disable the energy fix for offline driver
    heat_glob = 0._r8
#endif
! add (-) global mean total energy difference as heating
    ptend%s(:ncol,:pver) = heat_glob

! compute effective sensible heat flux
    do i = 1, ncol
       eshflx(i) = heat_glob * (state%pint(i,pver+1) - state%pint(i,1)) / gravit
    end do

    return
  end subroutine check_energy_fix


!===============================================================================
  subroutine check_tracers_init(state, tracerint)

!-----------------------------------------------------------------------
! Compute initial values of tracers integrals,
! zero cumulative tendencies
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------

    type(physics_state),   intent(in)    :: state
    type(check_tracers_data), intent(out)   :: tracerint

!---------------------------Local storage-------------------------------

    real(r8) :: tr(pcols)                          ! vertical integral of tracer
    real(r8) :: trpdel(pcols, pver)                ! pdel for tracer

    integer ncol                                   ! number of atmospheric columns
    integer  i,k,m                                 ! column, level,constituent indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: ixrain, ixsnow                      ! RAINQM and SNOWQM indices
    integer :: ixgrau                              ! GRAUQM index
!-----------------------------------------------------------------------

    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain,   abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow,   abort=.false.)
    call cnst_get_ind('GRAUQM', ixgrau,   abort=.false.)


    do m = 1,pcnst

       if ( any(m == (/ 1, ixcldliq, ixcldice, &
                           ixrain,   ixsnow, ixgrau /)) ) exit   ! dont process water substances
                                                                 ! they are checked in check_energy

       if (cnst_get_type_byind(m).eq.'dry') then
          trpdel(:ncol,:) = state%pdeldry(:ncol,:)
       else
          trpdel(:ncol,:) = state%pdel(:ncol,:)
       endif

       ! Compute vertical integrals of tracer
       tr = 0._r8
       do k = 1, pver
          do i = 1, ncol
             tr(i) = tr(i) + state%q(i,k,m)*trpdel(i,k)/gravit
          end do
       end do

       ! Compute vertical integrals of frozen static tracers and total water.
       do i = 1, ncol
          tracerint%tracer(i,m) = tr(i)
       end do

       ! zero cummulative boundary fluxes
       tracerint%tracer_tnd(:ncol,m) = 0._r8

       tracerint%count(m) = 0

    end do

    return
  end subroutine check_tracers_init

!===============================================================================
  subroutine check_tracers_chng(state, tracerint, name, nstep, ztodt, cflx)

!-----------------------------------------------------------------------
! Check that the tracers and water change matches the boundary fluxes
! these checks are not save when there are tracers transformations, as
! they only check to see whether a mass change in the column is
! associated with a flux
!-----------------------------------------------------------------------

    use cam_abortutils, only: endrun


    implicit none

!------------------------------Arguments--------------------------------

    type(physics_state)    , intent(in   ) :: state
    type(check_tracers_data), intent(inout) :: tracerint! tracers integrals and boundary fluxes
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: cflx(pcols,pcnst)       ! boundary flux of tracers       (kg/m2/s)

!---------------------------Local storage-------------------------------

    real(r8) :: tracer_inp(pcols,pcnst)                   ! total tracer of new (input) state
    real(r8) :: tracer_xpd(pcols,pcnst)                   ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tracer_dif(pcols,pcnst)                   ! tracer_inp - original tracer
    real(r8) :: tracer_tnd(pcols,pcnst)                   ! tendency from last process
    real(r8) :: tracer_rer(pcols,pcnst)                   ! relative error in tracer column

    real(r8) :: tr(pcols)                           ! vertical integral of tracer
    real(r8) :: trpdel(pcols, pver)                       ! pdel for tracer

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: ixrain, ixsnow                      ! RAINQM and SNOWQM indices
    integer :: ixgrau                              ! GRAUQM index
    integer :: m                            ! tracer index
    character(len=8) :: tracname   ! tracername
!-----------------------------------------------------------------------
!!$    if (.true.) return

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain,   abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow,   abort=.false.)
    call cnst_get_ind('GRAUQM', ixgrau,   abort=.false.)

    do m = 1,pcnst

       if ( any(m == (/ 1, ixcldliq, ixcldice, &
                           ixrain,   ixsnow, ixgrau /)) ) exit   ! dont process water substances
                                                                 ! they are checked in check_energy
       tracname = cnst_name(m)
       if (cnst_get_type_byind(m).eq.'dry') then
          trpdel(:ncol,:) = state%pdeldry(:ncol,:)
       else
          trpdel(:ncol,:) = state%pdel(:ncol,:)
       endif

       ! Compute vertical integrals tracers
       tr = 0._r8
       do k = 1, pver
          do i = 1, ncol
             tr(i) = tr(i) + state%q(i,k,m)*trpdel(i,k)/gravit
          end do
       end do

       ! Compute vertical integrals of tracer
       do i = 1, ncol
          tracer_inp(i,m) = tr(i)
       end do

       ! compute expected values and tendencies
       do i = 1, ncol
          ! change in tracers
          tracer_dif(i,m) = tracer_inp(i,m) - tracerint%tracer(i,m)

          ! expected tendencies from boundary fluxes for last process
          tracer_tnd(i,m) = cflx(i,m)

          ! cummulative tendencies from boundary fluxes
          tracerint%tracer_tnd(i,m) = tracerint%tracer_tnd(i,m) + tracer_tnd(i,m)

          ! expected new values from original values plus boundary fluxes
          tracer_xpd(i,m) = tracerint%tracer(i,m) + tracerint%tracer_tnd(i,m)*ztodt

          ! relative error, expected value - input value / original
          tracer_rer(i,m) = (tracer_xpd(i,m) - tracer_inp(i,m)) / tracerint%tracer(i,m)
       end do

!! final loop for error checking
!    do i = 1, ncol

!! error messages
!       if (abs(enrgy_rer(i)) > 1.E-14 .or. abs(water_rer(i)) > 1.E-14) then
!          tracerint%count = tracerint%count + 1
!          write(iulog,*) "significant conservations error after ", name,        &
!               " count", tracerint%count, " nstep", nstep, "chunk", lchnk, "col", i
!          write(iulog,*) enrgy_inp(i),enrgy_xpd(i),enrgy_dif(i),tracerint%enrgy_tnd(i)*ztodt,  &
!               enrgy_tnd(i)*ztodt,enrgy_rer(i)
!          write(iulog,*) water_inp(i),water_xpd(i),water_dif(i),tracerint%water_tnd(i)*ztodt,  &
!               water_tnd(i)*ztodt,water_rer(i)
!       end if
!    end do


       ! final loop for error checking
       if ( maxval(tracer_rer) > 1.E-14_r8 ) then
          write(iulog,*) "CHECK_TRACERS TRACER large rel error"
          write(iulog,*) tracer_rer
       endif

       do i = 1, ncol
          ! error messages
          if (abs(tracer_rer(i,m)) > 1.E-14_r8 ) then
             tracerint%count = tracerint%count + 1
             write(iulog,*) "CHECK_TRACERS TRACER significant conservation error after ", name,        &
                  " count", tracerint%count, " nstep", nstep, "chunk", lchnk, "col",i
             write(iulog,*)' process name, tracname, index ',  name, tracname, m
             write(iulog,*)" input integral              ",tracer_inp(i,m)
             write(iulog,*)" expected integral           ", tracer_xpd(i,m)
             write(iulog,*)" input - inital integral     ",tracer_dif(i,m)
             write(iulog,*)" cumulative tend      ",tracerint%tracer_tnd(i,m)*ztodt
             write(iulog,*)" process tend         ",tracer_tnd(i,m)*ztodt
             write(iulog,*)" relative error       ",tracer_rer(i,m)
             call endrun()
          end if
       end do
    end do

    return
  end subroutine check_tracers_chng

!#######################################################################

  subroutine calc_te_and_aam_budgets(state, outfld_name_suffix,vc)
    use physconst,       only: gravit,cpair,pi,rearth,omega
    use cam_thermo,      only: get_hydrostatic_energy,thermo_budget_num_vars,thermo_budget_vars, &
                               wvidx,wlidx,wiidx,seidx,keidx,moidx,mridx,ttidx,teidx,poidx
    use cam_history,     only: hist_fld_active, outfld
    use dyn_tests_utils, only: vc_physics, vc_height
    use cam_abortutils,  only: endrun
    use budgets,         only: budget_info_byname
    use cam_history_support, only: max_fieldname_len
!------------------------------Arguments--------------------------------

    type(physics_state), intent(inout) :: state
    character(len=*),    intent(in)    :: outfld_name_suffix ! suffix for "outfld"
    integer, optional,   intent(in)    :: vc                 ! vertical coordinate

!---------------------------Local storage-------------------------------
    real(r8) :: se(pcols)                          ! Dry Static energy (J/m2)
    real(r8) :: po(pcols)                          ! surface potential or potential energy (J/m2)
    real(r8) :: ke(pcols)                          ! kinetic energy    (J/m2)
    real(r8) :: wv(pcols)                          ! column integrated vapor       (kg/m2)
    real(r8) :: liq(pcols)                         ! column integrated liquid      (kg/m2)
    real(r8) :: ice(pcols)                         ! column integrated ice         (kg/m2)
    real(r8) :: tt(pcols)                          ! column integrated test tracer (kg/m2)
    real(r8) :: mr(pcols)                          ! column integrated wind axial angular momentum (kg*m2/s)
    real(r8) :: mo(pcols)                          ! column integrated mass axial angular momentum (kg*m2/s)
    real(r8) :: tt_tmp,mr_tmp,mo_tmp,cos_lat
    real(r8) :: mr_cnst, mo_cnst
    real(r8) :: cp_or_cv(pcols,pver)               ! cp for pressure-based vcoord and cv for height vcoord
    real(r8) :: temp(pcols,pver)                   ! temperature
    real(r8) :: scaling(pcols,pver)                ! scaling for conversion of temperature increment

    integer :: lchnk                               ! chunk identifier
    integer :: ncol                                ! number of atmospheric columns
    integer :: i,k                                 ! column, level indices
    integer :: vc_loc                              ! local vertical coordinate variable
    integer :: s_ind,b_ind                      ! budget array index
    integer :: ixtt                                ! test tracer index
    character(len=32) :: name_out1,name_out2,name_out3,name_out4,name_out5,name_out6,name_out7
    character(len=max_fieldname_len) :: name_out(thermo_budget_num_vars)

!-----------------------------------------------------------------------


    do i=1,thermo_budget_num_vars
       name_out(i)=trim(thermo_budget_vars(i))//'_'//trim(outfld_name_suffix)
    end do
    name_out1 = 'SE_'   //trim(outfld_name_suffix)
    name_out2 = 'KE_'   //trim(outfld_name_suffix)
    name_out3 = 'WV_'   //trim(outfld_name_suffix)
    name_out4 = 'WL_'   //trim(outfld_name_suffix)
    name_out5 = 'WI_'   //trim(outfld_name_suffix)
    name_out6 = 'TT_'   //trim(outfld_name_suffix)
    name_out7 = 'TE_'   //trim(outfld_name_suffix)

      lchnk = state%lchnk
      ncol  = state%ncol

      call budget_info_byname(trim(outfld_name_suffix),budget_ind=b_ind,state_ind=s_ind)

      if (present(vc)) then
        vc_loc = vc
      else
        vc_loc = vc_physics
      end if

      if (state%psetcols == pcols) then
        if (vc_loc == vc_height) then
          !
          ! compute cv if vertical coordinate is height: cv = cp - R
          !
          cp_or_cv(:,:) = cpairv(:,:,lchnk)-rairv(:,:,lchnk)!cv
        else
          cp_or_cv(:,:) = cpairv(:,:,lchnk)                 !cp
        end if
      else
        call endrun('calc_te_and_aam_budgets: energy diagnostics not implemented/tested for subcolumns')
      end if

      if (vc_loc == vc_height) then
        scaling(:,:) = cpairv(:,:,lchnk)/cp_or_cv(:,:) !cp/cv scaling for temperature increment under constant volume
      else
        scaling(:,:) = 1.0_r8
      end if
      ! scale accumulated temperature increment for constant volume (otherwise effectively do nothing)
      temp(1:ncol,:) = state%temp_ini(1:ncol,:)+scaling(1:ncol,:)*(state%T(1:ncol,:)- state%temp_ini(1:ncol,:))

      call get_hydrostatic_energy(state%q(1:ncol,1:pver,1:pcnst),&
           state%pdel(1:ncol,1:pver), cp_or_cv,                                        &
           state%u(1:ncol,1:pver), state%v(1:ncol,1:pver), temp(1:ncol,1:pver),        &
           vc_loc, ps = state%ps(1:ncol), phis = state%phis(1:ncol),                   &
           z_mid = state%z_ini(1:ncol,:), se = se, po = po, ke = ke, wv = wv, liq = liq, ice = ice)

      call cnst_get_ind('TT_LW' , ixtt    , abort=.false.)

      tt    = 0._r8
      if (ixtt > 1) then
        if (name_out6 == 'TT_pAM'.or.name_out6 == 'TT_zAM') then
          !
          ! after dme_adjust mixing ratios are all wet
          !
          do k = 1, pver
            do i = 1, ncol
              tt_tmp   = state%q(i,k,ixtt)*state%pdel(i,k)/gravit
              tt   (i) = tt(i)    + tt_tmp
            end do
          end do
        else
          do k = 1, pver
            do i = 1, ncol
              tt_tmp   = state%q(i,k,ixtt)*state%pdeldry(i,k)/gravit
              tt   (i) = tt(i)    + tt_tmp
            end do
          end do
        end if
      end if

      state%te_budgets(1:ncol,teidx,s_ind)=state%te_budgets(1:ncol,teidx,s_ind)+(se(1:ncol)+ke(1:ncol)+po(1:ncol))
      state%te_budgets(1:ncol,seidx,s_ind)=      state%te_budgets(1:ncol,seidx,s_ind)+se(1:ncol)
      state%te_budgets(1:ncol,poidx,s_ind)=      state%te_budgets(1:ncol,poidx,s_ind)+po(1:ncol)
      state%te_budgets(1:ncol,keidx,s_ind)=      state%te_budgets(1:ncol,keidx,s_ind)+ke(1:ncol)
      state%te_budgets(1:ncol,wvidx,s_ind)=      state%te_budgets(1:ncol,wvidx,s_ind)+wv(1:ncol)
      state%te_budgets(1:ncol,wlidx,s_ind)=      state%te_budgets(1:ncol,wlidx,s_ind)+liq(1:ncol)
      state%te_budgets(1:ncol,wiidx,s_ind)=      state%te_budgets(1:ncol,wiidx,s_ind)+ice(1:ncol)
      state%te_budgets(1:ncol,ttidx,s_ind)=      state%te_budgets(1:ncol,ttidx,s_ind)+tt(1:ncol)
      state%budget_cnt(b_ind)=state%budget_cnt(b_ind)+1
      ! Output energy diagnostics

      call outfld(name_out1  ,se+po     ,pcols   ,lchnk   )
      call outfld(name_out2  ,ke        ,pcols   ,lchnk   )
      call outfld(name_out3  ,wv        ,pcols   ,lchnk   )
      call outfld(name_out4  ,liq       ,pcols   ,lchnk   )
      call outfld(name_out5  ,ice       ,pcols   ,lchnk   )
      call outfld(name_out6  ,tt        ,pcols   ,lchnk   )
      call outfld(name_out7  ,se+ke+po  ,pcols   ,lchnk   )

!!$      call outfld(name_out(seidx)  ,se      , pcols   ,lchnk   )
!!$      call outfld(name_out(keidx)  ,ke      , pcols   ,lchnk   )
!!$      call outfld(name_out(wiidx)  ,wv      , pcols   ,lchnk   )
!!$      call outfld(name_out(wlidx)  ,liq     , pcols   ,lchnk   )
!!$      call outfld(name_out(wiidx)  ,ice     , pcols   ,lchnk   )
!!$      call outfld(name_out(ttidx)  ,tt      , pcols   ,lchnk   )
!!$      call outfld(name_out(teidx)  ,te      , pcols   ,lchnk   )

    !
    ! Axial angular momentum diagnostics
    !
    ! Code follows
    !
    ! Lauritzen et al., (2014): Held-Suarez simulations with the Community Atmosphere Model
    ! Spectral Element (CAM-SE) dynamical core: A global axial angularmomentum analysis using Eulerian
    ! and floating Lagrangian vertical coordinates. J. Adv. Model. Earth Syst. 6,129-140,
    ! doi:10.1002/2013MS000268
    !
    ! MR is equation (6) without \Delta A and sum over areas (areas are in units of radians**2)
    ! MO is equation (7) without \Delta A and sum over areas (areas are in units of radians**2)
    !
    name_out1 = 'MR_'   //trim(outfld_name_suffix)
    name_out2 = 'MO_'   //trim(outfld_name_suffix)

!!jt    if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2)) then
      lchnk = state%lchnk
      ncol  = state%ncol

      mr_cnst = rearth**3/gravit
      mo_cnst = omega*rearth**4/gravit

      mr = 0.0_r8
      mo = 0.0_r8
      do k = 1, pver
        do i = 1, ncol
          cos_lat = cos(state%lat(i))
          mr_tmp = mr_cnst*state%u(i,k)*state%pdel(i,k)*cos_lat
          mo_tmp = mo_cnst*state%pdel(i,k)*cos_lat**2

          mr(i) = mr(i) + mr_tmp
          mo(i) = mo(i) + mo_tmp
        end do
      end do
      state%te_budgets(1:ncol,moidx,s_ind)=mo(1:ncol)
      state%te_budgets(1:ncol,mridx,s_ind)=mr(1:ncol)
      call outfld(name_out(mridx)  ,mr, pcols,lchnk   )
      call outfld(name_out(moidx)  ,mo, pcols,lchnk   )

      call outfld(name_out1  ,mr, pcols,lchnk   )
      call outfld(name_out2  ,mo, pcols,lchnk   )
!!jt    end if
  end subroutine calc_te_and_aam_budgets


end module check_energy
