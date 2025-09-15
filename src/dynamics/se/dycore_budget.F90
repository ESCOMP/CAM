module dycore_budget
use shr_kind_mod, only: r8=>shr_kind_r8
implicit none

public :: print_budget
real(r8), parameter :: eps      = 1.0E-7_r8
real(r8), parameter :: eps_mass = 1.0E-12_r8

real(r8), save :: previous_dEdt_adiabatic_dycore    = 0.0_r8
real(r8), save :: previous_dEdt_dry_mass_adjust     = 0.0_r8
real(r8), save :: previous_dEdt_phys_dyn_coupl_err  = 0.0_r8
!=========================================================================================
contains
!=========================================================================================

subroutine print_budget(hstwr)

  use spmd_utils,             only: masterproc
  use cam_abortutils,         only: endrun  
  use cam_logfile,            only: iulog
  use cam_budget,             only: cam_budget_get_global, is_cam_budget, thermo_budget_histfile_num, thermo_budget_history
  use cam_thermo,             only: thermo_budget_vars_descriptor, thermo_budget_num_vars, thermo_budget_vars_massv, &
                                    teidx, seidx, keidx, poidx
  use dimensions_mod,         only: use_cslam
  use control_mod,            only: ftype

  ! arguments
  logical, intent(in) :: hstwr(:)

  ! Local variables
  character(len=*), parameter :: subname = 'dycore_budget:print_budgets:'
  !
  ! physics energy tendencies
  !
  integer  :: idx(4)
  real(r8) :: dEdt_param_physE(4)      ! dE/dt CAM physics using physics E formula (phAP-phBP)
  real(r8) :: dEdt_param_dynE(4)       ! dE/dt CAM physics using dycore E (dyAP-dyBP)

  real(r8) :: dEdt_efix_physE(4)       ! dE/dt energy fixer using physics E formula (phBP-phBF)
  real(r8) :: dEdt_efix_dynE(4)        ! dE/dt energy fixer using dycore E formula (dyBP-dyBF)

  real(r8) :: dEdt_dme_adjust_physE(4) ! dE/dt dry mass adjustment using physics E formula (phAM-phAP)
  real(r8) :: dEdt_dme_adjust_dynE(4)  ! dE/dt dry mass adjustment using dycore E (dyAM-dyAP)

  real(r8) :: dEdt_param_efix_physE(4) ! dE/dt CAM physics + energy fixer using physics E formula (phAP-phBF)
  real(r8) :: dEdt_param_efix_dynE(4)  ! dE/dt CAM physics + energy fixer using dycore E formula (dyAP-dyBF)

  real(r8) :: dEdt_phys_total_dynE(4)  ! dE/dt physics total using dycore E (dyAM-dyBF)
                                       ! physics total = parameterizations + efix + dry-mass adjustment
  !
  ! SE dycore specific energy tendencies
  !
  real(r8) :: dEdt_phys_total_in_dyn(4) ! dEdt of physics total in dynamical core
                                        ! physics total = parameterizations + efix + dry-mass adjustment
  real(r8) :: dEdt_dycore_phys          ! dEdt dycore (estimated in physics)
  !
  ! mass budgets physics
  !
  real(r8) :: dMdt_efix                 ! mass tendency energy fixer
  real(r8) :: dMdt_parameterizations    ! mass tendency physics paramterizations
  real(r8) :: dMdt_dme_adjust           ! mass tendency dry-mass adjustment
  real(r8) :: dMdt_phys_total           ! mass tendency physics total (energy fixer + parameterizations + dry-mass adjustment)
  !
  ! mass budgets dynamics
  !
  real(r8) :: dMdt_floating_dyn         ! mass tendency floating dynamics (dAL-dBL)
  real(r8) :: dMdt_vert_remap           ! mass tendency vertical remapping (dAR-dAD)
  real(r8) :: dMdt_del4_fric_heat       ! mass tendency del4 frictional heating (dAH-dCH)
  real(r8) :: dMdt_del4_tot             ! mass tendency del4 + del4 frictional heating (dAH-dBH)
  real(r8) :: dMdt_residual             ! mass tendency residual (time truncation errors) 
  real(r8) :: dMdt_phys_total_in_dyn    ! mass tendency physics total in dycore
  real(r8) :: dMdt_PDC                  ! mass tendency physics-dynamics coupling
  !
  ! energy budgets dynamics
  !
  real(r8) :: dEdt_floating_dyn         ! dE/dt floating dynamics (dAL-dBL)
  real(r8) :: dEdt_vert_remap           ! dE/dt vertical remapping (dAR-dAD)
  real(r8) :: dEdt_del4                 ! dE/dt del4 (dCH-dBH)
  real(r8) :: dEdt_del4_fric_heat       ! dE/dt del4 frictional heating (dAH-dCH)
  real(r8) :: dEdt_del4_tot             ! dE/dt del4 + del4 fricitional heating (dAH-dBH)
  real(r8) :: dEdt_del2_sponge          ! dE/dt del2 sponge (dAS-dBS)
  real(r8) :: dEdt_del2_del4_tot        ! dE/dt explicit diffusion total
  real(r8) :: dEdt_residual             ! dE/dt residual (dEdt_floating_dyn-dEdt_del2_del4_tot)
  real(r8) :: dEdt_dycore_dyn           ! dE/dt adiabatic dynamical core (calculated in dycore)
  !
  ! physics-dynamics coupling variables
  !
  real(r8) :: E_dBF(4)                  ! E of dynamics state at the end of dycore integration (on dycore deomposition)
  real(r8) :: E_dyBF(4)                 ! E of physics state using dycore E


  real(r8) :: diff, tmp                 ! dummy variables
  integer  :: m_cnst, i
  character(LEN=*), parameter :: fmt  = "(a40,a15,a1,F6.2,a1,F6.2,a1,E10.2,a5)"
  character(LEN=*), parameter :: fmtf = "(a48,F8.4,a6)"
  character(LEN=*), parameter :: fmtm = "(a48,E8.2,a9)"
  character(LEN=15)           :: str(4)
  character(LEN=5)            :: pf     ! pass or fail identifier
  !--------------------------------------------------------------------------------------
  
  if (masterproc .and. thermo_budget_history .and. hstwr(thermo_budget_histfile_num)) then
    idx(1) = teidx !total energy index
    idx(2) = seidx !enthaly index
    idx(3) = keidx !kinetic energy index
    idx(4) = poidx !surface potential energy index
    str(1) = "(total        )"
    str(2) = "(enthalpy     )"
    str(3) = "(kinetic      )"
    str(4) = "(srf potential)"
    do i=1,4
      !
      ! CAM physics energy tendencies
      !
      call cam_budget_get_global('phAP-phBP',idx(i),dEdt_param_physE(i))
      call cam_budget_get_global('phBP-phBF',idx(i),dEdt_efix_physE(i))
      call cam_budget_get_global('phAM-phAP',idx(i),dEdt_dme_adjust_physE(i))
      call cam_budget_get_global('phAP-phBF',idx(i),dEdt_param_efix_physE(i))
      !
      ! CAM physics energy tendencies using dycore energy formula scaling
      ! temperature tendencies for consistency with CAM physics
      !
      call cam_budget_get_global('dyAP-dyBP',idx(i),dEdt_param_dynE(i))
      call cam_budget_get_global('dyBP-dyBF',idx(i),dEdt_efix_dynE(i))
      call cam_budget_get_global('dyAM-dyAP',idx(i),dEdt_dme_adjust_dynE(i))
      call cam_budget_get_global('dyAP-dyBF',idx(i),dEdt_param_efix_dynE(i))
      call cam_budget_get_global('dyAM-dyBF',idx(i),dEdt_phys_total_dynE(i))
      call cam_budget_get_global('dyBF'     ,idx(i),E_dyBF(i))!state beginning physics
      !
      ! CAM physics energy tendencies in dynamical core
      !
      call cam_budget_get_global('dBD-dAF',idx(i),dEdt_phys_total_in_dyn(i))
      call cam_budget_get_global('dBF'    ,idx(i),E_dBF(i))  !state passed to physics
    end do

    call cam_budget_get_global('dAL-dBL',teidx,dEdt_floating_dyn)
    call cam_budget_get_global('dAR-dAD',teidx,dEdt_vert_remap)
    dEdt_dycore_dyn = dEdt_floating_dyn+dEdt_vert_remap

    call cam_budget_get_global('dCH-dBH',teidx,dEdt_del4)
    call cam_budget_get_global('dAH-dCH',teidx,dEdt_del4_fric_heat)
    call cam_budget_get_global('dAH-dBH',teidx,dEdt_del4_tot)
    call cam_budget_get_global('dAS-dBS',teidx,dEdt_del2_sponge)
    dEdt_del2_del4_tot      = dEdt_del4_tot+dEdt_del2_sponge
    dEdt_residual           = dEdt_floating_dyn-dEdt_del2_del4_tot

    write(iulog,*)" "
    write(iulog,*)"======================================================================"
    write(iulog,*)"Total energy diagnostics introduced in Lauritzen and Williamson (2019)"
    write(iulog,*)"(DOI:10.1029/2018MS001549)"
    write(iulog,*)"======================================================================"
    write(iulog,*)" "
    write(iulog,*)"Globally and vertically integrated total energy (E) diagnostics are"
    write(iulog,*)"computed at various points in the physics and dynamics loops to compute"
    write(iulog,*)"energy tendencies (dE/dt) and check for consistency (e.g., is E of"
    write(iulog,*)"state passed to physics computed using dycore state variables the same"
    write(iulog,*)"E of the state in the beginning of physics computed using the physics"
    write(iulog,*)"representation of the state)"
    write(iulog,*)" "
    write(iulog,*)"Energy stages in physics:"
    write(iulog,*)"-------------------------"
    write(iulog,*)" "
    write(iulog,*)"  xxBF: state passed to parameterizations, before energy fixer"
    write(iulog,*)"  xxBP: after energy fixer, before parameterizations"
    write(iulog,*)"  xxAP: after last phys_update in parameterizations and state "
    write(iulog,*)"        saved for energy fixer"
    write(iulog,*)"  xxAM: after dry mass adjustment"
    write(iulog,*)"  history files saved off here"
    write(iulog,*)" "
    write(iulog,*)"where xx='ph','dy' "
    write(iulog,*)" "
    write(iulog,*)"Suffix ph is CAM physics total energy"
    write(iulog,*)"(eq. 111 in Lauritzen et al. 2022; 10.1029/2022MS003117)"
    write(iulog,*)" "
    write(iulog,*)"Suffix dy is dycore energy computed in CAM physics using"
    write(iulog,*)"CAM physics state variables"
    write(iulog,*)" "
    write(iulog,*)" "
    write(iulog,*)"Energy stages in dynamics (specific to the SE dycore)"
    write(iulog,*)"-----------------------------------------------------"
    write(iulog,*)" "
    write(iulog,*)"suffix (d)"
    write(iulog,*)"dED: state from end of previous dynamics (= pBF + time sampling)"
    write(iulog,*)"   loop over vertical remapping and physics dribbling -------- (nsplit) -------"
    write(iulog,*)"            (dribbling and remapping always done together)                    |"
    write(iulog,*)"          dAF: state from previous remapping                                  |"
    write(iulog,*)"          dBD: state after physics dribble, before dynamics                   |"
    write(iulog,*)"          loop over vertical Lagrangian dynamics --------rsplit-------------  |"
    write(iulog,*)"              dynamics here                                                |  |"
    write(iulog,*)"              loop over hyperviscosity ----------hypervis_sub------------  |  |"
    write(iulog,*)"                 dBH   state before hyperviscosity                      |  |  |"
    write(iulog,*)"                 dCH   state after hyperviscosity                       |  |  |"
    write(iulog,*)"                 dAH   state after hyperviscosity momentum heating      |  |  |"
    write(iulog,*)"              end hyperviscosity loop -----------------------------------  |  |"
    write(iulog,*)"              dBS   state before del2 sponge                            |  |  |"
    write(iulog,*)"              dAS   state after del2+mom heating sponge                 |  |  |"
    write(iulog,*)"          end of vertical Lagrangian dynamics loop -------------------------  |"
    write(iulog,*)"      dAD  state after dynamics, before vertical remapping                    |"
    write(iulog,*)"      dAR     state after vertical remapping                                  |"
    write(iulog,*)"   end of remapping loop ------------------------------------------------------"
    write(iulog,*)"dBF  state passed to parameterizations = state after last remapping            "
    write(iulog,*)" "
    write(iulog,*)" "
    write(iulog,*)"FYI: all difference (diff) below are absolute normalized differences"
    write(iulog,*)" "
    write(iulog,*)"Consistency check 0:"
    write(iulog,*)"--------------------"
    write(iulog,*)" "
    write(iulog,*)"For energetic consistency we require that dE/dt [W/m^2] from energy "
    write(iulog,*)"fixer and all parameterizations computed using physics E and"
    write(iulog,*)"dycore in physics E are the same! Checking:"
    write(iulog,*)" "
    write(iulog,*)  "                                                        xx=ph   xx=dy  norm. diff."
    write(iulog,*)  "                                                        -----   -----  -----------"
    do i=1,4
      diff = abs_diff(dEdt_efix_physE(i),dEdt_efix_dynE(i),pf=pf)
      write(iulog,fmt)"dE/dt energy fixer          (xxBP-xxBF) ",str(i)," ",dEdt_efix_physE(i), " ", &
                       dEdt_efix_dynE(i)," ",diff,pf
      diff = abs_diff(dEdt_param_physE(i),dEdt_param_dynE(i),pf=pf)
      write(iulog,fmt)"dE/dt all parameterizations (xxAP-xxBP) ",str(i)," ",dEdt_param_physE(i)," ", &
                      dEdt_param_dynE(i)," ",diff,pf
      write(iulog,*) " "
      if (diff>eps) then
        write(iulog,*)"FAIL"
        call endrun(subname//"dE/dt's in physics inconsistent")
      end if
    end do
    write(iulog,*)" "
    write(iulog,*)" "
    write(iulog,*)"dE/dt from dry-mass adjustment will differ if dynamics and physics use"
    write(iulog,*)"different energy definitions! Checking:"
    write(iulog,*)" "
    write(iulog,*)  "                                                        xx=ph   xx=dy  diff"
    write(iulog,*)  "                                                        -----   -----  ----"
    do i=1,4
      diff = dEdt_dme_adjust_physE(i)-dEdt_dme_adjust_dynE(i)
      write(iulog,fmt)"dE/dt dry mass adjustment   (xxAM-xxAP) ",str(i)," ",dEdt_dme_adjust_physE(i)," ", &
                      dEdt_dme_adjust_dynE(i)," ",diff
    end do
    write(iulog,*)" "
    write(iulog,*)" "
    !
    ! these diagnostics only make sense time-step to time-step
    !
    write(iulog,*)" "
    write(iulog,*)"Some energy budget observations:"
    write(iulog,*)"--------------------------------"
    write(iulog,*)" "
    write(iulog,*)"Note that total energy fixer fixes:"
    write(iulog,*) " "
    write(iulog,*) "-dE/dt energy fixer(t=n) = dE/dt dry mass adjustment (t=n-1) +"
    write(iulog,*) "                      dE/dt adiabatic dycore (t=n-1)         +"
    write(iulog,*) "                      dE/dt physics-dynamics coupling errors (t=n-1)"
    write(iulog,*) " "
    write(iulog,*) "(equation 23 in Lauritzen and Williamson (2019))"
    write(iulog,*) " "

    tmp = previous_dEdt_phys_dyn_coupl_err+previous_dEdt_adiabatic_dycore+previous_dEdt_dry_mass_adjust
    diff = abs_diff(-dEdt_efix_dynE(1),tmp,pf)
    if (.not.use_cslam) then
      write(iulog,*) "Check if that is the case:", pf, diff
      write(iulog,*) " "
      if (abs(diff)>eps) then
        write(iulog,*) "dE/dt energy fixer(t=n)                        = ",dEdt_efix_dynE(1)
        write(iulog,*) "dE/dt dry mass adjustment (t=n-1)              = ",previous_dEdt_dry_mass_adjust
        write(iulog,*) "dE/dt adiabatic dycore (t=n-1)                 = ",previous_dEdt_adiabatic_dycore
        write(iulog,*) "dE/dt physics-dynamics coupling errors (t=n-1) = ",previous_dEdt_phys_dyn_coupl_err
      end if
    else
      previous_dEdt_phys_dyn_coupl_err = dEdt_efix_dynE(1)+previous_dEdt_dry_mass_adjust+previous_dEdt_adiabatic_dycore
      write(iulog,*) "dE/dt energy fixer(t=n)                        = ",dEdt_efix_dynE(1)
      write(iulog,*) "dE/dt dry mass adjustment (t=n-1)              = ",previous_dEdt_dry_mass_adjust
      write(iulog,*) "dE/dt adiabatic dycore (t=n-1)                 = ",previous_dEdt_adiabatic_dycore
      write(iulog,*) "dE/dt physics-dynamics coupling errors (t=n-1) = ",previous_dEdt_phys_dyn_coupl_err
      write(iulog,*) " "
      write(iulog,*) "Note: when running CSLAM the physics-dynamics coupling error is diagnosed"
      write(iulog,*) "      (using equation above) rather than explicitly computed"
      write(iulog,*) " "
      write(iulog,*) " "
      write(iulog,*) "Physics-dynamics coupling errors include: "
      write(iulog,*) " "
      write(iulog,*) " -dE/dt adiabatic dycore is computed on GLL grid;"
      write(iulog,*) " error in mapping to physics grid"
      write(iulog,*) " -dE/dt physics tendencies mapped to GLL grid"
      write(iulog,*) " (tracer tendencies mapped non-conservatively!)"
      write(iulog,*) " -dE/dt dynamics state mapped to GLL grid"
    end if
    write(iulog,*) ""
    if (.not.use_cslam) then
      dEdt_dycore_phys = -dEdt_efix_dynE(1)-previous_dEdt_phys_dyn_coupl_err-previous_dEdt_dry_mass_adjust
      write(iulog,*)               "Hence the dycore E dissipation estimated from energy fixer "
      write(iulog,'(A39,F6.2,A6)') "based on previous time-step values is ",dEdt_dycore_phys," W/M^2"
      write(iulog,*) " "
    end if
    write(iulog,*) " "
    write(iulog,*) "-------------------------------------------------------------------"
    write(iulog,*) " Consistency check 1: state passed to physics same as end dynamics?"
    write(iulog,*) "-------------------------------------------------------------------"
    write(iulog,*) " "
    write(iulog,*) "Is globally integrated total energy of state at the end of dynamics (dBF)"
    write(iulog,*) "and beginning of physics (using dynamics in physics energy; dyBF) the same?"
    write(iulog,*) ""
    if (.not.use_cslam) then
      if (abs(E_dyBF(1))>eps) then
        diff = abs_diff(E_dBF(1),E_dyBF(1))
        if (abs(diff)<eps) then
          write(iulog,'(A23,E8.3)')"yes. (dBF-dyBF)/dyBF =",diff
        else
          write(iulog,*)"Error in physics dynamics coupling!"
          write(iulog,*)" "
          do i=1,4
            write(iulog,*) str(i),":"
            write(iulog,*) "======"
            diff = abs_diff(E_dBF(i),E_dyBF(i),pf=pf)
            write(iulog,*) "diff, E_dBF, E_dyBF ",diff,E_dBF(i),E_dyBF(i)
            write(iulog,*) " "
          end do
          call endrun(subname//"Error in physics dynamics coupling")
        end if
      end if
    else
      write(iulog,*)" "
      write(iulog,*)"Since you are using a separate physics grid, the state in dynamics"
      write(iulog,*)"will not be the same on the physics grid since it is"
      write(iulog,*)"interpolated from the dynamics to the physics grid"
      write(iulog,*)" "
      do i=1,4
        write(iulog,*) str(i),":"
        write(iulog,*) "======"
        diff = abs_diff(E_dBF(i),E_dyBF(i),pf=pf)
        write(iulog,*) "diff, E_dBF, E_dyBF ",diff,E_dBF(i),E_dyBF(i)
        write(iulog,*) " "
      end do
    end if

    write(iulog,*)" "
    write(iulog,*)"-------------------------------------------------------------------------"
    write(iulog,*)" Consistency check 2: total energy increment in dynamics same as physics?"
    write(iulog,*)"-------------------------------------------------------------------------"
    write(iulog,*)" "
    if (.not.use_cslam) then
      previous_dEdt_phys_dyn_coupl_err = dEdt_phys_total_in_dyn(1)-dEdt_phys_total_dynE(1)
      diff = abs_diff(dEdt_phys_total_dynE(1),dEdt_phys_total_in_dyn(1),pf=pf)
      write(iulog,'(A40,E8.2,A7,A5)')" dE/dt physics-dynamics coupling errors       ",diff," W/M^2 ",pf
      if (abs(diff)>eps) then
        !
        ! if errors print details
        !
        if (ftype==1) then
          write(iulog,*) ""
          write(iulog,*) " You are using ftype==1 so physics-dynamics coupling errors should be round-off!"
          write(iulog,*) ""
          write(iulog,*) " Because of failure provide detailed diagnostics below:"
          write(iulog,*) ""
        else
          write(iulog,*) ""
          write(iulog,*) " Since ftype<>1 there are physics dynamics coupling errors"
          write(iulog,*) ""
          write(iulog,*) " Break-down below:"
          write(iulog,*) ""
        end if
        
        do i=1,4
          write(iulog,*) str(i),":"
          write(iulog,*) "======"
          diff = abs_diff(dEdt_phys_total_dynE(i),dEdt_phys_total_in_dyn(i),pf=pf)
          write(iulog,*) "dE/dt physics-dynamics coupling errors (diff) ",diff
          write(iulog,*) "dE/dt physics total in dynamics (dBD-dAF)     ",dEdt_phys_total_in_dyn(i)
          write(iulog,*) "dE/dt physics total in physics  (dyAM-dyBF)   ",dEdt_phys_total_dynE(i)
          write(iulog,*) " "
          write(iulog,*) "      physics total = parameterizations + efix + dry-mass adjustment"
          write(iulog,*) " "
        end do
!   Temporarily disable endrun until energy bias for consistancy check 2 is better understood.
!        if (ftype==1) then
!          call endrun(subname//"Physics-dynamics coupling error. See atm.log")
!        end if
      end if
    else
      write(iulog,'(a47,F6.2,a6)')" dE/dt physics tendency in dynamics (dBD-dAF)   ",dEdt_phys_total_in_dyn(1)," W/M^2"
      write(iulog,'(a47,F6.2,a6)')" dE/dt physics tendency in physics  (dyAM-dyBF) ",dEdt_phys_total_dynE(1)," W/M^2"
      write(iulog,*)" "
      write(iulog,*) " When runnig with a physics grid this consistency check does not make sense"
      write(iulog,*) " since it is computed on the GLL grid whereas we enforce energy conservation"
      write(iulog,*) " on the physics grid. To assess the errors of running dynamics on GLL"
      write(iulog,*) " grid, tracers on CSLAM grid and physics on physics grid we use the energy"
      write(iulog,*) " fixer check from above:"
      write(iulog,*) " "
      write(iulog,*) " dE/dt physics-dynamics coupling errors (t=n-1) =",previous_dEdt_phys_dyn_coupl_err
      write(iulog,*) ""
    end if
    write(iulog,*)" "
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" SE dycore energy tendencies"
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" "
    write(iulog,fmtf)"   dE/dt dycore                                  ",dEdt_dycore_dyn," W/M^2"
    write(iulog,*)" "
    write(iulog,*)"Adiabatic dynamics can be divided into quasi-horizontal and vertical remapping: "
    write(iulog,*)" "
    write(iulog,fmtf)"   dE/dt floating dynamics           (dAD-dBD)   ",dEdt_floating_dyn," W/M^2"
    write(iulog,fmtf)"   dE/dt vertical remapping          (dAR-dAD)   ",dEdt_vert_remap," W/M^2"

    write(iulog,*) " "
    write(iulog,*) "Breakdown of floating dynamics:"
    write(iulog,*) " "
    write(iulog,fmtf)"   dE/dt hypervis del4               (dCH-dBH)   ",dEdt_del4,          " W/M^2"
    write(iulog,fmtf)"   dE/dt hypervis frictional heating (dAH-dCH)   ",dEdt_del4_fric_heat," W/M^2"
    write(iulog,fmtf)"   dE/dt hypervis del4 total         (dAH-dBH)   ",dEdt_del4_tot, " W/M^2"
    write(iulog,fmtf)"   dE/dt hypervis sponge del2        (dAS-dBS)   ",dEdt_del2_sponge,   " W/M^2"
    write(iulog,fmtf)"   dE/dt explicit diffusion total                ",dEdt_del2_del4_tot,    " W/M^2"
    write(iulog,*) " "
    write(iulog,fmtf)"   dE/dt residual (time-truncation errors,...)   ",dEdt_residual,      " W/M^2"
    write(iulog,*)" "
    write(iulog,*)" "
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)"Tracer mass budgets"
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" "
    write(iulog,*)"Below the physics-dynamics coupling error is computed as    "
    write(iulog,*)"dMASS/dt physics tendency in dycore (dBD-dAF) minus"
    write(iulog,*)"dMASS/dt total physics              (pAM-pBF)"
    write(iulog,*)" "
    write(iulog,*)" "
    do m_cnst=1,thermo_budget_num_vars
      if (thermo_budget_vars_massv(m_cnst)) then
        write(iulog,*)thermo_budget_vars_descriptor(m_cnst)
        write(iulog,*)"------------------------------"
        call cam_budget_get_global('phBP-phBF',m_cnst,dMdt_efix)
        call cam_budget_get_global('phAM-phAP',m_cnst,dMdt_dme_adjust)
        call cam_budget_get_global('phAP-phBP',m_cnst,dMdt_parameterizations)
        call cam_budget_get_global('phAM-phBF',m_cnst,dMdt_phys_total)
        !
        ! total energy fixer should not affect mass - checking
        !
        if (abs(dMdt_efix)>eps_mass) then
          write(iulog,*) "dMASS/dt energy fixer        (pBP-pBF)           ",dMdt_efix," Pa/m^2/s"
          write(iulog,*) "ERROR: Mass not conserved in energy fixer. ABORT"      
          call endrun(subname//"Mass not conserved in energy fixer. See atm.log")
        endif
        !
        ! dry-mass adjustmnt should not affect mass - checking
        !
        if (abs(dMdt_dme_adjust)>eps_mass) then
          write(iulog,*)"dMASS/dt dry mass adjustment (pAM-pAP) ",dMdt_dme_adjust," Pa/m^2/s"
          write(iulog,*) "ERROR: Mass not conserved in dry mass adjustment. ABORT"
          call endrun(subname//"Mass not conserved in dry mass adjustment. See atm.log")
        end if
        !
        ! all of the mass-tendency should come from parameterization - checking
        !
        if (abs(dMdt_parameterizations-dMdt_phys_total)>eps_mass) then
          write(iulog,*) "Error: dMASS/dt parameterizations (pAP-pBP) .ne. dMASS/dt physics total (pAM-pBF)"
          write(iulog,*) "dMASS/dt parameterizations   (pAP-pBP) ",dMdt_parameterizations," Pa/m^2/s"
          write(iulog,*) "dMASS/dt physics total       (pAM-pBF) ",dMdt_phys_total," Pa/m^2/s"
          call endrun(subname//"mass change not only due to parameterizations. See atm.log")
        end if
        write(iulog,*)"  "
        !
        ! detailed mass budget in dynamical core
        !
        if (is_cam_budget('dAD').and.is_cam_budget('dBD').and.is_cam_budget('dAR').and.is_cam_budget('dCH')) then
          call cam_budget_get_global('dAL-dBL',m_cnst,dMdt_floating_dyn)
          call cam_budget_get_global('dAR-dAD',m_cnst,dMdt_vert_remap)
          tmp  = dMdt_floating_dyn+dMdt_vert_remap
          diff = abs_diff(tmp,0.0_r8,pf=pf)
          write(iulog,fmtm)"   dMASS/dt total adiabatic dynamics             ",diff,pf
          !
          ! check for mass-conservation in the adiabatic dynamical core -
          ! if not conserved provide detailed break-down
          !
          if (abs(diff)>eps_mass) then
            write(iulog,*) "Error: mass non-conservation in dynamical core"
            write(iulog,*) "(detailed budget below)"
            write(iulog,*) " "
            write(iulog,*)"dMASS/dt 2D dynamics            (dAL-dBL) ",dMdt_floating_dyn," Pa/m^2/s"
            write(iulog,*)"dE/dt vertical remapping        (dAR-dAD) ",dMdt_vert_remap
            write(iulog,*)" "
            write(iulog,*)"Breakdown of 2D dynamics:"
            write(iulog,*)" "
            call cam_budget_get_global('dAH-dCH',m_cnst,dMdt_del4_fric_heat)
            call cam_budget_get_global('dAH-dBH',m_cnst,dMdt_del4_tot)
            write(iulog,*)"dMASS/dt hypervis               (dAH-dBH) ",dMdt_del4_tot," Pa/m^2/s"
            write(iulog,*)"dMASS/dt frictional heating     (dAH-dCH) ",dMdt_del4_fric_heat," Pa/m^2/s"
            dMdt_residual = dMdt_floating_dyn-dMdt_del4_tot
            write(iulog,*)"dMASS/dt residual (time truncation errors)",dMdt_residual," Pa/m^2/s"
          end if
        end if
        if (is_cam_budget('dBD').and.is_cam_budget('dAF')) then
          !
          ! check if mass change in physics is the same as dynamical core
          !
          call cam_budget_get_global('dBD-dAF',m_cnst,dMdt_phys_total_in_dyn)
          dMdt_PDC = dMdt_phys_total-dMdt_phys_total_in_dyn
          write(iulog,fmtm)"   Mass physics-dynamics coupling error          ",dMdt_PDC," Pa/m^2/s"
          write(iulog,*)" "
          if (abs(dMdt_PDC)>eps_mass) then
            write(iulog,fmtm)"   dMASS/dt physics tendency in dycore (dBD-dAF) ",dMdt_phys_total_in_dyn," Pa/m^2/s"
            write(iulog,fmtm)"   dMASS/dt total physics                        ",dMdt_phys_total," Pa/m^2/s"
          end if
        end if
      end if
    end do
    !
    ! save adiabatic dycore dE/dt and dry-mass adjustment to avoid samping error
    !
    previous_dEdt_adiabatic_dycore = dEdt_dycore_dyn
    previous_dEdt_dry_mass_adjust  = dEdt_dme_adjust_dynE(1)
  end if
end subroutine print_budget
!=========================================================================================
function abs_diff(a,b,pf)
  real(r8),                   intent(in) :: a,b
  character(LEN=5), optional, intent(out):: pf
  real(r8)                               :: abs_diff
  if (abs(b)>eps) then
    abs_diff = abs((b-a)/b)
  else
    abs_diff = abs(b-a)
  end if
  If (present(pf)) then
    if (abs_diff>eps) then
      pf = ' FAIL'
    else
      pf = ' PASS'
    end if
  end if
end function abs_diff
end module dycore_budget
