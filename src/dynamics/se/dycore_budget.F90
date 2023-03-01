module dycore_budget
use shr_kind_mod, only: r8=>shr_kind_r8
implicit none

public :: print_budget
real(r8), parameter :: eps      = 1.0E-9_r8
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
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use budgets,                only: budget_get_global, is_budget, thermo_budget_histfile_num, thermo_budget_history
  use cam_thermo,             only: thermo_budget_vars_descriptor, thermo_budget_num_vars, thermo_budget_vars_massv, &
                                    teidx, seidx, keidx, poidx
  use dimensions_mod,         only: ntrac
  use control_mod,            only: ftype
  use cam_thermo,             only: teidx, seidx, keidx, poidx
  use cam_thermo,             only: thermo_budget_vars_descriptor, thermo_budget_num_vars, thermo_budget_vars_massv

  ! arguments
  logical, intent(in) :: hstwr(:)

  ! Local variables
  character(len=*), parameter :: subname = 'check_energy:print_budgets'

  integer,  dimension(4) :: idx
  real(r8), dimension(4) :: ph_param,ph_EFIX,ph_DMEA,ph_PARAM_AND_EFIX,ph_phys_total
  real(r8), dimension(4) :: dy_param,dy_EFIX,dy_DMEA,dy_param_and_efix,dy_phys_total
  real(r8), dimension(4) :: se_phys_total
  real(r8)               :: dycore, err, param, pefix, &
                            pdmea, phys_total, dyn_total, dyn_phys_total, &
                            rate_of_change_2D_dyn, rate_of_change_vertical_remapping, &
                            diffusion_del4, diffusion_fric, diffusion_del4_tot, diffusion_sponge, &
                            diffusion_total, twoDresidual, &
                            rate_of_change_heating_term_put_back_in, rate_of_change_hvis_sponge, &
                            dADIA, &
                            mass_change__2D_dyn,mass_change__vertical_remapping, &
                            mass_change__heating_term_put_back_in,mass_change__hypervis_total, &
                            error, mass_change__physics, dbd, daf, dar, dad, val

  real(r8) :: E_dBF(4), E_phBF, diff, tmp
  real(r8) :: E_dyBF(4)
  integer  :: m_cnst, i
  character(LEN=*), parameter :: fmt  = "(a40,a15,a1,F6.2,a1,F6.2,a1,E10.2,a5)"
  character(LEN=*), parameter :: fmt2 = "(a40,F6.2,a3)"
  character(LEN=15)           :: str(4)
  character(LEN=5)  :: pf! pass or fail identifier
  !--------------------------------------------------------------------------------------

  if (masterproc .and. thermo_budget_history .and. hstwr(thermo_budget_histfile_num)) then
    idx(1) = teidx !total energy index
    idx(2) = seidx !enthaly index
    idx(3) = keidx !kinetic energy index
    idx(4) = poidx !surface potential energy index
    str(1) = "(total)       )"
    str(2) = "(enthalpy     )"
    str(3) = "(kinetic      )"
    str(4) = "(srf potential)"
    do i=1,4
      !
      ! CAM physics energy tendencies
      !
      call budget_get_global('phAP-phBP',idx(i),ph_param(i))
      call budget_get_global('phBP-phBF',idx(i),ph_EFIX(i))
      call budget_get_global('phAM-phAP',idx(i),ph_dmea(i))
      call budget_get_global('phAP-phBF',idx(i),ph_param_and_efix(i))
      call budget_get_global('phAM-phBF',idx(i),ph_phys_total(i))
      !
      ! CAM physics energy tendencies using dycore energy formula scaling
      ! temperature tendencies for consistency with CAM physics
      !
      call budget_get_global('dyAP-dyBP',idx(i),dy_param(i))
      call budget_get_global('dyBP-dyBF',idx(i),dy_EFIX(i))
      call budget_get_global('dyAM-dyAP',idx(i),dy_dmea(i))
      call budget_get_global('dyAP-dyBF',idx(i),dy_param_and_efix(i))
      call budget_get_global('dyAM-dyBF',idx(i),dy_phys_total(i))
      call budget_get_global('dyBF'     ,idx(i),E_dyBF(i))!state beginning physics
      !
      ! CAM physics energy tendencies in dynamical core
      !
      call budget_get_global('dBD-dAF',idx(i),se_phys_total(i))
      call budget_get_global('dBF'    ,idx(i),E_dBF(i))  !state passed to physics
    end do

    call budget_get_global('dBF-dED',teidx,dyn_total)
    call budget_get_global('dAD-dBD',teidx,rate_of_change_2D_dyn)
    call budget_get_global('dAR-dAD',teidx,rate_of_change_vertical_remapping)
    dADIA = rate_of_change_2D_dyn+rate_of_change_vertical_remapping

    call budget_get_global('dCH-dBH',teidx,diffusion_del4)
    call budget_get_global('dAH-dCH',teidx,diffusion_fric)
    call budget_get_global('dAH-dBH',teidx,diffusion_del4_tot)
    call budget_get_global('dAS-dBS',teidx,diffusion_sponge)
    diffusion_total      = diffusion_del4_tot+diffusion_sponge

    rate_of_change_heating_term_put_back_in = diffusion_fric
    rate_of_change_hvis_sponge = diffusion_sponge

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
    write(iulog,*)"Energy stages in dynamics"
    write(iulog,*)"-------------------------"
    write(iulog,*)" "
    write(iulog,*)"suffix (dynamics)"
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
      diff = abs_diff(ph_EFIX(i),dy_EFIX(i),pf=pf)
      write(iulog,fmt)"dE/dt energy fixer          (xxBP-xxBF) ",str(i)," ",ph_EFIX(i), " ",dy_EFIX(i)," ",diff,pf
      diff = abs_diff(ph_param(i),dy_param(i),pf=pf)
      write(iulog,fmt)"dE/dt all parameterizations (xxAP-xxBP) ",str(i)," ",ph_param(i)," ",dy_param(i)," ",diff,pf
      write(iulog,*) " "
    end do
    if (diff>eps) then
      write(iulog,*)"FAIL"
      call endrun(subname//"dE/dts in physics inconsistent")
    end if
    write(iulog,*)" "
    write(iulog,*)" "
    write(iulog,*)"dE/dt from dry-mass adjustment will differ if dynamics and physics use"
    write(iulog,*)"different energy definitions! Checking:"
    write(iulog,*)" "
    write(iulog,*)  "                                                        xx=ph   xx=dy  diff"
    write(iulog,*)  "                                                        -----   -----  ----"
    do i=1,4
      diff = ph_dmea(i)-dy_dmea(i)
      write(iulog,fmt)"dE/dt dry mass adjustment   (xxAM-xxAP) ",str(i)," ",ph_dmea(i)," ",dy_dmea(i)," ",diff
      write(iulog,*) ""
      write(iulog,*) str(i),":"
      write(iulog,*) "======"
      write(iulog,*)"dE/dt dry mass adjustment   (phAM-phAP)"," ",ph_dmea(i)
      write(iulog,*)"dE/dt dry mass adjustment   (dyAM-dyAP)"," ",dy_dmea(i)
      write(iulog,*) " "
      write(iulog,*) " "
    end do
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
    diff = abs_diff(-dy_EFIX(1),tmp,pf)
    if (ntrac==0) then
      write(iulog,*) "Check if that is the case:", pf, diff
      write(iulog,*) " "
      if (abs(diff)>eps) then
        write(iulog,*) "dE/dt energy fixer(t=n)                        = ",dy_EFIX(1)
        write(iulog,*) "dE/dt dry mass adjustment (t=n-1)              = ",previous_dEdt_dry_mass_adjust
        write(iulog,*) "dE/dt adiabatic dycore (t=n-1)                 = ",previous_dEdt_adiabatic_dycore
        write(iulog,*) "dE/dt physics-dynamics coupling errors (t=n-1) = ",previous_dEdt_phys_dyn_coupl_err
        !      call endrun(subname//"Error in energy fixer budget")
      end if
    else
      previous_dEdt_phys_dyn_coupl_err = dy_EFIX(1)+previous_dEdt_dry_mass_adjust+previous_dEdt_adiabatic_dycore
      write(iulog,*) "dE/dt energy fixer(t=n)                        = ",dy_EFIX(1)
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
    if (ntrac==0) then
      dycore = -dy_EFIX(1)-previous_dEdt_phys_dyn_coupl_err-previous_dEdt_dry_mass_adjust
      write(iulog,*) "Hence the dycore E dissipation estimated from energy fixer "
      write(iulog,*) "based on previous time-step values is ",dycore," W/M^2"
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
    if (ntrac==0) then
      if (abs(E_dyBF(1))>eps) then
        diff = abs_diff(E_dBF(1),E_dyBF(1))
        if (abs(diff)<eps) then
          write(iulog,*)"yes. (dBF-dyBF)/dyBF =",diff
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
      write(iulog,*)"will not be the same on the physics grid"
      write(iulog,*)"interpolated from the physics to the dynamics "
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
    if (ntrac>0) then
      write(iulog,'(a46,F6.2,a6)')"dE/dt physics tendency in dynamics (dBD-dAF)   ",se_phys_total(1)," W/M^2"
      write(iulog,'(a46,F6.2,a6)')"dE/dt physics tendency in physics  (dyAM-dyBF) ",dy_phys_total(1)," W/M^2"
      write(iulog,*)" "
      write(iulog,*) " When runnig with a physics grid this consistency check does not make sense"
      write(iulog,*) " since it is computed on the GLL grid whereas we enforce energy conservation"
      write(iulog,*) " on the physics grid. To assess the errors of running dynamics on GLL"
      write(iulog,*) " grid, tracers on CSLAM grid and physics on physics grid we use the energy"
      write(iulog,*) " fixer check from above:"
      write(iulog,*) " "
      write(iulog,*) " dE/dt physics-dynamics coupling errors (t=n-1) =",previous_dEdt_phys_dyn_coupl_err
      write(iulog,*) ""
    else
      previous_dEdt_phys_dyn_coupl_err = se_phys_total(1)-dy_phys_total(1)
      diff = abs_diff(dy_phys_total(1),se_phys_total(1),pf=pf)
      write(iulog,*)"dE/dt physics-dynamics coupling errors       ",diff," W/M^2 "
      write(iulog,*) pf
      if (abs(diff)>eps) then
        !
        ! if errors print details
        !
        if (ftype==1) then
          write(iulog,*) ""
          write(iulog,*) "You are using ftype==1 so physics-dynamics coupling errors should be round-off!"
          write(iulog,*) ""
          write(iulog,*) "Because of failure provide detailed diagnostics below:"
          write(iulog,*) ""
        else
          write(iulog,*) ""
          write(iulog,*) "Since ftype<>1 there are physics dynamics coupling errors"
          write(iulog,*) ""
          write(iulog,*) "Break-down below:"
          write(iulog,*) ""
        end if
!        else
!          write(iulog,*)" "
!          write(iulog,*)"Since you are using a separate physics grid, the physics tendencies"
!          write(iulog,*)"in the dynamical core will not match due to the tendencies being"
!          write(iulog,*)"interpolated from the physics to the dynamics grid:"
!          write(iulog,*)" "
        do i=1,4
          write(iulog,*) str(i),":"
          write(iulog,*) "======"
          diff = abs_diff(dy_phys_total(i),se_phys_total(i),pf=pf)
          write(iulog,*) "dE/dt physics-dynamics coupling errors (diff) ",diff
          write(iulog,*) "dE/dt physics tendency in dynamics (dBD-dAF)  ",se_phys_total(i)
          write(iulog,*) "dE/dt physics tendency in physics  (pAM-pBF)  ",dy_phys_total(i)
          write(iulog,*) " "
        end do
      end if
    end if
    write(iulog,*)" "
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" SE dycore energy tendencies"
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" "
    !     write(iulog,*)"dE/dt dyn total (dycore+phys tendency   (dBF-dED) ",dyn_total," W/M^2"
    write(iulog,'(a46,F6.2,a6)')"dE/dt adiabatic dynamics                     ",dADIA," W/M^2"
    write(iulog,*)" "
    write(iulog,*)"Adiabatic dynamics can be divided into quasi-horizontal and vertical remapping: "
    write(iulog,*)" "
    write(iulog,'(a40,F6.2,a6)') "dE/dt 2D dynamics           (dAD-dBD)  ",rate_of_change_2D_dyn," W/M^2"
    write(iulog,'(a40,F6.2,a6)') "dE/dt vertical remapping    (dAR-dAD)  ",rate_of_change_vertical_remapping," W/M^2"

    write(iulog,*) " "
    write(iulog,*) "Breakdown of 2D dynamics:"
    write(iulog,*) " "
    write(iulog,'(a45,F6.2,a6)')"   dE/dt hypervis del4               (dCH-dBH) ",diffusion_del4," W/M^2"
    write(iulog,'(a45,F6.2,a6)')"   dE/dt hypervis frictional heating (dAH-dCH) ",diffusion_fric," W/M^2"
    write(iulog,'(a45,F6.2,a6)')"   dE/dt hypervis del4 total         (dAH-dBH) ",diffusion_del4_tot," W/M^2"
    write(iulog,'(a45,F6.2,a6)')"   dE/dt hypervis sponge total       (dAS-dBS) ",diffusion_sponge," W/M^2"
    write(iulog,'(a45,F6.2,a6)')"   dE/dt explicit diffusion total              ",diffusion_total," W/M^2"
    twoDresidual = rate_of_change_2D_dyn-diffusion_total
    write(iulog,'(a45,F6.2,a6)')"   dE/dt residual (time-truncation errors)     ",twoDresidual," W/M^2"
    write(iulog,*)" "
    write(iulog,*)" "
#ifdef xxx
    write(iulog,*)" "
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" CAM physics energy tendencies (using pressure coordinate)"
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" "
    write(iulog,'(a40,F6.2,a6)')"dE/dt energy fixer          (phBP-phBF) ",ph_EFIX," W/M^2"
    write(iulog,'(a40,F6.2,a6)')"dE/dt all parameterizations (phAP-phBP) ",ph_param," W/M^2"
    write(iulog,'(a40,F6.2,a6)')"dE/dt dry mass adjustment   (phAM-phAP) ",ph_DMEA," W/M^2"
    write(iulog,'(a40,F6.2,a6)')"dE/dt physics total         (phAM-phBF) ",ph_phys_total," W/M^2"
    write(iulog,*)" "
    write(iulog,*) " "
    write(iulog,*) "-dE/dt energy fixer = dE/dt dry mass adjustment              +"
    write(iulog,*) "                      dE/dt dycore                           +"
    write(iulog,*) "                      dE/dt physics-dynamics coupling errors +"
    write(iulog,*) "                      dE/dt energy formula differences       "
    write(iulog,*) " "
    write(iulog,*) "(equation 23 in Lauritzen and Williamson (2019))"
    write(iulog,*) " "
    dycore = -ph_EFIX-ph_DMEA
    dycore = -ph_EFIX-previous_dEdt_dry_mass_adjust
    write(iulog,*) ""
    write(iulog,*) "Dycore TE dissipation estimated from physics in pressure coordinate:"
    write(iulog,*) "(note: to avoid sampling error we need dE/dt from previous time-step)"
    write(iulog,*) ""
    write(iulog,*) "dE/dt adiabatic dycore estimated from physics (t=n-1) = "
    write(iulog,'(a58,F6.2,a6)') "-dE/dt energy fixer(t=n)-dE/dt dry-mass adjust(t=n-1) = ",dycore," W/M^2"
    write(iulog,*) ""
    write(iulog,'(a58,F6.2,a6)') "dE/dt adiabatic dycore computed in dycore (t=n-1)     = ",&
         previous_dEdt_adiabatic_dycore," W/M^2"
    write(iulog,'(a58,F6.2,a6)') "dE/dt dry-mass adjust  (t=n-1)                        = ",&
         previous_dEdt_dry_mass_adjust," W/M^2"
    write(iulog,*) ""
    if (abs(previous_dEdt_adiabatic_dycore)>eps) then
      diff = abs((dycore-previous_dEdt_adiabatic_dycore)/previous_dEdt_adiabatic_dycore)
      if (diff>eps) then
        write(iulog,*) "energy budget not closed: previous_dEdt_adiabatic_dycore <> dycore"
        write(iulog,*) "normalized difference is:",diff
        !        call endrun(subname//"physics energy budget consistency error 2")
      end if
    end if
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" Physics dynamics coupling errors"
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" "
    write(iulog,'(a46,F6.2,a6)')"dE/dt physics tendency in dynamics (dBD-dAF) ",se_phys_total_te," W/M^2"
    write(iulog,'(a46,F6.2,a6)')"dE/dt physics tendency in physics  (pAM-pBF) ",ph_phys_total," W/M^2"
    write(iulog,*)" "
    write(iulog,'(a46,F6.2,a6)')"dE/dt physics-dynamics coupling errors       ",ph_phys_total-se_phys_total_te," W/M^2"

    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" Consistency checks"
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" "
    !
    ! consistency check
    !
    if (abs(ph_param+ph_EFIX+ph_DMEA-ph_phys_total)>eps) then
      write(iulog,*) "Physics energy budget not adding up:"
      write(iulog,*) "(phBP-pBF)+(phAP-pBP)+(pAM-pAP) does not add up to (pAM-pBF)",\
      abs(ph_param+ph_EFIX+ph_DMEA-ph_phys_total)
      call endrun(subname//"physics energy budget consistency error")
    endif
    write(iulog,*) ""
    write(iulog,*) "Is globally integrated total energy of state at the end of dynamics (dBF)"
    write(iulog,*) "and beginning of physics (phBF) the same?"
    write(iulog,*) ""
    call budget_get_global('dBF' ,teidx,E_dBF(1))  !state passed to physics
    call budget_get_global('phBF',teidx,E_phBF)!state beginning physics
    !     if (abs(E_phBF)>eps) then
    diff = abs_diff(E_dBF(1),E_phBF)
    if (abs(diff)<eps) then
      write(iulog,*)"yes. (dBF-phBF)/phBF =",diff
      write(iulog,*)"E_dBF=",E_dBF(1),"; E_phBF=",E_phBF
    else
      write(iulog,*) "no. (dBF-phBF)/phBF =",diff,E_dBF(1),E_phBF
      write(iulog,*) "To run energy consistent version of SE use namelist"
      write(iulog,*) ""
      write(iulog,*) "se_ftype     = 1           !no dribbling of physics tendencies"
      write(iulog,*) "se_lcp_moist =  .false.    !no variable latent heats"
      write(iulog,*) "water_species_in_air = 'Q' !only water vapor energetically active"
      write(iulog,*) ""
    end if
    write(iulog,*) ""
    write(iulog,*) "Is globally integrated total energy of state at the end of dynamics (dBF)"
    write(iulog,*) "and beginning of physics dynamics energy (dyBF) the same?"
    write(iulog,*) ""
    diff = abs_diff(E_dBF(1),E_dyBF(1))
    if (abs(diff)<eps) then
      write(iulog,*)"yes. (dBF-dyBF)/dyBF =",diff
      write(iulog,*)"E_dBF=",E_dBF(1),"; E_dyBF=",E_dyBF(1)
    else
      write(iulog,*) "no. (dBF-dyBF)/dyBF =",diff,E_dBF(1),E_dyBF(1)
    end if
    !     end if
#endif
    do m_cnst=1,thermo_budget_num_vars
      if (thermo_budget_vars_massv(m_cnst)) then
        write(iulog,*)"------------------------------------------------------------"
        write(iulog,*)thermo_budget_vars_descriptor(m_cnst)//" budget"
        write(iulog,*)"------------------------------------------------------------"
        call budget_get_global('phBP-phBF',m_cnst,pEFIX)
        call budget_get_global('phAM-phAP',m_cnst,pDMEA)
        call budget_get_global('phAP-phBP',m_cnst,param)
        call budget_get_global('phAM-phBF',m_cnst,phys_total)
        if (abs(pEFIX)>eps_mass) then
          write(iulog,*) "dMASS/dt energy fixer        (pBP-pBF) ",pEFIX," Pa"
          write(iulog,*) "ERROR: Mass not conserved in energy fixer. ABORT"      
          call endrun(subname//"Mass not conserved in energy fixer. See atm.log")
        endif
        if (abs(pDMEA)>eps_mass) then
          write(iulog,*)"dMASS/dt dry mass adjustment (pAM-pAP) ",pDMEA," Pa"
          write(iulog,*) "ERROR: Mass not conserved in dry mass adjustment. ABORT"
          call endrun(subname//"Mass not conserved in dry mass adjustment. See atm.log")
        end if
        if (abs(param-phys_total)>eps_mass) then
          write(iulog,*) "Error: dMASS/dt parameterizations (pAP-pBP) .ne. dMASS/dt physics total (pAM-pBF)"
          write(iulog,*) "dMASS/dt parameterizations   (pAP-pBP) ",param," Pa"
          write(iulog,*) "dMASS/dt physics total       (pAM-pBF) ",phys_total," Pa"
          call endrun(subname//"mass change not only due to parameterizations. See atm.log")
        end if
        write(iulog,*)"dMASS/dt parameterizations   (pAP-pBP) ",param," Pa"
        write(iulog,*)"dMASS/dt physics total       (pAM-pBF) ",phys_total," Pa"
        write(iulog,*)"  "
        !
        ! detailed mass budget in dynamical core
        !
        if (is_budget('dAD').and.is_budget('dBD').and.is_budget('dAR').and.is_budget('dCH')) then
          call budget_get_global('dAD-dBD',m_cnst,mass_change__2D_dyn)
          call budget_get_global('dAR-dAD',m_cnst,mass_change__vertical_remapping)
          diff = mass_change__2D_dyn+mass_change__vertical_remapping
          write(iulog,*)"dMASS/dt total adiabatic dynamics       ",diff," Pa"
          if (abs(diff)>eps_mass) then
            write(iulog,*) "Error: mass non-conservation in dynamical core"
            write(iulog,*) "(detailed budget below)"
            write(iulog,*) " "
            write(iulog,*)"dMASS/dt 2D dynamics            (dAD-dBD) ",mass_change__2D_dyn," Pa"
            if (is_budget('dAR').and.is_budget('dAD')) then
              call budget_get_global('dAR',m_cnst,dar)
              call budget_get_global('dAD',m_cnst,dad)
              call budget_get_global('dAR-dAD',m_cnst,mass_change__vertical_remapping)
              write(iulog,*)"dE/dt vertical remapping        (dAR-dAD) ",mass_change__vertical_remapping
            end if
            write(iulog,*)" "
            write(iulog,*)"Breakdown of 2D dynamics:"
            write(iulog,*)" "
            call budget_get_global('dAH-dCH',m_cnst,mass_change__heating_term_put_back_in)
            call budget_get_global('dAH-dBH',m_cnst,mass_change__hypervis_total)
            write(iulog,*)"dMASS/dt hypervis               (dAH-dBH) ",mass_change__hypervis_total," Pa"
            write(iulog,*)"dMASS/dt frictional heating     (dAH-dCH) ",mass_change__heating_term_put_back_in," Pa"
            error = mass_change__2D_dyn-mass_change__hypervis_total
            write(iulog,*)"dMASS/dt residual (time truncation errors)",error," Pa"
          end if
        end if
        write(iulog,*)" "
        if (is_budget('dBD').and.is_budget('dAF')) then
          call budget_get_global('dBD',m_cnst,dbd)
          call budget_get_global('dAF',m_cnst,daf)
          call budget_get_global('dBD-dAF',m_cnst,mass_change__physics)
          write(iulog,*)"dMASS/dt physics tendency in dynamics (dBD-dAF) ",mass_change__physics," Pa"
          val = phys_total-mass_change__physics
          write(iulog,*) " "
          write(iulog,*) "Mass physics dynamics coupling error:",val
        end if
        write(iulog,*)""
      end if
    end do
    !
    ! save adiabatic dycore dE/dt and dry-mass adjustment to avoid samping error
    !
    previous_dEdt_adiabatic_dycore = dADIA
    previous_dEdt_dry_mass_adjust  = dy_DMEA(1)
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
