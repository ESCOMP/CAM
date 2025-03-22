module dycore_budget
use shr_kind_mod, only: r8=>shr_kind_r8  
implicit none

public :: print_budget
real(r8), parameter :: eps      = 1.0E-9_r8
real(r8), parameter :: eps_mass = 1.0E-12_r8
real(r8), save :: previous_dEdt_dry_mass_adjust          = 0.0_r8
real(r8), save :: previous_dEdt_phys_dyn_coupl_err_Agrid = 0.0_r8
!=========================================================================================
contains
!=========================================================================================

subroutine print_budget(hstwr)

  use cam_budget,     only: cam_budget_get_global, thermo_budget_histfile_num, thermo_budget_history
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use cam_thermo,     only: thermo_budget_vars_descriptor, thermo_budget_num_vars, thermo_budget_vars_massv
  use cam_thermo,     only: teidx, seidx, keidx, poidx

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
  ! dycore specific energy tendencies
  !
  real(r8) :: dEdt_phys_total_in_dyn(4) ! dEdt of physics total in dynamical core
                                        ! physics total = parameterizations + efix + dry-mass adjustment
  real(r8) :: dEdt_param_efix_in_dyn(4) ! dEdt CAM physics + energy fixer in dynamical core
  real(r8) :: dEdt_dme_adjust_in_dyn(4) ! dEdt of dme adjust in dynamical core
  real(r8) :: dEdt_dycore_and_pdc_estimated_from_efix ! dEdt dycore and PDC errors (estimated in physics)
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
  real(r8) :: dMdt_phys_total_in_dyn    ! mass tendency physics total in dycore
  real(r8) :: dMdt_PDC                  ! mass tendency physics-dynamics coupling
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
    str(2) = "(internal     )"
    str(3) = "(kinetic      )"
    str(4) = "(potential    )"
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
      call cam_budget_get_global('dAP-dBF',teidx,dEdt_param_efix_in_dyn(i))
      call cam_budget_get_global('dAM-dAP',teidx,dEdt_dme_adjust_in_dyn(i))
      call cam_budget_get_global('dAM-dBF',teidx,dEdt_param_efix_in_dyn(i))

      call cam_budget_get_global('dAM-dBF',idx(i),dEdt_phys_total_in_dyn(i))
      call cam_budget_get_global('dBF'    ,idx(i),E_dBF(i))  !state passed to physics
    end do
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
    write(iulog,*)"Energy stages in dynamics"
    write(iulog,*)"-------------------------"
    write(iulog,*)" "
    write(iulog,*)" dBF: dynamics state before physics (d_p_coupling)"
    write(iulog,*)" dAP: dynamics state with T,u,V increment but not incl water changes"
    write(iulog,*)" dAM: dynamics state with full physics increment (incl. water)"
    write(iulog,*)" "
    write(iulog,*)"Note that these energies are computed using the dynamical core"
    write(iulog,*)"state variables which may be different from the physics prognostic"
    write(iulog,*)"variables."
    write(iulog,*)" "
    write(iulog,*)" "
    write(iulog,*)"FYI : norm. diff = absolute normalized difference"
    write(iulog,*)"FYI2: diff       = difference (not normalized)"
    write(iulog,*)" "
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
      write(iulog,fmt)"dE/dt energy fixer          (xxBP-xxBF) ",str(i)," ",dEdt_efix_physE(i), " ",dEdt_efix_dynE(i)," ",diff,pf
      diff = abs_diff(dEdt_param_physE(i),dEdt_param_dynE(i),pf=pf)
      write(iulog,fmt)"dE/dt all parameterizations (xxAP-xxBP) ",str(i)," ",dEdt_param_physE(i)," ",dEdt_param_dynE(i)," ",diff,pf
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
      write(iulog,fmt)"dE/dt dry mass adjustment   (xxAM-xxAP) ",str(i)," ",dEdt_dme_adjust_physE(i), &
           " ",dEdt_dme_adjust_dynE(i)," ",diff
    end do
    write(iulog,*)" "
    write(iulog,*)"Compare to dry mass adjustment in dynamics (xx=d,dy):"
    write(iulog,*)  "                                                        xx=d    xx=dy  norm. diff"
    write(iulog,*)  "                                                        -----   -----  ----------"
    do i=1,4
      diff = abs_diff(dEdt_dme_adjust_in_dyn(i),dEdt_dme_adjust_dynE(i),pf=pf)
      write(iulog,fmt)"dE/dt dry mass adjustment   (xxAM-xxAP) ",str(i)," ",dEdt_dme_adjust_in_dyn(i),&
           " ",dEdt_dme_adjust_dynE(i)," ",diff,pf
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
    write(iulog,*)" Note that total energy fixer fixes:"
    write(iulog,*)" "
    write(iulog,*)" -dE/dt energy fixer(t=n) = dE/dt dry mass adjustment (t=n-1) +"
    write(iulog,*)"                            dE/dt adiabatic dycore (t=n-1)    +"
    write(iulog,*)"                            dE/dt physics-dynamics coupling errors (t=n-1)"
    write(iulog,*)" "
    write(iulog,*)" (equation 23 in Lauritzen and Williamson (2019))"
    write(iulog,*)" "
    write(iulog,*)" Technically this equation is only valid with instantaneous time-step to"
    write(iulog,*)" time-step output"
    write(iulog,*) " "
    write(iulog,*) " dE/dt energy fixer(t=n)           = ",dEdt_efix_dynE(1)
    write(iulog,*) " dE/dt dry mass adjustment (t=n-1) = ",previous_dEdt_dry_mass_adjust
    write(iulog,*) " dE/dt adiabatic dycore (t=n-1)    = unknown"
    write(iulog,*) " dE/dt PDC errors (A-grid) (t=n-1) = ",previous_dEdt_phys_dyn_coupl_err_Agrid
    write(iulog,*) " dE/dt PDC errors (other ) (t=n-1) = unknown"

    dEdt_dycore_and_pdc_estimated_from_efix = -dEdt_efix_dynE(1) - &
                                              previous_dEdt_phys_dyn_coupl_err_Agrid - &
                                              previous_dEdt_dry_mass_adjust
    write(iulog,*) " "
    write(iulog,*)               "Hence the dycore E dissipation and physics-dynamics coupling errors"
    write(iulog,*)               "associated with mapping wind tendencies to C-grid and dribbling    "
    write(iulog,*)               "tendencies in the dycore (PDC other), estimated from energy fixer "
    write(iulog,'(A39,F6.2,A6)') "based on previous time-step values is ",dEdt_dycore_and_pdc_estimated_from_efix," W/M^2"
    write(iulog,*) " "
    write(iulog,*) " "
    write(iulog,*) "-------------------------------------------------------------------"
    write(iulog,*) " Consistency check 1: state passed to physics same as end dynamics?"
    write(iulog,*) "-------------------------------------------------------------------"
    write(iulog,*) " "
    write(iulog,*) "Is globally integrated total energy of state at the end of dynamics (dBF)"
    write(iulog,*) "and beginning of physics (using dynamics in physics energy; dyBF) the same?"
    write(iulog,*) ""

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

    write(iulog,*)" "
    write(iulog,*)"----------------------------------------------------------------------------"
    write(iulog,*)" Consistency check 2: total energy increment on dynamics decomposition      "
    write(iulog,*)" on an A-grid (physics grid) the same as physics increment (also on A-grid)?"
    write(iulog,*)" Note that wind tendencies are mapped to C-grid in MPAS dycore and added    "
    write(iulog,*)" throughout the time-integration leading to additional physics-dynamics     "
    write(iulog,*)" coupling errors                                                            "
    write(iulog,*)"----------------------------------------------------------------------------"
    write(iulog,*)" "

    previous_dEdt_phys_dyn_coupl_err_Agrid = dEdt_phys_total_in_dyn(1)-dEdt_phys_total_dynE(1)
    diff = abs_diff(dEdt_phys_total_dynE(1),dEdt_phys_total_in_dyn(1),pf=pf)
    write(iulog,'(A50,E8.2,A7,A5)')" dE/dt physics-dynamics coupling errors (A-grid) ",diff," W/M^2 ",pf
    write(iulog,*)" "
    if (abs(diff)>eps) then
      do i=1,4
        write(iulog,*) str(i),":"
        write(iulog,*) "======"
        diff = abs_diff(dEdt_phys_total_dynE(i),dEdt_phys_total_in_dyn(i),pf=pf)
        write(iulog,*) "dE/dt physics-dynamics coupling errors (diff) ",diff
        write(iulog,*) "dE/dt physics total in dynamics (dAM-dBF)  ",dEdt_phys_total_in_dyn(i)
        write(iulog,*) "dE/dt physics total in physics  (pAM-pBF)  ",dEdt_phys_total_dynE(i)
        write(iulog,*) " "
        write(iulog,*) "      physics total = parameterizations + efix + dry-mass adjustment"
        write(iulog,*) " "
      end do
    end if
    write(iulog,*)" "
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" MPAS dycore energy tendencies"
    write(iulog,*)"------------------------------------------------------------"
    write(iulog,*)" "
    write(iulog,*)" Energy diagnostics have not been implemented in the MPAS"
    write(iulog,*)" dynamical core so a detailed budget is not available."
    write(iulog,*)" "
    write(iulog,*)" dE/dt adiabatic dynamical core must therefore be estimated"
    write(iulog,*)" from"
    write(iulog,*)" "
    write(iulog,*)" dE/dt adiabatic dycore (t=n-1) = "
    write(iulog,*)"    -dE/dt dry mass adjustment (t=n-1) +"
    write(iulog,*)"    -dE/dt energy fixer(t=n)"
    write(iulog,*)"    -dE/dt physics-dynamics coupling errors (t=n-1)"
    write(iulog,*)" "
    dEdt_dycore_and_pdc_estimated_from_efix = -dEdt_efix_dynE(1)-previous_dEdt_dry_mass_adjust
    write(iulog,'(A34,F6.2,A6)') "                                = ",dEdt_dycore_and_pdc_estimated_from_efix," W/M^2"
    write(iulog,*)" "
    write(iulog,*)" assuming no physics-dynamics coupling errors, that is,"
    write(iulog,*)" dE/dt physics-dynamics coupling errors (t=n-1) = 0"
    write(iulog,*)" "
    write(iulog,*)" For MPAS the physics-dynamics coupling errors include:"
    write(iulog,*)"  - `dribbling' temperature and wind tendencies during the"
    write(iulog,*)"    dynamical core time-integration"
    write(iulog,*)"  - mapping wind tendencies from A to C grid"
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
          write(iulog,*) "Error: dMASS/dt parameterizations (pAP-pBP) /= dMASS/dt physics total (pAM-pBF)"
          write(iulog,*) "dMASS/dt parameterizations   (pAP-pBP) ",dMdt_parameterizations," Pa/m^2/s"
          write(iulog,*) "dMASS/dt physics total       (pAM-pBF) ",dMdt_phys_total," Pa/m^2/s"
          call endrun(subname//"mass change not only due to parameterizations. See atm.log")
        end if
        write(iulog,*)"  "
        !
        ! check if mass change in physics is the same as dynamical core
        !
        call cam_budget_get_global('dAM-dBF',m_cnst,dMdt_phys_total_in_dyn)
        dMdt_PDC = dMdt_phys_total-dMdt_phys_total_in_dyn
        write(iulog,fmtm)"   Mass physics-dynamics coupling error          ",dMdt_PDC," Pa/m^2/s"
        write(iulog,*)" "
        if (abs(dMdt_PDC)>eps_mass) then
          write(iulog,fmtm)"   dMASS/dt physics tendency in dycore (dAM-dBF) ",dMdt_phys_total_in_dyn," Pa/m^2/s"
          write(iulog,fmtm)"   dMASS/dt total physics                        ",dMdt_phys_total," Pa/m^2/s"
        end if
      end if
    end do
    !
    ! save dry-mass adjustment to avoid sampling error
    !
    previous_dEdt_dry_mass_adjust = dEdt_dme_adjust_dynE(1)
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
   if (present(pf)) then
     if (abs_diff>eps) then
       pf = ' FAIL'
     else
       pf = ' PASS'
     end if
   end if
 end function abs_diff
end module dycore_budget

