module dycore_budget
use shr_kind_mod, only: r8=>shr_kind_r8  
implicit none

public :: print_budget
real(r8), parameter :: eps = 1.0E-11_r8
real(r8), save :: previous_dEdt_adiabatic_dycore = 0.0_r8
!=========================================================================================
contains
!=========================================================================================

subroutine print_budget(hstwr)

  use budgets,            only: budget_get_global, thermo_budget_histfile_num, thermo_budget_history
  use spmd_utils,         only: masterproc
  use cam_logfile,        only: iulog
  use cam_abortutils,     only: endrun
  use cam_thermo,             only: teidx, thermo_budget_vars_descriptor, thermo_budget_num_vars, thermo_budget_vars_massv

  ! arguments
  logical, intent(in) :: hstwr(:)

  ! Local variables
  real(r8),allocatable :: tmp(:,:)
  integer          :: i
  character(len=*), parameter :: subname = 'check_energy:print_budgets'

  real(r8)          :: ph_param,ph_EFIX,ph_dmea,ph_param_and_efix,ph_phys_total
  real(r8)          :: dy_param,dy_EFIX,dy_dmea,dy_param_and_efix,dy_phys_total
  real(r8)          :: mpas_param,mpas_dmea,mpas_phys_total, dycore, err, param, pefix, pdmea, param_mpas, phys_total
  real(r8)          :: E_dBF, E_dyBF
  real(r8)          :: diff
  integer           :: m_cnst 
  character(LEN=*), parameter :: fmt  = "(a40,F6.2,a1,F6.2,a1,E10.2,a4)"
  character(LEN=*), parameter :: fmt2 = "(a40,F6.2,a3)"
  character(LEN=5)  :: pf! pass or fail identifier
  !--------------------------------------------------------------------------------------

  if (masterproc .and. thermo_budget_history .and. hstwr(thermo_budget_histfile_num)) then
     call budget_get_global('phAP-phBP',teidx,ph_param)
     call budget_get_global('phBP-phBF',teidx,ph_EFIX)
     call budget_get_global('phAM-phAP',teidx,ph_dmea)
     call budget_get_global('phAP-phBF',teidx,ph_param_and_efix)
     call budget_get_global('phAM-phBF',teidx,ph_phys_total)
     
     call budget_get_global('dyAP-dyBP',teidx,dy_param)
     call budget_get_global('dyBP-dyBF',teidx,dy_EFIX)
     call budget_get_global('dyAM-dyAP',teidx,dy_dmea)
     call budget_get_global('dyAP-dyBF',teidx,dy_param_and_efix)
     call budget_get_global('dyAM-dyBF',teidx,dy_phys_total)
     
     call budget_get_global('dAP-dBF',teidx,mpas_param)
     call budget_get_global('dAM-dAP',teidx,mpas_dmea)
     call budget_get_global('dAM-dBF',teidx,mpas_phys_total)

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
     write(iulog,*)" dBF: dynamics state before physics (d_p_coupling)"
     write(iulog,*)" dAP: dynamics state with T,u,V increment but not incl water changes"
     write(iulog,*)" dAM: dynamics state with full physics increment (incl. water)"
     write(iulog,*)" "
     write(iulog,*)"Note that these energies are computed using the dynamical core"
     write(iulog,*)"state variables which may be different from the physics prognostic"
     write(iulog,*)"variables."
     write(iulog,*)" "
     write(iulog,*)" "
     write(iulog,*)"Consistency check 0:"
     write(iulog,*)"--------------------"
     write(iulog,*)" "
     write(iulog,*)"For energetic consistency we require that dE/dt [W/m^2] from energy "
     write(iulog,*)"fixer and all parameterizations computed using physics E and"
     write(iulog,*)"dycore in physics E are the same! Checking:"
     write(iulog,*)" "
     write(iulog,*)  "                                        xx=ph   xx=dy  norm. diff."
     write(iulog,*)  "                                        -----   -----  -----------"
     diff = abs_diff(ph_EFIX,dy_EFIX,pf=pf)
     write(iulog,fmt)"dE/dt energy fixer          (xxBP-xxBF) ",ph_EFIX,      " ",dy_EFIX," ",diff,pf

     diff = abs_diff(ph_param,dy_param,pf=pf)
     write(iulog,fmt)"dE/dt all parameterizations (xxAP-xxBP) ",ph_param,     " ",dy_param," ",diff,pf
     if (diff>eps) write(iulog,*)"FAIL"
     write(iulog,*)" "
     write(iulog,*)" "
     write(iulog,*)"dE/dt from dry-mass adjustment will differ if dynamics and physics use"
     write(iulog,*)"different energy definitions! Checking:"
     write(iulog,*)" "
     diff = ph_dmea-dy_dmea
     write(iulog,*)  "                                        xx=ph   xx=dy  difference"
     write(iulog,*)  "                                        -----   -----  -----------"
     write(iulog,fmt)"dE/dt dry mass adjustment   (xxAM-xxAP) ",ph_dmea,      " ",dy_dmea," ",diff
     write(iulog,*)" "
     write(iulog,*)" "
     write(iulog,*)"Some energy budget observations:"
     write(iulog,*)"--------------------------------"
     write(iulog,*)" "
     write(iulog,*)"Note that total energy fixer fixes:"
     write(iulog,*) " "     
     write(iulog,*) "-dE/dt energy fixer = dE/dt dry mass adjustment              +"
     write(iulog,*) "                      dE/dt dycore                           +"
     write(iulog,*) "                      dE/dt physics-dynamics coupling errors"
     write(iulog,*) " "
     write(iulog,*) "(equation 23 in Lauritzen and Williamson (2019))"
     write(iulog,*) " "  
     dycore = -dy_EFIX-dy_dmea
     write(iulog,*)"Hence the dycore E dissipation estimated from energy fixer is ",dycore," W/M^2"
     write(iulog,*)"(assuming no physics-dynamics coupling errors)"
     write(iulog,*)" "

!     dycore = -ph_EFIX-ph_dmea
!     dycore = -ph_EFIX-previous_dEdt_dry_mass_adjust
!     write(iulog,*) ""
!     write(iulog,*) "Dycore TE dissipation estimated from physics in pressure coordinate:"
!     write(iulog,*) "(note to avoid sampling error we need dE/dt from previous time-step)"
!     write(iulog,*) ""
!     write(iulog,*) "dE/dt adiabatic dycore estimated from physics (t=n-1) = "
!     write(iulog,'(a58,F6.2,a6)') "-dE/dt energy fixer(t=n)-dE/dt dry-mass adjust(t=n-1) = ",dycore," W/M^2"
!     write(iulog,*) ""
!     write(iulog,'(a58,F6.2,a6)') "dE/dt adiabatic dycore computed in dycore (t=n-1)     = ",&
!          previous_dEdt_adiabatic_dycore," W/M^2"
!     write(iulog,'(a58,F6.2,a6)') "dE/dt dry-mass adjust  (t=n-1)                        = ",&
!          previous_dEdt_dry_mass_adjust," W/M^2"
!     write(iulog,*) ""
!     if (abs(previous_dEdt_adiabatic_dycore)>eps) then
!       diff = abs((dycore-previous_dEdt_adiabatic_dycore)/previous_dEdt_adiabatic_dycore)
!       if (diff>eps) then
!         write(iulog,*) "energy budget not closed: previous_dEdt_adiabatic_dycore <> dycore"
!         write(iulog,*) "normalized difference is:",diff         
!         call endrun('dycore_budget module: physics energy budget consistency error 2')
!       end if
!     end if
     write(iulog,*) " "     
     write(iulog,*) "-------------------------------------------------------------------"
     write(iulog,*) " Consistency check 1: state passed to physics same as end dynamics?"
     write(iulog,*) "-------------------------------------------------------------------"
     write(iulog,*) " "     
     write(iulog,*) "Is globally integrated total energy of state at the end of dynamics (dBF)"
     write(iulog,*) "and beginning of physics (dyBF) the same?"
     write(iulog,*) ""     
     call budget_get_global('dBF',teidx,E_dBF)  !state passed to physics
     call budget_get_global('dyBF',teidx,E_dyBF)!state beginning physics
     if (abs(E_dyBF)>eps) then
       diff = abs_diff(E_dBF,E_dyBF)
       if (abs(diff)<eps) then
         write(iulog,*)"yes. (dBF-dyBF)/dyBF =",diff
       else
         write(iulog,*)"no. (dBF-dyBF)/dyBF =",diff
         write(iulog,*)"E_dBF=",E_dBF,"; E_dyBF=",E_dyBF
         write(iulog,*)"Error in physics dynamics coupling!"
!         call endrun('dycore_budget module: Error in physics dynamics coupling')
       end if
     end if
     write(iulog,*)" "
     write(iulog,*)"-------------------------------------------------------------------------"
     write(iulog,*)" Consistency check 2: total energy increment in dynamics same as physics?"
     write(iulog,*)"-------------------------------------------------------------------------"
     write(iulog,*)" "
     diff = abs_diff(mpas_param,dy_param+dy_EFIX,pf=pf)
     write(iulog,*)"Increment all parameterizations?"
     write(iulog,*)"            ((dAP-dBF)-(dyAP-dyBF))/(dyAP-dyBF) =",diff,pf
     diff = abs_diff(mpas_dmea,dy_dmea,pf=pf)
     write(iulog,*)"Increment dry-mass adj.?"
     write(iulog,*)"            ((dAM-dAP)-(dyAM-dyAP))/(dyAM-dyAP) =",diff,pf
     write(iulog,*)" "     
     write(iulog,*)" "     
     do m_cnst=1,thermo_budget_num_vars
        if (thermo_budget_vars_massv(m_cnst)) then 
           write(iulog,*)"------------------------------------------------------------"
           write(iulog,*)thermo_budget_vars_descriptor(m_cnst)//" budget"
           write(iulog,*)"------------------------------------------------------------"        
           call budget_get_global('phAP-phBP',m_cnst,param)
           call budget_get_global('phBP-phBF',m_cnst,pEFIX)
           call budget_get_global('phAM-phAP',m_cnst,pdmea)
           
           call budget_get_global('dAM-dBF',m_cnst,param_mpas)
           call budget_get_global('phAM-phBF',m_cnst,phys_total)
           
           write(iulog,fmt2)"dMASS/dt energy fixer               (pBP-pBF) ",pEFIX," Pa"
           write(iulog,fmt2)"dMASS/dt parameterizations          (pAP-pBP) ",param," Pa"
           write(iulog,fmt2)"dMASS/dt dry mass adjustment        (pAM-pAP) ",pdmea," Pa"
           write(iulog,fmt2)"dMass/dt physics total in MPAS      (dAM-dBF) ",param_mpas," Pa"
           err = (param_mpas-param)
           write(iulog,*)""
           write(iulog,*)"Is mass budget closed?    (pAP-pBP)-(dAM-dBF) ",err
           write(iulog,*)"-----------------------------------------------------------------"
           write(iulog,*)" "
           if (err>eps) write(iulog,*)" MASS BUDGET ERROR"
        end if
     end do
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

