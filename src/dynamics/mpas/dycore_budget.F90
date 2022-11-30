module dycore_budget
use shr_kind_mod, only: r8=>shr_kind_r8  
implicit none

public :: print_budget
real(r8), parameter :: eps = 1.0E-9_r8
real(r8), save :: previous_dEdt_adiabatic_dycore = 0.0_r8
!=========================================================================================
contains
!=========================================================================================

subroutine print_budget()

  use budgets,            only: budget_get_global
  use spmd_utils,         only: masterproc
  use cam_logfile,        only: iulog
  use cam_abortutils,     only: endrun
  ! Local variables
  integer :: b_ind,s_ind,is1,is2
  logical :: budget_outfld
  character(len=64)    :: name_out1,name_out2,name_out3,name_out4,name_out5,budget_name,name_out(9)
  character(len=3)     :: budget_pkgtype,budget_optype  ! budget type phy or dyn
  real(r8),allocatable :: tmp(:,:)
  real(r8), pointer :: te_budgets(:,:,:)! energy/mass budgets se,ke,wv,liq,ice
  integer, pointer :: budgets_cnt(:) ! budget counts for normalizating sum
  integer          :: i
  character(len=*), parameter :: subname = 'check_energy:print_budgets'

  real(r8)          :: ph_param,ph_EFIX,ph_dmea,ph_param_and_efix,ph_phys_total
  real(r8)          :: dy_param,dy_EFIX,dy_dmea,dy_param_and_efix,dy_phys_total
  real(r8)          :: mpas_param,mpas_dmea,mpas_phys_total, dycore, err, param, pefix, pdmea, param_mpas, phys_total
  real(r8)          :: E_dBF, E_dyBF
  real(r8)          :: diff
  integer           :: m_cnst  
  !--------------------------------------------------------------------------------------

  if (masterproc) then
     call budget_get_global('phAP-phBP',1,ph_param)
     call budget_get_global('phBP-phBF',1,ph_EFIX)
     call budget_get_global('phAM-phAP',1,ph_dmea)
     call budget_get_global('phAP-phBF',1,ph_param_and_efix)
     call budget_get_global('phAM-phBF',1,ph_phys_total)
     
     call budget_get_global('dyAP-dyBP',1,dy_param)
     call budget_get_global('dyBP-dyBF',1,dy_EFIX)
     call budget_get_global('dyAM-dyAP',1,dy_dmea)
     call budget_get_global('dyAP-dyBF',1,dy_param_and_efix)
     call budget_get_global('dyAM-dyBF',1,dy_phys_total)
     
     call budget_get_global('dAP-dBF',1,mpas_param)
     call budget_get_global('dAM-dAP',1,mpas_dmea)
     call budget_get_global('dAM-dBF',1,mpas_phys_total)

     write(iulog,*)" "
     write(iulog,*)" Total energy diagnostics introduced in Lauritzen and Williamson (2019)"     
     write(iulog,*)" (DOI:10.1029/2018MS001549)"
     write(iulog,*)" "    
     write(iulog,*)"------------------------------------------------------------"     
     write(iulog,*)"Physics time loop"
     write(iulog,*)"------------------------------------------------------------"     
     write(iulog,*)" "
     write(iulog,*)"phBF: state passed to parameterizations, before energy fixer"
     write(iulog,*)"phBP: after energy fixer, before parameterizations"
     write(iulog,*)"phAP: after last phys_update in parameterizations and state "
     write(iulog,*)"      saved for energy fixer"
     write(iulog,*)"phAM: after dry mass correction"
     write(iulog,*)"history files saved off here"     
     write(iulog,*)" "
     write(iulog,*)"------------------------------------------------------------"
     write(iulog,*)" CAM physics energy tendencies (using pressure coordinate)  "
     write(iulog,*)"------------------------------------------------------------"
     write(iulog,*)" "     
     write(iulog,'(a40,F6.2,a6)')"dE/dt energy fixer          (phBP-phBF) ",ph_EFIX," W/M^2"     
     write(iulog,'(a40,F6.2,a6)')"dE/dt all parameterizations (phAP-phBP) ",ph_param," W/M^2"
     write(iulog,'(a40,F6.2,a6)')"dE/dt dry mass adjustment   (phAM-phAP) ",ph_dmea," W/M^2"
     write(iulog,'(a40,F6.2,a6)')"dE/dt physics total         (phAM-phBF) ",ph_phys_total," W/M^2"
     write(iulog,*)" "   
     write(iulog,*)" "
     write(iulog,*)"------------------------------------------------------------"
     write(iulog,*)" CAM physics energy tendencies (using z coordinate)         "
     write(iulog,*)"------------------------------------------------------------"
     write(iulog,*)" "     
     write(iulog,'(a40,F6.2,a6)')"dE/dt energy fixer          (dyBP-dyBF) ",dy_EFIX," W/M^2"     
     write(iulog,'(a40,F6.2,a6)')"dE/dt all parameterizations (dyAP-dyBP) ",dy_param," W/M^2"
     write(iulog,'(a40,F6.2,a6)')"dE/dt dry mass adjustment   (dyAM-dyAP) ",dy_dmea," W/M^2"
     write(iulog,'(a40,F6.2,a6)')"dE/dt physics total         (dyAM-dyBF) ",dy_phys_total," W/M^2"
     write(iulog,*)" "   
     dycore = -dy_EFIX-dy_dmea
     write(iulog,*)"Dycore TE dissipation estimated from physics with dycore energy ",dycore," W/M^2"
     write(iulog,*)"(assuming no physics-dynamics coupling errors; -efix-dme_adjust)    "
     write(iulog,*)" "

     write(iulog,*) " "     
     write(iulog,*) "-dE/dt energy fixer = dE/dt dry mass adjustment              +"
     write(iulog,*) "                      dE/dt dycore                           +"
     write(iulog,*) "                      dE/dt physics-dynamics coupling errors +"
     write(iulog,*) "                      dE/dt energy formula differences        "
     write(iulog,*) " "
     write(iulog,*) "(equation 23 in Lauritzen and Williamson (2019))"
     write(iulog,*) " "  
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
     write(iulog,*) ""
     write(iulog,*) "Is globally integrated total energy of state at the end of dynamics (dBF)"
     write(iulog,*) "and beginning of physics (dyBF) the same?"
     write(iulog,*) ""     
     call budget_get_global('dBF',1,E_dBF)  !state passed to physics
     call budget_get_global('dyBF',1,E_dyBF)!state beginning physics
     if (abs(E_dyBF)>eps) then
       diff = abs_diff(E_dBF,E_dyBF)
       if (abs(diff)<eps) then
         write(iulog,*)"yes. (dBF-dyBF)/dyBF =",diff
         write(iulog,*)"E_dBF=",E_dBF,"; E_dyBF=",E_dyBF
       else
         write(iulog,*)"no. (dBF-dyBF)/dyBF =",diff
         write(iulog,*)"E_dBF=",E_dBF,"; E_dyBF=",E_dyBF
         write(iulog,*)"Error in physics dynamics coupling!"
         call endrun('dycore_budget module: Error in physics dynamics coupling')
       end if
     end if
     write(iulog,*)" "
     write(iulog,*)"------------------------------------------------------------"
     write(iulog,*)" MPAS energy tendencies                                     "
     write(iulog,*)"------------------------------------------------------------"
     write(iulog,*)" "
     write(iulog,*)              "dE/dt all parameterizations+ "
     write(iulog,'(a42,F6.2,a6)')"dE/dt energy fixer            (dAP-dBF) ",mpas_param," W/M^2"
     write(iulog,'(a42,F6.2,a6)')"dE/dt dry mass adjustment     (dAM-dAP) ",mpas_dmea," W/M^2"
     write(iulog,*)" "   
     write(iulog,*)"Are these values consistent with CAM physics dE/dt's?"
     write(iulog,*)" "
     diff = abs_diff(mpas_param,dy_param+dy_EFIX)
     write(iulog,*)"Physics tendency: ((dAP-dBF)-(dyAP-dyBF))/(dyAP-dyBF) =",diff
     if (abs(diff)>eps) then
!       call endrun('dycore_budget module: physics tendency in dynamics error')
     endif
     diff = abs_diff(mpas_dmea,dy_dmea)
     write(iulog,*)"Dry-mass adj.   : ((dAM-dAP)-(dyAM-dyAP))/(dyAM-dyAP) =",diff
     if (abs(diff)>eps) then
       write(iulog,*) "error: dry-mass adjustment in dynamics error"
!       call endrun('dycore_budget module: dry-mass adjustment in dynamics error')
     endif
     write(iulog,*)" "     
     do m_cnst=4,6
       write(iulog,*)"------------------------------------------------------------"
       if (m_cnst.eq.4) write(iulog,*)"Water vapor mass budget"
       if (m_cnst.eq.5) write(iulog,*)"Liquid water mass budget"
       if (m_cnst.eq.6) write(iulog,*)"Ice water mass budget"
       write(iulog,*)"------------------------------------------------------------"        
       call budget_get_global('phAP-phBP',m_cnst,param)
       call budget_get_global('phBP-phBF',m_cnst,pEFIX)
       call budget_get_global('phAM-phAP',m_cnst,pdmea)
       
       call budget_get_global('dAM-dBF',m_cnst,param_mpas)
       call budget_get_global('phAM-phBF',m_cnst,phys_total)
       
       write(iulog,*)"dMASS/dt energy fixer               (pBP-pBF) ",pEFIX," Pa"
       write(iulog,*)"dMASS/dt parameterizations          (pAP-pBP) ",param," Pa"
       write(iulog,*)"dMASS/dt dry mass adjustment        (pAM-pAP) ",pdmea," Pa"
       write(iulog,*)""
       write(iulog,*)""
       write(iulog,*)"dMass/dt physics total in MPAS      (dAM-dBF) ",param_mpas," Pa"
       err = (param_mpas-param)
       write(iulog,*)"Is mass budget closed?    (pAP-pBP)-(dAM-dBF) ",err
       write(iulog,*)"---------------------------------------------------------------------------------------------------"
       write(iulog,*)" "
     end do
   end if
 end subroutine print_budget
 !=========================================================================================
 function abs_diff(a,b)
   real(r8), intent(in)  :: a,b
  real(r8)              :: abs_diff
  if (abs(b)>eps) then
    abs_diff = abs((b-a)/b)
  else
    abs_diff = abs(b-a)
  end if
end function abs_diff
end module dycore_budget

