module dycore_budget

implicit none

public :: print_budget


!=========================================================================================
contains
!=========================================================================================

subroutine print_budget()

!!$  use budgets,                only : budget_num, budget_info, budget_ind_byname, budget_get_global
!!$  
!!$  ! Local variables
!!$  integer :: b_ind,s_ind,is1,is2
!!$  logical :: budget_outfld
!!$  character(len=64)    :: name_out1,name_out2,name_out3,name_out4,name_out5,budget_name,name_out(9)
!!$  character(len=3)     :: budget_pkgtype,budget_optype  ! budget type phy or dyn
!!$  real(r8),allocatable :: tmp(:,:)
!!$  real(r8), pointer :: te_budgets(:,:,:)! energy/mass budgets se,ke,wv,liq,ice
!!$  integer, pointer :: budgets_cnt(:) ! budget counts for normalizating sum
!!$  integer          :: i
!!$  character(len=*), parameter :: subname = 'check_energy:print_budgets'
!!$
!!$  real(r8)          :: ph_param,ph_EFIX,ph_DMEA,ph_param_and_efix,ph_phys_total
!!$  real(r8)          :: dy_param,dy_EFIX,dy_DMEA,dy_param_and_efix,dy_phys_total
!!$  real(r8)          :: mpas_param,mpas_dmea,mpas_phys_total, dycore, err, param, pefix, pdmea, param_mpas, phys_total
!!$  integer           :: m_cnst  
!!$  !--------------------------------------------------------------------------------------
!!$
!!$  if (masterproc) then
!!$     call budget_get_global('phAP-phBP',1,ph_param)
!!$     call budget_get_global('phBP-phBF',1,ph_EFIX)
!!$     call budget_get_global('phAM-phAP',1,ph_DMEA)
!!$     call budget_get_global('phAP-phBF',1,ph_param_and_efix)
!!$     call budget_get_global('phAM-phBF',1,ph_phys_total)
!!$     
!!$     call budget_get_global('dyAP-dyBP',1,dy_param)
!!$     call budget_get_global('dyBP-dyBF',1,dy_EFIX)
!!$     call budget_get_global('dyAM-dyAP',1,dy_DMEA)
!!$     call budget_get_global('dyAP-dyBF',1,dy_param_and_efix)
!!$     call budget_get_global('dyAM-dyBF',1,dy_phys_total)
!!$     
!!$     call budget_get_global('dAP-dBF',1,mpas_param)
!!$     call budget_get_global('dAM-dAP',1,mpas_dmea)
!!$     call budget_get_global('dAM-dBF',1,mpas_phys_total)
!!$     
!!$     
!!$     write(iulog,*)" "
!!$     write(iulog,*)"================================================================================="
!!$     write(iulog,*)"|                                                                               |"
!!$     write(iulog,*)"| ANALYSIS OF ENERGY DIAGNOSTICS IN PHYSICS                                      |"
!!$     write(iulog,*)"|                                                                               |"
!!$     write(iulog,*)"================================================================================="
!!$     write(iulog,*)" "
!!$     write(iulog,*)"-------------------------------------------------------"
!!$     write(iulog,*)" CAM physics energy increments (in pressure coordinate)"
!!$     write(iulog,*)"-------------------------------------------------------"
!!$     write(iulog,*)" "
!!$     write(iulog,*)"dE/dt parameterizations no efix (param) (pAP-pBP) ",ph_param," W/M^2"
!!$     write(iulog,*)"dE/dt energy fixer (efix)               (pBP-pBF) ",ph_EFIX," W/M^2"
!!$     write(iulog,*)"NOTE: energy fixer uses energy formula consistent with dycore (so this is not p-based for MPAS) "
!!$     write(iulog,*)"dE/dt parameterizations + efix          (pAP-pBF) ",ph_param_and_efix," W/M^2"
!!$     write(iulog,*)" "
!!$     write(iulog,*)"dE/dt dry mass adjustment (pwork)       (pAM-pAP) ",ph_DMEA," W/M^2"
!!$     write(iulog,*)"dE/dt physics total (phys)              (pAM-pBF) ",ph_phys_total," W/M^2"
!!$     write(iulog,*)" "
!!$     dycore = -ph_EFIX-ph_DMEA
!!$     write(iulog,*)"Dycore TE dissipation estimated from physics in pressure coordinate ",dycore," W/M^2"
!!$     write(iulog,*)"(assuming no physics-dynamics coupling errors)    "
!!$     write(iulog,*)" "
!!$     write(iulog,*)"-----------------------------------------------------------------------------------"
!!$     write(iulog,*)" CAM physics dynamical core consistent energy increments (for MPAS in z coordinate)"
!!$     write(iulog,*)"-----------------------------------------------------------------------------------"
!!$     write(iulog,*)" "
!!$     write(iulog,*)"dE/dt parameterizations no efix (param) (dyAP-dyBP) ",dy_param," W/M^2"
!!$     write(iulog,*)"dE/dt energy fixer (efix)               (dyBP-dyBF) ",dy_EFIX," W/M^2"
!!$     write(iulog,*)"dE/dt parameterizations + efix          (dyAP-dyBF) ",dy_param_and_efix," W/M^2"
!!$     write(iulog,*)" "
!!$     write(iulog,*)"dE/dt dry mass adjustment (pwork)       (dyAM-dyAP) ",dy_DMEA," W/M^2"
!!$     write(iulog,*)"dE/dt physics total (phys)              (dyAM-dyBF) ",dy_phys_total," W/M^2"
!!$     write(iulog,*)" "
!!$     dycore = -dy_EFIX-dy_DMEA
!!$     write(iulog,*)"Dycore TE dissipation estimated from physics with dycore energy ",dycore," W/M^2"
!!$     write(iulog,*)"(assuming no physics-dynamics coupling errors; -efix-dme_adjust)    "
!!$     write(iulog,*)" "
!!$     
!!$     
!!$     write(iulog,*)"================================================================================="
!!$     write(iulog,*)"|                                                                               |"
!!$     write(iulog,*)"| ANALYSIS OF ENERGY DIAGNOSTICS IN PHYSICS dp_coupling (MPAS)                  |"
!!$     write(iulog,*)"|                                                                               |"
!!$     write(iulog,*)"================================================================================="
!!$     write(iulog,*)" "
!!$     write(iulog,*)"  "
!!$     write(iulog,*)"dE/dt parameterizations + efix (total physics increment) in MPAS   "
!!$     write(iulog,*)"when adding as one increment - no dribbling (dAP-dBF) ",mpas_param," W/M^2"
!!$     err = ph_param_and_efix-mpas_param
!!$     write(iulog,*)"compare to same tendency in physics (MUST BE SMALL!)  ",err," W/M^2"
!!$     write(iulog,*)" "
!!$     write(iulog,*)"dE/dt dry mass adjustment in MPAS           (dAM-dAP) ",mpas_dmea," W/M^2"
!!$     err = dy_DMEA-mpas_dmea
!!$     write(iulog,*)"compare to same tendency in physics (MUST BE SMALL!)  ",err," W/M^2"
!!$     
!!$     do m_cnst=4,6
!!$        
!!$        if (m_cnst.eq.4) then 
!!$           
!!$           write(iulog,*)"Water vapor budget"
!!$           write(iulog,*)"------------------"
!!$        end if
!!$        if (m_cnst.eq.5) then
!!$           write(iulog,*)"Cloud liquid budget"
!!$           write(iulog,*)"------------------"
!!$        end if
!!$        if (m_cnst.eq.6) then 
!!$           write(iulog,*)"Cloud ice budget"
!!$           write(iulog,*)"------------------"
!!$        end if
!!$        write(iulog,*)""
!!$        
!!$        call budget_get_global('phAP-phBP',m_cnst,param)
!!$        call budget_get_global('phBP-phBF',m_cnst,pEFIX)
!!$        call budget_get_global('phAM-phAP',m_cnst,pDMEA)
!!$        
!!$        call budget_get_global('dAM-dBF',m_cnst,param_mpas)
!!$        call budget_get_global('phAM-phBF',m_cnst,phys_total)
!!$        
!!$        write(iulog,*)"dMASS/dt energy fixer                      (pBP-pBF) ",pEFIX," Pa"
!!$        write(iulog,*)"dMASS/dt parameterizations                 (pAP-pBP) ",param," Pa"
!!$        write(iulog,*)"dMASS/dt dry mass adjustment               (pAM-pAP) ",pDMEA," Pa"
!!$        write(iulog,*)""
!!$        write(iulog,*)""
!!$        write(iulog,*)"dMass/dt physics total in MPAS             (dAM-dBF) ",param_mpas," Pa"
!!$        err = (param_mpas-param)
!!$        write(iulog,*)"Is mass budget closed?          (pAP-pBP)- (dAM-dBF) ",err
!!$        write(iulog,*)"---------------------------------------------------------------------------------------------------"
!!$        write(iulog,*)" "
!!$     end do
!!$  end if
end subroutine print_budget
!=========================================================================================

end module dycore_budget

