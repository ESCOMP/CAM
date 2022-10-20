module dycore_budget
use shr_kind_mod, only: r8=>shr_kind_r8
implicit none

public :: print_budget
real(r8), parameter :: eps = 1.0E-12_r8

!=========================================================================================
contains
!=========================================================================================

subroutine print_budget()

  use spmd_utils,             only: masterproc
  use cam_abortutils,         only: endrun  
  use cam_logfile,            only: iulog
  use shr_kind_mod,           only: r8=>shr_kind_r8
  use budgets,                only: budget_get_global, is_budget
  use dimensions_mod,         only: lcp_moist,qsize
  use control_mod,            only: ftype
  ! Local variables
  integer          :: i
  character(len=*), parameter :: subname = 'check_energy:print_budgets'

  real(r8)          :: ph_param,ph_EFIX,ph_DMEA,ph_phys_total
  real(r8)          :: dy_param,dy_EFIX,dy_DMEA,dy_param_and_efix,dy_phys_total
  real(r8)          :: se_param,se_dmea,se_phys_total, dycore, err, param, pefix, &
                       pdmea, phys_total, dyn_total, dyn_phys_total, &
                       rate_of_change_2D_dyn, rate_of_change_vertical_remapping, &
                       diffusion_del4, diffusion_fric, diffusion_del4_tot, diffusion_sponge, &
                       diffusion_total, twoDresidual, rate_of_change_physics, &
                       rate_of_change_heating_term_put_back_in, rate_of_change_hvis_sponge, &
                       value_pdc, dADIA, ttt, fff, &
                       mass_change__2D_dyn,mass_change__vertical_remapping, &
                       mass_change__heating_term_put_back_in,mass_change__hypervis_total, &
                       error, mass_change__physics, dbd, daf, dar, dad, qneg, val,phbf,ded

  real(r8) :: E_dBF, E_phBF, diff
  

  integer           :: m_cnst, qsize_condensate_loading
  logical           :: te_consistent_version
  !--------------------------------------------------------------------------------------

  qsize_condensate_loading=qsize
  te_consistent_version=.false.
  if (qsize_condensate_loading.eq.1) then
    if (lcp_moist.eqv..false.) then
      write(iulog,*)"Using total energy consistent version: qsize_condensate_loading=1 and cp=cpdry"
      te_consistent_version=.true.
    else
      write(iulog,*)"WARNING: Total energy formulas for dynamics and physics are different:"
      write(iulog,*)"   Dynamics (cp includes water vapor; condensates not thermodynamically active)."
      write(iulog,*)"   Physics (cp=cp_dry in internal energy)."
    end if
  else
    write(iulog,*)"WARNING: Total energy formulaes for dynamics and physics are different"
    write(iulog,*)"in dynamics (cp and dp includes all water variables) and physics (cp=cp_dry in internal energy)."
  end if


  if (masterproc) then
     call budget_get_global('phAP-phBP',1,ph_param)
     call budget_get_global('phBP-phBF',1,ph_EFIX)
     call budget_get_global('phAM-phAP',1,ph_DMEA)
     call budget_get_global('phAM-phBF',1,ph_phys_total)
     
     call budget_get_global('dyAP-dyBP',1,dy_param)
     call budget_get_global('dyBP-dyBF',1,dy_EFIX)
     call budget_get_global('dyAM-dyAP',1,dy_DMEA)
     call budget_get_global('dyAP-dyBF',1,dy_param_and_efix)
     call budget_get_global('dyAM-dyBF',1,dy_phys_total)
     
!jt     call budget_get_global('dAP-dBP',1,se_param)
!jt     call budget_get_global('dAM-dAP',1,se_dmea)
!jt     call budget_get_global('dAM-dBF',1,se_phys_total)

     call budget_get_global('dBF-dED',1,dyn_total)
!jt     call budget_get_global('dAD-dAF',1,dyn_phys_total)
     call budget_get_global('dAD-dBD',1,rate_of_change_2D_dyn)
     call budget_get_global('dAR-dAD',1,rate_of_change_vertical_remapping)
     dADIA = rate_of_change_2D_dyn+rate_of_change_vertical_remapping

     call budget_get_global('dCH-dBH',1,diffusion_del4)
     call budget_get_global('dAH-dCH',1,diffusion_fric)
     call budget_get_global('dAH-dBH',1,diffusion_del4_tot)
     call budget_get_global('dAS-dBS',1,diffusion_sponge)
     diffusion_total      = diffusion_del4_tot+diffusion_sponge
     
     call budget_get_global('dBD-dAF',1,rate_of_change_physics)

     rate_of_change_heating_term_put_back_in = diffusion_fric
     rate_of_change_hvis_sponge = diffusion_sponge
     
     write(iulog,*)" "
     write(iulog,*)" Total energy diagnostics introduced in Lauritzen and Williamson (2019)"     
     write(iulog,*)" (DOI:10.1029/2018MS001549)"
     write(iulog,*)" "    
     
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
     write(iulog,*)" "
     
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
     write(iulog,*)" CAM physics energy tendencies (using pressure coordinate)"
     write(iulog,*)"------------------------------------------------------------"
     write(iulog,*)" "     
     write(iulog,*)"dE/dt energy fixer          (phBP-pBF) ",ph_EFIX," W/M^2"     
     write(iulog,*)"dE/dt all parameterizations (phAP-pBP) ",ph_param," W/M^2"
     write(iulog,*)"dE/dt dry mass adjustment   (pAM-pAP)  ",ph_DMEA," W/M^2"
     write(iulog,*)"dE/dt physics total         (pAM-pBF)  ",ph_phys_total," W/M^2"
     !
     ! consistency check
     !
     if (abs(ph_param+ph_EFIX+ph_DMEA-ph_phys_total)>eps) then
       write(iulog,*) "Physics energy budget not adding up:"
       write(iulog,*) "(phBP-pBF)+(phAP-pBP)+(pAM-pAP) does not add up to (pAM-pBF)",\
       abs(ph_param+ph_EFIX+ph_DMEA-ph_phys_total)
       call endrun('dycore_budget module: physics energy budget consistency error')
     endif
     write(iulog,*) " "     
     write(iulog,*) "-dE/dt energy fixer = dE/dt dry mass adjustment              +"
     write(iulog,*) "                      dE/dt dycore                           +"
     write(iulog,*) "                      dE/dt physics-dynamics coupling errors +"
     write(iulog,*) "                      dE/dt energy formula differences        "
     write(iulog,*) " "
     write(iulog,*) "(equation 23 in Lauritzen and Williamson (2019))"
     write(iulog,*) " "
     !
     ! check for energy formula difference
     !
     write(iulog,*) ""
     write(iulog,*) "Is globally integrated total energy of state at the end of dynamics (dBF)"
     write(iulog,*) "and beginning of physics (phBF) the same?"
     write(iulog,*) ""     
     call budget_get_global('dBF',1,E_dBF)  !state passed to physics
     call budget_get_global('phBF',1,E_phBF)!state beginning physics
     if (abs(E_phBF)>eps) then
       diff = abs_diff(E_dBF,E_phBF)
       if (abs(diff)<eps) then
         write(iulog,*)"yes. (dBF-phBF)/phBF =     (dBF-phBP)=",diff
       else
         write(iulog,*) "no. (dBF-phBF)/phBF =     (dBF-phBP)=",diff,E_dBF,E_phBF
         write(iulog,*) "To run energy consistent version of SE use namelist"
         write(iulog,*) ""
         write(iulog,*) "se_ftype     = 1           !no dribbling of physics tendencies"
         write(iulog,*) "se_lcp_moist =  .false.    !no variable latent heats"
         write(iulog,*) "water_species_in_air = 'Q' !only water vapor energetically active"
         write(iulog,*) ""         
       end if
     end if
     write(iulog,*)"dE/dt physics tendency in dynamics (dBD-dAF) ",rate_of_change_physics," W/M^2"
     write(iulog,*)"dE/dt physics tendency in physics  (pAM-pBF) ",ph_phys_total," W/M^2"
     write(iulog,*) ""
     write(iulog,*) "If there are no dribbling errors and no energy formula inconsistencies"
     write(iulog,*) "these should be the same:",abs_diff(rate_of_change_physics,ph_phys_total)
     
     dycore = -ph_EFIX-ph_DMEA
     write(iulog,*)"Dycore TE dissipation estimated from physics in pressure coordinate ",dycore," W/M^2"
     write(iulog,*)"(assuming no physics-dynamics coupling errors)    "

     write(iulog,*)"================================================================================="
     write(iulog,*)"|                                                                               |"
     write(iulog,*)"| ANALYSIS OF ENERGY DIAGNOSTICS IN DYNAMICS - specific for SE dycore           |"
     write(iulog,*)"|                                                                               |"
     write(iulog,*)"================================================================================="
     write(iulog,*)" "
     
     write(iulog,*)"dE/dt dyn total (dycore+phys tendency   (dBF-dED) ",dyn_total," W/M^2"
     write(iulog,*)"dE/dt total adiabatic dynamics (adiab)            ",dADIA," W/M^2"
     write(iulog,*)"dE/dt 2D dynamics (2D)                  (dAD-dBD) ",rate_of_change_2D_dyn," W/M^2"
     write(iulog,*)"dE/dt vertical remapping (remap)        (dAR-dAD) ",rate_of_change_vertical_remapping," W/M^2"

     write(iulog,*)" "
     write(iulog,*)"Breakdown of 2D dynamics:"
     write(iulog,*)" "
     
     write(iulog,*)"      dE/dt hypervis del4 (hvis)             (dCH-dBH) ",diffusion_del4," W/M^2"
     write(iulog,*)"      dE/dt hypervis frictional heating del4 (dAH-dCH) ",diffusion_fric," W/M^2"
     write(iulog,*)"      dE/dt hypervis del4 total (hvis)       (dAH-dBH) ",diffusion_del4_tot," W/M^2"
     write(iulog,*)"      dE/dt hypervis sponge total            (dAS-dBS) ",diffusion_sponge," W/M^2"
     write(iulog,*)"      dE/dt explicit diffusion total                   ",diffusion_total," W/M^2"

     twoDresidual = rate_of_change_2D_dyn-diffusion_total
     write(iulog,*)"      dE/dt residual (res)                        ",twoDresidual," W/M^2"
      write(iulog,*)""
      write(iulog,*)"================================================================================="
      write(iulog,*)"|                                                                               |"
      write(iulog,*)"| ANALYSIS OF ENERGY DIAGNOSTICS IN DYNAMICS-PHYSICS COMBINED                   |"
      write(iulog,*)"|                                                                               |"
      write(iulog,*)"================================================================================="
      write(iulog,*)""
      value_pdc = ph_phys_total-rate_of_change_physics
      if (te_consistent_version.eqv..true.) then
        write(iulog,*)"Your model is energy consistent (qsize_condensate_loading=1 and cpdry)"
        if (ftype .eq. 1) then
          write(iulog,*)""
          write(iulog,*)"You are using ftype=1 so PDC errors should be zero:"
          write(iulog,*)""


          write(iulog,*)"    dE/dt physics tendency in dynamics (dBD-dAF) should exactly match dE/dt physics total (pAM-pBF): ",value_pdc
          write(iulog,*)""
        else
          write(iulog,*)""
          write(iulog,*)"You are using ftype=0 or 2 so there are PDC errors (dribbling errors):"
          write(iulog,*)""
          write(iulog,*)"   Dribbling errors (pAM-pBF-(dBD-dAF))/dt: ",value_pdc
        end if
!jt        discr = "0       "
!jt        str_pdc = sprintf("%6.3g",10*value_pdc)
      else
        write(iulog,*)"Your model is energy inconsistent (qsize_condensate_loading<>1 and/or cp<>cpdry)"
        write(iulog,*)""
        write(iulog,*)"PDC errors can not be assesed trhough "
        write(iulog,*)""
        write(iulog,*)"   dE/dt physics tendency in dynamics (dBD-dAF) does not match dE/dt physics total (pAM-pBF) due to energy discrepancy:",value_pdc
	write(iulog,*)ph_phys_total,"  ",rate_of_change_physics
!jt        str_pdc = "undef"
      end if
      write(iulog,*)""
      write(iulog,*)"Some more consisitency/budget terms"
      write(iulog,*)"==================================="
      write(iulog,*)""
      write(iulog,*)"Energy fixer fixes dme_adjust (pDMEA), lack of energy conservation in adiabatic"
      write(iulog,*)"dynamical core (dADIA), energy discrepancy (EDIFF) and energy lost/gained in physics-dynamics coupling"
      write(iulog,*)""
      write(iulog,*)"dADIA                     ",dADIA," W/M^2"    
      write(iulog,*)"pDMEA                     ",ph_DMEA," W/M^2"    
      write(iulog,*)"physics-dynamics coupling ",value_pdc," W/M^2"
      write(iulog,*)""
!jt      str="dPDC+EDIFF"
      write(iulog,*)""
      write(iulog,*)"		-energy fixer = DME_adjust+adaib dycore+phys-dyn errors+discr"
      write(iulog,*)"            "
      ttt = -ph_DMEA-dADIA-value_pdc
!jt      discr = -99.0
      write(iulog,*)"          DME_adjust+adaib dycore+phys-dyn errors+discr = ",ttt
      write(iulog,*)"          Energy fixer                                  = ",ph_EFIX
      write(iulog,*)""
      fff = ttt-ph_EFIX
      write(iulog,*)"          Difference                                    = ",fff


      
      call budget_get_global('phBF',1,phbf)
      call budget_get_global('dED',1,ded)
      qneg=phbf-ded
      write(iulog,*)""
      write(iulog,*)" qneg: ",qneg
      write(iulog,*)""

    if (qsize.gt.0) then
      write(iulog,*)""
      write(iulog,*)""
      write(iulog,*)""
      write(iulog,*)"================================================================================="
      write(iulog,*)"|                                                                               |"
      write(iulog,*)"| ANALYSIS OF WATER VAPOR, CLOUD LIQUID AND CLOUD ICE BUDGETS                   |"
      write(iulog,*)"|                                                                               |"
      write(iulog,*)"================================================================================="
      write(iulog,*)""
    end if

!jt    do m_cnst=4,4+qsize-1
    do m_cnst=4,6
      if (m_cnst.eq.4) then 
        write(iulog,*)"Water vapor"
        write(iulog,*)"-----------"
      end if
      if (m_cnst.eq.5) then
        write(iulog,*)"Cloud liquid"
        write(iulog,*)"-----------"
      end if
      if (m_cnst.eq.6) then 
        write(iulog,*)"Cloud ice"
        write(iulog,*)"-----------"
      end if

      call budget_get_global('phBP-phBF',m_cnst,pEFIX)
      call budget_get_global('phAM-phAP',m_cnst,pDMEA)
      call budget_get_global('phAP-phBP',m_cnst,param)
!jt      call budget_get_global('dBF-dED',m_cnst,dyn_total)
      call budget_get_global('phAM-phBF',m_cnst,phys_total)

      write(iulog,*)"dMASS/dt energy fixer                      (pBP-pBF) ",pEFIX," Pa"
      write(iulog,*)"dMASS/dt parameterizations                 (pAP-pBP) ",param," Pa"
      write(iulog,*)"dMASS/dt dry mass adjustment               (pAM-pAP) ",pDMEA," Pa"
      write(iulog,*)" "
      val = pEFIX+pDMEA
      write(iulog,*)"=> dMASS/dt dynamical core (estimated from physics)  "
      write(iulog,*)"   dMASS/dt energy fixer + dMASS/dt dry mass adjustment ",val," Pa"

      write(iulog,*)"=> dMASS/dt physics total                   (pAM-pBF)",phys_total," Pa"
 

      write(iulog,*)"  "
      write(iulog,*)"  "
      write(iulog,*)"  "

      if (is_budget('dAD').and.is_budget('dBD').and.is_budget('dAR').and.is_budget('dCH')) then
        call budget_get_global('dAD-dBD',m_cnst,mass_change__2D_dyn)
        call budget_get_global('dAR-dAD',m_cnst,mass_change__vertical_remapping)
        dADIA = mass_change__2D_dyn+mass_change__vertical_remapping
        write(iulog,*)"dE/dt total adiabatic dynamics                    ",dADIA," Pa"
        write(iulog,*)"dE/dt 2D dynamics                       (dAD-dBD) ",mass_change__2D_dyn," Pa"
        write(iulog,*)" "
        write(iulog,*)"Breakdown of 2D dynamics:"
        write(iulog,*)" "
        call budget_get_global('dAH-dCH',m_cnst,mass_change__heating_term_put_back_in)
        call budget_get_global('dAH-dBH',m_cnst,mass_change__hypervis_total)
        write(iulog,*)"      dE/dt hypervis                    (dAH-dBH) ",mass_change__hypervis_total," Pa"
        write(iulog,*)"      dE/dt frictional heating          (dAH-dCH) ",mass_change__heating_term_put_back_in," Pa"
        error = mass_change__2D_dyn-mass_change__hypervis_total
        write(iulog,*)"      dE/dt residual (time truncation errors)     ",error," Pa"
      end if
      if (is_budget('dAR').and.is_budget('dAD')) then
        call budget_get_global('dAR',m_cnst,dar)
        call budget_get_global('dAD',m_cnst,dad)
        call budget_get_global('dAR-dAD',m_cnst,mass_change__vertical_remapping)
        write(iulog,*)"dE/dt vertical remapping                (dAR-dAD) ",mass_change__vertical_remapping," Pa","dar:",dar,"dad:",dad
      end if
      write(iulog,*)" "
      write(iulog,*)" "

      if (is_budget('dBD').and.is_budget('dAF')) then
        call budget_get_global('dBD',m_cnst,dbd)
        call budget_get_global('dAF',m_cnst,daf)
        call budget_get_global('dBD-dAF',m_cnst,mass_change__physics)
        write(iulog,*)"dE/dt physics tendency in dynamics      (dBD-dAF) ",mass_change__physics," Pa","dbd:",dbd,"daf:",daf
        val = phys_total-mass_change__physics
      end if
      if (is_budget('dBD').and.is_budget('dAF')) then
        if (ftype .eq. 1 .or.ftype .eq. 2) then
          write(iulog,*)" "
          write(iulog,*)"      Consistency check:"
          write(iulog,*)" "
          write(iulog,*)"      dE/dt physics tendency in dynamics (dBD-dAF) should exactly match dE/dt physics total (pAM-pBF):",val
          write(iulog,*)" "
        else
          write(iulog,*)"Dribbling errors (pAM-pBF-(dBD-dAF))",val
        end if
      end if
      write(iulog,*)""
      write(iulog,*)"================================================================================="
      write(iulog,*)""
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
