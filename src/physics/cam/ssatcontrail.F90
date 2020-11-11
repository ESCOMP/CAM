module ssatcontrail

    use shr_kind_mod,   only: r8 => shr_kind_r8
    use ppgrid,         only: pcols, pver
    use cam_history,    only: outfld
    use cam_logfile,    only: iulog
    use physics_types,  only: physics_state, physics_ptend, physics_tend
    use physics_types,  only: physics_ptend_sum, physics_update
    use physics_types,  only: physics_state_copy, physics_ptend_init
    use physconst,      only: cpair,mwdry,mwh2o, gravit, zvir, rair, pi, rearth
!    use phys_buffer,    only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx, pbuf_times
    use physics_buffer, only : pbuf_get_index, pbuf_get_field, physics_buffer_desc, pbuf_set_field,  pbuf_old_tim_idx
    use constituents, only: cnst_get_ind, pcnst    
    use aircraft_emit,    only: aircraft_cnt, spc_name_list
    use geopotential,     only: geopotential_dse
    use phys_grid, only    : get_wght_all_p
    use time_manager,       only: get_curr_date

    implicit none
    private
    save

    public ssatcontrail_d0

    ! Private data
!    real(r8), parameter :: ICIWC0 = 2.0e-6   ! ICIWC = 2 mg/m3, converted to kg/kg here by dividing rho_air
    real(r8), parameter :: rhoi = 500.0_r8             ! density of ice (500 kg/m3)
    real(r8), parameter :: radius = 3.75e-6            ! diameter of ice particle = 7.5 microns

  
contains

    subroutine ssatcontrail_d0(state1,pbuf,dtime,ptend_loc, tend)
    implicit none

    type(physics_state), intent(inout) :: state1
    type(physics_ptend), intent(inout) :: ptend_loc
    type(physics_tend) :: tend
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8),            intent(in)    :: dtime      ! time step
!------------------------Local storage------------------------------------------------------
    real(r8) :: Ma, Mh2o, epsi, Q, eta, p, G, T_contr, eslTc, eslT, RH_contr
    real(r8) :: w, esiT, ws, RH, ei
    integer  :: i,k 
    integer  :: lchnk,ncol
    real(r8) :: contrail(pcols,pver), pcontrail(pcols,pver)
    real(r8), pointer, dimension(:,:) ::  cld    ! cloud fraction
    real(r8), pointer, dimension(:,:) :: ac_H2O
    real(r8), pointer, dimension(:,:) :: ac_SLANT_DIST
    integer  :: itim, ifld
    integer :: ixcldice, ixcldliq                  ! indices for CLDICE and CLDLIQ
    integer :: ixnumice, ixnumliq
    real(r8):: zi, zm, rog
    logical :: has_aircraft_H2O
    real(r8) :: hkl, hkk, tv
    real(r8) :: particle_mass 
    real(r8) :: ICIWC0(pcols,pver), ICIWC, rho
    real(r8) :: qs
    real(r8) :: wght(pcols)
    real(r8) :: dz, ratio(pcols,pver)
    real(r8) :: dcld(pcols,pver) 
    real(r8) :: ac_Q, ac_Q1, ac_Q2
    real(r8) :: RHcts(pcols,pver)
    integer :: itype
    logical :: lq(pcnst)

    integer :: yr, mon, day, ncsec
    real(r8) :: curr_factor

!    ICIWC = ICIWC0/rhodair

    has_aircraft_H2O = .false.

!-----------------------------------------------------------------------------------------
! check if ac_H2O in namelist
! if not, then bypass this subroutine
!-----------------------------------------------------------------------------------------
    if(aircraft_cnt>0) then
     do i=1,aircraft_cnt
      if(trim(spc_name_list(i)) == 'ac_H2O') then
        has_aircraft_H2O = .true.
      endif
     enddo
    endif
    if(.not. has_aircraft_H2O)  return


    particle_mass = 4._r8/3._r8*pi*rhoi*radius**3   ! mass of ice particle
   
    rog = rair/gravit                               ! Rd/Cp


    contrail(:,:) = 0.0_r8
    pcontrail(:,:) = 0.0_r8
    ICIWC0(:,:) = 0.0_r8
    RHcts(:,:) = 0.0_r8

    lchnk = state1%lchnk
    ncol  = state1%ncol
    
    call get_wght_all_p(lchnk, ncol, wght)

    itim = pbuf_old_tim_idx()
    ifld = pbuf_get_index('CLD')
    call pbuf_get_field(pbuf, ifld, cld, (/1,1,itim/),(/pcols,pver,1/))

    ifld = pbuf_get_index('ac_H2O')
    call pbuf_get_field(pbuf,ifld,ac_H2O)
    ifld = pbuf_get_index('ac_SLANT_DIST')
    call pbuf_get_field(pbuf,ifld,ac_SLANT_DIST)
 
	
    ! Update constituents, all schemes use time split q: no tendency kept
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    ! Check for number concentration of cloud liquid and cloud ice (if not present)
    ! the indices will be set to -1)
    call cnst_get_ind('NUMICE', ixnumice, abort=.false.)
    call cnst_get_ind('NUMLIQ', ixnumliq, abort=.false.)

!------------------------------------------------------------------------------------------
    lq(:) = .FALSE.
    lq(1) = .TRUE.
    lq(ixcldice) = .TRUE.
    lq(ixnumice) = .TRUE.

    call physics_ptend_init(ptend_loc, state1%psetcols,'ssatcontrail',ls=.true.,lq = lq)
!-----------------------------------------------------------------------------------------

! adjust h2o to volume mixing ratio (mass adjustment and conversion from g/kg to kg/kg)
    Ma = mwdry   
    Mh2o = mwh2o

! contrail paramter (G = CF*p/epi)
! and Schumann 1996, reprinted by Ponater 2002, JGR (eq 6-8)

    epsi = Mh2o/Ma
    ei = 1.21_r8      ! water vapor emiision index (g) h2o per kg fuel (Schumann 96)?
    Q = 43.e6         ! specific combustion heat Schummann 1996, Q = 43 MJ/kg
    eta = 0.3_r8      ! propulsion effieciency (Ponater 2002)
    
    ratio(i,k) = 0._r8
    dcld(i,k) = 0._r8

    do i=1,ncol
      do k=1,pver
	
        p = (state1%pint(i,k)+state1%pint(i,k+1))/2.0_r8
       
        G = (ei*cpair*p)/(epsi*Q*(1.0_r8-eta))   ! eq 7, Ponater JGR 2002

        T_contr = -46.46_r8+9.43_r8*log(G-0.053_r8)+0.72_r8*log(G-0.053_r8)*log(G-0.053_r8) ! eq 8, Ponater JGR 2002
        T_contr = T_contr + 273.15_r8  ! convert to Kelvin
       
        ! compute saturation pressure        
        itype = 0
        call gffgch(T_contr,eslTc,itype)
        itype = 0
        call gffgch(state1%t(i,k),eslT,itype)

        RH_contr = (G*(state1%t(i,k)-T_contr)+eslTc)/eslT
        ! RH_contr ranges between 0 and 1
        if(RH_contr.gt.1.0_r8) RH_contr = 1.0_r8
        if(RH_contr.lt.0.0_r8) RH_contr = 0.0_r8
        
        w = state1%q(i,k,1)/(1.0_r8-state1%q(i,k,1))  ! mixing ratio from specific humidity        
        itype = 1
        call gffgch(state1%t(i,k),esiT,itype)   ! saturation water vapor with respect to ice        
        ws = epsi*esiT/(p-esiT)  ! saturation mixing ration with respect to ice
        qs = ws/(1.0_r8+ws)

        RH = w/ws  ! relative humidity with respect to ice
        if( RH.ge.1.0_r8 ) RHcts(i,k) = 1.0_r8

! Schumann, 2002: IWC(g/m3) = exp(6.97+0.103*T(C))*1e-3
!                 IWC(kg/m3) = exp(6.97+0.103*T(C))*1e-6 

        ICIWC0(i,k) = exp(6.97_r8+0.103_r8*(state1%t(i,k)-273.16_r8)) ! in mg/m3
        rho = p/(287._r8*state1%t(i,k))
        ICIWC = ICIWC0(i,k)/rho*1.0e-6


! persistent contrail condition
        if( (state1%t(i,k).lt.T_contr).and.(RH.gt.RH_contr).and.(RH.gt.1.0_r8).and.(ac_H2O(i,k).gt.0.0_r8) ) then

! if persistent contrail, H2O emitted from aircraft turns into cloud ice
          dz = state1%zi(i,k)-state1%zi(i,k+1) 
          ratio(i,k) = (ac_SLANT_DIST(i,k)*dtime*1.e4_r8)/(dz*rearth*rearth*wght(i))

          ac_Q = min(ac_H2O(i,k)*dtime + (state1%q(i,k,1)-qs)*ratio(i,k),ratio(i,k)*ICIWC)
          ptend_loc%q(i,k,ixcldice) = ac_Q/dtime

! take out water vapor from q
          ptend_loc%q(i,k,1) = -(ac_Q-ac_H2O(i,k)*dtime)/dtime

! modify cloud fraction
! by a prescribed ICIWC, we may deduce the new cloud fraction
        
         cld(i,k) = min(1._r8, cld(i,k)+ac_Q/ICIWC)
        
! modify cloud ice number concentration, 
! by assuming the particle size, the number of ice particles may be obtained
          ptend_loc%q(i,k,ixnumice) = ac_Q/particle_mass/dtime

        else
! if not persistent contrail, just add ac_H2O to state1%q(1) (vapor phase)

          ptend_loc%q(i,k,1) = ac_H2O(i,k)
          
        endif


      enddo

! modify dry static energy if water field is added to any grid cell
! this bypasses geopotential_t which assumes dry static energy conservation
! water vapor added to the system is assumed to increase dry static energy
! conservation of dry static energy by geopotential_t will lower temperature to compensate
 
      zi = 0.0_r8
      do k=pver,1,-1
         hkl = state1%lnpint(i,k+1)-state1%lnpint(i,k)
         hkk = 1._r8 - state1%pint(i,k) * hkl * state1%rpdel(i,k)

         tv = state1%t(i,k) * (1._r8 + zvir*(state1%q(i,k,1)+ptend_loc%q(i,k,1)*dtime))

         zm = zi + rog * tv * hkk       
         zi = zi + rog * tv * hkl

         ptend_loc%s(i,k) = (cpair*state1%t(i,k)+gravit*zm + state1%phis(i) - state1%s(i,k) )/dtime
      enddo

     enddo

    call physics_update(state1, ptend_loc, dtime, tend)

    end subroutine ssatcontrail_d0
   
subroutine gffgch(t       ,es      ,itype   )
   use shr_kind_mod, only: r8 => shr_kind_r8
   use physconst,    only: tmelt
   use abortutils,   only: endrun
   use cam_logfile,  only: iulog
    
   implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: t          ! Temperature
!
! Output arguments
!
   integer, intent(inout) :: itype   ! Flag for ice phase and associated transition

   real(r8), intent(out) :: es         ! Saturation vapor pressure
!
!---------------------------Local variables-----------------------------
   real(r8) e1         ! Intermediate scratch variable for es over water
   real(r8) e2         ! Intermediate scratch variable for es over water
   real(r8) eswtr      ! Saturation vapor pressure over water
   real(r8) f          ! Intermediate scratch variable for es over water
   real(r8) f1         ! Intermediate scratch variable for es over water
   real(r8) f2         ! Intermediate scratch variable for es over water
   real(r8) f3         ! Intermediate scratch variable for es over water
   real(r8) f4         ! Intermediate scratch variable for es over water
   real(r8) f5         ! Intermediate scratch variable for es over water
   real(r8) ps         ! Reference pressure (mb)
   real(r8) t0         ! Reference temperature (freezing point of water)
   real(r8) term1      ! Intermediate scratch variable for es over ice
   real(r8) term2      ! Intermediate scratch variable for es over ice
   real(r8) term3      ! Intermediate scratch variable for es over ice
   real(r8) tr         ! Transition range for es over water to es over ice
   real(r8) ts         ! Reference temperature (boiling point of water)
   real(r8) weight     ! Intermediate scratch variable for es transition
   integer itypo   ! Intermediate scratch variable for holding itype
!
!-----------------------------------------------------------------------
!
! Check on whether there is to be a transition region for es
!
   if (itype < 0) then
      tr    = abs(real(itype,r8))
      itypo = itype
      itype = 1
   else
      tr    = 0.0_r8
      itypo = itype
   end if
   if (tr > 40.0_r8) then
      write(iulog,900) tr
      call endrun ('GFFGCH')                ! Abnormal termination
   end if
!
   if(t < (tmelt - tr) .and. itype == 1) go to 10
!
! Water
!
   ps = 1013.246_r8
   ts = 373.16_r8
   e1 = 11.344_r8*(1.0_r8 - t/ts)
   e2 = -3.49149_r8*(ts/t - 1.0_r8)
   f1 = -7.90298_r8*(ts/t - 1.0_r8)
   f2 = 5.02808_r8*log10(ts/t)
   f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
   f4 = 8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
   f5 = log10(ps)
   f  = f1 + f2 + f3 + f4 + f5
   es = (10.0_r8**f)*100.0_r8
   eswtr = es
!
   if(t >= tmelt .or. itype == 0) go to 20
!
! Ice
!
10 continue
   t0    = tmelt
   term1 = 2.01889049_r8/(t0/t)
   term2 = 3.56654_r8*log(t0/t)
   term3 = 20.947031_r8*(t0/t)
   es    = 575.185606e10_r8*exp(-(term1 + term2 + term3))
!
   if (t < (tmelt - tr)) go to 20
!
! Weighted transition between water and ice
!
   weight = min((tmelt - t)/tr,1.0_r8)
   es = weight*es + (1.0_r8 - weight)*eswtr
!
20 continue
   itype = itypo
   return
!
900 format('GFFGCH: FATAL ERROR ******************************',/, &
           'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR', &
           ' PRESSURE, TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF', &
           ' 40.0 DEGREES C',/, ' TR = ',f7.2)
!
end subroutine gffgch


end module ssatcontrail
