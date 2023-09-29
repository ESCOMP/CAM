module ssatcontrail
! contrail parameterization
! see Chen et al., 2012: Global contrail coverage simulated
!     by CAM5 with the inventory of 2006 global aircraft emissions, JAMES
!     https://doi.org/10.1029/2011MS000105
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use ppgrid,         only: pcols, pver
    use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
    use physconst,      only: cpair,mwdry,mwh2o, gravit, zvir, rair, pi, rearth, tmelt
    use physics_buffer, only: pbuf_get_index, pbuf_get_field, physics_buffer_desc,  pbuf_old_tim_idx
    use constituents,   only: cnst_get_ind, pcnst
    use phys_grid,      only: get_wght_all_p
    use wv_saturation,  only: qsat_water, qsat_ice
    use aircraft_emit,  only: get_aircraft

    implicit none
    private
    save

    public ssatcontrail_d0

    ! Private data
    real(r8), parameter :: rhoi = 500.0_r8             ! density of ice (500 kg/m3)
    real(r8), parameter :: radius = 3.75e-6_r8         ! diameter of ice particle = 7.5 microns

 
contains

    subroutine ssatcontrail_d0(state1,pbuf,dtime,ptend_loc)
    implicit none

    type(physics_state), intent(in) :: state1
    type(physics_ptend), intent(inout) :: ptend_loc
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8),            intent(in)    :: dtime      ! time step
!------------------------Local storage------------------------------------------------------
    real(r8) :: Ma, Mh2o, epsi, Q, eta, p, G, T_contr, eslTc, eslT, RH_contr, qslTc, qslT
    real(r8) :: w, esiT, qsiT, ws, RH, ei
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
    logical :: has_aircraft_distance
    real(r8) :: hkl, hkk, tv
    real(r8) :: particle_mass
    real(r8) :: ICIWC0(pcols,pver), ICIWC, rho
    real(r8) :: qs
    real(r8) :: wght(pcols)
    real(r8) :: dz, ratio(pcols,pver)
    real(r8) :: dcld(pcols,pver)
    real(r8) :: ac_Q, ac_Q1, ac_Q2
    real(r8) :: RHcts(pcols,pver)
    logical :: lq(pcnst)

    integer :: yr, mon, day, ncsec
    real(r8) :: curr_factor
    integer :: aircraft_cnt
    character(len=16) :: spc_name_list(30)
 
    has_aircraft_H2O = .false.
    has_aircraft_distance = .false.

    ! Update constituents, all schemes use time split q: no tendency kept
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    ! Check for number concentration of cloud liquid and cloud ice (if not present)
    ! the indices will be set to -1)
    call cnst_get_ind('NUMICE', ixnumice, abort=.false.)
    call cnst_get_ind('NUMLIQ', ixnumliq, abort=.false.)

    call get_aircraft(aircraft_cnt, spc_name_list)
!-----------------------------------------------------------------------------------------
! check if ac_H2O in namelist
! if not, then bypass this subroutine
!-----------------------------------------------------------------------------------------
    if(aircraft_cnt>0) then
     do i=1,aircraft_cnt
      if(trim(spc_name_list(i)) == 'ac_H2O') then
        has_aircraft_H2O = .true.
      endif
      if(trim(spc_name_list(i)) == 'ac_SLANT_DIST' .or. trim(spc_name_list(i)) == 'ac_TRACK_DIST') then
        has_aircraft_distance = .true.
      endif
     enddo
    endif
!------------------------------------------------------------------------------------------
    lq(:) = .FALSE.
    lq(1) = .TRUE.
    lq(ixcldice) = .TRUE.
    lq(ixnumice) = .TRUE.

    call physics_ptend_init(ptend_loc, state1%psetcols,'ssatcontrail',ls=.true.,lq = lq)
!-----------------------------------------------------------------------------------------
    if(.not. has_aircraft_H2O)  return
    if(.not. has_aircraft_distance) return

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
 
! adjust h2o to volume mixing ratio (mass adjustment and conversion from g/kg to kg/kg)
    Ma = mwdry
    Mh2o = mwh2o

! contrail paramter (G = CF*p/epi)
! and Schumann 1996 DOI: 10.1127/metz/5/1996/4, reprinted by Ponater 2002, JGR (eq 6-8) DOI: 10.1029/2011MS000105

    epsi = Mh2o/Ma
    ei = 1.21_r8      ! water vapor emision index (g) h2o per kg fuel (Schumann 96)
    Q = 43.e6_r8      ! specific combustion heat Schummann 1996, Q = 43 MJ/kg
    eta = 0.3_r8      ! propulsion effieciency (Ponater 2002)
 
    ratio = 0._r8
    dcld = 0._r8

    do i=1,ncol
      do k=1,pver
        p = (state1%pint(i,k)+state1%pint(i,k+1))/2.0_r8
 
        G = (ei*cpair*p)/(epsi*Q*(1.0_r8-eta))   ! eq 7, Ponater JGR 2002

        if( G > 0.053_r8 ) then
            T_contr = -46.46_r8+9.43_r8*log(G-0.053_r8)+0.72_r8*log(G-0.053_r8)*log(G-0.053_r8) ! eq 6, Ponater JGR 2002
            T_contr = T_contr + tmelt  ! convert to Kelvin
 
            ! compute saturation pressure
            call qsat_water(T_contr, p, eslTc, qslTc)
            call qsat_water(state1%t(i,k), p, eslT, qslT)

            RH_contr = (G*(state1%t(i,k)-T_contr)+eslTc)/eslT
            ! RH_contr ranges between 0 and 1
            if(RH_contr>1.0_r8) RH_contr = 1.0_r8
            if(RH_contr<0.0_r8) RH_contr = 0.0_r8
 
            w = state1%q(i,k,1)/(1.0_r8-state1%q(i,k,1))  ! mixing ratio from specific humidity
            call qsat_ice(state1%t(i,k), p, esiT, qsiT)
            ws = epsi*esiT/(p-esiT)  ! saturation mixing ration with respect to ice
            qs = ws/(1.0_r8+ws)

            RH = w/ws  ! relative humidity with respect to ice
            if( RH>=1.0_r8 ) RHcts(i,k) = 1.0_r8

! Schumann, U. “Contrail Cirrus.” In Cirrus, edited by D. K. Lynch and others, 231–55. Oxford University Press, 2002
!                 IWC(g/m3) = exp(6.97+0.103*T(C))*1e-3
!                 IWC(kg/m3) = exp(6.97+0.103*T(C))*1e-6

            ICIWC0(i,k) = exp(6.97_r8+0.103_r8*(state1%t(i,k)-tmelt)) ! in mg/m3
            rho = p/(rair*state1%t(i,k))
            ICIWC = ICIWC0(i,k)/rho*1.0e-6_r8


! persistent contrail condition
            if( (state1%t(i,k)<T_contr).and.(RH>RH_contr).and.(RH>1.0_r8).and.(ac_H2O(i,k)>0.0_r8) ) then

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
       
       else
           ptend_loc%q(i,k,1) = ac_H2O(i,k)
       end if

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

    end subroutine ssatcontrail_d0

end module ssatcontrail
