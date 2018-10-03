module dpie_coupling
!
! Dynamics/Physics Ionosphere/Electrodynamics coupler.
! B. Foster, 2015.
!
  use shr_kind_mod    ,only: r8 => shr_kind_r8
  use cam_logfile     ,only: iulog
  use cam_history     ,only: outfld
  use cam_history     ,only: addfld, horiz_only
  use cam_history_support, only: fillvalue
  use cam_abortutils  ,only: endrun
  use spmd_utils      ,only: masterproc
  use savefield_waccm ,only: savefld_waccm
  use edyn_mpi        ,only: array_ptr_type
  use perf_mod        ,only: t_startf, t_stopf
  use amie_module     ,only: getamie
  use edyn_solve      ,only: phihm
  use edyn_params     ,only: dtr, rtd
  use edyn_mpi,        only: switch_model_format
  use aurora_params,   only: amie_period ! turns on overwrite of energy fields in aurora phys
  
  implicit none

  private
  public :: d_pie_init
  public :: d_pie_epotent  ! sets electric potential
  public :: d_pie_coupling ! handles coupling with edynamo and ion transport

  logical :: ionos_edyn_active, ionos_oplus_xport    ! if true, call oplus_xport for O+ transport
  integer :: nspltop  ! nsplit for oplus_xport

  logical :: debug = .false.

  real(r8) :: crad(2), crit(2)
  logical :: crit_user_set = .false.
  real(r8), parameter :: amie_default_crit(2) = (/ 35._r8, 40._r8 /)
  
contains
!----------------------------------------------------------------------
  subroutine d_pie_init( edyn_active_in, oplus_xport_in, oplus_nsplit_in, crit_colats_deg )

    logical, intent(in) :: edyn_active_in, oplus_xport_in
    integer, intent(in) :: oplus_nsplit_in
    real(r8),intent(in) :: crit_colats_deg(:)

    ionos_edyn_active = edyn_active_in
    ionos_oplus_xport = oplus_xport_in
    nspltop = oplus_nsplit_in

    crit_user_set = all( crit_colats_deg(:) > 0._r8 )
    if (crit_user_set) then
       crit(:) = crit_colats_deg(:)*dtr
    endif
    
    ! Dynamo inputs (called from dpie_coupling. Fields are in waccm format, in CGS units):
    call addfld ('DPIE_OMEGA',(/ 'lev' /), 'I', 'Pa/s    ','OMEGA input to DPIE coupling', gridname='fv_centers')
    call addfld ('DPIE_MBAR' ,(/ 'lev' /), 'I', '        ','MBAR Mean Mass from dpie_coupling', gridname='fv_centers')
    call addfld ('DPIE_TN   ',(/ 'lev' /), 'I', 'deg K   ','DPIE_TN'   , gridname='fv_centers')
    call addfld ('DPIE_UN   ',(/ 'lev' /), 'I', 'cm/s    ','DPIE_UN'   , gridname='fv_centers')
    call addfld ('DPIE_VN   ',(/ 'lev' /), 'I', 'cm/s    ','DPIE_VN'   , gridname='fv_centers')
    call addfld ('DPIE_WN   ',(/ 'lev' /), 'I', 'cm/s    ','DPIE_WN'   , gridname='fv_centers')
    call addfld ('DPIE_OM   ',(/ 'lev' /), 'I', 's-1     ','DPIE_OM'   , gridname='fv_centers')
    call addfld ('DPIE_ZHT  ',(/ 'lev' /), 'I', 'cm      ','DPIE_ZHT (geometric height,simple)', gridname='fv_centers')
    call addfld ('DPIE_ZGI  ',(/ 'lev' /), 'I', 'cm      ','DPIE_ZGI (geopotential height on interfaces)', gridname='fv_centers')
    call addfld ('DPIE_BARM ',(/ 'lev' /), 'I', '        ','DPIE_BARM' , gridname='fv_centers')
    call addfld ('DPIE_O2   ',(/ 'lev' /), 'I', 'mmr     ','DPIE_O2'   , gridname='fv_centers')
    call addfld ('DPIE_O    ',(/ 'lev' /), 'I', 'mmr     ','DPIE_O'    , gridname='fv_centers')
    call addfld ('DPIE_N2   ',(/ 'lev' /), 'I', 'mmr     ','DPIE_N2'   , gridname='fv_centers')
    call addfld ('DPIE_TE   ',(/ 'lev' /), 'I', 'deg K   ','DPIE_TE'   , gridname='fv_centers')
    call addfld ('DPIE_TI   ',(/ 'lev' /), 'I', 'deg K   ','DPIE_TI'   , gridname='fv_centers')

    call addfld ('DPIE_OPMMR' ,(/ 'lev' /), 'I', 'mmr'  ,'DPIE_OPMMR'  , gridname='fv_centers')
    call addfld ('DPIE_O2P',(/ 'lev' /), 'I', 'm^3','DPIE_O2P(dpie input)', gridname='fv_centers')
    call addfld ('DPIE_NOP',(/ 'lev' /), 'I', 'm^3','DPIE_NOP(dpie input)', gridname='fv_centers')
    call addfld ('DPIE_N2P',(/ 'lev' /), 'I', 'm^3','DPIE_N2P(dpie input)', gridname='fv_centers')

    call addfld ('OPLUS',   (/ 'lev' /), 'I', 'cm^3','O+ (oplus_xport output)', gridname='fv_centers')
    call addfld ('WACCM_UI'   ,(/ 'lev' /), 'I', 'm/s'  ,'WACCM_UI (dpie output)', gridname='fv_centers')
    call addfld ('WACCM_VI'   ,(/ 'lev' /), 'I', 'm/s'  ,'WACCM_VI (dpie output)', gridname='fv_centers')
    call addfld ('WACCM_WI'   ,(/ 'lev' /), 'I', 'm/s'  ,'WACCM_WI (dpie output)', gridname='fv_centers')
 
    call addfld ('HMF2'       , horiz_only , 'I', 'km'  ,'Height of the F2 Layer'      , gridname='fv_centers')
    call addfld ('NMF2'       , horiz_only , 'I', 'cm-3','Peak Density of the F2 Layer', gridname='fv_centers')

    call addfld ('Z3GM'   ,(/ 'lev' /), 'I', 'm'   ,'Geometric height'                        , gridname='fv_centers')
    call addfld ('Z3GMI  ',(/ 'lev' /), 'I', 'm' ,'Geometric height (Interfaces)', gridname='fv_centers')
    call addfld ('OpDens' ,(/ 'lev' /), 'I', 'cm^3','O+ Number Density'                       , gridname='fv_centers')
    call addfld ('EDens'  ,(/ 'lev' /), 'I', 'cm^3','e Number Density (sum of O2+,NO+,N2+,O+)', gridname='fv_centers')

  end subroutine d_pie_init

!-----------------------------------------------------------------------
  subroutine d_pie_epotent( highlat_potential_model, crit_out, i0,i1,j0,j1, efxg, kevg )
    use edyn_solve,  only: pfrac    ! NH fraction of potential (nmlonp1,nmlat0)
    use edyn_geogrid,only: nglblat=>nlat
    use time_manager,only: get_curr_date
    use heelis,      only: heelis_model
    use wei05sc,     only: weimer05  ! driver for weimer high-lat convection model
    use edyn_esmf,   only: edyn_esmf_update
    use solar_parms_data,only: solar_parms_advance
    use solar_wind_data, only: solar_wind_advance
    use solar_wind_data, only: bzimf=>solar_wind_bzimf, byimf=>solar_wind_byimf
    use solar_wind_data, only: swvel=>solar_wind_swvel, swden=>solar_wind_swden
    use edyn_mpi,        only: mlat0,mlat1,mlon0,omlon1
    use edyn_maggrid ,only: nmlonp1,nmlat
! Args:
!
    character(len=*),intent(in) :: highlat_potential_model
    real(r8),        intent(out):: crit_out(2) ! critical colatitudes (degrees)
    integer,optional,intent(in) :: &
      i0,                 & ! grid%ifirstxy
      i1,                 & ! grid%ilastxy
      j0,                 & ! grid%jfirstxy
      j1                    ! grid%jlastxy
    real(r8),optional,intent(out) :: efxg(:,:) ! energy flux from AMIE
    real(r8),optional,intent(out) :: kevg(:,:) ! characteristic mean energy from AMIE
    !
    ! local vars
    !

    real(r8) :: secs           ! time of day in seconds
    integer :: iyear,imo,iday,tod  ! tod is time-of-day in seconds
    real(r8) :: sunlons(nglblat)

    integer  :: iprint,amie_ibkg
    integer :: i, j, iamie
    type(array_ptr_type) :: ptrs(2)
    !
    ! AMIE fields (extra dimension added for longitude switch)
    !
    real(r8) :: amie_efxm(nmlonp1,nmlat), amie_kevm(nmlonp1,nmlat)  ! auroral energy flux and 
    real(r8) :: amie_phihm(nmlonp1,nmlat)
    real(r8),allocatable,target :: amie_efxg (:,:,:) ! AMIE energy flux
    real(r8),allocatable,target :: amie_kevg (:,:,:) ! AMIE characteristic mean energy
    
    call edyn_esmf_update()

    call get_curr_date(iyear,imo,iday,tod) ! tod is integer time-of-day in seconds
    secs = tod                             ! should promote from int to real(r8)

    ! update solar wind data (IMF, etc.)
    call solar_wind_advance()

    ! update kp -- phys timestep init happens later ...
    call solar_parms_advance()


    !
    ! Get sun's longitude at latitudes (geographic):
    !
    call sunloc(iday,secs,sunlons) ! sunlons(nglblat) is returned
    !
    ! Get high-latitude convection from empirical model (heelis or weimer).
    ! High-latitude potential phihm (edyn_solve) is defined for edynamo.
    !
    if (trim(highlat_potential_model) == 'heelis') then
       call heelis_model(sunlons) ! heelis.F90
    elseif (trim(highlat_potential_model) == 'weimer') then
       !
       call weimer05(byimf,bzimf,swvel,swden,sunlons)
       if (debug.and.masterproc) then
          write(iulog,"('dpie_coupling call weimer05: byimf,bzimf=',2f8.2,' swvel,swden=',2f8.2)") &
               byimf,bzimf,swvel,swden
       endif
    else
       call endrun('dpie_coupling: Unknown highlat_potential_model') 
    endif

    if (present(efxg)) then ! the presence of efxg indicate the user want to use prescribed potential
       iprint = 0
       amie_ibkg = 0
       iamie = 1
       if (masterproc) write(iulog,"('Calling getamie >>> iamie=',i2)") iamie
       
       call getamie(iyear,imo,iday,tod,sunlons(1),amie_ibkg,iprint,iamie, &
            amie_phihm,amie_efxm,amie_kevm,crad,efxg,kevg)

       if (masterproc) write(iulog,"('After Calling getamie >>> iamie=',i2)") iamie
       amie_period = iamie == 1

       do j=mlat0,mlat1
          call outfld('amie_phihm',amie_phihm(mlon0:omlon1,j),omlon1-mlon0+1,j)
          call outfld('amie_efxm',amie_efxm(mlon0:omlon1,j),omlon1-mlon0+1,j)
          call outfld('amie_kevm',amie_kevm(mlon0:omlon1,j),omlon1-mlon0+1,j)
       enddo

       if (amie_period) then

          phihm = amie_phihm

          ! Load AMIE fields into pointers for TIE-GCM to WACCM longitude swap
          !
          allocate(amie_efxg(1,i0:i1,j0:j1))
          allocate(amie_kevg(1,i0:i1,j0:j1))

          do i=i0,i1
             do j=j0,j1
                amie_efxg(1,i,j) = efxg(i-i0+1,j-j0+1)
                amie_kevg(1,i,j) = kevg(i-i0+1,j-j0+1)
             enddo
          enddo

          ptrs(1)%ptr => amie_efxg
          ptrs(2)%ptr => amie_kevg
          call switch_model_format(ptrs,1,1,i0,i1,j0,j1, 2)

          do i=i0,i1
             do j=j0,j1
                efxg(i-i0+1,j-j0+1) = amie_efxg(1,i,j)
                kevg(i-i0+1,j-j0+1) = amie_kevg(1,i,j)
             enddo
          enddo

          deallocate(amie_efxg)
          deallocate(amie_kevg)

       endif

       call savefld_waccm(efxg,'amie_efxg',1,i0,i1,j0,j1) 
       call savefld_waccm(kevg,'amie_kevg',1,i0,i1,j0,j1) 

    endif
 
    call calc_pfrac(sunlons(1),pfrac) ! returns pfrac for dynamo (edyn_solve)

    crit_out(:) = crit(:)*rtd ! degrees
  end subroutine d_pie_epotent

!-----------------------------------------------------------------------
  subroutine d_pie_coupling(omega,pe,zgi,zgpmid,u,v,tn,                 &
    sigma_ped,sigma_hall,te,ti,o2mmr,o1mmr,h1mmr,o2pmmr,                &
    nopmmr,n2pmmr,opmmr,opmmrtm1,ui,vi,wi,                              &
    rmassO2,rmassO1,rmassH,rmassN2,rmassO2p, rmassNOp,rmassN2p,rmassOp, &
    i0,i1,j0,j1)
!
! Call dynamo to calculate electric potential, electric field, and ion drifts.
! Then call oplus_xport to transport O+, which is passed back to physics.
!
! This routine is called from p_d_coupling (dynamics/fv/dp_coupling.F90) when
! nstep > 0.
!
    use edyn_geogrid, only: nlev, nilev
    use shr_const_mod,only:      &
      grav   => shr_const_g,     &   ! gravitational constant (m/s^2)
      kboltz => shr_const_boltz      ! Boltzmann constant (J/K/molecule)
    use time_manager, only: get_nstep
    use time_manager, only: get_curr_date
    use edynamo,      only: dynamo
    use edyn_mpi,      only: mp_geo_halos,mp_pole_halos
    use oplus,         only: oplus_xport
    use ref_pres,      only: pref_mid
!
! Args:
!
    integer,intent(in) :: &
      i0,                 & ! grid%ifirstxy
      i1,                 & ! grid%ilastxy
      j0,                 & ! grid%jfirstxy
      j1                    ! grid%jlastxy

    real(r8),intent(in) :: omega  (i0:i1,j0:j1,nlev)    ! pressure velocity on midpoints (Pa/s) (i,k,j)
    real(r8),intent(in) :: pe     (i0:i1,nilev,j0:j1)   ! interface pressure (Pa)  (note i,k,j dims)
    real(r8),intent(in) :: zgi    (i0:i1,j0:j1,nlev)    ! geopotential height (on interfaces) (m)
    real(r8),intent(in) :: zgpmid (i0:i1,j0:j1,nlev)    ! geopotential height (on midpoints) (m)
    real(r8),intent(in) :: u      (i0:i1,j0:j1,nlev)    ! U-wind (m/s)
    real(r8),intent(in) :: v      (i0:i1,j0:j1,nlev)    ! V-wind (m/s)
    real(r8),intent(in) :: tn     (i0:i1,j0:j1,nlev)    ! neutral temperature (K)
    real(r8),intent(in) :: sigma_ped (i0:i1,j0:j1,nlev) ! Pedersen conductivity
    real(r8),intent(in) :: sigma_hall(i0:i1,j0:j1,nlev) ! Hall conductivity
    real(r8),intent(in) :: te(i0:i1,j0:j1,nlev)         ! electron temperature
    real(r8),intent(in) :: ti(i0:i1,j0:j1,nlev)         ! ion temperature
    real(r8),intent(in) :: o2mmr(i0:i1,j0:j1,nlev)      ! O2 mass mixing ratio (for oplus)
    real(r8),intent(in) :: o1mmr(i0:i1,j0:j1,nlev)      ! O mass mixing ratio (for oplus)
    real(r8),intent(in) :: h1mmr(i0:i1,j0:j1,nlev)      ! H mass mixing ratio (for oplus)
    real(r8),intent(in) :: o2pmmr(i0:i1,j0:j1,nlev)     ! O2+ mass mixing ratio (for oplus)
    real(r8),intent(in) :: nopmmr(i0:i1,j0:j1,nlev)     ! NO+ mass mixing ratio (for oplus)
    real(r8),intent(in) :: n2pmmr(i0:i1,j0:j1,nlev)     ! N2+ mass mixing ratio (for oplus)
    real(r8),intent(inout) :: opmmr(i0:i1,j0:j1,nlev)   ! O+ mass mixing ratio (oplus_xport output)
    real(r8),intent(inout) :: opmmrtm1(i0:i1,j0:j1,nlev)   ! O+ previous time step (oplus_xport output)
    real(r8),intent(inout) :: ui(i0:i1,j0:j1,nlev)      ! zonal ion drift (edynamo or empirical)
    real(r8),intent(inout) :: vi(i0:i1,j0:j1,nlev)      ! meridional ion drift (edynamo or empirical)
    real(r8),intent(inout) :: wi(i0:i1,j0:j1,nlev)      ! vertical ion drift (edynamo or empirical)
    real(r8),intent(in) :: rmassO2                      ! O2 molecular weight kg/kmol
    real(r8),intent(in) :: rmassO1                      ! O atomic weight kg/kmol
    real(r8),intent(in) :: rmassH                       ! H atomic weight kg/kmol
    real(r8),intent(in) :: rmassN2                      ! N2 molecular weight kg/kmol
    real(r8),intent(in) :: rmassO2p                     ! O2+ molecular weight kg/kmol
    real(r8),intent(in) :: rmassNOp                     ! NO+ molecular weight kg/kmol
    real(r8),intent(in) :: rmassN2p                     ! N2+ molecular weight kg/kmol
    real(r8),intent(in) :: rmassOp                      ! O+ molecular weight kg/kmol
!
! Local:
!
    integer :: i,j,k
    integer :: kx                  ! Vertical index at peak of F2 layer electron density
    integer :: nstep
    integer :: nfields             ! Number of fields for multi-field calls
    integer :: iyear,imo,iday,tod  ! tod is time-of-day in seconds
    integer :: isplit              ! loop index

    real(r8) :: secs           ! time of day in seconds

    real(r8), parameter :: n2min = 1.e-6_r8  ! lower limit of N2 mixing ratios
    real(r8), parameter :: small = 1.e-25_r8 ! for fields not currently available
    real(r8) :: zht  (i0:i1,j0:j1,nlev) ! geometric height (m) (Simple method - interfaces)
    real(r8) :: zhtmid(i0:i1,j0:j1,nlev)! geometric height (m) (Simple method - midpoints)
    real(r8) :: wn   (i0:i1,j0:j1,nlev) ! vertical velocity (from omega)
    real(r8) :: mbar (i0:i1,j0:j1,nlev) ! mean molecular weight
    real(r8) :: n2mmr(i0:i1,j0:j1,nlev) ! N2 mass mixing ratio (for oplus)
    real(r8) :: pmid_inv(nlev)          ! inverted reference pressure at midpoints (Pa)
    real(r8) :: pmid(i0:i1,nlev,j0:j1)  ! pressure at midpoints (Pa) (global i,j)
    real(r8) :: re = 6.370e6_r8         ! earth radius (m)

    real(r8),dimension(i0:i1,j0:j1,nlev) :: & ! ion number densities (m^3)
      o2p,nop,n2p,op,ne, optm1

    real(r8),dimension(nlev,i0:i1,j0:j1) :: opmmr_kij
!
! Args for dynamo:
    real(r8),target :: edyn_tn   (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_un   (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_vn   (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_wn   (nlev,i0:i1,j0:j1) ! vertical wind (cm/s)
    real(r8),target :: edyn_zht  (nlev,i0:i1,j0:j1) ! geometric height (cm)
    real(r8),target :: edyn_mbar (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_ped  (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_hall (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_ui   (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_vi   (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_wi   (nlev,i0:i1,j0:j1)
!
! Additional fields needed by oplus_xport:
    real(r8),target :: edyn_te   (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_ti   (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_o2   (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_o1   (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_n2   (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_op   (nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_optm1(nlev,i0:i1,j0:j1)
    real(r8),target :: edyn_om   (nlev,i0:i1,j0:j1) ! omega vertical motion (1/s)
    real(r8),target :: edyn_zgi  (nlev,i0:i1,j0:j1) ! geopotential height (cm) (interfaces)
    real(r8),target :: op_out    (nlev,i0:i1,j0:j1) ! oplus_xport output
    real(r8),target :: opnm_out  (nlev,i0:i1,j0:j1) ! oplus_xport output at time n-1
    real(r8),target :: edyn_ne   (nlev,i0:i1,j0:j1) ! electron density diagnostic

    real(r8),target :: halo_tn  (nlev,i0-2:i1+2,j0-2:j1+2) ! neutral temperature (deg K)
    real(r8),target :: halo_te  (nlev,i0-2:i1+2,j0-2:j1+2) ! electron temperature (deg K)
    real(r8),target :: halo_ti  (nlev,i0-2:i1+2,j0-2:j1+2) ! ion temperature (deg K)
    real(r8),target :: halo_un  (nlev,i0-2:i1+2,j0-2:j1+2) ! neutral zonal wind (cm/s)
    real(r8),target :: halo_vn  (nlev,i0-2:i1+2,j0-2:j1+2) ! neutral meridional wind (cm/s)
    real(r8),target :: halo_om  (nlev,i0-2:i1+2,j0-2:j1+2) ! omega (1/s)
    real(r8),target :: halo_o2  (nlev,i0-2:i1+2,j0-2:j1+2) ! o2 (mmr)
    real(r8),target :: halo_o1  (nlev,i0-2:i1+2,j0-2:j1+2) ! o (mmr)
    real(r8),target :: halo_n2  (nlev,i0-2:i1+2,j0-2:j1+2) ! n2 (mmr)
    real(r8),target :: halo_mbar(nlev,i0-2:i1+2,j0-2:j1+2) ! mean molecular weight
    real(r8), allocatable :: polesign(:)
!
    real(r8) :: nmf2  (i0:i1,j0:j1) ! Electron number density at F2 peak (m-3 converted to cm-3)
    real(r8) :: hmf2  (i0:i1,j0:j1) ! Height of electron number density F2 peak (m converted to km)
    real(r8) :: &
      height(3),    &  ! Surrounding heights when locating electron density F2 peak
      nde(3)           ! Surround densities when locating electron density F2 peak
    real(r8) h12,h22,h32,deltx,atx,ax,btx,bx,ctx,cx ! Variables used for weighting when locating F2 peak
! 
    logical :: do_integrals
!
! Pointers for multiple-field calls:
    type(array_ptr_type),allocatable :: ptrs(:)

    call t_startf('d_pie_coupling')

    if (debug.and.masterproc) then

       nstep = get_nstep()
       call get_curr_date(iyear,imo,iday,tod) ! tod is integer time-of-day in seconds
       secs = tod ! integer to float

       write(iulog,"('Enter d_pie_coupling: nstep=',i8,' iyear,imo,iday=',3i5,' ut (hrs)=',f6.2)") &
            nstep,iyear,imo,iday,secs/3600._r8

       write(iulog,"('d_pie_coupling: nspltop = ',i3)") nspltop
    endif
!
! Get pressure at midpoints from pe (note pe is vertical dimension is nilev):
!
    do k=1,nlev
      pmid(i0:i1,k,j0:j1) = 0.5_r8*(pe(i0:i1,k,j0:j1)+pe(i0:i1,k+1,j0:j1))
    enddo
   
    !---------------------------------------------------------------
    ! Convert geopotential z to geometric height zht (m):
    !---------------------------------------------------------------

    zht(:,:,1:nlev) = zgi(:,:,1:nlev) * (1._r8 + zgi(:,:,1:nlev) / re) ! geometric height (interfaces)

    !---------------------------------------------------------------
    ! Need geometric height on midpoints for output
    !---------------------------------------------------------------

    zhtmid(:,:,1:nlev) = zgpmid(:,:,1:nlev) *(1._r8 + zgpmid(:,:,1:nlev) / re ) 
    
    !------------------------------------------------------------------------------------------
    ! Convert virtual potential temperature to temperature and compute mean molecular weight:
    !------------------------------------------------------------------------------------------
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          n2mmr(i,j,k) = max(1.0_r8-(o1mmr(i,j,k)+o2mmr(i,j,k)+h1mmr(i,j,k)),n2min)
          mbar(i,j,k) = 1.0_r8/(o1mmr(i,j,k)/rmassO1+o2mmr(i,j,k)/rmassO2 &
                               +h1mmr(i,j,k)/rmassH+n2mmr(i,j,k)/rmassN2)
        enddo
      enddo
    enddo
   
    !-----------------------------------------------------------------------------------------------
    ! Save analytically derived geometric height on interfaces and midpoints, omega (Pa/s) and mbar.
    !-----------------------------------------------------------------------------------------------
    do j=j0,j1
      call outfld('Z3GMI'    ,zht(i0:i1,j,1:nlev),i1-i0+1,j)
      call outfld('Z3GM'     ,zhtmid(i0:i1,j,1:nlev),i1-i0+1,j)
      call outfld('DPIE_OMEGA',omega(i0:i1,j,1:nlev),i1-i0+1,j)
      call outfld('DPIE_MBAR' ,mbar (i0:i1,j,1:nlev),i1-i0+1,j)
    enddo

    !---------------------------------------------------------------
    ! Calculate vertical neutral wind velocity wn(i,j,k).
    ! (omega is input Pa/s, grav is m/s^2, tn and mbar are calculated above)
    !---------------------------------------------------------------
    call calc_wn(tn,omega,pmid,mbar,grav,wn,i0,i1,j0,j1,nlev) ! wn is output (m/s)

    !---------------------------------------------------------------
    ! Convert from mmr to number densities (m^3):
    !---------------------------------------------------------------
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
! O2+, NO+, N2+, O+:
          o2p(i,j,k) = o2pmmr(i,j,k) * mbar(i,j,k) / rmassO2p * &
            pmid(i,k,j) / (kboltz * tn(i,j,k))
          nop(i,j,k) = nopmmr(i,j,k) * mbar(i,j,k) / rmassNOp * &
            pmid(i,k,j) / (kboltz * tn(i,j,k))
          n2p(i,j,k) = n2pmmr(i,j,k) * mbar(i,j,k) / rmassN2p * &
            pmid(i,k,j) / (kboltz * tn(i,j,k))
          op(i,j,k)  = opmmr(i,j,k)  * mbar(i,j,k) / rmassOp  * &
            pmid(i,k,j) / (kboltz * tn(i,j,k))
          optm1(i,j,k)  = opmmrtm1(i,j,k)  * mbar(i,j,k) / rmassOp  * &
            pmid(i,k,j) / (kboltz * tn(i,j,k))
        enddo
      enddo
    enddo ! k=1,nlev
!
! Save input ions to waccm history (m^3):
    do j=j0,j1
      call outfld('DPIE_O2P',o2p(i0:i1,j,1:nlev),i1-i0+1,j)
      call outfld('DPIE_NOP',nop(i0:i1,j,1:nlev),i1-i0+1,j)
      call outfld('DPIE_N2P',n2p(i0:i1,j,1:nlev),i1-i0+1,j)
      call outfld('OpDens'  ,op (i0:i1,j,1:nlev)/1.E6_r8,i1-i0+1,j)
      do k=1,nlev
        do i=i0,i1
          ne(i,j,k) = o2p(i,j,k)+nop(i,j,k)+n2p(i,j,k)+op(i,j,k)
        enddo
      enddo
      call outfld('EDens'   ,ne (i0:i1,j,1:nlev)/1.E6_r8,i1-i0+1,j)
    enddo ! j=j0,j1

    !-------------------------------------------------------------------------------
    !  Derive diagnostics nmF2 and hmF2 for output based on TIE-GCM algorithm
    !------------------------------------------------------------------------------- 
    jloop: do j=j0,j1
      iloop: do i=i0,i1

        kx = 0
        kloop: do k=2,nlev 
          if (ne(i,j,k) >= ne(i,j,k-1) .and. ne(i,j,k) >= ne(i,j,k+1)) then
            kx = k
            exit kloop
          endif
        enddo kloop

        if (kx==0) then
          hmf2(i,j) = fillvalue
          nmf2(i,j) = fillvalue 
          exit iloop
        endif

        height = (/zht(i,j,kx+1),zht(i,j,kx),zht(i,j,kx-1)/)
        nde = (/ne(i,j,kx+1),ne(i,j,kx),ne(i,j,kx-1)/)

        h12 = height(1)*height(1)
        h22 = height(2)*height(2)
        h32 = height(3)*height(3)

        deltx=h12*height(2)+h22*height(3)+h32*height(1)-h32*height(2)-h12*height(3)-h22*height(1)
        atx=nde(1)*height(2)+nde(2)*height(3)+nde(3)*height(1)-height(2)*nde(3)-height(3)*nde(1)-height(1)*nde(2)
        ax=atx/deltx

        btx=h12*nde(2)+h22*nde(3)+h32*nde(1)-h32*nde(2)-h12*nde(3)-h22*nde(1)
        bx=btx/deltx
        ctx=h12*height(2)*nde(3)+h22*height(3)*nde(1)+h32*height(1)*nde(2)-h32*height(2)*nde(1)- &
          h12*height(3)*nde(2)-h22*height(1)*nde(3)
        cx=ctx/deltx

        hmf2(i,j)=-(bx/(2._r8*ax)) * 1.E-03_r8
        nmf2(i,j)=-((bx*bx-4._r8*ax*cx)/(4._r8*ax)) * 1.E-06_r8

      enddo iloop ! i=i0,i1

      call outfld('HMF2',hmf2(i0:i1,j),i1-i0+1,j)
      call outfld('NMF2',nmf2(i0:i1,j),i1-i0+1,j)

    enddo jloop
!
! Save fields to waccm history:
! (must be transformed from (i,j,k) to (k,i,j))
!
    do j=j0,j1
      do i=i0,i1
        opmmr_kij(:,i,j) = opmmr(i,j,:)
      enddo
    enddo
   call savefld_waccm(opmmr_kij,'DPIE_OPMMR',nlev,i0,i1,j0,j1) ! mmr
!
! Prepare inputs to edynamo and oplus_xport:
!
    do k = 1,nlev
      edyn_tn   (k,i0:i1,j0:j1) = tn        (i0:i1,j0:j1,k)
      edyn_un   (k,i0:i1,j0:j1) = u         (i0:i1,j0:j1,k) * 100._r8 ! m/s -> cm/s
      edyn_vn   (k,i0:i1,j0:j1) = v         (i0:i1,j0:j1,k) * 100._r8 ! m/s -> cm/s
      edyn_wn   (k,i0:i1,j0:j1) = wn        (i0:i1,j0:j1,k) * 100._r8 ! m/s -> cm/s
      edyn_zgi  (k,i0:i1,j0:j1) = zgi       (i0:i1,j0:j1,k) * 100._r8 ! m -> cm
      edyn_zht  (k,i0:i1,j0:j1) = zht       (i0:i1,j0:j1,k) * 100._r8 ! m -> cm
      edyn_mbar (k,i0:i1,j0:j1) = mbar      (i0:i1,j0:j1,k)
      edyn_ped  (k,i0:i1,j0:j1) = sigma_ped (i0:i1,j0:j1,k)
      edyn_hall (k,i0:i1,j0:j1) = sigma_hall(i0:i1,j0:j1,k)
      edyn_ui   (k,i0:i1,j0:j1) = ui        (i0:i1,j0:j1,k) * 100._r8 ! zonal ion drift (m/s -> cm/s)
      edyn_vi   (k,i0:i1,j0:j1) = vi        (i0:i1,j0:j1,k) * 100._r8 ! meridional ion drift (m/s -> cm/s)
      edyn_wi   (k,i0:i1,j0:j1) = wi        (i0:i1,j0:j1,k) * 100._r8 ! vertical ion drift (m/s -> cm/s)
!
! Additional fields for oplus:
!
      edyn_te   (k,i0:i1,j0:j1) = te     (i0:i1,j0:j1,k)
      edyn_ti   (k,i0:i1,j0:j1) = ti     (i0:i1,j0:j1,k)
      edyn_o2   (k,i0:i1,j0:j1) = o2mmr  (i0:i1,j0:j1,k)
      edyn_o1   (k,i0:i1,j0:j1) = o1mmr  (i0:i1,j0:j1,k)
      edyn_n2   (k,i0:i1,j0:j1) = n2mmr  (i0:i1,j0:j1,k)
      edyn_om   (k,i0:i1,j0:j1) = -(omega(i0:i1,j0:j1,k) / pmid(i0:i1,k,j0:j1)) ! Pa/s -> 1/s
      edyn_op   (k,i0:i1,j0:j1) = op     (i0:i1,j0:j1,k) / 1.e6_r8  ! m^3 -> cm^3
      edyn_optm1(k,i0:i1,j0:j1) = optm1  (i0:i1,j0:j1,k) / 1.e6_r8  ! m^3 -> cm^3
    enddo
!
! At first timestep, allocate optm1 module data, and initialize local
! edyn_optm1 to op from physics. This will be input to oplus_xport.
! After oplus_xport, optm1 will be updated from local oplus_xport output.
! After first timestep, simply update edyn_optm1 from optm1.
! optm1 is m^3 for waccm, whereas edyn_optm1 is cm^3 for oplus_xport.
!
! At this point, everything is in waccm format. The locals edyn_op and
! edyn_optm1 will be converted to tiegcm format for the call to oplus_xport,
! then oplus_xport output (opnm_out) will be converted back to waccm format
! before using it to update optm1 module data.
!
!    if (nstep==1) then
!      optm1 = 0._r8
!      do k=1,nlev
!        edyn_optm1(k,i0:i1,j0:j1) = op(i0:i1,j0:j1,k) / 1.e6_r8  ! m^3 -> cm^3
!      enddo
!
! After the first step, edyn_optm1 input is updated from the module data
! (note edyn_optm1 will be converted to TIEGCM format before being
! passed in to oplus_xport)
!
!    else ! nstep > 1
!      do k=1,nlev
!        edyn_optm1(k,i0:i1,j0:j1) = optm1(i0:i1,j0:j1,k) / 1.e6_r8 ! m^3 -> cm^3
!      enddo
!    endif
!
! These are in WACCM format, and most are in CGS units (see above):
! (units are specified in addfld calls, edyn_init.F90)
!
    call savefld_waccm(edyn_tn   ,'DPIE_TN'  ,nlev,i0,i1,j0,j1)  ! deg K
    call savefld_waccm(edyn_un   ,'DPIE_UN'  ,nlev,i0,i1,j0,j1)  ! cm/s
    call savefld_waccm(edyn_vn   ,'DPIE_VN'  ,nlev,i0,i1,j0,j1)  ! cm/s
    call savefld_waccm(edyn_wn   ,'DPIE_WN'  ,nlev,i0,i1,j0,j1)  ! cm/s
    call savefld_waccm(edyn_om   ,'DPIE_OM'  ,nlev,i0,i1,j0,j1)  ! omega on midpoints (1/s)
    call savefld_waccm(edyn_zht  ,'DPIE_ZHT' ,nlev,i0,i1,j0,j1)  ! geometric height (cm)
    call savefld_waccm(edyn_zgi  ,'DPIE_ZGI' ,nlev,i0,i1,j0,j1) ! geopotential height on interfaces (cm)
    call savefld_waccm(edyn_mbar ,'DPIE_BARM',nlev,i0,i1,j0,j1)  ! mean mass
    call savefld_waccm(edyn_o2   ,'DPIE_O2'  ,nlev,i0,i1,j0,j1)  ! cm^3
    call savefld_waccm(edyn_o1   ,'DPIE_O'   ,nlev,i0,i1,j0,j1)  ! cm^3
    call savefld_waccm(edyn_n2   ,'DPIE_N2'  ,nlev,i0,i1,j0,j1)  ! cm^3
    call savefld_waccm(edyn_te   ,'DPIE_TE'  ,nlev,i0,i1,j0,j1)
    call savefld_waccm(edyn_ti   ,'DPIE_TI'  ,nlev,i0,i1,j0,j1)
!
! Save electron density to TIEGCM-format file (edynamo.nc):
! (ne(i,j,k) was calculated in m^3 above, save here in cm^3)
!
    do j=j0,j1
      do i=i0,i1
        do k=1,nlev 
          edyn_ne(k,i,j) = ne(i,j,k)*1.e-6_r8 ! m^3 -> cm^3
        enddo
      enddo
    enddo
!
! Convert input fields from "WACCM format" to "TIEGCM format" 
! (phase shift longitude data and invert the vertical dimension).
!
    if (ionos_edyn_active) then
       nfields = 21
       allocate(ptrs(nfields))
       !
       ! Fields needed for edynamo:
       ptrs(1)%ptr => edyn_tn    ; ptrs(2)%ptr => edyn_un   ; ptrs(3)%ptr => edyn_vn 
       ptrs(4)%ptr => edyn_wn    ; ptrs(5)%ptr => edyn_zht  ; ptrs(6)%ptr => edyn_zgi
       ptrs(7)%ptr => edyn_mbar  ; ptrs(8)%ptr => edyn_ped  ; ptrs(9)%ptr => edyn_hall
       !
       ! Additional fields needed for oplus (and Ne for diag):
       ptrs(10)%ptr => edyn_te  ; ptrs(11)%ptr => edyn_ti    ; ptrs(12)%ptr => edyn_o2   
       ptrs(13)%ptr => edyn_o1  ; ptrs(14)%ptr => edyn_n2    ; ptrs(15)%ptr => edyn_om   
       ptrs(16)%ptr => edyn_op  ; ptrs(17)%ptr => edyn_optm1 ; ptrs(18)%ptr => edyn_ne
       ptrs(19)%ptr => edyn_ui  ; ptrs(20)%ptr => edyn_vi    ; ptrs(21)%ptr => edyn_wi
       !
       ! Convert from WACCM to TIEGCM format:
       call switch_model_format(ptrs,1,nlev,i0,i1,j0,j1,nfields)
       deallocate(ptrs)
    endif
!
! Call electrodynamo (edynamo.F90)
! If using time3d conductances, tell dynamo to *not* do fieldline 
! integrations (i.e., do_integrals == false). In this case, edynamo
! conductances zigmxx,rim1,2 from time3d will be set by subroutine
! transform_glbin in time3d module.
!
    do_integrals = .true.
!
! If ionos_edyn_active=false, then empirical ion drifts were passed in from physics,
! otherwise dynamo calculates them here, and they will be passed to physics.
!
    if (ionos_edyn_active) then

       if (debug.and.masterproc) then
         write(iulog,"('dpie_coupling call dynamo... nstep=',i8)") nstep
         write(iulog,*) 'dpie_coupling: before dynamo MIN/MAX edyn_tn ', &
               MINVAL(edyn_tn(:,i0:i1,j0:j1)), MAXVAL(edyn_tn(:,i0:i1,j0:j1))   
         write(iulog,*) 'dpie_coupling: before dynamo MIN/MAX edyn_un ', &
              MINVAL(edyn_un(:,i0:i1,j0:j1)), MAXVAL(edyn_un(:,i0:i1,j0:j1))   
         write(iulog,*) 'dpie_coupling: before dynamo MIN/MAX edyn_vn ', &
              MINVAL(edyn_un(:,i0:i1,j0:j1)), MAXVAL(edyn_vn(:,i0:i1,j0:j1))   
         write(iulog,*) 'dpie_coupling: before dynamo MIN/MAX edyn_wn ', &
              MINVAL(edyn_wn(:,i0:i1,j0:j1)), MAXVAL(edyn_wn(:,i0:i1,j0:j1))   
         write(iulog,*) 'dpie_coupling: before dynamo MIN/MAX edynz_zgi ', &
              MINVAL(edyn_zgi(:,i0:i1,j0:j1)), MAXVAL(edyn_zgi(:,i0:i1,j0:j1))   
         write(iulog,*) 'dpie_coupling: before dynamo MIN/MAX edynz_ped ', &
              MINVAL(edyn_ped(:,i0:i1,j0:j1)), MAXVAL(edyn_ped(:,i0:i1,j0:j1))   
         write(iulog,*) 'dpie_coupling: before dynamo MIN/MAX edynz_hall ', &
              MINVAL(edyn_hall(:,i0:i1,j0:j1)), MAXVAL(edyn_hall(:,i0:i1,j0:j1))   
         write(iulog,*) 'dpie_coupling: before dynamo MIN/MAX edynz_op ', &
              MINVAL(edyn_op(:,i0:i1,j0:j1)), MAXVAL(edyn_op(:,i0:i1,j0:j1))   
         write(iulog,*) 'dpie_coupling: before dynamo MIN/MAX edynz_optm1 ', &
              MINVAL(edyn_optm1(:,i0:i1,j0:j1)), MAXVAL(edyn_optm1(:,i0:i1,j0:j1)) 
       endif

       call t_startf('dpie_ionos_dynamo')
       call dynamo(edyn_tn, edyn_un, edyn_vn, edyn_wn, edyn_zgi,   &
                   edyn_ped, edyn_hall, edyn_ui, edyn_vi, edyn_wi, &
                   1,nlev,i0,i1,j0,j1,do_integrals)
       call t_stopf ('dpie_ionos_dynamo')

       if (debug.and.masterproc) then
          write(iulog,"('dpie_coupling after dynamo: nstep=',i8)") nstep
          write(iulog,"('  ui min,max (cm/s)=',2es12.4)") minval(edyn_ui),maxval(edyn_ui)
          write(iulog,"('  vi min,max (cm/s)=',2es12.4)") minval(edyn_vi),maxval(edyn_vi)
          write(iulog,"('  wi min,max (cm/s)=',2es12.4)") minval(edyn_wi),maxval(edyn_wi)
       endif
    else
       if (debug.and.masterproc) then
          write(iulog,"('dpie_coupling (dynamo NOT called): nstep=',i8)") nstep
          write(iulog,"('  empirical ExB ui min,max (cm/s)=',2es12.4)") minval(ui),maxval(ui)
          write(iulog,"('  empirical ExB vi min,max (cm/s)=',2es12.4)") minval(vi),maxval(vi)
          write(iulog,"('  empirical ExB wi min,max (cm/s)=',2es12.4)") minval(wi),maxval(wi)
       endif
    endif
!
! Call O+ transport routine.  Now all inputs to oplus_xport should be in 
! tiegcm-format wrt longitude (-180->180), vertical (bot2top), and units (CGS).
! (Composition is mmr, ne is cm^3, winds are cm/s)
! Output op_out and opnm_out will be in cm^3, converted to mmr below.
!
    if (ionos_oplus_xport) then
      pmid_inv(1:nlev) = pref_mid(nlev:1:-1) ! invert ref pressure (Pa) as in tiegcm


!
! Transport O+ (all args in 'TIEGCM format')
! Subcycle oplus_xport nspltop times.
!
      if (debug.and.masterproc) &
         write(iulog,"('dpie_coupling before subcycling oplus_xport: nstep=',i8,' nspltop=',i3)") nstep,nspltop

      call t_startf('dpie_halo')
!$omp parallel do private(i, j, k)
      do k=1,nlev
         do j=j0,j1
            do i=i0,i1
               halo_tn(k,i,j)   = edyn_tn(k,i,j)
               halo_te(k,i,j)   = edyn_te(k,i,j)
               halo_ti(k,i,j)   = edyn_ti(k,i,j)
               halo_un(k,i,j)   = edyn_un(k,i,j)
               halo_vn(k,i,j)   = edyn_vn(k,i,j)
               halo_om(k,i,j)   = edyn_om(k,i,j)
               halo_o2(k,i,j)   = edyn_o2(k,i,j)
               halo_o1(k,i,j)   = edyn_o1(k,i,j)
               halo_n2(k,i,j)   = edyn_n2(k,i,j)
               halo_mbar(k,i,j) = edyn_mbar(k,i,j)
            enddo
         enddo
      enddo
      !
      ! Define halo points on inputs:
      ! WACCM has global longitude values at the poles (j=1,j=nlev)
      ! (they are constant for most, except the winds.)
      !
      ! Set two halo points in lat,lon:
      !
      nfields=10
      allocate(ptrs(nfields),polesign(nfields))
      ptrs(1)%ptr => halo_tn ; ptrs(2)%ptr => halo_te ; ptrs(3)%ptr => halo_ti
      ptrs(4)%ptr => halo_un ; ptrs(5)%ptr => halo_vn ; ptrs(6)%ptr => halo_om
      ptrs(7)%ptr => halo_o2 ; ptrs(8)%ptr => halo_o1 ; ptrs(9)%ptr => halo_n2
      ptrs(10)%ptr => halo_mbar

      polesign = 1._r8
      polesign(4:5) = -1._r8 ! un,vn
     
      call mp_geo_halos(ptrs,1,nlev,i0,i1,j0,j1,nfields)
      !
      ! Set latitude halo points over the poles (this does not change the poles).
      ! (the 2nd halo over the poles will not actually be used (assuming lat loops
      !  are lat=2,plat-1), because jp1,jm1 will be the pole itself, and jp2,jm2 
      !  will be the first halo over the pole)
      !
      call mp_pole_halos(ptrs,1,nlev,i0,i1,j0,j1,nfields,polesign)
      deallocate(ptrs,polesign)
      call t_stopf('dpie_halo')

      call t_startf('dpie_oplus_xport')
      do isplit=1,nspltop

         if (isplit > 1) then
            edyn_op = op_out
            edyn_optm1 = opnm_out
         endif

         call oplus_xport(halo_tn,halo_te,halo_ti,halo_un,halo_vn,halo_om,            &
              edyn_zgi,halo_o2,halo_o1,halo_n2,edyn_op,edyn_optm1,        &
              halo_mbar,edyn_ui,edyn_vi,edyn_wi,pmid_inv,    &
              op_out,opnm_out, &
              i0,i1,j0,j1,nspltop,isplit)

      enddo ! isplit=1,nspltop
      call t_stopf ('dpie_oplus_xport')

      if (debug.and.masterproc) then
         write(iulog,"('dpie_coupling after subcycling oplus_xport: nstep=',i8,' nspltop=',i3)") &
              nstep,nspltop
         write(iulog,"('  op_out   min,max (cm^3)=',2es12.4)") minval(op_out)  ,maxval(op_out)
         write(iulog,"('  opnm_out min,max (cm^3)=',2es12.4)") minval(opnm_out),maxval(opnm_out)
      endif

   endif ! ionos_oplus_xport
!
! Convert ion drifts and O+ output from TIEGCM to WACCM format:
!
   if (ionos_edyn_active) then
      nfields = 5 ! ui,vi,wi,op,opnm
      allocate(ptrs(nfields))
      ptrs(1)%ptr => edyn_ui ; ptrs(2)%ptr => edyn_vi ; ptrs(3)%ptr => edyn_wi
      ptrs(4)%ptr => op_out  ; ptrs(5)%ptr => opnm_out
      call switch_model_format(ptrs,1,nlev,i0,i1,j0,j1,nfields)
      deallocate(ptrs)
   endif
!
    if (ionos_oplus_xport) then
      call savefld_waccm(op_out,'OPLUS',nlev,i0,i1,j0,j1) ! cm^3
!
! Pass new O+ for current and previous time step back to physics (convert from cm^3 to m^3 and back to mmr).
!
      do k=1,nlev
        do j=j0,j1
          do i=i0,i1
            opmmr(i,j,k) = op_out(k,i,j)*1.e6_r8 * rmassOp / mbar(i,j,k) * &
              (kboltz * tn(i,j,k)) / pmid(i,k,j)
            op_out(k,i,j) = opmmr(i,j,k) ! for save to waccm hist in mmr
            opmmrtm1(i,j,k) = opnm_out(k,i,j)*1.e6_r8 * rmassOp / mbar(i,j,k) * &
              (kboltz * tn(i,j,k)) / pmid(i,k,j)
          enddo
        enddo
      enddo

    endif ! ionos_oplus_xport
!
! Convert ion drifts from cm/s to m/s for WACCM physics and history files.
!   real(r8),intent(inout) :: ui(i0:i1,j0:j1,nlev)      ! zonal ion drift (edynamo or empirical)
!
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          ui(i,j,k) = edyn_ui(k,i,j)/100._r8
          vi(i,j,k) = edyn_vi(k,i,j)/100._r8
          wi(i,j,k) = edyn_wi(k,i,j)/100._r8
        enddo
      enddo
    enddo
    call savefld_waccm(edyn_ui/100._r8,'WACCM_UI',nlev,i0,i1,j0,j1)
    call savefld_waccm(edyn_vi/100._r8,'WACCM_VI',nlev,i0,i1,j0,j1)
    call savefld_waccm(edyn_wi/100._r8,'WACCM_WI',nlev,i0,i1,j0,j1)

    call t_stopf('d_pie_coupling')

  end subroutine d_pie_coupling
!-----------------------------------------------------------------------
  subroutine calc_wn(tn,omega,pmid,mbar,grav,wn,i0,i1,j0,j1,nlev)
    use shr_const_mod,only : shr_const_rgas ! Universal gas constant
!
! Calculate neutral vertical wind on midpoints (m/s)
!
! Inputs:
    integer,intent(in) :: i0,i1,j0,j1,nlev
    real(r8),dimension(i0:i1,j0:j1,nlev),intent(in) :: &
      tn,   &  ! neutral temperature (deg K)
      omega,&  ! pressure velocity (Pa/s)
      mbar     ! mean molecular weight
    real(r8),dimension(i0:i1,nlev,j0:j1),intent(in) :: &
      pmid     ! pressure at midpoints (Pa)
    real(r8),intent(in) :: grav ! m/s^2
!
! Output:
    real(r8),intent(out) :: wn(i0:i1,j0:j1,nlev)    ! vertical velocity output (m/s)
!
! Local:
    integer :: i,j,k
    real(r8) :: scheight(i0:i1,j0:j1,nlev) ! dimensioned for vectorization

    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          scheight(i,j,k) = shr_const_rgas*tn(i,j,k)/(mbar(i,j,k)*grav)
          wn(i,j,k) = -omega(i,j,k)*scheight(i,j,k)/pmid(i,k,j)
        enddo
      enddo
    enddo
  end subroutine calc_wn
!-----------------------------------------------------------------------
  subroutine calc_pfrac(sunlon,pfrac)
!
! Calculate pfrac fractional presence of dynamo equation using critical
!  convection colatitudes crit(2).
!
    use edyn_maggrid  ,only: nmlonp1,ylonm,ylatm
    use edyn_solve    ,only: nmlat0
    use aurora_params ,only: offc, dskofc, theta0, aurora_params_set

    implicit none
!
! Args:
    real(r8),intent(in) :: sunlon  ! Sun's longitude in dipole coordinates
!
! Output: fractional presence of dynamo equation using critical colatitudes
!
    real(r8),intent(out) :: pfrac(nmlonp1,nmlat0) ! NH fraction of potential
!
! Local:
    integer :: j,i
    real(r8),dimension(nmlonp1,nmlat0) :: colatc
    real(r8) :: sinlat,coslat,aslonc,ofdc,cosofc,sinofc,crit1deg

    if (.not. crit_user_set) then
       if (amie_period) then
          crit(:) = amie_default_crit(:)*dtr
       else
          crit1deg = max(15._r8,0.5_r8*(theta0(1)+theta0(2))*rtd + 5._r8)
          crit1deg = min(30._r8,crit1deg)
          !
          ! Critical colatitudes:
          crit(1) = crit1deg*dtr
          crit(2) = crit(1) + 15._r8*dtr
       endif
    endif

    if (.not.aurora_params_set) then
       offc(:)   =  1._r8*dtr
       dskofc(:) =  0._r8
    endif

!    
! offc(2), dskofc(2) are for northern hemisphere aurora 
!
    ofdc = sqrt(offc(2)**2+dskofc(2)**2)
    cosofc = cos(ofdc)
    sinofc = sin(ofdc)
    aslonc = asin(dskofc(2)/ofdc)
!
! Define colatc with northern convection circle coordinates
!
    do j=1,nmlat0
      sinlat = sin(abs(ylatm(j+nmlat0-1)))
      coslat = cos(    ylatm(j+nmlat0-1))
      do i=1,nmlonp1
        colatc(i,j) = cos(ylonm(i)-sunlon+aslonc)
        colatc(i,j) = acos(cosofc*sinlat-sinofc*coslat*colatc(i,j))
      enddo ! i=1,nmlonp1
!
! Calculate fractional presence of dynamo equation at each northern
! hemisphere geomagnetic grid point. Output in pfrac(nmlonp1,nmlat0)
!
      do i=1,nmlonp1
        pfrac(i,j) = (colatc(i,j)-crit(1))/(crit(2)-crit(1))
        if (pfrac(i,j) < 0._r8) pfrac(i,j)  = 0._r8
        if (pfrac(i,j) >= 1._r8) pfrac(i,j) = 1._r8
      enddo ! i=1,nmlonp1
    enddo ! j=1,nmlat0
!   
  end subroutine calc_pfrac
!-----------------------------------------------------------------------
  subroutine sunloc(iday,secs,sunlons)
!
! Given day of year and ut, return sun's longitudes in dipole coordinates 
! in sunlons(nlat)
!
    use getapex       ,only: alonm     ! (nlonp1,0:nlatp1)
    use edyn_geogrid  ,only: nlon,nlat
    use edyn_params   ,only: pi
!
! Args:
    integer,intent(in)   :: iday          ! day of year
    real(r8),intent(in)  :: secs          ! ut in seconds
    real(r8),intent(out) :: sunlons(nlat) ! output
!
! Local:
    integer :: j,i,ii,isun,jsun
    real(r8) :: glats,glons,pisun,pjsun,sndlons,csdlons
    real(r8) :: dphi,dlamda
    real(r8) :: rlonm(nlon+4,nlat) ! (nlon+4,nlat)
    real(r8) :: r8_nlat, r8_nlon
    real(r8) :: r8_isun, r8_jsun

!
! Sun's geographic coordinates:
    r8_nlat = dble(nlat)
    r8_nlon = dble(nlon)
    glats   = asin(.398749_r8*sin(2._r8*pi*(iday-80)/365._r8))
    glons   = pi*(1._r8-2._r8*secs/86400._r8)
    dphi    = pi/r8_nlat
    dlamda  = 2._r8*pi/r8_nlon

    do j=1,nlat
       do i=1,nlon
          ii = i+2
          rlonm(ii,j) = alonm(i,j)
       enddo
       do i=1,2
          rlonm(i,j) = rlonm(i+nlon,j)
          rlonm(i+nlon+2,j) = rlonm(i+2,j)
       enddo
    enddo

    pisun = (glons+pi)/dlamda+1._r8
    pjsun = (glats+.5_r8*(pi-dphi))/dphi+1._r8
    isun = int(pisun)
    jsun = int(pjsun)
    r8_isun = dble(isun)
    r8_jsun = dble(jsun)
    pisun = pisun-r8_isun
    pjsun = pjsun-r8_jsun

    sndlons = (1._r8-pisun)*(1._r8-pjsun)*sin(rlonm(isun+2,jsun))+   &
               pisun*(1._r8-pjsun)       *sin(rlonm(isun+3,jsun))+   &
               pisun*pjsun               *sin(rlonm(isun+3,jsun+1))+ &
               (1._r8-pisun)*pjsun       *sin(rlonm(isun+2,jsun+1))
    csdlons = (1._r8-pisun)*(1._r8-pjsun)*cos(rlonm(isun+2,jsun))+   &
               pisun*(1._r8-pjsun)       *cos(rlonm(isun+3,jsun))+   &
               pisun*pjsun               *cos(rlonm(isun+3,jsun+1))+ &
               (1._r8-pisun)*pjsun       *cos(rlonm(isun+2,jsun+1))
    sunlons(1) = atan2(sndlons,csdlons)
    do j = 2,nlat
      sunlons(j) = sunlons(1)
    enddo

  end subroutine sunloc
!-----------------------------------------------------------------------
end module dpie_coupling
