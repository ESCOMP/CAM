module edynamo
!
!-----------------------------------------------------------------------
! Purpose:
! Electro-dynamo module
!-----------------------------------------------------------------------
!
  use shr_kind_mod ,only: r8 => shr_kind_r8 ! 8-byte reals
  use cam_logfile  ,only: iulog
  use cam_abortutils,only: endrun
  use spmd_utils     ,only: masterproc
#ifdef WACCMX_EDYN_ESMF
  use edyn_params  ,only: finit      ! initialization value 
  use edyn_maggrid ,only: nmlon,nmlonp1,nmlat,nmlath,nmlev
  use edyn_mpi     ,only: mlon0,mlon1,omlon1,mlat0,mlat1,mlev0,mlev1,mytid,&
                          lon0,lon1,lat0,lat1,lev0,lev1
  use edyn_solve   ,only: solve_edyn
  use time_manager, only: get_nstep ! for debug
  use cam_history,  only : outfld, hist_fld_active
  use savefield_waccm,only: savefld_waccm_switch
  use esmf, only : ESMF_KIND_R8, ESMF_Field        ! ESMF library module
#endif

  implicit none
  save
  private

#ifdef WACCMX_EDYN_ESMF
  integer :: nstep
!
! 3d pointers to fields regridded to magnetic subdomains (i,j,k):
! (mlon0:mlon1,mlat0:mlat1,nmlev)
!
  real(ESMF_KIND_R8),pointer,dimension(:,:,:) :: & ! 3d fields on mag grid
    ped_mag,     & ! pedersen conductivity on magnetic grid
    hal_mag,     & ! hall conductivity on magnetic grid
    zpot_mag,    & ! geopotential on magnetic grid
    scht_mag,    & ! scale height on magnetic grid
    adotv1_mag,  & ! ue1 (m/s)
    adotv2_mag     ! ue2 (m/s)
!
! 2d pointers to fields on magnetic subdomains (i,j):
! (mlon0:mlon1,mlat0:mlat1)
!
  real(ESMF_KIND_R8),pointer,dimension(:,:) :: &
    sini_mag,    & ! sin(I_m)
    adota1_mag,  & ! d(1)**2/D
    adota2_mag,  & ! d(2)**2/D
    a1dta2_mag,  & ! (d(1) dot d(2)) /D
    be3_mag        ! mag field strength (T)
!
! 2d coefficients and RHS terms for PDE on magnetic subdomains
! (including halo points).
! If use_time3d_integ==.true., these will be input from time3d
! (see use-association in time3d.F90)
!
  real(r8),allocatable,dimension(:,:) :: &
    zigm11,    & ! sigma11*cos(theta0)
    zigmc,     & ! sigmac
    zigm2,     & ! sigma2
    zigm22,    & ! sigma22/cos(theta0)
    rim1,rim2, & ! see description in comment below
    rhs,       & ! right-hand side of PDE
    phimsolv,  & ! solution direct from solver (nhem only)
    phim2d       ! solution with phihm and both nhem and shem
!
! 3d potential and electric field on mag subdomains (see sub pthreed):
! (mlon0:mlon1,mlat0:mlat1,mlev0:mlev1)
! Electric potential and field components are output fields of edynamo 
! (later, these can be output arguments of the main driver, sub dynamo)
!
  real(r8),allocatable,dimension(:,:,:) :: &
    phim3d,       & ! 3d electric potential
    ed13d,ed23d,  & ! 3d electric field for current calculations
    ephi3d,       & ! 3d eastward electric field
    elam3d,       & ! 3d equatorward electric field
    emz3d,        & ! 3d upward electric field
    zpotm3d         ! 3d geopotential (values at all levels)
!
! 3d ion drift velocities on geographic grid (output):
!
! real(r8),allocatable,dimension(:,:,:),save,target :: & ! (nlev,lon0:lon1,lat0:lat1)
!   ui,   & ! zonal ion drift
!   vi,   & ! meridional ion drift
!   wi      ! vertical ion drift
!
! 3d electric field on geographic subdomains (see sub pefield):
! (nlev,lon0-2,lon1+2,lat0:lat1)
  real(r8),allocatable,dimension(:,:,:) :: ex,ey,ez
!
! 3d electric potential on geographic subdomains (lon0:lon1,lat0:lat1,nlevp1)
! This will be regridded from phim3d for output to history files.
      real(r8),allocatable,dimension(:,:,:) :: phig3d ! (lon0:lon1,lat0:lat1,nlevp1)
      real(r8),allocatable,dimension(:,:,:) :: poten  ! (nlevp1,lon0:lon1,lat0:lat1)
!
! Fields at mag equator:
!
  real(r8),allocatable,dimension(:,:) ::   & ! (mlon0:mlon1,nmlev)
    ped_meq, hal_meq, adotv1_meq, adotv2_meq, zpot_meq
  real(r8),allocatable,dimension(:,:,:) :: & ! (mlon0:mlon1,nmlev,4)
    fmeq_out  
  real(r8),allocatable,dimension(:,:,:,:) :: & ! (mlon0:mlon1,mlat0:mlat1,nmlev,4)
    fmeq_in 
! 
! Global longitude values near mag equator and poles for complete_integrals and rhs.
! These are declared in module data because they are used by subs complete_integrals 
! and rhspde.  The nf2d 6 fields are: zigm11,zigm22,zigmc,zigm2,rim1,rim2, 
! order is important (see feq_jpm1 and fpole_jpm2)!
!
  integer,parameter :: nf2d=6            ! 6 2d fields
  real(r8) :: feq_jpm1(nmlonp1,2,nf2d)   ! 6 fields at 2 lats (eq-1, eq+1)
  real(r8) :: fpole_jpm2(nmlonp1,4,nf2d) ! fields at S pole+1,2 and N pole-1,2

  real(r8),parameter :: unitvm(nmlon)=1._r8
!
! ed1,ed2: 2d electric field output on mag grid:
! (use-associated by dpie_coupling)
!
  real(r8),allocatable,dimension(:,:) :: ed1,ed2 ! (mlon0-1:mlon1+1,mlat0-1:mlat1+1)
!
! Global inputs to time3d: Note dimension order switch: 
!   edynamo has subdomains (mlon,mlat), whereas time3d has global (nmlat,nmlonp1)
! These are use-associated by time3d, and are init to zero in edyn_init.
!
  real(r8),dimension(nmlat,nmlonp1) :: ed1_glb,ed2_glb
  logical :: do_integ           ! from input arg do_integrals
  logical :: debug=.false. ! set true for prints to stdout at each call

  public alloc_edyn,ed1,ed2,ed1_glb,ed2_glb
  public zigm11,zigmc,zigm2,zigm22,rim1,rim2
#endif
  public :: dynamo

  contains
!-----------------------------------------------------------------------
  subroutine dynamo(tn,un,vn,wn,zpot,ped,hall,ui,vi,wi, &
    lev0,lev1,lon0,lon1,lat0,lat1,do_integrals)
  use edyn_mpi,only:      &
    mp_mag_halos,         & ! set magnetic halo points
    mp_scatter_phim         ! scatter solution to slave tasks
  use edyn_solve,only: rim_glb   ! pde solver output (nmlonp1,nmlat,2)
!
! Main driver for edynamo.
! Note alloc_edyn and esmf_init are called from edyn_init.
!
! Args:
    integer,intent(in) :: & ! geographic subdomain
      lev0,lev1,  & ! first,last level indices (not distributed)
      lon0,lon1,  & ! first,last longitude indices of geographic subdomain
      lat0,lat1     ! first,last latitude indices of geographic subdomain
    logical,intent(in) :: do_integrals
!
! Inputs from neutral atmosphere (on geographic subdomain):
! (intent(inout) because they are passed to sub dynamo_input)
!
    real(r8),intent(inout),dimension(lev0:lev1,lon0:lon1,lat0:lat1) :: &
      tn,         & ! neutral temperature (deg K)
      un,         & ! neutral zonal wind velocity (cm/s)
      vn,         & ! neutral meridional wind velocity (cm/s)
      wn,         & ! neutral vertical wind velocity (cm/s)
      zpot,       & ! geopotential height (cm)
      ped,        & ! pedersen conductivity (S/m)
      hall,       & ! hall conductivity (S/m)
      ui,         & ! zonal ion drift (cm/s)
      vi,         & ! meridional ion drift (cm/s)
      wi            ! vertical ion drift (cm/s)
#ifdef WACCMX_EDYN_ESMF
    if (debug) then
       nstep = get_nstep()
       write(iulog,"('Enter dynamo: nstep=',i5,' do_integrals=',l1)") nstep,do_integrals
    endif

    do_integ = do_integrals ! do_integ is module data
!
! Regrid input fields from geographic to magnetic, and calculate 
! some additional fields. If conductances are passed in from
! time3d (.not.do_integrals), then we do not need these inputs.
!
    if (do_integrals) then
      call dynamo_input(tn,un,vn,wn,zpot,ped,hall,&
        lev0,lev1,lon0,lon1,lat0,lat1)
      if (debug) write(iulog,"('edynamo debug: after dynamo_input')")
    endif
!
! Fieldline integration:
!
! If *not* doing fieldline integrations, then global conductances
! were passed in to the driver from time3d, and transformed from
! (nmlat,nmlonp1) to (nmlonp1,nmlat), defining zigmxx and rim1,2
! for the solver.
!
    if (do_integrals) call fieldline_integrals
!
! Equatorial and polar values, hemisphere folding:
! (these will be time3d integrations if do_integrals==.false.)
!
    call complete_integrals
    if (debug) write(iulog,"('edynamo debug: after complete_integrals')")
!
! Calculate right-hand side on mag subdomains:
! (mag halos are needed in rim1,2 for rhs calculation)
!
    call mp_mag_halos(rim1,mlon0,mlon1,mlat0,mlat1,1) 
    call mp_mag_halos(rim2,mlon0,mlon1,mlat0,mlat1,1)
    call rhspde
    if (debug) write(iulog,"('edynamo debug: after rhspde')")
!
! Gather needed arrays to root task for the serial solver:
!
    call gather_edyn
    if (debug) write(iulog,"('edynamo debug: after gather_edyn')")
!
! Root task now sets up stencils and calls the PDE solver: 
!
    if (debug) write(iulog,"('edynamo debug: call solve_edyn (master only)')")
    if (mytid==0) then
      call solve_edyn
    endif
    if (debug) write(iulog,"('edynamo debug: after solve_edyn (master only)')")
!
! rim1 after solver is needed for highlat_poten. rim_glb is distributed
! to subdomains as rim1, and mag halos set. This will overwrite rim1 from
! fieldline_integrals, complete_integrals, etc.
!
    call mp_scatter_phim(rim_glb(:,:,1),rim1(mlon0:mlon1,mlat0:mlat1))
    if (debug) write(iulog,"('edynamo debug: after mp_scatter_phim')")

    call mp_mag_halos(rim1,mlon0,mlon1,mlat0,mlat1,1)
    if (debug) write(iulog,"('edynamo debug: after mp_mag_halos')")
!    
! Add high latitude potential from empirical model (heelis or weimer)
! to solution rim1, defining phim2d on mag subdomains.
!
    call highlat_poten 
    if (debug) write(iulog,"('edynamo debug: after highlat_poten')")
!
! Expand phim2d to phim3d, first setting mag halos in phim2d from
! hightlat_poten. phim3d will then be the final potential from pdynamo.
!
    call mp_mag_halos(phim2d,mlon0,mlon1,mlat0,mlat1,1)

    call pthreed
    if (debug) write(iulog,"('edynamo debug: after pthreed')")
!
! Convert electric field to geographic grid:
    call pefield
    if (debug) write(iulog,"('edynamo debug: after pefield')")

!
! Calculate ion drift velocities:
!

    call ionvel(zpot,ui,vi,wi)
    if (debug) write(iulog,"('edynamo debug: after ionvel')")
#else
    call endrun('ERROR: To use edymamo must build with cppdef WACCMX_EDYN_ESMF')
#endif
  end subroutine dynamo
!-----------------------------------------------------------------------
#ifdef WACCMX_EDYN_ESMF
  subroutine dynamo_input(tn,un,vn,wn,zpot,ped,hall,&
    lev0,lev1,lon0,lon1,lat0,lat1)
!     
! Input fields are in "TIEGCM format" and CGS units.
! Provide needed inputs to the dynamo by regridding the fields
! from geographic to magnetic.
!
  use edyn_params ,only: h0,kbotdyn
  use getapex     ,only: & ! (nlonp1,0:nlatp1)
      zb,                & ! downward component of magnetic field
      bmod                 ! magnitude of magnetic field (gauss?)
  use edyn_geogrid,only: nlev

  use edyn_mpi,only:   &
!   mp_periodic_f2d, & ! set 2d periodic points
!   mp_periodic_f3d, & ! set 3d periodic points
    mp_mageq           ! get global values at mag equator
    
  use edyn_esmf,only:   & ! use-associate grid definitions and subroutines
    geo_src_grid,     & ! geographic source grid (ESMF_Grid type)
    edyn_esmf_regrid,      & ! subroutine that calls ESMF to regrid a field
    edyn_esmf_set2d_geo,   & ! set values of a 2d ESMF field on geographic grid
    edyn_esmf_set3d_geo,   & ! set values of a 3d ESMF field on geographic grid
    edyn_esmf_get_3dfield, & ! retrieve values of a 3d ESMF field
    edyn_esmf_get_2dfield    ! retrieve values of a 2d ESMF field

  use edyn_esmf,only: & ! 3d ESMF fields on geographic grid
    geo_ped,        & ! pedersen conductivity
    geo_hal,        & ! hall conductivity
    geo_zpot,       & ! geopotential height
    geo_scht,       & ! scale height
    geo_adotv1,     & ! ue1 (m/s)
    geo_adotv2        ! ue2 (m/s)

  use edyn_esmf,only: & ! 2d ESMF fields on geographic grid
    geo_sini,       & ! sin(I_m)
    geo_adota1,     & ! d(1)**2/D
    geo_adota2,     & ! d(2)**2/D
    geo_a1dta2,     & ! (d(1) dot d(2)) /D
    geo_be3           ! mag field strength (T)

  use edyn_esmf,only: & ! 3d ESMF fields on geomagnetic grid
    mag_ped,        & ! pedersen conductivity
    mag_hal,        & ! hall conductivity
    mag_zpot,       & ! geopotential height
    mag_scht,       & ! scale height
    mag_adotv1,     & ! ue1 (m/s)
    mag_adotv2        ! ue2 (m/s)

  use edyn_esmf,only: & ! 3d fields on geographic grid (bundled?)
    nf_3dgeo,       & ! number of 3d geo fields
    f_3dgeo           ! array of nf_3dgeo pointers to 3d geo fields

  use edyn_esmf,only: & ! 2d ESMF fields on geomagnetic grid
    mag_sini,       & ! sin(I_m)
    mag_adota1,     & ! d(1)**2/D
    mag_adota2,     & ! d(2)**2/D
    mag_a1dta2,     & ! (d(1) dot d(2)) /D
    mag_be3           ! mag field strength (T)

  use edyn_esmf,only: edyn_esmf_update_step ! indicates ESMF updated the current time step with updated geo-mag coordinates
  use edyn_esmf,only: edyn_esmf_update_flag
!
! Args: Input fields on geographic grid:
!
    integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
    real(r8),dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(in) :: &
      tn,    & ! neutral temperature (deg K) 
      un,    & ! neutral zonal velocity (cm/s)
      vn,    & ! neutral meridional velocity (cm/s)
      wn       ! neutral vertical velocity (cm/s)

    real(r8),dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(inout) :: &
      zpot,  & ! geopotential height (cm)
      ped,   & ! pedersen conductivity (S/m)
      hall     ! hall conductivity (S/m)
!
! Local:
!
    integer :: j,i,k
    real(r8),dimension(lev0:lev1,lon0:lon1,lat0:lat1) :: &
      scheight,  & ! scale height (no longer necessary since wn calculated outside)
      adotv1,    & ! ue1 (m/s)
      adotv2       ! ue2 (m/s)
    real(r8),dimension(lon0:lon1,lat0:lat1) :: &
      sini,      & ! sin(I_m)
      adota1,    & ! d(1)**2/D
      adota2,    & ! d(2)**2/D
      a1dta2,    & ! (d(1) dot d(2)) /D
      be3          ! mag field strength (T)

!
! See nf_3dgeo above:
    real(r8),dimension(lev0:lev1,lon0:lon1,lat0:lat1,nf_3dgeo) :: f3d
    character(len=8) :: fnames(nf_3dgeo)
!
! For wc timing:
!   real(r8) :: starttime,endtime

    scheight = 0._r8

!   starttime = mpi_wtime()

    if (debug) write(iulog,"('Enter dynamo_input')")
!
! Save 3d input fields on geo grid to WACCM history:
   call savefld_waccm_switch(tn   ,'EDYN_TN'  ,nlev,lon0,lon1,lat0,lat1)
   call savefld_waccm_switch(un   ,'EDYN_UN'  ,nlev,lon0,lon1,lat0,lat1)
   call savefld_waccm_switch(vn   ,'EDYN_VN'  ,nlev,lon0,lon1,lat0,lat1)
   call savefld_waccm_switch(wn   ,'EDYN_WN'  ,nlev,lon0,lon1,lat0,lat1)
   call savefld_waccm_switch(zpot ,'EDYN_Z'   ,nlev,lon0,lon1,lat0,lat1)
   call savefld_waccm_switch(ped  ,'EDYN_PED' ,nlev,lon0,lon1,lat0,lat1)
   call savefld_waccm_switch(hall ,'EDYN_HALL',nlev,lon0,lon1,lat0,lat1)

!
    if (debug) write(iulog,"('dynamo_input after savefld_waccm calls')")
    if (debug) write(iulog,"('dynamo_input: kbotdyn=',i4)") kbotdyn
!
! Calculate some 2d and 3d fields:
    call calc_adotv(zpot,un,vn,wn,adotv1,adotv2,adota1,adota2, &
      a1dta2,be3,lev0,lev1,lon0,lon1,lat0,lat1)
    if (debug) write(iulog,"('dynamo_input after calc_adotv')")

    call savefld_waccm_switch(adotv1  ,'EDYN_ADOTV1',nlev,lon0,lon1,lat0,lat1)
    call savefld_waccm_switch(adotv2  ,'EDYN_ADOTV2',nlev,lon0,lon1,lat0,lat1)
!
! Calculate sini sin(I_m) (zb and bmod are from apex)
!
    do j=lat0,lat1
      do i=lon0,lon1
        sini(i,j) = zb(i,j)/bmod(i,j) ! sin(I_m)
      enddo
    enddo
!
! Set 3d field values on geographic source grid, including
! separate calculations at the poles.  This is consolidated
! into a single call, so mp_geopole_3d can be called by 
! esmf_set3d_geo once for all fields.
!
    fnames = (/'PED     ','HAL     ','ZPOT    ','SCHT    ',&
               'ADOTV1  ','ADOTV2  '/)

    f3d(:,:,:,1) = ped
    f3d(:,:,:,2) = hall
    f3d(:,:,:,3) = zpot
    f3d(:,:,:,4) = scheight
    f3d(:,:,:,5) = adotv1
    f3d(:,:,:,6) = adotv2

    f_3dgeo(1) = geo_ped
    f_3dgeo(2) = geo_hal
    f_3dgeo(3) = geo_zpot
    f_3dgeo(4) = geo_scht
    f_3dgeo(5) = geo_adotv1
    f_3dgeo(6) = geo_adotv2

    call edyn_esmf_set3d_geo(f_3dgeo,fnames,f3d,nf_3dgeo, &
      lev0,lev1,lon0,lon1,lat0,lat1)

    geo_ped    = f_3dgeo(1)
    geo_hal    = f_3dgeo(2)
    geo_zpot   = f_3dgeo(3)
    geo_scht   = f_3dgeo(4)
    geo_adotv1 = f_3dgeo(5)
    geo_adotv2 = f_3dgeo(6)

    ped     = f3d(:,:,:,1)
    hall    = f3d(:,:,:,2)
    zpot    = f3d(:,:,:,3)
    scheight= f3d(:,:,:,4)
    adotv1  = f3d(:,:,:,5)
    adotv2  = f3d(:,:,:,6)

    if (debug) write(iulog,"('dynamo_input after edyn_esmf_set3d_geo')")
!
! 2d fields need only be calculated in first timestep:
    if (edyn_esmf_update_step) then
!
! Set 2d field values on geographic grid:
! (esmf fields on source grid exclude periodic points)
!
      call edyn_esmf_set2d_geo(geo_sini,  geo_src_grid,'SINI    ',&
        sini,  lon0,lon1,lat0,lat1)
      call edyn_esmf_set2d_geo(geo_adota1,geo_src_grid,'ADOTA1  ',&
        adota1,lon0,lon1,lat0,lat1)
      call edyn_esmf_set2d_geo(geo_adota2,geo_src_grid,'ADOTA2  ',&
        adota2,lon0,lon1,lat0,lat1)
      call edyn_esmf_set2d_geo(geo_a1dta2,geo_src_grid,'A1DTA2  ',&
        a1dta2,lon0,lon1,lat0,lat1)
      call edyn_esmf_set2d_geo(geo_be3,   geo_src_grid,'BE3     ',&
        be3,   lon0,lon1,lat0,lat1)
      if (debug) write(iulog,"('dynamo_input after edyn_esmf_set2d_geo')")
    endif

! 
! Regrid 3d geo fields to mag grid:
    call edyn_esmf_regrid(geo_ped  ,mag_ped,  'geo2mag',3)
    call edyn_esmf_regrid(geo_hal  ,mag_hal,  'geo2mag',3)
    call edyn_esmf_regrid(geo_zpot ,mag_zpot, 'geo2mag',3)
    call edyn_esmf_regrid(geo_scht ,mag_scht, 'geo2mag',3)
    call edyn_esmf_regrid(geo_adotv1 ,mag_adotv1, 'geo2mag',3)
    call edyn_esmf_regrid(geo_adotv2 ,mag_adotv2, 'geo2mag',3)
    if (debug) write(iulog,"('dynamo_input after edyn_esmf_regrid')")
!
! Regrid time-independent 2d geo fields to mag grid:
    if (edyn_esmf_update_step) then
      call edyn_esmf_regrid(geo_sini   ,mag_sini  , 'geo2mag',2)
      call edyn_esmf_regrid(geo_adota1 ,mag_adota1, 'geo2mag',2)
      call edyn_esmf_regrid(geo_adota2 ,mag_adota2, 'geo2mag',2)
      call edyn_esmf_regrid(geo_a1dta2 ,mag_a1dta2, 'geo2mag',2)
      call edyn_esmf_regrid(geo_be3    ,mag_be3   , 'geo2mag',2)
    endif
!    
! Define edynamo module data pointers to the regridded mag fields.
! First arg of esmf_get_field is input esmf field (my_esmf module),
!   second arg is output data pointer (edynamo module)
! (These destination grid fields have periodic points allocated and set)
!     
! Get regridded 3d mag fields:
!
    call edyn_esmf_get_3dfield(mag_ped   ,ped_mag,   "PED     ")
    call edyn_esmf_get_3dfield(mag_hal   ,hal_mag,   "HAL     ")
    call edyn_esmf_get_3dfield(mag_zpot  ,zpot_mag,  "ZPOT    ")
    call edyn_esmf_get_3dfield(mag_scht  ,scht_mag,  "SCHT    ")
    call edyn_esmf_get_3dfield(mag_adotv1,adotv1_mag,"ADOTV1  ")
    call edyn_esmf_get_3dfield(mag_adotv2,adotv2_mag,"ADOTV2  ")
!
! Get regridded 2d mag fields (time-independent):
! First arg is input ESMF field, second is output pointer:
!
    if (edyn_esmf_update_step) then
      call edyn_esmf_get_2dfield(mag_sini  ,sini_mag  , "SINI    ")
      call edyn_esmf_get_2dfield(mag_adota1,adota1_mag, "ADOTA1  ")
      call edyn_esmf_get_2dfield(mag_adota2,adota2_mag, "ADOTA2  ")
      call edyn_esmf_get_2dfield(mag_a1dta2,a1dta2_mag, "A1A2M   ")
      call edyn_esmf_get_2dfield(mag_be3   ,be3_mag   , "BE3     ")
      call edyn_esmf_update_flag(.false.)
    endif
!
! fmeq_in are input fields on 3d mag subdomains.
! allocate(fmeq_in(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1,4)
!
    fmeq_in(:,:,:,1) = ped_mag(:,:,:)
    fmeq_in(:,:,:,2) = hal_mag(:,:,:)
    fmeq_in(:,:,:,3) = adotv1_mag(:,:,:)
    fmeq_in(:,:,:,4) = adotv2_mag(:,:,:)
!
! Tasks w/ mag equator send eq data(i,k) to other tasks in their tidi:
!
    call mp_mageq(fmeq_in,fmeq_out,4,mlon0,mlon1,mlat0,mlat1,nmlev)
!
! Output arrays now have mag equator data on longitude subdomain
!   and full column (mlon0:mlon1,nmlev)
! These will be used in fieldline_integrals.
!
    ped_meq(:,:)    = fmeq_out(:,:,1)
    hal_meq(:,:)    = fmeq_out(:,:,2)
    adotv1_meq(:,:) = fmeq_out(:,:,3)
    adotv2_meq(:,:) = fmeq_out(:,:,4)
!
! Save geopotential on magnetic grid in zpotm3d, then
! limit max zpot_mag to h0 for use in fieldline integrals
! and pthreed. This should set zpot_mag to constant h0
! below kbotdyn. It is not necessary to set poles of zpotm3d
! since sub pthreed does not reference the poles of zpotm3d.
!
    do k=mlev0,mlev1
      do j=mlat0,mlat1
        do i=mlon0,mlon1
          zpotm3d(i,j,k) = zpot_mag(i,j,k)
          if (zpot_mag(i,j,k) < h0) zpot_mag(i,j,k)=h0
        enddo
      enddo
    enddo
!
! Set 3d mag fields to zero below kbotdyn:
!
!   ped_mag(:,:,1:kbotdyn-1)    = finit
!   hal_mag(:,:,1:kbotdyn-1)    = finit
!   adotv1_mag(:,:,1:kbotdyn-1) = finit
!   adotv2_mag(:,:,1:kbotdyn-1) = finit

!   call savefld_waccm_switch(adota1_mag(mlon0:mlon1,mlat0:mlat1)   ,'ADOTA1_MAG'  ,1,mlon0,mlon1,mlat0,mlat1)

    do j=mlat0,mlat1
      call outfld('PED_MAG',ped_mag(mlon0:omlon1,j,mlev1:mlev0:-1),omlon1-mlon0+1,j)
      call outfld('HAL_MAG',hal_mag(mlon0:omlon1,j,mlev1:mlev0:-1),omlon1-mlon0+1,j)
      call outfld('ZPOT_MAG',zpot_mag(mlon0:omlon1,j,mlev1:mlev0:-1),omlon1-mlon0+1,j)
      call outfld('ADOTV1_MAG',adotv1_mag(mlon0:omlon1,j,mlev1:mlev0:-1),omlon1-mlon0+1,j)
      call outfld('ADOTV2_MAG',adotv2_mag(mlon0:omlon1,j,mlev1:mlev0:-1),omlon1-mlon0+1,j)
    enddo
!
! Save 3d input fields on geo grid to waccm files (switch to "waccm format"):
!   call savefld_waccm_switch(tn   ,'EDYN_TN'  ,nlev,lon0,lon1,lat0,lat1)
!   call savefld_waccm_switch(un   ,'EDYN_UN'  ,nlev,lon0,lon1,lat0,lat1)
!   call savefld_waccm_switch(vn   ,'EDYN_VN'  ,nlev,lon0,lon1,lat0,lat1)
!   call savefld_waccm_switch(wn   ,'EDYN_WN'  ,nlev,lon0,lon1,lat0,lat1)
!   call savefld_waccm_switch(zpot ,'EDYN_Z'   ,nlev,lon0,lon1,lat0,lat1)
!   call savefld_waccm_switch(ped  ,'EDYN_PED' ,nlev,lon0,lon1,lat0,lat1)
!   call savefld_waccm_switch(hall ,'EDYN_HALL',nlev,lon0,lon1,lat0,lat1)

!   call savefld_waccm_switch(scheight,'EDYN_SCHT'  ,nlev,lon0,lon1,lat0,lat1)
   call savefld_waccm_switch(adotv1  ,'EDYN_ADOTV1',nlev,lon0,lon1,lat0,lat1)
   call savefld_waccm_switch(adotv2  ,'EDYN_ADOTV2',nlev,lon0,lon1,lat0,lat1)
!
! Save 2d geo fields (lon0:lon1,lat0:lat1):
   call savefld_waccm_switch(sini  ,'EDYN_SINI'  ,1,lon0,lon1,lat0,lat1)
   call savefld_waccm_switch(adota1,'EDYN_ADOTA1',1,lon0,lon1,lat0,lat1)
   call savefld_waccm_switch(adota2,'EDYN_ADOTA2',1,lon0,lon1,lat0,lat1)
   call savefld_waccm_switch(a1dta2,'EDYN_A1DTA2',1,lon0,lon1,lat0,lat1)
   call savefld_waccm_switch(be3   ,'EDYN_BE3'   ,1,lon0,lon1,lat0,lat1)

!   endtime = mpi_wtime()
!   time_dynamo_input=time_dynamo_input+(endtime-starttime)
  end subroutine dynamo_input
!-----------------------------------------------------------------------
  subroutine calc_adotv(z,un,vn,wn,adotv1,adotv2,adota1,adota2,&
    a1dta2,be3,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Calculate adotv1,2, adota1,2, a1dta2 and be3.
!
  use edyn_params ,only: r0,h0
  use edyn_geogrid,only: jspole,jnpole
  use getapex,     only: &
    dvec,                & ! (nlonp1,nlat,3,2)
    dddarr,              & ! (nlonp1,nlat)
    be3arr,              & ! (nlonp1,nlat)
    alatm                  ! (nlonp1,0:nlatp1)
!
! Args:
    integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
    real(r8),dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(in) :: &
      z,    & ! geopotential height (cm)
      un,   & ! neutral zonal velocity (cm/s)
      vn      ! neutral meridional velocity (cm/s)
    real(r8),dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(in) :: &
      wn    ! vertical velocity (cm/s)
    real(r8),dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) :: &
      adotv1, adotv2
    real(r8),dimension(lon0:lon1,lat0:lat1),intent(out) :: &
      adota1, adota2, a1dta2, be3
!
! Local:
    integer :: k,i,j
    real(r8) :: r0or,rat,sinalat
    real(r8) :: clm2(lon0:lon1,lat0:lat1)
!
    adotv1 = finit 
    adotv2 = finit
    adota1 = finit
    adota2 = finit
    a1dta2 = finit
    be3 = finit

    do j=lat0,lat1
      if (j==jspole.or.j==jnpole) cycle
      do i=lon0,lon1
        sinalat = sin(alatm(i,j))          ! sin(lam)
        clm2(i,j) = 1._r8-sinalat*sinalat  ! cos^2(lam)
        be3(i,j) = 1.e-9_r8*be3arr(i,j)    ! be3 is in T (be3arr in nT)

        do k=lev0,lev1-1
!
! d_1 = (R_0/R)^1.5
          r0or = r0/(r0 + 0.5_r8*(z(k,i,j)+z(k+1,i,j))-h0)
          rat = 1.e-2_r8*r0or**1.5_r8 ! 1/100 conversion in cm
!
! A_1 dot V = fac( d_1(1) u + d_1(2) v + d_1(3) w
          adotv1(k,i,j) = rat*(        &
            dvec(i,j,1,1)*un(k,i,j)+ &
            dvec(i,j,2,1)*vn(k,i,j)+ &
            dvec(i,j,3,1)*wn(k,i,j))

!
! Note: clm2 is being used here to represent the squared cosine of the
!  quasi-dipole latitude, not of the M(90) latitude, since the wind
!  values are aligned vertically, not along the field line.
!           
          rat = rat*sqrt((4._r8-3._r8*clm2(i,j))/(4._r8-3._r8*r0or*clm2(i,j)))
! 
! A_2 dot V = fac( d_2(1) u + d_2(2) v + d_2(3) w
          adotv2(k,i,j) = rat*(        &
            dvec(i,j,1,2)*un(k,i,j)+ &
            dvec(i,j,2,2)*vn(k,i,j)+ &
            dvec(i,j,3,2)*wn(k,i,j))
        enddo ! k=lev0,lev1-1
! 
! Calculation of adota(n) = d(n)**2/D
!   a1dta2   = (d(1) dot d(2)) /D
!         
        adota1(i,j) = (dvec(i,j,1,1)**2 + dvec(i,j,2,1)**2 + &
                       dvec(i,j,3,1)**2)/dddarr(i,j)
        adota2(i,j) = (dvec(i,j,1,2)**2 + dvec(i,j,2,2)**2 + &
                       dvec(i,j,3,2)**2)/dddarr(i,j)
        a1dta2(i,j) = (dvec(i,j,1,1)*dvec(i,j,1,2) + &
                       dvec(i,j,2,1)*dvec(i,j,2,2) + &
                       dvec(i,j,3,1)*dvec(i,j,3,2))/dddarr(i,j)
      enddo ! i=lon0,lon1

   enddo ! j=lat0,lat1

   call savefld_waccm_switch(adota1   ,'ADOTA1'  ,1,lon0,lon1,lat0,lat1)

  end subroutine calc_adotv
!-----------------------------------------------------------------------
  subroutine alloc_edyn
    use edyn_geogrid,only: nlev
!
! Allocate and initialize arrays for parallel dynamo (module data)
! (called once per run)
!
    integer :: istat
    integer :: mlon00,mlon11,mlat00,mlat11
!
    mlon00=mlon0-1 ; mlon11=mlon1+1
    mlat00=mlat0-1 ; mlat11=mlat1+1
!
! 2d fields on mag subdomains (i,j):
! Certain fields are allocated with halos mlon0-1:mlon1+1,mlat0-1:mlat1+1
!
    allocate(zigm11(mlon00:mlon11,mlat00:mlat11),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: zigm11')
    zigm11 = finit
    allocate(zigmc(mlon00:mlon11,mlat00:mlat11) ,stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: zigmc')
    zigmc = finit
    allocate(zigm2(mlon00:mlon11,mlat00:mlat11) ,stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: zigm2')
    zigm2 = finit
    allocate(zigm22(mlon00:mlon11,mlat00:mlat11),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: zigm22')
    zigm22 = finit
    allocate(rhs(mlon00:mlon11,mlat00:mlat11)   ,stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: rhs')
    rhs = finit
    allocate(rim1(mlon00:mlon11,mlat00:mlat11)  ,stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: rim1')
    rim1 = finit
    allocate(rim2(mlon00:mlon11,mlat00:mlat11)  ,stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: rim2')
    rim2 = finit
    allocate(phimsolv(mlon00:mlon11,mlat00:mlat11),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: phimsolv')
    phimsolv = finit
    allocate(phim2d(mlon00:mlon11,mlat00:mlat11),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: phim2d')
    phim2d = finit
!
! 3d phim and electric field on mag subdomains:
    allocate(phim3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: phim3d')
    phim3d = finit
    allocate(ed13d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: ed13d')
    ed13d = finit 
    allocate(ed23d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: ed23d')
    ed23d = finit
    allocate(ephi3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: ephi3d')
    ephi3d = finit 
    allocate(elam3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: elam3d')
    elam3d = finit
    allocate(emz3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: emz3d')
    emz3d = finit
    allocate(zpotm3d(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: zpotm3d')
    zpotm3d = finit
!
! Fields at mag equator (subdomain longitudes and full column):
!
    allocate(ped_meq(mlon0:mlon1,mlev0:mlev1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: ped_meq')
    ped_meq = finit
    allocate(hal_meq(mlon0:mlon1,mlev0:mlev1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: hal_meq')
    hal_meq = finit
    allocate(adotv1_meq(mlon0:mlon1,mlev0:mlev1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: adotv1_meq')
    adotv1_meq = finit
    allocate(adotv2_meq(mlon0:mlon1,mlev0:mlev1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: adotv2_meq')
    adotv2_meq = finit
    allocate(zpot_meq(mlon0:mlon1,mlev0:mlev1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: zpot_meq')
    zpot_meq = finit
!
! Fields input to mp_mageq (4 fields at full mag subdomain i,j,k):
!
    allocate(fmeq_in(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1,4),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: fmeq_in')
    fmeq_in = finit
!
! Fields output by mp_mageq (4 fields at mag subdomain i,k)
!
    allocate(fmeq_out(mlon0:mlon1,mlev0:mlev1,4),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: fmeq_out')
    fmeq_out = finit
!
! 3d electric field on geographic subdomains, with halos:
!
    allocate(ex(nlev,lon0:lon1,lat0:lat1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: ex')
    allocate(ey(nlev,lon0:lon1,lat0:lat1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: ey')
    allocate(ez(nlev,lon0:lon1,lat0:lat1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: ez')
    ex=finit ; ey=finit ; ez=finit
!
! 3d electric potential on geographic subdomains (k,i,j):
      allocate(poten(nlev,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) call endrun('alloc_pdyn: poten')
      poten = finit
!
! 3d electric potential on geographic subdomains ((i,j,k) regridded from phim3d):
      allocate(phig3d(lon0:lon1,lat0:lat1,mlev0:mlev1),stat=istat)
      if (istat /= 0) call endrun('alloc_pdyn: phig3d')
      phig3d = finit
!
! 2d electric field components on mag grid (these may be input to time3d):
! real(r8),dimension(:,:) :: ed1,ed2 ! (mlon0-1:mlon1+1,mlat0-1:mlat1+1)
!
    allocate(ed1(mlon0-1:mlon1+1,mlat0-1:mlat1+1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: ed1')
    ed1 = finit

    allocate(ed2(mlon0-1:mlon1+1,mlat0-1:mlat1+1),stat=istat)
    if (istat /= 0) call endrun('alloc_edyn: ed2')
    ed2 = finit

  end subroutine alloc_edyn
!-----------------------------------------------------------------------
  subroutine fieldline_integrals
!
! Integrate along magnetic field lines, saving conductances and rims.
!
    use edyn_params,   only: r0,h0,finit,kbotdyn
    use edyn_maggrid,  only: ylatm
!
! Local:
    integer :: i,j,k
    real(r8) :: &
      sinlm,    &  ! sin(lam_m)
      clm2,     &  ! cos^2(lam_m)
      ra,       &  ! (R_E + H_0)/cos^2(lam_m)
      sqomrra,  &  ! sqrt(1/ R_0/R_A) = sin(lam_m)
      sqrra,    &  ! sqrt(R_0/R_A)
      afac,     &  ! 2*sqrt(R_A-R_0) (afac is NaN at mag poles)
      htfac        ! sqrt(R_A -3/4*R_0)

    real(r8) :: rora,del,omdel,sig1,sig2,ue1,ue2
    real(r8),dimension(mlon0:mlon1) :: aam
    real(r8),dimension(mlon0:mlon1,mlev0:mlev1) :: rrm, &
       rtramrm, htfunc, htfunc2
!
! Initialize coefficients:
!
    zigm11 = finit
    zigm22 = finit
    zigm2  = finit
    zigmc  = finit
    rim1   = finit
    rim2   = finit
!
! Subdomain latitude scan for field line integrations:
!
    do j=mlat0,mlat1
!
! Skip poles and equator:
      if (j==1.or.j==nmlat.or.j==nmlath) cycle

      sinlm   = sin(ylatm(j))
      clm2    = 1._r8 - sinlm*sinlm
      ra      = r0/clm2
      sqomrra = sqrt(1._r8-r0/ra)
      sqrra   = sqrt(r0/ra)
      afac    = 2._r8*sqrt(ra-r0)
      htfac   = sqrt(ra-0.75_r8*r0)
      do i=mlon0,mlon1
        aam(i) = afac/abs(sini_mag(i,j))
      enddo
!
! 2*sqrt( h_A - h_0 )/ |sin I_m | w.r to reference height A(h_R)
      do k=kbotdyn,mlev1
        do i=mlon0,mlon1
!
! rr = r0+z-z0 radius of magnetic point
! (Note zpot_mag min was set to h0 in dynamo_inputs)
          rrm(i,k) = r0+zpot_mag(i,j,k)-h0
!
! rtramr = ra-r if +ive, zero otherwise
          rtramrm(i,k) = max(0._r8,ra-rrm(i,k))
          rtramrm(i,k) = sqrt(rtramrm(i,k))
        enddo ! i=mlon0,mlon1
      enddo ! k=kbotdyn,mlev1
!
! Interpolate to midpoints:
! htfunc = factor by which to multiply AAM(I)*d(sqrt(ra-r)) = ds
      do k=kbotdyn,mlev1-1
        do i=mlon0,mlon1
          rrm(i,k)     = 0.5_r8*(rrm(i,k)+rrm(i,k+1))
          rtramrm(i,k) = rtramrm(i,k)-rtramrm(i,k+1)
          htfunc(i,k)  = sqrt(ra-0.75_r8*rrm(i,k))/htfac
          htfunc2(i,k) = htfunc(i,k)**2
        enddo ! i=mlon0,mlon1
      enddo ! k=kbotdyn,mlev1
!
! Compute integrals:
      do k=kbotdyn,mlev1-1
        do i=mlon0,mlon1
!
! (R_E+h)/(R_E+h_A) < 1 -> h_A > h
          rora = min(1._r8,rrm(i,k)/ra)
!
! (lam_m - lam) / lam_m =
! sqrt(1-r_0/r_A)sqrt(r/r_A) - sqrt(r_0/r_A)sqrt(1-r/r_A)
          del = (sqomrra*sqrt(rora)-sqrra*sqrt(1._r8-rora))/abs(ylatm(j))
          omdel = 1._r8 - del
!
! Interpolate conductivities and winds in latitude along field line, assuming
!   linear variation between foot of field line and magnetic equator.
!   (For field lines other than those near the magnetic equator, del is nearly
!   zero, so that the interpolated values are essentially the values for the
!   latitude of the foot of the field line; inaccuracy of the assumption of
!   linear variation is thus unimportant for these field lines.)
!
! Here, mag equator ped_meq, etc. are from mp_mageq, called from dynamo inputs:
          sig1 = omdel*ped_mag(i,j,k) + del*ped_meq(i,k)
          sig2 = omdel*hal_mag(i,j,k) + del*hal_meq(i,k)
          ue1  = omdel*adotv1_mag(i,j,k) + del*adotv1_meq(i,k)
          ue2  = omdel*adotv2_mag(i,j,k) + del*adotv2_meq(i,k)
!
! height varying factors: ds = aam*htfunc
!    d_1^2/D   = 1/htfunc * adota1_mag(i,j)
!    d_2^2/D   = htfunc   * adota2_mag(i,j)
!    d_1*d_2/D = 1        * a1dta2_mag(i,j)
!
! zigm11: int (sigma_p*d_1^2/D) ds : d_1^2/D
! zigm22: int (sigma_p*d_2^2/D) ds : d_2^2/D
!
          zigm11(i,j) = zigm11(i,j) + sig1*rtramrm(i,k)
          zigm22(i,j) = zigm22(i,j) + sig1*rtramrm(i,k)*htfunc2(i,k)
!
! zigmc: int (sigma_p*d_1*d_2/D) ds
! zigm2: int (sigma_h) ds
!
          zigmc(i,j) = zigmc(i,j) + sig1*rtramrm(i,k)*htfunc(i,k)
          zigm2(i,j) = zigm2(i,j) + sig2*rtramrm(i,k)*htfunc(i,k)
!
! rim1: int [sigma_p*d_1^2/D u_e2+(sigma_h-(sigma_p*d_1*d_2)/D) u_e1] ds
! rim2: int [(sigma_h+sigma_p*d_1*d_2/D) u_e2-sigma_p*d_2^2/D u_e1 ] ds
!
          rim1(i,j) = rim1(i,j) + (sig1*adota1_mag(i,j)*ue2 + &
           (sig2 - sig1*a1dta2_mag(i,j))*htfunc(i,k)*ue1)*rtramrm(i,k)
          rim2(i,j) = rim2(i,j) + (sig1*adota2_mag(i,j)*htfunc2(i,k)* &
            ue1 - (sig2 + sig1*a1dta2_mag(i,j))*htfunc(i,k)*ue2)*rtramrm(i,k)
        enddo ! i=mlon0,mlon1
      enddo ! k=kbotdyn,mlev1
!
! Complete calculation and place result in /coefm/ zigm's
!   rim's are in A/m multiply by 1/100 to convert from [cm] to [m]
!
! At this point,
!   zigm11 is int[sig_p*d_1^2/D] ds,   i.e. Sigma_(phi phi)/abs(sin Im)
!   zigm22 is int[sig_p*d_2^2/D] ds,   i.e. Sigma_(lam lam)*abs(sin Im)
!   zigmc  is int[sig_p*d_1*d_2/D] ds, i.e. Sigma_c
!   zigm2  is int[sigma_h] ds,         i.e. Sigma_h
!
!   rim1 is int[(sigma_h-sigma_p*d_1*d_2/D)u_e1 + sigma_p*d_1^2/D u_e2] *A(h_r)*
!                B_e3 ds, i.e.  K_(m phi)^D/abs(sin Im)
!   rim2 is int[(sigma_h+sigma_p*d_1*d_2/D)u_e2 - sigma_p*d_2^2/D u_e1] *A(h_r)*
!                B_e3 ds, K_(m lam)^D ( minus in northern hemisphere
!   Change sign of RIM(2) in S. hemisphere to be compatible with transf
! At this point, rim2 is +-K_(m lam)^D
!
      do i = mlon0,mlon1
        zigm11(i,j) = 1.e-2_r8*zigm11(i,j)*aam(i)*adota1_mag(i,j)
        zigm22(i,j) = 1.e-2_r8*zigm22(i,j)*aam(i)*adota2_mag(i,j)
        zigmc(i,j)  = 1.e-2_r8*zigmc (i,j)*aam(i)*a1dta2_mag(i,j)
        zigm2(i,j)  = 1.e-2_r8*zigm2 (i,j)*aam(i)
        rim1(i,j)   = 1.e-2_r8*rim1(i,j)*aam(i)*be3_mag(i,j)
        rim2(i,j)   = 1.e-2_r8*rim2(i,j)*aam(i)*be3_mag(i,j)
      enddo ! i = 1,nmlon
    enddo ! j=mlat0,mlat1 (without poles)

!   call savefld_waccm_switch(adota1_mag(mlon0:mlon1,mlat0:mlat1)   ,'adota1_mag_a'  ,1,mlon0,mlon1,mlat0,mlat1)

!   call savefld_waccm_switch(zigm11(mlon0:mlon1,mlat0:mlat1)   ,'ZIGM11_a'  ,1,mlon0,mlon1,mlat0,mlat1)

  end subroutine fieldline_integrals
!-----------------------------------------------------------------------
  subroutine complete_integrals
  use edyn_mpi,only: mlat0,mlat1,mlon0,mlon1,mp_mageq_jpm1,mp_magpole_2d,&
    mp_mag_foldhem,mp_mag_periodic_f2d
  use edyn_maggrid,only: rcos0s,dt1dts
!
! Field line integrals for each hemisphere have been calculated in
!   mag subdomains. Now, complete these arrays with equator and polar
!   values, and sum the integrals from the 2 hemispheres.
! This is done by obtaining global mag fields via mpi exchange
!   of mag subdomains, completing the global fields, and updating
!   the subdomains.
! The 6 2d fields are: zigm11,zigm22,zigmc,zigm2,rim1,rim2
!
! Local:
    integer :: i,j,ii,lonend
    real(r8) :: fmsub(mlon0:mlon1,mlat0:mlat1,nf2d)
    real(r8) :: corfac
    real(r8),parameter :: r8_nmlon = dble(nmlon)
!
! If do_integ==.false. (meaning use_time3d_integ=.true.), then these were passed in 
! from time3d (in this case, dynamo did not call fieldline_integrals). Otherwise, 
! they were calculated by this module, in sub fieldline_integrals and this routine.
! 

! 
! For equatorial values, we need latitudes eq+1 and eq-1:
! Local feq_jpm1(nmlonp1,2,6) is returned by mp_mageq_jpm1,
!   where the 2 dim contains lats nmlath-1, nmlath+1. These 
!   are global in lon, even tho each subd uses only its own i's.
! These mag equator values do not show up on plots because
!   of the small factor .06 and .125.
! The 6 2d fields are: zigm11,zigm22,zigmc,zigm2,rim1,rim2
! Must specify mlon0:mlon1,mlat0:mlat1 because these fields
!   were allocated to include single-point halos (these calls
!   exclude the halo points):
!
    fmsub(:,:,1) = zigm11(mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,2) = zigm22(mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,3) = zigmc (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,4) = zigm2 (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,5) = rim1  (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,6) = rim2  (mlon0:mlon1,mlat0:mlat1)

    call mp_mageq_jpm1(fmsub,mlon0,mlon1,mlat0,mlat1,nmlonp1,feq_jpm1,nf2d)
!
! From sub fieldline_integrals:
!   zigm11 is int[sig_p*d_1^2/D] ds,   i.e. Sigma_(phi phi)/abs(sin Im)
!   zigm22 is int[sig_p*d_2^2/D] ds,   i.e. Sigma_(lam lam)*abs(sin Im)
!   zigmc  is int[sig_p*d_1*d_2/D] ds, i.e. Sigma_c
!   zigm2  is int[sigma_h] ds,         i.e. Sigma_h
!
!   rim1 is int[(sigma_h-sigma_p*d_1*d_2/D)u_e1 + sigma_p*d_1^2/D u_e2] *A(h_r)*
!                B_e3 ds, i.e.  K_(m phi)^D/abs(sin Im)
!   rim2 is int[(sigma_h+sigma_p*d_1*d_2/D)u_e2 - sigma_p*d_2^2/D u_e1] *A(h_r)*
!                B_e3 ds, K_(m lam)^D ( minus in northern hemisphere
!   Change sign of RIM(2) in S. hemisphere to be compatible with transf
! At this point, rim2 is +-K_(m lam)^D
!
! Equatorial values:
! Assume that quantities primarily dependent on Pedersen conductivity
!   have field-line integrals 1/4 as large as the averages for next-higher
!   field lines; quantities primarily dependent on Hall conductivity
!   have field-line integrals 0.12 as large as the averages for next-higher
!   field lines.  Exact values chosen should not be important for potential
!   calculation, as long as they are physically reasonable and not too
!   different from adjacent values.
!
    do j=mlat0,mlat1
      if (j==nmlath) then ! mag equator

        zigm11(mlon0:mlon1,j) = .125_r8*(feq_jpm1(mlon0:mlon1,1,1)+ &
                                         feq_jpm1(mlon0:mlon1,2,1))
        zigm22(mlon0:mlon1,j) = .125_r8*(feq_jpm1(mlon0:mlon1,1,2)+ &
                                         feq_jpm1(mlon0:mlon1,2,2))
        zigmc (mlon0:mlon1,j) = .125_r8*(feq_jpm1(mlon0:mlon1,1,3)+ &
                                         feq_jpm1(mlon0:mlon1,2,3))
        zigm2 (mlon0:mlon1,j) = .060_r8*(feq_jpm1(mlon0:mlon1,1,4)+ &
                                         feq_jpm1(mlon0:mlon1,2,4))
        rim1  (mlon0:mlon1,j) = .060_r8*(feq_jpm1(mlon0:mlon1,1,5)+ &
                                         feq_jpm1(mlon0:mlon1,2,5))
        rim2  (mlon0:mlon1,j) = .060_r8*(feq_jpm1(mlon0:mlon1,1,6)+ &
                                         feq_jpm1(mlon0:mlon1,2,6))
!
! Include the boundary condition at the equator eq.(5.30) in
! Richmond (1995) Ionospheric Electrodynamics use. Mag. Apex Coord.
!   J.Geomag.Geoelectr. 47,191-212
! Sig_phiphi/abs(sin Im) = 0.5*Sig_cowling/abs(sin Im)
!   = 0.5/abs(sin Im)*(Sig_phiphi - Sig_philam*sig_lamphi/Sig_lamlam)
!   = 0.5/abs(sin Im)*(Sig_phiphi + (Sig_h-sig_c)*(Sig_h+sig_c)/Sig_lamlam)
!  rim(1) / |sin I_m| = I_1 = R/2*(K_mphi - Sig_philam/Sig_lamlam*K_mlam)
!
        do i=mlon0,mlon1
          zigm11(i,j) = zigm11(i,j)+(zigm2(i,j)-zigmc(i,j))* &
            (zigm2(i,j)+zigmc(i,j))/zigm22(i,j)
          rim1(i,j) = rim1(i,j) - (zigm2(i,j)-zigmc(i,j))/ &
            zigm22(i,j)*rim2(i,j)
        enddo ! i=mlon0,mlon1
      endif ! j at equator
    enddo ! j=mlat0,mlat1
!
    do j=mlat0,mlat1
      call outfld('EDYN_ZIGM11_PED',zigm11(mlon0:omlon1,j),omlon1-mlon0+1,j)
      call outfld('EDYN_ZIGM2_HAL',zigm2(mlon0:omlon1,j),omlon1-mlon0+1,j)
    enddo

!
! Using notation of Richmond (1995) on right-hand side below:
! Sigma_(phi phi) = zigm11*abs(sin I_m)
! Sigma_(lam lam) = zigm22/abs(sin I_m)
! Sigma_(phi lam) = +-(zigm2-zigmc)
! Sigma_(lam phi) = -+(zigm2+zigmc)
! K_(m phi)^D     = rim(1)*abs(sin I_m)
! K_(m lam)^D     = +-rim(2)
!
! Transforming PDE from original apex (theta_a) to new apex grid (theta_0)
!   which is equally spaced in magnetic latitude
! SCALE quantities to modified (0) magnetic latitude system, multiplying or
!   dividing by abs(sin I_m) [inverse contained in DT1DTS] as necessary.
! Sign of K_(m lam)^D in southern hemisphere remains reversed.
! for the mixed terms the transformation from the integration and differentiation
! canceled out (zigmc, zigm2)
! DT1DTS : d theta_0/ d theta_a / abs(sin I_m)
! RCOS0S : cos(theta_0)/ cos(theta_a)
!
! corfac: abs(I_m)*d theta_a/d theta_0 * cos(theta_0)/ cos(theta_a)
! zigm11: abs(I_m)*d theta_a/d theta_0 * cos(theta_0)/ cos(theta_a)
! zigm22: 1/abs(I_m)*d theta_0/d theta_a * cos(theta_a)/ cos(theta_0)
! rim(1): abs(I_m)*d theta_a/d theta_0
! rim(2): cos(theta_a)/ cos(theta_0)
!
    do j=mlat0,mlat1
      if (j==1.or.j==nmlat) cycle   ! skip poles
      corfac = rcos0s(j)/dt1dts(j)
      do i=mlon0,mlon1
        zigm11(i,j) = zigm11(i,j)*corfac
        zigm22(i,j) = zigm22(i,j)/corfac
        rim1(i,j) = rim1(i,j)/dt1dts(j)
        rim2(i,j) = rim2(i,j)/rcos0s(j)
      enddo
    enddo
! 
! For polar values, we need south pole plus 1 and 2 (j==2,3),
!   and north pole minus 1 and 2 (j==nmlat-1,nmlat-2). These 
!   are returned by sub mp_magpole_jpm2 (mpi.F):
! Must specify (mlon0:mlon1,mlat0:mlat1) because zigmxx and rims
!   are allocated to include halo cells. 
! 
    fmsub(:,:,1) = zigm11(mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,2) = zigm22(mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,3) = zigmc (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,4) = zigm2 (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,5) = rim1  (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,6) = rim2  (mlon0:mlon1,mlat0:mlat1)
! 
! mp_magpole_2d returns fpole_jpm2(nmlonp1,1->4,nf) as:
!   1: j = 2       (spole+1)
!   2: j = 3       (spole+2)
!   3: j = nmlat-1 (npole-1)
!   4: j = nmlat-2 (npole-2)
!
    call mp_magpole_2d(fmsub,mlon0,mlon1,mlat0,mlat1,nmlonp1, &
      1,nmlat,fpole_jpm2,nf2d)
!
! the PDE is divided by 1/ DT0DTS
! Sigma_(phi phi) = zigm11/ rcos0s * dt0dts
! Sigma_(lam lam) = zigm22 * rcos0s / dt0dts
! Sigma_(phi lam) = +-(zigm2-zigmc)
! Sigma_(lam phi) = -+(zigm2+zigmc)
! K_(m phi)^D     =   rim(1) * dt0dts
! K_(m lam)^D     = +-rim(2) * rcos0s
!
! Compute polar values for the conductances, 4th order interpolation:
!
    do j=mlat0,mlat1
!
! South pole:
      if (j==1) then ! south pole (use fpole_jpm2(nmlon,1->2,nf)
        zigm11(mlon0,j)=(4._r8*                         &
          dot_product(unitvm,fpole_jpm2(1:nmlon,1,1))-  &
          dot_product(unitvm,fpole_jpm2(1:nmlon,2,1)))/ &
          (3._r8*r8_nmlon)
        zigm22(mlon0,j)=(4._r8*                         &
          dot_product(unitvm,fpole_jpm2(1:nmlon,1,2))-  &
          dot_product(unitvm,fpole_jpm2(1:nmlon,2,2)))/ &
          (3._r8*r8_nmlon)
        zigmc (mlon0,j)=(4._r8*                         &
          dot_product(unitvm,fpole_jpm2(1:nmlon,1,3))-  &
          dot_product(unitvm,fpole_jpm2(1:nmlon,2,3)))/ &
          (3._r8*r8_nmlon)
        zigm2 (mlon0,j)=(4._r8*                         &
          dot_product(unitvm,fpole_jpm2(1:nmlon,1,4))-  &
          dot_product(unitvm,fpole_jpm2(1:nmlon,2,4)))/ &
          (3._r8*r8_nmlon)
!
! Extend south pole over longitude:
        do i=mlon0+1,mlon1
          zigm11(i,j) = zigm11(mlon0,j)
          zigm22(i,j) = zigm22(mlon0,j)
          zigmc (i,j) = zigmc (mlon0,j)
          zigm2 (i,j) = zigm2 (mlon0,j)
        enddo ! i=mlon0,mlon1
! 
! RHS vector (I_1,I_2): average over south pole:
! (use fpole_jpm2(i,1,nf), i.e. j==2, and lons across the pole)
        lonend = mlon1
        if (mlon1==nmlonp1) lonend = mlon1-1
        do i=mlon0,lonend
          ii = 1+mod(i-1+nmlon/2,nmlon)
          rim1(i,j) = 0.5_r8*(fpole_jpm2(i,1,5)-fpole_jpm2(ii,1,5))
          rim2(i,j) = 0.5_r8*(fpole_jpm2(i,1,6)-fpole_jpm2(ii,1,6))
        enddo
! 
! North pole:
      elseif (j==nmlat) then ! north pole (use fpole_jpm2(nmlon,3->4,1,nf)
        zigm11(mlon0,j)=(4._r8*                          &
          dot_product(unitvm,fpole_jpm2(1:nmlon,3,1))-   &
          dot_product(unitvm,fpole_jpm2(1:nmlon,4,1)))/  &
          (3._r8*r8_nmlon)
        zigm22(mlon0,j)=(4._r8*                          &
          dot_product(unitvm,fpole_jpm2(1:nmlon,3,2))-   &
          dot_product(unitvm,fpole_jpm2(1:nmlon,4,2)))/  &
          (3._r8*r8_nmlon)
        zigmc (mlon0,j)=(4._r8*                          &
          dot_product(unitvm,fpole_jpm2(1:nmlon,3,3))-   &
          dot_product(unitvm,fpole_jpm2(1:nmlon,4,3)))/  &
          (3._r8*r8_nmlon)
        zigm2 (mlon0,j)=(4._r8*                          &
          dot_product(unitvm,fpole_jpm2(1:nmlon,3,4))-   &
          dot_product(unitvm,fpole_jpm2(1:nmlon,4,4)))/  &
          (3._r8*r8_nmlon)
! 
! Extend north pole over longitude:
        do i=mlon0+1,mlon1
          zigm11(i,j) = zigm11(mlon0,j)
          zigm22(i,j) = zigm22(mlon0,j)
          zigmc (i,j) = zigmc (mlon0,j)
          zigm2 (i,j) = zigm2 (mlon0,j)
        enddo ! i=mlon0,mlon1
! 
! RHS vector (I_1,I_2): average over north pole:
! (use fpole_jpm2(i,3,nf), i.e. j==nmlat-1, and lons across the pole)
        lonend = mlon1
        if (mlon1==nmlonp1) lonend = mlon1-1
        do i=mlon0,lonend
          ii = 1+mod(i-1+nmlon/2,nmlon)
          rim1(i,j) = 0.5_r8*(fpole_jpm2(i,3,5)-fpole_jpm2(ii,3,5))
          rim2(i,j) = 0.5_r8*(fpole_jpm2(i,3,6)-fpole_jpm2(ii,3,6))
        enddo
      endif ! south or north pole
    enddo ! j=mlat0,mlat1
!
! Fold south hemisphere over onto north, and set periodic points:
    fmsub(:,:,1) = zigm11(mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,2) = zigm22(mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,3) = zigmc (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,4) = zigm2 (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,5) = rim1  (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,6) = rim2  (mlon0:mlon1,mlat0:mlat1)

    call mp_mag_foldhem(fmsub,mlon0,mlon1,mlat0,mlat1,nf2d)
    call mp_mag_periodic_f2d(fmsub,mlon0,mlon1,mlat0,mlat1,nf2d)

    zigm11(mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,1)
    zigm22(mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,2)
    zigmc (mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,3)
    zigm2 (mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,4)
    rim1  (mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,5)
    rim2  (mlon0:mlon1,mlat0:mlat1) = fmsub(:,:,6)
!
! Reverse sign of zigmc in northern hemisphere.
    do j=mlat0,mlat1
      if (j >= nmlath) then
        zigmc(mlon0:mlon1,j) = -zigmc(mlon0:mlon1,j)
      endif
      call outfld('EDYN_RIM1',rim1(mlon0:omlon1,j),omlon1-mlon0+1,j)
      call outfld('EDYN_RIM2',rim2(mlon0:omlon1,j),omlon1-mlon0+1,j)
    enddo

    if (debug.and.masterproc) then 
      write(iulog,"('complete_integrals: nstep=',i4)") nstep
      write(iulog,"('  zigm11 min,max=',2e12.4)") &
   	minval(zigm11(mlon0:mlon1,mlat0:mlat1)),maxval(zigm11(mlon0:mlon1,mlat0:mlat1))
      write(iulog,"('  zigm22 min,max=',2e12.4)") &
   	minval(zigm22(mlon0:mlon1,mlat0:mlat1)),maxval(zigm22(mlon0:mlon1,mlat0:mlat1))
      write(iulog,"('  zigmc  min,max=',2e12.4)") &
   	minval(zigmc (mlon0:mlon1,mlat0:mlat1)),maxval(zigmc (mlon0:mlon1,mlat0:mlat1))
      write(iulog,"('  zigm2  min,max=',2e12.4)") &
   	minval(zigm2 (mlon0:mlon1,mlat0:mlat1)),maxval(zigm2 (mlon0:mlon1,mlat0:mlat1))
      write(iulog,"('  rim1   min,max=',2e12.4)") &
   	minval(rim1  (mlon0:mlon1,mlat0:mlat1)),maxval(rim1  (mlon0:mlon1,mlat0:mlat1))
      write(iulog,"('  rim2   min,max=',2e12.4)") &
   	minval(rim2  (mlon0:mlon1,mlat0:mlat1)),maxval(rim2  (mlon0:mlon1,mlat0:mlat1))
    endif

    call savefld_waccm_switch(zigm11(mlon0:mlon1,mlat0:mlat1)   ,'EDYN_ZIGM11'  ,1,mlon0,mlon1,mlat0,mlat1)
    call savefld_waccm_switch(zigm22(mlon0:mlon1,mlat0:mlat1)   ,'EDYN_ZIGM22'  ,1,mlon0,mlon1,mlat0,mlat1)
    call savefld_waccm_switch(zigmc (mlon0:mlon1,mlat0:mlat1)   ,'EDYN_ZIGMC'   ,1,mlon0,mlon1,mlat0,mlat1)
    call savefld_waccm_switch(zigm2 (mlon0:mlon1,mlat0:mlat1)   ,'EDYN_ZIGM2'   ,1,mlon0,mlon1,mlat0,mlat1)
    call savefld_waccm_switch(rim1  (mlon0:mlon1,mlat0:mlat1)   ,'EDYN_RIM1'    ,1,mlon0,mlon1,mlat0,mlat1)
    call savefld_waccm_switch(rim2  (mlon0:mlon1,mlat0:mlat1)   ,'EDYN_RIM2'    ,1,mlon0,mlon1,mlat0,mlat1)

  end subroutine complete_integrals
!-----------------------------------------------------------------------
  subroutine rhspde
    use edyn_params   ,only: pi_dyn,r0
    use edyn_maggrid  ,only: dlatm,dlonm,rcos0s,dt1dts
!
! Calculate right-hand side from rim1,2 on mag subdomains.
! Use global longitude arrays for poles and equator obtained
!   by sub complete_integrals.
!
! Local:
    integer :: j,i
    real(r8),dimension(nmlat) :: tint1
    real(r8) :: &
      rim2_npm1(nmlonp1), & ! global rim2 at nmlat-1
      rim2_eqp1(nmlonp1), & ! global rim2 at meq+1
      rim2_meq (nmlonp1), & ! global rim2 at mag eq
      rim2_tmp (nmlonp1), & ! temp array
      rim1_meq (nmlonp1), & ! global rim1 at mag eq
      zigm2_meq(nmlonp1), & ! needed for rim1_meq
      zigmc_meq(nmlonp1), & ! needed for rim1_meq
      zigm22_meq(nmlonp1)   ! needed for rim1_meq
    real(r8),parameter :: r8_nmlon = dble(nmlon)

    do j=1,nmlat
      tint1(j) = cos(-pi_dyn/2._r8+(j-1)*dlatm)
    enddo
!
! Init rhs subdomains:
    rhs(:,:) = finit
!
! Will need rim2 at npole-1 and mag equator:
! rim2_npm1: global rim2 at nmlat-1:
    rim2_npm1(:) = fpole_jpm2(:,3,6)+fpole_jpm2(:,1,6)
!
! rim2_meq: global rim2 at mag equator:
    rim2_meq(:) = .060_r8*(feq_jpm1(:,1,6)+feq_jpm1(:,2,6))
    rim2_meq(:) = rim2_meq(:)/rcos0s(nmlath)
    rim2_meq(:) = rim2_meq(:)*2._r8 ! fold eq on itself
!
! Perform differentiation of rim(2) w.r.t. lam_0:
!  +/- (d [ K_(m lam)^D * cos(lam_m)]/ d lam_0 ) /cos ( lam_0) =
!  + (d [ K_(m lam)^D * cos(lam_m)]/ d lam_m ) /cos ( lam_m) / (RCOS0S*DT0DTS) =
!  +/- (d [ K_(m lam)^D(0) * cos(lam_0)]/ d lam_0 ) /cos ( lam_0) =
!
! Lat scan to define rhs subdomains:
    do j=mlat0,mlat1
!
! North Pole (redundant in longitude):
      if (j == nmlat) then       ! north pole
        do i=mlon0,mlon1
          rhs(i,j) = -2._r8/r8_nmlon*dot_product(unitvm,rim2_npm1(1:nmlon))/ &
            tint1(nmlat-1)
        enddo
!
! Include the boundary condition at the equator.
! rhs(equator)/R = d (K_mphi^DT(0) - sig_philam/sig_lamlam*K_mlam^DT(0)) / d phi_m
!                + d (cos lam_0 * K_mlam^DT(0))/ d lam_0
! from Cicely's notes:
! I_1 = 0.5*(K_(m phi)^DT(0) - Sig_(phi lam)/Sig_(lam lam)*K_(ml am)^DT(0))
! I_2 = K_(m lam)^DT(0)
! differentiate
! rhs = (I_1(i+1/2,j)-I_1(i-1/2,j))/dlonm +
!       (2*cos(lam_0)_(j+1/2)*I_2(i,j+1/2))/dlat_0
! (first calc global mag equator as in complete_integrals)
!
      elseif (j == nmlath) then  ! mag equator
        rim2_eqp1(:) = feq_jpm1(:,2,6)/rcos0s(j-1)+ &
                       feq_jpm1(:,1,6)/rcos0s(j+1)
        zigm22_meq(:)  = .125_r8*(feq_jpm1(:,1,2)+feq_jpm1(:,2,2))
        zigmc_meq (:)  = .125_r8*(feq_jpm1(:,1,3)+feq_jpm1(:,2,3))
        zigm2_meq (:)  = .060_r8*(feq_jpm1(:,1,4)+feq_jpm1(:,2,4))
        rim1_meq  (:)  = .060_r8*(feq_jpm1(:,1,5)+feq_jpm1(:,2,5))
        rim2_tmp  (:)  = .060_r8*(feq_jpm1(:,1,6)+feq_jpm1(:,2,6))
        do i=1,nmlonp1
          rim1_meq(i) = rim1_meq(i) - (zigm2_meq(i)-zigmc_meq(i))/ &
            zigm22_meq(i)*rim2_tmp(i)
        enddo
        rim1_meq(:) = rim1_meq(:)/dt1dts(j)
        rim1_meq(:) = rim1_meq(:)*2._r8 ! fold eq on itself

        do i=mlon0,mlon1
          if (i==1) then               ! western most lon
            rhs(i,j) = 0.5_r8/dlonm*(rim1(i+1,j)-rim1_meq(nmlon))
            rhs(i,j) = rhs(i,j)+1._r8/dlatm*(tint1(j)*rim2(i,j)+ &
              tint1(j+1)*rim2_eqp1(i))
          elseif (i==nmlonp1) then       ! eastern most lon
            rhs(i,j) = 0.5_r8/dlonm*(rim1_meq(1)-rim1(i-1,j))
            rhs(i,j) = rhs(i,j)+1._r8/dlatm*(tint1(j)*rim2(i,j)+ &
              tint1(j+1)*rim2_eqp1(i))
          else                         ! body of i subdomain
! Note that rim1 halos were set before calling this subroutine.
            rhs(i,j) = 0.5_r8/dlonm*(rim1(i+1,j)-rim1(i-1,j))
            rhs(i,j) = rhs(i,j)+1._r8/dlatm*(tint1(j)*rim2(i,j)+ &
               tint1(j+1)*rim2_eqp1(i))
            endif
          enddo ! i=mlon0,mlon1
!
! North hemisphere (not npole and not equator):
! (allow south hemisphere to remain 0)
! (use rim1 instead of tint33)
!
      elseif (j > nmlath) then ! north hem only (excluding npole)
        do i=mlon0,mlon1
          rhs(i,j) = 1._r8/(dlonm*tint1(j))*0.5_r8*(rim1(i+1,j)-rim1(i-1,j))
          if (j == nmlath+1) then
            rhs(i,j) = rhs(i,j)+1._r8/(dlatm*tint1(j))* &
             0.5_r8*(rim2(i,j+1)*tint1(j+1)-rim2_meq(i)*tint1(j-1))
          else
            rhs(i,j) = rhs(i,j)+1._r8/(dlatm*tint1(j))* &
             0.5_r8*(rim2(i,j+1)*tint1(j+1)-rim2(i,j-1)*tint1(j-1))
          endif
        enddo
      endif ! at poles or equator
    enddo ! j=mlat0,mlat1
!
! scale (multiply by earth radius in meter  = R0*1.E-2)
![( d K_(m phi)^D / d phi /(cos(theta_m)?) +
! (d [ K_(m lam)^D * cos(lam_m)]/ d lam_m ) /cos ( lam_m) ] * R / (RCOS0S*DT0DTS)
! ~ J_(Mr)*r^2*cos(theta_m)/cos(theta_0)/DT0DTS
! theta_m = theta_s
!
    do j=mlat0,mlat1
      do i=mlon0,mlon1
        rhs(i,j) = rhs(i,j)*r0*1.e-2_r8
      enddo
    enddo

  end subroutine rhspde
!-----------------------------------------------------------------------
  subroutine gather_edyn
!
! Gather needed global arrays to root task, so it can finish non-parallel
! part of dynamo (beginning after sub rhspde) as in original code
!
    use edyn_mpi,  only: mp_gather_edyn
    use edyn_solve,only: &  ! (nmlonp1,nmlat)
      zigm11_glb        ,&
      zigm22_glb        ,&
      zigmc_glb         ,&
      zigm2_glb         ,&
      rhs_glb         
    use edyn_solve ,only: rim_glb ! pde solver output (nmlonp1,nmlat,2)
!
! Local:
! 7 fields to gather: zigm11,zigm22,zigmc,zigm2,rim1,rim2,rhs
!
    integer,parameter :: nf = 7
    real(r8) :: fmsub(mlon0:mlon1,mlat0:mlat1,nf)
    real(r8) :: fmglb(nmlonp1,nmlat,nf)
    real(r8) :: rhs_nhem(nmlonp1,nmlat)
    integer :: i,j,jj
!
! These calls exclude halo points in zigm11, etc.
!
    fmsub(:,:,1) = zigm11(mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,2) = zigm22(mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,3) = zigmc (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,4) = zigm2 (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,5) = rim1  (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,6) = rim2  (mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,7) = rhs   (mlon0:mlon1,mlat0:mlat1)

    call mp_gather_edyn(fmsub,mlon0,mlon1,mlat0,mlat1, &
      fmglb,nmlonp1,nmlat,nf)
!
! Now root task can take over, and work with global arrays:
!
    if (mytid==0) then
      zigm11_glb(:,:) = fmglb(:,:,1)
      zigm22_glb(:,:) = fmglb(:,:,2)
      zigmc_glb(:,:)  = fmglb(:,:,3)
      zigm2_glb(:,:)  = fmglb(:,:,4)
      rim_glb(:,:,1)  = fmglb(:,:,5)
      rim_glb(:,:,2)  = fmglb(:,:,6)
      rhs_nhem(:,:)   = fmglb(:,:,7)
!
! Transfer local global rhs_nhem (from sub rhspde) to rhs_glb,
! so the latter has data in south hemisphere.
!
      rhs_glb= 0._r8 ! init
      do j=1,nmlat
!
! Transfer north pole to equator:
        if (j == nmlat) then
          do i=1,nmlonp1
            rhs_glb(i,nmlath) = rhs_nhem(i,j)
          enddo
! Transfer equator to south pole:
        elseif (j == nmlath) then
          do i=1,nmlonp1
            rhs_glb(i,1) = rhs_nhem(i,j)
          enddo
! Transfer north hem to south hem:
        elseif (j > nmlath) then ! 50 -> 96
          jj = j-nmlath+1        !  2 -> 48
          do i=1,nmlonp1
            rhs_glb(i,jj) = rhs_nhem(i,j)
          enddo
        endif
      enddo ! j=mlat0,mlat1

    endif ! mytid==0
  end subroutine gather_edyn
!-----------------------------------------------------------------------
  subroutine highlat_poten
  use edyn_solve,only: &
    phihm             ,& ! high-latitude potential (nmlonp1,nmlat)
    pfrac                ! NH fraction of potential (nmlonp1,nmlat0)
!
! Global PDE solution rim_glb(:,:,1) has been scattered to mag subdomains
! in rim1, and halos set (this overwrites previous rim1 from fieldline
! integrations).  Now add high latitude potential from empirical model
! (heelis or weimer), defining phim2d on mag subdomains. After this pthreed
! will expand phim2d to phim3d.
!
! Input:  rim1(mag subdomains)   ! Solution from mudpack solver (nhem only)
!         pfrac(nmlonp1,nmlat0)  ! NH fraction of potential
!         phihm(nmlonp1,nmlat)   ! potential in magnetic
! Output: phim2d(mag subdomains) ! solution with phihm in both nhem and shem
! Both rim1 and phim2d are dimensioned (mlon00:mlon11,mlat00:mlat11)
!
! Both phihm and pfrac have been set by either heelis or weimer.
! phihm is on 2d global mag grid, pfrac is in north hemisphere only
!
! Local:
    logical,parameter :: mod_heelis = .false.   ! true == modified
    integer :: i,j,jn,js
    real(r8) :: fac
!
! Add empirical model potential at high latitude:
!
    fac = 1.0_r8
    if (mod_heelis) fac = 0._r8  ! modified heelis
    do j=mlat0,mlat1
      if (j > nmlath) cycle   ! south only (including equator)
      jn = nmlat-j+1
      js = nmlath-j+1
      do i=mlon0,mlon1
        phim2d(i,j) = rim1(i,j)+fac*(1._r8-pfrac(i,js))*(phihm(i,j)- &
          phihm(i,jn))
      enddo
    enddo

    do j=mlat0,mlat1
      if (j <= nmlath) cycle ! north only (excluding equator)
      do i=mlon0,mlon1
        phim2d(i,j) = rim1(i,j)
      enddo
    enddo

    do j=mlat0,mlat1
      call outfld('PHIM2D',phim2d(mlon0:omlon1,j),omlon1-mlon0+1,j)
    enddo

  end subroutine highlat_poten
!-----------------------------------------------------------------------
  subroutine pthreed
!
! phim2d is now 2d electric potential solution on mag subdomains,
! with high-latitude potential added from empirical model (see subs
! heelis and highlat_poten), and mag halos set. Now expand phim2d in 
! vertical, defining phim3d. Also calculate electric field ed13d, ed23d 
! for later current calculations, and ephi3d, elam3d and emz3d for conversion
! to geographic grid (sub pefield), and subsequent calculation of ion drifts 
! by sub ionvel (not in edynamo).
!
  use edyn_params ,only: re,pi_dyn,r0,kbotdyn
  use edyn_maggrid,only: ylatm,dlatm,dlonm,rcos0s,dt1dts,dt0dts,table
  use edyn_mpi    ,only: &
    mp_mag_halos        ,& 
    mp_magpole_2d       ,&
    mp_mageq_jpm3       ,&
    mp_mag_jslot        ,&
    mp_magpoles         ,&
    mp_mag_periodic_f2d ,&
    ixfind
!
! Local:
    real(r8),parameter :: eps = 1.e-10_r8, unitvm(nmlon)=1._r8
    integer,parameter :: mxneed=nmlat+2
    integer :: i,j,k,n,mlon00,mlon11,mlat00,mlat11
    real(r8) :: csth0,cosltm,sym,pi,phims,phimn,real8
    real(r8),dimension(nmlonp1) :: thetam,pslot,qslot
    integer,dimension(nmlonp1) :: islot,jslot,ip1f,ip2f,ip3f

!   real(r8),dimension(mlon0-1:mlon1+1,mlat0-1:mlat1+1) :: ed1,ed2

    real(r8),dimension(mlon0-1:mlon1+1,mlat0-1:mlat1+1) :: ephi,elam
    real(r8) :: fpole2d_jpm2(nmlonp1,4,4) ! global lons at S pole+1,2 and N pole-1,2
    real(r8) :: fpoles(nmlonp1,2,1) ! global lons at poles (1 field only)
    real(r8) :: fmsub(mlon0:mlon1,mlat0:mlat1,4)
    real(r8) :: fmsub1(mlon0-1:mlon1+1,mlat0-1:mlat1+1,5)
    real(r8) :: feq_jpm3(nmlonp1,-3:3,1) ! global lons at equator +/- 3
    integer :: jneed(mxneed) ! lats needed from other tasks for interp
    integer :: njneed,icount
    real(r8),dimension(mlon0-1:mlon1+1,mxneed) :: &
      phineed,   & ! phim2d at needed latitudes
      ed1need,   & ! ed1 at needed latitudes
      ed2need,   & ! ed2 at needed latitudes
      ephineed,  & ! ephi at needed latitudes
      elamneed     ! elam at needed latitudes
    real(r8),dimension(mlon0-1:mlon1+1,mxneed,5) :: fmneed
    real(r8) :: phi0j0,phi1j0,phi0j1,phi1j1
    real(r8) :: ed1i0j0,ed1i1j0,ed1i0j1,ed1i1j1
    real(r8) :: ed2i0j0,ed2i1j0,ed2i0j1,ed2i1j1
    real(r8) :: ephi0j0,ephi1j0,ephi0j1,ephi1j1
    real(r8) :: elam0j0,elam1j0,elam0j1,elam1j1
    real(r8) :: fac_elam
!
    pi = pi_dyn
    mlon00=mlon0-1 ; mlon11=mlon1+1
    mlat00=mlat0-1 ; mlat11=mlat1+1
!
! Calculate ed1,ed2 components of electric field:
! phim2d has halos set, so when mlon0==1, i-1 should wrap
!   to value at i==nmlon, and when mlon1==nmlonp1, i+1 should
!   wrap to value at i==2.
!
    ed1 = 0._r8
    ephi = 0._r8
    do j=mlat0,mlat1
      if (j==1.or.j==nmlat) cycle
      csth0 = cos(-pi/2._r8+(j-1)*dlatm)
      do i=mlon0,mlon1
        ed1(i,j) = -(phim2d(i+1,j)-phim2d(i-1,j))/(2._r8*dlonm*csth0)* &
          rcos0s(j)/(r0*1.e-2_r8)
        ephi(i,j) = ed1(i,j)*(r0*1.e-2_r8)
      enddo ! i=mlon0,mlon1
    enddo ! j=mlat0,mlat1
!
! Southern hemisphere (excluding equator):
    ed2 = 0._r8
    elam = 0._r8
    do j=mlat0,mlat1
      if (j >= nmlath) cycle
      if (j==1.or.j==nmlat) cycle   ! skip poles
      do i=mlon0,mlon1
        ed2(i,j)= -(phim2d(i,j+1)-phim2d(i,j-1))/(2._r8*dlatm)*dt1dts(j)/ &
          (r0*1.e-2_r8)
        elam(i,j)= -(phim2d(i,j+1)-phim2d(i,j-1))/(2._r8*dlatm)*dt0dts(j)
      enddo ! i=mlon0,mlon1
    enddo ! j=mlat0,mlat1
!
! Northern hemisphere (excluding equator):
    do j=mlat0,mlat1
      if (j <= nmlath) cycle
      if (j==1.or.j==nmlat) cycle   ! skip poles
      do i=mlon0,mlon1
        ed2(i,j) = (phim2d(i,j+1)-phim2d(i,j-1))/(2._r8*dlatm)*dt1dts(j)/ &
          (r0*1.e-2_r8)
        elam(i,j)= -(phim2d(i,j+1)-phim2d(i,j-1))/(2._r8*dlatm)*dt0dts(j)
      enddo ! i=mlon0,mlon1
    enddo ! j=mlat0,mlat1

! Need ed1,2 at global longitudes at j==2 and j==nmlat-1:
! mp_magpole_2d: Return fpole_jpm2(nglblon,1->4,nf) as:
!   1: j = jspole+1 (spole+1)
!   2: j = jspole+2 (spole+2) (unused here)
!   3: j = jnpole-1 (npole-1)
!   4: j = jnpole-2 (npole-2) (unused here)
!
    fmsub(:,:,1) = ed1(mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,2) = ed2(mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,3) = ephi(mlon0:mlon1,mlat0:mlat1)
    fmsub(:,:,4) = elam(mlon0:mlon1,mlat0:mlat1)
    call mp_magpole_2d(fmsub,mlon0,mlon1,mlat0,mlat1,nmlonp1, &
      1,nmlat,fpole2d_jpm2,4)
!
! Poles: average over four surrounding points
    do i = 1,nmlonp1
      ip1f(i) = i + nmlon/4
      if (ip1f(i) > nmlonp1) ip1f(i) = ip1f(i) - nmlon
      ip2f(i) = i + nmlon/2
      if (ip2f(i) > nmlonp1) ip2f(i) = ip2f(i) - nmlon
      ip3f(i) = i + 3*nmlon/4
      if (ip3f(i) > nmlonp1) ip3f(i) = ip3f(i) - nmlon
    enddo ! i=1,nmlonp1
!
! S pole:
    do j=mlat0,mlat1
      if (j==1) then
        do i=mlon0,mlon1
          ed1(i,j)=.25_r8*(fpole2d_jpm2(i,1,1)-fpole2d_jpm2(ip2f(i),1,1)+ &
                     fpole2d_jpm2(ip1f(i),1,2)-fpole2d_jpm2(ip3f(i),1,2))
          ed2(i,j)=.25_r8*(fpole2d_jpm2(i,1,2)-fpole2d_jpm2(ip2f(i),1,2)- &
                     fpole2d_jpm2(ip1f(i),1,1)+fpole2d_jpm2(ip3f(i),1,1))

         ephi(i,j)=.25_r8*(fpole2d_jpm2(i,1,3)-fpole2d_jpm2(ip2f(i),1,3)+ &
                     fpole2d_jpm2(ip1f(i),1,4)-fpole2d_jpm2(ip3f(i),1,4))
         elam(i,j)=.25_r8*(fpole2d_jpm2(i,1,4)-fpole2d_jpm2(ip2f(i),1,4)- &
                     fpole2d_jpm2(ip1f(i),1,3)+fpole2d_jpm2(ip3f(i),1,3))
        enddo ! i=mlon0,mlon1
! N pole:
      elseif (j==nmlat) then
        do i=mlon0,mlon1
          ed1(i,j)=.25_r8*(fpole2d_jpm2(i,3,1)-fpole2d_jpm2(ip2f(i),3,1)+ &
                     fpole2d_jpm2(ip1f(i),3,2)-fpole2d_jpm2(ip3f(i),3,2))
          ed2(i,j)=.25_r8*(fpole2d_jpm2(i,3,2)-fpole2d_jpm2(ip2f(i),3,2)- &
                     fpole2d_jpm2(ip1f(i),3,1)+fpole2d_jpm2(ip3f(i),3,1))

         ephi(i,j)=.25_r8*(fpole2d_jpm2(i,3,3)-fpole2d_jpm2(ip2f(i),3,3)+ &
                     fpole2d_jpm2(ip1f(i),3,4)-fpole2d_jpm2(ip3f(i),3,4))
         elam(i,j)=.25_r8*(fpole2d_jpm2(i,3,4)-fpole2d_jpm2(ip2f(i),3,4)- &
                     fpole2d_jpm2(ip1f(i),3,3)+fpole2d_jpm2(ip3f(i),3,3))
        enddo ! i=mlon0,mlon1
      endif ! S or N pole
    enddo ! j=mlat0,mlat1
!
! Equator: derivative of quadratic polynomial (3 given points)
! For equator and equator +/- 1 of ed2, we need equator and
! equator +/- 3 of phim2d (note feq_jpm3(nmlonp1,-3:3,1)):
!
    fmsub(:,:,1) = phim2d(mlon0:mlon1,mlat0:mlat1)
    call mp_mageq_jpm3(fmsub(:,:,1),mlon0,mlon1,mlat0,mlat1,nmlonp1,feq_jpm3,1)
    do j=mlat0,mlat1
      if (j==nmlath-1) then     ! equator-1
        do i=mlon0,mlon1
          ed2(i,j) = (4._r8*feq_jpm3(i,-2,1)-feq_jpm3(i,-3,1)- &
            3._r8*feq_jpm3(i,-1,1))/(2._r8*dlatm)/(r0*1.e-2_r8)
        enddo
      elseif (j==nmlath) then   ! equator
        do i=mlon0,mlon1
          ed2(i,j) = (4._r8*feq_jpm3(i,1,1)-feq_jpm3(i,2,1)- &
            3._r8*feq_jpm3(i,0,1))/(2._r8*dlatm)/(r0*1.e-2_r8)
         elam(i,j) = (4._r8*feq_jpm3(i,1,1)-feq_jpm3(i,2,1)- &
            3._r8*feq_jpm3(i,0,1))/(2._r8*dlatm)
        enddo
      elseif (j==nmlath+1) then ! equator+1
        do i=mlon0,mlon1
          ed2(i,j) = (4._r8*feq_jpm3(i,2,1)-feq_jpm3(i,3,1)- &
            3._r8*feq_jpm3(i,1,1))/(2._r8*dlatm)/(r0*1.e-2_r8)
        enddo
      endif ! equator +/- 1
    enddo ! j=mlat0,mlat1
!
! Set halos for 3d calculations:
    fmsub1(:,:,1) = ed1
    fmsub1(:,:,2) = ed2
    fmsub1(:,:,3) = ephi
    fmsub1(:,:,4) = elam
    call mp_mag_halos(fmsub1,mlon0,mlon1,mlat0,mlat1,4)
    ed1 = fmsub1(:,:,1)
    ed2 = fmsub1(:,:,2)
    ephi = fmsub1(:,:,3)
    elam = fmsub1(:,:,4)

    do j=mlat0,mlat1
      call outfld('ED1',ed1(mlon0:omlon1,j),omlon1-mlon0+1,j)
      call outfld('ED2',ed2(mlon0:omlon1,j),omlon1-mlon0+1,j)
    enddo
!
! Determine latitudes needed for interpolation that fall
! outside a task's latitudinal subdomain:
!
    if (debug) write(iulog,*) "pthreed: kbotdyn ", kbotdyn
   
    njneed = 0    ! number of unique latitudes needed
    jneed(:) = -1 ! j-indices of needed latitudes
    do k=kbotdyn,nmlev
      do j=mlat0,mlat1
        if (j==1.or.j==nmlat) cycle ! exclude poles
        sym = 1._r8
        if (j < nmlath) sym = -1._r8
        cosltm = cos(ylatm(j))
        do i=mlon0,mlon1
          if (i==nmlonp1) cycle

          thetam(i)=(re+zpotm3d(i,j,kbotdyn))/(re+zpotm3d(i,j,k))
          thetam(i) = acos(sqrt(thetam(i))*cosltm*(1._r8-eps))

          pslot(i) = thetam(i)*180._r8/pi+1._r8
          islot(i) = pslot(i)
          real8 = dble(islot(i))
          pslot(i) = pslot(i)-real8

          thetam(i) = ((1._r8-pslot(i))*table(islot(i),2)+pslot(i)* &
            table(islot(i)+1,2))*sym ! thetam negative for south hem

          islot(i) = i
          pslot(i) = 0._r8
          qslot(i) = (thetam(i)+pi/2._r8)/dlatm+1._r8
          jslot(i) = qslot(i)
          real8 = dble(jslot(i))
          qslot(i) = qslot(i)-real8

! Save j index if outside subdomain w/ halos:
          if ((jslot(i) < mlat00 .or. jslot(i) > mlat11).and. &
             .not.any(jslot(i)==jneed)) then
            njneed = njneed+1
            if (njneed > mxneed) call endrun('njneed')
            jneed(njneed) = jslot(i)
          endif ! jslot is outside subdomain
!
! Save j+1 index if outside subdomain:
          if ((jslot(i)+1 < mlat00 .or. jslot(i)+1 > mlat11).and. &
               .not.any(jslot(i)+1==jneed)) then
            njneed = njneed+1
            if (njneed > mxneed) call endrun('njneed')
            jneed(njneed) = jslot(i)+1
          endif ! jslot(i)+1 is outside subdomain
        enddo ! i=mlon0,mlon1
      enddo ! j=mlat0,mlat1
    enddo ! k=kbotdyn,nmlev
! 
! Get phim2 at needed latitudes (note inclusion of phim2d halos).
!     real,intent(in) :: fin(mlon00:mlon11,mlat00:mlat11,nf) ! data at current subdomain
!     real,intent(out) :: fout(mlon00:mlon11,mxneed,nf) ! returned data at needed lats
!
    fmsub1(:,:,1) = phim2d
    fmsub1(:,:,2) = ed1
    fmsub1(:,:,3) = ed2
    fmsub1(:,:,4) = ephi
    fmsub1(:,:,5) = elam
    call mp_mag_jslot(fmsub1,mlon00,mlon11,mlat00,mlat11,fmneed,jneed,mxneed,5)
    phineed = fmneed(:,:,1)
    ed1need = fmneed(:,:,2)
    ed2need = fmneed(:,:,3)
    ephineed= fmneed(:,:,4)
    elamneed= fmneed(:,:,5)

    ephi3d = 0._r8
    elam3d = 0._r8
    emz3d  = 0._r8
    do k=kbotdyn,nmlev
      do j=mlat0,mlat1
        if (j==1.or.j==nmlat) cycle ! exclude poles
        sym = 1._r8
        if (j < nmlath) sym = -1._r8
        cosltm = cos(ylatm(j))
        do i=mlon0,mlon1
          if (i==nmlonp1) cycle

          thetam(i)=(re+zpotm3d(i,j,kbotdyn))/(re+zpotm3d(i,j,k))
          thetam(i) = acos(sqrt(thetam(i))*cosltm*(1._r8-eps))
          fac_elam = tan(ylatm(j))/tan(thetam(i)*sym)  ! tan(lambda_q)/tan(lambda_m)

          pslot(i) = thetam(i)*180._r8/pi+1._r8
          islot(i) = pslot(i)
          real8 = dble(islot(i))
          pslot(i) = pslot(i)-real8

          thetam(i) = ((1._r8-pslot(i))*table(islot(i),2)+pslot(i)* &
            table(islot(i)+1,2))*sym ! thetam negative for south hem

          islot(i) = i
          pslot(i) = 0._r8
          qslot(i) = (thetam(i)+pi/2._r8)/dlatm+1._r8
          jslot(i) = qslot(i)
          real8 = dble(jslot(i))
          qslot(i) = qslot(i)-real8
!
! Check for jslot in subdomain:
          if (jslot(i) >= mlat00.and.jslot(i) <= mlat11) then ! within subdomain
            phi0j0 = phim2d(islot(i)  ,jslot(i))
            phi1j0 = phim2d(islot(i)+1,jslot(i))
            ed1i0j0 =   ed1(islot(i)  ,jslot(i))
            ed1i1j0 =   ed1(islot(i)+1,jslot(i))
            ed2i0j0 =   ed2(islot(i)  ,jslot(i))
            ed2i1j0 =   ed2(islot(i)+1,jslot(i))
            ephi0j0 =  ephi(islot(i)  ,jslot(i))
            ephi1j0 =  ephi(islot(i)+1,jslot(i))
            elam0j0 =  elam(islot(i)  ,jslot(i))
            elam1j0 =  elam(islot(i)+1,jslot(i))
          else                                           ! jslot outside subdomain
            n = ixfind(jneed,mxneed,jslot(i),icount)
            if (n==0) then
              write(iulog,"('>>> pthreed: i=',i4,' j=',i4,' k=',i4)") i,j,k
              write(iulog,"('  Could not find jslot ',i4,' in jneed=',/,(10i4))") &
                i,j,k,jslot(i),jneed
              call endrun('jslot(i) not in jneed')
            endif
            phi0j0  = phineed(islot(i)  ,n)
            phi1j0  = phineed(islot(i)+1,n)
            ed1i0j0 = ed1need(islot(i)  ,n)
            ed1i1j0 = ed1need(islot(i)+1,n)
            ed2i0j0 = ed2need(islot(i)  ,n)
            ed2i1j0 = ed2need(islot(i)+1,n)
            ephi0j0 =ephineed(islot(i)  ,n)
            ephi1j0 =ephineed(islot(i)+1,n)
            elam0j0 =elamneed(islot(i)  ,n)
            elam1j0 =elamneed(islot(i)+1,n)
          endif
!
! Check for jslot+1 in subdomain:
          if (jslot(i)+1 >= mlat00.and.jslot(i)+1 <= mlat11) then ! within subdomain
            phi0j1 = phim2d(islot(i)  ,jslot(i)+1)
            phi1j1 = phim2d(islot(i)+1,jslot(i)+1)
            ed1i0j1 =   ed1(islot(i)  ,jslot(i)+1)
            ed1i1j1 =   ed1(islot(i)+1,jslot(i)+1)
            ed2i0j1 =   ed2(islot(i)  ,jslot(i)+1)
            ed2i1j1 =   ed2(islot(i)+1,jslot(i)+1)
            ephi0j1 =  ephi(islot(i)  ,jslot(i)+1)
            ephi1j1 =  ephi(islot(i)+1,jslot(i)+1)
            elam0j1 =  elam(islot(i)  ,jslot(i)+1)
            elam1j1 =  elam(islot(i)+1,jslot(i)+1)
          else                                           ! jslot+1 outside subdomain
            n = ixfind(jneed,mxneed,jslot(i)+1,icount)
            if (n==0) then
              write(iulog,"('>>> pthreed: i=',i4,' j=',i4,' k=',i4)") i,j,k
              write(iulog,"('  Could not find jslot+1 ',i4,' in jneed=',/,(10i4))") &
                i,j,k,jslot(i)+1,jneed
              call endrun('jslot(i)+1 not in jneed')
            endif
            phi0j1  = phineed(islot(i)  ,n)
            phi1j1  = phineed(islot(i)+1,n)
            ed1i0j1 = ed1need(islot(i)  ,n)
            ed1i1j1 = ed1need(islot(i)+1,n)
            ed2i0j1 = ed2need(islot(i)  ,n)
            ed2i1j1 = ed2need(islot(i)+1,n)
            ephi0j1 =ephineed(islot(i)  ,n)
            ephi1j1 =ephineed(islot(i)+1,n)
            elam0j1 =elamneed(islot(i)  ,n)
            elam1j1 =elamneed(islot(i)+1,n)
          endif
!
! Do the interpolation:
          phim3d(i,j,k) = (1._r8-qslot(i))*((1._r8-pslot(i))*           &
            phi0j0 + pslot(i) * phi1j0) + qslot(i)*((1._r8-pslot(i))*   &
            phi0j1 + pslot(i) * phi1j1)

          ed13d(i,j,k) = (1._r8-qslot(i))*((1._r8-pslot(i))*            &
            ed1i0j0 + pslot(i) * ed1i1j0) + qslot(i)*((1._r8-pslot(i))* &
            ed1i0j1 + pslot(i) * ed1i1j1)
          ed23d(i,j,k) = (1._r8-qslot(i))*((1._r8-pslot(i))*            &
            ed2i0j0 + pslot(i) * ed2i1j0) + qslot(i)*((1._r8-pslot(i))* &
            ed2i0j1 + pslot(i) * ed2i1j1)

          ephi3d(i,j,k) = (1._r8-qslot(i))*((1._r8-pslot(i))*           &
            ephi0j0 + pslot(i) * ephi1j0) + qslot(i)*((1._r8-pslot(i))* &
            ephi0j1 + pslot(i) * ephi1j1)
          elam3d(i,j,k) = (1._r8-qslot(i))*((1._r8-pslot(i))*           &
            elam0j0 + pslot(i) * elam1j0) + qslot(i)*((1._r8-pslot(i))* &
            elam0j1 + pslot(i) * elam1j1)
          elam3d(i,j,k) = elam3d(i,j,k)*fac_elam  ! add height variation
        enddo ! i=mlon0,mlon1
      enddo ! j=mlat0,mlat1
    enddo ! k=kbotdyn,nmlev

!
! Mag poles for phim: 
! mp_magpoles returns global longitudes at S,N poles in fpoles(nglblon,2,nf)
!
    call mp_magpoles(phim2d(mlon0:mlon1,mlat0:mlat1),   &
      mlon0,mlon1,mlat0,mlat1,nmlonp1,1,nmlat,fpoles,1)

    real8 = dble(nmlon)
    phims=dot_product(unitvm,fpoles(1:nmlon,1,1))/real8
    phimn=dot_product(unitvm,fpoles(1:nmlon,2,1))/real8

    do k=kbotdyn,nmlev
      do j=mlat0,mlat1
        if (j==1) then
          do i=mlon0,mlon1
            phim3d(i,j,k) = phims
            ed13d(i,j,k) = ed1(i,j)
            ed23d(i,j,k) = ed2(i,j)
            ephi3d(i,j,k) = ephi(i,j)
            elam3d(i,j,k) = ed2(i,j)*(r0*1.e-2_r8)
          enddo
        elseif (j==nmlat) then
          do i=mlon0,mlon1
            phim3d(i,j,k) = phimn
            ed13d(i,j,k) = ed1(i,j) 
            ed23d(i,j,k) = ed2(i,j)
            ephi3d(i,j,k) = ephi(i,j)
            elam3d(i,j,k) = -ed2(i,j)*(r0*1.e-2_r8)
          enddo
        endif ! poles
      enddo ! j=mlat0,mlat1
    enddo ! k=kbotdyn,nmlev
!
! Extend kbotdyn downward redundantly:
    do k=1,kbotdyn-1
      phim3d(:,:,k) = phim3d(:,:,kbotdyn)
      ephi3d(:,:,k) = ephi3d(:,:,kbotdyn)
      elam3d(:,:,k) = elam3d(:,:,kbotdyn)
    enddo
!
! Upward electric field:
    do k=kbotdyn,nmlev-1
      do j=mlat0,mlat1
        do i=mlon0,mlon1
          emz3d(i,j,k) = -(phim3d(i,j,k+1)-phim3d(i,j,k-1))
        enddo
      enddo
    enddo
! 
    do k=mlev0,mlev1
      call mp_mag_periodic_f2d(phim3d(:,:,k),mlon0,mlon1,mlat0,mlat1,1)
    enddo
!
    do j=mlat0,mlat1
      call outfld('EPHI3D',ephi3d(mlon0:omlon1,j,mlev1:mlev0:-1),omlon1-mlon0+1,j)
      call outfld('ELAM3D',elam3d(mlon0:omlon1,j,mlev1:mlev0:-1),omlon1-mlon0+1,j)
      call outfld('EMZ3D', emz3d(mlon0:omlon1,j,mlev1:mlev0:-1),omlon1-mlon0+1,j)

      call outfld('PHIM3D',phim3d(mlon0:omlon1,j,mlev1:mlev0:-1),omlon1-mlon0+1,j)
      call outfld('ED13D' ,ed13d (mlon0:omlon1,j,mlev1:mlev0:-1),omlon1-mlon0+1,j)
      call outfld('ED23D' ,ed23d (mlon0:omlon1,j,mlev1:mlev0:-1),omlon1-mlon0+1,j)
    enddo
  end subroutine pthreed
!-----------------------------------------------------------------------
  subroutine pefield
  use edyn_params ,only: pi
  use edyn_maggrid,only: dt0dts,dlatm,dlonm,rcos0s
  use edyn_geogrid,only: nlev
  use edyn_mpi    ,only: mp_magpole_3d,mp_mag_halos,mp_magpoles
  use edyn_esmf   ,only: mag_ephi3d,mag_elam3d,mag_emz3d,mag_phi3d,&
                         geo_ephi3d,geo_elam3d,geo_emz3d,geo_phi3d
!
! Local:
  integer :: i,ii,j,k
  real(r8) :: &
    phi3d(mlon0-1:mlon1+1,mlat0-1:mlat1+1,nmlev), & ! local phi w/ halos
    fpole3d_jpm2(nmlonp1,4,nmlev,1) ! global lons at S pole+1,2 and N pole-1,2
  real(r8) :: csth0,real8
  real(r8) :: fpoles(nmlonp1,2,nmlev) ! global lons at poles
  real(r8),dimension(lon0:lon1,lat0:lat1,nlev) :: exgeo,eygeo,ezgeo
!
! Copy phim3d to local phi3d, and set halo points:
    do j=mlat0,mlat1
      do i=mlon0,mlon1
        phi3d(i,j,:) = phim3d(i,j,:)
      enddo
    enddo
    call mp_mag_halos(phi3d,mlon0,mlon1,mlat0,mlat1,nmlev)
!
! Return fpole3d_jpm2(nglblon,1->4,nlev,nf) as:
!   1: j = jspole+1 (spole+1)
!   2: j = jspole+2 (spole+2) not used here
!   3: j = jnpole-1 (npole-1)
!   4: j = jnpole-2 (npole-2) not used here
!
    call mp_magpole_3d(phim3d(mlon0:mlon1,mlat0:mlat1,:),mlon0,&
      mlon1,mlat0,mlat1,nmlev,nmlonp1,1,nmlat,fpole3d_jpm2,1)
!
! Set j=0 and j=nmlat+1 of local phi3d. This overwrites the far
! north and south halo points set by mp_mag_halos above.
    do j=mlat0,mlat1
      if (j==1) then
        do i=mlon0,mlon1
          ii = 1+mod(i-1+nmlon/2,nmlon) ! over the south pole
          phi3d(i,j-1,:) = fpole3d_jpm2(ii,1,:,1)
        enddo
      elseif (j==nmlat) then
        do i=mlon0,mlon1
          ii = 1+mod(i-1+nmlon/2,nmlon) ! over the north pole
          phi3d(i,j+1,:) = fpole3d_jpm2(ii,3,:,1)
        enddo
      endif ! poles or not
    enddo ! j=mlat0,mlat1
!
! Meridional component of electric field:
    do j=mlat0,mlat1
      do i=mlon0,mlon1
        elam3d(i,j,:) = -(phi3d(i,j+1,:)-phi3d(i,j-1,:))/ &
          (2._r8*dlatm)*dt0dts(j)
      enddo
    enddo
!
! Zonal component of electric field:
    do j=mlat0,mlat1
      if (j==1.or.j==nmlat) cycle
      real8 = dble(j-1)
      csth0 = cos(-pi/2._r8+real8*dlatm)
      do i=mlon0,mlon1
        ephi3d(i,j,:) = -(phi3d(i+1,j,:)-phi3d(i-1,j,:))/ &
          (2._r8*dlonm*csth0)*rcos0s(j)
      enddo
    enddo
! 
! Polar values for ephi3d (need global lons at poles of elam3d):
    call mp_magpoles(elam3d,mlon0,mlon1,mlat0,mlat1,nmlonp1,1,nmlat,fpoles,nmlev)
    do j=mlat0,mlat1
      if (j==1) then         ! south pole
        do i=mlon0,mlon1
          ii = 1+mod(i-1+(nmlon/4),nmlon) ! over the south pole
          ephi3d(i,j,:) = fpoles(ii,1,:)
        enddo
      elseif (j==nmlat) then ! north pole
        do i=mlon0,mlon1
          ii = 1+mod(i-1+((3*nmlon)/4),nmlon) ! over the north pole
          ephi3d(i,j,:) = fpoles(ii,2,:)
        enddo
      endif ! poles or not
    enddo ! j=mlat0,mlat1
!       
! emz = d(phi)/dz
    do k=2,nmlev-1
      do j=mlat0,mlat1
        do i=mlon0,mlon1
          emz3d(i,j,k) = -(phim3d(i,j,k+1)-phi3d(i,j,k-1))
        enddo
      enddo
    enddo ! k=2,nmlev-1
!
! btf 6/18/14: mag2geo is not working due to error return rc=51 from 
!   ESMF_FieldSMM for 3d mag2geo (see sub esmf_regrid in edyn_esmf.F90)
!   (this is the call to do the data regridding, not the init call)
!
! Use ESMF to regrid the electric field to the geographic grid:
    call mag2geo_3d(ephi3d,exgeo ,mag_ephi3d,geo_ephi3d,'EPHI3D  ')
    call mag2geo_3d(elam3d,eygeo ,mag_elam3d,geo_elam3d,'ELAM3D  ')
    call mag2geo_3d(emz3d ,ezgeo ,mag_emz3d ,geo_emz3d ,'EMZ3D   ')
    call mag2geo_3d(phim3d,phig3d,mag_phi3d ,geo_phi3d ,'PHIM3D  ')
!
! Define ex,ey,ez on geographic subdomains for ionvel:
    do j=lat0,lat1
      do i=lon0,lon1
        ex(:,i,j) = exgeo(i,j,:)
        ey(:,i,j) = eygeo(i,j,:)
        ez(:,i,j) = ezgeo(i,j,:)
        poten(:,i,j) = phig3d(i,j,:)
      enddo
    enddo

! ex,ey,ez(nlev,lon0-2,lon1+2,lat0:lat1)
    if (debug) then
      write(iulog,"('pefield after mag2geo: ex=',2e12.4,' ey=',2e12.4,' ez=',2e12.4)") &
   	minval(ex(:,lon0:lon1,:)),maxval(ex(:,lon0:lon1,:)), &
   	minval(ey(:,lon0:lon1,:)),maxval(ey(:,lon0:lon1,:)), &
   	minval(ez(:,lon0:lon1,:)),maxval(ez(:,lon0:lon1,:))
    endif

    call savefld_waccm_switch(poten(1:nlev,lon0:lon1,lat0:lat1),'POTEN',&
         nlev,lon0,lon1,lat0,lat1)
    call savefld_waccm_switch(ex(1:nlev,lon0:lon1,lat0:lat1),'EX',&
         nlev,lon0,lon1,lat0,lat1)
    call savefld_waccm_switch(ey(1:nlev,lon0:lon1,lat0:lat1),'EY',&
         nlev,lon0,lon1,lat0,lat1)
    call savefld_waccm_switch(ez(1:nlev,lon0:lon1,lat0:lat1),'EZ',&
         nlev,lon0,lon1,lat0,lat1)

  end subroutine pefield
!-----------------------------------------------------------------------
  subroutine ionvel(z,ui,vi,wi)
!
! Calculate 3d ExB ion drifts from electric field (sub pefield)
! on geographic grid.
!
    use edyn_params  ,only: re
    use edyn_geogrid ,only: nlev
    use getapex      ,only: &
      rjac  ,& ! (nlon+1,jspole:jnpole,2,2)
      bmod  ,& ! magnitude of magnetic field (nlon+1,jspole:jnpole)
      xb,yb,zb ! north,east,down magnetic field (nlon+1,jspole:jnpole)
!
! Args:
    real(r8),intent(in),dimension(lev0:lev1,lon0:lon1,lat0:lat1) :: &
      z  ! geopotential from input (cm)
    real(r8),intent(out),dimension(lev0:lev1,lon0:lon1,lat0:lat1) :: &
      ui,vi,wi
!
! Local:
    integer :: i,ii,k,j
    real(r8),dimension(lev0:lev1,lon0:lon1) :: eex,eey,eez
    real(r8),dimension(lev0:lev1,lon0:lon1,lat0:lat1) :: rjac_out

! mag field diagnostics
    call savefld_waccm_switch(bmod(lon0:lon1,lat0:lat1),'BMOD',1,lon0,lon1,lat0,lat1)
    call savefld_waccm_switch(xb(lon0:lon1,lat0:lat1),'XB',1,lon0,lon1,lat0,lat1)
    call savefld_waccm_switch(yb(lon0:lon1,lat0:lat1),'YB',1,lon0,lon1,lat0,lat1)
    call savefld_waccm_switch(zb(lon0:lon1,lat0:lat1),'ZB',1,lon0,lon1,lat0,lat1)

!
! Scan geographic latitude subdomain:
!
    do j=lat0,lat1
      do i=lon0,lon1
        ii = i
        do k=lev0,lev1
          eex(k,i) = (rjac(ii,j,1,1)*ex(k,i,j)+ &
                      rjac(ii,j,2,1)*ey(k,i,j))/(re+z(k,i,j))
          eey(k,i) = (rjac(ii,j,1,2)*ex(k,i,j)+ &
                      rjac(ii,j,2,2)*ey(k,i,j))/(re+z(k,i,j))
        enddo ! k=lev0,lev1
      enddo
!
      do i=lon0,lon1
        do k=lev0+1,lev1-1
          eez(k,i) = ez(k,i,j)/(z(k+1,i,j)-z(k-1,i,j))
        enddo ! k=lev0+1,lev1-1
      enddo 
!
! Extrapolate for lower and upper boundaries:
      do i=lon0,lon1
        eez(lev0,i) = 2._r8*eez(2,i)-eez(3,i)
        eez(lev1,i) = 2._r8*eez(lev1-1,i)-eez(lev1-2,i)
      enddo

    if (debug.and.masterproc) then
      write(iulog,"('ionvel: j=',i4,' eex=',2e12.4,' eey=',2e12.4,' eez=',2e12.4)") &
        j,minval(eex),maxval(eex),minval(eey),maxval(eey),minval(eez),maxval(eez)
    endif

! 
! ion velocities = (e x b/b**2) (x 1.e6 for m/sec)
! ui = zonal, vi = meridional, wi = vertical
!
      do k=lev0,lev1
        do i=lon0,lon1
          ii = i
          ui(k,i,j) = -(eey(k,i)*zb(ii,j)+eez(k,i)*xb(ii,j))* &
            1.e6_r8/bmod(ii,j)**2
          vi(k,i,j) =  (eez(k,i)*yb(ii,j)+eex(k,i)*zb(ii,j))* &
            1.e6_r8/bmod(ii,j)**2
          wi(k,i,j) =  (eex(k,i)*xb(ii,j)-eey(k,i)*yb(ii,j))* &
            1.e6_r8/bmod(ii,j)**2
        enddo ! i=lon0,lon1
      enddo ! k=lev0,lev1

    if (debug.and.masterproc) then
      write(iulog,"('ionvel: j=',i4,' ui=',2e12.4,' vi=',2e12.4,' wi=',2e12.4)") &
        j,minval(ui),maxval(ui),minval(vi),maxval(vi),minval(wi),maxval(wi)
    endif
!
! Output ion drifts in cm/s for oplus_xport call from dpie_coupling:
      do i=lon0,lon1
        ui(:,i,j) = ui(:,i,j)*100._r8
        vi(:,i,j) = vi(:,i,j)*100._r8
        wi(:,i,j) = wi(:,i,j)*100._r8
      enddo
    enddo ! j=lat0,lat1 

    if (debug.and.masterproc) then
      write(iulog,"('ionvel: ion drifts on geo grid: ui=',2e12.4,' vi=',2e12.4,' wi=',2e12.4)") &
                  minval(ui),maxval(ui), minval(vi),maxval(vi), minval(wi),maxval(wi)
    endif

    if (hist_fld_active('RJAC11')) then
       do i=1,nlev
          rjac_out(i,lon0:lon1,lat0:lat1) = rjac(lon0:lon1,lat0:lat1,1,1)
       end do
       call savefld_waccm_switch(rjac_out,'RJAC11',nlev,lon0,lon1,lat0,lat1)
    endif

    if (hist_fld_active('RJAC12')) then
       do i=1,nlev
          rjac_out(i,lon0:lon1,lat0:lat1) = rjac(lon0:lon1,lat0:lat1,1,2)
       end do
       call savefld_waccm_switch(rjac_out,'RJAC12',nlev,lon0,lon1,lat0,lat1)
    endif

    if (hist_fld_active('RJAC21')) then
       do i=1,nlev
          rjac_out(i,lon0:lon1,lat0:lat1) = rjac(lon0:lon1,lat0:lat1,2,1)
       end do
       call savefld_waccm_switch(rjac_out,'RJAC21',nlev,lon0,lon1,lat0,lat1)
    endif

    if (hist_fld_active('RJAC22')) then
       do i=1,nlev
          rjac_out(i,lon0:lon1,lat0:lat1) = rjac(lon0:lon1,lat0:lat1,2,2)
       end do
       call savefld_waccm_switch(rjac_out,'RJAC22',nlev,lon0,lon1,lat0,lat1)
    endif

  end subroutine ionvel
!-----------------------------------------------------------------------
  subroutine mag2geo_3d(fmag,fgeo,ESMF_mag,ESMF_geo,fname)
!
! Convert field on geomagnetic grid fmag to geographic grid in fgeo.
!
  use edyn_esmf,only: edyn_esmf_set3d_mag,edyn_esmf_regrid,edyn_esmf_get_3dfield
  use edyn_geogrid,only: nlev
!
! Args:
! integer,intent(in) :: mlon0,mlon1,mlat0,mlat1,nmlev,lon0,lon1,lat0,lat1,nlev
  character(len=*) :: fname
  type(ESMF_Field),intent(inout) :: ESMF_mag, ESMF_geo
  real(r8),intent(in) :: fmag(mlon0:mlon1,mlat0:mlat1,nmlev)
  real(r8),intent(out) :: fgeo(lon0:lon1,lat0:lat1,nlev)
!
! Local:
  integer :: j
  character(len=8) :: fnames(1)
  type(ESMF_Field) :: magfields(1)
  real(r8),pointer,dimension(:,:,:) :: fptr

  fgeo = finit
  fnames(1) = fname
  magfields(1) = ESMF_mag
!
! Put fmag into ESMF mag field on mag source grid:
  call edyn_esmf_set3d_mag(magfields,fnames,fmag,1,1,nmlev,mlon0,mlon1,mlat0,mlat1)
!
! Regrid to geographic destination grid, defining ESMF_geo:
  call edyn_esmf_regrid(ESMF_mag,ESMF_geo,'mag2geo',3)
!
! Put regridded geo field into pointer:
  call edyn_esmf_get_3dfield(ESMF_geo,fptr,fname)
!
! Transfer from pointer to output arg:
  do j=lat0,lat1
    fgeo(:,j,:) = fptr(:,j,:)
  enddo
  end subroutine mag2geo_3d
#endif

!-----------------------------------------------------------------------
end module edynamo
