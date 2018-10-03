module oplus
!
! Horizontally transport the O+ ion, adapted for WACCM-X from TIEGCM. 
! Input O+ is received from WACCM physics/chemistry, transported O+ 
! (op_out and opnm_out) are passed back to chemistry.
!
! B. Foster (foster@ucar.edu), May, 2015.
!
  use shr_kind_mod   ,only: r8 => shr_kind_r8
  use cam_abortutils ,only: endrun
  use cam_logfile    ,only: iulog
  use savefield_waccm,only: savefld_waccm, savefld_waccm_switch     ! save field to waccm history
  use edyn_geogrid   ,only: dphi,dlamda,cs,zp,expz,p0 !, nlon, nlat, nlev
  use getapex        ,only: bx,by,bz,bmod2       ! (0:nlonp1,jspole-1:jnpole+1)
  use edyn_params    ,only: re
  use time_manager   ,only: get_step_size,is_first_step,is_first_restart_step
  use edyn_mpi       ,only: array_ptr_type
  use shr_const_mod  ,only: shr_const_g         ! gravitational constant (m/s^2)
  use spmd_utils     ,only: masterproc

  implicit none
  private
  public :: oplus_xport,  oplus_init
  public :: kbot

  real(r8) :: pi,rtd
!
! Constants in CGS:
!
  real(r8),parameter :: boltz = 1.38E-16_r8 ! boltzman's constant (erg/kelvin)
  real(r8),parameter :: gask  = 8.314e7_r8  ! gas constant (erg/mol)
!
! Collision factor (tuneable) (see also local colfac in iondrag.F90)
! FIX: Collision factor colfac is set locally in iondrag.F90 and here.
!      It should be in one location and shared between ionosphere and
!      dpie_coupling.
!
    real(r8),parameter :: colfac = 1.5_r8   ! see also iondrag.F90
!
! Reciprocal of molecular mass (multiply is cheaper than divide)
    real(r8),parameter :: rmassinv_o2=1._r8/32._r8, rmassinv_o1=1._r8/16._r8, &
                          rmassinv_n2=1._r8/28._r8
    real(r8),parameter :: rmass_op=16._r8

    real(r8) :: dzp                       ! delta zp (typically 0.5 from kbot to top)
    real(r8) :: grav_cm                   ! gravitational constant (cm/s^2)
    integer, protected :: kbot = -999     ! k-index corresponding to ~pbot
    real(r8),parameter :: pbot = 0.004_r8 ! Pa -- bottom of O+ transport (near 120 km)
!
! The shapiro constant .03 is used for spatial smoothing of oplus,
! (shapiro is tuneable, and maybe should be a function of timestep size).
! dtsmooth and dtsmooth_div2 are used in the time smoothing.
! To turn off all smoothing here, set shapiro=0. and dtsmooth = 1. 
!

    real(r8),parameter ::    & 
      dtsmooth = 0.95_r8,    &                ! for time smoother
      dtsmooth_div2 = 0.5_r8*(1._r8-dtsmooth)
 
    real(r8) :: adiff_limiter
    real(r8) :: shapiro_const
    logical  :: enforce_floor
    logical, parameter :: debug = .false.

  contains

!-----------------------------------------------------------------------
  subroutine oplus_init( adiff_limiter_in, shapiro_const_in, enforce_floor_in )

    use cam_history,  only : addfld, horiz_only
    use filter_module,only : filter_init

    real(r8), intent(in) :: adiff_limiter_in
    real(r8), intent(in) :: shapiro_const_in
    logical , intent(in) :: enforce_floor_in

    shapiro_const = shapiro_const_in
    enforce_floor = enforce_floor_in
    adiff_limiter = adiff_limiter_in
    
    call filter_init

    !
    ! Save fields from oplus module:
    !
    call addfld ('OPLUS_Z'    ,(/ 'lev' /), 'I', 'cm   ','OPLUS_Z'     , gridname='fv_centers')
    call addfld ('OPLUS_TN'   ,(/ 'lev' /), 'I', 'deg K','OPLUS_TN'    , gridname='fv_centers')
    call addfld ('OPLUS_TE'   ,(/ 'lev' /), 'I', 'deg K','OPLUS_TE'    , gridname='fv_centers')
    call addfld ('OPLUS_TI'   ,(/ 'lev' /), 'I', 'deg K','OPLUS_TI'    , gridname='fv_centers')
    call addfld ('OPLUS_UN'   ,(/ 'lev' /), 'I', 'cm/s' ,'OPLUS_UN'    , gridname='fv_centers')
    call addfld ('OPLUS_VN'   ,(/ 'lev' /), 'I', 'cm/s' ,'OPLUS_VN'    , gridname='fv_centers')
    call addfld ('OPLUS_OM'   ,(/ 'lev' /), 'I', 'Pa/s' ,'OPLUS_OM'    , gridname='fv_centers')
    call addfld ('OPLUS_O2'   ,(/ 'lev' /), 'I', 'mmr'  ,'OPLUS_O2'    , gridname='fv_centers')
    call addfld ('OPLUS_O1'   ,(/ 'lev' /), 'I', 'mmr'  ,'OPLUS_O1'    , gridname='fv_centers') 

    call addfld ('OPLUS_N2'   ,(/ 'lev' /), 'I', 'mmr'  ,'OPLUS_N2'    , gridname='fv_centers')
    call addfld ('OPLUS_OP'   ,(/ 'lev' /), 'I', 'cm^3' ,'OPLUS_OP'    , gridname='fv_centers')
    call addfld ('OPLUS_UI'   ,(/ 'lev' /), 'I', 'm/s'  ,'OPLUS_UI'    , gridname='fv_centers')
    call addfld ('OPLUS_VI'   ,(/ 'lev' /), 'I', 'm/s'  ,'OPLUS_VI'    , gridname='fv_centers')
    call addfld ('OPLUS_WI'   ,(/ 'lev' /), 'I', 'm/s'  ,'OPLUS_WI'    , gridname='fv_centers')
    call addfld ('OPLUS_MBAR' ,(/ 'lev' /), 'I', ' '    ,'OPLUS_MBAR'  , gridname='fv_centers')
    call addfld ('OPLUS_TR'   ,(/ 'lev' /), 'I', ' '    ,'OPLUS_TR'    , gridname='fv_centers')
    call addfld ('OPLUS_TP0'  ,(/ 'lev' /), 'I', ' '    ,'OPLUS_TP0'   , gridname='fv_centers')
    call addfld ('OPLUS_TP1'  ,(/ 'lev' /), 'I', ' '    ,'OPLUS_TP1'   , gridname='fv_centers')
    !       call addfld ('OPLUS_TP2'  ,(/ 'lev' /), 'I', ' '    ,'OPLUS_TP2'   , gridname='fv_centers')
    call addfld ('OPLUS_DJ'   ,(/ 'lev' /), 'I', ' '    ,'OPLUS_DJ'    , gridname='fv_centers')
    call addfld ('OPLUS_HJ'   ,(/ 'lev' /), 'I', ' '    ,'OPLUS_HJ'    , gridname='fv_centers')
    call addfld ('OPLUS_BVEL' ,(/ 'lev' /), 'I', ' '    ,'OPLUS_BVEL'  , gridname='fv_centers')
    call addfld ('OPLUS_DIFFJ',(/ 'lev' /), 'I', ' '    ,'OPLUS_DIFFJ' , gridname='fv_centers')
    call addfld ('OPLUS_OPNM' ,(/ 'lev' /), 'I', ' '    ,'OPLUS_OPNM'  , gridname='fv_centers')
    call addfld ('OPNM_SMOOTH',(/ 'lev' /), 'I', ' '    ,'OPNM_SMOOTH' , gridname='fv_centers')
    call addfld ('BDOTDH_OP'  ,(/ 'lev' /), 'I', ' '    ,'BDOTDH_OP'   , gridname='fv_centers')
    call addfld ('BDOTDH_OPJ' ,(/ 'lev' /), 'I', ' '    ,'BDOTDH_OPJ'  , gridname='fv_centers')
    call addfld ('BDOTDH_DIFF',(/ 'lev' /), 'I', ' '    ,'BDOTDH_DIFF' , gridname='fv_centers')
    call addfld ('BDZDVB_OP'  ,(/ 'lev' /), 'I', ' '    ,'BDZDVB_OP'   , gridname='fv_centers') 
    call addfld ('EXPLICIT0'  ,(/ 'lev' /), 'I', ' '    ,'EXPLICIT0'   , gridname='fv_centers')

    call addfld ('EXPLICITa'  ,(/ 'lev' /), 'I', ' '    ,'EXPLICITa'   , gridname='fv_centers') ! part a
    call addfld ('EXPLICITb'  ,(/ 'lev' /), 'I', ' '    ,'EXPLICITb'   , gridname='fv_centers') ! part b
    call addfld ('EXPLICIT1'  ,(/ 'lev' /), 'I', ' '    ,'EXPLICIT1'   , gridname='fv_centers') ! complete
    call addfld ('EXPLICIT'   ,(/ 'lev' /), 'I', ' '    ,'EXPLICIT'    , gridname='fv_centers') ! final w/ poles

    call addfld ('EXPLICIT2'  ,(/ 'lev' /), 'I', ' '    ,'EXPLICIT2'   , gridname='fv_centers')
    call addfld ('EXPLICIT3'  ,(/ 'lev' /), 'I', ' '    ,'EXPLICIT3'   , gridname='fv_centers')
    call addfld ('TPHDZ0'     ,(/ 'lev' /), 'I', ' '    ,'TPHDZ0'      , gridname='fv_centers')
    call addfld ('TPHDZ1'     ,(/ 'lev' /), 'I', ' '    ,'TPHDZ1'      , gridname='fv_centers')
    call addfld ('DIVBZ'      ,(/ 'lev' /), 'I', ' '    ,'DIVBZ'       , gridname='fv_centers')
    call addfld ('HDZMBZ'     ,(/ 'lev' /), 'I', ' '    ,'HDZMBZ'      , gridname='fv_centers')
    call addfld ('HDZPBZ'     ,(/ 'lev' /), 'I', ' '    ,'HDZPBZ'      , gridname='fv_centers')
    call addfld ('P_COEFF0'   ,(/ 'lev' /), 'I', ' '    ,'P_COEFF0'    , gridname='fv_centers')
    call addfld ('Q_COEFF0'   ,(/ 'lev' /), 'I', ' '    ,'Q_COEFF0'    , gridname='fv_centers')
    call addfld ('R_COEFF0'   ,(/ 'lev' /), 'I', ' '    ,'R_COEFF0'    , gridname='fv_centers')
    call addfld ('P_COEFF0a'  ,(/ 'lev' /), 'I', ' '    ,'P_COEFF0a'   , gridname='fv_centers')
    call addfld ('Q_COEFF0a'  ,(/ 'lev' /), 'I', ' '    ,'Q_COEFF0a'   , gridname='fv_centers')
    call addfld ('DJINT'      ,(/ 'lev' /), 'I', ' '    ,'DJINT'       , gridname='fv_centers')
    call addfld ('BDOTU'      ,(/ 'lev' /), 'I', ' '    ,'BDOTU'       , gridname='fv_centers')
    call addfld ('R_COEFF0a'  ,(/ 'lev' /), 'I', ' '    ,'R_COEFF0a'   , gridname='fv_centers')
    call addfld ('P_COEFF1'   ,(/ 'lev' /), 'I', ' '    ,'P_COEFF1'    , gridname='fv_centers')
    call addfld ('Q_COEFF1'   ,(/ 'lev' /), 'I', ' '    ,'Q_COEFF1'    , gridname='fv_centers')
    call addfld ('R_COEFF1'   ,(/ 'lev' /), 'I', ' '    ,'R_COEFF1'    , gridname='fv_centers')
    call addfld ('P_COEFF2'   ,(/ 'lev' /), 'I', ' '    ,'P_COEFF2'    , gridname='fv_centers')
    call addfld ('Q_COEFF2'   ,(/ 'lev' /), 'I', ' '    ,'Q_COEFF2'    , gridname='fv_centers')
    call addfld ('R_COEFF2'   ,(/ 'lev' /), 'I', ' '    ,'R_COEFF2'    , gridname='fv_centers')

    call addfld ('P_COEFF'    ,(/ 'lev' /), 'I', ' '    ,'P_COEFF'     , gridname='fv_centers') ! final w/ poles
    call addfld ('Q_COEFF'    ,(/ 'lev' /), 'I', ' '    ,'Q_COEFF'     , gridname='fv_centers') ! final w/ poles
    call addfld ('R_COEFF'    ,(/ 'lev' /), 'I', ' '    ,'R_COEFF'     , gridname='fv_centers') ! final w/ poles

    call addfld ('OP_SOLVE'   ,(/ 'lev' /), 'I', ' '    ,'OP_SOLVE'    , gridname='fv_centers')

    call addfld ('OP_OUT'     ,(/ 'lev' /), 'I', 'cm^3' ,'OPLUS (oplus_xport output)', gridname='fv_centers')
    call addfld ('OPNM_OUT'   ,(/ 'lev' /), 'I', 'cm^3' ,'OPNM_OUT'    , gridname='fv_centers')
    call addfld ('BMOD2'      ,(/ 'lev' /), 'I', ' '    ,'BMOD2'       , gridname='fv_centers')

    call addfld ('OPLUS_FLUX', horiz_only , 'I', ' ','OPLUS_FLUX', gridname='fv_centers')
    call addfld ('OPLUS_DIVB', horiz_only , 'I', ' ','OPLUS_DIVB', gridname='fv_centers')
    call addfld ('OPLUS_BX'  , horiz_only , 'I', ' ','OPLUS_BX'  , gridname='fv_centers')
    call addfld ('OPLUS_BY'  , horiz_only , 'I', ' ','OPLUS_BY'  , gridname='fv_centers')
    call addfld ('OPLUS_BZ'  , horiz_only , 'I', ' ','OPLUS_BZ'  , gridname='fv_centers')
    call addfld ('OPLUS_BMAG', horiz_only , 'I', ' ','OPLUS_BMAG', gridname='fv_centers')

  end subroutine oplus_init

!-----------------------------------------------------------------------
  subroutine oplus_xport(tn,te,ti,un,vn,om,zg,o2,o1,n2,op_in,opnm_in, &
                         mbar,ui,vi,wi,pmid,op_out,opnm_out, &
                         i0,i1,j0,j1,nspltop,ispltop )
!
! All input fields from dpie_coupling are in "TIEGCM" format, i.e., 
! longitude (-180->180), vertical (bot2top), and units (CGS).
!
    use edyn_mpi,only:  mp_geo_halos,mp_pole_halos,setpoles
    use edyn_geogrid,only : glat, nlat, nlev
    use trsolv_mod, only : trsolv
!
! Transport O+ ion.
! March-May, 2015 B.Foster: Adapted from TIEGCM (oplus.F) for WACCM-X.
!
! Notes:
! - waccmx_opt='ionosphere' must be set in user_nl_cam for te,ti inputs to have values
!
! Args:
!
    integer,intent(in) :: &
      i0,                 & ! grid%ifirstxy
      i1,                 & ! grid%ilastxy
      j0,                 & ! grid%jfirstxy
      j1                    ! grid%jlastxy
    integer,intent(in) :: nspltop,ispltop
!
! Input fields without halo points (lon +/-180, vertical bot2top, CGS units):
!
    real(r8),intent(in) :: tn   (nlev,i0-2:i1+2,j0-2:j1+2) ! neutral temperature (deg K)
    real(r8),intent(in) :: te   (nlev,i0-2:i1+2,j0-2:j1+2) ! electron temperature (deg K)
    real(r8),intent(in) :: ti   (nlev,i0-2:i1+2,j0-2:j1+2) ! ion temperature (deg K)
    real(r8),intent(in) :: un   (nlev,i0-2:i1+2,j0-2:j1+2) ! neutral zonal wind (cm/s)
    real(r8),intent(in) :: vn   (nlev,i0-2:i1+2,j0-2:j1+2) ! neutral meridional wind (cm/s)
    real(r8),intent(in) :: om   (nlev,i0-2:i1+2,j0-2:j1+2) ! omega (1/s)
    real(r8),intent(in) :: o2   (nlev,i0-2:i1+2,j0-2:j1+2) ! o2 (mmr)
    real(r8),intent(in) :: o1   (nlev,i0-2:i1+2,j0-2:j1+2) ! o (mmr)
    real(r8),intent(in) :: n2   (nlev,i0-2:i1+2,j0-2:j1+2) ! n2 (mmr)
    real(r8),intent(in) :: mbar (nlev,i0-2:i1+2,j0-2:j1+2) ! mean molecular weight

    real(r8),intent(in) :: op_in(nlev,i0:i1,j0:j1) ! O+ density (cm^3)
    real(r8),intent(in) :: opnm_in(nlev,i0:i1,j0:j1) ! O+ density (cm^3) at time-1
    real(r8),intent(in) :: zg   (nlev,i0:i1,j0:j1) ! geopotential height (cm)
!
! Ion drifts from edynamo (also in tiegcm-format):
!
    real(r8),intent(in) :: ui(nlev,i0:i1,j0:j1)   ! zonal ion drift
    real(r8),intent(in) :: vi(nlev,i0:i1,j0:j1)   ! meridional ion drift
    real(r8),intent(in) :: wi(nlev,i0:i1,j0:j1)   ! vertical ion drift
    real(r8),intent(in) :: pmid(nlev)             ! pressure at midpoints (Pa)
!
! Output:
!
    real(r8),intent(out) ::       &
      op_out  (nlev,i0:i1,j0:j1), & ! O+ output
      opnm_out(nlev,i0:i1,j0:j1)    ! O+ output at time n-1
!
! Local:
!
    integer :: i,j,k,lat,jm1,jp1,jm2,jp2,lat0,lat1
    real(r8),dimension(i0:i1,j0:j1) :: &
      opflux,   & ! upward number flux of O+ (returned by sub oplus_flux)
      dvb         ! divergence of B-field
!
! Local inputs with added halo points in lat,lon:
! 
    real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2),target :: op, opnm

    real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2),target :: & 
      tr          ,&   ! Reduced temperature (.5*(tn+ti))
      tp          ,&   ! Plasma temperature N(O+)*(te+ti)
      dj          ,&   ! diffusion coefficients
      bvel        ,&   ! bvel @ j   = (B.U)*N(O+)
      diffj       ,&   ! (D/(H*DZ)*2.*TP+M*G/R)*N(O+)
      bdotdh_op   ,&   ! (b(h)*del(h))*phi
      bdotdh_opj  ,&   ! (b(h)*del(h))*phi
      bdotdh_diff ,&   ! (b(h)*del(h))*phi
      opnm_smooth      ! O+ at time-1, smoothed

    real(r8),dimension(nlev,i0:i1,j0:j1) :: & ! for saving to histories
      diag0,diag1,diag2,diag3,diag4,diag5,diag6,diag7,diag8,diag9,&
      diag10,diag11,diag12,diag13,diag14,diag15,diag16,diag17,&
      diag18,diag19,diag20,diag21,diag22,diag23,diag24,diag25,&
      diag26,diag27
    real(r8),dimension(nlev,i0:i1,j0-1:j1+1) :: hj ! scale height
    real(r8) :: gmr,dtime,dtx2,dtx2inv
    real(r8),dimension(nlev,i0:i1) :: &
      bdzdvb_op,   &
      hdz,         &
      tp1,         &
      tphdz0,      &
      tphdz1,      &
      djint,       &
      divbz,       &
      hdzmbz,      &
      hdzpbz,      &
      bdotu
!
! Arguments for tridiagonal solver trsolv (no halos):
    real(r8),dimension(nlev,i0:i1,j0:j1) :: &
      explicit,explicit_a,explicit_b,p_coeff,q_coeff,r_coeff

    real(r8),dimension(i0:i1) :: ubca, ubcb ! O+ upper boundary
    real(r8),parameter :: one=1._r8
    logical :: calltrsolv
!
! Pointers for multiple-field calls (e.g., mp_geo_halos)
    integer :: nfields
    real(r8),allocatable :: polesign(:)
    type(array_ptr_type),allocatable :: ptrs(:)

    real(r8) :: zpmid(nlev), opfloor
    real(r8),parameter :: opmin=3000.0_r8
!
! Execute:
!
    dtime = get_step_size() ! step size in seconds
    dtime = dtime / dble(nspltop)
    dtx2 = 2._r8*dtime
    dtx2inv = 1._r8/dtx2

    if ((is_first_step().or.is_first_restart_step()).and.ispltop==1) then
      if (masterproc) write(iulog,"('oplus: shapiro=',es12.4,' dtsmooth=',es12.4,' dtsmooth_div2=',es12.4)") &
        shapiro_const,dtsmooth,dtsmooth_div2
      if (masterproc) write(iulog,"('oplus: shr_const_g=',f8.3)") shr_const_g
    endif

    !
    ! zp,expz are declared in edyn_geogrid.F90, and allocated in sub 
    ! set_geogrid (edyn_init.F90). pmid was passed in here (bot2top)
    ! from dpie_coupling.
    !
    ! kbot is the k-index at the bottom of O+ transport calculations,
    !   corresponding to pressure pbot. 
    !
    if ((is_first_step().or.is_first_restart_step()).and.ispltop==1) then
       kloop: do k=1,nlev
          if ( pmid(k) <= pbot) then
             kbot = k
             exit kloop
          end if
       enddo kloop
       do k=1,nlev
          zp(k) = -log(pmid(k)*10._r8/p0)
          expz(k) = exp(-zp(k))
       enddo
       if (debug.and.masterproc) then
          write(iulog,"('oplus: kbot=',i4,' pmid(kbot)=',es12.4,' zp(kbot)=',es12.4)") &
               kbot,pmid(kbot),zp(kbot)
       endif
    endif

    if (kbot < 1) then
       call endrun('oplus_xport: kbot is not set')
    endif

    dzp = zp(nlev)-zp(nlev-1)  ! use top 2 levels (typically dzp=0.5)

    if (debug.and.masterproc) then
      write(iulog,"('oplus: nlev=',i3,' zp (bot2top)   =',/,(6es12.3))") nlev,zp
      write(iulog,"('oplus: nlev=',i3,' expz (bot2top) =',/,(6es12.3))") nlev,expz
      write(iulog,"('oplus: nlev=',i3,' dzp  =',/,(6es12.3))") nlev,dzp
    endif
!
! Set subdomain blocks from input (composition is in mmr):
!
!$omp parallel do private(i, j, k)
    do k=1,nlev
      do j=j0,j1
        do i=i0,i1
          op(k,i,j)   = op_in(k,i,j)
          opnm(k,i,j) = opnm_in(k,i,j)
        enddo
      enddo
    enddo

!
! Define halo points on inputs:
! WACCM has global longitude values at the poles (j=1,j=nlev)
! (they are constant for most, except the winds.)
!
! Set two halo points in lat,lon:
!   real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2),target :: tn,te,etc.
!
    nfields = 2
    allocate(ptrs(nfields),polesign(nfields))

    ptrs(1)%ptr => op ; ptrs(2)%ptr => opnm
    polesign = 1._r8
!
! mp_geo_halos first arg:
!     type(array_ptr_type) :: fmsub(nf) ! (lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)
!
    call mp_geo_halos(ptrs,1,nlev,i0,i1,j0,j1,nfields)
!
! Set latitude halo points over the poles (this does not change the poles).
! (the 2nd halo over the poles will not actually be used (assuming lat loops
!  are lat=2,nlat-1), because jp1,jm1 will be the pole itself, and jp2,jm2 
!  will be the first halo over the pole)
!
! mp_pole_halos first arg:
!   type(array_ptr_type) :: f(nf) ! (nlev,i0-2:i1+2,j0-2:j1+2)

    call mp_pole_halos(ptrs,1,nlev,i0,i1,j0,j1,nfields,polesign)
    deallocate(ptrs,polesign)

!
! Use below to exclude the poles (lat=2,nlat-1) from latitude scans.
!
    lat0 = j0
    lat1 = j1
    if (j0 == 1)    lat0 = 2
    if (j1 == nlat) lat1 = nlat-1
!
! Save input fields to WACCM histories. Sub savefld_waccm_switch converts
! fields from tiegcm-format to waccm-format before saving to waccm histories.
!
    call savefld_waccm_switch(tn(:,i0:i1,j0:j1),'OPLUS_TN',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(te(:,i0:i1,j0:j1),'OPLUS_TE',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(ti(:,i0:i1,j0:j1),'OPLUS_TI',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(un(:,i0:i1,j0:j1),'OPLUS_UN',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(vn(:,i0:i1,j0:j1),'OPLUS_VN',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(om(:,i0:i1,j0:j1),'OPLUS_OM',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(zg(:,i0:i1,j0:j1),'OPLUS_Z' ,nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(o2(:,i0:i1,j0:j1),'OPLUS_O2',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(o1(:,i0:i1,j0:j1),'OPLUS_O1',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(n2(:,i0:i1,j0:j1),'OPLUS_N2',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(op(:,i0:i1,j0:j1),'OPLUS_OP',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(ui(:,i0:i1,j0:j1),'OPLUS_UI',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(vi(:,i0:i1,j0:j1),'OPLUS_VI',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(wi(:,i0:i1,j0:j1),'OPLUS_WI',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(mbar(:,i0:i1,j0:j1),'OPLUS_MBAR',nlev,i0,i1,j0,j1)
    call savefld_waccm_switch(opnm(:,i0:i1,j0:j1),'OPLUS_OPNM',nlev,i0,i1,j0,j1)
!
! Initialize output op_out with input op at 1:kbot-1, to retain values from 
! bottom of column up to kbot. This routine will change (transport) these 
! outputs only from kbot to the top (nlev).
!
    op_out   = 0._r8 
    opnm_out = 0._r8 
    op_out  (1:kbot-1,i0:i1,j0:j1) = op  (1:kbot-1,i0:i1,j0:j1)
    opnm_out(1:kbot-1,i0:i1,j0:j1) = opnm(1:kbot-1,i0:i1,j0:j1)
!
! Sub oplus_flux returns upward number flux of O+ in opflux
! Output opflux(i,j) is 2d lon x lat subdomain:
!
    call oplus_flux(opflux,i0,i1,j0,j1)
    call savefld_waccm(opflux(i0:i1,j0:j1),'OPLUS_FLUX',1,i0,i1,j0,j1)
!
! Divergence of B (mag field) is returned by divb in dvb(i0:i1,j0:j1)
!
    call divb(dvb,i0,i1,j0,j1)
    call savefld_waccm(dvb(i0:i1,j0:j1),'OPLUS_DIVB',1,i0,i1,j0,j1)
!
! The solver will be called only if calltrsolv=true. It is sometimes
! set false when skipping parts of the code for debug purposes.
! 
    calltrsolv = .true.

    tr  = 0._r8
    tp  = 0._r8
    dj  = 0._r8
    hj  = 0._r8
    bvel= 0._r8
    diffj = 0._r8
    opnm_smooth = 0._r8
    diag0 =0._r8
    grav_cm = shr_const_g * 100._r8 ! m/s^2 -> cm/s^2
!
!----------------------- Begin first latitude scan ---------------------
    do lat=lat0,lat1
      jm2 = lat-2
      jm1 = lat-1
      jp1 = lat+1
      jp2 = lat+2
!
! as of April, 2015, TIEGCM incorrectly uses te+ti instead of tn+ti
! This has not been fixed in TIEGCM, because fixing it causes a tuning 
! problem (ask Hanli and Wenbin). For WACCM, it is correct as below.
! (see also tp)
!
!$omp parallel do private(i,k)
      do i=i0,i1
! 
! Reduced temperature (tpj in tiegcm):
! 'OPLUS_TR' (has constants at poles)
!
        do k=kbot,nlev
          tr(k,i,jm1) = 0.5_r8*(tn(k,i,jm1)+ti(k,i,jm1))
          tr(k,i,lat) = 0.5_r8*(tn(k,i,lat)+ti(k,i,lat))
          tr(k,i,jp1) = 0.5_r8*(tn(k,i,jp1)+ti(k,i,jp1))
        enddo
      enddo ! i=i0,i1
!
! rrk returns ambipolar diffusion coefficients in d(jm1),dj(lat),djp1(jp1):
! 'OPLUS_DJ' (has constants at poles)
!
      call rrk(                                            &
        tn(kbot:nlev,i0:i1,jm1),mbar(kbot:nlev,i0:i1,jm1), &
        o2(kbot:nlev,i0:i1,jm1),o1  (kbot:nlev,i0:i1,jm1), &
        n2(kbot:nlev,i0:i1,jm1),tr  (kbot:nlev,i0:i1,jm1), &
        dj(kbot:nlev,i0:i1,jm1),i0,i1,kbot,nlev)

      call rrk(                                            &
        tn(kbot:nlev,i0:i1,lat),mbar(kbot:nlev,i0:i1,lat), &
        o2(kbot:nlev,i0:i1,lat),o1  (kbot:nlev,i0:i1,lat), &
        n2(kbot:nlev,i0:i1,lat),tr  (kbot:nlev,i0:i1,lat), &
        dj(kbot:nlev,i0:i1,lat),i0,i1,kbot,nlev)

      call rrk(                                            &
        tn(kbot:nlev,i0:i1,jp1),mbar(kbot:nlev,i0:i1,jp1), &
        o2(kbot:nlev,i0:i1,jp1),o1  (kbot:nlev,i0:i1,jp1), &
        n2(kbot:nlev,i0:i1,jp1),tr  (kbot:nlev,i0:i1,jp1), &
        dj(kbot:nlev,i0:i1,jp1),i0,i1,kbot,nlev)
!
! Plasma temperature:
! 'OPLUS_TP0' (tp will get poles from jm1 and jp1)
!
!$omp parallel do private(i,k)
      do i=i0,i1
        do k=kbot,nlev
          tp(k,i,jm1) = te(k,i,jm1)+ti(k,i,jm1)
          tp(k,i,lat) = te(k,i,lat)+ti(k,i,lat)
          tp(k,i,jp1) = te(k,i,jp1)+ti(k,i,jp1)
        enddo
      enddo
      diag0(kbot:nlev,i0:i1,lat) = tp(kbot:nlev,i0:i1,lat)
!
! Add poles to diag0:
      if (j0==1.and.lat==2)         diag0(kbot:nlev,i0:i1,j0) = tp(kbot:nlev,i0:i1,jm1)
      if (j1==nlat.and.lat==nlat-1) diag0(kbot:nlev,i0:i1,j1) = tp(kbot:nlev,i0:i1,jp1)
!
! Neutral scale height:
! 'OPLUS_HJ' (has constants at poles)
!
!$omp parallel do private(i,k)
      do i=i0,i1
        do k=kbot,nlev
          hj(k,i,jm1) = gask * tn(k,i,jm1) / (mbar(k,i,jm1) * grav_cm)
          hj(k,i,lat) = gask * tn(k,i,lat) / (mbar(k,i,lat) * grav_cm)
          hj(k,i,jp1) = gask * tn(k,i,jp1) / (mbar(k,i,jp1) * grav_cm)
        enddo
      enddo
!
! bvel @ jm1 = (B.U)*N(O+)    (J-1)
! bvel @ j   = (B.U)*N(O+)      (J)
! bvel @ jp1 = (B.U)*N(O+)    (J+1)
! 'OPLUS_BVEL' (has constants at poles)
!
! Note bx,by,bz were set globally for all tasks by sub magfield
! (getapex.F90)
!
!$omp parallel do private(i,k)
      do i=i0,i1
        do k=kbot,nlev
          bvel(k,i,jm1) = &
            (bx(i,jm1)*un(k,i,jm1)+by(i,jm1)*vn(k,i,jm1)+   &
            hj(k,i,jm1)*bz(i,jm1)*om(k,i,jm1))*op(k,i,jm1)
          bvel(k,i,lat) = &
            (bx(i,lat)*un(k,i,lat)+by(i,lat)*vn(k,i,lat)+   &
            hj(k,i,lat)*bz(i,lat)*om(k,i,lat))*op(k,i,lat)
          bvel(k,i,jp1) = &
            (bx(i,jp1)*un(k,i,jp1)+by(i,jp1)*vn(k,i,jp1)+   &
            hj(k,i,jp1)*bz(i,jp1)*om(k,i,jp1))*op(k,i,jp1)
        enddo ! k=kbot,nlev
      enddo ! i=lon0,lon1
!
! Ambipolar diffusion is returned in diffj:
! 'OPLUS_DIFFJ' (will have constants at poles after this lat scan)
!
      call diffus(tp(kbot:nlev,i0:i1,jm1),op(kbot:nlev,i0:i1,jm1),hj(kbot:nlev,:,jm1), &
        diffj(kbot:nlev,i0:i1,jm1),i0,i1,kbot,nlev,lat)
      call diffus(tp(kbot:nlev,i0:i1,lat),op(kbot:nlev,i0:i1,lat),hj(kbot:nlev,:,lat), &
        diffj(kbot:nlev,i0:i1,lat),i0,i1,kbot,nlev,lat)
      call diffus(tp(kbot:nlev,i0:i1,jp1),op(kbot:nlev,i0:i1,jp1),hj(kbot:nlev,:,jp1), &
        diffj(kbot:nlev,i0:i1,jp1),i0,i1,kbot,nlev,lat)
!
! 'OPLUS_TP1' (constants at the poles)
!
!$omp parallel do private(i,k)
      do i=i0,i1
        do k=kbot,nlev
          tp(k,i,jm2) = op(k,i,jm2)*(te(k,i,jm2)+ti(k,i,jm2))
          tp(k,i,jm1) = tp(k,i,jm1)*op(k,i,jm1)
          tp(k,i,lat) = tp(k,i,lat)*op(k,i,lat)
          tp(k,i,jp1) = tp(k,i,jp1)*op(k,i,jp1)
          tp(k,i,jp2) = op(k,i,jp2)*(te(k,i,jp2)+ti(k,i,jp2))
        enddo
      enddo
!
! Latidinal shapiro smoother: opnm is O+ at time n-1.
! opnm_smooth will be used in explicit terms below.
! Smooth in latitude:
! 'OPNM_SMOOTH' (zero at poles)
!
!$omp parallel do private(i,k)
      do i=i0,i1
        do k=kbot,nlev
          opnm_smooth(k,i,lat) = opnm(k,i,lat)-shapiro_const* &
                                (opnm(k,i,jp2)+opnm(k,i,jm2)-4._r8*         &
                                (opnm(k,i,jp1)+opnm(k,i,jm1))+6._r8*        &
                                 opnm(k,i,lat))
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
    enddo ! end first latitude scan (lat=lat0,lat1)
!
!------------------------- End first latitude scan ---------------------
!
! Set pole values for opnm_smooth. Do this before savefld calls, so plots will
! include the poles. All other fields in 1st lat scan got values at the poles 
! via jm1,jp1 above.
!
    call setpoles(opnm_smooth(kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1)
!
! Save to history file (exclude halo points)
!
   call savefld_waccm_switch(tr   (:,i0:i1,j0:j1),'OPLUS_TR'   ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(dj   (:,i0:i1,j0:j1),'OPLUS_DJ'   ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(hj   (:,i0:i1,j0:j1),'OPLUS_HJ'   ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(bvel (:,i0:i1,j0:j1),'OPLUS_BVEL' ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diffj(:,i0:i1,j0:j1),'OPLUS_DIFFJ',nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag0(:,i0:i1,j0:j1),'OPLUS_TP0'  ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(tp   (:,i0:i1,j0:j1),'OPLUS_TP1'  ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(opnm_smooth(:,i0:i1,j0:j1),'OPNM_SMOOTH',nlev,i0,i1,j0,j1)
!
! Set halo points where needed.
!
    nfields = 5
    allocate(ptrs(nfields),polesign(nfields))
    ptrs(1)%ptr => dj ; ptrs(2)%ptr => bvel ; ptrs(3)%ptr => diffj
    ptrs(4)%ptr => tp ; ptrs(5)%ptr => opnm_smooth
    polesign = 1._r8

    call mp_geo_halos (ptrs,1,nlev,i0,i1,j0,j1,5)
    call mp_pole_halos(ptrs,1,nlev,i0,i1,j0,j1,5,polesign)

    deallocate(ptrs,polesign)

!----------------------- Begin second latitude scan --------------------
    bdotdh_op  = 0._r8
    bdotdh_opj = 0._r8

    do lat=lat0,lat1
      jm2 = lat-2
      jm1 = lat-1
      jp1 = lat+1
      jp2 = lat+2
!
! bdotdh_op = (B(H).DEL(H))*(D/(H*DZ)*TP+M*G/R)*N(O+)
! then bdotdh_op = d*bz*bdotdh_op
! real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2) :: diffj
! real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2) :: bdotdh_op
! 'BDOTDH_OP' (zero at the poles)
!
      call bdotdh(                    &
        diffj(kbot:nlev,i0:i1,jm1),   &
        diffj(kbot:nlev,:,lat    ),   & ! includes longitude halos
        diffj(kbot:nlev,i0:i1,jp1),   &
        bdotdh_op(kbot:nlev,i0:i1,lat),i0,i1,kbot,nlev,lat)
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev
          bdotdh_op(k,i,lat) = dj(k,i,lat)*bz(i,lat)*bdotdh_op(k,i,lat) ! BDOTDH_OP
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
!
! bdotdh_opjm1 = (B(H).DEL(H))*2.*TP*N(O+)    (J-1)
! bdotdh_opj   = (B(H).DEL(H))*2.*TP*N(O+)      (J)
! bdotdh_opjp1 = (B(H).DEL(H))*2.*TP*N(O+)    (J+1)
! 'BDOTDH_OPJ' (has reasonable non-constant values at poles)
!
      call bdotdh(                &
        tp(kbot:nlev,i0:i1,jm2),  &
        tp(kbot:nlev,:,jm1),      &
        tp(kbot:nlev,i0:i1,lat),  &
        bdotdh_opj(kbot:nlev,i0:i1,jm1),i0,i1,kbot,nlev,jm1)
      call bdotdh(                &
        tp(kbot:nlev,i0:i1,jm1),  &
        tp(kbot:nlev,:,lat),      &
        tp(kbot:nlev,i0:i1,jp1),  &
        bdotdh_opj(kbot:nlev,i0:i1,lat),i0,i1,kbot,nlev,lat)
      call bdotdh(                &
        tp(kbot:nlev,i0:i1,lat),  &
        tp(kbot:nlev,:,jp1),      &
        tp(kbot:nlev,i0:i1,jp2),  &
        bdotdh_opj(kbot:nlev,i0:i1,jp1),i0,i1,kbot,nlev,jp1)
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev
          bdotdh_opj(k,i,jm1) = bdotdh_opj(k,i,jm1)*dj(k,i,jm1)
          bdotdh_opj(k,i,lat) = bdotdh_opj(k,i,lat)*dj(k,i,lat)
          bdotdh_opj(k,i,jp1) = bdotdh_opj(k,i,jp1)*dj(k,i,jp1)
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
    enddo ! lat=j0,j1 (end second lat scan)
!
!------------------------ End second latitude scan ---------------------
!
! bdotdh_opj already has non-constant polar values, but bdotdh_op poles are zero.
! Sub setpoles will set poles to the zonal average of the latitude below each pole.
!
! This may not be necessary, but do it for plotting: 
    call setpoles(bdotdh_op(kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1)

   call savefld_waccm_switch(bdotdh_op (:,i0:i1,j0:j1),'BDOTDH_OP' ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(bdotdh_opj(:,i0:i1,j0:j1),'BDOTDH_OPJ',nlev,i0,i1,j0,j1)
!
! Note mp_geo_halos will overwrite jm1,jp1 that was set above.
! bdotdh_opj needs longitude halos for the bdotdh call below.
!
! real(r8),dimension(nlev,i0-2:i1+2,j0-2:j1+2),target :: bdotdh_op,opj
!
    allocate(ptrs(1))
    ptrs(1)%ptr => bdotdh_opj
    call mp_geo_halos (ptrs,1,nlev,i0,i1,j0,j1,1)
    call mp_pole_halos(ptrs,1,nlev,i0,i1,j0,j1,1,(/1._r8/))
    deallocate(ptrs)
!
!----------------------- Begin third latitude scan ---------------------
!
    bdotdh_diff = 0._r8
    bdzdvb_op   = 0._r8
    explicit(1:nlev,i0:i1,j0:j1)    = 0._r8 ; explicit_a(1:nlev,i0:i1,j0:j1)=0._r8 ; explicit_b(1:nlev,i0:i1,j0:j1)=0._r8
    hdz         = 0._r8
    tphdz0      = 0._r8
    tphdz1      = 0._r8
    djint       = 0._r8
    divbz       = 0._r8
    hdzmbz      = 0._r8
    hdzpbz      = 0._r8
    p_coeff(1:nlev,i0:i1,j0:j1)     = 0._r8
    q_coeff(1:nlev,i0:i1,j0:j1)     = 0._r8
    r_coeff(1:nlev,i0:i1,j0:j1)     = 0._r8
    bdotu       = 0._r8

    diag1  = 0._r8 ; diag2 = 0._r8 ; diag3 = 0._r8 ; diag4 = 0._r8 ; diag5 = 0._r8
    diag6  = 0._r8 ; diag7 = 0._r8 ; diag8 = 0._r8 ; diag9 = 0._r8 ; diag10= 0._r8  
    diag11 = 0._r8 ; diag12= 0._r8 ; diag13= 0._r8 ; diag14= 0._r8 ; diag15= 0._r8  
    diag16 = 0._r8 ; diag17= 0._r8 ; diag18= 0._r8 ; diag19= 0._r8 ; diag20= 0._r8  
    diag21 = 0._r8 ; diag22= 0._r8 ; diag23= 0._r8 ; diag24= 0._r8 ; diag25= 0._r8
    diag26 = 0._r8 ; diag27= 0._r8

!
! gmr = G*M(O+)/(2.*R)
!
    gmr = grav_cm*rmass_op/(2._r8*gask)

!
! Globally, this loop is lat=2,nlat-1 (i.e., skipping the poles)
!
    do lat=lat0,lat1
      jm2 = lat-2
      jm1 = lat-1 ! this will be south pole for southern pes (j==1)
      jp1 = lat+1 ! this will be north pole for northern pes (j==nlat)
      jp2 = lat+2
!
! bdotdh_opj = (B(H).DEL(H))*D*(B(H).DEL(H))*2.*TP*N(O+)   (J)
! 'BDOTDH_DIFF' (zero at the poles)
!
      call bdotdh(                       &
        bdotdh_opj(kbot:nlev,i0:i1,jm1), &
        bdotdh_opj(kbot:nlev,:,lat),     & ! includes longitude halos
        bdotdh_opj(kbot:nlev,i0:i1,jp1), &
        bdotdh_diff(kbot:nlev,i0:i1,lat),i0,i1,kbot,nlev,lat) ! BDOTDH_DIFF
!
! bdzdvb_op = (BZ*D/(H*DZ)+DIV(*B))*S2
! bdzdvb returns bdzdvb_op(k,i).
! 'BDZDVB_OP' (zero at the poles)
!
! real(r8),dimension(i0:i1,j0:j1) :: dvb
! real(r8),dimension(nlev,i0:i1,j0-1:j1+1) :: hj ! scale height
! real(r8),dimension(nlev,i0:i1) :: bdzdvb_op

! subroutine bdzdvb(phi,dvb,h,ans,lev0,lev1,lon0,lon1,lat)
!   real(r8),intent(in) :: dvb(lon0:lon1)
!   real(r8),dimension(lev0:lev1,lon0:lon1),intent(in)    :: phi,h
!   real(r8),dimension(lev0:lev1,lon0:lon1),intent(out)   :: ans
!
      call bdzdvb(bdotdh_opj(kbot:nlev,i0:i1,lat),dvb(:,lat),hj(kbot:nlev,i0:i1,lat), &
        bdzdvb_op(kbot:nlev,i0:i1),kbot,nlev,i0,i1,lat)
      diag1(:,i0:i1,lat) = bdzdvb_op(:,i0:i1) ! BDZDVB_OP
!
! Collect explicit terms:
! 'EXPLICIT0' (this will have poles set after third lat scan, before
!              plotting. The poles will be constant in longitude, and
!              may differ structurally from adjacent latitudes. 
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev
          explicit(k,i,lat) = -one*(bdzdvb_op(k,i)+bdotdh_diff(k,i,lat)+ &
            bdotdh_op(k,i,lat))
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
      diag2(:,i0:i1,lat) = explicit(:,i0:i1,lat) ! EXPLICIT0
!
! Ion drifts are interpolated to midpoints (is this necessary in WACCM?).
!
! Need lon,lat halos for op, bvel, and bmod2
! op,bvel halos were set above, bmod2 was set in magfield (getapex.F90)
! (ui,vi,wi halos are not used here.)
!
! bmod2 halos are set in sub magfield (getapex.F90), including nlat-1,nlat,nlat+1,
!   and 1 halo point in longitude. Note bmod2 is global in lon and lat for all pe's.
! use getapex,only: bmod2  ! (0:nlonp1,jspole-1:jnpole+1)
!
! When looping lat=2,nlat-1, this explicit has zero pole values,
! but there are still problems at processor longitude boundaries,
! especially near the south pole:
! 'EXPLICIT1' (zero at the poles)
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev-1
!
! Original TIEGCM statement:
!         explicit(k,i) = explicit(k,i)+1._r8/(2._r8*re)*        &
!           (1._r8/(cs(lat)*dlamda)*(bx(i,lat)*                  &
!           (bvel(k,i+1,lat)-bvel(k,i-1,lat))+                   &
!           0.5_r8*(ui(k,i,lat)+ui(k+1,i,lat))*bmod2(i,lat)**2*  &
!           (op(k,i+1,lat)/bmod2(i+1,lat)**2-                    &
!            op(k,i-1,lat)/bmod2(i-1,lat)**2))+                  &
!
!           1._r8/dphi*(by(i,lat)*(bvel(k,i,jp1)-bvel(k,i,jm1))+ &
!           0.5_r8*(vi(k,i,lat)+vi(k+1,i,lat))*bmod2(i,lat)**2*  &
!           (op(k,i,jp1)/bmod2(i,jp1)**2-                        &
!            op(k,i,jm1)/bmod2(i,jm1)**2)))
!
! Break it into two pieces and put together for debug:
!
! 'EXPLICITa'
         explicit_a(k,i,lat) = (bx(i,lat)*                       &
            (bvel(k,i+1,lat)-bvel(k,i-1,lat))+                   &
            0.5_r8*(ui(k,i,lat)+ui(k+1,i,lat))*bmod2(i,lat)**2*  &
            (op(k,i+1,lat)/bmod2(i+1,lat)**2-                    &
             op(k,i-1,lat)/bmod2(i-1,lat)**2))                   
!
! 'EXPLICITb'
!
         explicit_b(k,i,lat) = &
            (by(i,lat)*(bvel(k,i,jp1)-bvel(k,i,jm1))+            &
            0.5_r8*(vi(k,i,lat)+vi(k+1,i,lat))*bmod2(i,lat)**2*  &
            (op(k,i,jp1)/bmod2(i,jp1)**2-                        &
             op(k,i,jm1)/bmod2(i,jm1)**2))
!
! 'EXPLICIT1'
! explicit will receive polar values after this latitude scan.
!
         explicit(k,i,lat) = explicit(k,i,lat)+1._r8/(2._r8*re)* &
           (1._r8/(cs(lat)*dlamda)*explicit_a(k,i,lat)+          &
           1._r8/dphi*explicit_b(k,i,lat))

!
! explicit is bad at i=1,72,73,144 near south pole (npole appears to be ok)
! This does not appear to adversely affect the final O+ output, and TIEGCM
! has the same high magnitudes, so am ignoring this for now. The high magnitudes
! are near the south pole, at processor longitude boundaries (implicating an error
! with longitude halo points).
!
         if (debug) then
            if (explicit(k,i,lat) < -300._r8 .or. explicit(k,i,lat) > 300._r8) then
               write(iulog,"('>>> bad explicit: k,i,lat=',3i4,' explicit=',es12.4)") &
                    k,i,lat,explicit(k,i,lat)
               write(iulog,"('  cs(lat)	       =',3es12.4)") cs(lat)
               write(iulog,"('  op(k,i-1:i+1,lat)  =',3es12.4)") op(k,i-1:i+1,lat)
               write(iulog,"('  op(k,i,jm1:jp1)    =',3es12.4)") op(k,i,jm1:jp1)
               write(iulog,"('  bvel(k,i-1:i+1,lat)=',3es12.4)") bvel(k,i-1:i+1,lat)
               write(iulog,"('  bvel(k,i,jm1:jp1)  =',3es12.4)") bvel(k,i,jm1:jp1)
               write(iulog,"('  bmod2(i-1:i+1,lat) =',3es12.4)") bmod2(i-1:i+1,lat)
               write(iulog,"('  bmod2(i,jm1:jp1)   =',3es12.4)") bmod2(i,jm1:jp1)
               write(iulog,"('  ui(k:k+1,i,lat)    =',2es12.4)") ui(k:k+1,i,lat)
               write(iulog,"('  vi(k:k+1,i,lat)    =',2es12.4)") vi(k:k+1,i,lat)
               write(iulog,"('  bx,by(i,lat)       =',2es12.4)") bx(i,lat),by(i,lat)
            endif
	 endif

        enddo ! k=kbot,nlev-1
      enddo ! i=i0,i1

!$omp parallel do private( k )
      do k=kbot,nlev
        diag25(k,i0:i1,lat) = bmod2(i0:i1,lat)    ! BMOD2 (redundant in vertical)
      enddo
      diag26(:,i0:i1,lat) = explicit_a(:,i0:i1,lat) ! EXPLICITa
      diag27(:,i0:i1,lat) = explicit_b(:,i0:i1,lat) ! EXPLICITb
      diag3 (:,i0:i1,lat) = explicit  (:,i0:i1,lat) ! EXPLICIT1

!$omp parallel do private( i )
      do i=i0,i1
        dvb(i,lat) = dvb(i,lat)/bz(i,lat)
      enddo ! i=i0,i1
 
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev
          hdz(k,i) = 1._r8/(hj(k,i,lat)*dzp)
          tp1(k,i) = 0.5_r8*(ti(k,i,lat)+te(k,i,lat))
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1

!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev-1
          tphdz1(k+1,i) = 2._r8*tp1(k+1,i)*(0.5_r8*(hdz(k,i)+hdz(k+1,i)))+gmr
          tphdz0(k+1,i) = 2._r8*tp1(k  ,i)*(0.5_r8*(hdz(k,i)+hdz(k+1,i)))-gmr
        enddo ! k=kbot,nlev-1
      enddo ! i=lon0,lon1
!
! Upper and lower boundaries:
! Both TPHDZ0 and TPHDZ1 are zero at the poles.
!
! 5/9/15: Appears to be a problem in TPHDZ0,1 near kbot, maybe is zero?
!
!$omp parallel do private( i )
      do i=i0,i1
        tphdz1(kbot,i) = 2._r8*tp1(kbot,i)*                             &
                         (1.5_r8*hdz(kbot,i)-0.5_r8*hdz(kbot+1,i))+gmr
        tphdz1(nlev,i) = 2._r8*(2._r8*tp1(nlev-1,i)-tp1(nlev-2,i))*     &
                         (1.5_r8*hdz(nlev-1,i)-0.5_r8*hdz(nlev-2,i))+gmr
        tphdz0(kbot,i) = 2._r8*(2._r8*tp1(kbot,i)-tp1(kbot+1,i))*       &
                         (1.5_r8*hdz(kbot,i)-0.5_r8*hdz(kbot+1,i))-gmr
        tphdz0(nlev,i) = 2._r8*tp1(nlev-1,i)*                           &
                         (1.5_r8*hdz(nlev-1,i)-0.5_r8*hdz(nlev-2,i))-gmr
      enddo ! i=i0,i1
      diag4(:,i0:i1,lat) = tphdz0(:,i0:i1) ! TPHDZ0
      diag5(:,i0:i1,lat) = tphdz1(:,i0:i1) ! TPHDZ1
!
! djint = dj diffusion at interfaces:
! 'DJINT' (zero at the poles - messes up the plots - may give
!          diag6 polar values after the lat scan, before plotting)
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev-1
          djint(k+1,i) = 0.5_r8*(dj(k,i,lat)+dj(k+1,i,lat))
        enddo 
        djint(kbot,i) = (1.5_r8*dj(kbot  ,i,lat)-0.5_r8*dj(kbot+1,i,lat))
        djint(nlev,i) = (1.5_r8*dj(nlev-1,i,lat)-0.5_r8*dj(nlev-2,i,lat))
      enddo ! i=i0,i1
      diag6(:,i0:i1,lat) = djint(:,i0:i1) ! DJINT
!
! divbz = (DIV(B)+(DH*D*BZ)/(D*BZ)
! 'DIVBZ' Field appears as a line following mins along magnetic equator (zero at poles)
! Field may be zero at kbot? Lat slices look strange between +/- 14 deg lat.
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev
          divbz(k,i) =                                                 &
            dvb(i,lat)+1._r8/(re*dj(k,i,lat)*bz(i,lat)**2)*(bx(i,lat)/ &
            cs(lat)*(dj(k,i+1,lat)*bz(i+1,lat)-dj(k,i-1,lat)*          &
            bz(i-1,lat))/(2._r8*dlamda)+by(i,lat)*(dj(k,i,jp1)*        &
            bz(i,jp1)-dj(k,i,jm1)*bz(i,jm1))/(2._r8*dphi))
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
      diag7(:,i0:i1,lat) = divbz(:,i0:i1) ! DIVBZ
!
! hdzmbz = (1./(H*DZ)-(DIV(B)+DH*D*BZ/(D*BZ))/(2*BZ))*BZ**2
! hdzpbz = (1./(H*DZ)+(DIV(B)+DH*D*BZ/(D*BZ))/(2*BZ))*BZ**2
! 'HDZMBZ' and 'HDZPBZ' are zero at the poles.
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev
          hdzmbz(k,i) = (hdz(k,i)-0.5_r8*divbz(k,i))*bz(i,lat)**2
          hdzpbz(k,i) = (hdz(k,i)+0.5_r8*divbz(k,i))*bz(i,lat)**2
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
      diag8(:,i0:i1,lat) = hdzmbz(:,i0:i1) ! HDZMBZ
      diag9(:,i0:i1,lat) = hdzpbz(:,i0:i1) ! HDZPBZ
!
! Sum O+ at time n-1 to explicit terms: N(O+)/(2*DT) (N-1)
! 'EXPLICIT2' (zero at the poles)
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev
          explicit(k,i,lat) = explicit(k,i,lat)-(opnm_smooth(k,i,lat)-shapiro_const* &
            (opnm_smooth(k,i+2,lat)+opnm_smooth(k,i-2,lat)-4._r8*      &
            (opnm_smooth(k,i+1,lat)+opnm_smooth(k,i-1,lat))+6._r8*     &
             opnm_smooth(k,i,lat)))*dtx2inv
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
      diag10(:,i0:i1,lat) = explicit(:,i0:i1,lat) ! EXPLICIT2
!
! Begin coefficients p_coeff, q_coeff, r_coeff
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev-1
          p_coeff(k,i,lat) =   hdzmbz(k,i)*djint(k  ,i)*tphdz0(k  ,i)
          q_coeff(k,i,lat) = -(hdzpbz(k,i)*djint(k+1,i)*tphdz0(k+1,i)+  &
                               hdzmbz(k,i)*djint(k  ,i)*tphdz1(k  ,i))
          r_coeff(k,i,lat) =   hdzpbz(k,i)*djint(k+1,i)*tphdz1(k+1,i)
        enddo ! k=kbot,nlev-1
      enddo ! i=i0,i1

      diag11(:,i0:i1,lat) = p_coeff(:,i0:i1,lat) ! P_COEFF0 (zero at poles)
      diag12(:,i0:i1,lat) = q_coeff(:,i0:i1,lat) ! Q_COEFF0 (zero at ubc)
      diag13(:,i0:i1,lat) = r_coeff(:,i0:i1,lat) ! R_COEFF0 (zero at ubc)
!
! bdotu = B.U
! Introducing neutral winds.
! Am not using 0.5*(om(k)+om(k+1)) here because waccm omega is on midpoints (?)
!           (tiegcm has 0.5*(w(k,i,j0)+w(k+1,i,j0)) )
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev
          bdotu(k,i) = bx(i,lat)*un(k,i,lat)+by(i,lat)*vn(k,i,lat)+ &
            hj(k,i,lat)*bz(i,lat)*om(k,i,lat)
        enddo ! k=kbot,nlev
      enddo ! i=i0,i1
      diag14(:,i0:i1,lat) = bdotu(:,i0:i1) ! BDOTU
!
! Continue coefficients with vertical ion drift:
! wi is converted from interfaces to midpoints (first use of wi).
! The p,q,r coeffs are still zero at top boundary k=nlev, and at poles.
!
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev-2

          p_coeff(k+1,i,lat) = p_coeff(k+1,i,lat)+(bz(i,lat)*bdotu(k,i)+  &
            0.5_r8*(wi(k+1,i,lat)+wi(k+2,i,lat)))*0.5_r8*hdz(k+1,i)

          q_coeff(k,i,lat) = q_coeff(k,i,lat)-0.5_r8*(wi(k,i,lat)+wi(k+1,i,lat))*6._r8/re

          r_coeff(k,i,lat) = r_coeff(k,i,lat)-(bz(i,lat)*bdotu(k+1,i)+  &
            0.5_r8*(wi(k,i,lat)+wi(k+1,i,lat)))*0.5_r8*hdz(k,i)

        enddo ! k=kbot,nlev-1
      enddo ! i=i0,i1

      diag22(:,i0:i1,lat) = p_coeff(:,i0:i1,lat) ! P_COEFF0a
      diag23(:,i0:i1,lat) = q_coeff(:,i0:i1,lat) ! Q_COEFF0a
      diag24(:,i0:i1,lat) = r_coeff(:,i0:i1,lat) ! R_COEFF0a
!
! Upper (nlev) and lower (kbot) boundaries of p,q,r_coeff:
! (convert wi to midpoints)
!
! tiegcm considers nlev-1 to be the top level. Do it tiegcm-style here,
! and then extrapolate to nlev.
!
!$omp parallel do private( i )
      do i=i0,i1
        p_coeff(kbot,i,lat) = p_coeff(kbot,i,lat)+(bz(i,lat)*  &  ! reset p_coeff lbc 
          (2._r8*bdotu(kbot,i)-bdotu(kbot+1,i))+               &
          0.5_r8*(wi(kbot,i,lat)+wi(kbot+1,i,lat)))*0.5_r8*hdz(kbot,i)

        q_coeff(nlev-1,i,lat) = q_coeff(nlev-1,i,lat)-         & 
          0.5_r8*(wi(nlev,i,lat)+wi(nlev-1,i,lat))*6._r8/re 

        r_coeff(nlev-1,i,lat) = r_coeff(nlev-1,i,lat)-(bz(i,lat)*  &  
          (2._r8*bdotu(nlev-1,i)-bdotu(nlev-2,i))+                 &  
          0.5_r8*(wi(nlev,i,lat)+wi(nlev-1,i,lat)))*0.5_r8*hdz(nlev-1,i)
      enddo ! i=i0,i1
!
! Extrapolate to top level (tiegcm does not do this):
!
      p_coeff(nlev,i0:i1,lat) = 1.5_r8*p_coeff(nlev-1,i0:i1,lat)- &
                                0.5_r8*p_coeff(nlev-2,i0:i1,lat)
      q_coeff(nlev,i0:i1,lat) = 1.5_r8*q_coeff(nlev-1,i0:i1,lat)- &
                                0.5_r8*q_coeff(nlev-2,i0:i1,lat)
      r_coeff(nlev,i0:i1,lat) = 1.5_r8*r_coeff(nlev-1,i0:i1,lat)- &
                                0.5_r8*r_coeff(nlev-2,i0:i1,lat)
!
! All P,Q,R are zero at the poles. Polar values will be set after third lat scan.
      diag15(:,i0:i1,lat) = p_coeff(:,i0:i1,lat) ! P_COEFF1 (zero at ubc and poles)
      diag17(:,i0:i1,lat) = r_coeff(:,i0:i1,lat) ! R_COEFF1 (ok at ubc, zero at poles)
!
! Additions to Q coefficients (includes q_coeff lbc,ubc):
!$omp parallel do private( i, k )
      do i=i0,i1
        do k=kbot,nlev
          q_coeff(k,i,lat) = q_coeff(k,i,lat)-bdotu(k,i)*dvb(i,lat)*bz(i,lat)-dtx2inv
        enddo ! k=kbot,nlev-1
      enddo ! i=i0,i1
!
! Plot Q_COEFF1 after ubc has been set.
      diag16(:,i0:i1,lat) = q_coeff(:,i0:i1,lat) ! Q_COEFF1 (ok at ubc, zero at poles)
!
! Upper boundary condition for O+:
!$omp parallel do private( i )
      do i=i0,i1
        ubca(i) = 0._r8
        ubcb(i) = -bz(i,lat)**2*djint(nlev,i)*tphdz0(nlev,i)-ubca(i)
        ubca(i) = -bz(i,lat)**2*djint(nlev,i)*tphdz1(nlev,i)+ubca(i)
!
! Q = Q+B/A*R
        q_coeff(nlev,i,lat) = q_coeff(nlev,i,lat)+ubcb(i)/ubca(i)* &
          r_coeff(nlev,i,lat)
!
! F = F -R/A*PHI
        explicit(nlev,i,lat) = explicit(nlev,i,lat)-opflux(i,lat)* &  ! explicit ubc
          r_coeff(nlev,i,lat)/ubca(i)
        r_coeff(nlev,i,lat) = 0._r8                   ! r_coeff ubc is reset to zero
      enddo ! i=i0,i1
!
! Ubc of EXPLICIT3 has a stripe along the mag equator, unlike the level below.
!
      diag18(:,i0:i1,lat) = explicit(:,i0:i1,lat) ! EXPLICIT3 (ubc ok, zero at poles)
      diag19(:,i0:i1,lat) = p_coeff(:,i0:i1,lat)  ! P_COEFF2  (zero at ubc, zero at poles)
      diag20(:,i0:i1,lat) = q_coeff(:,i0:i1,lat)  ! Q_COEFF2  (ubc ok, zero at poles)
      diag21(:,i0:i1,lat) = r_coeff(:,i0:i1,lat)  ! R_COEFF2  (zero at ubc, zero at poles)
!
! At this point, TIEGCM calculates "sources and sinks" xiop2p and xiop2d.
! Also calculates op_loss, which is subtracted from q_coeff.
! Then TIEGCM "Add source term to RHS (explicit terms)", and calculates
! lower boundary condition N(O+) = Q/L (q_coeff, explicit, p_coeff), and
! finally calls trsolv.
!
 300 continue
    enddo ! end third latitude scan (lat=lat0,lat1)
!
!------------------------ End third latitude scan ---------------------

!
! Set poles for selected diagnostics:
!
    call setpoles(diag26(kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1) ! EXPLICITa
    call setpoles(diag27(kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1) ! EXPLICITb
    call setpoles(diag2 (kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1) ! EXPLICIT0
    call setpoles(diag3 (kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1) ! EXPLICIT1
    call setpoles(diag6 (kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1) ! DJINT
!
! All tasks have global 2d bmod2.
! bmod2 was set by sub magfield (getapex.F90) 
!   allocate(bmod2(0:nlonp1,jspole-1:jnpole+1))
! Copy bmod2 poles to diagnostic array.
!
!$omp parallel do private( i, k )
    do i=i0,i1
      do k=kbot,nlev
        diag25(k,i,j0)   = bmod2(i,j0)
        diag25(k,i,j1)   = bmod2(i,j1)
      enddo
    enddo
   call savefld_waccm_switch(diag25,'BMOD2'     ,nlev,i0,i1,j0,j1)
!
! Assign polar values to coefficients for trsolv.
!
    call setpoles(explicit(kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1)
    call setpoles(p_coeff (kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1)
    call setpoles(q_coeff (kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1)
    call setpoles(r_coeff (kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1)
!
! Call solver, defining O+ output op_out:
!
! Its best not to call this unless the coefficients and explicit terms
! have been properly set in the third latitude scan above (e.g., during 
! "goto 300" debugging above, where the coeffs may not have been calculated).
!
    if (calltrsolv) then

!$omp parallel do private( lat )
      do lat=j0,j1

        call trsolv(p_coeff (kbot:nlev,i0:i1,lat), &
                    q_coeff (kbot:nlev,i0:i1,lat), &
                    r_coeff (kbot:nlev,i0:i1,lat), &
                    explicit(kbot:nlev,i0:i1,lat), &
                    op_out  (kbot:nlev,i0:i1,lat), &
                    kbot,nlev,kbot,nlev,i0,i1 )

      enddo

      call savefld_waccm_switch(op_out,'OP_SOLVE',nlev,i0,i1,j0,j1)

    else ! trsolv not called (debug only)
      op_out  (kbot:nlev,i0:i1,j0:j1) = op  (kbot:nlev,i0:i1,j0:j1)
      opnm_out(kbot:nlev,i0:i1,j0:j1) = opnm(kbot:nlev,i0:i1,j0:j1)
    endif ! calltrsolv
!
! Write fields from third latitude scan to waccm history:
!
   call savefld_waccm_switch(explicit,'EXPLICIT',nlev,i0,i1,j0,j1) ! non-zero at ubc
   call savefld_waccm_switch(p_coeff ,'P_COEFF' ,nlev,i0,i1,j0,j1) ! zero at ubc?
   call savefld_waccm_switch(q_coeff ,'Q_COEFF' ,nlev,i0,i1,j0,j1) ! non-zero at ubc
   call savefld_waccm_switch(r_coeff ,'R_COEFF' ,nlev,i0,i1,j0,j1) ! is set zero at ubc

   call savefld_waccm_switch(bdotdh_diff(:,i0:i1,j0:j1), 'BDOTDH_DIFF',nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag1 ,'BDZDVB_OP',nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag2 ,'EXPLICIT0',nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag26,'EXPLICITa',nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag27,'EXPLICITb',nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag3 ,'EXPLICIT1',nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag4 ,'TPHDZ0'   ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag5 ,'TPHDZ1'   ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag6 ,'DJINT'    ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag7 ,'DIVBZ'    ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag8 ,'HDZMBZ'   ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag9 ,'HDZPBZ'   ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag10,'EXPLICIT2',nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag11,'P_COEFF0' ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag12,'Q_COEFF0' ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag13,'R_COEFF0' ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag14,'BDOTU'    ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag15,'P_COEFF1' ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag16,'Q_COEFF1' ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag17,'R_COEFF1' ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag18,'EXPLICIT3',nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag19,'P_COEFF2' ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag20,'Q_COEFF2' ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag21,'R_COEFF2' ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag22,'P_COEFF0a',nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag23,'Q_COEFF0a',nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(diag24,'R_COEFF0a',nlev,i0,i1,j0,j1)
!
!------------------------------------------------------------------------
!
! Filter O+ output from solver:
! (TIMEGCM calls both filters, whereas TIEGCM calls only filter2)
!
!   call filter1_op(op_out(kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1)
!
    call filter2_op(op_out(kbot:nlev,i0:i1,j0:j1),kbot,nlev,i0,i1,j0,j1)
!
!----------------------- Begin fourth latitude scan ---------------------
!
!$omp parallel do private(lat, i, k, opfloor)
    do lat=j0,j1
      do i=i0,i1
        do k=kbot,nlev
          opnm_out(k,i,lat) = dtsmooth*op(k,i,lat)+dtsmooth_div2* &
            (opnm(k,i,lat)+op_out(k,i,lat))
        enddo
      enddo
!
! Insure non-negative O+ output:
      do i=i0,i1
        do k=kbot,nlev
          if (op_out  (k,i,lat) < 1.e-5_r8) op_out  (k,i,lat) = 1.e-5_r8
          if (opnm_out(k,i,lat) < 1.e-5_r8) opnm_out(k,i,lat) = 1.e-5_r8
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! Enforce O+ minimum if enforce_opfloor is true.
! Opfloor is Stan's "smooth floor" (product of two Gaussians, 
!   dependent on latitude and pressure level) (opmin=3000.0):
!
      if (enforce_floor) then
         zpmid(kbot:nlev) = log(50.e-6_r8/pmid(kbot:nlev)) ! tgcm levs -- maybe done once at init time
         do k=kbot,nlev
            opfloor = opmin*exp(-(glat(lat)/90.0_r8)**2/0.3_r8) &
                 *exp(-((zpmid(k)-4.25_r8)/zpmid(nlev))**2/0.1_r8)
            do i=i0,i1
               if (op_out(k,i,lat) < opfloor) then
                  op_out(k,i,lat) = opfloor
               endif ! opout < opfloor
            enddo ! i=lon0,lon1
         enddo ! k=lev0,lev1-1
      endif ! enforce_opfloor

    enddo ! lat=lat0,lat1

!
! Save O+ output to WACCM history (cm^3):
   call savefld_waccm_switch(op_out  (:,i0:i1,j0:j1),'OP_OUT'  ,nlev,i0,i1,j0,j1)
   call savefld_waccm_switch(opnm_out(:,i0:i1,j0:j1),'OPNM_OUT',nlev,i0,i1,j0,j1)
  end subroutine oplus_xport
!-----------------------------------------------------------------------
  subroutine oplus_flux(opflux,lon0,lon1,lat0,lat1)
!
! Calculate O+ number flux for sub oplus_xport.
! Flux is returned in opflux(lon0:lon1,lat0:lat1).
!
! alatm: geomagnetic latitude at each geographic grid point (radians)
  use getapex,only: alatm ! (nlonp1,jspole:jnpole)
!
! Args:
    integer,intent(in) :: lon0,lon1,lat0,lat1
    real(r8),intent(out) :: opflux(lon0:lon1,lat0:lat1)
!
! Local:
    integer :: i,j
    real(r8),dimension(lon0:lon1,lat0:lat1) :: chi  ! solar zenith angle
    real(r8),parameter :: &
      phid =  2.0e8_r8,   &
      phin = -2.0e8_r8,   &
!     phin = 0._r8,       &
      ppolar = 0._r8
    real(r8) :: a(lon0:lon1)    
    real(r8) :: fed(lon0:lon1)  
    real(r8) :: fen(lon0:lon1)  
!
! Set some paramaters:
    pi  = 4._r8*atan(1._r8)
    rtd = 45._r8/atan(1._r8)
!
! Sub get_zenith calls sub zenith (..cam/src/physics/cam/zenith.F90)
    call get_zenith(chi,lon0,lon1,lat0,lat1)
!
! Latitude scan:
    do j=lat0,lat1
!
! Longitude loop:
      do i=lon0,lon1
        if (abs(alatm(i,j))-pi/24._r8>=0._r8) then
          a(i) = 1._r8
        else
          a(i)=.5_r8*(1._r8+sin(pi*(abs(alatm(i,j))-pi/48._r8)/(pi/24._r8)))
          if (a(i) < 0.05_r8) a(i) = 0.05_r8
        endif
        fed(i) = phid*a(i)
        fen(i) = phin*a(i)
        if (chi(i,j)-0.5_r8*pi >= 0._r8) then
          opflux(i,j) = fen(i)
        else
          opflux(i,j) = fed(i)
        endif
        if ((chi(i,j)*rtd-80._r8)*(chi(i,j)*rtd-100._r8) < 0._r8) then
          opflux(i,j) = .5_r8*(fed(i)+fen(i))+.5_r8*(fed(i)-fen(i))* &
            cos(pi*(chi(i,j)*rtd-80._r8)/20._r8)
        endif
!
! Add ppolar if magnetic latitude >= 60 degrees:
! QUESTION: is the 60 deg here related to critical angles crit(2) in tiegcm?
! 3/15/15: opflux is comparable to tiegcm.
!
        if (abs(alatm(i,j))-pi/3._r8 >= 0._r8) &
          opflux(i,j) = opflux(i,j)+ppolar
      enddo ! i=lon0,lon1
    enddo ! j=lat0,lat1
!
  end subroutine oplus_flux
!-----------------------------------------------------------------------
  subroutine get_zenith(chi,i0,i1,j0,j1)
!
! Get solar zenith angle chi(i0:i1,j0:j1) (radians)
! Subroutine zenith returns cos(chi) at each (i,j)
! Note glon(i0:i1) from edyn_init is in -180 -> +180 (TIEGCM mode)
!
    use time_manager,only : get_curr_calday
    use edyn_geogrid,only : glon,glat
    use orbit,       only : zenith
!
! Args:
    integer,intent(in) :: i0,i1,j0,j1
    real(r8),intent(out) :: chi(i0:i1,j0:j1)
!
! Local:
    integer :: i,j
    real(r8) :: dtr,calday
    real(r8) :: cosZenAngR(1)

    dtr = pi/180._r8
    calday = get_curr_calday() ! fractional day of year
    do j=j0,j1
      do i=i0,i1
        call zenith(calday,(/dtr*glat(j)/),(/dtr*glon(i)/),cosZenAngR,1)
        chi(i,j) = acos(cosZenAngR(1))
      enddo 
    enddo
  end subroutine get_zenith
!-----------------------------------------------------------------------
  subroutine divb(dvb,i0,i1,j0,j1)
!
! Evaluate divergence of B, the unit magnetic field vector.
! (all processors have the full global 2d field)
!
! Args:
    integer,intent(in) :: i0,i1,j0,j1
    real(r8),intent(out) :: dvb(i0:i1,j0:j1) 
!
! Local:
    integer :: i,j,jm1,jp1
    real(r8),parameter :: re = 6.37122e8_r8  ! earth radius (cm)

    dvb = 0._r8

    call savefld_waccm(bx(i0:i1,j0:j1),'OPLUS_BX',1,i0,i1,j0,j1)
    call savefld_waccm(by(i0:i1,j0:j1),'OPLUS_BY',1,i0,i1,j0,j1)
    call savefld_waccm(bz(i0:i1,j0:j1),'OPLUS_BZ',1,i0,i1,j0,j1)
    call savefld_waccm(bmod2(i0:i1,j0:j1),'OPLUS_BMAG',1,i0,i1,j0,j1)
!
! Note re is in cm.
! (bx,by,bz are set by sub magfield (getapex.F90))
! (dphi,dlamda, and cs are set by sub set_geogrid (edyn_init.F90))
!
    do j=j0,j1
      jm1 = j-1
      jp1 = j+1
      do i=i0,i1
        dvb(i,j) = (((bx(i+1,j)-bx(i-1,j))/(2._r8*dlamda)+      &
          (cs(jp1)*by(i,jp1)-cs(jm1)*by(i,jm1))/(2._r8*dphi))/  &
          cs(j)+2._r8*bz(i,j))/re
      enddo ! i=i0,i1
    enddo ! j=j0,j1
  end subroutine divb
!-----------------------------------------------------------------------
  subroutine rrk(t,rms,ps1,ps2,n2,tr,ans,lon0,lon1,lev0,lev1)
!
! Returns ambipolar diffusion coefficient in ans.
!
! Args:
    integer,intent(in) :: lon0,lon1,lev0,lev1
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(in) :: &
      t,rms,ps1,ps2,n2,tr
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(out) :: ans
!
! Local:
!
    integer :: k,i
!
!$omp parallel do private(i,k)
    do i=lon0,lon1
      do k=lev0,lev1-1
 
        ans(k,i) = 1.42e17_r8*boltz*t(k,i)/(p0*expz(k)*.5_r8*(rms(k,i)+    &
          rms(k+1,i))*(ps2(k,i)*rmassinv_o1*sqrt(tr(k,i))*(1._r8-0.064_r8* &
          log10(tr(k,i)))**2*colfac+18.6_r8*n2(k,i)*rmassinv_n2+18.1_r8*   &
          ps1(k,i)*rmassinv_o2))

      enddo ! k=lev0,lev1
      ans(lev1,i) = ans(lev1-1,i) ! should not need to do this

    enddo ! i=lon0,lon1
!     
! Cap ambipolar diffusion coefficient in ans.
! 
    ! acceptable range for limiter 1.e8 to 1.e9 ...
    where( ans(:,:) > adiff_limiter )
      ans(:,:) = adiff_limiter
    endwhere

  end subroutine rrk
!-----------------------------------------------------------------------
  subroutine diffus(tp,en,hj,ans,i0,i1,lev0,lev1,lat)
!                                      kbot,nlev
! Evaluates ans = (d/(h*dz)*tp+m*g/r)*en
! Remember: "bot2top": lev0=kbot=bottom, lev1=nlev=top
!  
! Args:
    integer :: i0,i1,lev0,lev1,lat
    real(r8),dimension(lev0:lev1,i0:i1),intent(in) :: tp,en,hj
    real(r8),dimension(lev0:lev1,i0:i1),intent(out) :: ans
!
! Local:
    integer :: i,k
    real(r8) :: mgr

    mgr = rmass_op*grav_cm/gask

!$omp parallel do private(i,k)
    do i=i0,i1
      do k=lev0,lev1-2
        ans(k+1,i) = 1._r8/(2._r8*hj(k+1,i)*dzp)*(tp(k+2,i)*en(k+2,i)- &
          tp(k,i)*en(k,i))+mgr*en(k+1,i)
      enddo
      if (debug) then
        write(iulog,"('diffus: lat=',i4,' i=',i4,' ans(lev0:lev1-1,i)=',2es12.4)") &
          lat,i,minval(ans(lev0:lev1-1,i)),maxval(ans(lev0:lev1-1,i))
      endif
    enddo
!
! Upper and lower boundaries:
!
!$omp parallel do private(i)
    do i=i0,i1
!
! Upper boundary:
      ans(lev1,i) = 1._r8/(hj(lev1,i)*dzp)*(tp(lev1,i)*en(lev1,i)- &
        tp(lev1-1,i)*en(lev1-1,i))+mgr*en(lev1,i)
!
! Lower boundary:
      ans(lev0,i) = 1._r8/(hj(lev0,i)*dzp)*(tp(lev0+1,i)*en(lev0+1,i)- &
        tp(lev0,i)*en(lev0,i))+mgr*en(lev0,i)
    enddo
  end subroutine diffus
!-----------------------------------------------------------------------
  subroutine bdotdh(phijm1,phij,phijp1,ans,lon0,lon1,lev0,lev1,lat)
!
! Evaluates ans = (b(h)*del(h))*phi
!
! Args:
    integer,intent(in) :: lon0,lon1,lev0,lev1,lat
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(in) :: phijm1,phijp1
    real(r8),dimension(lev0:lev1,lon0-2:lon1+2),intent(inout) :: phij ! why intent(inout)?
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(out) :: ans
!
! Local:
    integer :: k,i
!
! Note phij longitude dimension is lon0-2:lon1+2 (only i-1 and i+1 are used).
! Halo longitudes i-1 and i+1 must have been set before this routine is
! called. ('by' is use-associated above)
!
!$omp parallel do private( i, k )
    do i=lon0,lon1
      do k=lev0,lev1
        ans(k,i) = 1._r8/re*(bx(i,lat)/(cs(lat)*2._r8*dlamda)* &
          (phij(k,i+1)-phij(k,i-1))+by(i,lat)*                 &
          (phijp1(k,i)-phijm1(k,i))/(2._r8*dphi))
      enddo ! k=lev0,lev1
    enddo ! i=lon0,lon1
!
  end subroutine bdotdh
!-----------------------------------------------------------------------
  subroutine bdzdvb(phi,dvb,h,ans,lev0,lev1,lon0,lon1,lat)
!
! Evaluates  ans = (bz*d/(h*dz)+divb)*phi
!
! Args:
    integer,intent(in) :: lev0,lev1,lon0,lon1,lat
    real(r8),intent(in) :: dvb(lon0:lon1)
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(in)    :: phi,h
    real(r8),dimension(lev0:lev1,lon0:lon1),intent(out)   :: ans
!
! Local:
    integer :: k,i
!
!$omp parallel do private( i, k )
    do i=lon0,lon1
      do k=lev0+1,lev1-1
        ans(k,i) = bz(i,lat)/(2._r8*h(k,i)*dzp)*(phi(k+1,i)-phi(k-1,i))+ &
          dvb(i)*phi(k,i)
      enddo ! k=lev0+1,lev1-1
    enddo ! i=lon0,lon1
!
! Upper and lower boundaries:
!$omp parallel do private( i )
    do i=lon0,lon1
      ans(lev1,i) = bz(i,lat)/(h(lev1,i)*dzp)*(phi(lev1,i)- &
        phi(lev1-1,i))+dvb(i)*phi(lev1,i)
      ans(lev0,i) = bz(i,lat)/(h(lev0,i)*dzp)* &
        (phi(lev0+1,i)-phi(lev0,i))+dvb(i)*phi(lev0,i)
    enddo ! i=lon0,lon1
  end subroutine bdzdvb
!-----------------------------------------------------------------------
  subroutine printpoles(f,klev,k0,k1,i0,i1,j0,j1,name)
    use edyn_geogrid,only : nlat
!
! Args:
    integer,intent(in) :: klev,k0,k1,i0,i1,j0,j1
    real(r8),intent(in) :: f(k0:k1,i0:i1,j0:j1)
    character(len=*),intent(in) :: name
!
! Print values at the poles at klev:
    if (j0==1) then
      if (debug.and.masterproc) write(iulog,"(/,'printpoles ',a,' spole: klev=',i4,' f(klev,i0:i1,j0)=',/,(8es12.4))") &
        name,klev,f(klev,i0:i1,j0)
    endif
    if (j1==nlat) then
      if (debug.and.masterproc) write(iulog,"(/,'printpoles ',a,' npole: klev=',i4,' f(klev,i0:i1,j1)=',/,(8es12.4))") &
        name,klev,f(klev,i0:i1,j1)
    endif

  end subroutine printpoles
!-----------------------------------------------------------------------
  subroutine filter1_op(f,k0,k1,i0,i1,j0,j1)
!
! Polar fft filter, option 1 (see filter.F90).
!
    use filter_module,only: filter1,kut1
    use edyn_mpi     ,only: mp_gatherlons_f3d,mytidi
    use edyn_mpi     ,only: mp_scatterlons_f3d
    use edyn_geogrid ,only: nlon
!
! Args:
    integer,intent(in) :: k0,k1,i0,i1,j0,j1
    real(r8),intent(inout) :: f(k0:k1,i0:i1,j0:j1)
!
! Local:
    integer :: i,j,k,nlevs
    real(r8) :: fik(nlon,k1-k0+1)
    type(array_ptr_type) :: fkij(1) ! fkij(1)%ptr(k1-k0+1,nlon,j0:j1)

    nlevs = k1-k0+1
!
! Define lons in fkij from current task subdomain:
!
    allocate(fkij(1)%ptr(nlevs,nlon,j0:j1))
    do j=j0,j1
      do i=i0,i1
        do k=k0,k1 ! kbot,nlev
          fkij(1)%ptr(k-k0+1,i,j) = f(k,i,j)
        enddo
      enddo
    enddo
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0
!   gather lons from other tasks in that row). This includes all latitudes.
!
    call mp_gatherlons_f3d(fkij,1,nlevs,i0,i1,j0,j1,1)
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
!
    if (mytidi==0) then
!
! Define 2d array with all longitudes for filter at each latitude:
!
      latscan: do j=j0,j1
        if (kut1(j) >= nlon/2) cycle latscan
        do i=1,nlon
          do k=k0,k1
            fik(i,k-k0+1) = fkij(1)%ptr(k-k0+1,i,j)
          enddo
        enddo 
!
! Remove wave numbers > kut(lat):
!
        call filter1(fik,1,nlevs,j)
!
! Return filtered array to fkij:
!
        do i=1,nlon
          do k=k0,k1
            fkij(1)%ptr(k-k0+1,i,j) = fik(i,k-k0+1)
          enddo
        enddo ! i=1,nlon
      enddo latscan ! j=j0,j1
    endif ! mytidi==0
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
    call mp_scatterlons_f3d(fkij,1,nlevs,i0,i1,j0,j1,1)
!
! Return filtered array to inout field at task subdomain:
!
    do j=j0,j1
      do i=i0,i1
        do k=k0,k1
          f(k,i,j) = fkij(1)%ptr(k-k0+1,i,j)
        enddo
      enddo
    enddo
    deallocate(fkij(1)%ptr)
  end subroutine filter1_op
!-----------------------------------------------------------------------
  subroutine filter2_op(f,k0,k1,i0,i1,j0,j1)
    use filter_module,only: filter2
    use edyn_mpi     ,only: mp_gatherlons_f3d,mytidi
    use edyn_mpi     ,only: mp_scatterlons_f3d
    use edyn_geogrid ,only: nlon
!
! Args:
    integer,intent(in) :: k0,k1,i0,i1,j0,j1
    real(r8),intent(inout) :: f(k0:k1,i0:i1,j0:j1)
!
! Local:
    integer :: i,j,k,nlevs
    real(r8) :: fik(nlon,k1-k0+1)
    type(array_ptr_type) :: fkij(1) ! fkij(1)%ptr(k1-k0+1,nlon,j0:j1)

    nlevs = k1-k0+1
!
! Define lons in fkij from current task subdomain:
!
    allocate(fkij(1)%ptr(nlevs,nlon,j0:j1))
!$omp parallel do private( i,j,k )
    do j=j0,j1
      do i=i0,i1
        do k=k0,k1
          fkij(1)%ptr(k-k0+1,i,j) = f(k,i,j)
        enddo
      enddo
    enddo
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0
!   gather lons from other tasks in that row). This includes all latitudes.
!
    call mp_gatherlons_f3d(fkij,1,nlevs,i0,i1,j0,j1,1)
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
!
    if (mytidi==0) then
!
! Define 2d array with all longitudes for filter at each latitude:
!
      do j=j0,j1
        do i=1,nlon
          do k=k0,k1
            fik(i,k-k0+1) = fkij(1)%ptr(k-k0+1,i,j)
          enddo
        enddo 
!
! Remove wave numbers > kut(lat):
!
        call filter2(fik,1,nlevs,j)
!
! Return filtered array to fkij:
!
        do i=1,nlon
          do k=k0,k1
            fkij(1)%ptr(k-k0+1,i,j) = fik(i,k-k0+1)
          enddo
        enddo ! i=1,nlon
      enddo ! j=j0,j1
    endif ! mytidi==0
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
    call mp_scatterlons_f3d(fkij,1,nlevs,i0,i1,j0,j1,1)
!
! Return filtered array to inout field at task subdomain:
    do j=j0,j1
      do i=i0,i1
        do k=k0,k1
          f(k,i,j) = fkij(1)%ptr(k-k0+1,i,j)
        enddo
      enddo
    enddo
    deallocate(fkij(1)%ptr)
  end subroutine filter2_op
!-----------------------------------------------------------------------
end module oplus
