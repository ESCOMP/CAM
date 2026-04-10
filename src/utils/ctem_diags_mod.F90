!-----------------------------------------------------------------------------
! For physics grid (and dycore) independent circulation diagnostics
! -- terms of the Transformed Eulerian Mean (TEM) equation
!
! This uses ESMF utilities to remap dynamical fields (U,V, etc) from the physics
! grid to a regular longitude / latitude grid where it is convenient to compute
! the zonal mean terms of the TEM equation
!-----------------------------------------------------------------------------
module ctem_diags_mod
  use shr_kind_mod, only: r8 => shr_kind_r8, cx => SHR_KIND_CX
  use ppgrid, only: begchunk, endchunk, pcols, pver, pverp
  use physics_types, only: physics_state
  use phys_grid, only: get_ncols_p
  use spmd_utils, only: masterproc
  use ref_pres, only: pref_mid
  use esmf_lonlat_grid_mod, only: beglon=>lon_beg, endlon=>lon_end, beglat=>lat_beg, endlat=>lat_end
  use cam_history,  only: addfld, outfld, horiz_only, add_default
  use cam_history_support, only : fillvalue
  use perf_mod, only: t_startf, t_stopf
  use cam_logfile, only: iulog
  use cam_abortutils, only: endrun

  implicit none

  private

  public :: ctem_diags_readnl
  public :: ctem_diags_reg
  public :: ctem_diags_init
  public :: ctem_diags_calc
  public :: ctem_diags_final

  integer :: ctem_diags_numlats = 0
  logical :: ctem_diags_active = .false.

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine ctem_diags_readnl(nlfile)
    use namelist_utils, only : find_group_name
    use spmd_utils, only : mpicom, masterprocid, mpi_integer, mpi_success
    use string_utils, only : to_lower

    character(len=*), intent(in) :: nlfile
    integer :: unitn, ierr
    character(len=cx) :: iomsg

    character(len=*), parameter :: prefix = 'ctem_diags_readnl: '

    namelist /ctem_diags_nl/ ctem_diags_numlats

    if (masterproc) then
       ! read namelist
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'ctem_diags_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, ctem_diags_nl, iostat=ierr, iomsg=iomsg)
          if (ierr /= 0) then
             call endrun(prefix//'ctem_diags_nl: ERROR reading namelist: '//trim(iomsg))
          end if
       end if
       close(unitn)
    end if

    call mpi_bcast(ctem_diags_numlats, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(prefix//'mpi_bcast error : ctem_diags_numlats')

    ctem_diags_active = ctem_diags_numlats > 0

    if (masterproc) then
       write(iulog,*) prefix//'ctem_diags_numlats: ', ctem_diags_numlats
       write(iulog,*) prefix//'ctem_diags_active : ', ctem_diags_active
    end if

  end subroutine ctem_diags_readnl

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine ctem_diags_reg()
    use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap, cam_grid_register
    use esmf_lonlat_grid_mod, only: glats, nlat, glons, nlon
    use esmf_lonlat_grid_mod, only: esmf_lonlat_grid_init
    use esmf_phys_mesh_mod, only: esmf_phys_mesh_init
    use esmf_phys2lonlat_mod, only: esmf_phys2lonlat_init

    integer, parameter :: zm_decomp  = 331 ! Must be unique within CAM
    integer, parameter :: reg_decomp = 332

    type(horiz_coord_t), pointer :: zmlon_coord
    type(horiz_coord_t), pointer :: zmlat_coord

    integer(iMap),       pointer :: grid_map(:,:)
    real(r8) :: zmlons(1)

    integer(iMap),       pointer :: coord_map(:) => null()
    type(horiz_coord_t), pointer :: lon_coord
    type(horiz_coord_t), pointer :: lat_coord
    integer :: i, j, ind, astat

    character(len=*), parameter :: subname = 'ctem_diags_reg: '

    if (.not.ctem_diags_active) return

    ! initialize grids and mapping
    call esmf_lonlat_grid_init(ctem_diags_numlats)
    call esmf_phys_mesh_init()
    call esmf_phys2lonlat_init()

    ! Zonal mean grid for history fields
    zmlons = 0._r8

    zmlat_coord => horiz_coord_create('zmlat', '', nlat, 'latitude', 'degrees_north', 1, nlat, glats)
    zmlon_coord => horiz_coord_create('zmlon', '', 1,    'longitude', 'degrees_east', 1, 1,   zmlons)

    ! grid decomposition map
    allocate(grid_map(4,endlat-beglat+1), stat=astat)
    if (astat/=0) then
       call endrun(subname//'not able to allocate grid_map array')
    end if

    ind = 0
    do j = beglat, endlat
       ind = ind + 1
       grid_map(1,ind) = 1
       grid_map(2,ind) = j
       if (beglon==1) then
          grid_map(3,ind) = 1
          grid_map(4,ind) = j
       else
          grid_map(3,ind) = 0
          grid_map(4,ind) = 0
       end if
    end do

    ! register the zonal average grid
    call cam_grid_register('ctem_zm', zm_decomp, zmlat_coord, zmlon_coord, grid_map, &
                           unstruct=.false., zonal_grid=.true.)

    nullify(grid_map)

    ! for the lon-lat grid
    allocate(grid_map(4, ((endlon - beglon + 1) * (endlat - beglat + 1))), stat=astat)
    if (astat/=0) then
       call endrun(subname//'not able to allocate grid_map array')
    end if

    ind = 0
    do i = beglat, endlat
       do j = beglon, endlon
          ind = ind + 1
          grid_map(1, ind) = j
          grid_map(2, ind) = i
          grid_map(3, ind) = j
          grid_map(4, ind) = i
       end do
    end do

    allocate(coord_map(endlat - beglat + 1), stat=astat)
    if (astat/=0) then
       call endrun(subname//'not able to allocate coord_map array')
    end if

    if (beglon==1) then
       coord_map = (/ (i, i = beglat, endlat) /)
    else
       coord_map = 0
    end if
    lat_coord => horiz_coord_create('reglat', '', nlat, 'latitude',  'degrees_north', beglat, endlat, &
                                    glats(beglat:endlat),  map=coord_map)

    nullify(coord_map)

    allocate(coord_map(endlon - beglon + 1), stat=astat)
    if (astat/=0) then
       call endrun(subname//'not able to allocate coord_map array')
    end if

    if (beglat==1) then
       coord_map = (/ (i, i = beglon, endlon) /)
    else
       coord_map = 0
    end if

    lon_coord => horiz_coord_create('reglon', '', nlon, 'longitude',  'degrees_east', beglon, endlon, &
                                    glons(beglon:endlon),  map=coord_map)

    nullify(coord_map)

    call cam_grid_register('ctem_lonlat', reg_decomp, lat_coord, lon_coord, grid_map, unstruct=.false.)

    nullify(grid_map)

  end subroutine ctem_diags_reg

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine ctem_diags_init()

    if (.not.ctem_diags_active) return

    ! fields on reg lon lat grid
    call addfld ('THtem', (/'lev'/), 'A','K',      'Potential temp', gridname='ctem_lonlat' )
    call addfld ('Utem',  (/'lev'/), 'A','m s-1',  'Zonal wind', gridname='ctem_lonlat' )
    call addfld ('Vtem',  (/'lev'/), 'A','m s-1',  'Meridional wind', gridname='ctem_lonlat' )
    call addfld ('Wtem',  (/'lev'/), 'A','m s-1',  'Vertical wind', gridname='ctem_lonlat' )
    call addfld ('VTHtem',(/'lev'/), 'A','K m s-1','Meridional Heat Flux:', gridname='ctem_lonlat')
    call addfld ('WTHtem',(/'lev'/), 'A','K m s-1','Vertical Heat Flux:', gridname='ctem_lonlat')
    call addfld ('UVtem', (/'lev'/), 'A','m2 s-2', 'Meridional Flux of Zonal Momentum', gridname='ctem_lonlat')
    call addfld ('UWtem', (/'lev'/), 'A','m2 s-2', 'Vertical Flux of Zonal Momentum', gridname='ctem_lonlat')

    ! fields on zonal mean grid
    call addfld ('Uzm',  (/'lev'/), 'A','m s-1',  'Zonal-Mean zonal wind', gridname='ctem_zm' )
    call addfld ('Vzm',  (/'lev'/), 'A','m s-1',  'Zonal-Mean meridional wind', gridname='ctem_zm' )
    call addfld ('Wzm',  (/'lev'/), 'A','m s-1',  'Zonal-Mean vertical wind', gridname='ctem_zm' )
    call addfld ('THzm', (/'lev'/), 'A','K',      'Zonal-Mean potential temp', gridname='ctem_zm' )
    call addfld ('VTHzm',(/'lev'/), 'A','K m s-1','Zonal-Mean meridional Heat Flux:', gridname='ctem_zm')
    call addfld ('WTHzm',(/'lev'/), 'A','K m s-1','Zonal-Mean vertical Heat Flux:', gridname='ctem_zm')
    call addfld ('UVzm', (/'lev'/), 'A','m2 s-2', 'Zonal-Mean meridional Flux of Zonal Momentum', gridname='ctem_zm')
    call addfld ('UWzm', (/'lev'/), 'A','m2 s-2', 'Zonal-Mean vertical Flux of Zonal Momentum', gridname='ctem_zm')

    call addfld ('PSzm',  horiz_only, 'A', 'Pa', 'Zonal-Mean surface pressure', gridname='ctem_zm' )
    call addfld ('PStem', horiz_only, 'A', 'Pa', 'Surface Pressure', gridname='ctem_lonlat')
    call addfld ('MSKtem',horiz_only, 'A', '1',  'TEM mask', gridname='ctem_lonlat' )

  end subroutine ctem_diags_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine ctem_diags_calc(phys_state)
    use air_composition, only: mbarv ! g/mole
    use shr_const_mod, only: rgas => shr_const_rgas ! J/K/kmole
    use shr_const_mod, only: grav => shr_const_g ! m/s2
    use esmf_phys2lonlat_mod, only: esmf_phys2lonlat_regrid
    use esmf_zonal_mean_mod, only: esmf_zonal_mean_calc, esmf_zonal_mean_wsums, esmf_zonal_mean_masked
    use interpolate_data, only: lininterp
    use esmf_phys2lonlat_mod, only: fields_bundle_t, nflds

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    real(r8), target :: u_phys(pver,pcols,begchunk:endchunk)
    real(r8), target :: v_phys(pver,pcols,begchunk:endchunk)
    real(r8), target :: w_phys(pver,pcols,begchunk:endchunk)
    real(r8), target :: t_phys(pver,pcols,begchunk:endchunk)

    real(r8), target :: uv_phys(pver,pcols,begchunk:endchunk)
    real(r8), target :: uw_phys(pver,pcols,begchunk:endchunk)
    real(r8), target :: vt_phys(pver,pcols,begchunk:endchunk)
    real(r8), target :: wt_phys(pver,pcols,begchunk:endchunk)

    real(r8) :: ps_phys(pcols,begchunk:endchunk)

    real(r8), target :: u_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8), target :: v_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8), target :: w_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8), target :: t_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8) :: ps_lonlat(beglon:endlon,beglat:endlat)
    real(r8) :: mskind1(beglon:endlon,beglat:endlat) ! vertical index where mountain masking begins

    real(r8), target  :: uv_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8), target  :: uw_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8), target  :: vt_lonlat(beglon:endlon,beglat:endlat,pver)
    real(r8), target  :: wt_lonlat(beglon:endlon,beglat:endlat,pver)

    real(r8) :: u_zm(beglat:endlat,pver)
    real(r8) :: v_zm(beglat:endlat,pver)
    real(r8) :: w_zm(beglat:endlat,pver)
    real(r8) :: t_zm(beglat:endlat,pver)

    real(r8) :: uv_zm(beglat:endlat,pver)
    real(r8) :: uw_zm(beglat:endlat,pver)
    real(r8) :: vt_zm(beglat:endlat,pver)
    real(r8) :: wt_zm(beglat:endlat,pver)

    real(r8) :: ps_zm(beglat:endlat)

    real(r8) :: vtp(beglon:endlon,beglat:endlat,pver)
    real(r8) :: wtp(beglon:endlon,beglat:endlat,pver)
    real(r8) :: uwp(beglon:endlon,beglat:endlat,pver)
    real(r8) :: uvp(beglon:endlon,beglat:endlat,pver)

    real(r8) :: vtp_zm(beglat:endlat,pver)
    real(r8) :: wtp_zm(beglat:endlat,pver)
    real(r8) :: uwp_zm(beglat:endlat,pver)
    real(r8) :: uvp_zm(beglat:endlat,pver)
    real(r8) :: wsums(beglat:endlat,pver)

    integer  :: lchnk, ncol, i, j, k

    real(r8) :: outtmp(beglon:endlon,pver)
    integer :: outcnt

    real(r8) :: wght(beglon:endlon,beglat:endlat,pver)

    type(fields_bundle_t) :: physflds(nflds)
    type(fields_bundle_t) :: lonlatflds(nflds)

    real(r8) :: theta(pver)
    real(r8), parameter :: hscale = 7000._r8 ! pressure scale height (meters)

    if (.not.ctem_diags_active) return

    call t_startf('ctem_diags_calc')

    call t_startf('ctem_diags_calc-setarrs')

    do lchnk = begchunk,endchunk
       ncol = phys_state(lchnk)%ncol
       do i = 1,ncol

          ! wind components -- vertically intepolate to ref press
          call lininterp( phys_state(lchnk)%u(i,:), phys_state(lchnk)%pmid(i,:), pver, &
                          u_phys(:,i,lchnk), pref_mid(:), pver )

          call lininterp( phys_state(lchnk)%v(i,:), phys_state(lchnk)%pmid(i,:), pver, &
                          v_phys(:,i,lchnk), pref_mid(:), pver )

          call lininterp( phys_state(lchnk)%omega(i,:), phys_state(lchnk)%pmid(i,:), pver, &
                          w_phys(:,i,lchnk), pref_mid(:), pver )

          ! omega -> vertical velocity
          w_phys(:,i,lchnk) = -hscale * w_phys(:,i,lchnk)/pref_mid(:)

          ! potential temperature -- vertically intepolate to ref press
          theta(:) =  phys_state(lchnk)%t(i,:) * phys_state(lchnk)%exner(i,:)
          call lininterp( theta(:), phys_state(lchnk)%pmid(i,:), pver, &
                          t_phys(:,i,lchnk), pref_mid(:), pver )

          ! surface pressure
          ps_phys(i,lchnk) = phys_state(lchnk)%ps(i)

          uv_phys(:,i,lchnk) = u_phys(:,i,lchnk)*v_phys(:,i,lchnk)
          uw_phys(:,i,lchnk) = u_phys(:,i,lchnk)*w_phys(:,i,lchnk)
          vt_phys(:,i,lchnk) = v_phys(:,i,lchnk)*t_phys(:,i,lchnk)
          wt_phys(:,i,lchnk) = w_phys(:,i,lchnk)*t_phys(:,i,lchnk)

       end do
    end do

    call t_stopf('ctem_diags_calc-setarrs')

    call t_startf('ctem_diags_calc-regrid')

    ! regrid to lon/lat grid

    ! set feild-bundle pointers for regridding utility
    physflds(1)%fld => u_phys ! regrid inputs
    physflds(2)%fld => v_phys
    physflds(3)%fld => w_phys
    physflds(4)%fld => t_phys
    physflds(5)%fld => uv_phys
    physflds(6)%fld => uw_phys
    physflds(7)%fld => vt_phys
    physflds(8)%fld => wt_phys

    lonlatflds(1)%fld => u_lonlat ! regrid outputs
    lonlatflds(2)%fld => v_lonlat
    lonlatflds(3)%fld => w_lonlat
    lonlatflds(4)%fld => t_lonlat
    lonlatflds(5)%fld => uv_lonlat
    lonlatflds(6)%fld => uw_lonlat
    lonlatflds(7)%fld => vt_lonlat
    lonlatflds(8)%fld => wt_lonlat

    ! regrid 3-D fields
    call esmf_phys2lonlat_regrid(physflds, lonlatflds)

    ! regrid 2-D field separately
    call esmf_phys2lonlat_regrid(ps_phys, ps_lonlat)

    call t_stopf('ctem_diags_calc-regrid')


    call t_startf('ctem_diags_calc-zonal_mean-uvwt')

    ! Mask out the mountains
    do i = beglon,endlon
       do j = beglat,endlat

          mskind1(i,j) = pver

          do k = 1,pver
             if (ps_lonlat(i,j)>=pref_mid(k)) then
                wght(i,j,k) = 1._r8
                mskind1(i,j) = k ! vertical index where masking begins
             else
                wght(i,j,k) = 0._r8
             end if
          end do
       end do
    end do

    ! calculate zonal means from interpolated fields
    ! mask out mountains from the zonal mean calculations
    wsums = esmf_zonal_mean_wsums(wght)

    ! compute zonal-mean fields
    call esmf_zonal_mean_masked(u_lonlat, wght, wsums, u_zm)
    call esmf_zonal_mean_masked(v_lonlat, wght, wsums, v_zm)
    call esmf_zonal_mean_masked(w_lonlat, wght, wsums, w_zm)
    call esmf_zonal_mean_masked(t_lonlat, wght, wsums, t_zm)
    call esmf_zonal_mean_masked(uv_lonlat, wght, wsums, uv_zm)
    call esmf_zonal_mean_masked(uw_lonlat, wght, wsums, uw_zm)
    call esmf_zonal_mean_masked(vt_lonlat, wght, wsums, vt_zm)
    call esmf_zonal_mean_masked(wt_lonlat, wght, wsums, wt_zm)

    call t_stopf('ctem_diags_calc-zonal_mean-uvwt')

    ! compute zonal-mean PS
    call esmf_zonal_mean_calc(ps_lonlat, ps_zm)

    call t_startf('ctem_diags_calc-calc_dev_flx')

    ! Calculate fluxes
    do k = 1,pver
       do j = beglat,endlat
          do i = beglon,endlon
             if (wght(i,j,k)>0._r8) then

                ! u'v' = (uv)zm - uzm*vzm

                uvp(i,j,k) = uv_zm(j,k) - (u_zm(j,k)*v_zm(j,k))
                uwp(i,j,k) = uw_zm(j,k) - (u_zm(j,k)*w_zm(j,k))
                vtp(i,j,k) = vt_zm(j,k) - (v_zm(j,k)*t_zm(j,k))
                wtp(i,j,k) = wt_zm(j,k) - (w_zm(j,k)*t_zm(j,k))

             else
                uvp(i,j,k) = fillvalue
                uwp(i,j,k) = fillvalue
                vtp(i,j,k) = fillvalue
                wtp(i,j,k) = fillvalue
             end if
          end do
       end do
    end do

    call t_stopf('ctem_diags_calc-calc_dev_flx')

    call t_startf('ctem_diags_calc-zonal_mean-p')

    ! compute zonal-mean flux terms
    call esmf_zonal_mean_masked(vtp, wght, wsums, vtp_zm)
    call esmf_zonal_mean_masked(wtp, wght, wsums, wtp_zm)
    call esmf_zonal_mean_masked(uwp, wght, wsums, uwp_zm)
    call esmf_zonal_mean_masked(uvp, wght, wsums, uvp_zm)

    call t_stopf('ctem_diags_calc-zonal_mean-p')

    call t_startf('ctem_diags_calc-output')

    outcnt = endlon-beglon+1

    ! output TEM diagnostics
    do j = beglat,endlat
       outtmp(beglon:endlon,1:pver) = t_lonlat(beglon:endlon,j,1:pver)
       call outfld('THtem',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = u_lonlat(beglon:endlon,j,1:pver)
       call outfld('Utem',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = v_lonlat(beglon:endlon,j,1:pver)
       call outfld('Vtem',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = w_lonlat(beglon:endlon,j,1:pver)
       call outfld('Wtem',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = vtp(beglon:endlon,j,1:pver)
       call outfld('VTHtem',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = wtp(beglon:endlon,j,1:pver)
       call outfld('WTHtem',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = uvp(beglon:endlon,j,1:pver)
       call outfld('UVtem',outtmp, outcnt, j)
       outtmp(beglon:endlon,1:pver) = uwp(beglon:endlon,j,1:pver)
       call outfld('UWtem',outtmp, outcnt, j)

       call outfld('PStem', ps_lonlat(:,j), outcnt, j)
       call outfld('MSKtem',mskind1(:,j), outcnt, j)

       call outfld('Uzm',  u_zm(j,:), 1,j)
       call outfld('Vzm',  v_zm(j,:), 1,j)
       call outfld('Wzm',  w_zm(j,:), 1,j)
       call outfld('THzm', t_zm(j,:), 1,j)
       call outfld('PSzm', ps_zm(j), 1,j)

       call outfld('VTHzm',vtp_zm(j,:),1,j)
       call outfld('WTHzm',wtp_zm(j,:),1,j)
       call outfld('UVzm', uvp_zm(j,:),1,j)
       call outfld('UWzm', uwp_zm(j,:),1,j)
    end do

    call t_stopf('ctem_diags_calc-output')

    call t_stopf('ctem_diags_calc')

  end subroutine ctem_diags_calc

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine ctem_diags_final()
    use esmf_phys2lonlat_mod, only: esmf_phys2lonlat_destroy
    use esmf_lonlat_grid_mod, only: esmf_lonlat_grid_destroy
    use esmf_phys_mesh_mod, only: esmf_phys_mesh_destroy

    if (.not.ctem_diags_active) return

    call esmf_phys2lonlat_destroy()
    call esmf_lonlat_grid_destroy()
    call esmf_phys_mesh_destroy()

  end subroutine ctem_diags_final

end module ctem_diags_mod
