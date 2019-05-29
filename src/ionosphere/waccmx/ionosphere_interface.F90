module ionosphere_interface

  use shr_kind_mod,        only: r8 => shr_kind_r8
  use phys_grid,           only: begchunk, endchunk, get_ncols_p
  use pmgrid,              only: plat, plon, plev
  use ppgrid,              only: pcols, pver

  use dpie_coupling,       only: d_pie_init
  use dpie_coupling,       only: d_pie_epotent
  use dpie_coupling,       only: d_pie_coupling         ! WACCM-X ionosphere/electrodynamics coupling
  use short_lived_species, only: slvd_index,slvd_pbf_ndx => pbf_idx ! Routines to access short lived species 

  use chem_mods,           only: adv_mass      ! Array holding mass values for short lived species
  use mo_chem_utls,        only: get_spc_ndx   ! Routine to get index of adv_mass array for short lived species
  use physics_buffer,      only: pbuf_get_chunk, pbuf_get_field, pbuf_get_index

  use cam_abortutils,      only: endrun
  use constituents,        only: cnst_get_ind, cnst_mw  !Needed to access constituent molecular weights
  use phys_grid,           only: get_lon_all_p, get_lat_all_p, transpose_block_to_chunk, transpose_chunk_to_block
  use phys_grid,           only: chunk_to_block_send_pters, chunk_to_block_recv_pters, block_to_chunk_send_pters, &
                                 block_to_chunk_recv_pters
  use physconst,           only: gravit
  use oplus,               only: oplus_init
  use edyn_init,           only: edynamo_init
  use pio,                 only: var_desc_t
  use spmd_dyn,            only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
  use dyn_internal_state,  only: get_dyn_state_grid
  use dynamics_vars,       only: t_fvdycore_grid
  use perf_mod
  use epotential_params,   only: epot_active, epot_crit_colats
  implicit none

  private

  public :: ionosphere_readnl
  public :: ionosphere_init
  public :: ionosphere_run1
  public :: ionosphere_run2
  public :: ionosphere_init_restart
  public :: ionosphere_write_restart
  public :: ionosphere_read_restart
  public :: ionosphere_final

  ! private data 

  ! this needs to persist from time-step to time-step and across restarts
  real(r8), allocatable :: opmmrtm1_blck(:,:,:)   ! O+ at previous time step(blocks)

  type(var_desc_t) :: Optm1_vdesc 
  integer :: index_ped, index_hall, index_te, index_ti
  integer :: index_ui, index_vi, index_wi

  integer :: ixo2=-1, ixo=-1, ixh=-1
  integer :: ixo2p=-1, ixnop=-1, ixn2p=-1, ixop=-1

  ! indices for accessing ions in pbuf when non-advected
  integer :: sIndxOp=-1, sIndxO2p=-1, sIndxNOp=-1, sIndxN2p=-1  

  real(r8) :: rmassO2    ! O2 molecular weight kg/kmol
  real(r8) :: rmassO1    ! O atomic weight kg/kmol
  real(r8) :: rmassH     ! H atomic weight kg/kmol
  real(r8) :: rmassN2    ! N2 molecular weight kg/kmol
  real(r8) :: rmassO2p   ! O2+ molecular weight kg/kmol
  real(r8) :: rmassNOp   ! NO+ molecular weight kg/kmol
  real(r8) :: rmassN2p   ! N2+ molecular weight kg/kmol
  real(r8) :: rmassOp    ! O+ molecular weight kg/kmol

  logical, public,  protected :: ionos_edyn_active = .true.   ! if true, edynamo will generate ion drifts
  logical, public,  protected :: ionos_xport_active = .true.  ! if true, call d_pie_coupling from dp_coupling
  !
  ! ionos_edyn_active = .true. will activate the edynamo which will generate ion drift velocities 
  !  used in oplus transport, otherwise empirical ion drifts calculated in exbdrift (physics) will be used.
  !
  logical, public,  protected :: ionos_oplus_xport = .true.    ! if true, call sub oplus (based on tiegcm oplus.F)
  integer, public,  protected :: ionos_xport_nsplit = 5        ! number of substeps for O+ transport per model time step

  real(r8), public, protected :: oplus_adiff_limiter = 1.5e+8_r8  ! limiter for ambipolar diffusion coefficient
  real(r8), public, protected :: oplus_shapiro_const = 0.03_r8    ! shapiro constant for spatial smoother
  logical,  public, protected :: oplus_enforce_floor = .true.     ! switch to apply Stan's  floor
  logical,  public, protected :: oplus_ring_polar_filter = .false. ! switch to apply ring polar filter

  character(len=256) :: wei05_coefs_file = 'NONE' !'wei05sc.nc'
  character(len=256) :: amienh_file  = 'NONE'
  character(len=256) :: amiesh_file  = 'NONE'

  character(len=16), public, protected :: ionos_epotential_model = 'none'
  logical,           public, protected :: ionos_epotential_amie = .false.
  integer ::  indxAMIEefxg=-1, indxAMIEkevg=-1

contains

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_readnl( nlfile )

    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use spmd_utils,     only: mpicom, masterprocid, mpi_real8, mpi_logical, mpi_integer, mpi_character
    use cam_logfile,    only: iulog
    use spmd_utils,     only: masterproc

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'ionosphere_readnl'

    namelist /ionosphere_nl/ ionos_xport_active, ionos_edyn_active, ionos_oplus_xport, ionos_xport_nsplit
    namelist /ionosphere_nl/ oplus_adiff_limiter, oplus_shapiro_const, oplus_enforce_floor, oplus_ring_polar_filter
    namelist /ionosphere_nl/ ionos_epotential_model, ionos_epotential_amie, wei05_coefs_file
    namelist /ionosphere_nl/ amienh_file, amiesh_file, wei05_coefs_file
    namelist /ionosphere_nl/ epot_crit_colats
    
    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'ionosphere_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, ionosphere_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(ionos_xport_active,  1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(ionos_edyn_active,   1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(ionos_oplus_xport,   1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(ionos_xport_nsplit,  1, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(oplus_adiff_limiter, 1, mpi_real8,   masterprocid, mpicom, ierr)
    call mpi_bcast(ionos_epotential_model, len(ionos_epotential_model), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(ionos_epotential_amie,1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(wei05_coefs_file, len(wei05_coefs_file), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(amienh_file, len(amienh_file), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(amiesh_file, len(amiesh_file), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(oplus_shapiro_const, 1, mpi_real8,   masterprocid, mpicom, ierr)
    call mpi_bcast(oplus_enforce_floor, 1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(oplus_ring_polar_filter,1, mpi_logical, masterprocid, mpicom, ierr)
    call mpi_bcast(epot_crit_colats,    2, mpi_real8,   masterprocid, mpicom, ierr)
    
    ! log the user settings
    if (masterproc) then
       write(iulog,*) 'ionosphere_readnl: ionos_xport_active     = ', ionos_xport_active
       write(iulog,*) 'ionosphere_readnl: ionos_edyn_active      = ', ionos_edyn_active
       write(iulog,*) 'ionosphere_readnl: ionos_oplus_xport      = ', ionos_oplus_xport
       write(iulog,*) 'ionosphere_readnl: ionos_xport_nsplit     = ', ionos_xport_nsplit
       write(iulog,*) 'ionosphere_readnl: ionos_epotential_model = ', trim(ionos_epotential_model)
       write(iulog,*) 'ionosphere_readnl: ionos_epotential_amie  = ', ionos_epotential_amie
       write(iulog,'(a,2(g12.4))') &
                     ' ionosphere_readnl: epot_crit_colats       = ', epot_crit_colats
       write(iulog,*) 'ionosphere_readnl: oplus_adiff_limiter    = ', oplus_adiff_limiter
       write(iulog,*) 'ionosphere_readnl: oplus_shapiro_const    = ', oplus_shapiro_const
       write(iulog,*) 'ionosphere_readnl: oplus_enforce_floor    = ', oplus_enforce_floor
       write(iulog,*) 'ionosphere_readnl: oplus_ring_polar_filter= ', oplus_ring_polar_filter
    endif
    epot_active = .true.
    
  end subroutine ionosphere_readnl

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_init()
    use physics_buffer, only: pbuf_add_field, dtype_r8
    use cam_history,    only: addfld, add_default, horiz_only
    use mo_apex,        only: mo_apex_init1
    use cam_control_mod,only: initial_run
    use dyn_grid,       only: get_horiz_grid_d
    use ref_pres,  only : & ! Hybrid level definitions:
      pref_mid,           & ! target alev(plev) midpoint levels coord
      pref_edge             ! target ailev(plevp) interface levels coord
    use amie_module,    only: init_amie
    use wei05sc,        only: weimer05_init
   
    ! local variables:
    type (t_fvdycore_grid), pointer :: grid
    integer :: sIndx

    integer :: mpicomm         ! MPI communicator
    integer :: ntaski, ntaskj  ! number of MPI tasks in lon,lat dimensions
    integer :: lat0,lat1       ! first and last latitude  indices
    integer :: lon0,lon1       ! first and last longitude indices
    integer :: lev0,lev1       ! first and last pressure indices
    real(r8), allocatable :: glon(:) ! global geo-graphic longitudes (degrees)
    real(r8), allocatable :: glat(:) ! global geo-graphic latitudes (degrees)

    if ( ionos_epotential_amie ) then
       call pbuf_add_field('AMIE_efxg', 'global', dtype_r8, (/pcols/), indxAMIEefxg)  ! Energy flux from AMIE
       call pbuf_add_field('AMIE_kevg', 'global', dtype_r8, (/pcols/), indxAMIEkevg)  ! Mean energy from AMIE  
    endif
    if (initial_run) then
       call ionosphere_read_ic()
    endif

    call mo_apex_init1()

    op_transport: if (ionos_xport_active) then

       grid => get_dyn_state_grid()

       index_ped  = pbuf_get_index('PedConduct')
       index_hall = pbuf_get_index('HallConduct')

       index_te   = pbuf_get_index('TElec')
       index_ti   = pbuf_get_index('TIon')
       !
       ! pbuf indices to empirical ion drifts, to be passed to oplus_xport,
       ! if ionos_edyn_active is false.
       !
       index_ui   = pbuf_get_index('UI')
       index_vi   = pbuf_get_index('VI')
       index_wi   = pbuf_get_index('WI')

       !-----------------------------------------------------------------------
       !  Get indices for neutrals to get mixing ratios from state%q and masses
       !-----------------------------------------------------------------------
       call cnst_get_ind('O2' ,ixo2 )
       call cnst_get_ind('O'  ,ixo )
       call cnst_get_ind('H'  ,ixh )
       !------------------------------------
       ! Get neutral molecular weights
       !------------------------------------
       rmassO2 = cnst_mw(ixo2)
       rmassO1 = cnst_mw(ixo)
       rmassH  = cnst_mw(ixh)
       rmassN2 = 28._r8

       call cnst_get_ind('Op',ixop, abort=.false.)
       if (ixop > 0) then
          rMassOp = cnst_mw(ixop)
       else
          sIndxOp  = slvd_index( 'Op' )
          if (sIndxOp > 0) then
             sIndx = get_spc_ndx( 'Op' )
             rmassOp = adv_mass(sIndx)
          else
             call endrun('ionosphere_init: Cannot find state or pbuf index for Op')
          endif
       endif

       call cnst_get_ind('O2p',ixo2p, abort=.false.)
       if (ixo2p > 0) then
          rMassO2p = cnst_mw(ixo2p)
       else
          sIndxO2p  = slvd_index( 'O2p' )
          if (sIndxO2p > 0) then
             sIndx = get_spc_ndx( 'O2p' )
             rmassO2p = adv_mass(sIndx)
          else
             call endrun('ionosphere_init: Cannot find state or pbuf index for O2p')
          endif
       endif

       call cnst_get_ind('NOp',ixnop, abort=.false.)
       if (ixnop > 0) then
          rMassNOp = cnst_mw(ixnop)
       else
          sIndxNOp  = slvd_index( 'NOp' )
          if (sIndxNOp > 0) then
             sIndx = get_spc_ndx( 'NOp' )
             rmassNOp = adv_mass(sIndx)
          else
             call endrun('ionosphere_init: Cannot find state or pbuf index for NOp')
          endif
       endif

       call cnst_get_ind('N2p',ixn2p, abort=.false.)
       if (ixn2p > 0) then
          rMassN2p = cnst_mw(ixn2p)
       else
          sIndxN2p  = slvd_index( 'N2p' )
          if (sIndxN2p > 0) then
             sIndx = get_spc_ndx( 'N2p' )
             rmassN2p = adv_mass(sIndx)
          else
             call endrun('ionosphere_init: Cannot find state or pbuf index for N2p')
          endif
       endif

       call d_pie_init( ionos_edyn_active, ionos_oplus_xport, ionos_xport_nsplit, epot_crit_colats )
 
       if ( grid%iam < grid%npes_xy ) then
          
          allocate(glon(plon))
          allocate(glat(plat))
          call get_horiz_grid_d( plon, lon_d_out=glon )
          call get_horiz_grid_d( plat, lat_d_out=glat )

          mpicomm = grid%commxy
          lon0 = grid%ifirstxy ; lon1 = grid%ilastxy
          lat0 = grid%jfirstxy ; lat1 = grid%jlastxy
          lev0 = 1             ; lev1 = grid%km
          ntaski = grid%nprxy_x
          ntaskj = grid%nprxy_y

          call edynamo_init( mpicomm, plon, plat, plev, lon0,lon1,lat0,lat1,lev0,lev1, ntaski,ntaskj, &
                             glon, glat, pref_mid,pref_edge )
          call ionosphere_alloc()
          call oplus_init( oplus_adiff_limiter, oplus_shapiro_const, oplus_enforce_floor, oplus_ring_polar_filter )

          deallocate(glon,glat)
       endif

       call addfld ('OpTM1&IC', (/ 'lev' /),'I','kg/kg','O+ at time step minus 1',gridname='fv_centers')
       call add_default ('OpTM1&IC',0, 'I')

    endif op_transport

    if (ionos_edyn_active) then
       call addfld ('UI',(/ 'lev' /),'I','m/s', 'UI Zonal ion drift from edynamo') 
       call addfld ('VI',(/ 'lev' /),'I','m/s', 'VI Meridional ion drift from edynamo')
       call addfld ('WI',(/ 'lev' /),'I','m/s', 'WI Vertical ion drift from edynamo')
       call addfld ('UI&IC', (/ 'lev' /), 'I','m/s', 'Zonal ion drift velocity')
       call addfld ('VI&IC', (/ 'lev' /), 'I','m/s', 'Meridional ion drift velocity')
       call addfld ('WI&IC', (/ 'lev' /), 'I','m/s', 'Vertical ion drift velocity')
       call add_default ('UI&IC', 0, ' ')
       call add_default ('VI&IC', 0, ' ')
       call add_default ('WI&IC', 0, ' ')
    endif
    if ( ionos_epotential_amie ) then
       call init_amie(amienh_file,amiesh_file)
       call addfld ('amie_efx_phys',horiz_only,'I','mW/m2', 'AMIE energy flux') 
       call addfld ('amie_kev_phys',horiz_only,'I','keV'  , 'AMIE mean energy')
    end if
    if ( trim(ionos_epotential_model) == 'weimer' ) then
       call weimer05_init(wei05_coefs_file)
    endif

  end subroutine ionosphere_init

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_run1(pbuf2d)
    use physics_buffer, only: physics_buffer_desc
    use cam_history,    only: outfld, write_inithist
    use phys_grid,      only: get_ncols_p

    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! local vars
    integer :: i, j, k, lchnk  ! indices
    integer :: ifirstxy, ilastxy, jfirstxy, jlastxy, km, idim
    real(r8), allocatable :: tmp(:,:)
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    type(t_fvdycore_grid), pointer :: grid

    real(r8), pointer :: pbuf_amie_efxg(:)     ! Pointer to access AMIE energy flux in pbuf
    real(r8), pointer :: pbuf_amie_kevg(:)     ! Pointer to access AMIE mean energy in pbuf
    
    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude in
    integer :: blksiz                ! number of columns in 2D block
    integer :: tsize                 ! amount of data per grid point passed to physics
    integer :: iam, astat
    integer :: ib, ic, jc,ncol
    integer, allocatable, dimension(:,:) :: bpter
                                     ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
    real(r8), allocatable, dimension(:) :: bbuffer, cbuffer
    real(r8), allocatable :: amie_efxg(:,:) ! energy flux from AMIE
    real(r8), allocatable :: amie_kevg(:,:) ! characteristic mean energy from AMIE

    grid => get_dyn_state_grid()
    iam = grid%iam

    ifirstxy     = grid%ifirstxy
    ilastxy      = grid%ilastxy
    jfirstxy     = grid%jfirstxy
    jlastxy      = grid%jlastxy
    km           = grid%km

    if( write_inithist() .and. ionos_xport_active ) then

       allocate( tmp(ifirstxy:ilastxy,km) )

       idim = ilastxy - ifirstxy + 1
       do j = jfirstxy, jlastxy
          do k = 1, km
             do i = ifirstxy, ilastxy
                tmp(i,k) = opmmrtm1_blck(i,j,k)
             enddo
          enddo
          call outfld ('OpTM1&IC', tmp, idim, j) 
       enddo

       deallocate( tmp )

    endif

    amie_active: if ( ionos_epotential_amie ) then
       allocate(amie_efxg(ifirstxy:ilastxy,jfirstxy:jlastxy))
       allocate(amie_kevg(ifirstxy:ilastxy,jfirstxy:jlastxy))

       ! data assimilated potential
       call d_pie_epotent( ionos_epotential_model, epot_crit_colats, &
                           i0=ifirstxy,i1=ilastxy,j0=jfirstxy,j1=jlastxy, &
                           efxg=amie_efxg,kevg=amie_kevg )

       ! transform to physics grid for aurora...

       ! blocks --> physics chunks

       blcks2phys_local: if (local_dp_map) then

          chnk_loop1 : do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)
             call get_lon_all_p(lchnk, ncol, lons)
             call get_lat_all_p(lchnk, ncol, lats)

             pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
             call pbuf_get_field(pbuf_chnk, indxAMIEefxg, pbuf_amie_efxg)
             call pbuf_get_field(pbuf_chnk, indxAMIEkevg, pbuf_amie_kevg)

             do i=1,ncol
                ic = lons(i)
                jc = lats(i)
                pbuf_amie_efxg(i) = amie_efxg(ic,jc)
                pbuf_amie_kevg(i) = amie_kevg(ic,jc)
             end do
             call outfld ( 'amie_efx_phys', pbuf_amie_efxg, pcols, lchnk )
             call outfld ( 'amie_kev_phys', pbuf_amie_kevg, pcols, lchnk )
          end do chnk_loop1

       else ! blcks2phys_local

          tsize = 2
          blksiz = (jlastxy-jfirstxy+1)*(ilastxy-ifirstxy+1)
          allocate( bpter(blksiz,0:km),stat=astat )
          allocate( bbuffer(tsize*block_buf_nrecs),stat=astat )
          allocate( cbuffer(tsize*chunk_buf_nrecs),stat=astat )

          if (iam < grid%npes_xy) then 
             call block_to_chunk_send_pters(iam+1,blksiz,pver+1,tsize,bpter)
          endif

          do j=jfirstxy,jlastxy
             do i=ifirstxy,ilastxy
                ib = (j-jfirstxy)*(ilastxy-ifirstxy+1) + (i-ifirstxy+1)
                bbuffer(bpter(ib,0)+0) = amie_efxg(i,j)
                bbuffer(bpter(ib,0)+1) = amie_kevg(i,j)
             end do
          end do

          call transpose_block_to_chunk(tsize, bbuffer, cbuffer)

          chnk_loop2: do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)

             pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
             call pbuf_get_field(pbuf_chnk, indxAMIEefxg, pbuf_amie_efxg)
             call pbuf_get_field(pbuf_chnk, indxAMIEkevg, pbuf_amie_kevg)
             call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)
             do i=1,ncol
                pbuf_amie_efxg(i) = cbuffer(cpter(i,0)+0)
                pbuf_amie_kevg(i) = cbuffer(cpter(i,0)+1)
             end do
             call outfld ( 'amie_efx_phys', pbuf_amie_efxg, pcols, lchnk )
             call outfld ( 'amie_kev_phys', pbuf_amie_kevg, pcols, lchnk )
          end do chnk_loop2

          deallocate(bpter)
          deallocate(bbuffer)
          deallocate(cbuffer)


       end if blcks2phys_local

       deallocate(amie_efxg,amie_kevg)

    else
       
       ! set cross tail potential before physics -- aurora uses weimer derived potential
       call d_pie_epotent( ionos_epotential_model, epot_crit_colats )

    end if amie_active

  end subroutine ionosphere_run1

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_run2( phys_state, dyn_in, pbuf2d )

    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc
    use dyn_comp,       only: dyn_import_t
    use cam_history,    only: outfld, write_inithist

    ! - pull some fields from pbuf and dyn_in
    ! - invoke ionosphere/electro-dynamics coupling
    ! - push some fields back to physics via pbuf...

    ! args
    type(physics_state),    intent(in) :: phys_state(begchunk:endchunk)
    type(dyn_import_t),  intent(inout) :: dyn_in  ! dynamics inputs
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local vars
    integer :: i,j,k, lchnk
    integer :: astat

    integer, allocatable, dimension(:,:) :: bpter
                                     ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
    real(r8), allocatable, dimension(:) :: bbuffer, cbuffer

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    real(r8), pointer :: sigma_ped_phys(:,:)  ! physics pointer to Pedersen Conductivity
    real(r8), pointer :: sigma_hall_phys(:,:) ! physics pointer fo Hall Conductivity
    real(r8), pointer :: te_phys(:,:)         ! te from pbuf
    real(r8), pointer :: ti_phys(:,:)         ! ti from pbuf
    real(r8), pointer :: mmrPO2p_phys(:,:)    ! Pointer to access O2+ in pbuf
    real(r8), pointer :: mmrPNOp_phys(:,:)    ! Pointer to access NO+ in pbuf
    real(r8), pointer :: mmrPN2p_phys(:,:)    ! Pointer to access N2+ in pbuf
    real(r8), pointer :: mmrPOp_phys(:,:)     ! Pointer to access O+ in pbuf
!
! Empirical ion drifts from exbdrift (to be converted to blocked for dpie_coupling):
    real(r8), pointer :: ui_phys(:,:)         ! zonal ion drift from pbuf
    real(r8), pointer :: vi_phys(:,:)         ! meridional ion drift from pbuf
    real(r8), pointer :: wi_phys(:,:)         ! vertical ion drift from pbuf

    real(r8), pointer :: o2pmmr_blck(:,:,:) => null()     ! O2+ (blocks)
    real(r8), pointer :: nopmmr_blck(:,:,:) => null()     ! NO+ (blocks)
    real(r8), pointer :: n2pmmr_blck(:,:,:) => null()     ! N2+ (blocks)
    real(r8), pointer :: opmmr_blck(:,:,:)  => null()     ! O+ (blocks)

    real(r8), pointer :: tracer(:,:,:,:)
    real(r8), pointer :: u3s(:,:,:)
    real(r8), pointer :: v3s(:,:,:)
    real(r8), pointer :: pexy(:,:,:)

    real(r8), pointer :: phis(:,:)            ! surface geopotential

    real(r8), pointer :: o2mmr_blck(:,:,:)
    real(r8), pointer :: o1mmr_blck(:,:,:)
    real(r8), pointer :: h1mmr_blck(:,:,:)

    integer :: ib, ic, jc, ifirstxy, ilastxy, jfirstxy, jlastxy, km, ncol

    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: nSIons                        ! number of ions set to non-advected
    integer :: ibuffOp,ibuffO2p,ibuffNOp, ibuffN2p ! Buffer indices for non-advected ions

    integer :: blksiz                 ! number of columns in 2D block
    integer :: tsize                  ! amount of data per grid point passed to physics
    integer :: iam

    real(r8), allocatable :: wuxy(:,:,:)
    real(r8), allocatable :: wvxy(:,:,:)
    real(r8), allocatable :: sigma_ped_blck (:,:,:)
    real(r8), allocatable :: sigma_hall_blck(:,:,:)
    real(r8), allocatable :: ti_blck(:,:,:)
    real(r8), allocatable :: te_blck(:,:,:)
    real(r8), allocatable :: zi_blck(:,:,:)
    real(r8), allocatable :: zm_blck(:,:,:)
    real(r8), allocatable :: ui_blck(:,:,:)
    real(r8), allocatable :: vi_blck(:,:,:)
    real(r8), allocatable :: wi_blck(:,:,:)
    real(r8), allocatable :: omega_blck(:,:,:)
    real(r8), allocatable :: tn_blck(:,:,:)

    type (t_fvdycore_grid), pointer :: grid

    ionos_cpl: if (ionos_xport_active) then 

       grid => get_dyn_state_grid()
       iam = grid%iam

       allocate( wuxy(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( wvxy(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( sigma_ped_blck (grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( sigma_hall_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( ti_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( te_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( zi_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( zm_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( ui_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( vi_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( wi_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( omega_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )
       allocate( tn_blck(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, grid%km) )

       ifirstxy = grid%ifirstxy
       ilastxy  = grid%ilastxy
       jfirstxy = grid%jfirstxy
       jlastxy  = grid%jlastxy
       km = grid%km

       phis   => dyn_in%phis

       tracer => dyn_in%tracer
       pexy   => dyn_in%pe

       u3s    => dyn_in%u3s
       v3s    => dyn_in%v3s

       if (iam < grid%npes_xy) then 
          call d2a3dijk( grid, u3s, v3s, wuxy, wvxy )
       endif

       if (sIndxOp>0) then
          allocate(opmmr_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
          if (astat /= 0) call endrun('ionos_intr_d_p_cplng: failed to allocate opmmr_blck')
       endif
       if (sIndxO2p>0) then
          allocate(o2pmmr_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
          if (astat /= 0) call endrun('ionos_intr_d_p_cplng: failed to allocate o2pmmr_blck')
       endif
       if (sIndxNOp>0) then
          allocate(nopmmr_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
          if (astat /= 0) call endrun('ionos_intr_d_p_cplng: failed to allocate nopmmr_blck')
       endif
       if (sIndxN2p>0) then
          allocate(n2pmmr_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
          if (astat /= 0) call endrun('ionos_intr_d_p_cplng: failed to allocate n2pmmr_blck')
       endif

       phys2blcks_local: if (local_dp_map) then

          do lchnk = begchunk,endchunk

             ncol = get_ncols_p(lchnk)
             call get_lon_all_p(lchnk, ncol, lons)
             call get_lat_all_p(lchnk, ncol, lats)
             pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

             ! Get Pedersen and Hall conductivities:
             call pbuf_get_field(pbuf_chnk, index_ped,  sigma_ped_phys)
             call pbuf_get_field(pbuf_chnk, index_hall, sigma_hall_phys)
             do k=1,km
                do i=1,ncol
                   sigma_ped_blck(lons(i),lats(i),k) = sigma_ped_phys(i,k)
                   sigma_hall_blck(lons(i),lats(i),k) = sigma_hall_phys(i,k)
                end do
             enddo

             ! Get ion and electron temperatures 
             call pbuf_get_field(pbuf_chnk, index_te, te_phys)
             call pbuf_get_field(pbuf_chnk, index_ti, ti_phys)
             do k=1,km
                do i=1,ncol
                   te_blck(lons(i),lats(i),k) = te_phys(i,k)
                   ti_blck(lons(i),lats(i),k) = ti_phys(i,k)
                end do
             enddo

             ! Get components of ion drift velocities
             call pbuf_get_field(pbuf_chnk, index_ui, ui_phys)
             call pbuf_get_field(pbuf_chnk, index_vi, vi_phys)
             call pbuf_get_field(pbuf_chnk, index_wi, wi_phys)
             do k=1,km
                do i=1,ncol
                   ui_blck(lons(i),lats(i),k) = ui_phys(i,k)
                   vi_blck(lons(i),lats(i),k) = vi_phys(i,k)
                   wi_blck(lons(i),lats(i),k) = wi_phys(i,k)
                   zi_blck(lons(i),lats(i),k)    = phys_state(lchnk)%zi(i,k)
                   zm_blck(lons(i),lats(i),k)    = phys_state(lchnk)%zm(i,k)
                   omega_blck(lons(i),lats(i),k) = phys_state(lchnk)%omega(i,k)
                   tn_blck(lons(i),lats(i),k)    = phys_state(lchnk)%t(i,k)
                enddo
             enddo

             !--------------------------------------------------------
             ! Get ions from physics buffer if non-transported
             !--------------------------------------------------------
             if (sIndxO2p > 0) then
                call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPO2p_phys, &
                     start=(/1,1,sIndxO2p/), kount=(/pcols,pver,1/) )
                do k=1,km
                   do i=1,ncol
                      o2pmmr_blck(lons(i),lats(i),k) = mmrPO2p_phys(i,k)
                   end do
                enddo
             endif
             if (sIndxNOp > 0) then
                call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPNOp_phys, &
                     start=(/1,1,sIndxNOp/), kount=(/pcols,pver,1/) )
                do k=1,km
                   do i=1,ncol
                      nopmmr_blck(lons(i),lats(i),k) = mmrPNOp_phys(i,k)
                   end do
                enddo
             endif
             if (sIndxN2p > 0) then
                call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPN2p_phys, &
                     start=(/1,1,sIndxN2p/), kount=(/pcols,pver,1/) )
                do k=1,km
                   do i=1,ncol
                      n2pmmr_blck(lons(i),lats(i),k) = mmrPN2p_phys(i,k)
                   end do
                enddo
             endif
             if (sIndxOp > 0) then
                call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPOp_phys, &
                     start=(/1,1,sIndxOp/), kount=(/pcols,pver,1/) )
                do k=1,km
                   do i=1,ncol
                      opmmr_blck(lons(i),lats(i),k) = mmrPOp_phys(i,k)
                   end do
                enddo
             endif

          enddo ! do lchnk = begchunk,endchunk

       else ! phys2blcks_local

          tsize = 11

          nSIons = 0
          if (sIndxOp > 0)  then 
             ibuffOp = tsize + nSIons
             nSIons = nSIons + 1
          endif
          if (sIndxO2p > 0) then
             ibuffO2p = tsize + nSIons
             nSIons = nSIons + 1
          endif
          if (sIndxNOp > 0) then
             ibuffNOp = tsize + nSIons
             nSIons = nSIons + 1
          endif
          if (sIndxN2p > 0) then
             ibuffN2p = tsize + nSIons
             nSIons = nSIons + 1
          endif
          tsize = tsize + nSIons

          blksiz = (jlastxy-jfirstxy+1)*(ilastxy-ifirstxy+1)
          allocate(bpter(blksiz,0:km))
          allocate(bbuffer(tsize*block_buf_nrecs))
          allocate(cbuffer(tsize*chunk_buf_nrecs))

          do lchnk = begchunk,endchunk
             ncol = get_ncols_p(lchnk)

             pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

             ! Get Pedersen and Hall conductivities:
             call pbuf_get_field(pbuf_chnk, index_ped,  sigma_ped_phys)
             call pbuf_get_field(pbuf_chnk, index_hall, sigma_hall_phys)

             ! Get ion and electron temperatures 
             call pbuf_get_field(pbuf_chnk, index_te,  te_phys)
             call pbuf_get_field(pbuf_chnk, index_ti,  ti_phys)

             ! Get components of ion drift velocities
             call pbuf_get_field(pbuf_chnk, index_ui,  ui_phys)
             call pbuf_get_field(pbuf_chnk, index_vi,  vi_phys)
             call pbuf_get_field(pbuf_chnk, index_wi,  wi_phys)
 
             !--------------------------------------------------------
             ! Get ions from physics buffer if non-transported
             !--------------------------------------------------------

             if (sIndxOp > 0)  call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPOp_phys,  &
                  start=(/1,1,sIndxOp/), kount=(/pcols,pver,1/) )
             if (sIndxO2p > 0) call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPO2p_phys, &
                  start=(/1,1,sIndxO2p/), kount=(/pcols,pver,1/) )     
             if (sIndxNOp > 0) call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPNOp_phys, &
                  start=(/1,1,sIndxNOp/), kount=(/pcols,pver,1/) )     
             if (sIndxN2p > 0) call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPN2p_phys, &
                  start=(/1,1,sIndxN2p/), kount=(/pcols,pver,1/) )

             call chunk_to_block_send_pters(lchnk,pcols,pver+1,tsize,cpter)

             do i=1,ncol
                cbuffer(cpter(i,0):cpter(i,0)+tsize-1) = 0.0_r8
             end do

             do k=1,km
                do i=1,ncol

                   cbuffer(cpter(i,k)+0) = sigma_ped_phys(i,k)
                   cbuffer(cpter(i,k)+1) = sigma_hall_phys(i,k)
                   cbuffer(cpter(i,k)+2) = te_phys(i,k)
                   cbuffer(cpter(i,k)+3) = ti_phys(i,k)
                   cbuffer(cpter(i,k)+4) = phys_state(lchnk)%zi(i,k)
                   cbuffer(cpter(i,k)+5) = phys_state(lchnk)%zm(i,k)
                   cbuffer(cpter(i,k)+6) = ui_phys(i,k)
                   cbuffer(cpter(i,k)+7) = vi_phys(i,k)
                   cbuffer(cpter(i,k)+8) = wi_phys(i,k)
                   cbuffer(cpter(i,k)+9) = phys_state(lchnk)%omega(i,k)
                   cbuffer(cpter(i,k)+10) = phys_state(lchnk)%t(i,k)

                   if (sIndxO2p > 0)cbuffer(cpter(i,k)+ibuffO2p) = mmrPO2p_phys(i,k)
                   if (sIndxNOp > 0)cbuffer(cpter(i,k)+ibuffNOp) = mmrPNOp_phys(i,k)
                   if (sIndxN2p > 0)cbuffer(cpter(i,k)+ibuffN2p) = mmrPN2p_phys(i,k)
                   if (sIndxOp > 0) cbuffer(cpter(i,k)+ibuffOp)  = mmrPOp_phys(i,k)

                end do

             end do

          end do

          call t_barrierf('sync_chk_to_blk', grid%commxy)
          call t_startf ('chunk_to_block')
          call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
          call t_stopf  ('chunk_to_block')

          if (iam < grid%npes_xy) then 
             call chunk_to_block_recv_pters(iam+1,blksiz,pver+1,tsize,bpter)
          endif

          do j=jfirstxy,jlastxy
             do k=1,km
                do i=ifirstxy,ilastxy
                   ib = (j-jfirstxy)*(ilastxy-ifirstxy+1) + (i-ifirstxy+1)

                   sigma_ped_blck(i,j,k)  = bbuffer(bpter(ib,k)+0)
                   sigma_hall_blck(i,j,k) = bbuffer(bpter(ib,k)+1)
                   te_blck(i,j,k)         = bbuffer(bpter(ib,k)+2)
                   ti_blck(i,j,k)         = bbuffer(bpter(ib,k)+3)
                   zi_blck(i,j,k)         = bbuffer(bpter(ib,k)+4)
                   zm_blck(i,j,k)         = bbuffer(bpter(ib,k)+5)
                   ui_blck(i,j,k)         = bbuffer(bpter(ib,k)+6)
                   vi_blck(i,j,k)         = bbuffer(bpter(ib,k)+7)
                   wi_blck(i,j,k)         = bbuffer(bpter(ib,k)+8)
                   omega_blck(i,j,k)      = bbuffer(bpter(ib,k)+9)
                   tn_blck(i,j,k)         = bbuffer(bpter(ib,k)+10)

                   if (sIndxO2p > 0) o2pmmr_blck(i,j,k) = bbuffer(bpter(ib,k)+ibuffO2p)
                   if (sIndxNOp > 0) nopmmr_blck(i,j,k) = bbuffer(bpter(ib,k)+ibuffNOp)
                   if (sIndxN2p > 0) n2pmmr_blck(i,j,k) = bbuffer(bpter(ib,k)+ibuffN2p)
                   if (sIndxOp > 0)  opmmr_blck(i,j,k)  = bbuffer(bpter(ib,k)+ibuffOp)

                enddo
             enddo
          enddo

          deallocate(bpter)
          deallocate(bbuffer)
          deallocate(cbuffer)

       endif phys2blcks_local

       !-------------------------------------------------------------------------------------------
       !  Set dpie_coupling input ions if they are advected ...
       !-------------------------------------------------------------------------------------------
       if (ixo2p > 0) then
          o2pmmr_blck => tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixo2p)
       endif
       if (ixnop > 0) then
          nopmmr_blck => tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixnop)
       endif
       if (ixn2p > 0) then
          n2pmmr_blck => tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixn2p)
       endif
       if (ixop > 0) then
          opmmr_blck => tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixop)
       endif

       !------------------------------------
       ! Get neutrals from advected tracers array
       !------------------------------------

       o2mmr_blck  => tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixo2)
       o1mmr_blck  => tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixo)
       h1mmr_blck  => tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixh)

       !
       !   Make geopotential height (m) for d_pie_coupling. 
       !
       do k=1,km
          do j=jfirstxy,jlastxy
             do i=ifirstxy,ilastxy
                zi_blck(i,j,k) = zi_blck(i,j,k)+phis(i,j)/gravit ! phis is redundant in k
                zm_blck(i,j,k) = zm_blck(i,j,k)+phis(i,j)/gravit ! phis is redundant in k
             enddo
          enddo
       enddo

       call t_startf('d_pie_coupling')

       if (iam < grid%npes_xy) then 
          ! waccmx ionosphere electro-dynamics -- transports O+ and provides updates to ion drift velocities
          call d_pie_coupling(omega_blck,pexy,zi_blck,zm_blck,wuxy,wvxy,tn_blck,                        &
               sigma_ped_blck,sigma_hall_blck,te_blck,ti_blck,                      &
               o2mmr_blck,o1mmr_blck,h1mmr_blck,o2pmmr_blck,nopmmr_blck,n2pmmr_blck, &
               opmmr_blck,opmmrtm1_blck,ui_blck,vi_blck,wi_blck,                                &
               rmassO2,rmassO1,rmassH,rmassN2,rmassO2p,rmassNOp,rmassN2p, rmassOp,       &
               ifirstxy,ilastxy, jfirstxy,jlastxy)
       endif

       call t_stopf ('d_pie_coupling')

       !
       !----------------------------------------
       !  Put data back in to state%q or pbuf
       !----------------------------------------
       if (ixop > 0) then
          tracer(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km,ixop) = opmmr_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km)           
       endif

       ! blocks --> physics chunks

       blcks2phys_local: if (local_dp_map) then

          chnk_loop1 : do lchnk = begchunk,endchunk
             ncol = phys_state(lchnk)%ncol
             call get_lon_all_p(lchnk, ncol, lons)
             call get_lat_all_p(lchnk, ncol, lats)

             pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

             call pbuf_get_field(pbuf_chnk, index_ui, ui_phys)
             call pbuf_get_field(pbuf_chnk, index_vi, vi_phys)
             call pbuf_get_field(pbuf_chnk, index_wi, wi_phys)
             if (sIndxOp > 0) then
                call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPOp_phys, &
                     start=(/1,1,sIndxOp/), kount=(/pcols,pver,1/) )
             endif
             do k=1,km
                do i=1,ncol
                   ic = lons(i)
                   jc = lats(i)
                   ui_phys(i,k) = ui_blck(ic,jc,k)
                   vi_phys(i,k) = vi_blck(ic,jc,k)
                   wi_phys(i,k) = wi_blck(ic,jc,k)
                   if (sIndxOp > 0) mmrPOp_phys(i,k) = opmmr_blck(ic,jc,k)
                end do
             end do

             if (ionos_edyn_active) then
                call outfld ( 'UI', ui_phys, pcols, lchnk )
                call outfld ( 'VI', vi_phys, pcols, lchnk )
                call outfld ( 'WI', wi_phys, pcols, lchnk )
                if (write_inithist()) then
                   call outfld ( 'UI&IC', ui_phys, pcols, lchnk )
                   call outfld ( 'VI&IC', vi_phys, pcols, lchnk )
                   call outfld ( 'WI&IC', wi_phys, pcols, lchnk )
                endif
             endif

          end do chnk_loop1

       else ! blcks2phys_local

          if (sIndxOp > 0) then
             tsize = 4 ! for ui,vi,wi,op
          else
             tsize = 3 ! for ui,vi,wi
          endif
          tsize=tsize+1

          blksiz = (jlastxy-jfirstxy+1)*(ilastxy-ifirstxy+1)
          allocate( bpter(blksiz,0:km),stat=astat )
          allocate( bbuffer(tsize*block_buf_nrecs),stat=astat )
          allocate( cbuffer(tsize*chunk_buf_nrecs),stat=astat )

          if (iam < grid%npes_xy) then 
             call block_to_chunk_send_pters(iam+1,blksiz,km+1,tsize,bpter)
          endif

          do j=jfirstxy,jlastxy
             do i=ifirstxy,ilastxy
                ib = (j-jfirstxy)*(ilastxy-ifirstxy+1) + (i-ifirstxy+1)

                do k=1,km

                   bbuffer(bpter(ib,k)) = ui_blck(i,j,k)
                   bbuffer(bpter(ib,k)+1) = vi_blck(i,j,k)
                   bbuffer(bpter(ib,k)+2) = wi_blck(i,j,k)
                   if (sIndxOp > 0) bbuffer(bpter(ib,k)+3) = opmmr_blck(i,j,k)

                end do
             end do
          end do

          call t_barrierf('sync_ionos_blk_to_chk', grid%commxy)
          call t_startf ('ionos_block_to_chunk')
          call transpose_block_to_chunk(tsize, bbuffer, cbuffer)
          call t_stopf  ('ionos_block_to_chunk')

          chnk_loop2: do lchnk = begchunk,endchunk
             ncol = phys_state(lchnk)%ncol

             pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

             call pbuf_get_field(pbuf_chnk, index_ui, ui_phys)
             call pbuf_get_field(pbuf_chnk, index_vi, vi_phys)
             call pbuf_get_field(pbuf_chnk, index_wi, wi_phys)
             if (sIndxOp > 0) then
                call pbuf_get_field(pbuf_chnk, slvd_pbf_ndx, mmrPOp_phys, &
                     start=(/1,1,sIndxOp/), kount=(/pcols,pver,1/) )
             endif

             call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

             do i=1,ncol

                do k=1,km
                   ui_phys(i,k) = cbuffer(cpter(i,k))
                   vi_phys(i,k) = cbuffer(cpter(i,k)+1)
                   wi_phys(i,k) = cbuffer(cpter(i,k)+2)
                   if (sIndxOp > 0) then
                      mmrPOp_phys(i,k) = cbuffer(cpter(i,k)+3)
                   endif
                end do ! k=1,km
             end do ! i=1,ncol

             if (ionos_edyn_active) then
                call outfld ( 'UI', ui_phys, pcols, lchnk )
                call outfld ( 'VI', vi_phys, pcols, lchnk )
                call outfld ( 'WI', wi_phys, pcols, lchnk )
                if (write_inithist()) then
                   call outfld ( 'UI&IC', ui_phys, pcols, lchnk )
                   call outfld ( 'VI&IC', vi_phys, pcols, lchnk )
                   call outfld ( 'WI&IC', wi_phys, pcols, lchnk )
                endif
             endif

          end do chnk_loop2

          deallocate(bpter)
          deallocate(bbuffer)
          deallocate(cbuffer)

       endif blcks2phys_local

       if (sIndxOp>0) then
          deallocate(opmmr_blck)
          nullify(opmmr_blck)
       endif
       if (sIndxO2p>0) then
          deallocate(o2pmmr_blck)
          nullify(o2pmmr_blck)
       endif
       if (sIndxNOp>0) then
          deallocate(nopmmr_blck)
          nullify(nopmmr_blck)
       endif
       if (sIndxN2p>0) then
          deallocate(n2pmmr_blck)
          nullify(n2pmmr_blck)
       endif

       deallocate( wuxy )
       deallocate( wvxy )
       deallocate( sigma_ped_blck )
       deallocate( sigma_hall_blck )
       deallocate( ti_blck )
       deallocate( te_blck )
       deallocate( zi_blck )
       deallocate( ui_blck )
       deallocate( vi_blck )
       deallocate( wi_blck )
       deallocate( omega_blck )
       deallocate( tn_blck )

    endif ionos_cpl

  end subroutine ionosphere_run2

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
 subroutine ionosphere_init_restart(File)
    use pio, only: file_desc_t, pio_double, pio_def_var
    use cam_pio_utils, only: cam_pio_def_dim
    use dyn_grid,      only: get_horiz_grid_dim_d

    type(File_desc_t),  intent(inout) :: File

    integer :: ierr,hdim1,hdim2, dimids(3)

    call get_horiz_grid_dim_d(hdim1, hdim2)

    call cam_pio_def_dim(File, 'lon',  hdim1, dimids(1),  existOK=.true.)
    call cam_pio_def_dim(File, 'lat',  hdim2, dimids(2),  existOK=.true.)
    call cam_pio_def_dim(File, 'lev',  pver,  dimids(3),  existOK=.true.)

    if (ionos_xport_active) then
       ierr = PIO_Def_Var(File, 'Optm1', pio_double, dimids, Optm1_vdesc)
    endif
  end subroutine ionosphere_init_restart

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_write_restart(File)
    use pio, only: io_desc_t, file_desc_t, pio_write_darray, pio_initdecomp, pio_double
    use cam_pio_utils, only: pio_subsystem
    use dyn_grid,      only: get_horiz_grid_dim_d

    type(File_desc_t), intent(inout) :: File

    type(io_desc_t) :: iodesc3d
    integer :: hdim1, hdim2
    integer, pointer :: ldof(:)
    integer :: ierr

    if (ionos_xport_active) then
       call get_horiz_grid_dim_d(hdim1, hdim2)
       ldof => get_restart_decomp(hdim1, hdim2, pver)
       call pio_initdecomp(pio_subsystem, pio_double, (/hdim1, hdim2, pver/), ldof, iodesc3d)
       deallocate(ldof)

       call pio_write_darray(File, Optm1_vdesc, iodesc3d, opmmrtm1_blck, ierr)
    endif

  end subroutine ionosphere_write_restart

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_read_restart(File)
    use pio, only: io_desc_t, file_desc_t, pio_inq_varid, pio_read_darray, pio_initdecomp, pio_double
    use cam_pio_utils, only: pio_subsystem
    use dyn_grid,      only: get_horiz_grid_dim_d

    type(file_desc_t), intent(inout) :: File

    integer :: ierr
    type(io_desc_t) :: iodesc3d
    integer :: hdim1, hdim2
    integer, pointer :: ldof(:)

    if (ionos_xport_active) then
       call ionosphere_alloc

       call get_horiz_grid_dim_d(hdim1, hdim2)
       ldof => get_restart_decomp(hdim1, hdim2, pver)
       call pio_initdecomp(pio_subsystem, pio_double, (/hdim1, hdim2, pver/), ldof, iodesc3d)
       deallocate(ldof)

       ierr = pio_inq_varid(File, 'Optm1', Optm1_vdesc)
       call pio_read_darray(File, Optm1_vdesc, iodesc3d, opmmrtm1_blck, ierr)
    endif

  end subroutine ionosphere_read_restart

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_final

#ifdef WACCMX_EDYN_ESMF
    use edyn_esmf, only: edyn_esmf_final

    call edyn_esmf_final()
#endif

    if (allocated(opmmrtm1_blck)) deallocate(opmmrtm1_blck)

  end subroutine ionosphere_final

!=========================================================================================
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_read_ic()

    use pio,          only: file_desc_t
    use ncdio_atm,    only: infld
    use cam_initfiles,      only: initial_file_get_id

    type(file_desc_t), pointer :: fh_ini    ! PIO filehandle

    type (t_fvdycore_grid), pointer :: grid
    integer :: ifirstxy,ilastxy,jfirstxy,jlastxy,km
    logical :: readvar

    if ( ionos_xport_active ) then
       call ionosphere_alloc()

       fh_ini   => initial_file_get_id()
       grid     => get_dyn_state_grid()
       ifirstxy =  grid%ifirstxy
       ilastxy  =  grid%ilastxy
       jfirstxy =  grid%jfirstxy
       jlastxy  =  grid%jlastxy
       km       =  grid%km

       ! try reading in OpTM1 from the IC file
       call infld('OpTM1', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
            1, km, opmmrtm1_blck, readvar, gridname='fv_centers')

       if (.not.readvar) then
          ! if OpTM1 is not included in the IC file then try using O+
          call infld('Op', fh_ini, 'lon', 'lat', 'lev', ifirstxy, ilastxy, jfirstxy, jlastxy, &
               1, km, opmmrtm1_blck, readvar, gridname='fv_centers')
       endif
    endif

  end subroutine ionosphere_read_ic

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_alloc

    type(T_FVDYCORE_GRID),pointer :: grid ! FV Dynamics grid
    integer :: ifirstxy, ilastxy, jfirstxy, jlastxy, km
    integer :: astat

    if (.not. allocated(opmmrtm1_blck)) then

       grid => get_dyn_state_grid()
       ifirstxy = grid%ifirstxy
       ilastxy  = grid%ilastxy
       jfirstxy = grid%jfirstxy
       jlastxy  = grid%jlastxy
       km = grid%km

       allocate(opmmrtm1_blck(ifirstxy:ilastxy,jfirstxy:jlastxy,km),stat=astat)
       if (astat /= 0) call endrun('ionosphere_init: failed to allocate opmmrtm1_blck')
       opmmrtm1_blck = 0._r8

    endif

  end subroutine ionosphere_alloc


  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
function get_restart_decomp(hdim1, hdim2, nlev) result(ldof)
   use dyn_grid, only: get_dyn_grid_parm

   ! Get the integer mapping of a variable in the dynamics decomp in memory.  
   ! The canonical ordering is as on the file. A 0 value indicates that the
   ! variable is not on the file (eg halo or boundary values)

   ! arguments
   integer, intent(in) :: hdim1, hdim2, nlev
   integer, pointer :: ldof(:)

   ! local variables
   integer :: i, k, j
   integer :: lcnt
   integer :: beglatxy, beglonxy, endlatxy, endlonxy
   !----------------------------------------------------------------------------

   beglonxy = get_dyn_grid_parm('beglonxy')
   endlonxy = get_dyn_grid_parm('endlonxy')
   beglatxy = get_dyn_grid_parm('beglatxy')
   endlatxy = get_dyn_grid_parm('endlatxy')

   lcnt = (endlatxy-beglatxy+1)*nlev*(endlonxy-beglonxy+1)
   allocate(ldof(lcnt))
   ldof(:) = 0	

   lcnt = 0
   do k = 1, nlev
      do j = beglatxy, endlatxy
         do i = beglonxy, endlonxy
            lcnt = lcnt + 1
            ldof(lcnt) = i + (j-(plat-hdim2+1))*hdim1+(k-1)*hdim1*hdim2
         end do
      end do
   end do

end function get_restart_decomp

!=========================================================================================


end module ionosphere_interface
