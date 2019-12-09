module atm_import_export

  use NUOPC           , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_IsConnected
  use NUOPC_Model     , only : NUOPC_ModelGet
  use ESMF            , only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF            , only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF            , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF            , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF            , only : operator(/=), operator(==)
  use shr_kind_mod    , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs, cx=>shr_kind_cx
  use shr_sys_mod     , only : shr_sys_abort
  use atm_shr_methods , only : chkerr
  use cam_logfile     , only : iulog
  use srf_field_check , only : set_active_Sl_ram1
  use srf_field_check , only : set_active_Sl_fv
  use srf_field_check , only : set_active_Sl_soilw
  use srf_field_check , only : set_active_Fall_flxdst1
  use srf_field_check , only : set_active_Fall_flxvoc
  use srf_field_check , only : set_active_Fall_flxfire
  use srf_field_check , only : set_active_Fall_fco2_lnd
  use srf_field_check , only : set_active_Faoo_fco2_ocn
  use srf_field_check , only : set_active_Faxa_nhx
  use srf_field_check , only : set_active_Faxa_noy

  implicit none
  private ! except

  public  :: advertise_fields
  public  :: realize_fields
  public  :: import_fields
  public  :: export_fields
  public  :: state_getfldptr

  private :: fldlist_add
  private :: fldlist_realize

  type fldlist_type
     character(len=128) :: stdname
     integer :: ungridded_lbound = 0
     integer :: ungridded_ubound = 0
  end type fldlist_type

  integer             , parameter         :: fldsMax = 100
  integer             , public, protected :: fldsToAtm_num = 0
  integer             , public, protected :: fldsFrAtm_num = 0
  type (fldlist_type) , public, protected :: fldsToAtm(fldsMax)
  type (fldlist_type) , public, protected :: fldsFrAtm(fldsMax)

  ! from lnd->atm
  character(len=cx)      :: carma_fields ! list of CARMA fields from lnd->atm
  integer                :: drydep_nflds ! number of dry deposition velocity fields lnd-> atm
  integer                :: megan_nflds  ! number of MEGAN voc fields from lnd-> atm
  integer                :: emis_nflds   ! number of fire emission fields from lnd-> atm

  ! from atm-<lnd
  integer, public        :: ndep_nflds   ! number  of nitrogen deposition fields from atm->lnd/ocn 

  integer                :: dbug_flag = 6    ! ESMF log output 
  integer, parameter     :: debug_import = 0 ! internal debug level
  integer, parameter     :: debug_export = 0 ! internal debug level
  character(*),parameter :: F01 = "('(cam_import_export) ',a,i8,2x,i8,2x,d21.14)"
  character(*),parameter :: F02 = "('(cam_import_export) ',a,i8,2x,i8,2x,i8,2x,d21.14)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, flds_scalar_name, rc)

    use spmd_utils        , only : masterproc
    use seq_drydep_mod    , only : seq_drydep_readnl, seq_drydep_init
    use shr_megan_mod     , only : shr_megan_readnl
    use shr_fire_emis_mod , only : shr_fire_emis_readnl
    use shr_carma_mod     , only : shr_carma_readnl
    use shr_ndep_mod      , only : shr_ndep_readnl

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)       :: importState
    type(ESMF_State)       :: exportState
    character(ESMF_MAXSTR) :: stdname
    character(ESMF_MAXSTR) :: cvalue
    character(len=2)       :: nec_str
    integer                :: n, num
    logical                :: flds_co2a      ! use case
    logical                :: flds_co2b      ! use case
    logical                :: flds_co2c      ! use case
    integer                :: ndep_nflds, megan_nflds, emis_nflds
    character(len=128)     :: fldname
    character(len=*), parameter :: subname='(atm_import_export:advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 10) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! determine necessary toggles for below
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2a
    call ESMF_LogWrite('flds_co2a = '// trim(cvalue), ESMF_LOGMSG_INFO)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2b
    call ESMF_LogWrite('flds_co2b = '// trim(cvalue), ESMF_LOGMSG_INFO)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2c
    call ESMF_LogWrite('flds_co2c = '// trim(cvalue), ESMF_LOGMSG_INFO)

    !--------------------------------
    ! Export fields
    !--------------------------------

    if (dbug_flag > 10) call ESMF_LogWrite(subname//' export fields', ESMF_LOGMSG_INFO)
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, trim(flds_scalar_name))

    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_topo'       )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_z'          )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_u'          )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_v'          )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_tbot'       )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_ptem'       )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_shum'       )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_pbot'       )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_dens'       )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_pslv'       )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_rainc'    )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_rainl'    )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_snowc'    )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_snowl'    )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_lwdn'     )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swndr'    )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swvdr'    )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swndf'    )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swvdf'    )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swnet'    )  ! only diagnostic

    ! from atm - black carbon deposition fluxes (3)
    ! (1) => bcphidry, (2) => bcphodry, (3) => bcphiwet
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_bcph', ungridded_lbound=1, ungridded_ubound=3)

    ! from atm - organic carbon deposition fluxes (3)
    ! (1) => ocphidry, (2) => ocphodry, (3) => ocphiwet
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_ocph', ungridded_lbound=1, ungridded_ubound=3)

    ! from atm - wet dust deposition frluxes (4 sizes)
    ! (1) => dstwet1, (2) => dstwet2, (3) => dstwet3, (4) => dstwet4 
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet', ungridded_lbound=1, ungridded_ubound=4)

    ! from atm - dry dust deposition frluxes (4 sizes)
    ! (1) => dstdry1, (2) => dstdry2, (3) => dstdry3, (4) => dstdry4 
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry', ungridded_lbound=1, ungridded_ubound=4)

    if (dbug_flag > 10) call ESMF_LogWrite(subname//' export fields co2', ESMF_LOGMSG_INFO)

    ! from atm co2 fields
    if (flds_co2a .or. flds_co2b .or. flds_co2c) then
       call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_co2prog' )
       call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_co2diag' )
    end if

    ! from atm - nitrogen deposition
    call shr_ndep_readnl("drv_flds_in", ndep_nflds)
    if (ndep_nflds > 0) then
       call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_ndep', ungridded_lbound=1, ungridded_ubound=ndep_nflds)
       call set_active_Faxa_nhx(.true.)
       call set_active_Faxa_noy(.true.)
    end if

    ! Now advertise above export fields
    if (dbug_flag > 10) call ESMF_LogWrite(subname//' advertise export fields ', ESMF_LOGMSG_INFO)
    do n = 1,fldsFrAtm_num
       call NUOPC_Advertise(exportState, standardName=fldsFrAtm(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !-----------------
    ! Import fields
    !-----------------

    if (dbug_flag > 10) call ESMF_LogWrite(subname//' Import Fields', ESMF_LOGMSG_INFO)

    call fldlist_add(fldsToAtm_num, fldsToAtm, trim(flds_scalar_name))

    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sx_anidr'  )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sx_avsdf'  )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sx_anidf'  )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sx_avsdr'  )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sl_lfrac'  )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Si_ifrac'  )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'So_ofrac'  )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sx_tref'   )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sx_qref'   )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sx_t'      )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'So_t'      )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sl_fv'     );  call set_active_Sl_fv(.true.)
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sl_ram1'   );  call set_active_Sl_ram1(.true.)
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sl_snowh'  )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Si_snowh'  )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'So_ssq'    )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'So_re'     )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sx_u10'    )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Faxx_taux' )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Faxx_tauy' )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Faxx_lat'  )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Faxx_sen'  )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Faxx_lwup' )
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Faxx_evap' )

    ! dust fluxes from land (4 sizes)
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Fall_flxdst', ungridded_lbound=1, ungridded_ubound=4)
    call set_active_Fall_flxdst1(.true.)

    ! co2 fields from land and ocean
    if (flds_co2b .or. flds_co2c) then
       call fldlist_add(fldsToAtm_num, fldsToAtm, 'Fall_fco2_lnd')
       call set_active_Fall_fco2_lnd(.true.)
    end if
    if (flds_co2c) then
       call fldlist_add(fldsToAtm_num, fldsToAtm, 'Faoo_fco2_ocn')
       call set_active_Faoo_fco2_ocn(.true.)
    end if

    ! dry deposition velocities from land - ALSO initialize drydep here
    call seq_drydep_readnl("drv_flds_in", drydep_nflds)
    if (drydep_nflds > 0) then
       call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sl_ddvel', ungridded_lbound=1, ungridded_ubound=drydep_nflds)
    end if
    call seq_drydep_init( )

    ! MEGAN VOC emissions fluxes from land
    call shr_megan_readnl('drv_flds_in', megan_nflds)
    if (megan_nflds > 0) then
       call fldlist_add(fldsToAtm_num, fldsToAtm, 'Fall_voc', ungridded_lbound=1, ungridded_ubound=megan_nflds)
       call set_active_Fall_flxvoc(.true.)
    end if

    ! fire emissions fluxes from land
    call shr_fire_emis_readnl('drv_flds_in', emis_nflds)
    if (emis_nflds > 0) then
       call fldlist_add(fldsToAtm_num, fldsToAtm, 'Fall_fire', ungridded_lbound=1, ungridded_ubound=emis_nflds)
       call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sl_fztop')
       call set_active_Fall_flxfire(.true.)
    end if

    ! CARMA volumetric soil water from land
    call shr_carma_readnl('drv_flds_in', carma_fields)
    if (carma_fields /= ' ') then
       call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sl_soilw') ! optional for carma
       call set_active_Sl_soilw(.true.) ! check for carma
    end if

    ! ------------------------------------------
    ! Now advertise above import fields
    ! ------------------------------------------
    do n = 1,fldsToAtm_num
       call NUOPC_Advertise(importState, standardName=fldsToAtm(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    if (dbug_flag > 10) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine advertise_fields

  !===============================================================================

  subroutine realize_fields(gcomp, Emesh, flds_scalar_name, flds_scalar_num, rc)

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_Mesh)     , intent(in)    :: Emesh
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num 
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    character(len=*), parameter :: subname='(atm_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrAtm, &
         numflds=fldsFrAtm_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':camExport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToAtm, &
         numflds=fldsToAtm_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':camImport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine realize_fields

  !===============================================================================

  subroutine import_fields( gcomp, cam_in, restart_init, rc)

    ! -----------------------------------------------------
    ! Set field pointers in import state and 
    ! copy from field pointer to chunk array data structure
    ! -----------------------------------------------------

    use spmd_utils        , only : masterproc
    use camsrfexch        , only : cam_in_t
    use phys_grid         , only : get_ncols_p
    use ppgrid            , only : begchunk, endchunk
    use shr_const_mod     , only : shr_const_stebol
    use co2_cycle         , only : c_i, co2_readFlux_ocn, co2_readFlux_fuel
    use co2_cycle         , only : co2_transport, co2_time_interp_ocn, co2_time_interp_fuel
    use co2_cycle         , only : data_flux_ocn, data_flux_fuel
    use physconst         , only : mwco2
    use time_manager      , only : is_first_step, get_nstep

    ! input/output variabes
    type(ESMF_GridComp)               :: gcomp
    type(cam_in_t)    , intent(inout) :: cam_in(begchunk:endchunk)
    logical, optional , intent(in)    :: restart_init
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_State)   :: importState
    integer            :: i,n,c,g, num  ! indices
    integer            :: ncols         ! number of columns
    integer            :: nstep
    logical            :: overwrite_flds
    logical            :: exists
    logical            :: exists_fco2_ocn
    logical            :: exists_fco2_lnd
    character(len=128) :: fldname
    real(r8), pointer  :: fldptr2d(:,:)
    real(r8), pointer  :: fldptr1d(:)
    real(r8), pointer  :: fldptr_lat(:)
    real(r8), pointer  :: fldptr_lwup(:)
    real(r8), pointer  :: fldptr_avsdr(:)
    real(r8), pointer  :: fldptr_anidr(:)
    real(r8), pointer  :: fldptr_avsdf(:)
    real(r8), pointer  :: fldptr_anidf(:)
    real(r8), pointer  :: fldptr_tsurf(:)
    real(r8), pointer  :: fldptr_tocn(:)
    real(r8), pointer  :: fldptr_tref(:)
    real(r8), pointer  :: fldptr_qref(:)
    real(r8), pointer  :: fldptr_u10(:)
    real(r8), pointer  :: fldptr_snowhland(:)
    real(r8), pointer  :: fldptr_snowhice(:)
    real(r8), pointer  :: fldptr_ifrac(:)
    real(r8), pointer  :: fldptr_ofrac(:)
    real(r8), pointer  :: fldptr_lfrac(:)
    real(r8), pointer  :: fldptr_taux(:)
    real(r8), pointer  :: fldptr_tauy(:)
    real(r8), pointer  :: fldptr_sen(:)
    real(r8), pointer  :: fldptr_evap(:)
    logical, save      :: first_time = .true.
    character(len=*), parameter :: subname='(atm_import_export:import_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 10) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get import state
    call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! don't overwrite fields if invoked during the initialization phase
    ! of a 'continue' or 'branch' run type with data from .rs file
    overwrite_flds = .true.
    if (present(restart_init)) overwrite_flds = .not. restart_init

    !--------------------------
    ! Required atmosphere input fields
    !--------------------------

    if (overwrite_flds) then
       call state_getfldptr(importState, 'Faxx_taux', fldptr=fldptr_taux, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Faxx_tauy', fldptr=fldptr_tauy, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Faxx_sen' , fldptr=fldptr_sen, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Faxx_evap', fldptr=fldptr_evap, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%wsx(i)    = -fldptr_taux(g)
             cam_in(c)%wsy(i)    = -fldptr_tauy(g)
             cam_in(c)%shf(i)    = -fldptr_sen(g)
             cam_in(c)%cflx(i,1) = -fldptr_evap(g)
             g = g + 1
          end do
       end do
    end if  ! end of overwrite_flds

    call state_getfldptr(importState, 'Faxx_lat', fldptr=fldptr_lat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Faxx_lwup', fldptr=fldptr_lwup, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sx_avsdr', fldptr=fldptr_avsdr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sx_anidr', fldptr=fldptr_anidr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sx_avsdf', fldptr=fldptr_avsdf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sx_anidf', fldptr=fldptr_anidf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sx_t', fldptr=fldptr_tsurf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'So_t', fldptr=fldptr_tocn, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sl_snowh', fldptr=fldptr_snowhland, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Si_snowh', fldptr=fldptr_snowhice,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sx_tref', fldptr=fldptr_tref, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sx_qref', fldptr=fldptr_qref, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sx_u10', fldptr=fldptr_u10, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Si_ifrac', fldptr=fldptr_ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'So_ofrac', fldptr=fldptr_ofrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sl_lfrac', fldptr=fldptr_lfrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    g = 1
    do c = begchunk,endchunk
       do i = 1,get_ncols_p(c)
          cam_in(c)%lhf(i)       = -fldptr_lat(g)
          cam_in(c)%lwup(i)      = -fldptr_lwup(g)
          cam_in(c)%asdir(i)     =  fldptr_avsdr(g)
          cam_in(c)%aldir(i)     =  fldptr_anidr(g)
          cam_in(c)%asdif(i)     =  fldptr_avsdf(g)
          cam_in(c)%aldif(i)     =  fldptr_anidf(g)
          cam_in(c)%ts(i)        =  fldptr_tsurf(g)
          cam_in(c)%sst(i)       =  fldptr_tocn(g)
          cam_in(c)%tref(i)      =  fldptr_tref(g)
          cam_in(c)%qref(i)      =  fldptr_qref(g)
          cam_in(c)%u10(i)       =  fldptr_u10(g)
          cam_in(c)%snowhland(i) =  fldptr_snowhland(g)
          cam_in(c)%snowhice(i)  =  fldptr_snowhice(g)
          cam_in(c)%icefrac(i)   =  fldptr_ifrac(g)
          cam_in(c)%ocnfrac(i)   =  fldptr_ofrac(g)
          cam_in(c)%landfrac(i)  =  fldptr_lfrac(g)
          g = g + 1
       end do
    end do

    ! Optional fields

    call state_getfldptr(importState, 'Sl_ram1', fldptr=fldptr1d, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          if ( associated(cam_in(c)%ram1) ) then
             do i = 1, get_ncols_p(c)
                cam_in(c)%ram1(i) = fldptr1d(g)
                g = g + 1
             end do
          end if
       end do
    end if

    call state_getfldptr(importState, 'Sl_fv', fldptr=fldptr1d, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          if ( associated(cam_in(c)%fv) ) then
             do i = 1,get_ncols_p(c)
                cam_in(c)%fv(i) = fldptr1d(g)
                g = g + 1
             end do
          end if
       end do
    end if

    ! For CARMA - soil water from land
    call state_getfldptr(importState, 'Sl_soilw', fldptr=fldptr1d, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          if ( associated(cam_in(c)%soilw)) then
             do i = 1,get_ncols_p(c)
                cam_in(c)%soilw(i) = fldptr1d(g)
                g = g+1
             end do
          end if
       end do
    end if

    ! dry deposition fluxes from land
    call state_getfldptr(importState, 'Fall_flxdst', fldptr2d=fldptr2d, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          if ( associated(cam_in(c)%dstflx) ) then
             do i = 1,get_ncols_p(c)
                do n = 1, size(fldptr2d, dim=1)
                   cam_in(c)%dstflx(i,n) = fldptr2d(n,g)
                end do
                g = g + 1
             end do
          end if
       end do
    end if

    ! MEGAN VOC emis fluxes from land
    call state_getfldptr(importState, 'Fall_voc', fldptr2d=fldptr2d, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c=begchunk,endchunk
          if ( associated(cam_in(c)%meganflx) ) then
             do i = 1,get_ncols_p(c)
                do n = 1, size(fldptr2d, dim=1)
                   cam_in(c)%meganflx(i,n) = fldptr2d(n,g)
                end do
                g = g + 1
             end do
          end if
       end do
    end if

    ! fire emission fluxes from land
    call state_getfldptr(importState, 'Fall_fire', fldptr2d=fldptr2d, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          if ( associated(cam_in(c)%fireflx) .and. associated(cam_in(c)%fireztop) ) then
             do i = 1,get_ncols_p(c)
                do n = 1, size(fldptr2d, dim=1)
                   cam_in(c)%fireflx(i,n) = fldptr2d(n,g)
                end do
                g = g + 1
             end do
          end if
       end do
    end if
    call state_getfldptr(importState, 'Sl_fztop', fldptr=fldptr1d, exists=exists, rc=rc)
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%fireztop(i) = fldptr1d(g)
             g = g + 1
          end do
       end do
    end if

    ! dry dep velocities
    call state_getfldptr(importState, 'Sl_ddvel', fldptr2d=fldptr2d, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             do n = 1, size(fldptr2d, dim=1)
                cam_in(c)%depvel(i,n) = fldptr2d(n,g)
             end do
             g = g + 1
          end do
       end do
    end if

    ! fields needed to calculate water isotopes to ocean evaporation processes
    call state_getfldptr(importState,  'So_ustar', fldptr=fldptr1d, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%ustar(i) = fldptr1d(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState,  'So_re', fldptr=fldptr1d, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%re(i)= fldptr1d(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState,  'So_ssq', fldptr=fldptr1d, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%ssq(i) = fldptr1d(g)
             g = g + 1
          end do
       end do
    end if

    ! bgc scenarios
    call state_getfldptr(importState,  'Fall_fco2_lnd', fldptr=fldptr1d, exists=exists_fco2_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists_fco2_lnd) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%fco2_lnd(i) = -fldptr1d(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState,  'Faoo_fco2_ocn', fldptr=fldptr1d, exists=exists_fco2_ocn, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists_fco2_ocn) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%fco2_ocn(i) = -fldptr1d(g)
             g = g + 1
          end do
       end do
    else
       ! Consistency check
       if (co2_readFlux_ocn) then
          call shr_sys_abort(subname // ':: co2_readFlux_ocn and x2a_Faoo_fco2_ocn cannot both be active')
       end if
    end if
    call state_getfldptr(importState,  'Faoo_dms_ocn', fldptr=fldptr1d, exists=exists, rc=rc)
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%fdms(i) = -fldptr1d(g)
             g = g + 1
          end do
       end do
    end if

    ! -----------------------------------
    ! Get total co2 flux from components,
    ! -----------------------------------

    ! Note - co2_transport determines if cam_in(c)%cflx(i,c_i(1:4)) is allocated

    if (co2_transport() .and. overwrite_flds) then

       ! Interpolate in time for flux data read in
       if (co2_readFlux_ocn) then
          call co2_time_interp_ocn
       end if
       if (co2_readFlux_fuel) then
          call co2_time_interp_fuel
       end if

       ! from ocn : data read in or from coupler or zero
       ! from fuel: data read in or zero
       ! from lnd : through coupler or zero
       ! all co2 fluxes in unit kgCO2/m2/s

       do c=begchunk,endchunk
          do i=1, get_ncols_p(c)

             ! co2 flux from ocn
             if (exists_fco2_ocn /= 0) then
                cam_in(c)%cflx(i,c_i(1)) = cam_in(c)%fco2_ocn(i)
             else if (co2_readFlux_ocn) then
                ! convert from molesCO2/m2/s to kgCO2/m2/s
                cam_in(c)%cflx(i,c_i(1)) = &
                     -data_flux_ocn%co2flx(i,c)*(1._r8- cam_in(c)%landfrac(i))*mwco2*1.0e-3_r8
             else
                cam_in(c)%cflx(i,c_i(1)) = 0._r8
             end if

             ! co2 flux from fossil fuel
             if (co2_readFlux_fuel) then
                cam_in(c)%cflx(i,c_i(2)) = data_flux_fuel%co2flx(i,c)
             else
                cam_in(c)%cflx(i,c_i(2)) = 0._r8
             end if

             ! co2 flux from land (cpl already multiplies flux by land fraction)
             if (exists_fco2_lnd) then
                cam_in(c)%cflx(i,c_i(3)) = cam_in(c)%fco2_lnd(i)
             else
                cam_in(c)%cflx(i,c_i(3)) = 0._r8
             end if

             ! merged co2 flux
             cam_in(c)%cflx(i,c_i(4)) = cam_in(c)%cflx(i,c_i(1)) + cam_in(c)%cflx(i,c_i(2)) + cam_in(c)%cflx(i,c_i(3))
          end do
       end do
    end if

    ! if first step, determine longwave up flux from the surface temperature
    if (first_time) then
       if (is_first_step()) then
          do c=begchunk, endchunk
             do i=1, get_ncols_p(c)
                cam_in(c)%lwup(i) = shr_const_stebol*(cam_in(c)%ts(i)**4)
             end do
          end do
       end if
       first_time = .false.
    end if

    !-----------------------------------------------------------------
    ! Debug import
    !-----------------------------------------------------------------

    if (debug_import > 0 .and. masterproc .and. get_nstep()<5) then
       nstep = get_nstep()
       g=1
       do c=begchunk, endchunk
          do i=1, get_ncols_p(c)
             write(iulog,F01)'import: nstep, g, Faxx_tauy = ',nstep,g,-cam_in(c)%wsy(i)
             write(iulog,F01)'import: nstep, g, Faxx_taux = ',nstep,g,-cam_in(c)%wsx(i)
             write(iulog,F01)'import: nstep, g, Faxx_shf  = ',nstep,g,-cam_in(c)%shf(i)
             write(iulog,F01)'import: nstep, g, Faxx_lhf  = ',nstep,g,-cam_in(c)%lhf(i)
             write(iulog,F01)'import: nstep, g, Faxx_evap = ',nstep,g,-cam_in(c)%cflx(i,1)
             write(iulog,F01)'import: nstep, g, Faxa_lwup = ',nstep,g,-cam_in(c)%lwup(i)
             write(iulog,F01)'import: nstep, g, Sx_asdir  = ',nstep,g, cam_in(c)%asdir(i)
             write(iulog,F01)'import: nstep, g, Sx_aldir  = ',nstep,g, cam_in(c)%aldir(i)
             write(iulog,F01)'import: nstep, g, Sx_asdif  = ',nstep,g, cam_in(c)%asdif(i)
             write(iulog,F01)'import: nstep, g, Sx_aldif  = ',nstep,g, cam_in(c)%aldif(i)
             write(iulog,F01)'import: nstep, g, Sx_t      = ',nstep,g, cam_in(c)%ts(i)
             write(iulog,F01)'import: nstep, g, So_t      = ',nstep,g, cam_in(c)%sst(i)
             write(iulog,F01)'import: nstep, g, Sl_snowh  = ',nstep,g, cam_in(c)%snowhland(i)
             write(iulog,F01)'import: nstep, g, Si_snowh  = ',nstep,g, cam_in(c)%snowhice(i)
             write(iulog,F01)'import: nstep, g, Si_ifrac  = ',nstep,g, cam_in(c)%icefrac(i)
             write(iulog,F01)'import: nstep, g, So_ofrac  = ',nstep,g, cam_in(c)%ocnfrac(i)
             write(iulog,F01)'import: nstep, g, Sl_lfrac  = ',nstep,g, cam_in(c)%landfrac(i)
             write(iulog,F01)'import: nstep, g, Sx_tref   = ',nstep,g, cam_in(c)%tref(i)
             write(iulog,F01)'import: nstep, g, Sx_qref   = ',nstep,g, cam_in(c)%qref(i)
             write(iulog,F01)'import: nstep, g, Sx_qu10   = ',nstep,g, cam_in(c)%u10(i)
             g = g + 1
          end do
       end do
    end if

    if (dbug_flag > 10) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine import_fields

  !===============================================================================

  subroutine export_fields( gcomp, cam_out, rc)

    ! -----------------------------------------------------
    ! Set field pointers in export set
    ! Copy from chunk array data structure into state fldptr
    ! -----------------------------------------------------

    use camsrfexch   , only : cam_out_t
    use phys_grid    , only : get_ncols_p
    use ppgrid       , only : begchunk, endchunk
    use time_manager , only : is_first_step, get_nstep
    use spmd_utils   , only : masterproc

    !-------------------------------
    ! Pack the export state
    !-------------------------------

    ! input/output variables
    type(ESMF_GridComp)           :: gcomp
    type(cam_out_t) , intent(in)  :: cam_out(begchunk:endchunk)
    integer         , intent(out) :: rc

    ! local variables
    type(ESMF_State)  :: exportState
    integer           :: i,m,c,n,g  ! indices
    integer           :: ncols      ! Number of columns
    integer           :: nstep
    logical           :: exists
    ! 2d pointers
    real(r8), pointer :: fldptr_ndep(:,:)
    real(r8), pointer :: fldptr_bcph(:,:)  , fldptr_ocph(:,:)
    real(r8), pointer :: fldptr_dstwet(:,:), fldptr_dstdry(:,:) 
    ! 1d pointers
    real(r8), pointer :: fldptr_soll(:)    , fldptr_sols(:)
    real(r8), pointer :: fldptr_solld(:)   , fldptr_solsd(:)
    real(r8), pointer :: fldptr_snowc(:)   , fldptr_snowl(:)
    real(r8), pointer :: fldptr_rainc(:)   , fldptr_rainl(:)
    real(r8), pointer :: fldptr_lwdn(:)    , fldptr_swnet(:)
    real(r8), pointer :: fldptr_topo(:)    , fldptr_zbot(:)
    real(r8), pointer :: fldptr_ubot(:)    , fldptr_vbot(:)
    real(r8), pointer :: fldptr_pbot(:)    , fldptr_tbot(:)
    real(r8), pointer :: fldptr_shum(:)    , fldptr_dens(:) 
    real(r8), pointer :: fldptr_ptem(:)    , fldptr_pslv(:) 
    real(r8), pointer :: fldptr_co2prog(:) , fldptr_co2diag(:)
    character(len=*), parameter :: subname='(atm_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get export state
    call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! required export state variables
    call state_getfldptr(exportState, 'Sa_topo', fldptr=fldptr_topo, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Sa_z'   , fldptr=fldptr_zbot, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Sa_u'   , fldptr=fldptr_ubot, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Sa_v'   , fldptr=fldptr_vbot, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Sa_tbot', fldptr=fldptr_tbot, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Sa_pbot', fldptr=fldptr_pbot, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Sa_shum', fldptr=fldptr_shum, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Sa_dens', fldptr=fldptr_dens, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Sa_ptem', fldptr=fldptr_ptem, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Sa_pslv', fldptr=fldptr_pslv, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    g = 1
    do c = begchunk,endchunk
       do i = 1,get_ncols_p(c)
          fldptr_topo(g) = cam_out(c)%topo(i)
          fldptr_zbot(g) = cam_out(c)%zbot(i)
          fldptr_ubot(g) = cam_out(c)%ubot(i)
          fldptr_vbot(g) = cam_out(c)%vbot(i)
          fldptr_pbot(g) = cam_out(c)%pbot(i)
          fldptr_tbot(g) = cam_out(c)%tbot(i)
          fldptr_shum(g) = cam_out(c)%qbot(i,1)
          fldptr_dens(g) = cam_out(c)%rho(i)
          fldptr_ptem(g) = cam_out(c)%thbot(i)
          fldptr_pslv(g) = cam_out(c)%psl(i)
          g = g + 1
       end do
    end do

    ! required export flux variables
    call state_getfldptr(exportState, 'Faxa_swnet', fldptr=fldptr_swnet,  rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Faxa_lwdn' , fldptr=fldptr_lwdn , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Faxa_rainc', fldptr=fldptr_rainc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Faxa_rainl', fldptr=fldptr_rainl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Faxa_snowc', fldptr=fldptr_snowc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Faxa_snowl', fldptr=fldptr_snowl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Faxa_swndr', fldptr=fldptr_soll, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Faxa_swvdr', fldptr=fldptr_sols, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Faxa_swndf', fldptr=fldptr_solld, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Faxa_swvdf', fldptr=fldptr_solsd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    g = 1
    do c = begchunk,endchunk
       do i = 1,get_ncols_p(c)
          fldptr_lwdn(g)  = cam_out(c)%flwds(i)
          fldptr_swnet(g) = cam_out(c)%netsw(i)
          fldptr_snowc(g) = cam_out(c)%precsc(i)*1000._r8
          fldptr_snowl(g) = cam_out(c)%precsl(i)*1000._r8
          fldptr_rainc(g) = (cam_out(c)%precc(i)-cam_out(c)%precsc(i))*1000._r8
          fldptr_rainl(g) = (cam_out(c)%precl(i)-cam_out(c)%precsl(i))*1000._r8
          fldptr_soll(g)  = cam_out(c)%soll(i)
          fldptr_sols(g)  = cam_out(c)%sols(i)
          fldptr_solld(g) = cam_out(c)%solld(i)
          fldptr_solsd(g) = cam_out(c)%solsd(i)
          g = g + 1
       end do
    end do

    ! aerosol deposition fluxes
    call state_getfldptr(exportState, 'Faxa_bcph', fldptr2d=fldptr_bcph, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Faxa_ocph', fldptr2d=fldptr_ocph, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Faxa_dstdry', fldptr2d=fldptr_dstdry, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Faxa_dstwet', fldptr2d=fldptr_dstwet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! (1) => bcphidry, (2) => bcphodry, (3) => bcphiwet
    ! (1) => ocphidry, (2) => ocphodry, (3) => ocphiwet
    g = 1
    do c = begchunk,endchunk
       do i = 1,get_ncols_p(c)
          fldptr_bcph(1,g)   = cam_out(c)%bcphidry(i)
          fldptr_bcph(2,g)   = cam_out(c)%bcphodry(i)
          fldptr_bcph(3,g)   = cam_out(c)%bcphiwet(i)
          fldptr_ocph(1,g)   = cam_out(c)%ocphidry(i)
          fldptr_ocph(2,g)   = cam_out(c)%ocphodry(i)
          fldptr_ocph(3,g)   = cam_out(c)%ocphiwet(i)
          fldptr_dstdry(1,g) = cam_out(c)%dstdry1(i)
          fldptr_dstdry(2,g) = cam_out(c)%dstdry2(i)
          fldptr_dstdry(3,g) = cam_out(c)%dstdry3(i)
          fldptr_dstdry(4,g) = cam_out(c)%dstdry4(i)
          fldptr_dstwet(1,g) = cam_out(c)%dstwet1(i)
          fldptr_dstwet(2,g) = cam_out(c)%dstwet2(i)
          fldptr_dstwet(3,g) = cam_out(c)%dstwet3(i)
          fldptr_dstwet(4,g) = cam_out(c)%dstwet4(i)
          g = g + 1
       end do
    end do

    call state_getfldptr(exportState, 'Sa_co2prog', fldptr=fldptr_co2prog, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr_co2prog(g) = cam_out(c)%co2prog(i) ! atm prognostic co2
             g = g + 1
          end do
       end do
    end if

    call state_getfldptr(exportState, 'Sa_co2diag', fldptr=fldptr_co2diag, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr_co2diag(g) = cam_out(c)%co2diag(i) ! atm diagnostic co2
             g = g + 1
          end do
       end do
    end if

    call state_getfldptr(exportState, 'Faxa_ndep', fldptr2d=fldptr_ndep, exists=exists, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       ! (1) => nhx, (2) => noy
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr_ndep(1,g) = cam_out(c)%nhx_nitrogen_flx(i)
             fldptr_ndep(2,g) = cam_out(c)%noy_nitrogen_flx(i)
             g = g + 1
          end do
       end do
    end if

    !-----------------------------------------------------------------
    ! Debug export
    !-----------------------------------------------------------------

    if (debug_export > 0 .and. masterproc .and. get_nstep()<5) then
       nstep = get_nstep()
       g=1
       do c=begchunk, endchunk
          do i=1, get_ncols_p(c)
             write(iulog,F01)'export: nstep, g, Sa_z          = ',nstep,g,cam_out(c)%zbot(i)
             write(iulog,F01)'export: nstep, g, Sa_topo       = ',nstep,g,cam_out(c)%topo(i)
             write(iulog,F01)'export: nstep, g, Sa_u          = ',nstep,g,cam_out(c)%ubot(i)
             write(iulog,F01)'export: nstep, g, Sa_v          = ',nstep,g,cam_out(c)%vbot(i)
             write(iulog,F01)'export: nstep, g, Sa_tbot       = ',nstep,g,cam_out(c)%tbot(i)
             write(iulog,F01)'export: nstep, g, Sa_ptem       = ',nstep,g,cam_out(c)%thbot(i)
             write(iulog,F01)'export: nstep, g, Sa_pbot       = ',nstep,g,cam_out(c)%pbot(i)
             write(iulog,F01)'export: nstep, g, Sa_shum       = ',nstep,g,cam_out(c)%qbot(i,1)
             write(iulog,F01)'export: nstep, g, Sa_dens       = ',nstep,g,cam_out(c)%rho(i)
             write(iulog,F01)'export: nstep, g, Faxa_swnet    = ',nstep,g,cam_out(c)%netsw(i)
             write(iulog,F01)'export: nstep, g, Faxa_lwdn     = ',nstep,g,cam_out(c)%flwds(i)
             write(iulog,F01)'export: nstep, g, Faxa_rainc    = ',nstep,g,(cam_out(c)%precc(i)-cam_out(c)%precsc(i))*1000._r8
             write(iulog,F01)'export: nstep, g, Faxa_rainl    = ',nstep,g,(cam_out(c)%precl(i)-cam_out(c)%precsl(i))*1000._r8
             write(iulog,F01)'export: nstep, g, Faxa_snowc    = ',nstep,g,cam_out(c)%precsc(i)*1000._r8
             write(iulog,F01)'export: nstep, g, Faxa_snowl    = ',nstep,g,cam_out(c)%precsl(i)*1000._r8
             write(iulog,F01)'export: nstep, g, Faxa_swndr    = ',nstep,g,cam_out(c)%soll(i)
             write(iulog,F01)'export: nstep, g, Faxa_swvdr    = ',nstep,g,cam_out(c)%sols(i)
             write(iulog,F01)'export: nstep, g, Faxa_swndf    = ',nstep,g,cam_out(c)%solld(i)
             write(iulog,F01)'export: nstep, g, Faxa_swvdf    = ',nstep,g,cam_out(c)%solsd(i)
             write(iulog,F01)'export: nstep, g, Faxa_bcphidry = ',nstep,g,cam_out(c)%bcphidry(i)
             write(iulog,F01)'export: nstep, g, Faxa_bcphodry = ',nstep,g,cam_out(c)%bcphodry(i)
             write(iulog,F01)'export: nstep, g, Faxa_bcphiwet = ',nstep,g,cam_out(c)%bcphiwet(i)
             write(iulog,F01)'export: nstep, g, Faxa_ocphidry = ',nstep,g,cam_out(c)%ocphidry(i)
             write(iulog,F01)'export: nstep, g, Faxa_ocphodry = ',nstep,g,cam_out(c)%ocphodry(i)
             write(iulog,F01)'export: nstep, g, Faxa_ocphidry = ',nstep,g,cam_out(c)%ocphiwet(i)
             write(iulog,F01)'export: nstep, g, Faxa_dstwet1  = ',nstep,g,cam_out(c)%dstwet1(i)
             write(iulog,F01)'export: nstep, g, Faxa_dstwet1  = ',nstep,g,cam_out(c)%dstdry1(i)
             write(iulog,F01)'export: nstep, g, Faxa_dstwet2  = ',nstep,g,cam_out(c)%dstwet2(i)
             write(iulog,F01)'export: nstep, g, Faxa_dstwet2  = ',nstep,g,cam_out(c)%dstdry2(i)
             write(iulog,F01)'export: nstep, g, Faxa_dstwet3  = ',nstep,g,cam_out(c)%dstwet3(i)
             write(iulog,F01)'export: nstep, g, Faxa_dstwet3  = ',nstep,g,cam_out(c)%dstdry3(i)
             write(iulog,F01)'export: nstep, g, Faxa_dstwet4  = ',nstep,g,cam_out(c)%dstwet4(i)
             write(iulog,F01)'export: nstep, g, Faxa_dstwet4  = ',nstep,g,cam_out(c)%dstdry4(i)
             write(iulog,F01)'export: nstep, g, Sa_co2prog    = ',nstep,g,cam_out(c)%co2prog(i)
             write(iulog,F01)'export: nstep, g, Sa_co2diag    = ',nstep,g,cam_out(c)%co2diag(i)
             g = g + 1
          end do
       end do
    end if

  end subroutine export_fields

  !===============================================================================

  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/otuput variables
    integer            , intent(inout) :: num
    type(fldlist_type) , intent(inout) :: fldlist(:)
    character(len=*)   , intent(in)    :: stdname
    integer, optional  , intent(in)    :: ungridded_lbound
    integer, optional  , intent(in)    :: ungridded_ubound

    ! local variables
    character(len=*), parameter :: subname='(atm_import_export:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       return
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fldlist_add

  !===============================================================================

  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    type(ESMF_State)    , intent(inout) :: state
    type(fldlist_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(atm_import_export:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                     ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                     gridToFieldMap=(/2/), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
             end if
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------

      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      ! input/output variables
      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(atm_import_export:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================

  subroutine state_getfldptr(State, fldname, fldptr, fldptr2d, exists, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    ! input/output variables
    type(ESMF_State)  , intent(in)    :: State
    character(len=*)  , intent(in)    :: fldname
    real(R8), optional, pointer       :: fldptr(:)
    real(R8), optional, pointer       :: fldptr2d(:,:)
    logical , optional, intent(out)   :: exists
    integer           , intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_StateItem_Flag)   :: itemFlag
    type(ESMF_Field)            :: lfield
    type(ESMF_Mesh)             :: lmesh
    integer                     :: nnodes, nelements
    character(len=*), parameter :: subname='(atm_import_export:state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (present(exists)) then
       ! if field exists then create output array - else do nothing
       if (itemflag == ESMF_STATEITEM_NOTFOUND) then
          exists = .false.
          RETURN
       else
          exists = .true.
       end if
    else
       if (itemflag == ESMF_STATEITEM_NOTFOUND) then
          call shr_sys_abort('variable '//trim(fldname)//' must be present ')
       end if
    end if
       
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, status=status, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    else
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (nnodes == 0 .and. nelements == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if

       if (present(fldptr)) then
          call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else if (present(fldptr2d)) then
          call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    endif  ! status

    if (dbug_flag > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif

  end subroutine state_getfldptr

end module atm_import_export
