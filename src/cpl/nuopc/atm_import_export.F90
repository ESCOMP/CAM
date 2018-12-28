module atm_import_export

  use ESMF                  , only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF                  , only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF                  , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF                  , only : operator(/=), operator(==)
  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_IsConnected
  use NUOPC_Model           , only : NUOPC_ModelGet
  use shr_kind_mod          , only : r8 => shr_kind_r8, cx=>shr_kind_cx, cxx=>shr_kind_cxx
  use shr_string_mod        , only : shr_string_listGetName, shr_string_listGetNum
  use shr_sys_mod           , only : shr_sys_abort
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkerr
  use shr_nuopc_scalars_mod , only : flds_scalar_name, flds_scalar_num
  use seq_drydep_mod        , only : seq_drydep_readnl, seq_drydep_init, n_drydep
  use shr_megan_mod         , only : shr_megan_readnl, shr_megan_mechcomps_n
  use shr_fire_emis_mod     , only : shr_fire_emis_readnl, shr_fire_emis_mechcomps_n, shr_fire_emis_ztop_token
  use shr_carma_mod         , only : shr_carma_readnl
  use shr_ndep_mod          , only : shr_ndep_readnl
  use cam_logfile           , only : iulog
  use srf_field_check       , only : set_active_Sl_ram1
  use srf_field_check       , only : set_active_Sl_fv
  use srf_field_check       , only : set_active_Sl_soilw
  use srf_field_check       , only : set_active_Fall_flxdst1
  use srf_field_check       , only : set_active_Fall_flxvoc
  use srf_field_check       , only : set_active_Fall_flxfire
  use srf_field_check       , only : set_active_Fall_fco2_lnd
  use srf_field_check       , only : set_active_Faoo_fco2_ocn
  use srf_field_check       , only : set_active_Faxa_nhx
  use srf_field_check       , only : set_active_Faxa_noy

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
  end type fldlist_type

  integer             , parameter         :: fldsMax = 100
  integer             , public, protected :: fldsToAtm_num = 0
  integer             , public, protected :: fldsFrAtm_num = 0
  type (fldlist_type) , public, protected :: fldsToAtm(fldsMax)
  type (fldlist_type) , public, protected :: fldsFrAtm(fldsMax)

  character(len=cxx)     :: drydep_fields       ! List of dry-deposition fields
  character(len=cxx)     :: megan_voc_fields    ! List of MEGAN VOC emission fields
  character(len=cxx)     :: fire_emis_fields    ! List of fire emission fields
  character(len=cx)      :: carma_fields        ! List of CARMA fields from lnd->atm
  character(len=cx)      :: ndep_fields         ! List of nitrogen deposition fields from atm->lnd/ocn

  integer                :: glc_nec
  integer, parameter     :: dbug = 10
  integer, parameter     :: debug_import = 0 ! internal debug level
  integer, parameter     :: debug_export = 0 ! internal debug level
  character(*),parameter :: F01 = "('(cam_import_export) ',a,i8,2x,i8,2x,d21.14)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, rc)

    use spmd_utils, only : masterproc, mpicom

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)       :: importState
    type(ESMF_State)       :: exportState
    character(ESMF_MAXSTR) :: stdname
    character(ESMF_MAXSTR) :: cvalue
    character(len=2)       :: nec_str
    integer                :: dbrc
    integer                :: n, num
    logical                :: flds_co2a      ! use case
    logical                :: flds_co2b      ! use case
    logical                :: flds_co2c      ! use case
    integer                :: ndep_nflds, megan_nflds, emis_nflds
    integer                :: seq_drydep_nflds
    integer                :: dbug_flag = 6
    character(len=128)     :: fldname
    character(len=*), parameter :: subname='(atm_import_export:advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! determine necessary toggles for below
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2a
    call ESMF_LogWrite('flds_co2a = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2b
    call ESMF_LogWrite('flds_co2b = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2c
    call ESMF_LogWrite('flds_co2c = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    !--------------------------------
    ! Export fields
    !--------------------------------
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' export fields', ESMF_LOGMSG_INFO, rc=dbrc)
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
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_bcphidry' )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_bcphodry' )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_bcphiwet' )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_ocphidry' )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_ocphodry' )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_ocphiwet' )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet1'  )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet2'  )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet3'  )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet4'  )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry1'  )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry2'  )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry3'  )
    call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry4'  )

    if (dbug_flag > 5) call ESMF_LogWrite(subname//' export fields co2', ESMF_LOGMSG_INFO, rc=dbrc)

    ! co2 fields
    if (flds_co2a .or. flds_co2b .or. flds_co2c) then
       call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_co2prog' )
       call fldlist_add(fldsFrAtm_num, fldsFrAtm, 'Sa_co2diag' )
    end if


    ! nitrogen deposition from atm
    call shr_ndep_readnl("drv_flds_in", ndep_fields, ndep_nflds)
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' export fields ndep', ESMF_LOGMSG_INFO, rc=dbrc)
    if (ndep_nflds > 0) then
       do n = 1, ndep_nflds
          call shr_string_listGetName(ndep_fields, n, fldname)
          call fldlist_add(fldsFrAtm_num, fldsFrAtm, trim(fldname))
          call set_active_Faxa_nhx(.true.)
          call set_active_Faxa_noy(.true.)
       end do
    end if

    ! Now advertise above export fields
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' advertise export fields ', ESMF_LOGMSG_INFO, rc=dbrc)
    do n = 1,fldsFrAtm_num
       call NUOPC_Advertise(exportState, standardName=fldsFrAtm(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !-----------------
    ! Import fields
    !-----------------
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' Import Fields', ESMF_LOGMSG_INFO, rc=dbrc)
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
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Fall_flxdst1')
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Fall_flxdst2')
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Fall_flxdst3')
    call fldlist_add(fldsToAtm_num, fldsToAtm, 'Fall_flxdst4')

    ! co2 fields from land and ocean
    if (flds_co2b .or. flds_co2c) then
       call fldlist_add(fldsToAtm_num, fldsToAtm, 'Fall_fco2_lnd')
       call set_active_Fall_fco2_lnd(.true.)
    end if
    if (flds_co2c) then
       call fldlist_add(fldsToAtm_num, fldsToAtm, 'Faoo_fco2_ocn')
       call set_active_Faoo_fco2_ocn(.true.)
    end if

    ! Dry Deposition fluxes from land - ALSO initialize drydep here
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' drydep readnl', ESMF_LOGMSG_INFO, rc=dbrc)
    call seq_drydep_readnl("drv_flds_in", drydep_fields, seq_drydep_nflds)
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' Import Fields (drydep)', ESMF_LOGMSG_INFO, rc=dbrc)
    if (seq_drydep_nflds > 0) then
       do n = 1, n_drydep
          call shr_string_listGetName(drydep_fields, n, fldname)
          call fldlist_add(fldsToAtm_num, fldsToAtm, trim(fldname))
       end do
       call seq_drydep_init( )
    endif
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' Import Fields megan', ESMF_LOGMSG_INFO, rc=dbrc)

    ! MEGAN VOC emissions fluxes from land
    call shr_megan_readnl('drv_flds_in', megan_voc_fields, megan_nflds)
    do n = 1, shr_megan_mechcomps_n
       call shr_string_listGetName(megan_voc_fields, n, fldname)
       call fldlist_add(fldsToAtm_num, fldsToAtm, trim(fldname))
       if (n == 1) call set_active_Fall_flxvoc(.true.)
    end do

    if (dbug_flag > 5) call ESMF_LogWrite(subname//' Import Fields fire', ESMF_LOGMSG_INFO, rc=dbrc)
    ! Fire emissions fluxes from land
    call shr_fire_emis_readnl('drv_flds_in', fire_emis_fields, emis_nflds)
    if (shr_fire_emis_mechcomps_n > 0) then
       do n = 1, shr_fire_emis_mechcomps_n
          call shr_string_listGetName(fire_emis_fields, n, fldname)
          call fldlist_add(fldsToAtm_num, fldsToAtm, trim(fldname))
       end do
       call fldlist_add(fldsToAtm_num, fldsToAtm, trim(shr_fire_emis_ztop_token))
       if (n == 1) call set_active_Fall_flxfire(.true.)
    end if
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' Import Fields carma', ESMF_LOGMSG_INFO, rc=dbrc)

    ! CARMA volumetric soil water from land
    call shr_carma_readnl('drv_flds_in', carma_fields)
    if (carma_fields /= ' ') then
       call fldlist_add(fldsToAtm_num, fldsToAtm, 'Sl_soilw') ! optional for carma
       call set_active_Sl_soilw(.true.) ! check for carma
    end if

    ! Now advertise above import fields
    do n = 1,fldsToAtm_num
       call NUOPC_Advertise(importState, standardName=fldsToAtm(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo
    if (dbug_flag > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine advertise_fields

  !===============================================================================

  subroutine realize_fields(gcomp, Emesh, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_Mesh)      :: Emesh
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    character(len=*), parameter :: subname='(lnd_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrAtm, &
         numflds=fldsFrAtm_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':clmExport',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToAtm, &
         numflds=fldsToAtm_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':clmImport',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine realize_fields

  !===============================================================================

  subroutine import_fields( gcomp, cam_in, restart_init, rc)

    !---------------------------------------------------------------------------
    ! Convert the input data from the mediator to the land model
    !---------------------------------------------------------------------------

    use camsrfexch        , only : cam_in_t
    use phys_grid         , only : get_ncols_p
    use ppgrid            , only : begchunk, endchunk
    use shr_const_mod     , only : shr_const_stebol
    use co2_cycle         , only : c_i, co2_readFlux_ocn, co2_readFlux_fuel
    use co2_cycle         , only : co2_transport, co2_time_interp_ocn, co2_time_interp_fuel
    use co2_cycle         , only : data_flux_ocn, data_flux_fuel
    use physconst         , only : mwco2
    use spmd_utils        , only : masterproc, mpicom
    use time_manager      , only : is_first_step, get_nstep
    use seq_drydep_mod    , only : seq_drydep_readnl, seq_drydep_init, n_drydep
    use shr_megan_mod     , only : shr_megan_readnl, shr_megan_mechcomps_n
    use shr_fire_emis_mod , only : shr_fire_emis_readnl, shr_fire_emis_mechcomps_n, shr_fire_emis_ztop_token
    use shr_carma_mod     , only : shr_carma_readnl
    use shr_ndep_mod      , only : shr_ndep_readnl

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
    real(r8), pointer  :: fldptr(:)
    logical            :: exists
    logical            :: exists_fco2_ocn
    logical            :: exists_fco2_lnd
    integer            :: dbrc
    character(len=128) :: fldname
    logical, save      :: first_time = .true.
    character(len=*), parameter :: subname='(atm_import_export:import_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! Get import state
    call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! don't overwrite fields if invoked during the initialization phase
    ! of a 'continue' or 'branch' run type with data from .rs file
    overwrite_flds = .true.
    if (present(restart_init)) overwrite_flds = .not. restart_init

    !--------------------------
    ! Required atmosphere input fields
    !--------------------------

    if (overwrite_flds) then
       call state_getfldptr(importState, 'Faxx_taux', fldptr, exists, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (exists) then
          g = 1
          do c = begchunk,endchunk
             do i = 1,get_ncols_p(c)
                cam_in(c)%wsx(i) = -fldptr(g)
                g = g + 1
             end do
          end do
       end if
       call state_getfldptr(importState, 'Faxx_tauy', fldptr, exists, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (exists) then
          g = 1
          do c = begchunk,endchunk
             do i = 1,get_ncols_p(c)
                cam_in(c)%wsy(i) = -fldptr(g)
                g = g + 1
             end do
          end do
       end if
       call state_getfldptr(importState, 'Faxx_sen' , fldptr, exists, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (exists) then
          g = 1
          do c = begchunk,endchunk
             do i = 1,get_ncols_p(c)
                cam_in(c)%shf(i) = -fldptr(g)
                g = g + 1
             end do
          end do
       end if
       call state_getfldptr(importState, 'Faxx_evap', fldptr, exists, rc=rc)
       if (exists) then
          g = 1
          do c = begchunk,endchunk
             do i = 1,get_ncols_p(c)
                cam_in(c)%cflx(i,1) = -fldptr(g)
                g = g + 1
             end do
          end do
       end if
    end if  ! end of overwrite_flds

    call state_getfldptr(importState, 'Faxx_lat', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%lhf(i) = -fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Faxx_lwup', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%lwup(i) = -fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Sx_avsdr', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%asdir(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Sx_anidr', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%aldir(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Sx_avsdf', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%asdif(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Sx_anidf', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%aldif(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Sx_t', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%ts(i)  = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'So_t', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%sst(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Sl_snowh', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%snowhland(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Si_snowh', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%snowhice(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Sx_tref', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%tref(i) =  fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Sx_qref', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%qref(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Sx_u10', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%u10(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Si_ifrac', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%icefrac(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'So_ofrac', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%ocnfrac(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Sl_lfrac', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%landfrac(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState, 'Sl_ram1', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          if ( associated(cam_in(c)%ram1) ) then
             do i = 1, get_ncols_p(c)
                cam_in(c)%ram1(i) = fldptr(g)
                g = g + 1
             end do
          end if
       end do
    end if
    call state_getfldptr(importState, 'Sl_fv', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          if ( associated(cam_in(c)%fv) ) then
             do i = 1,get_ncols_p(c)
                cam_in(c)%fv(i) = fldptr(g)
                g = g + 1
             end do
          end if
       end do
    end if
    call state_getfldptr(importState, 'Sl_soilw', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          if ( associated(cam_in(c)%soilw)) then
             do i = 1,get_ncols_p(c)
                cam_in(c)%soilw(i) = fldptr(g)
                g = g+1
             end do
          end if
       end do
    end if

    ! Dry Deposition fluxes from land - ALSO initialize drydep here
    call state_getfldptr(importState, 'Fall_flxdst1', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          if ( associated(cam_in(c)%dstflx) ) then
             do i = 1,get_ncols_p(c)
                cam_in(c)%dstflx(i,1) = fldptr(g)
                g = g + 1
             end do
          end if
       end do
    end if
    call state_getfldptr(importState, 'Fall_flxdst2', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          if ( associated(cam_in(c)%dstflx) ) then
             do i = 1,get_ncols_p(c)
                cam_in(c)%dstflx(i,2) = fldptr(g)
                g = g + 1
             end do
          end if
       end do
    end if
    call state_getfldptr(importState, 'Fall_flxdst3', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          if ( associated(cam_in(c)%dstflx) ) then
             do i = 1,get_ncols_p(c)
                cam_in(c)%dstflx(i,3) = fldptr(g)
                g = g + 1
             end do
          end if
       end do
    end if
    call state_getfldptr(importState, 'Fall_flxdst4', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          if ( associated(cam_in(c)%dstflx) ) then
             do i = 1,get_ncols_p(c)
                cam_in(c)%dstflx(i,4) = fldptr(g)
                g = g + 1
             end do
          end if
       end do
    end if
    call seq_drydep_init( )

    ! MEGAN VOC emis fluxes from land
    do num = 1, shr_megan_mechcomps_n
       call shr_string_listGetName(megan_voc_fields, num, fldname)
       call state_getfldptr(importState, trim(fldname), fldptr, exists, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (exists) then
          g = 1
          do c=begchunk,endchunk
             if ( associated(cam_in(c)%meganflx) ) then
                do i = 1,get_ncols_p(c)
                   cam_in(c)%meganflx(i,n) = fldptr(g)
                   g = g + 1
                end do
             end if
          end do
       end if
    end do

    ! Fire emission fluxes from land
    do num = 1, shr_string_listGetNum(fire_emis_fields)
       call shr_string_listGetName(fire_emis_fields, num, fldname)
       call state_getfldptr(importState, trim(fldname), fldptr, exists, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (exists) then
          g = 1
          do c = begchunk,endchunk
             if ( associated(cam_in(c)%fireflx) .and. associated(cam_in(c)%fireztop) ) then
                do i = 1,get_ncols_p(c)
                   cam_in(c)%fireflx(i,n) = fldptr(g)
                   g = g + 1
                end do
             end if
          end do
       end if
    end do
    call state_getfldptr(importState,  trim(shr_fire_emis_ztop_token), fldptr, exists, rc=rc)
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%fireztop(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if

    ! dry dep velocities
    do num = 1, n_drydep
       call shr_string_listGetName(drydep_fields, num, fldname)
       call state_getfldptr(importState, trim(fldname), fldptr, exists, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (exists) then
          g = 1
          do c = begchunk,endchunk
             do i = 1,get_ncols_p(c)
                cam_in(c)%depvel(i,n) = fldptr(g)
                g = g + 1
             end do
          end do
       end if
    end do

    ! fields needed to calculate water isotopes to ocean evaporation processes
    call state_getfldptr(importState,  'So_ustar', fldptr, exists, rc=rc)
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%ustar(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState,  'So_re', fldptr, exists, rc=rc)
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%re(i)= fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState,  'So_ssq', fldptr, exists, rc=rc)
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%ssq(i) = fldptr(g)
             g = g + 1
          end do
       end do
    end if

    ! bgc scenarios
    call state_getfldptr(importState,  'Fall_fco2_lnd', fldptr, exists_fco2_lnd, rc=rc)
    if (exists_fco2_lnd) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%fco2_lnd(i) = -fldptr(g)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(importState,  'Faoo_fco2_ocn', fldptr, exists_fco2_ocn, rc=rc)
    if (exists_fco2_ocn) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%fco2_ocn(i) = -fldptr(g)
             g = g + 1
          end do
       end do
    else
       ! Consistency check
       if (co2_readFlux_ocn) then
          call shr_sys_abort(subname // ':: co2_readFlux_ocn and x2a_Faoo_fco2_ocn cannot both be active')
       end if
    end if
    call state_getfldptr(importState,  'Faoo_dms_ocn', fldptr, exists, rc=rc)
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             cam_in(c)%fdms(i) = -fldptr(g)
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

  end subroutine import_fields

  !===============================================================================

  subroutine export_fields( gcomp, cam_out, rc)

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
    integer           :: avsize, avnat
    integer           :: i,m,c,n,g, num   ! indices
    integer           :: ncols            ! Number of columns
    integer           :: nstep
    real(r8), pointer :: fldptr(:)
    logical           :: exists
    character(len=*), parameter :: subname='(lnd_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get export state
    call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set field pointers in export set

    ! Copy from chunk array data structure into state fldptr

    call state_getfldptr(exportState, 'Sa_pslv', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%psl(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Sa_z', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%zbot(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Sa_topo', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%topo(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Sa_u', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%ubot(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Sa_v', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%vbot(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Sa_tbot', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%tbot(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Sa_ptem', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%thbot(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Sa_pbot', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%pbot(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Sa_shum', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%qbot(i,1)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Sa_dens', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%rho(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_swnet', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%netsw(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_lwdn', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%flwds(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_rainc', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = (cam_out(c)%precc(i)-cam_out(c)%precsc(i))*1000._r8
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_rainl', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = (cam_out(c)%precl(i)-cam_out(c)%precsl(i))*1000._r8
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_snowc', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%precsc(i)*1000._r8
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_snowl', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%precsl(i)*1000._r8
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_swndr', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%soll(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_swvdr', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%sols(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_swndf', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%solld(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_swvdf', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%solsd(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_bcphidry', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%bcphidry(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_bcphodry', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%bcphodry(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_bcphiwet', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%bcphiwet(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_ocphidry', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%ocphidry(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_ocphodry', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%ocphodry(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_ocphiwet', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%ocphiwet(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_dstwet1', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%dstwet1(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_dstdry1', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%dstdry1(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_dstwet2', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%dstwet2(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_dstdry2', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%dstdry2(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_dstwet3', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%dstwet3(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_dstdry3', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%dstdry3(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_dstwet4', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%dstwet4(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_dstdry4', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%dstdry4(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Sa_co2prog', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%co2prog(i) ! atm prognostic co2
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Sa_co2diag', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%co2diag(i) ! atm diagnostic co2
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_nhx', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%nhx_nitrogen_flx(i)
             g = g + 1
          end do
       end do
    end if
    call state_getfldptr(exportState, 'Faxa_noy', fldptr, exists, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (exists) then
       g = 1
       do c = begchunk,endchunk
          do i = 1,get_ncols_p(c)
             fldptr(g) = cam_out(c)%noy_nitrogen_flx(i)
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

  subroutine fldlist_add(num, fldlist, stdname)
    integer,                    intent(inout) :: num
    type(fldlist_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname

    ! local variables
    integer :: rc
    integer :: dbrc
    character(len=*), parameter :: subname='(dshr_nuopc_mod:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
       return
    endif
    fldlist(num)%stdname = trim(stdname)

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
    integer                :: dbrc
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(dshr_nuopc_mod:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the field
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
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

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(dshr_nuopc_mod:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================

  subroutine state_getfldptr(State, fldname, fldptr, exists, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------
    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    type(ESMF_State),  intent(in)    :: State
    character(len=*),  intent(in)    :: fldname
    real(R8),          pointer       :: fldptr(:)
    logical,           intent(out)   :: exists
    integer,           intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_StateItem_Flag)   :: itemFlag
    type(ESMF_Field)            :: lfield
    type(ESMF_Mesh)             :: lmesh
    integer                     :: dbrc
    integer                     :: nnodes, nelements
    character(len=*), parameter :: subname='(lnd_import_export:state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 10) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if field exists then create output array - else do nothing
    if (itemflag == ESMF_STATEITEM_NOTFOUND) then
       exists = .false.
       RETURN
    else
       exists = .true.
    end if

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, status=status, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    else
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (nnodes == 0 .and. nelements == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", ESMF_LOGMSG_INFO, rc=dbrc)
          rc = ESMF_FAILURE
          return
       end if

       call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif  ! status

    if (dbug > 10) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine state_getfldptr

end module atm_import_export
