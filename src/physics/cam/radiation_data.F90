!================================================================================================
! output/input data necessary to drive radiation offline
! Francis Vitt -- Created 15 Dec 2009
!================================================================================================
module radiation_data

  use shr_kind_mod,     only: r8=>shr_kind_r8
  use ppgrid,           only : pcols, pver, pverp, begchunk, endchunk
  use cam_history,      only: addfld, add_default, horiz_only, outfld
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_gas, rad_cnst_get_aer_mmr
  use radconstants,     only: nradgas, gaslist
  use cam_history_support, only: fieldname_len, fillvalue
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun

  use drv_input_data,only: drv_input_data_t
  use drv_input_data,only: drv_input_data_get, drv_input_2d_t, drv_input_3d_t, drv_input_4d_t, drv_input_2di_t
  use physics_types, only: physics_state
  use physics_buffer,only: physics_buffer_desc, pbuf_get_chunk

  implicit none
  private

  public :: rad_data_readnl
  public :: rad_data_register
  public :: rad_data_init
  public :: rad_data_write
  public :: rad_data_read
  public :: rad_data_enable

  integer :: volcgeom_ifld, volcgeom1_ifld, volcgeom2_ifld, volcgeom3_ifld
  integer :: cld_ifld,rel_ifld,rei_ifld
  integer :: dei_ifld,mu_ifld,lambdac_ifld,iciwp_ifld,iclwp_ifld
  integer :: des_ifld,icswp_ifld,cldfsnow_ifld
  integer :: cldemis_ifld, cldtau_ifld, cicewp_ifld, cliqwp_ifld, nmxrgn_ifld, pmxrgn_ifld
  integer :: qrs_ifld, qrl_ifld
  integer :: dgnumwet_ifld, qaerwat_ifld

  character(len=fieldname_len), parameter :: &
       lndfrc_fldn    = 'rad_lndfrc      ' , &
       icefrc_fldn    = 'rad_icefrc      ' , &
       snowh_fldn     = 'rad_snowh       ' , &
       asdir_fldn     = 'rad_asdir       ' , &
       asdif_fldn     = 'rad_asdif       ' , &
       aldir_fldn     = 'rad_aldir       ' , &
       aldif_fldn     = 'rad_aldif       ' , &
       coszen_fldn    = 'rad_coszen      ' , &
       asdir_pos_fldn = 'rad_asdir_pos   ' , &
       asdif_pos_fldn = 'rad_asdif_pos   ' , &
       aldir_pos_fldn = 'rad_aldir_pos   ' , &
       aldif_pos_fldn = 'rad_aldif_pos   ' , &
       lwup_fldn      = 'rad_lwup        ' , &
       ts_fldn        = 'rad_ts          ' , &
       temp_fldn      = 'rad_temp        ' , &
       pdel_fldn      = 'rad_pdel        ' , &
       pdeldry_fldn   = 'rad_pdeldry     ' , &
       pmid_fldn      = 'rad_pmid        ' , &
       watice_fldn    = 'rad_watice      ' , &
       watliq_fldn    = 'rad_watliq      ' , &
       watvap_fldn    = 'rad_watvap      ' , &
       zint_fldn      = 'rad_zint        ' , &
       pint_fldn      = 'rad_pint        ' , &
       cld_fldn       = 'rad_cld         ' , &
       cldemis_fldn   = 'rad_cldemis     ' , &
       cldtau_fldn    = 'rad_cldtau      ' , &
       cicewp_fldn    = 'rad_cicewp      ' , &
       cliqwp_fldn    = 'rad_cliqwp      ' , &
       nmxrgn_fldn    = 'rad_nmxrgn      ' , &
       pmxrgn_fldn    = 'rad_pmxrgn      ' , &
       cldfsnow_fldn  = 'rad_cldfsnow    ' , &
       rel_fldn       = 'rad_rel         ' , &
       rei_fldn       = 'rad_rei         ' , &
       dei_fldn       = 'rad_dei         ' , &
       des_fldn       = 'rad_des         ' , &
       mu_fldn        = 'rad_mu          ' , &
       lambdac_fldn   = 'rad_lambdac     ' , &
       iciwp_fldn     = 'rad_iciwp       ' , &
       iclwp_fldn     = 'rad_iclwp       ' , &
       icswp_fldn     = 'rad_icswp       ' , &
       qrs_fldn       = 'rad_qrs         ' , &
       qrl_fldn       = 'rad_qrl         ' , &
       volcgeom_fldn  = 'rad_volc_geom   ' , &
       volcgeom1_fldn = 'rad_volc_geom1  ' , &
       volcgeom2_fldn = 'rad_volc_geom2  ' , &
       volcgeom3_fldn = 'rad_volc_geom3  '

  ! for modal aerosols
  character(len=fieldname_len), allocatable :: dgnumwet_fldn(:)
  character(len=fieldname_len), allocatable :: qaerwat_fldn(:)
  integer :: nmodes=0

  ! rad constituents mixing ratios
  integer :: ngas, naer=0
  character(len=64), allocatable :: gasnames(:)
  character(len=64), allocatable :: aernames(:)
 
  ! control options  
  logical          :: rad_data_output = .false.
  integer          :: rad_data_histfile_num = 2
  character(len=1) :: rad_data_avgflag = 'I'

  ! MG microphys check
  logical :: mg_microphys

  logical :: fixed_dyn_heating = .false.

  ! for fixed dynamical heating ...
  logical :: do_fdh = .false.

  integer :: tcorr_idx = -1
  integer :: qrsin_idx = -1
  integer :: qrlin_idx = -1
  integer :: tropp_idx = -1

  logical :: enabled = .false.
  logical :: gmean_3modes = .false.

contains


!================================================================================================
!================================================================================================
  subroutine rad_data_register
    use physics_buffer,  only: pbuf_add_field, dtype_r8

    if ( do_fdh ) then
      call pbuf_add_field('tcorr', 'global',  dtype_r8, (/pcols,pver/), tcorr_idx)
      call pbuf_add_field('qrsin', 'physpkg', dtype_r8, (/pcols,pver/), qrsin_idx) 
      call pbuf_add_field('qrlin', 'physpkg', dtype_r8, (/pcols,pver/), qrlin_idx) 
      call pbuf_add_field('tropp', 'physpkg', dtype_r8, (/pcols/),      tropp_idx) 
    end if 

  end subroutine rad_data_register

!================================================================================================
!================================================================================================
  subroutine rad_data_readnl(nlfile)

    ! Read rad_data_nl namelist group.  Parse input.

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    logical :: rad_data_fdh = .false.

    ! Local variables
    integer :: unitn, ierr, i
    character(len=*), parameter :: subname = 'rad_data_readnl'

    namelist /rad_data_nl/ rad_data_output, rad_data_histfile_num, rad_data_avgflag, rad_data_fdh

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'rad_data_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, rad_data_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast (rad_data_output,       1,   mpilog ,  0, mpicom)
    call mpibcast (rad_data_fdh,          1,   mpilog ,  0, mpicom)
    call mpibcast (rad_data_histfile_num, 1,   mpiint ,  0, mpicom)
    call mpibcast (rad_data_avgflag,      1,   mpichar , 0, mpicom)
#endif
    do_fdh = rad_data_fdh
    enabled = rad_data_output

  end subroutine rad_data_readnl

  !================================================================================================
  !================================================================================================
  subroutine rad_data_enable()
    enabled = .true.
  end subroutine rad_data_enable

  !================================================================================================
  !================================================================================================
  subroutine rad_data_init( pbuf2d )
    use phys_control,     only: phys_getopts
    use physics_buffer, only: pbuf_get_index
    use physics_buffer, only: pbuf_set_field
    implicit none
    
    ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! local vars
    integer :: i, err
    integer :: m,l, nspec
    character(len=64) :: name
    character(len=32) :: aername =  ' '
    character(len=128):: long_name
    character(len=64) :: long_name_description
    character(len=16)  :: microp_scheme  ! microphysics scheme
    character(len=16)  :: rad_scheme
   
    if (.not.enabled) return

    call phys_getopts(microp_scheme_out=microp_scheme, radiation_scheme_out=rad_scheme)
    mg_microphys =  (trim(microp_scheme) == 'MG')

    volcgeom_ifld = pbuf_get_index('VOLC_RAD_GEOM',errcode=err) ! might need 3 more for 3-mode inputs
    volcgeom1_ifld = pbuf_get_index('VOLC_RAD_GEOM1',errcode=err) ! might need 3 more for 3-mode inputs
    volcgeom2_ifld = pbuf_get_index('VOLC_RAD_GEOM2',errcode=err) ! might need 3 more for 3-mode inputs
    volcgeom3_ifld = pbuf_get_index('VOLC_RAD_GEOM3',errcode=err) ! might need 3 more for 3-mode inputs

    gmean_3modes = volcgeom1_ifld > 0 .and. volcgeom2_ifld > 0 .and. volcgeom3_ifld > 0 .and. volcgeom_ifld < 1

    if (volcgeom_ifld > 0) then
       call addfld(volcgeom_fldn, (/ 'lev' /), 'I','m', 'combined volcanic aerosol geometric-mode radius' )
    endif
    if (gmean_3modes) then
       call addfld(volcgeom1_fldn, (/ 'lev' /), 'I','m', 'mode 1 volcanic aerosol geometric-mode radius' )
       call addfld(volcgeom2_fldn, (/ 'lev' /), 'I','m', 'mode 2 volcanic aerosol geometric-mode radius' )
       call addfld(volcgeom3_fldn, (/ 'lev' /), 'I','m', 'mode 3 volcanic aerosol geometric-mode radius' )
    endif

    cld_ifld    = pbuf_get_index('CLD')
    rel_ifld    = pbuf_get_index('REL')
    rei_ifld    = pbuf_get_index('REI')
    qrs_ifld    = pbuf_get_index('QRS')
    qrl_ifld    = pbuf_get_index('QRL')
    if (mg_microphys) then
       dei_ifld      = pbuf_get_index('DEI')
       des_ifld      = pbuf_get_index('DES')
       mu_ifld       = pbuf_get_index('MU')
       lambdac_ifld  = pbuf_get_index('LAMBDAC')
       iciwp_ifld    = pbuf_get_index('ICIWP')
       iclwp_ifld    = pbuf_get_index('ICLWP')
       icswp_ifld    = pbuf_get_index('ICSWP')
       cldfsnow_ifld = pbuf_get_index('CLDFSNOW')
    else
       cldemis_ifld= pbuf_get_index('CLDEMIS')
       cldtau_ifld = pbuf_get_index('CLDTAU')
       cicewp_ifld = pbuf_get_index('CICEWP')
       cliqwp_ifld = pbuf_get_index('CLIQWP')
       nmxrgn_ifld = pbuf_get_index('NMXRGN')
       pmxrgn_ifld = pbuf_get_index('PMXRGN')
    endif
    
    if ( do_fdh ) then
       call addfld('rad_TROP_P', horiz_only,  'A','Pa','Pressure at which tropopause is defined for radiation' )
       call addfld('rad_TCORR',  (/ 'lev' /), 'I','K', 'Fixed dynamical heating temperature correction ' )
       call pbuf_set_field(pbuf2d, tcorr_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, qrsin_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, qrlin_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, tropp_idx, -1.0_r8)
    endif

    call rad_cnst_get_info(0, ngas=ngas, naero=naer, nmodes=nmodes)

    ! The code to output the gases assumes that the rad_constituents module has
    ! ordered them in the same way that they are ordered in the "gaslist" array
    ! in module radconstants, and that there are nradgas of them.  This ordering 
    ! is performed in the internal init_lists routine in rad_constituents.
    if (ngas /= nradgas) then
       call endrun('rad_data_init: ERROR: ngas /= nradgas')
    end if

    allocate( gasnames(ngas) )
    call rad_cnst_get_info(0, gasnames=gasnames)

    if (naer > 0) then
       allocate( aernames(naer) )
       call rad_cnst_get_info(0, aernames=aernames)
    endif

    if (nmodes>0) then
       allocate(dgnumwet_fldn(nmodes),qaerwat_fldn(nmodes))
       dgnumwet_ifld = pbuf_get_index('DGNUMWET')
       qaerwat_ifld  = pbuf_get_index('QAERWAT')
       do m = 1, nmodes
          write(dgnumwet_fldn(m), 1001 ) m
          write(qaerwat_fldn(m), 1003 ) m
       enddo
    endif

    if (.not.rad_data_output) return

    call addfld (lndfrc_fldn,    horiz_only,   rad_data_avgflag, 'fraction', &
         'radiation input: land fraction')
    call addfld (icefrc_fldn,    horiz_only,   rad_data_avgflag, 'fraction', &
         'radiation input: ice fraction')
    call addfld (snowh_fldn,     horiz_only,   rad_data_avgflag, 'm',        &
         'radiation input: water equivalent snow depth')
    call addfld (asdir_fldn,     horiz_only,   rad_data_avgflag, '1',        &
         'radiation input: s     hort wave direct albedo',                    flag_xyfill=.true.)
    call addfld (asdif_fldn,     horiz_only,   rad_data_avgflag, '1',        &
         'radiation input: s     hort wave difuse albedo',                    flag_xyfill=.true.)
    call addfld (aldir_fldn,     horiz_only,   rad_data_avgflag, '1',        &
         'radiation input: long wave direct albedo',                          flag_xyfill=.true.)
    call addfld (aldif_fldn,     horiz_only,   rad_data_avgflag, '1',        &
         'radiation input: long wave difuse albedo',                          flag_xyfill=.true.)

    call addfld (coszen_fldn,    horiz_only,   rad_data_avgflag, '1',        &
         'radiation input: cosine solar zenith when positive',                flag_xyfill=.true.)
    call addfld (asdir_pos_fldn, horiz_only,   rad_data_avgflag, '1',        &
         'radiation input: s     hort wave direct albedo weighted by coszen', flag_xyfill=.true.)
    call addfld (asdif_pos_fldn, horiz_only,   rad_data_avgflag, '1',        &
         'radiation input: s     hort wave difuse albedo weighted by coszen', flag_xyfill=.true.)
    call addfld (aldir_pos_fldn, horiz_only,   rad_data_avgflag, '1',        &
         'radiation input: long wave direct albedo weighted by coszen',       flag_xyfill=.true.)
    call addfld (aldif_pos_fldn, horiz_only,   rad_data_avgflag, '1',        &
         'radiation input: long wave difuse albedo weighted by coszen',       flag_xyfill=.true.)
    
    call addfld (lwup_fldn,      horiz_only,   rad_data_avgflag, 'W/m2',     &
         'radiation input: long wave up radiation flux ')
    call addfld (ts_fldn,        horiz_only,   rad_data_avgflag, 'K',        &
         'radiation input: surface temperature')

    call addfld (temp_fldn,      (/ 'lev' /),  rad_data_avgflag, 'K',        &
         'radiation input: midpoint temperature')
    call addfld (pdel_fldn,      (/ 'lev' /),  rad_data_avgflag, 'Pa',       &
         'radiation input: pressure layer thickness')
    call addfld (pdeldry_fldn,   (/ 'lev' /),  rad_data_avgflag,'Pa',        &
         'radiation input: dry pressure layer thickness')
    call addfld (pmid_fldn,      (/ 'lev' /),  rad_data_avgflag, 'Pa',       &
         'radiation input: midpoint pressure')
    call addfld (watice_fldn,    (/ 'lev' /),  rad_data_avgflag, 'kg/kg',    &
         'radiation input: cloud ice')
    call addfld (watliq_fldn,    (/ 'lev' /),  rad_data_avgflag, 'kg/kg',    &
         'radiation input: cloud liquid water')
    call addfld (watvap_fldn,    (/ 'lev' /),  rad_data_avgflag, 'kg/kg',    &
         'radiation input: water vapor')

    call addfld (zint_fldn,      (/ 'ilev' /), rad_data_avgflag, 'km',       &
         'radiation input: interface height')
    call addfld (pint_fldn,      (/ 'ilev' /), rad_data_avgflag, 'Pa',       &
         'radiation input: interface pressure')

    call addfld (cld_fldn,       (/ 'lev' /),  rad_data_avgflag, 'fraction', &
         'radiation input: cloud fraction')
    call addfld (rel_fldn,       (/ 'lev' /),  rad_data_avgflag, 'micron',   &
         'radiation input: effective liquid drop radius')
    call addfld (rei_fldn,       (/ 'lev' /),  rad_data_avgflag, 'micron',   &
         'radiation input: effective ice partical radius')
    call addfld (qrs_fldn,       (/ 'lev' /),  rad_data_avgflag, 'J/s/kg',   &
         'radiation input: solar heating rate')
    call addfld (qrl_fldn,       (/ 'lev' /),  rad_data_avgflag, 'J/s/kg',   &
         'radiation input: longwave heating rate')

    if (mg_microphys) then
       call addfld (dei_fldn,      (/ 'lev' /), rad_data_avgflag, 'micron',   &
            'radiation input: effective ice partical diameter')
       call addfld (des_fldn,      (/ 'lev' /), rad_data_avgflag, 'micron',   &
            'radiation input: effective snow partical diameter')
       call addfld (mu_fldn,       (/ 'lev' /), rad_data_avgflag, ' ',        &
            'radiation input: ice gamma parameter for optics (radiation)')
       call addfld (lambdac_fldn,  (/ 'lev' /), rad_data_avgflag,' ',         &
            'radiation input: slope of droplet distribution for optics (radiation)')
       call addfld (iciwp_fldn,    (/ 'lev' /), rad_data_avgflag, 'kg/m2',    &
            'radiation input: In-cloud ice water path')
       call addfld (iclwp_fldn,    (/ 'lev' /), rad_data_avgflag, 'kg/m2',    &
            'radiation input: In-cloud liquid water path')
       call addfld (icswp_fldn,    (/ 'lev' /), rad_data_avgflag, 'kg/m2',    &
            'radiation input: In-cloud snow water path')
       call addfld (cldfsnow_fldn, (/ 'lev' /), rad_data_avgflag, 'fraction', &
            'radiation input: cloud liquid drops + snow')
    else

      call addfld (nmxrgn_fldn,    horiz_only,   rad_data_avgflag, ' ',       &
         'radiation input: ')
      call addfld (pmxrgn_fldn,    (/ 'ilev' /), rad_data_avgflag, 'Pa',      &
         'radiation input: ')
      call addfld (cldemis_fldn,   (/ 'lev' /), rad_data_avgflag, ' ',        &
         'radiation input: cloud property ')
      call addfld (cldtau_fldn,    (/ 'lev' /), rad_data_avgflag, ' ',        &
         'radiation input: cloud property ')
      call addfld (cicewp_fldn,    (/ 'lev' /), rad_data_avgflag, ' ',        &
         'radiation input: cloud property ')
      call addfld (cliqwp_fldn,    (/ 'lev' /), rad_data_avgflag, ' ',        &
         'radiation input: cloud property ')
    endif

    call add_default (lndfrc_fldn,    rad_data_histfile_num, ' ')
    call add_default (icefrc_fldn,    rad_data_histfile_num, ' ')
    call add_default (snowh_fldn,     rad_data_histfile_num, ' ')
    call add_default (asdir_fldn,     rad_data_histfile_num, ' ')
    call add_default (asdif_fldn,     rad_data_histfile_num, ' ')
    call add_default (aldir_fldn,     rad_data_histfile_num, ' ')
    call add_default (aldif_fldn,     rad_data_histfile_num, ' ')

    call add_default (coszen_fldn,    rad_data_histfile_num, ' ')
    call add_default (asdir_pos_fldn, rad_data_histfile_num, ' ')
    call add_default (asdif_pos_fldn, rad_data_histfile_num, ' ')
    call add_default (aldir_pos_fldn, rad_data_histfile_num, ' ')
    call add_default (aldif_pos_fldn, rad_data_histfile_num, ' ')

    call add_default (lwup_fldn,      rad_data_histfile_num, ' ')
    call add_default (ts_fldn,        rad_data_histfile_num, ' ')
    call add_default (temp_fldn,      rad_data_histfile_num, ' ')
    call add_default (pdel_fldn,      rad_data_histfile_num, ' ')
    call add_default (pdeldry_fldn,   rad_data_histfile_num, ' ')
    call add_default (pmid_fldn,      rad_data_histfile_num, ' ')
    call add_default (watice_fldn,    rad_data_histfile_num, ' ')
    call add_default (watliq_fldn,    rad_data_histfile_num, ' ')
    call add_default (watvap_fldn,    rad_data_histfile_num, ' ')
    call add_default (zint_fldn,      rad_data_histfile_num, ' ')
    call add_default (pint_fldn,      rad_data_histfile_num, ' ')

    call add_default (cld_fldn,       rad_data_histfile_num, ' ')
    call add_default (rel_fldn,       rad_data_histfile_num, ' ')
    call add_default (rei_fldn,       rad_data_histfile_num, ' ')
    call add_default (qrs_fldn,       rad_data_histfile_num, ' ')
    call add_default (qrl_fldn,       rad_data_histfile_num, ' ')
    
    if (mg_microphys) then
       call add_default (dei_fldn,       rad_data_histfile_num, ' ')
       call add_default (des_fldn,       rad_data_histfile_num, ' ')
       call add_default (mu_fldn,        rad_data_histfile_num, ' ')
       call add_default (lambdac_fldn,   rad_data_histfile_num, ' ')
       call add_default (iciwp_fldn,     rad_data_histfile_num, ' ')
       call add_default (iclwp_fldn,     rad_data_histfile_num, ' ')
       call add_default (icswp_fldn,     rad_data_histfile_num, ' ')
       call add_default (cldfsnow_fldn,  rad_data_histfile_num, ' ')
    else
       call add_default (cldemis_fldn,   rad_data_histfile_num, ' ')
       call add_default (cldtau_fldn,    rad_data_histfile_num, ' ')
       call add_default (cicewp_fldn,    rad_data_histfile_num, ' ')
       call add_default (cliqwp_fldn,    rad_data_histfile_num, ' ')
       call add_default (nmxrgn_fldn,    rad_data_histfile_num, ' ')
       call add_default (pmxrgn_fldn,    rad_data_histfile_num, ' ')
    endif

    ! stratospheric aersols geometric mean radius
    if (volcgeom_ifld > 0) then
       call add_default (volcgeom_fldn, rad_data_histfile_num, ' ')
    endif
    if (gmean_3modes) then
       call add_default (volcgeom1_fldn, rad_data_histfile_num, ' ')
       call add_default (volcgeom2_fldn, rad_data_histfile_num, ' ')
       call add_default (volcgeom3_fldn, rad_data_histfile_num, ' ')
    endif

    ! rad constituents

    long_name_description = ' mass mixing ratio used in rad climate calculation'

    do i = 1, ngas
       long_name = trim(gasnames(i))//trim(long_name_description)
       name = 'rad_'//gasnames(i)
       call addfld(trim(name), (/ 'lev' /), rad_data_avgflag, 'kg/kg', trim(long_name))
       call add_default (trim(name), rad_data_histfile_num, ' ')
    end do
    
    if (naer > 0) then

       do i = 1, naer
          long_name = trim(aernames(i))//trim(long_name_description)
          name = 'rad_'//aernames(i)
          call addfld(trim(name), (/ 'lev' /), rad_data_avgflag, 'kg/kg', trim(long_name))
          call add_default (trim(name), rad_data_histfile_num, ' ')
       end do
    endif
    if (nmodes>0) then

       ! for modal aerosol model
       do m = 1, nmodes
          write(long_name, 1002) m
          call addfld ( dgnumwet_fldn(m), (/ 'lev' /), rad_data_avgflag, '', trim(long_name) )
          call add_default (dgnumwet_fldn(m), rad_data_histfile_num, ' ')

          write(long_name, 1004) m
          call addfld ( qaerwat_fldn(m),  (/ 'lev' /), rad_data_avgflag, '', trim(long_name) )
          call add_default (qaerwat_fldn(m), rad_data_histfile_num, ' ')

          ! get mode info
          call rad_cnst_get_info(0, m, nspec=nspec)
          ! aerosol species loop
          do l = 1, nspec
             call rad_cnst_get_info(0,m,l, spec_name=aername)
             name = 'rad_'//trim(aername)
             call addfld(trim(name),      (/ 'lev' /), rad_data_avgflag, 'kg/kg', trim(long_name))
             call add_default (trim(name), rad_data_histfile_num, ' ')
          end do
       end do
    endif

 1001 format ( 'rad_dgnumwet_', I1 )
 1002 format ( 'radiation input: DGNUMWET for mode ', I1 )
 1003 format ( 'rad_qaerwat_', I1 )
 1004 format ( 'radiation input: QAERWAT for mode ', I1 )

  end subroutine rad_data_init

  !================================================================================================
  !================================================================================================
  subroutine rad_data_write(  pbuf, state, cam_in, coszen )

    use physics_types,    only: physics_state
    use camsrfexch,       only: cam_in_t     
    use constituents,     only: cnst_get_ind
    use physics_buffer,   only: pbuf_get_field, pbuf_old_tim_idx
    use drv_input_data,   only: drv_input_data_freq
    use physconst,        only: cpair

    implicit none
    type(physics_buffer_desc), pointer :: pbuf(:)
    
    type(physics_state), intent(in), target :: state
    type(cam_in_t),      intent(in) :: cam_in
    real(r8),            intent(in) :: coszen(pcols)

    ! Local variables
    integer :: i, k
    integer :: m,l, nspec
    character(len=32) :: name, aername
    real(r8), pointer :: mmr(:,:)

    integer :: lchnk, itim_old, ifld
    integer :: ixcldice              ! cloud ice water index
    integer :: ixcldliq              ! cloud liquid water index
    integer :: icol
    integer :: ncol
    
    ! surface albedoes weighted by (positive cosine zenith angle)
    real(r8):: coszrs_pos(pcols)    ! = max(coszrs,0)
    real(r8):: asdir_pos (pcols)    !
    real(r8):: asdif_pos (pcols)    !
    real(r8):: aldir_pos (pcols)    !
    real(r8):: aldif_pos (pcols)    !

    real(r8), pointer, dimension(:,:)  :: ptr
    integer , pointer, dimension(:)    :: iptr1d
    real(r8),          dimension(pcols) :: rptr1d

    ! for Fixed Dynamical Heating
    real(r8), pointer, dimension(:,:) :: tcorr
    real(r8), pointer, dimension(:,:) :: qrsin
    real(r8), pointer, dimension(:,:) :: qrlin
    real(r8), pointer, dimension(:)   :: tropp
    real(r8), pointer, dimension(:,:) :: qrs
    real(r8), pointer, dimension(:,:) :: qrl
    real(r8), pointer, dimension(:,:,:) :: dgnumwet_ptr
    real(r8), pointer, dimension(:,:,:) :: qaerwat_ptr

    if (.not.enabled) return
    call pbuf_get_field(pbuf, qrs_ifld, qrs )    
    call pbuf_get_field(pbuf, qrl_ifld, qrl )

    lchnk = state%lchnk
    ncol = state%ncol

   ! store temperature correction above tropopause in pbuf for Fixed Dynamical Heating
    if(do_fdh) then
       
       call pbuf_get_field(pbuf, tcorr_idx, tcorr)
       call pbuf_get_field(pbuf, qrsin_idx, qrsin)
       call pbuf_get_field(pbuf, qrlin_idx, qrlin)
       call pbuf_get_field(pbuf, tropp_idx, tropp)
       do k = 1, pver
          do i = 1, ncol
             if ( state%pmid(i,k) < tropp(i)) then
                tcorr(i,k) = tcorr(i,k) + drv_input_data_freq*(qrs(i,k) - qrsin(i,k) + qrl(i,k) - qrlin(i,k)) /cpair 
             endif
          enddo
       enddo

       call outfld('rad_TROP_P', tropp(:ncol), ncol, lchnk)
       call outfld('rad_TCORR', tcorr(:ncol,:), ncol, lchnk)

    end if

    if (.not.rad_data_output) return

    call outfld(qrs_fldn, qrs, pcols, lchnk )
    call outfld(qrl_fldn, qrl, pcols, lchnk )

    ! get index of (liquid+ice) cloud water
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

    do icol = 1, ncol
       coszrs_pos(icol)  = max(coszen(icol),0._r8)
    enddo
    asdir_pos(:ncol)  = cam_in%asdir(:ncol) * coszrs_pos(:ncol)
    asdif_pos(:ncol)  = cam_in%asdif(:ncol) * coszrs_pos(:ncol)
    aldir_pos(:ncol)  = cam_in%aldir(:ncol) * coszrs_pos(:ncol)
    aldif_pos(:ncol)  = cam_in%aldif(:ncol) * coszrs_pos(:ncol)

    call outfld(lndfrc_fldn, cam_in%landfrac,  pcols, lchnk)
    call outfld(icefrc_fldn, cam_in%icefrac,   pcols, lchnk)
    call outfld(snowh_fldn,  cam_in%snowhland, pcols, lchnk)
    call outfld(temp_fldn,   state%t,               pcols, lchnk   )
    call outfld(pdel_fldn,   state%pdel,            pcols, lchnk   )
    call outfld(pdeldry_fldn,state%pdeldry,         pcols, lchnk   )
    call outfld(watice_fldn, state%q(:,:,ixcldice), pcols, lchnk   )
    call outfld(watliq_fldn, state%q(:,:,ixcldliq), pcols, lchnk   )
    call outfld(watvap_fldn, state%q(:,:,1),        pcols, lchnk   )
    call outfld(zint_fldn,   state%zi,              pcols, lchnk   )
    call outfld(pint_fldn,   state%pint,            pcols, lchnk   )
    call outfld(pmid_fldn,   state%pmid,            pcols, lchnk   )

    call outfld(asdir_fldn, cam_in%asdir, pcols, lchnk   )
    call outfld(asdif_fldn, cam_in%asdif, pcols, lchnk   )
    call outfld(aldir_fldn, cam_in%aldir, pcols, lchnk   )
    call outfld(aldif_fldn, cam_in%aldif, pcols, lchnk   )

    call outfld(coszen_fldn, coszrs_pos, pcols, lchnk   )
    call outfld(asdir_pos_fldn, asdir_pos, pcols, lchnk   )
    call outfld(asdif_pos_fldn, asdif_pos, pcols, lchnk   )
    call outfld(aldir_pos_fldn, aldir_pos, pcols, lchnk   )
    call outfld(aldif_pos_fldn, aldif_pos, pcols, lchnk   )

    call outfld(lwup_fldn,  cam_in%lwup,  pcols, lchnk   )
    call outfld(ts_fldn,    cam_in%ts,    pcols, lchnk   )

    itim_old = pbuf_old_tim_idx()

    call pbuf_get_field(pbuf, cld_ifld,    ptr, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
    call outfld(cld_fldn,    ptr, pcols, lchnk )

    call pbuf_get_field(pbuf, rel_ifld,    ptr )
    call outfld(rel_fldn,    ptr, pcols, lchnk )

    call pbuf_get_field(pbuf, rei_ifld,    ptr )
    call outfld(rei_fldn,    ptr, pcols, lchnk )
    
    if (mg_microphys) then

       call pbuf_get_field(pbuf,  dei_ifld, ptr     )
       call outfld(dei_fldn,      ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  des_ifld, ptr     )
       call outfld(des_fldn,      ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  mu_ifld, ptr      )
       call outfld(mu_fldn,       ptr, pcols, lchnk   ) 

       call pbuf_get_field(pbuf,  lambdac_ifld, ptr )
       call outfld(lambdac_fldn,  ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  iciwp_ifld, ptr   )
       call outfld(iciwp_fldn,    ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  iclwp_ifld, ptr   )
       call outfld(iclwp_fldn,    ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  icswp_ifld, ptr   )
       call outfld(icswp_fldn,    ptr, pcols, lchnk   )       

       call pbuf_get_field(pbuf,  cldfsnow_ifld, ptr, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
       call outfld(cldfsnow_fldn, ptr, pcols, lchnk   )

    else

       call pbuf_get_field(pbuf, cldemis_ifld,    ptr )
       call outfld(cldemis_fldn,   ptr, pcols, lchnk )

       call pbuf_get_field(pbuf, cldtau_ifld,    ptr )
       call outfld(cldtau_fldn,    ptr, pcols, lchnk )

       call pbuf_get_field(pbuf, cicewp_ifld,    ptr )
       call outfld(cicewp_fldn,    ptr, pcols, lchnk )

       call pbuf_get_field(pbuf, cliqwp_ifld,    ptr )
       call outfld(cliqwp_fldn,    ptr, pcols, lchnk )

       call pbuf_get_field(pbuf, nmxrgn_ifld,    iptr1d )
       rptr1d = dble(iptr1d)
       call outfld(nmxrgn_fldn,    rptr1d, pcols, lchnk )

       call pbuf_get_field(pbuf, pmxrgn_ifld,    ptr )
       call outfld(pmxrgn_fldn,    ptr, pcols, lchnk )

    endif

    ! output mixing ratio of rad constituents 

    do i = 1, ngas
       name = 'rad_'//gasnames(i)
       call rad_cnst_get_gas(0, gaslist(i), state, pbuf, mmr)
       call outfld(trim(name), mmr, pcols, lchnk)
    end do

    if ( naer>0 ) then
       ! bulk aerosol model
       do i = 1, naer
          name = 'rad_'//aernames(i)
          call rad_cnst_get_aer_mmr(0, i, state, pbuf, mmr)
          call outfld(trim(name), mmr, pcols, lchnk)
       end do
    endif
    if (nmodes>0) then

       call pbuf_get_field(pbuf, dgnumwet_ifld, dgnumwet_ptr)
       call pbuf_get_field(pbuf, qaerwat_ifld, qaerwat_ptr)

       ! for modal aerosol model
       do m = 1, nmodes
          ptr => dgnumwet_ptr(:,:,m)
          call outfld( dgnumwet_fldn(m), ptr, pcols, lchnk )
          ptr => qaerwat_ptr(:,:,m)
          call outfld( qaerwat_fldn(m),  ptr, pcols, lchnk )

          ! get mode info
          call rad_cnst_get_info(0, m, nspec=nspec)
          ! aerosol species loop
          do l = 1, nspec
             call rad_cnst_get_info(0,m,l, spec_name=aername)
             call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, mmr)
             name = 'rad_'//aername
             call outfld(trim(name), mmr, pcols, lchnk)
          enddo
       enddo
    endif

    ! stratospheric aersols geometric mean radius
    if (volcgeom_ifld > 0) then
       call pbuf_get_field(pbuf, volcgeom_ifld, ptr)
       call outfld(volcgeom_fldn, ptr, pcols, lchnk)
    endif
    if (gmean_3modes) then
       call pbuf_get_field(pbuf, volcgeom1_ifld, ptr)
       call outfld(volcgeom1_fldn, ptr, pcols, lchnk)
       call pbuf_get_field(pbuf, volcgeom2_ifld, ptr)
       call outfld(volcgeom2_fldn, ptr, pcols, lchnk)
       call pbuf_get_field(pbuf, volcgeom3_ifld, ptr)
       call outfld(volcgeom3_fldn, ptr, pcols, lchnk)
    endif

  end subroutine rad_data_write

!=================================================================================
!=================================================================================
  subroutine rad_data_read(indata, phys_state, pbuf2d, cam_in, recno )

    use camsrfexch,       only: cam_in_t
    use physics_buffer,   only: pbuf_get_field, pbuf_old_tim_idx
    use constituents,     only: cnst_get_ind
    use tropopause,       only: tropopause_find, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE

    implicit none

    type(drv_input_data_t), intent(inout) :: indata
    type(physics_buffer_desc), pointer         :: pbuf2d(:,:)
    type(physics_state), target, intent(inout) :: phys_state(begchunk:endchunk)
    type(cam_in_t),      target, intent(inout) :: cam_in(begchunk:endchunk)
    integer,             intent(in)            :: recno

! local vars
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(drv_input_3d_t) :: cld_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: cldfsnow_ptrs(begchunk:endchunk)

    type(drv_input_3d_t) :: rel_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: rei_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: dei_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: des_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: mu_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: lambdac_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: iciwp_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: iclwp_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: icswp_ptrs(begchunk:endchunk)

    type(drv_input_3d_t) :: qrs_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: qrl_ptrs(begchunk:endchunk)

    type(drv_input_3d_t) :: watvap_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: watliq_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: watice_ptrs(begchunk:endchunk)

    type(drv_input_3d_t) :: zi_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: pint_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: lnpint_ptrs(begchunk:endchunk)

    type(drv_input_3d_t) :: temp_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: pdel_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: pdeldry_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: lnpmid_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: pmid_ptrs(begchunk:endchunk)

    type(drv_input_2d_t) :: lndfrac_ptrs(begchunk:endchunk)
    type(drv_input_2d_t) :: icefrac_ptrs(begchunk:endchunk)
    type(drv_input_2d_t) :: snowhland_ptrs(begchunk:endchunk)
    type(drv_input_2d_t) :: asdir_ptrs(begchunk:endchunk)
    type(drv_input_2d_t) :: asdif_ptrs(begchunk:endchunk)
    type(drv_input_2d_t) :: aldir_ptrs(begchunk:endchunk)
    type(drv_input_2d_t) :: aldif_ptrs(begchunk:endchunk)

    type(drv_input_2di_t):: nmxrgn_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: pmxrgn_ptrs(begchunk:endchunk)

    type(drv_input_3d_t) :: cldemis_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: cldtau_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: cicewp_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: cliqwp_ptrs(begchunk:endchunk)

    type(drv_input_2d_t) :: lwup_ptrs(begchunk:endchunk)
    type(drv_input_2d_t) :: ts_ptrs(begchunk:endchunk)

    type(drv_input_3d_t) :: volcgeom_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: volcgeom1_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: volcgeom2_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: volcgeom3_ptrs(begchunk:endchunk)

    integer :: i, k, c, ncol, itim

    integer :: ixcldice   ! cloud ice water index
    integer :: ixcldliq   ! cloud liquid water index

    ! for Fixed Dynamical Heating
    real(r8), pointer, dimension(:,:) :: tcorr
    real(r8), pointer, dimension(:,:) :: qrsin
    real(r8), pointer, dimension(:,:) :: qrlin
    real(r8), pointer, dimension(:)   :: tropp
    integer :: troplev(pcols)

    ! for modal aerosols
    type(drv_input_4d_t) :: dgnumwet_ptrs(begchunk:endchunk)
    type(drv_input_4d_t) :: qaerwat_ptrs(begchunk:endchunk)

    ! get index of (liquid+ice) cloud water
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

    ! phys buffer time index
    itim = pbuf_old_tim_idx()

    ! setup the data pointers
!$OMP PARALLEL DO PRIVATE (C,pbuf)
    do c=begchunk,endchunk
       pbuf => pbuf_get_chunk(pbuf2d, c)

       watvap_ptrs (c)%array => phys_state(c)%q(:,:,1)
       watliq_ptrs (c)%array => phys_state(c)%q(:,:,ixcldliq)
       watice_ptrs (c)%array => phys_state(c)%q(:,:,ixcldice)

       zi_ptrs     (c)%array => phys_state(c)%zi
       pint_ptrs   (c)%array => phys_state(c)%pint
       lnpint_ptrs (c)%array => phys_state(c)%lnpint

       temp_ptrs   (c)%array => phys_state(c)%t
       pdel_ptrs   (c)%array => phys_state(c)%pdel
       pdeldry_ptrs(c)%array => phys_state(c)%pdeldry
       lnpmid_ptrs (c)%array => phys_state(c)%lnpmid
       pmid_ptrs   (c)%array => phys_state(c)%pmid

       lndfrac_ptrs  (c)%array => cam_in(c)%landfrac
       icefrac_ptrs  (c)%array => cam_in(c)%icefrac
       snowhland_ptrs(c)%array => cam_in(c)%snowhland
       asdir_ptrs    (c)%array => cam_in(c)%asdir
       asdif_ptrs    (c)%array => cam_in(c)%asdif
       aldir_ptrs    (c)%array => cam_in(c)%aldir
       aldif_ptrs    (c)%array => cam_in(c)%aldif
       lwup_ptrs     (c)%array => cam_in(c)%lwup
       ts_ptrs       (c)%array => cam_in(c)%ts

       call pbuf_get_field(pbuf, cld_ifld, cld_ptrs   (c)%array, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
       call pbuf_get_field(pbuf, rel_ifld,   rel_ptrs(c)%array )
       call pbuf_get_field(pbuf, rei_ifld,   rei_ptrs(c)%array )

       call pbuf_get_field(pbuf, qrs_ifld,   qrs_ptrs(c)%array )
       call pbuf_get_field(pbuf, qrl_ifld,   qrl_ptrs(c)%array )

       if (mg_microphys) then
          call pbuf_get_field(pbuf, dei_ifld,   dei_ptrs(c)%array )
          call pbuf_get_field(pbuf, des_ifld,   des_ptrs(c)%array )
          call pbuf_get_field(pbuf, mu_ifld,    mu_ptrs     (c)%array )
          call pbuf_get_field(pbuf, lambdac_ifld, lambdac_ptrs(c)%array )
          call pbuf_get_field(pbuf, iciwp_ifld, iciwp_ptrs  (c)%array )
          call pbuf_get_field(pbuf, iclwp_ifld, iclwp_ptrs  (c)%array )
          call pbuf_get_field(pbuf, icswp_ifld, icswp_ptrs  (c)%array )
          call pbuf_get_field(pbuf, cldfsnow_ifld, cldfsnow_ptrs(c)%array, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
       else
          call pbuf_get_field(pbuf, nmxrgn_ifld,   nmxrgn_ptrs(c)%array )
          call pbuf_get_field(pbuf, pmxrgn_ifld,   pmxrgn_ptrs(c)%array )
          call pbuf_get_field(pbuf, cldemis_ifld,  cldemis_ptrs(c)%array )
          call pbuf_get_field(pbuf, cldtau_ifld,   cldtau_ptrs(c)%array )
          call pbuf_get_field(pbuf, cicewp_ifld,   cicewp_ptrs(c)%array )
          call pbuf_get_field(pbuf, cliqwp_ifld,   cliqwp_ptrs(c)%array )
       endif

       if (nmodes>0) then
          call pbuf_get_field(pbuf, dgnumwet_ifld, dgnumwet_ptrs(c)%array )
          call pbuf_get_field(pbuf, qaerwat_ifld, qaerwat_ptrs(c)%array )
       endif

       ! stratospheric aersols geometric mean radius
       if (volcgeom_ifld > 0) then
          call pbuf_get_field(pbuf, volcgeom_ifld, volcgeom_ptrs(c)%array )
       endif
       if (gmean_3modes) then
          call pbuf_get_field(pbuf, volcgeom1_ifld, volcgeom1_ptrs(c)%array )
          call pbuf_get_field(pbuf, volcgeom2_ifld, volcgeom2_ptrs(c)%array )
          call pbuf_get_field(pbuf, volcgeom3_ifld, volcgeom3_ptrs(c)%array )
       endif

    enddo

    
    if (nmodes>0) then
       call get_aeromode_data( indata, dgnumwet_fldn, recno, dgnumwet_ptrs )
       call get_aeromode_data( indata, qaerwat_fldn, recno, qaerwat_ptrs )
    endif

    ! get the 2D data

    call drv_input_data_get( indata, lndfrc_fldn, recno, lndfrac_ptrs )
    call drv_input_data_get( indata, icefrc_fldn, recno, icefrac_ptrs )
    call drv_input_data_get( indata, snowh_fldn,  recno, snowhland_ptrs )
    call drv_input_data_get( indata, asdir_fldn,  recno, asdir_ptrs )
    call drv_input_data_get( indata, asdif_fldn,  recno, asdif_ptrs )
    call drv_input_data_get( indata, aldir_fldn,  recno, aldir_ptrs )
    call drv_input_data_get( indata, aldif_fldn,  recno, aldif_ptrs )
    call drv_input_data_get( indata, lwup_fldn,   recno, lwup_ptrs )
    call drv_input_data_get( indata, ts_fldn,     recno, ts_ptrs )

    ! get the 3D data

    call drv_input_data_get( indata, cld_fldn,    'lev', pver, recno, cld_ptrs )
    call drv_input_data_get( indata, rel_fldn,    'lev', pver, recno, rel_ptrs )
    call drv_input_data_get( indata, rei_fldn,    'lev', pver, recno, rei_ptrs )
    call drv_input_data_get( indata, qrs_fldn,    'lev', pver, recno, qrs_ptrs )
    call drv_input_data_get( indata, qrl_fldn,    'lev', pver, recno, qrl_ptrs )

    if (mg_microphys) then
       call drv_input_data_get( indata, dei_fldn,    'lev', pver, recno, dei_ptrs )
       call drv_input_data_get( indata, des_fldn,    'lev', pver, recno, des_ptrs )
       call drv_input_data_get( indata, mu_fldn,     'lev', pver, recno, mu_ptrs )
       call drv_input_data_get( indata, lambdac_fldn,'lev', pver, recno, lambdac_ptrs )
       call drv_input_data_get( indata, iciwp_fldn,  'lev', pver, recno, iciwp_ptrs )
       call drv_input_data_get( indata, iclwp_fldn,  'lev', pver, recno, iclwp_ptrs )
       call drv_input_data_get( indata, icswp_fldn,  'lev', pver, recno, icswp_ptrs )
       call drv_input_data_get( indata, cldfsnow_fldn,'lev',pver, recno, cldfsnow_ptrs )
    else
       call drv_input_data_get( indata, nmxrgn_fldn, recno, nmxrgn_ptrs )
       call drv_input_data_get( indata, pmxrgn_fldn, 'ilev',pverp,recno, pmxrgn_ptrs )
       call drv_input_data_get( indata, cldemis_fldn,'lev', pver, recno, cldemis_ptrs )
       call drv_input_data_get( indata, cldtau_fldn, 'lev', pver, recno, cldtau_ptrs )
       call drv_input_data_get( indata, cicewp_fldn, 'lev', pver, recno, cicewp_ptrs )
       call drv_input_data_get( indata, cliqwp_fldn, 'lev', pver, recno, cliqwp_ptrs )
    endif

    ! stratospheric aersols geometric mean radius
    if (volcgeom_ifld > 0) then
       call drv_input_data_get( indata, volcgeom_fldn, 'lev', pver, recno, volcgeom_ptrs )
    endif
    if (gmean_3modes) then
       call drv_input_data_get( indata, volcgeom1_fldn, 'lev', pver, recno, volcgeom1_ptrs )
       call drv_input_data_get( indata, volcgeom2_fldn, 'lev', pver, recno, volcgeom2_ptrs )
       call drv_input_data_get( indata, volcgeom3_fldn, 'lev', pver, recno, volcgeom3_ptrs )
    endif

    call drv_input_data_get( indata, watvap_fldn, 'lev', pver, recno, watvap_ptrs )
    call drv_input_data_get( indata, watliq_fldn, 'lev', pver, recno, watliq_ptrs )
    call drv_input_data_get( indata, watice_fldn, 'lev', pver, recno, watice_ptrs )
    call drv_input_data_get( indata, temp_fldn,   'lev', pver, recno, temp_ptrs )
    call drv_input_data_get( indata, pdel_fldn,   'lev', pver, recno, pdel_ptrs )
    call drv_input_data_get( indata, pdeldry_fldn,'lev', pver, recno, pdeldry_ptrs )
    call drv_input_data_get( indata, pmid_fldn,   'lev', pver, recno, pmid_ptrs )

    call drv_input_data_get( indata, zint_fldn,   'ilev',pverp,recno, zi_ptrs )
    call drv_input_data_get( indata, pint_fldn,   'ilev',pverp,recno, pint_ptrs )

    do c=begchunk,endchunk
       ncol = phys_state(c)%ncol
       if (all(pmid_ptrs(c)%array(:ncol,:) > 0._r8 )) then
          do i=1,ncol
             lnpmid_ptrs(c)%array(i,:) = log(pmid_ptrs(c)%array(i,:))
             lnpint_ptrs(c)%array(i,:) = log(pint_ptrs(c)%array(i,:))
          enddo
       endif

       ! adjust temperatue input above tropopause for Fixed Dynamical Heating ....
       if (do_fdh) then
          pbuf => pbuf_get_chunk(pbuf2d, c)

          call pbuf_get_field(pbuf, tropp_idx, tropp)
          call pbuf_get_field(pbuf, tcorr_idx, tcorr)
          call pbuf_get_field(pbuf, qrsin_idx, qrsin)
          call pbuf_get_field(pbuf, qrlin_idx, qrlin)

          call tropopause_find(phys_state(c), troplev, tropP=tropp(:), primary=TROP_ALG_CLIMATE, &
                               backup=TROP_ALG_CLIMATE)

          qrsin(:,:) = qrs_ptrs(c)%array(:,:)
          qrlin(:,:) = qrl_ptrs(c)%array(:,:)

          do k = 1, pver
             do i = 1, ncol
                if ( phys_state(c)%pmid(i,k) < tropp(i) ) then
                   temp_ptrs(c)%array(i,k) = temp_ptrs(c)%array(i,k) + tcorr(i,k)
                endif
             enddo
          enddo
       endif

    enddo

    call get_rad_cnst_data(indata, phys_state, pbuf2d, recno)

  end subroutine rad_data_read

!================================================================================================
! Private routines
!================================================================================================


  !================================================================================================
  !================================================================================================
  subroutine get_aeromode_data(indata, infld_names, recno, chunk_ptrs)
    use drv_input_data, only: drv_input_data_read

    type(drv_input_data_t), intent(inout) :: indata
    character(len=*),    intent(in)    :: infld_names(nmodes)
    integer,             intent(in)    :: recno
    type(drv_input_4d_t),   intent(inout) :: chunk_ptrs(begchunk:endchunk)

    real(r8) :: data (pcols, pver, begchunk:endchunk)

    integer :: c, n

    do n = 1,nmodes

       data = drv_input_data_read( indata, infld_names(n), 'lev', pver, recno )
       do c=begchunk,endchunk
          chunk_ptrs(c)%array(:,:,n) = data(:,:,c)
       enddo

    enddo

  end subroutine get_aeromode_data

  !================================================================================================
  !================================================================================================
  subroutine get_rad_cnst_data(indata, state, pbuf2d, recno)

    use physics_types,    only: physics_state
    use physics_buffer,   only: physics_buffer_desc
    use rad_constituents, only: rad_cnst_get_info

    implicit none

    ! Arguments
    type(drv_input_data_t), intent(inout) :: indata
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_state), intent(inout) :: state(begchunk:endchunk)
    integer,             intent(in)    :: recno

    ! Local variables
    integer :: i, ncol
    integer :: m, l, nspec
    integer :: idx
    character(len=32) :: aername

    !-----------------------------------------------------------------------------

    ! read in mixing ratio of rad constituents 

    do i = 1,ngas
       call read_rad_gas_data(indata, gasnames(i), i, state, pbuf2d, recno )
    enddo

    do i = 1,naer
       call read_rad_aer_data(indata, aernames(i), i, state, pbuf2d, recno )
    enddo
    
    ! aerosol mode loop
    do m = 1, nmodes

       ! get mode info
       call rad_cnst_get_info(0, m, nspec=nspec)
       ! aerosol species loop
       do l = 1, nspec
          call rad_cnst_get_info(0,m,l, spec_name=aername)
          call read_rad_mam_data( indata, aername, m, l, state, pbuf2d, recno )
       end do
    end do

  end subroutine get_rad_cnst_data

  !=================================================================================
  !=================================================================================
  subroutine read_rad_gas_data(indata, name, idx, state, pbuf2d, recno )
    use rad_constituents, only: rad_cnst_get_gas
    use drv_input_data,   only: drv_input_data_read

    type(drv_input_data_t), intent(inout) :: indata
    character(len=32),   intent(in)    :: name
    integer,             intent(in)    :: idx
    type(physics_state),target, intent(inout) :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
    integer,             intent(in)    :: recno

    type(physics_buffer_desc), pointer, dimension(:) :: phys_buffer_chunk
    character(len=32) :: radname
    integer            :: c, ncol
    real(r8)          :: mass(pcols,pver,begchunk:endchunk)
    real(r8), pointer :: mmr_ptr(:,:)

    radname = 'rad_'//trim(name)
    mass = drv_input_data_read( indata, radname, 'lev', pver, recno )

    do c = begchunk,endchunk
       ncol = state(c)%ncol
       phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
       call rad_cnst_get_gas(0, gaslist(idx), state(c), phys_buffer_chunk, mmr_ptr)
       mmr_ptr(:ncol,:) = mass(:ncol,:,c)
    enddo

  end subroutine read_rad_gas_data

  !=================================================================================
  !=================================================================================
  subroutine read_rad_aer_data(indata, name, idx, state, pbuf2d, recno )
    use rad_constituents, only: rad_cnst_get_aer_mmr
    use drv_input_data,   only: drv_input_data_read

    type(drv_input_data_t), intent(inout) :: indata
    character(len=32),   intent(in)    :: name
    integer,             intent(in)    :: idx
    type(physics_state),target, intent(inout) :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
    integer,             intent(in)    :: recno

    type(physics_buffer_desc), pointer, dimension(:) :: phys_buffer_chunk
    character(len=32) :: radname
    integer :: c, ncol
    real(r8)          :: mass(pcols,pver,begchunk:endchunk)
    real(r8), pointer :: mmr_ptr(:,:)

    radname = 'rad_'//trim(name)
    mass = drv_input_data_read( indata, radname, 'lev', pver, recno )

    do c = begchunk,endchunk
       ncol = state(c)%ncol
       phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
       call rad_cnst_get_aer_mmr(0, idx, state(c), phys_buffer_chunk, mmr_ptr)
       mmr_ptr(:ncol,:) = mass(:ncol,:,c)
    enddo

  end subroutine read_rad_aer_data

  !=================================================================================
  !=================================================================================
  subroutine read_rad_mam_data(indata, name, mode_idx, spec_idx, state, pbuf2d, recno )
    use rad_constituents, only: rad_cnst_get_aer_mmr
    use drv_input_data,   only: drv_input_data_read

    type(drv_input_data_t), intent(inout) :: indata
    character(len=32),   intent(in)    :: name
    integer,             intent(in)    :: mode_idx
    integer,             intent(in)    :: spec_idx
    type(physics_state),target, intent(inout) :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
    integer,             intent(in)    :: recno

    type(physics_buffer_desc), pointer, dimension(:) :: phys_buffer_chunk
    character(len=32) :: radname
    integer :: c, ncol
    real(r8)          :: mass(pcols,pver,begchunk:endchunk)
    real(r8), pointer :: mmr_ptr(:,:)

    radname = 'rad_'//trim(name)
    mass = drv_input_data_read( indata, radname, 'lev', pver, recno )

    do c = begchunk,endchunk
       ncol = state(c)%ncol
       phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
       call rad_cnst_get_aer_mmr(0, mode_idx, spec_idx, 'a', state(c), phys_buffer_chunk, mmr_ptr)
       mmr_ptr(:ncol,:) = mass(:ncol,:,c)
    enddo

  end subroutine read_rad_mam_data

end module radiation_data
