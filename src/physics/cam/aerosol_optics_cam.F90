module aerosol_optics_cam
  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_kind_mod, only: cl => shr_kind_cl
  use cam_logfile,  only: iulog
  use radconstants, only: nswbands, nlwbands, idx_sw_diag, idx_uv_diag, idx_nir_diag
  use radconstants, only: get_lw_spectral_boundaries
  use phys_prop,    only: ot_length
  use physics_types,only: physics_state
  use physics_buffer,only: physics_buffer_desc
  use ppgrid, only: pcols, pver
  use physconst, only: rga, rair
  use cam_abortutils, only: endrun
  use spmd_utils, only: masterproc
  use rad_constituents,  only: n_diag, rad_cnst_get_call_list
  use cam_history,       only: addfld, add_default, outfld, horiz_only, fieldname_len
  use cam_history_support, only: fillvalue

  use tropopause, only : tropopause_findChemTrop

  use aerosol_properties_mod, only: aerosol_properties
  use modal_aerosol_properties_mod, only: modal_aerosol_properties

  use aerosol_state_mod,      only: aerosol_state
  use modal_aerosol_state_mod,only: modal_aerosol_state

  use aerosol_optics_mod,     only: aerosol_optics
  use refractive_aerosol_optics_mod, only: refractive_aerosol_optics

  implicit none

  private

  public :: aerosol_optics_cam_readnl
  public :: aerosol_optics_cam_init
  public :: aerosol_optics_cam_final
  public :: aerosol_optics_cam_sw
  public :: aerosol_optics_cam_lw

  type aero_props_t
     class(aerosol_properties), pointer :: obj => null()
  end type aero_props_t
  type aero_state_t
     class(aerosol_state), pointer :: obj => null()
  end type aero_state_t

  type(aero_props_t), allocatable :: aero_props(:) ! array of aerosol properties objects to allow for
                                                   ! multiple aerosol representations in the same sim
                                                   ! such as MAM and CARMA

  ! refractive index for water read in read_water_refindex
  complex(r8) :: crefwsw(nswbands) = -huge(1._r8) ! complex refractive index for water visible
  complex(r8) :: crefwlw(nlwbands) = -huge(1._r8) ! complex refractive index for water infrared
  character(len=cl) :: water_refindex_file = 'NONE' ! full pathname for water refractive index dataset

  logical :: modal_active = .false.
  integer :: num_aero_models = 0
  integer :: lw10um_indx = -1            ! wavelength index corresponding to 10 microns
  real(r8), parameter :: lw10um = 10._r8 ! microns

  character(len=4) :: diag(0:n_diag) = (/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

  type out_name
     character(len=fieldname_len), allocatable :: name(:) ! nbins
  end type out_name

  type(out_name), allocatable :: burden_fields(:) ! num_aero_models
  type(out_name), allocatable :: aodbin_fields(:)
  type(out_name), allocatable :: aoddust_fields(:)
  type(out_name), allocatable :: burdendn_fields(:) ! num_aero_models
  type(out_name), allocatable :: aodbindn_fields(:)
  type(out_name), allocatable :: aoddustdn_fields(:)

contains

  !===============================================================================
  subroutine aerosol_optics_cam_readnl(nlfile)
    use namelist_utils, only : find_group_name
    use spmd_utils,     only : mpicom, masterprocid, mpi_character, mpi_success

    character(len=*), intent(in)  :: nlfile  ! filepath for file containing namelist input

    integer                       :: unitn, ierr
    character(len=cl)             :: errmsg
    character(len=*), parameter   :: subname = 'aerosol_optics_cam_readnl'

    ! ===================
    ! Namelist definition
    ! ===================
    namelist /aerosol_optics_nl/ water_refindex_file

    ! =============
    ! Read namelist
    ! =============
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aerosol_optics_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aerosol_optics_nl, iostat=ierr)
          if (ierr /= 0) then
             write(errmsg,'(2a,i10)') subname,':: ERROR reading namelist, error code: ',ierr
             call endrun(errmsg)
          end if
       end if
       close(unitn)
    end if

    ! ============================
    ! Broadcast namelist variables
    ! ============================
    call mpi_bcast(water_refindex_file, len(water_refindex_file),  mpi_character, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname // ':: ERROR mpi_bcast '//trim(water_refindex_file))
    end if

    if (masterproc) then
       write(iulog,*) subname,': water_refindex_file = ',trim(water_refindex_file)
    end if

  end subroutine aerosol_optics_cam_readnl

  !===============================================================================
  subroutine aerosol_optics_cam_init
    use rad_constituents, only: rad_cnst_get_info
    use phys_control,     only: phys_getopts
    use ioFileMod,        only: getfil

    character(len=*), parameter :: prefix = 'aerosol_optics_cam_init: '
    integer :: nmodes=0, iaermod, istat, ilist, i

    logical :: call_list(0:n_diag)
    real(r8) :: lwavlen_lo(nlwbands), lwavlen_hi(nlwbands)
    integer :: m, n

    character(len=fieldname_len) :: fldname
    character(len=128) :: lngname
    logical :: history_aero_optics     ! output aerosol optics diagnostics
    logical :: history_amwg            ! output the variables used by the AMWG diag package
    logical :: history_dust            ! output dust diagnostics

    character(len=cl) :: locfile

    call phys_getopts(history_amwg_out        = history_amwg, &
                      history_aero_optics_out = history_aero_optics, &
                      history_dust_out        = history_dust )

    num_aero_models = 0

    call rad_cnst_get_info(0, nmodes=nmodes)
    modal_active = nmodes>0

    if (modal_active) then
       num_aero_models = num_aero_models+1 ! count aerosol models
    end if

    if (num_aero_models>0) then
       allocate(aero_props(num_aero_models), stat=istat)
       if (istat/=0) then
          call endrun(prefix//'array allocation error: aero_props')
       end if
    end if

    iaermod = 0

    if (modal_active) then
       iaermod = iaermod+1
       aero_props(iaermod)%obj => modal_aerosol_properties()
    end if

    if (water_refindex_file=='NONE') then
       call endrun(prefix//'water_refindex_file must be specified')
    else
       call getfil(water_refindex_file, locfile)
       call read_water_refindex(locfile)
    end if

    call get_lw_spectral_boundaries(lwavlen_lo, lwavlen_hi, units='um')
    do i = 1,nlwbands
       if ((lwavlen_lo(i)<=lw10um) .and. (lwavlen_hi(i)>=lw10um)) then
          lw10um_indx = i ! index corresponding to 10 microns
       end if
    end do
    call rad_cnst_get_call_list(call_list)

    do ilist = 0, n_diag
       if (call_list(ilist)) then
          call addfld ('EXTINCT'//diag(ilist),    (/ 'lev' /), 'A','/m',&
               'Aerosol extinction 550 nm, day only', flag_xyfill=.true.)
          call addfld ('EXTINCTUV'//diag(ilist),  (/ 'lev' /), 'A','/m',&
               'Aerosol extinction 350 nm, day only', flag_xyfill=.true.)
          call addfld ('EXTINCTNIR'//diag(ilist), (/ 'lev' /), 'A','/m',&
               'Aerosol extinction 1020 nm, day only', flag_xyfill=.true.)
          call addfld ('ABSORB'//diag(ilist),     (/ 'lev' /), 'A','/m',&
               'Aerosol absorption, day only', flag_xyfill=.true.)
          call addfld ('AODVIS'//diag(ilist),   horiz_only,  'A','  ', &
               'Aerosol optical depth 550 nm', flag_xyfill=.true.)
          call addfld ('AODVISst'//diag(ilist), horiz_only,  'A','  ', &
               'Stratospheric aerosol optical depth 550 nm, day only', flag_xyfill=.true.)
          call addfld ('AODNIRst'//diag(ilist), horiz_only,  'A','  ', &
               'Stratospheric aerosol optical depth 1020 nm, day only',flag_xyfill=.true.)
          call addfld ('AODUVst'//diag(ilist),  horiz_only,  'A','  ', &
               'Stratospheric aerosol optical depth 350 nm, day only', flag_xyfill=.true.)
          call addfld ('AODUV'//diag(ilist),      horiz_only,  'A','  ', &
               'Aerosol optical depth 350 nm, day only', flag_xyfill=.true.)
          call addfld ('AODNIR'//diag(ilist),     horiz_only,  'A','  ', &
               'Aerosol optical depth 1020 nm, day only',flag_xyfill=.true.)
          call addfld ('AODABS'//diag(ilist),     horiz_only,  'A','  ', &
               'Aerosol absorption optical depth 550 nm, day only', flag_xyfill=.true.)
          call addfld ('AODxASYM'//diag(ilist),   horiz_only,  'A','  ', &
               'Aerosol optical depth 550 * asymmetry factor, day only', flag_xyfill=.true.)
          call addfld ('EXTxASYM'//diag(ilist),   (/ 'lev' /), 'A','  ', &
               'extinction 550 nm * asymmetry factor, day only',  flag_xyfill=.true.)
          call addfld ('AODTOT'//diag(ilist), horiz_only, 'A','1',&
               'Aerosol optical depth summed over all sw wavelengths', flag_xyfill=.true.)

          call addfld ('EXTINCTdn'//diag(ilist),    (/ 'lev' /), 'A','/m',&
               'Aerosol extinction 550 nm, day only', flag_xyfill=.true.)
          call addfld ('EXTINCTUVdn'//diag(ilist),  (/ 'lev' /), 'A','/m',&
               'Aerosol extinction 350 nm, day only', flag_xyfill=.true.)
          call addfld ('EXTINCTNIRdn'//diag(ilist), (/ 'lev' /), 'A','/m',&
               'Aerosol extinction 1020 nm, day only', flag_xyfill=.true.)
          call addfld ('ABSORBdn'//diag(ilist),     (/ 'lev' /), 'A','/m',&
               'Aerosol absorption, day only', flag_xyfill=.true.)
          call addfld ('AODVISdn'//diag(ilist),   horiz_only,  'A','  ', &
               'Aerosol optical depth 550 nm', flag_xyfill=.true.)
          call addfld ('AODVISstdn'//diag(ilist), horiz_only,  'A','  ', &
               'Stratospheric aerosol optical depth 550 nm, day only', flag_xyfill=.true.)
          call addfld ('AODNIRstdn'//diag(ilist), horiz_only,  'A','  ', &
               'Stratospheric aerosol optical depth 1020 nm, day only', flag_xyfill=.true.)
          call addfld ('AODUVstdn'//diag(ilist),  horiz_only,  'A','  ', &
               'Stratospheric aerosol optical depth 350 nm, day only', flag_xyfill=.true.)
          call addfld ('AODUVdn'//diag(ilist),      horiz_only,  'A','  ', &
               'Aerosol optical depth 350 nm, day only', flag_xyfill=.true.)
          call addfld ('AODNIRdn'//diag(ilist),     horiz_only,  'A','  ', &
               'Aerosol optical depth 1020 nm, day only', flag_xyfill=.true.)
          call addfld ('AODABSdn'//diag(ilist),     horiz_only,  'A','  ', &
               'Aerosol absorption optical depth 550 nm, day only', flag_xyfill=.true.)
          call addfld ('AODxASYMdn'//diag(ilist),   horiz_only,  'A','  ', &
               'Aerosol optical depth 550 * asymmetry factor, day only', flag_xyfill=.true.)
          call addfld ('EXTxASYMdn'//diag(ilist),   (/ 'lev' /), 'A','  ', &
               'extinction 550 nm * asymmetry factor, day only',  flag_xyfill=.true.)
          call addfld ('AODTOTdn'//diag(ilist), horiz_only, 'A','1',&
               'Aerosol optical depth summed over all sw wavelengths, day only')

          if (lw10um_indx>0) then
             call addfld('AODABSLW'//diag(ilist), (/ 'lev' /), 'A','/m',&
                  'Aerosol long-wave absorption optical depth at 10 microns')
          end if
          call addfld ('TOTABSLW'//diag(ilist), (/ 'lev' /), 'A',' ', &
               'LW Aero total abs')

          if (ilist>0 .and. history_aero_optics) then
             call add_default ('EXTINCT'//diag(ilist), 1, ' ')
             call add_default ('ABSORB'//diag(ilist),  1, ' ')
             call add_default ('AODVIS'//diag(ilist),  1, ' ')
             call add_default ('AODVISst'//diag(ilist),  1, ' ')
             call add_default ('AODABS'//diag(ilist),  1, ' ')
          end if

       end if
    end do

    if (num_aero_models>0) then

       allocate(burden_fields(num_aero_models), stat=istat)
       if (istat/=0) then
          call endrun(prefix//'array allocation error: burden_fields')
       end if
       allocate(aodbin_fields(num_aero_models), stat=istat)
       if (istat/=0) then
          call endrun(prefix//'array allocation error: aodbin_fields')
       end if
       allocate(aoddust_fields(num_aero_models), stat=istat)
       if (istat/=0) then
          call endrun(prefix//'array allocation error: aoddust_fields')
       end if

       allocate(burdendn_fields(num_aero_models), stat=istat)
       if (istat/=0) then
          call endrun(prefix//'array allocation error: burdendn_fields')
       end if
       allocate(aodbindn_fields(num_aero_models), stat=istat)
       if (istat/=0) then
          call endrun(prefix//'array allocation error: aodbindn_fields')
       end if
       allocate(aoddustdn_fields(num_aero_models), stat=istat)
       if (istat/=0) then
          call endrun(prefix//'array allocation error: aoddustdn_fields')
       end if

       do n = 1,num_aero_models

          allocate(burden_fields(n)%name(aero_props(n)%obj%nbins()), stat=istat)
          if (istat/=0) then
             call endrun(prefix//'array allocation error: burden_fields(n)%name')
          end if
          allocate(aodbin_fields(n)%name(aero_props(n)%obj%nbins()), stat=istat)
          if (istat/=0) then
             call endrun(prefix//'array allocation error: aodbin_fields(n)%name')
          end if
          allocate(aoddust_fields(n)%name(aero_props(n)%obj%nbins()), stat=istat)
          if (istat/=0) then
             call endrun(prefix//'array allocation error: aoddust_fields(n)%name')
          end if

          allocate(burdendn_fields(n)%name(aero_props(n)%obj%nbins()), stat=istat)
          if (istat/=0) then
             call endrun(prefix//'array allocation error: burdendn_fields(n)%name')
          end if
          allocate(aodbindn_fields(n)%name(aero_props(n)%obj%nbins()), stat=istat)
          if (istat/=0) then
             call endrun(prefix//'array allocation error: aodbindn_fields(n)%name')
          end if
          allocate(aoddustdn_fields(n)%name(aero_props(n)%obj%nbins()), stat=istat)
          if (istat/=0) then
             call endrun(prefix//'array allocation error: aoddustdn_fields(n)%name')
          end if

          do m = 1, aero_props(n)%obj%nbins()

             write(fldname,'(a,i2.2)') 'BURDEN', m
             burden_fields(n)%name(m) = fldname
             write(lngname,'(a,i2.2)') 'Aerosol burden bin ', m
             call addfld (fldname, horiz_only, 'A', 'kg/m2', lngname, flag_xyfill=.true.)
             if (history_aero_optics) then
                call add_default (fldname, 1, ' ')
             end if

             fldname = 'AOD_'//trim(aero_props(n)%obj%bin_name(0,m))
             aodbin_fields(n)%name(m) = fldname
             lngname = 'Aerosol optical depth, day only, 550 nm, '//trim(aero_props(n)%obj%bin_name(0,m))
             call addfld (aodbin_fields(n)%name(m), horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
             if (history_aero_optics) then
                call add_default (fldname, 1, ' ')
             end if

             write(fldname,'(a,i2.2)') 'AODDUST', m
             aoddust_fields(n)%name(m) = fldname
             write(lngname,'(a,i2,a)') 'Aerosol optical depth, day only, 550 nm mode ',m,' from dust'
             call addfld (aoddust_fields(n)%name(m), horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
             if (history_aero_optics) then
                call add_default (fldname, 1, ' ')
             end if

             write(fldname,'(a,i2.2)') 'BURDENdn', m
             burdendn_fields(n)%name(m) = fldname
             write(lngname,'(a,i2)') 'Aerosol burden, day night, bin ', m
             call addfld (burdendn_fields(n)%name(m), horiz_only, 'A', 'kg/m2', lngname, flag_xyfill=.true.)
             if (history_aero_optics) then
                call add_default (fldname, 1, ' ')
             end if

             fldname = 'AODdn_'//trim(aero_props(n)%obj%bin_name(0,m))
             aodbindn_fields(n)%name(m) = fldname
             lngname = 'Aerosol optical depth 550 nm, day night, '//trim(aero_props(n)%obj%bin_name(0,m))
             call addfld (aodbindn_fields(n)%name(m), horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
             if (history_aero_optics) then
                call add_default (fldname, 1, ' ')
             end if

             write(fldname,'(a,i2.2)') 'AODdnDUST', m
             aoddustdn_fields(n)%name(m) = fldname
             write(lngname,'(a,i2,a)') 'Aerosol optical depth 550 nm, day night, bin ',m,' from dust'
             call addfld (aoddustdn_fields(n)%name(m), horiz_only, 'A', '  ', lngname, flag_xyfill=.true.)
             if (history_aero_optics) then
                call add_default (fldname, 1, ' ')
             end if

          end do

       end do

    end if

   call addfld ('AODDUST',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust, day only',         &
        flag_xyfill=.true.)
   call addfld ('AODSO4',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SO4, day only',          &
        flag_xyfill=.true.)
   call addfld ('AODPOM',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from POM, day only',          &
        flag_xyfill=.true.)
   call addfld ('AODSOA',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SOA, day only',          &
        flag_xyfill=.true.)
   call addfld ('AODBC',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from BC, day only',           &
        flag_xyfill=.true.)
   call addfld ('AODSS',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from seasalt, day only',      &
        flag_xyfill=.true.)
   call addfld ('AODABSBC',      horiz_only, 'A','  ',    'Aerosol absorption optical depth 550 nm from BC, day only',&
        flag_xyfill=.true.)
   call addfld ('BURDENDUST',    horiz_only, 'A','kg/m2', 'Dust aerosol burden, day only'        ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSO4',     horiz_only, 'A','kg/m2', 'Sulfate aerosol burden, day only'     ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENPOM',     horiz_only, 'A','kg/m2', 'POM aerosol burden, day only'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSOA',     horiz_only, 'A','kg/m2', 'SOA aerosol burden, day only'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENBC',      horiz_only, 'A','kg/m2', 'Black carbon aerosol burden, day only',                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSEASALT', horiz_only, 'A','kg/m2', 'Seasalt aerosol burden, day only'     ,                    &
        flag_xyfill=.true.)
   call addfld ('SSAVIS',        horiz_only, 'A','  ',    'Aerosol single-scatter albedo, day only',                  &
        flag_xyfill=.true.)

   call addfld ('AODDUSTdn',       horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from dust, day night',         &
        flag_xyfill=.true.)
   call addfld ('AODSO4dn',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SO4, day night',          &
        flag_xyfill=.true.)
   call addfld ('AODPOMdn',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from POM, day night',          &
        flag_xyfill=.true.)
   call addfld ('AODSOAdn',        horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from SOA, day night',          &
        flag_xyfill=.true.)
   call addfld ('AODBCdn',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from BC, day night',           &
        flag_xyfill=.true.)
   call addfld ('AODSSdn',         horiz_only, 'A','  ',    'Aerosol optical depth 550 nm from seasalt, day night',      &
        flag_xyfill=.true.)
   call addfld ('AODABSBCdn',      horiz_only, 'A','  ',    'Aerosol absorption optical depth 550 nm from BC, day night',&
        flag_xyfill=.true.)
   call addfld ('BURDENDUSTdn',    horiz_only, 'A','kg/m2', 'Dust aerosol burden, day night'        ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSO4dn',     horiz_only, 'A','kg/m2', 'Sulfate aerosol burden, day night'     ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENPOMdn',     horiz_only, 'A','kg/m2', 'POM aerosol burden, day night'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSOAdn',     horiz_only, 'A','kg/m2', 'SOA aerosol burden, day night'         ,                    &
        flag_xyfill=.true.)
   call addfld ('BURDENBCdn',      horiz_only, 'A','kg/m2', 'Black carbon aerosol burden, day night',                    &
        flag_xyfill=.true.)
   call addfld ('BURDENSEASALTdn', horiz_only, 'A','kg/m2', 'Seasalt aerosol burden, day night'     ,                    &
        flag_xyfill=.true.)
   call addfld ('SSAVISdn',        horiz_only, 'A','  ',    'Aerosol single-scatter albedo, day night',                  &
        flag_xyfill=.true.)

   if (history_amwg) then
      call add_default ('AODDUST01'     , 1, ' ')
      call add_default ('AODDUST03'     , 1, ' ')
      call add_default ('AODDUST'      , 1, ' ')
      call add_default ('AODVIS'       , 1, ' ')
   end if

   if (history_dust) then
      call add_default ('AODDUST01'     , 1, ' ')
      call add_default ('AODDUST02'     , 1, ' ')
      call add_default ('AODDUST03'     , 1, ' ')
   end if

   if (history_aero_optics) then
      call add_default ('AODDUST01'     , 1, ' ')
      call add_default ('AODDUST03'     , 1, ' ')
      call add_default ('ABSORB'       , 1, ' ')
      call add_default ('AODVIS'       , 1, ' ')
      call add_default ('AODUV'        , 1, ' ')
      call add_default ('AODNIR'       , 1, ' ')
      call add_default ('AODABS'       , 1, ' ')
      call add_default ('AODABSBC'     , 1, ' ')
      call add_default ('AODDUST'      , 1, ' ')
      call add_default ('AODSO4'       , 1, ' ')
      call add_default ('AODPOM'       , 1, ' ')
      call add_default ('AODSOA'       , 1, ' ')
      call add_default ('AODBC'        , 1, ' ')
      call add_default ('AODSS'        , 1, ' ')
      call add_default ('BURDEN01'      , 1, ' ')
      call add_default ('BURDEN02'      , 1, ' ')
      call add_default ('BURDEN03'      , 1, ' ')
      call add_default ('BURDENDUST'   , 1, ' ')
      call add_default ('BURDENSO4'    , 1, ' ')
      call add_default ('BURDENPOM'    , 1, ' ')
      call add_default ('BURDENSOA'    , 1, ' ')
      call add_default ('BURDENBC'     , 1, ' ')
      call add_default ('BURDENSEASALT', 1, ' ')
      call add_default ('SSAVIS'       , 1, ' ')
      call add_default ('EXTINCT'      , 1, ' ')
      call add_default ('AODxASYM'     , 1, ' ')
      call add_default ('EXTxASYM'     , 1, ' ')

      call add_default ('AODdnDUST01'     , 1, ' ')
      call add_default ('AODdnDUST03'     , 1, ' ')
      call add_default ('ABSORBdn'       , 1, ' ')
      call add_default ('AODVISdn'       , 1, ' ')
      call add_default ('AODUVdn'        , 1, ' ')
      call add_default ('AODNIRdn'       , 1, ' ')
      call add_default ('AODABSdn'       , 1, ' ')
      call add_default ('AODABSBCdn'     , 1, ' ')
      call add_default ('AODDUSTdn'      , 1, ' ')
      call add_default ('AODSO4dn'       , 1, ' ')
      call add_default ('AODPOMdn'       , 1, ' ')
      call add_default ('AODSOAdn'       , 1, ' ')
      call add_default ('AODBCdn'        , 1, ' ')
      call add_default ('AODSSdn'        , 1, ' ')
      call add_default ('BURDENdn01'      , 1, ' ')
      call add_default ('BURDENdn02'      , 1, ' ')
      call add_default ('BURDENdn03'      , 1, ' ')
      call add_default ('BURDENDUSTdn'   , 1, ' ')
      call add_default ('BURDENSO4dn'    , 1, ' ')
      call add_default ('BURDENPOMdn'    , 1, ' ')
      call add_default ('BURDENSOAdn'    , 1, ' ')
      call add_default ('BURDENBCdn'     , 1, ' ')
      call add_default ('BURDENSEASALTdn', 1, ' ')
      call add_default ('SSAVISdn'       , 1, ' ')
      call add_default ('EXTINCTdn'      , 1, ' ')
      call add_default ('AODxASYMdn'     , 1, ' ')
      call add_default ('EXTxASYMdn'     , 1, ' ')
   end if

  end subroutine aerosol_optics_cam_init

  !===============================================================================
  subroutine aerosol_optics_cam_final

    integer :: iaermod

    do iaermod = 1,num_aero_models
       if (associated(aero_props(iaermod)%obj)) then
          deallocate(aero_props(iaermod)%obj)
          nullify(aero_props(iaermod)%obj)
       end if
    end do

    if (allocated(aero_props)) then
       deallocate(aero_props)
    endif

  end subroutine aerosol_optics_cam_final

  !===============================================================================
  subroutine aerosol_optics_cam_sw(list_idx, state, pbuf, nnite, idxnite, tauxar, wa, ga, fa)

    ! calculates aerosol sw radiative properties

    integer,             intent(in) :: list_idx       ! index of the climate or a diagnostic list
    type(physics_state), intent(in), target :: state          ! state variables

    type(physics_buffer_desc), pointer :: pbuf(:)
    integer,             intent(in) :: nnite          ! number of night columns
    integer,             intent(in) :: idxnite(nnite) ! local column indices of night columns

    real(r8), intent(inout) :: tauxar(pcols,0:pver,nswbands) ! layer extinction optical depth
    real(r8), intent(inout) :: wa(pcols,0:pver,nswbands)     ! layer single-scatter albedo
    real(r8), intent(inout) :: ga(pcols,0:pver,nswbands)     ! asymmetry factor
    real(r8), intent(inout) :: fa(pcols,0:pver,nswbands)     ! forward scattered fraction

    character(len=*), parameter :: prefix = 'aerosol_optics_cam_sw: '

    integer :: ibin, nbins
    integer :: iwav, ilev
    integer :: icol, istat
    integer :: lchnk, ncol

    type(aero_state_t), allocatable :: aero_state(:) ! array of aerosol state objects to allow for
                                                     ! multiple aerosol representations in the same sim
                                                     ! such as MAM and CARMA

    class(aerosol_optics), pointer :: aero_optics

    real(r8) :: dopaer(pcols)
    real(r8) :: mass(pcols,pver)
    real(r8) :: air_density(pcols,pver)

    real(r8), allocatable :: pext(:)  ! parameterized specific extinction (m2/kg)
    real(r8), allocatable :: pabs(:)  ! parameterized specific absorption (m2/kg)
    real(r8), allocatable :: palb(:)  ! parameterized single scattering albedo
    real(r8), allocatable :: pasm(:)  ! parameterized asymmetry factor

    character(len=ot_length) :: opticstype
    integer :: iaermod

    real(r8) :: aodvis(pcols)              ! extinction optical depth in vis
    real(r8) :: aoduv(pcols)               ! extinction optical depth in uv
    real(r8) :: aodnir(pcols)              ! extinction optical depth in nir
    real(r8) :: absorb(pcols,pver)
    real(r8) :: aodabs(pcols)              ! absorption optical depth

    real(r8) :: aodabsbc(pcols)             ! absorption optical depth of BC

    real(r8) :: aodtot(pcols)

    real(r8) :: extinct(pcols,pver)
    real(r8) :: extinctnir(pcols,pver)
    real(r8) :: extinctuv(pcols,pver)

    real(r8) :: asymvis(pcols)              ! asymmetry factor * optical depth
    real(r8) :: asymext(pcols,pver)         ! asymmetry factor * extinction

    real(r8) :: wetvol(pcols,pver)
    real(r8) :: watervol(pcols,pver)

    real(r8) :: vol(pcols)
    real(r8) :: dustvol(pcols)

    real(r8) :: scatdust(pcols)
    real(r8) :: absdust(pcols)
    real(r8) :: dustaodbin(pcols)

    real(r8) :: scatbc(pcols)
    real(r8) :: absbc(pcols)

    real(r8) :: scatpom(pcols)
    real(r8) :: abspom(pcols)

    real(r8) :: scatsslt(pcols)
    real(r8) :: abssslt(pcols)

    real(r8) :: scatsoa(pcols)
    real(r8) :: abssoa(pcols)

    real(r8) :: scatsulf(pcols)
    real(r8) :: abssulf(pcols)

    real(r8) :: burden(pcols)
    real(r8) :: burdendust(pcols), burdenso4(pcols), burdenbc(pcols), &
                burdenpom(pcols), burdensoa(pcols), burdenseasalt(pcols)

    real(r8) :: hygrodust(pcols), hygrosulf(pcols), hygrobc(pcols), &
                hygropom(pcols), hygrosoa(pcols), hygrosslt(pcols)

    real(r8) :: aodbin(pcols)

    complex(r8), pointer :: specrefindex(:)     ! species refractive index

    class(aerosol_state),      pointer :: aerostate
    class(aerosol_properties), pointer :: aeroprops

    real(r8) :: specdens
    character(len=32) :: spectype            ! species type
    real(r8), pointer :: specmmr(:,:)
    real(r8)          :: hygro_aer           !

    real(r8) :: scath2o, absh2o, sumscat, sumabs, sumhygro

    real(r8) :: aodc                        ! aod of component

    ! total species AOD
    real(r8) :: dustaod(pcols), sulfaod(pcols), bcaod(pcols), &
                pomaod(pcols), soaaod(pcols), ssltaod(pcols)

    real(r8) :: aodvisst(pcols) ! stratospheric extinction optical depth
    real(r8) :: aoduvst(pcols)  ! stratospheric extinction optical depth in uv
    real(r8) :: aodnirst(pcols) ! stratospheric extinction optical depth in nir
    real(r8) :: ssavis(pcols)
    integer :: troplev(pcols)

    nullify(aero_optics)

    call tropopause_findChemTrop(state, troplev)

    lchnk = state%lchnk
    ncol  = state%ncol

    mass(:ncol,:)        = state%pdeldry(:ncol,:)*rga
    air_density(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))

    aodvis = 0._r8
    aodnir = 0._r8
    aoduv = 0._r8
    aodabs = 0._r8
    absorb = 0._r8
    aodtot = 0._r8
    tauxar = 0._r8
    extinct = 0._r8
    extinctnir = 0._r8
    extinctuv = 0._r8
    asymvis = 0.0_r8
    asymext = 0.0_r8
    ssavis = 0.0_r8
    aodvisst = 0.0_r8
    aoduvst = 0.0_r8
    aodnirst = 0.0_r8

    burdendust = 0.0_r8
    burdenso4 = 0.0_r8
    burdenbc = 0.0_r8
    burdenpom = 0.0_r8
    burdensoa = 0.0_r8
    burdenseasalt = 0.0_r8

    aodabsbc = 0.0_r8
    dustaod = 0.0_r8
    sulfaod = 0.0_r8
    pomaod = 0.0_r8
    soaaod = 0.0_r8
    bcaod = 0.0_r8
    ssltaod = 0.0_r8

    if (num_aero_models<1) return

    allocate(aero_state(num_aero_models), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: aero_state')
    end if

    iaermod = 0
    if (modal_active) then
       iaermod = iaermod+1
       aero_state(iaermod)%obj => modal_aerosol_state( state, pbuf )
    end if

    allocate(pext(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: pext')
    end if
    allocate(pabs(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: pabs')
    end if
    allocate(palb(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: palb')
    end if
    allocate(pasm(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: pasm')
    end if

    aeromodel: do iaermod = 1,num_aero_models

       aeroprops => aero_props(iaermod)%obj
       aerostate => aero_state(iaermod)%obj

       nbins=aeroprops%nbins(list_idx)

       binloop: do ibin = 1, nbins

          dustaodbin(:) = 0._r8
          burden(:) = 0._r8
          aodbin(:) = 0.0_r8

          call aeroprops%optics_params(list_idx, ibin, opticstype=opticstype)

          select case (trim(opticstype))
          case('modal') ! refractive method
             aero_optics=>refractive_aerosol_optics(aeroprops, aerostate, list_idx, ibin, &
                                                    ncol, pver, nswbands, nlwbands, crefwsw, crefwlw)
          case default
             call endrun(prefix//'optics method not recognized')
          end select

          if (associated(aero_optics)) then

             wetvol(:ncol,:pver) = aerostate%wet_volume(aeroprops, list_idx, ibin, ncol, pver)
             watervol(:ncol,:pver) = aerostate%water_volume(aeroprops, list_idx, ibin, ncol, pver)

             wavelength: do iwav = 1, nswbands

                vertical: do ilev = 1, pver

                   call aero_optics%sw_props(ncol, ilev, iwav, pext, pabs, palb, pasm )

                   call init_diags

                   column: do icol = 1,ncol
                      dopaer(icol) = pext(icol)*mass(icol,ilev)
                      tauxar(icol,ilev,iwav) = tauxar(icol,ilev,iwav) + dopaer(icol)
                      wa(icol,ilev,iwav) = wa(icol,ilev,iwav) + dopaer(icol)*palb(icol)
                      ga(icol,ilev,iwav) = ga(icol,ilev,iwav) + dopaer(icol)*palb(icol)*pasm(icol)
                      fa(icol,ilev,iwav) = fa(icol,ilev,iwav) + dopaer(icol)*palb(icol)*pasm(icol)*pasm(icol)

                      call update_diags

                   end do column

                end do vertical
             end do wavelength

          else
             call endrun(prefix//'aero_optics object pointer not associated')
          end if

          deallocate(aero_optics)
          nullify(aero_optics)

          call output_bin_diags

       end do binloop
    end do aeromodel

    call output_tot_diags

    deallocate(pext)
    deallocate(pabs)
    deallocate(palb)
    deallocate(pasm)

    do iaermod = 1,num_aero_models
       deallocate(aero_state(iaermod)%obj)
       nullify(aero_state(iaermod)%obj)
    end do

    deallocate(aero_state)

  contains

    !===============================================================================
    subroutine init_diags
      dustvol(:ncol)   = 0._r8
      scatdust(:ncol)  = 0._r8
      absdust(:ncol)   = 0._r8
      hygrodust(:ncol) = 0._r8
      scatsulf(:ncol)  = 0._r8
      abssulf(:ncol)   = 0._r8
      hygrosulf(:ncol) = 0._r8
      scatbc(:ncol)    = 0._r8
      absbc(:ncol)     = 0._r8
      hygrobc(:ncol)   = 0._r8
      scatpom(:ncol)   = 0._r8
      abspom(:ncol)    = 0._r8
      hygropom(:ncol)  = 0._r8
      scatsoa(:ncol)   = 0._r8
      abssoa(:ncol)    = 0._r8
      hygrosoa(:ncol)  = 0._r8
      scatsslt(:ncol)  = 0._r8
      abssslt(:ncol)   = 0._r8
      hygrosslt(:ncol) = 0._r8
    end subroutine init_diags

    !===============================================================================
    subroutine update_diags

      integer :: ispec

      if (iwav==idx_uv_diag) then
         aoduv(icol) = aoduv(icol) + dopaer(icol)
         extinctuv(icol,ilev) = extinctuv(icol,ilev) + dopaer(icol)*air_density(icol,ilev)/mass(icol,ilev)
         if (ilev<=troplev(icol)) then
            aoduvst(icol) = aoduvst(icol) + dopaer(icol)
         end if

      else if (iwav==idx_sw_diag) then ! vis
         aodvis(icol) = aodvis(icol) + dopaer(icol)
         aodabs(icol) = aodabs(icol) + pabs(icol)*mass(icol,ilev)
         extinct(icol,ilev) = extinct(icol,ilev) + dopaer(icol)*air_density(icol,ilev)/mass(icol,ilev)
         absorb(icol,ilev)  = absorb(icol,ilev) + pabs(icol)*air_density(icol,ilev)
         ssavis(icol)       = ssavis(icol) + dopaer(icol)*palb(icol)
         asymvis(icol)      = asymvis(icol) + dopaer(icol)*pasm(icol)
         asymext(icol,ilev) = asymext(icol,ilev) + dopaer(icol)*pasm(icol)*air_density(icol,ilev)/mass(icol,ilev)

         aodbin(icol) = aodbin(icol) + dopaer(icol)

         if (ilev<=troplev(icol)) then
            aodvisst(icol) = aodvisst(icol) + dopaer(icol)
         end if

         ! loop over species ...

         do ispec = 1, aeroprops%nspecies(list_idx,ibin)
            call aeroprops%get(ibin, ispec, list_ndx=list_idx, density=specdens, &
                 spectype=spectype, refindex_sw=specrefindex, hygro=hygro_aer)
            call aerostate%get_ambient_mmr(list_idx, ispec, ibin, specmmr)

            burden(icol) = burden(icol) + specmmr(icol,ilev)*mass(icol,ilev)

            vol(icol) = specmmr(icol,ilev)/specdens

            select case ( trim(spectype) )
            case('dust')
               dustvol(icol) = vol(icol)
               burdendust(icol) = burdendust(icol) + specmmr(icol,ilev)*mass(icol,ilev)
               scatdust(icol) = vol(icol) * specrefindex(iwav)%re
               absdust(icol)  =-vol(icol) * specrefindex(iwav)%im
               hygrodust(icol)= vol(icol)*hygro_aer
            case('black-c')
               burdenbc(icol) = burdenbc(icol) + specmmr(icol,ilev)*mass(icol,ilev)
               scatbc(icol) = vol(icol) * specrefindex(iwav)%re
               absbc(icol)  =-vol(icol) * specrefindex(iwav)%im
               hygrobc(icol)= vol(icol)*hygro_aer
            case('sulfate')
               burdenso4(icol) = burdenso4(icol) + specmmr(icol,ilev)*mass(icol,ilev)
               scatsulf(icol) = vol(icol) * specrefindex(iwav)%re
               abssulf(icol)  =-vol(icol) * specrefindex(iwav)%im
               hygrosulf(icol)= vol(icol)*hygro_aer
            case('p-organic')
               burdenpom(icol) = burdenpom(icol) + specmmr(icol,ilev)*mass(icol,ilev)
               scatpom(icol) = vol(icol) * specrefindex(iwav)%re
               abspom(icol)  =-vol(icol) * specrefindex(iwav)%im
               hygropom(icol)= vol(icol)*hygro_aer
            case('s-organic')
               burdensoa(icol) = burdensoa(icol) + specmmr(icol,ilev)*mass(icol,ilev)
               scatsoa(icol) = vol(icol) * specrefindex(iwav)%re
               abssoa(icol) = -vol(icol) * specrefindex(iwav)%im
               hygrosoa(icol)= vol(icol)*hygro_aer
            case('seasalt')
               burdenseasalt(icol) = burdenseasalt(icol) + specmmr(icol,ilev)*mass(icol,ilev)
               scatsslt(icol) = vol(icol) * specrefindex(iwav)%re
               abssslt(icol) = -vol(icol) * specrefindex(iwav)%im
               hygrosslt(icol)= vol(icol)*hygro_aer
            end select
         end do

         if (wetvol(icol,ilev)>1.e-40_r8 .and. vol(icol)>0._r8) then

            dustaodbin(icol) = dustaodbin(icol) + dopaer(icol)*dustvol(icol)/wetvol(icol,ilev)

            ! partition optical depth into contributions from each constituent
            ! assume contribution is proportional to refractive index X volume

            scath2o = watervol(icol,ilev)*crefwsw(iwav)%re
            absh2o = -watervol(icol,ilev)*crefwsw(iwav)%im
            sumscat = scatsulf(icol) + scatpom(icol) + scatsoa(icol) + scatbc(icol) + &
                 scatdust(icol) + scatsslt(icol) + scath2o
            sumabs  = abssulf(icol) + abspom(icol) + abssoa(icol) + absbc(icol) + &
                 absdust(icol) + abssslt(icol) + absh2o
            sumhygro = hygrosulf(icol) + hygropom(icol) + hygrosoa(icol) + hygrobc(icol) + &
                 hygrodust(icol) + hygrosslt(icol)

            scatdust(icol) = (scatdust(icol) + scath2o*hygrodust(icol)/sumhygro)/sumscat
            absdust(icol)  = (absdust(icol) + absh2o*hygrodust(icol)/sumhygro)/sumabs

            scatsulf(icol) = (scatsulf(icol) + scath2o*hygrosulf(icol)/sumhygro)/sumscat
            abssulf(icol)  = (abssulf(icol) + absh2o*hygrosulf(icol)/sumhygro)/sumabs

            scatpom(icol) = (scatpom(icol) + scath2o*hygropom(icol)/sumhygro)/sumscat
            abspom(icol)  = (abspom(icol) + absh2o*hygropom(icol)/sumhygro)/sumabs

            scatsoa(icol) = (scatsoa(icol) + scath2o*hygrosoa(icol)/sumhygro)/sumscat
            abssoa(icol)  = (abssoa(icol) + absh2o*hygrosoa(icol)/sumhygro)/sumabs

            scatbc(icol)= (scatbc(icol) + scath2o*hygrobc(icol)/sumhygro)/sumscat
            absbc(icol)  = (absbc(icol) +  absh2o*hygrobc(icol)/sumhygro)/sumabs

            scatsslt(icol) = (scatsslt(icol) + scath2o*hygrosslt(icol)/sumhygro)/sumscat
            abssslt(icol)  = (abssslt(icol) + absh2o*hygrosslt(icol)/sumhygro)/sumabs


            aodabsbc(icol) = aodabsbc(icol) + absbc(icol)*dopaer(icol)*(1.0_r8-palb(icol))



            aodc          = (absdust(icol)*(1.0_r8 - palb(icol)) + palb(icol)*scatdust(icol))*dopaer(icol)
            dustaod(icol) = dustaod(icol) + aodc

            aodc          = (abssulf(icol)*(1.0_r8 - palb(icol)) + palb(icol)*scatsulf(icol))*dopaer(icol)
            sulfaod(icol) = sulfaod(icol) + aodc

            aodc          = (abspom(icol)*(1.0_r8 - palb(icol)) + palb(icol)*scatpom(icol))*dopaer(icol)
            pomaod(icol)  = pomaod(icol) + aodc

            aodc          = (abssoa(icol)*(1.0_r8 - palb(icol)) + palb(icol)*scatsoa(icol))*dopaer(icol)
            soaaod(icol)  = soaaod(icol) + aodc

            aodc          = (absbc(icol)*(1.0_r8 - palb(icol)) + palb(icol)*scatbc(icol))*dopaer(icol)
            bcaod(icol)   = bcaod(icol) + aodc

            aodc          = (abssslt(icol)*(1.0_r8 - palb(icol)) + palb(icol)*scatsslt(icol))*dopaer(icol)
            ssltaod(icol) = ssltaod(icol) + aodc

         end if
      else if (iwav==idx_nir_diag) then
         aodnir(icol) = aodnir(icol) + dopaer(icol)
         extinctnir(icol,ilev) = extinctnir(icol,ilev) + dopaer(icol)*air_density(icol,ilev)/mass(icol,ilev)

         if (ilev<=troplev(icol)) then
            aodnirst(icol) = aodnirst(icol) + dopaer(icol)
         end if

      end if

      aodtot(icol) = aodtot(icol) + dopaer(icol)

    end subroutine update_diags

    !===============================================================================
    subroutine output_bin_diags

      integer :: icol

      if (list_idx == 0) then

         call outfld(burdendn_fields(iaermod)%name(ibin), burden, pcols, lchnk)
         call outfld(aoddustdn_fields(iaermod)%name(ibin), dustaodbin, pcols, lchnk)
         call outfld(aodbindn_fields(iaermod)%name(ibin),      aodbin, pcols, lchnk)

         do icol = 1, nnite
            burden(idxnite(icol))  = fillvalue
            aodbin(idxnite(icol)) = fillvalue
            dustaodbin(idxnite(icol)) = fillvalue
         end do

         call outfld(burden_fields(iaermod)%name(ibin), burden, pcols, lchnk)
         call outfld(aoddust_fields(iaermod)%name(ibin), dustaodbin, pcols, lchnk)
         call outfld(aodbin_fields(iaermod)%name(ibin),      aodbin, pcols, lchnk)

      endif

    end subroutine output_bin_diags

    !===============================================================================
    subroutine output_tot_diags

      integer :: icol

      call outfld('AODUVdn'//diag(list_idx),  aoduv,  pcols, lchnk)
      call outfld('AODVISdn'//diag(list_idx), aodvis, pcols, lchnk)
      call outfld('AODABSdn'//diag(list_idx),     aodabs,  pcols, lchnk)
      call outfld('AODNIRdn'//diag(list_idx), aodnir, pcols, lchnk)
      call outfld('AODTOTdn'//diag(list_idx), aodtot, pcols, lchnk)
      call outfld('EXTINCTUVdn'//diag(list_idx),  extinctuv,  pcols, lchnk)
      call outfld('EXTINCTNIRdn'//diag(list_idx), extinctnir, pcols, lchnk)
      call outfld('EXTINCTdn'//diag(list_idx),  extinct,  pcols, lchnk)
      call outfld('ABSORBdn'//diag(list_idx),   absorb,  pcols, lchnk)
      call outfld('EXTxASYMdn'//diag(list_idx), asymext, pcols, lchnk)
      call outfld('AODxASYMdn'//diag(list_idx), asymvis, pcols, lchnk)

      call outfld('AODVISstdn'//diag(list_idx), aodvisst,pcols, lchnk)
      call outfld('AODUVstdn'//diag(list_idx),  aoduvst, pcols, lchnk)
      call outfld('AODNIRstdn'//diag(list_idx), aodnirst,pcols, lchnk)

      do icol = 1, nnite
         aodvis(idxnite(icol)) = fillvalue
         aodnir(idxnite(icol)) = fillvalue
         aoduv(idxnite(icol)) = fillvalue
         aodabs(idxnite(icol)) = fillvalue
         aodtot(idxnite(icol)) = fillvalue
         extinct(idxnite(icol),:) = fillvalue
         extinctnir(idxnite(icol),:) = fillvalue
         extinctuv(idxnite(icol),:) = fillvalue
         absorb(idxnite(icol),:)  = fillvalue
         asymext(idxnite(icol),:) = fillvalue
         asymvis(idxnite(icol)) = fillvalue
         aodabs(idxnite(icol))    = fillvalue
         aodvisst(idxnite(icol))  = fillvalue
         aoduvst(idxnite(icol))   = fillvalue
         aodnirst(idxnite(icol))  = fillvalue
      end do

      call outfld('AODUV'//diag(list_idx),  aoduv,  pcols, lchnk)
      call outfld('AODVIS'//diag(list_idx), aodvis, pcols, lchnk)
      call outfld('AODABS'//diag(list_idx), aodabs,  pcols, lchnk)
      call outfld('AODNIR'//diag(list_idx), aodnir, pcols, lchnk)
      call outfld('AODTOT'//diag(list_idx), aodtot, pcols, lchnk)
      call outfld('EXTINCTUV'//diag(list_idx),  extinctuv,  pcols, lchnk)
      call outfld('EXTINCTNIR'//diag(list_idx), extinctnir, pcols, lchnk)
      call outfld('EXTINCT'//diag(list_idx),  extinct,  pcols, lchnk)
      call outfld('ABSORB'//diag(list_idx),   absorb,  pcols, lchnk)
      call outfld('EXTxASYM'//diag(list_idx), asymext, pcols, lchnk)
      call outfld('AODxASYM'//diag(list_idx), asymvis, pcols, lchnk)
      call outfld('AODVISst'//diag(list_idx), aodvisst,pcols, lchnk)
      call outfld('AODUVst'//diag(list_idx),  aoduvst, pcols, lchnk)
      call outfld('AODNIRst'//diag(list_idx), aodnirst,pcols, lchnk)

      ! These diagnostics are output only for climate list
      if (list_idx == 0) then
         do icol = 1, ncol
            if (aodvis(icol) > 1.e-10_r8) then
               ssavis(icol) = ssavis(icol)/aodvis(icol)
            else
               ssavis(icol) = 0.925_r8
            endif
         end do
         call outfld('SSAVISdn',        ssavis,        pcols, lchnk)

         call outfld('BURDENDUSTdn',    burdendust,    pcols, lchnk)
         call outfld('BURDENSO4dn' ,    burdenso4,     pcols, lchnk)
         call outfld('BURDENPOMdn' ,    burdenpom,     pcols, lchnk)
         call outfld('BURDENSOAdn' ,    burdensoa,     pcols, lchnk)
         call outfld('BURDENBCdn'  ,    burdenbc,      pcols, lchnk)
         call outfld('BURDENSEASALTdn', burdenseasalt, pcols, lchnk)

         call outfld('AODABSBCdn',      aodabsbc,      pcols, lchnk)

         call outfld('AODDUSTdn',       dustaod,       pcols, lchnk)
         call outfld('AODSO4dn',        sulfaod,       pcols, lchnk)
         call outfld('AODPOMdn',        pomaod,        pcols, lchnk)
         call outfld('AODSOAdn',        soaaod,        pcols, lchnk)
         call outfld('AODBCdn',         bcaod,         pcols, lchnk)
         call outfld('AODSSdn',         ssltaod,       pcols, lchnk)


         do icol = 1, nnite

            ssavis(idxnite(icol))     = fillvalue
            asymvis(idxnite(icol))    = fillvalue

            burdendust(idxnite(icol)) = fillvalue
            burdenso4(idxnite(icol))  = fillvalue
            burdenpom(idxnite(icol))  = fillvalue
            burdensoa(idxnite(icol))  = fillvalue
            burdenbc(idxnite(icol))   = fillvalue
            burdenseasalt(idxnite(icol)) = fillvalue
            aodabsbc(idxnite(icol)) = fillvalue

            dustaod(idxnite(icol))    = fillvalue
            sulfaod(idxnite(icol))     = fillvalue
            pomaod(idxnite(icol))     = fillvalue
            soaaod(idxnite(icol))     = fillvalue
            bcaod(idxnite(icol))      = fillvalue
            ssltaod(idxnite(icol)) = fillvalue

         end do

         call outfld('SSAVIS',        ssavis,        pcols, lchnk)
         call outfld('AODxASYM',      asymvis,       pcols, lchnk)
         call outfld('BURDENDUST',    burdendust,    pcols, lchnk)
         call outfld('BURDENSO4' ,    burdenso4,     pcols, lchnk)
         call outfld('BURDENPOM' ,    burdenpom,     pcols, lchnk)
         call outfld('BURDENSOA' ,    burdensoa,     pcols, lchnk)
         call outfld('BURDENBC'  ,    burdenbc,      pcols, lchnk)
         call outfld('BURDENSEASALT', burdenseasalt, pcols, lchnk)
         call outfld('AODABSBC',      aodabsbc,      pcols, lchnk)
         call outfld('AODDUST',       dustaod,       pcols, lchnk)
         call outfld('AODSO4',        sulfaod,       pcols, lchnk)
         call outfld('AODPOM',        pomaod,        pcols, lchnk)
         call outfld('AODSOA',        soaaod,        pcols, lchnk)
         call outfld('AODBC',         bcaod,         pcols, lchnk)
         call outfld('AODSS',         ssltaod,       pcols, lchnk)

      end if

    end subroutine output_tot_diags

  end subroutine aerosol_optics_cam_sw

  !===============================================================================
  subroutine aerosol_optics_cam_lw(list_idx, state, pbuf, tauxar)

    ! calculates aerosol lw radiative properties

    integer,             intent(in)  :: list_idx ! index of the climate or a diagnostic list
    type(physics_state), intent(in), target :: state    ! state variables

    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8), intent(inout) :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth


    real(r8) :: dopaer(pcols)
    real(r8) :: mass(pcols,pver)

    character(len=*), parameter :: prefix = 'aerosol_optics_cam_lw: '

    integer :: ibin, nbins
    integer :: iwav, ilev
    integer :: ncol, icol, istat

    type(aero_state_t), allocatable :: aero_state(:) ! array of aerosol state objects to allow for
                                                     ! multiple aerosol representations in the same sim
                                                     ! such as MAM and CARMA

    class(aerosol_optics), pointer :: aero_optics
    class(aerosol_state),      pointer :: aerostate
    class(aerosol_properties), pointer :: aeroprops

    real(r8), allocatable :: pabs(:)

    character(len=32) :: opticstype
    integer :: iaermod

    real(r8) :: lwabs(pcols,pver)
    lwabs = 0._r8
    tauxar = 0._r8

    nullify(aero_optics)

    allocate(aero_state(num_aero_models), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: aero_state')
    end if

    iaermod = 0
    if (modal_active) then
       iaermod = iaermod+1
       aero_state(iaermod)%obj => modal_aerosol_state( state, pbuf )
    end if

    ncol = state%ncol

    mass(:ncol,:) = state%pdeldry(:ncol,:)*rga

    allocate(pabs(ncol), stat=istat)
    if (istat/=0) then
       call endrun(prefix//'array allocation error: pabs')
    end if

    aeromodel: do iaermod = 1,num_aero_models

       aeroprops => aero_props(iaermod)%obj
       aerostate => aero_state(iaermod)%obj

       nbins=aero_props(iaermod)%obj%nbins(list_idx)

       binloop: do ibin = 1, nbins

          call aeroprops%optics_params(list_idx, ibin, opticstype=opticstype)

          select case (trim(opticstype))
          case('modal') ! refractive method
             aero_optics=>refractive_aerosol_optics(aeroprops, aerostate, list_idx, ibin, &
                                                    ncol, pver, nswbands, nlwbands, crefwsw, crefwlw)
          case default
             call endrun(prefix//'optics method not recognized')
          end select

          if (associated(aero_optics)) then

             wavelength: do iwav = 1, nlwbands

                vertical: do ilev = 1, pver
                   call aero_optics%lw_props(ncol, ilev, iwav, pabs )

                   column: do icol = 1, ncol
                      dopaer(icol) = pabs(icol)*mass(icol,ilev)
                      tauxar(icol,ilev,iwav) = tauxar(icol,ilev,iwav) + dopaer(icol)
                      lwabs(icol,ilev) = lwabs(icol,ilev) + pabs(icol)
                   end do column

                end do vertical

             end do wavelength

          else
             call endrun(prefix//'aero_optics object pointer not associated')
          end if

          deallocate(aero_optics)
          nullify(aero_optics)

       end do binloop
    end do aeromodel

    call outfld('TOTABSLW'//diag(list_idx),  lwabs(:,:), pcols, state%lchnk)

    if (lw10um_indx>0) then
       call outfld('AODABSLW'//diag(list_idx), tauxar(:,:,lw10um_indx), pcols, state%lchnk)
    end if

    deallocate(pabs)

    do iaermod = 1,num_aero_models
       deallocate(aero_state(iaermod)%obj)
       nullify(aero_state(iaermod)%obj)
    end do

    deallocate(aero_state)

  end subroutine aerosol_optics_cam_lw

  !===============================================================================
  ! Private routines
  !===============================================================================

  subroutine read_water_refindex(infilename)
    use cam_pio_utils, only: cam_pio_openfile
    use pio, only: file_desc_t, var_desc_t, pio_inq_dimlen, pio_inq_dimid, pio_inq_varid, &
                   pio_get_var, PIO_NOWRITE, pio_closefile, pio_noerr


    ! read water refractive index file and set module data

    character*(*), intent(in) :: infilename   ! modal optics filename

    ! Local variables

    integer            :: i, ierr
    type(file_desc_t)  :: ncid              ! pio file handle
    integer            :: did               ! dimension ids
    integer            :: dimlen            ! dimension lengths
    type(var_desc_t)   :: vid               ! variable ids
    real(r8) :: refrwsw(nswbands), refiwsw(nswbands) ! real, imaginary ref index for water visible
    real(r8) :: refrwlw(nlwbands), refiwlw(nlwbands) ! real, imaginary ref index for water infrared

    character(len=*), parameter :: prefix = 'read_water_refindex: '
    !----------------------------------------------------------------------------

    ! open file
    call cam_pio_openfile(ncid, infilename, PIO_NOWRITE)

    ! inquire dimensions.  Check that file values match parameter values.

    ierr = pio_inq_dimid(ncid, 'lw_band', did)
    if (ierr /= pio_noerr ) then
       call endrun(prefix//'pio_inq_dimid lw_band')
    end if
    ierr = pio_inq_dimlen(ncid, did, dimlen)
    if (ierr /= pio_noerr ) then
       call endrun(prefix//'pio_inq_dimlen lw_band')
    end if
    if (dimlen /= nlwbands) then
       write(iulog,*) 'lw_band len=', dimlen, ' from ', infilename, ' ne nlwbands=', nlwbands
       call endrun(prefix//'bad lw_band value')
    endif

    ierr = pio_inq_dimid(ncid, 'sw_band', did)
    if (ierr /= pio_noerr ) then
       call endrun(prefix//'pio_inq_dimid sw_band')
    end if
    ierr = pio_inq_dimlen(ncid, did, dimlen)
    if (ierr /= pio_noerr ) then
       call endrun(prefix//'pio_inq_dimlen sw_band')
    end if
    if (dimlen /= nswbands) then
       write(iulog,*) 'sw_band len=', dimlen, ' from ', infilename, ' ne nswbands=', nswbands
       call endrun(prefix//'bad sw_band value')
    endif

    ! read variables
    ierr = pio_inq_varid(ncid, 'refindex_real_water_sw', vid)
    if (ierr /= pio_noerr ) then
       call endrun(prefix//'pio_inq_varid refindex_real_water_sw')
    end if
    ierr = pio_get_var(ncid, vid, refrwsw)
    if (ierr /= pio_noerr ) then
       call endrun(prefix//'pio_get_var refrwsw')
    end if

    ierr = pio_inq_varid(ncid, 'refindex_im_water_sw', vid)
    if (ierr /= pio_noerr ) then
       call endrun(prefix//'pio_inq_varid refindex_im_water_sw')
    end if
    ierr = pio_get_var(ncid, vid, refiwsw)
    if (ierr /= pio_noerr ) then
       call endrun(prefix//'pio_get_var refiwsw')
    end if

    ierr = pio_inq_varid(ncid, 'refindex_real_water_lw', vid)
    if (ierr /= pio_noerr ) then
       call endrun(prefix//'pio_inq_varid refindex_real_water_lw')
    end if
    ierr = pio_get_var(ncid, vid, refrwlw)
    if (ierr /= pio_noerr ) then
       call endrun(prefix//'pio_get_var refrwlw')
    end if

    ierr = pio_inq_varid(ncid, 'refindex_im_water_lw', vid)
    if (ierr /= pio_noerr ) then
       call endrun(prefix//'pio_inq_varid refindex_im_water_lw')
    end if
    ierr = pio_get_var(ncid, vid, refiwlw)
    if (ierr /= pio_noerr ) then
       call endrun(prefix//'pio_get_var refiwlw')
    end if

    ! set complex representation of refractive indices as module data
    do i = 1, nswbands
       crefwsw(i) = cmplx(refrwsw(i), abs(refiwsw(i)), kind=r8)
    end do
    do i = 1, nlwbands
       crefwlw(i) = cmplx(refrwlw(i), abs(refiwlw(i)), kind=r8)
    end do

    call pio_closefile(ncid)

  end subroutine read_water_refindex

end module aerosol_optics_cam
