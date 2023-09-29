module upper_bc

!---------------------------------------------------------------------------------
! Module to compute the upper boundary conditions for temperature (dry static energy)
! and trace gases. Uses the MSIS model, and SNOE and TIME GCM and general prescribed UBC data.
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_kind_mod, only: cl => shr_kind_cl
  use shr_const_mod,only: grav   => shr_const_g,     &   ! gravitational constant (m/s^2)
                          kboltz => shr_const_boltz, &   ! Boltzmann constant
                          pi => shr_const_pi,        &   ! pi
                          rEarth => shr_const_rearth     ! Earth radius
  use ppgrid,       only: pcols, pver, pverp
  use constituents, only: pcnst
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc
  use ref_pres,     only: do_molec_diff, ptop_ref
  use shr_kind_mod, only: cx=>SHR_KIND_CX
  use cam_abortutils,only: endrun
  use cam_history,   only: addfld, horiz_only, outfld, fieldname_len

  use upper_bc_file, only: upper_bc_file_readnl, upper_bc_file_specified, upper_bc_file_adv, upper_bc_file_get
  use infnan,        only: nan, assignment(=)

  implicit none
  private
  save
!
! Public interfaces
!
  public :: ubc_readnl         ! read namelist options for UBCs
  public :: ubc_init           ! global initialization
  public :: ubc_timestep_init  ! time step initialization
  public :: ubc_get_vals       ! get ubc values for this step
  public :: ubc_get_flxs       ! get ub fluxes for this step
  public :: ubc_fixed_conc     ! returns true for constituents that have fixed UBC
  public :: ubc_fixed_temp     ! true if temperature at upper boundary is fixed

  character(len=64) :: ubc_specifier(pcnst) = 'NOTSET'

  character(len=16) :: ubc_flds(pcnst) = 'NOTSET'
  character(len=16) :: ubc_file_spfr(pcnst) = ' '
  character(len=32) :: ubc_source(pcnst) = ' '

  integer :: n_fixed_mmr=0
  integer :: n_fixed_vmr=0

  integer, allocatable :: fixed_mmr_ndx(:)
  integer, allocatable :: fixed_vmr_ndx(:)
  real(r8), allocatable :: fixed_mmr(:)
  real(r8), allocatable :: fixed_vmr(:)

  integer :: num_infile = 0
  integer :: num_fixed = 0
  character(len=2), parameter :: msis_flds(5) = &
       (/ 'H ','N ','O ','O2','T ' /)
  character(len=2), parameter :: tgcm_flds(1) = &
       (/ 'H2' /)
  character(len=2), parameter :: snoe_flds(1) = &
       (/ 'NO' /)

  logical, protected :: ubc_fixed_temp =.false.
  logical :: msis_active =.false.
  logical :: tgcm_active =.false.
  logical :: snoe_active =.false.

! Namelist variables
  character(len=cl) :: snoe_ubc_file = 'NONE'
  real(r8)          :: t_pert_ubc  = 0._r8
  real(r8)          :: no_xfac_ubc = 1._r8

  integer :: h_ndx=-1
  integer :: h_msis_ndx=-1, n_msis_ndx=-1, o_msis_ndx=-1, o2_msis_ndx=-1

  character(len=cl) :: tgcm_ubc_file = 'NONE'
  integer           :: tgcm_ubc_cycle_yr = 0
  integer           :: tgcm_ubc_fixed_ymd = 0
  integer           :: tgcm_ubc_fixed_tod = 0
  character(len=32) :: tgcm_ubc_data_type = 'CYCLICAL'

  logical :: apply_upper_bc = .false.

  integer, allocatable :: file_spc_ndx(:)
  integer, allocatable :: spc_ndx(:)
  character(len=fieldname_len), allocatable :: hist_names(:)

!================================================================================================
contains
!================================================================================================

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine ubc_readnl(nlfile)
    use namelist_utils, only : find_group_name
    use spmd_utils, only : mpicom, masterprocid, mpi_character, mpi_integer, mpi_real8
    use string_utils, only : to_lower

    character(len=*), intent(in) :: nlfile
    integer :: unitn, ierr, m, n, ndx_co, ndx_ar

    character(len=*), parameter :: prefix = 'ubc_readnl: '

    namelist /upper_bc_opts/ tgcm_ubc_file,tgcm_ubc_data_type,tgcm_ubc_cycle_yr,tgcm_ubc_fixed_ymd, &
                             tgcm_ubc_fixed_tod, snoe_ubc_file, no_xfac_ubc, t_pert_ubc
    namelist /upper_bc_opts/ ubc_specifier

    if (masterproc) then
       ! read namelist
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'upper_bc_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, upper_bc_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun(prefix//'upper_bc_opts: ERROR reading namelist')
          end if
       end if
       close(unitn)

       ! log the UBC options
       write(iulog,*) prefix//'tgcm_ubc_file = '//trim(tgcm_ubc_file)
       write(iulog,*) prefix//'tgcm_ubc_data_type = '//trim(tgcm_ubc_data_type)
       write(iulog,*) prefix//'tgcm_ubc_cycle_yr = ', tgcm_ubc_cycle_yr
       write(iulog,*) prefix//'tgcm_ubc_fixed_ymd = ', tgcm_ubc_fixed_ymd
       write(iulog,*) prefix//'tgcm_ubc_fixed_tod = ', tgcm_ubc_fixed_tod
       write(iulog,*) prefix//'snoe_ubc_file = '//trim(snoe_ubc_file)
       write(iulog,*) prefix//'t_pert_ubc = ', t_pert_ubc
       write(iulog,*) prefix//'no_xfac_ubc = ', no_xfac_ubc
       write(iulog,*) prefix//'ubc_specifier : '

       n=1
       m=1
       do while(ubc_specifier(n)/='NOTSET')
          write(iulog,'(i4,a)') n,'  '//trim(ubc_specifier(n))

          ndx_ar = index(ubc_specifier(n),'->')

          if (ndx_ar<1) then
             call endrun(prefix//'ubc_specifier "'//trim(ubc_specifier(n))//'" must include "->"')
          endif

          ubc_source(n) = trim(to_lower(adjustl(ubc_specifier(n)(ndx_ar+2:))))

          if (trim(ubc_source(n))=='ubc_file') then
             ubc_file_spfr(m) = trim(ubc_specifier(n)(:ndx_ar-1))
             m=m+1
          endif
          if (index(ubc_source(n),'mmr')>0) then
             n_fixed_mmr=n_fixed_mmr+1
          else if (index(ubc_source(n),'vmr')>0) then
             n_fixed_vmr=n_fixed_vmr+1
          end if

          ndx_co = index(ubc_specifier(n),':')

          if (ndx_co>0) then
             ubc_flds(n) = ubc_specifier(n)(:ndx_co-1)
          else
             ubc_flds(n) = ubc_specifier(n)(:ndx_ar-1)
          end if

          n=n+1
       end do
       num_fixed=n-1
       num_infile=m-1
    end if


    ! broadcast to all MPI tasks
    call mpi_bcast(num_fixed, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : num_fixed')
    call mpi_bcast(num_infile, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : num_infile')
    call mpi_bcast(n_fixed_mmr, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : n_fixed_mmr')
    call mpi_bcast(n_fixed_vmr, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : n_fixed_vmr')
    call mpi_bcast(tgcm_ubc_file, len(tgcm_ubc_file), mpi_character, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : tgcm_ubc_file')
    call mpi_bcast(tgcm_ubc_data_type, len(tgcm_ubc_data_type),mpi_character, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : tgcm_ubc_data_type')
    call mpi_bcast(tgcm_ubc_cycle_yr, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : tgcm_ubc_cycle_yr')
    call mpi_bcast(tgcm_ubc_fixed_ymd, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : tgcm_ubc_fixed_ymd')
    call mpi_bcast(tgcm_ubc_fixed_tod, 1, mpi_integer, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : tgcm_ubc_fixed_tod')
    call mpi_bcast(snoe_ubc_file, len(snoe_ubc_file), mpi_character, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : snoe_ubc_file')
    call mpi_bcast(t_pert_ubc, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : t_pert_ubc')
    call mpi_bcast(no_xfac_ubc,1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : no_xfac_ubc')
    call mpi_bcast(ubc_specifier, pcnst*len(ubc_specifier(1)), mpi_character,masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : ubc_specifier')
    call mpi_bcast(ubc_flds,      pcnst*len(ubc_flds(1)),      mpi_character,masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : ubc_flds')
    call mpi_bcast(ubc_file_spfr, pcnst*len(ubc_file_spfr(1)), mpi_character,masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : ubc_file_spfr')
    call mpi_bcast(ubc_source,    pcnst*len(ubc_source(1)),    mpi_character,masterprocid, mpicom, ierr)
    if (ierr /= 0) call endrun(prefix//'mpi_bcast error : ubc_source')

    apply_upper_bc = num_fixed>0

    if (apply_upper_bc) then
       call upper_bc_file_readnl(nlfile)
    end if

  end subroutine ubc_readnl

!===============================================================================

  subroutine ubc_init()
!-----------------------------------------------------------------------
! Initialization of time independent fields for the upper boundary condition
! Calls initialization routine for MSIS, TGCM and SNOE
!-----------------------------------------------------------------------
    use mo_tgcm_ubc, only: tgcm_ubc_inti
    use mo_snoe,     only: snoe_inti
    use mo_msis_ubc, only: msis_ubc_inti
    use constituents,only: cnst_get_ind
    use upper_bc_file, only: upper_bc_file_init

    !---------------------------Local workspace-----------------------------
    logical, parameter :: zonal_avg = .false.
    integer :: m, mm, ierr
    integer :: mmrndx, vmrndx, m_mmr, m_vmr

    real(r8) :: val
    character(len=32) :: str
    character(len=*), parameter :: prefix = 'ubc_init: '
    !-----------------------------------------------------------------------

    call cnst_get_ind('H', h_ndx, abort=.false.) ! for H fluxes UBC (WACCMX)

    if (.not.apply_upper_bc) return

    if (num_infile>0) then
       call upper_bc_file_init( ubc_file_spfr(:num_infile) )
    endif

    ! possible MSIS, TGCM, SNOE and ubc_file inputs

    mm=1

    allocate(hist_names(num_fixed), stat=ierr)
    if (ierr /= 0) call endrun(prefix//'allocate error : hist_names')
    allocate(spc_ndx(num_fixed), stat=ierr)
    if (ierr /= 0) call endrun(prefix//'allocate error : spc_ndx')
    spc_ndx=-1
    if (num_infile>0) then
       allocate(file_spc_ndx(num_infile), stat=ierr)
       if (ierr /= 0) call endrun(prefix//'allocate error : file_spc_ndx')
       file_spc_ndx=-1
    end if
    if (n_fixed_mmr>0) then
       allocate(fixed_mmr_ndx(n_fixed_mmr), stat=ierr)
       if (ierr /= 0) call endrun(prefix//'allocate error : fixed_mmr_ndx')
       allocate(fixed_mmr(n_fixed_mmr), stat=ierr)
       if (ierr /= 0) call endrun(prefix//'allocate error : fixed_mmr')
       fixed_mmr_ndx=-1
       fixed_mmr = nan
    end if
    if (n_fixed_vmr>0) then
       allocate(fixed_vmr_ndx(n_fixed_vmr), stat=ierr)
       if (ierr /= 0) call endrun(prefix//'allocate error : fixed_vmr_ndx')
       allocate(fixed_vmr(n_fixed_vmr), stat=ierr)
       if (ierr /= 0) call endrun(prefix//'allocate error : fixed_vmr')
       fixed_vmr_ndx=-1
       fixed_vmr = nan
    end if

    m_mmr = 0
    m_vmr = 0

    do m = 1,num_fixed
       hist_names(m) = trim(ubc_flds(m))//'_UBC'
       if (ubc_flds(m)=='T') then
          ubc_fixed_temp=.true.
          spc_ndx(m) = -1
          call addfld(hist_names(m), horiz_only, 'I', 'K', trim(ubc_flds(m))//' at upper boundary' )
       else
          call cnst_get_ind(ubc_flds(m), spc_ndx(m), abort=.true.)
          call addfld(hist_names(m), horiz_only, 'I', 'kg/kg', trim(ubc_flds(m))//' at upper boundary' )
       end if

       if (trim(ubc_source(m))=='msis') then
          if (do_molec_diff .and. any(msis_flds==ubc_flds(m))) then
             msis_active = .true.
             if (trim(ubc_flds(m))=='H') h_msis_ndx=spc_ndx(m)
             if (trim(ubc_flds(m))=='N') n_msis_ndx=spc_ndx(m)
             if (trim(ubc_flds(m))=='O') o_msis_ndx=spc_ndx(m)
             if (trim(ubc_flds(m))=='O2') o2_msis_ndx=spc_ndx(m)
          else
             call endrun(prefix//'MSIS is not allowed in this configuration')
          end if
       else if (trim(ubc_source(m))=='tgcm') then
          if (do_molec_diff .and. any(tgcm_flds==ubc_flds(m))) then
             tgcm_active = .true.
          else
             call endrun(prefix//'TGCM is not allowed in this configuration')
          end if
       else if (trim(ubc_source(m))=='snoe') then
          if (do_molec_diff .and. any(snoe_flds==ubc_flds(m))) then
             snoe_active = .true.
          else
             call endrun(prefix//'SNOE is not allowed in this configuration')
          end if
       else if (trim(ubc_source(m))=='ubc_file') then
          file_spc_ndx(mm) = spc_ndx(m)
          mm = mm+1
       else
          mmrndx = index(trim(ubc_source(m)),'mmr')
          vmrndx = index(trim(ubc_source(m)),'vmr')
          if (mmrndx>0 .and. vmrndx>0) then
             call endrun(prefix//'incorrect units in UBC source: '//trim(ubc_source(m)))
          end if
          if (mmrndx>0) then
             str = ubc_source(m)(:mmrndx-1)
             read(str,*) val
             m_mmr = m_mmr + 1
             fixed_mmr(m_mmr) = val
             fixed_mmr_ndx(m_mmr) = spc_ndx(m)
          else if (vmrndx>0) then
             str = ubc_source(m)(:vmrndx-1)
             read(str,*) val
             m_vmr = m_vmr + 1
             fixed_vmr(m_vmr) = val
             fixed_vmr_ndx(m_vmr) = spc_ndx(m)
          else
             call endrun(prefix//'unrecognized UBC source: '//trim(ubc_source(m)))
          end if
       end if
    end do

    if (tgcm_active) then
       !-----------------------------------------------------------------------
       !       ... initialize the tgcm upper boundary module
       !-----------------------------------------------------------------------
       call tgcm_ubc_inti( tgcm_ubc_file, tgcm_ubc_data_type, tgcm_ubc_cycle_yr, &
            tgcm_ubc_fixed_ymd, tgcm_ubc_fixed_tod)
       if (masterproc) write(iulog,*) 'ubc_init: after tgcm_ubc_inti'
    endif

    if (snoe_active) then
       !-----------------------------------------------------------------------
       !       ... initialize the snoe module
       !-----------------------------------------------------------------------
       call snoe_inti(snoe_ubc_file)
       if (masterproc) write(iulog,*) 'ubc_init: after snoe_inti'
    endif

    if (msis_active) then
       !-----------------------------------------------------------------------
       !       ... initialize the msis module
       !-----------------------------------------------------------------------
       call msis_ubc_inti( zonal_avg, n_msis_ndx,h_msis_ndx,o_msis_ndx,o2_msis_ndx )
       if (masterproc) write(iulog,*) 'ubc_init: after msis_ubc_inti'
    endif

  end subroutine ubc_init

!===============================================================================
!===============================================================================

  pure logical function ubc_fixed_conc(name)

    character(len=*), intent(in) :: name

    integer :: m

    ubc_fixed_conc = .false.

    do m = 1,num_fixed
       if ( trim(ubc_flds(m)) == trim(name) ) then
          ubc_fixed_conc = .true.
          return
       endif
    end do

  end function ubc_fixed_conc

!===============================================================================

  subroutine ubc_timestep_init(pbuf2d, state)
!-----------------------------------------------------------------------
! timestep dependent setting
!-----------------------------------------------------------------------

    use solar_parms_data, only: kp=>solar_parms_kp, ap=>solar_parms_ap, f107=>solar_parms_f107
    use solar_parms_data, only: f107a=>solar_parms_f107a, f107p=>solar_parms_f107p
    use mo_msis_ubc,      only: msis_timestep_init
    use mo_tgcm_ubc,      only: tgcm_timestep_init
    use mo_snoe,          only: snoe_timestep_init
    use physics_types,    only: physics_state
    use ppgrid,           only: begchunk, endchunk
    use physics_buffer,   only: physics_buffer_desc

    type(physics_state), intent(in) :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if (.not.apply_upper_bc) return

    if (num_infile>0) then
       call upper_bc_file_adv( pbuf2d, state )
    end if
    if (msis_active) then
       call msis_timestep_init( ap, f107p, f107a )
    end if
    if (tgcm_active) then
       call tgcm_timestep_init( pbuf2d, state )
    end if
    if (snoe_active) then
       call snoe_timestep_init( kp, f107 )
    end if

  end subroutine ubc_timestep_init

!===============================================================================

  subroutine ubc_get_vals (lchnk, ncol, pint, zi, ubc_temp, ubc_mmr)

!-----------------------------------------------------------------------
! interface routine for vertical diffusion and pbl scheme
!-----------------------------------------------------------------------
    use mo_msis_ubc,      only: get_msis_ubc
    use mo_snoe,          only: set_no_ubc, ndx_no
    use mo_tgcm_ubc,      only: set_tgcm_ubc
    use cam_abortutils,   only: endrun
    use air_composition,  only: rairv, mbarv ! gas constant, mean mass
    use constituents,     only: cnst_mw  ! Needed for ubc_flux

!------------------------------Arguments--------------------------------
    integer,  intent(in)  :: lchnk                 ! chunk identifier
    integer,  intent(in)  :: ncol                  ! number of atmospheric columns
    real(r8), intent(in)  :: pint(pcols,pverp)     ! interface pressures
    real(r8), intent(in)  :: zi(pcols,pverp)       ! interface geoptl height above sfc

    real(r8), intent(out) :: ubc_temp(pcols)      ! upper bndy temperature (K)
    real(r8), intent(out) :: ubc_mmr(pcols,pcnst)  ! upper bndy mixing ratios (kg/kg)

    !---------------------------Local storage-------------------------------
    integer :: m                                   ! constituent index
    real(r8) :: rho_top(pcols)                     ! density at top interface
    real(r8) :: z_top(pcols)                       ! height of top interface (km)
    real(r8) :: vals(pcols,num_infile)

    real(r8), parameter :: m2km = 1.e-3_r8         ! meter to km

    !-----------------------------------------------------------------------

    ubc_mmr(:,:) = 0._r8
    ubc_temp(:) = nan

    if (.not. apply_upper_bc) return

    ! UBC_FILE
    if (num_infile>0) then
       call upper_bc_file_get(lchnk, ncol, vals)
       do m = 1,num_infile
          if (file_spc_ndx(m)>0) then
             ubc_mmr(:ncol,file_spc_ndx(m)) = vals(:ncol,m)
          else
             ubc_temp(:ncol) = vals(:ncol,m)
          end if

       end do

    endif

    ! MSIS
    if (msis_active) then
       call get_msis_ubc( lchnk, ncol, ubc_temp, ubc_mmr )
       if( t_pert_ubc /= 0._r8 ) then
          ubc_temp(:ncol) = ubc_temp(:ncol) + t_pert_ubc
          if( any( ubc_temp(:ncol) < 0._r8 ) ) then
             write(iulog,*) 'ubc_get_vals: msis temp < 0 after applying offset = ',t_pert_ubc
             call endrun('ubc_get_vals: msis temp < 0 after applying t_pert_ubc')
          end if
       end if
    end if

    ! SNOE
    if (snoe_active) then

       rho_top(:ncol) = pint(:ncol,1) / (rairv(:ncol,1,lchnk)*ubc_temp(:ncol))
       z_top(:ncol)   = m2km * zi(:ncol,1)

       call set_no_ubc  ( lchnk, ncol, z_top, ubc_mmr, rho_top )
       if( ndx_no > 0 .and. no_xfac_ubc /= 1._r8 ) then
          ubc_mmr(:ncol,ndx_no) = no_xfac_ubc * ubc_mmr(:ncol,ndx_no)
       end if

    endif

    ! TIE-GCM
    if (tgcm_active) then
       call set_tgcm_ubc( lchnk, ncol, ubc_mmr )
    endif

    ! fixed values
    do m = 1,n_fixed_mmr
       ubc_mmr(:ncol,fixed_mmr_ndx(m)) = fixed_mmr(m)
    end do
    do m = 1,n_fixed_vmr
       ubc_mmr(:ncol,fixed_vmr_ndx(m)) = cnst_mw(fixed_vmr_ndx(m))*fixed_vmr(m)/mbarv(:ncol,1,lchnk)
    end do

    ! diagnostic output
    do m = 1,num_fixed
       if (ubc_flds(m)=='T') then
          call outfld(hist_names(m),ubc_temp(:ncol),ncol,lchnk)
       else
          call outfld(hist_names(m),ubc_mmr(:ncol,spc_ndx(m)),ncol,lchnk)
       end if
    end do

  end subroutine ubc_get_vals

!===============================================================================

  subroutine ubc_get_flxs (lchnk, ncol, pint, zi, t, q, phis, ubc_flux)

    use physconst,       only: avogad, rga
    use air_composition, only: rairv
    use constituents,    only: cnst_mw
!------------------------------Arguments--------------------------------
    integer,  intent(in)  :: lchnk                 ! chunk identifier
    integer,  intent(in)  :: ncol                  ! number of atmospheric columns
    real(r8), intent(in)  :: pint(pcols,pverp)     ! interface pressures
    real(r8), intent(in)  :: zi(pcols,pverp)       ! interface geoptl height above sfc
    real(r8), intent(in)  :: t(pcols,pver)         ! midpoint temperature
    real(r8), intent(in),target :: q(pcols,pver,pcnst)   ! contituent mixing ratios (kg/kg)
    real(r8), intent(in)  :: phis(pcols)           ! Surface geopotential (m2/s2)

    real(r8), intent(out) :: ubc_flux(pcols,pcnst) ! upper bndy flux (kg/s/m^2)

!---------------------------Local storage-------------------------------
    integer :: iCol                                ! column loop counter

    real(r8), parameter :: h_escape_flx_factor = 2.03e-13_r8 ! for hydrogen escape flux due to charge exchange
    ! adopted from TIME-GCM (R. G. Roble, pp. 1-21, AGU Geophys. Monogr. Ser 87, 1995) following
    ! Liu, S.C., and T. M. Donahue, Mesospheric hydrogen related to exospheric escape mechanisms, J. Atmos. Sci.,
    ! 31, 1466-1470, 1974. (Equation 4 there). DOI: 10.1175/1520-0469(1974)031<1466:Mhrtee>2.0.Co;2
    ! https://journals.ametsoc.org/view/journals/atsc/31/5/1520-0469_1974_031_1466_mhrtee_2_0_co_2.xml

    real(r8), parameter :: hfluxlimitfac = 0.72_r8 ! Hydrogen upper boundary flux limiting factor

    real(r8) :: nmbartop                           ! Top level density (rho)
    real(r8) :: zkt                                ! Factor for H Jean's escape flux calculation

    real(r8), pointer :: qh_top(:)                 ! Top level hydrogen mixing ratio (kg/kg)

    ubc_flux(:,:) = nan

    qh_top => q(:,1,h_ndx)

    do iCol = 1, ncol
       !--------------------------------------------------
       ! Get total density (rho) at top level
       !--------------------------------------------------
       nmbartop = 0.5_r8 * (pint(iCol,1) + pint(iCol,2)) / ( rairv(iCol,1,lchnk) * t(iCol,1) )

       !---------------------------------------------------------------------
       ! Calculate factor for Jean's escape flux once here, used twice below
       !---------------------------------------------------------------------
       zkt = (rEarth + ( 0.5_r8 * ( zi(iCol,1) + zi(iCol,2) ) + rga * phis(iCol) ) ) * &
            cnst_mw(h_ndx) / avogad * grav / ( kboltz * t(iCol,1) )

       ubc_flux(iCol,h_ndx) = hfluxlimitfac * SQRT(kboltz/(2.0_r8 * pi * cnst_mw(h_ndx) / avogad)) * &
            qh_top(iCol) * nmbartop * &
            SQRT(t(iCol,1)) * (1._r8 + zkt) * EXP(-zkt)

       ubc_flux(iCol,h_ndx) = ubc_flux(iCol,h_ndx) * &
            (h_escape_flx_factor * qh_top(iCol) * nmbartop / (cnst_mw(h_ndx) / avogad) * t(iCol,1))

    enddo

  end subroutine ubc_get_flxs

end module upper_bc
