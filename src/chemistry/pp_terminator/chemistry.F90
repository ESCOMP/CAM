!================================================================================================
! This is the "toy" chemistry module.
!================================================================================================

module chemistry
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use physics_types,       only: physics_state, physics_ptend, physics_ptend_init
  use ppgrid,              only: begchunk, endchunk, pcols
  use ppgrid,              only: pver
  use constituents,        only: pcnst, cnst_add
  use mo_gas_phase_chemdr, only: map2chm
  use mo_constants,        only: pi
  use shr_const_mod,       only: molw_dryair=>SHR_CONST_MWDAIR
  use mo_chem_utls,        only : get_spc_ndx
  use chem_mods,           only : gas_pcnst, adv_mass
  use mo_sim_dat, only: set_sim_dat
  implicit none
  private
  save
  !
  ! Public interfaces
  !
  public :: chem_is                        ! identify which chemistry is being used
  public :: chem_register                  ! register consituents
  public :: chem_is_active                 ! returns true if this package is active (ghg_chem=.true.)
  public :: chem_implements_cnst           ! returns true if consituent is implemented by this package
  public :: chem_init_cnst                 ! initialize mixing ratios if not read from initial file
  public :: chem_init                      ! initialize (history) variables
  public :: chem_timestep_init             ! time interpolate chemical loss frequencies
  public :: chem_timestep_tend             ! interface to tendency computation
  public :: chem_final
  public :: chem_write_restart
  public :: chem_read_restart
  public :: chem_init_restart
  public :: chem_readnl                    ! read chem namelist 
  public :: chem_reset_fluxes
  public :: chem_emissions

  interface chem_write_restart
     module procedure chem_write_restart_bin
     module procedure chem_write_restart_pio
  end interface
  interface chem_read_restart
     module procedure chem_read_restart_bin
     module procedure chem_read_restart_pio
  end interface

  ! Private data
  integer, parameter :: nspecies = 3
  
  integer :: idx_cl =-1
  integer :: idx_cl2=-1

  real(r8), parameter :: k1_lat_center =  pi* 20.0_r8/180._r8
  real(r8), parameter :: k1_lon_center =  pi*300.0_r8/180._r8

  character(len=8) :: species(nspecies) = (/ 'CL      ','CL2     ','RHO     '/)
  integer          :: indices(nspecies)

!================================================================================================
contains
!================================================================================================

  logical function chem_is (name)

    character(len=*), intent(in) :: name

    chem_is = .false.
    if (name == 'toy' ) then
       chem_is = .true.
    end if

  end function chem_is

!================================================================================================

  subroutine chem_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents for parameterized greenhouse gas chemistry
    ! 
    !-----------------------------------------------------------------------

    real(r8), parameter :: cptmp = 666._r8
    real(r8), parameter :: qmin = -1.e36_r8

    logical :: camout
    integer :: i, n
    
    do i = 1, nspecies
       camout = trim(species(i)) .eq. 'RHO'
       call cnst_add( species(i), adv_mass(i), cptmp, qmin, n, &  
                      readiv=.true.,mixtype='dry',cam_outfld=camout)
       indices(i) = n
       map2chm(n) = i
    enddo

    call set_sim_dat()

    idx_cl = get_spc_ndx('CL')
    idx_cl2 = get_spc_ndx('CL2')

  end subroutine chem_register

!================================================================================================

  subroutine chem_readnl(nlfile)

    ! args
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  end subroutine chem_readnl

!================================================================================================

  function chem_is_active()
    !-----------------------------------------------------------------------
    logical :: chem_is_active
    !-----------------------------------------------------------------------
    chem_is_active = .true.
  end function chem_is_active

!================================================================================================

  function chem_implements_cnst(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is implemented by this package
    ! 
    ! Author: B. Eaton
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: chem_implements_cnst        ! return value

    integer :: i
    
    chem_implements_cnst = .false.

    do i = 1, nspecies
       if (trim(species(i)) .eq. trim(name)) then
          chem_implements_cnst = .true.
          exit
       end if
    end do

  end function chem_implements_cnst

!===============================================================================
  
  subroutine chem_init(phys_state, pbuf2d)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: initialize parameterized greenhouse gas chemistry
    !          (declare history variables)
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer, only: physics_buffer_desc
    use cam_history,    only: addfld, add_default, horiz_only

    type(physics_state), intent(in):: phys_state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer :: i

    do i = 1, nspecies
       call addfld ( species(i), (/'lev'/), 'A', 'mole/mole', trim(species(i))//' mixing ratio' )
       call add_default ( species(i),   1, ' ')
    enddo

    call addfld( 'k1', horiz_only, 'A', '/s ', 'reaction rate ' )
    call addfld( 'k2', horiz_only, 'A', '/s ', 'reaction rate ' )
    call add_default ( 'k1',   1, ' ')
    call add_default ( 'k2',   1, ' ')

    call addfld ( 'CLy', (/'lev'/), 'A', 'mole/mole', 'CLy mixing ratio' )
    call add_default ( 'CLy',   1, ' ')
    !
    ! terminator diagnostics (DCMIP2016)
    !
    call addfld ( 'iCLy', horiz_only, 'A','mole/mole','Average mass-weighted column-integrated CLy mixing ratio')
    call add_default ( 'iCLy',   1, ' ')

    call addfld ( 'iCL', horiz_only, 'A','mole/mole','Average mass-weighted column-integrated CL mixing ratio')
    call add_default ( 'iCL',   1, ' ')

    call addfld ( 'iCL2', horiz_only, 'A','mole/mole','Average mass-weighted column-integrated CL2 mixing ratio')
    call add_default ( 'iCL2',   1, ' ')
  end subroutine chem_init

!===============================================================================

  subroutine chem_timestep_init(phys_state, pbuf2d)
    use physics_buffer,   only: physics_buffer_desc

    type(physics_state), intent(in):: phys_state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

  end subroutine chem_timestep_init

!===============================================================================

  subroutine chem_timestep_tend( state, ptend, cam_in, cam_out, dt, pbuf,  fh2o )

    use physics_buffer,   only: physics_buffer_desc
    use cam_history,      only: outfld
    use camsrfexch,       only: cam_in_t, cam_out_t
    use cam_history,      only: hist_fld_active
    !-----------------------------------------------------------------------
    !
    ! Arguments:
    !
    real(r8),            intent(in)    :: dt          ! time step
    type(physics_state), intent(in)    :: state       ! Physics state variables
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(in)    :: cam_out
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), optional,  intent(out)   :: fh2o(pcols) ! h2o flux to balance source from chemistry

    real(r8) :: a(pver),b(pver),c(pver),d(pver)
    
    real(r8) :: k1(pcols)
    real(r8) :: k2(pcols)

    real(r8) :: cly(pcols,pver)
    real(r8) :: cl (pcols,pver)
    real(r8) :: cl2(pcols,pver)
    real(r8) :: new_cl (pcols,pver)
    real(r8) :: new_cl2(pcols,pver)

    real(r8) :: integrated(pcols)

    integer :: ncol,lchnk, i
    logical :: lq(pcnst)
    integer :: m,n

    real(r8) :: r(pcols), det(pcols,pver), e(pcols,pver)
    real(r8) :: cl_f(pcols,pver),cl2_f(pcols,pver)
    real(r8) :: l(pcols,pver)

    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          lq(n) = .true.
       end if
    end do

    call physics_ptend_init(ptend, state%psetcols, 'chemistry', lq=lq)

    lchnk = state%lchnk
    ncol = state%ncol

    do i = 1,ncol
       k1(i) = max(0.0_r8, sin(state%lat(i))*sin(k1_lat_center) + &
                           cos(state%lat(i))*cos(k1_lat_center)*cos(state%lon(i)-k1_lon_center) )
    enddo
    k2(:) = 1._r8

    call outfld('k1', k1, pcols, lchnk )
    call outfld('k2', k2, pcols, lchnk )

    cl(:ncol,:)  = (molw_dryair/adv_mass(idx_cl ))*state%q(:ncol,:,indices(1))
    cl2(:ncol,:) = (molw_dryair/adv_mass(idx_cl2))*state%q(:ncol,:,indices(2))

    do i = 1,ncol

        r(i) = k1(i)/(4._r8*k2(i))
        cly(i,:) = cl(i,:) + 2._r8* cl2(i,:)

        det(i,:) = sqrt(r(i)*r(i) + 2._r8*r(i)*cly(i,:))
        e(i,:) = exp(-4._r8*k2(i)*det(i,:)*dt)

        where( abs(det(i,:) * k2(i) * dt) .gt. 1e-16_r8)
           l(i,:) = (1._r8 - e(i,:))/det(i,:)/dt
        elsewhere
           l(i,:) = 4._r8*k2(i)
        endwhere 

        cl_f(i,:) = -l(i,:)*(cl(i,:) - det(i,:) + r(i) )*(cl(i,:) + det(i,:) + r(i)) / ( 1._r8 +e(i,:) + dt*l(i,:)*(cl(i,:) + r(i)))
        cl2_f(i,:) = -cl_f(i,:) / 2._r8

    enddo

    ptend%q(:ncol,:,indices(1)) = (adv_mass(idx_cl )/molw_dryair)*( cl_f (:ncol,:) )
    ptend%q(:ncol,:,indices(2)) = (adv_mass(idx_cl2)/molw_dryair)*( cl2_f(:ncol,:) )

    cly(:ncol,:) = 2._r8*(cl2(:ncol,:) + dt * cl2_f(:ncol,:)) + (cl(:ncol,:) + dt * cl_f(:ncol,:))

    call outfld('CLy', cly, pcols, lchnk )
    call outfld ( species(1), cl (:ncol,:) + dt * cl_f (:ncol,:), ncol, lchnk )
    call outfld ( species(2), cl2(:ncol,:) + dt * cl2_f(:ncol,:), ncol, lchnk )
    !
    ! terminator diagnostics (DCMIP2016)
    !
    if ( hist_fld_active('iCL')) then
      integrated(:ncol) = SUM(state%pdeldry(:ncol,:)*cl(:ncol,:),DIM=2)/SUM(state%pdeldry(:ncol,:),DIM=2)
      call outfld('iCL', integrated, pcols, lchnk )
    end if
    if ( hist_fld_active('iCL2')) then
      integrated(:ncol) = SUM(state%pdeldry(:ncol,:)*cl2(:ncol,:),DIM=2)/SUM(state%pdeldry(:ncol,:),DIM=2)
      call outfld('iCL2', integrated, pcols, lchnk )
    end if
    if ( hist_fld_active('iCLy')) then
      integrated(:ncol) = SUM(state%pdeldry(:ncol,:)*cly(:ncol,:),DIM=2)/SUM(state%pdeldry(:ncol,:),DIM=2)
      call outfld('iCLy', integrated, pcols, lchnk )
    end if
    return
  end subroutine chem_timestep_tend

!===============================================================================
  subroutine chem_init_cnst(name, latvals, lonvals, mask, q)

!-----------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------
    character(len=*), intent(in)  :: name       !  constituent name
    real(r8),         intent(in)  :: latvals(:) ! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) ! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)     ! kg tracer/kg dry air (ncol, plev
    real(r8)            :: q_vmr(size(q, 1)) !  volume mixing ratio (ncol)
    real(r8)            :: det(size(q, 1))
    real(r8)            :: krat(size(q, 1))
    
    real(r8)            :: k1(size(q, 1))
    real(r8)            :: k2(size(q, 1))

    real(r8)            :: h
    integer             :: i, k
    integer             :: blksize

    real(r8), parameter :: init_vmr_cl2 = 2e-6_r8
    real(r8), parameter :: init_vmr_cl  = 0._r8

    blksize = size(q, 1)

    do i = 1, blksize
      k1(i) = max(0.0_r8, sin(latvals(i))*sin(k1_lat_center) + &
                          cos(latvals(i))*cos(k1_lat_center)*cos(lonvals(i)-k1_lon_center))
    end do
    k2(:) = 1._r8

    krat(:) = k1(:) / (4._r8 * k2(:))

    h = init_vmr_cl + 2._r8 * init_vmr_cl2
    
    det(:) = sqrt(krat(:) * krat(:) + 2._r8 * h * krat(:))

    if (trim(name) == trim(species(1)) ) then
      do k = 1, pver
        q_vmr(:) = (det(:)-krat(:))
        where(mask)
          q(:,k) = q_vmr(:) * adv_mass(idx_cl ) / molw_dryair
        end where
      end do
    else if (trim(name) == trim(species(2))) then
      do k = 1, pver
        q_vmr(:) = h / 2._r8 - (det(:) - krat(:)) / 2._r8
        where(mask)
          q(:,k) = q_vmr(:) * adv_mass(idx_cl2) / molw_dryair
        end where
      end do
    else if (trim(name) == trim(species(3))) then
      do k = 1, pver
        where(mask)
          q(:,k) = 1.0_r8
        end where
      end do
    end if

    return

  end subroutine chem_init_cnst

!===============================================================================
  subroutine chem_final
    return
  end subroutine chem_final
!===============================================================================
  subroutine chem_write_restart_bin( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    return
  end subroutine chem_write_restart_bin
!===============================================================================
  subroutine chem_read_restart_bin( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    return
  end subroutine chem_read_restart_bin
!===============================================================================
  subroutine chem_write_restart_pio( File )
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_write_restart_pio
!===============================================================================
  subroutine chem_read_restart_pio( File )
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_read_restart_pio
!===============================================================================
  subroutine chem_init_restart(File)
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_init_restart
!================================================================================
  subroutine chem_reset_fluxes( fptr, cam_in )
    use camsrfexch, only : cam_in_t     

    real(r8), pointer             :: fptr(:,:)        ! pointer into    array data
    type(cam_in_t), intent(inout) :: cam_in(begchunk:endchunk)

  end subroutine chem_reset_fluxes
!================================================================================
  subroutine chem_emissions( state, cam_in )
    use camsrfexch,       only: cam_in_t     

    ! Arguments:

    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t),         intent(inout) :: cam_in  ! import state

  end subroutine chem_emissions

end module chemistry
