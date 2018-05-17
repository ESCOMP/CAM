!========================================================================
! NOy at upper boundary for CAM-Chem
!========================================================================

module noy_ubc

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use spmd_utils,       only : masterproc
  use cam_abortutils,   only : endrun
  use cam_logfile,      only : iulog

  use tracer_data,      only : trfld,trfile,MAXTRCRS
  use cam_history,      only : addfld, horiz_only

  implicit none

  private
  public :: noy_ubc_init
  public :: noy_ubc_set
  public :: noy_ubc_advance
  public :: noy_ubc_readnl

  save

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  integer :: ub_nspecies
  character(len=16) :: ubc_name(MAXTRCRS)
  integer :: map(MAXTRCRS) = -1

  character(len=256) :: noy_ubc_filename  = 'NONE'
  character(len=256) :: noy_ubc_filelist  = ' '
  character(len=256) :: noy_ubc_datapath  = ' '
  character(len=32)  :: noy_ubc_datatype  = 'SERIAL'
  logical            :: noy_ubc_rmv_file  = .false.
  integer            :: noy_ubc_cycle_yr  = 0
  integer            :: noy_ubc_fixed_ymd = 0
  integer            :: noy_ubc_fixed_tod = 0

  real(r8)           :: fac_relax

  logical :: has_noy_ubc = .false.

contains
  
  !======================================================================
  !======================================================================
  subroutine noy_ubc_readnl(nlfile)
    
    use namelist_utils, only : find_group_name
    use units,          only : getunit, freeunit
    use spmd_utils,     only : mpicom, masterprocid, mpi_character, mpi_integer

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    
    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'noy_ubc_readnl'

    namelist /noy_ubc_nl/ &
         noy_ubc_filename, noy_ubc_filelist, noy_ubc_datapath, noy_ubc_datatype, &
         noy_ubc_cycle_yr, noy_ubc_fixed_ymd, noy_ubc_fixed_tod

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'noy_ubc_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, noy_ubc_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(noy_ubc_filename,  len(noy_ubc_filename), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(noy_ubc_filelist,  len(noy_ubc_filelist), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(noy_ubc_datapath,  len(noy_ubc_datapath), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(noy_ubc_datatype,  len(noy_ubc_datatype), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(noy_ubc_cycle_yr,  1, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(noy_ubc_fixed_ymd, 1, mpi_integer, masterprocid, mpicom, ierr)
    call mpi_bcast(noy_ubc_fixed_tod, 1, mpi_integer, masterprocid, mpicom, ierr)

    has_noy_ubc = len_trim(noy_ubc_filename) > 0 .and. noy_ubc_filename.ne.'NONE'

  end subroutine noy_ubc_readnl

  !======================================================================
  !======================================================================
  subroutine noy_ubc_init()

    !------------------------------------------------------------------
    !	... initialize upper boundary values
    !------------------------------------------------------------------
    use tracer_data,  only : trcdata_init
    use mo_chem_utls, only : get_spc_ndx

    !------------------------------------------------------------------
    !	... dummy args
    !------------------------------------------------------------------

    ! local vars
    integer :: vid, i,ii 

    integer,          parameter :: nubc = 4
    character(len=4), parameter :: species(nubc) = (/'NO  ','NO2 ','HNO3','N2O5'/)
    character(len=4)            :: specifier(nubc) = ' '

    if (.not.has_noy_ubc) return

    ii = 0

    do i = 1,nubc
       vid = get_spc_ndx(species(i))
       if( vid > 0 ) then
          ii = ii+1
          specifier(ii) = species(i) ! set specifier to the species that actually 
                                     ! are in the simulation so that the species mapping is correct
          map(ii) = vid
          ubc_name(ii) = trim(specifier(i))//'_ubc'
          call addfld( ubc_name(ii), horiz_only, 'I', 'mol/mol', 'upper boundary vmr' )

       end if
    enddo

    ub_nspecies = count( map(:)>0 )
    
    if (ub_nspecies > 0) then
       file%top_bndry = .true.
       allocate(file%in_pbuf(size(specifier)))
       file%in_pbuf(:) = .false.
       call trcdata_init( specifier, noy_ubc_filename, noy_ubc_filelist, noy_ubc_datapath, fields, file, &
            noy_ubc_rmv_file, noy_ubc_cycle_yr, noy_ubc_fixed_ymd, noy_ubc_fixed_tod, noy_ubc_datatype)
    endif

  end subroutine noy_ubc_init

  !======================================================================
  !======================================================================
  subroutine noy_ubc_advance(pbuf2d, state)

    use tracer_data,    only : advance_trcdata
    use physics_types,  only : physics_state
    use physics_buffer, only : physics_buffer_desc
    use time_manager,   only : get_step_size

    !--------------------------------------------------------------------
    !	... Advance ub values
    !--------------------------------------------------------------------
    implicit none

    ! args
    type(physics_state),    intent(in) :: state(:)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
!
    integer             :: dtime                         ! model time step (s)
    real(r8), parameter :: tau_relax = 864000._r8        ! 10 days

    if (.not.has_noy_ubc) return
!
! define relaxation factor
!
    dtime = get_step_size()
    fac_relax = 1._r8 - exp( -real(dtime)/tau_relax )
!
    if (ub_nspecies > 0) then
       call advance_trcdata( fields, file, state, pbuf2d )
    endif

  end subroutine noy_ubc_advance

  !======================================================================
  !	... Set the upper boundary values
  !======================================================================
  subroutine noy_ubc_set( lchnk, ncol, vmr )
    use cam_history,  only : outfld

    implicit none

    !--------------------------------------------------------------------
    !	... dummy args
    !--------------------------------------------------------------------
    integer,  intent(in)    :: lchnk            ! chunk id
    integer,  intent(in)    :: ncol              ! columns in chunk
    real(r8), intent(inout) :: vmr(:,:,:)

    integer  :: m,n,m1,m2,i
    real(r8) :: xno,xno2,xnox,rno,dtime
    real(r8) :: yno,yno2,ynox

    if (.not.has_noy_ubc) return
!
! only update model top layer (index=1)
!
    if (ub_nspecies > 0) then
       do m = 1,ub_nspecies
          if ( trim(fields(m)%fldnam) == 'NO' .or. trim(fields(m)%fldnam) == 'NO2' ) cycle
          n = map(m)
          vmr(:ncol,1,n) = fields(m)%data(:ncol,1,lchnk)
          call outfld( ubc_name(m), vmr(:ncol,1,n), ncol, lchnk )
       enddo
    endif
!
! special case for NO & NO2
!
    m1 = -99
    m2 = -99
    do m=1,ub_nspecies
       if ( trim(fields(m)%fldnam) == 'NO'  ) m1 = m
       if ( trim(fields(m)%fldnam) == 'NO2' ) m2 = m
    end do
    if ( m1 > 0 .and. m2 > 0 ) then
!
       do i=1,ncol
!
          xno  = vmr(i,1,map(m1))
          xno2 = vmr(i,1,map(m2))
          xnox = xno + xno2
          rno  = xno/xnox
!
          yno  = fields(m1)%data(i,1,lchnk)
          yno2 = fields(m2)%data(i,1,lchnk)
          ynox = yno + yno2
!
! relax model NOx towards the specified values
!
          xnox = xnox + (ynox - xnox) * fac_relax
!
! use original ratio to redistribute updated NOx between NO and NO2
!
          vmr(i,1,map(m1)) = rno * xnox
          vmr(i,1,map(m2)) = (1._r8-rno) * xnox
!
       end do
!
       call outfld( ubc_name(m1), vmr(:ncol,1,map(m1)), ncol, lchnk )
       call outfld( ubc_name(m2), vmr(:ncol,1,map(m2)), ncol, lchnk )
    end if
!
  end subroutine noy_ubc_set

end module noy_ubc
