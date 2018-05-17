!-------------------------------------------------------------------------------
! Energetic Particle Precipitation (EPP) forcings module
!  Manages ionization of the atmosphere due to energetic particles, which consists of
!  solar protons events (SPE), galactic cosmic rays(GCR), medium energy electrons (MEE)
!-------------------------------------------------------------------------------
module epp_ionization
  use shr_kind_mod,   only : r8 => shr_kind_r8, cs => shr_kind_cs, cl=> shr_kind_cl
  use spmd_utils,     only : masterproc
  use cam_abortutils, only : endrun
  use cam_logfile,    only : iulog
  use phys_grid,      only : pcols, pver, begchunk, endchunk, get_ncols_p
  use pio,            only : var_desc_t, file_desc_t
  use pio,            only : pio_get_var, pio_inq_varid, pio_get_att
  use pio,            only : pio_inq_varndims, pio_inq_vardimid, pio_inq_dimname, pio_inq_dimlen
  use pio,            only : PIO_NOWRITE
  use cam_pio_utils,  only : cam_pio_openfile
  use ioFileMod,      only : getfil
  use input_data_utils, only : time_coordinate

  implicit none
  private

  public :: epp_ionization_readnl  ! read namelist variables
  public :: epp_ionization_init    ! initialization
  public :: epp_ionization_adv     ! read and time/space interpolate the data
  public :: epp_ionization_ionpairs! ion pairs production rates
  public :: epp_ionization_setmag  ! update geomagnetic coordinates mapping
  public :: epp_ionization_active

  character(len=cl) :: epp_all_filepath = 'NONE'
  character(len=cs) :: epp_all_varname  = 'epp_ion_rates'
  character(len=cl) :: epp_mee_filepath = 'NONE'
  character(len=cs) :: epp_mee_varname  = 'iprm'
  character(len=cl) :: epp_spe_filepath = 'NONE'
  character(len=cs) :: epp_spe_varname  = 'iprp'
  character(len=cl) :: epp_gcr_filepath = 'NONE'
  character(len=cs) :: epp_gcr_varname  = 'iprg'

  logical, protected :: epp_ionization_active = .false.

  type input_obj_t
    type(file_desc_t) :: fid
    type(var_desc_t)  :: vid
    character(len=32) :: units
    integer :: nlevs = 0
    integer :: nglats = 0
    real(r8), allocatable :: press(:)
    real(r8), allocatable :: glats(:)
    real(r8), allocatable :: gwght(:,:)      ! (pcol, begchunk:endchunk)
    integer,  allocatable :: glatn(:,:)      ! (pcol, begchunk:endchunk)
    real(r8), allocatable :: indata(:,:,:,:) ! (pcol,nlevs,begchunk:endchunk,2) inputs at indexm and indexp
    type(time_coordinate) :: time_coord
  endtype input_obj_t

  type(input_obj_t), pointer :: epp_in => null()
  type(input_obj_t), pointer :: spe_in => null()
  type(input_obj_t), pointer :: mee_in => null()
  type(input_obj_t), pointer :: gcr_in => null()

contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine epp_ionization_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
    use spmd_utils,     only: mpicom, mpi_character, masterprocid

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'epp_ionization_readnl'

    namelist /epp_ionization_nl/ epp_all_filepath, epp_all_varname, &
         epp_mee_filepath, epp_mee_varname, epp_spe_filepath, epp_spe_varname, epp_gcr_filepath, epp_gcr_varname

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'epp_ionization_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, epp_ionization_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(epp_all_filepath, len(epp_all_filepath), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(epp_mee_filepath, len(epp_mee_filepath), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(epp_spe_filepath, len(epp_spe_filepath), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(epp_gcr_filepath, len(epp_gcr_filepath), mpi_character, masterprocid, mpicom, ierr)

    call mpi_bcast(epp_all_varname, len(epp_all_varname), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(epp_mee_varname, len(epp_mee_varname), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(epp_spe_varname, len(epp_spe_varname), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(epp_gcr_varname, len(epp_gcr_varname), mpi_character, masterprocid, mpicom, ierr)

    epp_ionization_active = epp_all_filepath /= 'NONE'
    epp_ionization_active = epp_mee_filepath /= 'NONE' .or. epp_ionization_active
    epp_ionization_active = epp_spe_filepath /= 'NONE' .or. epp_ionization_active
    epp_ionization_active = epp_gcr_filepath /= 'NONE' .or. epp_ionization_active

    if ( epp_ionization_active .and. masterproc ) then
       write(iulog,*) subname//':: epp_all_filepath = '//trim(epp_all_filepath)
       write(iulog,*) subname//':: epp_mee_filepath = '//trim(epp_mee_filepath)
       write(iulog,*) subname//':: epp_spe_filepath = '//trim(epp_spe_filepath)
       write(iulog,*) subname//':: epp_gcr_filepath = '//trim(epp_gcr_filepath)
    endif

  end subroutine epp_ionization_readnl

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine epp_ionization_init()
    use cam_history, only : addfld

    character(len=32) :: fldunits
    fldunits = ''
    
    if (epp_all_filepath /= 'NONE') then
       epp_in => create_input_obj(epp_all_filepath,epp_all_varname)
       fldunits = trim(epp_in%units)
    else
       if (epp_mee_filepath /= 'NONE') then
          mee_in => create_input_obj(epp_mee_filepath,epp_mee_varname)
          fldunits = trim(mee_in%units)
       endif
       if (epp_spe_filepath /= 'NONE') then
          spe_in => create_input_obj(epp_spe_filepath,epp_spe_varname)
          fldunits = trim(spe_in%units)
       endif
       if (epp_gcr_filepath /= 'NONE') then
          gcr_in => create_input_obj(epp_gcr_filepath,epp_gcr_varname)
          fldunits = trim(gcr_in%units)
       endif
    endif
    call addfld( 'EPPions', (/ 'lev' /), 'A', fldunits, 'EPP ionization data' )

  end subroutine epp_ionization_init

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine epp_ionization_setmag( maglat )
    real(r8), intent(in) :: maglat(pcols,begchunk:endchunk)

    if (.not.epp_ionization_active) return

    if ( associated(epp_in) ) then
       call set_wghts(maglat,epp_in)
    else
       if ( associated(mee_in) ) then
          call set_wghts(maglat,mee_in)
       endif
       if ( associated(spe_in) ) then
          call set_wghts(maglat,spe_in)
       endif
       if ( associated(gcr_in) ) then
          call set_wghts(maglat,gcr_in)
       endif
    endif

  end subroutine epp_ionization_setmag

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine epp_ionization_adv

    if (.not.epp_ionization_active) return

    if ( associated(epp_in) ) then
       call update_input(epp_in)
    else
       if ( associated(spe_in) ) then
          call update_input(spe_in)
       endif
       if ( associated(gcr_in) ) then
          call update_input(gcr_in)
       endif
       if ( associated(mee_in) ) then
          call update_input(mee_in)
       endif
    endif

  end subroutine epp_ionization_adv

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine epp_ionization_ionpairs( ncol, lchnk, pmid, temp, ionpairs )

    integer,  intent(in) :: ncol, lchnk
    real(r8), intent(in) :: pmid(:,:), temp(:,:)
    real(r8), intent(out) :: ionpairs(:,:) ! ion pair production rate

    ionpairs = 0._r8
    if (.not.epp_ionization_active) return

    if ( associated(epp_in) ) then
       ionpairs(:ncol,:) = ionpairs(:ncol,:) + interp_ionpairs( ncol, lchnk, pmid, temp, epp_in )
    else 
       if ( associated(spe_in) ) then
          ionpairs(:ncol,:) = ionpairs(:ncol,:) + interp_ionpairs( ncol, lchnk, pmid, temp, spe_in )
       endif
       if ( associated(gcr_in) ) then
          ionpairs(:ncol,:) = ionpairs(:ncol,:) + interp_ionpairs( ncol, lchnk, pmid, temp, gcr_in )
       endif
       if ( associated(mee_in) ) then
          ionpairs(:ncol,:) = ionpairs(:ncol,:) + interp_ionpairs( ncol, lchnk, pmid, temp, mee_in )
       endif
    endif

  end subroutine epp_ionization_ionpairs

  ! private methods
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine update_input( input )
    type(input_obj_t), pointer :: input

    if ( input%time_coord%read_more() ) then
       call input%time_coord%advance()
       call read_next_data( input )
    else
       call input%time_coord%advance()
    endif

  end subroutine update_input

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine read_next_data( input )
    type(input_obj_t), pointer :: input

    ! read data corresponding surrounding time indices
    if ( input%nglats > 0 ) then
      call read_2d_profile( input )
    else
      call read_1d_profile( input )
    endif

  end subroutine read_next_data

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  function interp_ionpairs( ncol, lchnk, pmid, temp, input ) result( ionpairs )
    use interpolate_data, only : lininterp
    use physconst,        only : rairv
    use cam_history,      only : outfld

    integer,  intent(in) :: ncol, lchnk
    real(r8), intent(in) :: pmid(:,:) ! Pa
    real(r8), intent(in) :: temp(:,:) ! K
    type(input_obj_t), pointer :: input
    real(r8) :: ionpairs(ncol,pver)
 
    real(r8) :: fctr1, fctr2
    real(r8) :: wrk(ncol,input%nlevs)
    real(r8) :: ions_diags(ncol,pver) ! for diagnostics
    integer :: i

    if (input%time_coord%time_interp) then
       ! time interpolate
       fctr1 = input%time_coord%wghts(1)
       fctr2 = input%time_coord%wghts(2)
       wrk(:ncol,:) = fctr1*input%indata(:ncol,:,lchnk,1) + fctr2*input%indata(:ncol,:,lchnk,2)
    else
       wrk(:ncol,:) = input%indata(:ncol,:,lchnk,1)
    endif

    ! vertical interpolate ...
    ! interpolate to model levels
    do i = 1,ncol

       ! interpolate over log pressure
       call lininterp( wrk(i,:input%nlevs), log(input%press(:input%nlevs)*1.e2_r8), input%nlevs, &
                       ionpairs(i,:pver), log(pmid(i,:pver)), pver )
       ions_diags(i,:pver) = ionpairs(i,:pver)
       
       if ( index(trim(input%units), 'g^-1') > 0 ) then
          ! convert to ionpairs/cm3/sec
          ionpairs(i,:pver) = ionpairs(i,:pver) *(1.e-3_r8*pmid(i,:pver)/(rairv(i,:pver,lchnk)*temp(i,:pver)))
       endif
    enddo

    call outfld( 'EPPions', ions_diags(:ncol,:), ncol, lchnk )

  end function interp_ionpairs

  !-----------------------------------------------------------------------------
  ! read 2D profile (geomag-lat vs press) and transfer to geographic grid 
  !-----------------------------------------------------------------------------
  subroutine read_2d_profile( input )

    type(input_obj_t), pointer :: input

    ! local vars
    real(r8) :: wrk2d( input%nglats, input%nlevs, 2 )
    integer  :: t, c, i, ntimes, ncols, ierr
    real(r8) :: wght1, wght2
    integer  :: gndx1, gndx2
    integer  :: cnt(3), strt(3)

    if (input%time_coord%time_interp) then
       ntimes = 2
    else
       ntimes = 1
    endif

    cnt(1) = input%nglats
    cnt(2) = input%nlevs
    cnt(3) = ntimes

    strt(:) = 1
    strt(3) = input%time_coord%indxs(1)

    ierr = pio_get_var( input%fid, input%vid, strt, cnt, wrk2d )

    do t = 1,ntimes
       do c=begchunk,endchunk
          ncols = get_ncols_p(c)
          do i = 1,ncols
             gndx1 = input%glatn(i,c)
             if (gndx1>0) then
                wght1 = input%gwght(i,c)
                gndx2 = gndx1+1
                if (gndx2.le.input%nglats) then
                  wght2 = 1._r8-wght1
                  input%indata(i,:,c,t) = wght1*wrk2d(gndx1,:,t) &
                                        + wght2*wrk2d(gndx2,:,t)
                else
                  input%indata(i,:,c,t) = wght1*wrk2d(gndx1,:,t)
                endif
             else
                input%indata(i,:,c,t) = 0._r8
             endif
          end do
       end do
    end do

  end subroutine read_2d_profile

  !-----------------------------------------------------------------------------
  ! read 1D vertical profile and transfer to geographic grid poleward of 60 degrees geomag-lat
  !-----------------------------------------------------------------------------
  subroutine read_1d_profile( input )

    type(input_obj_t), pointer :: input

    ! local vars
    real(r8) :: wrk( input%nlevs, 2 )
    integer  :: t, c, i, ntimes, ncols, ierr
    integer  :: cnt(2), strt(2)

    if (input%time_coord%time_interp) then
       ntimes = 2
    else
       ntimes = 1
    endif

    cnt(1) = input%nlevs
    cnt(2) = ntimes

    strt(:) = 1
    strt(2) = input%time_coord%indxs(1)

    ierr = pio_get_var( input%fid, input%vid, strt, cnt, wrk )

    do t = 1,ntimes
       do c=begchunk,endchunk
          ncols = get_ncols_p(c)
          do i = 1,ncols
             input%indata(i,:,c,t) = input%gwght(i,c)*wrk(:,t)
          end do
       end do
    end do

  end subroutine read_1d_profile

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  function create_input_obj( path, varname ) result(in_obj)
    use infnan,  only : nan, assignment(=)

    character(*), intent(in) :: path
    character(*), intent(in) :: varname
    type(input_obj_t), pointer :: in_obj

    character(len=cl) :: filen
    character(len=cl) :: data_units
    character(len=cs) :: dimname
    integer :: i, ierr
    integer, allocatable :: dimids(:)
    integer :: pres_did, pres_vid, glat_did, glat_vid, ndims

    if (path .eq. 'NONE') return

    allocate(in_obj)

    call in_obj%time_coord%initialize( path )

    call getfil( path, filen, 0 )
    call cam_pio_openfile( in_obj%fid, filen, PIO_NOWRITE )

    ierr = pio_inq_varid( in_obj%fid, varname, in_obj%vid )

    ierr = pio_get_att( in_obj%fid, in_obj%vid, 'units', data_units)
    in_obj%units = trim(data_units(1:32))

    ierr = pio_inq_varndims( in_obj%fid, in_obj%vid, ndims )
    allocate( dimids(ndims) )

    ierr = pio_inq_vardimid( in_obj%fid, in_obj%vid, dimids)
    pres_did = -1
    glat_did = -1
    do i = 1,ndims
       ierr = pio_inq_dimname( in_obj%fid, dimids(i), dimname )
       select case( trim(dimname(1:4)) )
         case ( 'pres', 'lev', 'plev' )
           pres_did = dimids(i)
           ierr = pio_inq_varid( in_obj%fid, dimname, pres_vid)
         case ( 'glat' )
           glat_did =  dimids(i)
           ierr = pio_inq_varid( in_obj%fid, dimname, glat_vid)
         case default
       end select
    end do

    deallocate( dimids )

    if (pres_did>0) then
       ierr = pio_inq_dimlen( in_obj%fid, pres_did, in_obj%nlevs )
       allocate( in_obj%press(in_obj%nlevs) )
       ierr = pio_get_var( in_obj%fid, pres_vid, in_obj%press )
    endif
    if (glat_did>0) then 
       ierr = pio_inq_dimlen( in_obj%fid, glat_did, in_obj%nglats )
       allocate( in_obj%glats(in_obj%nglats) )
       ierr = pio_get_var( in_obj%fid, glat_vid, in_obj%glats )
       allocate( in_obj%glatn(pcols,begchunk:endchunk) )
    endif
       
    allocate( in_obj%gwght(pcols,begchunk:endchunk) )

    if (in_obj%time_coord%time_interp) then
       allocate( in_obj%indata(pcols,in_obj%nlevs,begchunk:endchunk,2) )
    else
       allocate( in_obj%indata(pcols,in_obj%nlevs,begchunk:endchunk,1) )
    endif
    in_obj%indata = nan

  end function create_input_obj

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine set_wghts( maglat, input )

    real(r8), intent(in) :: maglat(pcols,begchunk:endchunk)
    type(input_obj_t), pointer :: input

    integer :: i, c, ncols, imag

    if (input%nglats>1) then ! read in general EPP 2D ionpairs production rates
       do c = begchunk,endchunk
          ncols = get_ncols_p(c)
          col_loop: do i = 1,ncols
             if ( maglat(i,c) .lt.  input%glats(1) ) then
                input%glatn(i,c) = 1
                input%gwght(i,c) = 1._r8
             elseif ( maglat(i,c) .gt. input%glats(input%nglats) ) then
                input%glatn(i,c) = input%nglats
                input%gwght(i,c) = 1._r8
             else
                mag_loop: do imag = 1,input%nglats-1
                   if ( maglat(i,c) .ge. input%glats(imag) .and. &
                        maglat(i,c) .lt. input%glats(imag+1) ) then
                      input%gwght(i,c) = (input%glats(imag+1)-maglat(i,c) ) &
                                       / (input%glats(imag+1)-input%glats(imag))
                      input%glatn(i,c) = imag
                      exit mag_loop
                   endif
                enddo mag_loop
             endif
          enddo col_loop
       enddo
    else ! read in 1D SPE ionpairs profile ...
       do c = begchunk,endchunk
          ncols = get_ncols_p(c)
          do i = 1,ncols
             if ( abs(maglat(i,c)) .ge. 60._r8 ) then ! poleward of 60 degrees
                input%gwght(i,c) = 1._r8
             else
                input%gwght(i,c) = 0._r8
             endif
          enddo
       enddo
    endif

    call read_next_data( input ) ! update the inputs when wghts are updated

  end subroutine set_wghts

end module epp_ionization
