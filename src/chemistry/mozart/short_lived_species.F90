!---------------------------------------------------------------------
! Manages the storage of non-transported short-lived chemical species
! in the physics buffer.
!
! Created by: Francis Vitt -- 20 Aug 2008
!---------------------------------------------------------------------
module short_lived_species

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods,    only : slvd_lst, nslvd, gas_pcnst
  use cam_logfile,  only : iulog
  use ppgrid,       only : pcols, pver, begchunk, endchunk
  use spmd_utils,   only : masterproc
  

  implicit none

  save
  private
  public :: map
  public :: register_short_lived_species
  public :: short_lived_species_initic
  public :: short_lived_species_writeic
  public :: initialize_short_lived_species
  public :: set_short_lived_species
  public :: get_short_lived_species
  public :: slvd_index
  public :: pbf_idx

  integer :: pbf_idx
  integer :: map(nslvd)

  character(len=16), parameter :: pbufname = 'ShortLivedSpecies'

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine register_short_lived_species
    use physics_buffer, only : pbuf_add_field, dtype_r8

    implicit none

    integer :: m

    if ( nslvd < 1 ) return

    call pbuf_add_field(pbufname,'global',dtype_r8,(/pcols,pver,nslvd/),pbf_idx)

  end subroutine register_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine short_lived_species_initic
#ifdef WACCMX_IONOS
    use cam_history, only : addfld, add_default

    integer :: m
    character(len=24) :: varname

    do m=1,nslvd
       varname = trim(slvd_lst(m))//'&IC'
       call addfld (varname, (/ 'lev' /),'I','kg/kg',trim(varname)//' not-transported species',gridname='physgrid')
       call add_default (varname,0, 'I')
    enddo
#endif
  end subroutine short_lived_species_initic

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine short_lived_species_writeic( lchnk, pbuf )
    use cam_history,    only : outfld, write_inithist
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field

    integer       , intent(in) :: lchnk  ! chunk identifier
    type(physics_buffer_desc), pointer :: pbuf(:)
#ifdef WACCMX_IONOS
    real(r8),pointer :: tmpptr(:,:)
    integer :: m
    character(len=24) :: varname
    
    if ( write_inithist() ) then
       do m=1,nslvd
          varname = trim(slvd_lst(m))//'&IC'
          call pbuf_get_field(pbuf, pbf_idx, tmpptr, start=(/1,1,m/), kount=(/ pcols,pver,1 /))
          call outfld(varname, tmpptr, pcols,lchnk)
       enddo
    endif
#endif
  end subroutine short_lived_species_writeic

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine initialize_short_lived_species(ncid_ini, pbuf2d)
    use cam_grid_support, only : cam_grid_check, cam_grid_id
    use cam_grid_support, only : cam_grid_get_dim_names
    use cam_abortutils,   only : endrun
    use mo_tracname,      only : solsym
    use ncdio_atm,        only : infld
    use pio,              only : file_desc_t
    use physics_buffer,   only : physics_buffer_desc, pbuf_set_field, pbuf_get_chunk, pbuf_get_field

    implicit none

    type(file_desc_t), intent(inout) :: ncid_ini
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer          :: m,n,lchnk
    integer          :: grid_id
    character(len=8) :: fieldname
    character(len=4) :: dim1name, dim2name
    logical          :: found
    real(r8),pointer :: tmpptr(:,:,:)   ! temporary pointer
    real(r8),pointer :: tmpptr2(:,:,:)   ! temporary pointer
    character(len=*), parameter :: subname='INITIALIZE_SHORT_LIVED_SPECIES'

    if ( nslvd < 1 ) return

    found = .false.

    grid_id = cam_grid_id('physgrid')
    if (.not. cam_grid_check(grid_id)) then
      call endrun(trim(subname)//': Internal error, no "physgrid" grid')
    end if
    call cam_grid_get_dim_names(grid_id, dim1name, dim2name)

    call pbuf_set_field(pbuf2d, pbf_idx, 0._r8)

    allocate(tmpptr(pcols,pver,begchunk:endchunk))

    do m=1,nslvd
       n = map(m)
       fieldname = solsym(n)
       call infld( fieldname,ncid_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
                   tmpptr, found, gridname='physgrid')

       if (.not.found) then
          tmpptr(:,:,:) = 1.e-36_r8
       endif

       call pbuf_set_field(pbuf2d, pbf_idx, tmpptr, start=(/1,1,m/),kount=(/pcols,pver,1/))
       
       if (masterproc) write(iulog,*)  fieldname, ' is set to short-lived'
  
    enddo

    deallocate(tmpptr)

  end subroutine initialize_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine set_short_lived_species( q, lchnk, ncol, pbuf )

    use physics_buffer, only : physics_buffer_desc, pbuf_set_field

    implicit none 

    real(r8), intent(in)               :: q(pcols,pver,gas_pcnst)
    integer,  intent(in)               :: lchnk, ncol
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: m,n

    if ( nslvd < 1 ) return

    do m=1,nslvd
       n = map(m)
       call pbuf_set_field(pbuf, pbf_idx, q(:,:,n), start=(/1,1,m/),kount=(/pcols,pver,1/))
    enddo

  end subroutine set_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine get_short_lived_species( q, lchnk, ncol, pbuf )
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field

    implicit none 

    real(r8), intent(inout)            :: q(pcols,pver,gas_pcnst)
    integer,  intent(in)               :: lchnk, ncol
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8),pointer                   :: tmpptr(:,:)


    integer :: m,n 

    if ( nslvd < 1 ) return

    do m=1,nslvd
       n = map(m)
       call pbuf_get_field(pbuf, pbf_idx, tmpptr, start=(/1,1,m/), kount=(/ pcols,pver,1 /))
       q(:ncol,:,n) = tmpptr(:ncol,:)
    enddo

  endsubroutine get_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  function slvd_index( name )
    implicit none

    character(len=*) :: name
    integer :: slvd_index

    integer :: m

    slvd_index = -1

    if ( nslvd < 1 ) return

    do m=1,nslvd
       if ( name == slvd_lst(m) ) then
          slvd_index = m
          return 
       endif
    enddo

  endfunction slvd_index

end module short_lived_species
