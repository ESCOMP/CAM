  ! Read the average temperature profile from the initial condition file.
  !
  ! NOTE: This needs to be in its own file to avoid circular references.
  subroutine carma_getT(T)
    use shr_kind_mod,     only: r8 => shr_kind_r8
    use cam_initfiles,    only: initial_file_get_id
    use pio,              only: file_desc_t
    use ppgrid,           only: pcols, pver, begchunk, endchunk
    use cam_abortutils,   only: endrun
    use cam_grid_support, only: cam_grid_check, cam_grid_id, cam_grid_get_dim_names
    use ncdio_atm,        only: infld
    use gmean_mod,        only: gmean

    real(r8), intent(out)   :: T(pver)      ! midpoint temperature (Pa)

    integer                    :: iz           ! vertical index
    type(file_desc_t), pointer :: ncid_ini
    logical                    :: found
    real(r8), pointer          :: init_t(:,:,:)
    integer  :: grid_id
    character(len=4) :: dim1name, dim2name
    character(len=*), parameter :: subname = 'carma_getT'

    ! For an initial run, if the file is missing, then create one using the average
    ! temperature from the initial condition file.
    ncid_ini => initial_file_get_id()

    allocate(init_t(pcols,pver,begchunk:endchunk))

    grid_id = cam_grid_id('physgrid')
    if (.not. cam_grid_check(grid_id)) then
      call endrun(subname//': Internal error, no "physgrid" grid')
    end if
    call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
    call infld('T', ncid_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, init_t, found, gridname='physgrid')

    if (.not. found) then
      call endrun(subname//': failed to find field T in IC file.')
    end if

    ! Just do a simple average. Could get gw and do a weighted average.
    do iz = 1, pver
      call gmean(init_t(:, iz, :), T(iz))
    end do

    deallocate(init_t)

  end subroutine carma_getT
