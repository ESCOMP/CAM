  ! Read the average H2SO4 profile from the initial condition file.
  !
  ! NOTE: This needs to be in its own file to avoid circular references.
  subroutine carma_getH2SO4(h2so4)
    use shr_kind_mod,     only: r8 => shr_kind_r8
    use cam_initfiles,    only: initial_file_get_id
    use pio,              only: file_desc_t
    use ppgrid,           only: pcols, pver, begchunk, endchunk
    use cam_abortutils,   only: endrun
    use cam_grid_support, only: cam_grid_check, cam_grid_id, cam_grid_get_dim_names
    use ncdio_atm,        only: infld
    use gmean_mod,        only: gmean

    real(r8), intent(out)   :: h2so4(pver)     ! midpoint h2so4 mmr (kg/kg)

    integer                    :: iz           ! vertical index
    type(file_desc_t), pointer :: ncid_ini
    logical                    :: found
    real(r8), pointer          :: init_h2so4(:,:,:)
    integer  :: grid_id
    character(len=4) :: dim1name, dim2name
    character(len=*), parameter :: subname = 'carma_getH2SO4'

    ! For an initial run, if the file is missing, then create one using the
    ! average concentration from the initial condition file.
    ncid_ini => initial_file_get_id()

    allocate(init_h2so4(pcols,pver,begchunk:endchunk))

    grid_id = cam_grid_id('physgrid')
    if (.not. cam_grid_check(grid_id)) then
      call endrun(subname//': Internal error, no "physgrid" grid')
    end if
    call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
    call infld('H2SO4', ncid_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, init_h2so4, found, gridname='physgrid')

    if (.not. found) then
      call endrun(subname//': failed to find field H2SO4 in IC file.')
    end if

    ! Just do a simple average. Could get gw and do a weighted average.
    do iz = 1, pver
      call gmean(init_h2so4(:, iz, :), h2so4(iz))
    end do

    deallocate(init_h2so4)

  end subroutine carma_getH2SO4
