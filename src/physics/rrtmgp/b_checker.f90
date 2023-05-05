module b_checker
!---------------------------------------------------------------------------------
!
! Instrumentation for debugging CAM interface to RRTMGP radiation parameterization.
!
!---------------------------------------------------------------------------------    
    use shr_kind_mod,          only: r8=>shr_kind_r8, cl=>shr_kind_cl
    use mo_gas_concentrations, only: ty_gas_concs
    use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
    use cam_abortutils,        only: endrun
    use cam_logfile,         only: iulog


    implicit none
    private
    public check_bounds_5d, check_bounds_4d, check_bounds_3d, check_bounds_2d, check_bounds_1d, check_bounds, &
           assert_shape_2dreal, assert_shape

    ! bpm -- interface for checking array bounds
    interface check_bounds
        module procedure check_bounds_1d, check_bounds_2d, check_bounds_3d, check_bounds_4d, check_bounds_5d, check_bounds_gas_concs, check_bounds_gas_optics
    end interface check_bounds

    interface assert_shape
        module procedure assert_shape_2dreal
    end interface assert_shape

    contains

    subroutine check_bounds_1d(arr, max_bound, min_bound, err_message)
        real(r8), intent(in) :: arr(:)
        real(r8), intent(in) :: max_bound, min_bound
        character(len=128), intent(out) :: err_message
        real(r8) :: mx, mn
        err_message=''
        mx = maxval(arr)
        mn = minval(arr)
        if (mn < min_bound) then
            err_message = "validate: array values too small "
        end if
        if (mx > max_bound ) then
            err_message = "validate: array values too large"
        end if
    end subroutine
 
  subroutine check_bounds_2d(arr, max_bound, min_bound, err_message)
    real(r8), intent(in) :: arr(:,:)
    real(r8), intent(in) :: max_bound, min_bound
    character(len=128), intent(out) :: err_message
    real(r8) :: mx, mn
    err_message = ''
    mx = maxval(arr)
    mn = minval(arr)
    if (mn < min_bound) then
        err_message = "validate: array values too small "
    end if
    if (mx > max_bound ) then
        err_message = "validate: array values too large"
    end if
  end subroutine
 
  subroutine check_bounds_3d(arr, max_bound, min_bound, err_message)
    real(r8), intent(in) :: arr(:,:,:)
    real(r8), intent(in) :: max_bound, min_bound
    character(len=128), intent(out) :: err_message
    real(r8) :: mx, mn
    err_message =  ''
    mx = maxval(arr)
    mn = minval(arr)
    if (mn < min_bound) then
        err_message = "validate: array values too small "
    end if
    if (mx > max_bound ) then
        err_message = "validate: array values too large"
    end if
 end subroutine
 
  subroutine check_bounds_4d(arr, max_bound, min_bound, err_message)
    real(r8), intent(in) :: arr(:,:,:,:)
    real(r8), intent(in) :: max_bound, min_bound
    character(len=128), intent(out) :: err_message
    real(r8) :: mx, mn
    err_message = ''
    mx = maxval(arr)
    mn = minval(arr)
    if (mn < min_bound) then
        err_message = "validate: array values too small "
    end if
    if (mx > max_bound ) then
        err_message = "validate: array values too large"
    end if
  end subroutine
 
  subroutine check_bounds_5d(arr, max_bound, min_bound, err_message)
    real(r8), intent(in) :: arr(:,:,:,:,:)
    real(r8), intent(in) :: max_bound, min_bound
    character(len=128), intent(out) :: err_message
    real(r8) :: mx, mn
    err_message = ''
    mx = maxval(arr)
    mn = minval(arr)
    if (mn < min_bound) then
        err_message = "validate: array values too small "
    end if
    if (mx > max_bound ) then
        err_message = "validate: array values too large"
    end if
  end subroutine

  subroutine check_bounds_gas_concs(ncol, nlay, gasconcs, err_message)
    integer,            intent(in)  :: ncol, nlay
    type(ty_gas_concs), intent(in)  :: gasconcs
    character(len=128), intent(out) :: err_message
    character(32), dimension(gasconcs%get_num_gases()) :: gc_gas_names
    integer :: i
    real(r8) :: vmr(ncol,nlay)
    gc_gas_names(:) = gasconcs%get_gas_names()
    do i = 1, gasconcs%get_num_gases()
        err_message = gasconcs%get_vmr(gc_gas_names(i), vmr)  ! gets values in vmr
        if (len_trim(err_message) > 0) then
            call endrun('check_bounds_gas_concs: error getting VMR for '//gc_gas_names(i)//' --> Error Message: '//trim(err_message))
        end if
        call check_bounds(vmr, 1.0_r8, 0.0_r8, err_message)
        if (len_trim(err_message) > 0) then
            err_message = 'check_bounds_gas_concs: VMR error for '//gc_gas_names(i)//' --> Error Message: '//trim(err_message)
        end if      
    end do
  end subroutine

  subroutine check_bounds_gas_optics(kdist, err_message)
    type(ty_gas_optics_rrtmgp), intent(in) :: kdist 
    character(len=128), intent(out) :: err_message
    write(iulog,*) '[check_bonds_gas_optics DRAFT] : kdist'
    ! write(iulog,*) 'number of gases: ',kdist%get_ngas()
    ! write(iulog,*) 'gas names: ',kdist%get_gases()
    ! write(iulog,*) 'kdist%source_is_external() = ',kdist%source_is_external()
    err_message = ""
  end subroutine


  subroutine assert_shape_2dreal(arr, shp, err_message)
    real(r8), intent(in)                :: arr(:,:) ! 2-D array to check
    integer, intent(in)             :: shp(2)  ! Expected shape
    character(len=*), intent(out) :: err_message
    character(len=512) :: err_append
    integer :: r ! rank of arr
    integer :: i
    r = RANK(arr)
    err_message = ''
    if (r .ne. SIZE(shp)) then
        err_message = 'Array is wrong rank (how could that happen?).' 
    end if
    if (len_trim(err_message) == 0) then
        do i = 1,r
            if (SIZE(arr, i) /= shp(i)) then
                write(err_append, "(a39,i3,a2)") 'Array size does not match on Dimension ', i, '._'
                err_message = trim(err_message) // trim(err_append)
            end if
        end do
    end if
end subroutine

end module b_checker
