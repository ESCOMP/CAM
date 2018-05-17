!===============================================================================
! Seasalt for Modal Aerosol Model
!===============================================================================
module seasalt_model
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl => shr_kind_cl
  use ppgrid,         only: pcols, pver
  use modal_aero_data,only: ntot_amode, nslt=>nSeaSalt

  implicit none
  private

  public :: seasalt_nbin
  public :: seasalt_nnum
  public :: seasalt_names
  public :: seasalt_indices
  public :: seasalt_init
  public :: seasalt_emis
  public :: seasalt_active

  integer, protected :: seasalt_nbin ! = nslt
  integer, protected :: seasalt_nnum ! = nnum

  character(len=6), protected, allocatable :: seasalt_names(:)
  integer, protected, allocatable :: seasalt_indices(:)

  logical :: seasalt_active = .false.

  real(r8):: emis_scale

contains
  
  !=============================================================================
  !=============================================================================
  subroutine seasalt_init(seasalt_emis_scale)
    use sslt_sections, only: sslt_sections_init
    use constituents,  only: cnst_get_ind
    use rad_constituents, only: rad_cnst_get_info

    real(r8), intent(in) :: seasalt_emis_scale
    integer :: m, l, nspec, ndx
    character(len=32) :: spec_name
    
    seasalt_nbin = nslt
    seasalt_nnum = nslt
    allocate(seasalt_names(2*nslt))
    allocate(seasalt_indices(2*nslt))

    ndx=0
    do m = 1, ntot_amode
       call rad_cnst_get_info(0, m, nspec=nspec)
       do l = 1, nspec
          call rad_cnst_get_info(0, m, l, spec_name=spec_name )
          if (spec_name(:3) == 'ncl') then
             ndx=ndx+1
             seasalt_names(ndx) = spec_name
             seasalt_names(nslt+ndx) = 'num_'//spec_name(5:)
             call cnst_get_ind(seasalt_names(     ndx), seasalt_indices(     ndx))
             call cnst_get_ind(seasalt_names(nslt+ndx), seasalt_indices(nslt+ndx))
          endif
       enddo
    enddo

    seasalt_active = any(seasalt_indices(:) > 0)
    if (.not.seasalt_active) return

    call sslt_sections_init()

    emis_scale = seasalt_emis_scale

  end subroutine seasalt_init

  !=============================================================================
  !=============================================================================
  subroutine seasalt_emis( u10cubed,  srf_temp, ocnfrc, ncol, cflx )

    use sslt_sections, only: nsections, fluxes, Dg, rdry
    use mo_constants,  only: dns_aer_sst=>seasalt_density, pi

    ! dummy arguments
    real(r8), intent(in) :: u10cubed(:)
    real(r8), intent(in) :: srf_temp(:)
    real(r8), intent(in) :: ocnfrc(:)
    integer,  intent(in) :: ncol
    real(r8), intent(inout) :: cflx(:,:)

    ! local vars
    integer  :: mn, mm, ibin, isec, i
    real(r8) :: fi(ncol,nsections)

    real(r8) :: sst_sz_range_lo (nslt)
    real(r8) :: sst_sz_range_hi (nslt)

    if (nslt==4) then
       sst_sz_range_lo (:) = (/ 0.08e-6_r8, 0.02e-6_r8, 0.3e-6_r8,  1.0e-6_r8 /) ! accu, aitken, fine, coarse
       sst_sz_range_hi (:) = (/ 0.3e-6_r8,  0.08e-6_r8, 1.0e-6_r8, 10.0e-6_r8 /)
    else if (nslt==3) then
       sst_sz_range_lo (:) =  (/ 0.08e-6_r8,  0.02e-6_r8,  1.0e-6_r8 /)  ! accu, aitken, coarse
       sst_sz_range_hi (:) =  (/ 1.0e-6_r8,   0.08e-6_r8, 10.0e-6_r8 /)
    endif

    fi(:ncol,:nsections) = fluxes( srf_temp, u10cubed, ncol )

    do ibin = 1,nslt
       mm = seasalt_indices(ibin)
       mn = seasalt_indices(nslt+ibin)
       
       if (mn>0) then
          do i=1, nsections
             if (Dg(i).ge.sst_sz_range_lo(ibin) .and. Dg(i).lt.sst_sz_range_hi(ibin)) then
                cflx(:ncol,mn)=cflx(:ncol,mn)+fi(:ncol,i)*ocnfrc(:ncol)*emis_scale  !++ ag: scale sea-salt
             endif
          enddo
       endif

       cflx(:ncol,mm)=0.0_r8
       do i=1, nsections
          if (Dg(i).ge.sst_sz_range_lo(ibin) .and. Dg(i).lt.sst_sz_range_hi(ibin)) then
             cflx(:ncol,mm)=cflx(:ncol,mm)+fi(:ncol,i)*ocnfrc(:ncol)*emis_scale  &   !++ ag: scale sea-salt
                  *4._r8/3._r8*pi*rdry(i)**3*dns_aer_sst  ! should use dry size, convert from number to mass flux (kg/m2/s)
          endif
       enddo

    enddo

  end subroutine seasalt_emis

end module seasalt_model
