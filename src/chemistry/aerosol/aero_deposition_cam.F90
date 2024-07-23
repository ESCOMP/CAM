module aero_deposition_cam

  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: cnst_get_ind, pcnst
  use camsrfexch,   only: cam_out_t
  use aerosol_properties_mod, only: aero_name_len
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private
  public :: aero_deposition_cam_init
  public :: aero_deposition_cam_setwet

  integer :: bcarbon_ndx( pcnst ) = -1
  integer :: bcarbon_cnt = 0
  integer :: ocarbon_ndx( pcnst ) = -1
  integer :: ocarbon_cnt = 0

  integer :: crse_dust_ndx( pcnst ) = -1
  integer :: fine_dust_ndx( pcnst ) = -1
  integer :: crse_dust_cnt = 0
  integer :: fine_dust_cnt = 0

contains

  !==============================================================================
  subroutine aero_deposition_cam_init(aero_props)

    class(aerosol_properties), intent(in) :: aero_props

    integer :: pcnt, scnt

    call get_indices( type='black-c',   indices=bcarbon_ndx, count=bcarbon_cnt )
    call get_indices( type='s-organic', indices=ocarbon_ndx, count=pcnt )
    call get_indices( type='p-organic', indices=ocarbon_ndx(pcnt+1:), count=scnt )
    ocarbon_cnt = pcnt+scnt

    ! fine dust has radius less than 1.25 microns
    call get_indices( type='dust', indices=fine_dust_ndx, count=fine_dust_cnt, max_radius=1.25e-6_r8 )
    call get_indices( type='dust', indices=crse_dust_ndx, count=crse_dust_cnt, min_radius=1.25e-6_r8 ) ! meters

  contains

    !==============================================================================
    subroutine get_indices( type, indices, count, min_radius, max_radius )

      character(len=*), intent(in) :: type
      integer, intent(out) :: indices(:)
      integer, intent(out) :: count
      real(r8), intent(in), optional :: min_radius ! meters
      real(r8), intent(in), optional :: max_radius ! meters

      integer :: ibin,ispc, ndx, nspec
      character(len=aero_name_len) :: spec_type, spec_name
      real(r8) :: minrad
      logical :: getndx

      count = 0
      indices(:) = -1


      do ibin = 1, aero_props%nbins()

         do ispc = 1, aero_props%nspecies(ibin)

            call aero_props%get(ibin,ispc, spectype=spec_type, specname=spec_name )

            if (spec_type==type) then

               getndx = .true.
               minrad = aero_props%min_mass_mean_rad(ibin,ispc) ! meters

               if (present(min_radius)) then
                  getndx = min_radius < minrad ! coarse
               elseif (present(max_radius)) then
                  getndx = max_radius > minrad .and. minrad > 0._r8 ! fine
               end if
               if (getndx) then
                  call cnst_get_ind(spec_name, ndx, abort=.false.)
                  if (ndx>0) then
                     count = count+1
                     indices(count) = ndx
                  endif
               endif

            endif

         enddo


      enddo

    end subroutine get_indices

  end subroutine aero_deposition_cam_init

  !==============================================================================
  subroutine aero_deposition_cam_setwet(aerdepwetis, aerdepwetcw, cam_out)


   ! Arguments:
   real(r8), intent(in) :: aerdepwetis(:,:)  ! aerosol wet deposition (interstitial)
   real(r8), intent(in) :: aerdepwetcw(:,:)  ! aerosol wet deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i, ispec
   integer :: ncol                      ! number of columns

   ncol = cam_out%ncol

   cam_out%bcphiwet(:) = 0._r8
   cam_out%ocphiwet(:) = 0._r8
   cam_out%dstwet1(:) = 0._r8
   cam_out%dstwet2(:) = 0._r8
   cam_out%dstwet3(:) = 0._r8
   cam_out%dstwet4(:) = 0._r8

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface,
   !        dry deposition fluxes are positive into surface.
   !        srf models want positive definite fluxes.
   do i = 1, ncol

      ! black carbon fluxes
      do ispec=1,bcarbon_cnt
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) &
                             - (aerdepwetis(i,bcarbon_ndx(ispec))+aerdepwetcw(i,bcarbon_ndx(ispec)))
      enddo

      ! organic carbon fluxes
      do ispec=1,ocarbon_cnt
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) &
                             - (aerdepwetis(i,ocarbon_ndx(ispec))+aerdepwetcw(i,ocarbon_ndx(ispec)))
      enddo

      !  Assign "fine" dust to bulk size bin 1:
      do ispec=1,fine_dust_cnt
         cam_out%dstwet1(i) = cam_out%dstwet1(i) &
                            -(aerdepwetis(i,fine_dust_ndx(ispec))+aerdepwetcw(i,fine_dust_ndx(ispec)))
      enddo

      !  Assign "coarse" dust to bulk size bin 3:
      do ispec=1,crse_dust_cnt
         cam_out%dstwet3(i) = cam_out%dstwet3(i) &
                            -(aerdepwetis(i,crse_dust_ndx(ispec))+aerdepwetcw(i,crse_dust_ndx(ispec)))
      enddo

      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphiwet(i) .lt. 0._r8) cam_out%bcphiwet(i) = 0._r8
      if (cam_out%ocphiwet(i) .lt. 0._r8) cam_out%ocphiwet(i) = 0._r8
      if (cam_out%dstwet1(i)  .lt. 0._r8) cam_out%dstwet1(i)  = 0._r8
      if (cam_out%dstwet3(i)  .lt. 0._r8) cam_out%dstwet3(i)  = 0._r8

   enddo

  end subroutine aero_deposition_cam_setwet

end module aero_deposition_cam
