module aero_deposition_cam
!------------------------------------------------------------------------------
! Purpose:
!
! Partition the contributions from aerosols of wet and dry
! deposition at the surface into the fields passed to the coupler.
!------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_kind_mod, only: cl => shr_kind_cl
  use constituents, only: cnst_get_ind, pcnst
  use camsrfexch,   only: cam_out_t
  use cam_abortutils,only: endrun
  use aerosol_properties_mod, only: aero_name_len
  use aerosol_properties_mod, only: aerosol_properties

  implicit none

  private

! Public interfaces

  public :: aero_deposition_cam_init
  public :: aero_deposition_cam_setwet
  public :: aero_deposition_cam_setdry

! Private module data

  integer :: bcphi_ndx( pcnst ) = -1
  integer :: bcphi_cnt = 0
  integer :: bcpho_ndx( pcnst ) = -1
  integer :: bcpho_cnt = 0
  integer :: ocphi_ndx( pcnst ) = -1
  integer :: ocphi_cnt = 0
  integer :: ocpho_ndx( pcnst ) = -1
  integer :: ocpho_cnt = 0

  class(aerosol_properties), pointer :: aero_props=>null()
  integer :: nele_tot=0            ! total number of aerosol elements

  ! bulk dust bins (meters)

  integer, parameter :: n_bulk_dst_bins = 4

  ! CAM4 bulk dust bin sizes (https://doi.org/10.1002/2013MS000279)
  real(r8), parameter :: bulk_dst_edges(n_bulk_dst_bins+1) = &
       (/0.1e-6_r8, 1.0e-6_r8, 2.5e-6_r8, 5.0e-6_r8, 10.e-6_r8/)

contains

  !============================================================================
  subroutine aero_deposition_cam_init(aero_props_in)

    class(aerosol_properties),target, intent(in) :: aero_props_in

    integer :: pcnt, scnt
    character(len=*), parameter :: subrname = 'aero_deposition_cam_init'

    ! construct the aerosol properties object
    aero_props => aero_props_in

    ! set the cam constituent indices and determine the counts
    ! for the specified aerosol types

    ! black carbons
    call get_indices( type='black-c',  hydrophilic=.true.,  indices=bcphi_ndx, count=bcphi_cnt )
    call get_indices( type='black-c',  hydrophilic=.false., indices=bcpho_ndx, count=bcpho_cnt )

    ! primary and secondary organics
    call get_indices( type='p-organic',hydrophilic=.true.,  indices=ocphi_ndx, count=pcnt )
    call get_indices( type='s-organic',hydrophilic=.true.,  indices=ocphi_ndx(pcnt+1:), count=scnt )
    ocphi_cnt = pcnt+scnt

    call get_indices( type='p-organic',hydrophilic=.false., indices=ocpho_ndx, count=pcnt )
    call get_indices( type='s-organic',hydrophilic=.false., indices=ocpho_ndx(pcnt+1:), count=scnt )
    ocpho_cnt = pcnt+scnt

    ! total number of aerosol elements
    nele_tot = aero_props%ncnst_tot()

  contains

    !==========================================================================
    ! returns CAM constituent indices of the aerosol tracers (and count)
    !==========================================================================
    subroutine get_indices( type, hydrophilic, indices, count)

      character(len=*), intent(in) :: type
      logical, intent(in ) :: hydrophilic
      integer, intent(out) :: indices(:)
      integer, intent(out) :: count

      integer :: ibin,ispc, ndx, nspec
      character(len=aero_name_len) :: spec_type, spec_name

      count = 0
      indices(:) = -1

      ! loop through aerosol bins / modes
      do ibin = 1, aero_props%nbins()

         ! check if the bin/mode is hydrophilic
         if ( aero_props%hydrophilic(ibin) .eqv. hydrophilic ) then
            do ispc = 1, aero_props%nspecies(ibin)

               call aero_props%get(ibin,ispc, spectype=spec_type, specname=spec_name)

               if (spec_type==type) then

                  ! get CAM constituent index
                  call cnst_get_ind(spec_name, ndx, abort=.false.)
                  if (ndx>0) then
                     count = count+1
                     indices(count) = ndx
                  endif

               endif

            enddo
         endif

      enddo

    end subroutine get_indices

  end subroutine aero_deposition_cam_init

  !============================================================================
  ! Set surface wet deposition fluxes passed to coupler.
  !============================================================================
  subroutine aero_deposition_cam_setwet(aerdepwetis, aerdepwetcw, cam_out)

   ! Arguments:
   real(r8), intent(in) :: aerdepwetis(:,:)  ! aerosol wet deposition (interstitial)
   real(r8), intent(in) :: aerdepwetcw(:,:)  ! aerosol wet deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i, ispec, ibin, mm, ndx
   integer :: ncol                      ! number of columns

   real(r8) :: dep_fluxes(nele_tot)
   real(r8) :: dst_fluxes(n_bulk_dst_bins)
   character(len=aero_name_len) :: specname, name_c
   integer :: errstat
   character(len=cl) :: errstr

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

      ! hydrophilic black carbon fluxes
      do ispec=1,bcphi_cnt
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) &
                             - (aerdepwetis(i,bcphi_ndx(ispec))+aerdepwetcw(i,bcphi_ndx(ispec)))
      enddo

      ! hydrophobic black carbon fluxes
      do ispec=1,bcpho_cnt
         cam_out%bcphiwet(i) = cam_out%bcphiwet(i) &
                             - (aerdepwetis(i,bcpho_ndx(ispec))+aerdepwetcw(i,bcpho_ndx(ispec)))
      enddo

      ! hydrophilic organic carbon fluxes
      do ispec=1,ocphi_cnt
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) &
                             - (aerdepwetis(i,ocphi_ndx(ispec))+aerdepwetcw(i,ocphi_ndx(ispec)))
      enddo

      ! hydrophobic organic carbon fluxes
      do ispec=1,ocpho_cnt
         cam_out%ocphiwet(i) = cam_out%ocphiwet(i) &
                             - (aerdepwetis(i,ocpho_ndx(ispec))+aerdepwetcw(i,ocpho_ndx(ispec)))
      enddo

      ! dust fluxes

      dep_fluxes = 0._r8
      dst_fluxes = 0._r8

      do ibin = 1,aero_props%nbins()
         do ispec = 0,aero_props%nmasses(ibin)
            if (ispec==0) then
               call aero_props%num_names(ibin, specname, name_c)
            else
               call aero_props%get(ibin,ispec, specname=specname)
            end if
            call cnst_get_ind(specname, ndx, abort=.false.)
            if (ndx>0) then
               mm = aero_props%indexer(ibin,ispec)
               dep_fluxes(mm) = - (aerdepwetis(i,ndx)+aerdepwetcw(i,ndx))
            end if
         end do
      end do

      ! rebin dust fluxes to bulk dust bins
      call aero_props%rebin_bulk_fluxes('dust', dep_fluxes, bulk_dst_edges, dst_fluxes, errstat, errstr)
      if (errstat/=0) then
         call endrun('aero_deposition_cam_setwet: '//trim(errstr))
      end if

      cam_out%dstwet1(i) = cam_out%dstwet1(i) + dst_fluxes(1)
      cam_out%dstwet2(i) = cam_out%dstwet2(i) + dst_fluxes(2)
      cam_out%dstwet3(i) = cam_out%dstwet3(i) + dst_fluxes(3)
      cam_out%dstwet4(i) = cam_out%dstwet4(i) + dst_fluxes(4)

      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphiwet(i) < 0._r8) cam_out%bcphiwet(i) = 0._r8
      if (cam_out%ocphiwet(i) < 0._r8) cam_out%ocphiwet(i) = 0._r8
      if (cam_out%dstwet1(i)  < 0._r8) cam_out%dstwet1(i)  = 0._r8
      if (cam_out%dstwet2(i)  < 0._r8) cam_out%dstwet2(i)  = 0._r8
      if (cam_out%dstwet3(i)  < 0._r8) cam_out%dstwet3(i)  = 0._r8
      if (cam_out%dstwet4(i)  < 0._r8) cam_out%dstwet4(i)  = 0._r8

   enddo

  end subroutine aero_deposition_cam_setwet

  !============================================================================
  ! Set surface dry deposition fluxes passed to coupler.
  !============================================================================
  subroutine aero_deposition_cam_setdry(aerdepdryis, aerdepdrycw, cam_out)

   ! Arguments:
   real(r8), intent(in) :: aerdepdryis(:,:)  ! aerosol dry deposition (interstitial)
   real(r8), intent(in) :: aerdepdrycw(:,:)  ! aerosol dry deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i, ispec, ibin, mm, ndx
   integer :: ncol                      ! number of columns

   real(r8) :: dep_fluxes(nele_tot)
   real(r8) :: dst_fluxes(n_bulk_dst_bins)
   character(len=aero_name_len) :: specname, name_c
   integer :: errstat
   character(len=cl) :: errstr

   ncol = cam_out%ncol

   cam_out%bcphidry(:) = 0._r8
   cam_out%ocphidry(:) = 0._r8
   cam_out%bcphodry(:) = 0._r8
   cam_out%ocphodry(:) = 0._r8
   cam_out%dstdry1(:) = 0._r8
   cam_out%dstdry2(:) = 0._r8
   cam_out%dstdry3(:) = 0._r8
   cam_out%dstdry4(:) = 0._r8

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface,
   !        dry deposition fluxes are positive into surface.
   !        srf models want positive definite fluxes.
   do i = 1, ncol

      ! hydrophilic black carbon fluxes
      do ispec=1,bcphi_cnt
         cam_out%bcphidry(i) = cam_out%bcphidry(i) &
                             + (aerdepdryis(i,bcphi_ndx(ispec))+aerdepdrycw(i,bcphi_ndx(ispec)))
      enddo

      ! hydrophobic black carbon fluxes
      do ispec=1,bcpho_cnt
         cam_out%bcphodry(i) = cam_out%bcphodry(i) &
                             + (aerdepdryis(i,bcpho_ndx(ispec))+aerdepdrycw(i,bcpho_ndx(ispec)))
      enddo

      ! hydrophilic organic carbon fluxes
      do ispec=1,ocphi_cnt
         cam_out%ocphidry(i) = cam_out%ocphidry(i) &
                             + (aerdepdryis(i,ocphi_ndx(ispec))+aerdepdrycw(i,ocphi_ndx(ispec)))
      enddo

      ! hydrophobic organic carbon fluxes
      do ispec=1,ocpho_cnt
         cam_out%ocphodry(i) = cam_out%ocphodry(i) &
                             + (aerdepdryis(i,ocpho_ndx(ispec))+aerdepdrycw(i,ocpho_ndx(ispec)))
      enddo

      ! dust fluxes

      dep_fluxes = 0._r8
      dst_fluxes = 0._r8

      do ibin = 1,aero_props%nbins()
         do ispec = 0,aero_props%nspecies(ibin)
            if (ispec==0) then
               call aero_props%num_names(ibin, specname, name_c)
            else
               call aero_props%get(ibin,ispec, specname=specname)
            end if
            call cnst_get_ind(specname, ndx, abort=.false.)
            if (ndx>0) then
               mm = aero_props%indexer(ibin,ispec)
               dep_fluxes(mm) = aerdepdryis(i,ndx)+aerdepdrycw(i,ndx)
            end if
         end do
      end do

      ! rebin dust fluxes to bulk dust bins
      call aero_props%rebin_bulk_fluxes('dust', dep_fluxes, bulk_dst_edges, dst_fluxes, errstat, errstr)
      if (errstat/=0) then
         call endrun('aero_deposition_cam_setdry: '//trim(errstr))
      end if

      cam_out%dstdry1(i) = cam_out%dstdry1(i) + dst_fluxes(1)
      cam_out%dstdry2(i) = cam_out%dstdry2(i) + dst_fluxes(2)
      cam_out%dstdry3(i) = cam_out%dstdry3(i) + dst_fluxes(3)
      cam_out%dstdry4(i) = cam_out%dstdry4(i) + dst_fluxes(4)

      ! in rare cases, integrated deposition tendency is upward
      if (cam_out%bcphidry(i) < 0._r8) cam_out%bcphidry(i) = 0._r8
      if (cam_out%ocphidry(i) < 0._r8) cam_out%ocphidry(i) = 0._r8
      if (cam_out%bcphodry(i) < 0._r8) cam_out%bcphodry(i) = 0._r8
      if (cam_out%ocphodry(i) < 0._r8) cam_out%ocphodry(i) = 0._r8
      if (cam_out%dstdry1(i)  < 0._r8) cam_out%dstdry1(i)  = 0._r8
      if (cam_out%dstdry2(i)  < 0._r8) cam_out%dstdry2(i)  = 0._r8
      if (cam_out%dstdry3(i)  < 0._r8) cam_out%dstdry3(i)  = 0._r8
      if (cam_out%dstdry4(i)  < 0._r8) cam_out%dstdry4(i)  = 0._r8

   enddo

  end subroutine aero_deposition_cam_setdry

end module aero_deposition_cam
