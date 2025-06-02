!------------------------------------------------------------------------------
! Provides methods for mapping from physics grid to regular longitude / latitude
! grid via ESMF regridding capabilities
!------------------------------------------------------------------------------
module esmf_phys2lonlat_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_logfile,  only: iulog
  use cam_abortutils, only: endrun
  use spmd_utils, only: masterproc
  use ppgrid, only: pver

  use ESMF, only: ESMF_RouteHandle, ESMF_Field, ESMF_ArraySpec, ESMF_ArraySpecSet
  use ESMF, only: ESMF_FieldCreate, ESMF_FieldRegridStore
  use ESMF, only: ESMF_FieldGet, ESMF_FieldRegrid
  use ESMF, only: ESMF_KIND_I4, ESMF_KIND_R8, ESMF_TYPEKIND_R8
  use ESMF, only: ESMF_REGRIDMETHOD_BILINEAR, ESMF_POLEMETHOD_ALLAVG, ESMF_EXTRAPMETHOD_NEAREST_IDAVG
  use ESMF, only: ESMF_TERMORDER_SRCSEQ, ESMF_MESHLOC_ELEMENT, ESMF_STAGGERLOC_CENTER
  use ESMF, only: ESMF_FieldDestroy, ESMF_RouteHandleDestroy
  use esmf_check_error_mod, only: check_esmf_error

  implicit none

  private

  public :: esmf_phys2lonlat_init
  public :: esmf_phys2lonlat_regrid
  public :: esmf_phys2lonlat_destroy
  public :: fields_bundle_t
  public :: nflds

  type(ESMF_RouteHandle) :: rh_phys2lonlat_3d
  type(ESMF_RouteHandle) :: rh_phys2lonlat_2d

  type(ESMF_Field) :: physfld_3d
  type(ESMF_Field) :: lonlatfld_3d

  type(ESMF_Field) :: physfld_2d
  type(ESMF_Field) :: lonlatfld_2d

  interface esmf_phys2lonlat_regrid
     module procedure esmf_phys2lonlat_regrid_2d
     module procedure esmf_phys2lonlat_regrid_3d
  end interface esmf_phys2lonlat_regrid

  type :: fields_bundle_t
     real(r8), pointer :: fld(:,:,:) => null()
  end type fields_bundle_t

  integer, parameter :: nflds = 5

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine esmf_phys2lonlat_init()
    use esmf_phys_mesh_mod, only: physics_grid_mesh
    use esmf_lonlat_grid_mod, only: lonlat_grid

    type(ESMF_ArraySpec) :: arrayspec
    integer(ESMF_KIND_I4), pointer :: factorIndexList(:,:)
    real(ESMF_KIND_R8),    pointer :: factorList(:)
    integer                        :: smm_srctermproc,  smm_pipelinedep, rc

    character(len=*), parameter :: subname  = 'esmf_phys2lonlat_init: '

    smm_srctermproc = 0
    smm_pipelinedep = 16

    ! create ESMF fields

    ! 3D phys fld
    call ESMF_ArraySpecSet(arrayspec, 3, ESMF_TYPEKIND_R8, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_ArraySpecSet 3D phys fld ERROR')

    physfld_3d = ESMF_FieldCreate(physics_grid_mesh, arrayspec, &
                                  gridToFieldMap=(/3/), meshloc=ESMF_MESHLOC_ELEMENT, &
                                  ungriddedLBound=(/1,1/), ungriddedUBound=(/pver,nflds/), rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldCreate 3D phys fld ERROR')

    ! 3D lon lat grid
    call ESMF_ArraySpecSet(arrayspec, 4, ESMF_TYPEKIND_R8, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_ArraySpecSet 3D lonlat fld ERROR')

    lonlatfld_3d = ESMF_FieldCreate( lonlat_grid, arrayspec, staggerloc=ESMF_STAGGERLOC_CENTER, &
                                     ungriddedLBound=(/1,1/), ungriddedUBound=(/pver,nflds/), rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldCreate 3D lonlat fld ERROR')

    ! 2D phys fld
    call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_ArraySpecSet 2D phys fld ERROR')

    physfld_2d = ESMF_FieldCreate(physics_grid_mesh, arrayspec, &
                                  meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldCreate 2D phys fld ERROR')

    ! 2D lon/lat grid
    call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_ArraySpecSet 2D lonlat fld ERROR')

    lonlatfld_2d = ESMF_FieldCreate( lonlat_grid, arrayspec, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldCreate 2D lonlat fld ERROR')

    call ESMF_FieldRegridStore(srcField=physfld_3d, dstField=lonlatfld_3d, &
         regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                          &
         polemethod=ESMF_POLEMETHOD_ALLAVG,                                &
         extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                     &
         routeHandle=rh_phys2lonlat_3d, factorIndexList=factorIndexList,   &
         factorList=factorList, srcTermProcessing=smm_srctermproc,         &
         pipelineDepth=smm_pipelinedep, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldRegridStore 3D routehandle ERROR')

    call ESMF_FieldRegridStore(srcField=physfld_2d, dstField=lonlatfld_2d, &
         regridMethod=ESMF_REGRIDMETHOD_BILINEAR,                          &
         polemethod=ESMF_POLEMETHOD_ALLAVG,                                &
         extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG,                     &
         routeHandle=rh_phys2lonlat_2d, factorIndexList=factorIndexList,   &
         factorList=factorList, srcTermProcessing=smm_srctermproc,         &
         pipelineDepth=smm_pipelinedep, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldRegridStore 3D routehandle ERROR')

  end subroutine esmf_phys2lonlat_init

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine esmf_phys2lonlat_regrid_3d(physflds, lonlatflds)
    use esmf_lonlat_grid_mod, only: lon_beg,lon_end,lat_beg,lat_end
    use ppgrid, only: pcols, pver, begchunk, endchunk
    use phys_grid, only: get_ncols_p

    type(fields_bundle_t) :: physflds(nflds)
    type(fields_bundle_t) :: lonlatflds(nflds)

    integer :: i, ichnk, ncol, ifld, ilev, icol, rc
    real(ESMF_KIND_R8), pointer :: physptr(:,:,:)
    real(ESMF_KIND_R8), pointer :: lonlatptr(:,:,:,:)

    character(len=*), parameter :: subname = 'esmf_phys2lonlat_regrid: '

    call ESMF_FieldGet(physfld_3d, localDe=0, farrayPtr=physptr, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldGet physptr')

    i = 0
    do ichnk = begchunk, endchunk
       ncol = get_ncols_p(ichnk)
       do icol = 1,ncol
          i = i+1
          do ifld = 1,nflds
             do ilev = 1,pver
                physptr(ilev,ifld,i) = physflds(ifld)%fld(ilev,icol,ichnk)
             end do
          end do
       end do
    end do

    call ESMF_FieldRegrid(physfld_3d, lonlatfld_3d, rh_phys2lonlat_3d, &
              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldRegrid physfld_3d->lonlatfld_3d')

    call ESMF_FieldGet(lonlatfld_3d, localDe=0, farrayPtr=lonlatptr, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldGet lonlatptr')

    do ifld = 1,nflds
       lonlatflds(ifld)%fld(lon_beg:lon_end,lat_beg:lat_end,1:pver) = lonlatptr(lon_beg:lon_end,lat_beg:lat_end,1:pver,ifld)
    end do

  end subroutine esmf_phys2lonlat_regrid_3d

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine esmf_phys2lonlat_regrid_2d(physarr, lonlatarr)
    use esmf_lonlat_grid_mod, only: lon_beg,lon_end,lat_beg,lat_end
    use ppgrid, only: pcols, pver, begchunk, endchunk
    use phys_grid, only: get_ncols_p

    real(r8),intent(in) :: physarr(pcols,begchunk:endchunk)
    real(r8),intent(out) :: lonlatarr(lon_beg:lon_end,lat_beg:lat_end)

    integer :: i, ichnk, ncol, icol, rc
    real(ESMF_KIND_R8), pointer :: physptr(:)
    real(ESMF_KIND_R8), pointer :: lonlatptr(:,:)

    character(len=*), parameter :: subname = 'esmf_phys2lonlat_regrid_2d: '

    call ESMF_FieldGet(physfld_2d, localDe=0, farrayPtr=physptr, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldGet physptr')

    i = 0
    do ichnk = begchunk, endchunk
       ncol = get_ncols_p(ichnk)
       do icol = 1,ncol
          i = i+1
          physptr(i) = physarr(icol,ichnk)
       end do
    end do

    call ESMF_FieldRegrid(physfld_2d, lonlatfld_2d, rh_phys2lonlat_2d, &
              termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldRegrid physfld_3d->lonlatfld_3d')

    call ESMF_FieldGet(lonlatfld_2d, localDe=0, farrayPtr=lonlatptr, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldGet lonlatptr')

    lonlatarr(lon_beg:lon_end,lat_beg:lat_end) = lonlatptr(lon_beg:lon_end,lat_beg:lat_end)

  end subroutine esmf_phys2lonlat_regrid_2d

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine esmf_phys2lonlat_destroy()

    integer :: rc
    character(len=*), parameter :: subname = 'esmf_phys2lonlat_destroy: '

    call ESMF_RouteHandleDestroy(rh_phys2lonlat_3d, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldDestroy rh_phys2lonlat_3d')

    call ESMF_RouteHandleDestroy(rh_phys2lonlat_2d, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldDestroy rh_phys2lonlat_2d')

    call ESMF_FieldDestroy(lonlatfld_3d, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldDestroy lonlatfld_3d')

    call ESMF_FieldDestroy(lonlatfld_2d, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldDestroy lonlatfld_2d')

    call ESMF_FieldDestroy(physfld_3d, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldDestroy physfld_3d')

    call ESMF_FieldDestroy(physfld_2d, rc=rc)
    call check_esmf_error(rc, subname//'ESMF_FieldDestroy physfld_2d')

  end subroutine esmf_phys2lonlat_destroy

end module esmf_phys2lonlat_mod
