!-------------------------------------------------------------------------------
! Utility module for mapping fields between CAM physics, oplus transport, and
! geomagnetic grids
!-------------------------------------------------------------------------------
module regridder
  use shr_kind_mod,only: r8 => shr_kind_r8 ! 8-byte reals
  use cam_abortutils, only: endrun

  use edyn_mpi, only: mlon0, mlon1, mlat0, mlat1, mlev0, mlev1
  use edyn_mpi, only: lon0, lon1, lat0, lat1, lev0, lev1

  use edyn_esmf, only: edyn_esmf_set3d_phys, edyn_esmf_regrid_phys2mag
  use edyn_esmf, only: edyn_esmf_regrid_phys2geo, edyn_esmf_get_3dfield
  use edyn_esmf, only: edyn_esmf_set2d_phys, edyn_esmf_get_2dfield, edyn_esmf_get_2dphysfield, edyn_esmf_set3d_geo
  use edyn_esmf, only: edyn_esmf_regrid_geo2mag, edyn_esmf_regrid_geo2phys
  use edyn_esmf, only: edyn_esmf_set2d_geo, edyn_esmf_set3d_mag, edyn_esmf_regrid_mag2geo
  use edyn_esmf, only: phys_3dfld, phys_2dfld
  use edyn_esmf, only: geo_3dfld, geo_2dfld
  use edyn_esmf, only: mag_des_3dfld, mag_des_2dfld
  use edyn_esmf, only: mag_src_3dfld, mag_src_2dfld
  use edyn_esmf, only: edyn_esmf_set2d_mag, edyn_esmf_regrid_mag2phys, edyn_esmf_get_1dfield

  implicit none

contains

  !-----------------------------------------------------------------------------
  ! map horizontal 2D fields from magnetic grid to physcis grid
  !-----------------------------------------------------------------------------
  subroutine regrid_mag2phys_2d(magfld, physfld, cols, cole)
    integer,  intent(in)  :: cols, cole
    real(r8), intent(in)  :: magfld(mlon0:mlon1,mlat0:mlat1)
    real(r8), intent(out) :: physfld(cols:cole)

    call edyn_esmf_set2d_mag( mag_src_2dfld, magfld, mlon0, mlon1, mlat0, mlat1 )
    call edyn_esmf_regrid_mag2phys( mag_src_2dfld, phys_2dfld, 2)
    call edyn_esmf_get_1dfield(phys_2dfld, physfld, cols, cole  )

  end subroutine regrid_mag2phys_2d

  !-----------------------------------------------------------------------------
  ! map 3D feilds from magnetic grid to oplus grid
  !-----------------------------------------------------------------------------
  subroutine regrid_mag2geo_3d(magfld,geofld)
    real(r8), intent(in)  :: magfld(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1)
    real(r8), intent(out) :: geofld(lon0:lon1,lat0:lat1,lev0:lev1)

    call edyn_esmf_set3d_mag( mag_src_3dfld, magfld, mlon0, mlon1, mlat0, mlat1, mlev0, mlev1 )
    call edyn_esmf_regrid_mag2geo(mag_src_3dfld, geo_3dfld, 3)
    call edyn_esmf_get_3dfield(geo_3dfld, geofld, lon0, lon1, lat0, lat1, lev0, lev1)

  end subroutine regrid_mag2geo_3d

  !-----------------------------------------------------------------------------
  ! map horizontal 2D fields from physcis grid to oplus grid
  !-----------------------------------------------------------------------------
  subroutine regrid_phys2geo_2d( physfld, geofld, cols, cole )
    integer,  intent(in)  :: cols, cole
    real(r8), intent(in)  :: physfld(cols:cole)
    real(r8), intent(out) :: geofld(lon0:lon1,lat0:lat1)

    call edyn_esmf_set2d_phys( phys_2dfld , physfld, cols, cole)
    call edyn_esmf_regrid_phys2geo(phys_2dfld, geo_2dfld, 2)
    call edyn_esmf_get_2dfield(geo_2dfld, geofld, lon0, lon1, lat0, lat1 )

  end subroutine regrid_phys2geo_2d

  !-----------------------------------------------------------------------------
  ! map 3D fields from physcis grid to oplus grid
  !-----------------------------------------------------------------------------
  subroutine regrid_phys2geo_3d( physfld, geofld, plev, cols, cole )
    integer,  intent(in)  :: plev, cols, cole
    real(r8), intent(in)  :: physfld(1:plev,cols:cole)
    real(r8), intent(out) :: geofld(lon0:lon1,lat0:lat1,lev0:lev1)

    call edyn_esmf_set3d_phys( phys_3dfld, physfld,  1, plev, cols, cole)
    call edyn_esmf_regrid_phys2geo(phys_3dfld, geo_3dfld, 3)
    call edyn_esmf_get_3dfield(geo_3dfld, geofld, lon0, lon1, lat0, lat1, lev0, lev1 )

  end subroutine regrid_phys2geo_3d

  !-----------------------------------------------------------------------------
  ! map 3D fields from oplus grid to magnetic grid
  !-----------------------------------------------------------------------------
  subroutine regrid_geo2mag_3d( geofld, magfld )
    real(r8), intent(in) :: geofld(lon0:lon1,lat0:lat1,lev0:lev1)
    real(r8), intent(out) :: magfld(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1)

    call edyn_esmf_set3d_geo( geo_3dfld, geofld, lon0, lon1, lat0, lat1, lev0, lev1 )
    call edyn_esmf_regrid_geo2mag(geo_3dfld, mag_des_3dfld, 3)
    call edyn_esmf_get_3dfield(mag_des_3dfld, magfld, mlon0, mlon1, mlat0, mlat1, mlev0, mlev1 )

  end subroutine regrid_geo2mag_3d

  !-----------------------------------------------------------------------------
  ! map horizontal 2D fields from oplus grid to magnetic grid
  !-----------------------------------------------------------------------------
  subroutine regrid_geo2mag_2d( geofld, magfld )
    real(r8), intent(in) :: geofld(lon0:lon1,lat0:lat1)
    real(r8), intent(out) :: magfld(mlon0:mlon1,mlat0:mlat1)

    call edyn_esmf_set2d_geo( geo_2dfld, geofld, lon0, lon1, lat0, lat1 )
    call edyn_esmf_regrid_geo2mag(geo_2dfld, mag_des_2dfld, 2)
    call edyn_esmf_get_2dfield(mag_des_2dfld, magfld, mlon0, mlon1, mlat0, mlat1 )

  end subroutine regrid_geo2mag_2d

  !-----------------------------------------------------------------------------
  ! map 3D fields from oplus grid to physics grid
  !-----------------------------------------------------------------------------
  subroutine regrid_geo2phys_3d( geofld, physfld, plev, cols, cole )
    integer,  intent(in)  :: plev, cols, cole
    real(r8), intent(in)  :: geofld(lon0:lon1,lat0:lat1,lev0:lev1)
    real(r8), intent(out) :: physfld(1:plev,cols:cole)


    call edyn_esmf_set3d_geo( geo_3dfld, geofld, lon0, lon1, lat0, lat1, lev0, lev1 )
    call edyn_esmf_regrid_geo2phys(geo_3dfld, phys_3dfld, 3)
    call edyn_esmf_get_2dphysfield(phys_3dfld, physfld, 1, plev, cols, cole  )

  end subroutine regrid_geo2phys_3d

  !-----------------------------------------------------------------------------
  ! map 3D fields from physics grid to magnetic
  !-----------------------------------------------------------------------------
  subroutine regrid_phys2mag_3d( physfld, magfld, plev, cols, cole )
    integer, intent(in)   :: plev, cols, cole
    real(r8), intent(in)  :: physfld(1:plev,cols:cole)
    real(r8), intent(out) :: magfld(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1)

    call edyn_esmf_set3d_phys( phys_3dfld, physfld,  1, plev, cols, cole)
    call edyn_esmf_regrid_phys2mag(phys_3dfld, mag_des_3dfld, 3)
    call edyn_esmf_get_3dfield(mag_des_3dfld, magfld, mlon0, mlon1, mlat0, mlat1, mlev0, mlev1 )

  end subroutine regrid_phys2mag_3d

end module regridder
