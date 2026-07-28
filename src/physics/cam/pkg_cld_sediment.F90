! CAM interface for CCPPized cloud_particle_sedimentation
module pkg_cld_sediment

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  save

  public :: cld_sediment_readnl

  ! namelist variables
  real(r8) :: cldsed_ice_stokes_fac = huge(1._r8)   ! factor applied to the ice fall velocity computed from
  ! stokes terminal velocity

!===============================================================================
contains
!===============================================================================

  subroutine cld_sediment_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use units, only: getunit, freeunit
    use mpishorthand
    use cam_abortutils, only: endrun
    use cam_logfile,  only: iulog
    use spmd_utils, only: masterproc
    use physconst, only: gravit, rhoh2o

    use cloud_particle_sedimentation, only: cloud_particle_sedimentation_init

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'cld_sediment_readnl'
    character(len=512) :: errmsg
    integer            :: errflg

    namelist /cldsed_nl/ cldsed_ice_stokes_fac
    !-----------------------------------------------------------------------------

    if (masterproc) then
      unitn = getunit()
      open (unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, 'cldsed_nl', status=ierr)
      if (ierr == 0) then
        read (unitn, cldsed_nl, iostat=ierr)
        if (ierr /= 0) then
          call endrun(subname//':: ERROR reading namelist')
        end if
      end if
      close (unitn)
      call freeunit(unitn)

    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(cldsed_ice_stokes_fac, 1, mpir8, 0, mpicom)
#endif

    ! Call CCPP-ized subroutine initialization
    call cloud_particle_sedimentation_init(&
      amIRoot = masterproc, iulog = iulog, &
      cldsed_ice_stokes_fac_in = cldsed_ice_stokes_fac, &
      rhoh2o = rhoh2o, gravit = gravit, &
      errmsg = errmsg, errflg = errflg)

  end subroutine cld_sediment_readnl
end module pkg_cld_sediment
