!-------------------------------------------------------------------------------
! Geometric mean radius parameterized optical properties for volcanic
! stratospheric aerosols
!-------------------------------------------------------------------------------
module volcrad_aerosol_optics_mod
  use shr_kind_mod, only: r8 => shr_kind_r8

  use aerosol_optics_mod, only: aerosol_optics
  use aerosol_properties_mod, only: aerosol_properties
  use aerosol_state_mod, only: aerosol_state

  implicit none

  private

  public :: volcrad_aerosol_optics

  type, extends(aerosol_optics) :: volcrad_aerosol_optics

     ! aerosol optics properties tables (from physprops files)
     real(r8), pointer :: r_sw_ext(:,:) => null()
     real(r8), pointer :: r_sw_scat(:,:) => null()
     real(r8), pointer :: r_sw_ascat(:,:) => null()
     real(r8), pointer :: r_lw_abs(:,:) => null()
     real(r8), pointer :: r_mu(:)

     ! from state
     real(r8), allocatable :: wmu(:,:) ! (-) weighting on left side values ! (pcols,pver)
     integer , allocatable :: kmu(:,:) ! index into rh mesh

     ! aerosol mass mixing ratio
     real(r8), pointer :: mmr(:,:)

   contains

     procedure :: sw_props
     procedure :: lw_props

     final :: destructor

  end type volcrad_aerosol_optics

  interface volcrad_aerosol_optics
     procedure :: constructor
  end interface volcrad_aerosol_optics

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(aero_props, aero_state, ibin, ncols, nlevs, geometric_radius) &
                result(newobj)
    class(aerosol_properties),intent(in) :: aero_props   ! aerosol_properties object
    class(aerosol_state),     intent(in) :: aero_state      ! aerosol_state object
    integer, intent(in) :: ibin   ! bin number
    integer, intent(in) :: ncols, nlevs
    real(r8),intent(in) :: geometric_radius(ncols,nlevs)

    type(volcrad_aerosol_optics), pointer :: newobj

    integer :: ierr, nmu, i, k
    real(r8) :: r_mu_min, r_mu_max, mutrunc, mu

    allocate(newobj, stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    allocate(newobj%wmu(ncols,nlevs), stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    allocate(newobj%kmu(ncols,nlevs), stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    ! optical properties tables
    call aero_props%optics_params(ibin, &
         r_sw_ext=newobj%r_sw_ext, &
         r_sw_scat=newobj%r_sw_scat, &
         r_sw_ascat=newobj%r_sw_ascat, &
         r_lw_abs=newobj%r_lw_abs, &
         r_mu=newobj%r_mu )

    call aero_state%get_ambient_mmr(species_ndx=1, bin_ndx=ibin, mmr=newobj%mmr)


! NOTE should try to use table_interp_mod utility !!!

    nmu = size(newobj%r_mu)
    r_mu_max = newobj%r_mu(nmu)
    r_mu_min = newobj%r_mu(1)

    do i = 1, ncols
       do k = 1, nlevs
          if(geometric_radius(i,k) > 0._r8) then
             mu = log(geometric_radius(i,k))
          else
             mu = 0._r8
          endif

          ASSOCIATE ( kmu=>newobj%kmu(i,k), wmu=>newobj%wmu(i,k), r_mu=>newobj%r_mu )

            mutrunc = max(min(mu,r_mu_max),r_mu_min)
            kmu = max(min(1 + (mutrunc-r_mu_min)/(r_mu_max-r_mu_min)*(nmu-1),nmu-1._r8),1._r8)
            wmu = max(min( (mutrunc -r_mu(kmu)) / (r_mu(kmu+1) - r_mu(kmu)) ,1._r8),0._r8)

          END ASSOCIATE

       end do
    end do

  end function constructor

  !------------------------------------------------------------------------------
  ! returns short wave aerosol optics properties
  !------------------------------------------------------------------------------
  subroutine sw_props(self, ncol, ilev, iwav, pext, pabs, palb, pasm)

    class(volcrad_aerosol_optics), intent(in) :: self

    integer, intent(in) :: ncol        ! number of columns
    integer, intent(in) :: ilev        ! vertical level index
    integer, intent(in) :: iwav        ! wave length index
    real(r8),intent(out) :: pext(ncol) ! parameterized specific extinction (m2/kg)
    real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)
    real(r8),intent(out) :: palb(ncol) ! parameterized single scattering albedo
    real(r8),intent(out) :: pasm(ncol) ! parameterized asymmetry factor

    real(r8) :: scat(ncol)
    real(r8) :: ascat(ncol)
    integer :: icol

    ! interpolate the properties tables
    do icol = 1, ncol

       pext(icol) = ((1._r8 - self%wmu(icol,ilev)) * self%r_sw_ext(iwav, self%kmu(icol,ilev)  ) + &
                             (self%wmu(icol,ilev)) * self%r_sw_ext(iwav, self%kmu(icol,ilev)+1))

       scat(icol) = ((1._r8 - self%wmu(icol,ilev)) * self%r_sw_scat(iwav, self%kmu(icol,ilev)  ) + &
                             (self%wmu(icol,ilev)) * self%r_sw_scat(iwav, self%kmu(icol,ilev)+1))

       ascat(icol) = ((1._r8 - self%wmu(icol,ilev)) * self%r_sw_ascat(iwav, self%kmu(icol,ilev)  ) + &
                              (self%wmu(icol,ilev)) * self%r_sw_ascat(iwav, self%kmu(icol,ilev)+1))


       palb(icol) = scat(icol) / pext(icol)

       if (scat(icol)>0._r8) then
          pasm(icol) = ascat(icol) / scat(icol)
       else
          pasm(icol) = 0._r8
       end if

       pext(icol) = pext(icol) * self%mmr(icol,ilev)

       pabs(icol) = pext(icol) * ( 1._r8 - palb(icol) )

    end do

  end subroutine sw_props

  !------------------------------------------------------------------------------
  ! returns long wave aerosol optics properties
  !------------------------------------------------------------------------------
  subroutine lw_props(self, ncol, ilev, iwav, pabs)

    class(volcrad_aerosol_optics), intent(in) :: self
    integer, intent(in) :: ncol        ! number of columns
    integer, intent(in) :: ilev        ! vertical level index
    integer, intent(in) :: iwav        ! wave length index
    real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)

    integer :: icol

    ! interpolate the properties tables
    do icol = 1, ncol
       pabs(icol) = ((1._r8 - self%wmu(icol,ilev)) * self%r_lw_abs(iwav, self%kmu(icol,ilev)  ) + &
                             (self%wmu(icol,ilev)) * self%r_lw_abs(iwav, self%kmu(icol,ilev)+1))
       pabs(icol) = pabs(icol) * self%mmr(icol,ilev)
    end do

  end subroutine lw_props

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)

    type(volcrad_aerosol_optics), intent(inout) :: self

    deallocate(self%wmu)
    deallocate(self%kmu)

  end subroutine destructor

end module volcrad_aerosol_optics_mod
