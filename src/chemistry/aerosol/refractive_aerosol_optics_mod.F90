module refractive_aerosol_optics_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use aerosol_optics_mod, only: aerosol_optics
  use physconst,  only: rhoh2o
  use aerosol_state_mod, only: aerosol_state
  use aerosol_properties_mod, only: aerosol_properties

  use table_interp_mod, only: table_interp, table_interp_wghts, table_interp_calcwghts

  implicit none

  private
  public :: refractive_aerosol_optics

  !> refractive_aerosol_optics
  !! Table look up implementation of aerosol_optics to parameterize aerosol radiative properties in terms of
  !! surface mode wet radius and wet refractive index using chebychev polynomials
  type, extends(aerosol_optics) :: refractive_aerosol_optics

     integer :: ibin, ilist
     class(aerosol_state), pointer :: aero_state      ! aerosol_state object
     class(aerosol_properties), pointer :: aero_props ! aerosol_properties object

     real(r8), allocatable :: watervol(:,:)   ! volume concentration of water in each mode (m3/kg)
     real(r8), allocatable :: wetvol(:,:)     ! volume concentration of wet mode (m3/kg)
     real(r8), allocatable :: cheb(:,:,:)     ! chebychev polynomials
     real(r8), allocatable :: radsurf(:,:)    ! aerosol surface mode radius
     real(r8), allocatable :: logradsurf(:,:) ! log(aerosol surface mode radius)

     ! refractive index for water read in read_water_refindex
     complex(r8), allocatable :: crefwsw(:) ! complex refractive index for water visible
     complex(r8), allocatable :: crefwlw(:) ! complex refractive index for water infrared

     real(r8), pointer :: extpsw(:,:,:,:) => null() ! specific extinction
     real(r8), pointer :: abspsw(:,:,:,:) => null() ! specific absorption
     real(r8), pointer :: asmpsw(:,:,:,:) => null() ! asymmetry factor
     real(r8), pointer :: absplw(:,:,:,:) => null() ! specific absorption

     real(r8), pointer :: refrtabsw(:,:) => null()  ! table of real refractive indices for aerosols
     real(r8), pointer :: refitabsw(:,:) => null()  ! table of imag refractive indices for aerosols
     real(r8), pointer :: refrtablw(:,:) => null()  ! table of real refractive indices for aerosols
     real(r8), pointer :: refitablw(:,:) => null()  ! table of imag refractive indices for aerosols

     ! Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
     ! in terms of refractive index and wet radius
     integer :: ncoef = -1  ! number of chebychev coeficients
     integer :: prefr = -1  ! number of real refractive indices
     integer :: prefi = -1  ! number of imaginary refractive indices

   contains

     procedure :: sw_props
     procedure :: lw_props

     final :: destructor

  end type refractive_aerosol_optics

  interface refractive_aerosol_optics
     procedure :: constructor
  end interface refractive_aerosol_optics

  ! radius limits (m)
  real(r8), parameter :: radmin = 0.01e-6_r8 ! min aerosol surface mode radius (m)
  real(r8), parameter :: radmax = 25.e-6_r8  ! max aerosol surface mode radius (m)
  real(r8), parameter :: xrmin=log(radmin)   ! min log(aerosol surface mode radius)
  real(r8), parameter :: xrmax=log(radmax)   ! max log(aerosol surface mode radius)

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor(aero_props, aero_state, ilist, ibin, ncol, nlev, nsw, nlw, crefwsw, crefwlw) &
       result(newobj)

    class(aerosol_properties),intent(in), target :: aero_props   ! aerosol_properties object
    class(aerosol_state),intent(in), target :: aero_state        ! aerosol_state object
    integer, intent(in) :: ilist  ! climate or a diagnostic list number
    integer, intent(in) :: ibin   ! bin number
    integer, intent(in) :: ncol   ! number of columns
    integer, intent(in) :: nlev   ! number of levels
    integer, intent(in) :: nsw    ! number of short wave lengths
    integer, intent(in) :: nlw    ! number of long wave lengths
    complex(r8), intent(in) :: crefwsw(nsw) ! complex refractive index for water visible
    complex(r8), intent(in) :: crefwlw(nlw) ! complex refractive index for water infrared

    type(refractive_aerosol_optics), pointer :: newobj

    integer :: ierr, icol, ilev, ispec, nspec
    real(r8) :: vol(ncol)             ! volume concentration of aerosol species (m3/kg)
    real(r8) :: dryvol(ncol)          ! volume concentration of aerosol mode (m3/kg)
    real(r8) :: specdens              ! species density (kg/m3)
    real(r8), pointer :: specmmr(:,:) ! species mass mixing ratio
    real(r8) :: logsigma              ! geometric standard deviation of number distribution

    real(r8) :: dgnumwet(ncol,nlev)   ! aerosol wet number mode diameter (m)
    real(r8) :: qaerwat(ncol,nlev)    ! aerosol water (g/g)

    real(r8), parameter :: rh2odens = 1._r8/rhoh2o

    allocate(newobj, stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    ! get mode properties
    call aero_props%optics_params(ilist, ibin, &
         refrtabsw=newobj%refrtabsw, refitabsw=newobj%refitabsw, &
         refrtablw=newobj%refrtablw, refitablw=newobj%refitablw,&
         extpsw=newobj%extpsw, abspsw=newobj%abspsw, asmpsw=newobj%asmpsw, &
         absplw=newobj%absplw, ncoef=newobj%ncoef, prefr=newobj%prefr, prefi=newobj%prefi)

    allocate(newobj%watervol(ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if
    allocate(newobj%wetvol(ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if
    allocate(newobj%cheb(newobj%ncoef,ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if
    allocate(newobj%radsurf(ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if
    allocate(newobj%logradsurf(ncol,nlev),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if

    allocate(newobj%crefwlw(nlw),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if
    newobj%crefwlw(:) = crefwlw(:)

    allocate(newobj%crefwsw(nsw),stat=ierr)
    if (ierr/=0) then
       nullify(newobj)
       return
    end if
    newobj%crefwsw(:) = crefwsw(:)

    call aero_state%water_uptake(aero_props, ilist, ibin,  ncol, nlev, dgnumwet, qaerwat)

    nspec = aero_props%nspecies(ilist,ibin)

    logsigma=aero_props%alogsig(ilist,ibin)

    ! calc size parameter for all columns
    call modal_size_parameters(newobj%ncoef, ncol, nlev, logsigma, dgnumwet, &
                               newobj%radsurf, newobj%logradsurf, newobj%cheb)

    do ilev = 1, nlev
       dryvol(:ncol) = 0._r8
       do ispec = 1, nspec
          call aero_state%get_ambient_mmr(ilist,ispec,ibin,specmmr)
          call aero_props%get(ibin, ispec, list_ndx=ilist, density=specdens)

          do icol = 1, ncol
             vol(icol) = specmmr(icol,ilev)/specdens
             dryvol(icol) = dryvol(icol) + vol(icol)

             newobj%watervol(icol,ilev) = qaerwat(icol,ilev)*rh2odens
             newobj%wetvol(icol,ilev) = newobj%watervol(icol,ilev) + dryvol(icol)
             if (newobj%watervol(icol,ilev) < 0._r8) then
                newobj%watervol(icol,ilev) = 0._r8
                newobj%wetvol(icol,ilev) = dryvol(icol)
             end if
          end do
       end do
    end do

    newobj%aero_state => aero_state
    newobj%aero_props => aero_props
    newobj%ilist = ilist
    newobj%ibin = ibin

  end function constructor

  !------------------------------------------------------------------------------
  ! returns short wave aerosol optics properties
  !------------------------------------------------------------------------------
  subroutine sw_props(self, ncol, ilev, iwav, pext, pabs, palb, pasm)

    class(refractive_aerosol_optics), intent(in) :: self
    integer, intent(in) :: ncol        ! number of columns
    integer, intent(in) :: ilev        ! vertical level index
    integer, intent(in) :: iwav        ! wave length index
    real(r8),intent(out) :: pext(ncol) ! parameterized specific extinction (m2/kg)
    real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)
    real(r8),intent(out) :: palb(ncol) ! parameterized asymmetry factor
    real(r8),intent(out) :: pasm(ncol) ! parameterized single scattering albedo

    real(r8) :: refr(ncol)      ! real part of refractive index
    real(r8) :: refi(ncol)      ! imaginary part of refractive index
    real(r8) :: cext(self%ncoef,ncol), cabs(self%ncoef,ncol), casm(self%ncoef,ncol)

    complex(r8) :: crefin(ncol) ! complex refractive index
    integer :: icol,icoef

    type(table_interp_wghts) :: wghtsr(ncol)
    type(table_interp_wghts) :: wghtsi(ncol)

    crefin(:ncol) = self%aero_state%refractive_index_sw(ncol, ilev, self%ilist, self%ibin, iwav, self%aero_props)

    do icol = 1, ncol
       crefin(icol) = crefin(icol) + self%watervol(icol,ilev)*self%crefwsw(iwav)
       crefin(icol) = crefin(icol)/max(self%wetvol(icol,ilev),1.e-60_r8)
       refr(icol) = real(crefin(icol))
       refi(icol) = abs(aimag(crefin(icol)))
    end do

    ! interpolate coefficients linear in refractive index

    wghtsr = table_interp_calcwghts( self%prefr, self%refrtabsw(:,iwav), ncol, refr(:ncol) )
    wghtsi = table_interp_calcwghts( self%prefi, self%refitabsw(:,iwav), ncol, refi(:ncol) )

    cext(:,:ncol)= table_interp( self%ncoef,ncol, self%prefr,self%prefi, wghtsr,wghtsi, self%extpsw(:,:,:,iwav))
    cabs(:,:ncol)= table_interp( self%ncoef,ncol, self%prefr,self%prefi, wghtsr,wghtsi, self%abspsw(:,:,:,iwav))
    casm(:,:ncol)= table_interp( self%ncoef,ncol, self%prefr,self%prefi, wghtsr,wghtsi, self%asmpsw(:,:,:,iwav))

    do icol = 1,ncol

       if (self%logradsurf(icol,ilev) <= xrmax) then
          pext(icol) = 0.5_r8*cext(1,icol)
          do icoef = 2, self%ncoef
             pext(icol) = pext(icol) + self%cheb(icoef,icol,ilev)*cext(icoef,icol)
          enddo
          pext(icol) = exp(pext(icol))
       else
          pext(icol) = 1.5_r8/(self%radsurf(icol,ilev)*rhoh2o) ! geometric optics
       endif

       ! convert from m2/kg water to m2/kg aerosol
       pext(icol) = pext(icol)*self%wetvol(icol,ilev)*rhoh2o
       pabs(icol) = 0.5_r8*cabs(1,icol)
       pasm(icol) = 0.5_r8*casm(1,icol)
       do icoef = 2, self%ncoef
          pabs(icol) = pabs(icol) + self%cheb(icoef,icol,ilev)*cabs(icoef,icol)
          pasm(icol) = pasm(icol) + self%cheb(icoef,icol,ilev)*casm(icoef,icol)
       enddo
       pabs(icol) = pabs(icol)*self%wetvol(icol,ilev)*rhoh2o
       pabs(icol) = max(0._r8,pabs(icol))
       pabs(icol) = min(pext(icol),pabs(icol))

       palb(icol) = 1._r8-pabs(icol)/max(pext(icol),1.e-40_r8)

    end do

  end subroutine sw_props

  !------------------------------------------------------------------------------
  ! returns long wave aerosol optics properties
  !------------------------------------------------------------------------------
  subroutine lw_props(self, ncol, ilev, iwav, pabs)

    class(refractive_aerosol_optics), intent(in) :: self
    integer, intent(in) :: ncol        ! number of columns
    integer, intent(in) :: ilev        ! vertical level index
    integer, intent(in) :: iwav        ! wave length index
    real(r8),intent(out) :: pabs(ncol) ! parameterized specific absorption (m2/kg)

    real(r8) :: refr(ncol)      ! real part of refractive index
    real(r8) :: refi(ncol)      ! imaginary part of refractive index
    real(r8) :: cabs(self%ncoef,ncol)

    complex(r8) :: crefin(ncol) ! complex refractive index
    integer :: icol, icoef

    type(table_interp_wghts) :: wghtsr(ncol)
    type(table_interp_wghts) :: wghtsi(ncol)

    crefin(:ncol) = self%aero_state%refractive_index_lw(ncol, ilev, self%ilist, self%ibin, iwav, self%aero_props)

    do icol = 1, ncol
       crefin(icol) = crefin(icol) + self%watervol(icol,ilev)*self%crefwlw(iwav)
       crefin(icol) = crefin(icol)/max(self%wetvol(icol,ilev), 1.e-40_r8)

       refr(icol) = real(crefin(icol))
       refi(icol) = aimag(crefin(icol))

    end do

    ! interpolate coefficients linear in refractive index

    wghtsr = table_interp_calcwghts( self%prefr, self%refrtablw(:,iwav), ncol, refr(:ncol) )
    wghtsi = table_interp_calcwghts( self%prefi, self%refitablw(:,iwav), ncol, refi(:ncol) )

    cabs(:,:ncol)= table_interp( self%ncoef,ncol, self%prefr,self%prefi, wghtsr,wghtsi, self%absplw(:,:,:,iwav))

    do icol = 1,ncol
       pabs(icol) = 0.5_r8*cabs(1,icol)
       do icoef = 2, self%ncoef
          pabs(icol) = pabs(icol) + self%cheb(icoef,icol,ilev)*cabs(icoef,icol)
       end do
       pabs(icol) = pabs(icol)*self%wetvol(icol,ilev)*rhoh2o
       pabs(icol) = max(0._r8,pabs(icol))
    end do

  end subroutine lw_props


  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)

    type(refractive_aerosol_optics), intent(inout) :: self

    deallocate(self%watervol)
    deallocate(self%wetvol)
    deallocate(self%cheb)
    deallocate(self%radsurf)
    deallocate(self%logradsurf)
    deallocate(self%crefwsw)
    deallocate(self%crefwlw)

    nullify(self%aero_state)
    nullify(self%aero_props)
    nullify(self%extpsw)
    nullify(self%abspsw)
    nullify(self%asmpsw)
    nullify(self%absplw)
    nullify(self%refrtabsw)
    nullify(self%refitabsw)
    nullify(self%refrtablw)
    nullify(self%refitablw)

  end subroutine destructor


  ! Private routines
  !===============================================================================

  !===============================================================================

  subroutine modal_size_parameters(ncoef,ncol,nlev, alnsg_amode, dgnumwet, radsurf, logradsurf, cheb)

    integer,  intent(in)  :: ncoef,ncol,nlev
    real(r8), intent(in)  :: alnsg_amode     ! geometric standard deviation of number distribution
    real(r8), intent(in)  :: dgnumwet(:,:)   ! aerosol wet number mode diameter (m)
    real(r8), intent(out) :: radsurf(:,:)    ! aerosol surface mode radius
    real(r8), intent(out) :: logradsurf(:,:) ! log(aerosol surface mode radius)
    real(r8), intent(out) :: cheb(:,:,:)

    integer  :: i, k, nc
    real(r8) :: explnsigma
    real(r8) :: xrad(ncol) ! normalized aerosol radius

    !-------------------------------------------------------------------------------

    explnsigma = exp(2.0_r8*alnsg_amode*alnsg_amode)

    do k = 1, nlev
       do i = 1, ncol
          ! convert from number mode diameter to surface area
          radsurf(i,k) = max(0.5_r8*dgnumwet(i,k)*explnsigma,radmin)
          logradsurf(i,k) = log(radsurf(i,k))
          ! normalize size parameter
          xrad(i) = max(logradsurf(i,k),xrmin)
          xrad(i) = min(xrad(i),xrmax)
          xrad(i) = (2._r8*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
          ! chebyshev polynomials
          cheb(1,i,k) = 1._r8
          cheb(2,i,k) = xrad(i)
          do nc = 3, ncoef
             cheb(nc,i,k) = 2._r8*xrad(i)*cheb(nc-1,i,k)-cheb(nc-2,i,k)
          end do
       end do
    end do

  end subroutine modal_size_parameters

end module refractive_aerosol_optics_mod
