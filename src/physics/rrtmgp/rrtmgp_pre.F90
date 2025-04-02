module rrtmgp_pre
 use ccpp_kinds,              only: kind_phys
 use ccpp_fluxes,             only: ty_fluxes_broadband_ccpp
 use ccpp_fluxes_byband,      only: ty_fluxes_byband_ccpp
 use atmos_phys_string_utils, only: to_lower

 public :: rrtmgp_pre_init
 public :: rrtmgp_pre_run
 public :: radiation_do_ccpp

CONTAINS

!> \section arg_table_rrtmgp_pre_init Argument Table
!! \htmlinclude rrtmgp_pre_init.html
!!
  subroutine rrtmgp_pre_init(nradgas, gaslist, available_gases, gaslist_lc, errmsg, errflg)
     integer,                             intent(in) :: nradgas          ! Number of radiatively active gases
     character(len=*),                    intent(in) :: gaslist          ! List of radiatively active gases
     type(ty_gas_concentrations_ccpp), intent(inout) :: available_gases  ! Gas concentrations object
     character(len=*),                   intent(out) :: gaslist_lc       ! Lowercase verison of radiatively active gas list
     character(len=512),                 intent(out) :: errmsg
     integer,                            intent(out) :: errflg

     ! Set error variables
     errmsg = ''
     errflg = 0

     ! Create lowercase version of the gaslist for RRTMGP.  The ty_gas_concs_ccpp objects
     ! work with CAM's uppercase names, but other objects that get input from the gas
     ! concs objects don't work.
     do i = 1, nradgas
        gaslist_lc(i) = to_lower(gaslist(i))
     end do

     errmsg = available_gases%gas_concs%init(gaslist_lc)
     if (len_trim(errmsg) /= 0) then
        errflg = 1
     end if

  end subroutine rrtmgp_pre_init

!> \section arg_table_rrtmgp_pre_run Argument Table
!! \htmlinclude rrtmgp_pre_run.html
!!
  subroutine rrtmgp_pre_run(coszrs, nstep, dtime, iradsw, iradlw, irad_always, ncol, &
                  nextsw_cday, idxday, nday, idxnite, nnite, dosw, dolw, nlay, nlwbands, &
                  nswbands, spectralflux, fsw, fswc, flw, flwc, errmsg, errflg)
     use time_manager,         only: get_curr_calday
     ! Inputs
     real(kind_phys), dimension(:),    intent(in) :: coszrs        ! Cosine solar zenith angle
     integer,                          intent(in) :: dtime         ! Timestep size [s]
     integer,                          intent(in) :: nstep         ! Timestep number
     integer,                          intent(in) :: iradsw        ! Freq. of shortwave radiation calc in time steps (positive) or hours (negative)
     integer,                          intent(in) :: iradlw        ! Freq. of longwave radiation calc in time steps (positive) or hours (negative)
     integer,                          intent(in) :: irad_always   ! Number of time steps to execute radiation continuously
     integer,                          intent(in) :: ncol          ! Number of columns
     integer,                          intent(in) :: nlay          ! Number of vertical layers
     integer,                          intent(in) :: nlwbands      ! Number of longwave bands
     integer,                          intent(in) :: nswbands      ! Number of shortwave bands
     logical,                          intent(in) :: spectralflux  ! Flag to calculate fluxes (up and down) per band
     ! Outputs
     class(ty_fluxes_broadband_ccpp), intent(out) :: fswc          ! Clear-sky shortwave flux object
     class(ty_fluxes_byband_ccpp),    intent(out) :: fsw           ! All-sky shortwave flux object
     class(ty_fluxes_broadband_ccpp), intent(out) :: flwc          ! Clear-sky longwave flux object
     class(ty_fluxes_byband_ccpp),    intent(out) :: flw           ! All-sky longwave flux object
     integer,                         intent(out) :: nday          ! Number of daylight columns
     integer,                         intent(out) :: nnite         ! Number of nighttime columns
     real(kind_phys),                 intent(out) :: nextsw_cday   ! The next calendar day during which radiation calculation will be performed
     integer, dimension(:),           intent(out) :: idxday        ! Indices of daylight columns
     integer, dimension(:),           intent(out) :: idxnite       ! Indices of nighttime columns
     logical,                         intent(out) :: dosw          ! Flag to do shortwave calculation
     logical,                         intent(out) :: dolw          ! Flag to do longwave calculation
     character(len=*),                intent(out) :: errmsg
     integer,                         intent(out) :: errflg

     ! Local variables
     integer :: idx
     integer :: offset
     integer :: nstep_next
     logical :: dosw_next
     real(kind_phys) :: caldayp1

     ! Set error variables
     errflg = 0
     errmsg = ''

     ! Gather night/day column indices.
     nday = 0
     nnite = 0
     do idx = 1, ncol
        if ( coszrs(idx) > 0.0_kind_phys ) then
           nday = nday + 1
           idxday(nday) = idx
        else
           nnite = nnite + 1
           idxnite(nnite) = idx
        end if
     end do

     ! Determine if we're going to do longwave and/or shortwave this timestep
     call radiation_do_ccpp('sw', nstep, iradsw, irad_always, dosw, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if
     call radiation_do_ccpp('lw', nstep, iradlw, irad_always, dolw, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if

     ! Get time of next radiation calculation - albedos will need to be
     ! calculated by each surface model at this time
     nextsw_cday = -1._kind_phys
     dosw_next = .false.
     offset = 0
     nstep_next = nstep
     do while (.not. dosw_next)
        nstep_next = nstep_next + 1
        offset = offset + dtime
        call radiation_do_ccpp('sw', nstep_next, iradsw, irad_always, dosw_next, errmsg, errflg)
        if (errflg /= 0) then
           return
        end if
        if (dosw_next) then
           nextsw_cday = get_curr_calday(offset=offset)
        end if
     end do
     if(nextsw_cday == -1._kind_phys) then
        errflg = 1
        errmsg = 'next calendar day with shortwave calculation not found'
        return
     end if

     ! determine if next radiation time-step not equal to next time-step
     if (nstep >= 1) then
        caldayp1 = get_curr_calday(offset=int(dtime))
        if (caldayp1 /= nextsw_cday) nextsw_cday = -1._kind_phys
     end if

     ! Allocate the flux arrays and init to zero.
     call initialize_rrtmgp_fluxes_byband(nday, nlay+1, nswbands, nswbands, spectralflux, fsw, errmsg, errflg, do_direct=.true.)
     if (errflg /= 0) then
        return
     end if
     call initialize_rrtmgp_fluxes_broadband(nday, nlay+1, nswbands, nswbands, spectralflux, fswc, errmsg, errflg, do_direct=.true.)
     if (errflg /= 0) then
        return
     end if
     call initialize_rrtmgp_fluxes_byband(ncol, nlay+1, nlwbands, nswbands, spectralflux, flw, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if
     call initialize_rrtmgp_fluxes_broadband(ncol, nlay+1, nlwbands, nswbands, spectralflux, flwc, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if

  end subroutine rrtmgp_pre_run

!================================================================================================

subroutine radiation_do_ccpp(op, nstep, irad, irad_always, radiation_do, errmsg, errflg)

   ! Return true if the specified operation is done this timestep.

   character(len=*),  intent(in) :: op             ! name of operation
   integer,           intent(in) :: nstep
   integer,           intent(in) :: irad
   integer,           intent(in) :: irad_always
   integer,          intent(out) :: errflg
   character(len=*), intent(out) :: errmsg
   logical,          intent(out) :: radiation_do   ! return value

   !-----------------------------------------------------------------------

   ! Set error variables
   errflg = 0
   errmsg = ''

   select case (op)
      case ('sw') ! do a shortwave heating calc this timestep?
         radiation_do = nstep == 0  .or.  irad == 1                     &
                     .or. (mod(nstep-1,irad) == 0  .and.  nstep /= 1)   &
                     .or. nstep <= irad_always
      case ('lw') ! do a longwave heating calc this timestep?
         radiation_do = nstep == 0  .or.  irad == 1                     &
                     .or. (mod(nstep-1,irad) == 0  .and.  nstep /= 1)   &
                     .or. nstep <= irad_always
      case default
         errflg = 1
         errmsg = 'radiation_do_ccpp: unknown operation:'//op
   end select

end subroutine radiation_do_ccpp

!=========================================================================================

subroutine initialize_rrtmgp_fluxes_broadband(ncol, nlevels, nbands, nswbands, spectralflux, fluxes, errmsg, errflg, do_direct)

   ! Allocate flux arrays and set values to zero.

   ! Arguments
   integer,                    intent(in)    :: ncol, nlevels, nbands, nswbands
   logical,                    intent(in)    :: spectralflux
   class(ty_fluxes_broadband_ccpp), intent(inout) :: fluxes
   logical, optional,          intent(in)    :: do_direct
   character(len=*),           intent(out)   :: errmsg
   integer,                    intent(out)   :: errflg

   ! Local variables
   logical :: do_direct_local
   character(len=256) :: alloc_errmsg
   character(len=*), parameter :: sub = 'initialize_rrtmgp_fluxes'
   !----------------------------------------------------------------------------

   if (present(do_direct)) then
      do_direct_local = .true.
   else
      do_direct_local = .false.
   end if

   ! Broadband fluxes
   allocate(fluxes%fluxes%flux_up(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_up". Message: ', &
              alloc_errmsg
      return
   end if
   allocate(fluxes%fluxes%flux_dn(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_dn". Message: ', &
              alloc_errmsg
      return
   end if
   allocate(fluxes%fluxes%flux_net(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_net". Message: ', &
              alloc_errmsg
      return
   end if
   if (do_direct_local) then
      allocate(fluxes%fluxes%flux_dn_dir(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
      if (errflg /= 0) then
         write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_dn_dir". Message: ', &
                 alloc_errmsg
         return
      end if
   end if

   ! Initialize
   call reset_fluxes_broadband(fluxes)

end subroutine initialize_rrtmgp_fluxes_broadband

!=========================================================================================

subroutine initialize_rrtmgp_fluxes_byband(ncol, nlevels, nbands, nswbands, spectralflux, fluxes, errmsg, errflg, do_direct)

   ! Allocate flux arrays and set values to zero.

   ! Arguments
   integer,                    intent(in)    :: ncol, nlevels, nbands, nswbands
   logical,                    intent(in)    :: spectralflux
   class(ty_fluxes_byband_ccpp),    intent(inout) :: fluxes
   logical, optional,          intent(in)    :: do_direct
   character(len=*),           intent(out)   :: errmsg
   integer,                    intent(out)   :: errflg

   ! Local variables
   logical :: do_direct_local
   character(len=256) :: alloc_errmsg
   character(len=*), parameter :: sub = 'initialize_rrtmgp_fluxes_byband'
   !----------------------------------------------------------------------------

   if (present(do_direct)) then
      do_direct_local = .true.
   else
      do_direct_local = .false.
   end if

   ! Broadband fluxes
   allocate(fluxes%fluxes%flux_up(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_up". Message: ', &
              alloc_errmsg
      return
   end if
   allocate(fluxes%fluxes%flux_dn(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_dn". Message: ', &
              alloc_errmsg
      return
   end if
   allocate(fluxes%fluxes%flux_net(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_net". Message: ', &
              alloc_errmsg
      return
   end if
   if (do_direct_local) then
      allocate(fluxes%fluxes%flux_dn_dir(ncol, nlevels), stat=errflg, errmsg=alloc_errmsg)
      if (errflg /= 0) then
         write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%flux_dn_dir". Message: ', &
                 alloc_errmsg
         return
      end if
   end if

   ! Fluxes by band always needed for SW.  Only allocate for LW
   ! when spectralflux is true.
   if (nbands == nswbands .or. spectralflux) then
      allocate(fluxes%fluxes%bnd_flux_up(ncol, nlevels, nbands), stat=errflg, errmsg=alloc_errmsg)
      if (errflg /= 0) then
         write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%bnd_flux_up". Message: ', &
                alloc_errmsg
         return
      end if
      allocate(fluxes%fluxes%bnd_flux_dn(ncol, nlevels, nbands), stat=errflg, errmsg=alloc_errmsg)
      if (errflg /= 0) then
         write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%bnd_flux_dn". Message: ', &
                 alloc_errmsg
         return
      end if
      allocate(fluxes%fluxes%bnd_flux_net(ncol, nlevels, nbands), stat=errflg, errmsg=alloc_errmsg)
      if (errflg /= 0) then
         write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%bnd_flux_net". Message: ', &
                 alloc_errmsg
         return
      end if
      if (do_direct_local) then
         allocate(fluxes%fluxes%bnd_flux_dn_dir(ncol, nlevels, nbands), stat=errflg, errmsg=alloc_errmsg)
         if (errflg /= 0) then
            write(errmsg, '(a,a,a)') sub, ': ERROR: failed to allocate "fluxes%fluxes%bnd_flux_dn_dir". Message: ', &
                    alloc_errmsg
            return
         end if
      end if
   end if

   ! Initialize
   call reset_fluxes_byband(fluxes)

end subroutine initialize_rrtmgp_fluxes_byband

!=========================================================================================

subroutine reset_fluxes_broadband(fluxes)

   ! Reset flux arrays to zero.

   class(ty_fluxes_broadband_ccpp), intent(inout) :: fluxes
   !----------------------------------------------------------------------------

   ! Reset broadband fluxes
   fluxes%fluxes%flux_up(:,:) = 0._kind_phys
   fluxes%fluxes%flux_dn(:,:) = 0._kind_phys
   fluxes%fluxes%flux_net(:,:) = 0._kind_phys
   if (associated(fluxes%fluxes%flux_dn_dir)) fluxes%fluxes%flux_dn_dir(:,:) = 0._kind_phys

end subroutine reset_fluxes_broadband

!=========================================================================================

subroutine reset_fluxes_byband(fluxes)

   ! Reset flux arrays to zero.

   class(ty_fluxes_byband_ccpp), intent(inout) :: fluxes
   !----------------------------------------------------------------------------

   ! Reset broadband fluxes
   fluxes%fluxes%flux_up(:,:) = 0._kind_phys
   fluxes%fluxes%flux_dn(:,:) = 0._kind_phys
   fluxes%fluxes%flux_net(:,:) = 0._kind_phys
   if (associated(fluxes%fluxes%flux_dn_dir)) fluxes%fluxes%flux_dn_dir(:,:) = 0._kind_phys

   ! Reset band-by-band fluxes
   if (associated(fluxes%fluxes%bnd_flux_up)) fluxes%fluxes%bnd_flux_up(:,:,:) = 0._kind_phys
   if (associated(fluxes%fluxes%bnd_flux_dn)) fluxes%fluxes%bnd_flux_dn(:,:,:) = 0._kind_phys
   if (associated(fluxes%fluxes%bnd_flux_net)) fluxes%fluxes%bnd_flux_net(:,:,:) = 0._kind_phys
   if (associated(fluxes%fluxes%bnd_flux_dn_dir)) fluxes%fluxes%bnd_flux_dn_dir(:,:,:) = 0._kind_phys

end subroutine reset_fluxes_byband

!=========================================================================================

end module rrtmgp_pre
