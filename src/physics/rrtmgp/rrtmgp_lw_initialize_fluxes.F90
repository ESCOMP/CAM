module rrtmgp_lw_initialize_fluxes

  public :: rrtmgp_lw_initialize_fluxes_run

contains
!> \section arg_table_rrtmgp_lw_initialize_fluxes_run Argument Table
!! \htmlinclude rrtmgp_lw_initialize_fluxes_run.html
!!
  subroutine rrtmgp_lw_initialize_fluxes_run(rrtmgp_phys_blksz, nlay, nlwbands, spectralflux, flux_allsky, flux_clrsky, &
                  errmsg, errflg)
   use mo_fluxes,        only: ty_fluxes_broadband
   use mo_fluxes_byband, only: ty_fluxes_byband
   ! Inputs
   integer, intent(in) :: rrtmgp_phys_blksz
   integer, intent(in) :: nlay
   integer, intent(in) :: nlwbands
   logical, intent(in) :: spectralflux
   ! Outputs
   class(ty_fluxes_broadband), intent(out) :: flux_clrsky
   class(ty_fluxes_broadband), intent(out) :: flux_allsky
   character(len=*),           intent(out) :: errmsg
   integer,                    intent(out) :: errflg

   ! Local variables
   character(len=256) :: alloc_errmsg
   integer            :: play

   play = nlay + 1

   ! Set error variables
   errmsg = ''
   errflg = 0

   ! Clearsky fluxes
   allocate(flux_clrsky%flux_up(rrtmgp_phys_blksz, play), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a)') 'rrtmgp_lw_initialize_fluxes_run: ERROR: failed to allocate "flux_clrsky%flux_up". Message: ', &
              alloc_errmsg
      return
   end if

   allocate(flux_clrsky%flux_dn(rrtmgp_phys_blksz, play), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a)') 'rrtmgp_lw_initialize_fluxes_run: ERROR: failed to allocate "flux_clrsky%flux_dn". Message: ', &
              alloc_errmsg
      return
   end if

   allocate(flux_clrsky%flux_net(rrtmgp_phys_blksz, play), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a)') 'rrtmgp_lw_initialize_fluxes_run: ERROR: failed to allocate "flux_clrsky%flux_net". Message: ', &
              alloc_errmsg
      return
   end if

   select type (flux_clrsky)
   type is (ty_fluxes_byband)
      ! Only allocate when spectralflux is true.
      if (spectralflux) then
         allocate(flux_clrsky%bnd_flux_up(rrtmgp_phys_blksz, play, nlwbands), stat=errflg, errmsg=alloc_errmsg)
         if (errflg /= 0) then
            write(errmsg, '(a,a)') 'rrtmgp_lw_initialize_fluxes_run: ERROR: failed to allocate "flux_clrsky%bnd_flux_up". Message: ', &
                    alloc_errmsg
            return
         end if

         allocate(flux_clrsky%bnd_flux_dn(rrtmgp_phys_blksz, play, nlwbands), stat=errflg, errmsg=alloc_errmsg)
         if (errflg /= 0) then
            write(errmsg, '(a,a)') 'rrtmgp_lw_initialize_fluxes_run: ERROR: failed to allocate "flux_clrsky%bnd_flux_dn". Message: ', &
                    alloc_errmsg
            return
         end if

         allocate(flux_clrsky%bnd_flux_net(rrtmgp_phys_blksz, play, nlwbands), stat=errflg, errmsg=alloc_errmsg)
         if (errflg /= 0) then
            write(errmsg, '(a,a)') 'rrtmgp_lw_initialize_fluxes_run: ERROR: failed to allocate "flux_clrsky%bnd_flux_net". Message: ', &
                    alloc_errmsg
            return
         end if
      end if
   end select

   ! Initialize
   call reset_fluxes(flux_clrsky)

   ! Allsky fluxes
   allocate(flux_allsky%flux_up(rrtmgp_phys_blksz, play), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a)') 'rrtmgp_lw_initialize_fluxes_run: ERROR: failed to allocate "flux_allsky%flux_up". Message: ', &
              alloc_errmsg
      return
   end if

   allocate(flux_allsky%flux_dn(rrtmgp_phys_blksz, play), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a)') 'rrtmgp_lw_initialize_fluxes_run: ERROR: failed to allocate "flux_allsky%flux_dn". Message: ', &
              alloc_errmsg
      return
   end if

   allocate(flux_allsky%flux_net(rrtmgp_phys_blksz, play), stat=errflg, errmsg=alloc_errmsg)
   if (errflg /= 0) then
      write(errmsg, '(a,a)') 'rrtmgp_lw_initialize_fluxes_run: ERROR: failed to allocate "flux_allsky%flux_net". Message: ', &
              alloc_errmsg
      return
   end if
!   if (do_direct_local) then
!      allocate(flux_allsky%flux_dn_dir(rrtmgp_phys_blksz, play), stat=errflg)
!      call handle_allocate_error(errflg, sub, 'flux_allsky%flux_dn_dir')
!   end if

   select type (flux_allsky)
   type is (ty_fluxes_byband)
      ! Fluxes by band always needed for SW.  Only allocate for LW
      ! when spectralflux is true.
      if (spectralflux) then
         allocate(flux_allsky%bnd_flux_up(rrtmgp_phys_blksz, play, nlwbands), stat=errflg, errmsg=alloc_errmsg)
         if (errflg /= 0) then
            write(errmsg, '(a,a)') 'rrtmgp_lw_initialize_fluxes_run: ERROR: failed to allocate "flux_allsky%bnd_flux_up". Message: ', &
                    alloc_errmsg
            return
         end if

         allocate(flux_allsky%bnd_flux_dn(rrtmgp_phys_blksz, play, nlwbands), stat=errflg, errmsg=alloc_errmsg)
         if (errflg /= 0) then
            write(errmsg, '(a,a)') 'rrtmgp_lw_initialize_fluxes_run: ERROR: failed to allocate "flux_allsky%bnd_flux_dn". Message: ', &
                    alloc_errmsg
            return
         end if

         allocate(flux_allsky%bnd_flux_net(rrtmgp_phys_blksz, play, nlwbands), stat=errflg, errmsg=alloc_errmsg)
         if (errflg /= 0) then
            write(errmsg, '(a,a)') 'rrtmgp_lw_initialize_fluxes_run: ERROR: failed to allocate "flux_allsky%bnd_flux_net". Message: ', &
                    alloc_errmsg
            return
         end if
   !      if (do_direct) then
   !         allocate(flux_allsky%bnd_flux_dn_dir(rrtmgp_phys_blksz, play, nlwbands), stat=errflg)
   !         call handle_allocate_error(errflg, sub, 'flux_allsky%bnd_flux_dn_dir')
   !         allocate(flux_allsky%bnd_flux_dn_dir(rrtmgp_phys_blksz, play, nbands), stat=errflg)
   !         call handle_allocate_error(errflg, sub, 'flux_allsky%bnd_flux_dn_dir')
   !      end if
      end if
   end select

   ! Initialize
   call reset_fluxes(flux_allsky)

  end subroutine rrtmgp_lw_initialize_fluxes_run

!=========================================================================================

  subroutine reset_fluxes(fluxes)
   use mo_fluxes,        only: ty_fluxes_broadband
   use mo_fluxes_byband, only: ty_fluxes_byband
   use ccpp_kinds,       only: kind_phys

   ! Reset flux arrays to zero.

   class(ty_fluxes_broadband), intent(inout) :: fluxes
   !----------------------------------------------------------------------------

   ! Reset broadband fluxes
   fluxes%flux_up(:,:) = 0._kind_phys
   fluxes%flux_dn(:,:) = 0._kind_phys
   fluxes%flux_net(:,:) = 0._kind_phys
   if (associated(fluxes%flux_dn_dir)) fluxes%flux_dn_dir(:,:) = 0._kind_phys

   select type (fluxes)
   type is (ty_fluxes_byband)
      ! Reset band-by-band fluxes
      if (associated(fluxes%bnd_flux_up)) fluxes%bnd_flux_up(:,:,:) = 0._kind_phys
      if (associated(fluxes%bnd_flux_dn)) fluxes%bnd_flux_dn(:,:,:) = 0._kind_phys
      if (associated(fluxes%bnd_flux_net)) fluxes%bnd_flux_net(:,:,:) = 0._kind_phys
      if (associated(fluxes%bnd_flux_dn_dir)) fluxes%bnd_flux_dn_dir(:,:,:) = 0._kind_phys
   end select

  end subroutine reset_fluxes

end module rrtmgp_lw_initialize_fluxes
