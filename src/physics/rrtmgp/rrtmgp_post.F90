module rrtmgp_post

  use ccpp_kinds, only: kind_phys
  use ccpp_optical_props,      only: ty_optical_props_1scl_ccpp, ty_optical_props_2str_ccpp
  use ccpp_source_functions,   only: ty_source_func_lw_ccpp
  use ccpp_fluxes,             only: ty_fluxes_broadband_ccpp
  use ccpp_fluxes_byband,      only: ty_fluxes_byband_ccpp

  public :: rrtmgp_post_run

contains
!> \section arg_table_rrtmgp_post_run Argument Table
!! \htmlinclude rrtmgp_post_run.html
!!
subroutine rrtmgp_post_run(ncol, qrs, qrl, pdel, atm_optics_sw, cloud_sw, aer_sw, &
                  fsw, fswc, sources_lw, cloud_lw, aer_lw, flw, flwc, errmsg, errflg)
   integer,                          intent(in)    :: ncol
   real(kind_phys), dimension(:,:),  intent(in)    :: pdel
   real(kind_phys), dimension(:,:),  intent(inout) :: qrs
   real(kind_phys), dimension(:,:),  intent(inout) :: qrl
   type(ty_optical_props_2str_ccpp), intent(inout) :: atm_optics_sw
   type(ty_optical_props_1scl_ccpp), intent(inout) :: aer_lw
   type(ty_optical_props_2str_ccpp), intent(inout) :: aer_sw
   type(ty_optical_props_1scl_ccpp), intent(inout) :: cloud_lw
   type(ty_optical_props_2str_ccpp), intent(inout) :: cloud_sw
   type(ty_fluxes_broadband_ccpp),   intent(inout) :: fswc
   type(ty_fluxes_broadband_ccpp),   intent(inout) :: flwc
   type(ty_fluxes_byband_ccpp),      intent(inout) :: fsw
   type(ty_fluxes_byband_ccpp),      intent(inout) :: flw
   type(ty_source_func_lw_ccpp),     intent(inout) :: sources_lw
   character(len=*),                 intent(out)   :: errmsg
   integer,                          intent(out)   :: errflg

   ! Set error varaibles
   errflg = 0
   errmsg = ''
   ! The radiative heating rates are carried in the physics buffer across timesteps
   ! as Q*dp (for energy conservation).
   qrs(:ncol,:) = qrs(:ncol,:) * pdel(:ncol,:)
   qrl(:ncol,:) = qrl(:ncol,:) * pdel(:ncol,:)
   call free_optics_sw(atm_optics_sw)
   call free_optics_sw(cloud_sw)
   call free_optics_sw(aer_sw)
   call free_fluxes_byband(fsw)
   call free_fluxes_broadband(fswc)

   call sources_lw%sources%finalize()
   call free_optics_lw(cloud_lw)
   call free_optics_lw(aer_lw)
   call free_fluxes_byband(flw)
   call free_fluxes_broadband(flwc)

end subroutine rrtmgp_post_run

  !=========================================================================================

subroutine free_optics_sw(optics)

   type(ty_optical_props_2str_ccpp), intent(inout) :: optics

   if (allocated(optics%optical_props%tau)) deallocate(optics%optical_props%tau)
   if (allocated(optics%optical_props%ssa)) deallocate(optics%optical_props%ssa)
   if (allocated(optics%optical_props%g)) deallocate(optics%optical_props%g)
   call optics%optical_props%finalize()

end subroutine free_optics_sw

!=========================================================================================

subroutine free_optics_lw(optics)

   type(ty_optical_props_1scl_ccpp), intent(inout) :: optics

   if (allocated(optics%optical_props%tau)) deallocate(optics%optical_props%tau)
   call optics%optical_props%finalize()

end subroutine free_optics_lw

!=========================================================================================

subroutine free_fluxes_broadband(fluxes)

   class(ty_fluxes_broadband_ccpp), intent(inout) :: fluxes

   if (associated(fluxes%fluxes%flux_up)) deallocate(fluxes%fluxes%flux_up)
   if (associated(fluxes%fluxes%flux_dn)) deallocate(fluxes%fluxes%flux_dn)
   if (associated(fluxes%fluxes%flux_net)) deallocate(fluxes%fluxes%flux_net)
   if (associated(fluxes%fluxes%flux_dn_dir)) deallocate(fluxes%fluxes%flux_dn_dir)

end subroutine free_fluxes_broadband

!=========================================================================================

subroutine free_fluxes_byband(fluxes)

   class(ty_fluxes_byband_ccpp), intent(inout) :: fluxes

   if (associated(fluxes%fluxes%flux_up)) deallocate(fluxes%fluxes%flux_up)
   if (associated(fluxes%fluxes%flux_dn)) deallocate(fluxes%fluxes%flux_dn)
   if (associated(fluxes%fluxes%flux_net)) deallocate(fluxes%fluxes%flux_net)
   if (associated(fluxes%fluxes%flux_dn_dir)) deallocate(fluxes%fluxes%flux_dn_dir)

   if (associated(fluxes%fluxes%bnd_flux_up)) deallocate(fluxes%fluxes%bnd_flux_up)
   if (associated(fluxes%fluxes%bnd_flux_dn)) deallocate(fluxes%fluxes%bnd_flux_dn)
   if (associated(fluxes%fluxes%bnd_flux_net)) deallocate(fluxes%fluxes%bnd_flux_net)
   if (associated(fluxes%fluxes%bnd_flux_dn_dir)) deallocate(fluxes%fluxes%bnd_flux_dn_dir)

end subroutine free_fluxes_byband

end module rrtmgp_post
