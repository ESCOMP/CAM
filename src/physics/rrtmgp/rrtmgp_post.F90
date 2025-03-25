module rrtmgp_post

  use ccpp_kinds, only: kind_phys
  use mo_optical_props,      only: ty_optical_props_1scl, ty_optical_props_2str
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_fluxes_byband,      only: ty_fluxes_byband

  public :: rrtmgp_post_run

contains
!> \section arg_table_rrtmgp_post_run Argument Table
!! \htmlinclude rrtmgp_post_run.html
!!
subroutine rrtmgp_post_run(ncol, qrs, qrl, pdel, atm_optics_sw, cloud_sw, aer_sw, &
                  fsw, fswc, sources_lw, cloud_lw, aer_lw, flw, flwc, errmsg, errflg)
   integer,                         intent(in)    :: ncol
   real(kind_phys), dimension(:,:), intent(in)    :: pdel
   real(kind_phys), dimension(:,:), intent(inout) :: qrs
   real(kind_phys), dimension(:,:), intent(inout) :: qrl
   type(ty_optical_props_2str),     intent(inout) :: atm_optics_sw
   type(ty_optical_props_1scl),     intent(inout) :: aer_lw
   type(ty_optical_props_2str),     intent(inout) :: aer_sw
   type(ty_optical_props_1scl),     intent(inout) :: cloud_lw
   type(ty_optical_props_2str),     intent(inout) :: cloud_sw
   type(ty_fluxes_broadband),       intent(inout) :: fswc
   type(ty_fluxes_broadband),       intent(inout) :: flwc
   type(ty_fluxes_byband),          intent(inout) :: fsw
   type(ty_fluxes_byband),          intent(inout) :: flw
   type(ty_source_func_lw),         intent(inout) :: sources_lw
   character(len=*),                intent(out)   :: errmsg
   integer,                         intent(out)   :: errflg

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
   call free_fluxes(fsw)
   call free_fluxes(fswc)

   call sources_lw%finalize()
   call free_optics_lw(cloud_lw)
   call free_optics_lw(aer_lw)
   call free_fluxes(flw)
   call free_fluxes(flwc)

end subroutine rrtmgp_post_run

  !=========================================================================================

subroutine free_optics_sw(optics)

   type(ty_optical_props_2str), intent(inout) :: optics

   if (allocated(optics%tau)) deallocate(optics%tau)
   if (allocated(optics%ssa)) deallocate(optics%ssa)
   if (allocated(optics%g)) deallocate(optics%g)
   call optics%finalize()

end subroutine free_optics_sw

!=========================================================================================

subroutine free_optics_lw(optics)

   type(ty_optical_props_1scl), intent(inout) :: optics

   if (allocated(optics%tau)) deallocate(optics%tau)
   call optics%finalize()

end subroutine free_optics_lw

!=========================================================================================

subroutine free_fluxes(fluxes)

   class(ty_fluxes_broadband), intent(inout) :: fluxes

   if (associated(fluxes%flux_up)) deallocate(fluxes%flux_up)
   if (associated(fluxes%flux_dn)) deallocate(fluxes%flux_dn)
   if (associated(fluxes%flux_net)) deallocate(fluxes%flux_net)
   if (associated(fluxes%flux_dn_dir)) deallocate(fluxes%flux_dn_dir)

   select type (fluxes)
   type is (ty_fluxes_byband)
      if (associated(fluxes%bnd_flux_up)) deallocate(fluxes%bnd_flux_up)
      if (associated(fluxes%bnd_flux_dn)) deallocate(fluxes%bnd_flux_dn)
      if (associated(fluxes%bnd_flux_net)) deallocate(fluxes%bnd_flux_net)
      if (associated(fluxes%bnd_flux_dn_dir)) deallocate(fluxes%bnd_flux_dn_dir)
   end select

end subroutine free_fluxes


end module rrtmgp_post
