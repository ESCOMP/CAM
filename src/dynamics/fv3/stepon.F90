module stepon

    ! MODULE: stepon -- FV Dynamics specific time-stepping

    use shr_kind_mod,  only: r8 => shr_kind_r8
    use physics_types, only: physics_state, physics_tend
    use ppgrid,        only: begchunk, endchunk
    use perf_mod,      only: t_startf, t_stopf, t_barrierf
    use spmd_utils,    only: iam, masterproc, mpicom
    use dyn_comp,      only: dyn_import_t, dyn_export_t
    use dyn_grid,      only: mytile
    use time_manager,  only: get_step_size
    use dimensions_mod, only: qsize_tracer_idx_cam2dyn

    implicit none
    private

    public stepon_init   ! Initialization
    public stepon_run1   ! run method phase 1
    public stepon_run2   ! run method phase 2
    public stepon_run3   ! run method phase 3
    public stepon_final  ! Finalization

!=======================================================================
contains
!=======================================================================

subroutine stepon_init(dyn_in, dyn_out)

    ! ROUTINE: stepon_init -- Time stepping initialization

  use cam_history,    only: addfld, add_default, horiz_only
  use constituents,   only: pcnst, cnst_name, cnst_longname
  
  type (dyn_import_t), intent(inout)   :: dyn_in             ! Dynamics import container
  type (dyn_export_t), intent(inout)   :: dyn_out            ! Dynamics export container
  
   ! local variables
   integer :: m_cnst,m_cnst_ffsl
   !----------------------------------------------------------------------------
   ! These fields on dynamics grid are output before the call to d_p_coupling.
   do m_cnst = 1, pcnst
     m_cnst_ffsl=qsize_tracer_idx_cam2dyn(m_cnst)
     call addfld(trim(cnst_name(m_cnst))//'_ffsl',  (/ 'lev' /), 'I', 'kg/kg',   &
          trim(cnst_longname(m_cnst)), gridname='FFSLHIST')
     call addfld(trim(cnst_name(m_cnst))//'_mass_ffsl',  (/ 'lev' /), 'I', 'kg/kg',   &
          trim(cnst_longname(m_cnst))//'*dp', gridname='FFSLHIST')
   end do
   call addfld('U_ffsl'     ,(/ 'lev' /), 'I', 'm/s ','U wind on A grid after dynamics',gridname='FFSLHIST')
   call addfld('V_ffsl'     ,(/ 'lev' /), 'I', 'm/s ','V wind on A grid after dynamics',gridname='FFSLHIST')
   call addfld('U_ffsl_ns'     ,(/ 'lev' /), 'I', 'm/s ','U wind on NS grid after dynamics',gridname='FFSLHIST_NS')
   call addfld('V_ffsl_ew'     ,(/ 'lev' /), 'I', 'm/s ','V wind on EW grid after dynamics',gridname='FFSLHIST_EW')
   call addfld('T_ffsl'     ,(/ 'lev' /), 'I', 'K '  ,'T on A grid grid after dynamics'     ,gridname='FFSLHIST')
   call addfld('PS_ffsl', horiz_only,  'I', 'Pa', 'Surface pressure on A grid after dynamics',gridname='FFSLHIST')
   call addfld('PHIS_ffsl', horiz_only,  'I', 'Pa', 'Geopotential height on A grid after dynamics',gridname='FFSLHIST')


   ! Fields for initial condition files
   call addfld('U&IC',   (/ 'lev' /),  'I', 'm/s', 'Zonal wind',     gridname='FFSLHIST' )
   call addfld('V&IC',   (/ 'lev' /),  'I', 'm/s', 'Meridional wind',gridname='FFSLHIST' )
   ! Don't need to register U&IC V&IC as vector components since we don't interpolate IC files
   call add_default('U&IC',0, 'I')
   call add_default('V&IC',0, 'I')
    
   call addfld('PS&IC', horiz_only,  'I', 'Pa', 'Surface pressure',gridname='FFSLHIST')
   call addfld('PHIS&IC', horiz_only,  'I', 'Pa', 'PHIS on ffsl grid',gridname='FFSLHIST')
   call addfld('T&IC',  (/ 'lev' /), 'I', 'K',  'Temperature',            gridname='FFSLHIST')
   call add_default('PS&IC',0, 'I')
   call add_default('PHIS&IC',0, 'I')
   call add_default('T&IC    ',0, 'I')

   do m_cnst = 1,pcnst
      call addfld(trim(cnst_name(m_cnst))//'&IC', (/ 'lev' /), 'I', 'kg/kg', &
                  trim(cnst_longname(m_cnst)), gridname='FFSLHIST')
      call add_default(trim(cnst_name(m_cnst))//'&IC', 0, 'I')
   end do

  
end subroutine stepon_init

!=======================================================================

subroutine stepon_run1(dtime_out, phys_state, phys_tend, pbuf2d, dyn_in, dyn_out)

    ! ROUTINE: stepon_run1 -- Phase 1 of dynamics run method.

    use physics_buffer, only: physics_buffer_desc
    use dp_coupling,    only: d_p_coupling

    real(r8),             intent(out)   :: dtime_out   ! Time-step
    type (physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type (physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
    type (physics_buffer_desc), pointer :: pbuf2d(:,:)
    type (dyn_import_t),  intent(inout) :: dyn_in  ! Dynamics import container
    type (dyn_export_t),  intent(inout) :: dyn_out ! Dynamics export container

    dtime_out = get_step_size()

    call diag_dyn_out(dyn_out,'')

    !----------------------------------------------------------
    ! Move data into phys_state structure.
    !----------------------------------------------------------

    call t_barrierf('sync_d_p_coupling', mpicom)
    call t_startf('d_p_coupling')
    call d_p_coupling(phys_state, phys_tend, pbuf2d, dyn_out)
    call t_stopf('d_p_coupling')

end subroutine stepon_run1

!=======================================================================

subroutine stepon_run2(phys_state, phys_tend, dyn_in, dyn_out)

    ! ROUTINE: stepon_run2 -- second phase run method

    use dp_coupling, only: p_d_coupling
    use dyn_comp,    only: calc_tot_energy_dynamics

    type (physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type (physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
    type (dyn_import_t),  intent(inout) :: dyn_in  ! Dynamics import container
    type (dyn_export_t),  intent(inout) :: dyn_out ! Dynamics export container

    ! copy from phys structures -> dynamics structures

    call t_barrierf('sync_p_d_coupling', mpicom)
#if ( defined CALC_ENERGY )
    call calc_tot_energy_dynamics(dyn_in%atm, 'dED')
#endif
    call t_startf('p_d_coupling')
    call p_d_coupling(phys_state, phys_tend, dyn_in)
    call t_stopf('p_d_coupling')

#if ( defined CALC_ENERGY )
    call calc_tot_energy_dynamics(dyn_in%atm, 'dBD')
#endif
end subroutine stepon_run2

!=======================================================================

subroutine stepon_run3(dtime, cam_out, phys_state, dyn_in, dyn_out)

    use camsrfexch, only: cam_out_t
    use dyn_comp,   only: dyn_run

    real(r8), intent(in) :: dtime            ! Time-step
    type (physics_state), intent(in):: phys_state(begchunk:endchunk)
    type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
    type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
    type (cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)

    call t_barrierf('sync_dyn_run', mpicom)
    call t_startf('dyn_run')
    call dyn_run(dyn_out)
    call t_stopf('dyn_run')

end subroutine stepon_run3

!=======================================================================

subroutine stepon_final(dyn_in, dyn_out)

    ! ROUTINE: stepon_final -- Dynamics finalization

    use dyn_comp,   only: dyn_final

    type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
    type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container

    call t_startf('dyn_final')
    call dyn_final(dyn_in, dyn_out)
    call t_stopf('dyn_final')

end subroutine stepon_final

!=======================================================================

subroutine diag_dyn_out(dyn_in,suffx)

  use cam_history,            only: write_inithist, outfld, hist_fld_active, fieldname_len
  use constituents,           only: cnst_name, pcnst
  use dyn_grid,               only: mytile
  use fv_arrays_mod,          only: fv_atmos_type
  use dimensions_mod,         only: nlev
  
  type (dyn_export_t),  intent(in) :: dyn_in
  character*(*)      ,  intent(in) :: suffx ! suffix for "outfld" names
  
  
  ! local variables
  integer              :: is,ie,js,je, j, m_cnst,m_cnst_ffsl
  integer              :: idim
  character(len=fieldname_len) :: tfname

  type (fv_atmos_type),  pointer :: Atm(:)
  
  !----------------------------------------------------------------------------
  
  Atm=>dyn_in%atm
  
  is = Atm(mytile)%bd%is
  ie = Atm(mytile)%bd%ie
  js = Atm(mytile)%bd%js
  je = Atm(mytile)%bd%je

  idim=ie-is+1
  ! Output tracer fields for analysis of advection schemes
  do m_cnst = 1, pcnst
     m_cnst_ffsl=qsize_tracer_idx_cam2dyn(m_cnst)
     tfname = trim(cnst_name(m_cnst))//'_ffsl'//trim(suffx)
     if (hist_fld_active(tfname)) then
        do j = js, je
           call outfld(tfname, RESHAPE(Atm(mytile)%q(is:ie, j, :, m_cnst_ffsl),(/idim,nlev/)), idim, j)
        end do
     end if
  end do

  ! Output tracer fields for analysis of advection schemes
  do m_cnst = 1, pcnst
     m_cnst_ffsl=qsize_tracer_idx_cam2dyn(m_cnst)
     tfname = trim(cnst_name(m_cnst))//'_mass_ffsl'//trim(suffx)
     if (hist_fld_active(tfname)) then
        do j = js, je
           call outfld(tfname,RESHAPE((Atm(mytile)%q(is:ie,j,:,m_cnst_ffsl)*Atm(mytile)%delp(is:ie,j,:)),(/idim,nlev/)),idim, j)
        end do
     end if
  end do
  
  if (hist_fld_active('U_ffsl'//trim(suffx)) .or. hist_fld_active('V_ffsl'//trim(suffx))) then
     do j = js, je
        call outfld('U_ffsl'//trim(suffx), RESHAPE(Atm(mytile)%ua(is:ie, j, :),(/idim,nlev/)), idim, j)
        call outfld('V_ffsl'//trim(suffx), RESHAPE(Atm(mytile)%va(is:ie, j, :),(/idim,nlev/)), idim, j)
     end do
  end if

  if (hist_fld_active('U_ffsl_ns'//trim(suffx))) then
     do j = js, je+1
        call outfld('U_ffsl_ns'//trim(suffx), RESHAPE(Atm(mytile)%u(is:ie, j, :),(/idim,nlev/)), idim, j)
     end do
  end if

  if (hist_fld_active('V_ffsl_ew'//trim(suffx))) then
     do j = js, je
        call outfld('V_ffsl_ew'//trim(suffx), RESHAPE(Atm(mytile)%v(is:ie+1, j, :),(/idim+1,nlev/)), idim+1, j)
     end do
  end if
  
  if (hist_fld_active('T_ffsl'//trim(suffx))) then
     do j = js, je
        call outfld('T_ffsl'//trim(suffx), RESHAPE(Atm(mytile)%pt(is:ie, j, :),(/idim,nlev/)), idim, j)
     end do
  end if

  if (hist_fld_active('PS_ffsl'//trim(suffx))) then
     do j = js, je
        call outfld('PS_ffsl'//trim(suffx), Atm(mytile)%ps(is:ie, j), idim, j)
     end do
  end if
  
  if (hist_fld_active('PHIS_ffsl'//trim(suffx))) then
     do j = js, je
        call outfld('PHIS_ffsl'//trim(suffx), Atm(mytile)%phis(is:ie, j), idim, j)
     end do
  end if
  
  if (write_inithist()) then
     
     do j = js, je
        call outfld('T&IC', RESHAPE(Atm(mytile)%pt(is:ie, j, :),(/idim,nlev/)), idim, j)
        call outfld('U&IC', RESHAPE(Atm(mytile)%ua(is:ie, j, :),(/idim,nlev/)), idim, j)
        call outfld('V&IC', RESHAPE(Atm(mytile)%va(is:ie, j, :),(/idim,nlev/)), idim, j)
        call outfld('PS&IC', Atm(mytile)%ps(is:ie, j), idim, j)
        call outfld('PHIS&IC', Atm(mytile)%phis(is:ie, j), idim, j)
        
        do m_cnst = 1, pcnst
           m_cnst_ffsl=qsize_tracer_idx_cam2dyn(m_cnst)
           call outfld(trim(cnst_name(m_cnst))//'&IC', RESHAPE(Atm(mytile)%q(is:ie, j, :, m_cnst_ffsl),(/idim,nlev/)), idim, j)
        end do
     end do
  end if  ! if (write_inithist)

end subroutine diag_dyn_out

end module stepon
