module MCSP

  use ccpp_kinds, only:  kind_phys

  implicit none

  save
  private                         ! Make default type private to the module

!
! PUBLIC: interfaces
!
  public MCSP_init
  public MCSP_run                ! MCSP scheme

  real(kind_phys) heat_coeff
  real(kind_phys) moisture_coeff
  real(kind_phys) uwind_coeff
  real(kind_phys) vwind_coeff
  real(kind_phys) storm_speed_pref
  real(kind_phys) conv_depth_min
  real(kind_phys) shear_min

contains

!===============================================================================
!> \section arg_table_MCSP_init Argument Table
!! \htmlinclude MCSP_init.html
!!
subroutine MCSP_init(MCSP_heat_coeff, MCSP_moisture_coeff, MCSP_uwind_coeff, MCSP_vwind_coeff, &
                     MCSP_storm_speed_pref, MCSP_conv_depth_min, mcsp_shear_min,masterproc,iulog)

   real(kind_phys),                              intent(in   ) :: MCSP_heat_coeff         ! heating coefficient for MCSP
   real(kind_phys),                              intent(in   ) :: MCSP_moisture_coeff     ! moisture coefficient for MCSP
   real(kind_phys),                              intent(in   ) :: MCSP_uwind_coeff        ! uwind coefficient for MCSP
   real(kind_phys),                              intent(in   ) :: MCSP_vwind_coeff        ! vwind coefficient for MCSP
   real(kind_phys),                              intent(in   ) :: MCSP_storm_speed_pref   ! pressure level for zonal wind in MCSP calculation [Pa]
   real(kind_phys),                              intent(in   ) :: MCSP_conv_depth_min     ! pressure thickness of convective heating [Pa]
   real(kind_phys),                              intent(in   ) :: mcsp_shear_min          ! min shear value for MCSP to be active
   logical, intent(in)                                         :: masterproc
   integer, intent(in)                                         :: iulog

   heat_coeff = MCSP_heat_coeff
   moisture_coeff = MCSP_moisture_coeff
   uwind_coeff = MCSP_uwind_coeff
   vwind_coeff = MCSP_vwind_coeff
   storm_speed_pref = MCSP_storm_speed_pref
   conv_depth_min = MCSP_conv_depth_min
   shear_min = mcsp_shear_min

   if(masterproc) then
      write(iulog,*) 'MCSP_hear_coeff', heat_coeff
      write(iulog,*) 'MCSP_moisture_coeff', moisture_coeff
      write(iulog,*) 'MCSP_uwind_coeff', uwind_coeff
      write(iulog,*) 'MCSP_vwind_coeff', vwind_coeff
      write(iulog,*) 'MCSP_storm_speed_pref', storm_speed_pref
      write(iulog,*) 'MCSP_conv_depth_min', conv_depth_min
      write(iulog,*) 'mcsp_shear_min', shear_min
   end if

end subroutine MCSP_init

!===============================================================================
!> \section arg_table_MCSP_run Argument Table
!! \htmlinclude MCSP_run.html
!!
subroutine MCSP_run( pcols, ncol, pver, pverp, cpair, pi,        &
                              ztodt, jctop,                      &
                              pmid, pint, pdel,                  &
                              s, q, u, v,                        & 
                              tend_t, tend_q,                    &
                              ptend_s, ptend_q, ptend_u, ptend_v,                 &
                              mcsp_dt_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out, &
                              mcsp_freq, mcsp_shear, conv_depth, mcsp_tend_s_max )

   !----------------------------------------------------------------------------
   ! Purpose: perform MCSP tendency calculations
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                               intent(in   ) :: pcols      ! number of atmospheric columns (max)
   integer,                               intent(in   ) :: ncol       ! number of atmospheric columns (actual)
   integer,                               intent(in   ) :: pver       ! number of mid-point vertical levels
   integer,                               intent(in   ) :: pverp      ! number of interface vertical levels
   real(kind_phys),                       intent(in   ) :: cpair      ! specific heat of dry air (J K-1 kg-1)
   real(kind_phys),                       intent(in   ) :: pi         ! 4*atan(1.0)
   real(kind_phys),                       intent(in   ) :: ztodt      ! 2x physics time step
   integer,  dimension(pcols),            intent(in   ) :: jctop      ! cloud top level indices
   real(kind_phys), dimension(pcols,pver),       intent(in   ) :: pmid ! physics state mid-point pressure
   real(kind_phys), dimension(pcols,pverp),      intent(in   ) :: pint ! physics state interface pressure
   real(kind_phys), dimension(pcols,pver),       intent(in   ) :: pdel ! physics state pressure thickness
   real(kind_phys), dimension(pcols,pver),       intent(in   ) :: s    ! physics state dry static energy
   real(kind_phys), dimension(pcols,pver),       intent(in   ) :: q    ! physics state specific humidity
   real(kind_phys), dimension(pcols,pver),       intent(in   ) :: u    ! physics state u momentum
   real(kind_phys), dimension(pcols,pver),       intent(in   ) :: v    ! physics state v momentum
   real(kind_phys), dimension(pcols,pver),       intent(in   ) :: tend_t ! input deep convective temperature tendency 
   real(kind_phys), dimension(pcols,pver),       intent(in   ) :: tend_q ! input deep convective tendency for specific humidity (qv)
   real(kind_phys), dimension(pcols,pver),       intent(  out) :: ptend_s    ! output tendency of DSE
   real(kind_phys), dimension(pcols,pver),       intent(  out) :: ptend_q    ! output tendency of qv
   real(kind_phys), dimension(pcols,pver),       intent(  out) :: ptend_u    ! output tendency of u-wind
   real(kind_phys), dimension(pcols,pver),       intent(  out) :: ptend_v    ! output tendency of v-wind
   real(kind_phys), dimension(pcols,pver),       intent(  out) :: mcsp_dt_out! final MCSP tendency for DSE
   real(kind_phys), dimension(pcols,pver),       intent(  out) :: mcsp_dq_out! final MCSP tendency for qv
   real(kind_phys), dimension(pcols,pver),       intent(  out) :: mcsp_du_out! final MCSP tendency for u wind
   real(kind_phys), dimension(pcols,pver),       intent(  out) :: mcsp_dv_out! final MCSP tendency for v wind
   real(kind_phys), dimension(pcols),            intent(  out) :: mcsp_freq  ! MSCP frequency for output
   real(kind_phys), dimension(pcols),            intent(  out) :: mcsp_shear ! shear used to check against threshold
   real(kind_phys), dimension(pcols),            intent(  out) :: conv_depth   ! pressure depth of deep convection
   real(kind_phys), dimension(pcols),            intent(  out) :: mcsp_tend_s_max ! max MCSP heating tendency
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i, k

   real(kind_phys) :: tend_k       ! temporary variable for kinetic energy tendency
   real(kind_phys) :: pdepth_mid_k ! temporary pressure depth used for vertical structure
   real(kind_phys) :: pdepth_total ! temporary pressure depth used for vertical structure

   real(kind_phys), dimension(pcols)      :: avg_tend_s   ! mass weighted column average DSE tendency from ZM
   real(kind_phys), dimension(pcols)      :: avg_tend_q   ! mass weighted column average qv tendency from ZM
   real(kind_phys), dimension(pcols)      :: pdel_sum        ! column integrated pressure thickness

   real(kind_phys), dimension(pcols,pver) :: mcsp_tend_s     ! MCSP tendency before energy fixer for DSE
   real(kind_phys), dimension(pcols,pver) :: mcsp_tend_q     ! MCSP tendency before energy fixer for qv
   real(kind_phys), dimension(pcols,pver) :: mcsp_tend_u     ! MCSP tendency before energy fixer for u wind
   real(kind_phys), dimension(pcols,pver) :: mcsp_tend_v     ! MCSP tendency before energy fixer for v wind

   real(kind_phys), dimension(pcols)      :: mcsp_avg_tend_s ! mass weighted column average MCSP tendency of DSE
   real(kind_phys), dimension(pcols)      :: mcsp_avg_tend_q ! mass weighted column average MCSP tendency of qv 
   real(kind_phys), dimension(pcols)      :: mcsp_avg_tend_k ! mass weighted column average MCSP tendency of kinetic energy

   logical :: do_mcsp_t = .false.   ! internal flag to enable tendency calculations
   logical :: do_mcsp_q = .false.   ! internal flag to enable tendency calculations
   logical :: do_mcsp_u = .false.   ! internal flag to enable tendency calculations
   logical :: do_mcsp_v = .false.   ! internal flag to enable tendency calculations

   !----------------------------------------------------------------------------
   ! initialize variables

   if (heat_coeff>0) do_mcsp_t = .true.
   if (moisture_coeff>0) do_mcsp_q = .true.
   if (uwind_coeff>0) do_mcsp_u = .true.
   if (vwind_coeff>0) do_mcsp_v = .true.

   ptend_s = 0.0_kind_phys
   ptend_q = 0.0_kind_phys
   ptend_u = 0.0_kind_phys
   ptend_v = 0.0_kind_phys

   avg_tend_s(1:ncol) = 0
   avg_tend_q(1:ncol) = 0

   pdel_sum(1:ncol) = 0

   mcsp_avg_tend_s(1:ncol) = 0
   mcsp_avg_tend_q(1:ncol) = 0
   mcsp_avg_tend_k(1:ncol) = 0

   mcsp_tend_s(1:ncol,1:pver) = 0
   mcsp_tend_q(1:ncol,1:pver) = 0
   mcsp_tend_u(1:ncol,1:pver) = 0
   mcsp_tend_v(1:ncol,1:pver) = 0
   mcsp_tend_s_max = 0.0_kind_phys

   if (heat_coeff>0 .OR. moisture_coeff>0 .OR. uwind_coeff>0 .OR. vwind_coeff>0) then

      !----------------------------------------------------------------------------
      ! calculate shear

      call mcsp_calculate_shear( ncol, pcols, pver, pmid, u, mcsp_shear )

      !----------------------------------------------------------------------------
      ! calculate mass weighted column average tendencies from deep convection

      conv_depth(1:ncol) = 0
      do i = 1,ncol
         if (jctop(i) .ne. pver) then
            ! integrate pressure and deep convective tendencies over column
            do k = jctop(i),pver
               avg_tend_s(i) = avg_tend_s(i) + tend_t(i,k) * pdel(i,k) * cpair
               avg_tend_q(i) = avg_tend_q(i) + tend_q(i,k) * pdel(i,k)
               pdel_sum(i) = pdel_sum(i) + pdel(i,k)
            end do
            ! normalize integrated deep convective tendencies by total mass
            avg_tend_s(i) = avg_tend_s(i) / pdel_sum(i)
            avg_tend_q(i) = avg_tend_q(i) / pdel_sum(i)
            ! calculate diagnostic deep convective depth
            conv_depth(i) = pint(i,pver+1) - pmid(i,jctop(i))
         else
            avg_tend_s(i) = 0
            avg_tend_q(i) = 0
            conv_depth(i)  = 0
         end if
      end do

      !----------------------------------------------------------------------------
      ! Note: To conserve total energy we need to account for the kinteic energy tendency
      ! which we can obtain from the velocity tendencies based on the following:
      !   KE_new = (u_new^2 + v_new^2)/2 
      !          = [ (u_old+du)^2 + (v_old+dv)^2 ]/2
      !          = [ ( u_old^2 + 2*u_old*du + du^2 ) + ( v_old^2 + 2*v_old*dv + dv^2 ) ]/2
      !          = ( u_old^2 + v_old^2 )/2 + ( 2*u_old*du + du^2 + 2*v_old*dv + dv^2 )/2
      !          = KE_old + [ 2*u_old*du + du^2 + 2*v_old*dv + dv^2 ] /2

      !----------------------------------------------------------------------------
      ! calculate MCSP tendencies

      do i = 1,ncol

         ! check that ZM produced tendencies over a depth that exceeds the threshold
         if ( conv_depth(i) >= conv_depth_min ) then
            ! check that ZM provided a non-zero column total heating
            if ( avg_tend_s(i) > 0 ) then
               ! check that there is sufficient wind shear to justify coherent organization
               if ( abs(mcsp_shear(i)).ge.shear_min ) then
                  do k = jctop(i),pver

                     ! See eq 7-8 of Moncrieff et al. (2017) - also eq (5) of Moncrieff & Liu (2006)
                     pdepth_mid_k = pint(i,pver+1) - pmid(i,k)
                     pdepth_total = pint(i,pver+1) - pmid(i,jctop(i))

                     ! specify the assumed vertical structure
                     if (do_mcsp_t) mcsp_tend_s(i,k) = -1*heat_coeff * sin(2.0_kind_phys*pi*(pdepth_mid_k/pdepth_total))
                     if (do_mcsp_q) mcsp_tend_q(i,k) = -1*moisture_coeff * sin(2.0_kind_phys*pi*(pdepth_mid_k/pdepth_total))
                     if (do_mcsp_u) mcsp_tend_u(i,k) =    uwind_coeff * (cos(pi*(pdepth_mid_k/pdepth_total)))
                     if (do_mcsp_v) mcsp_tend_v(i,k) =    vwind_coeff * (cos(pi*(pdepth_mid_k/pdepth_total)))

                     ! scale the vertical structure by the ZM heating/drying tendencies
                     if (do_mcsp_t) mcsp_tend_s(i,k) = avg_tend_s(i) * mcsp_tend_s(i,k)
                     if (do_mcsp_t) mcsp_tend_s_max(i) = avg_tend_s(i) * heat_coeff / cpair
                     if (do_mcsp_q) mcsp_tend_q(i,k) = avg_tend_q(i) * mcsp_tend_q(i,k)

                     ! integrate the DSE/qv tendencies for energy/mass fixer
                     if (do_mcsp_t) mcsp_avg_tend_s(i) = mcsp_avg_tend_s(i) + mcsp_tend_s(i,k) * pdel(i,k) / pdel_sum(i)
                     if (do_mcsp_q) mcsp_avg_tend_q(i) = mcsp_avg_tend_q(i) + mcsp_tend_q(i,k) * pdel(i,k) / pdel_sum(i)

                     ! integrate the change in kinetic energy (KE) for energy fixer
                     if (do_mcsp_u.or.do_mcsp_v) then
                        tend_k = ( 2.0_kind_phys*mcsp_tend_u(i,k)*ztodt*u(i,k) + mcsp_tend_u(i,k)*mcsp_tend_u(i,k)*ztodt*ztodt &
                                  +2.0_kind_phys*mcsp_tend_v(i,k)*ztodt*v(i,k) + mcsp_tend_v(i,k)*mcsp_tend_v(i,k)*ztodt*ztodt &
                          )/2.0_kind_phys/ztodt
                        mcsp_avg_tend_k(i) = mcsp_avg_tend_k(i) + tend_k*pdel(i,k) / pdel_sum(i)
                     end if

                  end do ! k = jctop(i),pver
               end if ! shear threshold
            end if ! zm_avg_tend_s(i) > 0
         end if ! convm_depth(i) >= conv_depth_min
      end do

   end if

   !----------------------------------------------------------------------------
   ! calculate final output tendencies

   mcsp_dt_out(1:ncol,1:pver) = 0
   mcsp_dq_out(1:ncol,1:pver) = 0
   mcsp_du_out(1:ncol,1:pver) = 0
   mcsp_dv_out(1:ncol,1:pver) = 0

   mcsp_freq(1:ncol) = 0

   do i = 1,ncol
      do k = jctop(i),pver

         ! update frequency if MCSP contributes any tendency in the column
         if ( abs(mcsp_tend_s(i,k)).gt.0 .or. abs(mcsp_tend_q(i,k)).gt.0 .or.&
              abs(mcsp_tend_u(i,k)).gt.0 .or. abs(mcsp_tend_v(i,k)).gt.0 ) then
            mcsp_freq(i) = 1
         end if

         ! subtract mass weighted average tendencies for energy/mass conservation
         mcsp_dt_out(i,k) = mcsp_tend_s(i,k) - mcsp_avg_tend_s(i)
         mcsp_dq_out(i,k) = mcsp_tend_q(i,k) - mcsp_avg_tend_q(i)
         mcsp_du_out(i,k) = mcsp_tend_u(i,k)
         mcsp_dv_out(i,k) = mcsp_tend_v(i,k)

         ! make sure kinetic energy correction is added to DSE tendency
         ! to conserve total energy whenever momentum tendencies are calculated
         if (do_mcsp_u.or.do_mcsp_v) then
            mcsp_dt_out(i,k) = mcsp_dt_out(i,k) - mcsp_avg_tend_k(i)
         end if

         ! update output tendencies
         if (do_mcsp_t) ptend_s(i,k) =  mcsp_dt_out(i,k)
         if (do_mcsp_q) ptend_q(i,k) =  mcsp_dq_out(i,k)
         if (do_mcsp_u) ptend_u(i,k) =  mcsp_du_out(i,k)
         if (do_mcsp_v) ptend_v(i,k) =  mcsp_dv_out(i,k)

         ! adjust units for diagnostic outputs
         if (do_mcsp_t) mcsp_dt_out(i,k) = mcsp_dt_out(i,k)/cpair

      end do
   end do

   !----------------------------------------------------------------------------
   return

end subroutine MCSP_run

!===================================================================================================

subroutine mcsp_calculate_shear( ncol, pcols, pver, pmid, u, mcsp_shear)
   !----------------------------------------------------------------------------
   ! Purpose: calculate shear for MCSP
   !----------------------------------------------------------------------------
   use interpolate_data, only: vertinterp
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                               intent(in   ) :: ncol       ! number of atmospheric columns (actual)
   integer,                               intent(in   ) :: pcols      ! number of atmospheric columns (max)
   integer,                               intent(in   ) :: pver       ! number of mid-point vertical levels
   real(kind_phys), dimension(pcols,pver),       intent(in   ) :: pmid ! physics state mid-point pressure
   real(kind_phys), dimension(pcols,pver),       intent(in   ) :: u    ! physics state u momentum
   real(kind_phys), dimension(pcols),            intent(  out) :: mcsp_shear
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i
   real(kind_phys), dimension(pcols) :: storm_u         ! u wind at storm reference level set by MCSP_storm_speed_pref
   real(kind_phys), dimension(pcols) :: storm_u_shear   ! u shear at storm reference level set by MCSP_storm_speed_pref
   !----------------------------------------------------------------------------
   ! Interpolate wind to pressure level specified by MCSP_storm_speed_pref
   call vertinterp( ncol, pcols, pver, pmid, storm_speed_pref, u, storm_u )

   !----------------------------------------------------------------------------
   ! calculate low-level shear
   do i = 1,ncol
      if (pmid(i,pver).gt.storm_speed_pref) then
         storm_u_shear(i) = storm_u(i)-u(i,pver)
      else
         storm_u_shear(i) = -999
      end if
      mcsp_shear(i) = storm_u_shear(i)
   end do

   !----------------------------------------------------------------------------
   return

end subroutine mcsp_calculate_shear

subroutine mcsp_finalize
end subroutine mcsp_finalize

!===================================================================================================

end module MCSP
