module zm_conv_mcsp
   !----------------------------------------------------------------------------
   ! Purpose: methods for mesoscale coherent structure parameterization (MCSP)
   !          This scheme essentially redistributes the ZM heating and drying
   !          tendencies vertically in order to mimic the effects of mesoscale
   !          organization. It can also add momentum tendencies, although this
   !          capability has not been extensively tested
   !----------------------------------------------------------------------------
   ! References:
   ! 
   !   Moncrieff, M. W., & Liu, C. (2006). Representing convective organization in
   !     prediction models by a hybrid strategy. J. Atmos. Sci., 63, 3404–3420.
   !     https://doi.org/10.1175/JAS3812.1
   !
   !   Chen, C.-C., Richter, J. H., Liu, C., Moncrieff, M. W., Tang, Q., Lin, W.,
   !     et al. (2021). Effects of organized convection parameterization on the MJO
   !     and precipitation in E3SMv1. Part I: Mesoscale heating. J. Adv.
   !     Mod. Earth Sys., 13, e2020MS002401, https://doi.org/10.1029/2020MS002401
   !
   !   Moncrieff, M. W., C. Liu, and P. Bogenschutz, 2017: Simulation, Modeling, 
   !     and Dynamically Based Parameterization of Organized Tropical Convection 
   !     for Global Climate Models. J. Atmos. Sci., 74, 1363–1380, 
   !     https://doi.org/10.1175/JAS-D-16-0166.1.
   !
   !   Moncrieff, M. W. (2019). Toward a Dynamical Foundation for Organized Convection
   !     Parameterization in GCMs. Geophys. Res. Lett., 46, 14103–14108.
   !     https://doi.org/10.1029/2019GL085316
   !
   !----------------------------------------------------------------------------

   use shr_kind_mod,     only: r8=>shr_kind_r8
   use cam_abortutils,   only: endrun
   use cam_logfile,      only: iulog
   use physconst,        only: cpair, pi


   implicit none
   private

   public :: zm_conv_mcsp_init ! Initialize MCSP output fields
   public :: zm_conv_mcsp_readnl ! Read MCSP namelist
   public :: zm_conv_mcsp_tend ! Perform MCSP tendency calculations
   public :: zm_conv_mcsp_hist ! Write diagnostic quantities to history files

   real(r8), parameter :: MCSP_storm_speed_pref = 600e2_r8 ! pressure level for winds in MCSP calculation [Pa]
   real(r8), parameter :: MCSP_conv_depth_min   = 500e2_r8 ! pressure thickness of convective heating [Pa]
   real(r8), parameter :: mcsp_shear_min        = 3.0_r8   ! min shear value for MCSP to be active
   real(r8), parameter :: mcsp_shear_max        = 200.0_r8 ! max shear value for MCSP to be active

   real(r8) :: zmconv_MCSP_heat_coeff      
   real(r8) :: zmconv_MCSP_moisture_coeff  
   real(r8) :: zmconv_MCSP_uwind_coeff      
   real(r8) :: zmconv_MCSP_vwind_coeff  

!===================================================================================================
contains
!===================================================================================================

subroutine zm_conv_mcsp_init()
   !----------------------------------------------------------------------------
   ! Purpose: initialize MCSP output fields
   !----------------------------------------------------------------------------
   use cam_history,     only: addfld, horiz_only
   use mpishorthand
   !----------------------------------------------------------------------------
   !++
   call addfld('MCSP_DT_max', horiz_only, 'A', 'K/s', 'MCSP max T tendency')
   !--
   call addfld('MCSP_DT', (/'lev'/), 'A', 'K/s',      'MCSP T tendency')
   call addfld('MCSP_DQ', (/'lev'/), 'A', 'kg/kg/s',  'MCSP qv tendency')
   call addfld('MCSP_DU', (/'lev'/), 'A', 'm/s/day',  'MCSP U wind tendency')
   call addfld('MCSP_DV', (/'lev'/), 'A', 'm/s/day',  'MCSP V wind tendency')

   call addfld('MCSP_freq',     horiz_only, 'A', '1',        'MCSP frequency of activation')
   call addfld('MCSP_shear',    horiz_only, 'A', 'm/s',      'MCSP vertical shear of zonal wind')
   call addfld('MCSP_zm_depth', horiz_only, 'A', 'Pa',  'ZM convection depth for MCSP')

end subroutine zm_conv_mcsp_init


!=========================================================================================

subroutine zm_conv_mcsp_readnl(nlfile)

   use spmd_utils,      only: mpicom, masterproc, masterprocid, mpi_real8, mpi_integer, mpi_logical
   use namelist_utils,  only: find_group_name

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'zm_conv_mcsp_readnl'

   namelist /zmconv_mcsp_nl/ zmconv_MCSP_heat_coeff, zmconv_MCSP_moisture_coeff, & 
                        zmconv_MCSP_uwind_coeff, zmconv_MCSP_vwind_coeff
   !-----------------------------------------------------------------------------

   if (masterproc) then
      open( newunit=unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'zmconv_mcsp_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, zmconv_mcsp_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)

   end if

   ! Broadcast namelist variables
   call mpi_bcast(zmconv_MCSP_heat_coeff,                1, mpi_integer, masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_mcsp_readnl: FATAL: mpi_bcast: zmconv_MCSP_heat_coeff")
   call mpi_bcast(zmconv_MCSP_moisture_coeff,            1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_mcsp_readnl: FATAL: mpi_bcast: zmconv_MCSP_moisture_coeff")
   call mpi_bcast(zmconv_MCSP_uwind_coeff,                1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_mcsp_readnl: FATAL: mpi_bcast: zmconv_MCSP_uwind_coeff")
   call mpi_bcast(zmconv_MCSP_vwind_coeff,                1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("zm_conv_mcsp_readnl: FATAL: mpi_bcast: zmconv_MCSP_vwind_coeff")

end subroutine zm_conv_mcsp_readnl

!===================================================================================================

subroutine zm_conv_mcsp_calculate_shear( pcols, ncol, pver, state_pmid, state_u, state_v, mcsp_shear)
   !----------------------------------------------------------------------------
   ! Purpose: calculate shear for MCSP
   !----------------------------------------------------------------------------
   use interpolate_data, only: vertinterp
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                               intent(in   ) :: pcols      ! number of atmospheric columns (max)
   integer,                               intent(in   ) :: ncol       ! number of atmospheric columns (actual)
   integer,                               intent(in   ) :: pver       ! number of mid-point vertical levels
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_pmid ! physics state mid-point pressure
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_u    ! physics state u momentum
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_v    ! physics state v momentum
   real(r8), dimension(pcols),            intent(  out) :: mcsp_shear
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i
   real(r8), dimension(pcols) :: storm_u         ! u wind at storm reference level set by MCSP_storm_speed_pref
   real(r8), dimension(pcols) :: storm_v         ! v wind at storm reference level set by MCSP_storm_speed_pref
   real(r8), dimension(pcols) :: storm_u_shear   ! u shear at storm reference level set by MCSP_storm_speed_pref
   real(r8), dimension(pcols) :: storm_v_shear   ! v shear at storm reference level set by MCSP_storm_speed_pref
   !----------------------------------------------------------------------------
   ! Interpolate wind to pressure level specified by MCSP_storm_speed_pref
   call vertinterp( ncol, pcols, pver, state_pmid, MCSP_storm_speed_pref, state_u, storm_u )
   call vertinterp( ncol, pcols, pver, state_pmid, MCSP_storm_speed_pref, state_v, storm_v )

   !----------------------------------------------------------------------------
   ! calculate low-level shear
   do i = 1,ncol
      if (state_pmid(i,pver).gt.MCSP_storm_speed_pref) then
         storm_u_shear(i) = storm_u(i)-state_u(i,pver)
         storm_v_shear(i) = storm_v(i)-state_v(i,pver)
      else
         storm_u_shear(i) = -999
         storm_v_shear(i) = -999
      end if
      mcsp_shear(i) = storm_u_shear(i)
   end do

   !----------------------------------------------------------------------------
   return

end subroutine zm_conv_mcsp_calculate_shear

!===================================================================================================

subroutine zm_conv_mcsp_tend( pcols, ncol, pver, pverp, &
                              ztodt, jctop,             &
                              state_pmid, state_pint, state_pdel, &
                              state_s, state_q, state_u, state_v, &
                              ptend_zm_s, ptend_zm_q, &
                              ptend_s, ptend_q, ptend_u, ptend_v, &
                              mcsp_dt_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out, &
                              mcsp_freq, mcsp_shear, zm_depth, mcsp_tend_s_max )
   !----------------------------------------------------------------------------
   ! Purpose: perform MCSP tendency calculations
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                               intent(in   ) :: pcols      ! number of atmospheric columns (max)
   integer,                               intent(in   ) :: ncol       ! number of atmospheric columns (actual)
   integer,                               intent(in   ) :: pver       ! number of mid-point vertical levels
   integer,                               intent(in   ) :: pverp      ! number of interface vertical levels
   real(r8),                              intent(in   ) :: ztodt      ! 2x physics time step
   integer,  dimension(pcols),            intent(in   ) :: jctop      ! cloud top level indices
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_pmid ! physics state mid-point pressure
   real(r8), dimension(pcols,pverp),      intent(in   ) :: state_pint ! physics state interface pressure
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_pdel ! physics state pressure thickness
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_s    ! physics state dry static energy
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_q    ! physics state specific humidity
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_u    ! physics state u momentum
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_v    ! physics state v momentum
   real(r8), dimension(pcols,pver),       intent(in   ) :: ptend_zm_s ! input ZM tendency for dry static energy (DSE)
   real(r8), dimension(pcols,pver),       intent(in   ) :: ptend_zm_q ! input ZM tendency for specific humidity (qv)
   real(r8), dimension(pcols,pver),       intent(inout) :: ptend_s    ! output tendency of DSE
   real(r8), dimension(pcols,pver),       intent(inout) :: ptend_q    ! output tendency of qv
   real(r8), dimension(pcols,pver),       intent(inout) :: ptend_u    ! output tendency of u-wind
   real(r8), dimension(pcols,pver),       intent(inout) :: ptend_v    ! output tendency of v-wind
   real(r8), dimension(pcols,pver),       intent(  out) :: mcsp_dt_out! final MCSP tendency for DSE
   real(r8), dimension(pcols,pver),       intent(  out) :: mcsp_dq_out! final MCSP tendency for qv
   real(r8), dimension(pcols,pver),       intent(  out) :: mcsp_du_out! final MCSP tendency for u wind
   real(r8), dimension(pcols,pver),       intent(  out) :: mcsp_dv_out! final MCSP tendency for v wind
   real(r8), dimension(pcols),            intent(  out) :: mcsp_freq  ! MSCP frequency for output
   real(r8), dimension(pcols),            intent(  out) :: mcsp_shear ! shear used to check against threshold
   real(r8), dimension(pcols),            intent(  out) :: zm_depth   ! pressure depth of ZM heating
   !++
   real(r8), dimension(pcols),            intent(  out) :: mcsp_tend_s_max ! max MCSP heating tendency
   !--
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i, k

   real(r8) :: tend_k       ! temporary variable for kinetic energy tendency
   real(r8) :: pdepth_mid_k ! temporary pressure depth used for vertical structure
   real(r8) :: pdepth_total ! temporary pressure depth used for vertical structure

   real(r8), dimension(pcols)      :: zm_avg_tend_s   ! mass weighted column average DSE tendency from ZM
   real(r8), dimension(pcols)      :: zm_avg_tend_q   ! mass weighted column average qv tendency from ZM
   real(r8), dimension(pcols)      :: pdel_sum        ! column integrated pressure thickness

   real(r8), dimension(pcols,pver) :: mcsp_tend_s     ! MCSP tendency before energy fixer for DSE
   real(r8), dimension(pcols,pver) :: mcsp_tend_q     ! MCSP tendency before energy fixer for qv
   real(r8), dimension(pcols,pver) :: mcsp_tend_u     ! MCSP tendency before energy fixer for u wind
   real(r8), dimension(pcols,pver) :: mcsp_tend_v     ! MCSP tendency before energy fixer for v wind

   real(r8), dimension(pcols)      :: mcsp_avg_tend_s ! mass weighted column average MCSP tendency of DSE
   real(r8), dimension(pcols)      :: mcsp_avg_tend_q ! mass weighted column average MCSP tendency of qv 
   real(r8), dimension(pcols)      :: mcsp_avg_tend_k ! mass weighted column average MCSP tendency of kinetic energy

   logical :: do_mcsp_t = .false.   ! internal flag to enable tendency calculations
   logical :: do_mcsp_q = .false.   ! internal flag to enable tendency calculations
   logical :: do_mcsp_u = .false.   ! internal flag to enable tendency calculations
   logical :: do_mcsp_v = .false.   ! internal flag to enable tendency calculations

   !----------------------------------------------------------------------------
   ! initialize variables

   if (zmconv_MCSP_heat_coeff>0) do_mcsp_t = .true.
   if (zmconv_MCSP_moisture_coeff>0) do_mcsp_q = .true.
   if (zmconv_MCSP_uwind_coeff>0) do_mcsp_u = .true.
   if (zmconv_MCSP_vwind_coeff>0) do_mcsp_v = .true.

   zm_avg_tend_s(1:ncol) = 0
   zm_avg_tend_q(1:ncol) = 0

   pdel_sum(1:ncol) = 0

   mcsp_avg_tend_s(1:ncol) = 0
   mcsp_avg_tend_q(1:ncol) = 0
   mcsp_avg_tend_k(1:ncol) = 0
   
   mcsp_tend_s(1:ncol,1:pver) = 0
   mcsp_tend_q(1:ncol,1:pver) = 0
   mcsp_tend_u(1:ncol,1:pver) = 0
   mcsp_tend_v(1:ncol,1:pver) = 0
   !++
   mcsp_tend_s_max = 0._r8
   !--

   if (zmconv_MCSP_heat_coeff>0 .OR. zmconv_MCSP_moisture_coeff>0 .OR. zmconv_MCSP_uwind_coeff>0 .OR. zmconv_MCSP_vwind_coeff>0) then

      !----------------------------------------------------------------------------
      ! calculate shear

      call zm_conv_mcsp_calculate_shear( pcols, ncol, pver, state_pmid, state_u, state_v, mcsp_shear )

      !----------------------------------------------------------------------------
      ! calculate mass weighted column average tendencies from ZM

      zm_depth(1:ncol) = 0
      do i = 1,ncol
         if (jctop(i) .ne. pver) then
            ! integrate pressure and ZM tendencies over column
            do k = jctop(i),pver
               zm_avg_tend_s(i) = zm_avg_tend_s(i) + ptend_zm_s(i,k) * state_pdel(i,k)
               zm_avg_tend_q(i) = zm_avg_tend_q(i) + ptend_zm_q(i,k) * state_pdel(i,k)
               pdel_sum(i) = pdel_sum(i) + state_pdel(i,k)
            end do
            ! normalize integrated ZM tendencies by total mass
            zm_avg_tend_s(i) = zm_avg_tend_s(i) / pdel_sum(i)
            zm_avg_tend_q(i) = zm_avg_tend_q(i) / pdel_sum(i)
            ! calculate diagnostic zm_depth
            zm_depth(i) = state_pint(i,pver+1) - state_pmid(i,jctop(i))
         else
            zm_avg_tend_s(i) = 0
            zm_avg_tend_q(i) = 0
            zm_depth(i)  = 0
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
         if ( zm_depth(i) >= MCSP_conv_depth_min ) then
            ! check that ZM provided a non-zero column total heating
            if ( zm_avg_tend_s(i) > 0 ) then
               ! check that there is sufficient wind shear to justify coherent organization
               if ( abs(mcsp_shear(i)).ge.mcsp_shear_min .and. &
                    abs(mcsp_shear(i)).lt.mcsp_shear_max ) then
                  do k = jctop(i),pver

                     ! See eq 7-8 of Moncrieff et al. (2017) - also eq (5) of Moncrieff & Liu (2006)
                     pdepth_mid_k = state_pint(i,pver+1) - state_pmid(i,k)
                     pdepth_total = state_pint(i,pver+1) - state_pmid(i,jctop(i))

                     ! specify the assumed vertical structure
                     if (do_mcsp_t) mcsp_tend_s(i,k) = -1*zmconv_MCSP_heat_coeff * sin(2.0_r8*pi*(pdepth_mid_k/pdepth_total))
                     if (do_mcsp_q) mcsp_tend_q(i,k) = -1*zmconv_MCSP_moisture_coeff * sin(2.0_r8*pi*(pdepth_mid_k/pdepth_total))
                     if (do_mcsp_u) mcsp_tend_u(i,k) =    zmconv_MCSP_uwind_coeff * (cos(pi*(pdepth_mid_k/pdepth_total)))
                     if (do_mcsp_v) mcsp_tend_v(i,k) =    zmconv_MCSP_vwind_coeff * (cos(pi*(pdepth_mid_k/pdepth_total)))

                     ! scale the vertical structure by the ZM heating/drying tendencies
                     if (do_mcsp_t) mcsp_tend_s(i,k) = zm_avg_tend_s(i) * mcsp_tend_s(i,k)
                     !++
                     if (do_mcsp_t) mcsp_tend_s_max(i) = zm_avg_tend_s(i) * zmconv_MCSP_heat_coeff / cpair
                     !--
                     if (do_mcsp_q) mcsp_tend_q(i,k) = zm_avg_tend_q(i) * mcsp_tend_q(i,k)

                     ! integrate the DSE/qv tendencies for energy/mass fixer
                     if (do_mcsp_t) mcsp_avg_tend_s(i) = mcsp_avg_tend_s(i) + mcsp_tend_s(i,k) * state_pdel(i,k) / pdel_sum(i)
                     if (do_mcsp_q) mcsp_avg_tend_q(i) = mcsp_avg_tend_q(i) + mcsp_tend_q(i,k) * state_pdel(i,k) / pdel_sum(i)

                     ! integrate the change in kinetic energy (KE) for energy fixer
                     if (do_mcsp_u.or.do_mcsp_v) then
                        tend_k = ( 2.0_r8*mcsp_tend_u(i,k)*ztodt*state_u(i,k) + mcsp_tend_u(i,k)*mcsp_tend_u(i,k)*ztodt*ztodt &
                                  +2.0_r8*mcsp_tend_v(i,k)*ztodt*state_v(i,k) + mcsp_tend_v(i,k)*mcsp_tend_v(i,k)*ztodt*ztodt )/2.0_r8/ztodt
                        mcsp_avg_tend_k(i) = mcsp_avg_tend_k(i) + tend_k*state_pdel(i,k) / pdel_sum(i)
                     end if

                  end do ! k = jctop(i),pver
               end if ! shear threshold
            end if ! zm_avg_tend_s(i) > 0
         end if ! zm_depth(i) >= MCSP_conv_depth_min
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
         if (do_mcsp_t) ptend_s(i,k) = ptend_s(i,k) + mcsp_dt_out(i,k)
         if (do_mcsp_q) ptend_q(i,k) = ptend_q(i,k) + mcsp_dq_out(i,k)
         if (do_mcsp_u) ptend_u(i,k) = ptend_u(i,k) + mcsp_du_out(i,k)
         if (do_mcsp_v) ptend_v(i,k) = ptend_v(i,k) + mcsp_dv_out(i,k)

         ! adjust units for diagnostic outputs
         if (do_mcsp_t) mcsp_dt_out(i,k) = mcsp_dt_out(i,k)/cpair

      end do
   end do

   !----------------------------------------------------------------------------
   return

end subroutine zm_conv_mcsp_tend

!===================================================================================================

subroutine zm_conv_mcsp_hist( lchnk, pcols, pver, &
                              mcsp_dt_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out, &
                              mcsp_freq, mcsp_shear, zm_depth, mcsp_tend_s_max )
   !----------------------------------------------------------------------------
   ! Purpose: write diagnostic quantities to history files
   !----------------------------------------------------------------------------
   use cam_history,      only: outfld
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in) :: lchnk           ! chunk identifier
   integer,                         intent(in) :: pcols           ! number of atmospheric columns (max)
   integer,                         intent(in) :: pver            ! number of mid-point vertical levels
   real(r8), dimension(pcols,pver), intent(in) :: mcsp_dt_out     ! final MCSP tendency for DSE
   real(r8), dimension(pcols,pver), intent(in) :: mcsp_dq_out     ! final MCSP tendency for qv
   real(r8), dimension(pcols,pver), intent(in) :: mcsp_du_out     ! final MCSP tendency for u wind
   real(r8), dimension(pcols,pver), intent(in) :: mcsp_dv_out     ! final MCSP tendency for v wind
   real(r8), dimension(pcols),      intent(in) :: mcsp_freq       ! MSCP frequency for output
   real(r8), dimension(pcols),      intent(in) :: mcsp_shear      ! shear used to check against threshold
   real(r8), dimension(pcols),      intent(in) :: zm_depth        ! pressure depth of ZM heating
   !++
   real(r8), dimension(pcols),      intent(in) :: mcsp_tend_s_max ! max MCSP heating tendency
   !--
   !----------------------------------------------------------------------------
   ! write out MCSP diagnostic history fields
   call outfld('MCSP_DT',       mcsp_dt_out, pcols, lchnk )
   call outfld('MCSP_DQ',       mcsp_dq_out, pcols, lchnk )
   call outfld('MCSP_DU',       mcsp_du_out, pcols, lchnk )
   call outfld('MCSP_DV',       mcsp_dv_out, pcols, lchnk )
   call outfld('MCSP_freq',     mcsp_freq,   pcols, lchnk )
   call outfld('MCSP_shear',    mcsp_shear,  pcols, lchnk )
   call outfld('MCSP_zm_depth', zm_depth,    pcols, lchnk )
   !++
   call outfld('MCSP_DT_max', mcsp_tend_s_max, pcols, lchnk)
   !--
   !----------------------------------------------------------------------------
   return

end subroutine zm_conv_mcsp_hist

!===================================================================================================

end module zm_conv_mcsp
