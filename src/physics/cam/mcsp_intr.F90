module mcsp_intr
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


   implicit none
   private

   public :: mcsp_register  ! Register fields wit the output buffer
   public :: mcsp_readnl    ! Read MCSP namelist
   public :: mcsp_intr_init ! initialize MCSP namelist variables
   public :: mcsp_tend      ! Perform MCSP tendency calculations

   real(r8) :: MCSP_storm_speed_pref                       ! pressure level for winds in MCSP calculation [Pa]
   real(r8) :: MCSP_conv_depth_min                         ! pressure thickness of convective heating [Pa]
   real(r8) :: mcsp_shear_min                              ! min shear value for MCSP to be active

   real(r8) :: MCSP_heat_coeff                             ! heating coefficient for MCSP 
   real(r8) :: MCSP_moisture_coeff                         ! moisture coefficient for MCSP  
   real(r8) :: MCSP_uwind_coeff                            ! uwind coefficient for MCSP
   real(r8) :: MCSP_vwind_coeff                            ! vwind coefficient for MCSP

!=========================================================================================
contains
!=========================================================================================
subroutine mcsp_register()
   !----------------------------------------------------------------------------
   ! Purpose: register MCSP output fields
   !----------------------------------------------------------------------------
   use cam_history,     only: addfld, horiz_only
   !----------------------------------------------------------------------------
   call addfld('MCSP_DT_max', horiz_only, 'A', 'K/s', 'MCSP max T tendency')
   call addfld('MCSP_DT', (/'lev'/), 'A', 'K/s',      'MCSP T tendency')
   call addfld('MCSP_DQ', (/'lev'/), 'A', 'kg/kg/s',  'MCSP qv tendency')
   call addfld('MCSP_DU', (/'lev'/), 'A', 'm/s/day',  'MCSP U wind tendency')
   call addfld('MCSP_DV', (/'lev'/), 'A', 'm/s/day',  'MCSP V wind tendency')

   call addfld('MCSP_freq',     horiz_only, 'A', '1',        'MCSP frequency of activation')
   call addfld('MCSP_shear',    horiz_only, 'A', 'm/s',      'MCSP vertical shear of zonal wind')
   call addfld('MCSP_conv_depth', horiz_only, 'A', 'Pa',  'ZM convection depth for MCSP')

end subroutine mcsp_register

!=========================================================================================

subroutine mcsp_readnl(nlfile)

   use spmd_utils,      only: mpicom, masterproc, masterprocid, mpi_real8, mpi_integer, mpi_logical
   use namelist_utils,  only: find_group_name

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'mcsp_readnl'

   namelist /mcsp_nl/ MCSP_heat_coeff, MCSP_moisture_coeff, & 
                      MCSP_uwind_coeff, MCSP_vwind_coeff,   &
                      MCSP_storm_speed_pref, MCSP_conv_depth_min, mcsp_shear_min
   !-----------------------------------------------------------------------------

   if (masterproc) then
      open( newunit=unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'mcsp_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, mcsp_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)

   end if

   ! Broadcast namelist variables
   call mpi_bcast(MCSP_heat_coeff,                1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("mcsp_readnl: FATAL: mpi_bcast: MCSP_heat_coeff")
   call mpi_bcast(MCSP_moisture_coeff,            1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("mcsp_readnl: FATAL: mpi_bcast: MCSP_moisture_coeff")
   call mpi_bcast(MCSP_uwind_coeff,                1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("mcsp_readnl: FATAL: mpi_bcast: MCSP_uwind_coeff")
   call mpi_bcast(MCSP_vwind_coeff,                1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("mcsp_readnl: FATAL: mpi_bcast: MCSP_vwind_coeff")
   call mpi_bcast(MCSP_storm_speed_pref,           1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("mcsp_readnl: FATAL: mpi_bcast: MCSP_storm_speed_pref")
   call mpi_bcast(MCSP_conv_depth_min,             1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("mcsp_readnl: FATAL: mpi_bcast: MCSP_conv_depth_min")
   call mpi_bcast(mcsp_shear_min,                  1, mpi_real8,   masterprocid, mpicom, ierr)
   if (ierr /= 0) call endrun("mcsp_readnl: FATAL: mpi_bcast: mcsp_shear_min")

end subroutine mcsp_readnl

!===================================================================================================

subroutine mcsp_intr_init()
  use MCSP,         only: MCSP_init
  use spmd_utils,   only: masterproc
  use cam_logfile,  only: iulog

  implicit none

  call MCSP_init(MCSP_heat_coeff, MCSP_moisture_coeff, MCSP_uwind_coeff, MCSP_vwind_coeff, &
                 MCSP_storm_speed_pref, MCSP_conv_depth_min, mcsp_shear_min, masterproc, iulog)

end subroutine mcsp_intr_init
        
!===================================================================================================

subroutine mcsp_tend( state, ptend, ztodt, jctop, tend_t, tend_q )

   use MCSP,           only: MCSP_run 
   use physconst,      only: cpair
   use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
   use ppgrid,         only: pver, pverp, pcols
   use physconst,      only: cpair, pi
   use constituents,   only: pcnst

   ! Arguments
   type(physics_state), intent( in) :: state   ! Physics state variables
   type(physics_ptend), intent(out) :: ptend   ! individual parameterization tendencies
   real(r8),                              intent(in   ) :: ztodt      ! 2x physics time step
   integer,  dimension(pcols),            intent(in   ) :: jctop      ! cloud top level indices
   real(r8), dimension(pcols,pver),       intent(in   ) :: tend_t ! input deep convective temperature tendency
   real(r8), dimension(pcols,pver),       intent(in   ) :: tend_q ! input deep convective tendency for specific humidity (qv)

   ! local variables
   integer                                              :: lchnk      ! chunk identifier
   integer                                              :: ncol       ! number of atmospheric columns (actual)
   real(r8), dimension(pcols,pver)                      :: mcsp_dt_out! final MCSP tendency for DSE
   real(r8), dimension(pcols,pver)                      :: mcsp_dq_out! final MCSP tendency for qv
   real(r8), dimension(pcols,pver)                      :: mcsp_du_out! final MCSP tendency for u wind
   real(r8), dimension(pcols,pver)                      :: mcsp_dv_out! final MCSP tendency for v wind
   real(r8), dimension(pcols)                           :: mcsp_freq  ! MSCP frequency for output
   real(r8), dimension(pcols)                           :: mcsp_shear ! shear used to check against threshold
   real(r8), dimension(pcols)                           :: conv_depth   ! pressure depth of deep convection
   real(r8), dimension(pcols)                           :: mcsp_tend_s_max ! max MCSP heating tendency
   logical                                              :: lq(pcnst)
   
   !----------------------------------------------------------------------------
   ! Perform MCSP tendency calculations
   !----------------------------------------------------------------------------
   lchnk = state%lchnk
   ncol = state%ncol
   lq(:) = .false.
   lq(1) = .true.

   call physics_ptend_init(ptend, state%psetcols, 'MCSP', ls = .true., lq = lq, lu = .true., lv = .true.)

   call MCSP_run( pcols, ncol, pver, pverp, cpair, pi,                    &
                              ztodt, jctop,                               &
                              state%pmid, state%pint, state%pdel,         &
                              state%s, state%q(:,:,1), state%u, state%v,  &
                              tend_t, tend_q,                             &
                              ptend%s, ptend%q(:,:,1), ptend%u, ptend%v,                               &
                              mcsp_dt_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out,                      &
                              mcsp_freq, mcsp_shear, conv_depth, mcsp_tend_s_max )

   !----------------------------------------------------------------------------
   ! outpuf MCSP variables
   !----------------------------------------------------------------------------
   call mcsp_hist( lchnk, pcols, pver, &
                              mcsp_dt_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out, &
                              mcsp_freq, mcsp_shear, conv_depth, mcsp_tend_s_max )

end subroutine mcsp_tend

!===================================================================================================

subroutine mcsp_hist( lchnk, pcols, pver, &
                              mcsp_dt_out, mcsp_dq_out, mcsp_du_out, mcsp_dv_out, &
                              mcsp_freq, mcsp_shear, conv_depth, mcsp_tend_s_max )
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
   real(r8), dimension(pcols),      intent(in) :: conv_depth      ! pressure depth of ZM heating
   real(r8), dimension(pcols),      intent(in) :: mcsp_tend_s_max ! max MCSP heating tendency
   !----------------------------------------------------------------------------
   ! write out MCSP diagnostic history fields
   call outfld('MCSP_DT',       mcsp_dt_out, pcols, lchnk )
   call outfld('MCSP_DQ',       mcsp_dq_out, pcols, lchnk )
   call outfld('MCSP_DU',       mcsp_du_out, pcols, lchnk )
   call outfld('MCSP_DV',       mcsp_dv_out, pcols, lchnk )
   call outfld('MCSP_freq',     mcsp_freq,   pcols, lchnk )
   call outfld('MCSP_shear',    mcsp_shear,  pcols, lchnk )
   call outfld('MCSP_conv_depth', conv_depth,    pcols, lchnk )
   call outfld('MCSP_DT_max',   mcsp_tend_s_max, pcols, lchnk)
   !----------------------------------------------------------------------------
   return

end subroutine mcsp_hist

!===================================================================================================

end module mcsp_intr


