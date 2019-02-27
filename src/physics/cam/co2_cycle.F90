module co2_cycle

!-------------------------------------------------------------------------------
!
! Purpose:
! Provides distributions of CO2_LND, CO2_OCN, CO2_FF, CO2
! Surface flux from CO2_LND and CO2_OCN can be provided by the flux coupler.
! Surface flux from CO2_FFF and CO2_OCN can be read from a file.
!
! Author: Jeff Lee, Keith Lindsay
!
!-------------------------------------------------------------------------------

   use shr_kind_mod,   only: r8 => shr_kind_r8, cl => shr_kind_cl
   use co2_data_flux,  only: co2_data_flux_type

   implicit none

   private

   ! Public interfaces
   public co2_cycle_readnl              ! read the namelist
   public co2_register                  ! register consituents
   public co2_transport                 ! turn on co2 tracers transport
   public co2_implements_cnst           ! returns true if consituent is implemented by this package
   public co2_init_cnst                 ! initialize mixing ratios if not read from initial file
   public co2_init                      ! initialize (history) variables
   public co2_time_interp_ocn           ! time interpolate co2 flux
   public co2_time_interp_fuel          ! time interpolate co2 flux
   public co2_cycle_set_ptend           ! set tendency from aircraft emissions
   public co2_cycle_set_cnst_type       ! set cnst_type for co2_cycle tracers

   ! Public data
   public data_flux_ocn                 ! data read in for co2 flux from ocn
   public data_flux_fuel                ! data read in for co2 flux from fuel
   public co2_readFlux_aircraft         ! true => read aircraft co2 flux from data file, namelist variable

   type(co2_data_flux_type) :: data_flux_ocn
   type(co2_data_flux_type) :: data_flux_fuel

   public c_i                           ! global index for new constituents
   public co2_readFlux_ocn              ! read ocn co2 flux from data file
   public co2_readFlux_fuel             ! read fuel co2 flux from data file

   ! Namelist variables
   logical :: co2_flag              = .false.         ! true => turn on co2 code, namelist variable
   logical :: co2_readFlux_ocn      = .false.         ! true => read ocn      co2 flux from data file, namelist variable
   logical :: co2_readFlux_fuel     = .false.         ! true => read fuel     co2 flux from data file, namelist variable
   logical :: co2_readFlux_aircraft = .false.         ! true => read aircraft co2 flux from data file, namelist variable
   character(len=cl) :: co2flux_ocn_file  = 'unset' ! co2 flux from ocn
   character(len=cl) :: co2flux_fuel_file = 'unset' ! co2 flux from fossil fuel

   !-------------------------------------------------------------------------------
   ! new constituents
   !-------------------------------------------------------------------------------

   integer, parameter :: ncnst=4                      ! number of constituents implemented

   character(len=7), dimension(ncnst), parameter :: & ! constituent names
      c_names = (/'CO2_OCN', 'CO2_FFF', 'CO2_LND', 'CO2    '/)

   integer :: co2_ocn_glo_ind ! global index of 'CO2_OCN'
   integer :: co2_fff_glo_ind ! global index of 'CO2_FFF'
   integer :: co2_lnd_glo_ind ! global index of 'CO2_LND'
   integer :: co2_glo_ind     ! global index of 'CO2'

   integer, dimension(ncnst) :: c_i                   ! global index

!===============================================================================

contains

!===============================================================================

subroutine co2_cycle_readnl(nlfile)

!-------------------------------------------------------------------------------
! Purpose: Read co2_cycle_nl namelist group.
!-------------------------------------------------------------------------------

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: masterproc
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_logical, mpi_character
   use cam_logfile,     only: iulog
   use cam_cpl_indices, only: index_x2a_Faoo_fco2_ocn
   use cam_abortutils,  only: endrun

   ! Arguments
   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'co2_cycle_readnl'

   namelist /co2_cycle_nl/ co2_flag, co2_readFlux_ocn, co2_readFlux_fuel, co2_readFlux_aircraft, &
                           co2flux_ocn_file, co2flux_fuel_file

   !----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'co2_cycle_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, co2_cycle_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   ! Broadcast namelist variables
   call mpi_bcast(co2_flag,                               1,   mpi_logical,   mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2_flag")
   call mpi_bcast(co2_readFlux_ocn,                       1,   mpi_logical,   mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2_readFlux_ocn")
   call mpi_bcast(co2_readFlux_fuel,                      1,   mpi_logical,   mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2_readFlux_fuel")
   call mpi_bcast(co2_readFlux_aircraft,                  1,   mpi_logical,   mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2_readFlux_aircraft")
   call mpi_bcast(co2flux_ocn_file,   len(co2flux_ocn_file),   mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2flux_ocn_file")
   call mpi_bcast(co2flux_fuel_file, len(co2flux_fuel_file),   mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: co2flux_fuel_file")

   ! Consistency check
   if (co2_readFlux_ocn .and. index_x2a_Faoo_fco2_ocn /= 0) then
      write(iulog,*)'error co2_readFlux_ocn and index_x2a_Faoo_fco2_ocn cannot both be active'
      call endrun(subname // ':: error co2_readFlux_ocn and index_x2a_Faoo_fco2_ocn cannot both be active')
   end if

end subroutine co2_cycle_readnl

!===============================================================================

subroutine co2_register

!-------------------------------------------------------------------------------
! Purpose: register advected constituents
!-------------------------------------------------------------------------------

   use physconst,      only: mwco2, cpair
   use constituents,   only: cnst_add

   ! Local variables
   real(r8), dimension(ncnst) :: &       
      c_mw,    &! molecular weights
      c_cp,    &! heat capacities
      c_qmin    ! minimum mmr

   integer  :: i

   !----------------------------------------------------------------------------

   if (.not. co2_flag) return

   c_mw   = (/     mwco2,     mwco2,     mwco2,     mwco2 /)
   c_cp   = (/     cpair,     cpair,     cpair,     cpair /)
   c_qmin = (/ -1.e36_r8, -1.e36_r8, -1.e36_r8, -1.e36_r8 /) ! disable qneg3

   ! register CO2 constiuents as dry tracers, set indices

   do i = 1, ncnst
      call cnst_add(c_names(i), c_mw(i), c_cp(i), c_qmin(i), c_i(i), longname=c_names(i), mixtype='dry')

      select case (trim(c_names(i)))
      case ('CO2_OCN')
         co2_ocn_glo_ind = c_i(i)
      case ('CO2_FFF')
         co2_fff_glo_ind = c_i(i)
      case ('CO2_LND')
         co2_lnd_glo_ind = c_i(i)
      case ('CO2')
         co2_glo_ind = c_i(i)
      end select
   end do

end subroutine co2_register

!===============================================================================

function co2_transport()

!-------------------------------------------------------------------------------
! Purpose: return true if this package is active
!-------------------------------------------------------------------------------

   ! Return value
   logical :: co2_transport

   !----------------------------------------------------------------------------

   co2_transport = co2_flag

end function co2_transport

!===============================================================================

function co2_implements_cnst(name)

!-------------------------------------------------------------------------------
! Purpose: return true if specified constituent is implemented by this package
!-------------------------------------------------------------------------------

   ! Return value
   logical :: co2_implements_cnst

   ! Arguments
   character(len=*), intent(in) :: name  ! constituent name

   ! Local variables
   integer :: m

   !----------------------------------------------------------------------------

   co2_implements_cnst = .false.

   if (.not. co2_flag) return

   do m = 1, ncnst
      if (name == c_names(m)) then
         co2_implements_cnst = .true.
         return
      end if
   end do

end function co2_implements_cnst

!===============================================================================

subroutine co2_init_cnst(name, latvals, lonvals, mask, q)

!-------------------------------------------------------------------------------
! Purpose:
! Set initial values of CO2_OCN, CO2_FFF, CO2_LND, CO2
! Need to be called from process_inidat in inidat.F90
! (or, initialize co2 in co2_timestep_init)
!-------------------------------------------------------------------------------

   use chem_surfvals,  only: chem_surfvals_get

   ! Arguments
   character(len=*), intent(in)  :: name       ! constituent name
   real(r8),         intent(in)  :: latvals(:) ! lat in degrees (ncol)
   real(r8),         intent(in)  :: lonvals(:) ! lon in degrees (ncol)
   logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
   real(r8),         intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol, plev)

   ! Local variables
   integer :: k

   !----------------------------------------------------------------------------

   if (.not. co2_flag) return

   do k = 1, size(q, 2)
      select case (name)
      case ('CO2_OCN')
         where(mask)
            q(:, k) = chem_surfvals_get('CO2MMR')
         end where
      case ('CO2_FFF')
         where(mask)
            q(:, k) = chem_surfvals_get('CO2MMR')
         end where
      case ('CO2_LND')
         where(mask)
            q(:, k) = chem_surfvals_get('CO2MMR')
         end where
      case ('CO2')
         where(mask)
            q(:, k) = chem_surfvals_get('CO2MMR')
         end where
      end select
   end do

end subroutine co2_init_cnst

!===============================================================================

subroutine co2_init

!-------------------------------------------------------------------------------
! Purpose: initialize co2,
!          declare history variables,
!          read co2 flux form ocn,  as data_flux_ocn
!          read co2 flux form fule, as data_flux_fuel
!-------------------------------------------------------------------------------

   use cam_history,    only: addfld, add_default, horiz_only
   use co2_data_flux,  only: co2_data_flux_init
   use constituents,   only: cnst_name, cnst_longname, sflxnam

   ! Local variables
   integer :: m, mm

   !----------------------------------------------------------------------------

   if (.not. co2_flag) return

    ! Add constituents and fluxes to history file
   do m = 1, ncnst
      mm = c_i(m)

      call addfld(trim(cnst_name(mm))//'_BOT', horiz_only,  'A', 'kg/kg',   trim(cnst_longname(mm))//', Bottom Layer')
      call addfld(cnst_name(mm),               (/ 'lev' /), 'A', 'kg/kg',   cnst_longname(mm))
      call addfld(sflxnam(mm),                 horiz_only,  'A', 'kg/m2/s', trim(cnst_name(mm))//' surface flux')

      call add_default(cnst_name(mm), 1, ' ')
      call add_default(sflxnam(mm),   1, ' ')

      ! The addfld call for the 'TM*' fields are made by default in the
      ! constituent_burden module.
      call add_default('TM'//trim(cnst_name(mm)), 1, ' ')
   end do

   ! Read flux data
   if (co2_readFlux_ocn) then
      call co2_data_flux_init ( co2flux_ocn_file,  'CO2_flux', data_flux_ocn )
   end if

   if (co2_readFlux_fuel) then
      call co2_data_flux_init ( co2flux_fuel_file, 'CO2_flux', data_flux_fuel )
   end if

   if (co2_readFlux_aircraft) then
      call addfld('TMac_CO2', horiz_only,'A', 'kg/m2/s', 'vertical integral of aircraft emission ac_CO2')
   end if

end subroutine co2_init

!===============================================================================

subroutine co2_time_interp_ocn

!-------------------------------------------------------------------------------
! Purpose: Time interpolate co2 flux to current time.
!          Read in new monthly data if necessary
!-------------------------------------------------------------------------------

   use time_manager,   only: is_first_step
   use co2_data_flux,  only: co2_data_flux_advance

   !----------------------------------------------------------------------------

   if (.not. co2_flag) return

   if (co2_readFlux_ocn)  then
      call co2_data_flux_advance ( data_flux_ocn )
   endif

end subroutine co2_time_interp_ocn

!===============================================================================

subroutine co2_time_interp_fuel

!-------------------------------------------------------------------------------
! Purpose: Time interpolate co2 flux to current time.
!          Read in new monthly data if necessary
!-------------------------------------------------------------------------------

   use time_manager,   only: is_first_step
   use co2_data_flux,  only: co2_data_flux_advance

   !----------------------------------------------------------------------------

   if (.not. co2_flag) return

   if (co2_readFlux_fuel) then
      call co2_data_flux_advance ( data_flux_fuel )
   endif

end subroutine co2_time_interp_fuel

!===============================================================================

subroutine co2_cycle_set_ptend(state, pbuf, ptend)

!-------------------------------------------------------------------------------
! Purpose:
! Set ptend, using aircraft CO2 emissions in ac_CO2 from pbuf
!-------------------------------------------------------------------------------

   use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
   use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
   use constituents,   only: pcnst
   use ppgrid,         only: pver
   use physconst,      only: gravit

   ! Arguments
   type(physics_state), intent(in)    :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_ptend), intent(out)   :: ptend     ! indivdual parameterization tendencies

   ! Local variables
   logical :: lq(pcnst)
   integer :: ifld, ncol, k
   real(r8), pointer :: ac_CO2(:,:)

   !----------------------------------------------------------------------------

   if (.not. co2_flag .or. .not. co2_readFlux_aircraft) then
      call physics_ptend_init(ptend, state%psetcols, 'none')
      return
   end if

   ! aircraft fluxes are added to 'CO2_FFF' and 'CO2' tendencies
   lq(:)               = .false.
   lq(co2_fff_glo_ind) = .true.
   lq(co2_glo_ind)     = .true.

   call physics_ptend_init(ptend, state%psetcols, 'co2_cycle_ac', lq=lq)

   ifld = pbuf_get_index('ac_CO2')
   call pbuf_get_field(pbuf, ifld, ac_CO2)

   ! [ac_CO2] = 'kg m-2 s-1'
   ! [ptend%q] = 'kg kg-1 s-1'
   ncol = state%ncol
   do k = 1, pver
      ptend%q(:ncol,k,co2_fff_glo_ind) = gravit * state%rpdeldry(:ncol,k) * ac_CO2(:ncol,k)
      ptend%q(:ncol,k,co2_glo_ind)     = gravit * state%rpdeldry(:ncol,k) * ac_CO2(:ncol,k)
   end do

end subroutine co2_cycle_set_ptend

!===============================================================================

subroutine co2_cycle_set_cnst_type(cnst_type_array, cnst_type_val)

!-------------------------------------------------------------------------------
! Purpose:
! set cnst_type for co2_cycle tracers
!-------------------------------------------------------------------------------

   ! Arguments
   character(len=*), intent(inout) :: cnst_type_array(:)
   character(len=*), intent(in) :: cnst_type_val

   ! Local variables
   integer :: m

   !----------------------------------------------------------------------------

   if (.not. co2_flag) return

   do m = 1, ncnst
      cnst_type_array(c_i(m)) = cnst_type_val
   end do

end subroutine co2_cycle_set_cnst_type

!===============================================================================

end module co2_cycle
