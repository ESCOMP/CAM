module co2_data_flux

!-------------------------------------------------------------------------------
! utilities for reading and interpolating co2 surface fluxes
!-------------------------------------------------------------------------------

   use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
   use input_data_utils, only: time_coordinate
   use cam_abortutils,   only: endrun

   implicit none
   private
   save

   ! Public interfaces
   public co2_data_flux_type

   public co2_data_flux_init
   public co2_data_flux_advance

!-------------------------------------------------------------------------------

   type :: co2_data_flux_type
      character(len=cl)     :: filename
      character(len=cl)     :: varname
      logical               :: initialized
      type(time_coordinate) :: time_coord
      real(r8), pointer     :: co2bdy(:,:,:)   ! bracketing data     (pcols,begchunk:endchunk,2)
      real(r8), pointer     :: co2flx(:,:)     ! Interpolated output (pcols,begchunk:endchunk)
   end type co2_data_flux_type

   ! dimension names for physics grid (physgrid)
   logical           :: dimnames_set = .false.
   character(len=8)  :: dim1name, dim2name

!===============================================================================
contains
!===============================================================================

subroutine co2_data_flux_init (input_file, varname, xin)

!-------------------------------------------------------------------------------
! Initialize co2_data_flux_type instance
!   including initial read of input and interpolation to the current timestep
!-------------------------------------------------------------------------------

   use ioFileMod,        only: getfil
   use ppgrid,           only: begchunk, endchunk, pcols
   use cam_grid_support, only: cam_grid_id, cam_grid_check
   use cam_grid_support, only: cam_grid_get_dim_names

   ! Arguments
   character(len=*),          intent(in)    :: input_file
   character(len=*),          intent(in)    :: varname
   type(co2_data_flux_type),  intent(inout) :: xin

   ! Local variables
   integer  :: grid_id
   real(r8) :: dtime
   character(len=*), parameter :: subname = 'co2_data_flux_init'
   !----------------------------------------------------------------------------

   if (.not. dimnames_set) then
      grid_id = cam_grid_id('physgrid')
      if (.not. cam_grid_check(grid_id)) then
         call endrun(subname // ': ERROR: no "physgrid" grid')
      endif
      call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
      dimnames_set = .true.
   end if

   call getfil(input_file, xin%filename)
   xin%varname = varname
   xin%initialized = .false.

   dtime = 1.0_r8 - 200.0_r8 / 86400.0_r8
   call xin%time_coord%initialize(input_file, delta_days=dtime)

   allocate( xin%co2bdy(pcols,begchunk:endchunk,2), &
             xin%co2flx(pcols,begchunk:endchunk)    )

   call co2_data_flux_advance(xin)

   xin%initialized = .true.

end subroutine co2_data_flux_init

!===============================================================================

subroutine co2_data_flux_advance (xin)

!-------------------------------------------------------------------------------
! Advance the contents of a co2_data_flux_type instance
!   including reading new data, if necessary
!-------------------------------------------------------------------------------

   use cam_pio_utils,    only: cam_pio_openfile
   use ncdio_atm,        only: infld
   use pio,              only: file_desc_t, pio_nowrite, pio_closefile
   use ppgrid,           only: begchunk, endchunk, pcols

   ! Arguments
   type(co2_data_flux_type),  intent(inout) :: xin

   ! Local variables
   character(len=*), parameter :: subname = 'co2_data_flux_advance'
   logical           :: read_data
   integer           :: indx2_pre_adv
   type(file_desc_t) :: fh_co2_data_flux
   logical           :: found

   !----------------------------------------------------------------------------

   read_data = xin%time_coord%read_more() .or. .not. xin%initialized

   indx2_pre_adv = xin%time_coord%indxs(2)

   call xin%time_coord%advance()

   if ( read_data ) then

      call cam_pio_openfile(fh_co2_data_flux, trim(xin%filename), PIO_NOWRITE)

      ! read time-level 1
      ! skip the read if the needed vals are present in time-level 2
      if (xin%initialized .and. xin%time_coord%indxs(1) == indx2_pre_adv) then
         xin%co2bdy(:,:,1) = xin%co2bdy(:,:,2)
      else
         call infld(trim(xin%varname), fh_co2_data_flux, dim1name, dim2name, &
              1, pcols, begchunk, endchunk, xin%co2bdy(:,:,1), found, &
              gridname='physgrid', timelevel=xin%time_coord%indxs(1))
         if (.not. found) then
            call endrun(subname // ': ERROR: ' // trim(xin%varname) // ' not found')
         endif
      endif

      ! read time-level 2
      call infld(trim(xin%varname), fh_co2_data_flux, dim1name, dim2name, &
           1, pcols, begchunk, endchunk, xin%co2bdy(:,:,2), found, &
           gridname='physgrid', timelevel=xin%time_coord%indxs(2))
      if (.not. found) then
         call endrun(subname // ': ERROR: ' // trim(xin%varname) // ' not found')
      endif

      call pio_closefile(fh_co2_data_flux)
   endif

   ! interpolate between time-levels
   ! If time:bounds is in the dataset, and the dataset calendar is compatible with CAM's,
   ! then the time_coordinate class will produce time_coord%wghts(2) == 0.0,
   ! generating fluxes that are piecewise constant in time.

   if (xin%time_coord%wghts(2) == 0.0_r8) then
      xin%co2flx(:,:) = xin%co2bdy(:,:,1)
   else
      xin%co2flx(:,:) = xin%co2bdy(:,:,1) + &
           xin%time_coord%wghts(2) * (xin%co2bdy(:,:,2) - xin%co2bdy(:,:,1))
   endif

end subroutine co2_data_flux_advance

!===============================================================================

end module co2_data_flux
