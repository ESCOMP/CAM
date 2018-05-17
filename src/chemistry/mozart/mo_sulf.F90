

      module mo_sulf
!---------------------------------------------------------------
!	... Annual cycle for sulfur
!---------------------------------------------------------------

      use shr_kind_mod,     only : r8 => shr_kind_r8

      use cam_abortutils,   only : endrun
      use cam_logfile,      only : iulog
      use tracer_data,      only : trfld,trfile
      use physics_types,    only : physics_state
      use ppgrid,           only : begchunk, endchunk
      use physics_buffer,   only : physics_buffer_desc
      use ppgrid,           only : pcols, pver

      use spmd_utils,       only : masterproc

      implicit none

      private
      public  :: sulf_inti, set_sulf_time, sulf_interp, sulf_readnl

      save

      type(trfld), pointer :: fields(:) => null()
      type(trfile) :: file

      logical :: read_sulf = .false.

      character(len=16)  :: fld_name = 'SULFATE'
      character(len=256) :: filename = 'NONE'
      character(len=256) :: filelist = ' '
      character(len=256) :: datapath = ' '
      character(len=32)  :: datatype = 'CYCLICAL'
      logical            :: rmv_file = .false.
      integer            :: cycle_yr  = 0
      integer            :: fixed_ymd = 0
      integer            :: fixed_tod = 0

      logical :: has_sulf_file = .false.

      contains 

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine sulf_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'sulf_readnl'

   character(len=16)  :: sulf_name
   character(len=256) :: sulf_file
   character(len=256) :: sulf_filelist
   character(len=256) :: sulf_datapath
   character(len=32)  :: sulf_type
   logical            :: sulf_rmfile
   integer            :: sulf_cycle_yr
   integer            :: sulf_fixed_ymd
   integer            :: sulf_fixed_tod

   namelist /sulf_nl/ &
      sulf_name,      &
      sulf_file,      &
      sulf_filelist,  &
      sulf_datapath,  &
      sulf_type,      &
      sulf_rmfile,    &
      sulf_cycle_yr,  &
      sulf_fixed_ymd, &
      sulf_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   sulf_name     = fld_name
   sulf_file     = filename
   sulf_filelist = filelist
   sulf_datapath = datapath
   sulf_type     = datatype
   sulf_rmfile   = rmv_file
   sulf_cycle_yr = cycle_yr
   sulf_fixed_ymd= fixed_ymd
   sulf_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'sulf_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, sulf_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(sulf_name,     len(sulf_name),     mpichar, 0, mpicom)
   call mpibcast(sulf_file,     len(sulf_file),     mpichar, 0, mpicom)
   call mpibcast(sulf_filelist, len(sulf_filelist), mpichar, 0, mpicom)
   call mpibcast(sulf_datapath, len(sulf_datapath), mpichar, 0, mpicom)
   call mpibcast(sulf_type,     len(sulf_type),     mpichar, 0, mpicom)
   call mpibcast(sulf_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(sulf_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(sulf_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(sulf_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   fld_name   = sulf_name
   filename   = sulf_file
   filelist   = sulf_filelist
   datapath   = sulf_datapath
   datatype   = sulf_type
   rmv_file   = sulf_rmfile
   cycle_yr   = sulf_cycle_yr
   fixed_ymd  = sulf_fixed_ymd
   fixed_tod  = sulf_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 .and. filename.ne.'NONE') has_sulf_file = .true.

end subroutine sulf_readnl

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine sulf_inti()
!-----------------------------------------------------------------------
! 	... Open netCDF file containing annual sulfur data.  Initialize
!           arrays with the data to be interpolated to the current time.
!
!           It is assumed that the time coordinate is increasing
!           and represents calendar days; range = [1.,366.).
!-----------------------------------------------------------------------
      use spmd_utils,    only : masterproc
      use mo_chem_utls,  only : get_spc_ndx, get_rxt_ndx
      use interpolate_data, only : lininterp_init, lininterp, lininterp_finish, interp_type
      use tracer_data,   only : trcdata_init
      use cam_history,   only : addfld

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: ndxs(5), so4_ndx

      character(len=8), parameter :: fld_names(1) = (/'SULFATE '/)

      ndxs(1) = get_rxt_ndx( 'usr_N2O5_aer' )
      ndxs(2) = get_rxt_ndx( 'usr_NO3_aer' )
      ndxs(3) = get_rxt_ndx( 'usr_NO2_aer' )
      ndxs(4) = get_rxt_ndx( 'usr_HO2_aer' )
      ndxs(5) = get_rxt_ndx( 'het1' )
      so4_ndx = get_spc_ndx( 'SO4' )
      if ( so4_ndx < 1 ) so4_ndx = get_spc_ndx( 'so4_a1' )

      read_sulf = any( ndxs > 0) .and. (so4_ndx < 0) .and. has_sulf_file

      if ( .not. read_sulf ) return

      allocate(file%in_pbuf(size(fld_names)))
      file%in_pbuf(:) = .false. 
      call trcdata_init( (/ fld_name /), filename, filelist, datapath, fields, file, &
           rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)

      call addfld('SULFATE', (/ 'lev' /), 'I','VMR', 'sulfate data' )

      end subroutine sulf_inti

      subroutine set_sulf_time( pbuf2d, state )
!--------------------------------------------------------------------
!	... Check and set time interpolation indicies
!--------------------------------------------------------------------
      use tracer_data,  only : advance_trcdata

      implicit none

!--------------------------------------------------------------------
!	... Dummy args
!--------------------------------------------------------------------
      type(physics_buffer_desc), pointer :: pbuf2d(:,:)
      type(physics_state), intent(in):: state(begchunk:endchunk)                 

      if ( .not. read_sulf ) return

      call advance_trcdata( fields, file, state, pbuf2d  )

      end subroutine set_sulf_time

      subroutine sulf_interp( ncol, lchnk, ccm_sulf )
!-----------------------------------------------------------------------
! 	... Time interpolate sulfatei to current time
!-----------------------------------------------------------------------
      use cam_history,  only : outfld

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)   :: ncol              ! columns in chunk
      integer, intent(in)   :: lchnk             ! chunk number
      real(r8), intent(out) :: ccm_sulf(:,:)     ! output sulfate

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------

      ccm_sulf(:,:) = 0._r8

      if ( .not. read_sulf ) return

      ccm_sulf(:ncol,:) = fields(1)%data(:ncol,:,lchnk)

      call outfld( 'SULFATE', ccm_sulf(:ncol,:), ncol, lchnk )

      end subroutine sulf_interp

      end module mo_sulf
