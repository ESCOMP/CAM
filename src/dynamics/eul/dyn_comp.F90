module dyn_comp
!-----------------------------------------------------------------------
!
! Eulerian dycore interface module
!
!-----------------------------------------------------------------------

use shr_kind_mod,    only: r8=>shr_kind_r8

use spmd_utils,      only: masterproc, npes, mpicom, mpir8

use physconst,       only: pi
use pmgrid,          only: plon, plat, plev, plevp, plnlv, beglat, endlat
use dyn_grid,        only: ptimelevels

use prognostics,     only: n3, ps, u3, v3, t3, q3, phis, pdeld, dpsm, dpsl, div, vort

use cam_control_mod, only: initial_run, ideal_phys, moist_physics, adiabatic
use phys_control,    only: phys_getopts
use constituents,    only: pcnst, cnst_name, cnst_longname, sflxnam, tendnam, &
                           fixcnam, tottnam, hadvnam, vadvnam, cnst_get_ind,  &
                           cnst_read_iv, qmin
use cam_initfiles,   only: initial_file_get_id, topo_file_get_id, pertlim
use inic_analytic,   only: analytic_ic_active, analytic_ic_set_ic
use cam_history,     only: addfld, add_default, horiz_only

use eul_control_mod, only: dif2, hdif_order, kmnhdn, hdif_coef, divdampn, eps, &
                           kmxhdc, eul_nsplit

use scamMod,         only: single_column, use_camiop, have_u, have_v, &
                           have_cldliq, have_cldice, loniop, latiop, scmlat, scmlon, &
                           qobs,tobs,scm_cambfb_mode

use cam_pio_utils,   only: clean_iodesc_list
use pio,             only: file_desc_t, pio_noerr, pio_inq_varid, pio_get_att, &
                           pio_inq_attlen, pio_inq_dimid, pio_inq_dimlen,      &
                           pio_get_var,var_desc_t, pio_seterrorhandling,       &
                           pio_bcast_error, pio_internal_error, pio_offset_kind

#if (defined SPMD)
use spmd_dyn,        only: spmd_readnl
#endif

use cam_logfile,     only: iulog
use cam_abortutils,  only: endrun

implicit none
private
save

public :: &
   dyn_import_t, &
   dyn_export_t, &
   dyn_readnl,   &
   dyn_register, &
   dyn_init

! these structures are not used in this dycore, but are included
! for interface compatibility.
type dyn_import_t
   integer :: placeholder
end type dyn_import_t

type dyn_export_t
   integer :: placeholder
end type dyn_export_t


real(r8), allocatable :: ps_tmp  (:,:  )
real(r8), allocatable :: phis_tmp(:,:  )
real(r8), allocatable :: q3_tmp  (:,:,:)
real(r8), allocatable :: t3_tmp  (:,:,:)
real(r8), allocatable :: arr3d_a (:,:,:)
real(r8), allocatable :: arr3d_b (:,:,:)

logical readvar            ! inquiry flag:  true => variable exists on netCDF file

!=========================================================================================
CONTAINS
!=========================================================================================

subroutine dyn_readnl(nlfile)

   ! Read dynamics namelist group.
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_real8

   ! args
   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! local vars
   integer :: unitn, ierr

   real(r8) :: eul_dif2_coef     ! del2 horizontal diffusion coeff.
   integer  :: eul_hdif_order    ! Order of horizontal diffusion operator
   integer  :: eul_hdif_kmnhdn   ! Nth order horizontal diffusion operator top level.
   real(r8) :: eul_hdif_coef     ! Nth order horizontal diffusion coefficient.
   real(r8) :: eul_divdampn      ! Number of days to invoke divergence damper
   real(r8) :: eul_tfilt_eps     ! Time filter coefficient. Defaults to 0.06.
   integer  :: eul_kmxhdc        ! Number of levels to apply Courant limiter

   namelist /dyn_eul_inparm/ eul_dif2_coef, eul_hdif_order, eul_hdif_kmnhdn, &
      eul_hdif_coef, eul_divdampn, eul_tfilt_eps, eul_kmxhdc, eul_nsplit

   character(len=*), parameter :: sub = 'dyn_readnl'
   !-----------------------------------------------------------------------------

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'dyn_eul_inparm', status=ierr)
      if (ierr == 0) then
         read(unitn, dyn_eul_inparm, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub//': ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(eul_dif2_coef, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_dif2_coef")

   call mpi_bcast(eul_hdif_order, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_hdif_order")

   call mpi_bcast(eul_hdif_kmnhdn, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_hdif_kmnhdn")

   call mpi_bcast(eul_hdif_coef, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_hdif_coef")

   call mpi_bcast(eul_divdampn, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_divdampn")

   call mpi_bcast(eul_tfilt_eps, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_tfilt_eps")

   call mpi_bcast(eul_kmxhdc, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_kmxhdc")

   call mpi_bcast(eul_nsplit, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_nsplit")

   dif2       = eul_dif2_coef
   hdif_order = eul_hdif_order
   kmnhdn     = eul_hdif_kmnhdn
   hdif_coef  = eul_hdif_coef
   divdampn   = eul_divdampn
   eps        = eul_tfilt_eps
   kmxhdc     = eul_kmxhdc

   ! Write namelist variables to logfile
   if (masterproc) then

      write(iulog,*) 'Eulerian Dycore Parameters:'


      ! Order of diffusion
      if (hdif_order < 2 .or. mod(hdif_order, 2) /= 0) then
         write(iulog,*) sub//': Order of diffusion must be greater than 0 and multiple of 2'
         write(iulog,*) 'hdif_order = ', hdif_order
         call endrun(sub//': ERROR: invalid eul_hdif_order specified')
      end if

      if (divdampn > 0._r8) then
         write(iulog,*) '  Divergence damper for spectral dycore invoked for days 0. to ',divdampn,' of this case'
      elseif (divdampn < 0._r8) then
         call endrun (sub//': divdampn must be non-negative')
      else
         write(iulog,*) '  Divergence damper for spectral dycore NOT invoked'
      endif

      if (kmxhdc >= plev .or. kmxhdc < 0) then
         call endrun (sub//':  ERROR:  KMXHDC must be between 0 and plev-1')
      end if

      write(iulog,9108) eps, hdif_order, kmnhdn, hdif_coef, kmxhdc, eul_nsplit

      if (kmnhdn > 1) then
         write(iulog,9109) dif2
      end if

   end if

#if (defined SPMD)
   call spmd_readnl(nlfile)
#endif

9108 format('   Time filter coefficient (EPS)                 ',f10.3,/,&
            '   Horizontal diffusion order (N)                ',i10/, &
            '   Top layer for Nth order horizontal diffusion  ',i10/, &
            '   Nth order horizontal diffusion coefficient    ',e10.3/, &
            '   Number of levels Courant limiter applied      ',i10/,   &
            '   Dynamics Subcycling                           ',i10)

9109 format('   DEL2 horizontal diffusion applied above Nth order diffusion',/,&
            '   DEL2 Horizontal diffusion coefficient (DIF2)  ',e10.3)


end subroutine dyn_readnl

!=========================================================================================

subroutine dyn_register()
end subroutine dyn_register

!=========================================================================================

subroutine dyn_init(dyn_in, dyn_out)

   use prognostics,          only: initialize_prognostics
   use scanslt,              only: scanslt_alloc

   use scamMod,              only: single_column
#if (defined SPMD)
   use spmd_dyn,             only: spmdbuf
#endif
#if (defined BFB_CAM_SCAM_IOP )
   use history_defaults,     only: initialize_iop_history
#endif

   ! Arguments are not used in this dycore, included for compatibility
   type(dyn_import_t), intent(out) :: dyn_in
   type(dyn_export_t), intent(out) :: dyn_out

   ! Local workspace
   integer :: m
   integer :: ixcldice, ixcldliq ! constituent indices for cloud liquid and ice water.
   logical :: history_amwg       ! output for AMWG diagnostics
   logical :: history_budget     ! output tendencies and state variables for CAM4
                                 ! temperature, water vapor, cloud ice and cloud
                                 ! liquid budgets.
   integer :: history_budget_histfile_num  ! output history file number for budget fields
   !----------------------------------------------------------------------------

   ! Initialize prognostics variables
   call initialize_prognostics
   call scanslt_alloc()

#if (defined SPMD)
   ! Allocate communication buffers for collective communications in realloc
   ! routines and in dp_coupling.  Call must come after phys_grid_init.
   call spmdbuf ()
#endif

   if (initial_run) then

#if (defined BFB_CAM_SCAM_IOP )
      call initialize_iop_history()
#endif
      call read_inidat()
      call clean_iodesc_list()
   end if

   call addfld ('ETADOT',(/ 'ilev' /),'A', '1/s','Vertical (eta) velocity', gridname='gauss_grid')
   call addfld ('U&IC',  (/ 'lev' /), 'I', 'm/s','Zonal wind',              gridname='gauss_grid' )
   call addfld ('V&IC',  (/ 'lev' /), 'I', 'm/s','Meridional wind',         gridname='gauss_grid' )
   call add_default ('U&IC',0, 'I')
   call add_default ('V&IC',0, 'I')

   call addfld ('PS&IC',horiz_only,'I',    'Pa','Surface pressure',         gridname='gauss_grid' )
   call addfld ('T&IC',(/ 'lev' /),'I', 'K','Temperature',                  gridname='gauss_grid' )
   call add_default ('PS&IC',0, 'I')
   call add_default ('T&IC',0, 'I')

   do m = 1, pcnst
      call addfld (trim(cnst_name(m))//'&IC',(/ 'lev' /),'I', 'kg/kg',cnst_longname(m), gridname='gauss_grid' )
      call add_default(trim(cnst_name(m))//'&IC',0, 'I')
      call addfld (hadvnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' horizontal advection tendency',  &
           gridname='gauss_grid')
      call addfld (vadvnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' vertical advection tendency',    &
           gridname='gauss_grid')
      call addfld (tendnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' total tendency',                 &
           gridname='gauss_grid')
      call addfld (tottnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' horz + vert + fixer tendency',   &
           gridname='gauss_grid')
      call addfld (fixcnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' tendency due to slt fixer',      &
           gridname='gauss_grid')
   end do

   call addfld ('DUH     ',(/ 'lev' /),'A', 'K/s     ','U horizontal diffusive heating', gridname='gauss_grid')
   call addfld ('DVH     ',(/ 'lev' /),'A', 'K/s     ','V horizontal diffusive heating', gridname='gauss_grid')
   call addfld ('DTH     ',(/ 'lev' /),'A', 'K/s     ','T horizontal diffusive heating', gridname='gauss_grid')

   call addfld ('ENGYCORR',(/ 'lev' /),'A', 'W/m2    ','Energy correction for over-all conservation', gridname='gauss_grid')
   call addfld ('TFIX    ',horiz_only ,'A', 'K/s     ','T fixer (T equivalent of Energy correction)', gridname='gauss_grid')

   call addfld ('FU      ',(/ 'lev' /),'A', 'm/s2    ','Zonal wind forcing term',         gridname='gauss_grid')
   call addfld ('FV      ',(/ 'lev' /),'A', 'm/s2    ','Meridional wind forcing term',    gridname='gauss_grid')
   call addfld ('UTEND   ',(/ 'lev' /),'A', 'm/s2    ','U tendency',                      gridname='gauss_grid')
   call addfld ('VTEND   ',(/ 'lev' /),'A', 'm/s2    ','V tendency',                      gridname='gauss_grid')
   call addfld ('TTEND   ',(/ 'lev' /),'A', 'K/s     ','T tendency',                      gridname='gauss_grid')
   call addfld ('LPSTEN  ',horiz_only ,'A', 'Pa/s    ','Surface pressure tendency',       gridname='gauss_grid')
   call addfld ('VAT     ',(/ 'lev' /),'A', 'K/s     ','Vertical advective tendency of T',gridname='gauss_grid')
   call addfld ('KTOOP   ',(/ 'lev' /),'A', 'K/s     ','(Kappa*T)*(omega/P)',             gridname='gauss_grid')

   call phys_getopts(history_amwg_out=history_amwg, &
        history_budget_out = history_budget, &
        history_budget_histfile_num_out = history_budget_histfile_num)

   if (history_amwg) then
      call add_default ('DTH     ', 1, ' ')
   end if

   if ( history_budget ) then
      if (.not.adiabatic) then
         call cnst_get_ind('CLDLIQ', ixcldliq)
         call cnst_get_ind('CLDICE', ixcldice)
      end if
      ! The following variables are not defined for single column
      if (.not. single_column) then
         call add_default(hadvnam(       1), history_budget_histfile_num, ' ')
         call add_default(vadvnam(       1), history_budget_histfile_num, ' ')
         if (.not.adiabatic) then
            call add_default(hadvnam(ixcldliq), history_budget_histfile_num, ' ')
            call add_default(hadvnam(ixcldice), history_budget_histfile_num, ' ')
            call add_default(vadvnam(ixcldliq), history_budget_histfile_num, ' ')
            call add_default(vadvnam(ixcldice), history_budget_histfile_num, ' ')
         end if
      end if
      call add_default(fixcnam(       1), history_budget_histfile_num, ' ')
      call add_default(tottnam(       1), history_budget_histfile_num, ' ')
      call add_default(tendnam(       1), history_budget_histfile_num, ' ')
      if (.not.adiabatic) then
         call add_default(fixcnam(ixcldliq), history_budget_histfile_num, ' ')
         call add_default(fixcnam(ixcldice), history_budget_histfile_num, ' ')
         call add_default(tottnam(ixcldliq), history_budget_histfile_num, ' ')
         call add_default(tottnam(ixcldice), history_budget_histfile_num, ' ')
         call add_default(tendnam(ixcldliq), history_budget_histfile_num, ' ')
         call add_default(tendnam(ixcldice), history_budget_histfile_num, ' ')
      end if
      call add_default('TTEND',           history_budget_histfile_num, ' ')
      call add_default('TFIX',            history_budget_histfile_num, ' ')
      call add_default('KTOOP',           history_budget_histfile_num, ' ')
      call add_default('VAT',             history_budget_histfile_num, ' ')
      call add_default('DTH',             history_budget_histfile_num, ' ')
   end if

end subroutine dyn_init

!=========================================================================================
! Private routines
!=========================================================================================

subroutine read_inidat()
   ! Read initial dataset and spectrally truncate as appropriate.
   ! Read and process the fields one at a time to minimize
   ! memory usage.

   use ppgrid,           only: begchunk, endchunk, pcols
   use phys_grid,        only: clat_p, clon_p
   use commap,           only: clat, clon
   use comspe,           only: alp, dalp

   use ncdio_atm,        only: infld
   use cam_pio_utils,    only: cam_pio_get_var
   use dyn_tests_utils,  only: vc_moist_pressure

   use iop,              only: setiopupdate,readiopdata

   ! Local variables

   integer i,c,m,n,lat                     ! indices
   integer ncol
   integer ixcldice, ixcldliq  ! indices into q3 array for cloud liq and cloud ice

   integer :: ierr, pio_errtype
   integer :: lonid, latid
   integer :: mlon, morec      ! lon/lat dimension lengths from IC file

   type(file_desc_t), pointer :: fh_ini, fh_topo

   real(r8), pointer, dimension(:,:,:)   :: convptr_2d
   real(r8), pointer, dimension(:,:,:,:) :: convptr_3d
   real(r8), pointer, dimension(:,:,:,:) :: cldptr
   real(r8), pointer, dimension(:,:    ) :: arr2d_tmp
   real(r8), pointer, dimension(:,:    ) :: arr2d
   character*16 fieldname                  ! field name

   real(r8) :: clat2d(plon,plat),clon2d(plon,plat)

   ! variables for analytic initial conditions
   integer,  allocatable       :: glob_ind(:)
   integer                     :: m_cnst(1)
   real(r8), allocatable       :: q4_tmp(:,:,:,:)

   integer londimid,dimlon,latdimid,dimlat,latvarid,lonvarid
   integer strt(3),cnt(3)
   character(len=3), parameter :: arraydims3(3) = (/ 'lon', 'lev', 'lat' /)
   character(len=3), parameter :: arraydims2(2) = (/ 'lon', 'lat' /)
   type(var_desc_t) :: varid
   real(r8), allocatable :: tmp2d(:,:)

   character(len=*), parameter :: sub='read_inidat'
   !----------------------------------------------------------------------------

   fh_ini  => initial_file_get_id()
   fh_topo => topo_file_get_id()

   allocate ( ps_tmp  (plon,plat     ) )
   allocate ( phis_tmp(plon,plat     ) )
   allocate ( q3_tmp  (plon,plev,plat) )
   allocate ( t3_tmp  (plon,plev,plat) )
   allocate ( arr3d_a (plon,plev,plat) )
   allocate ( arr3d_b (plon,plev,plat) )

   if (analytic_ic_active()) then
      allocate(glob_ind(plon * plat))
      m = 1
      do c = 1, plat
         do i = 1, plon
            ! Create a global column index
            glob_ind(m) = i + (c-1)*plon
            m = m + 1
         end do
      end do
      call analytic_ic_set_ic(vc_moist_pressure, clat(:), clon(:,1),            &
           glob_ind(:), U=arr3d_a, V=arr3d_b, T=t3_tmp, PS=ps_tmp, PHIS=phis_tmp)
      readvar = .false.
      call process_inidat('PS')
      call process_inidat('UV')
      call process_inidat('T')
      call process_inidat('PHIS')
      allocate(q4_tmp(plon,plev,plat,1))
      do m = 1, pcnst
         m_cnst(1) = m
         call analytic_ic_set_ic(vc_moist_pressure, clat(:), clon(:,1),          &
              glob_ind(:), Q=q4_tmp, m_cnst=m_cnst)
         arr3d_a(:,:,:) = q4_tmp(:,:,:,1)
         call process_inidat('CONSTS', m_cnst=m, fh=fh_ini)
      end do
      deallocate(q4_tmp)
      deallocate(glob_ind)
      deallocate ( arr3d_a  )
      deallocate ( arr3d_b  )
   else
      !---------------------
      ! Read required fields
      !---------------------

      call pio_seterrorhandling(fh_ini, PIO_BCAST_ERROR, pio_errtype)

      ierr = pio_inq_dimid(fh_ini, 'lon', lonid)
      ierr = pio_inq_dimid(fh_ini, 'lat', latid)
      ierr = pio_inq_dimlen(fh_ini, lonid, mlon)
      ierr = pio_inq_dimlen(fh_ini, latid, morec)
      if (.not. single_column .and. (mlon /= plon .or. morec /= plat)) then
         write(iulog,*) sub//': ERROR: model parameters do not match initial dataset parameters'
         write(iulog,*)'Model Parameters:    plon = ',plon,' plat = ',plat
         write(iulog,*)'Dataset Parameters:  dlon = ',mlon,' dlat = ',morec
         call endrun(sub//': ERROR: model parameters do not match initial dataset parameters')
      end if

      call pio_seterrorhandling(fh_ini, pio_errtype)
      !-----------
      ! 3-D fields
      !-----------

      fieldname = 'U'
      call cam_pio_get_var(fieldname, fh_ini, arraydims3, arr3d_a, found=readvar)
      if (.not. readvar) then
         call endrun(sub//': ERROR: reading '//trim(fieldname))
      end if

      fieldname = 'V'
      call cam_pio_get_var(fieldname, fh_ini, arraydims3, arr3d_b, found=readvar)
      if (.not. readvar) then
         call endrun(sub//': ERROR: reading '//trim(fieldname))
      end if

      call process_inidat('UV')

      fieldname = 'T'
      call cam_pio_get_var(fieldname, fh_ini, arraydims3, t3_tmp, found=readvar)
      if (.not. readvar) then
         call endrun(sub//': ERROR: reading '//trim(fieldname))
      end if

      call process_inidat('T')

      ! Constituents (read and process one at a time)

      do m = 1,pcnst

         readvar   = .false.
         fieldname = cnst_name(m)
         if (cnst_read_iv(m)) then
            call cam_pio_get_var(fieldname, fh_ini, arraydims3, arr3d_a, found=readvar)
         end if
         call process_inidat('CONSTS', m_cnst=m, fh=fh_ini)

      end do

      deallocate ( arr3d_a  )
      deallocate ( arr3d_b  )

      !-----------
      ! 2-D fields
      !-----------

      fieldname = 'PS'
      call cam_pio_get_var(fieldname, fh_ini, arraydims2, ps_tmp, found=readvar)
      if (.not. readvar) then
         call endrun(sub//': ERROR: reading '//trim(fieldname))
      end if
      call process_inidat('PS')
   end if

   ! PHIS processing.  This code allows an analytic specification of PHIS to be
   ! overridden by one from a specified topo file.
   fieldname = 'PHIS'
   readvar   = .false.
   if (associated(fh_topo)) then
      call cam_pio_get_var(fieldname, fh_topo, arraydims2, phis_tmp, found=readvar)
      if (.not. readvar) then
         call endrun(sub//': ERROR: reading '//trim(fieldname))
      end if
      call process_inidat('PHIS', fh=fh_topo)
   else if (.not. analytic_ic_active()) then
      phis_tmp(:,:) = 0._r8
      call process_inidat('PHIS')
   end if

   if (single_column) then
      ps(:,:,1) = ps_tmp(:,:)
   else
      ! Integrals of mass, moisture and geopotential height
      ! (fix mass of moisture as well)
      call global_int
   end if

   deallocate ( ps_tmp   )
   deallocate ( phis_tmp )

   if (single_column) then
      if ( scm_cambfb_mode ) then

         fieldname = 'CLAT1'
         call infld(fieldname, fh_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
            clat2d, readvar, gridname='physgrid')
         if (.not. readvar) then
            call endrun('CLAT not on iop initial file')
         else
            clat(:) = clat2d(1,:)
            clat_p(:)=clat(:)
         end if

         fieldname = 'CLON1'
         call infld(fieldname, fh_ini, 'lon', 'lat', 1, pcols, begchunk, endchunk, &
            clon2d, readvar, gridname='physgrid')
         if (.not. readvar) then
            call endrun('CLON not on iop initial file')
         else
            clon = clon2d
            clon_p(:)=clon(:,1)
         end if

         ! Get latdeg/londeg from initial file for bfb calculations
         ! needed for dyn_grid to determine bounding area and verticies
         ierr = pio_inq_dimid  (fh_ini, 'lon'  , londimid)
         ierr = pio_inq_dimlen (fh_ini, londimid, dimlon)
         ierr = pio_inq_dimid  (fh_ini, 'lat'  , latdimid)
         ierr = pio_inq_dimlen (fh_ini, latdimid, dimlat)
         strt(:)=1
         cnt(1)=dimlon
         cnt(2)=dimlat
         cnt(3)=1
         allocate(latiop(dimlat))
         allocate(loniop(dimlon))
         allocate(tmp2d(dimlon,dimlat))
         ierr = pio_inq_varid (fh_ini,'CLAT1', varid)
         ierr = pio_get_var(fh_ini,varid,strt,cnt,tmp2d)
         latiop(:)=tmp2d(1,:)
         ierr = pio_inq_varid (fh_ini,'CLON1', varid)
         ierr = pio_get_var(fh_ini,varid,strt,cnt,tmp2d)
         loniop(:)=tmp2d(:,1)
         deallocate(tmp2d)
      else

         ! Using a standard iop - make the default grid size is
         ! 4x4 degree square for mo_drydep deposition.(standard ARM IOP area)
         allocate(latiop(2))
         allocate(loniop(2))
         latiop(1)=(scmlat-2._r8)*pi/180_r8
         latiop(2)=(scmlat+2._r8)*pi/180_r8
         loniop(1)=(mod(scmlon-2.0_r8+360.0_r8,360.0_r8))*pi/180.0_r8
         loniop(2)=(mod(scmlon+2.0_r8+360.0_r8,360.0_r8))*pi/180.0_r8
         call setiopupdate()
         ! readiopdata will set all n1 level prognostics to iop value timestep 0
         call readiopdata(timelevel=1)
         ! set t3, and q3(n1) values from iop on timestep 0
         t3(1,:,1,1) = tobs
         q3(1,:,1,1,1) = qobs
      end if
   end if

   deallocate ( q3_tmp  )
   deallocate ( t3_tmp  )

   if (.not. single_column) then
      deallocate ( alp )
      deallocate ( dalp )
   end if

   call copytimelevels()

end subroutine read_inidat

!=========================================================================================

subroutine process_inidat(fieldname, m_cnst, fh)

! Post-process input fields

   use commap
   use comspe
   use spetru
   use dyn_grid,            only: get_horiz_grid_dim_d
   use const_init,          only: cnst_init_default
   use qneg_module,         only: qneg3

#if ( defined SPMD )
   use spmd_dyn,     only: compute_gsfactors
#endif

   ! arguments
   character(len=*),  intent(in)              :: fieldname ! fields to be processed
   integer,           intent(in),    optional :: m_cnst    ! constituent index
   type(file_desc_t), intent(inout), optional :: fh        ! pio file handle

   !---------------------------Local workspace-----------------------------

   integer i,j,k,n,lat,irow               ! grid and constituent indices
   integer :: nglon, nglat, rndm_seed_sz  ! For pertlim
   integer, allocatable :: rndm_seed(:)   ! For pertlim
   real(r8) pertval                       ! perturbation value
   integer  varid                         ! netCDF variable id
   integer  ret
   integer(pio_offset_kind) :: attlen                   ! netcdf return values
   logical  phis_hires                    ! true => PHIS came from hi res topo
   character*256 text
   character*256 trunits                  ! tracer untis

   real(r8), pointer, dimension(:,:,:) :: q_tmp
   real(r8), pointer, dimension(:,:,:) :: tmp3d_a, tmp3d_b, tmp3d_extend
   real(r8), pointer, dimension(:,:  ) :: tmp2d_a, tmp2d_b

#if ( defined BFB_CAM_SCAM_IOP )
   real(r8), allocatable :: ps_sav(:,:)
   real(r8), allocatable :: u3_sav(:,:,:)
   real(r8), allocatable :: v3_sav(:,:,:)
#endif

#if ( defined SPMD )
   integer :: numperlat                   ! number of values per latitude band
   integer :: numsend(0:npes-1)           ! number of items to be sent
   integer :: numrecv                     ! number of items to be received
   integer :: displs(0:npes-1)            ! displacement array
#endif
   character(len=*), parameter :: sub='process_inidat'
   !----------------------------------------------------------------------------

   select case (fieldname)

   !------------
   ! Process U/V
   !------------

   case ('UV')

      allocate ( tmp3d_a(plon,plev,plat) )
      allocate ( tmp3d_b(plon,plev,plat) )

      ! Spectral truncation

      if (single_column) then
         tmp3d_a(:,:,:) = 0._r8
         tmp3d_b(:,:,:) = 0._r8
      else
#if (( defined BFB_CAM_SCAM_IOP )  && ( ! defined DO_SPETRU ))
         allocate ( u3_sav (plon,plev,plat) )
         allocate ( v3_sav (plon,plev,plat) )
         u3_sav(:plon,:plev,:plat) = arr3d_a(:plon,:plev,:plat)
         v3_sav(:plon,:plev,:plat) = arr3d_b(:plon,:plev,:plat)
         call spetru_uv(u3_sav  ,v3_sav  ,vort=tmp3d_a, div=tmp3d_b)
         deallocate ( u3_sav )
         deallocate ( v3_sav )
#else
         call spetru_uv(arr3d_a ,arr3d_b ,vort=tmp3d_a, div=tmp3d_b)
#endif
      end if

#if ( defined SPMD )
      numperlat = plnlv
      call compute_gsfactors (numperlat, numrecv, numsend, displs)

      call mpiscatterv (arr3d_a ,numsend, displs, mpir8,u3  (:,:,beglat:endlat,1)  ,numrecv, mpir8,0,mpicom)
      call mpiscatterv (arr3d_b ,numsend, displs, mpir8,v3  (:,:,beglat:endlat,1)  ,numrecv, mpir8,0,mpicom)
      call mpiscatterv (tmp3d_a ,numsend, displs, mpir8,vort(:,:,beglat:endlat,1) ,numrecv, mpir8,0,mpicom)
      call mpiscatterv (tmp3d_b ,numsend, displs, mpir8,div (:,:,beglat:endlat,1) ,numrecv, mpir8,0,mpicom)
#else
      u3    (:,:,:,1) = arr3d_a(:plon,:plev,:plat)
      v3    (:,:,:,1) = arr3d_b(:plon,:plev,:plat)
      vort  (:,:,:,1) = tmp3d_a(:,:,:)
      div   (:,:,:,1) = tmp3d_b(:,:,:)
#endif
      deallocate ( tmp3d_a )
      deallocate ( tmp3d_b )

   !----------
   ! Process T
   !----------

   case ('T')

      ! Add random perturbation to temperature if required

      if (pertlim .ne. 0.0_r8) then
         if (masterproc) write(iulog,*) sub//': INFO: Adding random perturbation bounded by +/-', &
            pertlim,' to initial temperature field'

         call get_horiz_grid_dim_d(nglon, nglat)
         call random_seed(size=rndm_seed_sz)
         allocate(rndm_seed(rndm_seed_sz))

         do lat = 1, plat
            do i = 1, plon
               ! seed random_number generator based on global column index
               rndm_seed = i + (lat-1)*nglon
               call random_seed(put=rndm_seed)
               do k = 1, plev
                  call random_number (pertval)
                  pertval = 2._r8*pertlim*(0.5_r8 - pertval)
                  t3_tmp(i,k,lat) = t3_tmp(i,k,lat)*(1._r8 + pertval)
               end do
            end do
         end do
         deallocate(rndm_seed)
      end if

      ! Spectral truncation

      if (.not. single_column) then
#if ( ( ! defined BFB_CAM_SCAM_IOP ) || ( defined DO_SPETRU ) )
         call spetru_3d_scalar(t3_tmp)
#endif
      end if

#if ( defined SPMD )
      numperlat = plnlv
      call compute_gsfactors (numperlat, numrecv, numsend, displs)
      call mpiscatterv (t3_tmp  ,numsend, displs, mpir8,t3(:,:,beglat:endlat,1) ,numrecv, mpir8,0,mpicom)
#else
      t3    (:,:,:,1) = t3_tmp(:plon,:plev,:plat)
#endif

   !---------------------
   ! Process Constituents
   !---------------------

   case ('CONSTS')

      if (.not. present(m_cnst)) then
         call endrun(sub//': ERROR:  m_cnst needs to be present in the'// &
            ' argument list')
      end if

      allocate(tmp3d_extend(plon,plev,beglat:endlat))

      if (readvar) then
         ! Check that all tracer units are in mass mixing ratios
         ret = pio_inq_varid(fh, cnst_name(m_cnst), varid)
         ret = pio_get_att(fh, varid, 'units', trunits)
         if (trunits(1:5) .ne. 'KG/KG' .and. trunits(1:5) .ne. 'kg/kg') then
            call endrun(sub//': ERROR:  Units for tracer ' &
               //trim(cnst_name(m_cnst))//' must be in KG/KG')
         end if

      else if (.not. analytic_ic_active()) then

         ! Constituents not read from initial file are initialized by the
         ! package that implements them.  Note that the analytic IC code calls
         ! cnst_init_default internally

         if (m_cnst == 1 .and. moist_physics) then
            call endrun(sub//': ERROR:  Q must be on Initial File')
         end if

         call cnst_init_default(m_cnst, clat, clon(:,1), arr3d_a)
      end if

!$omp parallel do private(lat)
      do lat = 1,plat
         call qneg3(sub, lat, plon, plon, plev    , &
            m_cnst, m_cnst, qmin(m_cnst) ,arr3d_a(1,1,lat))
      end do

      ! if "Q", "CLDLIQ", or "CLDICE", save off for later use
      if (m_cnst == 1) q3_tmp(:plon,:,:) = arr3d_a(:plon,:,:)

#if ( defined SPMD )
      numperlat = plnlv
      call compute_gsfactors(numperlat, numrecv, numsend, displs)
      call mpiscatterv(arr3d_a, numsend, displs, mpir8, tmp3d_extend, numrecv, mpir8,0,mpicom)
      q3(:,:,m_cnst,:,1) = tmp3d_extend(:,:,beglat:endlat)
#else
      q3(:,:plev,m_cnst,:,1) = arr3d_a(:plon,:plev,:plat)
#endif
      deallocate ( tmp3d_extend )

   !-----------
   ! Process PS
   !-----------

   case ('PS')

      allocate ( tmp2d_a(plon,plat) )
      allocate ( tmp2d_b(plon,plat) )

      ! Spectral truncation

      if (single_column) then
         tmp2d_a(:,:) = 0._r8
         tmp2d_b(:,:) = 0._r8
      else
#if (( defined BFB_CAM_SCAM_IOP ) && ( ! defined DO_SPETRU ))
         allocate ( ps_sav(plon,plat) )
         ps_sav(:plon,:plat)=ps_tmp(:plon,:plat)
         call spetru_ps(ps_sav, tmp2d_a, tmp2d_b)
         deallocate ( ps_sav )
#else
         call spetru_ps(ps_tmp, tmp2d_a, tmp2d_b)
#endif
      end if

#if ( defined SPMD )
      numperlat = plon
      call compute_gsfactors (numperlat, numrecv, numsend, displs)
      call mpiscatterv (tmp2d_a ,numsend, displs, mpir8,dpsl ,numrecv, mpir8,0,mpicom)
      call mpiscatterv (tmp2d_b ,numsend, displs, mpir8,dpsm ,numrecv, mpir8,0,mpicom)
#else
      dpsl(:,:) = tmp2d_a(:,:)
      dpsm(:,:) = tmp2d_b(:,:)
#endif
      deallocate ( tmp2d_a )
      deallocate ( tmp2d_b )

      !-------------
      ! Process PHIS
      !-------------

   case ('PHIS')

      ! Check for presence of 'from_hires' attribute to decide whether to filter
      if (readvar) then
         ret = pio_inq_varid (fh, 'PHIS', varid)
         ! Allow pio to return errors in case from_hires doesn't exist
         call pio_seterrorhandling(fh, PIO_BCAST_ERROR)
         ret = pio_inq_attlen (fh, varid, 'from_hires', attlen)
         if (ret.eq.PIO_NOERR .and. attlen.gt.256) then
            call endrun(sub//': ERROR:  from_hires attribute length is too long')
         end if
         ret = pio_get_att(fh, varid, 'from_hires', text)

         if (ret.eq.PIO_NOERR .and. text(1:4).eq.'true') then
            phis_hires = .true.
            if(masterproc) write(iulog,*) sub//': INFO: Will filter input PHIS: attribute from_hires is true'
         else
            phis_hires = .false.
            if(masterproc) write(iulog,*) sub//': INFO: Will not filter input PHIS: attribute ', &
               'from_hires is either false or not present'
         end if
         call pio_seterrorhandling(fh, PIO_INTERNAL_ERROR)

      else
         phis_hires = .false.

      end if

      ! Spectral truncation

      if (.not. single_column) then
#if  (( ! defined BFB_CAM_SCAM_IOP ) || ( defined DO_SPETRU ))
         call spetru_phis(phis_tmp, phis_hires)
#endif
      end if

#if ( defined SPMD )
      numperlat = plon
      call compute_gsfactors (numperlat, numrecv, numsend, displs)
      call mpiscatterv (phis_tmp  ,numsend, displs, mpir8,phis ,numrecv, mpir8,0,mpicom)
#else
      phis = phis_tmp
#endif

   end select

end subroutine process_inidat

!=========================================================================================

subroutine global_int()

   ! Compute global integrals of mass, moisture and geopotential height
   ! and fix mass of atmosphere

   use commap
   use physconst,       only: gravit
#if ( defined SPMD )
   use mpishorthand
   use spmd_dyn,        only: compute_gsfactors
   use spmd_utils,      only: npes
#endif
   use hycoef,          only: hyai, ps0
   use eul_control_mod, only: pdela, qmass1, tmassf, fixmas, &
                              tmass0, zgsint, qmass2, qmassf
   use inic_analytic,   only: analytic_ic_active

   !---------------------------Local workspace-----------------------------

   integer i,k,lat,ihem,irow  ! grid indices
   real(r8) pdelb(plon,plev)  ! pressure diff between interfaces
                              ! using "B" part of hybrid grid only
   real(r8) pssum             ! surface pressure sum
   real(r8) dotproda          ! dot product
   real(r8) dotprodb          ! dot product
   real(r8) zgssum            ! partial sums of phis
   real(r8) hyad (plev)       ! del (A)
   real(r8) tmassf_tmp        ! Global mass integral
   real(r8) qmass1_tmp        ! Partial Global moisture mass integral
   real(r8) qmass2_tmp        ! Partial Global moisture mass integral
   real(r8) qmassf_tmp        ! Global moisture mass integral
   real(r8) zgsint_tmp        ! Geopotential integral

   integer platov2            ! plat/2 or plat (if in scm mode)
#if ( defined SPMD )
   integer :: numperlat         ! number of values per latitude band
   integer :: numsend(0:npes-1) ! number of items to be sent
   integer :: numrecv           ! number of items to be received
   integer :: displs(0:npes-1)  ! displacement array
#endif

   type(file_desc_t), pointer :: fh_topo

   character(len=*), parameter :: sub='global_int'
   !-----------------------------------------------------------------------

   fh_topo => topo_file_get_id()

   if (masterproc) then

      ! Initialize mass and moisture integrals for summation
      ! in a third calculation loop (assures bit-for-bit compare
      ! with non-random history tape).

      tmassf_tmp = 0._r8
      qmass1_tmp = 0._r8
      qmass2_tmp = 0._r8
      zgsint_tmp = 0._r8

      ! Compute pdel from "A" portion of hybrid vertical grid for later use in global integrals
      do k = 1,plev
         hyad(k) = hyai(k+1) - hyai(k)
      end do
      do k = 1,plev
         do i = 1,plon
            pdela(i,k) = hyad(k)*ps0
         end do
      end do

      ! Compute integrals of mass, moisture, and geopotential height
      if (single_column) then
         platov2 = 1
      else
         platov2 = plat/2
      endif
      do irow = 1,platov2
         do ihem = 1,2
            if (ihem.eq.1) then
               lat = irow
            else
               lat = plat - irow + 1
            end if

            ! Accumulate average mass of atmosphere
            call pdelb0 (ps_tmp(1,lat), pdelb, plon)
            pssum  = 0._r8
            do i = 1, plon
               pssum  = pssum  + ps_tmp  (i,lat)
            end do
            tmassf_tmp = tmassf_tmp + w(irow)*pssum/plon

            zgssum = 0._r8
            do i = 1, plon
               zgssum = zgssum + phis_tmp(i,lat)
            end do
            zgsint_tmp = zgsint_tmp + w(irow)*zgssum/plon

            ! Calculate global integrals needed for water vapor adjustment
            do k = 1,plev
               dotproda = 0._r8
               dotprodb = 0._r8
               do i = 1, plon
                  dotproda = dotproda + q3_tmp(i,k,lat)*pdela(i,k)
                  dotprodb = dotprodb + q3_tmp(i,k,lat)*pdelb(i,k)
               end do
               qmass1_tmp = qmass1_tmp + w(irow)*dotproda/plon
               qmass2_tmp = qmass2_tmp + w(irow)*dotprodb/plon
            end do
         end do
      end do                  ! end of latitude loop

      ! Normalize average mass, height
      tmassf_tmp = tmassf_tmp*.5_r8/gravit
      qmass1_tmp = qmass1_tmp*.5_r8/gravit
      qmass2_tmp = qmass2_tmp*.5_r8/gravit
      zgsint_tmp = zgsint_tmp*.5_r8/gravit
      qmassf_tmp = qmass1_tmp + qmass2_tmp

      if (analytic_ic_active()) then
         tmass0 = tmassf_tmp
      else
         ! Globally avgd sfc. partial pressure of dry air (i.e. global dry mass):
         tmass0 = 98222._r8/gravit
         if (.not. associated(fh_topo)) tmass0 = (101325._r8-245._r8)/gravit
         if (adiabatic)                 tmass0 =  tmassf_tmp
         if (ideal_phys )               tmass0 =  100000._r8/gravit
      end if

      if (masterproc) then
         write(iulog,*) sub//': INFO:'
         write(iulog,*) '  Mass of initial data before correction      = ', tmassf_tmp
         write(iulog,*) '  Dry mass will be held at                    = ', tmass0
         write(iulog,*) '  Mass of moisture after removal of negatives = ', qmassf_tmp
         write(iulog,*) '  Globally averaged geopotential height (m)   = ', zgsint_tmp
      end if

      if (analytic_ic_active()) then
         fixmas = 1._r8
      else
         ! Compute and apply an initial mass fix factor which preserves horizontal
         ! gradients of ln(ps).
         if (.not. moist_physics) then
            fixmas = tmass0/tmassf_tmp
         else
            fixmas = (tmass0 + qmass1_tmp)/(tmassf_tmp - qmass2_tmp)
         end if
         ps_tmp = ps_tmp*fixmas
      end if

      ! Global integerals
      tmassf = tmassf_tmp
      qmass1 = qmass1_tmp
      qmass2 = qmass2_tmp
      qmassf = qmassf_tmp
      zgsint = zgsint_tmp

   end if   ! end of if-masterproc

#if ( defined SPMD )
   call mpibcast (tmass0,1,mpir8,0,mpicom)
   call mpibcast (tmassf,1,mpir8,0,mpicom)
   call mpibcast (qmass1,1,mpir8,0,mpicom)
   call mpibcast (qmass2,1,mpir8,0,mpicom)
   call mpibcast (qmassf,1,mpir8,0,mpicom)
   call mpibcast (zgsint,1,mpir8,0,mpicom)

   numperlat = plon
   call compute_gsfactors(numperlat, numrecv, numsend, displs)
   call mpiscatterv(ps_tmp, numsend, displs, mpir8, ps(:,beglat:endlat,1), numrecv, &
                    mpir8, 0, mpicom)
#else
   ps(:,:,1) = ps_tmp(:,:)
#endif

end subroutine global_int

!=========================================================================================

subroutine copytimelevels()

   !---------------------------Local variables-----------------------------

   integer n,i,k,lat            ! index
   real(r8) pdel(plon,plev)     ! pressure arrays needed to calculate
   real(r8) pint(plon,plevp)    !     pdeld
   real(r8) pmid(plon,plev)     !

   ! If dry-type tracers are present, initialize pdeld
   ! First, set current time pressure arrays for model levels etc. to get pdel
   do lat = beglat, endlat
      call plevs0(plon, plon, plev, ps(:,lat,1), pint, pmid, pdel)
      do k = 1, plev
         do i = 1, plon
            pdeld(i,k,lat,1) = pdel(i,k)*(1._r8-q3(i,k,1,lat,1))
         end do
      end do
   end do

   ! Make all time levels of prognostics contain identical data.
   ! Fields to be convectively adjusted only *require* n3 time
   ! level since copy gets done in linems.
   do n = 2, ptimelevels
      ps(:,:,n)     = ps(:,:,1)
      u3(:,:,:,n)   = u3(:,:,:,1)
      v3(:,:,:,n)   = v3(:,:,:,1)
      t3(:,:,:,n)   = t3(:,:,:,1)
      q3(1:plon,:,:,:,n) = q3(1:plon,:,:,:,1)
      vort(:,:,:,n) = vort(:,:,:,1)
      div(:,:,:,n)  = div(:,:,:,1)
      pdeld(1:plon,:,:,n) = pdeld(1:plon,:,:,1)
   end do

end subroutine copytimelevels

!=========================================================================================

end module dyn_comp
