module physconst

   ! Physical constants.  Use csm_share values whenever available.
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use shr_const_mod,  only: shr_const_g
   use shr_const_mod,  only: shr_const_stebol
   use shr_const_mod,  only: shr_const_tkfrz
   use shr_const_mod,  only: shr_const_mwdair
   use shr_const_mod,  only: shr_const_rdair
   use shr_const_mod,  only: shr_const_mwwv
   use shr_const_mod,  only: shr_const_latice
   use shr_const_mod,  only: shr_const_latvap
   use shr_const_mod,  only: shr_const_cpdair
   use shr_const_mod,  only: shr_const_rhofw
   use shr_const_mod,  only: shr_const_cpwv
   use shr_const_mod,  only: shr_const_rgas
   use shr_const_mod,  only: shr_const_karman
   use shr_const_mod,  only: shr_const_pstd
   use shr_const_mod,  only: shr_const_rhodair
   use shr_const_mod,  only: shr_const_avogad
   use shr_const_mod,  only: shr_const_boltz
   use shr_const_mod,  only: shr_const_cpfw
   use shr_const_mod,  only: shr_const_rwv
   use shr_const_mod,  only: shr_const_zvir
   use shr_const_mod,  only: shr_const_pi
   use shr_const_mod,  only: shr_const_rearth
   use shr_const_mod,  only: shr_const_sday
   use shr_const_mod,  only: shr_const_cday
   use shr_const_mod,  only: shr_const_spval
   use shr_const_mod,  only: shr_const_omega
   use shr_const_mod,  only: shr_const_cpvir
   use shr_const_mod,  only: shr_const_tktrip
   use shr_const_mod,  only: shr_const_cpice
   use shr_flux_mod,   only: shr_flux_adjust_constants
   use cam_abortutils, only: endrun
   use constituents,   only: pcnst

   implicit none
   private
   save

   public :: physconst_readnl

   ! Constants based off share code or defined in physconst

   real(r8), public, parameter :: avogad      = shr_const_avogad     ! Avogadro's number (molecules kmole-1)
   real(r8), public, parameter :: boltz       = shr_const_boltz      ! Boltzman's constant (J K-1 molecule-1)
   real(r8), public, parameter :: cday        = shr_const_cday       ! sec in calendar day (seconds)
   real(r8), public, parameter :: cpliq       = shr_const_cpfw       ! specific heat of fresh h2o (J K-1 kg-1)
   real(r8), public, parameter :: cpice       = shr_const_cpice      ! specific heat of ice (J K-1 kg-1)
   real(r8), public, parameter :: karman      = shr_const_karman     ! Von Karman constant
   real(r8), public, parameter :: latice      = shr_const_latice     ! Latent heat of fusion (J kg-1)
   real(r8), public, parameter :: latvap      = shr_const_latvap     ! Latent heat of vaporization (J kg-1)
   real(r8), public, parameter :: pi          = shr_const_pi         ! 3.14...
#ifdef planet_mars
   real(r8), public, parameter :: pstd        = 6.0E1_r8             ! Standard pressure (Pascals)
#else
   real(r8), public, parameter :: pstd        = shr_const_pstd       ! Standard pressure (Pascals)
   real(r8), public, parameter :: tref        = 288._r8              ! Reference temperature (K)
   real(r8), public, parameter :: lapse_rate  = 0.0065_r8            ! reference lapse rate (K m-1)
#endif
   real(r8), public, parameter :: r_universal = shr_const_rgas       ! Universal gas constant (J K-1 kmol-1)
   real(r8), public, parameter :: rhoh2o      = shr_const_rhofw      ! Density of liquid water at STP (kg m-3)
   real(r8), public, parameter :: spval       = shr_const_spval      !special value
   real(r8), public, parameter :: stebol      = shr_const_stebol     ! Stefan-Boltzmann's constant (W m-2 K-4)
   real(r8), public, parameter :: h2otrip     = shr_const_tktrip     ! Triple point temperature of water (K)

   real(r8), public, parameter :: c0          = 2.99792458e8_r8      ! Speed of light in a vacuum (m s-1)
   real(r8), public, parameter :: planck      = 6.6260755e-34_r8     ! Planck's constant (J.s)
   real(r8), public, parameter :: amu         = 1.66053886e-27_r8    ! Atomic Mass Unit (kg)

   ! Molecular weights (g mol-1)
   real(r8), public, parameter :: mwco2       =  44._r8             ! molecular weight co2
   real(r8), public, parameter :: mwn2o       =  44._r8             ! molecular weight n2o
   real(r8), public, parameter :: mwch4       =  16._r8             ! molecular weight ch4
   real(r8), public, parameter :: mwf11       = 136._r8             ! molecular weight cfc11
   real(r8), public, parameter :: mwf12       = 120._r8             ! molecular weight cfc12
   real(r8), public, parameter :: mwo3        =  48._r8             ! molecular weight O3
   real(r8), public, parameter :: mwso2       =  64._r8             ! molecular weight so2
   real(r8), public, parameter :: mwso4       =  96._r8             ! molecular weight so4
   real(r8), public, parameter :: mwh2o2      =  34._r8             ! molecular weight h2o2
   real(r8), public, parameter :: mwdms       =  62._r8             ! molecular weight dms
   real(r8), public, parameter :: mwnh4       =  18._r8             ! molecular wieght nh4
   real(r8), public, protected :: mwh2o       =  shr_const_mwwv     ! molecular weight h2o
   real(r8), public, protected :: mwdry       =  shr_const_mwdair   ! molecular weight dry air

   ! modifiable physical constants for  other planets (including aquaplanet)
   real(r8), public, protected :: gravit  = shr_const_g            ! gravitational acceleration (m s-2)
   real(r8), public, protected :: sday    = shr_const_sday         ! sec in sidereal day (seconds)
   real(r8), public, protected :: cpwv    = shr_const_cpwv         ! specific heat of water vapor (J K-1 kg-1)
   real(r8), public, protected :: cpair   = shr_const_cpdair       ! specific heat of dry air (J K-1 kg-1)
   real(r8), public, protected :: rearth  = shr_const_rearth       ! radius of earth (m)
   real(r8), public, protected :: tmelt   = shr_const_tkfrz        ! Freezing point of water (K)

   !-----  Variables below here are derived from those above -----------------

   real(r8), public, protected :: rga        = 1._r8/shr_const_g         ! reciprocal of gravit (s2 m-1)
   real(r8), public, protected :: ra         = 1._r8/shr_const_rearth    ! reciprocal of earth radius (m-1)
   real(r8), public, protected :: omega      = shr_const_omega           ! earth rot (rad sec-1)
   real(r8), public, protected :: rh2o       = shr_const_rwv             ! Water vapor gas constant (J K-1 kg-1)
   real(r8), public, protected :: rair       = shr_const_rdair           ! Dry air gas constant     (J K-1 kg-1)
   real(r8), public, protected :: epsilo     = shr_const_mwwv/shr_const_mwdair   ! ratio of h2o to dry air molecular weights
   real(r8), public, protected :: zvir       = shr_const_zvir            ! (rh2o/rair) - 1
   real(r8), public, protected :: cpvir      = shr_const_cpvir           ! CPWV/CPDAIR - 1.0
   real(r8), public, protected :: rhodair    = shr_const_rhodair         ! density of dry air at STP (kg m-3)
   real(r8), public, protected :: cappa      = (shr_const_rgas/shr_const_mwdair)/shr_const_cpdair  ! R/Cp
   real(r8), public, protected :: ez                                     ! Coriolis expansion coeff -> omega/sqrt(0.375)
   real(r8), public, protected :: Cpd_on_Cpv = shr_const_cpdair/shr_const_cpwv

!==============================================================================
CONTAINS
!==============================================================================

   ! Read namelist variables.
   subroutine physconst_readnl(nlfile)
      use namelist_utils,  only: find_group_name
      use spmd_utils,      only: masterproc, mpicom, masterprocid
      use spmd_utils,      only: mpi_real8
      use cam_logfile,     only: iulog
      use dyn_tests_utils, only: vc_physics, vc_moist_pressure
      use dyn_tests_utils, only: string_vc, vc_str_lgth

      ! Dummy argument: filepath for file containing namelist input
      character(len=*), intent(in) :: nlfile

      ! Local variables
      integer                     :: unitn, ierr
      logical                     :: newg
      logical                     :: newsday
      logical                     :: newmwh2o
      logical                     :: newcpwv
      logical                     :: newmwdry
      logical                     :: newcpair
      logical                     :: newrearth
      logical                     :: newtmelt
      logical                     :: newomega
      integer,          parameter :: lsize = 76
      integer,          parameter :: fsize = 23
      character(len=*), parameter :: subname = 'physconst_readnl :: '
      character(len=vc_str_lgth)  :: str
      character(len=lsize)        :: banner
      character(len=lsize)        :: bline
      character(len=fsize)        :: field

      ! Physical constants needing to be reset
      !    (e.g., for aqua planet experiments)
      namelist /physconst_nl/  gravit, sday, mwh2o, cpwv, mwdry,              &
           cpair, rearth, tmelt, omega
      !-----------------------------------------------------------------------

      banner = repeat('*', lsize)
      bline = "***"//repeat(' ', lsize - 6)//"***"
2000  format("*** ",a,2("   ",E18.10),"  ***")
      if (masterproc) then
         open(newunit=unitn, file=trim(nlfile), status='old')
         call find_group_name(unitn, 'physconst_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, physconst_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname//'ERROR reading namelist, physconst_nl')
            end if
         end if
         close(unitn)
      end if

      ! Broadcast namelist variables
      call MPI_bcast(gravit, 1, mpi_real8, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: gravit")
      call MPI_bcast(sday,   1, mpi_real8, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: sday")
      call MPI_bcast(mwh2o,  1, mpi_real8, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: mwh20")
      call MPI_bcast(cpwv,   1, mpi_real8, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: cpwv")
      call MPI_bcast(mwdry,  1, mpi_real8, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: mwdry")
      call MPI_bcast(cpair,  1, mpi_real8, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: cpair")
      call MPI_bcast(rearth, 1, mpi_real8, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: rearth")
      call MPI_bcast(tmelt,  1, mpi_real8, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: tmelt")
      call MPI_bcast(omega,  1, mpi_real8, masterprocid, mpicom, ierr)
      if (ierr /= 0) call endrun(subname//": FATAL: mpi_bcast: omega")

      newg     =  gravit /= shr_const_g
      newsday  =  sday   /= shr_const_sday
      newmwh2o =  mwh2o  /= shr_const_mwwv
      newcpwv  =  cpwv   /= shr_const_cpwv
      newmwdry =  mwdry  /= shr_const_mwdair
      newcpair =  cpair  /= shr_const_cpdair
      newrearth=  rearth /= shr_const_rearth
      newtmelt =  tmelt  /= shr_const_tkfrz
      newomega =  omega  /= shr_const_omega

      if (newg .or. newsday .or. newmwh2o .or. newcpwv .or. newmwdry .or.     &
           newrearth .or. newtmelt .or. newomega) then
         if (masterproc) then
            write(iulog, *) banner
            write(iulog, *) '***    New Physical Constant Values set ',       &
                 'via namelist                     ***'
            write(iulog, *) bline
            write(iulog, *) '*** Physical Constant    Old Value                  New Value         ***'
            if (newg) then
               field = 'GRAVIT'
               write(iulog, 2000) field, shr_const_g, gravit
            end if
            if (newsday) then
               field = 'SDAY'
               write(iulog, 2000) field, shr_const_sday, sday
            end if
            if (newmwh2o) then
               field = 'MWH20'
               write(iulog, 2000) field, shr_const_mwwv, mwh2o
            end if
            if (newcpwv) then
               field = 'CPWV'
               write(iulog, 2000) field, shr_const_cpwv, cpwv
            end if
            if (newmwdry) then
               field = 'MWDRY'
               write(iulog, 2000) field, shr_const_mwdair, mwdry
            end if
            if (newcpair) then
               field = 'CPAIR'
               write(iulog, 2000) field, shr_const_cpdair, cpair
            end if
            if (newrearth) then
               field = 'REARTH'
               write(iulog, 2000) field, shr_const_rearth, rearth
            end if
            if (newtmelt) then
               field = 'TMELT'
               write(iulog, 2000) field, shr_const_tkfrz, tmelt
            end if
            if (newomega) then
               field = 'OMEGA'
               write(iulog, 2000) field, shr_const_omega, omega
            end if
            write(iulog,*) banner
         end if
         rga = 1._r8 / gravit
         ra  = 1._r8 / rearth
         if (.not. newomega) then
            omega = 2.0_r8 * pi / sday
         end if
         cpvir  = (cpwv / cpair) - 1._r8
         epsilo = mwh2o / mwdry

         !  defined rair and rh2o before any of the variables that use them
         rair = r_universal / mwdry
         rh2o = r_universal / mwh2o

         cappa       = rair / cpair
         rhodair     = pstd / (rair * tmelt)
         zvir        = (rh2o / rair) - 1.0_r8
         Cpd_on_Cpv  = cpair / cpwv

         ! Adjust constants in shr_flux_mod.
         call shr_flux_adjust_constants(zvir=zvir, cpvir=cpvir, gravit=gravit)
      end if

      ez = omega / sqrt(0.375_r8)
      !
      ! vertical coordinate info
      !
      vc_physics = vc_moist_pressure
      if (masterproc) then
         call string_vc(vc_physics, str)
         write(iulog, *) 'vertical coordinate physics : ', trim(str)
      end if

   end subroutine physconst_readnl

end module physconst
