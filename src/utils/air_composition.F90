module air_composition

   use shr_kind_mod,   only: r8 => shr_kind_r8
   use cam_abortutils, only: endrun

   implicit none
   private
   save

   public  :: air_composition_readnl
   public  :: air_composition_init

   private :: dry_air_species_info

   integer,  parameter :: unseti = -HUGE(1)
   real(r8), parameter :: unsetr = HUGE(1.0_r8)

   ! composition of air
   !
   integer, parameter :: num_names_max = 20 ! Should match namelist definition
   character(len=6)   :: dry_air_species(num_names_max)
   character(len=6)   :: water_species_in_air(num_names_max)

   integer, protected, public :: dry_air_species_num
   integer, protected, public :: water_species_in_air_num

   ! Thermodynamic variables
   integer,               protected, public :: thermodynamic_active_species_num = unseti
   integer,  allocatable, protected, public :: thermodynamic_active_species_idx(:)
   integer,  allocatable,            public :: thermodynamic_active_species_idx_dycore(:)
   real(r8), allocatable, protected, public :: thermodynamic_active_species_cp(:)
   real(r8), allocatable, protected, public :: thermodynamic_active_species_cv(:)
   real(r8), allocatable, protected, public :: thermodynamic_active_species_R(:)
   ! thermodynamic_active_species_mwi: inverse molecular weights dry air
   real(r8), allocatable, protected, public :: thermodynamic_active_species_mwi(:)
   ! thermodynamic_active_species_kv: molecular diffusion
   real(r8), allocatable, protected, public :: thermodynamic_active_species_kv(:)
   ! thermodynamic_active_species_kc: thermal conductivity
   real(r8), allocatable, protected, public :: thermodynamic_active_species_kc(:)
   !
   ! for energy computations liquid and ice species need to be identified
   !
   ! thermodynamic_active_species_liq_num: number of liquid water species
   integer,               protected, public :: thermodynamic_active_species_liq_num = unseti
   ! thermodynamic_active_species_ice_num: number of frozen water species
   integer,               protected, public :: thermodynamic_active_species_ice_num = unseti
   ! thermodynamic_active_species_liq_idx: index of liquid water species
   integer,  allocatable, protected, public :: thermodynamic_active_species_liq_idx(:)
   ! thermodynamic_active_species_liq_idx_dycore: index of liquid water species
   integer,  allocatable,            public :: thermodynamic_active_species_liq_idx_dycore(:)
   ! thermodynamic_active_species_ice_idx: index of ice water species
   integer,  allocatable, protected, public :: thermodynamic_active_species_ice_idx(:)
   ! thermodynamic_active_species_ice_idx_dycore: index of ice water species
   integer,  allocatable,            public :: thermodynamic_active_species_ice_idx_dycore(:)
   ! enthalpy_reference_state: choices: 'ice', 'liq', 'wv'
   character(len=3), public, protected :: enthalpy_reference_state = 'xxx'

   integer, protected, public :: wv_idx = -1 ! Water vapor index

   !------------- Variables for consistent themodynamics --------------------
   !

   ! standard dry air (constant composition)
   real(r8), public, protected :: mmro2 = unsetr  ! Mass mixing ratio of O2
   real(r8), public, protected :: mmrn2 = unsetr  ! Mass mixing ratio of N2
   real(r8), public, protected :: o2_mwi = unsetr ! Inverse mol. weight of O2
   real(r8), public, protected :: n2_mwi = unsetr ! Inverse mol. weight of N2
   real(r8), public, protected :: mbar = unsetr   ! Mean mass at mid level

CONTAINS

   ! Read namelist variables.
   subroutine air_composition_readnl(nlfile)
      use namelist_utils, only: find_group_name
      use spmd_utils,     only: masterproc, mpicom, masterprocid
      use spmd_utils,     only: mpi_character
      use cam_logfile,    only: iulog

      ! Dummy argument: filepath for file containing namelist input
      character(len=*), intent(in) :: nlfile

      ! Local variables
      integer                     :: unitn, ierr, indx
      integer,          parameter :: lsize = 76
      character(len=*), parameter :: subname = 'air_composition_readnl :: '
      character(len=lsize)        :: banner
      character(len=lsize)        :: bline

      ! Variable components of dry air and water species in air
      namelist /air_composition_nl/ dry_air_species, water_species_in_air
      !-----------------------------------------------------------------------

      banner = repeat('*', lsize)
      bline = "***"//repeat(' ', lsize - 6)//"***"
2000  format("*** ",a,2("   ",E18.10),"  ***")

      ! Read variable components of dry air and water species in air
      dry_air_species = (/ (' ', indx = 1, num_names_max) /)
      water_species_in_air = (/ (' ', indx = 1, num_names_max) /)

      if (masterproc) then
         open(newunit=unitn, file=trim(nlfile), status='old')
         call find_group_name(unitn, 'air_composition_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, air_composition_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname//'ERROR reading namelist, air_composition_nl')
            end if
         end if
         close(unitn)
      end if

      call mpi_bcast(dry_air_species, len(dry_air_species)*num_names_max,     &
           mpi_character, masterprocid, mpicom, ierr)
      call mpi_bcast(water_species_in_air,                                    &
           len(water_species_in_air)*num_names_max, mpi_character,            &
           masterprocid, mpicom, ierr)

      dry_air_species_num = 0
      water_species_in_air_num = 0
      do indx = 1, num_names_max
         if ( (LEN_TRIM(dry_air_species(indx)) > 0) .and.                     &
              (TRIM(dry_air_species(indx)) /= 'N2')) then
            dry_air_species_num = dry_air_species_num + 1
         end if
         if (LEN_TRIM(water_species_in_air(indx)) > 0) then
            water_species_in_air_num = water_species_in_air_num + 1
         end if
      end do

      ! Initialize number of thermodynamically active species
      thermodynamic_active_species_num =                                      &
           dry_air_species_num + water_species_in_air_num

      if (masterproc) then
         write(iulog, *) banner
         write(iulog, *) bline

         if (dry_air_species_num == 0) then
            write(iulog, *) " Thermodynamic properties of dry air are ",      &
                 "fixed at troposphere values"
         else
            write(iulog, *) " Thermodynamic properties of dry air are ",      &
                 "based on variable composition of the following species:"
            do indx = 1, dry_air_species_num
               write(iulog, *) '   ',  trim(dry_air_species(indx))
            end do
            write(iulog,*) ' '
         end if
         write(iulog,*) " Thermodynamic properties of moist air are ",        &
              "based on variable composition of the following water species:"
         do indx = 1, water_species_in_air_num
            write(iulog, *) '   ', trim(water_species_in_air(indx))
         end do
         write(iulog, *) bline
         write(iulog, *) banner
      end if

   end subroutine air_composition_readnl

   !===========================================================================

   subroutine air_composition_init()
      use string_utils, only: int2str
      use spmd_utils,   only: masterproc
      use cam_logfile,  only: iulog
      use physconst,    only: r_universal, cpair, rair, cpwv, rh2o, cpliq, cpice
      use constituents, only: cnst_get_ind, cnst_mw

      integer  :: icnst, ix, isize, ierr, idx
      integer  :: liq_num, ice_num
      integer  :: liq_idx(water_species_in_air_num)
      integer  :: ice_idx(water_species_in_air_num)
      logical  :: has_liq, has_ice
      real(r8) :: mw

      character(len=*), parameter :: subname = 'composition_init'
      character(len=*), parameter :: errstr = subname//": failed to allocate "

      !
      ! define cp and R for species in species_name
      !
      ! Last major species in namelist dry_air_species is derived from the
      !    other major species (since the sum of dry mixing ratios for
      !    major species of dry air add must add to one)
      !
      ! cv = R * dofx / 2;   cp = R * (1 + (dofx / 2))
      ! DOF == Degrees of Freedom
      ! dof1 = monatomic ideal gas, 3 translational DOF
      real(r8), parameter :: dof1 = 3._r8
      real(r8), parameter :: cv1 = 0.5_r8 * r_universal * dof1
      real(r8), parameter :: cp1 = 0.5_r8 * r_universal * (2._r8 + dof1)
      ! dof2 = diatomic ideal gas, 3 translational + 2 rotational = 5 DOF
      real(r8), parameter :: dof2 = 5._r8
      real(r8), parameter :: cv2 = 0.5_r8 * r_universal * dof2
      real(r8), parameter :: cp2 = 0.5_r8 * r_universal * (2._r8 + dof2)
      ! dof3 = polyatomic ideal gas, 3 translational + 3 rotational = 6 DOF
      real(r8), parameter :: dof3 = 6._r8
      real(r8), parameter :: cv3 = 0.5_r8 * r_universal * dof3
      real(r8), parameter :: cp3 = 0.5_r8 * r_universal * (2._r8 + dof3)

      liq_num = 0
      ice_num = 0
      has_liq = .false.
      has_ice = .false.
      ! standard dry air (constant composition)
      o2_mwi = 1._r8 / 32._r8
      n2_mwi = 1._r8 / 28._r8
      mmro2 = 0.235_r8
      mmrn2 = 0.765_r8
      mbar = 1._r8 / ((mmro2 * o2_mwi) + (mmrn2 * n2_mwi))

      ! init for variable composition dry air

      isize = dry_air_species_num + water_species_in_air_num
      allocate(thermodynamic_active_species_idx(isize), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"thermodynamic_active_species_idx")
      end if
      allocate(thermodynamic_active_species_idx_dycore(isize), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"thermodynamic_active_species_idx_dycore")
      end if
      allocate(thermodynamic_active_species_cp(0:isize), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"thermodynamic_active_species_cp")
      end if
      allocate(thermodynamic_active_species_cv(0:isize), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"thermodynamic_active_species_cv")
      end if
      allocate(thermodynamic_active_species_R(0:isize), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"thermodynamic_active_species_R")
      end if

      isize = dry_air_species_num
      allocate(thermodynamic_active_species_mwi(0:isize), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"thermodynamic_active_species_mwi")
      end if
      allocate(thermodynamic_active_species_kv(0:isize), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"thermodynamic_active_species_kv")
      end if
      allocate(thermodynamic_active_species_kc(0:isize), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"thermodynamic_active_species_kc")
      end if

      thermodynamic_active_species_idx        = -HUGE(1)
      thermodynamic_active_species_idx_dycore = -HUGE(1)
      thermodynamic_active_species_cp         = 0.0_r8
      thermodynamic_active_species_cv         = 0.0_r8
      thermodynamic_active_species_R          = 0.0_r8
      thermodynamic_active_species_mwi        = 0.0_r8
      thermodynamic_active_species_kv         = 0.0_r8
      thermodynamic_active_species_kc         = 0.0_r8
      !
      if (dry_air_species_num > 0) then
         !
         ! The last major species in dry_air_species is derived from the
         !    others and constants associated with it are initialized here
         !
         if (TRIM(dry_air_species(dry_air_species_num + 1)) == 'N2') then
            call dry_air_species_info('N2', ix, mw)
            mw = 2.0_r8 * cnst_mw(ix)
            icnst = 0 ! index for the derived tracer N2
            thermodynamic_active_species_cp(icnst) = cp2 / mw
            thermodynamic_active_species_cv(icnst) = cv2 / mw !N2
            thermodynamic_active_species_R  (icnst) = r_universal / mw
            thermodynamic_active_species_mwi(icnst) = 1.0_r8 / mw
            thermodynamic_active_species_kv(icnst)  = 3.42_r8
            thermodynamic_active_species_kc(icnst)  = 56._r8
            !
            ! if last major species is not N2 then add code here
            !
         else
            write(iulog, *) subname, ' derived major species not found: ',    &
                 dry_air_species(dry_air_species_num)
            call endrun(subname//': derived major species not found')
         end if
      else
         !
         ! dry air is not species dependent
         !
         icnst = 0
         thermodynamic_active_species_cp (icnst) = cpair
         thermodynamic_active_species_cv (icnst) = cpair - rair
         thermodynamic_active_species_R  (icnst) = rair
      end if
      !
      !************************************************************************
      !
      ! add prognostic components of dry air
      !
      !************************************************************************
      !
      icnst = 1
      do idx = 1, dry_air_species_num
         select case (TRIM(dry_air_species(idx)))
            !
            ! O
            !
         case('O')
            call cnst_get_ind('O' ,ix, abort=.false.)
            if (ix < 1) then
               write(iulog, *) subname, ' dry air component not found: ',     &
                    dry_air_species(idx)
               call endrun(subname//': dry air component not found')
            else
               mw = cnst_mw(ix)
               thermodynamic_active_species_idx(icnst) = ix
               thermodynamic_active_species_cp (icnst) = 0.5_r8*r_universal*(2._r8+dof1) / mw
               thermodynamic_active_species_cv (icnst) = 0.5_r8*r_universal*dof1 / mw
               thermodynamic_active_species_R  (icnst) = r_universal / mw
               thermodynamic_active_species_mwi(icnst) = 1.0_r8 / mw
               thermodynamic_active_species_kv(icnst)  = 3.9_r8
               thermodynamic_active_species_kc(icnst)  = 75.9_r8
               icnst = icnst + 1
            end if
            !
            ! O2
            !
         case('O2')
            call cnst_get_ind('O2' ,ix, abort=.false.)
            if (ix < 1) then
               write(iulog, *) subname, ' dry air component not found: ',     &
                    dry_air_species(idx)
               call endrun(subname//': dry air component not found')
            else
               mw = cnst_mw(ix)
               thermodynamic_active_species_idx(icnst) = ix
               thermodynamic_active_species_cp (icnst) = 0.5_r8*r_universal*(2._r8+dof2) / mw
               thermodynamic_active_species_cv (icnst) = 0.5_r8*r_universal*dof2 / mw
               thermodynamic_active_species_R  (icnst) = r_universal / mw
               thermodynamic_active_species_mwi(icnst) = 1.0_r8 / mw
               thermodynamic_active_species_kv(icnst)  = 4.03_r8
               thermodynamic_active_species_kc(icnst)  = 56._r8
               icnst = icnst + 1
            end if
            !
            ! H
            !
         case('H')
            call cnst_get_ind('H' ,ix, abort=.false.)
            if (ix < 1) then
               write(iulog, *) subname, ' dry air component not found: ',     &
                    dry_air_species(idx)
               call endrun(subname//': dry air component not found')
            else
               mw = cnst_mw(ix)
               thermodynamic_active_species_idx(icnst) = ix
               thermodynamic_active_species_cp (icnst) = 0.5_r8*r_universal*(2._r8+dof1) / mw
               thermodynamic_active_species_cv (icnst) = 0.5_r8*r_universal*dof1 / mw
               thermodynamic_active_species_R  (icnst) = r_universal / mw
               thermodynamic_active_species_mwi(icnst) = 1.0_r8 / mw
               thermodynamic_active_species_kv(icnst)  = 0.0_r8
               thermodynamic_active_species_kc(icnst)  = 0.0_r8
               icnst = icnst + 1
            end if
            !
            ! If support for more major species is to be included add code here
            !
         case default
            write(iulog, *) subname, ' dry air component not found: ',        &
                 dry_air_species(idx)
            call endrun(subname//': dry air component not found')
         end select

         if (masterproc) then
            write(iulog, *) "Dry air composition ",                           &
                 TRIM(dry_air_species(idx)),                                  &
                 icnst-1,thermodynamic_active_species_idx(icnst-1),           &
                 thermodynamic_active_species_mwi(icnst-1),                   &
                 thermodynamic_active_species_cp(icnst-1),                    &
                 thermodynamic_active_species_cv(icnst-1)
         end if
      end do
      isize = dry_air_species_num+1
      icnst = 0 ! N2
      if(isize > 0) then
         if(masterproc) then
            write(iulog, *) "Dry air composition ",                           &
                 TRIM(dry_air_species(idx)),                                  &
                 icnst, -1, thermodynamic_active_species_mwi(icnst),          &
                 thermodynamic_active_species_cp(icnst),                      &
                 thermodynamic_active_species_cv(icnst)
         end if
      end if
      !
      !************************************************************************
      !
      ! Add non-dry components of moist air (water vapor and condensates)
      !
      !************************************************************************
      !
      icnst = dry_air_species_num + 1
      do idx = 1, water_species_in_air_num
         select case (TRIM(water_species_in_air(idx)))
            !
            ! Q
            !
         case('Q')
            call cnst_get_ind('Q', ix, abort=.false.)
            if (ix < 1) then
               write(iulog, *) subname, ' moist air component not found: ',   &
                    water_species_in_air(idx)
               call endrun(subname//': moist air component not found')
            else
               wv_idx = ix
               mw = cnst_mw(ix)
               thermodynamic_active_species_idx(icnst) = ix
               thermodynamic_active_species_cp (icnst) = cpwv
               thermodynamic_active_species_cv (icnst) = 0.5_r8*r_universal*dof3 / mw
               thermodynamic_active_species_R  (icnst) = rh2o
               icnst = icnst + 1
            end if
            !
            ! CLDLIQ
            !
         case('CLDLIQ')
            call cnst_get_ind('CLDLIQ' ,ix, abort=.false.)
            if (ix < 1) then
               write(iulog, *) subname, ' moist air component not found: ',   &
                    water_species_in_air(idx)
               call endrun(subname//': moist air component not found')
            else
               thermodynamic_active_species_idx(icnst) = ix
               thermodynamic_active_species_cp (icnst) = cpliq
               thermodynamic_active_species_cv (icnst) = cpliq
               liq_num           = liq_num+1
               liq_idx (liq_num) = ix
               icnst = icnst + 1
               has_liq = .true.
            end if
            !
            ! CLDICE
            !
         case('CLDICE')
            call cnst_get_ind('CLDICE' ,ix, abort=.false.)
            if (ix < 1) then
               write(iulog, *) subname, ' moist air component not found: ',   &
                    water_species_in_air(idx)
               call endrun(subname//': moist air component not found')
            else
               thermodynamic_active_species_idx(icnst) = ix
               thermodynamic_active_species_cp (icnst) = cpice
               thermodynamic_active_species_cv (icnst) = cpice
               ice_num           = ice_num+1
               ice_idx(ice_num)  = ix
               icnst = icnst + 1
               has_ice = .true.
            end if
            !
            ! RAINQM
            !
         case('RAINQM')
            call cnst_get_ind('RAINQM' ,ix, abort=.false.)
            if (ix < 1) then
               write(iulog, *) subname, ' moist air component not found: ',   &
                    water_species_in_air(idx)
               call endrun(subname//': moist air component not found')
            else
               thermodynamic_active_species_idx(icnst) = ix
               thermodynamic_active_species_cp (icnst) = cpliq
               thermodynamic_active_species_cv (icnst) = cpliq
               liq_num           = liq_num+1
               liq_idx(liq_num)  = ix
               icnst = icnst + 1
               has_liq = .true.
            end if
            !
            ! SNOWQM
            !
         case('SNOWQM')
            call cnst_get_ind('SNOWQM' ,ix, abort=.false.)
            if (ix < 1) then
               write(iulog, *) subname, ' moist air component not found: ',   &
                    water_species_in_air(idx)
               call endrun(subname//': moist air component not found')
            else
               thermodynamic_active_species_idx(icnst) = ix
               thermodynamic_active_species_cp (icnst) = cpice
               thermodynamic_active_species_cv (icnst) = cpice
               ice_num           = ice_num+1
               ice_idx(ice_num)  = ix
               icnst = icnst + 1
               has_ice = .true.
            end if
            !
            ! GRAUQM
            !
         case('GRAUQM')
            call cnst_get_ind('GRAUQM' ,ix, abort=.false.)
            if (ix < 1) then
               write(iulog, *) subname, ' moist air component not found: ',   &
                    water_species_in_air(idx)
               call endrun(subname//': moist air component not found')
            else
               mw = cnst_mw(ix)
               thermodynamic_active_species_idx(icnst) = ix
               thermodynamic_active_species_cp (icnst) = cpice
               thermodynamic_active_species_cv (icnst) = cpice
               ice_num           = ice_num+1
               ice_idx(ice_num)  = ix
               icnst = icnst + 1
               has_ice = .true.
            end if
            !
            ! If support for more major species is to be included add code here
            !
         case default
            write(iulog, *) subname, ' moist air component not found: ',      &
                 water_species_in_air(idx)
            call endrun(subname//': moist air component not found')
         end select
         !
         !
         !
         if (masterproc) then
            write(iulog, *) "Thermodynamic active species ",                  &
                 TRIM(water_species_in_air(idx))
            write(iulog, *) "   global index                           : ",   &
                 icnst-1
            write(iulog, *) "   thermodynamic_active_species_idx       : ",   &
                 thermodynamic_active_species_idx(icnst-1)
            write(iulog, *) "   cp                                     : ",   &
                 thermodynamic_active_species_cp(icnst-1)
            write(iulog, *) "   cv                                     : ",   &
                 thermodynamic_active_species_cv(icnst-1)
            if (has_liq) then
               write(iulog, *) "   register phase (liquid or ice)         :", &
                    " liquid"
            end if
            if (has_ice) then
               write(iulog, *) "   register phase (liquid or ice)         :", &
                    " ice"
            end if
            write(iulog, *) "  "
         end if
         has_liq = .false.
         has_ice = .false.
      end do

      allocate(thermodynamic_active_species_liq_idx(liq_num), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"thermodynamic_active_species_liq_idx")
      end if
      allocate(thermodynamic_active_species_liq_idx_dycore(liq_num), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"thermodynamic_active_species_liq_idx_dycore")
      end if
      allocate(thermodynamic_active_species_ice_idx(ice_num), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"thermodynamic_active_species_ice_idx")
      end if
      allocate(thermodynamic_active_species_ice_idx_dycore(ice_num), stat=ierr)
      if (ierr /= 0) then
         call endrun(errstr//"thermodynamic_active_species_ice_idx_dycore")
      end if

      thermodynamic_active_species_liq_idx = liq_idx(1:liq_num)
      thermodynamic_active_species_liq_num = liq_num

      ! array initialized by the dycore
      thermodynamic_active_species_liq_idx_dycore = -99

      thermodynamic_active_species_ice_idx = ice_idx(1:ice_num)
      thermodynamic_active_species_ice_num = ice_num

      ! array initialized by the dycore
      thermodynamic_active_species_ice_idx_dycore = -99

      if (water_species_in_air_num /= 1 + liq_num+ice_num) then
         write(iulog, '(2a,2(i0,a))') subname,                                &
              "  water_species_in_air_num = ",                                &
              water_species_in_air_num, ", should be ",              &
              (1 + liq_num + ice_num), " (1 + liq_num + ice_num)"
         call endrun(subname//': water_species_in_air_num /= 1+liq_num+ice_num')
      end if
      enthalpy_reference_state = 'ice'
      if (masterproc) then
         write(iulog, *)   'Enthalpy reference state           : ',           &
              TRIM(enthalpy_reference_state)
      end if
   end subroutine air_composition_init

   !===========================================================================

   subroutine dry_air_species_info(name, index, molec_weight, caller)
      use cam_abortutils, only: endrun
      use cam_logfile,    only: iulog
      use constituents,   only: cnst_get_ind, cnst_mw
      ! Find the constituent index of <name> and return it in
      !    <index>. Return the constituent molecular weight in
      !    <molec_weight>

      ! Dummy arguments
      character(len=*),           intent(in)  :: name
      integer,                    intent(out) :: index
      real(r8),                   intent(out) :: molec_weight
      character(len=*), optional, intent(in)  :: caller
      ! Local parameter
      character(len=*), parameter :: subname = 'dry_air_species_info: '

      call cnst_get_ind(trim(name), index, abort=.false.)
      if (index < 1) then
         if (present(caller)) then
            write(iulog, *) trim(caller), ": dry air component not found, '", &
                 trim(name), "'"
            call endrun(trim(caller)//": dry air component not found, '"//    &
                 trim(name)//"'")
         else
            write(iulog, *) subname, "dry air component not found, '",        &
                 trim(name), "'"
            call endrun(subname//"dry air component not found, '"//           &
                 trim(name)//"'")
         end if
      else
         molec_weight = cnst_mw(index)
      end if

   end subroutine dry_air_species_info


end module air_composition
