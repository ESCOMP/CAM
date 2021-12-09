
module constituents

! Metadata manager for the advected constituents.

use shr_kind_mod,     only: r8 => shr_kind_r8
use shr_const_mod,    only: shr_const_rgas
use spmd_utils,       only: masterproc
use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog

implicit none
private
save

! Public interfaces
public :: &
   cnst_readnl,         &! read namelist
   cnst_add,            &! add a constituent to the list of advected constituents
   cnst_num_avail,      &! returns the number of available slots in the constituent array
   cnst_get_ind,        &! get the index of a constituent
   cnst_get_type_byind, &! return mixing ratio type of a constituent
   cnst_get_molec_byind,&! return molecular diffusion type of a constituent
   cnst_read_iv,        &! query whether constituent initial values are read from initial file
   cnst_chk_dim,        &! check that number of constituents added equals dimensions (pcnst)
   cnst_cam_outfld,     &! Returns true if default CAM output was specified in the cnst_add calls.
   cnst_set_spec_class, &! Sets the type of species class
   cnst_is_a_water_species,&! Returns true for constituents identified as water species
   cnst_set_convtran2    ! Override for convtran2 values set by the cnst_add routine

! Public data

integer, parameter, public :: pcnst  = PCNST      ! number of advected constituents (including water vapor)

character(len=16), public, protected :: cnst_name(pcnst)     ! constituent names
character(len=128),public, protected :: cnst_longname(pcnst) ! long name of constituents

! Namelist variables
logical, public, protected :: readtrace = .true.  ! true => obtain initial tracer data from IC file

integer, public, parameter :: cnst_spec_class_undefined   = 0
integer, public, parameter :: cnst_spec_class_cldphysics  = 1
integer, public, parameter :: cnst_spec_class_aerosol     = 2
integer, public, parameter :: cnst_spec_class_gas         = 3
integer, public, parameter :: cnst_spec_class_other       = 4

!
! Constants for each tracer

integer, public, protected :: cnst_species_class(pcnst) = cnst_spec_class_undefined  ! indicates species class &
                                                                                       ! (cldphysics, aerosol, gas )
real(r8),    public :: cnst_cp  (pcnst)          ! specific heat at constant pressure (J/kg/K)
real(r8),    public :: cnst_cv  (pcnst)          ! specific heat at constant volume (J/kg/K)
real(r8),    public :: cnst_mw  (pcnst)          ! molecular weight (kg/kmole)
character*3, public, protected :: cnst_type(pcnst)! wet or dry mixing ratio
character*5, public :: cnst_molec(pcnst)         ! major or minor species molecular diffusion
real(r8),    public :: cnst_rgas(pcnst)          ! gas constant ()
real(r8),    public :: qmin     (pcnst)          ! minimum permitted constituent concentration (kg/kg)
real(r8),    public :: qmincg   (pcnst)          ! for backward compatibility only
logical,     public :: cnst_fixed_ubc(pcnst) = .false.  ! upper bndy condition = fixed ?
logical,     public :: cnst_fixed_ubflx(pcnst) = .false.! upper boundary non-zero fixed constituent flux
logical, public, protected :: cnst_is_convtran1(pcnst) = .false.  ! do convective transport in phase 1
logical, public, protected :: cnst_is_convtran2(pcnst) = .false.  ! do convective transport in phase 2

!++bee - temporary... These names should be declared in the module that makes the addfld and outfld calls.
! Lists of tracer names and diagnostics
character(len=16), public :: apcnst    (pcnst)   ! constituents after physics  (FV core only)
character(len=16), public :: bpcnst    (pcnst)   ! constituents before physics (FV core only)
character(len=16), public :: hadvnam   (pcnst)   ! names of horizontal advection tendencies
character(len=16), public :: vadvnam   (pcnst)   ! names of vertical advection tendencies
character(len=16), public :: dcconnam  (pcnst)   ! names of convection tendencies
character(len=16), public :: fixcnam   (pcnst)   ! names of species slt fixer tendencies
character(len=16), public :: tendnam   (pcnst)   ! names of total tendencies of species
character(len=16), public :: ptendnam  (pcnst)   ! names of total physics tendencies of species
character(len=16), public :: dmetendnam(pcnst)   ! names of dme adjusted tracers (FV)
character(len=16), public :: sflxnam   (pcnst)   ! names of surface fluxes of species
character(len=16), public :: tottnam   (pcnst)   ! names for horz + vert + fixer tendencies

! Private data

integer :: padv = 0                      ! index pointer to last advected tracer
logical :: read_init_vals(pcnst)         ! true => read initial values from initial file
logical :: cam_outfld_(pcnst)            ! true  => default CAM output of constituents in kg/kg
                                         ! false => chemistry is responsible for making outfld
                                         !          calls for constituents

!==============================================================================================
CONTAINS
!==============================================================================================

subroutine cnst_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_logical


   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: sub = 'cnst_readnl'

   namelist /constituents_nl/ readtrace
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'constituents_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, constituents_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub//': FATAL: reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(readtrace, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: readtrace")

   if (masterproc) then
      write(iulog,*)'Summary of constituent module options:'
      write(iulog,*)'  Read constituent initial values from initial file by default: ', readtrace
   end if

end subroutine cnst_readnl

!=========================================================================================


subroutine cnst_add (name, mwc, cpc, qminc, &
                     ind, longname, readiv, mixtype, molectype, cam_outfld, &
                     fixed_ubc, fixed_ubflx, is_convtran1, is_convtran2, cnst_spec_class)

   ! Register a constituent.

   character(len=*), intent(in) :: &
      name      ! constituent name used as variable name in history file output (8 char max)
   real(r8),intent(in)    :: mwc    ! constituent molecular weight (kg/kmol)
   real(r8),intent(in)    :: cpc    ! constituent specific heat at constant pressure (J/kg/K)
   real(r8),intent(in)    :: qminc  ! minimum value of mass mixing ratio (kg/kg)
                                    ! normally 0., except water 1.E-12, for radiation.
   integer, intent(out)   :: ind    ! global constituent index (in q array)

   character(len=*), intent(in), optional :: &
      longname    ! value for long_name attribute in netcdf output (128 char max, defaults to name)
   logical,          intent(in), optional :: &
      readiv      ! true => read initial values from initial file (default: true)
   character(len=*), intent(in), optional :: &
      mixtype     ! mixing ratio type (dry, wet)
   character(len=*), intent(in), optional :: &
      molectype     ! molecular diffusion type (minor, major)
   logical,          intent(in), optional :: &
      cam_outfld  ! true => default CAM output of constituent in kg/kg
   logical,          intent(in), optional :: &
      fixed_ubc   ! true => const has a fixed upper bndy condition
   logical,          intent(in), optional :: &
      fixed_ubflx ! true => const has a non-zero fixed upper bndy flux value
   logical,          intent(in), optional :: &
      is_convtran1 ! true => convective transport in convtran1
   logical,          intent(in), optional :: &
      is_convtran2 ! true => convective transport in convtran2
   integer,          intent(in), optional :: &
      cnst_spec_class ! type of species class

   character(len=*), parameter :: sub='cnst_add'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   ! set tracer index and check validity
   padv = padv+1
   ind  = padv
   if (padv > pcnst) then
      write(errmsg, *) sub//': FATAL: advected tracer index greater than pcnst=', pcnst
      call endrun(errmsg)
   end if

   ! set tracer name and constants
   cnst_name(ind) = name
   if (present(longname)) then
      cnst_longname(ind) = longname
   else
      cnst_longname(ind) = name
   end if

   ! set whether to read initial values from initial file
   if (present(readiv)) then
      read_init_vals(ind) = readiv
   else
      read_init_vals(ind) = readtrace
   end if

   ! set constituent mixing ratio type
   if (present(mixtype)) then
      cnst_type(ind) = mixtype
   else
      cnst_type(ind) = 'wet'
   end if

   ! set constituent molecular diffusion type
   if (present(molectype)) then
      cnst_molec(ind) = molectype
   else
      cnst_molec(ind) = 'minor'
   end if

   ! set outfld type 
   ! (false: the module declaring the constituent is responsible for outfld calls)
   if (present(cam_outfld)) then
      cam_outfld_(ind) = cam_outfld
   else
      cam_outfld_(ind) = .true.
   end if

   ! set upper boundary condition type
   if (present(fixed_ubc)) then
      cnst_fixed_ubc(ind) = fixed_ubc
   else
      cnst_fixed_ubc(ind) = .false.
   end if

   ! set upper boundary flux type
   if (present(fixed_ubflx)) then
      cnst_fixed_ubflx(ind) = fixed_ubflx
   else
      cnst_fixed_ubflx(ind) = .false.
   end if

   ! Set flag for convective transport by first call to convtran (phase 1).
   if (present(is_convtran1)) then
      cnst_is_convtran1(ind) = is_convtran1
   else
      cnst_is_convtran1(ind) = .false.
   end if
   ! Set flag for convective transport after wetdep (phase 2).
   if (present(is_convtran2)) then
      cnst_is_convtran2(ind) = is_convtran2
   else
      ! The default is that all constituents except water vapor are transported in phase 2
      ! unless they were transported in phase 1 (typically the condensed water species)
      if (ind > 1) cnst_is_convtran2(ind) = .not. cnst_is_convtran1(ind)
   end if
   ! consistency check -- It is OK to completely turn off the deep scheme transport by setting
   ! both cnst_is_convtran1 and cnst_is_convtran2 to FALSE.  But it is an error to
   ! have both set TRUE.
   if (cnst_is_convtran1(ind) .and. cnst_is_convtran2(ind)) then
      call endrun(sub//': FATAL: cannot set both cnst_is_convtran1 and cnst_is_convtran2 to TRUE')
   end if

   ! Set type for species class
   if ( present(cnst_spec_class) ) then
      cnst_species_class(ind) = cnst_spec_class
   else
      cnst_species_class(ind) = cnst_spec_class_undefined
   end if

   cnst_cp  (ind) = cpc
   cnst_mw  (ind) = mwc
   qmin     (ind) = qminc
   if (ind == 1) then
      ! qmincg for water vapor set to zero
      qmincg(ind) = 0._r8
   else
      qmincg(ind) = qminc
   end if

   cnst_rgas(ind) = shr_const_rgas * mwc
   cnst_cv  (ind) = cpc - cnst_rgas(ind)

end subroutine cnst_add

!----------------------------------------------------------------------------------------------

subroutine cnst_set_convtran2(ind, is_convtran2)

   ! Allow user to override the value of cnst_is_convtran2 set by a previous cnst_add call.

   integer, intent(in) :: ind          ! global constituent index (in q array)
   logical, intent(in) :: is_convtran2 ! true => convect in convtran2

   character(len=*), parameter :: sub = 'cnst_set_convtran2'
   character(len=128)          :: errmsg
   !-----------------------------------------------------------------------

   ! check index
   if (ind <= 0 .or. ind > padv) then
      write(errmsg,*) sub//': FATAL: bad tracer index: padv, ind = ', padv, ind
      call endrun(errmsg)
   end if

   ! Set flag for convective transport after wetdep (phase 2).
   cnst_is_convtran2(ind) = is_convtran2

   ! consistency check -- It is OK to completely turn off the tracer convection by setting
   ! both cnst_is_convtran1 and cnst_is_convtran2 to FALSE.  But it is an error to
   ! have both set TRUE.
   if (cnst_is_convtran1(ind) .and. cnst_is_convtran2(ind)) then
      call endrun(sub//': FATAL: cannot set both cnst_is_convtran1 and cnst_is_convtran2 to TRUE')
   end if

end subroutine cnst_set_convtran2

!----------------------------------------------------------------------------------------------

subroutine cnst_set_spec_class(ind, cnst_spec_class_in)

   ! Allow user to override the value of cnst_spec_class set by a previous cnst_add call.

   integer, intent(in) :: ind                ! global constituent index (in q array)
   integer, intent(in) :: cnst_spec_class_in ! species class designator

   character(len=*), parameter :: subname = 'cnst_set_spec_class'
   !-----------------------------------------------------------------------

   ! check index
    if (ind <= 0 .or. ind > padv) then
       write(iulog,*) subname//': illegal tracer index: padv, ind = ', padv, ind
       call endrun(subname//': illegal tracer index')
    end if
    
    ! Check designator
    if (cnst_spec_class_in /= cnst_spec_class_undefined  .and. &
        cnst_spec_class_in /= cnst_spec_class_cldphysics .and. &
        cnst_spec_class_in /= cnst_spec_class_aerosol    .and. &
        cnst_spec_class_in /= cnst_spec_class_gas        .and. &
        cnst_spec_class_in /= cnst_spec_class_other ) then
          write(iulog,*) subname//': trying to use invalid cnst_spec_class designator', cnst_spec_class_in
          call endrun(subname//': invalid cnst_spec_class designator')
    end if

    ! Set flag for convective transport after wetdep (phase 2).
    cnst_species_class(ind) = cnst_spec_class_in

 end subroutine cnst_set_spec_class

!==============================================================================

function cnst_num_avail()

   ! return number of available slots in the constituent array

   integer cnst_num_avail

   cnst_num_avail = pcnst - padv

end function cnst_num_avail

!==============================================================================

subroutine cnst_get_ind (name, ind, abort)

   ! Get the index of a constituent.  Optional abort argument allows returning
   ! control to caller when constituent name is not found.  Default behavior is
   ! to call endrun when name is not found.

   !-----------------------------Arguments---------------------------------
   character(len=*),  intent(in)  :: name  ! constituent name
   integer,           intent(out) :: ind   ! global constituent index (in q array)
   logical, optional, intent(in)  :: abort ! optional flag controlling abort

   !---------------------------Local workspace-----------------------------
   integer :: m                                   ! tracer index
   logical :: abort_on_error
   character(len=*), parameter :: sub='cnst_get_ind'
   !-----------------------------------------------------------------------

   ! Find tracer name in list
   do m = 1, pcnst
      if (name == cnst_name(m)) then
         ind  = m
         return
      end if
   end do

   ! Unrecognized name
   abort_on_error = .true.
   if (present(abort)) abort_on_error = abort

   if (abort_on_error) then
      write(iulog, *) sub//': FATAL: name:', name,  ' not found in list:', cnst_name(:)
      call endrun(sub//': FATAL: name not found')
   end if

   ! error return
   ind = -1

end subroutine cnst_get_ind

!==============================================================================================

character*3 function cnst_get_type_byind(ind)

   ! Return the mixing ratio type of a constituent 

   !-----------------------------Arguments---------------------------------
   integer, intent(in)   :: ind    ! global constituent index (in q array)

   !---------------------------Local workspace-----------------------------
   character(len=*), parameter :: sub='cnst_get_type_byind'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   if (ind > 0 .and. ind <= pcnst) then
      cnst_get_type_byind = cnst_type(ind)
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value for constituent index=', ind
      call endrun(errmsg)
   end if

end function cnst_get_type_byind

!==============================================================================================

character*5 function cnst_get_molec_byind (ind)

   ! Return the molecular diffusion type of a constituent 

   !-----------------------------Arguments---------------------------------
   integer, intent(in)   :: ind    ! global constituent index (in q array)

   !---------------------------Local workspace-----------------------------
   character(len=*), parameter :: sub='cnst_get_molec_byind'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   if (ind > 0 .and. ind <= pcnst) then
      cnst_get_molec_byind = cnst_molec(ind)
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value for constituent index=', ind
      call endrun(errmsg)
   end if

end function cnst_get_molec_byind

!==============================================================================

function cnst_read_iv(m)

   ! Query whether to attempt to read constituent initial values from initial file.

   !-----------------------------Arguments---------------------------------
   integer, intent(in) :: m    ! constituent index

   logical :: cnst_read_iv     ! true => read initial values from inital file

   character(len=*), parameter :: sub='cnst_read_iv'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   if (m > 0 .and. m <= pcnst) then
      cnst_read_iv = read_init_vals(m)
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value for constiuent index=', m
      call endrun(errmsg)
   end if

end function cnst_read_iv

!==============================================================================

subroutine cnst_chk_dim

   ! Check that the number of registered constituents is pcnst
   ! Write constituent list to log file.

   integer :: i, m
   character(len=*), parameter :: sub='cnst_chk_dim'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   if (padv /= pcnst) then
      write(errmsg, *) sub//': FATAL: number of advected tracer (',padv, &
                       ') not equal to pcnst (', pcnst, ')'
      call endrun (errmsg)
   endif

   if (masterproc) then
      write(iulog,*) 'Advected constituent list:'
      do i = 1, pcnst
         write(iulog,'(2x,i4,2x,a8,2x,a128,2x,a3)') i, cnst_name(i), cnst_longname(i), &
                                                    cnst_type(i)
      end do
   end if

   ! Set names of advected tracer diagnostics
   do m=1,pcnst
      apcnst    (m)  = trim(cnst_name(m))//'AP'
      bpcnst    (m)  = trim(cnst_name(m))//'BP'
      hadvnam   (m)  = 'HA'//cnst_name(m)
      vadvnam   (m)  = 'VA'//cnst_name(m)
      fixcnam   (m)  = 'DF'//cnst_name(m)
      tendnam   (m)  = 'TE'//cnst_name(m)
      ptendnam  (m)  = 'PTE'//cnst_name(m)
      dmetendnam(m)  = 'DME'//cnst_name(m)
      tottnam   (m)  = 'TA'//cnst_name(m)
      sflxnam(m)     = 'SF'//cnst_name(m)
   end do

end subroutine cnst_chk_dim

!==============================================================================

function cnst_cam_outfld(m)

   ! Query whether default CAM outfld calls should be made.

   !----------------------------------------------------------------------- 
   integer, intent(in) :: m                ! constituent index

   logical             :: cnst_cam_outfld  ! true => use default CAM outfld calls

   character(len=*), parameter :: sub='cnst_cam_outfld'
   character(len=128) :: errmsg
   !-----------------------------------------------------------------------

   if (m > 0 .and. m <= pcnst) then
      cnst_cam_outfld = cam_outfld_(m)
   else
      ! index out of range
      write(errmsg,*) sub//': FATAL: bad value for constiuent index=', m
      call endrun(errmsg)
   end if

end function cnst_cam_outfld

!==============================================================================

pure logical function cnst_is_a_water_species(name)

   ! test whether the input name matches the name of a water species

   character(len=*), intent(in) :: name  
   !-------------------------------------------------------------------------

   cnst_is_a_water_species = .false.

   if (name == 'Q'      .or. &
       name == 'CLDLIQ' .or. &
       name == 'CLDICE' .or. &
       name == 'RAINQM' .or. &
       name == 'SNOWQM' .or. &
       name == 'GRAUQM'      ) cnst_is_a_water_species = .true.
      
end function cnst_is_a_water_species

!==============================================================================

end module constituents
