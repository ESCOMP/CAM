!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
!
! !MODULE: carma_aero_gasaerexch --- does carma aerosol gas-aerosol exchange for SOA
!
! !INTERFACE:
module carma_aero_gasaerexch

! !USES:
  use shr_kind_mod,    only:  r8 => shr_kind_r8
  use chem_mods,       only:  gas_pcnst
  use ref_pres,        only:  top_lev => clim_modal_aero_top_lev
  use ppgrid,          only:  pcols, pver
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_info_by_bin, rad_cnst_get_bin_props_by_idx, &
                              rad_cnst_get_info_by_bin_spec
  use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

  implicit none
  private
  public :: carma_aero_gasaerexch_sub
  public :: carma_aero_gasaerexch_init

  !PUBLIC DATA MEMBERS:

  ! description of bin aerosols
  integer, public, protected :: nspec_max = 0
  integer, public, protected :: nbins = 0
  integer, public, protected :: nsoa_vbs = 0
  integer, public, protected :: nsoa = 0
  integer, public, protected :: npoa = 0
  integer, public, protected, allocatable :: nspec(:)

  ! Misc private data
  character(len=32), allocatable :: fldname(:)    ! names for interstitial output fields
  character(len=32), allocatable :: fldname_cw(:)    ! names for cloud_borne output fields

  ! local indexing for bins
  integer, allocatable :: bin_idx(:,:) ! table for local indexing of modal aero number and mmr
  integer :: ncnst_tot                  ! total number of mode number conc + mode species

  real(r8) :: mw_soa = 250._r8
  integer :: fracvbs_idx = -1
  integer, allocatable :: dqdtsoa_idx(:,:)
  integer, allocatable  :: cnsoa(:)         ! true if soa gas is a species and carma soa in bin
  integer, allocatable  :: cnpoa(:)         ! true if soa gas is a species and carma soa in bin
  integer, allocatable  :: l_soag(:)         ! true if soa gas is a species and carma soa in bin

  logical, allocatable  :: do_soag_any(:)         ! true if soa gas is a species and carma soa in bin
! !DESCRIPTION: This module implements ...
!
! !REVISION HISTORY:
!
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! list private module data here

!EOC
!----------------------------------------------------------------------


contains

!----------------------------------------------------------------------

  subroutine carma_aero_gasaerexch_init

!-----------------------------------------------------------------------
!
! Purpose:
!      gas-aerosol exchange SOAG <-> soa
!
! Author: Simone Tilmes
!
!-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default, fieldname_len, horiz_only
    use constituents,   only: pcnst, cnst_name
    use phys_control,   only: phys_getopts
    use mo_chem_utls,   only: get_spc_ndx

!-----------------------------------------------------------------------
! arguments

!-----------------------------------------------------------------------
! local
    integer  :: j
    integer  :: i, ii
    integer  :: l
    integer  :: m
    integer  :: ns
    character(len=fieldname_len+3) :: fieldname
    character(len=32)              :: spectype
    character(len=32)              :: spec_name
    character(128)                 :: long_name
    character(8)                   :: unit
    character(len=2)               :: outsoa

    logical                        :: history_aerosol      ! Output the MAM aerosol tendencies
    !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out        = history_aerosol   )

    !
    ! get info about the bin aerosols
    ! get nbins

    call rad_cnst_get_info( 0, nbins=nbins)

    allocate( nspec(nbins) )
    allocate( cnsoa(nbins) )
    allocate( cnpoa(nbins) )

    do m = 1, nbins
       call rad_cnst_get_info_by_bin(0, m, nspec=nspec(m))
    end do

    nspec_max = maxval(nspec)

    ncnst_tot = nspec(1)
    do m = 2, nbins
       ncnst_tot = ncnst_tot + nspec(m)
    end do

    allocate(  bin_idx(nbins,nspec_max), &
               do_soag_any(nbins),       &
               fldname_cw(ncnst_tot),    &
               fldname(ncnst_tot) )

    ! Local indexing compresses the mode and number/mass indicies into one index.
    ! This indexing is used by the pointer arrays used to reference state and pbuf
    ! fields.
    ! for CARMA we add number = 0, total mass = 1, and mass from each constituence into mm.
    ii = 0
    do m = 1, nbins
       do l = 1, nspec(m)    ! do through nspec
          ii = ii + 1
          bin_idx(m,l) = ii
       end do
    end do

    ! SOAG / SOA / POM information
    ! Define number of VBS bins (nsoa) based on number of SOAG chemistry species

    nsoa_vbs = 0
    do i = 1, pcnst
       if (cnst_name(i)(:4) == 'SOAG') then
          nsoa_vbs = nsoa_vbs + 1
       end if
    end do
    allocate( l_soag(nsoa_vbs) )
    nsoa_vbs = 0
    do i = 1, pcnst
       if (cnst_name(i)(:4) == 'SOAG') then
          nsoa_vbs = nsoa_vbs + 1
          l_soag(nsoa_vbs) = get_spc_ndx(cnst_name(i))
       end if
    end do

    fracvbs_idx = pbuf_get_index('FRACVBS')

    ! identify number of SOA and POA in CARMA code (CARMA number cn)
    do m = 1, nbins
       cnsoa(m) = 0
       cnpoa(m) = 0
       do l = 1, nspec(m)
          call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
          if (trim(spectype) == 's-organic') then
             cnsoa(m) = cnsoa(m) + 1
          end if
          if (trim(spectype) == 'p-organic') then
             cnpoa(m) = cnpoa(m) + 1
          end if
       end do
    end do
    ! some bins don't contain soa or poa
    nsoa= maxval(cnsoa)
    npoa= maxval(cnpoa)

    allocate( dqdtsoa_idx(nbins,nsoa)                       )
    do m = 1, nbins
       ns = 0
       do l = 1, nspec(m)
          call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
          if (trim(spectype) == 's-organic') then
             call rad_cnst_get_info_by_bin_spec(0, m, l, spec_name=spec_name)
             ns = ns + 1
             dqdtsoa_idx(m,ns) = pbuf_get_index('DQDT_'//trim(spec_name))
          end if
       end do
    end do

    do m = 1, nbins
       do_soag_any(m) = cnsoa(m)>0
    end do

!---------define history fields for new cond/evap diagnostics----------------------------------------

    fieldname=trim('qcon_gaex')
    long_name = trim('3D fields for SOA condensation')
    unit = 'kg/kg/s'
    call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
    if ( history_aerosol ) then
       call add_default( fieldname,  1, ' ' )
    endif

    do j = 1, nsoa_vbs
       write (outsoa, "(I2.2)") j
       fieldname = 'qcon_gaex'//outsoa
       long_name = '3D fields for SOA condensation for VBS bin'//outsoa
       call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
       if ( history_aerosol ) then
          call add_default( fieldname,  1, ' ' )
       endif
       fieldname = 'qevap_gaex'//outsoa
       long_name = '3D fields for SOA evaporation for VBS bin'//outsoa
       call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
       if ( history_aerosol ) then
          call add_default( fieldname,  1, ' ' )
       endif
    end do

    fieldname=trim('qevap_gaex')
    long_name = trim('3D fields for SOA evaporation')
    call addfld(fieldname, (/'lev'/), 'A', unit, long_name )
    if ( history_aerosol ) then
       call add_default( fieldname,  1, ' ' )
    endif

!------------------------------------------------------------------------------

!  define history fields for basic gas-aer exchange
    do m = 1, nbins
       do l = 1, nspec(m)    ! do through nspec
          ii = bin_idx(m,l)
          if (l <= nspec(m) ) then   ! species
             call rad_cnst_get_info_by_bin_spec(0, m, l, spec_name=fldname(ii) )
             ! only write out SOA exchange here
             call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
             if (trim(spectype) == 's-organic') then
                fieldname= trim(fldname(ii)) // '_sfgaex1'
                long_name = trim(fldname(ii)) // ' gas-aerosol-exchange primary column tendency'
                unit = 'kg/m2/s'
                call addfld( fieldname, horiz_only, 'A', unit, long_name )
                if ( history_aerosol ) then
                   call add_default( fieldname, 1, ' ' )
                endif
             end if
          end if
       end do

       write(fieldname,'("WETRAD_bin",I2.2)') m
       write(long_name,'("bin ",I2.2," wet radius in carma_aero_gasaerexch")') m

       call addfld(fieldname, (/'lev'/), 'A', 'cm', long_name )
       if ( history_aerosol ) then
          call add_default( fieldname,  1, ' ' )
       endif

       write(fieldname,'("UPTKRATE_bin",I2.2)') m
       write(long_name,'("bin ",I2.2," up take rate in carma_aero_gasaerexch")') m

       call addfld(fieldname, (/'lev'/), 'A', 'sec-1', long_name )
       if ( history_aerosol ) then
          call add_default( fieldname,  1, ' ' )
       endif

       write(fieldname,'("NUMDENS_bin",I2.2)') m
       write(long_name,'("bin ",I2.2," number density carma_aero_gasaerexch")') m

       call addfld(fieldname, (/'lev'/), 'A', 'm-3', long_name )
       if ( history_aerosol ) then
          call add_default( fieldname,  1, ' ' )
       endif

    end do

    fieldname=trim('UPTKRATE')
    long_name = trim('total uptake rate in carma_aero_gasaerexch')
    call addfld(fieldname, (/'lev'/), 'A', 'sec-1', long_name )
    if ( history_aerosol ) then
       call add_default( fieldname,  1, ' ' )
    endif


  end subroutine carma_aero_gasaerexch_init


!----------------------------------------------------------------------

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  carma_aero_gasaerexch_sub --- ...
!
! !INTERFACE:
subroutine carma_aero_gasaerexch_sub(  state, &
                        pbuf, lchnk,    ncol,     nstep,         &
                        loffset,  deltat,   mbar,                &
                        t,        pmid,     pdel,                &
                        qh2o,               troplev,             &
                        q,                  raervmr,             &
                        wetr_n                                   )

  ! !USES:
  use cam_history,       only: outfld, fieldname_len
  use physconst,         only: gravit, mwdry
  use cam_abortutils,    only: endrun
  use time_manager,      only: is_first_step
  use carma_aerosol_state_mod, only: carma_aerosol_state
  use physics_types,     only: physics_state
  use physconst, only: mwdry, rair

! !PARAMETERS:
  type(physics_state), target, intent(in) :: state    ! Physics state variables
  type(physics_buffer_desc), pointer :: pbuf(:)

  integer,  intent(in)    :: lchnk                ! chunk identifier
  integer,  intent(in)    :: ncol                 ! number of atmospheric column
  integer,  intent(in)    :: nstep                ! model time-step number
  integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
  integer,  intent(in)    :: troplev(pcols)       ! tropopause vertical index
  real(r8), intent(in)    :: deltat               ! time step (s)
  real(r8), intent(in)    :: mbar(ncol,pver)      ! mean wet atmospheric mass ( amu )

  real(r8), intent(inout) :: q(ncol,pver,gas_pcnst) ! tracer mixing ratio (TMR) array
  ! *** MUST BE  #/kmol-air for number
  ! *** MUST BE mol/mol-air for mass
  ! *** NOTE ncol dimension
  real(r8), intent(in)    :: raervmr (ncol,pver,ncnst_tot) ! aerosol mixing rations (vmr)
  real(r8), intent(in)    :: t(pcols,pver)        ! temperature at model levels (K)
  real(r8), intent(in)    :: pmid(pcols,pver)     ! pressure at model levels (Pa)
  real(r8), intent(in)    :: pdel(pcols,pver)     ! pressure thickness of levels (Pa)
  real(r8), intent(in)    :: qh2o(pcols,pver)     ! water vapor mixing ratio (kg/kg)
  real(r8), intent(in)    :: wetr_n(pcols,pver,nbins) !wet geo. mean dia. (cm) of number distrib.

! !DESCRIPTION:
! this version does only do condensation for SOA for CARMA
!     method_soa=0 is no uptake
!     method_soa=1 is irreversible uptake done like h2so4 uptake
!     method_soa=2 is reversible uptake using subr carma_aero_soaexch
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
  integer, parameter :: ldiag1=-1, ldiag2=-1, ldiag3=-1, ldiag4=-1
  integer, parameter :: method_soa = 2

  real (r8), parameter :: mw_poa_host = 12.0_r8    ! molec wght of poa used in host code
  real (r8), parameter :: mw_soa_host = 250.0_r8   ! molec wght of soa used in host code

  integer :: i
  integer :: j, jsoa
  integer :: k
  integer :: l
  integer :: mm, m, n, nn, niter, niter_max

  character(len=fieldname_len+3) :: fieldname
  character(len=32)              :: spectype
  character(len=2)               :: outsoa

  real (r8) :: avg_uprt_soa(nsoa_vbs)
  real (r8) :: deltatxx
  real (r8) :: dqdt_soa_vbs(nbins,nsoa_vbs)
  real (r8) :: dqdt_soa_all(pcols,pver,nbins,nsoa)
  real (r8) :: dqdt_soag(nsoa_vbs)
  real (r8) :: fgain_soa(nbins,nsoa_vbs)
  real (r8) :: pdel_fac
  real (r8) :: num_bin(pcols,pver,nbins)
  real (r8) :: soa_vbs(pcols,pver,nbins,nsoa_vbs)
  real (r8) :: soa_c(pcols,pver,nbins,nsoa) ! SOA from CARMA
  real (r8) :: poa_c(pcols,pver,nbins,npoa) ! POA from CARMA
  real (r8) :: qold_poa(nbins,npoa)         ! POA from CARMA old
  real (r8) :: qold_soa(nbins,nsoa_vbs)     ! SOA on VBS bins  old
  real (r8) :: qnew_soa_vbs(nbins,nsoa_vbs) ! SOA on VBS bins  new
  real (r8) :: qnew_soa(nbins)              ! SOA new for combined VBS bin new for combined VBS binss
  real (r8) :: qold_soag(nsoa_vbs)
  real (r8) :: sum_dqdt_soa(nsoa_vbs)     !   sum_dqdt_soa = soa tendency from soa  gas uptake (mol/mol/s)
  real (r8) :: sum_uprt_soa(nsoa_vbs)     ! total soa uptake rate over all bin, for each soa vbs bin
  real (r8) :: uptkrate(pcols,pver,nbins)
  real (r8) :: uptkrate_all(pcols,pver)
  real (r8) :: uptkratebb(nbins)
  real (r8) :: uptkrate_soa(nbins,nsoa_vbs)
  ! gas-to-aerosol mass transfer rates (1/s)

  integer, parameter :: nsrflx = 1 ! only one dimension of qsrflx, no renaming or changes in size for CARMA currently
  real(r8) :: dqdt(ncol,pver,gas_pcnst)  ! TMR "delta q" array - NOTE dims
  real(r8) :: qsrflx(pcols,nbins,nsoa)
  ! process-specific column tracer tendencies
  ! (1=gas condensation)
  real(r8) :: qcon_vbs(pcols,pver,nsoa_vbs)
  real(r8) :: qevap_vbs(pcols,pver,nsoa_vbs)
  real(r8) :: qcon(pcols,pver)
  real(r8) :: qevap(pcols,pver)
  real(r8) :: total_soag
  real(r8) :: soag(nsoa_vbs)

  real(r8), pointer :: frac_vbs(:,:,:,:)   ! fraction of vbs SOA bins to total SOA
  real(r8), pointer :: dqdt_soa(:,:)

  real(r8) :: rhoair(pcols,pver)
  real(r8), pointer :: nmr(:,:)
  type(carma_aerosol_state), pointer :: aero_state

!----------------------------------------------------------------------
   aero_state => carma_aerosol_state(state, pbuf)

!  map CARMA soa to working soa(nbins,nsoa)

  call pbuf_get_field(pbuf, fracvbs_idx, frac_vbs)

  num_bin(:,:,:) = 0._r8
  soa_c(:,:,:,:) = 0._r8
  poa_c(:,:,:,:) = 0._r8

  rhoair(:ncol,:) = pmid(:ncol,:)/(rair*t(:ncol,:))   ! (kg-air/m3)

  do m = 1, nbins      ! main loop over aerosol bins
     if (do_soag_any(m)) then  ! only bins that contain soa
        n = 0
        nn = 0
        do l = 1, nspec(m)
           mm = bin_idx(m, l)
           call rad_cnst_get_bin_props_by_idx(0, m, l, spectype=spectype)
           if (trim(spectype) == 's-organic') then
              n = n + 1
              soa_c(:ncol,:,m,n) = raervmr(:ncol,:,mm)
           end if
           if (trim(spectype) == 'p-organic') then
              nn = nn + 1
              poa_c(:ncol,:,m,nn) = raervmr(:ncol,:,mm)
           end if
        end do
        if (npoa .gt. 1) then
           call endrun( 'carma_aero_gasaerexch_sub error: CARMA currently only supports 1 POA element' )
        end if

        if (nsoa_vbs.eq.nsoa) then
           soa_vbs(:ncol,:,:,:) = soa_c(:ncol,:,:,:)
        else
           if (nsoa.eq.1) then
              if (is_first_step()) then
                 !first time step initialization only
                 do k=top_lev,pver
                    do i=1,ncol
                       total_soag = 0.0_r8
                       do j = 1, nsoa_vbs
                          soag(j) = q(i,k,l_soag(j))
                          total_soag = total_soag + soag(j)
                       end do
                       if (total_soag .gt. 0.0_r8) then
                          do j= 1, nsoa_vbs
                             frac_vbs(i,k,m,j) = soag(j)/total_soag
                          end do
                       end if
                    end do
                 end do
              end if
              ! end first time step, after that use fraction from previous time step
              do k=top_lev,pver
                 do i=1,ncol
                    do j= 1, nsoa_vbs
                       soa_vbs(i,k,m,j) = frac_vbs(i,k,m,j)*soa_c(i,k,m,nsoa)
                    end do
                 end do
              end do
           else
              ! error message this code only works if SOAG and SOA CARMA have the same number of species,
              ! or if SOA CARMA has only one species.
              call endrun( 'carma_aero_gasaerexch_sub error in number of SOA species' )
           end if

        end if

        ! get bin number densities for all aerosols
        call aero_state%get_ambient_num(m,nmr) !  #/kg
        num_bin(:ncol,:,m) = nmr(:ncol,:)*rhoair(:ncol,:) ! #/m3

     end if
  end do


! SOA will be updated in CARMA

! zero out tendencies and other
  dqdt(:,:,:) = 0.0_r8
  qsrflx(:,:,:) = 0.0_r8

!-------Initialize evap/cond diagnostics (ncols x pver)-----------
  qcon_vbs(:,:,:) = 0.0_r8
  qevap_vbs(:,:,:) = 0.0_r8
  qcon(:,:) = 0.0_r8
  qevap(:,:) = 0.0_r8
!---------------------------------------------------
! compute gas-to-aerosol mass transfer rates
! check if only number is needed for this calculatuion!
  call gas_aer_uptkrates( ncol,       loffset,                &
                          num_bin,          t,          pmid,       &
                          wetr_n,                   uptkrate    )

  do m = 1, nbins

     write(fieldname,'("NUMDENS_bin",I2.2)') m
     call outfld(fieldname, num_bin(:ncol,:,m), ncol, lchnk )

     write(fieldname,'("WETRAD_bin",I2.2)') m
     call outfld(fieldname, wetr_n(:ncol,:,m), ncol, lchnk )

     write(fieldname,'("UPTKRATE_bin",I2.2)') m
     call outfld(fieldname, uptkrate(:ncol,:,m), ncol, lchnk )

     uptkrate_all(:ncol,:) = uptkrate_all(:ncol,:) + uptkrate(:ncol,:,m)
  end do

  fieldname = trim('UPTKRATE')
  call outfld(fieldname, uptkrate_all(:ncol,:), ncol, lchnk )

! use this for tendency calcs to avoid generating very small negative values
  deltatxx = deltat * (1.0_r8 + 1.0e-15_r8)

  dqdt_soa_all(:,:,:,:) = 0.0_r8
  do k=top_lev,pver
     do i=1,ncol
        sum_uprt_soa(:) = 0.0_r8
        uptkrate_soa(:,:) = 0.0_r8
        do n = 1, nbins
           if (do_soag_any(n)) then  ! only bins that contain soa
              uptkratebb(n) = uptkrate(i,k,n)
              if (npoa .gt. 0) then
                 do j = 1, npoa
                    qold_poa(n,j) = poa_c(i,k,n,j)
                 end do
              else
                 qold_poa(n,j) = 0.0_r8
              end if
              do jsoa = 1, nsoa_vbs
                 ! 0.81 factor is for gas diffusivity (soa/h2so4)
                 ! (differences in fuch-sutugin and accom coef ignored)
                 fgain_soa(n,jsoa) = uptkratebb(n)*0.81_r8
                 uptkrate_soa(n,jsoa) = fgain_soa(n,jsoa)
                 sum_uprt_soa(jsoa) = sum_uprt_soa(jsoa) + fgain_soa(n,jsoa)
                 qold_soa(n,jsoa) = soa_vbs(i,k,n,jsoa)
              end do
           else
              qold_poa(n,:) = 0.0_r8
              qold_soa(n,:) = 0.0_r8
              fgain_soa(n,:) = 0.0_r8
           end if
        end do ! n

        do jsoa = 1, nsoa_vbs
           if (sum_uprt_soa(jsoa) > 0.0_r8) then
              do n = 1, nbins
                 if (do_soag_any(n)) then  ! only bins that contain soa
                    fgain_soa(n,jsoa) = fgain_soa(n,jsoa) / sum_uprt_soa(jsoa)
                 end if
              end do
           end if
        end do

!   uptake amount (fraction of gas uptaken) over deltat
        do jsoa = 1, nsoa_vbs
           avg_uprt_soa(jsoa) = (1.0_r8 - exp(-deltatxx*sum_uprt_soa(jsoa)))/deltatxx
        end do

!   sum_dqdt_soa = soa_a tendency from soa   gas uptake (mol/mol/s)

        do jsoa = 1, nsoa_vbs
           sum_dqdt_soa(jsoa) = q(i,k,l_soag(jsoa)) * avg_uprt_soa(jsoa)
        end do

        if (method_soa > 1) then
!   compute TMR tendencies for soag and soa interstial aerosol
!   using soa parameterization
           niter_max = 1000
           dqdt_soa_vbs(:,:) = 0.0_r8
           dqdt_soag(:) = 0.0_r8
           do jsoa = 1, nsoa_vbs
              qold_soag(jsoa) = q(i,k,l_soag(jsoa))
           end do

           call carma_aero_soaexch( deltat, t(i,k), pmid(i,k), &
                niter, niter_max, nbins, nsoa_vbs, npoa, &
                mw_poa_host, mw_soa_host, &
                qold_soag, qold_soa, qold_poa, uptkrate_soa, &
                dqdt_soag, dqdt_soa_vbs )

           sum_dqdt_soa(:) = -dqdt_soag(:)

        else if ( method_soa .eq. 1) then
!   compute TMR tendencies for soa interstial aerosol
!   due to simple gas uptake

           do n = 1, nbins
              if (do_soag_any(n) ) then
                 do jsoa = 1, nsoa_vbs
                    dqdt_soa_vbs(n,jsoa) = fgain_soa(n,jsoa)*sum_dqdt_soa(jsoa)
                 end do
              end if
           end do

        end if

        !  update soa to calcuate fractions (state variables and pbuf is not updated for SOA, will be done in CARMA)
        pdel_fac = pdel(i,k)/gravit
        qnew_soa(:) =0.0_r8
        qnew_soa_vbs(:,:) =0.0_r8

        do n = 1, nbins
           if ( do_soag_any(n) ) then
              if (nsoa.eq.nsoa_vbs) then
                 do jsoa = 1, nsoa_vbs
                    qsrflx(i,n,jsoa) = qsrflx(i,n,jsoa) + dqdt_soa_vbs(n,jsoa)*pdel_fac
                    dqdt_soa_all(i,k,n,jsoa) = dqdt_soa_vbs(n,jsoa) !  sum up for different volatility bins
                 end do
              else if (nsoa.eq.1) then
                 do jsoa = 1, nsoa_vbs
                    !  sum up for different volatility bins
                    dqdt_soa_all(i,k,n,nsoa) = dqdt_soa_all(i,k,n,nsoa) + dqdt_soa_vbs(n,jsoa)
                 end do
                 do jsoa = 1, nsoa_vbs
                    qsrflx(i,n,nsoa) = qsrflx(i,n,nsoa) + dqdt_soa_vbs(n,jsoa)*pdel_fac
                    qnew_soa_vbs(n,jsoa) = qold_soa(n,jsoa) + dqdt_soa_vbs(n,jsoa)*deltat
                    qnew_soa(n) = qnew_soa(n) + qnew_soa_vbs(n,jsoa) ! derive new fraction of SOA bin contributions
                 end do
                 do jsoa = 1, nsoa_vbs
                    if (qnew_soa(n) .gt. 0.0_r8) then
                       frac_vbs(i,k,n,jsoa) = qnew_soa_vbs(n,jsoa) / qnew_soa(n)
                    end if
                 end do
              else
                 call endrun( 'carma_aero_gasaerexch_sub error' )
              end if

!------- Add code for condensation/evaporation diagnostics sum of all bin---
              do jsoa = 1, nsoa_vbs
                 if (dqdt_soa_vbs(n,jsoa).ge.0.0_r8) then
                    qcon_vbs(i,k,jsoa)=qcon_vbs(i,k,jsoa) + dqdt_soa_vbs(n,jsoa)*(mw_soa/mwdry)
                    qcon(i,k)=qcon(i,k)+dqdt_soa_vbs(n,jsoa)*(mw_soa/mwdry)
                 else if (dqdt_soa_vbs(n,jsoa).lt.0.0_r8) then
                    qevap_vbs(i,k,jsoa)=qevap_vbs(i,k,jsoa) + dqdt_soa_vbs(n,jsoa)*(mw_soa/mwdry)
                    qevap(i,k)=qevap(i,k)+dqdt_soa_vbs(n,jsoa)*(mw_soa/mwdry)
                 endif
              end do
!---------------------------------------------------------------------------------------------------------------------
           end if
        end do ! n

!   compute TMR tendencies for SAOG gas
!   due to simple gas uptake
        do jsoa = 1, nsoa
           dqdt(i,k,l_soag(jsoa)) = -sum_dqdt_soa(jsoa)
        end do

     end do   ! "i = 1, ncol"
  end do     ! "k = top_lev, pver"

!  This applies dqdt tendencies for SOAG only , soa is done in CARMA
!  apply the dqdt to update q
!
  do jsoa = 1, nsoa_vbs
     do k = top_lev, pver
        do i = 1, ncol
           q(i,k,l_soag(jsoa)) = max (q(i,k,l_soag(jsoa)) + dqdt(i,k,l_soag(jsoa))*deltat, 1.0e-40_r8)
        end do
     end do
  end do


  !-----Outfld for condensation/evaporation------------------------------
  call outfld(trim('qcon_gaex'), qcon(:,:), pcols, lchnk )
  call outfld(trim('qevap_gaex'), qevap(:,:), pcols, lchnk )
  do jsoa = 1, nsoa_vbs
     write (outsoa, "(I2.2)") jsoa
     call outfld(trim('qcon_gaex')//outsoa, qcon_vbs(:,:,jsoa), pcols, lchnk )
     call outfld(trim('qevap_gaex')//outsoa, qevap_vbs(:,:,jsoa), pcols, lchnk )
  end do
  !-----------------------------------------------------------------------
  !   do history file of column-tendency fields over SOA fields (as defined in CARMA) and set pointer
  do m = 1, nbins
     if (do_soag_any(m)) then
        j  = 0
        do l = 1, nspec(m)
           mm = bin_idx(m,l)
           call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
           if (trim(spectype) == 's-organic') then
              j = j + 1
              fieldname= trim(fldname(mm)) // '_sfgaex1'
              do i = 1, ncol
                 qsrflx(i,m,j) = qsrflx(i,m,j)*(mw_soa/mwdry)
              end do
              call outfld( fieldname, qsrflx(:,m,j), pcols, lchnk )

              !set pointer field
              call pbuf_get_field(pbuf, dqdtsoa_idx(m,j),  dqdt_soa )

              dqdt_soa(:ncol,:) = dqdt_soa_all(:ncol,:,m,j) *(mw_soa/mbar(:ncol,:))
           end if
        end do ! l = ...
     end if
  end do ! m = ...

end subroutine carma_aero_gasaerexch_sub


!----------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine gas_aer_uptkrates( ncol,       loffset,                &
                              num_bin,          t,          pmid,       &
                              wetr,                     uptkrate    )

!
!                         /
!   computes   uptkrate = | dx  dN/dx  gas_conden_rate(Dp(x))
!                         /
!   using Gauss-Hermite quadrature of order nghq=2
!
!       Dp = particle diameter (cm)
!       x = ln(Dp)
!       dN/dx = log-normal particle number density distribution
!       gas_conden_rate(Dp) = 2 * pi * gasdiffus * Dp * F(Kn,ac)
!           F(Kn,ac) = Fuchs-Sutugin correction factor
!           Kn = Knudsen number
!           ac = accomodation coefficient
!

  integer,  intent(in) :: ncol                 ! number of atmospheric column
  integer,  intent(in) :: loffset
  real(r8), intent(in) :: t(pcols,pver)        ! Temperature in Kelvin
  real(r8), intent(in) :: pmid(pcols,pver)     ! Air pressure in Pa
  real(r8), intent(in) :: wetr(pcols,pver,nbins)
  real(r8), intent(in) :: num_bin(pcols,pver,nbins)

  real(r8), intent(out) :: uptkrate(pcols,pver,nbins)
  ! gas-to-aerosol mass transfer rates (1/s)


! local
  integer, parameter :: nghq = 2
  integer :: i, k, n

  ! Can use sqrt here once Lahey is gone.
  real(r8), parameter :: tworootpi = 3.5449077_r8
  real(r8), parameter :: root2 = 1.4142135_r8
  real(r8), parameter :: beta = 2.0_r8

  real(r8) :: const
  real(r8) :: dp
  real(r8) :: gasdiffus, gasspeed
  real(r8) :: freepathx2, fuchs_sutugin
  real(r8) :: knudsen

  ! initialize to zero
  uptkrate(:,:,:) = 0.0_r8

! outermost loop over all bins
  do n = 1, nbins

! loops k and i
     do k=top_lev,pver
        do i=1,ncol
           if (wetr(i,k,n) .gt. 0.0_r8) then

!   gasdiffus = h2so4 gas diffusivity from mosaic code (m^2/s)
!               (pmid must be Pa)
              gasdiffus = 0.557e-4_r8 * (t(i,k)**1.75_r8) / pmid(i,k)
!   gasspeed = h2so4 gas mean molecular speed from mosaic code (m/s)
              gasspeed  = 1.470e1_r8 * sqrt(t(i,k))
!   freepathx2 = 2 * (h2so4 mean free path)  (m)
              freepathx2 = 6.0_r8*gasdiffus/gasspeed
              dp = wetr(i,k,n) * 1.e-2 ! meters
              const = tworootpi * num_bin(i,k,n) * 2.0_r8 * dp
              ! gas_conden_rate(Dp) = const *  gasdiffus *  F(Kn,ac)
              !   knudsen number
              knudsen = freepathx2/dp
              fuchs_sutugin = (0.4875_r8*(1._r8 + knudsen)) /   &
                   (knudsen*(1.184_r8 + knudsen) + 0.4875_r8)
              uptkrate(i,k,n) = const * gasdiffus * fuchs_sutugin

           else
              uptkrate(i,k,n) = 0.0_r8
           end if

        end do   ! "do i = 1, ncol"
     end do   ! "do k = 1, pver"

  end do   ! "do n = 1, nbins"

end subroutine gas_aer_uptkrates

!----------------------------------------------------------------------
subroutine carma_aero_soaexch( dtfull, temp, pres, &
     niter, niter_max, nbins, ntot_soaspec, ntot_poaspec,  &
     mw_poa_host, mw_soa_host, &
     g_soa_in, a_soa_in, a_poa_in, xferrate_in, &
     g_soa_tend, a_soa_tend )

!-----------------------------------------------------------------------
!
! Purpose:
!
! calculates condensation/evaporation of "soa gas"
! to/from multiple aerosol modes in 1 grid cell
!
! key assumptions
! (1) ambient equilibrium vapor pressure of soa gas
!     is given by p0_soa_298 and delh_vap_soa
! (2) equilibrium vapor pressure of soa gas at aerosol
!     particle surface is given by raoults law in the form
!     g_star = g0_soa*[a_soa/(a_soa + a_opoa)]
! (3) (oxidized poa)/(total poa) is equal to frac_opoa (constant)
!
!
! Author: R. Easter and R. Zaveri
! Additions to run with multiple BC, SOA and POM's: Shrivastava et al., 2015
!-----------------------------------------------------------------------

  use mo_constants, only: rgas ! Gas constant (J/K/mol)

  real(r8), intent(in)  :: dtfull                     ! full integration time step (s)
  real(r8), intent(in)  :: temp                       ! air temperature (K)
  real(r8), intent(in)  :: pres                       ! air pressure (Pa)
  integer,  intent(out) :: niter                      ! number of iterations performed
  integer,  intent(in)  :: niter_max                  ! max allowed number of iterations
  integer,  intent(in)  :: nbins                      ! number of bins
  integer,  intent(in)  :: ntot_poaspec               ! number of poa species
  integer,  intent(in)  :: ntot_soaspec               ! number of soa species
  real(r8), intent(in)  :: mw_poa_host                ! molec wght of poa used in host code
  real(r8), intent(in)  :: mw_soa_host                ! molec wght of soa used in host code
  real(r8), intent(in)  :: g_soa_in(ntot_soaspec)               ! initial soa gas mixrat (mol/mol at host mw)
  real(r8), intent(in)  :: a_soa_in(nbins,ntot_soaspec)    ! initial soa aerosol mixrat (mol/mol at host mw)
  real(r8), intent(in)  :: a_poa_in(nbins,ntot_poaspec)    ! initial poa aerosol mixrat (mol/mol at host mw)
  real(r8), intent(in)  :: xferrate_in(nbins,ntot_soaspec) ! gas-aerosol mass transfer rate (1/s)
  real(r8), intent(out) :: g_soa_tend(ntot_soaspec)             ! soa gas mixrat tendency (mol/mol/s at host mw)
  real(r8), intent(out) :: a_soa_tend(nbins,ntot_soaspec)  ! soa aerosol mixrat tendency (mol/mol/s at host mw)

  integer :: ll
  integer :: m

  logical :: skip_soamode(nbins)   ! true if this bin does not have soa

  real(r8), parameter :: a_min1 = 1.0e-40_r8
  real(r8), parameter :: g_min1 = 1.0e-40_r8
  real(r8), parameter :: alpha = 0.05_r8     ! parameter used in calc of time step
  real(r8), parameter :: dtsub_fixed = -1.0_r8  ! fixed sub-step for time integration (s)

  real(r8) :: a_ooa_sum_tmp(nbins)          ! total ooa (=soa+opoa) in a bin
  real(r8) :: a_opoa(nbins)                 ! oxidized-poa aerosol mixrat (mol/mol at actual mw)
  real(r8) :: a_soa(nbins,ntot_soaspec)     ! soa aerosol mixrat (mol/mol at actual mw)
  real(r8) :: a_soa_tmp(nbins,ntot_soaspec) ! temporary soa aerosol mixrat (mol/mol)
  real(r8) :: beta(nbins,ntot_soaspec)      ! dtcur*xferrate
  real(r8) :: delh_vap_soa(ntot_soaspec)           ! delh_vap_soa = heat of vaporization for gas soa (J/mol)
  real(r8) :: del_g_soa_tmp(ntot_soaspec)
  real(r8) :: dtcur                                ! current time step (s)
  real(r8) :: dtmax                                ! = (dtfull-tcur)
  real(r8) :: g0_soa(ntot_soaspec)                 ! ambient soa gas equilib mixrat (mol/mol at actual mw)
  real(r8) :: g_soa(ntot_soaspec)                  ! soa gas mixrat (mol/mol at actual mw)
  real(r8) :: g_star(nbins,ntot_soaspec)    ! soa gas mixrat that is in equilib
  ! with each aerosol mode (mol/mol)
  real(r8) :: mw_poa                               ! actual molec wght of poa
  real(r8) :: mw_soa                               ! actual molec wght of soa
  real(r8) :: opoa_frac(ntot_poaspec)              ! fraction of poa that is opoa
  real(r8) :: phi(nbins,ntot_soaspec)       ! "relative driving force"
  real(r8) :: p0_soa(ntot_soaspec)                 ! soa gas equilib vapor presssure (atm)
  real(r8) :: p0_soa_298(ntot_soaspec)             ! p0_soa_298 = soa gas equilib vapor presssure (atm) at 298 k
  real(r8) :: sat(nbins,ntot_soaspec)       ! sat(m,ll) = g0_soa(ll)/a_ooa_sum_tmp(m) = g_star(m,ll)/a_soa(m,ll)
  !    used by the numerical integration scheme -- it is not a saturation rato!
  real(r8) :: tcur                                 ! current integration time (from 0 s)
  real(r8) :: tmpa, tmpb, tmpf
  real(r8) :: tot_soa(ntot_soaspec)                ! g_soa + sum( a_soa(:) )
  real(r8) :: xferrate(nbins,ntot_soaspec)    ! gas-aerosol mass transfer rate (1/s)

! Changed by Manish Shrivastava
  opoa_frac(:) = 0.0_r8 !POA does not form solution with SOA for all runs; set opoa_frac=0.0_r8  by Manish Shrivastava
  mw_poa = 250.0_r8
  mw_soa = 250.0_r8

  ! New SOA properties added by Manish Shrivastava on 09/27/2012
  if (ntot_soaspec ==1) then
     p0_soa_298(:) = 1.0e-12_r8
     delh_vap_soa(:) = 156.0e3_r8
     opoa_frac(:) = 0.0_r8
  elseif (ntot_soaspec ==2) then
     ! same for anthropogenic and biomass burning species
     p0_soa_298 (1) = 1.0e-10_r8
     p0_soa_298 (2) = 1.0e-10_r8
     delh_vap_soa(:) = 156.0e3_r8
  elseif(ntot_soaspec ==5) then
     ! 5 volatility bins for each of the a combined SOA classes ( including biomass burning, fossil fuel, biogenic)
     p0_soa_298 (1) = 9.7831E-13_r8 !soaff0 C*=0.01ug/m3
     p0_soa_298 (2) = 9.7831E-12_r8 !soaff1 C*=0.10ug/m3
     p0_soa_298 (3) = 9.7831E-11_r8 !soaff2 C*=1.0ug/m3
     p0_soa_298 (4) = 9.7831E-10_r8 !soaff3 C*=10.0ug/m3
     p0_soa_298 (5) = 9.7831E-9_r8  !soaff4 C*=100.0ug/m3

     delh_vap_soa(1) = 153.0e3_r8
     delh_vap_soa(2) = 142.0e3_r8
     delh_vap_soa(3) = 131.0e3_r8
     delh_vap_soa(4) = 120.0e3_r8
     delh_vap_soa(5) = 109.0e3_r8
  elseif(ntot_soaspec ==15) then
     !
     ! 5 volatility bins for each of the 3 SOA classes ( biomass burning, fossil fuel, biogenic)
     ! SOA species 1-5 are for anthropogenic while 6-10 are for biomass burning SOA
     ! SOA species 11-15 are for biogenic SOA, based on Cappa et al., Reference needs to be updated
     ! For MW=250.0
     p0_soa_298 (1) = 9.7831E-13_r8 !soaff0 C*=0.01ug/m3
     p0_soa_298 (2) = 9.7831E-12_r8 !soaff1 C*=0.10ug/m3
     p0_soa_298 (3) = 9.7831E-11_r8 !soaff2 C*=1.0ug/m3
     p0_soa_298 (4) = 9.7831E-10_r8 !soaff3 C*=10.0ug/m3
     p0_soa_298 (5) = 9.7831E-9_r8  !soaff4 C*=100.0ug/m3
     p0_soa_298 (6) = 9.7831E-13_r8 !soabb0 C*=0.01ug/m3
     p0_soa_298 (7) = 9.7831E-12_r8 !soabb1 C*=0.10ug/m3
     p0_soa_298 (8) = 9.7831E-11_r8 !soabb2 C*=1.0ug/m3
     p0_soa_298 (9) = 9.7831E-10_r8 !soabb3 C*=10.0ug/m3
     p0_soa_298 (10) = 9.7831E-9_r8  !soabb4 C*=100.0ug/m3
     p0_soa_298 (11) = 9.7831E-13_r8 !soabg0 C*=0.01ug/m3
     p0_soa_298 (12) = 9.7831E-12_r8 !soabg1 C*=0.1ug/m3
     p0_soa_298 (13) = 9.7831E-11_r8 !soabg2 C*=1.0ug/m3
     p0_soa_298 (14) = 9.7831E-10_r8 !soabg3 C*=10.0ug/m3
     p0_soa_298 (15) = 9.7831E-9_r8  !soabg4 C*=100.0ug/m3

     !
     ! have to be adjusted to 15 species, following the numbers by Epstein et al., 2012
     !
     delh_vap_soa(1) = 153.0e3_r8
     delh_vap_soa(2) = 142.0e3_r8
     delh_vap_soa(3) = 131.0e3_r8
     delh_vap_soa(4) = 120.0e3_r8
     delh_vap_soa(5) = 109.0e3_r8
     delh_vap_soa(6) = 153.0e3_r8
     delh_vap_soa(7) = 142.0e3_r8
     delh_vap_soa(8) = 131.0e3_r8
     delh_vap_soa(9) = 120.0e3_r8
     delh_vap_soa(10) = 109.0e3_r8
     delh_vap_soa(11) = 153.0e3_r8
     delh_vap_soa(12) = 142.0e3_r8
     delh_vap_soa(13) = 131.0e3_r8
     delh_vap_soa(14) = 120.0e3_r8
     delh_vap_soa(15) = 109.0e3_r8
  endif

  !BSINGH - Initialized g_soa_tend and a_soa_tend to circumvent the undefined behavior (04/16/12)
  g_soa_tend(:)   = 0.0_r8
  a_soa_tend(:,:) = 0.0_r8
  xferrate(:,:) = 0.0_r8

  ! determine which modes have non-zero transfer rates
  !    and are involved in the soa gas-aerosol transfer
  ! for diameter = 1 nm and number = 1 #/cm3, xferrate ~= 1e-9 s-1
  do m = 1, nbins
     if (do_soag_any(m)) then
        skip_soamode(m) = .false.
        do ll = 1, ntot_soaspec
           xferrate(m,ll) = xferrate_in(m,ll)
        end do
     else
        skip_soamode(m) = .true.
     end if
  end do

  ! convert incoming mixing ratios from mol/mol at the "host-code" molec. weight (12.0 in cam5)
  !    to mol/mol at the "actual" molec. weight (currently assumed to be 250.0)
  ! also
  !    force things to be non-negative
  !    calc tot_soa(ll)
  !    calc a_opoa (always slightly >0)
  do ll = 1, ntot_soaspec
     tmpf = mw_soa_host/mw_soa
     g_soa(ll) = max( g_soa_in(ll), 0.0_r8 ) * tmpf
     tot_soa(ll) = g_soa(ll)
     do m = 1, nbins
        if ( skip_soamode(m) ) cycle
        a_soa(m,ll) = max( a_soa_in(m,ll), 0.0_r8 ) * tmpf
        tot_soa(ll) = tot_soa(ll) + a_soa(m,ll)
     end do
  end do


  tmpf = mw_poa_host/mw_poa
  do m = 1, nbins
     if ( skip_soamode(m) ) cycle
     a_opoa(m) = 0.0_r8
     !check since it seems like in the modal approach there is a bug, not summing up the values for each specie
     do ll = 1, ntot_poaspec
        tmpf = mw_poa_host/mw_poa
        a_opoa(m) = a_opoa(m) + opoa_frac(ll)*a_poa_in(m,ll)
        a_opoa(m) = max( a_opoa(m), 1.0e-40_r8 )  ! force to small non-zero value
     end do
  end do

  ! calc ambient equilibrium soa gas
  do ll = 1, ntot_soaspec
     p0_soa(ll) = p0_soa_298(ll) * &
          exp( -(delh_vap_soa(ll)/rgas)*((1.0_r8/temp)-(1.0_r8/298.0_r8)) )
     g0_soa(ll) = 1.01325e5_r8*p0_soa(ll)/pres
  end do

  ! IF mw of soa EQ 12 (as in the MAM3 default case), this has to be in
  ! should actully talk the mw from the chemistry mechanism and substitute with 12.0

  niter = 0
  tcur = 0.0_r8
  dtcur = 0.0_r8
  phi(:,:) = 0.0_r8
  g_star(:,:) = 0.0_r8

! integration loop -- does multiple substeps to reach dtfull
  time_loop: do while (tcur < dtfull-1.0e-3_r8 )

     niter = niter + 1
     if (niter > niter_max) exit

     tmpa = 0.0_r8  ! time integration parameter for all soa species
     do m = 1, nbins
        if ( skip_soamode(m) ) cycle
        a_ooa_sum_tmp(m) = sum( a_soa(m,1:ntot_soaspec) )
     end do
     do ll = 1, ntot_soaspec
        tmpb = 0.0_r8  ! time integration parameter for a single soa species
        do m = 1, nbins
           if ( skip_soamode(m) ) cycle
           sat(m,ll) = g0_soa(ll)/max( a_ooa_sum_tmp(m), a_min1 )
           g_star(m,ll) = sat(m,ll)*a_soa(m,ll)
           phi(m,ll) = (g_soa(ll) - g_star(m,ll))/max( g_soa(ll), g_star(m,ll), g_min1 )
           tmpb = tmpb + xferrate(m,ll)*abs(phi(m,ll))
        end do
        tmpa = max( tmpa, tmpb )
     end do

     if (dtsub_fixed > 0.0_r8) then
        dtcur = dtsub_fixed
        tcur = tcur + dtcur
     else
        dtmax = dtfull-tcur
        if (dtmax*tmpa <= alpha) then
! here alpha/tmpa >= dtmax, so this is final substep
           dtcur = dtmax
           tcur = dtfull
        else
           dtcur = alpha/tmpa
           tcur = tcur + dtcur
        end if
     end if

! step 1 - for modes where soa is condensing, estimate "new" a_soa(m,ll)
!    using an explicit calculation with "old" g_soa
!    and g_star(m,ll) calculated using "old" a_soa(m,ll)
! do this to get better estimate of "new" a_soa(m,ll) and sat(m,ll)
     do m = 1, nbins
        if ( skip_soamode(m) ) cycle
        do ll = 1, ntot_soaspec
           ! first ll loop calcs a_soa_tmp(m,ll) & a_ooa_sum_tmp
           a_soa_tmp(m,ll) = a_soa(m,ll)
           beta(m,ll) = dtcur*xferrate(m,ll)
           del_g_soa_tmp(ll) = g_soa(ll) - g_star(m,ll)
           if (del_g_soa_tmp(ll) > 0.0_r8) then
              a_soa_tmp(m,ll) = a_soa(m,ll) + beta(m,ll)*del_g_soa_tmp(ll)
           end if
        end do
        a_ooa_sum_tmp(m) =  sum( a_soa_tmp(m,1:ntot_soaspec) )
        do ll = 1, ntot_soaspec
           ! second ll loop calcs sat & g_star
           if (del_g_soa_tmp(ll) > 0.0_r8) then
              sat(m,ll) = g0_soa(ll)/max( a_ooa_sum_tmp(m), a_min1 )
              g_star(m,ll) = sat(m,ll)*a_soa_tmp(m,ll)   ! this just needed for diagnostics
           end if
        end do
     end do

! step 2 - implicit in g_soa and semi-implicit in a_soa,
!    with g_star(m,ll) calculated semi-implicitly
     do ll = 1, ntot_soaspec
        tmpa = 0.0_r8
        tmpb = 0.0_r8
        do m = 1, nbins
           if ( skip_soamode(m) ) cycle
           tmpa = tmpa + a_soa(m,ll)/(1.0_r8 + beta(m,ll)*sat(m,ll))
           tmpb = tmpb + beta(m,ll)/(1.0_r8 + beta(m,ll)*sat(m,ll))
        end do

        g_soa(ll) = (tot_soa(ll) - tmpa)/(1.0_r8 + tmpb)
        g_soa(ll) = max( 0.0_r8, g_soa(ll) )
        do m = 1, nbins
           if ( skip_soamode(m) ) cycle
           a_soa(m,ll) = (a_soa(m,ll) + beta(m,ll)*g_soa(ll))/   &
                (1.0_r8 + beta(m,ll)*sat(m,ll))
        end do
     end do

  end do time_loop

! calculate outgoing tendencies (at the host-code molec. weight)
! (a_soa & g_soa are at actual mw, but a_soa_in & g_soa_in are at host-code mw)
  do ll = 1, ntot_soaspec
     tmpf = mw_soa/mw_soa_host
     g_soa_tend(ll) = (g_soa(ll)*tmpf - g_soa_in(ll))/dtfull
     do m = 1, nbins
        if ( skip_soamode(m) ) cycle
        a_soa_tend(m,ll) = (a_soa(m,ll)*tmpf - a_soa_in(m,ll))/dtfull
     end do
  end do

end subroutine carma_aero_soaexch

!----------------------------------------------------------------------

end module carma_aero_gasaerexch
