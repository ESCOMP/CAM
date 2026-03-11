!-------------------------------------------------------------------------------
! Portable aerosol optics core:
!  Creates the aerosol_optics object
!  Calculates per-bin SW/LW aerosol optics.
!-------------------------------------------------------------------------------
module aerosol_optics_core
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: create_aerosol_optics_object
  public :: aerosol_optics_sw_bin
  public :: aerosol_optics_lw_bin

  ! Jasper Kok et al. (2017) Fig. 1d: 20-60 % higher mass extinction efficiency
  ! because dust is aspherical. Currently not captured by the spherical assumption
  ! in the optical calculation. Asphericity is strong for D > 1 um (coarse mode).
  real(r8), parameter, public :: dustaspherical_opts = 1.3_r8

contains

  !===============================================================================
  ! Dispatch to the appropriate concrete aerosol_optics constructor based on
  ! the opticstype string from aeroprops%optics_params.
  ! Returns a null pointer for unrecognized opticstype (caller handles error).
  !===============================================================================
  function create_aerosol_optics_object(aeroprops, aerostate, ibin, &
                                        ncol, nlev, nswbands, nlwbands, numrh, &
                                        relh, sulfwtpct, crefwsw, crefwlw, &
                                        geometric_radius) result(aero_optics)

    use phys_prop, only: ot_length

    use aerosol_properties_mod, only: aerosol_properties
    use aerosol_state_mod, only: aerosol_state
    use aerosol_optics_mod, only: aerosol_optics

    use refractive_aerosol_optics_mod, only: refractive_aerosol_optics
    use hygrocoreshell_aerosol_optics_mod, only: hygrocoreshell_aerosol_optics
    use hygrowghtpct_aerosol_optics_mod, only: hygrowghtpct_aerosol_optics
    use hygro_aerosol_optics_mod, only: hygro_aerosol_optics
    use hygroscopic_aerosol_optics_mod, only: hygroscopic_aerosol_optics
    use insoluble_aerosol_optics_mod, only: insoluble_aerosol_optics
    use volcrad_aerosol_optics_mod, only: volcrad_aerosol_optics

    class(aerosol_properties), intent(in), target :: aeroprops
    class(aerosol_state),      intent(in), target :: aerostate
    integer,     intent(in)                   :: ibin
    integer,     intent(in)                   :: ncol
    integer,     intent(in)                   :: nlev
    integer,     intent(in)                   :: nswbands
    integer,     intent(in)                   :: nlwbands
    integer,     intent(in)                   :: numrh
    real(r8),    intent(in)                   :: relh(:, :)
    real(r8),    intent(in)                   :: sulfwtpct(:, :)
    complex(r8), intent(in)                   :: crefwsw(:)
    complex(r8), intent(in)                   :: crefwlw(:)
    real(r8),    intent(in), optional, target :: geometric_radius(:, :)

    class(aerosol_optics), pointer :: aero_optics

    character(len=ot_length) :: opticstype

    nullify (aero_optics)

    call aeroprops%optics_params(bin_ndx=ibin, opticstype=opticstype)

    select case (trim(opticstype))
    case ('modal') ! refractive method
      aero_optics => refractive_aerosol_optics(aeroprops, aerostate, ibin, &
                                               ncol, nlev, nswbands, nlwbands, crefwsw, crefwlw)
    case ('hygroscopic_coreshell')
      aero_optics => hygrocoreshell_aerosol_optics(aeroprops, aerostate, &
                                                   ibin, ncol, nlev, relh)
    case ('hygroscopic_wtp')
      aero_optics => hygrowghtpct_aerosol_optics(aeroprops, aerostate, &
                                                 ibin, ncol, nlev, sulfwtpct)
    case ('hygro')
      ! Short-wave hygroscopic aerosol, Long-wave non-hygroscopic
      ! aerosol optical properties
      aero_optics => hygro_aerosol_optics(aeroprops, aerostate, &
                                          ibin, ncol, nlev, numrh, relh)
    case ('hygroscopic')
      ! Short-wave and long-wave hygroscopic aerosol properties
      aero_optics => hygroscopic_aerosol_optics(aeroprops, aerostate, ibin, &
                                                ncol, nlev, numrh, relh)
    case ('nonhygro', 'insoluble')
      aero_optics => insoluble_aerosol_optics(aeroprops, aerostate, ibin)

    case ('volcanic_radius', 'volcanic_radius1', 'volcanic_radius2', 'volcanic_radius3', 'volcanic_radius5')
      if (present(geometric_radius)) then
        aero_optics => volcrad_aerosol_optics(aeroprops, aerostate, &
                                              ibin, ncol, nlev, geometric_radius)
      end if

    end select

  end function create_aerosol_optics_object

  !===============================================================================
  ! Per-bin SW aerosol optics including dust asphericity correction.
  !
  ! Returns per-bin extinction optical depth (tau_bin), single-scatter albedo
  ! (ssa_bin), and asymmetry parameter (asm_bin). For coarse dust modes,
  ! tau_bin at the visible band (idx_sw_diag) is modified by the asphericity
  ! correction (1.3x for dust-attributed AOD).
  !===============================================================================
  subroutine aerosol_optics_sw_bin(aeroprops, aerostate, ibin, &
                                   ncol, nlev, top_lev, nswbands, nlwbands, numrh, &
                                   idx_sw_diag, &
                                   relh, sulfwtpct, mass, crefwsw, crefwlw, &
                                   geometric_radius, &
                                   tau_bin, ssa_bin, asm_bin, &
                                   pabs_vis, dopaer0_vis, &
                                   errmsg, errflg)
    use aerosol_properties_mod, only: aerosol_properties
    use aerosol_properties_mod, only: aero_name_len
    use aerosol_state_mod,      only: aerosol_state
    use aerosol_optics_mod,     only: aerosol_optics

    class(aerosol_properties), intent(in), target :: aeroprops
    class(aerosol_state),      intent(in), target :: aerostate
    integer,     intent(in) :: ibin
    integer,     intent(in) :: ncol, nlev, top_lev
    integer,     intent(in) :: nswbands
    integer,     intent(in) :: nlwbands
    integer,     intent(in) :: numrh
    integer,     intent(in) :: idx_sw_diag
    real(r8),    intent(in) :: relh(:, :)
    real(r8),    intent(in) :: sulfwtpct(:, :)
    real(r8),    intent(in) :: mass(:, :)                    ! layer mass (pdeldry*rga)
    complex(r8), intent(in) :: crefwsw(:)
    complex(r8), intent(in) :: crefwlw(:)
    real(r8),    intent(in), optional, target :: geometric_radius(:, :)

    real(r8),    intent(out) :: tau_bin(:, :, :)  ! (ncol,nlev,nswbands) extinction OD
    real(r8),    intent(out) :: ssa_bin(:, :, :)  ! (ncol,nlev,nswbands) single scatter albedo
    real(r8),    intent(out) :: asm_bin(:, :, :)  ! (ncol,nlev,nswbands) asymmetry parameter

    ! Diagnostic outputs for BFB absorption diagnostics in CAM
    real(r8),    intent(out) :: pabs_vis(:, :)     ! (ncol,nlev) specific absorption at vis band
    real(r8),    intent(out) :: dopaer0_vis(:, :)  ! (ncol,nlev) pre-asphericity tau at vis band

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    class(aerosol_optics), pointer :: aero_optics
    real(r8) :: pext(ncol), pabs(ncol), palb(ncol), pasm(ncol)
    integer :: iwav, ilev, icol, ispec

    ! For asphericity computation
    logical :: coarse_dust_mode
    character(len=aero_name_len) :: modetype
    real(r8) :: wetvol(ncol, nlev), watervol(ncol, nlev)
    real(r8) :: vol(ncol)
    real(r8) :: scatdust(ncol), absdust(ncol), hygrodust(ncol)
    real(r8) :: scatbc(ncol), absbc(ncol), hygrobc(ncol)
    real(r8) :: scatpom(ncol), abspom(ncol), hygropom(ncol)
    real(r8) :: scatsoa(ncol), abssoa(ncol), hygrosoa(ncol)
    real(r8) :: scatsulf(ncol), abssulf(ncol), hygrosulf(ncol)
    real(r8) :: scatsslt(ncol), abssslt(ncol), hygrosslt(ncol)
    real(r8) :: scath2o, absh2o, sumscat, sumabs, sumhygro
    real(r8) :: aodc, dustaod0
    real(r8) :: specdens
    character(len=32) :: spectype
    real(r8), pointer :: specmmr(:, :)
    real(r8) :: hygro_aer
    complex(r8), pointer :: specrefindex(:)

    errflg = 0
    errmsg = ''

    tau_bin = 0._r8
    ssa_bin = 0._r8
    asm_bin = 0._r8
    pabs_vis = 0._r8
    dopaer0_vis = 0._r8

    ! Create aerosol optics object
    aero_optics => create_aerosol_optics_object(aeroprops, aerostate, ibin, &
                                                ncol, nlev, nswbands, nlwbands, numrh, &
                                                relh, sulfwtpct, crefwsw, crefwlw, &
                                                geometric_radius)

    if (.not. associated(aero_optics)) then
      errflg = 1
      errmsg = 'unrecognized aerosol optics type, could not create object'
    end if

    ! Determine if this is a coarse dust mode (MAM only)
    coarse_dust_mode = .false.
    if (aeroprops%model_is('MAM')) then
      modetype = aeroprops%bin_name(bin_ndx=ibin)
      coarse_dust_mode = (modetype == 'coarse' .or. modetype == 'coarse_dust')
    end if

    ! Main optics loop
    do iwav = 1, nswbands
      do ilev = top_lev, nlev
        call aero_optics%sw_props(ncol, ilev, iwav, pext, pabs, palb, pasm)

        do icol = 1, ncol
          tau_bin(icol, ilev, iwav) = pext(icol)*mass(icol, ilev)
          ssa_bin(icol, ilev, iwav) = palb(icol)
          asm_bin(icol, ilev, iwav) = pasm(icol)
        end do

        ! Save specific absorption at visible band for BFB diagnostics
        if (iwav == idx_sw_diag) then
          pabs_vis(1:ncol, ilev) = pabs(1:ncol)
        end if
      end do
    end do

    ! Save pre-asphericity tau at visible band for BFB diagnostics
    if (idx_sw_diag > 0) then
      dopaer0_vis(1:ncol, top_lev:nlev) = tau_bin(1:ncol, top_lev:nlev, idx_sw_diag)
    end if

    ! Apply asphericity correction at visible band for coarse dust
    ! dmleung 20 Oct 2025: coarse-mode dust is aspherical, with ~30 % enhanced
    ! extinction compared with spherical coarse-mode dust.
    ! ref: Fig. 1d of Jasper F. Kok et al. (2017)
    if (coarse_dust_mode .and. idx_sw_diag > 0) then
      wetvol(:ncol, :nlev) = aerostate%wet_volume(aeroprops, ibin, ncol, nlev)
      watervol(:ncol, :nlev) = aerostate%water_volume(aeroprops, ibin, ncol, nlev)

      do ilev = top_lev, nlev
        scatdust(:ncol) = 0._r8
        absdust(:ncol) = 0._r8
        hygrodust(:ncol) = 0._r8
        scatsulf(:ncol) = 0._r8
        abssulf(:ncol) = 0._r8
        hygrosulf(:ncol) = 0._r8
        scatbc(:ncol) = 0._r8
        absbc(:ncol) = 0._r8
        hygrobc(:ncol) = 0._r8
        scatpom(:ncol) = 0._r8
        abspom(:ncol) = 0._r8
        hygropom(:ncol) = 0._r8
        scatsoa(:ncol) = 0._r8
        abssoa(:ncol) = 0._r8
        hygrosoa(:ncol) = 0._r8
        scatsslt(:ncol) = 0._r8
        abssslt(:ncol) = 0._r8
        hygrosslt(:ncol) = 0._r8

        ! loop over species ...
        do ispec = 1, aeroprops%nspecies(ibin)
          call aeroprops%get(ibin, ispec, density=specdens, &
                             spectype=spectype, refindex_sw=specrefindex, hygro=hygro_aer)
          call aerostate%get_ambient_mmr(species_ndx=ispec, bin_ndx=ibin, mmr=specmmr)

          do icol = 1, ncol
            vol(icol) = specmmr(icol, ilev)/specdens

            select case (trim(spectype))
            case ('dust')
              if (associated(specrefindex)) then
                scatdust(icol) = vol(icol)*specrefindex(idx_sw_diag)%re
                absdust(icol) = -vol(icol)*specrefindex(idx_sw_diag)%im
              end if
              hygrodust(icol) = vol(icol)*hygro_aer
            case ('black-c')
              if (associated(specrefindex)) then
                scatbc(icol) = vol(icol)*specrefindex(idx_sw_diag)%re
                absbc(icol) = -vol(icol)*specrefindex(idx_sw_diag)%im
              end if
              hygrobc(icol) = vol(icol)*hygro_aer
            case ('sulfate')
              if (associated(specrefindex)) then
                scatsulf(icol) = vol(icol)*specrefindex(idx_sw_diag)%re
                abssulf(icol) = -vol(icol)*specrefindex(idx_sw_diag)%im
              end if
              hygrosulf(icol) = vol(icol)*hygro_aer
            case ('p-organic')
              if (associated(specrefindex)) then
                scatpom(icol) = vol(icol)*specrefindex(idx_sw_diag)%re
                abspom(icol) = -vol(icol)*specrefindex(idx_sw_diag)%im
              end if
              hygropom(icol) = vol(icol)*hygro_aer
            case ('s-organic')
              if (associated(specrefindex)) then
                scatsoa(icol) = vol(icol)*specrefindex(idx_sw_diag)%re
                abssoa(icol) = -vol(icol)*specrefindex(idx_sw_diag)%im
              end if
              hygrosoa(icol) = vol(icol)*hygro_aer
            case ('seasalt')
              if (associated(specrefindex)) then
                scatsslt(icol) = vol(icol)*specrefindex(idx_sw_diag)%re
                abssslt(icol) = -vol(icol)*specrefindex(idx_sw_diag)%im
              end if
              hygrosslt(icol) = vol(icol)*hygro_aer
            end select
          end do
        end do

        do icol = 1, ncol
          if (wetvol(icol, ilev) > 1.e-40_r8 .and. vol(icol) > 0._r8) then

            ! partition optical depth into contributions from each constituent
            ! assume contribution is proportional to refractive index X volume

            scath2o = watervol(icol, ilev)*crefwsw(idx_sw_diag)%re
            absh2o = -watervol(icol, ilev)*crefwsw(idx_sw_diag)%im
            sumscat = scatsulf(icol) + scatpom(icol) + scatsoa(icol) + scatbc(icol) + &
                      scatdust(icol) + scatsslt(icol) + scath2o
            sumabs = abssulf(icol) + abspom(icol) + abssoa(icol) + absbc(icol) + &
                     absdust(icol) + abssslt(icol) + absh2o
            sumhygro = hygrosulf(icol) + hygropom(icol) + hygrosoa(icol) + hygrobc(icol) + &
                       hygrodust(icol) + hygrosslt(icol)

            scatdust(icol) = (scatdust(icol) + scath2o*hygrodust(icol)/sumhygro)/sumscat
            absdust(icol) = (absdust(icol) + absh2o*hygrodust(icol)/sumhygro)/sumabs

            aodc = (absdust(icol)*(1.0_r8 - ssa_bin(icol, ilev, idx_sw_diag)) &
                    + ssa_bin(icol, ilev, idx_sw_diag)*scatdust(icol)) &
                      *tau_bin(icol, ilev, idx_sw_diag)

            ! dustaod0 is the single-level spherical dust AOD
            dustaod0 = aodc

            ! scale up dust AOD by 30 %
            tau_bin(icol, ilev, idx_sw_diag) = tau_bin(icol, ilev, idx_sw_diag) &
                                               - dustaod0 + dustaod0*dustaspherical_opts

          end if
        end do
      end do ! ilev
    end if ! if coarse_dust_mode && idx_sw_diag > 0

    deallocate (aero_optics)
  end subroutine aerosol_optics_sw_bin

  !===============================================================================
  ! Per-bin LW aerosol optics. Returns absorption optical depth (tau_lw_bin)
  ! and raw specific absorption (absorp_bin) for diagnostic use.
  !===============================================================================
  subroutine aerosol_optics_lw_bin(aeroprops, aerostate, ibin, &
                                   ncol, nlev, nswbands, nlwbands, numrh, &
                                   relh, sulfwtpct, mass, crefwsw, crefwlw, &
                                   geometric_radius, &
                                   tau_lw_bin, absorp_bin, &
                                   errmsg, errflg)
    use aerosol_properties_mod, only: aerosol_properties
    use aerosol_state_mod, only: aerosol_state
    use aerosol_optics_mod, only: aerosol_optics

    class(aerosol_properties), intent(in), target :: aeroprops
    class(aerosol_state),      intent(in), target :: aerostate
    integer,     intent(in)         :: ibin
    integer,     intent(in)         :: ncol, nlev
    integer,     intent(in)         :: nswbands
    integer,     intent(in)         :: nlwbands
    integer,     intent(in)         :: numrh
    real(r8),    intent(in)         :: relh(:, :)
    real(r8),    intent(in)         :: sulfwtpct(:, :)
    real(r8),    intent(in)         :: mass(:, :)
    complex(r8), intent(in)         :: crefwsw(:)
    complex(r8), intent(in)         :: crefwlw(:)
    real(r8),    intent(in), optional, target :: geometric_radius(:, :)
    real(r8),    intent(out)        :: tau_lw_bin(:, :, :)   ! (ncol,nlev,nlwbands) absorption OD
    real(r8),    intent(out)        :: absorp_bin(:, :, :)   ! (ncol,nlev,nlwbands) raw specific absorption (pabs)

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    class(aerosol_optics), pointer :: aero_optics
    real(r8) :: pabs(ncol)
    integer  :: iwav, ilev, icol

    errflg = 0
    errmsg = ''

    tau_lw_bin = 0._r8
    absorp_bin = 0._r8

    ! Create aerosol optics object
    aero_optics => create_aerosol_optics_object(aeroprops, aerostate, ibin, &
                                                ncol, nlev, nswbands, nlwbands, numrh, &
                                                relh, sulfwtpct, crefwsw, crefwlw, &
                                                geometric_radius)

    if (.not. associated(aero_optics)) then
      errflg = 1
      errmsg = 'unrecognized aerosol optics type, could not create object'
    end if

    do iwav = 1, nlwbands
      do ilev = 1, nlev
        call aero_optics%lw_props(ncol, ilev, iwav, pabs)

        do icol = 1, ncol
          tau_lw_bin(icol, ilev, iwav) = pabs(icol)*mass(icol, ilev)
          absorp_bin(icol, ilev, iwav) = pabs(icol)
        end do
      end do
    end do

    deallocate (aero_optics)

  end subroutine aerosol_optics_lw_bin

end module aerosol_optics_core
