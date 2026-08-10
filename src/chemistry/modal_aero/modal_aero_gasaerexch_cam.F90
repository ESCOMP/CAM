! CAM wrapper for modal_aero_gasaerexch.
! Resolves CAM-specific species indices (applying loffset) and mode metadata,
! calls the portable modal_aero_gasaerexch_init, and registers history fields.
module modal_aero_gasaerexch_cam
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: modal_aero_gasaerexch_cam_init
contains

  ! Resolve species indices from CAM constituent arrays,
  ! call portable modal_aero_gasaerexch_init,
  ! and register history fields (addfld) for gas-aerosol exchange diagnostics.
  subroutine modal_aero_gasaerexch_cam_init()
    use modal_aero_gasaerexch, only: modal_aero_gasaerexch_init
    use modal_aero_data, only: ntot_amode, nsoa, npoa, nspec_max, &
                               nspec_amode, modeptr_pcarbon, modeptr_accum, &
                               alnsg_amode, sigmag_amode, specmw_amode, specdens_amode, spechygro, &
                               lptr_so4_a_amode, lptr_nh4_a_amode, &
                               lptr2_soa_a_amode, lptr2_soa_g_amode, lptr2_pom_a_amode, &
                               numptr_amode, lmassptr_amode, cnst_name_cw
    use modal_aero_rename_cam, only: npair_renamexf, nspecfrm_renamexf, &
                                 lspecfrma_renamexf, lspectooa_renamexf, &
                                 lspecfrmc_renamexf, lspectooc_renamexf
    use radiative_aerosol, only: rad_aer_get_info
    use constituents, only: pcnst, cnst_name, cnst_get_ind
    use physconst, only: rair, mwdry, r_universal
    use cam_history, only: addfld, add_default, fieldname_len, horiz_only
    use cam_abortutils, only: endrun
    use spmd_utils, only: masterproc
    use phys_control, only: phys_getopts, cam_chempkg_is
    use cam_logfile, only: iulog

    ! indices are in pcnst (q) space
    ! run phase receives loffset for gas_pcnst (vmr) space

    ! local
    integer                          :: ipair, iq, iqfrm, iqtoo
    integer                          :: jac, jsoa, j
    integer                          :: l, lsfrm, lstoo
    integer                          :: mfrm, mtoo
    integer                          :: n, nspec
    integer                          :: nchfrm, nchfrmskip, nchtoo, nchtooskip

    logical                          :: dotend(pcnst), dotendqqcw(pcnst)

    character(len=fieldname_len)     :: tmpnamea
    character(len=fieldname_len + 3) :: fieldname
    character(128)                   :: long_name
    character(8)                     :: unit

    logical                          :: history_aerosol

    ! Resolved species indices (pcnst-space, for portable _init)
    integer                          :: idx_h2so4, idx_nh3, idx_msa
    integer                          :: idx_soag(nsoa)
    integer                          :: idx_so4_a(ntot_amode), idx_nh4_a(ntot_amode)
    integer                          :: idx_soa_a(ntot_amode, nsoa)
    integer                          :: idx_pom_a(ntot_amode, npoa)
    integer                          :: idx_num(ntot_amode)
    integer                          :: idx_mass(nspec_max, ntot_amode)

    ! pcage resolved arrays (pcnst-space)
    integer                          :: nspecfrm_pcage_resolved
    integer                          :: lspecfrm_pcage_resolved(nspec_max)
    integer                          :: lspectoo_pcage_resolved(nspec_max)

    ! SOA/POA molecular weights from host
    real(r8)                         :: mw_soa_host_resolved(nsoa)
    real(r8)                         :: mw_poa_host_resolved(npoa)

    character(len=32)                :: spec_type
    character(len=256) :: errmsg
    integer            :: errflg

    call phys_getopts(history_aerosol_out=history_aerosol)

    !-----------------------------------------------------------------------
    ! Part A: Resolve arguments and call portable _init
    !-----------------------------------------------------------------------

    ! --- Gas-phase species indices (pcnst-space) ---
    call cnst_get_ind('H2SO4', idx_h2so4, .false.)
    if ((idx_h2so4 <= 0) .or. (idx_h2so4 > pcnst)) then
      write (*, '(/a/a,i7)') &
        '*** modal_aero_gasaerexch_cam_init -- cannot find H2SO4 species', &
        '    idx_h2so4=', idx_h2so4
      call endrun('modal_aero_gasaerexch_cam_init error: H2SO4 not found')
    end if

    call cnst_get_ind('NH3', idx_nh3, .false.)
    if (.not. ((idx_nh3 > 0) .and. (idx_nh3 <= pcnst))) idx_nh3 = 0

    if (.not. cam_chempkg_is('geoschem_mam4')) then
      call cnst_get_ind('MSA', idx_msa, .false.)
    else
      idx_msa = 0
    end if
    if (.not. ((idx_msa > 0) .and. (idx_msa <= pcnst))) idx_msa = 0

    ! --- SOA gas-phase species indices (pcnst-space) ---
    do jsoa = 1, nsoa
      l = lptr2_soa_g_amode(jsoa)
      if ((l > 0) .and. (l <= pcnst)) then
        idx_soag(jsoa) = l
      else
        idx_soag(jsoa) = 0
      end if
    end do

    ! --- Aerosol species indices (per mode, pcnst-space) ---
    do n = 1, ntot_amode
      l = lptr_so4_a_amode(n)
      if ((l > 0) .and. (l <= pcnst)) then
        idx_so4_a(n) = l
      else
        idx_so4_a(n) = 0
      end if

      l = lptr_nh4_a_amode(n)
      if ((l > 0) .and. (l <= pcnst)) then
        idx_nh4_a(n) = l
      else
        idx_nh4_a(n) = 0
      end if

      do jsoa = 1, nsoa
        l = lptr2_soa_a_amode(n, jsoa)
        if ((l > 0) .and. (l <= pcnst)) then
          idx_soa_a(n, jsoa) = l
        else
          idx_soa_a(n, jsoa) = 0
        end if
      end do

      do j = 1, npoa
        l = lptr2_pom_a_amode(n, j)
        if ((l > 0) .and. (l <= pcnst)) then
          idx_pom_a(n, j) = l
        else
          idx_pom_a(n, j) = 0
        end if
      end do

      l = numptr_amode(n)
      if ((l > 0) .and. (l <= pcnst)) then
        idx_num(n) = l
      else
        idx_num(n) = 0
      end if

      do l = 1, nspec_amode(n)
        idx_mass(l, n) = lmassptr_amode(l, n)
      end do
    end do

    ! --- Resolve pcage species matching ---
    ! Copied from original modal_aero_gasaerexch_init
    ! Indices stored in pcnst-space.
    nspecfrm_pcage_resolved = 0
    lspecfrm_pcage_resolved(:) = 0
    lspectoo_pcage_resolved(:) = 0

    if ((modeptr_pcarbon > 0) .and. (modeptr_accum > 0)) then
      l = lptr_so4_a_amode(modeptr_accum)
      if ((l >= 1) .and. (l <= pcnst)) then

        mfrm = modeptr_pcarbon
        mtoo = modeptr_accum

        if (mfrm < 10) then
          nchfrmskip = 1
        else if (mfrm < 100) then
          nchfrmskip = 2
        else
          nchfrmskip = 3
        end if
        if (mtoo < 10) then
          nchtooskip = 1
        else if (mtoo < 100) then
          nchtooskip = 2
        else
          nchtooskip = 3
        end if
        nspec = 0

        aa_iqfrm: do iqfrm = -1, nspec_amode(mfrm)

          if (iqfrm == -1) then
            lsfrm = numptr_amode(mfrm)
            lstoo = numptr_amode(mtoo)
          else if (iqfrm == 0) then
            ! bypass transfer of aerosol water due to primary-carbon aging
            cycle aa_iqfrm
          else
            lsfrm = lmassptr_amode(iqfrm, mfrm)
            lstoo = 0
          end if
          if ((lsfrm < 1) .or. (lsfrm > pcnst)) cycle aa_iqfrm

          if (lsfrm > 0 .and. iqfrm > 0) then
            nchfrm = len(trim(cnst_name(lsfrm))) - nchfrmskip

            ! find "too" species having same cnst_name (except for last 1/2/3 characters which are the mode index)
            do iqtoo = 1, nspec_amode(mtoo)
              lstoo = lmassptr_amode(iqtoo, mtoo)
              nchtoo = len(trim(cnst_name(lstoo))) - nchtooskip
              if (cnst_name(lsfrm) (1:nchfrm) == cnst_name(lstoo) (1:nchtoo)) then
                exit
              else
                lstoo = 0
              end if
            end do
          end if

          if ((lstoo < 1) .or. (lstoo > pcnst)) lstoo = 0
          nspec = nspec + 1
          lspecfrm_pcage_resolved(nspec) = lsfrm
          lspectoo_pcage_resolved(nspec) = lstoo
        end do aa_iqfrm

        nspecfrm_pcage_resolved = nspec

        ! output results
        if (masterproc) then
          write (iulog, 9310)
          write (iulog, 9320) 1, mfrm, mtoo

          do iq = 1, nspecfrm_pcage_resolved
            lsfrm = lspecfrm_pcage_resolved(iq)
            lstoo = lspectoo_pcage_resolved(iq)
            if (lstoo > 0) then
              write (iulog, 9330) lsfrm, cnst_name(lsfrm), &
                lstoo, cnst_name(lstoo)
            else
              write (iulog, 9340) lsfrm, cnst_name(lsfrm)
            end if
          end do

          write (iulog, *)
        end if ! ( masterproc )

9310    format(/'subr. modal_aero_gasaerexch_cam_init - primary carbon aging pointers')
9320    format('pair', i3, 5x, 'mode', i3, ' ---> mode', i3)
9330    format(5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a)
9340    format(5x, 'spec', i3, '=', a, ' ---> LOSS')

      end if ! lptr_so4_a_amode(modeptr_accum) valid
    end if ! modeptr_pcarbon > 0 .and. modeptr_accum > 0

    ! --- Resolve SOA/POA molecular weights using rad_aer_get_info ---
    mw_soa_host_resolved(:) = 0.0_r8
    mw_poa_host_resolved(:) = 0.0_r8

    do n = 1, ntot_amode
      do l = 1, nspec_amode(n)
        call rad_aer_get_info(0, n, l, spec_type=spec_type)
        select case (spec_type)
        case ('s-organic')
          mw_soa_host_resolved(:) = specmw_amode(l, n)
        case ('p-organic')
          mw_poa_host_resolved(:) = specmw_amode(l, n)
        end select
      end do
    end do

    ! --- Call portable init ---
    call modal_aero_gasaerexch_init( &
      ntot_amode         = ntot_amode,                &
      nsoa               = nsoa,                      &
      npoa               = npoa,                      &
      nspec_max          = nspec_max,                 &
      nspec_amode        = nspec_amode,               &
      modeptr_pcarbon    = modeptr_pcarbon,           &
      modeptr_accum      = modeptr_accum,             &
      alnsg_amode        = alnsg_amode,               &
      sigmag_amode       = sigmag_amode,              &
      specmw_amode       = specmw_amode,              &
      specdens_amode     = specdens_amode,            &
      spechygro          = spechygro,                 &
      idx_h2so4          = idx_h2so4,                 &
      idx_nh3            = idx_nh3,                   &
      idx_msa            = idx_msa,                   &
      idx_soag           = idx_soag,                  &
      idx_so4_a          = idx_so4_a,                 &
      idx_nh4_a          = idx_nh4_a,                 &
      idx_soa_a          = idx_soa_a,                 &
      idx_pom_a          = idx_pom_a,                 &
      idx_num            = idx_num,                   &
      idx_mass           = idx_mass,                  &
      pcnst_in           = pcnst,                     &
      nspecfrm_pcage_in  = nspecfrm_pcage_resolved,   &
      lspecfrm_pcage_in  = lspecfrm_pcage_resolved,   &
      lspectoo_pcage_in  = lspectoo_pcage_resolved,   &
      mw_soa_host        = mw_soa_host_resolved,      &
      mw_poa_host        = mw_poa_host_resolved,      &
      rair               = rair,                      &
      mwdry              = mwdry,                     &
      r_universal        = r_universal,               &
      errmsg             = errmsg,                    &
      errflg             = errflg)

    if (errflg /= 0) then
      call endrun('modal_aero_gasaerexch_cam_init: '//trim(errmsg))
    end if

    !-----------------------------------------------------------------------
    ! Part B: History field registration (addfld calls)
    !-----------------------------------------------------------------------

    ! --- Tendency flags for _sfgaex1 fields ---
    ! Determine which constituents get gas-aerosol exchange tendency output
    dotend(:) = .false.

    ! H2SO4 (required). below indices are in pcnst indexing:
    dotend(idx_h2so4) = .true.

    ! NH3 (optional)
    if ((idx_nh3 > 0) .and. (idx_nh3 <= pcnst)) dotend(idx_nh3) = .true.

    ! MSA (optional)
    if ((idx_msa > 0) .and. (idx_msa <= pcnst)) dotend(idx_msa) = .true.

    ! SOA gases
    do jsoa = 1, nsoa
      l = lptr2_soa_g_amode(jsoa)
      if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
    end do

    ! Aerosol species in modes
    do n = 1, ntot_amode
      l = lptr_so4_a_amode(n)
      if ((l > 0) .and. (l <= pcnst)) then
        dotend(l) = .true.
        if (idx_nh3 > 0) then
          l = lptr_nh4_a_amode(n)
          if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
        end if
      end if
      do jsoa = 1, nsoa
        if (idx_soag(jsoa) > 0) then
          l = lptr2_soa_a_amode(n, jsoa)
          if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
        end if
      end do
    end do

    ! pcage species
    if (nspecfrm_pcage_resolved > 0) then
      do iq = 1, nspecfrm_pcage_resolved
        ! lspec*_pcage_resolved are already pcnst-space (resolved above), so
        ! they index dotend directly -- matching the original gasaerexch_init.
        lsfrm = lspecfrm_pcage_resolved(iq)
        lstoo = lspectoo_pcage_resolved(iq)
        if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
          dotend(lsfrm) = .true.
          if ((lstoo > 0) .and. (lstoo <= pcnst)) then
            dotend(lstoo) = .true.
          end if
        end if
      end do
    end if

    ! --- SOA condensation/evaporation diagnostics ---
    fieldname = trim('qconff_gaex')
    long_name = trim('3D fields for Fossil SOA condensation')
    unit = 'kg/kg/s'
    call addfld(fieldname, (/'lev'/), 'A', unit, long_name)
    if (history_aerosol) then
      call add_default(fieldname, 1, ' ')
    end if
    if (masterproc) write (*, '(3(a,3x))') 'qconff addfld', fieldname, unit

    fieldname = trim('qevapff_gaex')
    long_name = trim('3D fields for Fossil SOA evaporation')
    call addfld(fieldname, (/'lev'/), 'A', unit, long_name)
    if (history_aerosol) then
      call add_default(fieldname, 1, ' ')
    end if
    if (masterproc) write (*, '(3(a,3x))') 'qevapff addfld', fieldname, unit

    fieldname = trim('qconbb_gaex')
    long_name = trim('3D fields for Biomass SOA condensation')
    call addfld(fieldname, (/'lev'/), 'A', unit, long_name)
    if (history_aerosol) then
      call add_default(fieldname, 1, ' ')
    end if
    if (masterproc) write (*, '(3(a,3x))') 'qconbb addfld', fieldname, unit

    fieldname = trim('qevapbb_gaex')
    long_name = trim('3D fields for Biomass SOA evaporation')
    call addfld(fieldname, (/'lev'/), 'A', unit, long_name)
    if (history_aerosol) then
      call add_default(fieldname, 1, ' ')
    end if
    if (masterproc) write (*, '(3(a,3x))') 'qevapbb addfld', fieldname, unit

    fieldname = trim('qconbg_gaex')
    long_name = trim('3D fields for Biogenic SOA condensation')
    call addfld(fieldname, (/'lev'/), 'A', unit, long_name)
    if (history_aerosol) then
      call add_default(fieldname, 1, ' ')
    end if
    if (masterproc) write (*, '(3(a,3x))') 'qconbg addfld', fieldname, unit

    fieldname = trim('qevapbg_gaex')
    long_name = trim('3D fields for Biogenic SOA evaporation')
    call addfld(fieldname, (/'lev'/), 'A', unit, long_name)
    if (history_aerosol) then
      call add_default(fieldname, 1, ' ')
    end if
    if (masterproc) write (*, '(3(a,3x))') 'qevapbg addfld', fieldname, unit

    fieldname = trim('qcon_gaex')
    long_name = trim('3D fields for SOA condensation')
    call addfld(fieldname, (/'lev'/), 'A', unit, long_name)
    if (history_aerosol) then
      call add_default(fieldname, 1, ' ')
    end if
    if (masterproc) write (*, '(3(a,3x))') 'qcon addfld', fieldname, unit

    fieldname = trim('qevap_gaex')
    long_name = trim('3D fields for Biogenic SOA evaporation')
    call addfld(fieldname, (/'lev'/), 'A', unit, long_name)
    if (history_aerosol) then
      call add_default(fieldname, 1, ' ')
    end if
    if (masterproc) write (*, '(3(a,3x))') 'qevap addfld', fieldname, unit

    ! --- Per-species _sfgaex1 fields (gas-aerosol exchange primary tendency) ---
    do l = 1, pcnst
      if (.not. dotend(l)) cycle

      tmpnamea = cnst_name(l)
      fieldname = trim(tmpnamea)//'_sfgaex1'
      long_name = trim(tmpnamea)//' gas-aerosol-exchange primary column tendency'
      unit = 'kg/m2/s'
      call addfld(fieldname, horiz_only, 'A', unit, long_name)
      if (history_aerosol) then
        call add_default(fieldname, 1, ' ')
      end if
      if (masterproc) write (*, '(3(a,3x))') 'gasaerexch addfld', fieldname, unit

    end do   ! l = ...

    ! --- Per-species _sfgaex2 fields (renaming tendency) ---
    dotend(:) = .false.
    dotendqqcw(:) = .false.
    do ipair = 1, npair_renamexf
      do iq = 1, nspecfrm_renamexf(ipair)
        lsfrm = lspecfrma_renamexf(iq, ipair)
        lstoo = lspectooa_renamexf(iq, ipair)
        if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
          dotend(lsfrm) = .true.
          if ((lstoo > 0) .and. (lstoo <= pcnst)) then
            dotend(lstoo) = .true.
          end if
        end if

        lsfrm = lspecfrmc_renamexf(iq, ipair)
        lstoo = lspectooc_renamexf(iq, ipair)
        if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
          dotendqqcw(lsfrm) = .true.
          if ((lstoo > 0) .and. (lstoo <= pcnst)) then
            dotendqqcw(lstoo) = .true.
          end if
        end if
      end do ! iq = ...
    end do ! ipair = ...

    do l = 1, pcnst
    do jac = 1, 2
      if (jac == 1) then
        if (.not. dotend(l)) cycle
        tmpnamea = cnst_name(l)
      else
        if (.not. dotendqqcw(l)) cycle
        tmpnamea = cnst_name_cw(l)
      end if

      fieldname = trim(tmpnamea)//'_sfgaex2'
      long_name = trim(tmpnamea)//' gas-aerosol-exchange renaming column tendency'
      unit = 'kg/m2/s'
      if ((tmpnamea(1:3) == 'num') .or. &
          (tmpnamea(1:3) == 'NUM')) unit = '#/m2/s'
      call addfld(fieldname, horiz_only, 'A', unit, long_name)
      if (history_aerosol) then
        call add_default(fieldname, 1, ' ')
      end if
      if (masterproc) write (*, '(3(a,3x))') 'gasaerexch addfld', fieldname, unit
    end do   ! jac = ...
    end do   ! l = ...

  end subroutine modal_aero_gasaerexch_cam_init
end module modal_aero_gasaerexch_cam
