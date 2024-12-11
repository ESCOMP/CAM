module carma_aerosol_properties_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: pi
  use aerosol_properties_mod, only: aerosol_properties, aero_name_len
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_bin_props_by_idx, &
                              rad_cnst_get_info_by_bin, rad_cnst_get_info_by_bin_spec, rad_cnst_get_bin_props
  use infnan, only: nan, assignment(=)

  implicit none

  private

  public :: carma_aerosol_properties

  type, extends(aerosol_properties) :: carma_aerosol_properties
     private
     integer, allocatable :: ibl(:)
   contains
     procedure :: number_transported
     procedure :: get
     procedure :: amcube
     procedure :: actfracs
     procedure :: num_names
     procedure :: mmr_names
     procedure :: amb_num_name
     procedure :: amb_mmr_name
     procedure :: species_type
     procedure :: icenuc_updates_num
     procedure :: icenuc_updates_mmr
     procedure :: apply_number_limits
     procedure :: hetfrz_species
     procedure :: optics_params
     procedure :: nbins_rlist
     procedure :: nspecies_per_bin_rlist
     procedure :: alogsig_rlist
     procedure :: soluble
     procedure :: min_mass_mean_rad
     procedure :: bin_name
     procedure :: scav_diam
     procedure :: resuspension_resize
     procedure :: rebin_bulk_fluxes
     procedure :: hydrophilic

     final :: destructor
  end type carma_aerosol_properties

  interface carma_aerosol_properties
     procedure :: constructor
  end interface carma_aerosol_properties

  real(r8), parameter :: onethird = 1._r8/3._r8

contains

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  function constructor() result(newobj)

    type(carma_aerosol_properties), pointer :: newobj

    integer :: l, m, nbins, ncnst_tot
    integer,allocatable :: nspecies(:)
    integer,allocatable :: nmasses(:)
    real(r8),allocatable :: alogsig(:)
    real(r8),allocatable :: f1(:)
    real(r8),allocatable :: f2(:)
    integer :: ierr

    integer, pointer :: ibl(:)
    integer :: ii, imx, imx_num, imx_mmr, ipr, ipr_num, ipr_mmr
    character(len=32) :: spectype
    character(len=32) :: bin_name
    character(len=32) :: bin_name_l    ! bin name of the larger bin

    integer, allocatable :: imx_bl(:)     ! index used to map pure sulfate bin to mixed sulfate bin
    integer, allocatable :: imx_mmr_bl(:) ! index used to map pure sulfate bin to mixed sulfate bin for mmr
    integer, allocatable :: imx_num_bl(:) ! index used to map pure sulfate bin to mixed sulfate bin for num

    allocate(newobj,stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    call rad_cnst_get_info( 0, nbins=nbins)

    allocate( nspecies(nbins),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate( nmasses(nbins),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate( alogsig(nbins),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate( f1(nbins),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    allocate( f2(nbins),stat=ierr )
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if

    ncnst_tot = 0

    do m = 1, nbins
       call rad_cnst_get_info_by_bin(0, m, nspec=nspecies(m))
       ncnst_tot = ncnst_tot + nspecies(m) + 1
       nmasses(m) = nspecies(m)
    end do

    alogsig(:) = log(2._r8)  !!!! ???? IS THIS RIGHT ???? !!!
    f1 = 1._r8
    f2 = 1._r8

    call newobj%initialize(nbins,ncnst_tot,nspecies,nmasses,alogsig,f1,f2,ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    deallocate(nspecies)
    deallocate(nmasses)
    deallocate(alogsig)
    deallocate(f1)
    deallocate(f2)

    allocate(newobj%ibl(ncnst_tot),stat=ierr)
    if( ierr /= 0 ) then
       nullify(newobj)
       return
    end if
    ibl => newobj%ibl

    ibl = -1

    allocate(imx_num_bl(nbins))
    allocate(imx_mmr_bl(nbins))
    allocate(imx_bl(nbins))

    imx = 0
    imx_mmr = 0
    imx_num = 0
    ipr = 0
    ipr_mmr = 0
    ipr_num = 0

    do m = 1,nbins
       bin_name = newobj%bin_name(0,m)
       bin_name_l = ' '
       if (m<nbins) then
          bin_name_l = newobj%bin_name(0,m+1)
       end if

       do l = 0,newobj%nspecies(m)
          ii = newobj%indexer(m,l)
          ibl(ii) = ii

          ! derive index  array for larger bin, for evaporation into larger bi
          if (l>0 .and. l<=newobj%nspecies(m)) then
             call newobj%species_type(m,l,spectype)
          else
             spectype = 'other'
          end if
          ! identification is required for pure and mixed aerosols, mixed aeroosols are moved to
          ! larger bin, pure aerosols are moved to mixed sulfate


          if (index(bin_name,'MXAER')>0 .and. index(bin_name_l,'MXAER')>0) then
             ! for mixed aerosols
             ! find larger bin
             ibl(ii) = newobj%indexer(m+1,l)
             ! define mixed aerosol sulfate index to be used for pure sulfate only
             if (trim(spectype) == 'sulfate') then
                imx = imx + 1
                imx_bl(imx) = ibl(ii)
             end if
             if (l == newobj%nspecies(m)+1) then !  only for mmr
                imx_mmr = imx_mmr + 1
                ibl(ii) = newobj%indexer(m+1,l)
                imx_mmr_bl(imx_mmr) = ibl(ii)
             end if
             if (l == 0) then !  only for num
                imx_num = imx_num + 1
                ibl(ii) =  newobj%indexer(m+1,l)
                imx_num_bl(imx_num) = ibl(ii)
             end if
          end if ! MXAER

          if (index(bin_name,'PRSUL')>0 .and. index(bin_name_l,'PRSUL')>0) then
             ! assuming  PRSULF and  MXSULF have the same number of bins
             if (trim(spectype) == 'sulfate') then
                ipr = ipr +1
                ibl(ii) = imx_bl(ipr)
             end if
             if (l == newobj%nspecies(m)+1) then ! only for mmr reset counter to only go from 1-20 bins
                ipr_mmr = ipr_mmr + 1
                ibl(ii) = imx_mmr_bl(ipr_mmr)
             end if
             if (l == 0 ) then ! only for num reset counter to only go from 1-20 bins
                ipr_num = ipr_num + 1
                ibl(ii) = imx_num_bl(ipr_num)
             end if
          end if
          if (ibl(ii).eq.0) then
             ibl(ii) = ii
          end if
       end do
    end do

    deallocate(imx_mmr_bl, imx_num_bl, imx_bl)

  end function constructor

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  subroutine destructor(self)
    type(carma_aerosol_properties), intent(inout) :: self

    call self%final()

  end subroutine destructor

  !------------------------------------------------------------------------------
  ! returns number of transported aerosol constituents
  !------------------------------------------------------------------------------
  integer function number_transported(self)
    class(carma_aerosol_properties), intent(in) :: self
    ! to be implemented later
    number_transported = -1
  end function number_transported

  !------------------------------------------------------------------------
  ! returns aerosol properties:
  !  density
  !  hygroscopicity
  !  species type
  !  species name
  !  short wave species refractive indices
  !  long wave species refractive indices
  !  species morphology
  !------------------------------------------------------------------------
  subroutine get(self, bin_ndx, species_ndx, list_ndx, density, hygro, &
                 spectype, specname, specmorph, refindex_sw, refindex_lw)

    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    integer, intent(in) :: species_ndx         ! species index
    integer, optional, intent(in) :: list_ndx  ! climate or a diagnostic list number
    real(r8), optional, intent(out) :: density ! density (kg/m3)
    real(r8), optional, intent(out) :: hygro   ! hygroscopicity
    character(len=*), optional, intent(out) :: spectype  ! species type
    character(len=*), optional, intent(out) :: specname  ! species name
    character(len=*), optional, intent(out) :: specmorph ! species morphology
    complex(r8), pointer, optional, intent(out) :: refindex_sw(:) ! short wave species refractive indices
    complex(r8), pointer, optional, intent(out) :: refindex_lw(:) ! long wave species refractive indices

    integer :: ilist

    if (present(list_ndx)) then
       ilist = list_ndx
    else
       ilist = 0
    end if

    if (present(density)) then
       call rad_cnst_get_bin_props_by_idx(ilist, bin_ndx, species_ndx, density_aer=density)
    end if
    if (present(hygro)) then
       call rad_cnst_get_bin_props_by_idx(ilist, bin_ndx, species_ndx, hygro_aer=hygro)
    end if
    if (present(spectype)) then
       call rad_cnst_get_bin_props_by_idx(ilist, bin_ndx, species_ndx, spectype=spectype)
    end if
    if (present(refindex_sw)) then
       call rad_cnst_get_bin_props_by_idx(ilist, bin_ndx, species_ndx, refindex_aer_sw=refindex_sw)
    end if
    if (present(refindex_lw)) then
       call rad_cnst_get_bin_props_by_idx(ilist, bin_ndx, species_ndx, refindex_aer_lw=refindex_lw)
    end if
    if (present(specmorph)) then
       call rad_cnst_get_bin_props_by_idx(ilist, bin_ndx, species_ndx, specmorph=specmorph)
    end if
    if (present(specname)) then
       if (species_ndx>self%nspecies(bin_ndx)) then
          call rad_cnst_get_info_by_bin(0, bin_ndx,  mmr_name=specname)
       else
          call rad_cnst_get_info_by_bin_spec(ilist, bin_ndx, species_ndx, spec_name=specname)
       end if
    end if

  end subroutine get

  !------------------------------------------------------------------------
  ! returns optics type and table parameters
  !------------------------------------------------------------------------
  subroutine optics_params(self, list_ndx, bin_ndx, opticstype, extpsw, abspsw, asmpsw, absplw, &
       refrtabsw, refitabsw, refrtablw, refitablw, ncoef, prefr, prefi, sw_hygro_ext_wtp, &
       sw_hygro_ssa_wtp, sw_hygro_asm_wtp, lw_hygro_ext_wtp, wgtpct, nwtp, &
       sw_hygro_coreshell_ext, sw_hygro_coreshell_ssa, sw_hygro_coreshell_asm, lw_hygro_coreshell_ext, &
       corefrac, bcdust, kap, relh, nfrac, nbcdust, nkap, nrelh )

    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx             ! bin index
    integer, intent(in) :: list_ndx            ! rad climate/diags list

    character(len=*), optional, intent(out) :: opticstype

    ! refactive index table parameters
    real(r8),  optional, pointer     :: extpsw(:,:,:,:) ! short wave specific extinction
    real(r8),  optional, pointer     :: abspsw(:,:,:,:) ! short wave specific absorption
    real(r8),  optional, pointer     :: asmpsw(:,:,:,:) ! short wave asymmetry factor
    real(r8),  optional, pointer     :: absplw(:,:,:,:) ! long wave specific absorption
    real(r8),  optional, pointer     :: refrtabsw(:,:)  ! table of short wave real refractive indices for aerosols
    real(r8),  optional, pointer     :: refitabsw(:,:)  ! table of short wave imaginary refractive indices for aerosols
    real(r8),  optional, pointer     :: refrtablw(:,:)  ! table of long wave real refractive indices for aerosols
    real(r8),  optional, pointer     :: refitablw(:,:)  ! table of long wave imaginary refractive indices for aerosols
    integer,   optional, intent(out) :: ncoef  ! number of chebychev polynomials
    integer,   optional, intent(out) :: prefr  ! number of real refractive indices in table
    integer,   optional, intent(out) :: prefi  ! number of imaginary refractive indices in table

    ! hygrowghtpct table parameters
    real(r8),  optional, pointer     :: sw_hygro_ext_wtp(:,:) ! short wave extinction table
    real(r8),  optional, pointer     :: sw_hygro_ssa_wtp(:,:) ! short wave single-scatter albedo table
    real(r8),  optional, pointer     :: sw_hygro_asm_wtp(:,:) ! short wave asymmetry table
    real(r8),  optional, pointer     :: lw_hygro_ext_wtp(:,:) ! long wave absorption table
    real(r8),  optional, pointer     :: wgtpct(:)   ! weight precent of H2SO4/H2O solution
    integer,   optional, intent(out) :: nwtp        ! number of weight precent values

    ! hygrocoreshell table parameters
    real(r8),  optional, pointer     :: sw_hygro_coreshell_ext(:,:,:,:,:) ! short wave extinction table
    real(r8),  optional, pointer     :: sw_hygro_coreshell_ssa(:,:,:,:,:) ! short wave single-scatter albedo table
    real(r8),  optional, pointer     :: sw_hygro_coreshell_asm(:,:,:,:,:) ! short wave asymmetry table
    real(r8),  optional, pointer     :: lw_hygro_coreshell_ext(:,:,:,:,:) ! long wave absorption table
    real(r8),  optional, pointer     :: corefrac(:) ! core fraction dimension values
    real(r8),  optional, pointer     :: bcdust(:)   ! bc/(bc + dust) fraction dimension values
    real(r8),  optional, pointer     :: kap(:)      ! hygroscopicity dimension values
    real(r8),  optional, pointer     :: relh(:)     ! relative humidity dimension values
    integer,   optional, intent(out) :: nfrac       ! core fraction dimension size
    integer,   optional, intent(out) :: nbcdust     ! bc/(bc + dust) fraction dimension size
    integer,   optional, intent(out) :: nkap        ! hygroscopicity dimension size
    integer,   optional, intent(out) :: nrelh       ! relative humidity dimension size

    if (present(extpsw)) then
       nullify(extpsw)
    end if
    if (present(abspsw)) then
       nullify(abspsw)
    end if
    if (present(asmpsw)) then
       nullify(asmpsw)
    end if
    if (present(absplw)) then
       nullify(absplw)
    end if
    if (present(refrtabsw)) then
       nullify(refrtabsw)
    end if
    if (present(refitabsw)) then
       nullify(refitabsw)
    end if
    if (present(refrtablw)) then
       nullify(refrtablw)
    end if
    if (present(refitablw)) then
       nullify(refitablw)
    end if
    if (present(ncoef)) then
       ncoef = huge(1)
    end if
    if (present(prefr)) then
       prefr = huge(1)
    end if
    if (present(prefi)) then
       prefi = huge(1)
    end if

    call rad_cnst_get_bin_props(list_ndx,bin_ndx, &
                                opticstype=opticstype, &
                                sw_hygro_ext_wtp=sw_hygro_ext_wtp, &
                                sw_hygro_ssa_wtp=sw_hygro_ssa_wtp, &
                                sw_hygro_asm_wtp=sw_hygro_asm_wtp, &
                                lw_hygro_ext_wtp=lw_hygro_ext_wtp, &
                                wgtpct=wgtpct, &
                                nwtp=nwtp, &
                                sw_hygro_coreshell_ext=sw_hygro_coreshell_ext, &
                                sw_hygro_coreshell_ssa=sw_hygro_coreshell_ssa, &
                                sw_hygro_coreshell_asm=sw_hygro_coreshell_asm, &
                                lw_hygro_coreshell_ext=lw_hygro_coreshell_ext, &
                                corefrac=corefrac, &
                                bcdust=bcdust, &
                                kap=kap, &
                                relh=relh, &
                                nbcdust=nbcdust, &
                                nkap=nkap, &
                                nrelh=nrelh, &
                                nfrac=nfrac )

  end subroutine optics_params

  !------------------------------------------------------------------------------
  ! returns radius^3 (m3) of a given bin number
  !------------------------------------------------------------------------------
  pure elemental real(r8) function amcube(self, bin_ndx, volconc, numconc)

    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number
    real(r8), intent(in) :: volconc ! volume conc (m3/m3)
    real(r8), intent(in) :: numconc ! number conc (1/m3)

    amcube = 3._r8/(4._r8*pi)*volconc/numconc

  end function amcube

  !------------------------------------------------------------------------------
  ! returns mass and number activation fractions
  !------------------------------------------------------------------------------
  subroutine actfracs(self, bin_ndx, smc, smax, fn, fm )
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx   ! bin index
    real(r8),intent(in) :: smc       ! critical supersaturation for particles of bin radius
    real(r8),intent(in) :: smax      ! maximum supersaturation for multiple competing aerosols
    real(r8),intent(out) :: fn       ! activation fraction for aerosol number
    real(r8),intent(out) :: fm       ! activation fraction for aerosol mass

    fn = 0._r8
    fm = 0._r8

    if (smc < smax) then
       fn = 1._r8
       fm = 1._r8
    end if

  end subroutine actfracs

  !------------------------------------------------------------------------
  ! returns constituents names of aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine num_names(self, bin_ndx, name_a, name_c)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), intent(out) :: name_a ! constituent name of ambient aerosol number dens
    character(len=*), intent(out) :: name_c ! constituent name of cloud-borne aerosol number dens

    call rad_cnst_get_info_by_bin(0, bin_ndx, num_name=name_a, num_name_cw=name_c)

  end subroutine num_names

  !------------------------------------------------------------------------
  ! returns constituents names of aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine mmr_names(self, bin_ndx, species_ndx, name_a, name_c)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: name_a ! constituent name of ambient aerosol MMR
    character(len=*), intent(out) :: name_c ! constituent name of cloud-borne aerosol MMR

    if (species_ndx>0) then
       call rad_cnst_get_info_by_bin_spec(0, bin_ndx, species_ndx, spec_name=name_a, spec_name_cw=name_c)
    else
       call rad_cnst_get_info_by_bin(0, bin_ndx,  mmr_name=name_a, mmr_name_cw=name_c)
    end if

  end subroutine mmr_names

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol number mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_num_name(self, bin_ndx, name)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    character(len=*), intent(out) :: name   ! constituent name of ambient aerosol number dens

    call rad_cnst_get_info_by_bin(0, bin_ndx, num_name=name)

  end subroutine amb_num_name

  !------------------------------------------------------------------------
  ! returns constituent name of ambient aerosol mass mixing ratios
  !------------------------------------------------------------------------
  subroutine amb_mmr_name(self, bin_ndx, species_ndx, name)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: name   ! constituent name of ambient aerosol MMR

    if (species_ndx>0) then
       call rad_cnst_get_info_by_bin_spec(0, bin_ndx, species_ndx, spec_name=name)
    else
       call rad_cnst_get_info_by_bin(0, bin_ndx,  mmr_name=name)
    end if

  end subroutine amb_mmr_name

  !------------------------------------------------------------------------
  ! returns species type
  !------------------------------------------------------------------------
  subroutine species_type(self, bin_ndx, species_ndx, spectype)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number
    character(len=*), intent(out) :: spectype ! species type

    call rad_cnst_get_info_by_bin_spec(0, bin_ndx, species_ndx, spec_type=spectype)

  end subroutine species_type

  !------------------------------------------------------------------------------
  ! returns TRUE if Ice Nucleation tendencies are applied to given aerosol bin number
  !------------------------------------------------------------------------------
  function icenuc_updates_num(self, bin_ndx) result(res)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number

    logical :: res

    character(len=aero_name_len) :: spectype
    integer :: spc_ndx

    res = .false.

    do spc_ndx = 1, self%nspecies(bin_ndx)
       call self%species_type( bin_ndx, spc_ndx, spectype)
       if (trim(spectype)=='dust') res = .true.
       if (trim(spectype)=='sulfate') res = .true.
    end do

  end function icenuc_updates_num

  !------------------------------------------------------------------------------
  ! returns TRUE if Ice Nucleation tendencies are applied to a given species within a bin
  !------------------------------------------------------------------------------
  function icenuc_updates_mmr(self, bin_ndx, species_ndx) result(res)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number

    logical :: res

    character(len=aero_name_len) :: spectype

    res = .false.

    if (species_ndx==0) then
       res = self%icenuc_updates_num(bin_ndx)
    else
       call self%species_type( bin_ndx, species_ndx, spectype)
       if (trim(spectype)=='dust') res = .true.
       if (trim(spectype)=='sulfate') res = .true.
    end if

  end function icenuc_updates_mmr

  !------------------------------------------------------------------------------
  ! apply max / min to number concentration
  !------------------------------------------------------------------------------
  subroutine apply_number_limits( self, naerosol, vaerosol, istart, istop, m )
    class(carma_aerosol_properties), intent(in) :: self
    real(r8), intent(inout) :: naerosol(:)  ! number conc (1/m3)
    real(r8), intent(in)    :: vaerosol(:)  ! volume conc (m3/m3)
    integer,  intent(in) :: istart          ! start column index (1 <= istart <= istop <= pcols)
    integer,  intent(in) :: istop           ! stop column index
    integer,  intent(in) :: m               ! mode or bin index

  end subroutine apply_number_limits

  !------------------------------------------------------------------------------
  ! returns TRUE if species `spc_ndx` in aerosol subset `bin_ndx` contributes to
  ! the particles' ability to act as heterogeneous freezing nuclei
  !------------------------------------------------------------------------------
  function hetfrz_species(self, bin_ndx, spc_ndx) result(res)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number
    integer, intent(in) :: spc_ndx  ! species number

    logical :: res

    character(len=aero_name_len) :: species_type

    res = .false.

    call self%species_type(bin_ndx, spc_ndx, species_type)
    if ( trim(species_type)=='black-c' .or. trim(species_type)=='dust' ) then
       res = .true.
    end if

  end function hetfrz_species

  !------------------------------------------------------------------------------
  ! returns TRUE if soluble
  !------------------------------------------------------------------------------
  logical function soluble(self,bin_ndx)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number

    soluble = .true.

  end function soluble

  !------------------------------------------------------------------------------
  ! returns minimum mass mean radius (meters)
  !------------------------------------------------------------------------------
  function min_mass_mean_rad(self,bin_ndx,species_ndx) result(minrad)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx           ! bin number
    integer, intent(in) :: species_ndx       ! species number

    real(r8) :: minrad  ! meters

    minrad = 0.0_r8

  end function min_mass_mean_rad

  !------------------------------------------------------------------------------
  ! returns the total number of bins for a given radiation list index
  !------------------------------------------------------------------------------
  function nbins_rlist(self, list_ndx)  result(res)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx  ! radiation list number

    integer :: res

    call rad_cnst_get_info(list_ndx, nbins=res)

  end function nbins_rlist

  !------------------------------------------------------------------------------
  ! returns number of species in a bin for a given radiation list index
  !------------------------------------------------------------------------------
  function nspecies_per_bin_rlist(self, list_ndx,  bin_ndx)  result(res)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx ! radiation list number
    integer, intent(in) :: bin_ndx  ! bin number

    integer :: res

    call rad_cnst_get_info_by_bin(list_ndx, bin_ndx, nspec=res)

  end function nspecies_per_bin_rlist

  !------------------------------------------------------------------------------
  ! returns the natural log of geometric standard deviation of the number
  ! distribution for radiation list number and aerosol bin
  !------------------------------------------------------------------------------
  function alogsig_rlist(self, list_ndx,  bin_ndx)  result(res)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx ! radiation list number
    integer, intent(in) :: bin_ndx  ! bin number

    real(r8) :: res

    res = self%alogsig(bin_ndx)

  end function alogsig_rlist

  !------------------------------------------------------------------------------
  ! returns name for a given radiation list number and aerosol bin
  !------------------------------------------------------------------------------
  function bin_name(self, list_ndx,  bin_ndx) result(name)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: list_ndx ! radiation list number
    integer, intent(in) :: bin_ndx  ! bin number

    character(len=32) name

    call rad_cnst_get_info_by_bin(list_ndx, bin_ndx, bin_name=name)

  end function bin_name

  !------------------------------------------------------------------------------
  ! returns scavenging diameter (cm) for a given aerosol bin number
  !------------------------------------------------------------------------------
  function scav_diam(self, bin_ndx) result(diam)

    use carma_intr, only: carma_get_bin_rmass
    use carma_intr, only: carma_get_group_by_name

    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx  ! bin number

    real(r8) :: diam ! cm

    real(r8) :: mass   ! the bin mass (g)
    real(r8) :: rho    ! density (kg/m3)
    integer :: ispec
    character(len=32) :: spectype

    character(len=aero_name_len) :: bin_name, shortname
    integer :: igroup, ibin, rc, nchr

    call rad_cnst_get_info_by_bin(0, bin_ndx, bin_name=bin_name)

    nchr = len_trim(bin_name)-2
    shortname = bin_name(:nchr)

    call carma_get_group_by_name(shortname, igroup, rc)

    read(bin_name(nchr+1:),*) ibin

    call carma_get_bin_rmass(igroup, ibin, mass, rc)

    do ispec = 1, self%nspecies(bin_ndx)
       call self%species_type(bin_ndx,ispec, spectype)
       if (trim(spectype) == 'sulfate') then
          call self%get(bin_ndx,ispec,density=rho)
       end if
    end do

    ! specdens kg/m3 to g/cm3, convert from radius to diameter
    diam = 2._r8*((0.75*mass / pi  / (1.0e-3_r8*rho))**onethird)

  end function scav_diam

  !------------------------------------------------------------------------------
  ! adjust aerosol concentration tendencies to create larger sizes of aerosols
  ! during resuspension
  !------------------------------------------------------------------------------
  subroutine resuspension_resize(self, dcondt)
    class(carma_aerosol_properties), intent(in) :: self
    real(r8), intent(inout) :: dcondt(:)

    integer :: m

    ! move dcondt_prevap to larger bin
    do m = 1, self%ncnst_tot()
       if (self%ibl(m) /= m) then
          dcondt(self%ibl(m)) = dcondt(self%ibl(m)) + dcondt(m)
          dcondt(m) = 0._r8
       end if
    end do

  end subroutine resuspension_resize

  !------------------------------------------------------------------------------
  ! returns dust deposition fluxes rebinned to specified diameter limits
  !------------------------------------------------------------------------------
  subroutine rebin_bulk_fluxes(self, bulk_type, dep_fluxes, diam_edges, bulk_fluxes, &
                               error_code, error_string)

    class(carma_aerosol_properties), intent(in) :: self
    character(len=*),intent(in) :: bulk_type ! aerosol type to rebin
    real(r8), intent(in) :: dep_fluxes(:) ! kg/m2/sec
    real(r8), intent(in) :: diam_edges(:) ! meters
    real(r8), intent(out) :: bulk_fluxes(:) ! kg/m2/sec
       integer,  intent(out) :: error_code            ! error code (0 if no error)
       character(len=*), intent(out) :: error_string  ! error string

    real(r8) :: mflx, mflx_tot
    real(r8) :: rho, mass, frac, diam
    integer :: i, m,l,mm
    integer :: n_bulk_bins
    character(len=aero_name_len) :: spectype
    logical :: type_not_found

    error_code = 0
    error_string = ' '

    n_bulk_bins = size(bulk_fluxes)

    bulk_fluxes(:) = 0._r8
    type_not_found = .true.

    bin_loop: do m = 1,self%nbins()

       mflx_tot = 0._r8
       mflx = 0._r8

       species: do l = 1,self%nmasses(m)
          mm = self%indexer(m,l)

          if (l>self%nspecies(m)) then
             ! use mass flux for the entire bin (concentration element) if available
             ! -- override the total summed below
             mflx_tot = dep_fluxes(mm)
          else
             ! this sums up the total assuming all species are transported
             mflx_tot = mflx_tot + dep_fluxes(mm)

             call self%get(m,l,spectype=spectype)

             if (spectype==bulk_type) then
                ! get mass flux and density of the specified type
                mflx = dep_fluxes(mm)
                call self%get(m,l,density=rho) ! kg/m3
                type_not_found = .false.
             end if
          end if
       end do species

       if (mflx>0._r8 .and. mflx_tot>0._r8) then
          ! mass flux fraction
          frac = mflx/mflx_tot

          ! mass of the specified aerosol type
          mass = frac * bin_mass(m) ! kg

          ! diameter in meters
          diam = 2._r8*((0.75_r8*mass/pi/rho)**onethird)

          ! add the flux to the corresponding bulk bin
          blk_loop: do i = 1,n_bulk_bins-1
             if (diam>diam_edges(i) .and. diam<=diam_edges(i+1)) then
                bulk_fluxes(i) = bulk_fluxes(i) + mflx
                exit blk_loop
             end if
          end do blk_loop
       endif

    end do bin_loop

    if (type_not_found) then
       bulk_fluxes(:) = nan
       error_code = 1
       write(error_string,*) 'aerosol_properties::rebin_bulk_fluxes ERROR : ',trim(bulk_type),' not found'
    end if

  contains

    !---------------------------------------------------------------
    ! get mass of the specified bin in kg -- could be done at init time ...
    !---------------------------------------------------------------
    real(r8) function bin_mass(bin_ndx) ! (kg)
      use carma_intr, only: carma_get_bin_rmass, carma_get_group_by_name

      integer, intent(in) :: bin_ndx

      character(len=aero_name_len) :: bin_name, shortname
      integer :: ibin, igroup, rc, nchr
      real(r8) :: rmass

      call rad_cnst_get_info_by_bin(0, bin_ndx, bin_name=bin_name)

      nchr = len_trim(bin_name)-2
      shortname = bin_name(:nchr)

      call carma_get_group_by_name(shortname, igroup, rc)

      read(bin_name(nchr+1:),*) ibin

      call carma_get_bin_rmass(igroup, ibin, rmass, rc)
      bin_mass = rmass * 1.e-3_r8 ! g->kg

    end function bin_mass

  end subroutine rebin_bulk_fluxes

  !------------------------------------------------------------------------------
  ! Returns TRUE if bin is hydrophilic, otherwise FALSE
  !------------------------------------------------------------------------------
  logical function hydrophilic(self, bin_ndx)
    class(carma_aerosol_properties), intent(in) :: self
    integer, intent(in) :: bin_ndx ! bin number

    hydrophilic = .true.

  end function hydrophilic

end module carma_aerosol_properties_mod
