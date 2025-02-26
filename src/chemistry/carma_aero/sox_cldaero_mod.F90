!----------------------------------------------------------------------------------
! CARMA implementation
!----------------------------------------------------------------------------------
module sox_cldaero_mod

  use shr_kind_mod,    only : r8 => shr_kind_r8
  use cam_abortutils,  only : endrun
  use ppgrid,          only : pcols, pver
  use mo_chem_utls,    only : get_spc_ndx
  use cldaero_mod,     only : cldaero_conc_t, cldaero_allocate, cldaero_deallocate
  use cam_logfile,     only : iulog
  use chem_mods,       only : adv_mass
  use physconst,       only : gravit
  use phys_control,    only : phys_getopts
  use cldaero_mod,     only : cldaero_uptakerate
  use chem_mods,       only : gas_pcnst
  use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_info_by_bin, rad_cnst_get_bin_props_by_idx

  implicit none
  private

  public :: sox_cldaero_init
  public :: sox_cldaero_create_obj
  public :: sox_cldaero_update
  public :: sox_cldaero_destroy_obj

  integer :: id_msa, id_h2so4, id_so2, id_h2o2, id_nh3

  real(r8), parameter :: small_value = 1.e-20_r8

  ! description of bin aerosols
  integer, public, protected :: nspec_max = 0
  integer, public, protected :: nbins = 0
  integer, public, protected, allocatable :: nspec(:)

  ! local indexing for bins
  integer, allocatable :: bin_idx(:,:) ! table for local indexing of modal aero number and mmr
  integer :: ncnst_tot                  ! total number of mode number conc + mode species

contains

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

  subroutine sox_cldaero_init

    integer :: l, m, ii
    logical :: history_aerosol      ! Output the MAM aerosol tendencies

    id_msa = get_spc_ndx( 'MSA' )
    id_h2so4 = get_spc_ndx( 'H2SO4' )
    id_so2 = get_spc_ndx( 'SO2' )
    id_h2o2 = get_spc_ndx( 'H2O2' )
    id_nh3 = get_spc_ndx( 'NH3' )

    if (id_h2so4<1 .or. id_so2<1 .or. id_h2o2<1) then
      call endrun('sox_cldaero_init:MAM mech does not include necessary species' &
                  //' -- should not invoke sox_cldaero_mod ')
    endif

   call phys_getopts( history_aerosol_out        = history_aerosol   )
    !
    !   add to history
    !

    ! get info about the modal aerosols
    ! get nbins

    call rad_cnst_get_info( 0, nbins=nbins)

    allocate( nspec(nbins) )

    do m = 1, nbins
       call rad_cnst_get_info_by_bin(0, m, nspec=nspec(m))
    end do
    ! add plus one to include number, total mmr and nspec
    nspec_max = maxval(nspec)

    ncnst_tot = nspec(1)
    do m = 2, nbins
      ncnst_tot = ncnst_tot + nspec(m)
    end do

   allocate(  bin_idx(nbins,nspec_max) )


   ! Local indexing compresses the mode and number/mass indicies into one index.
   ! This indexing is used by the pointer arrays used to reference state and pbuf
   ! fields.
   ! for CARMA we add number = 0, total mass = 1, and mass from each constituence into mm.
   ii = 0
   do m = 1, nbins
      do l = 1, nspec(m)    ! loop through species
         ii = ii + 1
         bin_idx(m,l) = ii
      end do
   end do


  end subroutine sox_cldaero_init

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
  function sox_cldaero_create_obj(cldfrc, qcw, lwc, cfact, ncol, loffset) result( conc_obj )

    real(r8), intent(in) :: cldfrc(:,:)
    real(r8), intent(in) :: qcw(:,:,:)
    real(r8), intent(in) :: lwc(:,:)
    real(r8), intent(in) :: cfact(:,:)
    integer,  intent(in) :: ncol
    integer,  intent(in) :: loffset

    real(r8) :: so4mmr(pcols,pver)

    type(cldaero_conc_t), pointer :: conc_obj

    character(len=32) :: spectype

    integer :: l,m
    integer :: i,k,mm

    ! local indexing for bins
    !integer, allocatable :: bin_idx(:,:) ! table for local indexing of modal aero number and mmr


    conc_obj => cldaero_allocate()

    do k = 1,pver
       do i = 1,ncol
          if( cldfrc(i,k) >0._r8) then
             conc_obj%xlwc(i,k) = lwc(i,k) *cfact(i,k) ! cloud water L(water)/L(air)
             conc_obj%xlwc(i,k) = conc_obj%xlwc(i,k) / cldfrc(i,k) ! liquid water in the cloudy fraction of cell
          else
             conc_obj%xlwc(i,k) = 0._r8
          endif
       enddo
    enddo

    conc_obj%no3c(:,:) = 0._r8
    conc_obj%nh4c(:,:) = 0._r8
    conc_obj%so4c(:,:) = 0._r8

    so4mmr(:,:) = 0._r8
    do k = 1,pver
       do i = 1,ncol
          do m = 1, nbins
            do l = 1, nspec(m)
                  mm = bin_idx(m, l)
                 call rad_cnst_get_bin_props_by_idx(0, m, l,spectype=spectype)
                 if (trim(spectype) == 'sulfate') then
                       so4mmr(i,k) =  so4mmr(i,k) +  qcw(i,k,mm)
                 end if
            end do
          end do
       end do
    end do
    conc_obj%so4c = so4mmr

  end function sox_cldaero_create_obj


!----------------------------------------------------------------------------------
! Update the mixing ratios
!----------------------------------------------------------------------------------
  subroutine sox_cldaero_update(  &
       state, ncol, lchnk, loffset, dtime, mbar, pdel, press, tfld, cldnum, cldfrc, cfact, xlwc, &
       delso4_hprxn, xh2so4, xso4, xso4_init, nh3g, hno3g, xnh3, xhno3, xnh4c,  xno3c, xmsa, xso2, xh2o2, qcw, qin, &
       aqso4, aqh2so4, aqso4_h2o2, aqso4_o3, aqso4_h2o2_3d, aqso4_o3_3d)

    use aerosol_properties_mod, only: aero_name_len
    use physics_types, only: physics_state
    use carma_intr, only: carma_get_group_by_name, carma_get_dry_radius

    ! args

    type(physics_state), intent(in) :: state     ! Physics state variables

    integer,  intent(in) :: ncol
    integer,  intent(in) :: lchnk ! chunk id
    integer,  intent(in) :: loffset

    real(r8), intent(in) :: dtime ! time step (sec)

    real(r8), intent(in) :: mbar(:,:) ! mean wet atmospheric mass ( amu )
    real(r8), intent(in) :: pdel(:,:)
    real(r8), intent(in) :: press(:,:)
    real(r8), intent(in) :: tfld(:,:)

    real(r8), intent(in) :: cldnum(:,:)
    real(r8), intent(in) :: cldfrc(:,:)
    real(r8), intent(in) :: cfact(:,:)
    real(r8), intent(in) :: xlwc(:,:)

    real(r8), intent(in) :: delso4_hprxn(:,:)
    real(r8), intent(in) :: xh2so4(:,:)
    real(r8), intent(in) :: xso4(:,:)
    real(r8), intent(in) :: xso4_init(:,:)
    real(r8), intent(in) :: nh3g(:,:)
    real(r8), intent(in) :: hno3g(:,:)
    real(r8), intent(in) :: xnh3(:,:)
    real(r8), intent(in) :: xhno3(:,:)
    real(r8), intent(in) :: xnh4c(:,:)
    real(r8), intent(in) :: xmsa(:,:)
    real(r8), intent(in) :: xso2(:,:)
    real(r8), intent(in) :: xh2o2(:,:)
    real(r8), intent(in) :: xno3c(:,:)

    real(r8), intent(inout) :: qcw(:,:,:) ! cloud-borne aerosol (vmr)  vmrcw(ncol,pver,ncnst_tot)
    real(r8), intent(inout) :: qin(:,:,:) ! xported species ( vmr )

    real(r8), intent(out) :: aqso4(:,:)                   ! aqueous phase chemistry
    real(r8), intent(out) :: aqh2so4(:,:)                 ! aqueous phase chemistry
    real(r8), intent(out) :: aqso4_h2o2(:)                ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
    real(r8), intent(out) :: aqso4_o3(:)                  ! SO4 aqueous phase chemistry due to O3 (kg/m2)
    real(r8), intent(out), optional :: aqso4_h2o2_3d(:,:)                ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
    real(r8), intent(out), optional :: aqso4_o3_3d(:,:)                  ! SO4 aqueous phase chemistry due to O3 (kg/m2)

    ! local vars ...
    real(r8) :: dryr(pcols,pver)   ! CARMA dry radius in cm
    real(r8) :: rho(pcols,pver)   !
    real(r8) :: dryr_n(nbins,ncol,pver)   ! CARMA dry radius in cm
    real(r8) :: dqdt_aqso4(ncol,pver,ncnst_tot), &
         dqdt_aqh2so4(ncol,pver,ncnst_tot), &
         dqdt_aqhprxn(ncol,pver), dqdt_aqo3rxn(ncol,pver)

    real(r8) :: faqgain_so4(nbins)
    real(r8) :: wt_mass(nbins)

    real(r8) :: delso4_o3rxn, &
         dso4dt_aqrxn, dso4dt_hprxn, &
         dso4dt_gasuptk, dmsadt_gasuptk_toso4, &
         dqdt_aq, dqdt_wr, dqdt

    real(r8) :: fwetrem, uptkrate

    integer :: l, n, mm
    integer :: ntot_msa_c

    integer :: i,k
    real(r8) :: xl
    real(r8) :: wt_sum
    real(r8) :: specmw_so4_amode

    character(len=32) :: spectype

    character(len=*), parameter :: subname = 'sox_cldaero_update'
    character(len=aero_name_len) :: bin_name, shortname
    integer :: igroup, ibin, rc, nchr

    ! make sure dqdt is zero initially, for budgets
    dqdt_aqso4(:,:,:) = 0.0_r8
    dqdt_aqh2so4(:,:,:) = 0.0_r8
    dqdt_aqhprxn(:,:) = 0.0_r8
    dqdt_aqo3rxn(:,:) = 0.0_r8
    dryr_n(:,:,:) = 0.0_r8

    ntot_msa_c = 0.0_r8
    aqso4 = 0.0_r8
    aqh2so4 = 0.0_r8
    aqso4_h2o2 = 0.0_r8
    aqso4_o3 = 0.0_r8

    do n = 1, nbins
       call rad_cnst_get_info_by_bin(0, n, nspec=nspec(n), bin_name=bin_name)


       nchr = len_trim(bin_name)-2
       shortname = bin_name(:nchr)

       call carma_get_group_by_name(shortname, igroup, rc)
       if (rc/=0) then
          call endrun(subname//': ERROR in carma_get_group_by_name')
       end if

       read(bin_name(nchr+1:),*) ibin

       call carma_get_dry_radius(state, igroup, ibin, dryr, rho, rc)
       if (rc/=0) then
          call endrun(subname//': ERROR in carma_get_dry_radius')
       end if

       dryr(:ncol,:) = dryr(:ncol,:)*1.e2_r8 ! cm

       if (index(bin_name,'MXAER')>0) then
          dryr_n(n,:ncol,:) = dryr(:ncol,:)
       end if
    end do

    lev_loop: do k = 1,pver
       col_loop: do i = 1,ncol
          cloud: if (cldfrc(i,k) >= 1.0e-5_r8) then
             xl = xlwc(i,k)

             if (xl .ge. 1.e-8_r8) then !! when cloud is present

                delso4_o3rxn = xso4(i,k) - xso4_init(i,k)

                ! the factors are proportional to the activated particle MR for each
                ! bin, which is the MR of cloud drops "associated with" the mode
                ! thus we are assuming the cloud drop size is independent of the
                ! associated aerosol mode properties (i.e., drops associated with
                ! Aitken and coarse sea-salt particles are same size)
                !
                ! qnum_c(n) = activated particle number MR for mode n (these are just
                ! used for partitioning among modes, so don't need to divide by cldfrc)

                !faqgain_so4(n) = fraction of total so4_c gain going to mode n
                wt_sum = 0._r8
                wt_mass(:) = 0._r8
                faqgain_so4(:) = 0.0_r8
                do n = 1, nbins
                   if (dryr_n(n,i,k) > 0._r8) then
                      wt_mass(n) = delso4_o3rxn / dryr_n(n,i,k) / dryr_n(n,i,k)
                      wt_sum = wt_sum +  wt_mass(n)
                   end if
                end do
                do n = 1, nbins
                   if (wt_mass(n) > 0._r8) then
                      faqgain_so4(n) = wt_mass(n)/wt_sum
                   end if
                end do

                uptkrate = cldaero_uptakerate( xl, cldnum(i,k), cfact(i,k), cldfrc(i,k), tfld(i,k),  press(i,k) )
                ! average uptake rate over dtime
                uptkrate = (1.0_r8 - exp(-min(100._r8,dtime*uptkrate))) / dtime

                dso4dt_gasuptk = xh2so4(i,k) * uptkrate

                ! if no modes have msa aerosol, then "rename" scavenged msa gas to so4
                dmsadt_gasuptk_toso4 = 0.0_r8

                !-----------------------------------------------------------------------
                ! now compute TMR tendencies
                ! this includes the above aqueous so2 chemistry AND
                ! the uptake of highly soluble aerosol precursor gases (h2so4, msa, ...)
                ! AND the wetremoval of dissolved, unreacted so2 and h2o2

                dso4dt_aqrxn = (delso4_o3rxn + delso4_hprxn(i,k)) / dtime
                dso4dt_hprxn = delso4_hprxn(i,k) / dtime
                !write(iulog,*) 'dso4dt_aqrxn ',dso4dt_aqrxn

                ! fwetrem = fraction of in-cloud-water material that is wet removed
                ! fwetrem = max( 0.0_r8, (1.0_r8-exp(-min(100._r8,dtime*clwlrat(i,k)))) )
                fwetrem = 0.0_r8 ! don't have so4 & msa wet removal here

                ! compute TMR tendencies for so4, not done currently for msa aerosol-in-cloud-water
                do n = 1, nbins
                   do l = 1, nspec(n)
                      mm = bin_idx(n, l)
                       call rad_cnst_get_bin_props_by_idx(0, n, l,spectype=spectype)
                       if (trim(spectype) == 'sulfate') then
                         if (faqgain_so4(n) .gt. 0.0_r8) then
                          dqdt_aqso4(i,k,mm) = faqgain_so4(n)*dso4dt_aqrxn*cldfrc(i,k)

                          dqdt_aqh2so4(i,k,mm) = faqgain_so4(n)* &
                            (dso4dt_gasuptk + dmsadt_gasuptk_toso4)*cldfrc(i,k)
                          dqdt_aq = dqdt_aqso4(i,k,mm) + dqdt_aqh2so4(i,k,mm)
                          dqdt_wr = -fwetrem*dqdt_aq
                          dqdt= dqdt_aq + dqdt_wr
                          !write(iulog,*) 'qcw(i,k,mm) before ', m, qcw(i,k,mm)
                          qcw(i,k,mm) = qcw(i,k,mm) + dqdt*dtime
                          !write(iulog,*) 'qcw(i,k,mm) after', m, qcw(i,k,mm)
                         end if
                       end if
                   end do
                 end do


                ! For gas species, tendency includes
                ! reactive uptake to cloud water that essentially transforms the gas to
                ! a different species. Wet removal associated with this is applied
                ! to the "new" species (e.g., so4_c) rather than to the gas.
                ! wet removal of the unreacted gas that is dissolved in cloud water.
                ! Need to multiply both these parts by cldfrc

                ! h2so4 (g) & msa (g)
                qin(i,k,id_h2so4) = qin(i,k,id_h2so4) - dso4dt_gasuptk * dtime * cldfrc(i,k)

                ! so2 -- the first order loss rate for so2 is frso2_c*clwlrat(i,k)
                ! fwetrem = max( 0.0_r8, (1.0_r8-exp(-min(100._r8,dtime*frso2_c*clwlrat(i,k)))) )
                fwetrem = 0.0_r8 ! don't include so2 wet removal here

                dqdt_wr = -fwetrem*xso2(i,k)/dtime*cldfrc(i,k)
                dqdt_aq = -dso4dt_aqrxn*cldfrc(i,k)
                dqdt = dqdt_aq + dqdt_wr
                dqdt = dqdt_aq
                qin(i,k,id_so2) = qin(i,k,id_so2) + dqdt * dtime
                qin(i,k,id_so2) =  MAX( qin(i,k,id_so2),    small_value )

                ! h2o2 -- the first order loss rate for h2o2 is frh2o2_c*clwlrat(i,k)
                ! fwetrem = max( 0.0_r8, (1.0_r8-exp(-min(100._r8,dtime*frh2o2_c*clwlrat(i,k)))) )
                fwetrem = 0.0_r8 ! don't include h2o2 wet removal here

                dqdt_wr = -fwetrem*xh2o2(i,k)/dtime*cldfrc(i,k)
                dqdt_aq = -dso4dt_hprxn*cldfrc(i,k)
                dqdt = dqdt_aq + dqdt_wr
                dqdt = dqdt_aq
                qin(i,k,id_h2o2) = qin(i,k,id_h2o2) + dqdt * dtime
                qin(i,k,id_h2o2) =  MAX( qin(i,k,id_h2o2),    small_value )

                ! for SO4 from H2O2/O3 budgets
                dqdt_aqhprxn(i,k) = dso4dt_hprxn*cldfrc(i,k)
                dqdt_aqo3rxn(i,k) = (dso4dt_aqrxn - dso4dt_hprxn)*cldfrc(i,k)

            endif !! when cloud is present
          endif cloud
       enddo col_loop
    enddo lev_loop

    !==============================================================
    ! ... Update the mixing ratios
    !==============================================================

    ! diagnostics

    specmw_so4_amode = 96.0_r8
      do n = 1, nbins
        ! while looking through all species, only dqdt_aqso4 from sulfates  is gt zero
        do l = 1, nspec(n)
           mm = bin_idx(n, l)
           aqso4(:,n)=0._r8
            do k=1,pver
               do i=1,ncol
                  aqso4(i,n)=aqso4(i,n)+dqdt_aqso4(i,k,mm)*specmw_so4_amode/mbar(i,k) &
                       *pdel(i,k)/gravit ! kg/m2/s
               enddo
            enddo

            aqh2so4(:,n)=0._r8
            do k=1,pver
               do i=1,ncol
                  aqh2so4(i,n)=aqh2so4(i,n)+dqdt_aqh2so4(i,k,mm)*specmw_so4_amode/mbar(i,k) &
                       *pdel(i,k)/gravit ! kg/m2/s
               enddo
            enddo
         end do
      end do

    aqso4_h2o2(:) = 0._r8
    do k=1,pver
       do i=1,ncol
           aqso4_h2o2(i)=aqso4_h2o2(i)+dqdt_aqhprxn(i,k)*specmw_so4_amode/mbar(i,k) &
                   *pdel(i,k)/gravit ! kg SO4 /m2/s
       enddo
    enddo

    if (present(aqso4_h2o2_3d)) then
        aqso4_h2o2_3d(:,:) = 0._r8
        do k=1,pver
           do i=1,ncol
              aqso4_h2o2_3d(i,k)=dqdt_aqhprxn(i,k)*specmw_so4_amode/mbar(i,k) &
                                 *pdel(i,k)/gravit ! kg SO4 /m2/s
           enddo
        enddo
    end if

    aqso4_o3(:)=0._r8
    do k=1,pver
        do i=1,ncol
           aqso4_o3(i)=aqso4_o3(i)+dqdt_aqo3rxn(i,k)*specmw_so4_amode/mbar(i,k) &
                   *pdel(i,k)/gravit ! kg SO4 /m2/s
        enddo
    enddo

    if (present(aqso4_o3_3d)) then
        aqso4_o3_3d(:,:)=0._r8
        do k=1,pver
           do i=1,ncol
              aqso4_o3_3d(i,k)=dqdt_aqo3rxn(i,k)*specmw_so4_amode/mbar(i,k) &
                               *pdel(i,k)/gravit ! kg SO4 /m2/s
           enddo
        enddo
    end if

  end subroutine sox_cldaero_update

  !----------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------
  subroutine sox_cldaero_destroy_obj( conc_obj )
    type(cldaero_conc_t), pointer :: conc_obj

    call cldaero_deallocate( conc_obj )

  end subroutine sox_cldaero_destroy_obj

end module sox_cldaero_mod
