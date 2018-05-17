!---------------------------------------------------------------------
!	... compute chemical potential heating
!---------------------------------------------------------------------

module mo_cph

  use shr_kind_mod,only : r8 => shr_kind_r8
  use chem_mods,   only : ncph=>enthalpy_cnt, exotherm=>cph_enthalpy, cph_rid

  implicit none

  save

  !==============================================================
  !... Doug Kinnison, dkin@ucar.edu
  !
  !... Enthalpy Data are taken from Atkinson et al., 
  !    Evaluated kinetic and photochemical data for atmospheric
  !    chemistry: Volume I, Atmos. Chem. Phys., 4, 1461-1738.
  !... Heats of formation at 0K.
  !... Units: kJ mol-1
  !
  !... Exception to the Atkinson et al. reference  (@0K unless noted)
  !    (4), (5), (8), (9), (10, (11), (14), (15), (27), (28) at 298K  
  !    (7)  h + o2 -> oh + o2 is multiplied by 0.6 (Mlynczak) to represent
  !         AG loss of excited OH.
  !    (25) n2d + o2 -> no + o1d taken from Roble, UMLT, Johnson and Killeen.
  !    (26) n2d + o  -> n  + o   taken from Roble, UMLT, Johnson and Killeen.
  !    (30-41) Taken from Roble, UMLT, Johnson and Killeen Ed., Geophys. Mono. 87
  !==============================================================
  !... Enthalpy Data are specified in preprocessor input mechanism file
  !... F Vitt -- 29 Oct 2013
  !==============================================================

  private
  public :: cph, init_cph

  logical :: has_cph
  character(len=24) :: fldnames(ncph)
  logical, parameter :: debug = .false.

contains

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  subroutine init_cph

    use mo_chem_utls,   only : get_rxt_ndx, get_spc_ndx
    use cam_history,    only : addfld, add_default
    use chem_mods,      only : rxt_tag_lst, rxt_tag_map, rxt_tag_cnt
    use cam_abortutils, only : endrun

    implicit none

    character(len=64) :: longname
    integer :: i, n, tagndx

    has_cph = ncph > 0

    if (.not.has_cph) return

    if ( any(exotherm(:) == 0._r8) ) then
       call endrun('init_cph: Enthalpies for chemical heating must be specified in mechanism file')
    endif
    
    do i = 1,ncph

       findtagndx: do n = 1,rxt_tag_cnt
          if ( rxt_tag_map(n) == cph_rid(i) ) then
             tagndx = n
             exit findtagndx
          endif
       enddo findtagndx

       if (debug) then
          if ( i< 10 ) then
             write(fldnames(i),fmt='("CPH",i1)') i
          else if (i<100) then
             write(fldnames(i),fmt='("CPH",i2)') i
          else if (i<1000) then
             write(fldnames(i),fmt='("CPH",i3)') i
          endif
       else
          fldnames(i) = 'CPH_'//trim(rxt_tag_lst(tagndx))
       endif

       write( longname, fmt='(f12.6)') exotherm(i)
       longname = trim(adjustl(longname))//' kcal/mol chem pot heating rate for rxtn '//trim(rxt_tag_lst(tagndx))
       call addfld( fldnames(i), (/ 'lev' /), 'I', 'K/s', trim(longname) )
       if (debug) then
          call add_default( fldnames(i), 10, ' ' )
       endif
    enddo

    call addfld( 'QCP', (/ 'lev' /), 'I',   'K/s', 'chem pot heating rate' )

  end subroutine init_cph
  
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  subroutine cph( cph_tot, vmr, rxt, cp, mbar, kbot, ncol, lchnk )

    !-----------------------------------------------------------------------
    !      	... forms the chemical potential heating rates
    !-----------------------------------------------------------------------

    use chem_mods,     only : gas_pcnst, rxntot
    use ppgrid,        only : pver
    use cam_history,   only : outfld
    use mo_rxt_rates_conv, only : set_rates

    implicit none

    !-----------------------------------------------------------------------
    !     	... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in)   ::  ncol                                ! columns in chunck
    integer, intent(in)   ::  lchnk                               ! chunk index
    integer, intent(in)   ::  kbot                                ! bottom vert index
    real(r8), intent(in)  ::  rxt(ncol,pver,rxntot)               ! rxt rates (1/cm^3/s)
    real(r8), intent(in)  ::  cp(ncol,pver)                       ! specific heat capacity (J/K/kg)
    real(r8), intent(in)  ::  mbar(ncol,pver)                     ! atm mean mass (g/mole)
    real(r8), intent(in)  ::  vmr(ncol,pver,gas_pcnst)            ! concentrations (mol/mol)
    real(r8), intent(out) ::  cph_tot(ncol,pver)                  ! total heating (K/s)

    !-----------------------------------------------------------------------
    !     	... local variables
    !-----------------------------------------------------------------------
    integer  ::  i, k
    real(r8) ::  tmp(ncol,pver)
    real(r8) ::  cph_rate(ncol,pver,ncph)
    real(r8) ::  rxt_rates(ncol,pver,rxntot)

    if (.not.has_cph) return

    ! get the reaction rates from rate constants ...
    rxt_rates(:ncol,:,:) = rxt(:ncol,:,:)
    call set_rates( rxt_rates, vmr, ncol )

    ! compute corresponding chem heating rates ...
    cph_rate(:,:,:) = 0._r8
    tmp(:ncol,:) =  1._r8 / (1.e-6_r8*cp(:ncol,:)*mbar(:ncol,:))
    do i = 1,ncph
       cph_rate(:ncol,:,i) = tmp(:ncol,:) * rxt_rates(:ncol,:,cph_rid(i)) * exotherm(i)
    enddo

    ! compute total heating rate ...
    cph_tot(:,:) = 0._r8
    do k = 1,kbot
       do i = 1,ncol
          cph_tot(i,k) = sum( cph_rate(i,k,:) )
       end do
    end do

    ! output diagnostics
    do i = 1,ncph
       if ( exotherm(i)>0._r8) then
          call outfld( fldnames(i), cph_rate(:,:,i), ncol, lchnk )
       endif
    enddo

    call outfld( 'QCP', cph_tot(:,:), ncol, lchnk )

  end subroutine cph

end module mo_cph
