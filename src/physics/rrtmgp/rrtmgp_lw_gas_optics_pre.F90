!> \file rrtmgp_lw_gas_optics_pre.F90
!!

!> This module contains an init routine to initialize the k-distribution data
!! and functions needed to compute the longwave gaseous optical properties in RRTMGP.
!! It also contains a run routine to compute gas optics during the radiation subcycle
module rrtmgp_lw_gas_optics_pre
  use machine,                 only: kind_phys
  use ccpp_gas_concentrations, only: ty_gas_concs_ccpp

  implicit none

  public :: rrtmgp_lw_gas_optics_pre_run
contains

!> \section arg_table_rrtmgp_lw_gas_optics_pre_run Argument Table
!! \htmlinclude rrtmgp_lw_gas_optics_pre_run.html
!!
  subroutine rrtmgp_lw_gas_optics_pre_run(icall, rad_const_array, pmid, pint, nlay, ncol, gaslist, idxday, &
                  pverp, ktoprad, ktopcam, dolw, nradgas, gas_concs, errmsg, errflg)

    ! Set gas vmr for the gases in the radconstants module's gaslist.

    ! The memory management for the gas_concs object is internal. The arrays passed to it
    ! are copied to the internally allocated memory.

    integer,                     intent(in) :: icall      ! index of climate/diagnostic radiation call
    character(len=*),            intent(in) :: gaslist(:)
    integer,                     intent(in) :: nlay               ! number of layers in radiation calculation
    integer,                     intent(in) :: ncol      ! number of columns, ncol for LW, nday for SW
    integer,                     intent(in) :: pverp
    integer,                     intent(in) :: idxday(:)          ! indices of daylight columns in a chunk
    integer,                     intent(in) :: ktoprad
    integer,                     intent(in) :: ktopcam
    integer,                     intent(in) :: nradgas
    logical,                     intent(in) :: dolw
    real(kind_phys),             intent(in) :: pmid(:,:)
    real(kind_phys),             intent(in) :: pint(:,:)
    real(kind_phys),             intent(in) :: rad_const_array(:,:,:) ! array of radiatively-active constituent vmrs
                                                                      !  last index corresponds to index in gaslist

    type(ty_gas_concs_ccpp),     intent(inout) :: gas_concs  ! the result is VRM inside gas_concs
    character(len=*),            intent(out)   :: errmsg
    integer,                     intent(out)   :: errflg

    ! Local variables
    integer :: i, gas_idx, idx(ncol)
    integer :: istat
    real(kind_phys)              :: gas_mmr(ncol, nlay)
    real(kind_phys)              :: gas_vmr(ncol, nlay)
    real(kind_phys)              :: mmr(ncol, nlay)
    real(kind_phys) :: massratio
    character(len=256) :: alloc_errmsg

    ! For ozone profile above model
    real(kind_phys) :: P_top, P_int, P_mid, alpha, beta, a, b, chi_mid, chi_0, chi_eff

    character(len=*), parameter :: sub = 'rrtmgp_lw_gas_optics_pre_run'
    !----------------------------------------------------------------------------

    ! Set error variables
    errmsg = ''
    errflg = 0

    if (.not. dolw) then
       return
    end if

    ! set the column indices; just count for longwave
    do i = 1, ncol
       idx(i) = i
    end do

    do gas_idx = 1, nradgas

       ! grab mass mixing ratio of gas
       gas_mmr = rad_const_array(:,:,gas_idx)

       do i = 1, ncol
          mmr(i,ktoprad:) = gas_mmr(idx(i),ktopcam:)
       end do

       ! If an extra layer is being used, copy mmr from the top layer of CAM to the extra layer.
       if (nlay == pverp) then
          mmr(:,1) = mmr(:,2)
       end if

       ! special case: H2O is specific humidity, not mixing ratio. Use r = q/(1-q):
       if (gaslist(gas_idx) == 'H2O') then 
          mmr = mmr / (1._kind_phys - mmr)
       end if  

       ! convert MMR to VMR, multipy by ratio of dry air molar mas to gas molar mass.
       call get_molar_mass_ratio(gaslist(gas_idx), massratio, errmsg, errflg)
       if (errflg /= 0) then
          return
       end if
       gas_vmr = mmr * massratio

       ! special case: Setting O3 in the extra layer:
       ! 
       ! For the purpose of attenuating solar fluxes above the CAM model top, we assume that ozone 
       ! mixing decreases linearly in each column from the value in the top layer of CAM to zero at 
       ! the pressure level set by P_top. P_top has been set to 50 Pa (0.5 hPa) based on model tuning 
       ! to produce temperatures at the top of CAM that are most consistent with WACCM at similar pressure levels. 

       if ((gaslist(gas_idx) == 'O3') .and. (nlay == pverp)) then
          P_top = 50.0_kind_phys
          do i = 1, ncol
             P_int = pint(idx(i),1) ! pressure (Pa) at upper interface of CAM
             P_mid = pmid(idx(i),1) ! pressure (Pa) at midpoint of top layer of CAM
             alpha = log(P_int/P_top)
             beta =  log(P_mid/P_int)/log(P_mid/P_top)

             a =  ( (1._kind_phys + alpha) * exp(-alpha) - 1._kind_phys ) / alpha
             b =  1._kind_phys - exp(-alpha)
   
             if (alpha .gt. 0) then             ! only apply where top level is below 80 km
                chi_mid = gas_vmr(i,1)          ! molar mixing ratio of O3 at midpoint of top layer
                chi_0 = chi_mid /  (1._kind_phys + beta)
                chi_eff = chi_0 * (a + b)
                gas_vmr(i,1) = chi_eff
             end if
          end do
       end if

       errmsg = gas_concs%gas_concs%set_vmr(gaslist(gas_idx), gas_vmr)
       if (len_trim(errmsg) > 0) then
          errflg = 1
          return
       end if

!       deallocate(gas_vmr)
!       deallocate(mmr)
    end do

  end subroutine rrtmgp_lw_gas_optics_pre_run

!=========================================================================================

  subroutine get_molar_mass_ratio(gas_name, massratio, errmsg, errflg)

    ! return the molar mass ratio of dry air to gas based on gas_name

    character(len=*), intent(in)  :: gas_name
    real(kind_phys),  intent(out) :: massratio
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! local variables
    real(kind_phys), parameter :: amdw = 1.607793_kind_phys    ! Molecular weight of dry air / water vapor
    real(kind_phys), parameter :: amdc = 0.658114_kind_phys    ! Molecular weight of dry air / carbon dioxide
    real(kind_phys), parameter :: amdo = 0.603428_kind_phys    ! Molecular weight of dry air / ozone
    real(kind_phys), parameter :: amdm = 1.805423_kind_phys    ! Molecular weight of dry air / methane
    real(kind_phys), parameter :: amdn = 0.658090_kind_phys    ! Molecular weight of dry air / nitrous oxide
    real(kind_phys), parameter :: amdo2 = 0.905140_kind_phys   ! Molecular weight of dry air / oxygen
    real(kind_phys), parameter :: amdc1 = 0.210852_kind_phys   ! Molecular weight of dry air / CFC11
    real(kind_phys), parameter :: amdc2 = 0.239546_kind_phys   ! Molecular weight of dry air / CFC12

    character(len=*), parameter :: sub='get_molar_mass_ratio'
    !----------------------------------------------------------------------------
    ! Set error variables
    errmsg = ''
    errflg = 0

    select case (trim(gas_name)) 
       case ('H2O') 
          massratio = amdw
       case ('CO2')
          massratio = amdc
       case ('O3')
          massratio = amdo
       case ('CH4')
          massratio = amdm
       case ('N2O')
          massratio = amdn
       case ('O2')
          massratio = amdo2
       case ('CFC11')
          massratio = amdc1
       case ('CFC12')
          massratio = amdc2
       case default
          write(errmsg, '(a,a,a)') sub, ': Invalid gas: ', trim(gas_name)
          errflg = 1
    end select

end subroutine get_molar_mass_ratio


end module rrtmgp_lw_gas_optics_pre
