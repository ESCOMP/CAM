module mo_lightning
  !----------------------------------------------------------------------
  ! ... the lightning module
  !----------------------------------------------------------------------

  use shr_kind_mod,      only : r8 => shr_kind_r8
  use ppgrid,            only : begchunk, endchunk, pcols, pver
  use phys_grid,         only : ngcols_p
  use cam_abortutils,    only : endrun
  use cam_logfile,       only : iulog
  use spmd_utils,        only : masterproc, mpicom

  implicit none

  private
  public  :: lightning_inti
  public  :: lightning_no_prod
  public  :: prod_no

  save

  real(r8) :: csrf
  real(r8) :: factor = 0.1_r8              ! user-controlled scaling factor to achieve arbitrary no prod.
  real(r8) :: geo_factor                   ! grid cell area factor
  real(r8) :: vdist(16,3)                  ! vertical distribution of lightning
  real(r8), allocatable :: prod_no(:,:,:)
  real(r8), allocatable :: glob_prod_no_col(:,:)
  real(r8), allocatable :: flash_freq(:,:)
  integer :: no_ndx,xno_ndx
  logical :: has_no_lightning_prod = .false.

contains

  subroutine lightning_inti( lght_no_prd_factor )
    !----------------------------------------------------------------------
    !       ... initialize the lightning module
    !----------------------------------------------------------------------
    use mo_constants,  only : pi
    use ioFileMod,     only : getfil
    !use mo_chem_utls,  only : get_spc_ndx

    use cam_history,   only : addfld, add_default, horiz_only
    use dyn_grid,      only : get_dyn_grid_parm
    use phys_control,  only : phys_getopts

    implicit none

    !----------------------------------------------------------------------
    !	... dummy args
    !----------------------------------------------------------------------
    real(r8), intent(in) :: lght_no_prd_factor        ! lightning no production factor

    !!----------------------------------------------------------------------
    !!	... local variables
    !!----------------------------------------------------------------------
    !integer  :: astat
    !integer  :: ncid
    !integer  :: dimid
    !integer  :: vid
    !integer  :: gndx
    !integer  :: jl, ju
    !integer  :: nlat, nlon
    !integer  :: plon, plat
    !real(r8), allocatable :: lats(:)
    !real(r8), allocatable :: lons(:)
    !real(r8), allocatable :: landmask(:,:)
    !character(len=256) :: locfn
    !logical :: history_cesm_forcing

    !call phys_getopts( history_cesm_forcing_out = history_cesm_forcing )

    !no_ndx = get_spc_ndx('NO')
    !xno_ndx = get_spc_ndx('XNO')

    !has_no_lightning_prod = no_ndx>0 .or. xno_ndx>0
    !if (.not.has_no_lightning_prod) return

    !
    !if( lght_no_prd_factor /= 1._r8 ) then
    !   factor = factor*lght_no_prd_factor
    !end if


    !if (masterproc) write(iulog,*) 'lght_inti: lightning no production scaling factor = ',factor

    !!----------------------------------------------------------------------
    !!       ... vdist(kk,itype) = % of lightning nox between (kk-1) and (kk)
    !!           km for profile itype
    !!----------------------------------------------------------------------
    !vdist(:,1) = (/  3.0_r8, 3.0_r8, 3.0_r8, 3.0_r8, 3.4_r8, 3.5_r8, 3.6_r8, 4.0_r8, &       ! midlat cont
    !                 5.0_r8, 7.0_r8, 9.0_r8, 14.0_r8, 16.0_r8, 14.0_r8, 8.0_r8, 0.5_r8 /)
    !vdist(:,2) = (/  2.5_r8, 2.5_r8, 2.5_r8, 2.5_r8, 2.5_r8, 2.5_r8, 2.5_r8, 6.1_r8, &       ! trop marine
    !                 17.0_r8, 15.4_r8, 14.5_r8, 13.0_r8, 12.5_r8, 2.8_r8, 0.9_r8, 0.3_r8 /)
    !vdist(:,3) = (/  2.0_r8, 2.0_r8, 2.0_r8, 1.5_r8, 1.5_r8, 1.5_r8, 3.0_r8, 5.8_r8, &       ! trop cont
    !                 7.6_r8, 9.6_r8, 11.0_r8, 14.0_r8, 14.0_r8, 14.0_r8, 8.2_r8, 2.3_r8 /)

    !allocate( prod_no(pcols,pver,begchunk:endchunk),stat=astat )
    !if( astat /= 0 ) then
    !   write(iulog,*) 'lght_inti: failed to allocate prod_no; error = ',astat
    !   call endrun
    !end if
    !allocate( flash_freq(pcols,begchunk:endchunk),stat=astat )
    !if( astat /= 0 ) then
    !   write(iulog,*) 'lght_inti: failed to allocate flash_freq; error = ',astat
    !   call endrun
    !end if
    !allocate( glob_prod_no_col(pcols,begchunk:endchunk),stat=astat )
    !if( astat /= 0 ) then
    !   write(iulog,*) 'lght_inti: failed to allocate glob_prod_no_col; error = ',astat
    !   call endrun
    !end if
    !prod_no(:,:,:)   = 0._r8
    !flash_freq(:,:)  = 0._r8
    !geo_factor = ngcols_p/(4._r8*pi)


    !call addfld( 'LNO_COL_PROD', horiz_only,  'I', 'TG N/YR', 'lighting column NO source' )
    !call addfld( 'LNO_PROD',     (/ 'lev' /), 'I', '/cm3/s',  'lighting insitu NO source' )
    !call addfld( 'FLASHFRQ',     horiz_only,  'I', '1/MIN',   'lighting flash rate' )        ! flash frequency in grid box per minute (PPP)
    !call addfld( 'FLASHENGY',    horiz_only,  'I', '   ',     'lighting flash rate' )          ! flash frequency in grid box per minute (PPP)
    !call addfld( 'CLDHGT',       horiz_only,  'I', 'KM',      'cloud top height' )              ! cloud top height
    !call addfld( 'DCHGZONE',     horiz_only,  'I', 'KM',      'depth of discharge zone' )       ! depth of discharge zone
    !call addfld( 'CGIC',         horiz_only,  'I', 'RATIO',   'ratio of cloud-ground/intracloud discharges' ) ! ratio of cloud-ground/intracloud discharges

    !if ( history_cesm_forcing ) then
    !   call add_default('LNO_COL_PROD',1,' ')
    !endif

  end subroutine lightning_inti

  subroutine lightning_no_prod( state, pbuf2d,  cam_in )
    !----------------------------------------------------------------------
    !	... set no production from lightning
    !----------------------------------------------------------------------
    use physics_types,    only : physics_state
    
    use physics_buffer,   only : pbuf_get_index, physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    use physconst,        only : rga
    use phys_grid,        only : get_rlat_all_p, get_lat_all_p, get_lon_all_p, get_wght_all_p
    use cam_history,      only : outfld
    use camsrfexch,       only : cam_in_t
    use shr_reprosum_mod, only : shr_reprosum_calc
    !use mo_constants,  only : rearth, d2r
    implicit none

    !----------------------------------------------------------------------
    !	... dummy args
    !----------------------------------------------------------------------
    type(physics_state), intent(in) :: state(begchunk:endchunk) ! physics state
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(cam_in_t), intent(in) :: cam_in(begchunk:endchunk) ! physics state

    !----------------------------------------------------------------------
    !	... local variables
    !----------------------------------------------------------------------

    !----------------------------------------------------------------------
    ! 	... parameters to determine cg/ic ratio [price and rind, 1993]
    !----------------------------------------------------------------------

    if (.not.has_no_lightning_prod) return

    ! < === INSERT CALCULATION HERE === >

    !!--------------------------------------------------------------------------------
    !!       ... output lightning no production to history file
    !!--------------------------------------------------------------------------------
    !do c = begchunk,endchunk
    !   lchnk = state(c)%lchnk
    !   call outfld( 'LNO_PROD',     prod_no(:,:,c),        pcols, lchnk )
    !   call outfld( 'LNO_COL_PROD', glob_prod_no_col(:,c), pcols, lchnk )
    !   call outfld( 'FLASHFRQ',     flash_freq(:,c),       pcols, lchnk )
    !   call outfld( 'FLASHENGY',    flash_energy(:,c),     pcols, lchnk )
    !   call outfld( 'CLDHGT',       cldhgt(:,c),           pcols, lchnk )
    !   call outfld( 'DCHGZONE',     dchgzone(:,c),         pcols, lchnk )
    !   call outfld( 'CGIC',         cgic(:,c),             pcols, lchnk )
    !enddo

  end subroutine lightning_no_prod

end module mo_lightning
