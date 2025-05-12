module mo_lightning
  !----------------------------------------------------------------------
  ! ... the lightning module
  !----------------------------------------------------------------------

  use shr_kind_mod,      only : r8 => shr_kind_r8
  use ppgrid,            only : begchunk, endchunk, pcols, pver
  use phys_grid,         only : ngcols_p => num_global_phys_cols
  use cam_abortutils,    only : endrun
  use cam_logfile,       only : iulog
  use spmd_utils,        only : masterproc, mpicom

  use physics_buffer, only : pbuf_get_index, physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
  use physics_buffer, only : pbuf_add_field, pbuf_set_field, dtype_r8

  implicit none

  private

  public  :: lightning_readnl
  public  :: lightning_register
  public  :: lightning_init
  public  :: lightning_no_prod
  public  :: prod_no

  real(r8),protected, allocatable :: prod_no(:,:,:)

  real(r8) :: factor = 0.1_r8              ! user-controlled scaling factor to achieve arbitrary no prod.
  real(r8) :: geo_factor = -huge(1._r8)    ! grid cell area factor
  real(r8), allocatable :: vdist(:,:)      ! vertical distribution of lightning

  logical :: calc_nox_prod = .false.
  logical :: calc_lightning = .false.

  integer :: flsh_frq_ndx = -1
  integer :: cldtop_ndx = -1, cldbot_ndx = -1

  ! namelist parameter
  real(r8) :: lght_no_prd_factor = -huge(1._r8)

contains

  !-------------------------------------------------------------------------
  ! Read namelist options
  !-------------------------------------------------------------------------
  subroutine lightning_readnl(nlfile)
    use namelist_utils, only : find_group_name
    use spmd_utils,     only : mpicom, masterprocid, mpi_real8, mpi_success

    character(len=*), intent(in)  :: nlfile  ! filepath for file containing namelist input

    integer                       :: unitn, ierr
    character(len=*), parameter   :: subname = 'lightning_readnl'

    ! ===================
    ! Namelist definition
    ! ===================
    namelist /lightning_nl/ lght_no_prd_factor

    ! =============
    ! Read namelist
    ! =============
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'lightning_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, lightning_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! ============================
    ! Broadcast namelist variables
    ! ============================
    call mpi_bcast(lght_no_prd_factor, 1, mpi_real8, masterprocid, mpicom, ierr)
    if (ierr/=mpi_success) then
       call endrun(subname//': MPI_BCAST ERROR: lght_no_prd_factor')
    end if

    if (masterproc) then
       write(iulog,*) subname,' lght_no_prd_factor: ',lght_no_prd_factor
    end if

    if( lght_no_prd_factor /= 1._r8 ) then
       factor = factor*lght_no_prd_factor
    end if

  end subroutine lightning_readnl

  !-------------------------------------------------------------------------
  ! register phys buffer field for cloud to ground lightning flash frequency
  ! to pass to the mediator for land model
  !-------------------------------------------------------------------------
  subroutine lightning_register()

    call pbuf_add_field('LGHT_FLASH_FREQ','global',dtype_r8,(/pcols/),flsh_frq_ndx) ! per minute

  end subroutine lightning_register

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  subroutine lightning_init( pbuf2d )
    !----------------------------------------------------------------------
    !       ... initialize the lightning module
    !----------------------------------------------------------------------
    use mo_constants,  only : pi

    use cam_history,   only : addfld, add_default, horiz_only
    use phys_control,  only : phys_getopts
    use time_manager,  only : is_first_step

    !----------------------------------------------------------------------
    !	... dummy args
    !----------------------------------------------------------------------
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !----------------------------------------------------------------------
    !	... local variables
    !----------------------------------------------------------------------
    integer  :: astat, err
    logical :: history_cesm_forcing
    character(len=*),parameter :: prefix = 'lightning_init: '

    cldtop_ndx = pbuf_get_index('CLDTOP',errcode=err)
    cldbot_ndx = pbuf_get_index('CLDBOT',errcode=err)
    calc_lightning = cldtop_ndx>0 .and. cldbot_ndx>0

    if (.not.calc_lightning) return

    calc_nox_prod = lght_no_prd_factor>0._r8

    if (calc_nox_prod) then

       if (masterproc) write(iulog,*) prefix,'lightning no production scaling factor = ',factor

       !----------------------------------------------------------------------
       !       ... vdist(kk,itype) = % of lightning nox between (kk-1) and (kk)
       !           km for profile itype
       !----------------------------------------------------------------------
       allocate(vdist(16,3),stat=astat)
       if( astat /= 0 ) then
          write(iulog,*) prefix,'failed to allocate vdist; error = ',astat
          call endrun(prefix//'failed to allocate vdist')
       end if
       vdist(:,1) = (/  3.0_r8, 3.0_r8, 3.0_r8, 3.0_r8, 3.4_r8, 3.5_r8, 3.6_r8, 4.0_r8, &       ! midlat cont
            5.0_r8, 7.0_r8, 9.0_r8, 14.0_r8, 16.0_r8, 14.0_r8, 8.0_r8, 0.5_r8 /)
       vdist(:,2) = (/  2.5_r8, 2.5_r8, 2.5_r8, 2.5_r8, 2.5_r8, 2.5_r8, 2.5_r8, 6.1_r8, &       ! trop marine
            17.0_r8, 15.4_r8, 14.5_r8, 13.0_r8, 12.5_r8, 2.8_r8, 0.9_r8, 0.3_r8 /)
       vdist(:,3) = (/  2.0_r8, 2.0_r8, 2.0_r8, 1.5_r8, 1.5_r8, 1.5_r8, 3.0_r8, 5.8_r8, &       ! trop cont
            7.6_r8, 9.6_r8, 11.0_r8, 14.0_r8, 14.0_r8, 14.0_r8, 8.2_r8, 2.3_r8 /)

       allocate( prod_no(pcols,pver,begchunk:endchunk),stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) prefix, 'failed to allocate prod_no; error = ',astat
          call endrun(prefix//'failed to allocate prod_no')
       end if
       geo_factor = ngcols_p/(4._r8*pi)

       call addfld( 'LNO_COL_PROD', horiz_only,  'I', 'Tg N yr-1', 'lightning column NO source' )
       call addfld( 'LNO_PROD',     (/ 'lev' /), 'I', 'molecules/cm3/s', 'lightning insitu NO source' )
       call addfld( 'FLASHENGY',    horiz_only,  'I', 'J', 'lightning flash energy' ) ! flash energy

       call phys_getopts( history_cesm_forcing_out = history_cesm_forcing )
       if ( history_cesm_forcing ) then
          call add_default('LNO_COL_PROD',1,' ')
       endif

       if (is_first_step()) then
          call pbuf_set_field(pbuf2d, flsh_frq_ndx, 0.0_r8)
       endif

    endif

    call addfld( 'FLASHFRQ', horiz_only,  'I', 'min-1', 'lightning flash rate' )    ! flash frequency in grid box per minute (PPP)
    call addfld( 'CLDHGT',   horiz_only,  'I', 'km',    'cloud top height' )        ! cloud top height
    call addfld( 'DCHGZONE', horiz_only,  'I', 'km',    'depth of discharge zone' ) ! depth of discharge zone
    call addfld( 'CGIC',     horiz_only,  'I', '1',     'ratio of cloud-ground/intracloud discharges' ) ! ratio of cloud-ground/intracloud discharges
    call addfld( 'LGHTNG_CLD2GRND', horiz_only,  'I', 'min-1', 'clound-to-ground lightning flash rate') ! clound to ground flash frequency

  end subroutine lightning_init

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  subroutine lightning_no_prod( state, pbuf2d,  cam_in )
    !----------------------------------------------------------------------
    !	... set no production from lightning
    !----------------------------------------------------------------------
    use physics_types,    only : physics_state
    use physconst,        only : rga
    use phys_grid,        only : get_rlat_all_p, get_wght_all_p
    use cam_history,      only : outfld
    use camsrfexch,       only : cam_in_t
    use shr_reprosum_mod, only : shr_reprosum_calc
    use mo_constants,     only : rearth, d2r

    !----------------------------------------------------------------------
    !	... dummy args
    !----------------------------------------------------------------------
    type(physics_state), intent(in) :: state(begchunk:endchunk) ! physics state
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(cam_in_t), intent(in) :: cam_in(begchunk:endchunk) ! physics state

    !----------------------------------------------------------------------
    !	... local variables
    !----------------------------------------------------------------------
    real(r8), parameter :: land   = 1._r8
    real(r8), parameter :: secpyr = 365._r8 * 8.64e4_r8

    integer :: cldtind             ! level index for cloud top
    integer :: cldbind             ! level index for cloud base > 273k
    integer :: k, kk, zlow_ind, zhigh_ind, itype
    real(r8) :: glob_flashfreq     ! global flash frequency [s-1]
    real(r8) :: glob_noprod        ! global rate of no production [as tgn/yr]
    real(r8) :: frac_sum           ! work variable
    real(r8) :: zlow
    real(r8) :: zhigh
    real(r8) :: zlow_scal
    real(r8) :: zhigh_scal
    real(r8) :: fraction
    real(r8) :: dchgz
    real(r8) :: dchgzone(pcols,begchunk:endchunk)           ! depth of discharge zone [km]
    real(r8) :: cldhgt(pcols,begchunk:endchunk)             ! cloud top height [km]
    real(r8) :: cgic(pcols,begchunk:endchunk)               ! cloud-ground/intracloud discharge ratio
    real(r8) :: flash_energy(pcols,begchunk:endchunk)       ! energy of flashes per second
    real(r8) :: prod_no_col(pcols,begchunk:endchunk)        ! global no production rate for diagnostics
    real(r8) :: wrk, wrk1, wrk2(1)
    integer  :: icol                                    ! column index
    integer  :: ncol                                    ! columns per chunk
    integer  :: lchnk                                   ! chunk index
    real(r8),pointer :: cldtop(:)                       ! cloud top level index
    real(r8),pointer :: cldbot(:)                       ! cloud bottom level index
    real(r8) :: zmid(pcols,pver)                        ! geopot height above surface at midpoints (m)
    real(r8) :: zint(pcols,pver+1,begchunk:endchunk)    ! geopot height above surface at interfaces (m)
    real(r8) :: zsurf(pcols)                            ! geopot height above surface at interfaces (m)
    real(r8) :: rlats(pcols)                            ! column latitudes in chunks
    real(r8) :: wght(pcols)

    real(r8) :: glob_prod_no_col(pcols,begchunk:endchunk)
    real(r8) :: flash_freq(pcols,begchunk:endchunk)

    !----------------------------------------------------------------------
    ! 	... parameters to determine cg/ic ratio [price and rind, 1993]
    !----------------------------------------------------------------------
    real(r8), parameter  :: ca = .021_r8
    real(r8), parameter  :: cb = -.648_r8
    real(r8), parameter  :: cc = 7.49_r8
    real(r8), parameter  :: cd = -36.54_r8
    real(r8), parameter  :: ce = 64.09_r8
    real(r8), parameter  :: t0 = 273._r8
    real(r8), parameter  :: m2km  = 1.e-3_r8
    real(r8), parameter  :: km2cm = 1.e5_r8
    real(r8), parameter  :: lat25 = 25._r8*d2r      ! 25 degrees latitude in radians

    real(r8) :: flash_freq_land, flash_freq_ocn
    real(r8), pointer :: cld2grnd_flash_freq(:)

    if (.not.calc_lightning) return

    nullify(cld2grnd_flash_freq)

    !----------------------------------------------------------------------
    !	... initialization
    !----------------------------------------------------------------------

    flash_freq(:,:)       = 0._r8
    cldhgt(:,:)           = 0._r8
    dchgzone(:,:)         = 0._r8
    cgic(:,:)             = 0._r8
    flash_energy(:,:)     = 0._r8

    if (calc_nox_prod) then
       prod_no(:,:,:)     = 0._r8
       prod_no_col(:,:)   = 0._r8
       glob_prod_no_col(:,:) = 0._r8
    end if

    !--------------------------------------------------------------------------------
    !	... estimate flash frequency and resulting no emissions
    !           [price, penner, prather, 1997 (jgr)]
    !    lightning only occurs in convective clouds with a discharge zone, i.e.
    !    an altitude range where liquid water, ice crystals, and graupel coexist.
    !    we test this by examining the temperature at the cloud base.
    !    it is assumed that only one thunderstorm occurs per grid box, and its
    !    flash frequency is determined by the maximum cloud top height (not the
    !    depth of the discharge zone). this is somewhat speculative but yields
    !    reasonable results.
    !
    !       the cg/ic ratio is determined by an empirical formula from price and
    !    rind [1993]. the average energy of a cg flash is estimated as 6.7e9 j,
    !    and the average energy of a ic flash is assumed to be 1/10 of that value.
    !       the no production rate is assumed proportional to the discharge energy
    !    with 1e17 n atoms per j. the total number of n atoms is then distributed
    !    over the complete column of grid boxes.
    !--------------------------------------------------------------------------------
    Chunk_loop : do lchnk = begchunk,endchunk
       ncol  = state(lchnk)%ncol
       call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk), flsh_frq_ndx, cld2grnd_flash_freq )
       call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk), cldtop_ndx, cldtop )
       call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk), cldbot_ndx, cldbot )
       zsurf(:ncol) = state(lchnk)%phis(:ncol)*rga
       call get_wght_all_p(lchnk, pcols, wght)

       do k = 1,pver
          zmid(:ncol,k)   = state(lchnk)%zm(:ncol,k) + zsurf(:ncol)
          zint(:ncol,k,lchnk) = state(lchnk)%zi(:ncol,k) + zsurf(:ncol)
       end do
       zint(:ncol,pver+1,lchnk) = state(lchnk)%zi(:ncol,pver+1) + zsurf(:ncol)

       cld2grnd_flash_freq(:) = 0.0_r8

       col_loop : do icol = 1,ncol
          !--------------------------------------------------------------------------------
          ! 	... find cloud top and bottom level above 273k
          !--------------------------------------------------------------------------------
          cldtind = nint( cldtop(icol) )
          cldbind = nint( cldbot(icol) )
          do
             if( cldbind <= cldtind .or. state(lchnk)%t(icol,cldbind) < t0 ) then
                exit
             end if
             cldbind = cldbind - 1
          end do
          cloud_layer : if( cldtind < pver .and. cldtind > 0 .and. cldtind < cldbind ) then
             !--------------------------------------------------------------------------------
             !       ... compute cloud top height and depth of charging zone
             !--------------------------------------------------------------------------------
             cldhgt(icol,lchnk)   = m2km*max( 0._r8,zint(icol,cldtind,lchnk) )
             dchgz = cldhgt(icol,lchnk) - m2km*zmid(icol,cldbind)
             dchgzone(icol,lchnk) = dchgz
             !--------------------------------------------------------------------------------
             !       ... compute flash frequency for given cloud top height
             !           (flashes storm^-1 min^-1)
             !--------------------------------------------------------------------------------
             flash_freq_land = 3.44e-5_r8 * cldhgt(icol,lchnk)**4.9_r8
             flash_freq_ocn  = 6.40e-4_r8 * cldhgt(icol,lchnk)**1.7_r8
             flash_freq(icol,lchnk) = cam_in(lchnk)%landfrac(icol)*flash_freq_land + &
                               cam_in(lchnk)%ocnfrac(icol) *flash_freq_ocn

             !--------------------------------------------------------------------------------
             !   cgic = proportion of cloud-to-ground flashes
             ! NOx from lightning 1. Global distribution based on lightning physics, C Price et al
             ! JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 102, NO. D5, PAGES 5929-5941, MARCH 20, 1997
             ! (https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/96JD03504)
             ! eq 14
             !--------------------------------------------------------------------------------
             cgic(icol,lchnk) = 1._r8/((((ca*dchgz + cb)*dchgz + cc) *dchgz + cd)*dchgz + ce)
             if( dchgz < 5.5_r8 ) then
                cgic(icol,lchnk) = 0._r8
             else if( dchgz > 14._r8 ) then
                cgic(icol,lchnk) = .02_r8
             end if

             cld2grnd_flash_freq(icol) = cam_in(lchnk)%landfrac(icol)*flash_freq_land*cgic(icol,lchnk) ! cld-to-grnd flash frq (per min)

             if (calc_nox_prod) then
                !--------------------------------------------------------------------------------
                !       ... compute flash energy (cg*6.7e9 + ic*6.7e8)
                !           and convert to total energy per second
                !           set ic = cg
                !--------------------------------------------------------------------------------
                flash_energy(icol,lchnk) = 6.7e9_r8 * flash_freq(icol,lchnk)/60._r8
                !--------------------------------------------------------------------------------
                !       ... LKE Aug 23, 2005: scale production to account for different grid
                !           box sizes. This requires a reduction in the overall fudge factor
                !           (e.g., from 1.2 to 0.5)
                !--------------------------------------------------------------------------------
                flash_energy(icol,lchnk) =  flash_energy(icol,lchnk) * wght(icol) * geo_factor
                !--------------------------------------------------------------------------------
                ! 	... compute number of n atoms produced per second
                !           and convert to n atoms per second per cm2 and apply fudge factor
                !--------------------------------------------------------------------------------
                prod_no_col(icol,lchnk) = 1.e17_r8*flash_energy(icol,lchnk)/(1.e4_r8*rearth*rearth*wght(icol)) * factor

                !--------------------------------------------------------------------------------
                ! 	... compute global no production rate in tgn/yr:
                !           tgn per second: * 14.00674 * 1.65979e-24 * 1.e-12
                !             nb: 1.65979e-24 = 1/avo
                !           tgn per year: * secpyr
                !--------------------------------------------------------------------------------
                glob_prod_no_col(icol,lchnk) = 1.e17_r8*flash_energy(icol,lchnk) &
                     * 14.00674_r8 * 1.65979e-24_r8 * 1.e-12_r8 * secpyr * factor
             end if
          end if cloud_layer
       end do Col_loop

       call outfld( 'LGHTNG_CLD2GRND', cld2grnd_flash_freq, pcols, lchnk )
    end do Chunk_loop

    do lchnk = begchunk,endchunk
       call outfld( 'FLASHFRQ',     flash_freq(:,lchnk),       pcols, lchnk )
       call outfld( 'CGIC',         cgic(:,lchnk),             pcols, lchnk )
       call outfld( 'CLDHGT',       cldhgt(:,lchnk),           pcols, lchnk )
       call outfld( 'DCHGZONE',     dchgzone(:,lchnk),         pcols, lchnk )
    enddo

    if (.not.calc_nox_prod) return

    !--------------------------------------------------------------------------------
    ! 	... Accumulate global total, convert to flashes per second
    ! 	... Accumulate global NO production rate
    !--------------------------------------------------------------------------------
    kk = pcols*(endchunk-begchunk+1)
    call shr_reprosum_calc( flash_freq, wrk2,kk,kk,1, commid=mpicom)
    glob_flashfreq=wrk2(1)/60._r8
    call shr_reprosum_calc( glob_prod_no_col, wrk2,kk,kk,1, commid=mpicom)
    glob_noprod = wrk2(1)
    if( masterproc ) then
       write(iulog,*) ' '
       write(iulog,'(''Global flash freq (/s), lightning NOx (TgN/y) = '',2f10.4)') &
            glob_flashfreq, glob_noprod
    end if

    if( glob_noprod > 0._r8 ) then
       !--------------------------------------------------------------------------------
       !	... Distribute production up to cloud top [Pickering et al., 1998 (JGR)]
       !--------------------------------------------------------------------------------
       do lchnk = begchunk,endchunk
          call get_rlat_all_p(lchnk, pcols, rlats)
          ncol = state(lchnk)%ncol
          call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk), cldtop_ndx, cldtop )
          do icol = 1,ncol
             cldtind = nint( cldtop(icol) )
             if( prod_no_col(icol,lchnk) > 0._r8 ) then
                if( cldhgt(icol,lchnk) > 0._r8 ) then
                   if( abs( rlats(icol) ) > lat25 ) then
                      itype = 1                                                     ! midlatitude continental
                   else if( nint( cam_in(lchnk)%landfrac(icol) ) == land ) then
                      itype = 3                                                     ! tropical continental
                   else
                      itype = 2                                                     ! topical marine
                   end if
                   frac_sum = 0._r8
                   do k = cldtind,pver
                      zlow       = zint(icol,k+1,lchnk) * m2km                      ! lower interface height (km)
                      zlow_scal  = zlow * 16._r8/cldhgt(icol,lchnk)                 ! scale to 16 km convection height
                      zlow_ind   = max( 1,INT(zlow_scal)+1 )                        ! lowest vdist index to include in layer
                      zhigh      = zint(icol,k,lchnk) * m2km                        ! upper interface height (km)
                      zhigh_scal = zhigh * 16._r8/cldhgt(icol,lchnk)                ! height (km) scaled to 16km convection height
                      zhigh_ind  = max( 1,MIN( 16,INT(zhigh_scal)+1 ) )             ! highest vdist index to include in layer
                      do kk = zlow_ind,zhigh_ind
                         wrk  = kk
                         wrk1 = kk - 1
                         fraction = min( zhigh_scal,wrk ) &                         ! fraction of vdist in this model layer
                              - max( zlow_scal,wrk1 )
                         fraction = max( 0._r8, min( 1._r8,fraction ) )
                         frac_sum = frac_sum + fraction*vdist(kk,itype)
                         prod_no(icol,k,lchnk) = prod_no(icol,k,lchnk) &            ! sum the fraction of column NOx in layer k
                              + fraction*vdist(kk,itype)*.01_r8
                      end do
                      prod_no(icol,k,lchnk) = prod_no_col(icol,lchnk) * prod_no(icol,k,lchnk) & ! multiply fraction by column amount
                                            / (km2cm*(zhigh - zlow))                ! and convert to atom N cm^-3 s^-1
                   end do
                end if
             end if
          end do
       end do
    end if

    !--------------------------------------------------------------------------------
    !       ... output lightning no production to history file
    !--------------------------------------------------------------------------------
    do lchnk = begchunk,endchunk
       call outfld( 'LNO_PROD',     prod_no(:,:,lchnk),        pcols, lchnk )
       call outfld( 'LNO_COL_PROD', glob_prod_no_col(:,lchnk), pcols, lchnk )
       call outfld( 'FLASHENGY',    flash_energy(:,lchnk),     pcols, lchnk )
    enddo

  end subroutine lightning_no_prod

end module mo_lightning
