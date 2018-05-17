!================================================================================
! aurora unit test driver
!================================================================================
module unit_driver

  use shr_kind_mod, only: r8=>SHR_KIND_R8
  use ppgrid,       only: pcols, pver, begchunk, endchunk

  implicit none
  private
  save

  public :: unit_driver_init
  public :: unit_driver_run

contains

  !==============================================================================
  !==============================================================================
  subroutine unit_driver_init
    use cam_history,  only : addfld, horiz_only

    call addfld( 'Tempr', (/ 'lev' /), 'I', '/s ', ' ')
    call addfld( 'O2vmr', (/ 'lev' /), 'I', '/s ', ' ')
    call addfld( 'O1vmr', (/ 'lev' /), 'I', '/s ', ' ')
    call addfld( 'qo2p', (/ 'lev' /), 'I', '/s ', ' ')
    call addfld( 'qop',  (/ 'lev' /), 'I', '/s ', ' ')
    call addfld( 'qn2p', (/ 'lev' /), 'I', '/s ', ' ')
    call addfld( 'qnp',  (/ 'lev' /), 'I', '/s ', ' ')
    
    call addfld('BNORTH',horiz_only,  'I','GAUSS', 'Northward component of magnetic field')
    call addfld('BEAST', horiz_only,  'I','GAUSS', 'Eastward component of magnetic field')
    call addfld('BDOWN', horiz_only,  'I','GAUSS', 'Downward component of magnetic field')
    call addfld('BMAG',  horiz_only,  'I','GAUSS', 'Magnetic field magnitude')

    call addfld('D1VEC', horiz_only,  'I', ' ', ' ')
    call addfld('D2VEC', horiz_only,  'I', ' ', ' ')

  end subroutine unit_driver_init

  !==============================================================================
  !==============================================================================
  subroutine unit_driver_run(indata, phys_state, pbuf2d, cam_out, cam_in, recno)
    use physics_types,  only: physics_state
    use ppgrid,         only: begchunk, endchunk
    use camsrfexch,     only: cam_out_t, cam_in_t     
    use drv_input_data, only: drv_input_data_t
    use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk
    use shr_const_mod,  only: shr_const_cpdair 
    use time_manager,   only: get_curr_calday
    use physconst,      only: mwdry

    use mo_aurora,      only: aurora_timestep_init, aurora
    use phys_grid,      only: get_rlat_all_p, get_rlon_all_p
    use ref_pres,       only: pref_mid
    use drv_input_data, only: drv_input_data_read
    use cam_history,    only: outfld
    use mo_apex,        only: bnorth, beast, bdown, bmag
    use mo_apex,        only: d1vec, d2vec

    type(drv_input_data_t), intent(inout) :: indata
    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer,             intent(in)    :: recno

    real(r8) :: mbar(pcols,pver,begchunk:endchunk) 
    real(r8) :: tfld(pcols,pver,begchunk:endchunk)
    real(r8) :: o1vmr(pcols,pver,begchunk:endchunk)
    real(r8) :: o2vmr(pcols,pver,begchunk:endchunk)
    real(r8) :: o1mmr(pcols,pver)
    real(r8) :: o2mmr(pcols,pver)
    real(r8) :: press(pcols,pver)

    real(r8) :: rlats(pcols)
    real(r8) :: rlons(pcols)

    real(r8) :: qo2p(pcols,pver)
    real(r8) :: qop(pcols,pver)
    real(r8) :: qn2p(pcols,pver)
    real(r8) :: qnp(pcols,pver)


    real(r8) :: calday

    integer :: c, ncol, k
    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8),parameter :: mwo1 = 16._r8
    real(r8),parameter :: mwo2 = 32._r8

    real(r8) :: aur_hrate(pcols,pver)                          ! chunk auroral heating rate
    real(r8) :: cpair(pcols,pver)                              ! specific heat capacity (J/K/kg)

    tfld = drv_input_data_read(indata, 'T', 'lev', pver, recno, abort=.true.)
    o1vmr = drv_input_data_read(indata, 'O', 'lev', pver, recno, abort=.true.)
    o2vmr = drv_input_data_read(indata, 'O2', 'lev', pver, recno, abort=.true.)

    mbar = mwdry  

    call  aurora_timestep_init

    do k=1,pver
      press(:,k) = pref_mid(k)
    enddo

    !-----------------------------------------------------------------------
    ! get current calendar day of year
    !-----------------------------------------------------------------------
    calday = get_curr_calday( )

    cpair(:,:) = shr_const_cpdair

!$OMP PARALLEL DO PRIVATE ( c, pbuf, ncol, rlats, rlons,qo2p,qop,qn2p,qnp  )
    do c=begchunk,endchunk

      ncol = phys_state(c)%ncol
      
      o1mmr(:ncol,:) = o1vmr(:ncol,:,c)*mwo1/mwdry
      o2mmr(:ncol,:) = o2vmr(:ncol,:,c)*mwo2/mwdry

      call get_rlat_all_p( c, ncol, rlats(:ncol) )
      call get_rlon_all_p( c, ncol, rlons(:ncol) )

      pbuf => pbuf_get_chunk(pbuf2d, c)

      call aurora( tfld(:pcols,:,c), o2mmr(:ncol,:), o1mmr(:ncol,:), mbar(:ncol,:,c), rlats(:ncol), &
                   qo2p(:ncol,:),qop(:ncol,:),qn2p(:ncol,:),qnp(:ncol,:), &
                   press(:,:), c, calday, ncol, rlons(:ncol), pbuf )

      call outfld( 'Tempr', tfld(:ncol,:,c), ncol, c )
      call outfld( 'O2vmr', o2vmr(:ncol,:,c), ncol, c )
      call outfld( 'O1vmr', o1vmr(:ncol,:,c), ncol, c )
      call outfld( 'qo2p', qo2p(:ncol,:), ncol, c )
      call outfld( 'qop',  qop (:ncol,:), ncol, c )
      call outfld( 'qn2p', qn2p(:ncol,:), ncol, c )
      call outfld( 'qnp',  qnp (:ncol,:), ncol, c )

      call outfld('BNORTH', bnorth(:ncol,c), ncol, c )
      call outfld('BEAST',  beast(:ncol,c),  ncol, c )
      call outfld('BDOWN',  bdown(:ncol,c),  ncol, c )
      call outfld('BMAG',   bmag(:ncol,c),   ncol, c )
      call outfld('D1VEC',  d1vec(1,:ncol,c),  ncol, c )
      call outfld('D2VEC',  d2vec(1,:ncol,c),  ncol, c )

      call aurora( tfld(:pcols,:,c), o2mmr(:ncol,:), o1mmr(:ncol,:), mbar(:ncol,:,c), rlats(:ncol), &
                   aur_hrate(:ncol,:), cpair(:ncol,:), &
                   press(:,:), c, calday, ncol, rlons(:ncol) )

      call outfld( 'QRS_AUR', aur_hrate(:ncol,:), ncol, c )

    end do

  end subroutine unit_driver_run

end module unit_driver
