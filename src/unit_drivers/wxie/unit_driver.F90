!================================================================================
! stub unit driver
!================================================================================
module unit_driver

  use shr_kind_mod, only: r8=>SHR_KIND_R8

  use short_lived_species, only: slvd_index, slvd_pbf_ndx => pbf_idx ! Routines to access short lived species

  implicit none
  private
  save

  public :: unit_driver_run
  public :: unit_driver_init
  public :: unit_driver_reg

  integer :: index_ped,index_hall,index_te,index_ti,index_ui,index_vi,index_wi
  integer :: index_op, index_nop, index_n2p, index_o2p
  integer :: index_o2, index_o, index_h

contains

  !================================================================================
  !================================================================================
  subroutine unit_driver_reg
  end subroutine unit_driver_reg

  !==============================================================================
  !==============================================================================
  subroutine unit_driver_run(indata, phys_state, pbuf2d, cam_out, cam_in, recno)
    use physics_types,    only: physics_state
    use ppgrid,           only: begchunk, endchunk, pver, pverp, pcols
    use camsrfexch,       only: cam_out_t, cam_in_t     
    use physics_buffer,   only: physics_buffer_desc
    use drv_input_data,   only: drv_input_data_t, drv_input_2d_t, drv_input_3d_t, drv_input_data_get
    use physics_buffer,   only: pbuf_get_field, pbuf_get_chunk

    use ionosphere_interface,only: ionosphere_run1, ionosphere_run2
    use mo_apex, only : mo_apex_init
    
    type(drv_input_data_t), intent(inout) :: indata
    type(physics_state), target, intent(inout) :: phys_state(begchunk:endchunk)
    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer,             intent(in)    :: recno

    type(drv_input_3d_t) :: ped_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: hall_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: te_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: ti_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: ui_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: vi_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: wi_ptrs(begchunk:endchunk)

    type(drv_input_3d_t) :: op_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: nop_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: n2p_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: o2p_ptrs(begchunk:endchunk)
    
    type(drv_input_3d_t) :: o2_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: o_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: h_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: t_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: u_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: v_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: zm_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: zi_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: omega_ptrs(begchunk:endchunk)
    type(drv_input_3d_t) :: pmid_ptrs(begchunk:endchunk)

    type(drv_input_2d_t) :: phis_ptrs(begchunk:endchunk)

    type(physics_buffer_desc), pointer :: pbuf(:)
    integer :: c

    call mo_apex_init(phys_state)
    
    do c=begchunk,endchunk
       pbuf => pbuf_get_chunk(pbuf2d, c)
       
       call pbuf_get_field(pbuf, index_hall, hall_ptrs(c)%array )
       call pbuf_get_field(pbuf, index_ped,  ped_ptrs(c)%array )
       call pbuf_get_field(pbuf, index_te,   te_ptrs(c)%array )
       call pbuf_get_field(pbuf, index_ti,   ti_ptrs(c)%array )
       call pbuf_get_field(pbuf, index_ui,   ui_ptrs(c)%array )
       call pbuf_get_field(pbuf, index_vi,   vi_ptrs(c)%array )
       call pbuf_get_field(pbuf, index_wi,   wi_ptrs(c)%array )
       call pbuf_get_field(pbuf, slvd_pbf_ndx,  op_ptrs(c)%array, start=(/1,1,index_op /), kount=(/pcols,pver,1/) )
       call pbuf_get_field(pbuf, slvd_pbf_ndx, nop_ptrs(c)%array, start=(/1,1,index_nop/), kount=(/pcols,pver,1/) )
       call pbuf_get_field(pbuf, slvd_pbf_ndx, n2p_ptrs(c)%array, start=(/1,1,index_n2p/), kount=(/pcols,pver,1/) )
       call pbuf_get_field(pbuf, slvd_pbf_ndx, o2p_ptrs(c)%array, start=(/1,1,index_o2p/), kount=(/pcols,pver,1/) )

       o2_ptrs(c)%array => phys_state(c)%q(:,:,index_o2)
       o_ptrs(c)%array => phys_state(c)%q(:,:,index_o)
       h_ptrs(c)%array => phys_state(c)%q(:,:,index_h)
       
       phis_ptrs(c)%array => phys_state(c)%phis
       t_ptrs(c)%array => phys_state(c)%t
       u_ptrs(c)%array => phys_state(c)%u
       v_ptrs(c)%array => phys_state(c)%v
       zm_ptrs(c)%array => phys_state(c)%zm
       zi_ptrs(c)%array => phys_state(c)%zi
       omega_ptrs(c)%array => phys_state(c)%omega
       pmid_ptrs(c)%array => phys_state(c)%pmid

    end do

    call drv_input_data_get( indata, 'DPIE_OPmmr', 'lev', pver, recno, op_ptrs)
    call drv_input_data_get( indata, 'DPIE_NOPmmr', 'lev', pver, recno, nop_ptrs)
    call drv_input_data_get( indata, 'DPIE_N2Pmmr', 'lev', pver, recno, n2p_ptrs)
    call drv_input_data_get( indata, 'DPIE_O2Pmmr', 'lev', pver, recno, o2p_ptrs)
    
    call drv_input_data_get( indata, 'DPIE_O', 'lev', pver, recno, o_ptrs )
    call drv_input_data_get( indata, 'DPIE_O2', 'lev', pver, recno, o2_ptrs )
    call drv_input_data_get( indata, 'DPIE_Hmmr', 'lev', pver, recno, h_ptrs )
    call drv_input_data_get( indata, 'T', 'lev', pver, recno, t_ptrs )
    call drv_input_data_get( indata, 'U', 'lev', pver, recno, u_ptrs )
    call drv_input_data_get( indata, 'V', 'lev', pver, recno, v_ptrs )
    call drv_input_data_get( indata, 'Zm', 'lev', pver, recno, zm_ptrs )
    call drv_input_data_get( indata, 'Zi', 'ilev', pverp, recno, zi_ptrs )
    call drv_input_data_get( indata, 'OMEGA', 'lev', pver, recno, omega_ptrs )
    call drv_input_data_get( indata, 'PMID', 'lev', pver, recno, pmid_ptrs )
    call drv_input_data_get( indata, 'SIGMAHAL', 'lev', pver, recno, hall_ptrs )
    call drv_input_data_get( indata, 'SIGMAPED', 'lev', pver, recno, ped_ptrs )
    call drv_input_data_get( indata, 'TElec', 'lev', pver, recno, te_ptrs )
    call drv_input_data_get( indata, 'TIon',  'lev', pver, recno, ti_ptrs )
    call drv_input_data_get( indata, 'UI',  'lev', pver, recno, ui_ptrs )
    call drv_input_data_get( indata, 'VI',  'lev', pver, recno, vi_ptrs )
    call drv_input_data_get( indata, 'WI',  'lev', pver, recno, wi_ptrs )

    call drv_input_data_get( indata, 'PHIS',      recno, phis_ptrs )

    call ionosphere_run1(pbuf2d)

    call ionosphere_run2( phys_state, pbuf2d )

  end subroutine unit_driver_run

  !==============================================================================
  !==============================================================================
  subroutine unit_driver_init
    use physics_buffer, only: pbuf_get_index
    use constituents,   only: cnst_get_ind

    index_ped  = pbuf_get_index('PedConduct')
    index_hall = pbuf_get_index('HallConduct')

    index_te   = pbuf_get_index('TElec')
    index_ti   = pbuf_get_index('TIon')
    index_ui   = pbuf_get_index('UI')
    index_vi   = pbuf_get_index('VI')
    index_wi   = pbuf_get_index('WI')

    index_op  = slvd_index( 'Op' )
    index_nop = slvd_index( 'NOp' )
    index_n2p = slvd_index( 'N2p' )
    index_o2p = slvd_index( 'O2p' )

    call cnst_get_ind('O', index_o)
    call cnst_get_ind('O2', index_o2)
    call cnst_get_ind('H', index_h)
    
  end subroutine unit_driver_init

end module unit_driver
