!================================================================================
! radiation unit driver
!================================================================================
module unit_driver

  use shr_kind_mod,     only: r8=>SHR_KIND_R8
  use cam_abortutils,   only: endrun
  use spmd_utils,       only: masterproc
  use ppgrid,           only: pcols, pver, pverp, begchunk, endchunk
  use cam_logfile,      only: iulog

  implicit none
  private
  save

  public :: unit_driver_run
  public :: unit_driver_init
  public :: unit_driver_reg

contains

!================================================================================
!================================================================================
  subroutine unit_driver_reg
    use radiation_data, only: rad_data_enable

    call rad_data_enable()

  end subroutine unit_driver_reg

!================================================================================
!================================================================================
  subroutine unit_driver_run( indata, phys_state, pbuf2d, cam_out, cam_in, recno)
    use physics_types,    only: physics_state
    use physics_types,    only: physics_ptend
    use ppgrid,           only: begchunk, endchunk
    use time_manager,     only: timemgr_set_date_time
    use drv_input_data,   only: drv_input_data_t
    use camsrfexch,       only: cam_out_t, cam_in_t
    use physics_buffer,   only: physics_buffer_desc, pbuf_get_chunk
    use radiation_data,   only: rad_data_read
    use solar_data,       only: solar_data_advance

    implicit none

    type(drv_input_data_t), intent(inout) :: indata
    type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer,             intent(in)    :: recno

    integer :: c                                 ! chunk index

    type(physics_buffer_desc), pointer :: pbuf(:)
    
    ! Solar irradiance
    call solar_data_advance()

    ! get data needed to drive radiation ...
    call rad_data_read( indata, phys_state, pbuf2d, cam_in, recno=recno )

!$OMP PARALLEL DO PRIVATE ( c, pbuf )
    do c=begchunk,endchunk

       pbuf => pbuf_get_chunk(pbuf2d, c)
       call unit_driver_exec ( phys_state(c), pbuf, cam_out(c), cam_in(c) )

    end do

  end subroutine unit_driver_run

!================================================================================
!================================================================================
  subroutine unit_driver_exec ( state, pbuf, cam_out, cam_in )
    use physics_types,    only: physics_state
    use physics_types,    only: physics_ptend, physics_ptend_dealloc, physics_update
    use camsrfexch,       only: cam_out_t, cam_in_t
    use physics_buffer,   only: physics_buffer_desc
    use radiation,        only: radiation_tend

    type(physics_state), intent(inout) :: state
    type(cam_out_t),     intent(inout) :: cam_out
    type(cam_in_t),      intent(inout) :: cam_in
    type(physics_buffer_desc), pointer :: pbuf(:)

  ! local vars
    type(physics_ptend) :: ptend
    real(r8) :: net_flx(pcols)

    call radiation_tend( &
       state, ptend, pbuf, cam_out, cam_in, net_flx)

    call physics_ptend_dealloc(ptend)

  end subroutine unit_driver_exec 

!================================================================================
!================================================================================
  subroutine unit_driver_init
    use radheat,        only: radheat_disable_waccm

    implicit none

    call radheat_disable_waccm()

  end subroutine unit_driver_init

end module unit_driver
