module ionosphere_interface

 ! Dummy interface -- actual ionosphere interface exist in src/ionosphere/modelname

  implicit none

  private

  public :: ionosphere_readnl
  public :: ionosphere_init
  public :: ionosphere_run1
  public :: ionosphere_run2
  public :: ionosphere_init_restart
  public :: ionosphere_write_restart
  public :: ionosphere_read_restart
  public :: ionosphere_final

contains

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_readnl( nlfile )

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  end subroutine ionosphere_readnl

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_init()
    
  end subroutine ionosphere_init

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_run1(pbuf2d)
    use physics_buffer, only: physics_buffer_desc

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

  end subroutine ionosphere_run1

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_run2( phys_state, dyn_in, pbuf2d )

    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc
    use phys_grid,      only: begchunk, endchunk
    use dyn_comp,       only: dyn_import_t

    ! args
    type(physics_state),    intent(in) :: phys_state(begchunk:endchunk)
    type(dyn_import_t),     intent(in) :: dyn_in  ! dynamics import 

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

  end subroutine ionosphere_run2

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_init_restart(File)
    use pio, only: file_desc_t

    type(File_desc_t),  intent(inout) :: File

  end subroutine ionosphere_init_restart

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_write_restart(File)
    use pio, only: file_desc_t

    type(File_desc_t), intent(inout) :: File

  end subroutine ionosphere_write_restart

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_read_restart(File)
    use pio, only: file_desc_t

    type(file_desc_t), intent(inout) :: File

  end subroutine ionosphere_read_restart

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  subroutine ionosphere_final

  end subroutine ionosphere_final

end module ionosphere_interface
