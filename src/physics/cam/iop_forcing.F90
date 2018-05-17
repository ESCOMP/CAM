
  module iop_forcing

  ! ----------------------------------------------- !
  ! User-defined manipulation of forcings for SCAM. !
  !                                                 ! 
  ! Author : Sungsu Park and John Truesdale         !
  !                                                 !
  ! ----------------------------------------------- !

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  private  

  public scam_use_iop_srf
  public scam_set_iop_Tg

  contains

  subroutine scam_use_iop_srf( cam_in )
  ! ------------------------------------------------------- !
  ! Use the SCAM-IOP specified surface LHFLX/SHFLX/ustar/Tg !
  ! instead of using internally-computed values.            !
  !                                                         !
  ! Author : John Truesdale and Sungsu Park                 !
  ! ------------------------------------------------------- !    
    use ppgrid,           only: begchunk, endchunk
    use camsrfexch,       only: cam_in_t    
    use physconst,        only: stebol, latvap
    use scamMod
    use cam_abortutils,   only: endrun

    implicit none
    save

    type(cam_in_t), intent(INOUT) :: cam_in(begchunk:endchunk)
    integer                       :: c    ! Chunk index
    integer                       :: ncol ! Number of columns

    if( scm_iop_lhflxshflxTg .and. scm_iop_Tg ) then
        call endrun( 'scam_use_iop_srf : scm_iop_lhflxshflxTg and scm_iop_Tg must not be specified at the same time.')
    end if

    if( scm_iop_lhflxshflxTg ) then
        do c = begchunk, endchunk
           ncol = cam_in(c)%ncol
           if( have_lhflx ) then
               cam_in(c)%lhf(1)    = lhflxobs(1)
               cam_in(c)%cflx(1,1) = lhflxobs(1) / latvap
           endif
           if( have_shflx ) cam_in(c)%shf(1) = shflxobs(1)
           if( have_Tg ) then
               cam_in(c)%ts(1)   = tground(1)
               cam_in(c)%lwup(1) = stebol * tground(1)**4
           endif
        end do
    endif

    if( scm_iop_Tg .or. scm_crm_mode) then
        do c = begchunk, endchunk
           ncol = cam_in(c)%ncol
           if( have_Tg ) then
               cam_in(c)%ts(1) = tground(1)
               cam_in(c)%lwup(1) = stebol * tground(1)**4
           endif
        end do
    endif

  end subroutine scam_use_iop_srf

  subroutine scam_set_iop_Tg( cam_out )
  ! ----------------------------- !
  ! USE the SCAM-IOP specified Tg !
  ! ----------------------------- !
    use ppgrid,           only: begchunk, endchunk
    use camsrfexch,       only: cam_out_t    
    use scamMod
    use cam_abortutils,   only: endrun

    implicit none
    save

    type(cam_out_t), intent(INOUT) :: cam_out(begchunk:endchunk)
    integer                        :: c    ! Chunk index

    if( scm_iop_lhflxshflxTg .and. scm_iop_Tg ) then
        call endrun( 'scam_use_iop_srf : scm_iop_lhflxshflxTg and scm_iop_Tg must not be specified at the same time.')
    end if
    if( scm_iop_Tg .or. scm_crm_mode ) then
        do c = begchunk, endchunk
           cam_out(c)%tbot(1) = tground(1)
        end do
    endif

  end subroutine scam_set_iop_Tg

  end module iop_forcing
