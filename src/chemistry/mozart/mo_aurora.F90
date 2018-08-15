!-----------------------------------------------------------------------
! Stub version 
! 
!-----------------------------------------------------------------------

module mo_aurora

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pver,pcols

  implicit none

  interface aurora
     module procedure aurora_prod
     module procedure aurora_hrate
  end interface aurora

contains


  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine aurora_register

  endsubroutine aurora_register

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine aurora_inti(pbuf2d)
    use physics_buffer, only: physics_buffer_desc
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  end subroutine aurora_inti

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine aurora_timestep_init

  end subroutine aurora_timestep_init

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine aurora_prod( tn, o2, o1, mbar, rlats, &
                          qo2p, qop, qn2p, qnp, pmid, &
                          lchnk, calday,  ncol, rlons, pbuf )

    use physics_buffer,only: physics_buffer_desc

    !-----------------------------------------------------------------------
    !   ... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in) ::  &
         ncol, &                           ! column count
         lchnk                             ! chunk index
    real(r8), intent(in) :: &
         calday                           ! calendar day of year
    real(r8), intent(in) :: &
         tn(pcols,pver), &                ! neutral gas temperature (K)
         o2(ncol,pver), &                 ! O2 concentration (kg/kg)
         o1(ncol,pver), &                 ! O concentration (kg/kg)
         mbar(ncol,pver)                  ! mean molecular weight (g/mole)
    real(r8), intent(in) :: &
         pmid(pcols,pver)                 ! midpoint pressure (Pa)
    real(r8), intent(in) :: &
         rlats(ncol), &                   ! column latitudes (radians)
         rlons(ncol)
    real(r8), intent(out) :: &
         qo2p(ncol,pver), &               ! o2+ production
         qop(ncol,pver), &                ! o+ production
         qn2p(ncol,pver), &               ! n2+ production
         qnp(ncol,pver)                   ! n+ production

    type(physics_buffer_desc),pointer :: pbuf(:)

  end subroutine aurora_prod

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine aurora_hrate( tn, mbar, rlats, &
                           aur_hrate, cpair, pmid, lchnk, calday, &
                           ncol, rlons, pbuf  )
    use physics_buffer,only: physics_buffer_desc
    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in) ::  &
         ncol, &                           ! column count
         lchnk                             ! chunk index
    real(r8), intent(in) :: &
         calday                           ! calendar day of year
    real(r8), intent(in) :: &
         tn(pcols,pver), &                ! neutral gas temperature (K)
         mbar(ncol,pver)                  ! mean molecular weight (g/mole)
    real(r8), intent(in) :: &
         cpair(ncol,pver)                 ! specific heat capacity (J/K/kg)
    real(r8), intent(in) :: &
         pmid(pcols,pver)                 ! midpoint pressure (Pa)
    real(r8), intent(in) :: &
         rlats(ncol), &                   ! column latitudes (radians)
         rlons(ncol)
    real(r8), intent(out) :: &
         aur_hrate(ncol,pver)             ! auroral heating rate
    type(physics_buffer_desc),pointer :: pbuf(:)

  end subroutine aurora_hrate

end module mo_aurora
