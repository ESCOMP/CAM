module beljaars_drag

  implicit none
  private      
  save

  public init_blj                             ! Initialization
  public compute_blj                          ! Full routine

  ! ------------ !
  ! Private data !
  ! ------------ !

  integer,  parameter :: r8 = selected_real_kind(12) ! 8 byte real

  real(r8), parameter :: horomin= 1._r8       ! Minimum value of subgrid orographic height for mountain stress [ m ]
  real(r8), parameter :: z0max  = 100._r8     ! Maximum value of z_0 for orography [ m ]
  real(r8), parameter :: dv2min = 0.01_r8     ! Minimum shear squared [ m2/s2 ]
  real(r8)            :: orocnst              ! Converts from standard deviation to height [ no unit ]
  real(r8)            :: z0fac                ! Factor determining z_0 from orographic standard deviation [ no unit ] 
  real(r8)            :: karman               ! von Karman constant
  real(r8)            :: gravit               ! Acceleration due to gravity
  real(r8)            :: rair                 ! Gas constant for dry air

contains

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine init_blj( kind, gravit_in, rair_in , errstring )

    integer,  intent(in) :: kind
    real(r8), intent(in) :: gravit_in, rair_in

    character(len=*), intent(out) :: errstring

    errstring = ' '

    if ( kind /= r8 ) then
       errstring = 'inconsistent KIND of reals passed to init_blj'
       return
    endif

   gravit   = gravit_in
   rair     = rair_in
    
  end subroutine init_blj

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine compute_blj( pcols    , pver    , ncol    ,                     &
                          u        , v       , t       , pmid    , delp    , &
                          zm       , sgh     , drag    , taux    , tauy    , & 
                          landfrac )

    !------------------------------------------------------------------------------ !
    ! Beljaars Sub-Grid Orographic (SGO) Form drag parameterization                 !  
    !                                                                               !
    ! Returns drag profile and integrated stress associated with subgrid mountains  !
    ! with horizontal length scales nominally below 3km.  Similar to TMS but        !
    ! drag is distributed in the vertical (Beljaars et al., 2003, QJRMS).           !
    !                                                                               !
    ! First cut follows TMS.     J. Bacmeister, March 2016                          !
    !------------------------------------------------------------------------------ !

    ! ---------------------- !
    ! Input-Output Arguments ! 
    ! ---------------------- !

    integer,  intent(in)  :: pcols                 ! Number of columns dimensioned
    integer,  intent(in)  :: pver                  ! Number of model layers
    integer,  intent(in)  :: ncol                  ! Number of columns actually used

    real(r8), intent(in)  :: u(pcols,pver)         ! Layer mid-point zonal wind [ m/s ]
    real(r8), intent(in)  :: v(pcols,pver)         ! Layer mid-point meridional wind [ m/s ]
    real(r8), intent(in)  :: t(pcols,pver)         ! Layer mid-point temperature [ K ]
    real(r8), intent(in)  :: pmid(pcols,pver)      ! Layer mid-point pressure [ Pa ]
    real(r8), intent(in)  :: delp(pcols,pver)      ! Layer thickness [ Pa ]
    real(r8), intent(in)  :: zm(pcols,pver)        ! Layer mid-point height [ m ]
    real(r8), intent(in)  :: sgh(pcols)            ! Standard deviation of orography [ m ]
    real(r8), intent(in)  :: landfrac(pcols)       ! Land fraction [ fraction ]
    
    real(r8), intent(out) :: drag(pcols,pver)      ! SGO drag profile [ kg/s/m2 ]
    real(r8), intent(out) :: taux(pcols)           ! Surface zonal      wind stress [ N/m2 ]
    real(r8), intent(out) :: tauy(pcols)           ! Surface meridional wind stress [ N/m2 ]

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer  :: i,k                                ! Loop indices
    integer  :: kb, kt                             ! Bottom and top of source region
   
    real(r8) :: vmag                               ! Velocity magnitude [ m /s ]

    real(r8) :: alpha,beta,Cmd,Ccorr,n1,n2,k1,kflt,k2,IH
    real(r8) :: a1(pcols),a2(pcols)

    alpha =  12._r8
    beta  =  1._r8
    n1    = -1.9_r8
    n2    = -2.8_r8

    Cmd   = 0.005_r8
    Ccorr = 0.6_r8 * 5._r8

    kflt  = 0.00035_r8  ! m-1
    k1    = 0.003_r8    ! m-1
    IH    = 0.00102_r8  ! m-1

    a1(1:ncol)    = (sgh(1:ncol)*sgh(1:ncol)) / ( IH* (kflt**n1) )
    a2(1:ncol)    = a1(1:ncol) * k1**(n1-n2)


    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
       
    do k = 1, pver
    do i = 1, ncol
       Vmag      = SQRT( u(i,k)**2 + v(i,k)**2)
       drag(i,k) = -alpha * beta * Cmd * Ccorr * Vmag * 2.109_r8 *  & 
                   EXP ( -(zm(i,k)/1500._r8 )*SQRT(zm(i,k)/1500._r8) ) * ( zm(i,k)**(-1.2_r8) )  &
                   * a2(i)
    end do
    end do
    

    !---------------------------------!
    ! Diagnose effective surface drag !
    ! in X and Y by integrating in    !
    ! the vertical                    !
    !---------------------------------!
    ! FIXME: uses 'state' u and v. 
    ! Should updated u and v's be used?

    taux=0._r8
    tauy=0._r8
    do k = 1, pver
    do i = 1, ncol
       taux(i)  = taux(i) + drag(i,k)*u(i,k)*delp(i,k)/gravit
       tauy(i)  = tauy(i) + drag(i,k)*v(i,k)*delp(i,k)/gravit
    end do
    end do

    return
  end subroutine compute_blj

end module beljaars_drag
