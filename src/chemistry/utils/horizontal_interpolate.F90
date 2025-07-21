module horizontal_interpolate

!   
! Modules Used:
!
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use shr_const_mod,      only: SHR_CONST_PI
  use cam_abortutils,     only: endrun
  use scamMod,            only: single_column
  use cam_logfile,        only: iulog
  implicit none
  private
  save

  public :: xy_interp_init, xy_interp

contains

  subroutine check_invariant(invariant_condition, message)
    logical, intent(in)  :: invariant_condition
    character(len=*), intent(in) :: message
    if ( .not.(invariant_condition) ) then
       call endrun("Invariant check failed: " // message)
    endif
  end subroutine check_invariant


  real(r8) function normalize_lon_right(left, right)
    ! Normalize the right side of a longitude interval such that
    ! norm_right > left and (norm_right - left) is in (0, 360]
    real(r8), intent(in)  :: left, right

    normalize_lon_right = right

    do while (normalize_lon_right <= left)
       normalize_lon_right = normalize_lon_right + 360.0_r8
    end do
    do while (normalize_lon_right - 360.0_r8 > left)
       normalize_lon_right = normalize_lon_right - 360.0_r8
    end do
    return
  end function normalize_lon_right


  real(r8) function lon_length(left, right)
    ! Compute a longitude interval length, accouting for wrapping
    ! around the globe.  Input values are normalized so that the
    ! return value is always in the range (0, 360].
    real(r8), intent(in)  :: left, right

    lon_length = normalize_lon_right(left, right) - left
    return
  end function lon_length


  real(r8) function normalize_lon_value(lon)
    ! Normalize a lon value to be in [0, 360)
    real(r8), intent(in)  :: lon
    real(r8) :: norm_lon

    norm_lon = lon
    do while (norm_lon < 0)
       norm_lon = norm_lon + 360.0_r8
    end do
    do while (norm_lon >= 360.0_r8)
       norm_lon = norm_lon - 360.0_r8
    end do

    normalize_lon_value = norm_lon
    return
  end function normalize_lon_value


  real(r8) function calculate_lon_overlap(input_left, input_right, sim_left, sim_right)
    ! Return the length of the overlap between the input and
    ! simulation longitude ranges.  Values are normalized before
    ! calculation to ensure correct handling of ranges which wrap over
    ! zero.
    real(r8), intent(in)  :: input_left, input_right, sim_left, sim_right
    real(r8) norm_input_left, norm_input_right, norm_sim_left, norm_sim_right
    real(r8) overlap_left, overlap_right

    ! Normalzie so norm_sim_left is in [0, 360) and norm_sim_left < norm_sim_right
    norm_sim_left = normalize_lon_value(sim_left)
    norm_sim_right = normalize_lon_right(norm_sim_left, sim_right)

    ! Normalize the input values to ensure that norm_input_left < norm_sim_left
    norm_input_left = normalize_lon_value(input_left) - 360.0_r8  ! now in [-360, 0)
    norm_input_right = normalize_lon_right(norm_input_left, input_right)

    ! if norm_input is strictly to the left of norm_sim, slide up by 360
    do while (norm_input_right <= norm_sim_left)
       norm_input_left = norm_input_left + 360.0_r8
       norm_input_right = norm_input_right + 360.0_r8
    end do

    ! Compute overlap
    overlap_left = merge(norm_input_left, norm_sim_left, norm_input_left > norm_sim_left)
    overlap_right = merge(norm_input_right, norm_sim_right, norm_input_right < norm_sim_right)
    if (overlap_left < overlap_right) then
       calculate_lon_overlap = overlap_right - overlap_left
    else
       calculate_lon_overlap = 0
    endif
    return
  end function calculate_lon_overlap


  real(r8) function lon_weight(input_left, input_right, sim_left, sim_right, use_flight_distance)
    ! Compute how much input data in the (input_left, input_right)
    ! longitude band contributes to simulation data in the (sim_left,
    ! sim_right) band.
    !
    ! use_flight_distance indicates how the input data should be
    ! interpreted.  If true, input data is assumed to be a total value
    ! for the input grid cell.  If false, input data is assumed to be
    ! a mixing ratio.
    real(r8), intent(in)  :: input_left, input_right, sim_left, sim_right
    logical,  intent(in)  :: use_flight_distance   !.true. = flight distance, .false. =  all mixing ratios
    real(r8)  :: overlap_len

    ! Sanity check that inputs aren't too large or too small.  For
    ! really huge floating point values, adding / subtracting 360 will
    ! encounter roundoff errors, and could cause while loops to go
    ! forever.  Values much outside the [-360, 360] range are probably
    ! an error anyway, so abort if we encounter anything suspicious here.
    call check_invariant(((-2000.0 < input_left).and.(input_left < 2000.0)), "input_left is out of range")
    call check_invariant(((-2000.0 < input_right).and.(input_right < 2000.0)), "input_right is out of range")
    call check_invariant(((-2000.0 < sim_left).and.(sim_left < 2000.0)), "sim_left is out of range")
    call check_invariant(((-2000.0 < sim_right).and.(sim_right < 2000.0)), "sim_right is out of range")

    overlap_len = calculate_lon_overlap(input_left, input_right, sim_left, sim_right)

    if (overlap_len == 0) then
       ! No overlap; weight is zero
       lon_weight = 0

    elseif (use_flight_distance) then
       ! Data values are total within the grid cell.  Hence, the
       ! weight is just the fraction of the area of the original grid
       ! cell which overlaps the new cell.
       lon_weight = overlap_len / lon_length(input_left, input_right)
    else
       ! Data values are mixing ratios.  To compute how much this
       ! input grid cell contributes to the mixing ratio in the
       ! simulation grid cell, we multiply by the fraction of the sim
       ! grid cell occupied by the overlap.  Fortunately this is just
       ! the ratio of the longitude ranges, since the grid cells are
       ! rectangular in lat,lon coordinates.
       lon_weight = overlap_len / lon_length(sim_left, sim_right)
    endif

    ! Floating point precision can wind up with lon_weight slightly
    ! out of the [0, 1] range, so fix up the values if that happens.
    ! If things are farther off, then crash with an error.
    call check_invariant((-0.0000001_r8 <= lon_weight) .and. (lon_weight <= 1.0000001_r8), "lon_weight must be in [0, 1]")
    if (lon_weight < 0.0_r8) then
       lon_weight = 0.0_r8
    endif
    if (lon_weight > 1.0_r8) then
       lon_weight = 1.0_r8
    endif

    return
  end function lon_weight


  real(r8) function lat_band_weight(low, high)
    ! Retuns the unit-less weight of the latitude band on the globe
    ! from [low, high].  The area of a latitudinal slice of height h
    ! is just 2*pi*r*h.  Since we only need to get ratios of areas,
    ! this just computes the difference in the y-coordinates from low
    ! to high.
    real(r8), intent(in)  :: low, high

    lat_band_weight = sin(high * SHR_CONST_PI / 180.0_r8) - sin(low * SHR_CONST_PI / 180.0_r8)
    return
  end function lat_band_weight


  real(r8) function normalize_lat(lat)
    ! Truncate a latitude value to be within [-90, 90].  Values much
    ! outside that range are probably an error in the calling code, so
    ! exit if we see those.
    real(r8), intent(in)  :: lat
    call check_invariant(((-91.0_r8 < lat).and.(lat < 91.0_r8)), "Lat value is out of expected range.")
    if (lat < -90.0_r8) then
       normalize_lat = -90.0_r8
    else if (lat > 90.0_r8) then
       normalize_lat = 90.0_r8
    else
       normalize_lat = lat
    endif
  end function normalize_lat


  real(r8) function lat_weight(input_bot, input_top, sim_bot, sim_top, use_flight_distance)
    ! Compute how much input data in the (input_left, input_right)
    ! longitude band contributes to simulation data in the (sim_left,
    ! sim_right) band.
    !
    ! use_flight_distance indicates how the input data should be
    ! interpreted.  If true, input data is assumed to be a total value
    ! for the input grid cell.  If false, input data is assumed to be
    ! a mixing ratio.
    real(r8), intent(in)  :: input_bot, input_top, sim_bot, sim_top
    logical,  intent(in)  :: use_flight_distance   !.true. = flight distance, .false. =  all mixing ratios
    real(r8)  :: norm_input_bot, norm_input_top, norm_sim_bot, norm_sim_top
    real(r8)  :: overlap_bot, overlap_top

    ! Make sure that inputs are in [-90, 90] and that bot < top
    norm_input_bot = normalize_lat(input_bot)
    norm_input_top = normalize_lat(input_top)
    norm_sim_bot = normalize_lat(sim_bot)
    norm_sim_top = normalize_lat(sim_top)
    call check_invariant(norm_input_bot < norm_input_top, "must have input_bot < input_top")
    call check_invariant(norm_sim_bot < norm_sim_top, "must have sim_bot < sim_top")

    overlap_bot = merge(norm_input_bot, norm_sim_bot, norm_input_bot > norm_sim_bot)
    overlap_top = merge(norm_input_top, norm_sim_top, norm_input_top < norm_sim_top)

    if ( (norm_input_bot == norm_input_top) .or. &
         (norm_sim_bot == norm_sim_top) .or. &
         (overlap_top <= overlap_bot) ) then
       ! No overlap
       lat_weight = 0
    elseif (use_flight_distance) then
       ! Input values are a total for the grid cell, so the weight is
       ! just the fraction of the input cell which overlaps the
       ! sim cell.
       lat_weight = lat_band_weight(overlap_bot, overlap_top) / lat_band_weight(norm_input_bot, norm_input_top)
       call check_invariant((0 <= lat_weight) .and. (lat_weight <= 1), "dist: lat_weight must be in [0, 1]")
    else
       ! Input values are a mixing ratio.  The amount added in the
       ! overlap is spread evenly over the sim cell.
       lat_weight = lat_band_weight(overlap_bot, overlap_top) / lat_band_weight(norm_sim_bot, norm_sim_top)
       call check_invariant((0 <= lat_weight) .and. (lat_weight <= 1), "mixrat: lat_weight must be in [0, 1]")
    endif
    return
  end function lat_weight


  subroutine xy_interp_init(num_input_lons, num_input_lats, input_lon_radians, input_lat_radians, &
       num_sim_lons, num_sim_lats, weight_x, weight_y, use_flight_distance)
!------------------------------------------------------------------------------------------------------------
! This program computes weighting functions to map a variable of (num_input_lons,num_input_lats) resolution to (num_sim_lons,num_sim_lats) resolution
! weight_x(num_sim_lons,num_input_lons) is the weighting function for zonal interpolation
! weight_y(num_sim_lats,num_input_lats) is the weighting function for meridional interpolation
!
! Author: Chih-Chieh (Jack) Chen  -- May 2010
!         Rob von Behren (jrvb@google.com)  -- Oct 2024
!
!------------------------------------------------------------------------------------------------------------
  integer,  intent(in)  :: num_input_lons, num_input_lats, num_sim_lons, num_sim_lats
  logical,  intent(in)  :: use_flight_distance   !.true. = flight distance, .false. =  all mixing ratios
  real(r8), intent(in)  :: input_lon_radians(num_input_lons), input_lat_radians(num_input_lats)
  real(r8), intent(out) :: weight_x(num_sim_lons,num_input_lons), weight_y(num_sim_lats,num_input_lats)

  real(r8) :: input_lon(num_input_lons), input_lat(num_input_lats)
  real(r8) :: sim_lon(num_sim_lons), sim_lat(num_sim_lats)
  real(r8) :: input_lon_edge(num_input_lons+1), sim_lon_edge(num_sim_lons+1)
  real(r8) :: input_lat_edge(num_input_lats+1), sim_lat_edge(num_sim_lats+1)
  real(r8) :: x1_west, x1_east, x2_west, x2_east
  real(r8) :: y1_south, y1_north, y2_south, y2_north
  integer  :: i1, j1, i2, j2, i

  weight_x(:,:) = 0.0_r8
  weight_y(:,:) = 0.0_r8

! input_lon_radians & input_lat_radians are longitude & latitude on the source mesh in radians
! convert input_lon, input_lat from radians to degrees
  input_lon(:) = input_lon_radians(:)/SHR_CONST_PI*180.0_r8
  input_lat(:) = input_lat_radians(:)/SHR_CONST_PI*180.0_r8

  ! Set up sim_lon, sim_lat (target mesh), in CAM convention.  The
  ! (lon,lat) pairs are the center points of the grid cells.
  do i2=1,num_sim_lons
     sim_lon(i2) = (float(i2)-1.0_r8)*360.0_r8/float(num_sim_lons)
  enddo
  do j2=1,num_sim_lats
     sim_lat(j2) = -90.0_r8+(float(j2)-1.0_r8)*180.0_r8/(float(num_sim_lats)-1.0_r8)
  enddo
  ! make sure the highest value is exactly +90, since the
  ! multiplication above could give a value that is slightly off due
  ! to rounding.
  sim_lat(num_sim_lats) = 90.0_r8

  ! Calculate the grid cell edges from the midpoints above.  We set
  ! things up so it is easy to find the boundary for grid cell i, j as:
  !    lat bounds:   (lat_edge(i), lat_edge(i+1))
  !    lon boudns:   (lon_edge(j), lon_edge(j+1))
  ! For latitudes, we ensure that the smallest and largest edges are
  ! at -90 and 90, respectively.  For longitudes, the values wrap so
  ! lon_edge(1) == lon_edge(num+1) - 360

  ! Input longitude edges
  do i1=2, num_input_lons
     input_lon_edge(i1) = (input_lon(i1-1) + input_lon(i1)) / 2.0_r8
  enddo
  input_lon_edge(1) = (input_lon(num_input_lons) - 360_r8 + input_lon(1)) / 2.0_r8
  input_lon_edge(num_input_lons+1) = input_lon_edge(1) + 360_r8

  do i2=2,num_sim_lons
     sim_lon_edge(i2) = (sim_lon(i2-1) + sim_lon(i2)) / 2.0_r8
  enddo
  sim_lon_edge(1) = (sim_lon(num_sim_lons) - 360_r8 + sim_lon(1)) / 2.0_r8
  sim_lon_edge(num_sim_lons+1) = sim_lon_edge(1) + 360_r8

  ! set up staggered lattiudes (cell edges in y)
  input_lat_edge(1)=-90.0_r8
  do j1=2,num_input_lats
	input_lat_edge(j1) = (input_lat(j1-1)+input_lat(j1))/2.0_r8
  enddo
  input_lat_edge(num_input_lats+1)=90.0_r8

  sim_lat_edge(1) = -90.0_r8
  do j2=2,num_sim_lats
     sim_lat_edge(j2) = (sim_lat(j2-1) + sim_lat(j2)) / 2.0_r8
  enddo
  sim_lat_edge(num_sim_lats+1) = 90.0_r8

  ! Compute the weight for all (input_lon, sim_lon) pairs
  do i2=1, num_sim_lons
     do i1=1, num_input_lons
        weight_x(i2, i1) = lon_weight(input_lon_edge(i1), input_lon_edge(i1+1), sim_lon_edge(i2), sim_lon_edge(i2+1), use_flight_distance)
     enddo
  enddo

  ! Compute the weight for all (input_lat, sim_lat) pairs
  do j2=1, num_sim_lats
     do j1=1, num_input_lats
        weight_y(j2, j1) = lat_weight(input_lat_edge(j1), input_lat_edge(j1+1), sim_lat_edge(j2), sim_lat_edge(j2+1), use_flight_distance)
     enddo
  enddo

  end subroutine xy_interp_init


  subroutine xy_interp(im1,jm1,km1,im2,jm2,pcols,ncols,weight_x,weight_y,var_src,var_trg,lons,lats,count_x,count_y,index_x,index_y)
!-------------------------------------------------------------------------------------------------------------
! This program interpolates var_src(im1,jm1,km1) to var_trg(im2,jm2,km1) based on weighting functions weight_x & weight_y.
!-------------------------------------------------------------------------------------------------------------
  implicit none
  integer,  intent(in)  :: im1   ! source number of longitudes
  integer,  intent(in)  :: jm1   ! source number of latitudes
  integer,  intent(in)  :: km1   ! source/target number of levels
  integer,  intent(in)  :: im2   ! target number of longitudes
  integer,  intent(in)  :: jm2   ! target number of latitudes
  integer,  intent(in)  :: pcols
  integer,  intent(in)  :: ncols
  real(r8), intent(in)  :: weight_x(im2,im1), weight_y(jm2,jm1)
  real(r8), intent(in)  :: var_src(im1,jm1,km1)
  integer,  intent(in)  :: lons(pcols), lats(pcols)
  integer,  intent(in)  :: count_x(im2), count_y(jm2)
  integer,  intent(in)  :: index_x(im2,im1), index_y(jm2,jm1)
  real(r8), intent(out) :: var_trg(pcols,km1)
  integer  :: n, i1, j1, k1, i2, j2, ii, jj
  real(r8) :: sum_x

  var_trg(:,:) = 0.0_r8
 

  do k1=1,km1
     do n=1,ncols
! interpolate in x
        do jj=1,count_y(lats(n))
           sum_x = 0.0_r8
           do ii=1,count_x(lons(n))
              sum_x = sum_x + var_src(index_x(lons(n),ii),index_y(lats(n),jj),k1)* &
                             weight_x(lons(n),index_x(lons(n),ii))
           enddo
           var_trg(n,k1) = var_trg(n,k1)+sum_x*weight_y(lats(n),index_y(lats(n),jj))
        enddo
     enddo
  enddo

  end subroutine xy_interp



end module horizontal_interpolate
