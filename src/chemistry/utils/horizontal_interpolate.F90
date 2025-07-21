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

  ! TODO(jrvb): These globals are only used by old_xy_interp_init.  Remove for the final PR
  real(r8) :: gw1(1000), gw2(1000)

  public :: xy_interp_init, xy_interp, xy_interp_run_unit_tests

contains

  ! TODO(jrvb): remove for final PR
  logical function env_var_equals(varname, checkval)
    use shr_sys_mod,      only: shr_sys_getenv

    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: checkval

    character(len=80)  :: envval
    integer :: rcode             ! shr_sys_getenv return code

    envval = ' '
    call shr_sys_getenv (varname, envval, rcode)

    if (envval .eq. checkval) then
       env_var_equals =  .true.
    else
       env_var_equals =  .false.
    endif

  end function env_var_equals


  subroutine old_xy_interp_init(im1,jm1,lon0,lat0,im2,jm2,weight_x,weight_y,use_flight_distance)
!------------------------------------------------------------------------------------------------------------
! This program computes weighting functions to map a variable of (im1,jm1) resolution to (im2,jm2) resolution
! weight_x(im2,im1) is the weighting function for zonal interpolation
! weight_y(jm2,jm1) is the weighting function for meridional interpolation
! 
! Author: Chih-Chieh (Jack) Chen  -- May 2010
!
!------------------------------------------------------------------------------------------------------------
  implicit none
  integer,  intent(in)  :: im1, jm1, im2, jm2
  logical,  intent(in)  :: use_flight_distance   !.true. = flight distance, .false. =  all mixing ratios
  real(r8), intent(in)  :: lon0(im1), lat0(jm1)
  real(r8), intent(out) :: weight_x(im2,im1), weight_y(jm2,jm1)

  real(r8) :: lon1(im1), lat1(jm1)
  real(r8) :: lon2(im2), lat2(jm2)
  real(r8) :: slon1(im1+1), slon2(im2+1), slat1(jm1+1), slat2(jm2+1)  
  real(r8) :: x1_west, x1_east, x2_west, x2_east
  real(r8) :: y1_south, y1_north, y2_south, y2_north
  integer  :: i1, j1, i2, j2

  ! NOTE(jrvb): the code below assumes input lat,lon values are sorted
  ! and bounds are within a typical range, but these invarianbts are
  ! not checked.

  weight_x(:,:) = 0.0_r8
  weight_y(:,:) = 0.0_r8

! lon0 & lat0 are longitude & latitude on the source mesh in radians
! convert lon1, lat1 from radians to degrees
  lon1(:) = lon0(:)/SHR_CONST_PI*180.0_r8
  lat1(:) = lat0(:)/SHR_CONST_PI*180.0_r8

! set up lon2, lat2 (target mesh), in CAM convention
  do i2=1,im2
   lon2(i2) = (float(i2)-1.0_r8)*360.0_r8/float(im2)
  enddo
  do j2=1,jm2
   lat2(j2) = -90.0_r8+(float(j2)-1.0_r8)*180.0_r8/(float(jm2)-1.0_r8)
  enddo


! set up staggered longitudes (cell edges in x)

  ! NOTE(jrvb): Subtracting half the size of a different cell as is
  ! done for slon(1) and slon(im1+1) will get the wrong value if the
  ! grid is not evenly spaced.  It would be more robust to calculate
  ! based on wrapping the end values and then averaging.
  do i1=2,im1
	slon1(i1) = (lon1(i1-1)+lon1(i1))/2.0_r8
  enddo
  slon1(1) = lon1(1)-(lon1(2)-lon1(1))/2.0_r8
  slon1(im1+1) = lon1(im1)+(lon1(im1)-lon1(im1-1))/2.0_r8
  ! NOTE(jrvb): for lon1(1) == 0, this is now negative
  ! NOTE(jrvb): this looks like it might break if slon1(1) is not close to slon(im1+1) - 360

  ! NOTE(jrvb): similar wrapping issues here
  do i2=2,im2
	slon2(i2) = (lon2(i2-1)+lon2(i2))/2.0_r8
  enddo
  slon2(1) = lon2(1)-(lon2(2)-lon2(1))/2.0_r8
  slon2(im2+1) = lon2(im2)+(lon2(im2)-lon2(im2-1))/2.0_r8

! set up staggered lattiudes (cell edges in y)
  slat1(1)=-90.0_r8
  do j1=2,jm1
	slat1(j1) = (lat1(j1-1)+lat1(j1))/2.0_r8
  enddo
  slat1(jm1+1)=90.0_r8

  slat2(1)=-90.0_r8
  do j2=2,jm2
	slat2(j2)=(lat2(j2-1)+lat2(j2))/2.0_r8
  enddo	
  slat2(jm2+1)=90.0_r8

! compute Guassian weight for two meshes (discrete form of cos(lat).)
  do j1=1,jm1
   	gw1(j1) = sin(slat1(j1+1)/180.0_r8*SHR_CONST_PI)-sin(slat1(j1)/180.0_r8*SHR_CONST_PI)
  enddo

  do j2=1,jm2
	gw2(j2) = sin(slat2(j2+1)/180.0_r8*SHR_CONST_PI)-sin(slat2(j2)/180.0_r8*SHR_CONST_PI)
  enddo


! add 360 to slon1 and slon2 

  ! NOTE(jrvb): This doesn't seem useful.  It is possible this was
  ! done with the thought that avoiding negative values would fix
  ! wrapping issues over the prime meridian below, but this doesn't do
  ! the trick.
  slon1(:) = slon1(:)+360.0_r8
  slon2(:) = slon2(:)+360.0_r8

 do i2=1,im2

! target grid east-west boundaries
  x2_west=slon2(i2)
  x2_east=slon2(i2+1)

  do i1=1,im1

! source grid east-west boundaries
   x1_west=slon1(i1)
   x1_east=slon1(i1+1)

! check if there is any overlap between the source grid and the target grid
! if no overlap, then weighting is zero
! there are three scenarios overlaps can take place 

   ! NOTE(jrvb): longitude wrap is not handled correctly below.  Eg,
   !    x1 (west, east) = (350, 10)   (x1_west was -10, added 360)
   !    x2 (west, east) = (2, 3)
   ! x2 is wholely contained in x1, but the routines here don't catch it
   if( (x1_west>=x2_west).and.(x1_east<=x2_east) ) then
! case 1: 
!                x1_west             x1_east
!                  |-------------------|
!            |---------------------------------|
!          x2_west                           x2_east
    if(use_flight_distance) then
     weight_x(i2,i1) = 1.0_r8
    else
     weight_x(i2,i1) =  (x1_east-x1_west)/(x2_east-x2_west)
    endif
   elseif ( (x1_west>=x2_west).and.(x1_west<x2_east) ) then
! case 2: 
!                x1_west                          x1_east
!                  |--------------------------------|
!            |---------------------------------|
!          x2_west                           x2_east
   if(use_flight_distance) then
     weight_x(i2,i1) = (x2_east-x1_west)/(x1_east-x1_west)
   else
     weight_x(i2,i1) = (x2_east-x1_west)/(x2_east-x2_west)
   endif
   elseif ( (x1_east>x2_west).and.(x1_east<=x2_east) ) then
! case 3: 
!       x1_west                          x1_east
!         |--------------------------------|
!                |---------------------------------|
!              x2_west                           x2_east
   if(use_flight_distance) then
     weight_x(i2,i1) = (x1_east-x2_west)/(x1_east-x1_west)
   else
     weight_x(i2,i1) = (x1_east-x2_west)/(x2_east-x2_west)
   endif
   elseif ( (x1_east>x2_east).and.(x1_west<x2_west) ) then
! case 4: 
!       x1_west                          x1_east
!         |--------------------------------|
!                |--------------------|
!              x2_west              x2_east
     weight_x(i2,i1) = 1._r8
      ! NOTE(jrvb): This looks broken for the use_flight_distance=TRUE
      ! case.  Only the portion of x1 which overlaps x2 should be
      ! added since the value is a total not a mixing ratio.
   endif

   enddo	
  enddo


! consider end points
      if(slon1(im1+1)>slon2(im2+1)) then
! case 1:
!           slon1(im1)                slon1(im1+1) <--- end point
!              |-------------------------|
!           |----------------|......................|
!        slon2(im2)         slon2(im2+1)        slon2(2)  (note: slon2(im2+1) = slon2(1))
         if(use_flight_distance) then
            weight_x(1,im1)= weight_x(1,im1)+(slon1(im1+1)-slon2(im2+1))/(slon1(im1+1)-slon1(im1))
         else
            weight_x(1,im1)= weight_x(1,im1)+(slon1(im1+1)-slon2(im2+1))/(slon2(2)-slon2(1))
         endif
      endif	

      if(slon1(im1+1)<slon2(im2+1)) then
! case 1:
!           slon1(im1)                slon1(im1+1)                  slon1(2)    (note: slon1(im1+1) = slon1(1))
!              |-------------------------|.............................|
!                   |-------------------------------|
!               slon2(im2)                        slon2(im2+1) <--- end point
         if(use_flight_distance) then
            weight_x(im2,1) = weight_x(im2,1)+(slon2(1)-slon1(1))/(slon1(2)-slon1(1))
         else
            weight_x(im2,1) = weight_x(im2,1)+(slon2(1)-slon1(1))/(slon2(2)-slon2(1)) 
         endif
      endif



      do j2=1,jm2
! target grid north-south boundaries
         y2_south=slat2(j2)
         y2_north=slat2(j2+1)

         do j1=1,jm1

! source grid north-south boundaries
            y1_south=slat1(j1)
            y1_north=slat1(j1+1)
 
! check if there is any overlap between the source grid and the target grid
! if no overlap, then weighting is zero
! there are three scenarios overlaps can take place 
! note: there is Guassian weight to consider in the meridional direction!

            if( (y1_south>=y2_south).and.(y1_north<=y2_north) ) then
! case 1: 
!                y1_south             y1_north
!                  |-------------------|
!            |---------------------------------|
!          y2_south                           y2_north
                if(use_flight_distance) then
                   weight_y(j2,j1) =  1.0_r8
                else
                   weight_y(j2,j1) =  gw1(j1)/gw2(j2)
                endif
            elseif ( (y1_south>=y2_south).and.(y1_south<y2_north) ) then
! case 2: 
!                y1_south                          y1_north
!                  |--------------------------------|
!            |---------------------------------|
!          y2_south                           y2_north
                if(use_flight_distance) then
                   weight_y(j2,j1) = (y2_north-y1_south)/(y1_north-y1_south)
                else
                   ! NOTE(jrvb): This looks wrong.  Assuming the gw*
                   ! values are the surface areas of the bands around
                   ! the sphere, we don't want the ratio of these two.
                   ! We want the ratio of the overlap area to the y2
                   ! area.  The linear factor here will not give us
                   ! that.
                   weight_y(j2,j1) = (y2_north-y1_south)/(y1_north-y1_south)*gw1(j1)/gw2(j2)
                endif
            elseif ( (y1_north>y2_south).and.(y1_north<=y2_north) ) then
! case 3: 
!       y1_south                          y1_north
!         |--------------------------------|
!                |---------------------------------|
!              y2_south                           y2_north
                if(use_flight_distance) then
                   weight_y(j2,j1) = (y1_north-y2_south)/(y1_north-y1_south)
                else
                   weight_y(j2,j1) = (y1_north-y2_south)/(y1_north-y1_south)*gw1(j1)/gw2(j2)
                endif
            elseif ( (y1_north>y2_north).and.(y1_south<y2_south) ) then
! case 4: 
!       y1_south                          y1_north
!         |--------------------------------|
!                |---------------------|
!              y2_south             y2_north
                if(use_flight_distance) then
                   ! NOTE(jrvb): This is wrong.  weight should be relative to overlap size, since not 100% of the values here belong
                   weight_y(j2,j1) = 1._r8
                else
                   ! NOTE(jrvb): the division here is inverted. The
                   ! weight should be the area of the overlap divided
                   ! by the area of the original, or gw2(j2)/gw1(j1)
                   ! in this case.
                   weight_y(j2,j1) = gw1(j1)/gw2(j2)
                endif
            endif

          enddo
        enddo

  end subroutine old_xy_interp_init

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

    ! NOTE: the calling routines sometimes feed in values that are almost identical, eg
    !    input_bot = 0.0,  input_top=2.853526808928757E-315
    ! For some reason noralize_lat turns the second into zero, maybe
    ! due to some floating point type conversion?  These are weird and
    ! broken-seeming inputs anyway, but hopefully allowing for zero
    ! intervals will resolve this?
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

  if (env_var_equals('CESM_USE_OLD_HORIZONTAL_INTERP', '1')) then
     write(iulog,*) 'CESM_USE_OLD_HORIZONTAL_INTERP=1: using old interp'
     call old_xy_interp_init(num_input_lons, num_input_lats, input_lon_radians, input_lat_radians, &
          num_sim_lons, num_sim_lats, weight_x, weight_y, use_flight_distance)
     return
  endif

  ! TODO(jrvb): code below definitely assumes input lat,lon values are
  ! sorted and bounds are as expected.  Add invariant checks here.

  weight_x(:,:) = 0.0_r8
  weight_y(:,:) = 0.0_r8

! input_lon_radians & input_lat_radians are longitude & latitude on the source mesh in radians
! convert input_lon, input_lat from radians to degrees

  ! NOTE(jrvb): it looks like the intent is to support arbitrary
  ! spacing of lat,lon on the input side.  Presumably there is an
  ! unstated assumption about things being sorted.
  input_lon(:) = input_lon_radians(:)/SHR_CONST_PI*180.0_r8
  input_lat(:) = input_lat_radians(:)/SHR_CONST_PI*180.0_r8

  ! Set up sim_lon, sim_lat (target mesh), in CAM convention.  The
  ! (lon,lat) pairs are the center points of the grid cells.
  ! TODO(jrvb): Why not do these in radians?
  do i2=1,num_sim_lons
     sim_lon(i2) = (float(i2)-1.0_r8)*360.0_r8/float(num_sim_lons)
  enddo
  do j2=1,num_sim_lats
     sim_lat(j2) = -90.0_r8+(float(j2)-1.0_r8)*180.0_r8/(float(num_sim_lats)-1.0_r8)
  enddo
  ! make sure the highest value is exactly +90, since the
  ! multiplication above will could slightly off due to rounding.
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
  ! NOTE(jrvb): this will correctly set up 1/2 height cells at the poles, as used in the CESM grid
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


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Test functions.  These should be moved to a separate unit test main prog.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_almost_equal(val1, val2, message)
    real(r8), intent(in)  :: val1, val2
    character(len=*), intent(in) :: message
    character(len=1000) :: out

    if ((val2 < 0.99999 * val1) .or. (val2 > 1.00001 * val1)) then
       ! call endrun(val1//' !~ '//val2//'  '//message)
       write(out, *) 'CHECK FAILED: ', val1, ' !~ ', val2, '  ', message
       call endrun(out)
       ! write(iulog, *) trim(out)
    endif
  end subroutine check_almost_equal

  ! real(r8) function normalize_lon_right(left, right)
  subroutine test_normalize_lon_right()
    call check_almost_equal(normalize_lon_right(10.0_r8, 20.0_r8), 20.0_r8, "test_normalize_lon_right:  right 1")
    call check_almost_equal(normalize_lon_right(-10.0_r8, 20.0_r8), 20.0_r8, "test_normalize_lon_right:  right 2")
    call check_almost_equal(normalize_lon_right(10.0_r8, 380.0_r8), 20.0_r8, "test_normalize_lon_right:  right 3")
    call check_almost_equal(normalize_lon_right(10.0_r8, 5.0_r8), 365.0_r8, "test_normalize_lon_right:  left 1")
    call check_almost_equal(normalize_lon_right(10.0_r8, -355.0_r8), 365.0_r8, "test_normalize_lon_right:  left 2")
    call check_almost_equal(normalize_lon_right(10.0_r8, -30.0_r8), 330.0_r8, "test_normalize_lon_right:  left 3")
    call check_almost_equal(normalize_lon_right(10.0_r8, 10.0_r8), 370.0_r8, "test_normalize_lon_right:  equal 1")
    call check_almost_equal(normalize_lon_right(10.0_r8, 370.0_r8), 370.0_r8, "test_normalize_lon_right:  equal 2")
    call check_almost_equal(normalize_lon_right(10.0_r8, -350.0_r8), 370.0_r8, "test_normalize_lon_right:  equal 3")
    write(iulog,*) 'test_normalize_lon_length: succeeded'
  end subroutine test_normalize_lon_right


  subroutine test_lon_length()
    call check_almost_equal(lon_length(0.0_r8, 10.0_r8), 10.0_r8, "test_lon_length:  basic")
    call check_almost_equal(lon_length(-20.0_r8, 10.0_r8), 30.0_r8, "test_lon_length:  negative value")
    call check_almost_equal(lon_length(340.0_r8, 10.0_r8), 30.0_r8, "test_lon_length:  wrapped")
    call check_almost_equal(lon_length(0.0_r8, 0.0_r8), 360.0_r8, "test_lon_length:  equal 1")
    call check_almost_equal(lon_length(-10.0_r8, 710.0_r8), 360.0_r8, "test_lon_length:  equal 2")
    write(iulog,*) 'test_lon_length: succeeded'
  end subroutine test_lon_length


  ! real(r8) function normalize_lon_value(lon)
  subroutine test_normalize_lon_value()
    call check_almost_equal(normalize_lon_value(-10.0_r8), 350.0_r8, "test_normalize_lon_value:  negative 1")
    call check_almost_equal(normalize_lon_value(-710.0_r8), 10.0_r8, "test_normalize_lon_value:  negative 2")
    call check_almost_equal(normalize_lon_value(0.0_r8), 0.0_r8, "test_normalize_lon_value:  zero 1")
    call check_almost_equal(normalize_lon_value(720.0_r8), 0.0_r8, "test_normalize_lon_value:  zero 2")
    call check_almost_equal(normalize_lon_value(20.0_r8), 20.0_r8, "test_normalize_lon_value:  positive 1")
    call check_almost_equal(normalize_lon_value(380.0_r8), 20.0_r8, "test_normalize_lon_value:  positive 2")
    write(iulog,*) 'test_normalize_lon_value: succeeded'
  end subroutine test_normalize_lon_value


  subroutine test_calculate_lon_overlap()
    call check_almost_equal(calculate_lon_overlap(10.0_r8, 20.0_r8, 30.0_r8, 40.0_r8), 0.0_r8, "test_calculate_lon_overlap:  disjoint")
    call check_almost_equal(calculate_lon_overlap(10.0_r8, 40.0_r8, 20.0_r8, 30.0_r8), 10.0_r8, "test_calculate_lon_overlap:  a contains b")
    call check_almost_equal(calculate_lon_overlap(20.0_r8, 30.0_r8, 10.0_r8, 40.0_r8), 10.0_r8, "test_calculate_lon_overlap:  b contains a")
    call check_almost_equal(calculate_lon_overlap(20.0_r8, 30.0_r8, 25.0_r8, 35.0_r8), 5.0_r8, "test_calculate_lon_overlap:  a left of b")
    call check_almost_equal(calculate_lon_overlap(20.0_r8, 30.0_r8, 15.0_r8, 25.0_r8), 5.0_r8, "test_calculate_lon_overlap:  a right of b")
    call check_almost_equal(calculate_lon_overlap(-5.0_r8, 5.0_r8, 358.0_r8, 359.0_r8), 1.0_r8, "test_calculate_lon_overlap:  wrapped 1")
    call check_almost_equal(calculate_lon_overlap(-5.0_r8, 5.0_r8, -362.0_r8, 370.0_r8), 7.0_r8, "test_calculate_lon_overlap:  wrapped 2")
    call check_almost_equal(calculate_lon_overlap(5.0_r8, 30.0_r8, -10.0_r8, 10.0_r8), 5.0_r8, "test_calculate_lon_overlap:  double wrap 1")
    write(iulog,*) 'test_calculate_lon_overlap: succeeded'
  end subroutine test_calculate_lon_overlap


  ! real(r8) function lon_weight(input_left, input_right, sim_left, sim_right, use_flight_distance)
  subroutine test_lon_weight()
    call check_almost_equal(lon_weight(-5.0_r8, 5.0_r8, 10.0_r8, 20.0_r8, .false.), 0.0_r8, "test_lon_weight:  zero 1")
    call check_almost_equal(lon_weight(-5.0_r8, 5.0_r8, 5.0_r8, 20.0_r8, .false.), 0.0_r8, "test_lon_weight:  zero 2")
    call check_almost_equal(lon_weight(-5.0_r8, 5.0_r8, -15.0_r8, -5.0_r8, .false.), 0.0_r8, "test_lon_weight:  zero 3")
    call check_almost_equal(lon_weight(-5.0_r8, 5.0_r8, -355.0_r8, 20.0_r8, .false.), 0.0_r8, "test_lon_weight:  zero 4")
    call check_almost_equal(lon_weight(-5.0_r8, 5.0_r8, -15.0_r8, -365.0_r8, .false.), 0.0_r8, "test_lon_weight:  zero 5")

    ! data is a mixing ratio; the weight should be the overlap / sim cell size
    call check_almost_equal(lon_weight(10.0_r8, 40.0_r8, 10.0_r8, 40.0_r8, .false.), 1.0_r8, "test_lon_weight:  mixing ratio; a == b")
    call check_almost_equal(lon_weight(20.0_r8, 30.0_r8, 15.0_r8, 35.0_r8, .false.), 0.5_r8, "test_lon_weight:  mixing ratio; a in b, 1")
    call check_almost_equal(lon_weight(-2.0_r8, -1.0_r8, -5.0_r8, 5.0_r8, .false.), 0.1_r8, "test_lon_weight:  mixing ratio; a in b, 2")
    call check_almost_equal(lon_weight(15.0_r8, 35.0_r8, 20.0_r8, 30.0_r8, .false.), 1.0_r8, "test_lon_weight:  mixing ratio; b in a, 1")
    call check_almost_equal(lon_weight(-5.0_r8, 5.0_r8, -2.0_r8, -1.0_r8, .false.), 1.0_r8, "test_lon_weight:  mixing ratio; b in a, 2")
    call check_almost_equal(lon_weight(10.0_r8, 30.0_r8, 25.0_r8, 35.0_r8, .false.), 0.5_r8, "test_lon_weight:  mixing ratio; a on left, 1")
    call check_almost_equal(lon_weight(-10.0_r8, -2.5_r8, -5.0_r8, 5.0_r8, .false.), 0.25_r8, "test_lon_weight:  mixing ratio; a on left, 2")
    call check_almost_equal(lon_weight(25.0_r8, 50.0_r8, 20.0_r8, 30.0_r8, .false.), 0.5_r8, "test_lon_weight:  mixing ratio; a on right, 1")
    call check_almost_equal(lon_weight(5.0_r8, 30.0_r8, -10.0_r8, 10.0_r8, .false.), 0.25_r8, "test_lon_weight:  mixing ratio; a on right, 2")

    ! data is a total for the cell; the weight should be the overlap / input cell size
    call check_almost_equal(lon_weight(10.0_r8, 40.0_r8, 10.0_r8, 40.0_r8, .true.), 1.0_r8, "test_lon_weight:  total; a == b")
    call check_almost_equal(lon_weight(20.0_r8, 30.0_r8, 10.0_r8, 40.0_r8, .true.), 1.0_r8, "test_lon_weight:  total; a in b, 1")
    call check_almost_equal(lon_weight(-2.0_r8, -1.0_r8, -5.0_r8, 5.0_r8, .true.), 1.0_r8, "test_lon_weight:  total; a in b, 2")
    call check_almost_equal(lon_weight(15.0_r8, 35.0_r8, 20.0_r8, 30.0_r8, .true.), 0.5_r8, "test_lon_weight:  total; b in a, 1")
    call check_almost_equal(lon_weight(-5.0_r8, 5.0_r8, -2.0_r8, -1.0_r8, .true.), 0.1_r8, "test_lon_weight:  total; b in a, 2")
    call check_almost_equal(lon_weight(10.0_r8, 30.0_r8, 25.0_r8, 35.0_r8, .true.), 0.25_r8, "test_lon_weight:  total; a on left, 1")
    call check_almost_equal(lon_weight(-7.5_r8, -2.5_r8, -5.0_r8, 5.0_r8, .true.), 0.5_r8, "test_lon_weight:  total; a on left, 2")
    call check_almost_equal(lon_weight(25.0_r8, 35.0_r8, 5.0_r8, 30.0_r8, .true.), 0.5_r8, "test_lon_weight:  total; a on right, 1")
    call check_almost_equal(lon_weight(5.0_r8, 25.0_r8, -30.0_r8, 10.0_r8, .true.), 0.25_r8, "test_lon_weight:  total; a on right, 2")

    write(iulog,*) 'test_lon_weight: succeeded'
  end subroutine test_lon_weight

  ! real(r8) function lat_band_weight(low, high)
  subroutine test_lat_band_weight()
    ! TODO(jrvb): probably nix this one; not really anything to test
    write(iulog,*) 'test_lat_band_weight: succeeded'
  end subroutine test_lat_band_weight

  ! real(r8) function lat_weight(input_bot, input_top, sim_bot, sim_top, use_flight_distance)
  subroutine test_lat_weight()
    call check_almost_equal(lat_weight(-5.0_r8, 5.0_r8, 10.0_r8, 20.0_r8, .false.), 0.0_r8, "test_lat_weight:  zero 1f")
    call check_almost_equal(lat_weight(-5.0_r8, 5.0_r8, 10.0_r8, 20.0_r8, .true.), 0.0_r8, "test_lat_weight:  zero 1t")
    call check_almost_equal(lat_weight(-5.0_r8, 5.0_r8, 5.0_r8, 20.0_r8, .false.), 0.0_r8, "test_lat_weight:  zero 2f")
    call check_almost_equal(lat_weight(-5.0_r8, 5.0_r8, 5.0_r8, 20.0_r8, .true.), 0.0_r8, "test_lat_weight:  zero 2t")
    call check_almost_equal(lat_weight(-5.0_r8, 5.0_r8, -15.0_r8, -5.0_r8, .false.), 0.0_r8, "test_lat_weight:  zero 3f")
    call check_almost_equal(lat_weight(-5.0_r8, 5.0_r8, -15.0_r8, -5.0_r8, .true.), 0.0_r8, "test_lat_weight:  zero 3t")
    call check_almost_equal(lat_weight(85.0_r8, 90.0_r8, 75.4_r8, 76.3_r8, .false.), 0.0_r8, "test_lat_weight:  zero 4f")
    call check_almost_equal(lat_weight(85.0_r8, 90.0_r8, 75.4_r8, 76.3_r8, .true.), 0.0_r8, "test_lat_weight:  zero 4t")

    ! check top and bottom
    call check_almost_equal(lat_weight(85.0_r8, 90.0_r8, 89.5_r8, 90.0_r8, .false.), 1.0_r8, "test_lat_weight:  b in a, f, top")
    call check_almost_equal(lat_weight(85.0_r8, 90.0_r8, 88.6_r8, 89.5_r8, .false.), 1.0_r8, "test_lat_weight:  b in a, f, near top")
    call check_almost_equal(lat_weight(-90.0_r8, -85.0_r8, -90.0_r8, -89.5_r8, .false.), 1.0_r8, "test_lat_weight:  b in a, f, bottom")
    call check_almost_equal(lat_weight(-90.0_r8, -85.0_r8, -89.5_r8, -88.6_r8, .false.), 1.0_r8, "test_lat_weight:  b in a, f, near bottom")
    call check_almost_equal(lat_weight(85.0_r8, 90.0_r8, 89.5_r8, 90.0_r8, .true.), 0.01000628511133214_r8, "test_lat_weight:  b in a, t, top")
    call check_almost_equal(lat_weight(85.0_r8, 90.0_r8, 88.6_r8, 89.5_r8, .true.), 0.06843958489161171_r8, "test_lat_weight:  b in a, t, near top")
    call check_almost_equal(lat_weight(-90.0_r8, -85.0_r8, -90.0_r8, -89.5_r8, .true.), 0.01000628511133214_r8, "test_lat_weight:  b in a, t, bottom")
    call check_almost_equal(lat_weight(-90.0_r8, -85.0_r8, -89.5_r8, -88.6_r8, .true.), 0.06843958489161171_r8, "test_lat_weight:  b in a, t, near bottom")

    ! data is a mixing ratio; the weight should be the overlap / sim cell size
    call check_almost_equal(lat_weight(10.0_r8, 40.0_r8, 10.0_r8, 40.0_r8, .false.), 1.0_r8, "test_lat_weight:  mixing ratio; a == b")
    call check_almost_equal(lat_weight(20.0_r8, 30.0_r8, 15.0_r8, 35.0_r8, .false.), 0.5019099187716736_r8, "test_lat_weight:  mixing ratio; a in b, 1")
    call check_almost_equal(lat_weight(-2.0_r8, -1.0_r8, -5.0_r8, 5.0_r8, .false.), 0.10009145533721157_r8, "test_lat_weight:  mixing ratio; a in b, 2")
    call check_almost_equal(lat_weight(15.0_r8, 35.0_r8, 20.0_r8, 30.0_r8, .false.), 1.0_r8, "test_lat_weight:  mixing ratio; b in a, 1")
    call check_almost_equal(lat_weight(-5.0_r8, 5.0_r8, -2.0_r8, -1.0_r8, .false.), 1.0_r8, "test_lat_weight:  mixing ratio; b in a, 2")
    call check_almost_equal(lat_weight(10.0_r8, 30.0_r8, 25.0_r8, 35.0_r8, .false.), 0.5126038285706509_r8, "test_lat_weight:  mixing ratio; a on left, 1")
    call check_almost_equal(lat_weight(-10.0_r8, -2.5_r8, -5.0_r8, 5.0_r8, .false.), 0.24976182870917_r8, "test_lat_weight:  mixing ratio; a on left, 2")
    call check_almost_equal(lat_weight(25.0_r8, 50.0_r8, 20.0_r8, 30.0_r8, .false.), 0.48982028397974603_r8, "test_lat_weight:  mixing ratio; a on right, 1")
    call check_almost_equal(lat_weight(5.0_r8, 30.0_r8, -10.0_r8, 10.0_r8, .false.), 0.24904504061416316_r8, "test_lat_weight:  mixing ratio; a on right, 2")

    ! data is a total for the cell; the weight should be the overlap / input cell size
    call check_almost_equal(lat_weight(10.0_r8, 40.0_r8, 10.0_r8, 40.0_r8, .true.), 1.0_r8, "test_lat_weight:  total; a == b")
    call check_almost_equal(lat_weight(20.0_r8, 30.0_r8, 10.0_r8, 40.0_r8, .true.), 1.0_r8, "test_lat_weight:  total; a in b, 1")
    call check_almost_equal(lat_weight(-2.0_r8, -1.0_r8, -5.0_r8, 5.0_r8, .true.), 1.0_r8, "test_lat_weight:  total; a in b, 2")
    call check_almost_equal(lat_weight(15.0_r8, 35.0_r8, 20.0_r8, 30.0_r8, .true.), 0.5019099187716736_r8, "test_lat_weight:  total; b in a, 1")
    call check_almost_equal(lat_weight(-5.0_r8, 5.0_r8, -2.0_r8, -1.0_r8, .true.), 0.10009145533721157_r8, "test_lat_weight:  total; b in a, 2")
    call check_almost_equal(lat_weight(10.0_r8, 30.0_r8, 25.0_r8, 35.0_r8, .true.), 0.23711140236969752_r8, "test_lat_weight:  total; a on left, 1")
    call check_almost_equal(lat_weight(-7.5_r8, -2.5_r8, -5.0_r8, 5.0_r8, .true.), 0.5009545047145969_r8, "test_lat_weight:  total; a on left, 2")
    call check_almost_equal(lat_weight(25.0_r8, 35.0_r8, 5.0_r8, 30.0_r8, .true.), 0.5126038285706509_r8, "test_lat_weight:  total; a on right, 1")
    call check_almost_equal(lat_weight(5.0_r8, 25.0_r8, -30.0_r8, 10.0_r8, .true.), 0.2578303983971048_r8, "test_lat_weight:  total; a on right, 2")

    ! ncall check_almost_equal(lat_weight(-90.0_r8, -85.0_r8, --85.8_r8, -84.8_r8, .true.), 0.2578303983971048_r8, "test_lat_weight:  total; a on right, 2")
    ! real(r8) function lat_weight(input_bot, input_top, sim_bot, sim_top, use_flight_distance)

    write(iulog,*) 'test_lat_weight: succeeded'
  end subroutine test_lat_weight

  ! subroutine xy_interp_init(num_input_lons, num_input_lats, input_lon_radians, input_lat_radians, &
  !      num_sim_lons, num_sim_lats, weight_x, weight_y, use_flight_distance)
  subroutine test_xy_interp_init()
    write(iulog,*) 'test_xy_interp_init: succeeded'
  end subroutine test_xy_interp_init

  ! TODO(jrvb): remove for final PR
  subroutine xy_interp_run_unit_tests()
    ! Run all tests
    call test_normalize_lon_right()
    call test_lon_length()
    call test_normalize_lon_value()
    call test_calculate_lon_overlap()
    call test_lon_weight()
    ! call test_lat_band_weight()
    call test_lat_weight()
    ! call test_xy_interp_init()
    ! call test_xy_interp()

    ! exit, as a reminder to remove the calls to this code.
    write(iulog,*) '############################################################'
    write(iulog,*) '##           horizontal unit tests passed!!               ##'
    write(iulog,*) '############################################################'
  end subroutine xy_interp_run_unit_tests


end module horizontal_interpolate
