! NCLFORTSTART
      subroutine mask_ice(plat, plon  , ice, landfrac)
!
!--------1---------2---------3---------4---------5--------6---------7--
!
! Use land mask to overwrite land fraction with ice fraction to 
! prepare for interpolation of ice fraction to new resolution
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer plat      ! latitude  dimension of input/output fields
      integer plon      ! longitude dimension of input/output fields
!
      real*8 ice      (plon,plat ) ! input/output ice fraction
      real*8 landfrac (plon,plat ) ! land mask
!
!-----------------------------------------------------------------------
!

! NCLEND

!---------------------------Local workspace-----------------------------
!
      integer i,ii,j
      integer icount
      real*8  tw                   ! on lat circle, ocn point west of land
      integer iw                   ! index of western ocn point
      integer iw_tmp               ! index of western ocn point
      real*8  te                   ! on lat circle, ocn point east of land
      integer ie                   ! index of eastern ocn point
      integer ie_tmp               ! index of eastern ocn point
      integer offshore_pts         ! how many points offshore to sample ice fraction
      real*8  land_threshold       ! minimum fraction to be considered land
      real*8  wgt_w, wgt_e         ! interpolation weights
!
      logical found_ocn_w          ! true when ocn point is found on lat 
!                                  ! circle west of current land point
      logical found_ocn_e          ! true when ocn point is found on lat 
!                                  ! circle east of current land point
      logical found_lnd            ! true when lnd point is found

      land_threshold = 0.0
      offshore_pts   = 3
!
      if(offshore_pts .lt. 1 .or. offshore_pts .gt. 10) then
         write(6,*) 'number of points offshore (in east/west direction)'
         write(6,*) 'to sample ice must be at least 1 but not more'
         write(6,*) 'than 10'
         write(6,*) 'offshore_pts = ', offshore_pts
         call abort
      end if
!
! Overwrite land fraction with interpolated ice fraction
!
      do j = 1, plat
!
         found_ocn_w = .false.
         found_ocn_e = .false.
         iw          = -999
         ie          = -999
!
         do i = 1,plon
!
! If found ocean, set tw/iw and move on to next point eastward
!
            if(landfrac(i,j) .le. land_threshold) then
               found_ocn_w = .true.
               iw          = i
!
! else, current point is land and will be overwritten
!
            else
!
! If "i=1", search westward over entire latitude circle for first
! ocean point west of current land point.
!
               if (i .eq. 1) then
                  found_ocn_w = .false.
                  ii = plon + 1
                  do while (ii .gt. 1 .and. .not. found_ocn_w)
                     ii = ii - 1
                     if(landfrac(ii,j) .le. land_threshold) then
                        iw          = ii-plon
                        found_ocn_w = .true.
                     end if
                  end do
               end if
!
! Now search eastward over entire latitude circle for first ocean point
! encountered east of current land point
!
               if(ie .lt. i .and. found_ocn_w) then
                  found_ocn_e = .false.
                  ii = i
                  ie = i
                  do while (ie .le. plon+i .and. .not. found_ocn_e)
                     ie = ie + 1
                     ii = ii + 1
                     if(ii .gt. plon) ii = 1
                     if(landfrac(ii,j) .le. land_threshold) then
                        found_ocn_e = .true.
                     end if
                  end do
               end if
!
! If no ocean points were found in this latitude circle, assume we
! are at Antarctica and set all ice fractions to 1.
!
               if      (.not. found_ocn_w .and. .not. found_ocn_e) then
                  ice(i,j) = 1.
!
! Else, interpolate ice fraction over land point
!
               else if (      found_ocn_w .and.       found_ocn_e) then
                  iw_tmp = iw
                  ie_tmp = ie
!
! First, decide how many points offshore (in both directions)  to
! sample ice fraction for interpolation over land
!
                  if (offshore_pts .gt. 1) then
!
! Eastern shore:  look east to make sure we are staying within half-way of
!                 the next land mass to the east
!
                     ii        = ie_tmp
                     icount    = 1
                     found_lnd = .false.
                     do while (icount .lt. 2*offshore_pts .and. &
                                                      .not. found_lnd)
                        icount = icount + 1
                        ii     = ii     + 1
                        if(ii .gt. plon) ii = ii - plon
                        if(landfrac(ii,j) .gt. land_threshold) then
                           found_lnd = .true.
                        end if
                     end do
                     ie_tmp = ie_tmp + icount/2 - 1
!
! Western shore:  look west to make sure we are staying within half-way of
!                 the next land mass to the west
!
                     ii        = iw_tmp
                     icount    = 1
                     found_lnd = .false.
                     do while (icount .lt. 2*offshore_pts .and. &
                                                      .not. found_lnd)
                        icount = icount + 1
                        ii     = ii     - 1
                        if(ii .lt.    1) ii = ii + plon
                        if(landfrac(ii,j) .gt. land_threshold) then
                           found_lnd = .true.
                        end if
                     end do
                     iw_tmp = iw_tmp - icount/2 + 1
                  end if
!
                  if(i .le. iw_tmp .or. i .ge. ie_tmp) then
                     write(6,*) 'Error:  Should never reach this'
                     write(6,*) 'branch.  Current land point should'
                     write(6,*) 'be between the east and west ocn'
                     write(6,*) 'points.'
                     write(6,*) 'iw, i, ie = ', iw_tmp, i, ie_tmp
                     call abort
                  end if
!
! get E/W ice fraction from offshore points
!
                  ii = iw_tmp
                  if(ii .lt.    1) ii = ii + plon
                  tw       = ice(ii,j)
                  ii = ie_tmp
                  if(ii .gt. plon) ii = ii - plon
                  te       = ice(ii,j)
!
                  wgt_w = float(ie_tmp -      i)/float(ie_tmp-iw_tmp)
                  wgt_e = float(i      - iw_tmp)/float(ie_tmp-iw_tmp)
                  if(wgt_w .gt. wgt_e) then
                     ice(i,j) = tw
                  else
                     ice(i,j) = te
                  end if
!                  ice(i,j) = tw*wgt_w + te*wgt_e
               else if (.not. found_ocn_w .and.       found_ocn_e) then
                  write(6,*) 'Error:  should never reach this branch.'
                  write(6,*) 'If even just one ocean point exists on '
                  write(6,*) 'a latitude circle, it should be found in'
                  write(6,*) 'both directions'
                  call abort
               end if
!
            end if
         end do
      end do
!
! Sanity checks
!
      do j = 1, plat
         do i = 1,plon
            if(ice     (i,j) .gt.   1. .or. ice     (i,j) .lt. 0.) then
               write(6,*) 'ice fraction out of bounds'
               write(6,*) 'ice fraction  = ', ice (i,j)
               write(6,*) 'at i,j        = ', i,j
               call abort
            end if
         end do
      end do
!
      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine mask_sst(plat, plon  , sst, landfrac, icefrac )
!
!--------1---------2---------3---------4---------5--------6---------7--
!
! Use land and ice masks to overwrite land fraction with SSTs to 
! prepare for interpolation of SSTs to new resolution
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer plat      ! latitude  dimension of input/output fields
      integer plon      ! longitude dimension of input/output fields
!
      real*8 sst      (plon,plat ) ! input:  sst (K); output: masked sst 
!                                  ! converted to degrees C
      real*8 landfrac (plon,plat ) ! land mask
      real*8 icefrac  (plon,plat ) ! ice  mask
!
!-----------------------------------------------------------------------
!

! NCLEND

!---------------------------Local workspace-----------------------------
!
      integer i,ii,j
      integer icount
      real*8 tfreez                ! freezing point of fresh water(K)
!                                  ! (C/K conversion factor)
      real*8 tfreezsw              ! freezing point of salt water (K)
      real*8  tw                   ! on lat circle, sst point west of land
      integer iw                   ! index of western sst point
      integer iw_tmp               ! index of western sst point
      real*8  te                   ! on lat circle, sst point east of land
      integer ie                   ! index of eastern sst point
      integer ie_tmp               ! index of eastern sst point
      integer offshore_pts         ! how many points offshore to sample SST's
      real*8  land_threshold       ! minimum fraction to be considered land
      real*8  wgt_w, wgt_e         ! interpolation weights
!
      logical found_ocn_w          ! true when ocn point is found on lat 
!                                  ! circle west of current land point
      logical found_ocn_e          ! true when ocn point is found on lat 
!                                  ! circle east of current land point
      logical found_lnd            ! true when lnd point is found

      tfreez         = 273.15
      tfreezsw       = tfreez - 1.8
      land_threshold = 0.0
      offshore_pts   = 1
!
      if(offshore_pts .lt. 1 .or. offshore_pts .gt. 10) then
         write(6,*) 'number of points offshore (in east/west direction)'
         write(6,*) 'to sample SSTs must be at least 1 but not more'
         write(6,*) 'than 10'
         write(6,*) 'offshore_pts = ', offshore_pts
         call abort
      end if
!
! For each grid box that has seaice, set sst to "tfreezsw"
!
      do j = 1, plat
         do i = 1,plon
            if(icefrac(i,j) .gt. 0.) then
               sst(i,j) = tfreezsw
            end if
         end do
      end do
!
! Now, set all T's to a minimum of "tfreezsw"
! (no, this is not being redundant)
!
      do j = 1, plat
         do i = 1,plon
            if(sst(i,j) .lt. tfreezsw) then
               sst(i,j) = tfreezsw
            end if
         end do
      end do
!
! Overwrite land Ts with interpolated SST's
!
      do j = 1, plat
!
         found_ocn_w = .false.
         found_ocn_e = .false.
         iw          = -999
         ie          = -999
!
         do i = 1,plon
!
! If found ocean, set tw/iw and move on to next point eastward
!
            if(landfrac(i,j) .le. land_threshold) then
               found_ocn_w = .true.
               iw          = i
!
! else, current point is land and will be overwritten
!
            else
!
! If "i=1", search westward over entire latitude circle for first
! ocean point west of current land point.
!
               if (i .eq. 1) then
                  found_ocn_w = .false.
                  ii = plon + 1
                  do while (ii .gt. 1 .and. .not. found_ocn_w)
                     ii = ii - 1
                     if(landfrac(ii,j) .le. land_threshold) then
                        iw          = ii-plon
                        found_ocn_w = .true.
                     end if
                  end do
               end if
!
! Now search eastward over entire latitude circle for first ocean point
! encountered east of current land point
!
               if(ie .lt. i .and. found_ocn_w) then
                  found_ocn_e = .false.
                  ii = i
                  ie = i
                  do while (ie .le. plon+i .and. .not. found_ocn_e)
                     ie = ie + 1
                     ii = ii + 1
                     if(ii .gt. plon) ii = 1
                     if(landfrac(ii,j) .le. land_threshold) then
                        found_ocn_e = .true.
                     end if
                  end do
               end if
!
! If no ocean points were found in this latitude circle, assume we
! are at Antarctica and set all SST's to seaice values
!
               if      (.not. found_ocn_w .and. .not. found_ocn_e) then
                  sst(i,j) = tfreezsw
!
! Else, interpolate SST's over land point
!
               else if (      found_ocn_w .and.       found_ocn_e) then
                  iw_tmp = iw
                  ie_tmp = ie
!
! First, decide how many points offshore (in both directions)  to
! sample SST's for interpolation over land
!
                  if (offshore_pts .gt. 1) then
!
! Eastern shore:  look east to make sure we are staying within half-way of
!                 the next land mass to the east
!
                     ii        = ie_tmp
                     icount    = 1
                     found_lnd = .false.
                     do while (icount .lt. 2*offshore_pts .and. &
                                                      .not. found_lnd)
                        icount = icount + 1
                        ii     = ii     + 1
                        if(ii .gt. plon) ii = ii - plon
                        if(landfrac(ii,j) .gt. land_threshold) then
                           found_lnd = .true.
                        end if
                     end do
                     ie_tmp = ie_tmp + icount/2 - 1
!
! Western shore:  look west to make sure we are staying within half-way of
!                 the next land mass to the west
!
                     ii        = iw_tmp
                     icount    = 1
                     found_lnd = .false.
                     do while (icount .lt. 2*offshore_pts .and. &
                                                      .not. found_lnd)
                        icount = icount + 1
                        ii     = ii     - 1
                        if(ii .lt.    1) ii = ii + plon
                        if(landfrac(ii,j) .gt. land_threshold) then
                           found_lnd = .true.
                        end if
                     end do
                     iw_tmp = iw_tmp - icount/2 + 1
                  end if
!
                  if(i .le. iw_tmp .or. i .ge. ie_tmp) then
                     write(6,*) 'Error:  Should never reach this'
                     write(6,*) 'branch.  Current land point should'
                     write(6,*) 'be between the east and west SST'
                     write(6,*) 'points.'
                     write(6,*) 'iw, i, ie = ', iw_tmp, i, ie_tmp
                     call abort
                  end if
!
! get E/W SSTs from offshore points
!
                  ii = iw_tmp
                  if(ii .lt.    1) ii = ii + plon
                  tw       = sst(ii,j)
                  ii = ie_tmp
                  if(ii .gt. plon) ii = ii - plon
                  te       = sst(ii,j)
!
                  wgt_w = float(ie_tmp -      i)/float(ie_tmp-iw_tmp)
                  wgt_e = float(i      - iw_tmp)/float(ie_tmp-iw_tmp)
                  if(wgt_w .gt. wgt_e) then
                     sst(i,j) = tw
                  else
                     sst(i,j) = te
                  end if
!                  sst(i,j) = tw*wgt_w + te*wgt_e
               else if (.not. found_ocn_w .and.       found_ocn_e) then
                  write(6,*) 'Error:  should never reach this branch.'
                  write(6,*) 'If even just one ocean point exists on '
                  write(6,*) 'a latitude circle, it should be found in'
                  write(6,*) 'both directions'
                  call abort
               end if
!
            end if
         end do
      end do
!
! Sanity checks
!
      do j = 1, plat
         do i = 1,plon
            if(sst(i,j) .gt. 340. .or. sst(i,j) .lt. 160.) then
               write(6,*) 'sst out of bounds'
               write(6,*) 'sst           = ', sst     (i,j)
               write(6,*) 'at i,j        = ', i,j
               call abort
            end if
         end do
      end do
!
      do j = 1, plat
         do i = 1,plon
            sst(i,j) = sst(i,j) - tfreez
         end do
      end do
!
      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine jra_25_press_full_levels(plevp1, plev   ,plat   , &
                                          plon  , psi    ,psm    )
!
!--------1---------2---------3---------4---------5--------6---------7--
!
! Interpolate JRA-25 full-level pressures from interface pressures
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer plevp1    ! vertical  dimension of interface  pressure
!                       !                              field (input )
      integer plev      ! vertical  dimension of full-level pressure
!                       !                              field (output)
      integer plat      ! latitude  dimension of input/output fields
      integer plon      ! longitude dimension of input/output fields
!
      real*8 psi (plon,plat,plevp1) ! input  interface  pressure field
      real*8 psm (plon,plat,plev  ) ! output full-level pressure field
!
!-----------------------------------------------------------------------
!

! NCLEND

!---------------------------Local workspace-----------------------------
!
      integer i,j,k
      real*8 one, two
!
!-----------------------------------------------------------------------
!
      one = 1.
      two = 2.
!
      do k = 2,plev
         do j = 1,plat
            do i = 1,plon
               psm(i,j,k) = ( psi(i,j,k+1) * log( psi(i,j,k+1) ) ) - &
                            ( psi(i,j,k  ) * log( psi(i,j,k  ) ) )
               psm(i,j,k) =   psm(i,j,k  )/(psi(i,j,k+1) - psi(i,j,k))
               psm(i,j,k) =   psm(i,j,k  ) - one
               psm(i,j,k) = exp( psm(i,j,k) )
            end do
         end do
      end do

      do j = 1,plat
         do i = 1,plon
            psm(i,j,1) = psi(i,j,2)/two
         end do
      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine binning(plev   ,plato   ,plono   ,plat    ,plon , &
                         xx     ,yy      ,clat    ,clon    ,gw   , &
                         clato  ,clono   ,gwo     ,bin_factor    , &
                         dyn_flag, dyn_flago      )
!
!--------1---------2---------3---------4---------5--------6---------7--
!
! Grid-Box Binning
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer plev      ! vertical dimension of input/output field
      integer plato     ! latitude dimension of output field
      integer plono     ! longitude dimension of output field
      integer plat      ! latitude dimension of input field
      integer plon      ! longitude dimension of input field
!
      real*8 xx     (plon ,plat ,plev) ! input analysis field
      real*8 yy     (plono,plato,plev) ! horizontally interpolated
!                                      ! (output) field
      real*8 clat (plat )              ! Input latitude in degrees
!                                      ! starting from southern-most lat
      real*8 clon (plon )              ! Input longitude in degrees
!                                      ! starting at 0 deg and moving
!                                      ! eastward
      real*8 gw   (plat )              ! Input Gaussian wgts (if relevant grid)
      real*8 clato(plato)              ! Output latitude in degrees
!                                      ! starting from southern-most lat
      real*8 clono(plono)              ! Output longitude in degrees
!                                      ! starting at 0 deg and moving
!                                      ! eastward
      real*8 gwo  (plato)              ! Output Gaussian wgts (if relevant grid)
      integer dyn_flag                 ! Dynamics flag of input grid:   FV=0
      integer dyn_flago                ! Dynamics flag of output grid:  FV=0
      real*8 bin_factor                ! bin-box area expansion/contraction factor relative to
!                                      ! output grid-box area.
!
!-----------------------------------------------------------------------
!

! NCLEND

!---------------------------Local workspace-----------------------------
!
      integer max_segs                          ! Max # of box segments
      parameter (max_segs = 10000 )               
      integer i, j, ii, jj, k, platp1, platop1  ! Indices
      integer nx, ny, nx_max, ny_max            ! Indices
      integer plon2, plonhalf
      integer i_in(max_segs),i_out(max_segs)
      integer j_in(max_segs),j_out(max_segs)
      real*8 sum
      real*8 xx_loc(plon*2, plat, plev)
      real*8 pi, pio180, pio2, one, factor
      real*8 flat  (plat)   ,flon  (plon*2  )
      real*8 flato (plato  ),flono (plono   )
      real*8 flati (plat+1), floni (plon*2+1)
      real*8 flatoi(plato+1),flonoi(plono+1  )
      real*8 tmps, tmpn, tmp(plono,plato)
      real*8 edge_w (plon*2), edge_e (plon*2), edge_s (plat ), &
                                               edge_n (plat )
      real*8 edgeo_w(plono ), edgeo_e(plono ), edgeo_s(plato), &
                                               edgeo_n(plato)
      real*8 sin_s (plat ),sin_n (plat )
      real*8 sino_s(plato),sino_n(plato)
      real*8 dx(max_segs), dy(max_segs)
      real*8 distmin, dist, two, zero
!
!-----------------------------------------------------------------------
!
      zero     = 0.
      platp1   = plat  + 1
      platop1  = plato + 1
      plon2    = plon*2
      plonhalf = plon/2
      one      = 1.
      pi       = 4.*atan(one)
      pio180   = pi/180.
      pio2     = pi/2.
      two      = 2.
!
! Sanity checks
!
      if(bin_factor .lt. 0.05) then
         write(6,*) 'ERROR ("BINNING"):  binning factor out of range'
         write(6,*) 'bin_factor = ', bin_factor
         call abort
      end if
!
! Copy input data to wrap-around array (wrap half-way around globe at each end of x-direction)
!
      do k = 1,plev
         do j = 1,plat
            ii = plonhalf
            do i = 1,plon2
               ii = ii + 1
               if(ii .gt. plon) ii = 1
               xx_loc(i,j,k) = xx(ii,j,k)
            end do
         end do
      end do
!
! Convert to radians (wrap half-way around globe at each end of x-direction
! of input grid)
!
      ii   = plonhalf
      do i = 1,plon2
         ii = ii + 1
         if(ii .gt. plon    )      ii = 1
         if(i  .le. plonhalf)      flon(i) = clon(ii)*pio180-4*pio2
         if(i  .gt. plonhalf+plon) flon(i) = clon(ii)*pio180+4*pio2
         if(i  .gt. plonhalf .and. i .le. plonhalf+plon) &
                                   flon(i) = clon(ii)*pio180
      end do

      do j = 1,plat
         flat (j) = clat (j)*pio180
      end do

      do i = 1,plono
         flono(i) = clono(i)*pio180
      end do

      do j = 1,plato
         flato(j) = clato(j)*pio180
      end do
!
! Compute box edges
!
      floni(       1) = 0.5*(flon(1) + flon(plon*2) - pio2*8)
      floni(plon*2+1) = 0.5*(flon(1) + flon(plon*2) + pio2*8)
      do i = 2,plon*2
         floni(i) = 0.5*(flon(i-1) + flon(i))
      end do

      sum = 0.
      flati(1     ) = -pio2
      flati(platp1) =  pio2
      do j = 1,plat
         if(dyn_flag .eq. 1) then
            sum = sum + gw (j)
            if(sum .gt. 2. .and. j .lt. plat) then
               write(6,*) 'ERROR ("BINNING"): Something wrong with'
               write(6,*) 'input Gaussian weights'
               call abort
            elseif(sum .gt. 2. .and. j .eq. plat) then
               sum = min (sum, two)
            end if
            flati(j+1) = asin( sum-one )
         elseif(dyn_flag .eq. 0) then
            if(j .lt. plat) flati(j+1) = 0.5*(flat(j) + flat(j+1))
         else
            write(6,*) 'ERROR:  "BINNING" should not reach this branch'
            call abort
         end if
      end do
      
      flonoi(      1) = 0.5*(flono(1) + flono(plono) - pio2*4)
      flonoi(plono+1) = 0.5*(flono(1) + flono(plono) + pio2*4)
      do i = 2,plono
         flonoi(i) = 0.5*(flono(i-1) + flono(i))
      end do

      sum = 0.
      flatoi(1      ) = -pio2
      flatoi(platop1) =  pio2
      sum = 0.
      do j = 1,plato
         if(dyn_flago .eq. 1) then
            sum = sum + gwo (j)
            if(sum .gt. 2. .and. j .lt. plato) then
               write(6,*) 'ERROR ("BINNING"): Something wrong with'
               write(6,*) 'output Gaussian weights'
               call abort
            elseif(sum .gt. 2. .and. j .eq. plato) then
               sum = min (sum, two)
            end if
            flatoi(j+1) = asin( sum-one )
         elseif(dyn_flago .eq. 0) then
            if(j .lt. plato) flatoi(j+1) = 0.5*(flato(j) + flato(j+1))
         else
            write(6,*) 'ERROR:  "BINNING" should not reach this branch'
            call abort
         end if
      end do
!
! Copy grid interfaces to "edge" arrays
!
      do i = 1,plon*2
        edge_w(i) = floni(i  )
        edge_e(i) = floni(i+1)
      end do

      do j = 1,plat
        edge_s(j) = flati(j  )
        edge_n(j) = flati(j+1)
        sin_s (j) = sin(edge_s(j))
        sin_n (j) = sin(edge_n(j))
      end do
!
! Expand/contract bin box area for each output grid box by "bin_factor"
!
      factor = sqrt(bin_factor)

      do i = 1,plono
        edgeo_w(i) = flono(i) - ( flono (i  ) - flonoi(i) )*factor
        edgeo_e(i) = flono(i) + ( flonoi(i+1) - flono (i) )*factor
      end do

      do j = 1,plato
        tmps       = flato(j) - ( flato (j  ) - flatoi(j) )*factor
        tmpn       = flato(j) + ( flatoi(j+1) - flato (j) )*factor
        edgeo_s(j) = max( tmps, -pio2) - max( ( tmpn - pio2), zero)
        edgeo_n(j) = min( tmpn,  pio2) + max( (-pio2 - tmps), zero)
        sino_s (j) = sin(edgeo_s(j))
        sino_n (j) = sin(edgeo_n(j))
      end do
!
! Make vector of box segments in x-direction
!
      nx = 0
      do i = 1,plono
         do ii = 1,plon*2
            if(edge_e (ii) .gt. edgeo_w( i) .and. &
               edgeo_e( i) .gt. edge_w (ii) ) then
               nx = nx + 1
               if(nx .gt. max_segs) then
                  write(6,*) 'ERROR  ("BINNING"):  number of box'
                  write(6,*) 'segments greater than "max_segs"'
                  call abort
               end if
               i_in (nx) = ii
               i_out(nx) = i
               dx   (nx) = min(min(min(edge_e(ii)-edge_w(ii), &
                                       edgeo_e(i)-edgeo_w(i) ), &
                                       edge_e(ii)-edgeo_w(i) ), &
                                       edgeo_e(i)-edge_w(ii) )
            end if
            if(edge_w (ii) .ge. edgeo_e( i)) exit
         end do
      end do
!
! Make vector of box segments in y-direction
!
      ny = 0
      do j = 1,plato
         do jj = 1,plat
            if(edge_n (jj) .gt. edgeo_s( j) .and. &
               edgeo_n( j) .gt. edge_s (jj) ) then
               ny = ny + 1
               if(ny .gt. max_segs) then
                  write(6,*) 'ERROR  ("BINNING"):  number of box'
                  write(6,*) 'segments greater than "max_segs"'
                  call abort
               end if
               j_in (ny) = jj
               j_out(ny) = j
               distmin   = edge_n(jj)-edge_s(jj)
               dy(ny)    = sin_n (jj)-sin_s (jj)
               dist      = edgeo_n(j)-edgeo_s(j)
               if(dist .lt. distmin) then
                  distmin = dist
                  dy(ny)  = sino_n(j)-sino_s(j)
               end if
               dist      = edge_n(jj)-edgeo_s(j)
               if(dist .lt. distmin) then
                  distmin = dist
                  dy(ny)  = sin_n(jj)-sino_s(j)
               end if
               dist      = edgeo_n(j)-edge_s(jj)
               if(dist .lt. distmin) then
                  distmin = dist
                  dy(ny)  = sino_n(j)-sin_s(jj)
               end if
            end if
            if(edge_s (jj) .ge. edgeo_n( j)) exit
         end do
      end do

      nx_max = nx
      ny_max = ny
!
! Begin weighted binning
!
      do k = 1,plev
         do j = 1,plato
            do i = 1,plono
               yy(i,j,k) = 0.
            end do
         end do
      end do

      do k = 1,plev
         do ny = 1,ny_max
            j  = j_out(ny)
            jj = j_in (ny)
            do nx = 1,nx_max
               i  = i_out(nx)
               ii = i_in (nx)
               yy(i,j,k) = yy(i,j,k) + xx_loc(ii,jj,k)*dx(nx)*dy(ny)
            end do
         end do
      end do
!
! Normalize
!
      do j = 1,plato
         do i = 1,plono
            tmp(i,j) = (edgeo_e(i) - edgeo_w(i))*(sino_n(j) - sino_s(j))
         end do
      end do

      do k = 1,plev
         do j = 1,plato
            do i = 1,plono
               yy(i,j,k) = yy(i,j,k)/tmp(i,j)
            end do
         end do
      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine cubic_opt1(plev    ,plato   ,plono   ,plat    ,plon , &
                            platm2  ,xx      ,yy      ,pext    ,xx_exts, &
                            clat    ,clon    ,clato   ,clono   ,limdr)
!
!--------1---------2---------3---------4---------5--------6---------7--
!
! Horizontal Cubic Interpolation Driver
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer plev      ! vertical dimension of input/output field
      integer plato     ! latitude dimension of output field
      integer plono     ! longitude dimension of output field
      integer plat      ! latitude dimension of input field
      integer plon      ! longitude dimension of input field
      integer platm2    ! input latitude dimension (less possible pole points)
      integer pext      ! # of latitude extensions

!
      real*8 xx     (plon ,plat ,plev) ! input analysis field
      real*8 yy     (plono,plato,plev) ! horizontally interpolated
!                                      ! (output) field
      real*8 xx_exts(plon ,pext ,plev) ! input latitude extensions
      real*8 clat (plat )              ! Input latitude in degrees
!                                      ! starting from southern-most lat
      real*8 clon (plon )              ! Input longitude in degrees
!                                      ! starting at 0 deg and moving
!                                      ! eastward
      real*8 clato(plato)              ! Output latitude in degrees
!                                      ! starting from southern-most lat
      real*8 clono(plono)              ! Output longitude in degrees
!                                      ! starting at 0 deg and moving
!                                      ! eastward
      logical limdr                    ! Flag to use SCM0 derivative estimate limiter

!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
!
      integer i, j, k       ! Indices
      integer nxpt          ! Extra point in local 4-pt interpolation
!                           ! stencil
      integer jintmx        ! number of extra latitudes at southern and
!                           ! northern extremes of input latitude grid
!                           ! (created for interpolation purposes)
      integer platd         ! Total latitude dimension of input grid
!                           ! including extensions
      integer plond         ! Total longitude dimension of input grid
!                           ! including extensions
!
      logical npole         ! flag set if northern-most input latitude is a pole point
      logical spole         ! flag set if southern-most input latitude is a pole point
      logical lpole         ! flag set if both input latitudes are pole points

      parameter (nxpt   = 1 )          ! set "nxpt"
      parameter (jintmx = 2 )          ! set "jintmx"

      real*8 xxm2 (plon ,platm2 ,plev) ! input analysis field (less possible pole points)

!
!-----------------------------------------------------------------------
!
      npole  = .false.
      spole  = .false.
      lpole  = .false.
      if(clat(   1) .lt. -89.9999) spole = .true.
      if(clat(plat) .gt.  89.9999) npole = .true.
      if(platm2 .ne. plat-2 ) then
         write(6,*) 'Error in CUBIC_OPT1:  '
         write(6,*) '"platm2" must be same as "plat-2" '
         call abort
      end if
      if(spole .and. .not. npole .or. .not. spole .and. npole ) then
         write(6,*) 'Error in CUBIC_OPT1:  '
         write(6,*) 'input data must either include BOTH pole points ', &
                    'or none'
         call abort
      end if

      platd = plat + 2*nxpt + 2*jintmx
      plond = plon + 1 + 2*nxpt

!
! If input data has pole points, strip them out for now
!
      if(spole .and. npole) then
         platd = platm2 + 2*nxpt + 2*jintmx
         lpole = .true.
         do k = 1,plev
            do j = 1,platm2
               do i = 1,plon
                  xxm2(i,j,k) = xx(i,j+1,k)
               end do
            end do
         end do
      end if
!
! Call cubic interpolator
!
      if (lpole) then
         call cubic_slav(nxpt    ,jintmx  ,platd   ,plond   , &
                         plev    ,plato   ,plono   ,platm2  ,plon , &
                         xxm2    ,yy      ,pext    ,xx_exts ,clat(2) , &
                         clon    ,clato   ,clono   ,limdr   )
      else
         call cubic_slav(nxpt    ,jintmx  ,platd   ,plond   , &
                         plev    ,plato   ,plono   ,plat    ,plon , &
                         xx      ,yy      ,pext    ,xx_exts ,clat , &
                         clon    ,clato   ,clono   ,limdr   )
      end if
!
      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      subroutine cubic_slav(nxpt    ,jintmx  ,platd   ,plond   , &
                            plev    ,plato   ,plono   ,plat    ,plon , &
                            xx      ,yy      ,pext    ,xx_exts ,clat , &
                            clon    ,clato   ,clono   ,limdr   )
!
!-----------------------------------------------------------------------
!
! Horizontal Cubic Interpolation
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer nxpt          ! Extra point in local 4-pt interpolation
!                           ! stencil
      integer jintmx        ! number of extra latitudes at southern and
!                           ! northern extremes of input latitude grid
!                           ! (created for interpolation purposes)
      integer platd         ! Total latitude dimension of input grid
!                           ! including extensions
      integer plond         ! Total longitude dimension of input grid
!                           ! including extensions
!
      integer plev      ! vertical dimension of input/output field
      integer plato     ! latitude dimension of output field
      integer plono     ! longitude dimension of output field
      integer plat      ! latitude dimension of input field
      integer plon      ! longitude dimension of input field
      integer pext      ! # of latitude extensions
!
      real*8 xx   (plon ,plat ,plev) ! input analysis field
      real*8 yy   (plono,plato,plev) ! horizontally interpolated
!                                    ! (output) field
      real*8 xx_exts(plon ,pext ,plev) ! input latitude extensions
      real*8 clat (plat )            ! Input latitude in degrees
!                                    ! starting from southern-most lat
      real*8 clon (plon )            ! Input longitude in degrees
!                                    ! starting at 0 deg and moving
!                                    ! eastward
      real*8 clato(plato)            ! Output latitude in degrees
!                                    ! starting from southern-most lat
      real*8 clono(plono)            ! Input longitude in degrees
!                                    ! starting at 0 deg and moving
!                                    ! eastward
      logical limdr                  ! Flag to use SCM0 derivative estimate limiter
!
!---------------------------Local workspace-----------------------------
!
      real*8 xxx   (plond ,platd ,plev) ! input field on extended grid
      real*8 dphi  (platd)              ! latitude intervals (radians)
      real*8 rdphi                      ! reciprocal of dphi for interpolant
                                        ! y-interval
      real*8 rdx                        ! recip of del-x
      real*8 phi   (platd)              ! latitudes on extended grid
!                                       ! (radians)
      real*8 lam   (plond)              ! longitudes on extended grid
!                                       ! radians
      real*8 dlam  (plond)              ! delta-lam
      real*8 phidp (plato)              ! latitudes of output field
!                                       ! (radians)
      real*8 lamdp (plono)              ! longitudes of output field
!                                       ! (radians)
      real*8 lbasdx(4,2,plono)          ! longitudinal basis functions
!                                       ! for interpolation
      real*8 lbasdy(4,2,plato)          ! latitudinal basis functions
!                                       ! for interpolation

      real*8 pi                         ! 3.1415926...
      real*8 pio180                     ! pi/180.
      real*8 pi2                        ! pi*2
      real*8 fintx(4)
      real*8 hs (plato)
      real*8 hn (plato)
      real*8 dhs(plato)
      real*8 dhn(plato)
      real*8 hl (plono)
      real*8 hr (plono)
      real*8 dhl(plono)
      real*8 dhr(plono)
!
      integer jdp  (plato)
      integer idp  (plono)
      integer istart                  ! index for first analysis long.
      integer istop                   ! index for last  analysis long.
      integer jstart                  ! index for first analysis lat.
      integer jstop                   ! index for last  analysis lat.
      integer i, j, k, n              ! Indices
      integer ii, jj
!
!-----------------------------------------------------------------------
!
! Define constants
!
      pi     = 4.*atan(1.)
      pio180 = pi/180.
      pi2    = pi*2.
      istart = nxpt   + 1
      jstart = nxpt   + 1 + jintmx
      istop  = istart - 1 + plon
      jstop  = jstart - 1 + plat
!
! Make sure latitude dimension of "xx_ext" is same as
! total number of extended latitudes
!
      if( (jstart-1)*2 .ne. pext) then
         write(6,*) 'Error in CUBIC_SLAV:  '
         write(6,*) '   dimension of "xx_ext not same as total number'
         write(6,*) '   of extened latitudes'
         call abort
      end if
!
! Define lats/lons (in radians) for input and output grids
!
      do j = 1,plat
        phi  (jstart-1+j) = clat (j)*pio180
      end do
      do i = 1,plon
        lam  (istart-1+i) = clon (i)*pio180
      end do
      do j = 1,plato
        phidp(j)          = clato(j)*pio180
      end do
      do i = 1,plono
        lamdp(i)          = clono(i)*pio180
        if(lamdp(i) .lt. 0.) lamdp(i) = lamdp(i) + pi2
      end do
!
! North and south poles.
!
      phi(jstart-1) = -pi/2.0
      phi(jstop +1) =  pi/2.0
!
! Extend Gauss latitudes below south pole so that the spacing above
! the pole is symmetric, and phi is decreasing, i.e., phi < -pi/2
!
      if( jstart .gt. 2 )then
        do j = 1,jstart-2
          phi(j) = -pi - phi(2*jstart-2-j)
        end do
      end if
!
! Analogously for Northern Hemisphere
!
      if( platd .gt. jstop+1 )then
        do j = jstop+2,platd
          phi(j) = pi - phi(2*jstop+2-j)
        end do
      end if
!
! Fill East/West extensions of lam array
!
      do i = 1,istart-1
        lam(i) = lam(i+istop-istart+1) - pi2
      end do
      do i = istop+1,plond
        lam(i) = lam(i-istop+istart-1) + pi2
      end do
!
! Compute delta values in X/Y
!
      do i = 1,plond-1
         dlam(i) = lam(i+1) - lam(i)
      end do
      do j = 1,platd-1
         dphi(j) = phi(j+1) - phi(j)
      end do
!
! Copy input to extended grid and fill extensions
!
      do k = 1,plev
        do j = 1,plat
          do i = 1,plon
            xxx(istart-1+i,jstart-1+j,k) = xx(i,j,k)
          end do
        end do
!
! fill N/S poles
!
        do i = 1,plon
           xxx(istart-1+i,jstart-1,k) = xx_exts(i,pext/2  ,k)
           xxx(istart-1+i,jstop +1,k) = xx_exts(i,pext/2+1,k)
        end do
!
! fill beyond N/S pole latitudes
!
        if( jstart-2 .ge. 1 )then
           do j=1,jstart-2
              do i = 1,plon
                 xxx(istart-1+i,j,k) = xx_exts(i,j,k)
              end do
           end do
        end if

        if( jstop+2 .le. platd ) then
           jj = 1
           do j=jstop+2,platd
              jj = jj + 1
              do i = 1,plon
                 xxx(istart-1+i,j,k) = xx_exts(i,pext/2+jj,k)
              end do
           end do
        end if
      end do
!
! Fill E/W points
!
      call extx (nxpt    ,platd   ,plond   ,plev    , &
                 plon    ,istart  ,istop   ,jstart  ,jstop   , &
                 xxx     )
!
! Interpolation weights for x-direction
!
      do i = 1,plono
         call lcdbas(plond, lam, dlam, lamdp(i), idp(i), &
                     lbasdx(1,1,i), lbasdx(1,2,i) , &
                     hl(i), hr(i) , dhl(i), dhr(i) )
      end do
!
! Interpolation weights for y-direction
!
      do j = 1,plato
         call lcdbas(platd, phi, dphi, phidp(j), jdp(j), &
                     lbasdy(1,1,j), lbasdy(1,2,j), &
                     hs(j), hn(j) , dhs(j), dhn(j) )
      end do
!
! Hermite Cubic interpolation (using 4X4 input stencil for each output grid point)
!
      do k = 1,plev
        do jj = 1,plato

          rdphi = 1./dphi(jdp(jj))

          do ii = 1,plono
            rdx = 1./dlam(idp(ii))
!
! x-interpolation (over 4 adjacent latitudes in stencil)
!
            do n = 1,4
               j = jdp(jj)-2+n
               call lcdint(xxx(idp(ii)-1,j,k),lbasdx(1,1,ii), rdx, &
                           limdr, hl(ii), hr(ii), dhl(ii), dhr(ii), &
                           fintx(n))
            end do
!
! y-interpolation
!
            call lcdint(fintx, lbasdy(1,1,jj), rdphi, limdr, &
                        hs(jj), hn(jj), dhs(jj), dhn(jj), &
                        yy(ii,jj,k))
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      subroutine extx (nxpt    ,platd   ,plond   ,plev    , &
                       plon    ,istart  ,istop   ,jstart  ,jstop   , &
                       fb      )

!-----------------------------------------------------------------------
! 
! Purpose: 
! Copy data to the longitude extensions of the extended array
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------

      implicit none

!------------------------------Arguments--------------------------------
!
      integer nxpt          ! Extra point in local 4-pt interpolation
!                           ! stencil
      integer platd         ! Total latitude dimension of input grid
!                           ! including extensions
      integer plond         ! Total longitude dimension of input grid
!                           ! including extensions
!
      integer plev                 ! vertical dimension of input/output
!                                  ! field
      integer plon                 ! longitude dimension of input field
      integer istart               ! index for first analysis long.
      integer istop                ! index for last  analysis long.
      integer jstart               ! index for first analysis lat.
      integer jstop                ! index for last  analysis lat.
      real*8  fb(plond,platd,plev) ! input field on extended grid
!
!---------------------------Local variables-----------------------------
!
      integer i                 ! longitude index
      integer j                 ! latitude  index
      integer k                 ! vertical  index
      integer i2pi              ! start of eastern long. extension
!
!-----------------------------------------------------------------------
!
! Fill west edge points.
!
      if(nxpt .ge. 1) then
        do k=1,plev
          do j=1,platd
            do i=1,nxpt
              fb(i,j,k) = fb(i+plon,j,k)
            end do
          end do
        end do
      end if
!
! Fill east edge points
!
      i2pi = nxpt + plon + 1
      do k=1,plev
        do j=1,platd
          do i=i2pi,plond
            fb(i,j,k) = fb(i-plon,j,k)
          end do
        end do
      end do
!
      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine extys(plon    ,plev    ,plat    ,extent_dim, clat    , &
                       fb      ,fb_extents)

!-----------------------------------------------------------------------
! 
! Purpose: 
! Fill latitude extensions of a scalar extended array
! 
! Method: 
! This is done in 2 steps:
!   1) interpolate to the pole points; use the mean field value on the
!      Gaussian latitude closest to the pole.
!   2) add latitude lines beyond the poles.
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------

      implicit none

!------------------------------Arguments--------------------------------
      integer plon               ! longitude dimension of input field
      integer plev               ! vertical dimension of input/output
      integer plat               ! Total latitude dimension of input grid
      integer extent_dim         ! # of latitude extensions
      real*8  clat (plat )       ! Input latitudes in degrees
      real*8  fb(plon,plat,plev) ! input field on extended grid
      real*8  fb_extents(plon,extent_dim,plev) ! latitude extensions

!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local variables-----------------------------
!
      real*8  zave              ! accumulator for zonal averaging
      integer i,j,k             ! indices
      integer plon2             ! half the number of real longitudes
      integer platn, plats      ! indices
      logical npole             ! flag set if northern-most input latitude is a pole point
      logical spole             ! flag set if southern-most input latitude is a pole point
!
!-----------------------------------------------------------------------
!
! This code Hardwired to 6 latitude extensions
!
      if(extent_dim .ne. 6) then
         write(6,*) 'Error in EXTYS:  latitude extensions hardwired'
         write(6,*) 'to 6'
         call abort
      end if
!
! Check if input dataset has pole points
!
      npole  = .false.
      spole  = .false.
      if(clat(   1) .lt. -89.9999) spole = .true.
      if(clat(plat) .gt.  89.9999) npole = .true.

      plon2 = plon/2
!
! Fill north pole line.
!
      if (npole) then
         do k=1,plev
            do i=1,plon
               fb_extents(i,4,k) = fb(i,plat,k)
            end do
         end do
      else
         do k=1,plev
            zave = 0.0
            do i = 1,plon
               zave = zave + fb(i,plat,k)
            end do
            zave = zave/plon
            do i=1,plon
               fb_extents(i,4,k) = zave
            end do
         end do
      end if
!
! Fill northern lines beyond pole line.
!
      if (npole) then
         platn = plat-1
         plats = plat-2
      else
         platn = plat
         plats = plat-1
      end if

      do k=1,plev
         do i=1,plon2
            fb_extents(      i,5,k) = fb(plon2+i,platn,k)
            fb_extents(plon2+i,5,k) = fb(      i,platn,k)
            fb_extents(      i,6,k) = fb(plon2+i,plats,k)
            fb_extents(plon2+i,6,k) = fb(      i,plats,k)
         end do
      end do
!
! Fill south pole line.
!
      if (spole) then
         do k=1,plev
            do i=1,plon
               fb_extents(i,3,k) = fb(i,1,k)
            end do
         end do
      else
         do k=1,plev
            zave = 0.0
            do i = 1,plon
               zave = zave + fb(i,1,k)
            end do
            zave = zave/plon
            do i=1,plon
               fb_extents(i,3,k) = zave
            end do
         end do
      end if
!
! Fill southern lines beyond pole line.
!
      if (spole) then
         platn = 3
         plats = 2
      else
         platn = 2
         plats = 1
      end if

      do k=1,plev
         do i=1,plon2
            fb_extents(      i,1,k) = fb(plon2+i,platn,k)
            fb_extents(plon2+i,1,k) = fb(      i,platn,k)
            fb_extents(      i,2,k) = fb(plon2+i,plats,k)
            fb_extents(plon2+i,2,k) = fb(      i,plats,k)
         end do
      end do
!
      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine extyv(plon    ,plev    ,plat    ,extent_dim,clon   , &
                       clat    , &
                       fbu     ,fbv     ,fbu_extents, fbv_extents   )

!-----------------------------------------------------------------------
! 
! Purpose: 
! Fill latitude extensions of a vector quantity
! 
! Method: 
! This is done in 2 steps:
!   1) interpolate to the pole points; project the orthogonal wave 1
!      of U and V to the pole; use the Gaussian latitude closest to the pole.
!   2) add latitude lines beyond the poles.
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------

      implicit none

!------------------------------Arguments--------------------------------
      integer plon               ! longitude dimension of input field
      integer plev               ! vertical dimension of input/output
      integer plat               ! Total latitude dimension of input grid
      integer extent_dim         ! # of latitude extensions
      real*8  clon (plon )       ! Input longitude in degrees
      real*8  clat (plat )       ! Input latitude  in degrees
      real*8  fbu(plon,plat,plev) ! input field on extended grid
      real*8  fbv(plon,plat,plev) ! input field on extended grid
      real*8  fbu_extents(plon,extent_dim,plev) ! latitude extensions
      real*8  fbv_extents(plon,extent_dim,plev) ! latitude extensions

!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local variables-----------------------------
!
      real*8  zave              ! accumulator for zonal averaging
      integer i,j,k             ! indices
      integer plon2             ! half the number of real longitudes
      real*8  pi
      real*8  zavecv
      real*8  zavesv
      real*8  zavecu
      real*8  zavesu
      real*8  zavcus  
      real*8  zaucvs
      real*8  sinlam(plon)      ! Sin of lamda at all points in lat circle
      real*8  coslam(plon)      ! Cos of lamda at all points in lat circle
      logical npole             ! flag set if northern-most input latitude is a pole point
      logical spole             ! flag set if southern-most input latitude is a pole point
!
!-----------------------------------------------------------------------
!
! This code Hardwired to 6 latitude extensions
!
      if(extent_dim .ne. 6) then
         write(6,*) 'Error in EXTYV:  latitude extensions hardwired'
         write(6,*) 'to 6'
         call abort
      end if
!
! Check if input dataset has pole points
!
      npole  = .false.
      spole  = .false.
      if(clat(   1) .lt. -89.9999) spole = .true.
      if(clat(plat) .gt.  89.9999) npole = .true.

      plon2 = plon/2
      pi = 4.*atan(1.)
      do i = 1,plon
         sinlam(i) = sin( clon(i)*pi/180. )
         coslam(i) = cos( clon(i)*pi/180. )
      end do
!
! Fill north pole line.
!
      if (npole) then
         do k = 1,plev
            do i = 1,plon
               fbv_extents(i,4,k) = fbv(i,plat,k)
               fbu_extents(i,4,k) = fbu(i,plat,k)
            end do
         end do
      else
         do k = 1,plev
            zavecv = 0.0
            zavesv = 0.0
            zavecu = 0.0
            zavesu = 0.0
            do i = 1,plon
               zavecv = zavecv + fbv(i,plat,k  )*coslam(i)
               zavesv = zavesv + fbv(i,plat,k  )*sinlam(i)
               zavecu = zavecu + fbu(i,plat,k  )*coslam(i)
               zavesu = zavesu + fbu(i,plat,k  )*sinlam(i)
            end do
            zavcus = (zavecv + zavesu)/plon
            zaucvs = (zavecu - zavesv)/plon
            do i = 1,plon
               fbv_extents(i,4,k) = zavcus*coslam(i) - zaucvs*sinlam(i)
               fbu_extents(i,4,k) = zaucvs*coslam(i) + zavcus*sinlam(i)
            end do
         end do
      end if
!
! Fill northern lines beyond pole line.
!
      do k=1,plev
         do i=1,plon2
            fbu_extents(      i,5,k) = fbu(plon2+i,plat  ,k)
            fbu_extents(plon2+i,5,k) = fbu(      i,plat  ,k)
            fbv_extents(      i,5,k) = fbv(plon2+i,plat  ,k)
            fbv_extents(plon2+i,5,k) = fbv(      i,plat  ,k)
            fbu_extents(      i,6,k) = fbu(plon2+i,plat-1,k)
            fbu_extents(plon2+i,6,k) = fbu(      i,plat-1,k)
            fbv_extents(      i,6,k) = fbv(plon2+i,plat-1,k)
            fbv_extents(plon2+i,6,k) = fbv(      i,plat-1,k)
         end do
      end do
!
! Fill south pole line.
!
      if (spole) then
         do k = 1,plev
            do i = 1,plon
               fbv_extents(i,3,k) = fbv(i,1,k)
               fbu_extents(i,3,k) = fbu(i,1,k)
            end do
         end do
      else
         do k = 1,plev
            zavecv = 0.0
            zavesv = 0.0
            zavecu = 0.0
            zavesu = 0.0
            do i = 1,plon
               zavecv = zavecv + fbv(i,1,k  )*coslam(i)
               zavesv = zavesv + fbv(i,1,k  )*sinlam(i)
               zavecu = zavecu + fbu(i,1,k  )*coslam(i)
               zavesu = zavesu + fbu(i,1,k  )*sinlam(i)
            end do
            zavcus = (zavecv - zavesu)/plon
            zaucvs = (zavecu + zavesv)/plon
            do i = 1,plon
               fbv_extents(i,3,k) = zavcus*coslam(i) + zaucvs*sinlam(i)
               fbu_extents(i,3,k) = zaucvs*coslam(i) - zavcus*sinlam(i)
            end do
         end do
      end if
!
! Fill southern lines beyond pole line.
!
      do k=1,plev
         do i=1,plon2
            fbu_extents(      i,2,k) = fbu(plon2+i,1,k)
            fbu_extents(plon2+i,2,k) = fbu(      i,1,k)
            fbv_extents(      i,2,k) = fbv(plon2+i,1,k)
            fbv_extents(plon2+i,2,k) = fbv(      i,1,k)
            fbu_extents(      i,1,k) = fbu(plon2+i,2,k)
            fbu_extents(plon2+i,1,k) = fbu(      i,2,k)
            fbv_extents(      i,1,k) = fbv(plon2+i,2,k)
            fbv_extents(plon2+i,1,k) = fbv(      i,2,k)
         end do
      end do
!
      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      subroutine lcdbas(idim, grd, dgrd, dp, idp, dbas2, dbas3, &
                        hb  , ht , dhb , dht    )
!
!-----------------------------------------------------------------------
! 
! Purpose: 
! Evaluate cubic basis functions for a 4-point stencil
! 
! Method: 
! 
! Author: J. Olson
!
!-----------------------------------------------------------------------

      implicit none

!------------------------------Arguments--------------------------------
!
      integer idim       ! dimension of input grid
      real*8  grd (idim) ! input grid
      real*8  dgrd(idim) ! input grid intervals
      real*8  dp         ! target grid point
      integer idp        ! pointer into grid for target grid pt
      real*8  dbas2(4)   ! derivatives at grid point 2.
      real*8  dbas3(4)   ! derivatives at grid point 3.
      real*8  hb         ! Interpolation weight for bot field value
      real*8  ht         ! Interpolation weight for top field value
      real*8  dhb        ! Interpolation weight for bot derivative
      real*8  dht        ! Interpolation weight for top derivative
!
!---------------------------Local variables-----------------------------
!
      integer i
      real*8 dx
      real*8 xb
      real*8 xt
      real*8 x1                   !  |
      real*8 x2                   !  |- grid values
      real*8 x3                   !  |
      real*8 x4                   !  |
      real*8 x1mx2                !  |
      real*8 x1mx3                !  |
      real*8 x1mx4                !  |- differences of grid values
      real*8 x2mx3                !  |
      real*8 x2mx4                !  |
      real*8 x3mx4                !  |
!
!-----------------------------------------------------------------------
!
      idp = -1
      do i = 2,idim-2
        if(dp .ge. grd(i) .and. dp .lt. grd(i+1) ) then
          idp = i
        end if
      end do
      if(idp .lt. 0) then
        write(6,*) 'Error in LCDBAS:  target grid point outside of'
        write(6,*) 'input grid'
        call abort
      end if
!
      x1 = grd(idp-1)
      x2 = grd(idp  )
      x3 = grd(idp+1)
      x4 = grd(idp+2)

      x1mx2 = x1 - x2
      x1mx3 = x1 - x3
      x1mx4 = x1 - x4
      x2mx3 = x2 - x3
      x2mx4 = x2 - x4
      x3mx4 = x3 - x4

      dbas2(1) =   x2mx3 * x2mx4 / ( x1mx2 * x1mx3 * x1mx4 )
      dbas2(2) =   -1./x1mx2 + 1./x2mx3 + 1./x2mx4
      dbas2(3) = - x1mx2 * x2mx4 / ( x1mx3 * x2mx3 * x3mx4 )
      dbas2(4) =   x1mx2 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )

      dbas3(1) = - x2mx3 * x3mx4 / ( x1mx2 * x1mx3 * x1mx4 )
      dbas3(2) =   x1mx3 * x3mx4 / ( x1mx2 * x2mx3 * x2mx4 )
      dbas3(3) =   -1./x1mx3 - 1./x2mx3 + 1./x3mx4
      dbas3(4) = - x1mx3 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )

      dx       = dgrd(idp)
      xb       = ( grd(idp+1) - dp )/dx
      xt       = 1. - xb
      hb       = ( 3.0 - 2.0*xb )*xb**2
      ht       = ( 3.0 - 2.0*xt )*xt**2
      dhb      = -dx*( xb - 1.  )*xb**2
      dht      =  dx*( xt - 1.  )*xt**2

      return
!
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      subroutine lcdint(x, lbasdx, rdx, limdr, hb, ht, dhb, dht, &
                        y        )
!
!-----------------------------------------------------------------------
! 
! Purpose: 
! Interpolate data on 4-pt stencil using cubic basis functions
! 
! Method: 
! 
! Author: J. Olson
!
!-----------------------------------------------------------------------

      implicit none

!------------------------------Arguments--------------------------------
!
      real*8  x(4)        ! input field on 4-pt stencil
      real*8  lbasdx(4,2) ! basis functions
      real*8  rdx         ! grid interval
      logical limdr       ! Flag to invoke SCM0 derivative limiting
      real*8  hb          ! Interpolation weight for bot field value
      real*8  ht          ! Interpolation weight for top field value
      real*8  dhb         ! Interpolation weight for bot derivative
      real*8  dht         ! Interpolation weight for top derivative
      real*8  y           ! Interpolated value
!
!---------------------------Local variables-----------------------------
!
      real*8 fbot, ftop
      real*8 deli, tmp1, tmp2
      real*8 fac
!
!-----------------------------------------------------------------------
!
      fbot = lbasdx(1,1)*x(1) &
           + lbasdx(2,1)*x(2) &
           + lbasdx(3,1)*x(3) &
           + lbasdx(4,1)*x(4)
      ftop = lbasdx(1,2)*x(1) &
           + lbasdx(2,2)*x(2) &
           + lbasdx(3,2)*x(3) &
           + lbasdx(4,2)*x(4)
!
! Apply SCM0 limiter to derivative estimates.
!
      if(limdr) then
         fac  = 3.*(1. - 10.*epsilon(fac))
         deli = ( x(3) - x(2) )*rdx
         tmp1 = fac*deli
         tmp2 = abs( tmp1 )
         if( deli*fbot   .le. 0.0  ) fbot = 0.
         if( deli*ftop   .le. 0.0  ) ftop = 0.
         if( abs( fbot ) .gt. tmp2 ) fbot = tmp1
         if( abs( ftop ) .gt. tmp2 ) ftop = tmp1
      end if

      y = x(2)*hb + fbot*dhb + x(3)*ht + ftop*dht

      return
!
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine tsadj(plat    ,plon    ,phis_old,phis_new,ts      )
!
!-----------------------------------------------------------------------
!
! Adjust Ts based on difference between old and new phis.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     INPUTS
!
      integer plat !  latitude dimension
      integer plon !  longitude dimension
!
      real*8 phis_old(plon,plat) ! analysis phis (e.g., ECMWF)
      real*8 phis_new(plon,plat) ! model phis
!
!     INPUT/OUTPUT
!
      real*8 ts      (plon,plat) ! Surface Temp
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
!
      real*8 dtdz
      real*8 gravit
      real*8 del_z
      integer i, j, k             ! Indices
!
      dtdz    = -0.0065           ! -6.5 deg/km
      gravit  = 9.80616           ! acceleration of gravity ~ m/s^2

      do j = 1,plat
        do i = 1,plon

          del_z = ( phis_new(i,j) - phis_old(i,j) )/gravit
          ts(i,j) = ts(i,j) + dtdz*del_z

        end do
      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine vert_quad_opt1(plevi , plevip1 , plev , plat  ,plon, &
                       t_old ,pressi_m ,pressi_i ,presso_m, phis_old, &
                       ps_old  , t_new , loglin )
!
!-----------------------------------------------------------------------
!
! Quadratic interpolation (designed for Temperature interpolation)
!
!                         (if "loglin" == 1: in P
!                          if "loglin" /= 1: in ln(P) )
!
!  Above input top           :  quadratic using top, levels 1 and 2 (top
!                               defined as 1.e-10 Pa for now).  Top
!                               value set to value at level 1
!  Between levels 1 and "bot":  quadratic interp using 3 closest levels
!  Between levels "bot"
!                 and surface:  linear interpolation using "Tbot"and
!                               "Tsurf".
!  Below surface             :  You don"t wanna know
!
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     INPUTS
!
      integer plevi     ! vertical dimension of input analysis fields
      integer plevip1   ! "plevi+1 (vert dimension of input interfaces)
      integer plev      ! vertical dimension of model fields
      integer plat      ! latitude dimension
      integer plon      ! longitude dimension
!
      real*8 t_old   (plevi  ,plon,plat)   ! analysis tempertatures
      real*8 pressi_m (plevi  ,plon,plat)  ! analysis pressures at all levels
      real*8 pressi_i (plevip1,plon,plat)  ! analysis interface pressures
      real*8 presso_m (plev   ,plon,plat)  ! model pressures (based on adjusted PS)
      real*8 phis_old(plon,     plat)      ! analysis phis
      real*8 ps_old  (plon,     plat)      ! analysis surface pressure
      real*8 ps_new  (plon,     plat)      ! "adjusted" model surface pressure
      integer loglin                       ! interpolation flag
!
!     OUTPUT
!
      real*8 t_new   (plev   ,plon,plat)   ! Interpolated Temperatures
!
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
!
      real*8 tsurf
      real*8 t0
      real*8 t_ref1
      real*8 t_ref2
      real*8 t_ref3
      real*8 t_ref3_top
      real*8 t_ref3_bot
      real*8 z_ref_top
      real*8 z_ref_bot
      real*8 tbot
      real*8 pbot
      real*8 psurf
      real*8 dtdz
      real*8 lapse
      real*8 boltz
      real*8 avogad
      real*8 mwdair
      real*8 rgas
      real*8 rdair
      real*8 gravit
      real*8 x
      real*8 threshold
      real*8 tmp
      real*8 z
      real*8 z_min
      real*8 z_incr
      real*8 hkk
      real*8 x1
      real*8 x2
      real*8 x3
      real*8 p1
      real*8 p2
      real*8 p3
      real*8 px
      real*8 pt                       ! top input pressure (ghost)
      real*8 xt                       ! top input value (linear extrap)
      real*8 tmp1
      real*8 tmp2
      real*8 tmp3
      real*8 tmpt
      real*8 beta
      real*8 zero
      integer i, j, k, kk, kkp1, kkp2 ! Indices
      integer k_bot
!
      zero    = 0.
      dtdz    = -0.0065           ! -6.5 deg/km
      gravit  = 9.80616           ! acceleration of gravity ~ m/s^2
      boltz   = 1.38065e-23       ! boltzmann's constant ~ J/k/molecule
      avogad  = 6.02214e26        ! avogadro's number ~ molecules/kmole
      mwdair  = 28.966            ! molecular weight dry air ~ kg/kmole

      rgas    = avogad*boltz      ! universal gas constant ~ J/k/kmole
      rdair   = rgas/mwdair       ! constant for dry air   ~ J/k/kg

      t_ref1    = 290.5
      t_ref2    = 255.0
      t_ref3    = 298.0
      z_ref_bot = 2000.
      z_ref_top = 2500.

      threshold = 0.001
      pt = 1.e-10
      if(loglin .ne. 1) pt = log(pt)
      
      do j = 1,plat
        do i = 1,plon
!
! Tbot and Pbot are determined from the first model level that is at
! least 150m above the surface
!  ... will be used later
!
            z_min = 150.
            z     = 0.

            do k = plevi,1,-1
              k_bot  = k
              hkk    = 0.5*( pressi_i(k+1,i,j) - pressi_i(k,i,j) )/ &
                                                         pressi_m(k,i,j)
              z_incr = (rdair/gravit)*t_old(k,i,j)*hkk
              z      = z + z_incr
              if(z .gt. z_min) go to 10
              z      = z + z_incr
            end do

            write(6,*) 'Error:  could not find model level above ',z_min
            call abort

   10       continue
            lapse = -dtdz

            tbot  = t_old  (k_bot,i,j)
            pbot  = pressi_m(k_bot,i,j)
            tmp   = lapse*(rdair/gravit)*(ps_old(i,j)/pbot - 1.)
            tsurf = tbot*(1. + tmp)

          kk   = 1
          kkp1 = kk + 1
          kkp2 = kk + 2
!
! Find bracketting input levels if possible
!
          do k = 1,plev
   20       continue
            if(pressi_m(kkp1,i,j) .le. presso_m(k,i,j) ) then
              beta = presso_m(k,i,j) - pressi_m(kkp1,i,j)
              beta = beta/( pressi_m(kkp2,i,j) - pressi_m(kkp1,i,j) )
              if(beta .ge. 0.5) then
                kk   = kk + 1
                kkp1 = kk + 1
                kkp2 = kk + 2
                if(kkp2 .ge. plevi+1) then
                  kk   = plevi - 2
                  kkp1 = kk + 1
                  kkp2 = kk + 2
                  go to 30
                endif
                go to 20
              endif
            endif

   30       continue
            
            x1    = t_old   (kk   ,i,j)
            x2    = t_old   (kkp1 ,i,j)
            x3    = t_old   (kkp2 ,i,j)
            p1    = pressi_m(kk   ,i,j)
            p2    = pressi_m(kkp1 ,i,j)
            p3    = pressi_m(kkp2 ,i,j)
            pbot  = pressi_m(k_bot,i,j)
            psurf = ps_old  (      i,j)

            px    = presso_m(k    ,i,j)
!
! Convert to log(P) for log(P) interpolation
!
            if(loglin .ne. 1) then
              p1    = log(p1)
              p2    = log(p2)
              p3    = log(p3)
              pbot  = log(pbot)
              psurf = log(psurf)
              px    = log(px)
            endif
!
! If above 1st analysis level: quadratic interp
!
            if    (px .lt. p1 ) then
              xt         = x1

              tmpt       = ( (px-p1)*(px-p2) )/( (pt-p1)*(pt-p2) )
              tmp1       = ( (px-pt)*(px-p2) )/( (p1-pt)*(p1-p2) )
              tmp2       = ( (px-pt)*(px-p1) )/( (p2-pt)*(p2-p1) )

              t_new(k,i,j) = xt*tmpt + x1*tmp1 + x2*tmp2
!PFC-diag+++
!!           if(t_new(k,i,j).lt.161.) t_new(k,i,j)=161.
              t_new(k,i,j) = x1
!PFC-diag===

!
! Elseif between "pbot" and analysis surface: linear interp
!
            elseif(px .ge. pbot .and. px .le. psurf ) then
              t_new(k,i,j) = ( tbot*(psurf - px) + tsurf*(px - pbot) ) &
                                                         /(psurf - pbot)
!
! Elseif below analysis surface: special case - heuristic equations
!
            elseif(px .gt. psurf ) then
              t_ref3_bot =      tsurf - z_ref_bot*dtdz
              t_ref3_top = min( tsurf - z_ref_top*dtdz, t_ref3 )
              z = phis_old(i,j)/gravit

              if(z .ge. z_ref_bot) then
                if(z .ge. z_ref_top) then
                  t0 = min( tsurf - z*dtdz, t_ref3 )
                else
                  t0 = ( t_ref3_bot*(z_ref_top - z        )   + &
                         t_ref3_top*(z         - z_ref_bot) ) / &
                                    (z_ref_top - z_ref_bot)
                endif
                lapse = max ( (t0 - tsurf)/z , zero)
              else
                lapse = -dtdz
              endif

              x   = lapse*rdair/gravit*log( presso_m(k    ,i,j)/ &
                                            ps_old  (      i,j) )
              t_new(k,i,j) = tsurf*(1. + x + x*x/2. + x*x*x/6.)
!
! Should never happen
!
            elseif(px .ge. p3 ) then
              write(6,*) 'Error:  extrapolation below input levels'
              write(6,*) 'not allowed'
              call abort
!
! Else between 1st analysis level and "pbot":  quadratic interp
!
            else
              tmp1       = ( (px-p2)*(px-p3) )/( (p1-p2)*(p1-p3) )
              tmp2       = ( (px-p1)*(px-p3) )/( (p2-p1)*(p2-p3) )
              tmp3       = ( (px-p1)*(px-p2) )/( (p3-p1)*(p3-p2) )

              t_new(k,i,j) = x1*tmp1 + x2*tmp2 + x3*tmp3

            endif
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine vert_int_opt2(plat    ,plon    ,plevi   ,plevo   , &
                               pressi  ,presso  ,xxi     ,xxo     , &
                               loglin  )
!
!-----------------------------------------------------------------------
!
! Designed for vertical interpolation of U/V
!
! Linear and Quadratic interpolation (if "loglin" == 1: in P
!                                     if "loglin" /= 1: in ln(P) )
!
!  Above input top        :  quadratic using top, levels 1 and 2 (top
!                            defined as 1.e-10 Pa for now).  Top value
!                            determined from linear extrapolation from
!                            levels 1 and 2
!  Between levels 1 and 2 :  quadratic interp using levels 1,2, & 3
!  Between levels 2 and K :  linear interpolation using adjacent levels
!  Below level K          :  set equal to level K
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     INPUTS
!
      integer plat  ! latitude dimension
      integer plon  ! longitude dimension
      integer plevi ! vertical dimension of analysis fields
      integer plevo ! vertical dimension of model fields
!
      real*8 pressi  (plevi,plon,plat) ! analysis pressures
      real*8 presso  (plevo,plon,plat) ! model pressures (based on adjusted PS)
      real*8 xxi     (plevi,plon,plat) ! input analysis field
      integer loglin                   ! interpolation flag
!
!     OUTPUTS
!
      real*8 xxo     (plevo,plon,plat) ! model field
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
!
      real*8 x1
      real*8 x2
      real*8 x3
      real*8 p1
      real*8 p2
      real*8 p3
      real*8 px
      real*8 pt                       ! top input pressure (ghost)
      real*8 xt                       ! top input value (linear extrap)
      real*8 tmp1
      real*8 tmp2
      real*8 tmp3
      real*8 tmpt
      integer i, j, k, kk, kkp1, kkp2 ! Indices
!
!-----------------------------------------------------------------------
!
      pt = 1.e-10
      if(loglin .ne. 1) pt = log(pt)
!
      do j = 1,plat
        do i = 1,plon

          kk   = 1
          kkp1 = kk + 1
          kkp2 = kk + 2
!
! Find bracketting analysis pressure levels
!
          do k = 1,plevo
   10       continue
            if(pressi(kkp1,i,j) .le. presso(k,i,j) ) then
              kk   = kk + 1
              kkp1 = kk + 1
              kkp2 = kk + 2
              if(kkp1 .eq. plevi+1) then
                kk   = plevi - 1
                kkp1 = kk + 1
                kkp2 = kk + 2
                go to 20
              endif
              go to 10
            endif

   20       continue
            
            x1 = xxi   (kk  ,i,j)
            x2 = xxi   (kkp1,i,j)
            p1 = pressi(kk  ,i,j)
            p2 = pressi(kkp1,i,j)
            px = presso(k   ,i,j)

            if(kkp2 .le. plevi) then
              x3 = xxi   (kkp2,i,j)
              p3 = pressi(kkp2,i,j)
            else
              x3 = 1.e+36
              p3 = 1.e+36
            endif
!
! Convert to log(P) for log(P) interpolation
!
            if(loglin .ne. 1) then
              p1 = log(p1)
              p2 = log(p2)
              p3 = log(p3)
              px = log(px)
            endif
!
! If above 1st analysis level:  quadratic interp
!
            if    (px .lt. p1 ) then
              xt         = ( x1*(p2 - pt) - x2*(p1 - pt) )/(p2 - p1)

              tmpt       = ( (px-p1)*(px-p2) )/( (pt-p1)*(pt-p2) )
              tmp1       = ( (px-pt)*(px-p2) )/( (p1-pt)*(p1-p2) )
              tmp2       = ( (px-pt)*(px-p1) )/( (p2-pt)*(p2-p1) )

              xxo(k,i,j) = xt*tmpt + x1*tmp1 + x2*tmp2

!PFC-daig+++
              xxo(k,i,j) = x1
!PFC-daig===

!
! Elseif below bottom analysis level:  output = bottome analysis field value
!
            elseif(px .ge. p2 ) then
              xxo(k,i,j) = x2
!
! Elseif between 1st and 2nd analysis levels:  quadratic interp
!
            elseif(kk .eq. 1 ) then
              tmp1       = ( (px-p2)*(px-p3) )/( (p1-p2)*(p1-p3) )
              tmp2       = ( (px-p1)*(px-p3) )/( (p2-p1)*(p2-p3) )
              tmp3       = ( (px-p1)*(px-p2) )/( (p3-p1)*(p3-p2) )

              xxo(k,i,j) = x1*tmp1 + x2*tmp2 + x3*tmp3
!
! Else, Linear interpolation
!
            else
              xxo(k,i,j) = ( x1*(p2 - px) + x2*(px - p1) )/(p2 - p1)

            endif
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine vert_int_opt1(plat    ,plon    ,plevi   ,plevo   , &
                               pressi  ,presso  ,xxi     ,xxo     , &
                               loglin  )
!
!-----------------------------------------------------------------------
!
! Designed for moisture fields like q, cloud water, cloud ice, cloud frac, etc.
!
! Linearly interpolate (if "loglin" == 1: in P
!                       if "loglin" /= 1: in ln(P) )
!
!  Above input top        :  set equal to level 1
!  Between levels 1 and K :  linear interpolation using adjacent levels
!  Below level K          :  set equal to level K
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     INPUTS
!
      integer plat   ! latitude dimension
      integer plon   ! longitude dimension
      integer plevi  ! vertical dimension of analysis fields
      integer plevo  ! vertical dimension of model fields
!
      real*8 pressi  (plevi,plon,plat)  ! analysis pressures
      real*8 presso  (plevo,plon,plat)  ! model pressures (based on adjusted PS)
      real*8 xxi     (plevi,plon,plat)  ! analysis field
      integer loglin                    ! interpolation flag
!
!     OUTPUTS
!
      real*8 xxo     (plevo,plon,plat)  ! model field
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
!
      real*8 p1
      real*8 p2
      real*8 px
      integer i, j, k, kk, kkp1         ! Indices
!
!-----------------------------------------------------------------------
!
      do j = 1,plat
        do i = 1,plon
!
! Find bracketting analysis pressure levels
!
          kk   = 1
          kkp1 = kk + 1

          do k = 1,plevo
   10       continue
            if(pressi(kkp1,i,j) .le. presso(k,i,j) ) then
              kk   = kk + 1
              kkp1 = kk + 1
              if(kkp1 .eq. plevi+1) then
                kk   = plevi - 1
                kkp1 = kk + 1
                go to 20
              endif
              go to 10
            endif

   20       continue
!
! If above 1st analysis level:  output = top analysis field value
!
            if    (presso(k,i,j) .lt. pressi(kk,i,j) ) then
              xxo(k,i,j) = xxi(kk  ,i,j)
!
! If below bottom analysis level:  output = bottom analysis field value
!
            elseif(presso(k,i,j) .ge. pressi(kkp1,i,j) ) then
              xxo(k,i,j) = xxi(kkp1,i,j)
!
! Else, Linear interpolation
!
            else
              p1 = pressi(kk  ,i,j)
              p2 = pressi(kkp1,i,j)
              px = presso(k   ,i,j)
              if(loglin .ne. 1) then
                p1 = log(p1)
                p2 = log(p2)
                px = log(px)
              endif
              xxo(k,i,j) = xxi(kk  ,i,j)*(p2 - px) +  &
                           xxi(kkp1,i,j)*(px - p1)
              xxo(k,i,j) = xxo(k   ,i,j)/(p2 - p1)
            endif
          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine myminmax(plev    ,plat    ,plon    ,x       ,fmin   , &
                          fmax    )
!
!-----------------------------------------------------------------------
!
! Bracket "x" between fmin and fmax
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer plev
      integer plat
      integer plon
!
      real*8 x       (plon,plat,plev)
      real*8 fmin
      real*8 fmax
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
!
      integer i, j, k           ! Indices

      do k = 1,plev
        do j = 1,plat
          do i = 1,plon

            x(i,j,k) = min(max(x(i,j,k),fmin),fmax)

          end do
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine psadj(plev    ,plevp1  ,plat    ,plon    ,t       , &
                       press_m ,press_i ,phis_old,phis_new,ps_old  , &
                       ps_new  )
!
!-----------------------------------------------------------------------
!
! Adjust Ps based on difference between "analysis" phis and model phis.
! Also uses T and P arrays
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     INPUTS
!
      integer plev      ! vertical dimension
      integer plevp1    ! "plev+1"
      integer plat      ! latitude dimension
      integer plon      ! longitude dimension
!
      real*8 t       (plon,plat,plev)    ! analysis Temperatures
      real*8 press_m (plon,plat,plev)    ! analysis pressures
      real*8 press_i (plon,plat,plevp1)  ! analysis pressures (interfaces)
      real*8 phis_old(plon,     plat)    ! analysis phis
      real*8 phis_new(plon,     plat)    ! model phis
      real*8 ps_old  (plon,     plat)    ! analysis Ps (horizontally
!                                        ! interpolated to model grid)
!
!     OUTPUTS
!
      real*8 ps_new  (plon,     plat)    ! adjusted model Ps
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
!
      real*8 tsurf
      real*8 t_ref1
      real*8 t_ref2
      real*8 t0
      real*8 tbot
      real*8 pbot
      real*8 dtdz
      real*8 lapse
      real*8 boltz
      real*8 avogad
      real*8 mwdair
      real*8 rgas
      real*8 rdair
      real*8 gravit
      real*8 del_phis
      real*8 x
      real*8 threshold
      real*8 tmp
      real*8 z
      real*8 z_min
      real*8 z_incr
      real*8 hkk
      integer i, j, k, kk         ! Indices
!
      dtdz    = -0.0065           ! -6.5 deg/km
      gravit  = 9.80616           ! acceleration of gravity ~ m/s^2
      boltz   = 1.38065e-23       ! boltzmann's constant ~ J/k/molecule
      avogad  = 6.02214e26        ! avogadro's number ~ molecules/kmole
      mwdair  = 28.966            ! molecular weight dry air ~ kg/kmole

      rgas    = avogad*boltz      ! universal gas constant ~ J/k/kmole
      rdair   = rgas/mwdair       ! constant for dry air   ~ J/k/kg

      t_ref1    = 290.5
      t_ref2    = 255.0
      threshold = 0.001
      
      do j = 1,plat
        do i = 1,plon

          del_phis = phis_old(i,j) - phis_new(i,j)
!
! If difference between analysis and model phis is negligible,
! then set model Ps = analysis
!
          if(abs(del_phis) .le. threshold) then
            ps_new(i,j) = ps_old(i,j)
!
! Else, go nuts...
!
          else
!
! Tbot and Pbot are determined from the first model level that is at
! least 150m above the surface
!
            z_min = 150.
            z     = 0.

            do k = plev,1,-1
              kk     = k
              hkk    = 0.5*( press_i(i,j,k+1) - press_i(i,j,k) )/ &
                                                          press_m(i,j,k)
              z_incr = (rdair/gravit)*t(i,j,k)*hkk
              z      = z + z_incr
              if(z .gt. z_min) go to 10
              z      = z + z_incr
            end do

            write(6,*) 'Error:  could not find model level above ',z_min
            call abort

   10       continue
            lapse = -dtdz
            k     = kk
!
! Define Tbot & Pbot
!
            tbot  = t      (i,j,k)
            pbot  = press_m(i,j,k)
            tmp   = lapse*(rdair/gravit)*(ps_old(i,j)/pbot - 1.)
            tsurf = tbot*(1. + tmp)
            t0    = tsurf + lapse*phis_old(i,j)/gravit
!
! Based on heuristic equations:
!
            if     (t0 .gt. t_ref1 .and. tsurf .le. t_ref1) then
              lapse = (t_ref1 - tsurf)*gravit/phis_old(i,j)
            elseif (t0 .gt. t_ref1 .and. tsurf .gt. t_ref1) then
              lapse = 0.
              tsurf = (t_ref1 + tsurf)*0.5
            endif

            if(tsurf .lt. t_ref2) then
              lapse = -dtdz
              tsurf = (t_ref2 + tsurf)*0.5
            end if              

            x   = lapse*del_phis/(gravit*tsurf)
            tmp = 1. - x/2. + x**2./3.
            tmp = del_phis/(rdair*tsurf)*tmp
            ps_new(i,j) = ps_old(i,j)*exp(tmp)

          endif
        end do
      end do

      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine q2rh(plat    ,plev    ,plon    ,q       ,t       , &
                      press   )
!
!-----------------------------------------------------------------------
!
! Compute RH from T, Pressure, and Specific Humidity
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer plat
      integer plev
      integer plon
!
      real*8 q    (plon,plat,plev)
      real*8 t    (plon,plat,plev)
      real*8 press(plon,plat,plev)
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------- Commons ----------------------------------
!
! Common block and statement functions for saturation vapor pressure
! look-up procedure, J. J. Hack, February 1990
!
      integer plenest  ! length of saturation vapor pressure table
      parameter (plenest=250)
!
! Table of saturation vapor pressure values es from tmin degrees
! to tmax+1 degrees k in one degree increments.  ttrice defines the
! transition region where es is a combination of ice & water values
!
      common/comes/estbl(plenest) ,tmin  ,tmax  ,ttrice ,pcf(6) , &
                   icephs
!
      real*8 estbl      ! table values of saturation vapor pressure
      real*8 tmin       ! min temperature (K) for table
      real*8 tmax       ! max temperature (K) for table
      real*8 ttrice     ! transition range from es over H2O to over ice
      real*8 pcf        ! polynomial coeffs -> es transition h2o to ice
      logical icephs  ! false => saturation vapor press over water only
!
! Dummy variables for statement functions
!
      real*8 td         ! dummy variable for function evaluation
      real*8 tlim       ! intermediate var for es look-up with estbl4
      real*8 estblf     ! statement function es look-up
      real*8 estbl4     ! statement function es look-up
!
!---------------------------Local workspace-----------------------------
!
      real*8 es
      real*8 qs
      real*8 epsilo 
      real*8 latvap 
      real*8 latice 
      real*8 rh2o   
      real*8 cpair  
      real*8 omeps              ! 1 - 0.622
      real*8 qmin,qmax
      real*8 one
      integer countneg,countgt1
      integer i, j, k           ! Indices
      integer imin, jmin, kmin
      integer imax, jmax, kmax
      integer itype
!
!-------------------------- Statement Functions ------------------------
!
! Statement functions used in saturation vapor pressure table lookup
! there are two ways to use these three statement functions.
! For compilers that do a simple in-line expansion:
! => ttemp = tlim(t)
!    es    = estbl4(ttemp)
!
! For compilers that provide real optimization:
! => es    = estblf(t)
!
      tlim(td) = max(min(td,tmax),tmin)
!
      estblf(td) =  (tmin + int(tlim(td)-tmin) - tlim(td) + 1.0) &
                  *estbl(int(tlim(td)-tmin)+1) &
                  -(tmin + int(tlim(td)-tmin) - tlim(td)      ) &
                  *estbl(int(tlim(td)-tmin)+2)
!
      estbl4(td) =  (tmin+int(td-tmin)+1.0-td)*estbl(int(td-tmin)+1) &
                  + ( td-(tmin+int(td-tmin)) )*estbl(int(td-tmin)+2)
!
!-----------------------------------------------------------------------
!
      one    = 1.0
      epsilo = 0.622
      latvap = 2.5104e06
      latice = 3.336e5
      rh2o   = 4.61e2
      cpair  = 1004.64
      omeps  = 1.0 - epsilo
!
! Build es table
!
      call esinti (epsilo ,latvap ,latice ,rh2o ,cpair )
!
      qmin  =  1.e+36
      qmax  = -1.e+36
      countneg = 0
      countgt1 = 0
!
      do k = 1,plev
        do j = 1,plat
          do i = 1,plon
!
! Saturation specific humidity
!
            es = estblf( t(i,j,k) )
            qs = epsilo*es/(press(i,j,k) - omeps*es)
!   
! The following check is to avoid the generation of negative values
! that can occur in the upper stratosphere and mesosphere
!
            qs = min(one,qs)
!         
            if (qs .lt. 0.0) then
              qs = 1.0
            end if
!
! Compute RH
!
            q(i,j,k) = q(i,j,k)/qs
!
            if (q(i,j,k) .lt. qmin) then
              qmin = q(i,j,k)
              imin = i
              kmin = k
              jmin = j
            endif
            if (q(i,j,k) .gt. qmax) then
              qmax = q(i,j,k)
              imax = i
              kmax = k
              jmax = j
            endif
!
            if (q(i,j,k) .lt. 0.0) then
              countneg = countneg + 1
              q(i,j,k) = 0.
            endif
!
            if (q(i,j,k) .gt. 1.) then
              countgt1 = countgt1 + 1
              q(i,j,k) = 1.
            endif
!
          end do
        end do
      end do
!
      write(6,*) ' '
      write(6,*) ' '
      write(6,1000) qmax*100.,imax,jmax,kmax
      write(6,2000) qmin*100.,imin,jmin,kmin
      if(countneg .gt. 0) write(6,3000) countneg, plon*plev*plat, &
                                  float(countneg)/(plon*plev*plat)*100.
      if(countgt1 .gt. 0) write(6,4000) countgt1, plon*plev*plat,  &
                                  float(countgt1)/(plon*plev*plat)*100.
      write(6,*) ' '
      write(6,*) ' '
!
 1000 format(' Maximum RH = ',f12.5,'% at i,j,k = ',3i5)
 2000 format(' Minimum RH = ',f12.5,'% at i,j,k = ',3i5)
 3000 format(' WARNING in Q2RH:  ',i10,' points (out of ',i10, &
             ' (',f12.5,'%) ) were negative.  All were set to 0')
 4000 format(' WARNING in Q2RH:  ',i10,' points (out of ',i10, &
             ' (',f12.5,'%) ) were  .gt. 1.   All were set to 1')
!
      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine rh2q(plat    ,plev    ,plon    ,q       ,t       , &
                      press   )
!
!-----------------------------------------------------------------------
!
! Compute Q from T, Pressure, and RH
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer plat
      integer plev
      integer plon
!
      real*8 q    (plon,plat,plev)
      real*8 t    (plon,plat,plev)
      real*8 press(plon,plat,plev)
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------- Commons ----------------------------------
!
! Common block and statement functions for saturation vapor pressure
! look-up procedure, J. J. Hack, February 1990
!
      integer plenest  ! length of saturation vapor pressure table
      parameter (plenest=250)
!
! Table of saturation vapor pressure values es from tmin degrees
! to tmax+1 degrees k in one degree increments.  ttrice defines the
! transition region where es is a combination of ice & water values
!
      common/comes/estbl(plenest) ,tmin  ,tmax  ,ttrice ,pcf(6) , &
                   icephs
!
      real*8 estbl      ! table values of saturation vapor pressure
      real*8 tmin       ! min temperature (K) for table
      real*8 tmax       ! max temperature (K) for table
      real*8 ttrice     ! transition range from es over H2O to over ice
      real*8 pcf        ! polynomial coeffs -> es transition h2o to ice
      logical icephs  ! false => saturation vapor press over water only
!
! Dummy variables for statement functions
!
      real*8 td         ! dummy variable for function evaluation
      real*8 tlim       ! intermediate var for es look-up with estbl4
      real*8 estblf     ! statement function es look-up
      real*8 estbl4     ! statement function es look-up
!
!---------------------------Local workspace-----------------------------
!
      real*8 es
      real*8 qs
      real*8 epsilo 
      real*8 latvap 
      real*8 latice 
      real*8 rh2o   
      real*8 cpair  
      real*8 omeps              ! 1 - 0.622
      real*8 one
      integer i, j, k           ! Indices
      integer itype
!
!-------------------------- Statement Functions ------------------------
!
! Statement functions used in saturation vapor pressure table lookup
! there are two ways to use these three statement functions.
! For compilers that do a simple in-line expansion:
! => ttemp = tlim(t)
!    es    = estbl4(ttemp)
!
! For compilers that provide real optimization:
! => es    = estblf(t)
!
      tlim(td) = max(min(td,tmax),tmin)
!
      estblf(td) =  (tmin + int(tlim(td)-tmin) - tlim(td) + 1.0) &
                  *estbl(int(tlim(td)-tmin)+1) &
                  -(tmin + int(tlim(td)-tmin) - tlim(td)      ) &
                  *estbl(int(tlim(td)-tmin)+2)
!
      estbl4(td) =  (tmin+int(td-tmin)+1.0-td)*estbl(int(td-tmin)+1) &
                  + ( td-(tmin+int(td-tmin)) )*estbl(int(td-tmin)+2)
!
!-----------------------------------------------------------------------
!
      one    = 1.0
      epsilo = 0.622
      latvap = 2.5104e06
      latice = 3.336e5
      rh2o   = 4.61e2
      cpair  = 1004.64
      omeps  = 1.0 - epsilo
!
! Build es table
!
      call esinti (epsilo ,latvap ,latice ,rh2o ,cpair )
!
      do k = 1,plev
        do j = 1,plat
          do i = 1,plon
!
! Saturation specific humidity
!
            es = estblf( t(i,j,k) )
            qs = epsilo*es/(press(i,j,k) - omeps*es)
!   
! The following check is to avoid the generation of negative values
! that can occur in the upper stratosphere and mesosphere
!
            qs = min(one,qs)
!         
            if (qs .lt. 0.0) then
              qs = 1.0
            end if
!
! Compute Q from RH and Qs
!
            q(i,j,k) = q(i,j,k)*qs
!
          end do
        end do
      end do
!
      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      subroutine esinti(epslon  ,latvap  ,latice  ,rh2o    ,cpair   )
!-----------------------------------------------------------------------
!
! Initialize es lookup tables 
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Hack
! Standardized:      L. Buja, Jun 1992, Feb 1996
! Reviewed:          J. Hack, G. Taylor, Aug 1992
!                    J. Hack, Feb 1996
!
!-----------------------------------------------------------------------
      implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
      real*8 epslon        ! Ratio of h2o to dry air molecular weights 
      real*8 latvap        ! Latent heat of vaporization
      real*8 latice        ! Latent heat of fusion
      real*8 rh2o          ! Gas constant for water vapor
      real*8 cpair         ! Specific heat of dry air
!
!---------------------------Local workspace-----------------------------
!
      real*8 tmn           ! Minimum temperature entry in table
      real*8 tmx           ! Maximum temperature entry in table
      real*8 trice         ! Trans range from es over h2o to es over ice
      logical ip           ! Ice phase (true or false)
!
!-----------------------------------------------------------------------
!
! Specify control parameters first
!
      tmn   = 173.16
      tmx   = 375.16
      trice =  20.00
      ip    = .true.
!
! Call gestbl to build saturation vapor pressure table.
!
      call gestbl(tmn     ,tmx     ,trice   ,ip      ,epslon  , &
                  latvap  ,latice  ,rh2o    ,cpair   )
!
      return
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      subroutine gestbl(tmn     ,tmx     ,trice   ,ip      ,epsil   , &
                        latvap  ,latice  ,rh2o    ,cpair   )
!-----------------------------------------------------------------------
!
! Builds saturation vapor pressure table for later lookup procedure.
! Uses Goff & Gratch (1946) relationships to generate the table
! according to a set of free parameters defined below.  Auxiliary
! routines are also included for making rapid estimates (well with 1%)
! of both es and d(es)/dt for the particular table configuration.
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Hack
! Standardized:      L. Buja, Jun 1992,  Feb 1996
! Reviewed:          J. Hack, G. Taylor, Aug 1992
!                    J. Hack, Aug 1992
!
!-----------------------------------------------------------------------
      implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
      real*8 tmn          ! Minimum temperature entry in es lookup table
      real*8 tmx          ! Maximum temperature entry in es lookup table
      real*8 epsil        ! Ratio of h2o to dry air molecular weights
      real*8 trice        ! Transition range from es over range to es
!                         ! over ice
      real*8 latvap       ! Latent heat of vaporization
      real*8 latice       ! Latent heat of fusion
      real*8 rh2o         ! Gas constant for water vapor
      real*8 cpair        ! Specific heat of dry air
!
!---------------------------Local variables-----------------------------
!
      real*8 t            ! Temperature
      integer n           ! Increment counter
      integer lentbl      ! Calculated length of lookup table
      integer itype       ! Ice phase: 0 -> no ice phase
                          !            1 -> ice phase, no transition
                          !           -x -> ice phase, x degree transitn
      logical ip          ! Ice phase logical flag
!
!---------------------------- Commons ----------------------------------
!
! Common block and statement functions for saturation vapor pressure
! look-up procedure, J. J. Hack, February 1990
!
      integer plenest  ! length of saturation vapor pressure table
      parameter (plenest=250)
!
! Table of saturation vapor pressure values es from tmin degrees
! to tmax+1 degrees k in one degree increments.  ttrice defines the
! transition region where es is a combination of ice & water values
!
      common/comes/estbl(plenest) ,tmin  ,tmax  ,ttrice ,pcf(6) , &
                   icephs
!
      real*8 estbl      ! table values of saturation vapor pressure
      real*8 tmin       ! min temperature (K) for table
      real*8 tmax       ! max temperature (K) for table
      real*8 ttrice     ! transition range from es over H2O to over ice
      real*8 pcf        ! polynomial coeffs -> es transition h2o to ice
      logical icephs  ! false => saturation vapor press over water only
!
!-----------------------------------------------------------------------
!
! Set es table parameters
!
      tmin   = tmn       ! Minimum temperature entry in table
      tmax   = tmx       ! Maximum temperature entry in table
      ttrice = trice     ! Trans. range from es over h2o to es over ice
      icephs = ip        ! Ice phase (true or false)
!
! Set physical constants required for es calculation
!
      lentbl = tmax-tmin+2.000001
      if (lentbl .gt. plenest) then
         write(6,9000) tmax, tmin, plenest
         call abort     ! Abnormal termination
      end if
!
! Begin building es table.
! Check whether ice phase requested.
! If so, set appropriate transition range for temperature
!
      if (icephs) then
         if(ttrice.ne.0.0) then
            itype = -ttrice
         else
            itype = 1
         end if
      else
         itype = 0
      end if
!
      t = tmin - 1.0
      do n=1,lentbl
         t = t + 1.0
         call gffgch(t,estbl(n),itype)
      end do
!
      do n=lentbl+1,plenest
         estbl(n) = -99999.0
      end do
!
! Table complete -- Set coefficients for polynomial approximation of
! difference between saturation vapor press over water and saturation
! pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
! is valid in the range -40 < t < 0 (degrees C).
!
!                  --- Degree 5 approximation ---
!
      pcf(1) =  5.04469588506e-01
      pcf(2) = -5.47288442819e+00
      pcf(3) = -3.67471858735e-01
      pcf(4) = -8.95963532403e-03
      pcf(5) = -7.78053686625e-05
!
!                  --- Degree 6 approximation ---
!
!-----pcf(1) =  7.63285250063e-02
!-----pcf(2) = -5.86048427932e+00
!-----pcf(3) = -4.38660831780e-01
!-----pcf(4) = -1.37898276415e-02
!-----pcf(5) = -2.14444472424e-04
!-----pcf(6) = -1.36639103771e-06
!
      return
!
 9000 format('GESTBL: FATAL ERROR *********************************',/, &
           ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH', &
           ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/, &
           ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)
! 
      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      subroutine gffgch(t       ,es      ,itype   )
!-----------------------------------------------------------------------
!
! Computes saturation vapor pressure over water and/or over ice using
! Goff & Gratch (1946) relationships.  T (temperature), and itype are
! input parameters, while es (saturation vapor pressure) is an output
! parameter.  The input parameter itype serves two purposes: a value of
! zero indicates that saturation vapor pressures over water are to be
! returned (regardless of temperature), while a value of one indicates
! that saturation vapor pressures over ice should be returned when t is
! less than 273.16 degrees k.  If itype is negative, its absolute value
! is interpreted to define a temperature transition region below 273.16
! degrees k in which the returned saturation vapor pressure is a
! weighted average of the respective ice and water value.  That is, in
! the temperature range 0 => -itype degrees c, the saturation vapor
! pressures are assumed to be a weighted average of the vapor pressure
! over supercooled water and ice (all water at 0 c; all ice at -itype
! c).  Maximum transition range => 40 c
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Hack
! Standardized:      L. Buja, Jun 1992,  Feb 1996
! Reviewed:          J. Hack, G. Taylor, Aug 1992
!                    J. Hack, Feb 1996 
!
!-----------------------------------------------------------------------
      implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
      real*8 t          ! Temperature
      integer itype   ! Flag for ice phase and associated transition
!
! Output arguments
!
      real*8 es         ! Saturation vapor pressure
!
!---------------------------Local variables-----------------------------
!
      real*8 e1        ! Intermediate scratch variable for es over water
      real*8 e2        ! Intermediate scratch variable for es over water
      real*8 eswtr     ! Saturation vapor pressure over water
      real*8 f         ! Intermediate scratch variable for es over water
      real*8 f1        ! Intermediate scratch variable for es over water
      real*8 f2        ! Intermediate scratch variable for es over water
      real*8 f3        ! Intermediate scratch variable for es over water
      real*8 f4        ! Intermediate scratch variable for es over water
      real*8 f5        ! Intermediate scratch variable for es over water
      real*8 ps        ! Reference pressure (mb)
      real*8 t0        ! Reference temperature (freezing point of water)
      real*8 term1     ! Intermediate scratch variable for es over ice
      real*8 term2     ! Intermediate scratch variable for es over ice
      real*8 term3     ! Intermediate scratch variable for es over ice
      real*8 tr        ! Transition range for es over liq to es over ice
      real*8 ts        ! Reference temperature (boiling point of water)
      real*8 weight    ! Intermediate scratch variable for es transition
      real*8 tfreez, one
      integer itypo    ! Intermediate scratch variable for holding itype
!
!-----------------------------------------------------------------------
!      
      tfreez = 273.16
      one    = 1.0
!
! Check on whether there is to be a transition region for es
!
      if (itype.lt.0) then
        tr    = abs(float(itype))
        itypo = itype
        itype = 1
      else
        tr    = 0.0
        itypo = itype
      end if
      if (tr .gt. 40.0) then
        write(6,900) tr
        call abort                 ! Abnormal termination
      end if
!
      if(t .lt. (tfreez - tr) .and. itype.eq.1) go to 10
!
! Water
!
      ps = 1013.246
      ts = 373.16
      e1 = 11.344*(1.0 - t/ts)
      e2 = -3.49149*(ts/t - 1.0)
      f1 = -7.90298*(ts/t - 1.0)
      f2 = 5.02808*log10(ts/t)
      f3 = -1.3816*(10.0**e1 - 1.0)/10000000.0
      f4 = 8.1328*(10.0**e2 - 1.0)/1000.0
      f5 = log10(ps)
      f  = f1 + f2 + f3 + f4 + f5
      es = (10.0**f)*100.0
      eswtr = es
!
      if(t.ge.tfreez .or. itype.eq.0) go to 20
!
! Ice
!
   10 continue
      t0    = tfreez
      term1 = 2.01889049/(t0/t)
      term2 = 3.56654*log(t0/t)
      term3 = 20.947031*(t0/t)
      es    = 575.185606e10*exp(-(term1 + term2 + term3))
!
      if (t.lt.(tfreez - tr)) go to 20
!
! Weighted transition between water and ice
!
      weight = min((tfreez - t)/tr,one)
      es = weight*es + (1.0 - weight)*eswtr
!
   20 continue
      itype = itypo
      return
!
  900 format('GFFGCH: FATAL ERROR ******************************',/, &
             'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR', &
             ' PRESSURE, TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF', &
             ' 40.0 DEGREES C',/, ' TR = ',f7.2)
!
      end
 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine mass_fixer(plev    , plevp1  ,plat    ,plon    ,q    , &
                            hyai    , hybi    ,gw      ,gravit  ,ps0  , &
                            tmass0  , ps      )
!
!-----------------------------------------------------------------------
!
! Adjust atmospheric mass based upon Q.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     INPUTS
!
      integer plev   !  level dimension
      integer plevp1 !  level dimension + 1
      integer plat   !  latitude dimension
      integer plon   !  longitude dimension
!
      real*8 hyai(plevp1)
      real*8 hybi(plevp1)
      real*8 gw  (plat)
      real*8 q   (plon,plat,plev)   ! Specific humidity
      real*8 gravit                 ! acceleration of gravity ~ m/s^2
      real*8 ps0                    ! Ref. Surface pressure (10**5 Pa)
      real*8 tmass0                 ! Dry mass of Ref. atmosphere
!
!     INPUT/OUTPUT
!
      real*8 ps  (plon,plat     )   ! Surface Pressure
!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
!
      real*8 tmassf               ! unfixed mass of atmosphere
      real*8 qmass1               ! "a" contribution to mass of q
      real*8 qmass2               ! "b" contribution to mass of q
      real*8 pdela(plon,plev)     ! "a" contribution to del-P
      real*8 pdelb(plon,plev)     ! "b" contribution to del-P
      real*8 pssum                ! !
      real*8 dotproda             ! ! - accumulators
      real*8 dotprodb             ! !
      real*8 fixmas               ! Mass fix coefficient
      integer i, lat, k           ! Indices

      tmassf = 0.
      qmass1 = 0.
      qmass2 = 0.
!
! Compute pdel from "A" portion of hybrid vertical grid for later use in global integrals
!
      do k = 1,plev
         do i = 1,plon
            pdela(i,k) = (hyai(k+1) - hyai(k))*ps0
         end do
      end do
!
! Compute integrals of mass, moisture, and geopotential height
!
      do lat = 1,plat
!              
! Accumulate average mass of atmosphere
!
         do k = 1,plev
            do i = 1,plon
               pdelb(i,k) = (hybi(k+1) - hybi(k))*ps(i,lat)
            end do
         end do

         pssum  = 0.
         do i = 1,plon
            pssum  = pssum  + ps(i,lat)
         end do
         tmassf = tmassf + gw(lat)*pssum/plon
!
! Calculate global integrals needed for water vapor adjustment
!
         do k = 1,plev
            dotproda = 0.
            dotprodb = 0.
            do i = 1,plon
               dotproda = dotproda + q(i,lat,k)*pdela(i,k)
               dotprodb = dotprodb + q(i,lat,k)*pdelb(i,k)
            end do
            qmass1 = qmass1 + gw(lat)*dotproda/plon
            qmass2 = qmass2 + gw(lat)*dotprodb/plon
         end do

      end do                  ! end of latitude loop
!
! Normalize average mass, height
!
      tmassf = tmassf*.5/gravit
      qmass1 = qmass1*.5/gravit
      qmass2 = qmass2*.5/gravit
!
! Compute and apply an initial mass fix factor which preserves horizontal
! gradients of ln(ps).
!
      fixmas = (tmass0 + qmass1)/(tmassf - qmass2)
!
      do lat = 1,plat
         do i = 1,plon
            ps(i,lat) = ps(i,lat)*fixmas
         end do
      end do
!
      write(6,1000) fixmas
1000  format("             FIXMAS = ", f18.14/)

      return

      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! NCLFORTSTART
      subroutine read_namelist(NAMENUMBER,I_MYTMPDIR,I_MYOUTDIR,I_INPUTDIR, I_TMP_clean, &
                               I_REF_DATE,I_CASE,I_DYCORE,I_PRECISION,I_VORT_DIV_TO_UV,  &
                               I_SST_MASK,I_ICE_MASK,I_OUTPUT_PHIS,I_REGRID_ALL,         &
                               I_ADJUST_STATE_FROM_TOPO,I_fname_phis_output,             &
                               I_ftype_phis_output,I_fname_grid_info,I_fields,           &
                               I_source_files,I_fname_phis_input,I_ftype_phis_input,     &
                               I_fname0,I_fname1,I_fname2,I_fname3,I_fname4,I_fname5,    &
                               I_ftype0,I_ftype1,I_ftype2,I_ftype3,I_ftype4,I_ftype5,    &
                               I_fdate0,I_fdate1,I_fdate2,I_fdate3,I_fdate4,I_fdate5)
!
!-----------------------------------------------------------------------
!
! Read Processing Namelist
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!     INPUTS
!
      integer     NAMENUMBER
!
!     OUTPUT
!
      character*(*) I_MYTMPDIR
      character*(*) I_MYOUTDIR
      character*(*) I_INPUTDIR
      character*(*) I_TMP_clean
      character*(*) I_REF_DATE
      character*(*) I_CASE
      character*(*) I_DYCORE
      character*(*) I_PRECISION
      character*(*) I_VORT_DIV_TO_UV
      character*(*) I_SST_MASK
      character*(*) I_ICE_MASK
      character*(*) I_OUTPUT_PHIS
      character*(*) I_REGRID_ALL
      character*(*) I_ADJUST_STATE_FROM_TOPO
      character*(*) I_fname_phis_output
      character*(*) I_ftype_phis_output
      character*(*) I_fname_grid_info
      character*(*) I_fields
      character*(*) I_source_files
      character*(*) I_fname_phis_input
      character*(*) I_ftype_phis_input
      character*(*) I_fname0,I_fname1,I_fname2,I_fname3,I_fname4,I_fname5
      character*(*) I_ftype0,I_ftype1,I_ftype2,I_ftype3,I_ftype4,I_ftype5
      character*(*) I_fdate0,I_fdate1,I_fdate2,I_fdate3,I_fdate4,I_fdate5

      integer,parameter::MAXLEN=256
      character(Len=MAXLEN) MYTMPDIR
      character(Len=MAXLEN) MYOUTDIR
      character(Len=MAXLEN) INPUTDIR
      character(Len=MAXLEN) ESMF_interp
      character(Len=MAXLEN) ESMF_pole
      character(Len=MAXLEN) ESMF_clean
      character(Len=MAXLEN) TMP_clean
      character(Len=MAXLEN) REF_DATE
      character(Len=MAXLEN) CASE
      character(Len=MAXLEN) DYCORE
      character(Len=MAXLEN) PRECISION
      character(Len=MAXLEN) VORT_DIV_TO_UV
      character(Len=MAXLEN) SST_MASK
      character(Len=MAXLEN) ICE_MASK
      character(Len=MAXLEN) OUTPUT_PHIS
      character(Len=MAXLEN) REGRID_ALL
      character(Len=MAXLEN) ADJUST_STATE_FROM_TOPO
      character(Len=MAXLEN) MASS_FIX
      character(Len=MAXLEN) fname_phis_output
      character(Len=MAXLEN) ftype_phis_output
      character(Len=MAXLEN) fname_grid_info
      character(Len=MAXLEN) fields
      character(Len=MAXLEN) source_files
      character(Len=MAXLEN) fname_phis_input
      character(Len=MAXLEN) ftype_phis_input
      character(Len=MAXLEN) fname0,fname1,fname2,fname3,fname4,fname5
      character(Len=MAXLEN) ftype0,ftype1,ftype2,ftype3,ftype4,ftype5
      character(Len=MAXLEN) fdate0,fdate1,fdate2,fdate3,fdate4,fdate5

!
!-----------------------------------------------------------------------

! NCLEND

!---------------------------Local workspace-----------------------------
      character*2 CNUM
      integer     fname_phis_in

      namelist /makeic_nl/  MYTMPDIR,MYOUTDIR,INPUTDIR,TMP_clean,          &
                            REF_DATE,CASE,DYCORE,PRECISION,VORT_DIV_TO_UV, &
                            SST_MASK,ICE_MASK,OUTPUT_PHIS,REGRID_ALL,      &
                            ADJUST_STATE_FROM_TOPO,fname_phis_output,      &
                            ftype_phis_output,fname_grid_info,fields,      &
                            source_files, fname_phis_in,                   &
                            fname0,fname1,fname2,fname3,fname4,fname5,     &
                            ftype0,ftype1,ftype2,ftype3,ftype4,ftype5,     &
                            fdate0,fdate1,fdate2,fdate3,fdate4,fdate5

! open namelist file and read in values
!-------------------------------------------
      write(CNUM(1:2),'(i2.2)') NAMENUMBER
      open(unit=31,file='./Config/Config_makeIC-'//CNUM//'.nl',status='old')
      read(31,makeic_nl)
      close(unit=31)
!
! Trim the input string to length
!-----------------------------------
      I_MYTMPDIR              =trim(MYTMPDIR)
      I_MYOUTDIR              =trim(MYOUTDIR)
      I_INPUTDIR              =trim(INPUTDIR)
      I_TMP_clean             =trim(TMP_clean)
      I_REF_DATE              =trim(REF_DATE)
      I_CASE                  =trim(CASE)
      I_DYCORE                =trim(DYCORE)
      I_PRECISION             =trim(PRECISION)
      I_VORT_DIV_TO_UV        =trim(VORT_DIV_TO_UV)
      I_SST_MASK              =trim(SST_MASK)
      I_ICE_MASK              =trim(ICE_MASK)
      I_OUTPUT_PHIS           =trim(OUTPUT_PHIS)
      I_REGRID_ALL            =trim(REGRID_ALL)
      I_ADJUST_STATE_FROM_TOPO=trim(ADJUST_STATE_FROM_TOPO)
      I_fname_phis_output     =trim(fname_phis_output)
      I_ftype_phis_output     =trim(ftype_phis_output)
      I_fname_grid_info       =trim(fname_grid_info)
      I_fields                =trim(fields)
      I_source_files          =trim(source_files)
      I_fname_phis_input      =trim(fname_phis_input)
      I_ftype_phis_input      =trim(ftype_phis_input)
      I_fname0                =trim(fname0)
      I_fname1                =trim(fname1)
      I_fname2                =trim(fname2)
      I_fname3                =trim(fname3)
      I_fname4                =trim(fname4)
      I_fname5                =trim(fname5)
      I_ftype0                =trim(ftype0)
      I_ftype1                =trim(ftype1)
      I_ftype2                =trim(ftype2)
      I_ftype3                =trim(ftype3)
      I_ftype4                =trim(ftype4)
      I_ftype5                =trim(ftype5)
      I_fdate0                =trim(fdate0)
      I_fdate1                =trim(fdate1)
      I_fdate2                =trim(fdate2)
      I_fdate3                =trim(fdate3)
      I_fdate4                =trim(fdate4)
      I_fdate5                =trim(fdate5)
!
! Set the PHIS input file
!-----------------------------
      if(fname_phis_in.eq.0) then
       I_fname_phis_input=I_fname0
       I_ftype_phis_input=I_ftype0
      elseif(fname_phis_in.eq.1) then
       I_fname_phis_input=I_fname1
       I_ftype_phis_input=I_ftype1
      elseif(fname_phis_in.eq.2) then
       I_fname_phis_input=I_fname2
       I_ftype_phis_input=I_ftype2
      elseif(fname_phis_in.eq.3) then
       I_fname_phis_input=I_fname3
       I_ftype_phis_input=I_ftype3
      elseif(fname_phis_in.eq.4) then
       I_fname_phis_input=I_fname4
       I_ftype_phis_input=I_ftype4
      elseif(fname_phis_in.eq.5) then
       I_fname_phis_input=I_fname5
       I_ftype_phis_input=I_ftype5
      else

      endif

      return

      end

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
