module spmd_dyn

  ! Purpose: SPMD implementation of CAM FV3 dynamics.
  
  implicit none
  private
  
  ! These variables are not used locally, but are set and used in phys_grid.
  ! They probably should be moved there.
  logical, public :: local_dp_map=.true.  ! flag indicates that mapping between dynamics
                                          !  and physics decompositions does not require
                                          !  interprocess communication
  integer, public :: block_buf_nrecs      ! number of local grid points (lon,lat,lev)
                                          !  in dynamics decomposition (including level 0)
  integer, public :: chunk_buf_nrecs      ! number of local grid points (lon,lat,lev)
                                          !  in physics decomposition (including level 0)
                                          !  assigned in phys_grid.F90
end module spmd_dyn
