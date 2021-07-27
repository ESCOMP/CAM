module spmd_dyn

    implicit none

    logical, public :: local_dp_map = .true.    ! flag indicates that mapping between dynamics
                                                ! and physics decompositions does not require
                                                ! interprocess communication
    integer, public :: block_buf_nrecs
    integer, public :: chunk_buf_nrecs

end module spmd_dyn
