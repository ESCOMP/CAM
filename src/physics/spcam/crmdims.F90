module crmdims
#ifdef CRM
    use shr_kind_mod,    only: r8 => shr_kind_r8
    implicit none

       integer, parameter :: nclubbvars = 17

       integer, parameter ::  crm_nx=SPCAM_NX, crm_ny=SPCAM_NY, crm_nz=SPCAM_NZ
       real(r8), parameter :: crm_dx=SPCAM_DX, crm_dy=SPCAM_DX, crm_dt=SPCAM_DT
#endif
end module crmdims
