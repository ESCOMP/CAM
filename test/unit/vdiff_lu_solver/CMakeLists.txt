# Local pFUnit files.
set(pf_sources
  test_fd_solver.pf test_fv_solver.pf)

# Sources to test.
set(sources_needed
  coords_1d.F90 linear_1d_operators.F90 vdiff_lu_solver.F90)
extract_sources("${sources_needed}" "${cam_sources}" test_sources)

# Do source preprocessing and add the executable.
create_pFUnit_test(vdiff_lu_solver vdiff_lu_solver_exe "${pf_sources}"
  "${test_sources}")

target_link_libraries(vdiff_lu_solver_exe csm_share)
