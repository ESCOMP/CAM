# Local pFUnit files.
set(pf_sources
  test_coords_1d.pf)

# Sources to test.
set(sources_needed
  coords_1d.F90)
extract_sources("${sources_needed}" "${cam_sources}" test_sources)

# Do source preprocessing and add the executable.
create_pFUnit_test(coords_1d coords_1d_exe "${pf_sources}"
  "${test_sources}")

target_link_libraries(coords_1d_exe csm_share)
