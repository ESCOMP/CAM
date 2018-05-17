# Local pFUnit files.
set(pf_sources
  test_MGPacker.pf test_MGFieldPostProc.pf test_MGPostProc.pf)

# Sources to test.
set(sources_needed
  micro_mg_data.F90)
extract_sources("${sources_needed}" "${cam_sources}" test_sources)

# Do source preprocessing and add the executable.
create_pFUnit_test(micro_mg_data micro_mg_data_exe "${pf_sources}"
  "${test_sources}")

target_link_libraries(micro_mg_data_exe csm_share)
