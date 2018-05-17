# Local pFUnit files.
set(pf_sources
  test_mg_processes.pf test_mg_size_utils.pf)

# Sources to test.
set(sources_needed
  micro_mg_utils.F90)
extract_sources("${sources_needed}" "${cam_sources}" test_sources)

# BUG(wjs, 2016-03-21, bugz 2301) micro_mg_utils unit tests currently give a lot
# of build errors; for now, commenting them out. We'd like the users to be able
# to see that the test was skipped but didn't fail, but CTest has no mechanism
# for this.  Instead, hack around it by adding a test that will always pass but
# is likely to stand out as a skipped test.
add_test(SKIPPED_micro_mg_utils true)

# # Do source preprocessing and add the executable.
# create_pFUnit_test(micro_mg_utils micro_mg_utils_exe "${pf_sources}"
#   "${test_sources}")

# target_link_libraries(micro_mg_utils_exe csm_share)
