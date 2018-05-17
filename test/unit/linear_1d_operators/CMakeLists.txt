# Local pFUnit files.
set(pf_sources
  test_diagonal.pf test_derivatives.pf test_arithmetic.pf)

# Sources to test.
set(sources_needed
  coords_1d.F90 linear_1d_operators.F90)
extract_sources("${sources_needed}" "${cam_sources}" test_sources)

# NOTE(wjs, 2016-04-29) The linear_1d_operators test gives segfaults with gnu
# (gfortran 5.3.0 on my Mac); for now, commenting it out. We'd like the users to
# be able to see that the test was skipped but didn't fail, but CTest has no
# mechanism for this.  Instead, hack around it by adding a test that will always
# pass but is likely to stand out as a skipped test.
add_test(SKIPPED_linear_1d_operators true)

# # Do source preprocessing and add the executable.
# create_pFUnit_test(linear_1d_operators linear_1d_operators_exe
#   "${pf_sources}" "${test_sources}")

# target_link_libraries(linear_1d_operators_exe csm_share)
