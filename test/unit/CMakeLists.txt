cmake_minimum_required(VERSION 2.8)

# Unit test boilerplate.
list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})
include(CIME_initial_setup)

project(CAM_tests Fortran C)

# More unit test boilerplate.
include(CIME_utils)

# Useful paths.
set(CAMROOT ../../)
set(CESMROOT ${CAMROOT}../../)
set(SHAREROOT ${CESMROOT}/cime/src/share)

# Add share code to source list.
include_directories(${SHAREROOT}/include ${SHAREROOT}/util)
add_subdirectory(${SHAREROOT}/util csm_share)
add_subdirectory(${SHAREROOT}/unit_test_stubs/util csm_share_stubs)

# CAM physics.
add_subdirectory(${CAMROOT}src/physics/cam physics_cam)
add_subdirectory(${CAMROOT}src/utils utils_cam)

# Build part of csm_share as a library.
add_library(csm_share ${share_sources})
declare_generated_dependencies(csm_share "${share_genf90_sources}")
get_target_property(includes csm_share INCLUDE_DIRECTORIES)
if(NOT includes)
  unset(includes)
endif()
list(APPEND includes "${PFUNIT_INCLUDE_DIRS}")
set_target_properties(csm_share PROPERTIES
  INCLUDE_DIRECTORIES "${includes}")

# Pick up csm_share module files.
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Test subdirectories.
add_subdirectory(coords_1d)
add_subdirectory(linear_1d_operators)
add_subdirectory(micro_mg_data)
add_subdirectory(micro_mg_utils)
add_subdirectory(vdiff_lu_solver)
