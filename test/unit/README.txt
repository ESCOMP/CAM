# To run all CAM unit tests on caldera, run the following command:
#
# Also note that, on caldera, this requires 'module load all-python-libs'
#
# The creation of a temporary directory ensures that you are doing a completely
# clean build of the unit tests. (The use of the --clean flag to run_tests.py
# cleans most, but not all of the files created by the unit test build.) For
# rerunning the tests after an incremental change, you can instead use an
# existing build directory.
#
# The specification of -DUSE_CONTIGUOUS is a workaround for https://github.com/ESMCI/cime/issues/1250
# It can be removed once that issue is resolved

../../../../cime/scripts/fortran_unit_testing/run_tests.py --build-dir `mktemp -d --tmpdir=. unit_tests.XXXXXXXX` --cmake-args ' -DCPPDEFS=-DUSE_CONTIGUOUS=contiguous,'

