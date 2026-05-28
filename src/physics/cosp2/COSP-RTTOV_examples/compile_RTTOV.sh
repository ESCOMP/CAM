# Run this script after unpacking the RTTOVv13.2 tarball into a directory named "RTTOV_src/".
# e.g. tar -xf ${RTTOV_TARBALL} -C ./RTTOV_src/
# Run from ./RTTOV_src

# Create a new ifx arch file
rm -f "RTTOV_src/build/arch/ifx"
echo '
FC=ifx
FC77=ifx
CC=gcc
LDFLAGS_ARCH=
CFLAGS_ARCH=
FFLAGS_ARCH=-fPIC -O3 -fp-model source
AR=ar r

# Loop unrolling causes ifort v13 and later to take a long time to compile these subroutines
FFLAGS_ARCH_rttov_opdep_9_ad=-fPIC -O3 -unroll0 -fp-model source
FFLAGS_ARCH_rttov_opdep_9_k=-fPIC -O3 -unroll0 -fp-model source
FFLAGS_ARCH_rttov_opdep_13_ad=-fPIC -O3 -unroll0 -fp-model source
FFLAGS_ARCH_rttov_opdep_13_k=-fPIC -O3 -unroll0 -fp-model source

# -fp-model source ensures more consistent floating point results

F2PY=f2py --fcompiler=intelem
F2PYFLAGS_ARCH="-fPIC"
F2PYLDFLAGS_ARCH= 
' >> "RTTOV_src/build/arch/ifx"

# set myarch variable
myarch="ifx"
clean="y"

# set installdir path
installdir="../RTTOV_build/"

# Compile with hdf5 and netcdf
hdf5=1
netcdf=1
# Set HDF5_PREFIX, FFLAGS_HDF5, and LDFLAGS_HDF5 in Makefile.inc
# Set NETCDF_PREFIX, FFLAGS_NETCDF and LDFLAGS_NETCDF in Makefile.inc
# Use NCAR HDF5 and NETCDF paths if on an NCAR machine, otherwise the user must set these environment variables before running this script.
echo '

HDF5_PREFIX  = ${NCAR_ROOT_HDF5}
FFLAGS_HDF5  = -D_RTTOV_HDF $(FFLAG_MOD)$(HDF5_PREFIX)/include
LDFLAGS_HDF5 = -L$(HDF5_PREFIX)/lib -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz

NETCDF_PREFIX  = ${NCAR_ROOT_NETCDF}
FFLAGS_NETCDF  = -D_RTTOV_NETCDF -I$(NETCDF_PREFIX)/include
LDFLAGS_NETCDF = -L$(NETCDF_PREFIX)/lib -lnetcdff

' >> "RTTOV_src/build/Makefile.inc"

# Do not compile with LAPACK, python interface, or GUI.
lapack=0
f2py=0
gui=0

# Prepare make command
cd "RTTOV_src/src/"
cmd_mkfile="../build/Makefile.PL RTTOV_HDF=${hdf5} RTTOV_F2PY=${f2py} RTTOV_USER_LAPACK=${lapack}"
cmd_clean="make ARCH=$myarch INSTALLDIR=$installdir clean $makeflags"
cmd_build="make ARCH=$myarch INSTALLDIR=$installdir $makeflags"

echo "Regenerating Makefiles using:"
echo "$ $cmd_mkfile"
echo
echo "Compiling RTTOV using:"
if [[ $clean = "y" ]]; then
    echo "$ $cmd_clean"
fi
echo "$ $cmd_build"


# Compile RTTOV
echo
echo "Regenerating Makefiles..."
$cmd_mkfile

echo
echo "Compiling RTTOV..."
if [[ $clean = "y" ]]; then
    $cmd_clean
fi

$cmd_build
if [[ $? -eq 0 ]]; then
    echo
    echo "RTTOV compiled successfully"
    echo
fi

# Download example file, all coefficients can also be downloaded by executing "rtcoef_rttov13/rttov_coef_download.sh" in the RTTOV source code. 
echo "Downloading coefficient files for test..."
wget -np -l1 "https://nwp-saf.eumetsat.int/downloads/rtcoef_rttov13/rttov13pred101L/rtcoef_eos_2_airs_l1c_7gas.H5" -P"../../RTTOV_coefs/"
wget -np -l1 "https://nwp-saf.eumetsat.int/downloads/rtcoef_rttov13/cldaer_visir/sccldcoef_eos_2_airs_l1c.H5" -P"../../RTTOV_coefs/"
echo "...complete"
