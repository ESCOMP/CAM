#! /bin/csh -f

#-----------------------------------
# Written in Jun 2012 by Jerry Olson
#-----------------------------------

set echo

# Set temporary directories for input and output analysis files

setenv MYTMPDIR               ./TMP
setenv MYTMPDIRO              ./OUTPUT

# Build fortran library

#/usr/local/bin/WRAPIT         MAKEIC.stub MAKEIC.f90
WRAPIT         MAKEIC.stub MAKEIC.f90

# Reference date to put on output file (YYYYMMDD)

setenv REF_DATE               20070901

# Output file format

setenv CASE                   ERAI_f09_L30  # Case name that will be appended to name of output file
setenv DYCORE                 fv        # Dycore ("eul" or "fv" are the current choices)
setenv PRECISION              float     # "double" or "float" are the current choices of output precision
setenv PTRM                    -1       # "M" spectral truncation (for "eul" dycore only; ignored for other dycores; "-1" = no trunc)
setenv PTRN                    -1       # "N" spectral truncation (for "eul" dycore only; ignored for other dycores; "-1" = no trunc)
setenv PTRK                    -1       # "K" spectral truncation (for "eul" dycore only; ignored for other dycores; "-1" = no trunc)
setenv PLAT                   192       # Number of latitudes  on output IC file
setenv PLON                   288       # Number of longitudes on output IC file
setenv PLEV                    30       # Number of vert levs  on output IC file
                                        # (if PLEV = 0, no vertical levels will be generated in file)

# File from which to pull hyai, hybi, hyam, hybm info to define OUPUT levels (must be a CAM file or a file with level info in CAM format)

setenv FNAME_lev_info         /glade/p/cesm/cseg/inputdata/atm/cam/inic/fv/cami-mam3_0000-01-01_0.9x1.25_L30_c100618.nc

# List of full input file pathnames (disk or HPSS) from which to pull fields to be regridded
# (up to 6 files)

setenv FNAME0                 /glade/p/rda/data/ds627.0/ei.oper.an.ml/200709/ei.oper.an.ml.regn128sc.2007090100
setenv FNAME1                 none
setenv FNAME2                 none
setenv FNAME3                 none
setenv FNAME4                 none
setenv FNAME5                 none

# Regrid ALL input fields, if the first input file is a CAM file (otherwise, just regrid the fields listed below)

setenv REGRID_ALL             False

# Time slice to pull from each file (YYYYMMDDSSSSS or time index (0, 1, 2, 3, etc.))

setenv FDATE                  2007090100000

# List of CAM fields to be regridded (must contain, at minimum, U, V [or US and VS, if fv dycore], T, Q, and PS fields )

setenv FIELDS                 U,V,T,Q,PS

# Input analysis file index in which each field can be found

setenv SOURCE_FILES           0,0,0,0,0

## Input file type (The "FTYPE" list maps to the above list of filenames)

##---------------------------------------------------------------------------------
## Current input file types:    CAM
##                              YOTC_PS_Z
##                              YOTC_sfc
##                              YOTC_sfc_fcst
##                              YOTC_sh
##                              ECMWF_gg
##                              ECMWF_sh
##                              NASA_MERRA
##                              NASA_MERRA_PREVOCA
##                              JRA_25
##                              Era_Interim_627.0_sc
##                              ERA40_ds117.2
##---------------------------------------------------------------------------------

setenv FTYPE                  Era_Interim_627.0_sc

# Adjust PS and near-surface T based on differences between input and output topography (PHIS)

setenv ADJUST_STATE_FROM_TOPO True

# File from which to pull input topography (PHIS) for use in T and Ps adjustment near surface

setenv FNAME_phis_input       $FNAME0
#setenv FTYPE_phis_input       CAM
setenv FTYPE_phis_input       Era_Interim_627.0_sc

# CAM File from which to pull output topography (PHIS; must already be on output grid) for use in T and Ps adjustment near surface

setenv FNAME_phis_output      /glade/p/cesm/cseg/inputdata/atm/cam/topo/USGS-gtopo30_0.9x1.25_remap_c051027.nc
setenv FTYPE_phis_output      FV_TOPOGRAPHY

# Processing options

setenv VORT_DIV_TO_UV         True      # U/V determined from vort/div input
setenv SST_MASK               False     # Use landfrac and icefrac masks to isolate and interpolate SSTs from Ts
                                        # (ignored if "SST_cpl" is not being output)
setenv ICE_MASK               False     # Use landfrac to isolate and interpolate ice fraction
                                        # (ignored if "ice_cov" is not being output)
setenv OUTPUT_PHIS            True      # Copy output PHIS to the output initial file.

setenv FNAME                  $FNAME0,$FNAME1,$FNAME2,$FNAME3,$FNAME4,$FNAME5

ncl < ./makeIC.ncl

exit
