#!/bin/csh -f
#
# Makes base case that can later be spawned to different
# lats and lons. Start with ARM SGP coords cuz we're comfortable
# there, lots of data etc. ...


if ( "$#argv" != 3) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 lat"
  echo "  -arg 2 lon"
  echo "  -arg 3 case string"
  exit
endif

set n = 1
set case_lat = "$argv[$n]"
set n = 2
set case_lon = "$argv[$n]"
set n = 3
set loc_string = "$argv[$n]"

set COMPSET=FSCAM

set src=cam6_3_041.dtsens

#set mach=izumi
#set queue=short
#set srcpath=/home/$USER/src
#set scratchdir=/scratch/cluster/$USER

set mach=cheyenne
set queue=share
set srcpath=/glade/u/home/$USER/src
set scratchdir=/glade/scratch/$USER

set case_year = 2010
set case_mon = 05
set case_day = 01

set case_date = $case_year$case_mon$case_day
set case_sdate = $case_year"-"$case_mon"-"$case_day

echo $case_date
echo $case_sdate

set laa = `echo $case_lat | cut -d '.' -f 1`
echo $laa
set loo = `echo $case_lon | cut -d '.' -f 1`
echo $loo

# set basecase name
set CASE="${src}_${COMPSET}_L58dev_CAMFORC_${loc_string}_${case_date}_c`date '+%y%m%d'`_test0"

# create new basecase
${srcpath}/${src}/cime/scripts/create_newcase --case ${scratchdir}/${CASE} --compset ${COMPSET} --res T42_T42 --user-mods-dir ${srcpath}/${src}/cime_config/usermods_dirs/scam_STUB --walltime 00:30:00 --mach ${mach} --pecount 1 --compiler intel --driver mct --queue ${queue} --run-unsupported

cd ${scratchdir}/${CASE}

#sed -i 's/intel\/18.0.3/intel\/20.0.1/' ./env_mach_specific.xml
#sed -i 's/intel\/mvapich2-2.3rc2-intel-18.0.3/intel\/mvapich2-2.1-qlc/' ./env_mach_specific.xml
./case.setup 

#./xmlchange DEBUG=TRUE
./xmlchange DOUT_S=FALSE

# Append to CAM configure options
./xmlchange --append CAM_CONFIG_OPTS='-phys cam_dev -nlev 58'

# ATM_NCPL should be at least 192 to accomodate
# high wind cases in SH winter
./xmlchange ATM_NCPL=96

# Default to 123 days of runtime
# i.e., 123*96=11808
./xmlchange STOP_N=11807
./xmlchange START_TOD=00000
./xmlchange STOP_OPTION=nsteps

echo "scm_use_ana_iop = .true.">>user_nl_cam

echo "cld_macmic_num_steps=3">>user_nl_cam
#echo "clubb_timestep=150.D0">>user_nl_cam
#echo "clubb_gamma_coef = 0.27D0">>user_nl_cam
#echo "clubb_c14 = 1.6D0">>user_nl_cam
#echo "clubb_l_trapezoidal_rule_zm = .false.">>user_nl_cam
#echo "clubb_l_trapezoidal_rule_zt = .false.">>user_nl_cam

#Set case specific variables
./xmlchange PTS_LAT=${case_lat}  
./xmlchange PTS_LON=${case_lon}
./xmlchange RUN_STARTDATE=${case_sdate}

cp ${srcpath}/${src}/cime_config/usermods_dirs/scam_STUB/scripts/STUB_iop.nc ./

ncap2 --overwrite -s "bdate=${case_date}" STUB_iop.nc STUB_iop.nc
ncap2 --overwrite -s "lat[lat]=${case_lat}" STUB_iop.nc STUB_iop.nc
ncap2 --overwrite -s "lon[lon]=${case_lon}" STUB_iop.nc STUB_iop.nc

pwd

echo "READY TO BUILD/SUBMIT "${CASE}
./case.build
./case.submit
exit
