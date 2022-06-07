#!/bin/csh -f
#
# Makes base case that can later be spawned to different
# lats and lons. Start with ARM SGP coords cuz we're comfortable
# there, lots of data etc. ...
# Example:
# $$> ./make_basecase.auto.csh -49.48 286.25 58 SAndes_x4


if ( "$#argv" != 4) then
  echo "Wrong number of arguments specified:"
  echo "  -arg 1 lat"
  echo "  -arg 2 lon"
  echo "  -arg 3 nlev"
  echo "  -arg 4 case string"
  exit
endif

set n = 1
set case_lat = "$argv[$n]"
set n = 2
set case_lon = "$argv[$n]"
set n = 3
set case_nlev = "$argv[$n]"
set n = 4
set loc_string = "$argv[$n]"

set srcpath=/project/amp/juliob/scam/
set scratchdir=/scratch/cluster/$USER
#set COMPSET=FSCAM
set COMPSET=2000_CAM60%SCAM_CLM50%SP_CICE5%PRES_DOCN%DOM_SROF_SGLC_SWAV
set case_year = 2010
#set case_mon = 10
#set case_day = 01
set case_mon = 07
set case_day = 01
#set case_year = 1999
#set case_mon = 02
#set case_day = 23
 
set case_date = $case_year$case_mon$case_day
set case_sdate = $case_year"-"$case_mon"-"$case_day

echo $case_date
echo $case_sdate

set laa = `echo $case_lat | cut -d '.' -f 1`
echo $laa
set loo = `echo $case_lon | cut -d '.' -f 1`
echo $loo

# set basecase name
set CASE="${loc_string}""_L""${case_nlev}"

# DEBUG this
#-----------------------------------
# create new basecase
#./create_newcase --debug --case  ../../cases/${CASE} --compset ${COMPSET} --res T42_T42 --user-mods-dir ../../cime_config/usermods_dirs/scam_STUB --walltime 01:00:00 --mach izumi --pecount 1 --compiler intel --run-unsupported

# NO debug
#-----------------------------------
./create_newcase --case  ../../cases/${CASE} --compset ${COMPSET} --res T42_T42 --driver mct --user-mods-dir ../../cime_config/usermods_dirs/scam_STUB --walltime 01:00:00 --mach izumi --pecount 1 --compiler intel --run-unsupported

cd ../../cases/${CASE}

#sed -i 's/intel\/18.0.3/intel\/20.0.1/' ./env_mach_specific.xml
#sed -i 's/intel\/mvapich2-2.3rc2-intel-18.0.3/intel\/mvapich2-2.1-qlc/' ./env_mach_specific.xml

./case.setup 

# DEBUG this
#-----------------------------------
#./xmlchange DEBUG=TRUE


# Archiving
#------------------------
./xmlchange DOUT_S_ROOT='/project/amp/${USER}/scam/archive/${CASE}'

# Append chnages to CAM configure options
#------------------------------------------
./xmlchange --append CAM_CONFIG_OPTS="-phys cam6 -nlev ${case_nlev}"
#./xmlchange --append CAM_CONFIG_OPTS="-phys cam_dev -nlev ${case_nlev}"
#./xmlchange CAM_CONFIG_OPTS="-dyn eul -scam -phys cam_dev -nlev ${case_nlev}"


# ATM_NCPL should be at least 192 to accomodate
# high wind cases in SH winter
#----------------------------------------------
./xmlchange ATM_NCPL=192

# Default to 123 days of runtime
# i.e., 123*96=11808
./xmlchange STOP_N=5952
./xmlchange START_TOD=00000
./xmlchange STOP_OPTION=nsteps

echo "scm_use_ana_iop = .true.">>user_nl_cam

echo "cld_macmic_num_steps=6">>user_nl_cam
#echo "deep_scheme = 'off'">>user_nl_cam

#echo "clubb_timestep=150.D0">>user_nl_cam
#echo "clubb_gamma_coef = 0.27D0">>user_nl_cam
#echo "clubb_c14 = 1.6D0">>user_nl_cam
#echo "clubb_l_trapezoidal_rule_zm = .false.">>user_nl_cam
#echo "clubb_l_trapezoidal_rule_zt = .false.">>user_nl_cam

echo "clubb_mf_nup = 100">>user_nl_cam
echo "clubb_mf_L0 = 50.D0">>user_nl_cam
echo "clubb_mf_Lopt = 3">>user_nl_cam
echo "clubb_mf_a0 = 1.D0">>user_nl_cam
echo "clubb_mf_b0 = 0.5D0">>user_nl_cam
echo "clubb_mf_alphturb = 3.D0">>user_nl_cam

echo "do_clubb_mf = .true.">>user_nl_cam
echo "do_clubb_mf_diag = .true.">>user_nl_cam
#echo "zmconv_num_cin = 1">>user_nl_cam
echo "use_gw_front = .false.">>user_nl_cam
echo "use_gw_convect_dp = .false.">>user_nl_cam

#Set case specific variables
./xmlchange PTS_LAT=${case_lat}  
./xmlchange PTS_LON=${case_lon}
./xmlchange RUN_STARTDATE=${case_sdate}
cp ../../cime_config/usermods_dirs/scam_STUB/scripts/STUB_iop.nc ./

ncap2 --overwrite -s "bdate=${case_date}" STUB_iop.nc STUB_iop.nc
ncap2 --overwrite -s "lat[lat]=${case_lat}" STUB_iop.nc STUB_iop.nc
ncap2 --overwrite -s "lon[lon]=${case_lon}" STUB_iop.nc STUB_iop.nc


pwd

echo "READY TO BUILD/SUBMIT "${CASE}
#./case.build
#./case.submit
exit
