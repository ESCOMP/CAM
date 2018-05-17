#!/bin/sh 
#

if [ $# -ne 3 ]; then
    echo "TCB_ccsm.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TCB_ccsm.$1.$2

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCB_ccsm.sh: CESM configure and build test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCB_ccsm.sh: CESM configure and build test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCB_ccsm.sh: this CESM configure and build test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CAM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CAM_TESTDIR}/${test_name} ${CAM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
        if [ $3 = "run_only" ]; then
            echo "TCB_ccsm.sh: CAM_RBOPTIONS set to run_only, this test needs to be built.  Try build_only or run_and_build"
            exit 2
        fi
    fi
fi

blddir=${CAM_TESTDIR}/${test_name}
if [ -d ${blddir} ]; then
    rm -rf ${blddir}
fi
mkdir -p ${blddir} 
if [ $? -ne 0 ]; then
    echo "TCB_ccsm.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${blddir}

echo "TCB_ccsm.sh: building ccsm executable; output in ${CAM_TESTDIR}/${test_name}/test.log" 
if [ -d ${CAM_TESTDIR}/case.$1.$2 ]; then
    rm -rf ${CAM_TESTDIR}/case.$1.$2
fi

# determine if chemistry preprocessor needs to be invoked
compset=${2%+*}
usrmech=${2#*+}

echo ${CCSM_MACH}
if [[ -n "$CCSM_MPILIB" ]]; then
    echo ${CCSM_MPILIB}
    mpiopt="-mpilib $CCSM_MPILIB"
else
    mpiopt=""
fi
echo "${CAM_ROOT}/cime/scripts/create_newcase --case ${CAM_TESTDIR}/case.$1.$2 \
    --res $1 --compset $compset  ${mpiopt} --run-unsupported "

${CAM_ROOT}/cime/scripts/create_newcase --case ${CAM_TESTDIR}/case.$1.$2 \
    --res $1 --compset $compset ${mpiopt} --run-unsupported  >test.log 2>&1
rc=$?
if [ $rc -eq 0 ]; then
    echo "TCB_ccsm.sh: create_newcase was successful" 
else
    echo "TCB_ccsm.sh: create_newcase failed, error from create_newcase= $rc" 
    echo "TCB_ccsm.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

cd ${CAM_TESTDIR}/case.$1.$2
echo "./xmlchange EXEROOT=${CAM_TESTDIR}/case.$1.$2/bld"
./xmlchange EXEROOT=${CAM_TESTDIR}/case.$1.$2/bld

echo "./xmlchange RUNDIR=${CAM_TESTDIR}/case.$1.$2/run"
./xmlchange RUNDIR=${CAM_TESTDIR}/case.$1.$2/run

# chemistry preprocessor
if [ $usrmech != $2 ]; then
   string1=`grep CAM_CONFIG_OPTS env_build.xml`
   string2=`echo $string1 | cut -d "=" -f 3`
   cfgstring=`echo $string2 | cut -d "\"" -f 2`
   echo "./xmlchange CAM_CONFIG_OPTS=""$cfgstring -usr_mech_infile ${CAM_SCRIPTDIR}/config_files/$usrmech -build_chem_proc"" "
   ./xmlchange CAM_CONFIG_OPTS="$cfgstring -usr_mech_infile ${CAM_SCRIPTDIR}/config_files/$usrmech -build_chem_proc" 
fi

#
# Override CESM pe layout.
#

for comp in ATM LND ICE OCN CPL GLC ROF WAV ESP; do
    echo "./xmlchange NTASKS_${comp}=$CAM_TASKS"
    ./xmlchange NTASKS_${comp}=$CAM_TASKS
#    echo "./xmlchange DEBUG=TRUE"
#    ./xmlchange DEBUG=TRUE
    echo "./xmlchange NTHRDS_${comp}=$CAM_THREADS"
    ./xmlchange NTHRDS_${comp}=$CAM_THREADS
done


echo "./case.setup"
./case.setup >> ${CAM_TESTDIR}/${test_name}/test.log 2>&1
rc=$?
if [ $rc -eq 0 ]; then
    echo "TCB_ccsm.sh: CESM configure was successful" 
else
    echo "TCB_ccsm.sh: CESM configure failed, error from configure= $rc" 
    echo "TCB_ccsm.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 5
fi


buildscript=`ls *.build`
echo "./$buildscript"
./$buildscript >> ${CAM_TESTDIR}/${test_name}/test.log 2>&1
rc=$?
if [ $rc -eq 0 ]; then
    echo "TCB_ccsm.sh: CESM build was successful" 
else
    echo "TCB_ccsm.sh: CESM build failed, error from build= $rc" 
    echo "TCB_ccsm.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
    echo "FAIL.job${JOBID}" > TestStatus
    exit 6
fi

cd ${blddir}
echo "TCB_ccsm.sh: CESM configure and build test passed"
echo "PASS" > TestStatus
if [ $CAM_RETAIN_FILES != "TRUE" ]; then
    echo "TCB_ccsm.sh: removing some unneeded files to save disc space" 
    for dir in atm glc ice lnd ocn rof wav cpl mct pio cesm csm_share lib
    do
        rm -rf ${CAM_TESTDIR}/case.$1.$2/$dir
    done
fi

exit 0
