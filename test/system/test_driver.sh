#!/bin/sh
#
# test_driver.sh:  driver for the testing of CAM with standalone scripts
#
# usage on hobart, izumi, leehill, cheyenne
# ./test_driver.sh
#
# **more details in the CAM testing user's guide, accessible
#   from the CAM developers web page

##
## General syntax help function
## Usage: help <exit status>
##

help () {
  local hname="Usage: `basename ${0}` [ OPTION [ OPTION ... ] ]"
  local hprefix="`echo ${hname} | tr '[!-~]' ' '`"
  hprefix="  "
  echo "${hname}  "
  echo "${hprefix} [ -b ] (support baseline scripts for cam5_2_12 and earlier)"
  echo "${hprefix} [ -e ] (email summary to $USER)"
  echo "${hprefix} [ -f ] (force batch submission -- avoids user prompt)"
  echo "${hprefix} [ -h ] (displays this help message)"
  echo "${hprefix} [ -i ] (interactive usage)"
  echo "${hprefix} [ -j ] (number of jobs for gmake)"
  echo "${hprefix} [ --baseline-dir <directory> ] (directory for saving baselines of cime tests)"
  echo "${hprefix} [ --no-baseline] (baselines of cime tests are not saved)"
  echo "${hprefix} [ --xml-driver <driver_name> ] (mct or nuopc)"
  echo "${hprefix} [ --cesm <test_name(s)> ] (default aux_cam)"
  echo "${hprefix} [ --rerun-cesm <test_id> ] (rerun the cesm tests with the --use-existing-flag)"
  echo "${hprefix} [ --namelists-only ] (Only perform namelist actions for tests.  Incompatible with --rerun-cesm.)"
  echo "${hprefix} [ --batch ] (Allow cime tests to run in parallel.)"
  echo ""
  echo "${hprefix} **pass environment variables by preceding above commands with:"
  echo "${hprefix}   'env var1=setting var2=setting '"
  echo ""
  echo "Supported ENVIRONMENT variables"
  echo "BL_TESTDIR:        Default = none (used to set baseline compare dir)"
  echo "CAM_ACCOUNT:       Default = none"
  echo "CAM_BATCHQ:        Default = machine dependent"
  echo "CAM_FC:            Default = machine dependent"
  echo "CAM_INPUT_TESTS:   Default = tests_pretag_<machine>[_<compiler>]"
  echo "CAM_RESTART_TASKS: Default = 64"
  echo "CAM_RETAIN_FILES:  Default = FALSE"
  echo "CAM_ROOT:          Default = set relative to CAM_SCRIPTDIR"
  echo "NB: If script is not called as ./`basename ${0}`, CAM_ROOT must be specified"
  echo "CAM_SCRIPTDIR:     Default = <current directory>"
  echo "CAM_TAG:           Default = none (used to set CESM baseline dir)"
  echo "CAM_TASKS:         Default = (depends on system)"
  echo "CAM_TESTDIR:       Default = <user_scratch_dir>/test-driver.<jobid>"
  echo ""
  echo "Less common ENVIRONMENT variables"
  echo "CALDERA_BATCHQ:    Default = caldera"
  echo "CAM_RBOPTIONS:     Default = build_only"
  echo "CAM_SOFF:          Default = none (stop of first test fail if TRUE)"
  echo "CIME_MODEL:        Default = none (should be set to cesm)"
  echo "EMAIL:             Default = $USER@ucar.edu"
  echo "SUMMARY_FILE:      Default = `pwd -P`/cam_test_summaries}"
  exit $1
}

##
## Error output function (should be handed a string)
##
perr() {
  echo -e "\nERROR: ${@}\n"
  help 1
}

## These variables may be overridden from the user's environment
EMAIL=${EMAIL:-"${USER}@ucar.edu"}

# These variables may be modified by script switches (./test_driver.sh -h)
cam_email_summary=false
cesm_test_suite="aux_cam"
force=false
gmake_j=0
interactive=false
use_existing=''
namelists_only=false
batch=false

# Understand where we are and where the CAM root and CIME reside
if [ -n "${CAM_ROOT}" ]; then
   # The user specified a CAM_ROOT, make sure it exists and that this
   # script was called from its test directory.
   if [ -d "${CAM_ROOT}" ]; then
     test_dir="${CAM_ROOT}/test/system"
     if [ -f "${test_dir}/test_driver.sh" ]; then
       # Check against this script.
       script_dir="$( cd $( dirname $0 ); pwd -P )"
       script_file="${script_dir}/test_driver.sh"
       if [ "${test_dir}/test_driver.sh" != "${script_file}" ]; then
         perr "CAM_ROOT test dir is ${test_dir} but script is ${script_file}"
       fi # Else, everything is fine
     else
       perr "No test_driver.sh found in ${test_dir}"
     fi
   else
     perr "CAM_ROOT, '${CAM_ROOT}', does not exist"
   fi
else
  # The user did not specify CAM_ROOT, find it relative to this script
  test_dir="$( dirname $0 )"
  export CAM_ROOT="$(dirname $(dirname $( cd ${test_dir}; pwd -P )))"
fi
# Now, find CIME_ROOT, first try a CAM standalone checkout
if [ -d "${CAM_ROOT}/cime/scripts" ]; then
  export CIME_ROOT="${CAM_ROOT}/cime"
else
  export CIME_ROOT="$( dirname $( dirname ${CAM_ROOT} ) )/cime"
  if [ ! -d "${CIME_ROOT}/scripts" ]; then
    perr "No CIME found from CAM_ROOT = '${CAM_ROOT}'"
  fi
fi

# Initialize variables which may not be set
submit_script_cime=''

while [ "${1:0:1}" == "-" ]; do
    case $1 in

        --baseline-dir )
            if [ $# -lt 2 ]; then
                perr "${1} requires a directory name)"
            fi
            baseline_dir="${2}"
            shift
            ;;

        --no-baseline ) no_baseline=false
             ;;

        -b ) export CAM_BASEBACK="YES"
             ;;

        --cesm )
            if [ $# -lt 2 ]; then
                perr "${1} requires a CESM test name or test suite name (e.g., aux_cam)"
            fi
            if [ "${2:0:1}" == "-" ]; then
                perr "Invalid CESM test name, '${2}'"
            fi
            cesm_test_suite="${2}"
            shift
            ;;

        -e ) cam_email_summary=true
             ;;

	-f ) force=true
             if  $interactive ; then
               echo "test_driver.sh: FATAL ERROR: -i and -f were set"
               exit 1
             fi
             ;;

	-h | --help )
             help 0
             ;;

	-i ) interactive=true
             if  $force ; then
               echo "test_driver.sh: FATAL ERROR: -i and -f were set"
               exit 1
             fi
             ;;

	-j ) shift; gmake_j=$1
             ;;

        --rerun-cesm )
            if [ $# -lt 2 ]; then
                perr "${1} requires a test_id from a previous run)"
            fi
            use_existing="${2}"
            shift
            if $namelists_only ; then
              echo "test_driver.sh: FATAL ERROR: --rerun-cesm and --namelists-only were set"
              exit 1
            fi
            ;;

        --xml-driver )
            if [ $# -lt 2 ]; then
                perr "${1} specify mct or nuopc)"
            fi
            xml_driver="${2}"
            shift
            ;;

        --namelists-only )
            namelists_only=true
            if [ "${use_existing}" != "" ]; then
              echo "test_driver.sh: FATAL ERROR: --namelists-only and --rerun-cesm were set"
              exit 1
            fi
            ;;

        --batch )
            batch=true
            ;;

    esac
    shift
done

# Currently, we don't support non-options, should we?
if [ $# -gt 0 ]; then
  perr "Unrecognized arguments: '$*'"
fi

#will attach timestamp onto end of script name to prevent overwriting
start_date="`date --iso-8601=seconds`"
cur_time=`date '+%H%M%S'`
date_str="`date '+%Y%m%d%H%M%S'`"

hostname=`hostname`

case $hostname in

    ##cheyenne
    ch* | r* )
    submit_script_cime="`pwd -P`/test_driver_cheyenne_cime_${cur_time}.sh"

    if [ -z "$CAM_ACCOUNT" ]; then
        echo "ERROR: Must set the environment variable CAM_ACCOUNT"
        exit 2
    fi

    if [ -z "$CAM_BATCHQ" ]; then
        export CAM_BATCHQ="regular"
    fi

    # wallclock for run job
    wallclock_limit="5:00:00"

    if [ $gmake_j = 0 ]; then
        gmake_j=36
    fi

    # run tests on 2 nodes using 18 tasks/node, 2 threads/task
    CAM_TASKS=36
    CAM_THREADS=2

    # change parallel configuration on 2 nodes using 32 tasks, 1 threads/task
    CAM_RESTART_TASKS=32
    CAM_RESTART_THREADS=1

    mach_workspace="/glade/scratch"

    # Check for CESM baseline directory
    if [ -n "${BL_TESTDIR}" ] && [ ! -d "${BL_TESTDIR}" ]; then
        echo "CESM_BASELINE ${BL_TESTDIR} not found.  Check BL_TESTDIR for correct tag name."
        exit
    fi

#-------------------------------------------

cat > ${submit_script_cime} << EOF
#!/bin/bash
#
#PBS -N cime-tests
#PBS -q $CAM_BATCHQ
#PBS -A $CAM_ACCOUNT
#PBS -l walltime=$wallclock_limit
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -j oe
#PBS -l inception=login

EOF

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;
    ##derecho
    derecho* | dec* )
    submit_script_cime="`pwd -P`/test_driver_derecho_cime_${cur_time}.sh"

    if [ -z "$CAM_ACCOUNT" ]; then
        echo "ERROR: Must set the environment variable CAM_ACCOUNT"
        exit 2
    fi

    if [ -z "$CAM_BATCHQ" ]; then
        export CAM_BATCHQ="main"
    fi

    # wallclock for run job
    wallclock_limit="5:00:00"

    if [ $gmake_j = 0 ]; then
        gmake_j=128
    fi

    # run tests on 1 node using 64 tasks/node, 2 threads/task
    # These settings are ignored on cheyenne and derecho.
    # PE layouts come from config_pes.xml.
    CAM_TASKS=64
    CAM_THREADS=2

    # change parallel configuration on 1 nodes using 32 tasks, 1 threads/task
    # These settings are ignored on cheyenne and derecho.
    # PE layouts come from config_pes.xml.
    CAM_RESTART_TASKS=32
    CAM_RESTART_THREADS=1

    mach_workspace="/glade/derecho/scratch"

    # Check for CESM baseline directory
    if [ -n "${BL_TESTDIR}" ] && [ ! -d "${BL_TESTDIR}" ]; then
        echo "CESM_BASELINE ${BL_TESTDIR} not found.  Check BL_TESTDIR for correct tag name."
        exit 3
    fi

#-------------------------------------------

cat > ${submit_script_cime} << EOF
#!/bin/bash
#
#PBS -N cime-tests
#PBS -q $CAM_BATCHQ
#PBS -A $CAM_ACCOUNT
#PBS -l walltime=$wallclock_limit
#PBS -l select=1:ncpus=128:mpiprocs=128
#PBS -j oe

EOF

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;


    ##hobart
    hob* | h[[:digit:]]* )
    submit_script_cime="`pwd -P`/test_driver_hobart_cime_${cur_time}.sh"
    export PATH=/cluster/torque/bin:${PATH}

    # Default setting is 12 hr in the long queue; the short queue only
    # allows 1 hr runs.
    wallclock_limit="12:00:00"
    gmake_j=24
    if [ -z "$CAM_BATCHQ" ]; then
        export CAM_BATCHQ="long"
    elif [[ "$CAM_BATCHQ" == short ]]; then
        wallclock_limit="1:00:00"
    fi

    if [ $gmake_j = 0 ]; then
        gmake_j=24
    fi

    if [ -z "$CAM_TASKS" ]; then
        CAM_TASKS=24
    fi
    if [ -z "$CAM_RESTART_TASKS" ]; then
        CAM_RESTART_TASKS=$(( $CAM_TASKS / 2))
    fi

    mach_workspace="/scratch/cluster"

    # Check for CESM baseline directory
    if  [ -n "{$BL_TESTDIR}" ] && [ ! -d "${BL_TESTDIR}" ]; then
        echo "CESM_BASELINE ${BL_TESTDIR} not found.  Check BL_TESTDIR for correct tag name."
        exit
    fi

#-------------------------------------------

cat > ${submit_script_cime} << EOF
#!/bin/bash
#
# Name of the queue (CHANGE THIS if needed)
#PBS -q $CAM_BATCHQ
# Number of nodes (CHANGE THIS if needed)
#PBS -l walltime=$wallclock_limit,nodes=1:ppn=24
# output file base name
#PBS -N cime-tests
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V
# End of options

EOF

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##izumi
    izu* | i[[:digit:]]* )

    submit_script_cime="`pwd -P`/test_driver_izumi_cime_${cur_time}.sh"
    export PATH=/cluster/torque/bin:${PATH}

    # Default setting is 12 hr in the long queue; the short queue only
    # allows 1 hr runs.
    wallclock_limit="12:00:00"
    gmake_j=24
    if [ -z "$CAM_BATCHQ" ]; then
        export CAM_BATCHQ="long"
    elif [[ "$CAM_BATCHQ" == short ]]; then
        wallclock_limit="1:00:00"
    fi

    if [ $gmake_j = 0 ]; then
        gmake_j=24
    fi

    if [ -z "$CAM_TASKS" ]; then
        CAM_TASKS=24
    fi
    if [ -z "$CAM_RESTART_TASKS" ]; then
        CAM_RESTART_TASKS=$(( $CAM_TASKS / 2))
    fi

    mach_workspace="/scratch/cluster"

    # Check for CESM baseline directory
    if  [ -n "{$BL_TESTDIR}" ] && [ ! -d "${BL_TESTDIR}" ]; then
        echo "CESM_BASELINE ${BL_TESTDIR} not found.  Check BL_TESTDIR for correct tag name."
        exit
    fi

#-------------------------------------------

cat > ${submit_script_cime} << EOF
#!/bin/bash
#
# Name of the queue (CHANGE THIS if needed)
#PBS -q $CAM_BATCHQ
# Number of nodes (CHANGE THIS if needed)
#PBS -l walltime=$wallclock_limit,nodes=1:ppn=24
# output file base name
#PBS -N cime-tests
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V
# End of options

EOF

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##casper
    casper* | crhtc* )
    submit_script_cime="`pwd -P`/test_driver_casper_cime_${cur_time}.sh"

    if [ -z "$CAM_ACCOUNT" ]; then
        echo "ERROR: Must set the environment variable CAM_ACCOUNT"
        exit 2
    fi

    if [ -z "$CAM_BATCHQ" ]; then
        export CAM_BATCHQ="casper"
    fi

    # wallclock for run job
    wallclock_limit="00:59:00"

    if [ $gmake_j = 0 ]; then
        gmake_j=36
    fi

    # run tests on 1 nodes using 18 tasks/node, 2 threads/task
    CAM_TASKS=18
    CAM_THREADS=2

    # change parallel configuration on 1 nodes using 32 tasks, 1 threads/task
    CAM_RESTART_TASKS=32
    CAM_RESTART_THREADS=1

    mach_workspace="/glade/scratch"

    # Check for CESM baseline directory
    if [ -n "${BL_TESTDIR}" ] && [ ! -d "${BL_TESTDIR}" ]; then
        echo "CESM_BASELINE ${BL_TESTDIR} not found.  Check BL_TESTDIR for correct tag name."
        exit
    fi

#-------------------------------------------

cat > ${submit_script_cime} << EOF
#!/bin/bash
#
#PBS -N cime-tests
#PBS -q $CAM_BATCHQ
#PBS -A $CAM_ACCOUNT
#PBS -l walltime=$wallclock_limit
#PBS -l select=1:ncpus=36:mpiprocs=36:mem=300GB
#PBS -j oe
#PBS -V
EOF

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    * ) echo "ERROR: machine $hostname not currently supported"; exit 1 ;;
esac


##vvvvvvvvvvvvvvvvvvvvvv start CAM aux test suite vvvvvvvvvvvvvvvvvvvvvvvvvvvv

cesm_test_mach=""
comp=""
if [ "${hostname:0:4}" == "chey" ]; then
  cesm_test_mach="cheyenne"
fi
if [ "${hostname:0:5}" == "derec" ] || [ "${hostname:0:3}" == "dec" ]; then
  cesm_test_mach="derecho"
fi
if [ "${hostname:0:6}" == "hobart" ]; then
  cesm_test_mach="hobart"
fi
if [ "${hostname:0:5}" == "izumi" ]; then
  cesm_test_mach="izumi"
fi
if [ "${hostname:0:6}" == "casper" ] || [ "${hostname:0:5}" == "crhtc" ]; then
  cesm_test_mach="casper"
fi
if [ -n "${CAM_FC}" ]; then
  comp="_${CAM_FC,,}"
fi

if [ "${cesm_test_suite}" != "none" -a -n "${cesm_test_mach}" ]; then
  if [ "${hostname:0:5}" != "izumi" ] && [ "${hostname:0:7}" != "derecho" ]; then
    module load python
  fi


  for cesm_test in ${cesm_test_suite}; do
    # Force derecho to run the cheyenne testlist.
    # After the transition to derecho is completed, this if statement can be removed and 
    # just the else needs to remain.
    if [ "${cesm_test_mach}" == "derecho" ]; then  
      testargs="--xml-category ${cesm_test} --xml-machine cheyenne --mach ${cesm_test_mach} --retry 2"
    else
      testargs="--xml-category ${cesm_test} --xml-machine ${cesm_test_mach} --retry 2"
    fi

    if [ -n "${use_existing}" ]; then
      test_id="${use_existing}"
    else
      test_id=${cesm_test}${comp}"_"${date_str}
    fi
    currdir="`pwd -P`"
    logfile="${currdir}/${test_id}.log"
    # Create an empty logfile so that other tasks can append to it
    if [ -f "${logfile}" ]; then
      rm -f ${logfile}
    fi
    touch ${logfile}
    script_dir="${CIME_ROOT}/scripts"
    if [ ! -d "${script_dir}" ]; then
      echo "ERROR: CIME scripts dir not found at ${script_dir}"
      exit 1
    fi
    if [ ! -x "${script_dir}/create_test" ]; then
      echo "ERROR: create_test script dir not found in ${script_dir}"
      exit 1
    fi

    ## If this is a Nag test, run the r8 and git tests
    if [ "${comp}" == "_nag" ]; then
      sepstr="################################################################"
      echo "${sepstr}" | tee -a ${logfile}
      ark_file="/fs/cgd/csm/tools/addrealkind/addrealkind"
      tr8_script="${CAM_ROOT}/test/system/TR8.sh"
      export ADDREALKIND_EXE="${ark_file}"; ${tr8_script} | tee -a ${logfile}
      res=${PIPESTATUS[0]}
      if [ $res -eq 0 ]; then
        echo "TR8 test PASS" | tee -a ${logfile}
      else
        echo "TR8 test FAIL, rc = $res" | tee -a ${logfile}
      fi
      echo "${sepstr}" | tee -a ${logfile}
      ${CAM_ROOT}/test/system/TGIT.sh | tee -a ${logfile}
      res=${PIPESTATUS[0]}
      if [ $res -eq 0 ]; then
        echo "TGIT test PASS" | tee -a ${logfile}
      else
        echo "TGIT test FAIL, rc = $res" | tee -a ${logfile}
      fi
      echo "${sepstr}" | tee -a ${logfile}
    fi

    ## Setup CESM work directory
    if [ "${hostname:0:6}" == "casper" ] || [ "${hostname:0:5}" == "crhtc" ]; then
       ## Would fail to compile on Casper with long folder name
       cesm_testdir=$mach_workspace/$LOGNAME/$cesm_test
    else
       cesm_testdir=$mach_workspace/$LOGNAME/$test_id
    fi

    if [ -e ${cesm_testdir} ]; then
      if [ -n "${use_existing}" ]; then
        echo " Using existing tests in ${cesm_testdir}"
      else
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
        echo "!! ERROR: ${cesm_testdir} already exists and << --rerun-cesm >> was not specified "
        echo "!!        Either remove ${cesm_testdir} or specify << --rerun-cesm >> "
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
        exit 1
      fi
    else
      mkdir $cesm_testdir
    fi

    if [ -n "${CAM_FC}" ]; then
      testargs="${testargs} --xml-compiler ${CAM_FC,,}"
    else
      testargs="${testargs} --xml-compiler intel"
    fi
    case $hostname in
        # cheyenne
        chey* | r* )
          testargs="${testargs} --queue ${CAM_BATCHQ} --test-root ${cesm_testdir} --output-root ${cesm_testdir}"
          ;;
        # derecho
        derec* | dec* )
          testargs="${testargs} --queue ${CAM_BATCHQ} --test-root ${cesm_testdir} --output-root ${cesm_testdir}"
          ;;
        # casper
        casper* | crhtc* )
          testargs="${testargs} --queue ${CAM_BATCHQ} --test-root ${cesm_testdir} --output-root ${cesm_testdir}"
          ;;
        *)
          if  $batch; then
            testargs="${testargs} --queue ${CAM_BATCHQ} --test-root ${cesm_testdir} --output-root ${cesm_testdir}"
          else
            testargs="${testargs} --test-root ${cesm_testdir} --output-root ${cesm_testdir}"
            testargs="${testargs} --no-batch"
          fi
    esac
    if [ -n "${CAM_ACCOUNT}" ]; then
      testargs="${testargs} --project ${CAM_ACCOUNT}"
    fi
    testargs="${testargs} --test-id ${test_id}"
    if [ -n "${BL_TESTDIR}" ]; then
      testargs="${testargs} --compare ${BL_TESTDIR} "
    fi
    if [ -n "${use_existing}" ]; then
      testargs="${testargs} --use-existing -o "
    fi
    if  $namelists_only ; then
      testargs="${testargs} --namelists-only "
    fi
    #  Check for a change in BL_TESTDIR                                                                   #
    if [ -n "${BL_TESTDIR}" ] && [ "${use_existing}" != "" ]; then
        #Check if BL_TESTDIR changed
        cmd="query_testlists --xml-category $cesm_test --xml-machine  ${cesm_test_mach}"
        if [ -n "${CAM_FC}" ]; then
            cmd="${cmd} --xml-compiler ${CAM_FC,,}"
        else
            cmd="${cmd} --xml-compiler intel"
        fi
        cmd="${CIME_ROOT}/scripts/"$cmd
        cime_testlist=`$cmd`
        for i in $(echo $cime_testlist | tr " " "\n")
        do
          if [[ $i =~ ${cesm_test_mach} ]]; then
            orig_baseline=`cd $cesm_testdir/$i*$test_id && ./xmlquery BASELINE_NAME_CMP --value`
            if [ $orig_baseline != ${BL_TESTDIR} ]; then
              echo "Changing BL_TESTDIR for $i."
                `cd $cesm_testdir/$i*$test_id && ./xmlchange BASELINE_NAME_CMP=$BL_TESTDIR`
              if [[ $i == ERI* ]]; then #Need to do special stuff to get ERI to rerun with new baseline.
                `cd $cesm_testdir/$i*$test_id && sed -i '/RUN/c\FAIL '$i' RUN' TestStatus`
                result=`cd $cesm_testdir/$i*$test_id && pwd && ./.case.test --reset -s`
              else
                `cd $cesm_testdir/$i*$test_id && sed -i '/RUN/c\PEND '$i' RUN' TestStatus`
              fi
            else
              echo "Checking for changed BL_TESTDIR for $i."
            fi
          fi
        done
    fi

    if [ "$no_baseline" != false ]; then
       if [ -n "${baseline_dir}" ]; then
         testargs="${testargs} --generate ${baseline_dir}"
       else
        testargs="${testargs} --generate ${cesm_testdir}/baselines"
      fi
    fi

    if [ -n "${xml_driver}" ]; then
      testargs="${testargs} --xml-driver ${xml_driver}"
    fi

    echo ""
    echo "CESM test results will be in: ${cesm_testdir}" | tee -a ${logfile}
    echo "Running ./create_test ${testargs}"             | tee -a ${logfile}

    if [ "${hostname:0:2}" == "ch" ]; then
      echo "cd ${script_dir}" >> ${submit_script_cime}
      echo "module load python" >> ${submit_script_cime}
      echo './create_test' ${testargs} >> ${submit_script_cime}
      chmod u+x ${submit_script_cime}
      qsub ${submit_script_cime}
    fi

    if [ "${hostname:0:2}" == "de" ]; then
      echo "cd ${script_dir}" >> ${submit_script_cime}
      echo './create_test' ${testargs} >> ${submit_script_cime}
      chmod u+x ${submit_script_cime}
      qsub ${submit_script_cime}
   fi

    if [ "${hostname:0:6}" == "hobart" ]; then
      echo "cd ${script_dir}" >> ${submit_script_cime}
      echo './create_test' ${testargs} >> ${submit_script_cime}
      if [ "${submit_script}" != "${submit_script_cime}" ]; then
        chmod u+x ${submit_script_cime}
        qsub ${submit_script_cime}
      fi
    fi

    if [ "${hostname:0:5}" == "izumi" ]; then
      echo "cd ${script_dir}" >> ${submit_script_cime}
      echo './create_test' ${testargs} >> ${submit_script_cime}
      if [ "${submit_script}" != "${submit_script_cime}" ]; then
        chmod u+x ${submit_script_cime}
        qsub ${submit_script_cime}
      fi
    fi

    if [ "${hostname:0:6}" == "casper" ] || [ "${hostname:0:5}" == "crhtc" ]; then
      echo "cd ${script_dir}" >> ${submit_script_cime}
      echo "module load python" >> ${submit_script_cime}
      echo './create_test' ${testargs} >> ${submit_script_cime}
      chmod u+x ${submit_script_cime}
      qsub ${submit_script_cime}
    fi

  done
fi

##^^^^^^^^^^^^^^^^^^^^^^ start CAM aux test suite ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

exit 0
