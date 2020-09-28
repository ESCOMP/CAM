#!/bin/sh -f
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
  echo "${hprefix} [ --archive-cime <directory> ] (directory for archiving baselines of cime tests)"
  echo "${hprefix} [ --cesm <test_name(s)> ] (default aux_cam)"
  echo "${hprefix} [ --no-cesm ] (do not run any CESM test or test suite)"
  echo "${hprefix} [ --no-cam ] (do not run CAM regression tests"
  echo "${hprefix} [ --rerun-cesm <test_id> ] (rerun the cesm tests with the --use-existing-flag)"
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
  echo "CAM_TASKS:         Default = 32"
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
SUMMARY_FILE="${SUMMARY_FILE:-`pwd -P`/cam_test_summaries}"

# These variables may be modified by script switches (./test_driver.sh -h)
cam_email_summary=false
cesm_test_suite="aux_cam"
force=false
gmake_j=0
interactive=false
run_cam_regression=true
use_existing=''

# Initialize variables which may not be set
submit_script=''
submit_script_cb=''
submit_script_cime=''

while [ "${1:0:1}" == "-" ]; do
    case $1 in

        --archive-cime )
            if [ $# -lt 2 ]; then
                perr "${1} requires a directory name)"
            fi
            if [ -z "$BL_TESTDIR" ]; then
                echo "\$BL_TESTDIR needs to be set when using --archive-cime."
                exit 1
            fi
            archive_dir="${BL_TESTDIR%/*}/${2}"
            shift
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

        --no-cam )
            run_cam_regression=false
            ;;

        --no-cesm )
            cesm_test_suite="none"
            ;;

        --rerun-cesm )
            if [ $# -lt 2 ]; then
                perr "${1} requires a test_id from a previous run)"
            fi
            use_existing="${2}"
            shift
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

hostname=`hostname`

case $hostname in

    ##cheyenne
    ch* | r* )
    submit_script="`pwd -P`/test_driver_cheyenne_${cur_time}.sh"
    submit_script_cb="`pwd -P`/test_driver_cheyenne_cb_${cur_time}.sh"
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

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv

cat > ${submit_script_cb} << EOF
#!/bin/sh
#
#PBS -N test_dr
#PBS -q $CAM_BATCHQ
#PBS -A $CAM_ACCOUNT
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -j oe
#PBS -l inception=login

export TMPDIR=/glade/scratch/$USER

if [ -n "\$PBS_JOBID" ]; then    #batch job
   export JOBID=\`echo \${PBS_JOBID} | cut -f1 -d'.'\`
   initdir=`pwd -P`
   interactive=false
else
   interactive=true
fi

export CAM_RBOPTIONS="build_only"

## create_newcase looks for account number in ACCOUNT environment variable
export ACCOUNT=$CAM_ACCOUNT

# tasks and threads need to be set in the cb script because TCB_ccsm.sh uses
# them to set the pe_layout
export CAM_THREADS=$CAM_THREADS
export CAM_TASKS=$CAM_TASKS

source /glade/u/apps/ch/opt/lmod/7.5.3/lmod/lmod/init/sh

module load intel/19.0.5
module load mkl
module list

export INC_NETCDF=\${NCAR_ROOT_NETCDF}/include
export LIB_NETCDF=\${NCAR_ROOT_NETCDF}/lib

export CFG_STRING="-cc mpicc -fc mpif90 -fc_type intel -ldflags -mkl=cluster"
export MAKE_CMD="gmake -j $gmake_j"
export CCSM_MACH="cheyenne"
export MACH_WORKSPACE="$mach_workspace"
dataroot=${CESMDATAROOT}
echo_arg="-e"
input_file="tests_pretag_cheyenne"

EOF

#-------------------------------------------

cat > ${submit_script} << EOF
#!/bin/sh
#
#PBS -N test_dr
#PBS -q $CAM_BATCHQ
#PBS -A $CAM_ACCOUNT
#PBS -l walltime=$wallclock_limit
#PBS -l select=2:ncpus=36:mpiprocs=18:ompthreads=2
#PBS -j oe

export TMPDIR=/glade/scratch/$USER

if [ -n "\$PBS_JOBID" ]; then    #batch job
   export JOBID=\`echo \${PBS_JOBID} | cut -f1 -d'.'\`
   initdir=`pwd -P`
   interactive=false
else
   interactive=true
fi

export CAM_RBOPTIONS="run_only"
ulimit -c unlimited

##omp threads
export OMP_STACKSIZE=256M
export CAM_THREADS=$CAM_THREADS
export CAM_RESTART_THREADS=$CAM_RESTART_THREADS

##mpi tasks
export CAM_TASKS=$CAM_TASKS
export CAM_RESTART_TASKS=$CAM_RESTART_TASKS

##Cheyenne hacks to avoid MPI_LAUNCH_TIMEOUT
MPI_IB_CONGESTED=1
MPI_LAUNCH_TIMEOUT=40

source /glade/u/apps/ch/opt/lmod/7.5.3/lmod/lmod/init/sh

module load intel/19.0.5
module load mkl
module list

export CCSM_MACH="cheyenne"
export MACH_WORKSPACE="$mach_workspace"
export CPRNC_EXE=${CESMDATAROOT}/tools/cime/tools/cprnc/cprnc.cheyenne

dataroot=${CESMDATAROOT}

echo_arg="-e"

input_file="tests_pretag_cheyenne"

EOF

#-------------------------------------------

cat > ${submit_script_cime} << EOF
#!/bin/bash
#
#PBS -N cime-tests
#PBS -q $CAM_BATCHQ
#PBS -A $CAM_ACCOUNT
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -j oe
#PBS -l inception=login

EOF

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##hobart
    hob* | h[[:digit:]]* )
    submit_script="`pwd -P`/test_driver_hobart_${cur_time}.sh"
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

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ${submit_script} << EOF
#!/bin/sh
#

# Name of the queue (CHANGE THIS if needed)
#PBS -q $CAM_BATCHQ
# Number of nodes (CHANGE THIS if needed)
#PBS -l walltime=$wallclock_limit,nodes=1:ppn=24
# output file base name
#PBS -N test_dr
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V
# End of options

# Make sure core dumps are created
ulimit -c unlimited

if [ -n "\$PBS_JOBID" ]; then    #batch job
    export JOBID=\`echo \${PBS_JOBID} | cut -f1 -d'.'\`
    initdir=\${PBS_O_WORKDIR}
fi

if [ "\$PBS_ENVIRONMENT" = "PBS_BATCH" ]; then
    interactive=false
else
    interactive=true
fi

if [ -z "$CAM_RBOPTIONS" ]; then
    export CAM_RBOPTIONS="run_and_build"
fi

##omp threads
export CAM_THREADS=1
export CAM_RESTART_THREADS=2

##mpi tasks
export CAM_TASKS=$CAM_TASKS
export CAM_RESTART_TASKS=$CAM_RESTART_TASKS

export P4_GLOBMEMSIZE=500000000
source /usr/share/Modules/init/sh
module purge

export LAPACK_LIBDIR=/usr/lib64

if [ "\$CAM_FC" = "INTEL" ]; then
    module load compiler/intel/14.0.2
    export CFG_STRING=" -cc mpicc -fc_type intel -fc mpif90 -cppdefs -DNO_MPI2 -cppdefs -DNO_MPIMOD "
    export INC_NETCDF=\${NETCDF_PATH}/include
    export LIB_NETCDF=\${NETCDF_PATH}/lib
    input_file="tests_pretag_hobart_nag"
    export CCSM_MACH="hobart_intel"
elif [ "\$CAM_FC" = "NAG" ]; then
    module load compiler/nag/6.2

    export CFG_STRING="-cc mpicc -fc mpif90 -fc_type nag "
    export INC_NETCDF=\${NETCDF_PATH}/include
    export LIB_NETCDF=\${NETCDF_PATH}/lib
    input_file="tests_pretag_hobart_nag"
    export CCSM_MACH="hobart_nag"
else
    module load compiler/pgi/18.1
    export CFG_STRING=" -cc mpicc -fc_type pgi -fc mpif90 -cppdefs -DNO_MPI2 -cppdefs -DNO_MPIMOD "
    export INC_NETCDF=\${NETCDF_PATH}/include
    export LIB_NETCDF=\${NETCDF_PATH}/lib
    input_file="tests_pretag_hobart_pgi"
    export CCSM_MACH="hobart_pgi"
fi
export MAKE_CMD="gmake --output-sync -j $gmake_j"
export MACH_WORKSPACE="$mach_workspace"
export CPRNC_EXE=/fs/cgd/csm/tools/bin/cprnc
export ADDREALKIND_EXE=/fs/cgd/csm/tools/addrealkind/addrealkind
dataroot="/fs/cgd/csm"
echo_arg="-e"

EOF

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

# To lower number of nodes required for regression testing on Hobart,
# run CIME test suites sequentially after standalone regression tests
submit_script_cime="${submit_script}"
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##izumi
    izu* | i[[:digit:]]* )
    submit_script="`pwd -P`/test_driver_izumi${cur_time}.sh"
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

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv
cat > ${submit_script} << EOF
#!/bin/sh
#

# Name of the queue (CHANGE THIS if needed)
#PBS -q $CAM_BATCHQ
# Number of nodes (CHANGE THIS if needed)
#PBS -l walltime=$wallclock_limit,nodes=1:ppn=24
# output file base name
#PBS -N test_dr
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V
# End of options

# Make sure core dumps are created
ulimit -c unlimited

if [ -n "\$PBS_JOBID" ]; then    #batch job
    export JOBID=\`echo \${PBS_JOBID} | cut -f1 -d'.'\`
    initdir=\${PBS_O_WORKDIR}
fi

if [ "\$PBS_ENVIRONMENT" = "PBS_BATCH" ]; then
    interactive=false
else
    interactive=true
fi

if [ -z "$CAM_RBOPTIONS" ]; then
    export CAM_RBOPTIONS="run_and_build"
fi

##omp threads
export CAM_THREADS=1
export CAM_RESTART_THREADS=2

##mpi tasks
export CAM_TASKS=$CAM_TASKS
export CAM_RESTART_TASKS=$CAM_RESTART_TASKS

export P4_GLOBMEMSIZE=500000000
source /usr/share/Modules/init/sh
module purge

export LAPACK_LIBDIR=/usr/lib64

if [ "\$CAM_FC" = "INTEL" ]; then
    module load compiler/intel/14.0.2
    export CFG_STRING=" -cc mpicc -fc_type intel -fc mpif90 -cppdefs -DNO_MPI2 -cppdefs -DNO_MPIMOD "
    export INC_NETCDF=\${NETCDF_PATH}/include
    export LIB_NETCDF=\${NETCDF_PATH}/lib
    input_file="tests_pretag_izumi_nag"
    export CCSM_MACH="izumi_intel"
elif [ "\$CAM_FC" = "NAG" ]; then
    module load compiler/nag/6.2-8.1.0

    export CFG_STRING="-cc mpicc -fc mpif90 -fc_type nag "
    export INC_NETCDF=\${NETCDF_PATH}/include
    export LIB_NETCDF=\${NETCDF_PATH}/lib
    input_file="tests_pretag_izumi_nag"
    export CCSM_MACH="izumi_nag"
else
    module load compiler/pgi/20.1
    export CFG_STRING=" -cc mpicc -fc_type pgi -fc mpif90 -cppdefs -DNO_MPI2 -cppdefs -DNO_MPIMOD "
    export INC_NETCDF=\${NETCDF_PATH}/include
    export LIB_NETCDF=\${NETCDF_PATH}/lib
    input_file="tests_pretag_izumi_pgi"
    export CCSM_MACH="izumi_pgi"
fi
export MAKE_CMD="gmake --output-sync -j $gmake_j"
export MACH_WORKSPACE="$mach_workspace"
export CPRNC_EXE=/fs/cgd/csm/tools/cime/tools/cprnc/cprnc
export ADDREALKIND_EXE=/fs/cgd/csm/tools/addrealkind/addrealkind
dataroot="/fs/cgd/csm"
echo_arg="-e"

EOF

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

# To lower number of nodes required for regression testing on izumi,
# run CIME test suites sequentially after standalone regression tests
submit_script_cime="${submit_script}"
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;

    ##leehill - always run with CAM_FC=PGI and -i
    le* )
    submit_script="`pwd -P`/test_driver_leehill_${cur_time}.sh"
    export PATH=/cluster/torque/bin:${PATH}

    # Default setting is 12 hr in the long queue; the short queue only
    # allows 1 hr runs.
    wallclock_limit="12:00:00"
    if [ -z "$CAM_BATCHQ" ]; then
        export CAM_BATCHQ="long"
    elif [[ "$CAM_BATCHQ" == short ]]; then
        wallclock_limit="1:00:00"
    fi

    if [ $gmake_j = 0 ]; then
        gmake_j=8
    fi

    mach_workspace="/scratch/cluster"

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv

cat > ${submit_script} << EOF
#!/bin/sh
#

# Make sure core dumps are created
ulimit -c unlimited

if [ -n "\$PBS_JOBID" ]; then    #batch job
    export JOBID=\`echo \${PBS_JOBID} | cut -f1 -d'.'\`
    initdir=\${PBS_O_WORKDIR}
fi

if [ "\$PBS_ENVIRONMENT" = "PBS_BATCH" ]; then
    interactive=false
else
    interactive=true
fi

echo "interactive = ${interactive}"

if [ -z "$CAM_RBOPTIONS" ]; then
    export CAM_RBOPTIONS="run_and_build"
fi

export INTEL=/usr/local/intel-cluster
export NAG=/usr/local/nag
export PGI=/usr/local/pgi
export LD_LIBRARY_PATH=\${PGI}/linux86-64/lib:/cluster/torque/lib:\${INTEL}/lib/intel64
echo \${LD_LIBRARY_PATH}
export P4_GLOBMEMSIZE=500000000

# Only PGI is supported on leehill at this time
#if [ "\$CAM_FC" = "PGI" ]; then
    export LAPACK_LIBDIR=\${PGI}/linux86-64/lib
    export INC_NETCDF=/usr/local/netcdf-pgi/include
    export LIB_NETCDF=/usr/local/netcdf-pgi/lib
    export PATH=\${PGI}/linux86-64/bin:\${PATH}:\${LIB_NETCDF}
    export LD_LIBRARY_PATH=\${PGI}/linux86-64/libso:\${LIB_NETCDF}:\${LD_LIBRARY_PATH}
    export CFG_STRING=" -cc pgcc -fc_type pgi -fc pgf90 "
    export input_file="tests_pretag_leehill"
    export CCSM_MACH="leehill_pgi"
#fi
export MAKE_CMD="gmake -j $gmake_j"
export MACH_WORKSPACE="$mach_workspace"
export CPRNC_EXE=/fs/cgd/csm/tools/cprnc_64/cprnc
export ADDREALKIND_EXE=/fs/cgd/csm/tools/addrealkind/addrealkind
dataroot="/fs/cgd/csm"
echo_arg="-e"

EOF

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^
    ;;


    * ) echo "ERROR: machine $hostname not currently supported"; exit 1 ;;
esac

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv writing to batch script vvvvvvvvvvvvvvvvvvv

for file in ${submit_script} ${submit_script_cb}
do
cat >> ${file} << EOF

##check if interactive job
if  \$interactive ; then

    echo "test_driver.sh: interactive run - setting JOBID to \$\$"
    export JOBID=\$\$
    if [ \$0 = "test_driver.sh" ]; then
	initdir="."
    else
	initdir=\${0%/*}
    fi
fi

##establish script dir and cam_root
if [ -f \${initdir}/test_driver.sh ]; then
    export CAM_SCRIPTDIR=\`cd \${initdir}; pwd \`
    if [ -d "\${CAM_SCRIPTDIR}/../../components" ]; then
        export CAM_ROOT=\`cd \${CAM_SCRIPTDIR}/../.. ; pwd \`
    else
        export CAM_ROOT=\`cd \${CAM_SCRIPTDIR}/../../../.. ; pwd \`
    fi
else
    if [ -n "\${CAM_ROOT}" ] && [ -f \${CAM_ROOT}/components/cam/test/system/test_driver.sh ]; then
	export CAM_SCRIPTDIR=\`cd \${CAM_ROOT}/components/cam/test/system; pwd \`
    else
	if [ -n "\${CAM_ROOT}"  -a  -f "\${CAM_ROOT}/test/system/test_driver.sh" ]; then
            export CAM_SCRIPTDIR=\`cd \${CAM_ROOT}/test/system; pwd \`
        else
            echo "ERROR: unable to determine script directory "
	    echo "       if initiating batch job from directory other than the one containing test_driver.sh, "
	    echo "       you must set the environment variable CAM_ROOT to the full path of directory containing "
            echo "       <components>. "
	    exit 3
        fi
    fi
fi

##output files
cam_log=\${initdir}/td.\${JOBID}.log
if [ -f \$cam_log ]; then
    rm \$cam_log
fi
cam_status=\${initdir}/td.\${JOBID}.status
if [ -f \$cam_status ]; then
    rm \$cam_status
fi

##setup test work directory
if [ -z "\$CAM_TESTDIR" ]; then
    export CAM_TESTDIR=$mach_workspace/\$LOGNAME/test-driver.\${JOBID}
    if [ -d \$CAM_TESTDIR ]; then
        rm -rf \$CAM_TESTDIR
    fi
fi
if [ ! -d \$CAM_TESTDIR ]; then
    mkdir -p \$CAM_TESTDIR
    if [ \$? -ne 0 ]; then
	echo "ERROR: unable to create work directory \$CAM_TESTDIR"
	exit 4
    fi
fi

##set our own environment vars
export CSMDATA=\${dataroot}/inputdata
export MPI_TYPE_MAX=100000

##process other env vars possibly coming in
if [ -z "\$CAM_RETAIN_FILES" ]; then
    export CAM_RETAIN_FILES=FALSE
fi
if [ -n "\${CAM_INPUT_TESTS}" ]; then
    input_file=\$CAM_INPUT_TESTS
else
    input_file=\${CAM_SCRIPTDIR}/\${input_file}
fi

if [ ! -f \${input_file} ]; then
    echo "ERROR: unable to locate input file \${input_file}"
    exit 5
fi

if  \$interactive ; then
    echo "reading tests from \${input_file}"
else
    echo "reading tests from \${input_file}" >> \${cam_log}
fi

num_tests=\`wc -w < \${input_file}\`
echo "STATUS OF CAM TESTING UNDER JOB \${JOBID};  scheduled to run \$num_tests tests from:" >> \${cam_status}
echo "\$input_file" >> \${cam_status}
echo "" >> \${cam_status}
echo "Start testing at \$(date)" >> \${cam_status}
echo "On node \$(hostname)" >> \${cam_status}
uptime >> \${cam_status}
echo "" >> \${cam_status}
if  \$interactive ; then
    echo "see \${cam_log} for more detailed output" >> \${cam_status}
fi
echo "" >> \${cam_status}

test_list=""
while read input_line; do
    test_list="\${test_list}\${input_line} "
done < \${input_file}

##initialize flags, counter
skipped_tests="NO"
pending_tests="NO"
count=0

##loop through the tests of input file
for test_id in \${test_list}; do
    count=\`expr \$count + 1\`
    while [ \${#count} -lt 3 ]; do
        count="0\${count}"
    done

    master_line=\`grep \$test_id \${CAM_SCRIPTDIR}/input_tests_master\`
    status_out=""
    for arg in \${master_line}; do
        status_out="\${status_out}\${arg} "
    done


    test_cmd=\${status_out#* }

    status_out="\${count} \${status_out}"

    if  \$interactive ; then
        echo ""
        echo "************************************************************"
        echo "\${status_out}"
        echo "************************************************************"
    else
        echo "" >> \${cam_log}
        echo "************************************************************"\
            >> \${cam_log}
        echo "\$status_out" >> \${cam_log}
        echo "************************************************************"\
            >> \${cam_log}
    fi

    while [ \${#status_out} -lt 95 ]; do
        status_out="\${status_out}."
    done

    echo \$echo_arg "\$status_out\c" >> \${cam_status}



    if  \$interactive ; then
        \${CAM_SCRIPTDIR}/\${test_cmd} \$CAM_RBOPTIONS
        rc=\$?
    else
        \${CAM_SCRIPTDIR}/\${test_cmd} \$CAM_RBOPTIONS >> \${cam_log} 2>&1
        rc=\$?
    fi
    if [ \$rc -eq 0 ]; then
        if [ \${CAM_RBOPTIONS} = "build_only" ]; then
          echo "BUILT at \$(date)" >> \${cam_status}
        else
          echo "PASS at \$(date)" >> \${cam_status}
        fi
    elif [ \$rc -eq 255 ]; then
        echo "SKIPPED* at \$(date)" >> \${cam_status}
        skipped_tests="YES"
    elif [ \$rc -eq 254 ]; then
        echo "PENDING** at \$(date)" >> \${cam_status}
        pending_tests="YES"
    else
        echo "FAIL! rc= \$rc at \$(date)" >> \${cam_status}
        if  \$interactive ; then
	    if [ "\$CAM_SOFF" != "FALSE" ]; then
		echo "stopping on first failure"
		echo "stopping on first failure" >> \${cam_status}
		exit 6
	    fi
	else
	    if [ "\$CAM_SOFF" == "TRUE" ]; then
		echo "stopping on first failure" >> \${cam_status}
		echo "stopping on first failure" >> \${cam_log}
		exit 6
	    fi
	fi
    fi
done


echo "end of input" >> \${cam_status}
if  \$interactive ; then
    echo "end of input"
else
    echo "end of input" >> \${cam_log}
fi

if [ "\$skipped_tests" = "YES" ]; then
    echo "*  please verify that any skipped tests are not required of your cam commit" >> \${cam_status}
fi
if [ "\$pending_tests" = "YES" ]; then
    echo "** tests that are pending must be checked manually for a successful completion" >> \${cam_status}
    if  \$interactive ; then
	echo "   see the test's output in \${cam_log} " >> \${cam_status}
	echo "   for the location of test results" >> \${cam_status}
    fi
fi
EOF

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ writing to batch script ^^^^^^^^^^^^^^^^^^^

done

## Make sure we have a place to store test summaries
if [ ! -f "${SUMMARY_FILE}" ]; then
  touch "${SUMMARY_FILE}"
fi

##vvvvvvvvvvvvvvvvvvvvvvvvvvvvv add summary to output or email vvvvvvvvvvvvvvv

cat >> ${submit_script} << EOF
banner="========================================"
subj="CAM regression test summary from \$CCSM_MACH"
np="Num PASS = \$( grep PASS \${cam_status} | wc -l )"
nf="Num FAIL = \$( grep FAIL \${cam_status} | wc -l )"
js="$CAM_FC job started at ${start_date}"
je="$CAM_FC job finished at \$( date --iso-8601=seconds )"
echo "${banner}" | tee -a ${SUMMARY_FILE}
echo "\${js}" | tee -a ${SUMMARY_FILE}
echo "" | tee -a ${SUMMARY_FILE}
echo "\${np}" | tee -a ${SUMMARY_FILE}
echo "\${nf}" | tee -a ${SUMMARY_FILE}
if [ "${nf}" != "Num FAIL = 0" ]; then
  grep FAIL \${cam_status} | tee -a ${SUMMARY_FILE}
fi
echo "" | tee -a ${SUMMARY_FILE}
echo "\${je}" | tee -a ${SUMMARY_FILE}
echo "" | tee -a ${SUMMARY_FILE}
EOF

if $cam_email_summary; then
  cat >> ${submit_script} << EOF
echo -e "\${js}\n\n\${np}\n\${nf}\n\n${je}" | mail -s "\${subj}" ${EMAIL}
EOF
fi
if [ "${submit_script}" != "${submit_script_cime}" ]; then
  echo "exit 0" >> ${submit_script}
fi

# If setting up a separate "configure and build" script then add command to
# submit the run script after builds are done.
if [ -n "${submit_script_cb}" ]; then

    case $hostname in
        # cheyenne
        chey* | r* )
            batch_queue_submit='qsub '
	    ;;
        *)
            echo "ERROR: machine $hostname not currently supported for batch builds"
            exit 1
            ;;
    esac

cat >> ${submit_script_cb} << EOF
echo "$batch_queue_submit ${submit_script}" >> \${cam_log}
$batch_queue_submit ${submit_script} >> \${cam_log} 2>&1
exit 0
EOF

fi

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ add summary to output or email ^^^^^^^^^^^^^^^

for file in ${submit_script} ${submit_script_cb}
do
  chmod a+x $file
done

if $interactive ; then
    if [ ${submit_script_cb} ]; then
      ${submit_script_cb}
    else
      ${submit_script}
    fi
    exit 0
fi

if ! $force ; then
    echo ""
    echo "**********************"
    echo "$submit_script has been created and will be submitted to the batch queue..."
    echo "(ret) to continue, (a) to abort"
    read ans
    case $ans in
	[aA]* )
	echo "aborting...type ./test_driver.sh -h for help message"
	exit 0
	;;
    esac
fi

##vvvvvvvvvvvvvvvvvvvvvv start CAM aux test suite vvvvvvvvvvvvvvvvvvvvvvvvvvvv

cesm_test_mach=""
comp=""
if [ "${hostname:0:4}" == "chey" ]; then
  cesm_test_mach="cheyenne"
fi
if [ "${hostname:0:6}" == "hobart" ]; then
  cesm_test_mach="hobart"
fi
if [ "${hostname:0:5}" == "izumi" ]; then
  cesm_test_mach="izumi"
fi
if [ -n "${CAM_FC}" ]; then
  comp="_${CAM_FC,,}"
fi

if [ "${cesm_test_suite}" != "none" -a -n "${cesm_test_mach}" ]; then
  if [ "${hostname:0:6}" != "hobart" ]; then
    module load python
  fi
  if [ "${hostname:0:5}" != "izumi" ]; then
    module load python
  fi

  for cesm_test in ${cesm_test_suite}; do
    testargs="--xml-category ${cesm_test} --xml-machine ${cesm_test_mach}"

    if [ -n "${use_existing}" ]; then
      test_id="${use_existing}"
    else
      idstr="`date '+%Y%m%d%H%M%S'`"
      test_id=${cesm_test}${comp}"_"${idstr}
    fi
    currdir="`pwd -P`"
    logfile="${currdir}/${test_id}.log"
    tdir="$( cd $( dirname $0 ); pwd -P )"
    trial_dir="$( dirname $( dirname $( dirname $( dirname ${tdir} ) ) ) )"
    if [ -d "${trial_dir}/cime/scripts" ]; then
      root_dir=$trial_dir
    else
      root_dir="$( dirname $( dirname ${tdir} ) )"
    fi

    script_dir="${root_dir}/cime/scripts"
    if [ ! -d "${script_dir}" ]; then
      echo "ERROR: CIME scripts dir not found at ${script_dir}"
      exit 1
    fi
    if [ ! -x "${script_dir}/create_test" ]; then
      echo "ERROR: create_test script dir not found in ${script_dir}"
      exit 1
    fi

    ##setup CESM work directory
    cesm_testdir=$mach_workspace/$LOGNAME/$test_id

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
    testargs="${testargs} --queue ${CAM_BATCHQ} --test-root ${cesm_testdir} --output-root ${cesm_testdir}"
    if [ -n "${CAM_ACCOUNT}" ]; then
      testargs="${testargs} --project ${CAM_ACCOUNT}"
    fi
    testargs="${testargs} --test-id ${test_id}"
    if [ -n "${BL_TESTDIR}" ]; then
      testargs="${testargs} --compare ${BL_TESTDIR} "
    fi
    if [ -n "${use_existing}" ]; then
      testargs="${testargs} --use-existing"
    fi
    if [ -n "${archive_dir}" ]; then
      testargs="${testargs} --generate ${archive_dir}"
    fi

    echo ""
    echo "CESM test results will be in: ${cesm_testdir}" | tee ${logfile}
    echo "Running ./create_test ${testargs}"             | tee -a ${logfile}

    if [ "${hostname:0:2}" == "ch" ]; then
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

  done
fi

##^^^^^^^^^^^^^^^^^^^^^^ start CAM aux test suite ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

if $run_cam_regression; then

  echo ""
  echo "submitting cam regression tests..."

  case $hostname in
    ##cheyenne
    ch* | r* )
	qsub ${submit_script_cb}
	;;

    ##hobart
    hob* | h[[:digit:]]* )
#        qsub ${submit_script}
    echo "************************************************************************************************************"
    echo "****************   IMPORTANT *******************************************************************************"
    echo "************************************************************************************************************"
    echo "hobart testing is disabled - use izumi or find this message in test_driver.sh and uncomment the qsub command"
    echo "************************************************************************************************************"
    echo "************************************************************************************************************"
        ;;

    ##izumi
    izu* | i[[:digit:]]* )
        qsub ${submit_script}
        ;;

    ##leehill
    le* )
        ${submit_script}
        ;;

  esac

fi

exit 0
