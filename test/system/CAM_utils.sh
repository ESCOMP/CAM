#!/bin/sh 
#
# define functions used by the test scripts

get_run_mode()
{
    # determine the CAM run mode based on arguments in the configure options file
    # arg $1 is name of config options file

    if [ $# -ne 1 ]; then
        echo "get_run_mode(): incorrect number of input arguments"
        exit 1
    fi

    if [ ! -f ${CAM_SCRIPTDIR}/config_files/$1 ]; then
        echo "get_run_mode(): configure options file ${CAM_SCRIPTDIR}/config_files/$1 not found"
        exit 2
    fi

    # search config options file for parallelization info
    spmd=1
    if grep -ic NOSPMD ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
        spmd=0
    fi

    smp=1
    if grep -ic NOSMP ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
        smp=0
    fi

    if [[ "$spmd" -eq "0" && "$smp" -eq "0" ]]; then
        # serial
        run_mode="serial"
    elif [[ "$spmd" -eq "0" && "$smp" -eq "1" ]]; then
        # openmp only
        run_mode="omp"
    elif [[ "$spmd" -eq "1" && "$smp" -eq "0" ]]; then
        # mpi only
        run_mode="mpi"
    else
        # hybrid
        run_mode="hybrid"
    fi

    echo -n $run_mode
}
