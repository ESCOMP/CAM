#!/bin/sh

# Set up run command for serial, openmp, mpi, or hybrid configurations.
# Input args:
# $1 - name of config file
# $2 - number of tasks
# $3 - number of threads

# utility functions
. $CAM_SCRIPTDIR/CAM_utils.sh

if [ $# -ne 3 ]; then
    echo "CAM_runcmnd.sh: incorrect number of input arguments"
    exit 1
fi

run_mode=`get_run_mode $1`

# set up threading via use of OMP_NUM_THREADS environment variable
if [ $run_mode = serial ]; then
    cmnd=""

elif [ $run_mode = omp ] || [ $run_mode = hybrid ]; then
    cmnd="env OMP_NUM_THREADS=$3 "

elif [ $run_mode = mpi ]; then
    cmnd="env OMP_NUM_THREADS=1 "
fi

# add MPI job launcher to command

if [ $run_mode = mpi ] || [ $run_mode = hybrid ]; then

    hostname=`hostname`
    case $hostname in

        # cheyenne
        ch* | r* )

	    cmnd="${cmnd} mpiexec_mpt -np $2 omplace -vv "
 	    # cmnd="${cmnd} ddt --connect mpiexec_mpt -np $ntasks omplace -vv "
            ;;

        # hobart and leehill
        hob* | h[[:digit:]]* | le* | izu* | i[[:digit:]]* )

            cmnd="${cmnd} mpiexec -n $2 "
            ;;

        * ) 
            echo "CAM_runcmnd.sh: unable to construct run command for unsupported machine $hostname "
            exit 3
            ;;
    esac

fi

#store command in temporary file for calling script to access
echo ${cmnd} > cam_run_command.txt
exit 0
