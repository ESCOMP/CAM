#!/bin/sh
#
### Job name
#
#PBS -N scam_run

### Declare job non-rerunable
#
#PBS -r n

### Output files - sort to top of directory.
#
#PBS -e scam_run.err
#PBS -o scam_run.log

# Mail to user
#
#PBS -m ae

### Queue name (short, medium, long, verylong)
#
#PBS -q medium
#
# Number of nodes, number of processors
#
# nodes = physical host
# ppn   = processors per node (i.e., number of cores)
#
#PBS -l nodes=1:ppn=48

#
# This job's working directory
#
echo `date`
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

# May be necessary for some OMP jobs.
#
#export KPM_STACKSIZE=50m

echo "Environment:"
echo "--------------"
echo ""

# Print out some job information for debugging.
#
echo Running $PROGNAME on host `hostname`
echo Time is `date`
echo Directory is `pwd`

# Configure the run environment.
#
module load compiler/intel/default

#Until I found out otherwise use dumb command
../bld/cesm.exe

# Submit like this:            
#/usr/local/torque/bin/qsub ens_run.sh



# Lots of stuff I don't understand ... ...
#------------------------------------------------
# Convert the host file to use IB
#
#/cluster/bin/make_ib_hosts.sh

# Get the number of procs by counting the nodes file,
# which was generated from the #PBS -l line above.
#
#NPROCS=`wc -l < $PBS_NODEFILE`

#echo "Node File:"
#echo "----------"
#cat  "$PBS_NODEFILE"
#echo ""

# Run the parallel MPI executable 
#
#echo "`date` mpiexec - Start"

#mpiexec -v -np $NPROCS ../bld/cesm.exe

#echo ""
#echo "`date` MPIRUN - END"          


exit 0


