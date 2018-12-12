#!/bin/tcsh
#
## LSF batch script to run an MPI application
#
#BSUB -P PROJECTNUMBER                         # project code
#BSUB -W 24:00                                 # wall-clock time (hrs:mins)
#BSUB -n 1                                     # number of tasks in job         
#BSUB -J          Gen_ERAI_fv09_001.04         # job name
#BSUB -o ./OUTPUT/Gen_ERAI_fv09_001.04.%J.out  # output file name in which %J is replaced by the job ID
#BSUB -e ./OUTPUT/Gen_ERAI_fv09_001.04.%J.err  # error file name in which %J is replaced by the job ID
#BSUB -q geyser                                # queue

#run the executable
module load ncl
./Gen_ERAI_fv09_001.01.csh  >& ./TMP/TMP_001.01/DUMP
