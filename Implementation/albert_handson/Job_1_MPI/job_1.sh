#!/bin/sh
#
# Your job name
#$ -N MPI_JOB
#
# Use current working directory
#$ -cwd
#
# Run job through bash shell
#$ -S /bin/bash
 
# If modules are needed, source modules environment:
. /etc/profile.d/modules.sh
 
# Add any modules you might require:
module add shared openmpi/gcc/64/1.8.8
 
# The following output will show in the output file
echo "Got $NSLOTS processors."

# Run your application
mpirun -np $1 1_nbody_mpi Output/$1 ../NBodyInput/input_$2.txt 1000 0.1 20.0 No > output_$2_$1_`date +%s%N`.txt
