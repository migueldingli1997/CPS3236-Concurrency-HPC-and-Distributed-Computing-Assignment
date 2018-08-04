#!/bin/sh
#
# Your job name
#$ -N HYBRID_JOB
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

export OMP_NUM_THREADS=$2

# Run your application 
mpirun -npernode 1 --bind-to none 3_nbody_hybrid Output/"$1_$2" ../NBodyInput/input_$3.txt 1000 0.1 20.0 No > output_$3_$1_$2_`date +%s%N`.txt
