declare -a numPart=() # Number of particles e.g. (64 1024 4096 16384)
declare -a numNode=() # Number of MPI procs e.g. (4 3 2 1)
iters=4

mkdir -p Output/1
mkdir -p Output/2
mkdir -p Output/3
mkdir -p Output/4

################################################################################

for P in "${numPart[@]}"
do
	for N in "${numNode[@]}"
	do
	    echo "Submitting $P $N jobs"
	    for i in `seq 1 $iters`
	    do
	        sleep .5
	        qsub -pe openmpi $N -o Output/$N -e Output/$N job_1.sh $N $P
	    done
	done
done

################################################################################
