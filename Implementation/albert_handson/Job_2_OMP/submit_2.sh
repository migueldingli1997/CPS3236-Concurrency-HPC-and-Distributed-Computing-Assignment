declare -a numPart=() # Number of particles e.g. (64 1024 4096 16384)
declare -a numProc=() # Processes per node  e.g. (12 6 4 2 1)
iters=4

mkdir -p Output/1
mkdir -p Output/2
mkdir -p Output/4
mkdir -p Output/6
mkdir -p Output/12

################################################################################

for P in "${numPart[@]}"
do
	for N in "${numProc[@]}"
	do
        echo "Submitting $P $N jobs"
	    for i in `seq 1 $iters`
	    do
	        sleep .5
	        if [ $N -eq 1 ]
			then
				qsub -pe openmpi $N -o Output/$N -e Output/$N job_2.sh $N $P
			else
				qsub -pe openmpi_"$N"x1 $N -o Output/$N -e Output/$N job_2.sh $N $P
			fi
	    done
	done
done

################################################################################
