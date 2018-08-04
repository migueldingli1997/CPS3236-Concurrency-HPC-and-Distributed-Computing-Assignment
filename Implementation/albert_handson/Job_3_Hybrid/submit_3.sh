declare -a numPart=() # Number of particles e.g .(64 1024 4096 16384)
declare -a numProc=() # Processors per node e.g. (12 6 4 2)
declare -a numNode=() # Number of nodes     e.g. (4 3 2 1)
iters=4

mkdir -p Output/1_2
mkdir -p Output/1_4
mkdir -p Output/1_6
mkdir -p Output/1_12

mkdir -p Output/2_2
mkdir -p Output/2_4
mkdir -p Output/2_6
mkdir -p Output/2_12

mkdir -p Output/3_2
mkdir -p Output/3_4
mkdir -p Output/3_6
mkdir -p Output/3_12

mkdir -p Output/4_2
mkdir -p Output/4_4
mkdir -p Output/4_6
mkdir -p Output/4_12

################################################################################

for P in "${numPart[@]}"
do
	for N2 in "${numProc[@]}"
	do
		for N1 in "${numNode[@]}"
    	do
  		    echo "Submitting $P $N1 $N2 jobs"
	        for i in `seq 1 $iters`
	        do
	            sleep .5
	            qsub -pe openmpi_"$N2"x1 "$(($N1 * $N2))" -o Output/"$N1"_"$N2" -e Output/"$N1"_"$N2" job_3.sh $N1 $N2 $P
	        done
	    done
	done
done

################################################################################
