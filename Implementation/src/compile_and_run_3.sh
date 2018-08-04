SRC="3_nbody_hybrid.cpp"    # source code file
EXEC="3_nbody"              # executable

OUTF="Output_3" # output folder
INPF=$1         # input file (e.g. NBodyInput/input_64.txt)
ITER=1000       # number of iterations
DELT=0.01       # delta_t
GTRM=20.0       # g_term
OUTPUT="yes"    # whether to output particle positions

mpic++ -std=c++11 -o $EXEC $SRC -fopenmp    # compile
rm -r $OUTF                                 # clear previous output, if necessary
mkdir $OUTF                                 # create output folder, if necessary

export OMP_NUM_THREADS=2                                   # number of OpenMP threads
mpirun -np 2 $EXEC $OUTF $INPF $ITER $DELT $GTRM $OUTPUT   # run
rm $EXEC                                                   # remove executable
