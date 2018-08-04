SRC="1_nbody_mpi.cpp"   # source code file
EXEC="1_nbody"          # executable

OUTF="Output_1" # output folder
INPF=$1         # input file (e.g. NBodyInput/input_64.txt)
ITER=1000       # number of iterations
DELT=0.01       # delta_t
GTRM=20.0       # g_term
OUTPUT="yes"    # whether to output particle positions

mpic++ -std=c++11 -o $EXEC $SRC     # compile
rm -r $OUTF                         # clear previous output, if necessary
mkdir $OUTF                         # create output folder, if necessary

mpirun -np 4 $EXEC $OUTF $INPF $ITER $DELT $GTRM $OUTPUT    # run
rm $EXEC                                                    # remove executable
