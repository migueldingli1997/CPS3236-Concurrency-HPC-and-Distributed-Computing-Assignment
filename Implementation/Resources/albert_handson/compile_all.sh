f1="Job_1_MPI/"
f2="Job_2_OMP/"
f3="Job_3_Hybrid/"

#------------------------

cd $f1
mpic++ -std=c++11 -o 1_nbody_mpi 1_nbody_mpi.cpp
cd ../

#------------------------

cd $f2
g++ -std=c++11 -o 2_nbody_omp 2_nbody_omp.cpp -fopenmp
cd ../

#------------------------

cd $f3
mpic++ -std=c++11 -o 3_nbody_hybrid 3_nbody_hybrid.cpp -fopenmp
cd ../

#------------------------
