In the Resources folder:

    - albert_handson: Indicates exactly the arrangement used while executing on Albert.

                    The folders Job_1_MPI, Job_2_OMP, and Job_3_Hybrid each contain
                    the source code, and the job.sh and submit.sh scripts for the
                    MPI, OpenMP, and Hybrid implementations, respectively. Each folder
                    also contains the vector2.h which is identical in each folder.

                    The folder NBodyInput contains the four provided input files.

                    The compile_all.sh script is the script used to compile all of
                    of the source code files in each folder into four executables.

    - src: Contains the source code files presented in the documentation. These are
                    identical to the ones found in albert_handson. For each source,
                    file, a compile_and_run script is  provided, which only requires
                    an input-file argument.

                    Example: ./compile_and_run_1.sh NBodyInput/input_1024.txt

    - Performance Measurement Data: Contains three excel sheets with the data obtained by
                    executing the three different solutions, as well as an averages for
                    each group of four readings.

    - Outputs for nbody_1024.txt: Contains the nbody_x.txt files obtained by running each
                    of the four solutions on a laptop.