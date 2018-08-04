//----------------------------------------------------------------------------------------------
//	Filename:	nbody.cpp
//	Author:		Keith Bugeja
//----------------------------------------------------------------------------------------------
//  CPS3236 assignment for academic year 2017/2018:
//	Sample naive [O(n^2)] implementation for the n-Body problem.
//----------------------------------------------------------------------------------------------

#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>

#include "vector2.h"

#include "omp.h"
#include "mpi.h"
MPI_Datatype MPI_PARTICLE;  // main particle datatype
MPI_Datatype MPI_REDUCEDP;  // particle datatype w/o positions
int numtasks, rank;         // MPI number of tasks and task rank
int currIndex;              // index to keep track of current bodies subset

/*
 * Constant definitions for field dimensions, and particle masses
 */
const int fieldWidth = 1000;
const int fieldHalfWidth = fieldWidth >> 1;
const int fieldHeight = 1000;
const int fieldHalfHeight = fieldHeight >> 1;

const float minBodyMass = 2.5f;
const float maxBodyMassVariance = 5.f;

/*
 * Particle structure
 */
struct Particle
{
    Vector2 Position;
    Vector2 Velocity;
    float	Mass;

    Particle(void)
        : Position( ((float)rand()) / RAND_MAX * fieldWidth - fieldHalfWidth,
                    ((float)rand()) / RAND_MAX * fieldHeight - fieldHalfHeight)
        , Velocity( 0.f, 0.f )
        , Mass ( ((float)rand()) / RAND_MAX * maxBodyMassVariance + minBodyMass )
    { }
};

/*
 * Compute forces of particles exerted on one another
 */
void ComputeForces(std::vector<Particle> &p_bodies1, float p_gravitationalTerm, float p_deltaT, std::vector<Particle> &p_bodies2)
{
    #pragma omp parallel for default(shared)
    for (size_t j = 0; j < p_bodies1.size(); ++j)
    {
        Particle &p1 = p_bodies1[j];
        float f1 = 0.0f, f2 = 0.0f;

        #pragma omp parallel for reduction(+:f1,f2) default(shared)
        for (size_t k = 0; k < p_bodies2.size(); ++k)
        {
            if (currIndex == rank && k == j) continue;

            Particle &p2 = p_bodies2[k];

            // Compute direction vector
            Vector2 direction = p2.Position - p1.Position;

            // Limit distance term to avoid singularities
            const float distance = std::max<float>( 0.5f * (p2.Mass + p1.Mass), fabs(direction.Length()) );

            // Accumulate force
            Vector2 temp = direction / (distance * distance * distance) * p2.Mass;
            f1 += temp.Element[0];
            f2 += temp.Element[1];
        }

        // Compute acceleration for body 
        Vector2 acceleration = Vector2(f1, f2) * p_gravitationalTerm;

        // Integrate velocity (m/s)
        p1.Velocity += acceleration * p_deltaT;
    }
}

/*
 * Update particle positions
 */
void MoveBodies(std::vector<Particle> &p_bodies, float p_deltaT)
{
    #pragma omp parallel for
    for (size_t j = 0; j < p_bodies.size(); ++j)
    {
        p_bodies[j].Position += p_bodies[j].Velocity * p_deltaT;
    }
}

/*
 * Commit particle masses and positions to file in CSV format
 */
void PersistPositions(const std::string &p_strFilename, std::vector<Particle> &p_bodies, bool append)
{
    std::cout << "Writing to file: " << p_strFilename << std::endl;

    std::ofstream output(p_strFilename.c_str(), (append ? std::ofstream::app : std::ofstream::out));

    if (output.is_open())
    {
        for (int j = 0; j < p_bodies.size(); j++)
        {
            output << p_bodies[j].Mass << ", " <<
                   p_bodies[j].Position.Element[0] << ", " <<
                   p_bodies[j].Position.Element[1] << std::endl;
        }

        output.close();
    }
    else
        std::cerr << "Unable to persist data to file:" << p_strFilename << std::endl;

}

int main(int argc, char **argv)
{
    // Time start
    const double start_time = MPI_Wtime();

    // Check number of arguments
    const int expNoOfArg = 6;
    if (argc != expNoOfArg + 1) {
        std::cout << "Expected: ./program OUT_FOLDER IN_FILE MAX_ITERS DELTA_T G_TERM OUTPUT_FILES(yes/no)" << std::endl;
        return 1;
    }

    // Vectors of particles
    std::vector<Particle> bodies;       // main bodies
    std::vector<Particle> subBodies;    // main bodies subset
    std::vector<Particle> recvBuffer;   // buffer for subset to receive
    std::vector<Particle> sendBuffer;   // buffer for subset to send
    
    // Main Variables
    const int maxIteration = std::stoi(argv[3]);            // example: 1000
    const float deltaT = std::stof(argv[4]);                // example: 0.01
    const float gTerm = std::stof(argv[5]);                 // example: 20.0
    const bool outputFiles = strcmp(argv[6], "yes") == 0;   // example: yes
    int particleCount;                                      // depends on input
    
    // Communication request variables
    MPI_Request sendReq;
    MPI_Request recvReq;

    // File Input/Output
    std::ifstream fileInput(argv[2]);
    std::stringstream fileOutput;

    // Initialise MPI environment
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // MPI_PARTICLE structured type and commit
    MPI_Datatype oldTypes1[1] = {MPI_FLOAT};
    int blockCounts1[1] = {5};
    MPI_Aint offsets1[1] = { offsetof(Particle, Position) };
    MPI_Type_create_struct(1, blockCounts1, offsets1, oldTypes1, &MPI_PARTICLE);
    MPI_Type_commit(&MPI_PARTICLE);
    
    // MPI_REDUCEDP structured type and commit
    MPI_Datatype oldTypes2[2] = {MPI_FLOAT, MPI_FLOAT};
    int blockCounts2[2] = {2, 1};
    MPI_Aint offsets2[2] = { offsetof(Particle, Position), offsetof(Particle, Mass) };
    MPI_Type_create_struct(2, blockCounts2, offsets2, oldTypes2, &MPI_REDUCEDP);
    MPI_Type_commit(&MPI_REDUCEDP);
    
    // Setting particle count and bodies vector
    if (rank == 0) {
        float mass, posX, posY;
        char comma1, comma2;
        
        while (fileInput >> mass >> comma1 >> posX >> comma2 >> posY)
        {
            Particle toAdd = Particle();
            toAdd.Mass = mass;
            toAdd.Position = Vector2(posX, posY);
            bodies.push_back(toAdd);
        }
        fileInput.close();
        particleCount = (int) bodies.size();
    }
    MPI_Bcast(&particleCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Setting mainSize and displs arrays
    int mainSize[numtasks]; // sizes of main subsets held by each node (a.k.a send counts)
    int displs[numtasks];   // displacements
    const int baseAmount = particleCount / numtasks; // minimum particles per MPI task
    const int remainder = particleCount % numtasks; // remaining particles to be distributed
    int displacement = 0;
    for (int i = 0; i < numtasks; i++) {
        mainSize[i] = baseAmount + (i < remainder ? 1 : 0);
        displs[i] = displacement;
        displacement += mainSize[i];
    }

    // Distribute the subsets of bodies
    subBodies.resize(mainSize[rank]);
    MPI_Scatterv(&bodies[0], mainSize, displs, MPI_PARTICLE, &subBodies[0], mainSize[rank], MPI_PARTICLE, 0, MPI_COMM_WORLD);

    // Rank of source and destination
    const int src = (rank + 1) % numtasks;              // source node
    const int dst = (rank - 1 + numtasks) % numtasks;   // destination node

    // Main iterations loop
    for (int iteration = 0; iteration < maxIteration; ++iteration)
    {
        // Initial subset and its index
        recvBuffer.resize(mainSize[rank]);  // resize to fit subBodies
        recvBuffer = subBodies;             // first subset held will be subBodies
        currIndex = rank;                   // index of subset currently held
        
        // Send-receive-compute-persist loop
        int i;
        for (i = 0; i < numtasks - 1; i++)
        {
            const int tempIndex = (currIndex + 1) % numtasks;   // compute next subset's index
            sendBuffer.resize(mainSize[currIndex]);             // sendBuffer must be able to hold current subset
            sendBuffer = recvBuffer;                            // current subset copied to sendBuffer
            recvBuffer.resize(mainSize[tempIndex]);             // recvBuffer must be able to hold next subset
            
            // Non-blocking exchange of bodies subsets
            MPI_Isend(&sendBuffer[0], mainSize[currIndex], MPI_REDUCEDP, dst, 0, MPI_COMM_WORLD, &sendReq);
            MPI_Irecv(&recvBuffer[0], mainSize[tempIndex], MPI_REDUCEDP, src, 0, MPI_COMM_WORLD, &recvReq);
            
            // Compute and persist
            ComputeForces(subBodies, gTerm, deltaT, sendBuffer);
            if (rank == 0 && iteration != 0 && outputFiles) {
                PersistPositions(fileOutput.str(), sendBuffer, i != 0);
            }
            
            // Finalize
            currIndex = tempIndex;
            MPI_Wait(&sendReq, MPI_STATUS_IGNORE);
            MPI_Wait(&recvReq, MPI_STATUS_IGNORE);
        }
        
        // Compute-persist (no need to send/receive last subset)
        ComputeForces(subBodies, gTerm, deltaT, recvBuffer);
        if (rank == 0 && iteration != 0 && outputFiles) {
            PersistPositions(fileOutput.str(), recvBuffer, i != 0);
        }
        
        // Move bodies in main subset
        MoveBodies(subBodies, deltaT);
        
        // Set file output stream for next iteration (or just print iteration)
        if (rank == 0) {
            if (outputFiles) {
                fileOutput.str(std::string());
                fileOutput << argv[1] << "/nbody_" << iteration << ".txt";
            } else {
                printf("%d\n", iteration);
            }
        }
    }
    
    // Gather bodies for final persist
    if (outputFiles) {
        MPI_Gatherv(&subBodies[0], mainSize[rank], MPI_PARTICLE, &bodies[0], mainSize, displs, MPI_PARTICLE, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            PersistPositions(fileOutput.str(), bodies, false);
        }
    }

    // Time end
    double end_time = MPI_Wtime() - start_time;
    double max_end_time;
    MPI_Reduce(&end_time, &max_end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Execution time: %.2f ms\n", max_end_time * 1e3);
    }

    // Free derived data type and close down MPI environment
    MPI_Type_free(&MPI_PARTICLE);
    MPI_Finalize();

    return 0;
}
