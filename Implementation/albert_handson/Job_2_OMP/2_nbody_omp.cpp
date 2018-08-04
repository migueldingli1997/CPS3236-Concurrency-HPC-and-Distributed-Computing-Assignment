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
void ComputeForces(std::vector<Particle> &p_bodies, float p_gravitationalTerm, float p_deltaT)
{
    #pragma omp parallel for default(shared)
    for (size_t j = 0; j < p_bodies.size(); ++j)
    {
        Particle &p1 = p_bodies[j];
        float f1 = 0.0f, f2 = 0.0f;

        #pragma omp parallel for reduction(+:f1,f2) default(shared)
        for (size_t k = 0; k < p_bodies.size(); ++k)
        {
            if (k == j) continue;

            Particle &p2 = p_bodies[k];

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
void PersistPositions(const std::string &p_strFilename, std::vector<Particle> &p_bodies)
{
    std::cout << "Writing to file: " << p_strFilename << std::endl;

    std::ofstream output(p_strFilename.c_str());

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
    const double start_time = omp_get_wtime();

    // Check number of arguments
    const int expNoOfArg = 6;
    if (argc != expNoOfArg + 1) {
        std::cout << "Expected: ./program OUT_FOLDER IN_FILE MAX_ITERS DELTA_T G_TERM OUTPUT_FILES(yes/no)" << std::endl;
        return 1;
    }

    // Main Particles
    std::vector<Particle> bodies;
    
    // Main Variables
    const int maxIteration = std::stoi(argv[3]);            // example: 1000
    const float deltaT = std::stof(argv[4]);                // example: 0.01
    const float gTerm = std::stof(argv[5]);                 // example: 20.0
    const bool outputFiles = strcmp(argv[6], "yes") == 0;   // example: yes

    // File Input/Output
    std::ifstream fileInput(argv[2]);
    std::stringstream fileOutput;

    // Setting particle count and bodies vector
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

    // Main iterations loop
    for (int iteration = 0; iteration < maxIteration; ++iteration)
    {
        ComputeForces(bodies, gTerm, deltaT);
        MoveBodies(bodies, deltaT);

        if (outputFiles) {
            fileOutput.str(std::string());
            fileOutput << argv[1] << "/nbody_" << iteration << ".txt";
            PersistPositions(fileOutput.str(), bodies);
        } else {
            printf("%d\n", iteration);
        }
    }

    // Time end
    printf("Execution time: %.2f ms\n", (omp_get_wtime() - start_time) * 1e3);

    return 0;
}
