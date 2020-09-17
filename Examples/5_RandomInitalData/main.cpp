#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <random>
#include <vector>

#include "DimensionDefinitions.hpp"
#include "geodesic_shooter.hpp"
#include "render.hpp"
#include "rk4.hpp"
#include "schwarzschild.hpp"
#include "tensor.hpp"

int main(void)
{
    // ==========================================
    // ========== Shooting some test geod========
    // ==========================================

    // Setting up inital data
    const int seed = 231;
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> velocity_rnd(0.1, 0.4);


    int numtasks, rank, size, size_per_task ;
    const int number_of_particles = 400;

    // ============= MPI INIT ================
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    size_per_task = (int)ceil(number_of_particles/ (double)numtasks);
    // ============= GEODESIC INIT ================
    const double x = 15;
    const double y = 0.0;
    const double z = 0.0;
    const double t = 0.0;
    const double vx = 0.0;
    const double vy = 0.0;
    const double vz = 0.0;
    const double vt = -1.00;

    const double start_time = 0;
    const double end_time = 10000;
    const double dt = 0.5;
    double intialdata[8];

    for(int index = size_per_task*rank ; index < size_per_task*(rank+1) ;index++){
        intialdata[0] = x;
        intialdata[1] = y;
        intialdata[2] = z;
        intialdata[3] = t;
        intialdata[4] = vx;
        intialdata[5] = velocity_rnd(generator);
        intialdata[6] = vz;
        intialdata[7] = -vt;

        geodesic_shooter<Black_Hole> pewpew;

        pewpew.single_shot(intialdata, index, end_time, start_time, dt);

    }

    MPI_Finalize();

    return 0;
}
