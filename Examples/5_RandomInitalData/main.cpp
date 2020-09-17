#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
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

    const double x = 15;
    const double y = -12.5;
    const double z = 0.0;
    const double t = 0.0;
    const double vx = -1.0;
    const double vy = 0.0;
    const double vz = 0.0;
    const double vt = -1.0;

    const double start_time = 0;
    const bool null_geodesic = true;
    const double end_time = 150;
    const double dt = 0.5;

    std::mt19937 generator(123);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    double intialdata[8];
    intialdata[0] = x;
    intialdata[1] = y;
    intialdata[2] = z;
    intialdata[3] = t;
    intialdata[4] = vx;
    intialdata[5] = vy;
    intialdata[6] = vz;
    intialdata[7] = -vt;

    geodesic_shooter<Black_Hole> pewpew;

    pewpew.single_shot(intialdata, 0, end_time, start_time, dt);

    return 0;
}
