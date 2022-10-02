#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
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
    // Position of first geodesic
    const double center_x = 15;
    const double center_y = -12.5;
    const double center_z = 0.0;
    const double start_time = 0.0;
    // Speed of geodesic
    const double velocity_x = -1.0;
    const double velocity_y = 0.0;
    const double velocity_z = 0.0;
    const double lapse = -1.0;
    // The norm will be changed such that the geodesic is null
    const bool null_geodesic = true;

    // All the next geodesics are shifted in y direction by
    // the amount shift_y
    const double shift_y = 0.25;
    const int numberofgeodesics = 100;

    // Cutoff time
    const double end_time = 150;
    // This is only telling how often we are writing out data (timesteps are
    // chosen dynamically )
    const double dt = 0.5;

    // Precision of arguments for the solver
    const double epsabs = 1e-6;
    const double epsrel = 1e-6;
    const double hstart = 1e-6;
    // Number of maximal steps
    const int nmax = 1000;

    const Vec3 initial_data(center_x, center_y, center_z, start_time,
                            velocity_x, velocity_y, velocity_z, lapse, epsabs,
                            epsrel, hstart, nmax);

    geodesic_shooter<Black_Hole> pewpew;

    pewpew.shoot(initial_data, shift_y, numberofgeodesics, null_geodesic,
                 end_time, start_time, dt);

    return 0;
}
