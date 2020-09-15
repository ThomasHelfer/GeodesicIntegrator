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
#include "isometric.hpp"
#include "maximalslicing.hpp"

int main(void)
{
    // ==========================================
    // ========== Shooting some test geod========
    // ==========================================

    // Setting up inital data

    const double center_x = 0.0;
    const double center_y = 2;
    const double center_z = 0.0;
    const double start_time = 0.0;
    const double velocity_x = 0.0;
    const double velocity_y = 1.0;
    const double velocity_z = 0.0;

    const double shift_y = -0.1;
    const int numberofgeodesics = 20;

    const double lapse = -1.0;
    const bool null_geodesic = true;
    const double end_time = 150;
    const double dt = 0.01;

    const Vec3 initial_data(center_x, center_y, center_z, start_time,
                            velocity_x, velocity_y, velocity_z, lapse);

    geodesic_shooter<Black_Hole_maximal_slicing_isometric> pewpew;

    pewpew.shoot(initial_data, shift_y, numberofgeodesics, null_geodesic, end_time, start_time,
                 dt);

    return 0;
}
