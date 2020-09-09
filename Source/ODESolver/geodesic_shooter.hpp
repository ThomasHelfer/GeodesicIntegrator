#ifndef GEODESIC_SHOOTER_HPP
#define GEODESIC_SHOOTER_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Rk4vec.hpp"
#include "ODESolverCore.hpp"

template <typename data_t> class geodesic_shooter
{
    public:

    void shoot(Vec3 center, double shift = 1 / 5., int shoot = 10,
               bool set_geodesic_null = true, const double time_end = 150.0,
               const double time_start = 0.0, const double dt = 0.1);
};

#include "geodesic_shooter.impl.hpp"

#endif
