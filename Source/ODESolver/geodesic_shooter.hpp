#ifndef GEODESIC_SHOOTER_HPP
#define GEODESIC_SHOOTER_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Rk4vec.hpp"

template <typename data_t> class geodesic_shooter
{
  public:
    void shoot(Vec3 center, double shift = 1 / 5., int shoot = 10,
               bool set_geodesic_null = true, const double TIME_MAXIMUM = 150.0,
               const double T_START = 0.0, const double DT = 0.1);
};

#include "geodesic_shooter.impl.hpp"

#endif
