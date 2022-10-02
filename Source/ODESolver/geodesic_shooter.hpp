#ifndef GEODESIC_SHOOTER_HPP
#define GEODESIC_SHOOTER_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "ODESolverCore.hpp"
#include "Rk4vec.hpp"
#include "TensorAlgebra.hpp"

template <typename data_t> class geodesic_shooter
{
  public:
    void shoot(Vec3 center, double shift = 1 / 5., int numberofgeodesics = 10,
               bool set_geodesic_null = true, const double time_end = 150.0,
               const double time_start = 0.0, const double dt = 0.1,
               const double epsabs = 1e-6, const double epsrel = 1e-6,
               const double hstart = 1e-6, const int nmax = 10000);

    void single_shot(double y[], const int index = 0,
                     const double time_end = 150.0,
                     const double time_start = 0.0, const double dt = 0.1,
                     const double epsabs = 1e-6, const double epsrel = 1e-6,
                     const double hstart = 1e-6, const int nmax = 10000);
};

#include "geodesic_shooter.impl.hpp"

#endif
