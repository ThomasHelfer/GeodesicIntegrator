#ifndef GEODESIC_SHOOTER_HPP
#define GEODESIC_SHOOTER_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

#include "Rk4vec.hpp"

template<typename data_t>
class geodesic_shooter{
        public:
        void shoot(Vec3 center, double shift = 1/5., int shoot = 10 , bool set_null = true,  
		   const double TIME_MAXIMUM = 150.0, const double T_START = 0.0, const double DT = 0.1);
};

#endif
