#ifndef SCHWARZSCHILD_HPP
#define SCHWARZSCHILD_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Rk4vec.hpp"

template <typename data_t> class render_black_hole {
public:
  void picture(int *red, int *green, int *blue, Vec3 center, double max_x,
               double max_y, const double alpha = 0, const int start_ind = 0,
               int end_ind = (H * H), const double TIME_MAXIMUM = 150.0,
               const double DT = 0.1, const double T_START = 0);

  void render_circle(int *red, int *green, int *blue, double max_x,
                     double max_y, const int start_ind = 0,
                     int end_ind = (H * H), double Radius = 2,
                     double thickness = 0.1);

  void render(int *red, int *green, int *blue, std::string file_name);
};

#endif
