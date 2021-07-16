#ifndef OSCILLOTON_HPP
#define OSCILLOTON_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "DimensionDefinitions.hpp"
#include "Rk4vec.hpp"
#include "TensorAlgebra.hpp"
#include "tensor.hpp"

// These are some definitions relevant for reading in Oscilloton data

#define column_max 8
#define row_max 2600
#define spacing 0.01

class Oscilloton
{
  public:
    static double a202[row_max][8]; // a_j
    static double c202[row_max][8]; // c_j
    static double m_omega;

    Oscilloton();

    static int eval_diff_eqn(double t, const double y[], double f[],
                             void *params);



    static double get_a202(double rr, int component);

    static double get_c202(double rr, int component);

    static tensor<2, double> get_metric(double M, double x, double y, double z,
                                        double t);

    static tensor<3, double> get_metric_deriv(double M, double x, double y,
                                              double z, double t);

};

#endif
