#ifndef ISOMETRIC_HPP
#define ISOMETRIC_HPP

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "DimensionDefinitions.hpp"
#include "Rk4vec.hpp"
#include "TensorAlgebra.hpp"
#include "tensor.hpp"

class Black_Hole_isometric
{

  public:
    static tensor<2, double> get_metric(double M, double x, double y, double z);

    static tensor<3, double> get_metric_deriv(double M, double x, double y,
                                              double z);

    static tensor<3, double> get_chris(tensor<2, double> g_UU,
                                       tensor<3, double> dg);

    static int eval_diff_eqn(double t, const double y[], double f[],
                             void *params);

};

#endif
