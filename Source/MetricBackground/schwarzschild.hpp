#ifndef SCHWARZSCHILD_HPP
#define SCHWARZSCHILD_HPP

#include "DimensionDefinitions.hpp"
#include "Rk4vec.hpp"
#include "tensor.hpp"

using namespace std;

class Black_Hole
{

  private:
    static tensor<2, double> get_metric(double M, double x, double y, double z);

    static tensor<3, double> get_metric_deriv(double M, double x, double y,
                                              double z);

    // Calculate inverse, only true for this metric, not a general 4x4 inversion
    static tensor<2, double> calculate_spatial_inverse(tensor<2, double> g);

    static tensor<3, double> get_chris(tensor<2, double> g_UU,
                                       tensor<3, double> dg);

  public:
    static Vec3 eval_diff_eqn(double t, Vec3 v);

    double calculate_norm(Vec3 v, double M = 1);
    // Change the vt component to set norm to any value
    Vec3 set_norm(Vec3 v, double norm_val = 0, double M = 1);
};

#endif
