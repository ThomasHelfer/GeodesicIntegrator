#ifndef TENSORALGEBRA_HPP_
#define TENSORALGEBRA_HPP_

#include "tensor.hpp"

namespace TensorAlgebra
{
/// Computes the inverse of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
template <class data_t>
inline tensor<2, data_t> compute_inverse(const tensor<2, data_t, 4> &matrix)
{
    tensor<2, data_t> h_UU;

    const double A2323 =
        matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2];
    const double A1323 =
        matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1];
    const double A1223 =
        matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1];
    const double A0323 =
        matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0];
    const double A0223 =
        matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0];
    const double A0123 =
        matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0];
    const double A2313 =
        matrix[1][2] * matrix[3][3] - matrix[1][3] * matrix[3][2];
    const double A1313 =
        matrix[1][1] * matrix[3][3] - matrix[1][3] * matrix[3][1];
    const double A1213 =
        matrix[1][1] * matrix[3][2] - matrix[1][2] * matrix[3][1];
    const double A2312 =
        matrix[1][2] * matrix[2][3] - matrix[1][3] * matrix[2][2];
    const double A1312 =
        matrix[1][1] * matrix[2][3] - matrix[1][3] * matrix[2][1];
    const double A1212 =
        matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1];
    const double A0313 =
        matrix[1][0] * matrix[3][3] - matrix[1][3] * matrix[3][0];
    const double A0213 =
        matrix[1][0] * matrix[3][2] - matrix[1][2] * matrix[3][0];
    const double A0312 =
        matrix[1][0] * matrix[2][3] - matrix[1][3] * matrix[2][0];
    const double A0212 =
        matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0];
    const double A0113 =
        matrix[1][0] * matrix[3][1] - matrix[1][1] * matrix[3][0];
    const double A0112 =
        matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0];

    const double det =
        matrix[0][0] * (matrix[1][1] * A2323 - matrix[1][2] * A1323 +
                        matrix[1][3] * A1223) -
        matrix[0][1] * (matrix[1][0] * A2323 - matrix[1][2] * A0323 +
                        matrix[1][3] * A0223) +
        matrix[0][2] * (matrix[1][0] * A1323 - matrix[1][1] * A0323 +
                        matrix[1][3] * A0123) -
        matrix[0][3] * (matrix[1][0] * A1223 - matrix[1][1] * A0223 +
                        matrix[1][2] * A0123);
    const double invdet = 1 / det;

    h_UU[0][0] = invdet * (matrix[1][1] * A2323 - matrix[1][2] * A1323 +
                           matrix[1][3] * A1223);
    h_UU[0][1] = invdet * -(matrix[0][1] * A2323 - matrix[0][2] * A1323 +
                            matrix[0][3] * A1223);
    h_UU[0][2] = invdet * (matrix[0][1] * A2313 - matrix[0][2] * A1313 +
                           matrix[0][3] * A1213);
    h_UU[0][3] = invdet * -(matrix[0][1] * A2312 - matrix[0][2] * A1312 +
                            matrix[0][3] * A1212);
    h_UU[1][0] = invdet * -(matrix[1][0] * A2323 - matrix[1][2] * A0323 +
                            matrix[1][3] * A0223);
    h_UU[1][1] = invdet * (matrix[0][0] * A2323 - matrix[0][2] * A0323 +
                           matrix[0][3] * A0223);
    h_UU[1][2] = invdet * -(matrix[0][0] * A2313 - matrix[0][2] * A0313 +
                            matrix[0][3] * A0213);
    h_UU[1][3] = invdet * (matrix[0][0] * A2312 - matrix[0][2] * A0312 +
                           matrix[0][3] * A0212);
    h_UU[2][0] = invdet * (matrix[1][0] * A1323 - matrix[1][1] * A0323 +
                           matrix[1][3] * A0123);
    h_UU[2][1] = invdet * -(matrix[0][0] * A1323 - matrix[0][1] * A0323 +
                            matrix[0][3] * A0123);
    h_UU[2][2] = invdet * (matrix[0][0] * A1313 - matrix[0][1] * A0313 +
                           matrix[0][3] * A0113);
    h_UU[2][3] = invdet * -(matrix[0][0] * A1312 - matrix[0][1] * A0312 +
                            matrix[0][3] * A0112);
    h_UU[3][0] = invdet * -(matrix[1][0] * A1223 - matrix[1][1] * A0223 +
                            matrix[1][2] * A0123);
    h_UU[3][1] = invdet * (matrix[0][0] * A1223 - matrix[0][1] * A0223 +
                           matrix[0][2] * A0123);
    h_UU[3][2] = invdet * -(matrix[0][0] * A1213 - matrix[0][1] * A0213 +
                            matrix[0][2] * A0113);
    h_UU[3][3] = invdet * (matrix[0][0] * A1212 - matrix[0][1] * A0212 +
                           matrix[0][2] * A0112);

    return h_UU;
}

inline double calculate_norm( const double vx, const  double vy, const  double vz,double vt, const tensor<2,double> g)
{
    double norm = 0;
    tensor<1, double> dx = {vx, vy, vz, vt};

    FOR2(i, j) { norm += g[i][j] * dx[i] * dx[j]; }
    return norm;
}

// Change the vt component to set norm to any value
template<class data_t>
inline Vec3 set_norm(Vec3 v, const double norm_val, const double M=1)
{
    data_t metric;
    tensor<1, double> dx = {v.vx, v.vy, v.vz, v.vt};
    tensor<2, double> g = metric.get_metric(M, v.x, v.y, v.z);
    double norm = calculate_norm(v.vx,v.vy,v.vz,v.vt,g);

    double b = 0;
    for (int i = 0; i < 3; i++)
    {
        b += g[i][3] * dx[i] + g[3][i] * dx[i];
    }

    double discriminant = b * b;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            discriminant += -(4 * g[3][3] * g[i][j] * dx[i] * dx[j]);
        }
    }

    v.vt = (-b - sqrt(discriminant)) / (2.0 * g[3][3]);

    return v;
}

} // namespace TensorAlgebra

#endif
