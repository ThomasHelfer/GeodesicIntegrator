#include "isometric.hpp"

tensor<2, double> Black_Hole_isometric::get_metric(double M, double x, double y,
                                                   double z)
{

    tensor<2, double> g;
    double r = sqrt(x * x + y * y + z * z);
    const double psi = 1.0 + M / (2.0 * r);
    FOR2(i, j) g[i][j] = 0;
    for (int i = 0; i < 3; i++)
        g[i][i] = psi * psi * psi * psi;

    const double alpha = (1.0 - M / (2. * r)) / (1.0 + M / (2. * r));

    g[3][3] = -alpha * alpha;

    return g;
}

tensor<3, double> Black_Hole_isometric::get_metric_deriv(double M, double x,
                                                         double y, double z)
{

    tensor<3, double> dg;
    tensor<2, double> g;
    tensor<2, double> g_dx;
    tensor<2, double> g_dy;
    tensor<2, double> g_dz;
    double h = 1e-8;

    g = get_metric(M, x, y, z);
    g_dx = get_metric(M, x - h, y, z);
    g_dy = get_metric(M, x, y - h, z);
    g_dz = get_metric(M, x, y, z - h);

    FOR2(i, j)
    {
        dg[i][j][0] = (g[i][j] - g_dx[i][j]) / h;
        dg[i][j][1] = (g[i][j] - g_dy[i][j]) / h;
        dg[i][j][2] = (g[i][j] - g_dz[i][j]) / h;
        dg[i][j][3] = 0;
    }
    return dg;
}

tensor<3, double> Black_Hole_isometric::get_chris(tensor<2, double> g_UU,
                                                  tensor<3, double> dg)
{

    // Init
    tensor<3, double> chris_LLL; // Christoffel index low low low
    tensor<3, double> chris_ULL; // Christoffel index high low low

    FOR3(i, j, k)
    {
        chris_LLL[i][j][k] = 0;
        chris_ULL[i][j][k] = 0;
    }
    // Calculation of Christoffel symbols
    FOR3(i, j, k)
    {
        chris_LLL[i][j][k] = 0.5 * (dg[j][i][k] + dg[k][i][j] - dg[j][k][i]);
    }
    FOR3(i, j, k)
    {
        chris_ULL[i][j][k] = 0;
        FOR1(l) { chris_ULL[i][j][k] += g_UU[i][l] * chris_LLL[l][j][k]; }
    }

    return chris_ULL;
}

int Black_Hole_isometric::eval_diff_eqn(double t, const double y[], double f[],
                                        void *params)
{

    double M = 1;
    tensor<2, double> g; // Metix Index low low
    tensor<2, double> g_UU;
    tensor<3, double> dg;        //
    tensor<3, double> chris_ULL; // Christoffel index high low low
    tensor<2, double> jacobian = {};

    Vec3 v(y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);

    tensor<1, double> dx = {v.vx, v.vy, v.vz, v.vt};
    tensor<1, double> ddx;

    FOR1(i) { ddx[i] = 0; }
    FOR2(i, j)
    {
        g[i][j] = 0;
        g_UU[i][j] = 0;
        jacobian[i][j] = 0;
    }
    FOR3(i, j, k) { dg[i][j][k] = 0; }

    double rr2 = pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2);
    double rho2 = pow(v.x, 2) + pow(v.y, 2);

    double rr = sqrt(rr2);
    double rho = sqrt(rho2);
    // sinus(theta)
    double sintheta = rho / rr;
    double costheta = v.z / rr;
    // cos(phi)
    double cosphi = v.x / rho;
    // sin(phi)
    double sinphi = v.y / rho;
    double eps = 1e-7;

    // Freezing out geodesics that are too close to Horizon (Metric is singular
    // at horizon)
    /*    if (rr < 2 * M + eps)
        {
            FOR1(i)
            {
                ddx[i] = 0;
                dx[i] = 0;
            }
            Vec3 out(dx[0], dx[1], dx[2], dx[3], ddx[0], ddx[1], ddx[2],
       ddx[3]); return out;
        }
    */
    // ====================

    g = get_metric(M, v.x, v.y, v.z);

    g_UU = TensorAlgebra::compute_inverse(g);

    dg = get_metric_deriv(M, v.x, v.y, v.z);

    //=========================

    chris_ULL = get_chris(g_UU, dg);

    //=========================

    FOR1(i)
    {
        FOR2(k, l) { ddx[i] += -chris_ULL[i][k][l] * dx[k] * dx[l]; }
    }

    Vec3 out(dx[0], dx[1], dx[2], dx[3], ddx[0], ddx[1], ddx[2], ddx[3]);

    out.write_to_array(f);

    return GSL_SUCCESS;
}

double Black_Hole_isometric::calculate_norm(Vec3 v, double M)
{
    double norm = 0;
    tensor<1, double> dx = {v.vx, v.vy, v.vz, v.vt};
    tensor<2, double> g = get_metric(M, v.x, v.y, v.z);

    FOR2(i, j) { norm += g[i][j] * dx[i] * dx[j]; }
    return norm;
}
// Change the vt component to set norm to any value
Vec3 Black_Hole_isometric::set_norm(Vec3 v, double norm_val, double M)
{
    double norm = calculate_norm(v);
    tensor<1, double> dx = {v.vx, v.vy, v.vz, v.vt};
    tensor<2, double> g = get_metric(M, v.x, v.y, v.z);

    double tmp = norm - g[3][3] * dx[3] * dx[3];
    v.vt = sqrt(1.0 / g[3][3] * (norm_val - tmp));
    return v;
}
