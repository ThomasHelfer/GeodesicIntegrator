#include "maximalslicing.hpp"

// Mass of the black hole
#define M 1

tensor<2, double> Black_Hole_maximal_slicing_isometric::get_metric(double x,
                                                                   double y,
                                                                   double z,
                                                                   double t)
{
    /*Definition used from https://arxiv.org/pdf/gr-qc/0701037.pdf*/

    tensor<2, double> g;
    double beta_u[3];
    double beta[3];
    double gamma[3][3];
    double normal[3];
    int i, j;
    double betasquared = 0;

    const double R = sqrt(x * x + y * y + z * z);
    const double r = R * (1.0 - M / R - M * M / (2 * R * R));

    /*    const double r = (2.0*R + M + sqrt(4.0*R*R+4.0*M*R+3*M*M))/4.0 *
           pow(((4.0 + 3.0*sqrt(2))*(2.0*R - 3*M ) )/(8.0*R + 6.0*M +
       3*sqrt( 8.0 * R*R + 8.0 * M * R + 6.0 * M*M )),1.0/sqrt(2.0) ) ;
    */
    normal[0] = x / R;
    normal[1] = y / R;
    normal[2] = z / R;

    const double psi = sqrt(R / r);

    const double alpha =
        sqrt(1.0 - 2.0 * M / R + 27. * M * M * M * M / (16.0 * R * R * R * R));

    for (i = 0; i < 3; i++)
    {
        beta[i] = 0;
        for (j = 0; j < 3; j++)
        {
            gamma[i][j] = 0;
        }
    }

    for (i = 0; i < 3; i++)
        gamma[i][i] = psi * psi * psi * psi;

    for (i = 0; i < 3; i++)
        beta_u[i] = 3.0 * sqrt(3.0) * M * M / 4.0 * r / (R * R * R) * normal[i];

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            beta[i] += gamma[i][j] * beta_u[j];
        }
    }

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            betasquared += gamma[i][j] * beta_u[i] * beta_u[j];
        }
    }

    for (i = 0; i < 3; i++)
    {
        g[i][3] = beta[i];
        g[3][i] = beta[i];
    }

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            g[i][j] = gamma[i][j];
        }
    }

    g[3][3] = -(alpha * alpha - betasquared);

    return g;
}

tensor<3, double>
Black_Hole_maximal_slicing_isometric::get_metric_deriv(double x, double y,
                                                       double z, double t)
{

    tensor<3, double> dg;
    tensor<2, double> g;
    tensor<2, double> g_dx;
    tensor<2, double> g_dy;
    tensor<2, double> g_dz;
    tensor<2, double> g_dt;
    double h = 1e-4;

    g = get_metric(x, y, z, t);
    g_dx = get_metric(x - h, y, z, t);
    g_dy = get_metric(x, y - h, z, t);
    g_dz = get_metric(x, y, z - h, t);
    g_dt = get_metric(x, y, z, t - h);

    FOR2(i, j)
    {
        dg[i][j][0] = (g[i][j] - g_dx[i][j]) / h;
        dg[i][j][1] = (g[i][j] - g_dy[i][j]) / h;
        dg[i][j][2] = (g[i][j] - g_dz[i][j]) / h;
        dg[i][j][3] = (g[i][j] - g_dt[i][j]) / h;
    }
    return dg;
}

int Black_Hole_maximal_slicing_isometric::eval_diff_eqn(double t,
                                                        const double y[],
                                                        double f[],
                                                        void *params)
{

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

    g = get_metric(v.x, v.y, v.z, v.t);

    g_UU = TensorAlgebra::compute_inverse(g);

    dg = get_metric_deriv(v.x, v.y, v.z, v.t);

    //=========================

    chris_ULL = TensorAlgebra::get_chris(g_UU, dg);

    //=========================

    FOR1(i)
    {
        FOR2(k, l) { ddx[i] += -chris_ULL[i][k][l] * dx[k] * dx[l]; }
    }

    Vec3 out(dx[0], dx[1], dx[2], dx[3], ddx[0], ddx[1], ddx[2], ddx[3]);

    out.write_to_array(f);

    return GSL_SUCCESS;
}
