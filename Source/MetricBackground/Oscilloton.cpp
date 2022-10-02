#include "Oscilloton.hpp"

double Oscilloton::a202[row_max][8];
double Oscilloton::c202[row_max][8];
double Oscilloton::m_omega;

tensor<2, double> Oscilloton::get_metric(double x, double y, double z, double t)
{

    tensor<2, double> jacobian;
    tensor<2, double> g;
    tensor<2, double> g_spher;

    FOR2(i, j)
    {
        g[i][j] = 0;
        g_spher[i][j] = 0;
        jacobian[i][j] = 0;
    }

    double rr2 = pow(x, 2) + pow(y, 2) + pow(z, 2);
    double rho2 = pow(x, 2) + pow(y, 2);

    double rr = sqrt(rr2);
    double rho = sqrt(rho2);
    // sinus(theta)
    double sintheta = rho / rr;
    double costheta = z / rr;
    // cos(phi)
    double cosphi = x / rho;
    // sin(phi)
    double sinphi = y / rho;

    double A = 0;
    double C = 0;

    if (rr < spacing * row_max)
    {
        for (int i = 0; i < 7; i++)
        {
            A += get_a202(rr, i + 1) * cos((2 * i) * m_omega * t);
            C += get_c202(rr, i + 1) * cos((2 * i) * m_omega * t);
        }
    }
    else
    {
        A = 1;
        C = 1;
    }

    g_spher[0][0] = A;
    // Define theta theta component
    g_spher[1][1] = rr2;
    // Define phi phi component
    g_spher[2][2] = rr2 * pow(sintheta, 2);
    // Devine tt component
    g_spher[3][3] = -A / C;

    jacobian[0][0] = x / rr;
    jacobian[1][0] = cosphi * z / rr2;
    jacobian[2][0] = -y / rho2;
    jacobian[0][1] = y / rr;
    jacobian[1][1] = sinphi * z / rr2;
    jacobian[2][1] = x / rho2;
    jacobian[0][2] = z / rr;
    jacobian[1][2] = -rho / rr2;
    jacobian[2][2] = 0.0;
    jacobian[3][3] = 1.0;

    // ====================

    FOR2(i, j)
    {
        FOR2(k, l)
        {
            g[i][j] += g_spher[k][l] * jacobian[k][i] * jacobian[l][j];
        }
    }

    return g;
}

tensor<3, double> Oscilloton::get_metric_deriv(double x, double y, double z,
                                               double t)
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

Oscilloton::Oscilloton()
{
    std::ifstream ifs("Oscilloton_data/general_a202.dat");
    std::ifstream ifs2("Oscilloton_data/general_c202.dat");

    for (int i = 0; i < row_max; ++i)
    {
        for (int j = 0; j < column_max; j++)
        {
            ifs >> a202[i][j];
            ifs2 >> c202[i][j];
        }
    }
    m_omega = sqrt(c202[row_max - 1][1]) / a202[row_max - 1][1];
}

int Oscilloton::eval_diff_eqn(double t, const double y[], double f[],
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

double Oscilloton::get_a202(double rr, int component)
{

    int indxL = static_cast<int>(floor(rr / spacing));
    int indxH = static_cast<int>(ceil(rr / spacing));

    double data_L = a202[indxL][component];
    double data_H = a202[indxH][component];

    double out = data_L + (rr / spacing - indxL) * (data_H - data_L);

    return out;
}

double Oscilloton::get_c202(double rr, int component)
{

    int indxL = static_cast<int>(floor(rr / spacing));
    int indxH = static_cast<int>(ceil(rr / spacing));

    double data_L = c202[indxL][component];
    double data_H = c202[indxH][component];

    double out = data_L + (rr / spacing - indxL) * (data_H - data_L);

    return out;
}
