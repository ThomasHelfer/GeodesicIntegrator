#include "minkowski.hpp"



tensor<2, double, DIM> Minkowski::get_metric(double x, double y, double z,
                                                   double t,double theta)
{

    tensor<2, double, DIM> g;

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            {g[i][j] = 0;}


    for (int i = 0; i < DIM; i++)
            {g[i][i] = 1;}
        
    g[DIM-1][DIM-1] = -1;

    return g;
}

tensor<3, double,DIM> Minkowski::get_metric_deriv(double x, double y,
                                                         double z, double t,double theta)
{

    tensor<3, double,DIM> dg;
    tensor<2, double,DIM> g;
    tensor<2, double,DIM> g_dx;
    tensor<2, double,DIM> g_dy;
    tensor<2, double,DIM> g_dz;
    tensor<2, double,DIM> g_dt;
    tensor<2, double,DIM> g_dtheta;
    double h = 1e-4;

    g = get_metric(x, y, z, t,theta);
    g_dx = get_metric(x - h, y, z, t,theta);
    g_dy = get_metric(x, y - h, z, t,theta);
    g_dz = get_metric(x, y, z - h, t,theta);
    g_dt = get_metric(x, y, z, t - h,theta);
    g_dtheta = get_metric(x, y, z, t,theta-h);

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < DIM; k++)
    {
        dg[i][j][0] = (g[i][j] - g_dx[i][j]) / h;
        dg[i][j][1] = (g[i][j] - g_dy[i][j]) / h;
        dg[i][j][2] = (g[i][j] - g_dz[i][j]) / h;
        dg[i][j][3] = (g[i][j] - g_dt[i][j]) / h;
        dg[i][j][4] = (g[i][j] - g_dtheta[i][j]) / h;
    }
    return dg;
}

int Minkowski::eval_diff_eqn(double t, const double y[], double f[],
                                        void *params)
{

    tensor<2, double,DIM> g; // Metix Index low low
    tensor<2, double,DIM> g_UU;
    tensor<3, double,DIM> dg;        //
    tensor<3, double,DIM> chris_ULL; // Christoffel index high low low
    tensor<2, double,DIM> jacobian = {};

    Vec3 v(y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);

    tensor<1, double> dx = {v.vx, v.vy, v.vz, v.vt};
    tensor<1, double> ddx;
   
    for (int i = 0; i < DIM; i++) { ddx[i] = 0; }
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
    {
        g[i][j] = 0;
        g_UU[i][j] = 0;
        jacobian[i][j] = 0;
    }
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            for (int k = 0; k < DIM; k++)
                { dg[i][j][k] = 0; }
    
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
    
    g_UU = invert_a_matrix(g);
    
    dg = get_metric_deriv(v.x, v.y, v.z, v.t);
    
    //=========================

    
    chris_ULL = TensorAlgebra::get_chris<DIM>(g_UU, dg);
    
    //=========================
    
    FOR1(i)
    {
        FOR2(k, l) { ddx[i] += /*-chris_ULL[i][k][l]*/0 * dx[k] * dx[l]; }
    }

    Vec3 out(dx[0], dx[1], dx[2], dx[3], ddx[0], ddx[1], ddx[2], ddx[3]);

    out.write_to_array(f);
    
    return GSL_SUCCESS;
}

// This is a terrible hack 
tensor<2, double,DIM> Minkowski::invert_a_matrix( tensor<2, double,DIM> g){
    gsl_matrix *mat = gsl_matrix_alloc(DIM, DIM);
    gsl_matrix *invmat = gsl_matrix_alloc(DIM, DIM);


    tensor<2, double,DIM> ginv;

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            gsl_matrix_set(mat, i, j, g[i][j]);
    
    invmat = invert_a_matrix(mat);

    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < DIM; j++)
            ginv[i][j] = gsl_matrix_get(invmat,i,j);

    return ginv;
}

gsl_matrix* Minkowski::invert_a_matrix(gsl_matrix *matrix)
{
    gsl_permutation *p = gsl_permutation_alloc(DIM);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(DIM, DIM);
    gsl_linalg_LU_invert(matrix, p, inv);

    gsl_permutation_free(p);

    return inv;
}