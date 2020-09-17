#include "ODESolverCore.hpp"

int ODE_Solver(gsl_odeiv2_system sys, double y[], const double time_start,
               const double time_end, const double hstart, const double epsabs,
               const double epsrel, const int nmax)
{

    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);

    int i;
    double t = time_start;

    gsl_odeiv2_driver_set_nmax(d, nmax);

    int status = gsl_odeiv2_driver_apply(d, &t, time_end, y);

    if (status != GSL_SUCCESS)
    {
        printf("error, return value=%d\n", status);
        return status;
    }

    gsl_odeiv2_driver_free(d);

    return 0;
}
