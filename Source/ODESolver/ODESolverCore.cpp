#include "ODESolverCore.hpp"

void ODE_Solver(gsl_odeiv2_system sys, double y[], const double time_start,
                const double time_end, const int NumberOutputs, const double hstart,
                const double epsabs, const double epsrel)
{

    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rk8pd, hstart, epsabs, epsrel);

    int i;
    double t = time_start;

    for (i = 1; i <= NumberOutputs; i++)
    {
        double ti = i * time_end /(double)NumberOutputs;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

        if (status != GSL_SUCCESS)
        {
            printf("error, return value=%d\n", status);
            break;
        }

        printf("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }

    gsl_odeiv2_driver_free(d);
}
