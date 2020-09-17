#ifndef ODESOLVERCORE_HPP
#define ODESOLVERCORE_HPP

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

int ODE_Solver(gsl_odeiv2_system sys, double y[], const double time_start = 0,
               const double time_end = 100,
               const double hstart = 1e-6, const double epsabs = 1e-6,
               const double epsrel = 0, const int nmax = 10000);

#endif
