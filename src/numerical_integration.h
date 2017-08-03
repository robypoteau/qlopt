#ifndef NUMERICAL_INTEGRATION_H
#define NUMERICAL_INTEGRATION_H

#include <misc.h>
#include <spline.h>
#include <nonlinear_odes.h>
#include <gsl/gsl_integration.h>
//#include <mp_nonlinear_odes.h>

mat rungekutta4(string fname, const vec& time, const vec& u, const vec& yNot);
//mp_mat mp_rungekutta4(string fname, const mp_vec& time, const mp_vec& u, const mp_vec& yNot);
mat qLinearRungeKutta4(string fname, const vec& time, const vec& u, const vec& yNot, const mat& xNminus);
double simpson(const vec& t, const vec& x);
mat qlinear(sys fhandle, const double& t, const vec& x, const vec& u, thesis::spline* Xn, int n);
double gsl_integration(const vec& t, const vec& x);
#endif
