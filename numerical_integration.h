#ifndef NUMERICAL_INTEGRATION_H
#define NUMERICAL_INTEGRATION_H

#include "misc.h"
#include "spline.h"

mat rungekutta4(sys fhandle, const vec& time, const vec& u, const vec& yNot);
mat qLinearRungeKutta4(sys fhandle, const vec& time, const vec& u, const vec& yNot, const mat& xNminus);
double simpson(const vec& t, const vec& x);
mat qlinear(sys fhandle, const double& t, const vec& x, const vec& u, thesis::spline* Xn, int n);

//mat rungekutta4(mat (*fhandle)(const double&, const vec&, const vec&, const mat&, const vec&), const vec& time, const vec& u, const vec& yNot, const mat& xNminus);

#endif