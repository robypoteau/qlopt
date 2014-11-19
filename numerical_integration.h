#ifndef NUMERICAL_INTEGRATION_H
#define NUMERICAL_INTEGRATION_H

#include <unsupported/Eigen/MPRealSupport>
#include "misc.h"
#include "spline.h"
	
using namespace mpfr;
using namespace Eigen;
using namespace std;

typedef mat (*sys)(const mpreal& , const vec&, const vec&);

mat rungekutta4(mat (*fhandle)(const mpreal&, const vec&, const vec&), const vec& time, const vec& u, const vec& yNot);
mat rungekutta4(mat (*fhandle)(const mpreal&, const vec&, const vec&, const mat&, const vec&), const vec& time, const vec& u, const vec& yNot, const mat& xNminus);
mat qLinearRungeKutta4(sys, const vec& time, const vec& u, const vec& yNot, const mat& xNminus);
mpreal simpson(const vec& t, const vec& x);
mat qlinear(sys, const mpreal& t, const vec& x, const vec& u, const mat& xn, const vec& time);

#endif