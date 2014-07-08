#ifndef NUMERICAL_INTEGRATION_H
#define NUMERICAL_INTEGRATION_H

#include <unsupported/Eigen/MPRealSupport>
#include "misc.h"
	
using namespace mpfr;
using namespace Eigen;
using namespace std;

mat rungekutta4(mat (*fhandle)(const mpreal&, const vec&, const vec&), const vec& time, const vec& u, const vec& yNot);
mat rungekutta4(mat (*fhandle)(const mpreal&, const vec&, const vec&, const mat&, const int&), const vec& time, const vec& u, const vec& yNot, const mat& xNminus);
mpreal simpson(const vec& t, const vec& x);

//mat rungekutta4(mat (*sys)(vec t, mat x, vec u, vec xNminus1), vec time, vec yNot);
#endif