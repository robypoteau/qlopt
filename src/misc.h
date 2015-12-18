#ifndef MISC_H
#define MISC_H

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
	
typedef Matrix<double, Dynamic, 1> vec;
typedef Matrix<double, Dynamic, Dynamic> mat;
typedef mat (*sys)(const double& t, const vec& x, const vec& u);
typedef struct{
	sys ode;
	vec *time;
	vec *initial_cond;
	vec *initial_params;
	vec *actual_params;
	mat *nth_soln;
	mat *measurements;
} soln_env;
#endif