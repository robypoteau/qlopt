#ifndef MISC_H
#define MISC_H

#include <Eigen/Dense>
#include <mpreal.h>
#include <unsupported/Eigen/MPRealSupport>

using namespace Eigen;
using namespace std;
using mpfr::mpreal;
	
typedef Matrix<double, Dynamic, 1> vec;
typedef Matrix<double, Dynamic, Dynamic> mat;

typedef Matrix<mpreal, Dynamic, 1> mp_vec;
typedef Matrix<mpreal, Dynamic, Dynamic> mp_mat;

typedef mat (*sys)(const double& t, const vec& x, const vec& u);
typedef mp_mat (*mp_sys)(const mpreal& t, const mp_vec& x, const mp_vec& u);

typedef struct{
	string *ode;
	vec *time;
	vec *initial_cond;
	vec *initial_params;
	vec *actual_params;
	mat *nth_soln;
	//mat *np_measurements;
	mat *measurements;
} soln_env;
#endif