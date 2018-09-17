#ifndef MISC_H
#define MISC_H

#include <eigen3/Eigen/Dense>
//#include <mpreal.h>
//#include <eigen3/unsupported/Eigen/MPRealSupport>

using namespace Eigen;
using namespace std;
//using mpfr::mpreal;

//#define BITS 128

typedef Matrix<double, Dynamic, 1> vec;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> mat;

//typedef Matrix<mpreal, Dynamic, 1> mp_vec;
//typedef Matrix<mpreal, Dynamic, Dynamic, RowMajor> mp_mat;

//typedef double (*functype)(double t, void *params);
//typedef mp_mat (*mp_sys)(const mpreal& t, const mp_vec& x, const mp_vec& u);

/*typedef struct{
	string *ode;
	vec *time;
	vec *initial_cond;
	vec *initial_params;
	vec *u_guess;
	mat *nth_soln;
	mat *measurements;
	double *lambda;

	//remove from struct, no longer needed
	//mp_mat *mp_nth_soln;
	//mp_mat *mp_measurements;
} soln_env;*/
#endif
