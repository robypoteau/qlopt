#ifndef THESIS_FUNCTIONS_H
#define THESIS_FUNCTIONS_H

#include <misc.h>
#include <dbg.h>

#include <float.h>
#include <math.h>
#include <algorithm>    // std::max
#include <vector>
#include <stack>          // std::stack
#include <queue>

//#include <glpk.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

#include <numerical_integration.h>
#include <latex_output.h>
#include <nonlinear_odes.h>
#include <nonlinear_odes.h>
#include <spline.h>

typedef struct output{
	int fevals;
	int iterations;
	vec du;
} output;

mat findA(const vec& t, const mat& U, int m);
vec findP(const vec& t, const mat& U, const vec& dx, int m);
double findO(const vec& t, const vec& dx);
double findGamma(double initialGuess, void * params);
double innerProd(const vec& u1, const vec& u2, const vec& time);
output findActualParam(soln_env *env, int regs, const int numdivs);
mat reshape(const mat& U, int n, int m);
double norm(const mat& M);
mat inverse(const mat& M);
void vecToGslVec(const vec& v, gsl_vector *gslv);
void matToGslMat(const mat& m, gsl_matrix *gslm);
vec gslVecToVec(gsl_vector *gslv);
mat gslMatToMat(gsl_matrix *gslm);
bool allpositive(const vec& x);
double cond(const mat& A);
mat ichol(const mat& A);
mat corrMat(const mat& M);
double rcond(const mat& M);
double matnorm1(const mat& M);

//vec fnnls(const mat& A, const vec& b);
//vec dulp(const mat& A, const vec& b, const vec& u);
#endif
