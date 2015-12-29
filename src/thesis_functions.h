#ifndef THESIS_FUNCTIONS_H
#define THESIS_FUNCTIONS_H

#include <misc.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <numerical_integration.h>

mat findA(const vec& t, const mat& U, int m);
vec findP(const vec& t, const mat& U, const vec& dx, int m);
double innerProd(const vec& u1, const vec& u2, const vec& time);
vec findActualParam(soln_env *env, bool regs);
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
 vec dulp(const mat& A, const vec& b, const vec& u);
//vec regularization(soln_env *env, const vec& du);
#endif