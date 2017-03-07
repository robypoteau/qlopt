#ifndef REGULARIZATION_H
#define REGULARIZATION_H
#include <misc.h>
#include <math.h>
#include <numerical_integration.h>
#include <thesis_functions.h>

void regparamexp1a(soln_env *env, vec u, int brk);
void regparamexp1b(soln_env *env, vec u, int brk);
void reg_guess_plots(soln_env *env, vec u, vec u_guess, int brk);
vec reg_guess(soln_env *env, vec u_guess, double gamma);
vec reg_guess2(soln_env *env, vec u_guess);
double func(double a, mat A, vec P, double O, mat B, vec u_guess, vec uNot);
double findLambda(mat A, vec P, double O, mat B, vec u_guess, vec uNot);

#endif
