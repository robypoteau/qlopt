#ifndef REGULARIZATION_H
#define REGULARIZATION_H
#include <misc.h>
#include <math.h>
#include <numerical_integration.h>
#include <thesis_functions.h>
#include <gsl/gsl_min.h>

struct curve_params {
    mat A;
    vec P;
    double O;
};

struct ymsmt_params {
    string *fhandle;
    vec times;
	mat A;
	vec P;
	mat msmt;
	vec uNot;
    vec yNot;
};

void regparamexp1a(soln_env *env, vec u, double gamma, int brk);
void regparamexp1b(soln_env *env, vec u, int brk);
void reg_guess_plots(soln_env *env, vec u, vec u_guess, double gamma, int brk);
vec reg_guess(soln_env *env, vec u_guess, double gamma);
vec reg_guess2(soln_env *env, vec u_guess);
double func(double a, mat A, vec P, double O, mat B, vec u_guess, vec uNot);
double findLambda(mat A, vec P, double O, mat B, vec u_guess, vec uNot);
vec reg1(soln_env *env, double gamma);
vec reg2(soln_env *env, vec u_guess);
double curvature (double x, void * params);
double alpha(mat A, vec P, vec O);
void curvature_test_plots(soln_env *env, int brk);
void curvetest(mat A, vec P, double O);
void alpha_plot(functype fhandle, void *params);
void ymsmt_test_plots(soln_env *env, double gamma, int brk);
void alpha_plot(functype fhandle, void *params);
double ymsmt_plot_function(double a, void * params);
vec dtregs(double gamma, vec u, mat A, vec P );
double findGamma(mat A, vec P, vec uNot, vec u);
vec reg_plot_and_finda(soln_env *env, vec u);
vec g_reg_plot_and_finda(soln_env *env, vec u, vec ug);
double g_findGamma(mat A, vec P, vec uNot, vec u, vec ug);

#endif
