#ifndef NONLINEAR_ODES_H_INCLUDED
#define NONLINEAR_ODES_H_INCLUDED

#include "misc.h"

typedef mat (*odeFuncPtr)(const mpreal& , const vec&, const vec&);
typedef mat (*qLinFuncPtr)(const mpreal& , const vec&, const vec&, const mat&, const int&);

mat lotka_volterra(const mpreal& time, const vec& x, const vec& u);
mat lotka_volterra_linearization(const mpreal& time, const vec& x, const vec& u, const mat& xn, const int& i);

#endif


/*
mat lotka_volterra(const vec& time, const vec& x, const vec& u);
double lotka_volterra_linearization(double* time, double* x, double* u, double* xn);

double angiogenesis();*/