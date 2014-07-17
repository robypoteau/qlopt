#ifndef NONLINEAR_ODES_H
#define NONLINEAR_ODES_H

#include "misc.h"
#include "spline.h"

typedef mat (*odeFuncPtr)(const mpreal& , const vec&, const vec&);
typedef mat (*qLinFuncPtr)(const mpreal& , const vec&, const vec&, const mat&, const vec&);

mat lotka_volterra(const mpreal& t, const vec& x, const vec& u);
mat lotka_volterra_linearization(const mpreal& t, const vec& x, const vec& u, const mat& xn, const vec& time);

#endif

/*
double angiogenesis();*/