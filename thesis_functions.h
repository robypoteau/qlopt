#ifndef THESIS_FUNCTIONS_H
#define THESIS_FUNCTIONS_H

#include "misc.h"

mat findA(const vec& t, const mat& U, int m);
vec findP(const vec& t, const mat& U, const vec& dx, int m);
mpreal innerProd(const vec& u1, const vec& u2, const vec& time);

#endif