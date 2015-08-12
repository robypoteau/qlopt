#ifndef LATEX_OUTPUT_H
#define LATEX_OUTPUT_H

#include "misc.h"
#include "dbg.h"

void latexOutput(const mat& xn, const vec& u, int p, string buf);
void longlatexOutput(const mat& otpt);
void shortlatexOutput(const mat& otpt);
vec colWiseStdDev(const mat& M);
vec colWiseMean(const mat& M);
void tableheader(int n);
void tablefooter();
#endif