#ifndef LATEX_OUTPUT_H
#define LATEX_OUTPUT_H

#include <misc.h>
#include <dbg.h>
#include <math.h>

void latexOutput(const mat& xn, const vec& u, int p, string buf);
void timelatexOutput(const vec& t, string buf, int n, int p);
void longlatexOutput(const mat& otpt);
void shortlatexOutput(const mat& otpt);
void R(double dt, const mat& otpt, int n);
void M(const mat& otpt, int n);
void shortNormalizedLatexOutput(const mat& M);
vec colWiseStdDev(const mat& M);
vec colWiseMean(const mat& M);
void tableheader(int n);
void tablefooter();
#endif