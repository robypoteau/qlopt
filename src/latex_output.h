#ifndef LATEX_OUTPUT_H
#define LATEX_OUTPUT_H

#include <misc.h>
#include <limits>
#include <math.h>
#include <iostream>
#include <fstream>

void plotOutput(vec x, vec y);
void pythonplot(vec x, vec y);
void latexplot(vec x, vec y);
void convertVec(vec x);
void parameterOutput(const mat& uvals, const vec& alpha);
void tableheader(int n);
void tablefooter();
#endif
