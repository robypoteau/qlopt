#ifndef REGULARIZATION_H
#define REGULARIZATION_H

#include <misc.h>
#include <odesolver.h>
#include <odeWrapper.h>
#include <qlopt.h>

namespace thesis{
    double findAlpha1(const mat& A, const vec& P, const double Om, const double Nd);
    double findAlpha2(const mat& A, const vec& P, const double Om, const double Nd);
    double findAlpha3(const mat& A, const vec& P, const double Om, const vec& uguess, const double Nd);
    double findGamma(const mat& A, const vec& P, const vec& uNot, const vec& u);
    double findGamma2(const mat& A, const vec& P, const vec& uNot, const vec& uguess, const vec& u);
    double findAlpha5(const mat& A, const vec& P, const vec& u0, OdeWrapper ow, const vector<mat>& msmt, const vector<vec>& input, vec y0, vec ts, const vector<vector<thesis::spline>>& spl_pairs);
    //double findGamma(double initialGuess, void * params);
}

#endif
