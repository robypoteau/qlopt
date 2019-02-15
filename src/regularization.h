#ifndef REGULARIZATION_H
#define REGULARIZATION_H

#include <misc.h>
#include <odesolver.h>
#include <odeWrapper.h>
#include <qlopt.h>

namespace thesis{

    double findAlpha(mat A, vec P, double Om, double Nd);
}

#endif
