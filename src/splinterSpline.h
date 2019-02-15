#ifndef SPLINTER_SPLINE_H
#define SPLINTER_SPLINE_H

#include <misc.h>

#include <iostream>

#include <SPLINTER/datatable.h>
#include <SPLINTER/bspline.h>
#include <SPLINTER/bsplinebuilder.h>

using std::cout;
using std::endl;

using namespace SPLINTER;

namespace thesis {
    vec splinterSpline(vec x, vec y, double alpha);
}

#endif
