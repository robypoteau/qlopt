#ifndef NOISE_REMOVAL_H
#define NOISE_REMOVAL_H
// External headers
#include <eigen3/Eigen/QR>

// My headers
#include <misc.h>

namespace thesis {
    vec lsNoiseRemoval(const vec& y, double noise);
    double regs(const vec& y, double noise);
}

#endif
