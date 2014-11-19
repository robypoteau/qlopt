#ifndef MISC_H
#define MISC_H

#include <mpreal.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MPRealSupport>
#include <string>
	
using namespace mpfr;
using namespace Eigen;
using namespace std;
	
typedef Matrix<mpreal, Dynamic, 1> vec;
typedef Matrix<mpreal, Dynamic, Dynamic> mat;
#endif