#ifndef MISC_H
#define MISC_H


#include <mpreal.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MPRealSupport>

using namespace mpfr;
using namespace Eigen;
using namespace std;

//typedef mpfr::mpreal mpreal; 	
typedef Matrix<mpreal, Dynamic, 1> vec;
typedef Matrix<mpreal, Dynamic, Dynamic> mat;

/*
#include <Eigen/Dense>
#include <math.h>

using namespace Eigen;
using namespace std;

typedef double mpreal; 	
typedef Matrix<mpreal, Dynamic, 1> vec;
typedef Matrix<mpreal, Dynamic, Dynamic> mat;*/
#endif