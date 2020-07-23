#ifndef MISC_H
#define MISC_H

//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <nvector/nvector_serial.h>  // access to serial N_Vector

#define NV_CONTENT_S(v) ( (N_VectorContent_Serial)(v->content) )
#define NV_OWN_DATA_S(v) ( NV_CONTENT_S(v)->own_data )
#define NV_DATA_S(v) ( NV_CONTENT_S(v)->data )
#define NV_LENGTH_S(v) ( NV_CONTENT_S(v)->length )
#define Ith(v,i) ( NV_DATA_S(v)[i] )

//#include <mpreal.h>
//#include <eigen3/unsupported/Eigen/MPRealSupport>

using namespace Eigen;
using namespace std;
//using mpfr::mpreal;
//#define BITS 128

typedef Matrix<double, Dynamic, 1> vec;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> mat;
typedef vector< double > vector_type;
//typedef boost::numeric::ublas::matrix< double > matrix_type;

typedef vec (*odefunction)(
    const double& t,
    const vec& x,
    const vec& u,
    const vec& control);
#endif
