#ifndef TSQR_H
#define TSQR_H

#include <misc.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_linalg.h>

namespace thesis{
	class tallskinnyqr {
		private:
			size_t num_obs, order;
			
			gsl_multilarge_linear_workspace * w;
			gsl_matrix *gslA;
			gsl_vector *gslx, *gslb;
			
			double rnorm, snorm, rcond, lambda;
			
			gsl_vector *reg_param;
			gsl_vector *rho;
			gsl_vector *eta;


			void vecToGslVec(const vec& v, gsl_vector *gslv);
			void matToGslMat(const mat& m, gsl_matrix *gslm);
			vec gslVecToVec(gsl_vector *gslv);
			mat gslMatToMat(gsl_matrix *gslm);
			
		public:
			tallskinnyqr(){}
			~tallskinnyqr();
			tallskinnyqr(const size_t n, const size_t p, const size_t nlcurve);
			void init(const size_t n, const size_t p, const size_t nlcurve);
			void update(const mat& A, const vec& P);
			vec solve();
	};
}

#endif
