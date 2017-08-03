#ifndef TSQR_H
#define TSQR_H

#include <misc.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_linalg.h>

namespace thesis{
	class tsqr {
		private:
			size_t num_obs;
			size_t order, divs;
			gsl_matrix *X;
			gsl_vector *gslx, *gsly, *c;
			gsl_multilarge_linear_workspace * w;
			double rnorm, snorm, rcond, lambda;

			void vecToGslVec(const vec& v, gsl_vector *gslv);
			void generateX(gsl_matrix *M, gsl_vector *gslv);
			void generate_xi(double xi, gsl_vector *gslv);

		public:
			tsqr(){}
			~tsqr();
			tsqr(const size_t n, const size_t p);
			void init(const size_t n, const size_t p);
			void update(const vec& x, const vec& y);
			double interpolate(double xi);
	};
}

#endif
