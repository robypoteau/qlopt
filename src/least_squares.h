#ifndef LEAST_SQUARES_H
#define LEAST_SQUARES_H

#include <misc.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>

namespace thesis{
	class lsquares {
		private:
			size_t num_obs;
			size_t order;
			gsl_matrix *X, *cov;
			gsl_vector *c,  *gslx, *gsly, *tau, *res;
			gsl_multifit_linear_workspace *mws;
			
			void vecToGslVec(const vec& v, gsl_vector *gslv);
			void generateX(gsl_matrix *M, gsl_vector *gslv);
			void generate_xi(double xi, gsl_vector *gslv);
			
		public:
			lsquares(){}
			~lsquares();
			lsquares(const size_t n, const size_t p);
			void init(const size_t n, const size_t p);
			void update(const vec& x, const vec& y);
			double interpolate(double xi);
	};
}

#endif