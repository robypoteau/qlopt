#ifndef LEAST_SQUARES_H
#define LEAST_SQUARES_H

#include <misc.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>

namespace thesis{
	class lsquares {
		private:
			gsl_matrix *X, *cov;
			gsl_vector *c, *y;
			
			gsl_multifit_linear_workspace *mw;
			
		public:
			lsquares(){}
			lsquares(const size_t num_obs, const size_t order);
			void init(const size_t num_obs, const size_t order);
			void update(const vec& x, const vec& y);
			double interpolate();
	};
}

#endif