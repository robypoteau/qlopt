#ifndef BSPLINE_H
#define BSPLINE_H

#include <misc.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>

namespace thesis{
	class bspline{
		private:
			gsl_bspline_workspace* bws;
			gsl_vector *gslx, *gsly;
			size_t ncoeffs;
			size_t n;
			
			gsl_matrix *X, *cov;
			gsl_vector *B;
			gsl_vector *c, *w;
			
			gsl_multifit_linear_workspace *mw;
			
			void vecToGslVec(const vec& v, gsl_vector* gslv);
			
		public:
			bspline(){}
			bspline(const size_t order, const size_t ncoeffs, const size_t vlen);
			void init(const size_t order, const size_t ncoeffs, const size_t vlen);
			~bspline();
			void update(const vec& x, const vec& y);
			void findB();
			void findCoeff();
			double interpolate(double xi);
	};
}

#endif