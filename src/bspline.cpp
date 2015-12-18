#include <bspline.h>
#include <dbg.h>

namespace thesis{

	bspline::bspline(const size_t order, const size_t ncoeffs, const size_t vlen){
		init(order, ncoeffs, vlen);
	}
	
	void bspline::init(const size_t order, const size_t num_coeffs, const size_t vlen){
		ncoeffs = num_coeffs;
		size_t nbreak = ncoeffs  - order + 2;
		
		bws  = gsl_bspline_alloc(order, nbreak);
		B    = gsl_vector_alloc(ncoeffs);
		
		n    = vlen;
		gslx = gsl_vector_alloc(n);
		gsly = gsl_vector_alloc(n);
		X    = gsl_matrix_alloc(n, ncoeffs);
		c    = gsl_vector_alloc(ncoeffs);
		w    = gsl_vector_alloc(n);
		cov  = gsl_matrix_alloc(ncoeffs, ncoeffs);
		mw   = gsl_multifit_linear_alloc(n, ncoeffs);
	}
			
	bspline::~bspline(){
		gsl_bspline_free(bws);
		gsl_vector_free(gslx);
		gsl_vector_free(gsly);
		
		gsl_vector_free(B);
		gsl_matrix_free(X);
		gsl_vector_free(c);
		gsl_vector_free(w);
		gsl_matrix_free(cov);
		gsl_multifit_linear_free(mw);
	}
	
	void bspline::update(const vec& x, const vec& y){
			check(x.size() == y.size(), "Mismatched Vector Sizes");
			
			vecToGslVec(x, gslx);
			vecToGslVec(y, gsly);
			
			for (size_t i=0; i<n; i++){
				gsl_vector_set(w, i, 1.0);///(.01*y(i)*y(i)));
			}
			gsl_bspline_knots_uniform(x(0),x(n-1), bws);
			findB();
			findCoeff();
	}
	
	void bspline::findB(){
		for(size_t i = 0; i < n; ++i)
		{
			double xi = gsl_vector_get(gslx, i);
		
			/* compute B_j(xi) for all j */
			gsl_bspline_eval(xi, B, bws);
		
			/* fill in row i of X */
			for(size_t j = 0; j < ncoeffs; ++j)
			{
				check(j<ncoeffs, "Index out of bounds");
				double Bj = gsl_vector_get(B, j);
				gsl_matrix_set(X, i, j, Bj);
			}
		}
	}
	
	void bspline::findCoeff(){
		/* do the fit */
		double chisq;
		gsl_multifit_wlinear(X, w, gsly, c, cov, &chisq, mw);
	}
	
	double bspline::interpolate(double xi){
		/* output the smoothed curve */
		double yi, yerr;
		
		gsl_bspline_eval(xi, B, bws);
		gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
		return yi;
	}
	
	void bspline::vecToGslVec(const vec& v, gsl_vector* gslv)
	{
		check(v.size() == gslv->size, "Mismatched Vector Sizes");	
		for (size_t i=0; i<gslv->size; i++){
			gsl_vector_set(gslv, i, v(i));
		}
	}
}
