#include <expo_tsqr.h>

namespace thesis{
	expo_tsqr::expo_tsqr(const size_t n)
	{
		init(n);
	}
			
	void expo_tsqr::init(const size_t n)
	{
		num_obs = n;
		w = gsl_multilarge_linear_alloc(gsl_multilarge_linear_tsqr, order);
		X    = gsl_matrix_alloc(num_obs, order);
		gslx = gsl_vector_alloc(num_obs);
		gsly = gsl_vector_alloc(num_obs);
		c    = gsl_vector_alloc(order);
		lambda = 0.0;
	}

	expo_tsqr::~expo_tsqr()
	{
		gsl_multilarge_linear_free(w);
		gsl_matrix_free(X);
		gsl_vector_free(gslx);
		gsl_vector_free(gsly);
		gsl_vector_free(c);
	}
	
	void expo_tsqr::update(const vec& x, const vec& y)
	{
		vecToGslVec(x, gslx);
		generate_logy(y, gsly);
		generateX(X, gslx);
		gsl_multilarge_linear_reset(w);
		gsl_multilarge_linear_accumulate(X, gsly, w);
		gsl_multilarge_linear_solve (lambda, c, &rnorm, &snorm, w);
	}
	
	void expo_tsqr::vecToGslVec(const vec& v, gsl_vector* gslv)
	{	
		for (size_t i=0; i<gslv->size; i++){
			gsl_vector_set(gslv, i, v(i));
		}
	}

	void expo_tsqr::generateX(gsl_matrix *M, gsl_vector* gslv)
	{
		gsl_vector *temp = gsl_vector_alloc(order);
		for (size_t i=0; i<M->size1; i++){
			generate_xi(gsl_vector_get(gslv, i), temp);
			gsl_matrix_set_row(M, i, temp);
		}
	}
	
	void expo_tsqr::generate_xi(double xi, gsl_vector *gslv)
	{
		//gsl_vector_set(gslv, 0, 1);
		gsl_vector_set(gslv, 0, log(xi));
		gsl_vector_set(gslv, 1, xi);
	}

	void expo_tsqr::generate_logy(const vec& y, gsl_vector *gslv)
	{
		for (size_t i=0; i<num_obs; i++){
			gsl_vector_set(gslv, i, log(y(i)));
		}
	}
	
	double expo_tsqr::interpolate(double xi)
	{
		/*
		double logA = gsl_vector_get(c,0);
		double B = gsl_vector_get(c,1);
		double C = gsl_vector_get(c,2);
		return exp(logA)*pow(xi,B)*exp(C*xi);
		*/
		double B = gsl_vector_get(c,0);
		double C = gsl_vector_get(c,1);
		return pow(xi,B)*exp(C*xi);
	}
}