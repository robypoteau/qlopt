#include <least_squares.h>

namespace thesis{
	lsquares::lsquares(const size_t num_obs, const size_t order)
	{
		init(num_obs, order);
	}
			
	void lsquares::init(const size_t num_obs, const size_t order)
	{
		X = gsl_matrix_alloc(num_obs, order);
		y = gsl_vector_alloc (num_obs);
		c = gsl_vector_alloc (num_obs);
		cov = gsl_matrix_alloc (order, order);
	}

	lsquares::~lsquares()
	{
		free(X);
		free(cov);
		free(y);
		free(c);
	}
	
	void lsquares::update(const vec& x, const vec& y)
	{
	
	}

	double lsquares::interpolate()
	{
		return 0.0;
	}
}