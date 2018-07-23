#include <tallskinnyqr.h>

namespace thesis{
	tallskinnyqr::tallskinnyqr(const size_t n, const size_t p, 
		const size_t nlcurve)
	{
		init(n, p, nlcurve);
	}

	void tallskinnyqr::init(const size_t n, const size_t p,
		const size_t nlcurve)
	{
		num_obs = n;
		order = p;
		w = gsl_multilarge_linear_alloc(gsl_multilarge_linear_tsqr, p);
		gslA = gsl_matrix_alloc(num_obs, order);
		gslx = gsl_vector_alloc(order);
		gslb = gsl_vector_alloc(num_obs);
		lambda = 0.0;
		
		reg_param = gsl_vector_alloc(nlcurve);
		rho = gsl_vector_alloc(nlcurve);
		eta = gsl_vector_alloc(nlcurve);
	}

	tallskinnyqr::~tallskinnyqr()
	{
		gsl_multilarge_linear_free(w);
		gsl_matrix_free(gslA);
		gsl_vector_free(gslx);
		gsl_vector_free(gslb);
	}

	void tallskinnyqr::update(const mat& A, const vec& P)
	{
		vecToGslVec(P, gslb);
		matToGslMat(A, gslA);
		gsl_multilarge_linear_reset(w);
		gsl_multilarge_linear_accumulate(gslA, gslb, w);
		
	}
	vec tallskinnyqr::solve()
	{
		gsl_multilarge_linear_lcurve(reg_param, rho, eta, w);
		gsl_multilarge_linear_solve(lambda, gslx, &rnorm, &snorm, w);
		gsl_multilarge_linear_rcond(&rcond, w);
		std::cout << "rcond = " << rcond << endl << endl;
		return gslVecToVec(gslx);
	}
	vec tallskinnyqr::rsolve(double alpha)
	{
		gsl_multilarge_linear_solve(alpha, gslx, &rnorm, &snorm, w);
		gsl_multilarge_linear_rcond(&rcond, w);
		std::cout << "rcond = " << rcond << endl << endl;
		return gslVecToVec(gslx);
	}
	void tallskinnyqr::vecToGslVec(const vec& v, gsl_vector* gslv)
	{
		for (size_t i=0; i<gslv->size; i++){
			gsl_vector_set(gslv, i, v(i));
		}
	}
	void tallskinnyqr::matToGslMat(const mat& m, gsl_matrix *gslm)
	{
		//gslm = gsl_matrix_alloc(m.rows(),m.cols());
		for(int i=0; i<m.rows(); i++){
			for(int j=0; j<m.cols(); j++){
				gsl_matrix_set(gslm, i, j, m(i,j));
			}
		}
	}
	vec tallskinnyqr::gslVecToVec(gsl_vector* gslv)
	{
		vec v(gslv->size);
		for (size_t i=0; i<gslv->size; i++){
			v(i) = gsl_vector_get(gslv, i);
		}
		return v;
	}

	mat tallskinnyqr::gslMatToMat(gsl_matrix *gslm){
		mat m(gslm->size1, gslm->size2);
		for (size_t i=0; i<gslm->size1; i++){
			for (size_t j=0; j<gslm->size2; j++){
				m(i,j) = gsl_matrix_get(gslm, i,j);
			}
		}
		return m;
	}

}
