#include <thesis_functions.h>
#include <algorithm>    // std::max

struct parameters{
	gsl_multilarge_linear_workspace *w;
	gsl_matrix * qr;
	gsl_vector * b;

	//mat *A;
	//vec *P;
	//double *O;
	//vec *du;
} params;

mat findA(const vec& t, const mat& U, int m)
{
	mat A(m,m);

	for(int i = 0; i<m; i++){
		for(int j = 0; j<m; j++){
			if(i <= j){
				A(i,j) = innerProd(U.row(i), U.row(j), t);
				A(j,i) = A(i,j);
			}
		}
	}
	return A;
}

vec findP(const vec& t, const mat& U, const vec& dx, int m)
{
	vec P(m);

	for(int i = 0; i<m; i++)
	{
		P(i) = innerProd(U.row(i), dx, t);
	}

	return P;
}

double findO(const vec& t, const vec& dx)
{
	return innerProd(dx, dx, t);
}

double innerProd(const vec& u1, const vec& u2, const vec& time)
{
	int mid = time.size();
	int n = u1.size()/mid;

	vec a = u1.array() * u2.array();

	vec aij = a.head(mid);

	if(mid != n*mid){
		for(int k=1; k<n; k++){
			aij += a.segment(k*mid, mid);
		}
	}

	return simpson(time, aij);
}

vec findActualParam(soln_env *env, bool regs=false, const int numdivs = 1)
{
	int n = (*env->nth_soln).rows();
	int m = (*env->initial_params).size();
	int lt = (*env->time).size();

	int divs = (int) (lt/numdivs+1);
	double TOL = .00001;

	const mat measurements = *env->nth_soln;

	vec uNot(m);
	uNot = *env->initial_params;
	mat bob;
	mat U(n*m, n*lt);
	mat A(m,m), AT(m,m);
	mat B = ((*env->time)(lt-1) - (*env->time)(0))*mat::Identity(m, m);
	mat I = mat::Identity(m, m);
	for(int i=0; i<m-1; i++){
		//B(i+1,i) = -1;
	}
	mat BT = B.transpose();
	vec P(m);
	vec du(m), u1(m), u2(m);
	//double gamma;

	//params.A = &A;
	//params.P = &P;

 	double rnorm, snorm, lambda=0.05;
	gsl_matrix *qr = gsl_matrix_alloc(m,m);
	gsl_vector *b = gsl_vector_alloc(m);
	gsl_vector *x = gsl_vector_alloc(m);
	gsl_multilarge_linear_workspace *w \
		= gsl_multilarge_linear_alloc(gsl_multilarge_linear_tsqr, m);

	params.w = w;
	params.qr = qr;
	params.b = b;

	//int signum;
	//double now, last, temp = gamma;

	int LIMIT = 500;
	if(regs){
		//TOL = .0001;
		for(int j = 0; j<numdivs; j++)
		{
			if(j == numdivs-1){
				for(int i = 0; i<LIMIT; i++)
				{
					bob = qLinearRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, measurements);
					U = reshape(bob.bottomRows(n*m), m, n*lt);
					*env->nth_soln = bob.topRows(n);
					A = findA(*env->time, U, m);
					P = findP(*env->time, U, reshape(measurements - *env->nth_soln, 1, n*lt).row(0), m);
					//O = findO(*env->time, reshape(measurements - *env->nth_soln, 1, n*lt).row(0));
					//cout << "cond(A) = "<< cond(A) <<"\nDeterminant(A) = " << A.determinant() << endl; cout << "rank(A) = " << A.fullPivHouseholderQr().rank() <<endl;

					/*matToGslMat(A, qr);
					vecToGslVec(P, b);

					gsl_multilarge_linear_reset(w);
					gsl_multilarge_linear_accumulate(qr, b, w);
					gsl_multilarge_linear_solve (lambda, x, &rnorm, &snorm, w);
					*/

					// cout << "No regs norm Au-P " << (norm(A*gslVecToVec(x)-P)) << endl;
					// gamma = norm(A*gslVecToVec(x)-P);
					// if(cond(A) > 1e-8 || A.fullPivHouseholderQr().rank() < m){
					// 	lambda = 1.0;
					// 	for(int k = 0; k<55; k++){
					// 		gsl_multilarge_linear_solve (lambda, x, &rnorm, &snorm, w);
					// 		u1 = gslVecToVec(x);
					//
					// 		cout << "Norm Au-P  " << (norm(A*u1-P)) << endl;
					//
					// 		if(norm(A*u1-P) < std::max(1.1*gamma,.0005)){
					// 			break;
					// 		}
					// 		lambda = lambda/1.5;
					// 	}
					// 	du = u1;
					// }else{
					// 	gsl_multilarge_linear_solve (0.0, x, &rnorm, &snorm, w);
					// }
					du = inverse(A + lambda*B.transpose()*B)*P;
					//du = gslVecToVec(x);

					//gamma = (double) (du.transpose()*((AT*A)+lambda*lambda*B)*du)/(P.transpose()*A*du);
					//cout << "gamma = "  << (du.transpose()*((AT*A)+lambda*lambda*B)*du)/((P.transpose()*A*du)) << endl;
					//cout << "du = " << du.transpose() << endl;
					//cout << "g*du = " << gamma*du.transpose() << endl;
					//cout << "minimzed du = " << findGamma(1, &params) << endl;
					uNot += du;
					//cout << " u = " << uNot.transpose() << endl;
					latexOutput(*env->nth_soln, uNot, i+1, " &");
					if(du.norm() < TOL || std::isnan(du.norm())){
						break;
					} else if (i >= LIMIT-1){
						note("u = ");
						note(uNot);
						if(du.norm() < TOL){
							break;
						}
						else{
							log_err("Function did not converge.");
							uNot(1) = NAN;
						}
					}
				}
			}else{
				for(int i = 0; i<LIMIT; i++)
				{
					bob = qLinearRungeKutta4(*env->ode, (*env->time).head((j+1)*divs), \
						uNot, *env->initial_cond, measurements.leftCols((j+1)*divs));
					U = reshape(bob.bottomRows(n*m), m, n*(j+1)*divs);
					A = findA((*env->time).head((j+1)*divs), U, m);
					P = findP((*env->time).head((j+1)*divs), U, reshape(measurements.leftCols((j+1)*divs) - bob.topRows(n), 1, n*(j+1)*divs).row(0), m);
					//O = findO(*env->time, reshape(measurements - *env->nth_soln, 1, n*lt).row(0));

					cout <<" cond(A) = "<< cond(A) <<"\nDeterminant(A) = " << A.determinant() << endl;
					cout << "rank(A) = " << A.fullPivHouseholderQr().rank() <<endl;

					matToGslMat(A, qr);
					vecToGslVec(P, b);

					gsl_multilarge_linear_reset(w);
					gsl_multilarge_linear_accumulate(qr, b, w);

					// if(cond(A) > 1e-8 || A.fullPivHouseholderQr().rank() < m){
					// 	lambda = 1.0;
					// 	for(int k = 0; k<15; k++){
					// 		gsl_multilarge_linear_solve (lambda, x, &rnorm, &snorm, w);
					// 		u1 = gslVecToVec(x);
					//
					// 		if(norm(A*u1-P) < TOL){
					// 			break;
					// 		}
					// 		lambda = lambda/2;
					// 	}
					// 	du = u1;
					// }

					gsl_multilarge_linear_solve (lambda, x, &rnorm, &snorm, w);
					du = gslVecToVec(x);

					uNot += du;

					latexOutput(*env->nth_soln, uNot, i+1, " &");
					if(du.norm() < TOL || std::isnan(du.norm())){
						break;
					} else if (i >= LIMIT-1){
						log_err("Function did not converge.");
						note("u = ");
						note(uNot);
						uNot(1) = NAN;
					}
				}
			}
		}
		/*for(int i = 0; i<LIMIT; i++)
		{
			bob = qLinearRungeKutta4(*env->ode, *env->time, uNot, *env->initial_cond, measurements);

			U = reshape(bob.bottomRows(n*m), m, n*lt);
			*env->nth_soln = bob.topRows(n);

			A = findA(*env->time, U, m);
			P = findP(*env->time, U, reshape(measurements - *env->nth_soln, 1, n*lt).row(0), m);

			AT = A.transpose();

			cout <<" cond(A) = "<< cond(A) <<"\nDeterminant(A) = " << A.determinant() << endl;
			cout << "rank(A) = " << A.fullPivHouseholderQr().rank() <<endl;

			/*
			matToGslMat(A, gslA);
			matToGslMat(A, qr);
			vecToGslVec(P, b);

			if(A.fullPivHouseholderQr().rank() < m){
				log_err("A_N Matrix is Singular");
				gsl_linalg_QR_decomp(qr, tau);
				gsl_linalg_QR_lssolve(qr, tau, b, x, residual);
				//uNot(1) = NAN;
				//break;
			}else{
				gsl_linalg_LU_decomp (qr, perm, &signum);
				gsl_linalg_LU_solve (qr, perm, b, x);
				gsl_linalg_LU_refine (gslA, qr, perm, b, x, residual);
			}
			du = gslVecToVec(x);
			*/
			/*
			gamma = 4.0;
			du = inverse(AT*A + gamma*gamma*BT*B)*(AT*P);
			last = norm(A*du-P) + norm(gamma*(B*du));	note(last);

			while(gamma > 0.0){
				if(std::isnan(du.norm())){
					break;
				}
				du = inverse(AT*A + gamma*gamma*BT*B)*(AT*P);
				now = norm(A*du-P) + norm(gamma*(B*du));	//note(now);
				if(now < last){
					temp = gamma;
					last = now;
				}
				gamma -= 0.01;
			}
			gamma = temp;
			du = inverse(AT*A + gamma*gamma*BT*B)*(AT*P);
			now = norm(A*du-P) + norm(gamma*(B*du));
			cout << "FINAL gamma = " << gamma << "\nFINAL now = " << now << endl;

			uNot += du;
			latexOutput(*env->nth_soln, uNot, i+1, " &");
			if(du.norm() < 0.0001 || std::isnan(du.norm())){
				break;
			} else if (i >= LIMIT-1){
				log_err("Function did not converge.");
				note("u = ");
				note(uNot);
				uNot(1) = NAN;
				break;
				//exit(1);
			} else if (uNot.norm() > 1000){
				log_err("uNot is increasing without bound.");
				uNot(1) = NAN;
				break;
			}
		}*/
	}else{
		for(int j = 0; j<numdivs; j++)
		{
			if(j == numdivs-1){
				for(int i = 0; i<LIMIT; i++)
				{
					bob = qLinearRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, measurements);
					U = reshape(bob.bottomRows(n*m), m, n*lt);
					*env->nth_soln = bob.topRows(n);
					A = findA(*env->time, U, m);
					P = findP(*env->time, U, reshape(measurements - *env->nth_soln, 1, n*lt).row(0), m);
					//matToGslMat(A, qr);
					//vecToGslVec(P, b);

					du = A.inverse()*P;
					uNot += du;
					latexOutput(*env->nth_soln, uNot, i+1, " &");
					cout << i << endl;
					if(du.norm() < 0.00001 || std::isnan(du.norm())){
						break;
					} else if (i >= LIMIT-1){
						log_err("Function did not converge.");
						note("u = ");
						note(uNot);
						uNot(1) = NAN;
					}
				}
			}else{
				for(int i = 0; i<LIMIT; i++)
				{
					bob = qLinearRungeKutta4(*env->ode, (*env->time).head((j+1)*divs), \
						uNot, *env->initial_cond, measurements.leftCols((j+1)*divs));
					U = reshape(bob.bottomRows(n*m), m, n*(j+1)*divs);
					A = findA((*env->time).head((j+1)*divs), U, m);
					P = findP((*env->time).head((j+1)*divs), U, reshape(measurements.leftCols((j+1)*divs) - bob.topRows(n), 1, n*(j+1)*divs).row(0), m);

					//cout <<" cond(A) = "<< cond(A) <<"\nDeterminant(A) = " << A.determinant() << endl;
					//cout << "rank(A) = " << A.fullPivHouseholderQr().rank() <<endl;

					du = A.inverse()*P;
					uNot += du;
					latexOutput(*env->nth_soln, uNot, i+1, " &");
					if(du.norm() < 0.00001 || std::isnan(du.norm())){
						break;
					} else if (i >= LIMIT-1){
						log_err("Function did not converge.");
						note("u = ");
						note(uNot);
						uNot(1);
						break;
					}
				}
			}
		}
	}
	return uNot;
}

mat reshape(const mat& U, int n, int m)
{
	mat newU(n,m);
	newU.fill(0);
	int olt = U.row(0).size(); 	//old time(row) length
	int on = m/olt;				//on col length

	for(int i = 0; i<n; i++){
		for(int j = 0; j<on; j++){
			newU.block(i, j*olt, 1, olt) = U.row(i*on + j);
		}
	}

	return newU;
}

double norm(const mat& M){
	return M.norm();//std::sqrt(M.cwiseAbs2().sum())  ;
}

mat inverse(const mat& M){
	return M.inverse();
}

void vecToGslVec(const vec& v, gsl_vector* gslv)
{
	for (size_t i=0; i<gslv->size; i++){
		gsl_vector_set(gslv, i, v(i));
	}
}
void matToGslMat(const mat& m, gsl_matrix *gslm)
{
	//gslm = gsl_matrix_alloc(m.rows(),m.cols());
	for(int i=0; i<m.rows(); i++){
		for(int j=0; j<m.cols(); j++){
			gsl_matrix_set(gslm, i, j, m(i,j));
		}
	}
}
vec gslVecToVec(gsl_vector* gslv)
{
	vec v(gslv->size);
	for (size_t i=0; i<gslv->size; i++){
		v(i) = gsl_vector_get(gslv, i);
	}
	return v;
}

mat gslMatToMat(gsl_matrix *gslm){
	mat m(gslm->size1, gslm->size2);
	for (size_t i=0; i<gslm->size1; i++){
		for (size_t j=0; j<gslm->size2; j++){
			m(i,j) = gsl_matrix_get(gslm, i,j);
		}
	}
	return m;
}

bool allpositive(const vec& x){
	bool allpos = true;
	for(int i=0; i < x.size(); i++){
		if(x(i)>0){
		}
		else{
			allpos = false;
			break;
		}
	}
	return allpos;
}

double cond(const mat& A){
	return A.norm()*A.inverse().norm();
}

/*double fn1 (double x, void * params)
{
	struct parameters * pParams = (struct parameters *) params;
	mat A = *(pParams->A);
	vec P = *(pParams->P);
	double O = *(pParams->O);
	vec du = *(pParams->du);

	return O - 2*x*P.transpose()*du + x*x*du.transpose()*A*du;
}*/

double fn1 (double x, void * params)
{
	double rnorm, snorm;
	struct parameters * pParams = (struct parameters *) params;
	gsl_vector * du = gsl_vector_alloc((pParams->b)->size);
	gsl_multilarge_linear_reset(pParams->w);
	gsl_multilarge_linear_accumulate(pParams->qr, pParams->b, pParams->w);
	gsl_multilarge_linear_solve (x, du, &rnorm, &snorm, pParams->w);

	return rnorm;
}

double findGamma(double initialGuess, void * params)
{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m = initialGuess;
	double a = 0.0, b = 10000.0;
	gsl_function F;
	F.function = &fn1;
	F.params = params;
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc (T);
	gsl_min_fminimizer_set (s, &F, m, a, b);

	do
	{
		iter++;
		status = gsl_min_fminimizer_iterate (s);
		m = gsl_min_fminimizer_x_minimum (s);
		a = gsl_min_fminimizer_x_lower (s);
		b = gsl_min_fminimizer_x_upper (s);
		status = gsl_min_test_interval (a, b, 0.001, 0.0);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	gsl_min_fminimizer_free (s);
	return m;
}

/*
mat ichol(const mat& A){
	gsl_matrix *L;
	L = gsl_matrix_alloc(A.rows(), A.cols());
	matToGslMat(A,L);
	gsl_linalg_cholesky_decomp(L);
	for (int i=0; i<A.rows(); i++){
		for (int j=0; j<A.cols(); j++){
			if(A(i,j)<1E-13){
				gsl_matrix_set(L, i, j, 0);
			}
		}
	}
	return gslMatToMat(L);
}*/

 /* vec dulp(const mat& A, const vec& b, const vec& u){
	int m = b.size();
	int ia[m*m+1], ja[m*m+1];
    double ar[m*m+1];
	//double x[1+m]; //x[0] is the soln to obj func, x[1] -> x[m] are variables
	vec v(m);
	int result;

	glp_prob *lp;
	lp = glp_create_prob();
	glp_add_rows(lp, m);
	glp_add_cols(lp, m);
	for (int j=1; j<m+1; j++){
		glp_set_row_bnds(lp, j, GLP_FX, b(j-1), b(j-1)); // RHS of the equal-sign
		glp_set_col_bnds(lp, j, GLP_FR, -u(j-1), 2.0); //variables
		glp_set_obj_coef(lp, j, 1.0);
	}

	for (int i=0; i<m; i++){
		for (int j=0; j<m; j++){
			ia[m*i+j+1] = i+1, ja[m*i+j+1] = j+1, ar[m*i+j+1] = A(i,j);
		}
	}
	glp_load_matrix(lp, m*m, ia, ja, ar);
	result = glp_simplex(lp, NULL);
	if(result == 0){
		for (int i=0; i<m; i++){
			v(i) = glp_get_col_prim(lp, i+1);
		}
	}else{
		for (int i=0; i<m; i++){
			v(i) = NAN;
		}
	}

	glp_delete_prob(lp);
    glp_free_env();
	return v;
} */


/* du1 = dulp(A, P, uNot);
if(std::isnan(du.norm())){
	du = du1;
}else{
	matToGslMat(A, qr);
	gsl_linalg_QR_decomp(qr, tau);
	vecToGslVec(P, b);
	//gsl_linalg_QR_solve(qr, tau, b, x);
	gsl_linalg_QR_lssolve(qr, tau, b, x, residual);
	du = gslVecToVec(x);*/
//}
