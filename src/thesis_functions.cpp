#include <dbg.h>
#include <thesis_functions.h>
#include <glpk.h>
#include <nonlinear_odes.h>
#include <regs.h>

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
vec findActualParam(soln_env *env, bool regs=false)
{
	int n = (*env->measurements).rows(); 
	int m = (*env->initial_params).size(); 
	int lt = (*env->time).size(); 		
	
	vec uNot(m);
	uNot = *env->initial_params; 
	mat bob(n*(m+1), lt);
	mat U(n*m, n*lt);
	mat A(m,m);
	mat B = mat::Identity(m, m);
	vec P(m);
	vec du(m);
	
	gsl_matrix *qr; qr = gsl_matrix_alloc(m,m);
	gsl_vector *tau; tau = gsl_vector_alloc(m);
	gsl_vector *b; b = gsl_vector_alloc(m);
	gsl_vector *x; x = gsl_vector_alloc(m);
	gsl_vector *residual; residual = gsl_vector_alloc(m);
	
	int LIMIT = 280;
	
	if(regs){
		uNot = regularization(env);
	}else{
		for(int i = 0; i<LIMIT; i++)
		{
			bob = qLinearRungeKutta4(*env->ode, *env->time, uNot, *env->initial_cond, *env->nth_soln);
	
			U = reshape(bob.bottomRows(n*m), m, n*lt);
			*env->nth_soln = bob.topRows(n);
		
			A = findA(*env->time, U, m);
			P = findP(*env->time, U, reshape(*env->measurements - *env->nth_soln, 1, n*lt).row(0), m);
			
			
				/* du1 = dulp(A, P, uNot);
				if(std::isnan(du.norm())){
					du = du1;
				}else{ */
					matToGslMat(A, qr);
					gsl_linalg_QR_decomp(qr, tau);
					vecToGslVec(P, b);
					gsl_linalg_QR_lssolve(qr, tau, b, x, residual);
					du = gslVecToVec(x);
				//}
			
			uNot += du; 
			if(du.norm() < 0.00001 || std::isnan(du.norm())){
				break;
			} else if (i >= LIMIT-1){
				log_err("Function did not converge.");
				note("u = ");
				note(uNot);
				exit(1);
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
	return M.norm();
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
}  */