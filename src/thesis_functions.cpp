#include <dbg.h>
#include <thesis_functions.h>

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
	double gamma;
	
	gsl_matrix *qr; qr = gsl_matrix_alloc(m,m);
	gsl_vector *tau; tau = gsl_vector_alloc(m);
	gsl_vector *b; b = gsl_vector_alloc(m);
	gsl_vector *x; x = gsl_vector_alloc(m);
	gsl_vector *residual; residual = gsl_vector_alloc(m);
	
	int LIMIT = 220;
	for(int i = 0; i<LIMIT; i++)
	{
		bob = qLinearRungeKutta4(env->ode, *env->time, uNot, *env->initial_cond, *env->nth_soln);

		U = reshape(bob.bottomRows(n*m), m, n*lt);
		*env->nth_soln = bob.topRows(n);
	
		A = findA(*env->time, U, m);
		P = findP(*env->time, U, reshape(*env->measurements - *env->nth_soln, 1, n*lt).row(0), m);
		
		
		if((regs) /*&& cond(A) > 1E+6) || isnan(cond(A))*/){ //create a reg function that accepts different types of reg function
		gamma = 1.0;
			do{
				gamma *= .5;
				du = inverse(A.transpose()*A + gamma*gamma*B.transpose()*B)*A.transpose()*P;
			}while(norm(A*du-P) > 0.1);
				
		}else{
			matToGslMat(A, qr);
			gsl_linalg_QR_decomp(qr, tau);
			vecToGslVec(P, b);
			gsl_linalg_QR_lssolve(qr, tau, b, x, residual);
			du = gslVecToVec(x);
			//du = A.inverse()*P;
		}
		
		uNot += du; 
		if(du.norm() < 0.00001 || isnan(du.norm())){
			break;
		} else if (i >= LIMIT){
			log_err("Function did not converge.");
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