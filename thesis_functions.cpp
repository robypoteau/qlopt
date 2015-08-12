#include "dbg.h"
#include "thesis_functions.h"
#include "numerical_integration.h"

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
	
	for(int i = 0; i<150; i++)
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
			du = A.inverse()*P;
		}
		
		uNot += du; 
		if(du.norm() < 0.00001 || isnan(du.norm() || i==149)){
			break;
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