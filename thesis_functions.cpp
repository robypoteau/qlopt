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
vec findActualParam(soln_env *env, bool regs)
{
	int n = (*env->measurements).rows();
	int m = (*env->initial_params).size();
	int lt = (*env->time).size();
	
	vec uNot;
	uNot << *env->initial_params; 
	//check_mem(uNot);
	//check(uNot(0) == 2,"Not initializing");
	mat bob((*env->initial_cond).size(), lt);
	mat U(n*m, n*lt);
	mat A(m,m);
	mat B = mat::Identity(m, m);
	vec P(m);
	vec du(m);
	
	bob = qLinearRungeKutta4(*env->ode, *env->time, *env->initial_params, *env->initial_cond, *env->nth_soln);

	U = reshape(bob.bottomRows(n*m), m, n*lt);
	*env->nth_soln = bob.topRows(n);
	
	*env->nth_soln = *env->measurements;
	
	double gamma = 1.0;
	
	for(int i = 0; i<150; i++)
	{
		A = findA(*env->time, U, m);
		P = findP(*env->time, U, reshape(*env->measurements - *env->nth_soln, 1, n*lt).row(0), m);
		
		if((regs) /*&& cond(A) > 1E+6) || isnan(cond(A))*/){ //create a reg function that accepts different types of reg function
			do{
				gamma *= .5;
				du = inverse(A.transpose()*A + gamma*gamma*B.transpose()*B)*A.transpose()*P;
			}while(norm(A*du-P) > 0.1);
				
		}else{
			du = A.inverse()*P;
		}
		
		uNot += du;
		du(0) = du.norm();
		if(du(0) < 0.00001 || isnan(du(0) || i==149)){
			//latexOutput(*env->nth_soln = msmt,*env->initial_params, -1, ",");
			//cout << endl;
			break;
		}
	}
	
	return du;

error:
	exit(1);
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