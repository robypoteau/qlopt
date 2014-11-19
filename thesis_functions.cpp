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

mpreal innerProd(const vec& u1, const vec& u2, const vec& time)
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