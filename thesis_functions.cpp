#include "thesis_functions.h"
#include "numerical_integration.h"

mat findA(const vec& t, const mat& U, int m)
{
	mat A(m,m);
	A.fill(0);
	
	for(int i = 0; i<m; i++){
		for(int j = 0; j<m; j++){
			A(i,j) = innerProd(U.row(i), U.row(j), t);
			//refactor can be further optimized - since Aij = Aji
		}
	}
	return A;
}

vec findP(const vec& t, const mat& U, const vec& dx, int m)
{
	vec P(m);
	P.fill(0);
	
	for(int i = 0; i<m; i++)
	{
		P(i) = innerProd(U.row(i), dx, t);
	}
	
	return P;
}

mpreal innerProd(const vec& u1, const vec& u2, const vec& time)
{
	mpreal ans = 0;
	int n = u1.size()/time.size();
	
	vec a = u1.array() * u2.array();
	
	int mid = u1.size()/n;
	
	vec aij = a.head(mid);
	
	if(mid != n*mid){
		for(int k=0; k<(n-1); k++){
			aij += a.segment(k*mid, mid);
		}
	}
	
	ans = simpson(time, aij);
	
	return ans;
}