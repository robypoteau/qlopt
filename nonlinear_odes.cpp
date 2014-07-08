#include "nonlinear_odes.h"

mat lotka_volterra(const mpreal& time, const vec& x, const vec& u)
{
    mat result(2,1); 
	result(0) = u(0)*x(0) - u(1)*x(0)*x(1);
	result(1) = u(1)*x(0)*x(1) - u(2)*x(1);
	
	return result;
}

mat lotka_volterra_linearization(const mpreal& time, const vec& x, const vec& u, const mat& xn, const int& i)
{
	mat result(8,1);
	result << (u(0)-u(1)*xn(1,i))*x(0) - u(1)*xn(0,i)*x(1) + u(1)*xn(0,i)*xn(1,i),
             u(1)*xn(1,i)*x(0) + (u(1)*xn(0,i) - u(2))*x(1) - u(1)*xn(0,i)*xn(1,i),
            (u(0) - u(1)*xn(1,i))*x(2) + x(0) - u(1)*xn(0,i)*x(3),
             u(1)*xn(1,i)*x(2) + (u(1)*xn(0,i) - u(2))*x(3),
            (u(0) - u(1)*xn(1,i))*x(4) + xn(1,i)*(xn(0,i) - x(0)) - xn(0,i)*x(1) - u(1)*xn(0,i)*x(5),
             u(1)*xn(1,i)*x(4) + xn(1,i)*x(0) + (u(1)*xn(0,i) - u(2))*x(5) + xn(0,i)*(x(1) - xn(1,i)),
            (u(0) - u(1)*xn(1,i))*x(6) - u(1)*xn(0,i)*x(7),
             u(1)*xn(1,i)*x(6) + (u(1)*xn(0,i) - u(2))*x(7) - x(1);
	
	return result;
}

/*mat lotka_volterra(const vec& time, const vec& x, const vec& u)
{// might not need this overloaded version.
	int N = time.size();
	mat result(2, N);
    for(int i = 0; i<N; i++)
	{
		result(0, i) = u(0)*x(0) - u(1)*x(0)*x(1);
		result(1, i) = u(1)*x(0)*x(1) - u(2)*x(1);
	}
	return result;
}


double angiogenesis()
{
	return 0;
}*/
