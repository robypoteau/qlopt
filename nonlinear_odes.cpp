#include "nonlinear_odes.h"

mat lotka_volterra(const mpreal& t, const vec& x, const vec& u)
{
    mat result(2,1); 
	result(0) = u(0)*x(0) - u(1)*x(0)*x(1);
	result(1) = u(1)*x(0)*x(1) - u(2)*x(1);
	
	return result;
}

mat lotka_volterra_linearization(const mpreal& t, const vec& x, const vec& u, const mat& xn, const vec& time)
{
	thesis::spline Xn(time, xn.row(0));
	thesis::spline Yn(time, xn.row(1));
	
	mpreal xnt = Xn.interpolate(t);
	mpreal ynt = Yn.interpolate(t);
	
	mat result(8,1);
	result << (u(0)-u(1)*ynt)*x(0) - u(1)*xnt*x(1) + u(1)*xnt*ynt,
             u(1)*ynt*x(0) + (u(1)*xnt - u(2))*x(1) - u(1)*xnt*ynt,
            (u(0) - u(1)*ynt)*x(2) + x(0) - u(1)*xnt*x(3),
             u(1)*ynt*x(2) + (u(1)*xnt - u(2))*x(3),
            (u(0) - u(1)*ynt)*x(4) + ynt*(xnt - x(0)) - xnt*x(1) - u(1)*xnt*x(5),
             u(1)*ynt*x(4) + ynt*x(0) + (u(1)*xnt - u(2))*x(5) + xnt*(x(1) - ynt),
            (u(0) - u(1)*ynt)*x(6) - u(1)*xnt*x(7),
             u(1)*ynt*x(6) + (u(1)*xnt - u(2))*x(7) - x(1);
	
	return result;
}

double angiogenesis()
{
	return 0;
}