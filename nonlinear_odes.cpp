#include "nonlinear_odes.h"

mat lotka_volterra(const mpreal& time, const vec& x, const vec& u)
{
    mat result(2,1); 
	result(0) = u(0)*x(0) - u(1)*x(0)*x(1);
	result(1) = u(1)*x(0)*x(1) - u(2)*x(1);
	
	return result;
}

mat lotka_volterra_linearization(const mpreal& t, const vec& x, const vec& u, const mat& xn, const int& i /*const vec& time*/)
{
	//splinef Xn(time, x);
	//splinef Yn(time, x);
	
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

double angiogenesis()
{
	return 0;
}

mpreal spline(const vec& xVals, const vec& yVals, const mpreal& xInst)
{
	mpreal ans = 0;
	
	
	
	
	/*vec x = xVals;
	int end = xVals.size() - 1;
	
	try{
		if(xInst > xVals(end) && xInst < xVals(0)){
			throw "ERROR: xInst is Out of bounds.\n";
		}
	}
	catch(char *str){
		cout << str;
	}
	
	//find upper and lower of xInst
	x = x.array()-xInst
	int i = indexOf(x, min(x));
	
	mpreal t = (xInst - x(i))/(x(i+1)- x(i))
	mpreal a = ;
	mpreal b = ;
	//apply the functions*/
	
	
	return ans;
}
/*
mpreal min(const vec& x){
	mpreal min = x(0);
	for(int i = 1; i<x.size(); i++){
		if(x(i) < x(i-1)){
			min = x(i);
		}
	}
	
	return min;
}

// for strictly increasing or decreasing functions
int indexOf(const vec& x, mpreal xi){
	for(int i = 0; i<x.size(); i++){
		if(x(i) == xi){
			return i;
		}
	}
}*/
