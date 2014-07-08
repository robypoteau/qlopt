#include "numerical_integration.h"

using namespace mpfr;
using namespace Eigen;
using namespace std;

mat rungekutta4(mat (*fhandle)(const mpreal&, const vec&, const vec&), const vec& time, const vec& u, const vec& yNot)
{
	int N = time.size() - 1;
	mpreal a = time(0);
	mpreal b = time(N);
	
	//number of equations
	int m = yNot.size();
	
	//timesteps
	mpreal h = (b-a)/N;
	
	//Init stuff
	mat w(m, N+1);
	w.fill(0);
	w.col(0) = yNot;
	
	vec k1(m), k2(m), k3(m), k4(m);
	
	for (int i = 1; i<N+1; i++)
	{
       k1 = h*fhandle(time(i), w.col(i-1), u);
       k2 = h*fhandle(time(i) + h/2, w.col(i-1) + k1/2, u);
       k3 = h*fhandle(time(i) + h/2, w.col(i-1) + k2/2, u);
       k4 = h*fhandle(time(i) + h, w.col(i-1) + k3, u);
	   
	   w.col(i) = w.col(i-1) + (k1 + 2*k2 + 2*k3 + k4)/6;
	}
	
	return w;
}

mat rungekutta4(mat (*fhandle)(const mpreal&, const vec&, const vec&, const mat&, const int&), const vec& time, const vec& u, const vec& yNot, const mat& xNminus)
{
	int N = time.size() - 1;
	mpreal a = time(0);
	mpreal b = time(N);
	
	//number of equations
	int m = yNot.size();
	
	//timesteps
	mpreal h = (b-a)/N;
	
	//Init stuff
	mat w(m, N+1);
	w.fill(0);
	w.col(0) = yNot;
	
	vec k1(m), k2(m), k3(m), k4(m);

	for (int i = 1; i<N+1; i++)
	{
		k1 = h*fhandle(time(i), w.col(i-1), u, xNminus, i);
		k2 = h*fhandle(time(i) + h/2, w.col(i-1) + k1/2, u, xNminus, i);
		k3 = h*fhandle(time(i) + h/2, w.col(i-1) + k2/2, u, xNminus, i);
		k4 = h*fhandle(time(i) + h, w.col(i-1) + k3, u, xNminus, i);
	   
		w.col(i) = w.col(i-1) + (k1 + 2*k2 + 2*k3 + k4)/6;
	}
	
	return w;
}

mpreal simpson(const vec& t, const vec& x)
{
	int N = t.size();
	int end = N-1;
	mpreal h = (t(end)- t(0))/N;
	
    mpreal area = x(0) + x(end);
    for(int i = 1; i<end; i++)
	{
        if((i+1)%2 == 0)
		{
            area += 4*x(i);
        }
		else
		{
            area += 2*x(i);
		}
    }
	
    area = h/3 * area;
	return area;
}