#ifndef NUMERICAL_INTEGRATION_H
#define NUMERICAL_INTEGRATION_H

#include <misc.h>
#include <spline.h>
#include <nonlinear_odes.h>
#include <gsl/gsl_integration.h>
#include <boost/numeric/odeint.hpp>
#include <vector>

template<typename statetype>
class nonlinear_system_wrapper
{
	sys ode;
	vec param;
	unsigned int counter;
	vec stToVec(statetype x)
	{
		vec y(x.size());
		int i=0;
		for (auto xval : x)
		{ // access by value, the type of i is int
        	y(i++)=xval;
        }
        return x;
	}

	statetype vecToSt(vec x)
	{
		int n = x.size();
		statetype y(n);
		
		for (int i=0; i<n; i++)
		{
			y[i]=x(i);	
        }
        return y;
	}
	
  public:
	nonlinear_system_wrapper(sys ode)
	{
		this->ode = ode;
		counter = 0;
	}
	
	void setParam(vec param)
	{
		this->param = param;
	}
	
	unsigned int getCounter(){
		return counter;
	}
	
	void operator()( const statetype &x , statetype &dxdt , const double t )
	{
		vec xv = stToVec(x);
		dxdt = vecToSt(sys(t, xv, this->param));
		counter++;
	}
};

template<typename statetype>
class linearized_system_wrapper
{
	qsys ode;
	vec param;
	thesis::spline* xn;
	
	unsigned int counter;
	vec stToVec(statetype x)
	{
		vec y(x.size());
		int i=0;
		for (auto xval : x)
		{ // access by value, the type of i is int
        	y(i++)=xval;
        }
        return x;
	}

	statetype vecToSt(vec x)
	{
		int n = x.size();
		statetype y(n);
		
		for (int i=0; i<n; i++)
		{
			y[i]=x(i);	
        }
        return y;
	}
	
  public:
	linearized_system_wrapper(qsys ode)
	{
		this->ode = ode;
		counter = 0;
	}
	
	void setParam(vec param, thesis::spline* xn)
	{
		this->param = param;
	}
	
	unsigned int getCounter(){
		return counter;
	}
	
	void operator()( const statetype &x , statetype &dxdt , const double t )
	{
		vec xv = stToVec(x);
		dxdt = vecToSt(sys(t, xv, this->param));
		counter++;
	}
};


mat rungekutta4(string fname, const vec& time, const vec& u, const vec& yNot);
mat qLinearRungeKutta4(string fname, const vec& time, const vec& u, const vec& yNot, std::vector<thesis::spline> Xn);
mat qlRungeKutta4(string fname, const vec& time, const vec& u, const vec& yNot, std::vector<thesis::spline> Xn);
mat qlOdeInt(string fname, const vec& time, const vec& u, const vec& yNot, std::vector<thesis::spline> Xn);
double simpson(const vec& t, const vec& x);
mat qlinear(sys fhandle, const double& t, const vec& x, const vec& u, std::vector<thesis::spline> Xn);
double gsl_integration(const vec& t, const vec& x);
#endif
