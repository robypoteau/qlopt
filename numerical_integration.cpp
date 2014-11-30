#include "numerical_integration.h"
#include "spline.h"
#include <iostream>

mat rungekutta4(mat (*fhandle)(const mpreal&, const vec&, const vec&), const vec& time, const vec& u, const vec& yNot){
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
	   
	   w.col(i) = w.col(i-1) + (k1 + 2*(k2 + k3) + k4)/6;
	}
	
	return w;
}

mat rungekutta4(mat (*fhandle)(const mpreal&, const vec&, const vec&, const mat&, const vec&), const vec& time, const vec& u, const vec& yNot, const mat& xNminus){
	int N = time.size();
	mpreal a = time(0);
	mpreal b = time(N-1);
	
	//number of equations
	int m = yNot.size();
	
	//timesteps
	mpreal h = (b-a)/(N-1);
	
	//Init stuff
	mat w(m, N);
	w.fill(0);
	w.col(0) = yNot;
	
	vec k1(m), k2(m), k3(m), k4(m);

	for (int i = 1; i<N; i++)
	{
		k1 = h*fhandle(time(i), w.col(i-1), u, xNminus, time);
		k2 = h*fhandle(time(i) + h/2, w.col(i-1) + k1/2, u, xNminus, time);
		k3 = h*fhandle(time(i) + h/2, w.col(i-1) + k2/2, u, xNminus, time);
		k4 = h*fhandle(time(i) + h, w.col(i-1) + k3, u, xNminus, time);
	   
		w.col(i) = w.col(i-1) + (k1 + 2*(k2 + k3) + k4)/6;
	}
	
	return w;
}

mat qLinearRungeKutta4(mat (*sys)(const mpreal&, const vec&, const vec&), const vec& time, const vec& u, const vec& yNot, const mat& xNminus)
{
	int N = time.size();
	int n = xNminus.col(1).size();
	
	//number of equations
	int m = yNot.size();
	
	//timesteps
	mpreal h = (time(N-1) - time(0))/(N-1);
	
	//Init stuff
	mat w(m, N);
	w.fill(0);
	w.col(0) = yNot;
	
	vec k1(m), k2(m), k3(m), k4(m);
	
	thesis::spline Xn[n];
	for(int i = 0; i < n; i++){
		Xn[i].update(time, xNminus.row(i));
	}
	
	for (int i = 1; i<N; i++)
	{
		k1 = h*qlinear(sys, time(i), w.col(i-1), u, Xn, n);
		k2 = h*qlinear(sys, time(i) + h/2, w.col(i-1) + k1/2, u, Xn, n);
		k3 = h*qlinear(sys, time(i) + h/2, w.col(i-1) + k2/2, u, Xn, n);
		k4 = h*qlinear(sys, time(i) + h, w.col(i-1) + k3, u, Xn, n);
	   
		w.col(i) = w.col(i-1) + (k1 + 2*(k2 + k3) + k4)/6;
	}
	
	return w;
}

mpreal simpson(const vec& t, const vec& x)
{
	thesis::spline sim(t,x);
	mpreal temp;
	
	int N = 1001;
	int end = t.size()-1;
	mpreal h = (t(end)- t(0))/(N-1);
	
    mpreal area = x(0) + x(end);
    for(int i = 1; i<N-1; i++)
	{
		temp = sim.interpolate(t(0)+h*i);
        if((i+1)%2 == 0)
		{
            area += 2*temp;
        }
		else
		{
            area += 4*temp;
		}
    }
	
    area = h/3 * area;
	return area;
}

mat der(const mat& dx, const mpreal& dt){
	mat ans(dx.size(),1);
	ans << dx/dt;
	return ans;
}

mat qlinear(mat (*sys)(const mpreal&, const vec&, const vec&), const mpreal& t, const vec& x, const vec& u, thesis::spline* Xn, int n)
{
	int m = u.size();
	mpreal step = 2.2E-16;
	
	vec xn1(n); // this is x_N-1
	vec dxn(n);
	
	for(int i=0; i<n; i++)
	{
		xn1(i) = Xn[i].interpolate(t);
		dxn(i) = x(i) - xn1(i); // x_N - x_{N-1}
	}
	
	mat dx(n,n);// for x derivative
	dx << mat::Identity(n,n)*step;
	
	//First few lines of linearization
	mat fx(n,1);
	fx = sys(t, xn1, u);
	mat ans(n+n*m, 1);
	ans << mat::Zero(n+n*m, 1);
	ans.block(0, 0, n, 1) = fx;
	for(int j=0; j<n; j++)
	{
		ans.block(0, 0, n, 1) += der(sys(t, xn1+dx.col(j), u) - fx, step)*dxn(j);
	}
	
	mat dfdx(n,1);
	mat dfdu(n,1);
	
	// the Un part of t the linearization
	mat dun(m, m);// for u derivative
	dun << mat::Identity(m,m)*step;
	
	int ind;
	for(int j=0; j<m; j++)
	{
		ind = (j+1)*n;
		dfdu = sys(t, xn1, u+dun.col(j));
		ans.block(ind, 0, n, 1) = der(dfdu - fx, step); // df/du
		for(int k=0; k<n; k++){
			dfdx = sys(t, xn1+dx.col(k), u);
			ans.block(ind, 0, n, 1) += der(dfdx  - fx, step)*x(ind+k); //J*Un
			ans.block(ind, 0, n, 1) += der(sys(t, xn1+dx.col(k), u+dun.col(j)) - dfdu - dfdx + fx, step*step)*dxn(k); //phi_ij
		}
	}
	
	return ans;
}