#include <numerical_integration.h>

mat rungekutta4(sys fhandle, const vec& time, const vec& u, const vec& yNot){
	int N = time.size();
		
	//number of equations
	int m = yNot.size();
	
	//timesteps
	double h = time(1)-time(0);
	
	//Init stuff
	mat w(m, N);
	w.fill(0);
	w.col(0) = yNot;
	
	vec k1(m), k2(m), k3(m), k4(m);
	
	for (int i = 0; i<N-1; i++)
	{
       k1 = h*fhandle(time(i), 		 w.col(i), 		  u);
       k2 = h*fhandle(time(i) + h/2, w.col(i) + k1/2, u);
       k3 = h*fhandle(time(i) + h/2, w.col(i) + k2/2, u);
       k4 = h*fhandle(time(i) + h, 	 w.col(i) + k3,   u);
	   
	   w.col(i+1) = w.col(i) + (k1 + 2*(k2 + k3) + k4)/6;
	}
	
	return w;
}

mat qLinearRungeKutta4(sys fhandle, const vec& time, const vec& u, const vec& yNot, const mat& xNminus)
{
	int N = time.size();
	int n = xNminus.col(1).size();
	
	//number of equations
	int m = yNot.size();
	
	//timestep
	double h = time(1)-time(0);
	
	//init stuff
	mat w(m, N);
	w.fill(0);
	w.col(0) = yNot;
	
	vec k1(m), k2(m), k3(m), k4(m);
	
	thesis::spline Xn[n];
	for(int i = 0; i < n; i++){
		Xn[i].update(time, xNminus.row(i));
	}
	
	for (int i = 0; i<N-1; i++)
	{
		
		k1 = h*qlinear(fhandle, time(i),       w.col(i),        u, Xn, n);
		k2 = h*qlinear(fhandle, time(i) + h/2, w.col(i) + k1/2, u, Xn, n);
		k3 = h*qlinear(fhandle, time(i) + h/2, w.col(i) + k2/2, u, Xn, n);
		k4 = h*qlinear(fhandle, time(i) + h,   w.col(i) + k3,   u, Xn, n);
		
		//cout <<"(" << i <<")\nk1"<< k1 << "\nk2:" << k2 << "\nk3:" << k3 << "\nk4:" << k4 <<endl;
		w.col(i+1) = w.col(i) + (k1 + 2*(k2 + k3) + k4)/6;
	}
	
	return w;
}

double simpson(const vec& t, const vec& x)
{
	thesis::spline sim(t,x);
	double temp;
	
	int N = 1001;
	int end = t.size()-1;
	double h = (t(end)- t(0))/(N-1);
	
    double area = x(0) + x(end);
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

mat der(const mat& dx, const double& dt){
	mat ans(dx.size(),1);
	ans << dx/dt;
	return ans;
}

mat qlinear(sys fhandle, const double& t, const vec& x, const vec& u, thesis::spline* Xn, int n)
{
	int m = u.size();
	double step = 1E-8;
	
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
	fx = fhandle(t, xn1, u);
	mat ans(n+n*m, 1);
	ans << mat::Zero(n+n*m, 1);
	ans.block(0, 0, n, 1) = fx;
	for(int j=0; j<n; j++)
	{
		ans.block(0, 0, n, 1) += der(fhandle(t, xn1+dx.col(j), u) - fx, step)*dxn(j);
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
		dfdu = fhandle(t, xn1, u+dun.col(j));
		ans.block(ind, 0, n, 1) = der(dfdu - fx, step); // df/du
		for(int k=0; k<n; k++){
			dfdx = fhandle(t, xn1+dx.col(k), u);
			ans.block(ind, 0, n, 1) += der(dfdx  - fx, step)*x(ind+k); //J*Un
			ans.block(ind, 0, n, 1) += der(fhandle(t, xn1+dx.col(k), u+dun.col(j)) - dfdu - dfdx + fx, step*step)*dxn(k); //phi_ij
		}
	}
	
	return ans;
}

/*
mat rungekutta4(mat (*fhandle)(const double&, const vec&, const vec&, const mat&, const vec&), const vec& time, const vec& u, const vec& yNot, const mat& xNminus){
	int N = time.size();
	double a = time(0);
	double b = time(N-1);
	
	//number of equations
	int m = yNot.size();
	
	//timesteps
	double h = time(1)-time(0);
	
	//Init stuff
	mat w(m, N);
	w.fill(0);
	w.col(0) = yNot;
	
	vec k1(m), k2(m), k3(m), k4(m);

	for (int i = 1; i<N; i++)
	{
		k1 = h*fhandle(time(i-1),		w.col(i-1), u, xNminus, time);
		k2 = h*fhandle(time(i-1) + h/2, w.col(i-1) + k1/2, u, xNminus, time);
		k3 = h*fhandle(time(i-1) + h/2, w.col(i-1) + k2/2, u, xNminus, time);
		k4 = h*fhandle(time(i-1) + h, 	w.col(i-1) + k3, u, xNminus, time);
	   
		w.col(i) = w.col(i-1) + (k1 + 2*(k2 + k3) + k4)/6;
	}
	
	return w;
}
*/