#include <numerical_integration.h>

mat rungekutta4(string fname, const vec& time, const vec& u, const vec& yNot){
	thesis::nonlinearOdes no;
	sys fhandle = no.odeFuncMap[fname];

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

/*mp_mat mp_rungekutta4(string fname, const mp_vec& time, const mp_vec& u, const mp_vec& yNot){
	thesis::mpnonlinearOdes no;
	mp_sys fhandle = no.mpOdeFuncMap[fname];

	int N = time.size();

	//number of equations
	int m = yNot.size();

	//timesteps
	mpreal h = time(1)-time(0);

	//Init stuff
	mp_mat w(m, N);
	w.fill(0);
	w.col(0) = yNot;

	mp_vec k1(m), k2(m), k3(m), k4(m);

	for (int i = 0; i<N-1; i++)
	{
       k1 = h*fhandle(time(i), 		 w.col(i), 		  u);
       k2 = h*fhandle(time(i) + h/2, w.col(i) + k1/2, u);
       k3 = h*fhandle(time(i) + h/2, w.col(i) + k2/2, u);
       k4 = h*fhandle(time(i) + h, 	 w.col(i) + k3,   u);

	   w.col(i+1) = w.col(i) + (k1 + 2*(k2 + k3) + k4)/6;
	}

	return w;
}*/

mat qLinearRungeKutta4(string fname, const vec& time, const vec& u, const vec& yNot, std::vector<thesis::spline> Xn)
{
	thesis::nonlinearOdes no;
	sys fhandle = no.odeFuncMap[fname];

	int N = time.size();
	int n = Xn.size();

	//number of equations
	int m = yNot.size();

	//timestep
	double h = time(1)-time(0);

	//init stuff
	mat w(m, N);
	w.fill(0);
	w.col(0) = yNot;

	vec k1(m), k2(m), k3(m), k4(m);

	for (int i = 0; i<N-1; i++)
	{

		k1 = h*qlinear(fhandle, time(i),       w.col(i),        u, Xn);
		k2 = h*qlinear(fhandle, time(i) + h/2, w.col(i) + k1/2, u, Xn);
		k3 = h*qlinear(fhandle, time(i) + h/2, w.col(i) + k2/2, u, Xn);
		k4 = h*qlinear(fhandle, time(i) + h,   w.col(i) + k3,   u, Xn);

		//cout <<"(" << i <<")\nk1"<< k1 << "\nk2:" << k2 << "\nk3:" << k3 << "\nk4:" << k4 <<endl;
		w.col(i+1) = w.col(i) + (k1 + 2*(k2 + k3) + k4)/6;
	}

	return w;
}

mat qlOdeInt(string fname, const vec& time, const vec& u, const vec& yNot, std::vector<thesis::spline> xNminus)
{
	thesis::nonlinearOdes no;
	qsys fhandle = no.qLinFuncMap[fname+"_linearization"];
	//linearized_system_wrapper(fhandle);
	int N = time.size();
	//int n = xNminus.size();

	//number of equations
	int m = yNot.size();
	
	
}

mat qlRungeKutta4(string fname, const vec& time, const vec& u, const vec& yNot, std::vector<thesis::spline> xNminus)
{
	thesis::nonlinearOdes no;
	qsys fhandle = no.qLinFuncMap[fname+"_linearization"];

	int N = time.size();

	//number of equations
	int m = yNot.size();

	//timestep
	double h;

	//init stuff
	mat w(m, N);
	w.fill(0);
	w.col(0) = yNot;

	vec k1(m), k2(m), k3(m), k4(m);

	for (int i = 0; i<N-1; i++)
	{
		h  = time(i+1)-time(i);
		k1 = h*fhandle(time(i-1),		w.col(i-1), u, xNminus);
		k2 = h*fhandle(time(i-1) + h/2, w.col(i-1) + k1/2, u, xNminus);
		k3 = h*fhandle(time(i-1) + h/2, w.col(i-1) + k2/2, u, xNminus);
		k4 = h*fhandle(time(i-1) + h, 	w.col(i-1) + k3, u, xNminus);

		//cout <<"(" << i <<")\nk1"<< k1 << "\nk2:" << k2 << "\nk3:" << k3 << "\nk4:" << k4 <<endl;
		w.col(i+1) = w.col(i) + (k1 + 2*(k2 + k3) + k4)/6;
	}

	return w;
}

/*double simpson(const vec& t, const vec& x)
{
	double temp;

	int N = t.size();
	int end = N-1;

    double area = x(0) + x(end);
    for(int i = 1; i<N-1; i++)
	{
		temp = t(i);
        if((i+1)%2 == 0)
		{
            area += 2*temp;
        }
		else
		{
            area += 4*temp;
		}
    }
	return (t(end)- t(0))/(N-1)/3 * area;
}*/

double simpson(const vec& t, const vec& x)
{
	thesis::spline sim(t,x);
	double temp;

	int N = 3001;
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

struct curve_params { vec t; vec x;};

double f (double t, void * params) {
	struct curve_params * p = (struct curve_params *) params;
	vec tvec = (p->t);
	vec x = (p->x);
	thesis::spline sim(tvec,x);
	return sim.interpolate(t);
}

double gsl_integration(const vec& t, const vec& x)
{
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	struct curve_params params = {t,x};
	double result, error;
	gsl_function F;
	F.function = &f;
	F.params = &params;

	gsl_integration_qags (&F, t(0), t(t.size()-1), 0, 1e-7, 1000, w, &result, &error);
	gsl_integration_workspace_free (w);

	return result;
}

mat der(const mat& dx, const double& dt){
	mat ans(dx.size(),1);
	ans << dx/dt;
	return ans;
}

mat jac(sys f, double t, const mat& x, const mat& u, const double& h){
	int n = x.size();
	mat dx(n,n);
	dx << mat::Identity(n,n)*h;
	mat fprime(n,n);
	for(int j=0; j<n; j++)
	{
		fprime.col(j) = -f(t, x+2*dx.col(j), u) + 8*f(t, x+dx.col(j), u) - 8*f(t, x-dx.col(j), u) + f(t, x-2*dx.col(j), u);
	}

	return fprime/(12*h);
}

mat qlinear(sys fhandle, const double& t, const vec& x, const vec& u, std::vector<thesis::spline> Xn)
{
	int m = u.size();
	int n = Xn.size();
	double step = 5E-6;

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
	ans.block(0, 0, n, 1) = fx + jac(fhandle, t, xn1, u, step)*dxn;

	/*
	ans.block(0, 0, n, 1) = fx;
	for(int j=0; j<n; j++)
	{
		ans.block(0, 0, n, 1) += der(fhandle(t, xn1+dx.col(j), u) - fx, step)*dxn(j);
	}
	*/

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
			//ans.block(ind, 0, n, 1) += der(fhandle(t, xn1+dx.col(k), u+dun.col(j)) - dfdu - dfdx + fx, step*step)*dxn(k); //phi_ij
			ans.block(ind, 0, n, 1) += der(
				  fhandle(t, xn1+dx.col(k), u+dun.col(j))
				- fhandle(t, xn1+dx.col(k), u-dun.col(j))
				- fhandle(t, xn1-dx.col(k), u+dun.col(j))
				+ fhandle(t, xn1-dx.col(k), u-dun.col(j)), 4*step*step)*dxn(k); //phi_ij
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
