#include <qlopt.h>

namespace thesis{
	outputStruct qlopt(odefunction fun, 
		const vec& t, 
		const vec& u0,
		const vec& uguess,
		const vec& y0,
		const vector<vec>& input, 
		const vector<mat>& data, 
		const inputStruct& params)
	{
		outputStruct results;
		
		//Use values from input struct and data info
		size_t n = params.gen.numOfStates;
		size_t m = u0.size();
		size_t lt = (int)((params.dat.endTime-params.dat.initialTime)
                      /params.dat.timeIncrement) + 1;//t.size();
		std::cout << "lt  = " << lt << endl << endl;
		
		size_t ds = data.size();
		vec ts(lt);
		ts(0) = params.dat.initialTime;
		for(size_t j=1; j<lt; j++)
		{
			ts(j) += ts(j-1) + params.dat.timeIncrement;
		}

		//TODO construct another more time dense things
		vec ly0(n*(m+1));
		ly0.fill(0);
		ly0.head(n) = y0;

		mat bob, temp, U, A(m,m);
		mat I = mat::Identity(m, m);
		vec P(m), du(m);
		double O, objval, objval2, alpha = params.reg.alpha;
		OdeWrapper odewrapper(fun);
		
		//Spline all the data sets
		std::vector<std::vector<thesis::spline>> spl_pairs(ds);
		for(size_t j=0; j<ds; j++){	
			for(size_t i=0; i<n; i++){
				spl_pairs[j].push_back(thesis::spline(t, data[j].row(i)));
				//std::cout << spl_pairs[j][i].interpolate(1) << endl;
			}
		}
		results.ufinal = u0;
		results.uvals = u0;
		results.iterations = 0;
		std::cout << "u  = " << results.uvals << endl << endl;
		//TODO improve U mat stuff
		//std::vector<std::vector<vec>> Us(n, std::vector<vec>(m));
		
		cout.precision(7);
		for(size_t k = 0; k < 1; k++)
		{
			for(size_t i = 0; i < ds; i++)
			{
				odewrapper.setControl(input[i]);
			
				for(size_t j=0; j<params.tol.maxiter; j++)
				{
					results.iterations++;
					if(params.gen.finitediff){ 
						bob = qloptRungeKutta4(odewrapper, ts, results.ufinal, ly0, spl_pairs[i]);
					}else{
						//bob = qlRungeKutta4(*env->ode, (*env->time), uNot, y0, ext_data);
					}
			
					U = reshape(bob.bottomRows(n*m), m, n*lt);
					temp = reshape(data[i] - bob.topRows(n), 1, n*lt).row(0).transpose();
			
					A = findA(ts, U, m);
					P = findP(ts, U, temp, m);
					O = findO(ts, temp);
					
					if(std::isnan(P.norm()) || std::isnan(A.norm())){
						std::cerr << "Termination: du is NaN." << endl; 
						exit(0);
					}

					alpha  = findGamma(A, P, u0, uguess);
					du = inverse(A + alpha*I)*P;
				    cout << "alpha = " << alpha << endl;
				
					std::cout << "iteration = " << results.iterations << endl << endl;
					std::cout << "du = " << du.transpose() << endl << endl;
					
					results.ufinal += du;
					results.uvals.conservativeResize(NoChange, results.uvals.cols()+1); 
					results.uvals.col(results.uvals.cols()-1) = results.ufinal;
					
					std::cout << "u  = " << results.ufinal << endl << endl;
					//std::cout << "u  = " << results.uvals << endl << endl;
			
					// Check the termination conditions
					if (std::isnan(du.norm())){
						std::cerr << "Termination: value for parameter is NaN." << endl;
						exit(0);
					}else if (j >= params.tol.maxiter-1){
						std::cout << "Termination: max iterations reached." << endl;
					}else if(du.norm()/u0.norm() < params.tol.relparam){
						std::cout << "Termination: relative parameter value tolerance." << endl;
						std::cout << "du/u = " << du.norm()/u0.norm() << endl;
						break;
					}else if(du.norm() < params.tol.absparam){
						std::cout << "Termination: absolute parameter value tolerance." << endl;
						std::cout << "du = "<< du.norm() << endl;
						break;
					}/*else if(O < params.tol.absparam){
						std::cout << "Termination: absolute normed difference tolerance." << endl;
						std::cout << "|x- x_N|^2 = "<< O << endl;
						break;
					}*/
					//TODO to avoid calculating the Obj func val use O with reps ||x-xn||^2
					if(j > 0){
						objval2 = O 
							- 2*P.segment(i*m, m).transpose()*du 
							+ du.transpose()*A.middleRows(i*m, m)*du;
					
						if(abs(objval - objval2) < params.tol.absobj){
							std::cout << "Termination: absolute objective function value tolerance." << endl;
							std::cout << "|J_new - J_old| = " << abs(objval - objval2) << endl;
							break;
						}
						objval = objval2;
					}else{
						objval = O 
							- 2*P.segment(i*m, m).transpose()*du 
							+ du.transpose()*A.middleRows(i*m, m)*du;
					}
				}
				std::cout << "|u_n - u*| = " << norm(results.ufinal.transpose() - uguess.transpose()) << endl;
			}
		}
		results.uvals.conservativeResize(NoChange, results.uvals.cols()+1); 
		results.uvals.col(results.uvals.cols()-1) = uguess;
		std::cout << results.ufinal.transpose() << endl;
		parameterOutput(results.uvals);
		return results;
	}

	double findGamma(mat A, vec P, vec uNot, vec u){
		int m = P.size();
		mat I = mat::Identity(m, m);
		double gamma, dg, dvs = 9,
			hold = pow(10,-14),
			nu,
			nup = norm(uNot + inverse(A + hold*I)*P - u);

		vec du(m), total(m);
		mat B(m,m);
		
		for(int j = -14; j<0; j++){
		    gamma = pow(10,j);
		    dg = .1; //(pow(10,j)-pow(10,j-1))/dvs;
		    for(int k = 0; k<(dvs); k++){
		        B = A + gamma*I;
		        du = B.inverse()*P;

		        total = uNot + du - u;
				nu = total.norm();
		        cout << "("
	                   << gamma
	                   << ","
	                   << nu
	                   << ")"
	                   << endl;
	            if(nup > nu){
					nup = nu;
					hold = gamma;
				}
				dg += .1;
				gamma = pow(10,j+dg);
		    }
		}
		cout << "alpha = " << hold << endl;
		return hold;
	}

	mat qloptRungeKutta4(OdeWrapper& fhandle, const vec& time, 
		const vec& u, const vec& yNot, std::vector<thesis::spline>& Xn)
	{
		size_t N = time.size();

		//number of equations
		size_t m = yNot.size();

		//timestep
		double h = 0.0;

		//init stuff
		mat w(m, N);
		w.fill(0);
		w.col(0) = yNot;

		vec k1(m), k2(m), k3(m), k4(m);

		for (size_t i = 0; i<N-1; i++)
		{
			h = time(i+1) - time(i);
			k1 = h*qlinear(fhandle, time(i),       w.col(i),        u, Xn);
			k2 = h*qlinear(fhandle, time(i) + h/2, w.col(i) + k1/2, u, Xn);
			k3 = h*qlinear(fhandle, time(i) + h/2, w.col(i) + k2/2, u, Xn);
			k4 = h*qlinear(fhandle, time(i) + h,   w.col(i) + k3,   u, Xn);

			//cout <<"(" << i <<")\nk1"<< k1 << "\nk2:" << k2 << "\nk3:" << k3 << "\nk4:" << k4 <<endl;
			w.col(i+1) = w.col(i) + (k1 + 2*(k2 + k3) + k4)/6;
		}

		return w;
	}
	
	mat qlinear(OdeWrapper& fhandle, const double& t, const vec& x, 
		const vec& u, std::vector<thesis::spline>& Xn)
	{
		size_t m = u.size();
		size_t n = Xn.size();
		double step = sqrt(2.2E-16);

		vec xn1(n); // this is x_N-1
		vec dxn(n);

		for(size_t i=0; i<n; i++)
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

		mat dfdx(n,1);
		mat dfdu(n,1);

		// the Un part of t the linearization
		mat dun(m, m);// for u derivative
		dun << mat::Identity(m,m)*step;

		size_t ind;
		for(size_t j=0; j<m; j++)
		{
			ind = (j+1)*n;
			dfdu = fhandle(t, xn1, u+dun.col(j)*u(j));
			ans.block(ind, 0, n, 1) = der(dfdu - fx, step*u(j)); // df/du
			for(size_t k=0; k<n; k++){
				dfdx = fhandle(t, xn1+dx.col(k)*xn1(k), u);
				ans.block(ind, 0, n, 1) += der(dfdx  - fx, step*xn1(k))*x(ind+k);//J*Un

				ans.block(ind, 0, n, 1) += der(
					  fhandle(t, xn1+dx.col(k)*xn1(k), u+dun.col(j)*u(j))
					- fhandle(t, xn1+dx.col(k)*xn1(k), u-dun.col(j)*u(j))
					- fhandle(t, xn1-dx.col(k)*xn1(k), u+dun.col(j)*u(j))
					+ fhandle(t, xn1-dx.col(k)*xn1(k), u-dun.col(j)*u(j)), 
					4*step*xn1(k)*step*u(j))*dxn(k); //phi_ij
			}
		}

		return ans;
	}
	
	mat der(const mat& dx, const double& dt)
	{
		mat ans(dx.size(),1);
		ans << dx/dt;
		return ans;
	}

	mat jac(OdeWrapper& f, double t, const mat& x, const mat& u, 
		const double& h)
	{
		int n = x.size();
		mat dx(n,n);
		dx << mat::Identity(n,n)*h;
		mat fprime(n,n);
		for(int j=0; j<n; j++)
		{
			fprime.col(j) = (-f(t, x+2*dx.col(j)*x(j), u) + 8*f(t, x+dx.col(j)*x(j), u) - 8*f(t, x-dx.col(j)*x(j), u) + f(t, x-2*dx.col(j)*x(j), u))/x(j);
		}

		return fprime/(12*h);
	}
	
	mat reshape(const mat& U, int n, int m)
	{
		mat newU(n,m);
		newU.fill(0);
		int olt = U.row(0).size(); 	//old time(row) length
		int on = m/olt;				//on col length

		for(int i = 0; i<n; i++){
			for(int j = 0; j<on; j++){
				newU.block(i, j*olt, 1, olt) = U.row(i*on + j);
			}
		}

		return newU;
	}

	mat findA(const vec& t, const mat& U, const size_t& m)
	{
		mat A(m,m);

		for(size_t i = 0; i<m; i++){
			for(size_t j = 0; j<m; j++){
				if(i <= j){
					A(i,j) = innerProd(U.row(i), U.row(j), t);
					A(j,i) = A(i,j);
				}
			}
		}
		return A;
	}

	vec findP(const vec& t, const mat& U, const vec& dx, const size_t& m)
	{
		vec P(m);

		for(size_t i = 0; i<m; i++)
		{
			P(i) = innerProd(U.row(i), dx, t);
		}

		return P;
	}

	double findO(const vec& t, const vec& dx)
	{
		return innerProd(dx, dx, t);
	}

	double innerProd(const vec& u1, const vec& u2, const vec& time)
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
		//return gsl_integration(time, aij);
	}
	
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
	
	double norm(const mat& M)
	{
		return M.norm();
	}

	mat inverse(const mat& M)
	{
		return M.inverse();
	}

	mat corrMat(const mat& M)
	{
		int n = M.rows();
		mat C(n,n);
	
		for(int i = 0; i<n; i++){
			for(int j = 0; j<n; j++){
				C(i,j) = M(i,j)/sqrt(M(i,i)*M(j,j));
			}	
		}
	
		return C;
	}

	double rcond(const mat& M)
	{
		return 1.0/(matnorm1(M)*matnorm1(M.inverse()));
	}
	
	double matnorm1(const mat& M){
		return M.cwiseAbs().colwise().sum().maxCoeff();	
	}
}
