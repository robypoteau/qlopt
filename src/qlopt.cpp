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

		mat bob, temp, U, A((ds+1)*m,m), A1(m,m);
		mat I = mat::Identity(m, m);
		vec P((ds+1)*m), P1(m), du(m);
		double O, objval, objval2, alpha = params.reg.alpha;
		double gamma = 0.0, dg = 0.0, told, tnew;
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
		std::cout << "u  = " << results.ufinal.transpose() << endl << endl;
			
		//TODO improve U mat stuff
		//std::vector<std::vector<vec>> Us(n, std::vector<vec>(m));
		
		tallskinnyqr tsqr(ds*m, m, 300);
		A.fill(0.0); P.fill(0.0);
		cout.precision(10);
		for(size_t j=0; j<params.tol.maxiter; j++)
		{		
			//TODO Solution to quasilinear
			A1.fill(0.0), P1.fill(0.0), O = 0.0;
			for(size_t i = 0; i < ds; i++)
			{
				odewrapper.setControl(input[i]);
				if(params.gen.finitediff){ 
					bob = qloptRungeKutta4(odewrapper, ts, results.ufinal, ly0, spl_pairs[i]);
				}else{
					//bob = qlRungeKutta4(*env->ode, (*env->time), uNot, y0, ext_data);
				}
				
				U = reshape(bob.bottomRows(n*m), m, n*lt);
				temp = reshape(data[i] - bob.topRows(n), 1, n*lt).row(0).transpose();
				
				/*for(size_t q = 0; q < n ; q++){
					for(size_t w = 0; w < m; w++){ 
						Us[q][w] = bob.bottomRows(n*m).row(q*w);
					}
				}*/
				
				if(params.reg.type == 0){
					A.middleRows(i*m, m) = findA(ts, U, m);
					P.segment(i*m, m) = findP(ts, U, temp, m);
				}else if(params.reg.type == 1){
					A.middleRows(i*m, m) = findA(ts, U, m);
					P.segment(i*m, m) = findP(ts, U, temp, m);
				}else if(params.reg.type == 2){
					A1 += findA(ts, U, m);
					P1 += findP(ts, U, temp, m);
				}else if(params.reg.type == 3){
					A1 += findA(ts, U, m);
					P1 += findP(ts, U, temp, m);
				}else if(params.reg.type == 4){
					A.middleRows(i*m, m) = findA(ts, U, m);
					P.segment(i*m, m) = findP(ts, U, temp, m);
				}else if(params.reg.type == 5){
          A.middleRows(i*m, m) = findA(ts, U, m);
					P.segment(i*m, m) = findP(ts, U, temp, m);
          A1 += findA(ts, U, m);
					P1 += findP(ts, U, temp, m);
				}else if(params.reg.type == 6){
					A1 += findA(ts, U, m);
					P1 += findP(ts, U, temp, m);
				}
				
				O += findO(ts, temp);
			}
			//cout << P << endl;
			if(std::isnan(P.norm()) || std::isnan(P1.norm())){
				std::cerr << "Termination: du is NaN." << endl; 
				exit(0);
			}

			if(params.reg.type == 0)
			{
				std::cout << "rcond(A) = " << rcond(A.transpose()*A) << endl;
				std::cout << "rank(A) = " << A.fullPivHouseholderQr().rank() << endl;
				du = A.fullPivHouseholderQr().solve(P);
			}
			else if(params.reg.type == 1)
			{
				tsqr.update(A,P);
				du = tsqr.solve();
			}
			else if(params.reg.type == 2)
			{
				tsqr.update(A1,P1);
				du = tsqr.solve();
			}
			else if(params.reg.type == 3)
			{
				A1 += alpha*I;
				du = A1.inverse()*P1;
			}
			else if(params.reg.type == 4)
			{
        		alpha = 0.0;
				A.middleRows(ds*m, m) = alpha*I;
				//std::cout << "rcond(A) = " << rcond(A) << endl;
				std::cout << "rcond(A^T*A) = " << rcond(A.transpose()*A) << endl;
				std::cout << "rank(A) = " << A.fullPivHouseholderQr().rank() << endl;
				//std::cout << "corr(A) = " << corrMat(A) << endl;
				du = A.fullPivHouseholderQr().solve(P);

				//if(rcond(A.transpose()*A) < 1E-6){
				    told = norm(results.ufinal + du - uguess);
					for(int j=-7; j<1; j++){
					  gamma = pow(10,j);
					  dg = (pow(10,j+1)-pow(10,j))/9; //10 divisions
					  for(int k = 0; k<9; k++)
					    {
					      gamma += dg;

					      A.middleRows(ds*m, m) = gamma*I;
					      du = A.fullPivHouseholderQr().solve(P);

					      tnew = norm(results.ufinal + du - uguess);
					      if(tnew < told){
					        alpha = gamma;
					        told = tnew;
					        //goto BREAK;
					      }
					      cout << "("
					           << gamma
					           << ","
					           << norm(results.ufinal + du - uguess)
					           << ","
					           << rcond(A.transpose()*A)
					           << ","
					           << A.fullPivHouseholderQr().rank()
					           << ")"
					           << endl;
					    }
					}
				//}
				//BREAK:
				du = A.fullPivHouseholderQr().solve(P);
				cout << "alpha = " << alpha << endl;
				//exit(0);
			}
			else if(params.reg.type == 5) //tsqr
			{
				alpha  = findGamma(A.middleRows(0, m), P.segment(0,m), u0, uguess);
				du = inverse(A.middleRows(0, m) + alpha*I)*P.segment(0,m);
                cout << "alpha = " << alpha << endl;
			}
			else if(params.reg.type == 6)
			{
				du = inverse(A + pow(10,-12)*I)*P;
				told = norm(results.ufinal + du - 	uguess);
				
				for(int j=-12; j<2; j++){
					gamma = pow(10,j);
					dg = (pow(10,j+1)-pow(10,j))/9; //10 divisions
					for(int k = 0; k<10; k++)
					{
						du = inverse(A + gamma*I)*P;
						tnew = norm(results.ufinal + du - uguess);
						if(tnew - told < 1E-2){
							alpha = gamma;
              told = tnew;
						}
						//cout << gamma << "," << 
						//	norm(results.ufinal + du - uguess) << endl;
						gamma += dg;
					}
				}
				
				du = inverse(A + alpha*I)*P;
				cout << "alpha = " << alpha << endl;
        		cout << "dx = " << O << endl;
				//exit(0);
			}
			else
			{
				std::cerr << "Termination: No such regularization method." << endl;
				exit(0);
			}
			//std::cout << "A = \n" << A1 << endl;
			//std::cout << "P = \n" << P1 << endl;
			
			
			std::cout << "du = " << du.transpose() << endl << endl;
			results.ufinal += du;
			std::cout << "u  = " << results.ufinal.transpose() << endl << endl;
			
			results.iterations++;
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
			}else if(O < params.tol.absparam){
				std::cout << "Termination: absolute normed difference tolerance." << endl;
				std::cout << "|x- x_N|^2 = "<< O << endl;
				break;
			}
			//TODO to avoid calculating the Obj func val use O with reps ||x-xn||^2
			if(j > 0){
				objval2 = O;
				for(size_t i = 0; i < ds; i++){
					objval2 = objval2
						- 2*P.segment(i*m, m).transpose()*du 
						+ du.transpose()*A.middleRows(i*m, m)*du;
				}
				
				if(abs(objval - objval2) < params.tol.absobj){
					std::cout << "Termination: absolute objective function value tolerance." << endl;
					std::cout << "|J_new - J_old| = " << abs(objval - objval2) << endl;
					break;
				}
				objval = objval2;
			}else{
				objval = O;
				for(size_t i = 0; i < ds; i++){
					objval = objval 
						- 2*P.segment(i*m, m).transpose()*du 
						+ du.transpose()*A.middleRows(i*m, m)*du;
				}
			}
		}
		std::cout << results.ufinal.transpose() << endl;
		std::cout << results.ufinal.transpose() - uguess.transpose() << endl;
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

		for(int j = -7; j<2; j++){
		    gamma = 1*pow(10,j);
		    dg = (1*pow(10,j)-1*pow(10,j-1))/dvs;
		    for(int k = 0; k<(dvs+1); k++){
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
				gamma += dg;
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

		for(int i = 0; i<m; i++){
			for(int j = 0; j<m; j++){
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

		for(int i = 0; i<m; i++)
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
