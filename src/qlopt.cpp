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

		// Get starting timepoint
		auto start = high_resolution_clock::now();

		//Use values from input struct and data info
		size_t n = params.gen.numOfStates;
		size_t m = u0.size();
		size_t lt = (int)((params.dat.endTime-params.dat.initialTime)
                  /params.dat.timeIncrement) + 1;
		std::cout << "lt  = " << lt << endl << endl;

		size_t ds = data.size();
		/**/vec ts(lt);
		ts(0) = params.dat.initialTime;
		for(size_t j=1; j<lt; j++)
		{
			ts(j) += ts(j-1) + params.dat.timeIncrement;
		}

		//TODO construct another more time dense things
		vec ly0(n*(m+1));
		ly0.fill(0);
		ly0.head(n) = y0;

		mat bob, robert, msmt, temp, U, A((ds+1)*m,m);
		mat I = mat::Identity(m, m);
		vec P((ds+1)*m), du(m);
		double O, objval, objval2, alpha = params.reg.alpha;
		double alpha2 = 0.0;

		OdeWrapper odewrapper(fun, input[0], u0);

		//Spline all the data sets
		std::vector<std::vector<thesis::spline>> spl_pairs(ds);
		for(size_t j=0; j<ds; j++){
			for(size_t i=0; i<n; i++){
				spl_pairs[j].push_back(thesis::spline(t, data[j].row(i)));
				//std::cout << spl_pairs[j][i].interpolate(80) << endl;
			}
		}

		/*for(size_t j=0; j<ds; j++){
			cout << "**** Data Set "<< j+1 << " ****\n";
			for(size_t i=0; i<n; i++){
					cout << "state " << i+1 << ": ";
				for(size_t k=0; k<lt; k++){
					std::cout << spl_pairs[j][i].interpolate(ts(k)) << "\t";
				}
				cout << endl <<endl;
			}
		}exit(0);*/

		results.ufinal = u0;
		results.uvals = u0;
		results.alpha(0);
		results.objval(0);
		results.iterations = 0;
		std::cout << "u  = " << results.uvals.transpose() << endl << endl;
		//TODO improve U mat stuff
		//std::vector<std::vector<vec>> Us(n, std::vector<vec>(m));

		cout.precision(7);
		A.fill(0.0); P.fill(0.0);
		for(size_t k = 0; k < params.gen.divisions; k++)
		{
			for(size_t j=0; j<params.tol.maxiter; j++)
			{
				results.iterations++;
				std::cout << "iteration = " << results.iterations << endl << endl;

				O = 0.0;

				for(size_t i = 0; i < ds; i++)
				{
					odewrapper.setControl(input[i]);
					odewrapper.setParameter(results.ufinal);
					odewrapper.setPreviousIteration(spl_pairs[i]);
					if(params.gen.finitediff)
					{
						robert = OdeIntWrapper(odewrapper, ly0, ts);
						bob = robert.bottomRows(n+n*m);
						lt = bob.cols();
						ts = robert.topRows(1).transpose();
					}else{
						//bob = qlRungeKutta4(*env->ode, (*env->time), uNot, y0, ext_data);
					}

					U = reshape(bob.bottomRows(n*m), m, n*lt);
					msmt.resize(n,lt);
					for(size_t j=0; j<n; j++){
						for(size_t k=0; k<lt; k++){
							msmt(j,k) = spl_pairs[i][j].interpolate(ts(k));
						}
					}
					temp = reshape(msmt - bob.topRows(n), 1, n*lt).row(0).transpose();

					//cout << U << endl << endl;
					//cout << temp.transpose() << endl << endl;

					A.middleRows(i*m, m) = findA(ts, U, m);
					//cout << "Condition numbers:" << endl;
					//cout << (A.middleRows(i*m, m)) << endl;

					P.segment(i*m, m) = findP(ts, U, temp, m);
					//cout << (P.segment(i*m, m)) << endl;

					O += findO(ts, temp);

					if(std::isnan(P.norm()) || std::isnan(A.norm())){
						std::cerr << "Termination: du is NaN." << endl;
						exit(1);
					}

					switch(params.reg.type)
					{
						case 0: alpha = 0.0;
								break;

						case 1: alpha = params.reg.alpha;
								break;

						case 2: alpha = alpha2 = params.reg.alpha;
								break;

						case 3: alpha = findAlpha(A, P);
								break;

						case 4: alpha = params.reg.alpha*O;
								break;

						case 5: alpha = findAlpha2(A, P, params.reg.alpha);
								break;

						default: cerr << "Chose a regularization option 0-5." << endl;
								exit(1);
					}
				}

				A.bottomRows(m) = alpha*I;
				du = A.colPivHouseholderQr().solve(P);
				results.alpha.conservativeResize(results.alpha.size()+1);
				results.alpha(j) = alpha;
				results.objval.conservativeResize(results.objval.size()+1);
				results.objval(j) = O;
				results.deltau.conservativeResize(results.deltau.size()+1);
				results.deltau(j) = du.norm();

				/*cout << "\tO = " << O << endl << endl;
				cout << "\tO/m = " << O/m << endl << endl;
				*/cout << "\talpha = " << alpha << endl << endl;
				cout << "\t||du|| = " << norm(du) << endl << endl;
				cout << "\talpha/||du|| = " << alpha/norm(du) << endl << endl;
				std::cout << "\tdu = " << du.transpose() << endl << endl;

				results.ufinal += du;
				results.uvals.conservativeResize(NoChange, results.uvals.cols()+1);
				results.uvals.col(results.uvals.cols()-1) = results.ufinal;

				std::cout << "\tu  = " << results.ufinal.transpose() << endl << endl;

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
				}

				if(j > 0){
					objval2 = O;
						for(size_t i = 0; i < ds; i++)
						{
							objval2 += du.transpose()*(A.middleRows(i*m, m)*du
								- 2*P.segment(i*m, m));
						}
					if(abs(objval - objval2) < params.tol.absobj){
						std::cout << "Termination: absolute objective function value tolerance." << endl;
						std::cout << "|J_new - J_old| = " << abs(objval - objval2) << endl;
						break;
					}
					objval = objval2;
				}else{
					objval = O;
					for(size_t i = 0; i < ds; i++)
					{
						objval += du.transpose()*(A.middleRows(i*m, m)*du
							- 2*P.segment(i*m, m));
					}
				}
			}
		}
		results.uvals.conservativeResize(NoChange, results.uvals.cols()+1);
		// Get ending timepoint
    	auto end = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(end - start);

		cout << "Computational Time:" << duration.count() << " ms" << endl;
		return results;
	}

	double findAlpha(mat A, vec P){
		int m = A.cols();
		int ds = A.rows()/m;

		double alpha = 1.0E2, atemp;



		mat I;
		I = mat::Identity(m,m);
		A.bottomRows(m) = alpha*I;

		vec du(m);
		du = A.colPivHouseholderQr().solve(P);

		double obj1=0.0, obj2=0.0;
		for(int i=0; i<ds-1; i++){
			obj1 += du.transpose()*(A.middleRows(i*m, m)*du
				- 2*P.segment(i*m, m));
		}

		for(int i = -6; i<2; i++){
			atemp = pow(10,i);
			A.bottomRows(m) = atemp*I;
			du = A.colPivHouseholderQr().solve(P);
			for(int i=0; i<ds-1; i++){
				obj2 += du.transpose()*(A.middleRows(i*m, m)*du
					- 2*P.segment(i*m, m));
			}
			if(obj1 > obj2){
				alpha = atemp;
				obj1 = obj2;
			}
			obj2 = 0.0;
		}
		return alpha;
	}

	double findAlpha2(mat A, vec P, const double max){
		int m = A.cols();
		double alpha;
		mat I;
		I = mat::Identity(m, m);
		vec du(m);

		for(int i = -12; i<3; i++){ // -7 -> 1
			alpha = pow(10,i);
			A.bottomRows(m) = alpha*I;
			du = A.colPivHouseholderQr().solve(P);

			if(du.norm() < max){
				break;
			}
		}
		return alpha;
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

		int N = 1000;
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
