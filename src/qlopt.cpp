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

		mat bob, robert, temp, U, A((ds+1)*m,m);
		std::vector<mat> msmt(ds);
		mat I = mat::Identity(m, m);
		vec P((ds+1)*m), du(m);
		double O, objval, objval2, alpha = params.reg.alpha;
		double alpha2 = 0.0;

		OdeWrapper odewrapper(fun, input[0], u0);

		//Spline all the data sets
		std::vector<std::vector<thesis::spline>> spl_pairs(ds);
		for(size_t i=0; i<ds; i++){
		    msmt[i].resize(n,lt);
		    for(size_t j=0; j<n; j++){
		        //spl_pairs[j].push_back(thesis::spline(t, splinterSpline(t, data[j].row(i), params.noise.regParam) ) );
		        spl_pairs[i].push_back(thesis::spline(t, data[i].row(j)));
		        for(size_t k=0; k<lt; k++){
		            msmt[i](j,k) = spl_pairs[i][j].interpolate(ts(k));
		        }
		    }
		}

		results.ufinal = u0;
		results.uvals = u0;
		results.alpha(0);
		results.objval(0);
		results.iterations = 0;
		std::cout << "u0 = " << results.uvals.transpose() << endl << endl;
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
				odewrapper.setParameter(results.ufinal);
				for(size_t i = 0; i < ds; i++)
				{
					odewrapper.setControl(input[i]);
					odewrapper.setPreviousIteration(spl_pairs[i]);
					if(params.gen.finitediff)
					{
						robert = OdeIntWrapper(odewrapper, ly0, ts);
						bob = robert.bottomRows(n+n*m);
						lt = bob.cols();
						ts = robert.topRows(1).transpose();
					}else{
					}

					U = reshape(bob.bottomRows(n*m), m, n*lt);
					temp = reshape(msmt[i] - bob.topRows(n), 1, n*lt).row(0).transpose();
					A.middleRows(i*m, m) = findA(ts, U, m);
					P.segment(i*m, m) = findP(ts, U, temp, m);
					O += findO(ts, temp);

					if(std::isnan(P.norm()) || std::isnan(A.norm())){
						std::cerr << "Termination: du is NaN." << endl;
						exit(1);
					}
				}

				switch(params.reg.type)
				{
					case 0: alpha = 0.0;
							break;

					case 1: alpha = params.reg.alpha;
							break;

					case 2: alpha = alpha2 = params.reg.alpha;
							break;

					case 3: //alpha = findAlpha(A, P, O, lt*ds);
							alpha = gcv(A, P, results.ufinal, odewrapper, msmt, input, y0, ts, spl_pairs);
							break;

					case 4: //if(O > 0.2)
								alpha = params.reg.alpha*pow(O,2);
							//else alpha = params.reg.alpha*O;
							break;

					case 5: alpha = findAlpha2(A, P, params.reg.alpha);
							break;

					case 6: alpha = findGamma(A, P, O, results.ufinal, uguess);
							break;

					default: cerr << "Chose a regularization option 0-6." << endl;
							exit(1);
				}

				A.bottomRows(m) = alpha*I;
				du = A.colPivHouseholderQr().solve(P);

				results.alpha.conservativeResize(results.alpha.size()+1);
				results.alpha(j) = alpha;
				results.omegaval.conservativeResize(results.omegaval.size()+1);
				results.omegaval(j) = sqrt(O);
				results.deltau.conservativeResize(results.deltau.size()+1);
				results.deltau(j) = du.norm();

				results.ufinal += du;
				results.uvals.conservativeResize(NoChange, results.uvals.cols()+1);
				results.uvals.col(results.uvals.cols()-1) = results.ufinal;

				objval = O + params.reg.alpha*du.transpose()*du ;
				for(size_t i = 0; i < ds; i++)
				{
					objval += du.transpose()*(A.middleRows(i*m, m)*du
						- 2*P.segment(i*m, m));
				}
				results.objval.conservativeResize(results.objval.size()+1);
				results.objval(j) = objval;

				std::cout << "\tdu = " << du.transpose() << endl << endl;
				std::cout << "\tu  = " << results.ufinal.transpose() << endl << endl;

				cout << "\t||x-x_N|| = " << sqrt(O) << endl << endl;
				cout << "\t||x-x_N||/ds = " << sqrt(O)/ds << endl << endl;
				cout << "\talpha = " << alpha << endl << endl;
				cout << "\t||du|| = " << norm(du) << endl << endl;
				std::cout << "\t||du||/||u||= " << du.norm()/results.ufinal.norm() << endl;
				std::cout << "\tJ(du) = " << objval << endl;

				// Check the termination conditions
				if (std::isnan(du.norm())){
					std::cerr << "Termination: value for parameter is NaN." << endl;
					exit(0);
				}else if (j >= params.tol.maxiter-1){
					std::cout << "Termination: max iterations reached." << endl;
				}else if (sqrt(O) < params.tol.normdiff){
					std::cout << "Termination: Normed difference of data and model." << endl;
					break;
				}else if (objval < params.tol.objval){
					std::cout << "Termination: Objective function value." << endl;
					break;
				}

				/*else if(du.norm()/u0.norm() < params.tol.relparam){
					std::cout << "Termination: relative parameter value tolerance." << endl;
					std::cout << "du/u = " << du.norm()/results.ufinal.norm() << endl;
					break;
				}else if(du.norm() < params.tol.absparam){
					std::cout << "Termination: absolute parameter value tolerance." << endl;
					std::cout << "du = "<< du.norm() << endl;
					break;
				}*/

				/*if(j > 0){
					objval2 = O + params.reg.alpha*du.transpose()*du ;
					for(size_t i = 0; i < ds; i++)
					{
						objval2 += du.transpose()*(A.middleRows(i*m, m)*du
							- 2*P.segment(i*m, m));
					}

					results.objval.conservativeResize(results.objval.size()+1);
					results.objval(j) = objval2;
					std::cout << "Objective function value." << endl;
					std::cout << "J_new = " << objval2 << endl;

					if(abs(objval - objval2) < params.tol.absobj){
						std::cout << "Termination: absolute objective function value tolerance." << endl;
						std::cout << "|J_new - J_old| = " << abs(objval - objval2) << endl;
						break;
					}

					else if(abs(objval - objval2)/objval2 < params.tol.relobj){
						std::cout << "Termination: relative objective function value tolerance." << endl;
						std::cout << "|J_new - J_old|/J_new = " << abs(objval - objval2)/objval2 << endl;
						break;
					}

					std::cout << "Absolute objective function value tolerance." << endl;
					std::cout << "|J_new - J_old| = " << abs(objval - objval2) << endl;
					std::cout << "Relative objective function value tolerance." << endl;
					std::cout << "|J_new - J_old|/J_new = " << abs(objval - objval2)/objval2 << endl;
					objval = objval2;
				}else{}*/
			}
		}
		results.uvals.conservativeResize(NoChange, results.uvals.cols()+1);
		// Get ending timepoint
    	auto end = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(end - start);

		cout << "Computational Time:" << duration.count()/1E6 << " s" << endl;
		return results;
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

	double findGamma(mat A, vec P, double O, vec uNot, vec u){
		int m = P.size();
		mat I = mat::Identity(m, m);
		double gamma, dg = .1, dvs = 10,
			hold = pow(10,-9),//*pow(O,2),
			nu;
		A.bottomRows(m) = hold*I;
		double nup = norm(uNot + A.colPivHouseholderQr().solve(P) - u);

		vec du(m), total(m);
		mat B(m,m);

		for(int j = -8; j<4; j++){
			gamma = pow(10,j);//*pow(O,2);
		    for(int k = 0; k<dvs; k++){
				A.bottomRows(m) = gamma*I;
				du = A.colPivHouseholderQr().solve(P);

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
				//dg += .1;
				gamma = pow(10,j+(k+1)*dg);//*pow(O,2);
		    }
		}
		cout << "alpha = " << hold << endl;
		//cout << "c = " << hold/pow(O,2) << endl << endl;
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
		double ans = 0.0;

		int mid = time.size();
		int n = u1.size()/mid;

		vec a = u1.array() * u2.array();

		vec aij = a.head(mid);

		if(mid != n*mid){
			for(int k=1; k<n; k++){
				aij += a.segment(k*mid, mid);
			}
		}

		for(int k=1; k<aij.size(); k++){
			ans += aij(k);
		}

		return ans*(time(1)-time(0));
		//return simpson(time, aij);
		//return gsl_integration(time, aij);
	}

	double simpson(const vec& t, const vec& x)
	{
		thesis::spline sim(t,x);
		double temp;

		int N = 5000;
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
