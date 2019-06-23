#include <qlopt.h>

namespace thesis{
	outputStruct qlopt(odefunction fun,
		const vec& t,
		const vec& u0,
		const vec& uguess,
		const vec& y0,
		const vector<vec>& input,
		const vector<mat>& data,
		const inputStruct& params,
		const vec& u)
	{
		outputStruct results;

		// Get starting timepoint
		auto start = high_resolution_clock::now();

		//Use values from input struct and data info
		size_t n = y0.size();
		size_t m = u0.size();
		size_t lt = t.size();
		size_t ds = data.size();

		//TODO construct another more time dense things
		vec ly0 = vec::Zero(n*(m+1));
		ly0.head(n) = y0;

		mat bob, robert, temp, U, A((ds+1)*m,m);
		mat SumA, Aplus;
		vec SumP, Pplus;
		mat I = mat::Identity(m, m);
		vec P((ds+1)*m), du(m), intmd;
		double O, objval, objval2, Oaddl, alpha = params.reg.alpha, tol;

		OdeWrapper odewrapper(fun, input[0], u0);

		//Spline all the data sets
		std::vector<std::vector<thesis::spline>> spl_pairs(ds);
		for(size_t i=0; i<ds; i++){
		    for(size_t j=0; j<n; j++){
		        //spl_pairs[j].push_back(thesis::spline(t, splinterSpline(t, data[j].row(i), params.noise.regParam) ) );
		        spl_pairs[i].push_back(thesis::spline(t, data[i].row(j)));
		    }
		}

		/******************************************************
			Enter first value for all the result containers
		******************************************************/

		results.ufinal = u0;
		results.uvals = u0;
		results.alpha(0);
		results.objval(0);
		results.omegaval(0);
		results.iterations = 0;

		std::cout << "u0 = " << results.uvals.transpose() << endl << endl;

		cout.precision(9);

		for(size_t k = 0; k < params.gen.divisions; k++)
		{
			for(size_t j=0; j < params.tol.maxiter; j++)
			{
				results.iterations++;
				std::cout << "iteration = " << results.iterations << endl << endl;

				/******************************
					Construct As, Ps and Os
				******************************/
				O = 0.0;
				odewrapper.setParameter(results.ufinal);

				for(size_t i = 0; i < ds; i++)
				{
					odewrapper.setControl(input[i]);
					odewrapper.setPreviousIteration(spl_pairs[i]);

					bob =  OdeIntWrapper(odewrapper, y0, t, tol);

					U = reshape(bob.bottomRows(n*m), m, n*lt);
					temp = reshape(data[i] - bob.topRows(n), 1, n*lt).row(0).transpose();
					A.middleRows(i*m, m) = findA(t, U, m);
					//cout << A.middleRows(i*m, m)  << endl << endl;
					P.segment(i*m, m) = findP(t, U, temp, m);
					O += findO(t, temp);

					if(std::isnan(P.norm()) || std::isnan(A.norm())){
						std::cerr << "Termination: du is NaN." << endl;
						exit(1);
					}
				}

				/*****************************************
					Determine regularization parameter
					And find the corresponding du
				*****************************************/

				switch(params.reg.type)
				{
					case 0: alpha = 0.0;
							A.bottomRows(m) = alpha*I;
							du = A.colPivHouseholderQr().solve(P);
							Oaddl = alpha*du.norm();
							break;

					case 1: alpha = params.reg.alpha;
							A.bottomRows(m) = alpha*I;
							du = A.colPivHouseholderQr().solve(P);
							Oaddl = alpha*du.norm();
							break;

					case 2: alpha  = params.reg.alpha*pow(O,2);
							A.bottomRows(m) = alpha*I;
							SumA = A.topRows(m);
							SumP = P.head(m);
							for(size_t  i = 1; i <= ds; i++)
							{
								SumA = SumA + A.middleRows(i*m, m);
								SumP = SumP + P.segment(i*m,m);
							}
							Aplus = SumA + alpha*I;
				            Pplus = SumP + alpha*(uguess-results.ufinal);
							du = Aplus.colPivHouseholderQr().solve(Pplus);

							intmd = uguess + du - results.ufinal;
							Oaddl = ds*alpha*intmd.transpose()*intmd;
							break;

					case 3: alpha  = findAlpha2(A, P, O, lt*ds);
							A.bottomRows(m) = alpha*I;
							du = A.colPivHouseholderQr().solve(P);
							Oaddl = alpha*du.norm();
							break;

					case 4: alpha = params.reg.alpha*pow(O,2);
							//cout << "before alpha" << results.ufinal.transpose() << endl;
							//cout << endl << uguess << endl;
							findGamma(A, P, results.ufinal, u);
							A.bottomRows(m) = alpha*I;
							du = A.colPivHouseholderQr().solve(P);
							Oaddl = alpha*du.norm();
							break;

					case 5: //alpha  = findAlpha3(A, P, O, uguess-results.ufinal, lt*ds);
							alpha = params.reg.alpha*pow(O,2);
							A.bottomRows(m) = alpha*I;
							SumA = A.topRows(m);
							SumP = P.head(m);
							for(size_t  i = 1; i <= ds; i++)
							{
								SumA = SumA + A.middleRows(i*m, m);
								SumP = SumP + P.segment(i*m,m);
							}
							Aplus = SumA + alpha*I;
				            Pplus = SumP + alpha*(uguess-results.ufinal);
							du = Aplus.colPivHouseholderQr().solve(Pplus);

							intmd = uguess + du - results.ufinal;
							Oaddl = ds*alpha*intmd.transpose()*intmd;
							break;

					case 6: alpha = findGamma(A, P, results.ufinal, u);
							A.bottomRows(m) = alpha*I;
							du = A.colPivHouseholderQr().solve(P);
							Oaddl = alpha*du.norm();
							break;

					case 7: alpha = findGamma2(A, P, results.ufinal, uguess, u);
							A.bottomRows(m) = alpha*I;
							SumA = A.topRows(m);
							SumP = P.head(m);
							for(size_t  i = 1; i <= ds; i++)
							{
								SumA = SumA + A.middleRows(i*m, m);
								SumP = SumP + P.segment(i*m,m);
							}
							Aplus = SumA + alpha*I;
				            Pplus = SumP + alpha*(uguess-results.ufinal);
							du = Aplus.colPivHouseholderQr().solve(Pplus);

							intmd = uguess + du - results.ufinal;
							Oaddl = ds*alpha*intmd.transpose()*intmd;
							break;

					case 8: alpha = findAlpha5(A, P, results.ufinal, odewrapper, data,
										input, ly0, t, spl_pairs);
							A.bottomRows(m) = alpha*I;
							du = A.colPivHouseholderQr().solve(P);
							Oaddl = alpha*du.norm();
							break;

					default: cerr << "Chose a regularization option 0-8." << endl;
							exit(1);
				}

				/****************************************
					Collect Results
				****************************************/

				results.alpha.conservativeResize(results.alpha.size()+1);
				results.alpha(j) = alpha;
				results.omegaval.conservativeResize(results.omegaval.size()+1);
				results.omegaval(j) = sqrt(O);
				results.deltau.conservativeResize(results.deltau.size()+1);
				results.deltau(j) = du.norm();

				results.ufinal += du;
				results.uvals.conservativeResize(NoChange, results.uvals.cols()+1);
				results.uvals.col(results.uvals.cols()-1) = results.ufinal;

				objval = O;
				for(size_t i = 0; i < ds; i++)
				{
					objval += du.transpose()*(A.middleRows(i*m, m)*du
						- 2*P.segment(i*m, m));
				}
				objval = objval + Oaddl;
				results.objval.conservativeResize(results.objval.size()+1);
				results.objval(j) = objval;

				/****************************************
					CLI output
				****************************************/

				std::cout << "\tdu = " << du.transpose() << endl << endl;
				std::cout << "\tu  = " << results.ufinal.transpose() << endl << endl;
				std::cout << "\talpha = " << alpha << endl << endl;

				std::cout << "Normed difference of data and model." << endl;
				std::cout << "\t\\Sigma||x-x_N|| = " << sqrt(O) << endl << endl;
				std::cout << "\t\\Sigma||x-x_N||/ds = " << sqrt(O)/ds << endl << endl;
				std::cout << "Absolute parameter value tolerance." << endl;
				std::cout << "\t||du|| = " << norm(du) << endl << endl;
				std::cout << "Relative parameter value tolerance." << endl;
				std::cout << "\t||du||/||u||= " << du.norm()/results.ufinal.norm() << endl;
				std::cout << "Objective function value." << endl;
				std::cout << "\tJ(du) = " << objval << endl;

				/****************************************
					Check the termination conditions
				****************************************/

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
				}else if(du.norm()/u0.norm() < params.tol.relparam){
					std::cout << "Termination: relative parameter value tolerance." << endl;
					break;
				}else if(du.norm() < params.tol.absparam){
					std::cout << "Termination: absolute parameter value tolerance." << endl;
					break;
				}

				if(j > 0){
					std::cout << "Absolute objective function value tolerance." << endl;
					std::cout << "\t|J_new - J_old| = " << abs(objval - objval2) << endl;
					std::cout << "Relative objective function value tolerance." << endl;
					std::cout << "\t|J_new - J_old|/J_new = " << abs(objval - objval2)/objval << endl;

					if(abs(objval - objval2) < params.tol.absobj){
						std::cout << "Termination: absolute objective function value tolerance." << endl;
						break;
					}else if(abs(objval - objval2)/objval < params.tol.relobj){
						std::cout << "Termination: relative objective function value tolerance." << endl;
						break;
					}
					objval2 = objval;
				}else{
					objval2 = objval;
				}
			}
		}
		results.uvals.conservativeResize(NoChange, results.uvals.cols()+1);
		// Get ending timepoint
    	auto end = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(end - start);

		cout << "Computational Time:" << duration.count()/1E6 << " s" << endl;
		return results;
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

		//return ans*(time(1)-time(0));
		return simpson(time, aij);
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
