#include <thesis_functions.h>

int thesis::nonlinearOdes::counter = 0;

mat findA(const vec& t, const mat& U, int m)
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

vec findP(const vec& t, const mat& U, const vec& dx, int m)
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

output findActualParam(soln_env *env, const int regs=0, const int numdivs = 1)
{
	int n = (*env->nth_soln).rows();
	int m = (*env->initial_params).size();
	int lt = (*env->time).size();
	
	int divs = (int) (lt/numdivs+1);
	double TOL = 1E-7;
	double lambda = *env->lambda;

	const mat measurements = *env->nth_soln;

	vec uNot(m);
	uNot = *env->initial_params;
	mat bob;
	mat U(n*m, n*lt);
	mat A(m,m), AT(m,m);
	mat I = mat::Identity(m, m);
	/*for(int i=0; i<m-1; i++){
        I(i+1,i) = -1;
    }
	I = I.transpose()*I;*/
	
	vec P(m);
	vec du(m), u1(m), u2(m);	

 	double O, objval, objval2;
	
	thesis::nonlinearOdes::setCounter(0);
	int LIMIT = 500;
	bool FD = true;

	std::vector<thesis::spline> ext_data;
	for(int i = 0; i < n; i++){
		ext_data.push_back(thesis::spline((*env->time), measurements.row(i))) ;
	}
	
	int i;
	output results = {0, 0, uNot};
	vec temp;
	for(int j = 0; j<numdivs; j++)
	{
		if(j == numdivs-1)
		{
			for(i = 0; i<LIMIT; i++)
			{	
				if(FD){ 
					bob = qLinearRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, ext_data);
				}else{
					bob = qlRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, ext_data);
				}
				
				U = reshape(bob.bottomRows(n*m), m, n*lt);
				*env->nth_soln = bob.topRows(n);
				//cout << "bob = \n" << bob << endl;
				A = findA(*env->time, U, m);
				//cout << "C = \n" << corrMat(A) << endl;
				cout << "rcond(A) = " << rcond(A) << endl;
				if(rcond(A) < 2.2E-10)
					cout << "A is near to being singular." << endl;
					
				temp = reshape(measurements - *env->nth_soln, 1, n*lt).row(0);
				P = findP(*env->time, U, temp, m);
				O = findO(*env->time, temp);
				
				//cout <<"\nDeterminant(A) = " << A.determinant() << endl;
				///cout << "rank(A) = " << A.fullPivHouseholderQr().rank() <<endl;

				
				if(regs == 0)
				{
					du = A.inverse()*P;
				}
				else if(regs == 1)
				{
					du = inverse(A + lambda*I)*P;
				}
				else if(regs == 2)
				{
					du = inverse(A + lambda*I)*(P + lambda*I*((*env->u_guess) - uNot));
					
					cout << "rcond(A + lambda*I) = " << rcond(A + lambda*I) << endl;
					if(rcond(A + lambda*I) < 2.2E-10)
						cout << "A + lambda*I is near to being singular." << endl;
				}
				else if(regs == 3)
				{
					//du = fnnls(A, P+A*uNot);
					//cout << "du = \n" << uNot.transpose() + du.transpose() << endl;
				}
				else
				{
					log_err("Termination: No such regularization method.");
					exit(0);
				}
				
				if(regs < 3){
					uNot += du;
				}else{
					uNot = du;
				}
				
				//latexOutput(*env->nth_soln, uNot, i+1, " &");
				results.iterations++;
				if (std::isnan(du.norm())){
					log_err("Termination: value for parameter is NaN.");
					exit(0);
				}else if(du.norm() < TOL){
					log_info("Termination: absolute parameter value tolerance.");
					log_info(du.norm());
					break;
				}else if (i >= LIMIT-1){
					log_info("Termination: max iterations reached.");
				}else if(du.norm()/uNot.norm() < TOL){
					log_info("Termination: relative parameter value tolerance.");
					log_info(du.norm()/uNot.norm());
					break;
				}
				if(i > 0){
					objval2 = O - 2*P.transpose()*du + du.transpose()*A*du;
					if(abs(objval - objval2) < TOL){
						log_info("Termination: absolute objective function value tolerance.");
						log_info(abs(objval - objval2));
						break;
					}
					objval = objval2;
				}else{
					objval = O - 2*P.transpose()*du + du.transpose()*A*du;
				}
			}
		}else{
			for(int i = 0; i<LIMIT; i++)
			{
				if(FD){ 
					bob = qLinearRungeKutta4(*env->ode, (*env->time).head((j+1)*divs), \
						uNot, *env->initial_cond, ext_data);
				}else{
					bob = qlRungeKutta4(*env->ode, (*env->time).head((j+1)*divs), \
						uNot, *env->initial_cond,  ext_data);
				}
				
				U = reshape(bob.bottomRows(n*m), m, n*(j+1)*divs);
				A = findA((*env->time).head((j+1)*divs), U, m);
				P = findP((*env->time).head((j+1)*divs), U, reshape(measurements.leftCols((j+1)*divs) - bob.topRows(n), 1, n*(j+1)*divs).row(0), m);
				O = findO(*env->time, reshape(measurements - *env->nth_soln, 1, n*lt).row(0));

				//cout <<"\nDeterminant(A) = " << A.determinant() << endl;
				//cout << "rank(A) = " << A.fullPivHouseholderQr().rank() <<endl;
				
				if(regs == 0)
				{
					du = A.inverse()*P;
				}
				else if(regs == 1)
				{
					du = inverse(A + lambda*I)*P;
				}
				else if(regs == 2)
				{
					du = inverse(A + lambda*I)*(P + lambda*I*((*env->u_guess) - uNot));
					
					cout << "rcond(A + lambda*I) = " << rcond(A + lambda*I) << endl;
					if(rcond(A + lambda*I) < 2.2E-10)
						cout << "A + lambda*I is near to being singular." << endl;
				}
				else if(regs == 3)
				{
					//du = fnnls(A, P+A*uNot);
					//cout << "du = \n" << uNot.transpose() + du.transpose() << endl;
				}
				else
				{
					log_err("Termination: No such regularization method.");
					exit(0);
				}
				
				if(regs < 3){
					uNot += du;
				}else{
					uNot = du;
				}
				//latexOutput(*env->nth_soln, uNot, i+1, " &");
				results.iterations++;
				if (std::isnan(du.norm())){
					log_err("Termination: value for parameter is NaN.");
					exit(0);
				}else if(du.norm() < TOL){
					log_info("Termination: absolute parameter value tolerance.");
					log_info(du.norm());
					break;
				}else if (i >= LIMIT-1){
					log_info("Termination: max iterations reached.");
					log_info(i+1);
				}else if(du.norm()/uNot.norm() < TOL){
					log_info("Termination: relative parameter value tolerance.");
					log_info(du.norm()/uNot.norm());
					break;
				}
				if(i > 0){
					objval2 = O - 2*P.transpose()*du + du.transpose()*A*du;
					if(abs(objval - objval2) < TOL){
						log_info("Termination: absolute objective function value tolerance.");
						log_info(abs(objval - objval2));
						break;
					}
					objval = objval2;
				}else{
					objval = O - 2*P.transpose()*du + du.transpose()*A*du;
				}
			}
		}
	}
	cout << lt << endl;
	if(FD){ 
		bob = qLinearRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, ext_data);
	}else{
		bob = qlRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, ext_data);
	}
	O = sqrt(findO(*env->time, reshape(measurements - *env->nth_soln, 1, n*lt).row(0)));
	cout << "("	<< lt << "," << O << ")" << endl;
	//cout << "number of function evaluations: " << thesis::nonlinearOdes::getCounter() << endl;
	results.fevals = thesis::nonlinearOdes::getCounter();
	results.du = uNot;
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

double norm(const mat& M){
	return M.norm();//std::sqrt(M.cwiseAbs2().sum())  ;
}

mat inverse(const mat& M){
	return M.inverse();
}

mat corrMat(const mat& M){
	int n = M.rows();
	mat C(n,n);
	
	for(int i = 0; i<n; i++){
		for(int j = 0; j<n; j++){
			C(i,j) = M(i,j)/sqrt(M(i,i)*M(j,j));
		}	
	}
	
	return C;

}

double rcond(const mat& M){
	return 1.0/(matnorm1(M)*matnorm1(M.inverse()));
}

double matnorm1(const mat& M){
	return M.cwiseAbs().colwise().sum().maxCoeff();	
}

/*vec fnnls(const mat& A, const vec& B)
{
	size_t n = B.size();
	vec x(n), s;
	x.fill(0.0);
	s = x;
	std::stack<int> P;
	std::queue<int> R;
	for(int i=0; i<n; i++)
		R.push(i);
	vec w;
	w = B-A*x;
	
	int j;
	while(R.size() > 0 || w.maxCoeff() > 1E-7){
		j = std::max_element(w.begin(), w.end())
		
		
		s.head(P.back()) = A.block()*B.head(P.back())
		if(s.head(P.back()).minCoeff() <= 0){
			
		}
	}
	
	
	
	return x;
}

double argmax(const vec& v){
	
	return
}*/

