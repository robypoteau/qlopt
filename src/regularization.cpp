#include <regularization.h>

//Used to find plots for the ||u_n + \delta u_n+1 - u*|| vs \aplha
//When the problem requires regularization
void regparamexp1a(soln_env *env, vec u, int brk){
	int n = (*env->nth_soln).rows();
	int m = (*env->initial_params).size();
	int lt = (*env->time).size();
	const mat measurements = *env->nth_soln;

	vec uNot(m);
	uNot = *env->initial_params;
	mat bob;
	mat U(n*m, n*lt);
	mat A(m,m), AT(m,m);
	mat B = mat::Identity(m, m);
	for(int i=0; i<m-1; i++){
		//B(i+1,i) = -1;
	}
	mat BT = B.transpose();
	vec P(m);
	vec du(m);
	double TOL = 0.00001;
	double gamma = 10.0;
	int LIMIT = 200;
	cout.precision(10);

	double inc = 0.1;

	for(int i = 0; i<LIMIT; i++)
	{
		bob = qLinearRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, measurements);
		U = reshape(bob.bottomRows(n*m), m, n*lt);
		*env->nth_soln = bob.topRows(n);
		A = findA(*env->time, U, m);
		P = findP(*env->time, U, reshape(measurements - *env->nth_soln, 1, n*lt).row(0), m);

		//matToGslMat(A + a*BT*B, qr);
		//vecToGslVec(P, b);

		if(i != brk){
			//gamma *= .5;
			//du = inverse(A.transpose()*A + gamma*gamma*B.transpose()*B)*A.transpose()*P;
			du = inverse(A + gamma*B.transpose()*B)*P;
			//cout << du.transpose() << endl;
		}else{
			gamma = inc+.0003;
			do{
				du = inverse(A + gamma*B.transpose()*B)*P;
				cout << gamma << "," << norm(uNot + du - u) << endl;
				gamma += inc;
			}while(gamma <= 50);
			exit(0);
		}
		uNot += du;
		if(du.norm() < TOL || std::isnan(du.norm())){
			break;
		}
	}
	log_err("Function stopped before break");
	exit(0);
}

//Used to find plots for the ||u_n + \delta u_n+1 - u*|| vs \aplha
//When the problem does not require regularization
void regparamexp1b(soln_env *env, vec u, int brk){
	int n = (*env->nth_soln).rows();
	int m = (*env->initial_params).size();
	int lt = (*env->time).size();
	const mat measurements = *env->nth_soln;

	vec uNot(m);
	uNot = *env->initial_params;
	mat bob;
	mat U(n*m, n*lt);
	mat A(m,m), AT(m,m);
	mat B = mat::Identity(m, m);
	for(int i=0; i<m-1; i++){
		B(i+1,i) = -1;
	}
	mat BT = B.transpose();
	vec P(m);
	vec du(m);
	double TOL = 0.00001;
	double gamma;
	int LIMIT = 200;
	cout.precision(10);

	for(int i = 0; i<LIMIT; i++)
	{
		bob = qLinearRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, measurements);
		U = reshape(bob.bottomRows(n*m), m, n*lt);
		*env->nth_soln = bob.topRows(n);
		A = findA(*env->time, U, m);
		P = findP(*env->time, U, reshape(measurements - *env->nth_soln, 1, n*lt).row(0), m);

		//matToGslMat(A + a*BT*B, qr);
		//vecToGslVec(P, b);

		if(i != brk){
			du = A.inverse()*P;
			cout << du.transpose() << endl;
		}else{
			gamma = 0.001;
			do{
				du = inverse(A + gamma*B.transpose()*B)*P;
				cout << gamma << "," << norm(uNot + du - u) << endl;
				gamma += 0.001;
			}while(gamma <= 5.0);
			exit(0);
		}
		uNot += du;
		if(du.norm() < TOL || std::isnan(du.norm())){
			break;
		}
	}
	log_err("Function stopped before break.");
	exit(0);
}

//Used to find plots for the ||u_n + \delta u_n+1 - u*|| vs \aplha
// When we try to increase the range of convergance using
// + a^2 ||u - u_0|| regularization
void reg_guess_plots(soln_env *env, vec u, vec u_guess, double gamma, int brk){
	int n = (*env->nth_soln).rows();
	int m = (*env->initial_params).size();
	int lt = (*env->time).size();
	const mat measurements = *env->nth_soln;

	vec uNot(m);
	uNot = *env->initial_params;
	mat bob;
	mat U(n*m, n*lt);
	mat A(m,m), AT(m,m);
	mat B = mat::Identity(m, m);
	vec P(m);
	vec du(m);
	double TOL = 0.00001, g;
	int LIMIT = 200;
	cout.precision(10);

	for(int i = 0; i<LIMIT; i++)
	{
		bob = qLinearRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, measurements);
		U = reshape(bob.bottomRows(n*m), m, n*lt);
		*env->nth_soln = bob.topRows(n);
		A = findA(*env->time, U, m);
		P = findP(*env->time, U, reshape(measurements - *env->nth_soln, 1, n*lt).row(0), m);
		if(i != brk){
			du = inverse(A + gamma*B)*(P + gamma*(u_guess - uNot));
		}else{
			g = 0.0000;
			do{
				du = inverse(A + gamma*B)*(P + gamma*(u_guess - uNot));
				//cout << gamma << "," << norm(uNot + du - u) << endl;
				//cout << norm(P-A*du) << "," << norm(uNot + du - u) << endl;
				cout << norm(P-A*du) << "," << norm(du) << endl;
				g += 0.000001;
			}while(g <= 2.0);
			exit(0);
		}
		uNot += du;
		if(du.norm() < TOL || std::isnan(du.norm())){
			break;
		}
	}
	log_err("Function stopped before break");
	exit(0);
}

// This function is to test if we can increase the range of convergance using
// + a^2 ||u - u_0|| regularization
vec reg_guess(soln_env *env, vec u_guess, double gamma){
    int n = (*env->nth_soln).rows();
    int m = (*env->initial_params).size();
    int lt = (*env->time).size();
    const mat measurements = *env->nth_soln;

    vec uNot(m);
    uNot = *env->initial_params;
    mat bob;
    mat U(n*m, n*lt);
    mat A(m,m), AT(m,m);
    mat B = mat::Identity(m, m);
    vec P(m);
    vec du(m);
    double TOL = 0.00001;
    int LIMIT = 1500;
    cout.precision(10);

	latexOutput(*env->nth_soln, uNot, 0, " &");
    for(int i = 0; i<LIMIT; i++)
    {
        bob = qLinearRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, measurements);
        U = reshape(bob.bottomRows(n*m), m, n*lt);
        *env->nth_soln = bob.topRows(n);
        A = findA(*env->time, U, m);
        P = findP(*env->time, U, reshape(measurements - *env->nth_soln, 1, n*lt).row(0), m);
		du = inverse(A + gamma*B)*(P + gamma*(u_guess - uNot));

		uNot += du;
        if(du.norm() < TOL || std::isnan(du.norm())){
			latexOutput(*env->nth_soln, uNot, i+1, " &");
			cout << gamma << " &" << endl;
            break;
        }
		latexOutput(*env->nth_soln, uNot, i+1, " &");
		cout << gamma << " &" << endl;
    }
	*env->nth_soln = measurements;
	return uNot;
}

// This function is to test if we can increase the range of convergance using
// + a^2 ||u - u_0|| regularization and the discrepancy principle to help with
// finding the regularization parameter
vec reg_guess2(soln_env *env, vec u_guess){
    int n = (*env->nth_soln).rows();
    int m = (*env->initial_params).size();
    int lt = (*env->time).size();
    const mat measurements = *env->nth_soln;

    vec uNot(m);
    uNot = *env->initial_params;
    mat bob;
    mat U(n*m, n*lt);
    mat A(m,m), AT(m,m);
    mat B = mat::Identity(m, m);
    for(int i=0; i<m-1; i++){
        //B(i+1,i) = -1;
    }
    mat BT = B.transpose();
    vec P(m);
    vec du(m);
    double TOL = 0.00001, g, O;
    int LIMIT = 1500;
    cout.precision(10);

	latexOutput(*env->nth_soln, uNot, 0, " &");
    for(int i = 0; i<LIMIT; i++)
    {
        bob = qLinearRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, measurements);
        U = reshape(bob.bottomRows(n*m), m, n*lt);
        *env->nth_soln = bob.topRows(n);
        A = findA(*env->time, U, m);
        P = findP(*env->time, U, reshape(measurements - *env->nth_soln, 1, n*lt).row(0), m);
		O = findO(reshape(measurements - *env->nth_soln, 1, n*lt).row(0), *env->time);
		g = findLambda(A, P, O, B, u_guess, uNot);
		uNot += inverse(A + g*B.transpose()*B)*(P + g*(u_guess - uNot));
        if(du.norm() < TOL || std::isnan(du.norm())){
			latexOutput(*env->nth_soln, uNot, i+1, " \\\\");
            break;
        }
		latexOutput(*env->nth_soln, uNot, i+1, " &");
    }
	return uNot;
}

double func(double a, mat A, vec P, double O, mat B, vec u_guess, vec uNot){
	vec du = inverse(A + a*B.transpose()*B)*(P + a*(u_guess - uNot));
	return std::sqrt(O - 2*P.dot(du) + du.dot(A*du)) - (1.21)*std::sqrt(O);
}

double findLambda(mat A, vec P, double O, mat B, vec u_guess, vec uNot){
	double a1 = 1.0001, a2 = 1.00011, a3;
	double f1 = func(a1, A, P, O, B, u_guess, uNot);
	double f2 = func(a2, A, P, O, B, u_guess, uNot);
	do{
		a3 = a2 - f2*(a2-a1)/(f2-f1);
		a1 = a2;
		a2 = a3;
	}while(abs(a2-a1) < .0000001);
	if(std::isnan(a3)){
		log_err("Regularization parameter is nan.");
	}
	log_err(a3);
	return a3;
}

// This function is to test if we can increase the range of convergance using
// + a^2 ||u|| regularization
vec reg1(soln_env *env, double gamma){
    int n = (*env->nth_soln).rows();
    int m = (*env->initial_params).size();
    int lt = (*env->time).size();
    const mat measurements = *env->nth_soln;

    vec uNot(m);
    uNot = *env->initial_params;
    mat bob;
    mat U(n*m, n*lt);
    mat A(m,m);
    mat B = mat::Identity(m, m);
    vec P(m);
    vec du(m);
    double TOL = 0.00001;
    int LIMIT = 1500;
    cout.precision(7);
	latexOutput(*env->nth_soln, uNot, 0, " &");
    for(int i = 0; i<LIMIT; i++)
    {
        bob = qLinearRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, measurements);
        U = reshape(bob.bottomRows(n*m), m, n*lt);
        *env->nth_soln = bob.topRows(n);
        A = findA(*env->time, U, m);
        P = findP(*env->time, U, reshape(measurements - *env->nth_soln, 1, n*lt).row(0), m);
		du = inverse(A + gamma*B)*(P);

		uNot += du;
        if(du.norm() < TOL || std::isnan(du.norm())){
			latexOutput(*env->nth_soln, uNot, i+1, " &");
			cout << gamma << " &" << endl;
            break;
        }
		latexOutput(*env->nth_soln, uNot, i+1, " &");
		cout << gamma << " &" << endl;
    }
	return uNot;
}

// This function is to test if we can increase the range of convergance using
// + a^2 ||u|| regularization
vec reg2(soln_env *env, vec u_guess){
    int n = (*env->nth_soln).rows();
    int m = (*env->initial_params).size();
    int lt = (*env->time).size();
    const mat measurements = *env->nth_soln;

    vec uNot(m);
    uNot = *env->initial_params;
    mat bob;
    mat U(n*m, n*lt);
    mat A(m,m), AT(m,m);
    mat B = mat::Identity(m, m);
    vec P(m);
    vec du(m);
    double TOL = 0.00001, g, O;
    int LIMIT = 1500;
    cout.precision(10);

	latexOutput(*env->nth_soln, uNot, 0, " &");
    for(int i = 0; i<LIMIT; i++)
    {
        bob = qLinearRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, measurements);
        U = reshape(bob.bottomRows(n*m), m, n*lt);
        *env->nth_soln = bob.topRows(n);
        A = findA(*env->time, U, m);
        P = findP(*env->time, U, reshape(measurements - *env->nth_soln, 1, n*lt).row(0), m);
		O = findO(reshape(measurements - *env->nth_soln, 1, n*lt).row(0), *env->time);
		g = findLambda(A, P, O, B, u_guess, uNot);
		uNot += inverse(A + g*B)*P;
        if(du.norm() < TOL || std::isnan(du.norm())){
			latexOutput(*env->nth_soln, uNot, i+1, " \\\\");
            break;
        }
		latexOutput(*env->nth_soln, uNot, i+1, " &");
    }
	return uNot;
}

struct curve_params { mat A; vec P; double O; };

double curvature(double a, void * params)
{
	struct curve_params * p = (struct curve_params *) params;
	mat A = (p->A);
	vec P = (p->P);
	double O = (p->O);
	int m = P.size();
	mat I = mat::Identity(m, m);
	vec x = inverse(A + a*I)*P;
	vec z = inverse(A + a*I)*(A*x - P);
	double nx = x.squaredNorm(),
			ny = abs(O - 2*P.dot(x) + x.dot(A*x)),
			nxprime = 4/sqrt(a)*x.dot(z);

	return -2*(nx/nxprime)*ny*((a*nxprime*ny + 2*sqrt(a)*nx*ny +
		a*a*nx*nxprime)/pow(a*nx*nx + ny*ny, 3/2));
}

double alpha(mat A, vec P, double O)
{
	struct curve_params params = {A,P,O};

	int status, iter = 0, max_iter = 150;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;

	gsl_function F;
	F.function = &curvature;
	F.params = &params;
	double m = 5.0e-7, a = 1.0e-13, b = 1;

	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc(T);
	gsl_min_fminimizer_set(s, &F, m, a, b);

	do
	{
		iter++;
		status = gsl_min_fminimizer_iterate(s);
		m = gsl_min_fminimizer_x_minimum(s);
		a = gsl_min_fminimizer_x_lower(s);
		b = gsl_min_fminimizer_x_upper(s);
		status = gsl_min_test_interval(a, b, 1.0e-10, 0.0);
		cout << "m = " << m << endl;
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	cout << "alpha values: "<< m << endl;
	gsl_min_fminimizer_free (s);
	return m;
}

void curvetest(mat A, vec P, double O)
{
	struct curve_params params = {A,P,O};
	double g = 1e-13;
	do{
		cout << g << "," << curvature(g, &params) << endl;
		g += 0.000001;
	}while(g <= 1.0e-1);
	exit(0);
}

//Used to find plots for the ||u_n + \delta u_n+1 - u*|| vs \aplha
// When we try to increase the range of convergance using
// + a^2 ||u - u_0|| regularization
void curvature_test_plots(soln_env *env, int brk){
	int n = (*env->nth_soln).rows();
	int m = (*env->initial_params).size();
	int lt = (*env->time).size();
	const mat measurements = *env->nth_soln;

	vec uNot(m);
	uNot = *env->initial_params;
	mat bob;
	mat U(n*m, n*lt);
	mat A(m,m), AT(m,m);
	mat B = mat::Identity(m, m);
	vec P(m);
	vec du(m);
	double TOL = 0.00001, gamma= 0.0005, O;
	int LIMIT = 200;
	cout.precision(10);

	for(int i = 0; i<LIMIT; i++)
	{
		bob = qLinearRungeKutta4(*env->ode, (*env->time), uNot, *env->initial_cond, measurements);
		U = reshape(bob.bottomRows(n*m), m, n*lt);
		*env->nth_soln = bob.topRows(n);
		A = findA(*env->time, U, m);
		P = findP(*env->time, U, reshape(measurements - *env->nth_soln, 1, n*lt).row(0), m);
		if(i != brk){
			du = inverse(A + gamma*B)*(P);
		}else{
			O = findO(reshape(measurements - *env->nth_soln, 1, n*lt).row(0), *env->time);
			curvetest(A, P, O);
		}
		uNot += du;
		if(du.norm() < TOL || std::isnan(du.norm())){
			break;
		}
	}
	log_err("Function stopped before break");
	exit(0);
}
