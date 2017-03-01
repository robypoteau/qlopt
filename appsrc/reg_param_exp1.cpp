#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <Eigen/LU>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <misc.h>
#include <nonlinear_odes.h>
#include <numerical_integration.h>
#include <thesis_functions.h>
#include <input.h>
#include <bspline.h>
#include <latex_output.h>
#include <mp_spline.h>

using namespace thesis;

bool isNegative(const vec& x);
bool isLessThanOne(const vec& x);
mat noise(const mat& M, double noise, int seed);

int main(int argc, char *argv[])
{
	const int bits = 128;	//const int bytes = bits/8;
	mpreal::set_default_prec(bits);

	time_t begin;
	begin = time(NULL);

	nonlinearOdes no;
	input in(argc, argv);

	//Gather input from the console
	string system = in.getSystem();

	vec times;
	times = in.getTimeData();

	vec u;
	u = in.getU();

	vec yNot;
	yNot = in.getYNot();

	vec uNot;
	uNot = in.getUNot();

	vec t;
	t = in.getInterval();
	mp_mat mp_t = t.cast<mpreal>();
	int OUT_ARR_SIZE = in.getNumberOfIterations();
	//End Console Input


	// Init key parameters
	int lt = t.size();
	int m  = uNot.size();
	int n  = yNot.size();

	// Create measurement and add noise if necessary
	mat measure;
	mat measure2;
	measure = rungekutta4(no.odeFuncMap[system], times, u, yNot);

	mp_mat mp_measure;
	mp_mat mp_measure2;
	//mp_measure = mp_rungekutta4(system, times.cast<mpreal>(), u.cast<mpreal>(),\
	//	yNot.cast<mpreal>());

	spline spl_msmtRow[n];
	mp_spline mp_spl_msmtRow[n];
	size_t ncoeffs = 12;
	size_t order = 4;
	//check(0 < lt-order-1, "number of coeffs is negative");
	bspline msmtRows(order, ncoeffs, lt);

	//size_t order = in.getNcoeffs();
    //lsquares lsq_msmt(times.size(), order);
    //tsqr lsq_msmt(times.size(), order);
    //nnls lsq_msmt(times.size(), order);
    //logittsqr lsq_msmt(times.size(), order);
	//expo_tsqr lsq_msmt(times);

	mat msmt(n,lt);
	mp_mat mp_msmt(n,times.size());

	vec du(m);

	vec lyNot(n*(m+1));
	lyNot.fill(0);

	soln_env* env;
	env = (soln_env*) malloc(sizeof(string*) + 4*sizeof(vec*) + 2*sizeof(mat*) \
		+ 2*sizeof(mp_mat*));
	env->ode = &system;
	env->time = &t;
	env->initial_cond = &lyNot;
	env->initial_params = &uNot;
	env->nth_soln = &msmt;
	env->measurements = &measure;
	env->mp_nth_soln = &mp_msmt;
	env->mp_measurements = &mp_measure;

	cout.precision(7);
	mat output(m,OUT_ARR_SIZE+3);
	mat output2(n*(OUT_ARR_SIZE+1),lt);
	output.col(0) = u;
	output2.topRows(n) = measure;

	timelatexOutput(t, " &", measure.rows(), u.size());
	latexOutput(measure, u, 0, " &");
	int badrun = 0;

	for(int q=0; q<OUT_ARR_SIZE; q++)
	{
		if(in.isNoisy() == true){
			measure2 = noise(measure, in.getNoise(), q+badrun);
			if(in.useBSpline() == true){
				for(int i=0; i<n; i++){
					msmtRows.update(times, measure2.row(i));
					//lsq_msmt.update(times, measure2.row(i));
					for(int j = 0; j<lt; j++){
						msmt(i,j) = msmtRows.interpolate(t(j));
						//lsq_msmt.interpolate(t(j));
					}
				}
			}else{
				for(int i=0; i<n; i++){
					spl_msmtRow[i].update(times, measure2.row(i));
					for(int j = 0; j<lt; j++){
						msmt(i,j) = spl_msmtRow[i].interpolate(t(j));
					}
				}
			}
		}
		latexOutput(msmt, u, -11111, " &");
		//cout << endl << msmt << endl;
        log_info(norm(measure.leftCols(msmt.cols())-msmt));
		lyNot.head(n) = msmt.col(0);

		du = findActualParam(env, in.isRegularized(), in.getNumDivs());

		if(std::isnan(du.norm())){
			q -= 1;
			//latexOutput(msmt, du, q+1, " ....bad run");
			badrun += 1;
			cout << "badrun "<< badrun << endl;
		}else{
			output.col(q+1) = du;
			output2.middleRows((q+1)*n,n) = msmt;
			if(q != OUT_ARR_SIZE-1){
			cout <<"goodrun "<< q << ", ";// << endl;
				//latexOutput(msmt, du, q+1, " &");
			}else{
				cout <<"goodrun "<< q << endl;
				latexOutput(msmt, du, q+1, " \\\\");
			}
		}
	}
	//longlatexOutput(output);
	cout << "(bad runs: " << badrun << ")" <<endl;
	shortlatexOutput(output);
	//shortNormalizedLatexOutput(output);
	R(t(1)-t(0), output2, n);
	//M(output2, n); Mi(output2, n);
	time_t end;
	end = time(NULL);
	cout << end - begin << endl;
	return 0;
}

bool isNegative(const vec& x){
	bool value = false;
	for(int i=0; i < x.size(); i++){
		if(x(i) < 0){
			value = true;
			break;
		}
	}
	return value;
}

bool isLessThanOne(const vec& x){
	bool value = false;
	for(int i=0; i < x.size()-1; i++){
		if(x(i) < .9){
			value = true;
			break;
		}
	}
	if(x(x.size()-1) < 1)
			value = true;

	return value;
}

mat noise(const mat& M, double noise, int seed){
	mat oM = M;

	const gsl_rng_type* T;
	gsl_rng* r;
	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	gsl_rng_set(r, time(NULL)+seed);
	double temp;

	for(int j=0; j<oM.cols(); j++){
		temp = noise*gsl_ran_ugaussian(r);
		for(int i=0; i<oM.rows(); i++){
			oM(i,j) += temp;
		}
	}
	gsl_rng_free (r);

	return oM;
}
