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
//#include <least_squares.h>
//#include <tsqr.h>
//#include <nnls.h>
//#include <logittsqr.h>
//#include <expo_tsqr.h>
#include <latex_output.h>

using namespace thesis;

bool isNegative(const vec& x);
bool isLessThanOne(const vec& x);
mat noise(const mat& M, double noise, int seed);

int main(int argc, char *argv[])
{
	time_t begin;
	begin = time(NULL);

	nonlinearOdes no;
	input in(argc, argv);

	//Gather input from the console
	string system = in.getSystem();

	vec times;
	times = in.getTimeData();

	//vec times(21);
	//times << 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20;

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

	//mat measure(2,16);
	//measure << 35,65,78,82,65,26,15,10,1,2,3,22,75,95,78,20,
	//			0.287,0.38,0.762,1.467,1.904,3.029,3.178,1.806,0.664,0.349,0.288,0.575,0.931,1.265,1.692,1.861;

	mat measure;
	//measure << 28,20,15,15,25,35,65,78,82,65,26,15,10,1,2,3,22,75,95,78,20,
	//			4.242,4.664,1.889,0.722,0.317,0.287,0.38,0.762,1.467,1.904,3.029,3.178,1.806,0.664,0.349,0.288,0.575,0.931,1.265,1.692,1.861;

	//measure << 30,47.2,70.2,77.4,36.3,20.6,18.1,21.4,22,25.4,27.1,40.3,57,76.6,52.3,19.5,11.2,7.6,14.6,16.2,24.7,
	//			4,6.1,9.8,35.2,59.4,41.7,19,13,8.3,9.1,7.4,8,12.3,19.5,45.7,51.1,29.7,15.8,9.7,10.1,8.6;

	//mat measure;
	mat measure2;
	measure = rungekutta4(no.odeFuncMap[system], times, u, yNot);

	spline spl_msmtRow[n];
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

	vec du(m);

	vec lyNot(n*(m+1));
	lyNot.fill(0);

	soln_env* env;
	env = (soln_env*) malloc(sizeof(string*) + 4*sizeof(vec*) + 2*sizeof(mat*));
	env->ode = &system;
	env->time = &t;
	env->initial_cond = &lyNot;
	env->initial_params = &uNot;
	env->nth_soln = &msmt;
	env->measurements = &measure;

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
				spl_msmtRow[i].update(times, measure.row(i));
				for(int j = 0; j<lt; j++){
					msmt(i,j) = spl_msmtRow[i].interpolate(t(j));
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
