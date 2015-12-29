#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <Eigen/LU>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <dbg.h>
#include <misc.h>
#include <nonlinear_odes.h>
#include <numerical_integration.h>
#include <thesis_functions.h>
#include <input.h>
#include <bspline.h>
#include <least_squares.h>
#include <latex_output.h>

#define NCOEFFS 10

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
	
	vec u;
	u = in.getU();
	
	vec yNot;
	yNot = in.getYNot();
		
	vec uNot;
	uNot = in.getUNot();
	
	vec t;
	t = in.getInterval();
	
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
	cout.precision(7);
	
	spline spl_msmtRow[n];	
	//size_t order = 4;
	//check(0 < lt-order-1, "number of coeffs is negative");
	//bspline msmtRows(order, NCOEFFS-order-1, lt);
	
	size_t order = 7;
	lsquares lsq_msmt(lt, order);
	mat msmt(n,lt);
	
	vec du(m);

	vec lyNot(n*(m+1));
	lyNot.fill(0); 
	
	soln_env* env;
	env = (soln_env*) malloc(sizeof(sys) + 4*sizeof(vec*) + 2*sizeof(mat*));
	env->ode = no.odeFuncMap[system];
	env->time = &t;
	env->initial_cond = &lyNot;
	env->initial_params = &uNot;
	env->actual_params = &u;
	env->nth_soln = &msmt;
	env->measurements = &measure;

	mat output(m,OUT_ARR_SIZE+3);
	output.col(0) = u;
	
	timelatexOutput(t, " &", measure.rows(), u.size());
	latexOutput(measure, u, 0, " &");
	
	for(int q=0; q<OUT_ARR_SIZE; q++)
	{
		if(in.isNoisy() == true){
			measure2 = noise(measure, in.getNoise(), q);
			if(in.useBSpline() == true){
				for(int i=0; i<n; i++){
					lsq_msmt.update(times, measure2.row(i));
					for(int j = 0; j<lt; j++){
						msmt(i,j) = lsq_msmt.interpolate(t(j));
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
		
		else{
			for(int i=0; i<n; i++){
				spl_msmtRow[i].update(times, measure.row(i));
				for(int j = 0; j<lt; j++){
					msmt(i,j) = spl_msmtRow[i].interpolate(t(j));
				}
			}
		}
		lyNot.head(n) = msmt.col(0);
		du = findActualParam(env, in.isRegularized());
		if(isnan(du.norm())){
			q -= 1;
			latexOutput(msmt, du, q+1, " ....bad run");
		}else{
			output.col(q+1) = du;
			if(q != OUT_ARR_SIZE-1){
				latexOutput(msmt, du, q+1, " &");
			}else{
				latexOutput(msmt, du, q+1, " \\\\");	
			}
		}
	}
	
	longlatexOutput(output);	
	shortlatexOutput(output);
	shortNormalizedLatexOutput(output);
	
	for(int i=0; i<n; i++){
		spl_msmtRow[i].update(times, measure.row(i));
	
		for(int j = 0; j<lt; j++){
			msmt(i,j) = spl_msmtRow[i].interpolate(t(j));
		}
	}
	
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
	
	for(int i=0; i<oM.rows(); i++){
		for(int j=0; j<oM.cols(); j++){
			oM(i,j) += noise*gsl_ran_ugaussian(r);
		}
	}
	gsl_rng_free (r);
	
	return oM;
}