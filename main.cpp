#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <Eigen/LU>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "dbg.h"
#include "misc.h"
#include "nonlinear_odes.h"
#include "numerical_integration.h"
#include "thesis_functions.h"
#include "input.h"
#include "bspline.h"
#include "latex_output.h"

#define OUT_ARR_SIZE 5

using namespace thesis;

double cond(const mat& A);
bool isNegative(const vec& x);
bool isLessThanOne(const vec& x);
mat noise(const mat& M, double noise);

int main(int argc, char *argv[])
{
	time_t begin;
	begin = time(NULL);
	
	nonlinearOdes no;
	input in(argc, argv);
	
	//Gather input from the console
	bool reg = in.isRegularized();
	bool noisy = in.isNoisy();
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
	
	spline msmtRow[n];	
	size_t order = 4;
	check(0 < lt-order-1, "number of coeffs is negative")
	bspline msmtRows(order, lt-order-1, lt);
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
	
	for(int q=0; q<OUT_ARR_SIZE; q++)
	{
		if(noisy == true){
			measure2 = noise(measure, in.getNoise());
			for(int i=0; i<n; i++){
				msmtRows.update(times, measure2.row(i));
				for(int j = 0; j<lt; j++){
					msmt(i,j) = msmtRows.interpolate(t(j));
				}
			}
		}
		
		else{
			for(int i=0; i<n; i++){
				msmtRow[i].update(times, measure.row(i));
				for(int j = 0; j<lt; j++){
					msmt(i,j) = msmtRow[i].interpolate(t(j));
				}
			}
		}
		lyNot.head(n) = msmt.col(0);
		du = findActualParam(env, reg);
		output.col(q+1) = du;
	}
	
	longlatexOutput(output);	
	shortlatexOutput(output);
	
	for(int i=0; i<n; i++){
		msmtRow[i].update(times, measure.row(i));
	
		for(int j = 0; j<lt; j++){
			msmt(i,j) = msmtRow[i].interpolate(t(j));
		}
	}
	time_t end;
	end = time(NULL);
	cout << end - begin << endl;
	return 0;
}

double cond(const mat& A){
	return A.norm()*A.inverse().norm();
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

mat noise(const mat& M, double noise){
	mat oM = M;

	const gsl_rng_type* T;
	gsl_rng* r;
	gsl_rng_env_setup();
	
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	
	gsl_rng_set (r, time(NULL));
	
	for(int i=0; i<oM.rows(); i++){
		for(int j=0; j<oM.cols(); j++){
			oM(i,j) += noise*gsl_ran_ugaussian(r);
		}
	}
	gsl_rng_free (r);
	
	return oM;
}