#include <iostream>
#include <stdlib.h>
#include <math.h>

#include <Eigen/LU>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <misc.h>
#include <nonlinear_odes.h>
#include <numerical_integration.h>
#include <regularization.h>
#include <input.h>
#include <bspline.h>

using namespace thesis;

int main(int argc, char *argv[])
{
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
	//End Console Input


	// Init key parameters
	int lt = t.size();
	int m  = uNot.size();
	int n  = yNot.size();

	// Create measurement and add noise if necessary
	mat measure;
	mat measure2;
	measure = rungekutta4(system, times, u, yNot);
	spline spl_msmtRow[n];
	size_t ncoeffs = 12;
	size_t order = 4;
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

	timelatexOutput(t, " &", measure.rows(), u.size());
	cout << "\\lambda &" << endl;

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
			spl_msmtRow[i].update(times, measure.row(i));
			for(int j = 0; j<lt; j++){
				msmt(i,j) = spl_msmtRow[i].interpolate(t(j));
			}
		}
	}
	lyNot.head(n) = msmt.col(0);
	//getNoise is being used to acquire a regularization parameter
	reg2(env, in.getUGuess());
	latexOutput(measure, u, 0, " \\\\");
	cout << " \\\\" << endl;
	return 0;
}
