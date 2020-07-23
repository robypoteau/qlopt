#ifndef QLOPT_H
#define QLOPT_H
// Built-in headers
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <chrono>

// External headers
#include <eigen3/Eigen/QR>
#include <eigen3/Eigen/SVD>

// My headers
#include <misc.h>
#include <spline.h>
#include <latex_output.h>
#include <odesolver.h>
#include <odeWrapper.h>
#include <regularization.h>
//#include <splinterSpline.h>

using namespace std::chrono;

namespace thesis{
	//TODO Make constructor to input values
	typedef struct input_struct {
		input_struct(){
			tolerance();
			regularization();
			data();
			initialValueProblem();
			general();
			smoothing();
		}
		struct tolerance {
			tolerance(){
				absobj = 1E-7;
				relobj = 1E-7;
				absparam = 1E-7;
				relparam = 1E-7;
				normdiff = 1E-7;
				objval = 1E-7;
				maxiter = 150;
			};
			double absobj;
			double absparam;
			double relobj;
			double relparam;
			double normdiff;
			double objval;
			unsigned int maxiter;
		} tol;
		struct data {
			data(){
				//spacing = "uniform";
				initialTime = 0.0;
				endTime = 0.0;
				timeIncrement = 0.0;
				numOfDataSets = 1;
			};
			//std::string spacing;
			double initialTime;
			double endTime;
			double timeIncrement;
			unsigned int numOfDataSets;
		} dat;
		struct regularization {
			regularization() : type(0), alpha(0.0)
				{};
			unsigned int type;
			double alpha;
		} reg;
		struct initialValueProblem {
			initialValueProblem(){
				solver="rk4";//"boost_rk4";//"cvodes";
			};
			std::string solver;
		} ivp ;
		struct general {
			general(){
				numOfStates = 1;
				numOfParams = 1;
				divisions = 1;
				finitediff = true;
			};
			unsigned int numOfStates;
			unsigned int numOfParams;
			unsigned int divisions;
			bool finitediff;
		} gen;
		struct smoothing {
			smoothing() {
				regParam = 0.03;
			}
			double regParam;
		} noise;
	} inputStruct;

	typedef struct output_struct{
		unsigned int numfevals;
		unsigned int iterations;
		double comptime;
		vec ufinal;
		mat xvals;
		mat uvals;
		vec omegaval;
        vec objval;
        vec alpha;
        vec deltau;


		output_struct(){
			numfevals = 0;
			iterations = 0;
			comptime = 0.0;
		}
		output_struct(const thesis::output_struct& op){
			numfevals = op.numfevals;
			iterations = op.iterations;
			comptime = op.comptime;
		}
	} outputStruct;

	outputStruct qlopt(odefunction fun,
		const vec& t,
		const vec& u0,
		const vec& uguess,
		const vec& y0,
		const vector<vec>& input,
		const vector<mat>& data,
		const inputStruct& params,
		const vec& u);

	mat findA(const vec& t, const mat& U, const size_t& m);
	vec findP(const vec& t, const mat& U, const vec& dx, const size_t& m);
	double findO(const vec& t, const vec& dx);
	double innerProd(const vec& u1, const vec& u2, const vec& time);
	double simpson(const vec& t, const vec& x);
	mat reshape(const mat& U, int n, int m);
	double norm(const mat& M);
	mat inverse(const mat& M);
	mat corrMat(const mat& M);
	double rcond(const mat& M);
	double matnorm1(const mat& M);
}
#endif
