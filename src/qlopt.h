#ifndef QLOPT_H
#define QLOPT_H
// Built-in headers
#include <vector>
#include <iostream>
#include <string>
#include <cmath>

// External headers
 #include <Eigen/QR>

// My headers
#include <misc.h>
#include <spline.h>
#include <tallskinnyqr.h>
#include <latex_output.h>

namespace thesis{
	typedef vec (*odefunction)(const double& t, 
		const vec& x, const vec& u, const vec& control);

	//TODO Make constructor to input values
	typedef struct input_struct {
		input_struct(){
			tolerance();
			regularization();
			data();
			initialValueProblem();
			general();			
		}
		struct tolerance {
			tolerance(){
				absparam =  1E-7;
				relparam =  1E-7;
				absobj =  1E-7;
				relobj =  1E-7;
				maxiter = 500;
				//unsigned int maxfunceval = 3E2;
			};
			double absparam;
			double relparam;
			double absobj;
			double relobj;
			unsigned int maxiter;
			//double maxfunceval;
		} tol ;
		struct data {
			data(){
				spacing = "uniform";
				initialTime = 0.0;
				endTime = 0.0;
				timeIncrement = 0.0;
				numOfDataSets = 1;
			};
			std::string spacing;
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
				solver="boost_rk4";//"cvodes";			
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
	} inputStruct;

	typedef struct output_struct{
		unsigned int numfevals;
		unsigned int iterations;
		double comptime;
		vec ufinal;
		mat xvals;
		mat uvals;
				
		
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

	class OdeWrapper
	{
	private:
		vec control;
		odefunction fun; 
	public:
		OdeWrapper(odefunction of){fun = of;}
		void setControl(vec input){control = input;}
		vec operator() (double t, vec  x, vec u)
		{
			//TODO check that the control has been set.
			return fun(t, x, u, control);
		}
	};

	outputStruct qlopt(odefunction fun, 
		const vec& t, 
		const vec& u0,
		const vec& uguess, 
		const vec& y0,
		const vector<vec>& input, 
		const vector<mat>& data, 
		const inputStruct& params);
	
	double findGamma(mat A, vec P, vec uNot, vec u);
		
	mat qloptRungeKutta4(OdeWrapper& fhandle, const vec& time, 
		const vec& u, const vec& yNot,  std::vector<thesis::spline>& Xn);
		
	mat qlinear(OdeWrapper& fhandle, const double& t, const vec& x, 
		const vec& u, std::vector<thesis::spline>& Xn);
			
	mat der(const mat& dx, const double& dt);
	
	mat jac(OdeWrapper& f, double t, const mat& x, const mat& u, 
		const double& h);
	
	mat findA(const vec& t, const mat& U, const size_t& m);
	vec findP(const vec& t, const mat& U, const vec& dx, const size_t& m);
	double findO(const vec& t, const vec& dx);
	double findGamma(double initialGuess, void * params);
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
