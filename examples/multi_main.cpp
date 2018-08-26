#include <qlopt.h>

//The namespace for the objects found in qlopt.h
using namespace thesis;

// Construct the ODE system it should have the same structure/arguments
vec lotka (const double& t, const vec&  x, const vec& u, const vec& input)
{
	(void) t;
	mat result(2,1);
	result(0) = u(0)*x(0) - u(1)*x(0)*x(1) + input(0);
	result(1) = u(1)*x(0)*x(1) - u(2)*x(1) + input(1);

	return result;
}

int main(int argc, char *argv[])
{
	(void) argc;
	(void) argv;
	/*
		Construct the input parameters for qlopt method. The struct has default
		values, but some values still need to be adjusted. All the elements of 
		the struct will be displayed and some changed.
	*/
	inputStruct params;
	//Tolerance parameter.
	//params.tol.absparam = 1E-7; //Change optional, default value given
	//params.tol.relparam = 1E-7; //Change optional, default value given
	//params.tol.absobj = 1E-7; 	//Change optional, default value given
	//params.tol.relobj = 1E-7; 	//Change optional, default value given
	//params.tol.maxiter = 500; 	//Change optional, default value given
	
	//Data parameters.
	params.dat.spacing = "uniform";	//Options: "uniform", "nonuniform"
	params.dat.initialTime = 0.000; 	//Should be set to proper value
	params.dat.endTime = 0.011;		//Should be set to proper value
	params.dat.timeIncrement = 0.001;	//Should be set to proper value
	params.dat.numOfDataSets = 2;
	
	//Regularization parameters.
	params.reg.type = 1; 	// 0 - none, 
							// 1 - Type 1, 
							// 2 - Type 2
	params.reg.alpha = 0.0;						
	
	//Initial value problem parameters.
	//params.ivp.solver = "rk4";	// "rk4", "boost_rk4", "cvodes", ...
	
	//General parameters.
	params.gen.numOfStates = 2;		//Should be set to proper value
	params.gen.numOfParams = 3;		//Should be set to proper value
	params.gen.divisions = 1;		//Change optional, default value given
	params.gen.finitediff = true;	//Change optional, default value given
	
	mat d1(2,12), d2(2,12);
	d1 << 35,35.013,35.026,35.04,35.053,35.066,35.079,35.092,35.105,35.119,35.132,35.145,
		  4,3.9999,3.9998,3.9998,3.9997,3.9996,3.9995,3.9995,3.9994,3.9993,3.9993,3.9992;
		  
	//d2 << 88.06,68.51,32.19,12.64,21.49,30.35,2.18,152.65,148.36,85.81,41.41,
	//	  34.38,29.59,21.30,13.69,7.650,4.080,4.09,14.330,38.220,60.78,70.77;
		  
	//Inputs and data
	std::vector<vec> input(2);
	std::vector<mat> data(2);
	vec t(12);
	vec u0(3);
  vec u(3);
  vec y0(2);

	input[0] = vec::Zero(2);
	input[1] = vec::Zero(2); //input[1](0) = 68.48; input[1](1) = 4.29;

	t << 0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011;
	data[0] = d1;
	data[1] = d1;
	u0 << .5,.5,.5;
  u << 1,2,1;
	y0 << 35, 4;

	/*
		This structure contains the many outputs of the method. Which the user
		can manipulate to get the desired graphics and numerical summaries.
	*/
	qlopt(lotka, t, u0, u, y0, input, data, params);
  //qlopt(benchmark, t2, u0, u, y0, input, data, params);
	return 0;
}



