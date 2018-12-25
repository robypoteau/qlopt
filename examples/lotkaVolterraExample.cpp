#include <qlopt.h>

//The namespace for the objects found in qlopt.h
using namespace thesis;

// Construct the ODE system for your examples it should have the same structure/arguments
vec lotka (const double& t, const vec&  x, const vec& u, const vec& input)
{
	(void) t;
	(void) input;
	mat result(2,1);
	result(0) = u(0)*x(0) - u(1)*x(0)*x(1);
	result(1) = u(1)*x(0)*x(1) - u(2)*x(1);

	return result;
}

mat rungekutta4(odefunction fhandle, const vec& t, const vec& u,
                const vec& y0, const vec& input);

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
	//params.tol.absparam = 1E-7; 	//Change optional, default value given
	//params.tol.relparam = 1E-7; 	//Change optional, default value given
	//params.tol.absobj = 1E-7; 	//Change optional, default value given
	//params.tol.relobj = 1E-7; 	//Change optional, default value given
	//params.tol.maxiter = 150; 	//Change optional, default value given

	//Data parameters. Time information used in qlopt.
	params.dat.initialTime = 0.0; 	//Should be set to proper value
	params.dat.endTime = 1.0;			//Should be set to proper value
	params.dat.timeIncrement = 0.01;	//Should be set to proper value
	params.dat.numOfDataSets = 1;

	//Regularization parameters.
	params.reg.type = 4; 	// 0 - none,
							// 1 - Type 1 Tikhonov using ||delta u_{N} - 0||
							// 2 - Type 2 Tikhonov using ||u_{N+1} - u_{N}||
							// 3 - Brute force search for alpha
                            // 4 -
	params.reg.alpha = 0.000001;

	//General parameters.
	params.gen.numOfStates = 2;		//Should be set to proper value
	params.gen.numOfParams = 3;		//Should be set to proper value
	params.gen.divisions = 1;		//Change optional, default value given
	params.gen.finitediff = true;	//Finite difference only

	//Inputs and data
	vec t;
	std::vector<vec> input(params.dat.numOfDataSets, vec::Zero(1));
	std::vector<mat> data(params.dat.numOfDataSets);

	//Simulated data
	vec u(3);
  	vec u0(3);
  	vec uguess(3);
  	vec y0(2);

  	u << 1,2,1;
	u0 << 3,4,3;
	uguess << 0,0,0; // The guess value for Type II Regularization
	y0 << 35, 4;

	t = vec::LinSpaced(101,0.0,1.0);
	for(size_t i = 0; i<params.dat.numOfDataSets; i++)
    {
		data[i] = rungekutta4(lotka, t, u, y0, input[i]);
    }

	/*
		This structure contains the many outputs of the method. Which the user
		can manipulate to get the desired graphics and numerical summaries.
		It contains u values at each iteration.
	*/
	outputStruct results;

	results = qlopt(lotka, t, u0, uguess, y0, input, data, params);

	results.uvals.col(results.uvals.cols()-1) = u;

	//Using the results from qlopt to construct a latex table
	std::cout << "\nu =" << results.ufinal.transpose() << endl;
	cout << "\nLatex Output:" << endl << endl;
	parameterOutput(results.uvals);

	return 0;
}

mat rungekutta4(odefunction fhandle, const vec& t, const vec& u,
                const vec& y0, const vec& input)
{
  const unsigned int n = y0.size();
  const unsigned int N = t.size();
  double h = 0.0;

  mat w(n, N);
  w.fill(0.0);
	w.col(0) = y0;

	vec k1(n), k2(n), k3(n), k4(n);

	for (size_t i = 0; i<N-1; i++)
    {
      h = t(i+1)-t(i);

      k1 = h*fhandle(t(i), 		 w.col(i), 		  u, input);
      k2 = h*fhandle(t(i) + h/2, w.col(i) + k1/2, u, input);
      k3 = h*fhandle(t(i) + h/2, w.col(i) + k2/2, u, input);
      k4 = h*fhandle(t(i) + h, 	 w.col(i) + k3,   u, input);

      w.col(i+1) = w.col(i) + (k1 + 2*(k2 + k3) + k4)/6;
    }

	return w;
}
