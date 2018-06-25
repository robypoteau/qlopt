#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/LU>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <misc.h>
#include <nonlinear_odes.h>
#include <numerical_integration.h>
#include <thesis_functions.h>
#include <input.h>
#include <bspline.h>
#include <spline.h>
//#include <least_squares.h>
#include <tsqr.h>
//#include <nnls.h>
//#include <logittsqr.h>
//#include <expo_tsqr.h>
#include <latex_output.h>

using namespace thesis;

mat noise(const mat& M, double noise, gsl_rng* r);
std::vector<std::string> split(const std::string &s, char delim);
mat convertData(string filname);

class lotka
{
private:
	vec input;
public:
	lotka(vec control){input = control;}
	mat operator() (double t, vec  x, vec u)
	{
		mat result(2,1);
		result(0) = u(0)*x(0) - u(1)*x(0)*x(1) + input(0);
		result(1) = u(1)*x(0)*x(1) - u(2)*x(1) + input(1);
		
		return result;
	}
};

int main(int argc, char *argv[])
{
	nonlinearOdes no;
	//Gather input from the console
	input in(argc, argv);
	string system = in.getSystem();

	vec u;
	u = in.getU();

	vec yNot;
	yNot = in.getYNot();

	vec uNot;
	uNot = in.getUNot();
	
	vec uGuess;
	if(in.getRegs()==2){
		uGuess = in.getUGuess();
	}else{
		uGuess = uNot;
	}

	vec t;
	t = in.getInterval();
	
	double lambda;
	if(in.getRegs() > 0){
		lambda = in.getLambda();
	}else{
		lambda = 0.0;
	}
	
	size_t OUT_ARR_SIZE = in.getNumberOfIterations();
	//End Console Input


	/*
		This block of code is to get data from a 
		CSV file and convert it into a time vector
		and data matrix.
	*/
	
	// Init key parameter
	size_t n;
	size_t lt;
	vec timeData;
	mat featureData;
	
	if(in.isSimulatedData())
	{
		timeData = in.getTimeData();
		featureData = rungekutta4(system, timeData, u, yNot);
		n  = yNot.size();
	}
	else
	{
		string filename = "data/" + system + "_data.csv"; 
	
		std::ifstream csvfile;
		std::string line;
		std::vector<std::string> dataAsStr, headers, dataRow;
		csvfile.open(filename);
		if(csvfile.is_open())
		{
			while(std::getline(csvfile, line))
			{	
				dataAsStr.push_back(line);
			}
			csvfile.close();
		}
		else
		{
			std::cerr << "ERROR: Unable to open file"; 
		}
	
		headers = split(dataAsStr[0], ',');
		dataAsStr.erase(dataAsStr.begin());
	
		n = headers.size() - 1;
		lt = dataAsStr.size();
		timeData.setZero(lt);
		featureData.setZero(n, lt);
	
		for(size_t i=0; i<lt; i++)
		{
			dataRow = split(dataAsStr[i], ',');
			timeData(i) = std::stod(dataRow[0]);
			for(size_t j = 0; j<n; j++)
			{
				featureData(j,i) = std::stod(dataRow[j+1]);
			}
		}
	}
	
	lt = t.size();
	size_t m  = uNot.size();
	//int n  = yNot.size();
	
	// Create measurement and add noise if necessary
	spline fdsa;
	mat measure(n,lt);
	for(size_t i=0; i<n; i++)
	{
		fdsa.update(timeData, featureData.row(i));
		for(size_t j=0; j<lt; j++)
		{	
			measure(i,j) = fdsa.interpolate(t(j));
		}
	}

	mat measure2;	
	//spline spl_msmtRow[n];
	//size_t ncoeffs = 12;
	//size_t order = 4;
	//check(0 < lt-order-1, "number of coeffs is negative");
	//bspline msmtRows(order, ncoeffs, lt);

	size_t lsorder = in.getNcoeffs();
    //lsquares lsq_msmt(times.size(), order);
    tsqr lsq_msmt(timeData.size(), lsorder);
    //nnls lsq_msmt(times.size(), order);
    //logittsqr lsq_msmt(times.size(), order);
	//expo_tsqr lsq_msmt(times);

	mat msmt(n,lt);

	vec du(m);

	vec lyNot(n*(m+1));
	lyNot.fill(0);

	soln_env* env;
	env = (soln_env*) malloc(sizeof(string*) + 4*sizeof(vec*) + 2*sizeof(mat*) + sizeof(double));
	env->ode = &system;
	env->time = &t;
	env->initial_cond = &lyNot;
	env->initial_params = &uNot;
	env->u_guess = &uGuess;
	env->nth_soln = &msmt;
	env->measurements = &measure;
	env->lambda = &lambda;

	cout.precision(7);
	mat duoutput(m,OUT_ARR_SIZE);
	mat output2(n*(OUT_ARR_SIZE+1),lt);
	vec ctimeOutput(OUT_ARR_SIZE);
	vec iterOutput(OUT_ARR_SIZE);
	vec fevalOutput(OUT_ARR_SIZE);
	output2.topRows(n) = measure;

	//timelatexOutput(t, " &", measure.rows(), u.size());
	//latexOutput(measure, u, 0, " &");
	int badrun = 0;
	std::chrono::duration<double, milli> elapsed;
	output results;
	
	//Random number generator
	const gsl_rng_type* T;
	gsl_rng* r;
	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	//end rand

	gsl_rng_set(r, time(NULL));
	
	
	for(size_t q=0; q<OUT_ARR_SIZE; q++)
	{
		if(in.isNoisy() == true)
		{
			measure2 = noise(featureData, in.getNoise(), r);
			for(size_t i=0; i<n; i++)
			{
				//msmtRows.update(timeData, measure2.row(i));
				lsq_msmt.update(timeData, measure2.row(i));
				for(size_t j = 0; j<lt; j++)
				{
					//msmt(i,j) = msmtRows.interpolate(t(j));
					msmt(i,j) = lsq_msmt.interpolate(t(j));
				}
			}
			log_info(norm(measure.leftCols(msmt.cols())-msmt));
		}
		else
		{
			msmt = measure;
		}
		
		//latexOutput(msmt, u, -11111, " &");
		lyNot.head(n) = msmt.col(0);
		
    	std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
		results = findActualParam(env, in.getRegs(), in.getNumDivs());
    	std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
		elapsed = end-start;
		du = results.du;
		
		if(std::isnan(du.norm())){
			q -= 1;
			//latexOutput(msmt, du, q+1, " ....bad run");
			badrun += 1;
			cout << "badrun "<< badrun << endl;
		}else{
			//Output Results
			ctimeOutput(q) = (elapsed.count()/1000.0);
			fevalOutput(q) = (results.fevals);
			iterOutput(q) = (results.iterations);
			duoutput.col(q) = du;
			output2.middleRows((q+1)*n,n) = msmt;
			
			
			if(q != OUT_ARR_SIZE-1){
			cout <<"goodrun "<< q << ", ";// << endl;
				//latexOutput(msmt, du, q+1, " &");
			}else{
				cout <<"goodrun "<< q << endl;
				//latexOutput(msmt, du, q+1, " \\\\");
			}
		}
	}
	gsl_rng_free (r);
	//longlatexOutput(output);
	//cout << "(bad runs: " << badrun << ")" <<endl;
	//shortlatexOutput(duoutput);
	//shortNormalizedLatexOutput(output);
	//R(t(1)-t(0), output2, n);
	//M(output2, n); Mi(output2, n);
	//cout.precision(14);
	vec dumid;
	cout << "\nSolutions for " << OUT_ARR_SIZE << " run(s):\n" << duoutput << endl << endl; 
	cout << "True Solution: \n" << u << endl << endl;
	for(size_t i =0; i<m; i++){
		for(size_t j =0; j<OUT_ARR_SIZE; j++){
			duoutput(i,j) = abs(duoutput(i,j) - u(i))/u(i);
		}
	}
	dumid = duoutput.colwise().sum()/m;
	cout << "accuracy of parameter (relative mean error): " << dumid.sum()/dumid.size() << "\n";
	cout << "accuracy of parameter (relative max  error): " << dumid.sum()/dumid.size() << "\n";
	cout << "mean iterations: "  <<  iterOutput.sum()/iterOutput.size()  << "\n";	
	cout << "mean computational time (local search): " << ctimeOutput.sum()/ctimeOutput.size() << "s\n";
	
	return 0;
}

/*
	Noise model: x + noise*x*e
	x - true signal
	noise - a percentage
	e ~ N(0,1)
*/
mat noise(const mat& M, double noise, gsl_rng* r){
	mat oM = noise/100*M;
	
	for(int j=0; j<oM.cols(); j++){
		for(int i=0; i<oM.rows(); i++){
			oM(i,j) = M(i,j) + oM(i,j)*gsl_ran_ugaussian(r);;
		}
	}

	return oM;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> tokens;
    while (std::getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

mat convertData(string filename)
{
	/*
		This block of code is to get data from a 
		CSV file and convert it into a time vector
		and data matrix.
		
		TODO Convert this into a function for readability
	*/
	std::ifstream csvfile;
	std::string line;
	std::vector<std::string> dataAsStr, headers, dataRow;
	csvfile.open(filename);
	if(csvfile.is_open())
	{
		while(std::getline(csvfile, line))
		{	
			dataAsStr.push_back(line);
		}
    	csvfile.close();
    }
	else
	{
		std::cerr << "ERROR: Unable to open file"; 
	}
	
	headers = split(dataAsStr[0], ',');
	dataAsStr.erase(dataAsStr.begin());
	
	size_t n = headers.size() - 1;
	size_t lt = dataAsStr.size();
	vec timeData(lt);
	mat featureData(n, lt);
	
	for(size_t i=0; i<lt; i++)
	{
		dataRow = split(dataAsStr[i], ',');
		timeData(i) = std::stod(dataRow[0]);
		for(size_t j = 0; j<n; j++)
		{
			featureData(j,i) = std::stod(dataRow[j+1]);
		}
	}
}


