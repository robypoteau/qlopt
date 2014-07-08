#ifndef MISC_H
#define MISC_H

#include <mpreal.h>
#include <Eigen/Dense>
#include <string>
	
using namespace mpfr;
using namespace Eigen;
using namespace std;
	
typedef Matrix<mpreal, Dynamic, 1> vec;
typedef Matrix<mpreal, Dynamic, Dynamic> mat;
	


	/*struct solverInput{
		string funcName;
		vec time;
		vec msmt;
		vec uNot;
		vec yNot;
	} ;
	
	solverInput getAll();
	string getFuncName();
	vec getTimeInterval(mpreal a, mpreal b, int n);
	vec getMsmt();
	vec getInitParam();
	vec getInitVal();*/

#endif