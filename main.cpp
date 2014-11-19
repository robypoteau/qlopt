#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <Eigen/LU>
#include <Eigen/QR>
#include <unsupported/Eigen/MPRealSupport>

#include "misc.h"
#include "nonlinear_odes.h"
#include "numerical_integration.h"
#include "thesis_functions.h"
#include "input.h"
#include "regs.h"

using namespace std;
using namespace mpfr;
using namespace Eigen;
using namespace thesis;

mat reshape(const mat& U, int n, int m);
void latexOutput(const mat& xn, const vec& u, int p, string buf);
mpreal cond(const mat& A);

int main(int argc, char *argv[])
{
	time_t begin;
	begin = time(NULL);
	
	const int bits = 128;	//const int bytes = bits/8;
	mpreal::set_default_prec(bits);
	
	nonlinearOdes no;
	input in(argc, argv);
	
	//Gather input from the console
	bool reg = in.isRegularized();
	
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
	//End input
	
	int lt = t.size();  
	int m = uNot.size();
	int n = yNot.size();
	
	mat measure;
	measure = rungekutta4(no.odeFuncMap[system], times, u, yNot);
	
	spline msmtRow[n];
	for(int i=0; i<n; i++){
		msmtRow[i].update(times, measure.row(i));
	}
	
	mat msmt(n,lt);
	for(int i = 0; i<n; i++){
		for(int j = 0; j<lt; j++){
			msmt(i,j) = msmtRow[i].interpolate(t(j));
		}
	}	
	//refactor into a function		
	vec zeros(n*m); 
	zeros.fill(0);
	vec lyNot((m+1)*n);
	//lyNot << msmt(0,0), msmt(0,1), zeros;// update where yNot gets its values
	lyNot << yNot, zeros;// update where yNot gets its values

	mat xNminus(n, lt);
	xNminus << msmt; 
	
	mat bob(lyNot.size(), lt);
	mat U(zeros.size()/* n*m */, n*lt);
	mat A(m,m);
	vec P(m);
	vec du(m);
	
	cout.precision(5);
	
	int interval = (int)(lt/10+.5);
	if(interval == 0){
		interval = 1;
	}
	
	int l = 0;
	//cout << "N: &" << endl;
	cout << "time: &" << endl;
	while(l<lt){
		cout << t(l) << " & " << endl;
		l += interval;
	}
	cout << "time: &" << endl;
	l=0;
	while(l<lt){
		cout << t(l) << " & " << endl;
		l += interval;
	}
	cout << "\t &" << endl;
	for(int k=0; k<m; k++){
		cout << "$u_{" << k+1 << "}$ & " << endl;
	}
	
	regs r;
	mat next(n, lt);
	mat last(n, lt);
	mpreal gamma;
	mpreal p, p0;
	
	//Run the rest of the iterations

	for(int i = 0; i<170; i++)
	{
		bob = qLinearRungeKutta4(no.odeFuncMap[system], t, uNot, lyNot, xNminus);
		//bob = rungekutta4(no.qLinFuncMap[system + "_linearization"], t, uNot, lyNot, xNminus); 
		//cout << bob << endl << endl;
		
		U = reshape(bob.bottomRows(zeros.size()/* n*m */), m, n*lt);
		xNminus = bob.topRows(n);
		
		A = findA(t, U, m);
		cout << "cond(A) = " << cond(A) << "\ndet(A) = " << A.determinant() << "\nrank(A) = " << A.fullPivHouseholderQr().rank() <<endl;
		cout << A << endl;
		
		P = findP(t, U, reshape(msmt - xNminus, 1, n*lt).row(0), m);
		//cout << P << endl;
		cout << "deltau\n" << A.inverse()*P << endl;
		
		if((reg && cond(A) > 1E+20) /*|| isnan(cond(A))*/){
				gamma = 1E-24;
				r.update(A, P);
				//next 
				last = xNminus;
				do{
					//last = next;
					gamma *=.9;
					du = r.regularization(gamma);
					bob = qLinearRungeKutta4(no.odeFuncMap[system], t, uNot, lyNot, xNminus);
					//cout << uNot+du << endl;
					//bob = rungekutta4(no.qLinFuncMap[system + "_linearization"], t, du, lyNot, xNminus);
					U = reshape(bob.bottomRows(zeros.size()), m, n*lt);
					next = bob.topRows(n);
					p0 = r.distanceMeasure(msmt, last);
					p = r.distanceMeasure(msmt, next);
					cout << p << ",";
				}while((p0 < p));
				cout << "gamma:\n" << gamma <<endl;
				cout << "regu:\n" << du <<endl;
		}else{
			//du = A.fullPivHouseholderQr().solve(P);  //A.lu().solve(P);
			du = A.inverse()*P;
		}
		uNot += du;

		latexOutput(xNminus, uNot, i+1, " & "); //cout << "n = " << i <<":\n" << xNminus << "\nParameter Estimates\n"<< uNot.transpose() << endl << endl;
		du(0) = du.norm();
		
		if(du(0) < 0.00001 || isnan(du(0)) ){
			break;
		}
	}
	latexOutput(msmt, u, -1, " \\\\ ");
	
	time_t end;
	end = time(NULL);
	cout << end - begin << endl;
	return 0;
}

mat reshape(const mat& U, int n, int m)
{
	mat newU(n,m);
	newU.fill(0);
	int olt = U.row(0).size(); 	//old time(row) length
	int on = m/olt;				//on col length
	
	for(int i = 0; i<n; i++){
		for(int j = 0; j<on; j++){
			newU.block(i, j*olt, 1, olt) = U.row(i*on + j);
		}
	}
	
	return newU;
}

void latexOutput(const mat& xn, const vec& u, int p, string buf){
	int n = xn.rows();
	int m = u.size();
	int N = xn.cols();
	
	int interval = (int)(N/10+.5);
	if(interval == 0){
		interval = 1;
	}
	int j;
	for(int i=0; i<n; i++){
		j = 0;
		cout << "$\\vec{x}_{"<< p <<"}(t,u)$ " << buf << endl;
		while(j<N){
			cout << xn(i, j) << buf << endl;
			j+=interval;
		}
	}
	cout << "\t "<< buf << endl;
	for(int k=0; k<m; k++){
		cout << u(k) << buf << endl;
	}
}

mpreal cond(const mat& A){
	return A.norm()*A.inverse().norm();
}

/**/