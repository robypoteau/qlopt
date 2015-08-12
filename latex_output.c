#include "latex_output.h"

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

void longlatexOutput(const mat& M){
	mat otpt;
	otpt = M;
	
	int x = 3;
	int r = otpt.rows();
	int c = otpt.cols();
	
	check(c > x-1, "Not enough columns");
	
	if(c > 3){
		otpt.col(c-2) = colWiseMean(otpt.block(0,1,r,c-x)).col(0);
		otpt.col(c-1) = colWiseStdDev(otpt.block(0,1,r,c-x)).col(0);
	}
	
	tableheader(c-1);
	
	cout << "Actual &";
	for(int j=1; j<c-(x-1); j++){
		cout << "\t&";
	}
	cout << "Mean & Deviation \\\\" << endl;
	
	for(int i=0; i<r; i++){
		for(int j=0; j<c-1; j++){
			cout << otpt(i,j) << " & " ;
		}
		cout << otpt(i,c-1) << "\\\\" << endl;
	}
	
	tablefooter();
}

void shortlatexOutput(const mat& otpt){
	int c = otpt.cols();
	int r = otpt.rows();
	int x = 3;
	mat M(r,x);
	
	M << otpt.col(0), \
		colWiseMean(otpt.block(0,1,r,c-x)), \
		colWiseStdDev(otpt.block(0,1,r,c-x));
	
	longlatexOutput(M);	
}

vec colWiseStdDev(const mat& M){
	vec vect(M.rows());
	vect.fill(0);
	
	for(int i=0; i<M.cols()-1; i++){
		vect = vect.array() + (M.col(i) - M.col(M.cols()-1)).array().square();
	}
	
	return (vect/(M.cols()-1)).cwiseSqrt();
}

vec colWiseMean(const mat& M){
	return M.rowwise().sum()/(M.cols());
}

void tableheader(int n){
	cout << "\\begin{table}\n" \ 
		<< "\\centering\n" \
	    << "\t\\begin{tabular}{|l|";
	for(int i=0; i<n; i++){
		cout <<  "c |";
	}
	  cout << "}\n" \
		   << "\t\t\\hline" << endl;
}
void tablefooter(){
	cout << "\t\t\\hline" << endl;
	cout << "\t\\end{tabular}" << endl;
	cout << "\\caption{ Table Caption }\n" \
		 <<	"\\label{}\n" \
		 << "\\end{table}" ;
}