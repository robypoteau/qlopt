#include <latex_output.h>

void latexOutput(const mat& xn, const vec& u, int p, string buf){
	int n = xn.rows();
	int m = u.size();
	int lt = xn.cols();

	int interval = (int)(lt/10+.5);
	if(interval == 0){
		interval = 1;
	}
	int j;
	for(int i=0; i<n; i++){
		j = 0;
		cout << "$\\vec{x}_{"<< p <<"}(t,u)$ " << buf << endl;
		while(j<lt){
			cout << xn(i, j) << buf << endl;
			j+=interval;
		}
	}
	cout << "\t "<< buf << endl;
	for(int k=0; k<m; k++){
		cout << u(k) << buf << endl;
	}
}

void timelatexOutput(const vec& t, string buf, int n, int p){
	int lt = t.size();

	int interval = (int)(lt/10+.5);
	if(interval == 0){
		interval = 1;
	}

	int j;
	for(int i=0; i<n; i++){
		j = 0;
		cout << "time " << buf << endl;
		while(j<lt){
			cout << t(j) << buf << endl;
			j+=interval;
		}
	}
	for(int k=0; k<p; k++){
		cout << "\t " << buf << endl;
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
	cout << "Mean & Deviation \\\\ \\hline" << endl;

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

void R(double dt, const mat& otpt, int n){
	mat M = otpt;
	int r = otpt.rows();
	int N = r/n;
//note("N should equal -k param+1");
//note(N);
	double temp, temp2, ans = 0, mean, vmin=numeric_limits<double>::max(), vmax=0;
	vec v;
	vec v2(N-1);

	for(int i=1;i<N;i++){
		M.middleRows(i*n,n) -= M.topRows(n);
	}

	M = M.array().square();
	M *= dt;
	v = M.rowwise().sum();

	for(int i=1;i<N;i++){
		temp = 0;
		for(int j=0;j<n;j++){
			temp += v(i*n+j);
		}
		temp2 = sqrt(temp);
		v2(i-1) = temp2;
		ans += temp2;
		vmax = max(vmax,temp2);
		vmin = min(vmin,temp2);
	}
	mean = ans/N;
	ans = 0;
	for(int i=1;i<N;i++){
		temp = 0;
		for(int j=0;j<n;j++){
			temp += v(i*n+j);
		}
		ans += pow(sqrt(temp)-mean,2);
	}
	cout << "Error Data Vector\n" << v2.transpose() << endl;
	cout << "L2N= " << mean << "\nstd = " << sqrt(ans/(N-1)) << endl;
	cout << "Max = " << vmax << "\nMin = " << vmin << endl;
}
void M(const mat& otpt, int n){
	mat M = otpt;
	int r = otpt.rows();
	int c = otpt.cols();
	int dim = r/n;
	double temp;
	//note("Output\n");
	//note(M);
	for(int i=1;i<(r/n);i++){
		M.middleRows(i*n,n) -= M.topRows(n);
	}
	M.topRows(n) -= M.topRows(n);
	//note("Approximations Rows - Actual Data:\n");
	//note(M);
	M = M.array().abs();
	temp = M(dim,0);
	note(temp);
	for(int i=dim;i<r;i++){
		for(int j=0;j<c;j++){
			if(temp < M(i,j)){
				temp = M(i,j);
			}
		}
	}
	//note("\n");
	//note(M);
	cout << "Max-Norm = " << temp << endl;
}

void Mi(const mat& otpt, int n){
	mat M = otpt;
	int r = otpt.rows();
	int c = otpt.cols();
	int dim = r/n;
	double temp;
	for(int i=1;i<dim;i++){
		M.middleRows(i*n,n) -= M.topRows(n);
	}
	M.topRows(n) -= M.topRows(n);
	M = M.array().abs();
	temp = M(dim,0);
	for(int i=dim;i<r;i++){
		for(int j=0;j<c;j++){
			if(temp > M(i,j)){
				temp = M(i,j);
			}
		}
	}
	cout << "Min-Norm = " << temp << endl;
}

void shortNormalizedLatexOutput(const mat& M){
	int r = M.rows();
	int c = M.cols();
	int x = 3;

	mat otpt(r,c); otpt = M;
	mat N(r,x);

	for(int i=1; i<c; i++){
		otpt.col(i) -= otpt.col(0);
	}

	N << otpt.col(0), \
		colWiseMean(otpt.block(0,1,r,c-x)), \
		colWiseStdDev(otpt.block(0,1,r,c-x));

	longlatexOutput(N);
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
	cout << "\\begin{table}\n"
		<< "\\centering\n"
	    << "\t\\begin{tabularx}{|X|";
	for(int i=0; i<n; i++){
		cout <<  "C|";
	}
	  cout << "}\n" \
		   << "\t\t\\hline" << endl;
}
void tablefooter(){
	cout << "\t\t\\hline" << endl;
	cout << "\t\\end{tabularx}" << endl;
	cout << "\\caption{ Table Caption }\n" \
		 <<	"\\label{}\n" \
		 << "\\end{table}" << endl;
}
