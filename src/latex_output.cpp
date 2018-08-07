#include <latex_output.h>

void parameterOutput(const mat& uvals)
{
	
	size_t m = uvals.rows(), 
		divs = 3;
	size_t nmax = uvals.cols()-1,
		n = ceil(nmax/divs);
	
	tableheader(divs+1);
	cout << "\t\t";
	for(size_t i=0; i<nmax-n; i+=n){
		cout << "$\\vecu_{"<< i <<"}$ " << "& ";
	}
	cout << "$\\vecu_{"<< nmax-1 <<"}$ "  << "& ";
	cout << "$\\vecu_$ " << "\\\\ \\midrule\n";
	
	for(size_t j=0; j<m; j++)
	{
		cout << "\t\t";
		for(size_t i=0; i<nmax-n; i+=n)
		{
			cout << uvals(j,i) << " & ";
		}
		cout << uvals(j,nmax-1) << " & ";
		cout << uvals(j,nmax) << " \\\\\n";
	}
	tablefooter();
}

void latexOutput(const mat& xn, const vec& u, int p, string buf){
	int n = xn.rows();
	int m = u.size();
	int lt = xn.cols();

	int interval = (int)(lt/5+.5);
	if(interval == 0){
		interval = 1;
	}
	int j;
	for(int i=0; i<n; i++){
		j = 0;
		cout << "$x_{"<< p <<"}$ " << buf << endl;
		while(j<lt){
			cout << xn(i, j) << buf << endl;
			j+=interval;
		}
	}
	cout << "$\\vecu_{"<< p <<"}$ " << buf << endl;
	for(int k=0; k<m; k++){
		cout << u(k) << buf << endl;
	}
}

void timelatexOutput(const vec& t, string buf, int n, int p){
	int lt = t.size();

	int interval = (int)(lt/5+.5);
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
	for(int k=0; k<p+1; k++){
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
		<< "\\sisetup{\n"
		<< "\tround-mode = places,\n"
		<< "\toutput-exponent-marker = \\text{e},\n"
		<< "}\n"
		<< "\t\\begin{tabularx}{\\linewidth}{T";
	for(int i=0; i<n; i++){
		cout <<  "S";
	}
	  cout << "}\n" \
		   << "\t\t\\toprule" << endl;
}
void tablefooter(){
	cout << "\t\t\\bottomrule" << endl;
	cout << "\t\\end{tabularx}" << endl;
	cout << "\\caption{ Table Caption }\n" \
		 <<	"\\label{}\n" \
		 << "\\end{table}" << endl;
}

void output_uNot_u_fig(string name, mat A, vec P, vec uNot, vec u, int N){
    ofstream outfile;
    string filename = "tiks/figs/" + name + std::to_string(N+1) + ".tex";
    outfile.open(filename);

    double gamma, dg, dvs = 100;
    int m = P.size();
    vec du(m), total(m);
    mat I = mat::Identity(m, m);
    mat B(m,m);

    outfile << "\\begin{figure}" << endl;
    outfile << "\\begin{tikzpicture}" << endl;
    outfile << "\\begin{semilogxaxis}[" << endl;
    outfile << "xlabel=$\\alpha$," << endl;
    outfile << "ylabel=$| \\vec{u}_{N-1} + \\Delta \\vec{u}_N  - \\vec{u}|$," << endl;
    outfile << "grid=major," << endl;
    outfile << "legend style={font=\\tiny}" << endl;
    outfile << "]" << endl;
    outfile << "\\addplot[color=red] coordinates {" << endl;

	total = uNot + A.inverse()*P - u;
	outfile << "(" << 0.0 << "," << total.norm() << ")" << endl;

    for(int j = -10; j<3; j++){
        gamma = 1*pow(10,j);
        dg = (1*pow(10,j)-1*pow(10,j-1))/dvs;
        for(int k = 0; k<(dvs+1); k++){
            B = A + gamma*I;

            du = B.inverse()*P;
            total = uNot + du - u;
            outfile << "(" << gamma << "," << total.norm() << ")" << endl;
            gamma += dg;
        }
    }
    outfile << "};" << endl;
    outfile << "\\end{semilogxaxis}" << endl;
    outfile << "\\end{tikzpicture}" << endl;
    outfile << "\\caption{Graph of the $\\alpha$ vs. $| \\vec{u}_{N-1} + \\Delta \\vec{u}_N  - \\vec{u}|$ at iteration $" << N+1 << "$.} " << endl;
    outfile << "\\label{fig:" << name << N+1 << "}" << endl;
    outfile << "\\end{figure}" << endl;

    outfile.close();
}

void g_output_uNot_u_fig(string name, mat A, vec P, vec uNot, vec u, vec ug, int N){
    ofstream outfile;
    string filename = "tiks/figs/" + name + std::to_string(N+1) + ".tex";
    outfile.open(filename);

    double gamma, dg, dvs = 100;
    int m = P.size();
    vec du(m), total(m);
    mat I = mat::Identity(m, m);
    mat B(m,m);

    outfile << "\\begin{figure}" << endl;
    outfile << "\\begin{tikzpicture}" << endl;
    outfile << "\\begin{semilogxaxis}[" << endl;
    outfile << "xlabel=$\\alpha$," << endl;
    outfile << "ylabel=$| \\vec{u}_{N-1} + \\Delta \\vec{u}_N  - \\vec{u}|$," << endl;
    outfile << "grid=major," << endl;
    outfile << "legend style={font=\\tiny}" << endl;
    outfile << "]" << endl;
    outfile << "\\addplot[color=red] coordinates {" << endl;
    for(int j = -14; j<2; j++){
        gamma = 1*pow(10,j);
        dg = (1*pow(10,j)-1*pow(10,j-1))/dvs;
        for(int k = 0; k<(dvs+1); k++){
            B = A + gamma*I;
            du = B.inverse()*(P+gamma*(ug - uNot));
            total = uNot + du - u;
            outfile << "(" << gamma << "," << total.norm() << ")" << endl;
            gamma += dg;
        }
    }
    outfile << "};" << endl;
    outfile << "\\end{semilogxaxis}" << endl;
    outfile << "\\end{tikzpicture}" << endl;
    outfile << "\\caption{Graph of the $\\alpha$ vs. $| \\vec{u}_{N-1} + \\Delta \\vec{u}_N  - \\vec{u}|$ at iteration $" << N+1 << "$.} " << endl;
    outfile << "\\label{fig:" << name << N+1 << "}" << endl;
    outfile << "\\end{figure}" << endl;

    outfile.close();
}
