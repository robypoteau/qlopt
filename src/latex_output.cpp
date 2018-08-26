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

