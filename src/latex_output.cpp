#include <latex_output.h>

void plotOutput(vec x, vec y)
{
	latexplot(x, y);
	pythonplot(x, y);
}

void pythonplot(vec x, vec y)
{
	//plot
	cout << "from matplotlib import pyplot as plt" << endl;
	cout << "x = ";
	convertVec(x);
	cout << "\ny = ";
	convertVec(y);
	cout << "\nplt.plot(x, y, marker=\'o\')" << endl;

	//legend and other text
	cout << "plt.xlabel(\'x label\')\n" \
		<< "plt.ylabel(\'y label\')\n" \
		<< "plt.title(\'Graph Title\')\n" \
		<< "plt.legend()" << endl;
	cout << "plt.show()" << endl;
}

void latexplot(vec x, vec y)
{
	cout << "\\begin{figure}\n" \
		<< "\\begin{tikzpicture}\n" \
		<< "\\begin{axis}[\n" \
		<< "\txlabel=X Label,\n" \
		<< "\tylabel=Y Label,\n" \
		<< "\tgrid=major,\n" \
		<< "\tlegend style={font=\\tiny}\n"
		<< "]\n"
		<< "\\addplot[color=red,grid=major] coordinates {\n	";
	for(int i=0; i<x.size(); i++)
	{
		cout << "("
			<< x(i)
			<< ","
			<< y(i)
			<< ")\n";
	}
	cout << "};\n" \
		<< "\\end{axis}\n" \
		<< "\\end{tikzpicture}\n" \
		<< "\\caption{Caption.}\n" \
		<< "\\label{fig:CHANGE_LABEL}\n" \
		<< "\\end{figure}" << endl;
}

void convertVec(vec x)
{
	cout << "[";
	for(int i = 0; i < x.size()-1; i++)
	{
		cout << x(i) << ",";
	}
	cout << x(x.size()-1);
	cout << "]";
}

void parameterOutput(const mat& uvals, const vec& alpha)
{

	size_t m = uvals.rows(),
		divs = 4;
	size_t nmax = uvals.cols()-1,
		n = ceil((double)nmax/(double)divs);

	tableheader(divs);
	cout << "\t\t$\\vecu & $";
	for(size_t i=0; i<nmax-n; i+=n){
		cout << "$\\vecu_{"<< i <<"}$ " << "& ";
	}
	cout << "$\\vecu_{"<< nmax-1 <<"}$ "  << "& ";
	cout << "$\\vecu$ " << "\\\\ \\midrule\n";

	for(size_t j=0; j<m; j++)
	{
		cout << "\t\t\\(p_{" << j+1 << "}\\) & ";
		for(size_t i=0; i<nmax-n; i+=n)
		{
			cout << uvals(j,i) << " & ";
		}
		cout << uvals(j,nmax-1) << " & ";
		cout << uvals(j,nmax) << " \\\\\n";
	}
	cout << "\t\t$\\alpha$ & & ";
	for(size_t i=n; i<nmax-n; i+=n)
	{
		cout << alpha(i-1) << " & ";
	}
	cout << alpha(nmax-2) << " & ";
	cout << " \\\\\n";

	tablefooter();
}

void tableheader(int n){
	cout << "\\begin{table}\n"
		<< "\\centering\n"
		<< "\\begin{adjustbox}{clip=false}\n"
		<< "\t\\begin{tabular}{l";
	for(int i=0; i<n; i++){
		cout <<  "S";
	}
	  cout << "}\n" \
		   << "\t\t\\toprule" << endl;
}

void tablefooter(){
	cout << "\t\t\\bottomrule" << endl;
	cout << "\t\\end{tabular}" << endl;
	cout << "\\end{adjustbox}" << endl;
	cout << "\\caption{ Table Caption }\n" \
		 <<	"\\label{tab:tablename}\n" \
		 << "\\end{table}" << endl;
}
