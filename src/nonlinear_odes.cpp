#include <nonlinear_odes.h>
#include <complex>

namespace thesis{

	nonlinearOdes::nonlinearOdes(){
		odeFuncMap["lotka_volterra"] = nonlinearOdes::lotka_volterra;
		odeFuncMap["lotka4"] = nonlinearOdes::lotka4;
		odeFuncMap["general_lv"] = nonlinearOdes::general_lv;
		odeFuncMap["lorenz"] = nonlinearOdes::lorenz;
		qLinFuncMap["lotka_volterra_linearization"] = nonlinearOdes::lotka_volterra_linearization;
		odeFuncMap["pielou"] = nonlinearOdes::pielou;
		//qLinFuncMap["pielou_linearization"] = nonlinearOdes::pielou_linearization;
		odeFuncMap["angiogenesis"] = nonlinearOdes::angiogenesis;
		odeFuncMap["cancer"] = nonlinearOdes::cancer;
		//qLinFuncMap["angiogenesis_linearization"] = nonlinearOdes::angiogenesis_linearization;
		//odeFuncMap["coral"] = nonlinearOdes::coral;
		odeFuncMap["coral5"] = nonlinearOdes::coral5;
		//odeFuncMap["coral_pw"] = nonlinearOdes::coral_pw;
		//odeFuncMap["coral_two"] = nonlinearOdes::coral_two;
		//odeFuncMap["coral_four"] = nonlinearOdes::coral_four;
		//qLinFuncMap["coral_linearization"] = nonlinearOdes::coral_linearization;
		odeFuncMap["bistable_switch"] = nonlinearOdes::bistable_switch;
		odeFuncMap["bistable_switch_two"] = nonlinearOdes::bistable_switch_two;
		qLinFuncMap["bistable_switch_linearization"] = nonlinearOdes::bistable_switch_linearization;
		odeFuncMap["eight_part"] = nonlinearOdes::eight_part;
		odeFuncMap["eight_part_spc"] = nonlinearOdes::eight_part_spc;
		//qLinFuncMap["eight_part_linearization"] = nonlinearOdes::eight_part_linearization;
		odeFuncMap["toggle_switch"] = nonlinearOdes::toggle_switch;
		odeFuncMap["toggle_switch_config1"] = nonlinearOdes::toggle_switch_config1;
		odeFuncMap["toggle_switch_config2"] = nonlinearOdes::toggle_switch_config2;
		odeFuncMap["toggle_switch_config3"] = nonlinearOdes::toggle_switch_config3;
		odeFuncMap["toggle_switch_config4"] = nonlinearOdes::toggle_switch_config4;
		odeFuncMap["repressilator"] = nonlinearOdes::repressilator;
		odeFuncMap["general_repressilator"] = nonlinearOdes::general_repressilator;
	}


	mat nonlinearOdes::lotka_volterra(const double& t, const vec& x, const vec& u)
	{
		(void) t;
		mat result(2,1);
		result(0) = u(0)*x(0) - u(1)*x(0)*x(1);
		result(1) = u(1)*x(0)*x(1) - u(2)*x(1);
		
		return result;
	}

	mat nonlinearOdes::lotka4(const double& t, const vec& x, const vec& u)
	{
		mat result(2,1);
		result(0) = u(0)*x(0) - u(1)*x(0)*x(1);
		result(1) = u(2)*x(0)*x(1) - u(3)*x(1);

		return result;
	}

	mat nonlinearOdes::general_lv(const double& t, const vec& x, const vec& u)
	{
		size_t n = x.size();
		check(n*(n+1) == u.size(), "In general_lv the parameter vector u is the wrong dimension.");

		mat result(n,1);
		for(size_t i=0; i<n; i++)
		{
			result(i) = u(n*i)*x(i);
			for(size_t j=0; j<n; j++)
			{
				result(i) += u(n*i + j+1)*x(i)*x(j);
			}
		}

		return result;
	}

	mat nonlinearOdes::lorenz(const double& t, const vec& x, const vec& u)
	{
		mat result(3,1);
		result(0) = u(0)*(x(1) - x(0));
		result(1) = x(0)*(u(1) - x(2)) - x(1);
		result(2) = x(0)*x(1) - u(2)*x(2);
		nonlinearOdes::addCounter();
		return result;
	}
	
	mat nonlinearOdes::lotka_volterra_linearization(const double& t, const vec& x, const vec& u, std::vector<thesis::spline> xn)
	{
		double xnt = xn[0].interpolate(t);
		double ynt = xn[1].interpolate(t);

		mat result(8,1);
		result << (u(0)-u(1)*ynt)*x(0) - u(1)*xnt*x(1) + u(1)*xnt*ynt,
				 u(1)*ynt*x(0) + (u(1)*xnt - u(2))*x(1) - u(1)*xnt*ynt,
				(u(0) - u(1)*ynt)*x(2) + x(0) - u(1)*xnt*x(3),
				 u(1)*ynt*x(2) + (u(1)*xnt - u(2))*x(3),
				(u(0) - u(1)*ynt)*x(4) + ynt*(xnt - x(0)) - xnt*x(1) - u(1)*xnt*x(5),
				 u(1)*ynt*x(4) + ynt*x(0) + (u(1)*xnt - u(2))*x(5) + xnt*(x(1) - ynt),
				(u(0) - u(1)*ynt)*x(6) - u(1)*xnt*x(7),
				 u(1)*ynt*x(6) + (u(1)*xnt - u(2))*x(7) - x(1);

		return result;
	}

	mat nonlinearOdes::pielou(const double& t, const vec& x, const vec& u)
	{
		double k = 250;

		mat result(2,1);
		result(0) = u(0)*(1-x(0)/k)*x(0) - u(1)*x(0)*x(1);
		result(1) = -u(2)*x(1) + u(3)*x(0)*x(1);

		return result;
	}

	/*mat nonlinearOdes::pielou_linearization(const double& t, const vec& x, const vec& u, const mat& xn, const vec& time)
	{
		double r = u(0);
		double b = u(1);
		double c = u(2);
		double d = u(3);

		double k = 250;

		thesis::spline Xn(time, xn.row(0));
		thesis::spline Yn(time, xn.row(1));

		double x1 = Xn.interpolate(t);
		double x2 = Yn.interpolate(t);

		mat result(10,1);
		result << (r*(x1*x1+k*x(1-1)-2*x1*x(1-1))-b*k*(x(1-1)*x2+x1*(-x2+x(2-1))))/k,
			-(c*x(2-1))+d*(x(1-1)*x2+x1*(-x2+x(2-1))),
			((k-2*x1)*x(1-1)+(k*r-2*r*x1-b*k*x2)*x(3-1)+x1*(x1-b*k*x(4-1)))/k,
			d*x2*x(3-1)+(-c+d*x1)*x(4-1),
			x1*x2-x2*x(1-1)-x1*x(2-1)+r*x(5-1)-(2*r*x1*x(5-1))/k-b*x2*x(5-1)-b*x1*x(6-1),
			d*x2*x(5-1)+(-c+d*x1)*x(6-1),
			((-2*r*x1+k*(r-b*x2))*x(7-1)-b*k*x1*x(8-1))/k,
			-x(2-1)+d*x2*x(7-1)+(-c+d*x1)*x(8-1),
			((-2*r*x1+k*(r-b*x2))*x(9-1)-b*k*x1*x(10-1))/k,
			-(x1*x2)+x2*x(1-1)+x1*x(2-1)+d*x2*x(9-1)-c*x(10-1)+d*x1*x(10-1);

		return result;
	}*/

	mat nonlinearOdes::angiogenesis(const double& t, const vec& x, const vec& u)
	{
		double m = u(0); // m>0
		double n = 1; //u(1); // n = .25, = 1
		double a = u(2); // a>=0
		double c = u(3); // c(t) = c = 0 => no treatment, c = const implies fixed treatment
		double w = u(4); // w>=0
		double g = u(5); // g>=0

		/*
		double m = u(0); // m>0
		double n = .25; // n = .25, = 1
		double a = u(1); // a>=0
		double c = u(2); // c(t) = c = 0 => no treatment, c = const implies fixed treatment
		double w = u(3); // w>=0
		double g = u(4); // g>=0
		*/

		double f;
		if(n == 0){
			f = -m*log(x(0)/x(1));
		}else{
			f = m*x(0)/n*(1-pow(x(0)/x(1),n));
		}

		mat result(2,1);
		result(0) = -a*c*x(1) + f;
		result(1) = -a*c*x(1) + w*x(0) - g*pow(x(0)*x(0),1/3)*x(1);

		return result;
	}

	/*mat nonlinearOdes::angiogenesis_linearization(const double& t, const vec& x, const vec& u, const mat& xn, const vec& time)
	{
		double m = u(0);
		double n = u(1);
		double a = u(2);
		double c = u(3);
		double w = u(4);
		double g = u(5);

		thesis::spline Xn(time, xn.row(0));
		thesis::spline Yn(time, xn.row(1));

		double x1 = Xn.interpolate(t);
		double x2 = Yn.interpolate(t);

		mat result(14,1);
		result << m*x1*pow(x2,-1 - n)*(-((2 + n)*pow(x1,n)*x(1-1)*x2) - x1*pow(x2,1 + n) + 2*x(1-1)*pow(x2,1 + n) + pow(x1,1 + n)*(x2 + n*x(2-1))),
		w*x(1-1) + (-3*a*c*pow(x1,1/3)*x(2-1) + g*(2*x1*x2 - 2*x(1-1)*x2 - 3*x1*x(2-1)))/(3*pow(x1,1/3)),
		x1*pow(x2,-1 - n)*(-(x1*pow(x2,1 + n)) - (2 + n)*pow(x1,n)*x2*x(1-1) + 2*pow(x2,1 + n)*x(1-1) + pow(x1,1 + n)*(x2 + n*x(2-1)) + m*(-(x2*((2 + n)*pow(x1,n) - 2*pow(x2,n))*x(3-1)) + n*pow(x1,1 + n)*x(4-1))),
		((3*w*pow(x1,1/3) - 2*g*x2)*x(3-1) - 3*(a*c + g*pow(x1,2/3))*pow(x1,1/3)*x(4-1))/(3*pow(x1,1/3)),
		m*x1*pow(x2,-1 - n)*(pow(x1,1 + n)*x2*log(x1) - pow(x1,1 + n)*x2*log(x2) - pow(x1,n)*x2*(1 + (2 + n)*log(x1) - (2 + n)*log(x2))*x(1-1) + pow(x1,1 + n)*(1 + n*log(x1) - n*log(x2))*x(2-1) - 2*pow(x1,n)*x2*x(5-1) - n*pow(x1,n)*x2*x(5-1) + 2*pow(x2,1 + n)*x(5-1) + n*pow(x1,1 + n)*x(6-1)),
		((3*w*pow(x1,1/3) - 2*g*x2)*x(5-1) - 3*(a*c + g*pow(x1,2/3))*pow(x1,1/3)*x(6-1))/(3*pow(x1,1/3)),
		m*x1*pow(x2,-1 - n)*(-(x2*((2 + n)*pow(x1,n) - 2*pow(x2,n))*x(7-1)) + n*pow(x1,1 + n)*x(8-1)),
		(-3*c*pow(x1,1/3)*x(2-1) + (3*w*pow(x1,1/3) - 2*g*x2)*x(7-1) - 3*(a*c + g*pow(x1,2/3))*pow(x1,1/3)*x(8-1))/(3*pow(x1,1/3)),
		m*x1*pow(x2,-1 - n)*(-(x2*((2 + n)*pow(x1,n) - 2*pow(x2,n))*x(9-1)) + n*pow(x1,1 + n)*x(10-1)),
		(-3*a*pow(x1,1/3)*x(2-1) + (3*w*pow(x1,1/3) - 2*g*x2)*x(9-1) - 3*(a*c + g*pow(x1,2/3))*pow(x1,1/3)*x(10-1))/(3*pow(x1,1/3)),
		m*x1*pow(x2,-1 - n)*(-(x2*((2 + n)*pow(x1,n) - 2*pow(x2,n))*x(11-1)) + n*pow(x1,1 + n)*x(12-1)),
		(3*pow(x1,1/3)*x(1-1) + (3*w*pow(x1,1/3) - 2*g*x2)*x(11-1) - 3*(a*c + g*pow(x1,2/3))*pow(x1,1/3)*x(12-1))/(3*pow(x1,1/3)),
		m*x1*pow(x2,-1 - n)*(-(x2*((2 + n)*pow(x1,n) - 2*pow(x2,n))*x(13-1)) + n*pow(x1,1 + n)*x(14-1)),
		w*x(13-1) + (2*x1*x2 - 2*x2*x(1-1) - 3*x1*x(2-1) - 2*g*x2*x(13-1) - 3*a*c*pow(x1,1/3)*x(14-1) - 3*g*x1*x(14-1))/(3*pow(x1,1/3));

		return result;
	}*/

	mat nonlinearOdes::cancer(const double& t, const vec& x, const vec& u)
	{
		double d   = u(0); // >=0
		double tau = u(1); // >=0
		double g   = u(2); // >=0
		double a   = u(3); // >=0, a <= g/2
		double k   = u(4); // >=0
		double R   = 1.0;
		mat result(2,1);
		result(0) = d*R - x(0)/tau - g*x(0)*x(0);
		result(1) = -(a*R + k*x(0))*x(1);

		return result;
	}

	/* mat nonlinearOdes::coral(const double& t, const vec& x, const vec& u)
	{
		//input
		double m = u(0);
		double r = u(1);
		double k1 = u(2);
		double p = u(3);
		double k2 = u(4);
		double k3 = u(5);

		//Pre-set parameter
		double a2 = 1.5;
		double b2 = .0013;
		double N = 1750;
		double I0 = 2000;
		double pi = atan(1)*4;

		mat result(2,1);
		result(0) = m*x(1-1) + r*x(1-1)*(1 - x(1-1)/N) - k1*x(2-1)*x(1-1) - p*x(1-1);
		result(1) = k2*x(2-1)+ k3*a2*b2*pow(2,b2*sin(pi*t/12))*I0*log(2)/12*pi*cos(pi*t/12);

		return result;
	} */

	mat nonlinearOdes::coral5(const double& t, const vec& x, const vec& u)
	{
		//input
		double m = u(0);
		double r = u(1);
		double k1 = u(2);
		//double p = u(0);
		double k2 = u(3);
		double k3 = u(4);

		//Pre-set parameter
		double a2 = 1.5;
		double b2 = .0013;
		double N = 1750;
		double I0 = 2000;
		double pi = atan(1)*4;

		mat result(2,1);
		result(0) = m*x(1-1) + r*x(1-1)*(1 - x(1-1)/N) - k1*x(2-1)*x(1-1);// - p*x(1-1);
		result(1) = k2*x(2-1)+ k3*a2*b2*pow(2,b2*sin(pi*t/12))*I0*log(2)/12*pi*cos(pi*t/12);

		return result;
	}

	/* mat nonlinearOdes::coral_pw(const double& t, const vec& x, const vec& u)
	{
		//input
		double m = u(0);
		double r = u(1);
		double k1 = u(2);
		double p = u(3);
		double k2 = u(4);
		double k3 = u(5);

		//Pre-set parameter
		double a2 = 1.5;
		double b2 = .0013;
		double N = 1750;
		double I0 = 2000;
		double pi = atan(1)*4;

		mat result(2,1);
		result(0) = m*x(1-1) + r*x(1-1)*(1 - x(1-1)/N) - k1*x(2-1)*x(1-1) - p*x(1-1);
		result(1) = k2*x(2-1);

		/*double modulus = mod(t,24);//= t/24 - floor(t/24);
		if(modulus >= 0 && modulus < 12){
			result(1) += k3*a2*b2*pow(2,b2*sin(pi*t/12))*I0*log(2)/12*pi*cos(pi*t/12);
		}*/
		//return result;
	//}*/

	/* mat nonlinearOdes::coral_two(const double& t, const vec& x, const vec& u)
	{
		//input
		double m = u(0);
		double r = u(1); // 0.00125; //u(1);
		double k1 = u(2); // 0.000001;// u(2);
		double p = u(3); //u(3);
		double k2 = u(4); // 0.015;//u(4);
		double k3 = u(5); // .0025;//u(5);

		//Pre-set parameter
		double a2 = u(6); // 1.5;
		double b2 = u(7); //.0013;
		double N = u(8); // 1750;
		double I0 = u(9); //2000;
		double pi = atan(1)*4;

		mat result(2,1);
		result(0) = m*x(1-1) + r*x(1-1)*(1 - x(1-1)/N) - k1*x(2-1)*x(1-1) - p*x(1-1);
		result(1) = k2*x(2-1);

		/*double modulus = mod(t,24);//= t/24 - floor(t/24);
		if(modulus >= 0 && modulus < 13){
			result(1) += k3*a2*b2*pow(2,b2*sin(pi*t/12))*I0*log(2)/12*pi*cos(pi*t/12);
		}*/
		/*return result;
	}

		mat nonlinearOdes::coral_four(const double& t, const vec& x, const vec& u)
	{
		//input
		double m = .02;
		double r = u(0);
		double k1 = u(1);
		double p = .02;
		double k2 = u(2);
		double k3 = u(3);

		//Pre-set parameter
		double a2 = 1.5;
		double b2 = .0013;
		double N = 1750;
		double I0 = 2000;
		double pi = atan(1)*4;

		mat result(2,1);
		result(0) = m*x(1-1) + r*x(1-1)*(1 - x(1-1)/N) - k1*x(2-1)*x(1-1) - p*x(1-1);
		result(1) = k2*x(2-1)+ k3*a2*b2*pow(2,b2*sin(pi*t/12))*I0*log(2)/12*pi*cos(pi*t/12);

		return result;
	} */

	/*mat nonlinearOdes::coral_linearization(const double& t, const vec& x, const vec& u, const mat& xn, const vec& time)
	{
		//input
		double m = u(0);
		double r = u(1);
		double k1 = u(2);
		double p = u(3);
		double k2 = u(4);
		double k3 = u(5);

		//Pre-set parameter
		double a2 = 1.5;
		double b2 = .0013;
		double I0 = 2000;
		double N = 1750;

		thesis::spline Xn(time, xn.row(0));
		thesis::spline Yn(time, xn.row(1));

		double x1 = Xn.interpolate(t);
		double x2 = Yn.interpolate(t);

		double pi = atan(1)*4;

		mat result(14,1);
		result <<
				m*x1 + r*x1*(1 - x1/N) - k1*x1*x2 - p*x2 + (m + r - 2*r*x1/N - k1*x2)*(x(1-1)-x1) - (k1*x1 + p)*(x(2-1) - x2),
				k2*x2 + a2*b2*k3*pow(2,b2*I0*sin(pi*t/12))*log(2)*pi*I0/12*cos(pi*t/12) + k2*(x(2-1) - x2),
				x1 + (x(1-1)-x1) + (m + r - 2*r*x1/N - k1*x2)*(x(3-1)) - (k1*x1 + p)*(x(4-1)),
				k2*x(4-1),
				x1*(1 - x1/N) + (x(1-1)-x1) + (m + r - 2*r*x1/N - k1*x2)*(x(5-1)) - (k1*x1 + p)*(x(6-1)),
				k2*x(6-1),
				-x1*x2 - x2*(x(1-1)-x1) + (m + r - 2*r*x1/N - k1*x2)*(x(7-1)) - (k1*x1 + p)*(x(8-1)),
				k2*x(8-1),
				(m + r - 2*r*x1/N - k1*x2)*x(9-1) - x(2-1) - (k1*x1 + p)*(x(10-1)),
				k2*x(10-1),
				(m + r - 2*r*x1/N - k1*x2)*x(11-1) - (k1*x1 + p)*x(12-1),
				x(2-1) + k2*x(12-1),
				(m + r - 2*r*x1/N - k1*x2)*x(13-1) - (k1*x1 + p)*x(14-1),
				a2*b2*pow(2,b2*I0*sin(pi*t/12))*log(2)*pi*I0/12*cos(pi*t/12)+ k2*x(14-1);
		return result;
	}*/

	mat nonlinearOdes::bistable_switch(const double& t, const vec& x, const vec& p)
	{
		double alpha = p(0);
		double u     = p(1);
		double n     = p(2);

		mat result(2,1);
		result(0) = alpha/(1 + pow(u*x(2-1),n)) - x(1-1);
		result(1) = alpha/(1 + pow(x(1-1),n))- x(2-1);
		
		nonlinearOdes::addCounter();
		return result;
	}

	mat nonlinearOdes::bistable_switch_two(const double& t, const vec& x, const vec& p)
	{
		double alpha = p(0);
		double u     = 3.2;
		double n     = p(1);

		mat result(2,1);
		result(0) = alpha/(1 + pow(u*x(1),n)) - x(0);
		result(1) = alpha/(1 + pow(x(0),n))- x(1);
				nonlinearOdes::addCounter();
		return result;
	}

	mat nonlinearOdes::bistable_switch_linearization(const double& t, const vec& x, const vec& p, std::vector<thesis::spline> xn)
	{
		double alpha = p(0);
		double u 	 = p(1);
		double n 	 = p(2);

		double A = xn[0].interpolate(t);
		double B = xn[1].interpolate(t);

		double g;
		double d;

		g = 1 + pow(u*B,n);
		d = 1 + pow(A,n);
		
		mat result(8,1);
		result << alpha/g - x(1-1) - alpha*n*pow(u,n)*pow(B,n-1)/pow(g,2)*(x(2-1)-B),
				  alpha/d - x(2-1) - alpha*n*pow(A,n-1)/pow(d,2)*(x(1-1)-A),
			// wrt alpha
			1/g - x(3-1) - n*pow(u,n)*pow(B,n-1)/pow(g,2)*(x(2-1)-B) - alpha*n*pow(u,n)*pow(B,n-1)/pow(g,2)*(x(4-1)),
			1/d - x(4-1) - n*pow(A,n-1)/pow(d,2)*(x(1-1)-A) - alpha*n*pow(A,n-1)/pow(d,2)*(x(3-1)),
			//wrt u
			-alpha*n*pow(u,n-1)*pow(B,n)/pow(g,2) - x(5-1) - (alpha*n*n*pow(u*B,n-1)*g - 2*alpha*n*pow(u*B,2*n-1))/pow(g,3)*(x(2-1)-B) - alpha*n*pow(u,n)*pow(B,n-1)/pow(g,2)*(x(6-1)),
			-x(6-1) - alpha*n*pow(A,n-1)/pow(d,2)*(x(5-1)),
			// wrt n
			-alpha*pow(u*B,n)*log(u*B)/pow(g,2) - x(7-1) - alpha*n*pow(u,n)*pow(B,n-1)/pow(g,2)*((1/n + log(u*B) - 2*pow(u*B,n)*log(u*B)/g)*(x(2-1)-B) + x(8-1)),
			-alpha*pow(A,n)*log(A)/pow(d,2) - x(8-1) - alpha*n*pow(A,n-1)/pow(d,2)*((1/n + log(A) - 2*pow(A,n)*log(A)/d)*(x(1-1)-A) + x(7-1));
		
		nonlinearOdes::addCounter();
		return result;
	}
	
	mat nonlinearOdes::eight_part(const double& t, const vec& x, const vec& pp)
	{
		double P = .05;
		double S = 10;

		vec p(36);
		if(pp.size() < 36){
			p << 1.0,1.0,2.0,1.0,2.0,1.0,1.0,1.0,2.0,1.0,2.0,1.0,1.0,1.0,2.0,1.0,2.0,1.0,0.1,1.0,0.1,0.1,1.0,0.1,0.1,1.0,0.1,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0;
			for(int i=0; i<pp.size(); i++){
				p(i) = pp(i);
			}
		}else{
			p << pp;
		}
		mat result(8,1);
		result << p(1-1)/(1 + pow(P/p(2-1), p(3-1))   + pow(p(4-1)/S      , p(5-1)))  - p(6-1)*x(1-1),
				  p(7-1)/(1 + pow(P/p(8-1), p(9-1))   + pow(p(10-1)/x(7-1), p(11-1))) - p(12-1)*x(2-1),
				  p(13-1)/(1 + pow(P/p(14-1),p(15-1)) + pow(p(16-1)/x(8-1), p(17-1))) - p(18-1)*x(3-1),
				  p(19-1)*x(1-1)/(p(20-1) + x(1-1)) - p(21-1)*x(4-1),
			 	  p(22-1)*x(2-1)/(p(23-1) + x(2-1)) - p(24-1)*x(5-1),
				  p(25-1)*x(3-1)/(p(26-1) + x(3-1)) - p(27-1)*x(6-1),
				  p(28-1)*x(4-1)*(S - x(7-1))/(p(29-1)*(1 + S/p(29-1) + x(7-1)/p(30-1))) - p(31-1)*x(5-1)*(x(7-1) - x(8-1))/(p(32-1)*(1 + x(7-1)/p(32-1) + x(8-1)/p(33-1))),
				  p(31-1)*x(5-1)*(x(7-1) - x(8-1))/(p(32-1)*(1 + x(7-1)/p(32-1) + x(8-1)/p(33-1))) -  p(34-1)*x(6-1)*(x(8-1) - P)/(p(35-1)*(1 + x(8-1)/p(35-1) + P/p(36-1)));

		return result;
	}

	 mat nonlinearOdes::eight_part_spc(const double& t, const vec& x, const vec& pp)
	{
		double P = 1;
		double S = .1;

		vec p(36);
		p << 1.0,1.0,2.0,1.0,2.0,1.0,1.0,1.0,2.0,1.0,2.0,1.0,1.0,1.0,2.0,1.0,2.0,1.0,0.1,1.0,0.1,0.1,1.0,0.1,0.1,1.0,0.1,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0;

		p(5) = pp(0);
		p(10) = pp(1);
		p(11) = pp(2);
		for(int i=16; i<27; i++){
			p(i) = pp(i-13);
		}
		for(int i=30; i<34; i++){
			p(i) = pp(i-16);
		}
		mat result(8,1);
		result << p(1-1)/(1 + pow(P/p(2-1), p(3-1))   + pow(p(4-1)/S      , p(5-1)))  - p(6-1)*x(1-1),
				  p(7-1)/(1 + pow(P/p(8-1), p(9-1))   + pow(p(10-1)/x(7-1), p(11-1))) - p(12-1)*x(2-1),
				  p(13-1)/(1 + pow(P/p(14-1),p(15-1)) + pow(p(16-1)/x(8-1), p(17-1))) - p(18-1)*x(3-1),
				  p(19-1)*x(1-1)/(p(20-1) + x(1-1)) - p(21-1)*x(4-1),
			 	  p(22-1)*x(2-1)/(p(23-1) + x(2-1)) - p(24-1)*x(5-1),
				  p(25-1)*x(3-1)/(p(26-1) + x(3-1)) - p(27-1)*x(6-1),
				  p(28-1)*x(4-1)*(S - x(7-1))/(p(29-1)*(1 + S/p(29-1) + x(7-1)/p(30-1))) - p(31-1)*x(5-1)*(x(7-1) - x(8-1))/(p(32-1)*(1 + x(7-1)/p(32-1) + x(8-1)/p(33-1))),
				  p(31-1)*x(5-1)*(x(7-1) - x(8-1))/(p(32-1)*(1 + x(7-1)/p(32-1) + x(8-1)/p(33-1))) -  p(34-1)*x(6-1)*(x(8-1) - P)/(p(35-1)*(1 + x(8-1)/p(35-1) + P/p(36-1)));

		return result;
	} /**/

	/*mat nonlinearOdes::eight_part(const double& t, const vec& x, const vec& pp)
	{
		double P = 1;
		double S = .1;

		vec p(36);
		if(pp.size() < 36){
			p << 1.0,1.0,2.0,1.0,2.0,1.0,1.0,1.0,2.0,1.0,2.0,1.0,1.0,1.0,2.0,1.0,2.0,1.0,0.1,1.0,0.1,0.1,1.0,0.1,0.1,1.0,0.1,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0;
			for(int i=0; i<pp.size(); i++){
				//if(i != 1)
					p(i) = pp(i);
			}
		}else{
			p << pp;
		}
		mat result(8,1);
		result << p(1-1)/(1 + pow(P/p(2-1), p(3-1)) + pow(p(4-1)/S      , p(5-1)))  - p(6-1)*x(1-1),
				  p(7-1)/(1 + pow(P/p(8-1), p(9-1)) + pow(p(10-1)/x(7-1), p(11-1))) - p(12-1)*x(2-1),
				 p(13-1)/(1 + pow(P/p(14-1),p(15-1)) + pow(p(16-1)/x(8-1), p(17-1))) - p(18-1)*x(3-1),
				  p(19-1)*x(1-1)/(p(20-1) + x(1-1)) - p(21-1)*x(4-1),
			 	  p(22-1)*x(2-1)/(p(23-1) + x(2-1)) - p(24-1)*x(5-1),
				  p(25-1)*x(3-1)/(p(26-1) + x(3-1)) - p(27-1)*x(6-1),
				  p(28-1)*x(4-1)*(S - x(7-1))/(p(29-1)*(1 + S/p(29-1) + x(7-1)/p(30-1))) - p(31-1)*x(5-1)*(x(7-1) - x(8-1))/(p(32-1)*(1 + x(7-1)/p(32-1) + x(8-1)/p(33-1))),
				  p(31-1)*x(5-1)*(x(7-1) - x(8-1))/(p(32-1)*(1 + x(7-1)/p(32-1) + x(8-1)/p(33-1))) - p(34-1)*x(6-1)*(x(8-1) - P)/(p(35-1)*(1 + x(8-1)/p(35-1) + P/p(36-1)));

		return result;
	}*/

	/*mat nonlinearOdes::eight_part_linearization(const double& t, const vec& x, const vec& p, const mat& xn, const vec& time)
	{
		double P = 1;
		double S = 1;

		thesis::spline Xn1(time, xn.row(0));
		thesis::spline Xn2(time, xn.row(1));
		thesis::spline Xn3(time, xn.row(2));
		thesis::spline Xn4(time, xn.row(3));
		thesis::spline Xn5(time, xn.row(4));
		thesis::spline Xn6(time, xn.row(5));
		thesis::spline Xn7(time, xn.row(6));
		thesis::spline Xn8(time, xn.row(7));

		double x1 = Xn1.interpolate(t);
		double x2 = Xn2.interpolate(t);
		double x3 = Xn3.interpolate(t);
		double x4 = Xn4.interpolate(t);
		double x5 = Xn5.interpolate(t);
		double x6 = Xn6.interpolate(t);
		double x7 = Xn7.interpolate(t);
		double x8 = Xn8.interpolate(t);

		mat result(296,1);
		result << p(1-1)/(1+pow(P/p(2-1), p(3-1))+pow(p(4-1)/S, p(5-1)))-p(6-1)*x(1-1),
		-(p(12-1)*x(2-1))+(p(7-1)*((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))-p(11-1)*pow(p(10-1)/x7,p(11-1)))*x7+p(11-1)*pow(p(10-1)/x7,p(11-1))*x(7-1)))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(3-1))+(p(13-1)*((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))-p(17-1)*pow(p(16-1)/x8,p(17-1)))*x8+p(17-1)*pow(p(16-1)/x8,p(17-1))*x(8-1)))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*pow(x1,2)+p(19-1)*p(20-1)*x(1-1)-p(21-1)*pow(p(20-1)+x1,2)*x(4-1))/pow(p(20-1)+x1,2),
		(p(22-1)*pow(x2,2)+p(22-1)*p(23-1)*x(2-1)-p(24-1)*pow(p(23-1)+x2,2)*x(5-1))/pow(p(23-1)+x2,2),
		(p(25-1)*pow(x3,2)+p(25-1)*p(26-1)*x(3-1)-p(27-1)*pow(p(26-1)+x3,2)*x(6-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*x4*(S-x7))/(p(30-1)*S+p(29-1)*(p(30-1)+x7))-(p(31-1)*p(33-1)*x5*(x7-x8))/(p(33-1)*x7+p(32-1)*(p(33-1)+x8))-(p(28-1)*p(30-1)*(S-x7)*(x4-x(4-1)))/(p(30-1)*S+p(29-1)*(p(30-1)+x7))+(p(31-1)*p(33-1)*(x7-x8)*(x5-x(5-1)))/(p(33-1)*x7+p(32-1)*(p(33-1)+x8))+(-((p(28-1)*p(30-1)*(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4)/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2))-(p(31-1)*p(33-1)*x5*(p(33-1)*x8+p(32-1)*(p(33-1)+x8)))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2))*(-x7+x(7-1))-(p(31-1)*p(33-1)*x5*(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*(x8-x(8-1)))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(31-1)*p(33-1)*x5*(x7-x8))/(p(33-1)*x7+p(32-1)*(p(33-1)+x8))+(p(34-1)*p(36-1)*x6*(P-x8))/(P*p(35-1)+p(36-1)*(p(35-1)+x8))-(p(31-1)*p(33-1)*(x7-x8)*(x5-x(5-1)))/(p(33-1)*x7+p(32-1)*(p(33-1)+x8))-(p(34-1)*p(36-1)*(P-x8)*(x6-x(6-1)))/(P*p(35-1)+p(36-1)*(p(35-1)+x8))-(p(31-1)*p(33-1)*x5*(p(33-1)*x8+p(32-1)*(p(33-1)+x8))*(x7-x(7-1)))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2)+(-((p(31-1)*p(33-1)*x5*(p(33-1)*x7+p(32-1)*(p(33-1)+x7)))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2))-(p(34-1)*p(36-1)*(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6)/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2))*(-x8+x(8-1)),
		(1/(1+pow(P/p(2-1), p(3-1))+pow(p(4-1)/S, p(5-1))))-p(6-1)*x(9-1),
		-(p(12-1)*x(10-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(15-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(11-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(16-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(9-1)-p(21-1)*pow(p(20-1)+x1,2)*x(12-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(10-1)-p(24-1)*pow(p(23-1)+x2,2)*x(13-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(11-1)-p(27-1)*pow(p(26-1)+x3,2)*x(14-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(12-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(15-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(13-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(15-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(16-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(14-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(16-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(13-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(15-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(16-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(1-1)*pow(P/p(2-1), p(3-1))*p(3-1))/(p(2-1)*pow((1+pow(P/p(2-1), p(3-1))+pow(p(4-1)/S, p(5-1))),2))-p(6-1)*x(17-1),
		-(p(12-1)*x(18-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(23-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(19-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(24-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(17-1)-p(21-1)*pow(p(20-1)+x1,2)*x(20-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(18-1)-p(24-1)*pow(p(23-1)+x2,2)*x(21-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(19-1)-p(27-1)*pow(p(26-1)+x3,2)*x(22-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(20-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(23-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(21-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(23-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(24-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(22-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(24-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(21-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(23-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(24-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-((p(1-1)*pow(P/p(2-1), p(3-1))*log(P/p(2-1)))/pow((1+pow(P/p(2-1), p(3-1))+pow(p(4-1)/S, p(5-1))),2))-p(6-1)*x(25-1),
		-(p(12-1)*x(26-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(31-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(27-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(32-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(25-1)-p(21-1)*pow(p(20-1)+x1,2)*x(28-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(26-1)-p(24-1)*pow(p(23-1)+x2,2)*x(29-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(27-1)-p(27-1)*pow(p(26-1)+x3,2)*x(30-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(28-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(31-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(29-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(31-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(32-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(30-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(32-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(29-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(31-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(32-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-((p(1-1)*p(5-1)*pow(p(4-1)/S, p(5-1)))/(p(4-1)*pow((1+pow(P/p(2-1), p(3-1))+pow(p(4-1)/S, p(5-1))),2)))-p(6-1)*x(33-1),
		-(p(12-1)*x(34-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(39-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(35-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(40-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(33-1)-p(21-1)*pow(p(20-1)+x1,2)*x(36-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(34-1)-p(24-1)*pow(p(23-1)+x2,2)*x(37-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(35-1)-p(27-1)*pow(p(26-1)+x3,2)*x(38-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(36-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(39-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(37-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(39-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(40-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(38-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(40-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(37-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(39-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(40-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-((p(1-1)*pow(p(4-1)/S, p(5-1))*log(p(4-1)/S))/pow((1+pow(P/p(2-1), p(3-1))+pow(p(4-1)/S, p(5-1))),2))-p(6-1)*x(41-1),
		-(p(12-1)*x(42-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(47-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(43-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(48-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(41-1)-p(21-1)*pow(p(20-1)+x1,2)*x(44-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(42-1)-p(24-1)*pow(p(23-1)+x2,2)*x(45-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(43-1)-p(27-1)*pow(p(26-1)+x3,2)*x(46-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(44-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(47-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(45-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(47-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(48-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(46-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(48-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(45-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(47-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(48-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-x1-p(6-1)*x(49-1),
		-(p(12-1)*x(50-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(55-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(51-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(56-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(49-1)-p(21-1)*pow(p(20-1)+x1,2)*x(52-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(50-1)-p(24-1)*pow(p(23-1)+x2,2)*x(53-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(51-1)-p(27-1)*pow(p(26-1)+x3,2)*x(54-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(52-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(55-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(53-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(55-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(56-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(54-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(56-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(53-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(55-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(56-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(57-1)),
		(p(11-1)*pow(p(10-1)/x7,p(11-1))*x(7-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7)-p(12-1)*x(58-1)+((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))-p(11-1)*pow(p(10-1)/x7,p(11-1)))*x7+p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(63-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(59-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(64-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(57-1)-p(21-1)*pow(p(20-1)+x1,2)*x(60-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(58-1)-p(24-1)*pow(p(23-1)+x2,2)*x(61-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(59-1)-p(27-1)*pow(p(26-1)+x3,2)*x(62-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(60-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(63-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(61-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(63-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(64-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(62-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(64-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(61-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(63-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(64-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(65-1)),
		(2*p(7-1)*pow(P/p(8-1), p(9-1))*p(9-1)*((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))-p(11-1)*pow(p(10-1)/x7,p(11-1)))*x7+p(11-1)*pow(p(10-1)/x7,p(11-1))*x(7-1)))/(p(8-1)*pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),3)*x7)-p(12-1)*x(66-1)+(p(7-1)*(-((pow(P/p(8-1), p(9-1))*p(9-1)*x7)/p(8-1))+p(11-1)*pow(p(10-1)/x7,p(11-1))*x(71-1)))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(67-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(72-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(65-1)-p(21-1)*pow(p(20-1)+x1,2)*x(68-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(66-1)-p(24-1)*pow(p(23-1)+x2,2)*x(69-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(67-1)-p(27-1)*pow(p(26-1)+x3,2)*x(70-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(68-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(71-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(69-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(71-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(72-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(70-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(72-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(69-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(71-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(72-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(73-1)),
		(-2*p(7-1)*pow(P/p(8-1), p(9-1))*log(P/p(8-1))*((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))-p(11-1)*pow(p(10-1)/x7,p(11-1)))*x7+p(11-1)*pow(p(10-1)/x7,p(11-1))*x(7-1)))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),3)*x7)-p(12-1)*x(74-1)+(p(7-1)*(pow(P/p(8-1), p(9-1))*x7*log(P/p(8-1))+p(11-1)*pow(p(10-1)/x7,p(11-1))*x(79-1)))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(75-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(80-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(73-1)-p(21-1)*pow(p(20-1)+x1,2)*x(76-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(74-1)-p(24-1)*pow(p(23-1)+x2,2)*x(77-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(75-1)-p(27-1)*pow(p(26-1)+x3,2)*x(78-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(76-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(79-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(77-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(79-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(80-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(78-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(80-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(77-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(79-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(80-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(81-1)),
		(pow(p(11-1),2)*p(7-1)*(1+pow(P/p(8-1), p(9-1))-pow(p(10-1)/x7,p(11-1)))*pow(p(10-1)/x7,-1+p(11-1))*x(7-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),3)*pow(x7,2))-p(12-1)*x(82-1)+(p(11-1)*p(7-1)*pow(p(10-1)/x7,-1+p(11-1))*(-((1+pow(P/p(8-1), p(9-1))+p(11-1)*(1+pow(P/p(8-1), p(9-1))-pow(p(10-1)/x7,p(11-1)))+pow(p(10-1)/x7,p(11-1)))*x7)+p(10-1)*(1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1)))*x(87-1)))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),3)*pow(x7,2)),
		-(p(18-1)*x(83-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(88-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(81-1)-p(21-1)*pow(p(20-1)+x1,2)*x(84-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(82-1)-p(24-1)*pow(p(23-1)+x2,2)*x(85-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(83-1)-p(27-1)*pow(p(26-1)+x3,2)*x(86-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(84-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(87-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(85-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(87-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(88-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(86-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(88-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(85-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(87-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(88-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(89-1)),
		(-2*p(7-1)*pow(p(10-1)/x7,p(11-1))*log(p(10-1)/x7)*((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))-p(11-1)*pow(p(10-1)/x7,p(11-1)))*x7+p(11-1)*pow(p(10-1)/x7,p(11-1))*x(7-1)))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),3)*x7)-p(12-1)*x(90-1)-(p(7-1)*pow(p(10-1)/x7,p(11-1))*(x7-x7*log(p(10-1)/x7)+p(11-1)*x7*log(p(10-1)/x7)-(1+p(11-1)*log(p(10-1)/x7))*x(7-1)-p(11-1)*x(95-1)))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(91-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(96-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(89-1)-p(21-1)*pow(p(20-1)+x1,2)*x(92-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(90-1)-p(24-1)*pow(p(23-1)+x2,2)*x(93-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(91-1)-p(27-1)*pow(p(26-1)+x3,2)*x(94-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(92-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(95-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(93-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(95-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(96-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(94-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(96-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(93-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(95-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(96-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(97-1)),
		-x2-p(12-1)*x(98-1)+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(103-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(99-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(104-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(97-1)-p(21-1)*pow(p(20-1)+x1,2)*x(100-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(98-1)-p(24-1)*pow(p(23-1)+x2,2)*x(101-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(99-1)-p(27-1)*pow(p(26-1)+x3,2)*x(102-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(100-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(103-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(101-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(103-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(104-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(102-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(104-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(101-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(103-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(104-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(105-1)),
		-(p(12-1)*x(106-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(111-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		(p(17-1)*pow(p(16-1)/x8,p(17-1))*x(8-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8)-p(18-1)*x(107-1)+((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))-p(17-1)*pow(p(16-1)/x8,p(17-1)))*x8+p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(112-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(105-1)-p(21-1)*pow(p(20-1)+x1,2)*x(108-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(106-1)-p(24-1)*pow(p(23-1)+x2,2)*x(109-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(107-1)-p(27-1)*pow(p(26-1)+x3,2)*x(110-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(108-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(111-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(109-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(111-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(112-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(110-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(112-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(109-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(111-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(112-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(113-1)),
		-(p(12-1)*x(114-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(119-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		(2*p(13-1)*pow(P/p(14-1), p(15-1))*p(15-1)*((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))-p(17-1)*pow(p(16-1)/x8,p(17-1)))*x8+p(17-1)*pow(p(16-1)/x8,p(17-1))*x(8-1)))/(p(14-1)*pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),3)*x8)-p(18-1)*x(115-1)+(p(13-1)*(-((pow(P/p(14-1), p(15-1))*p(15-1)*x8)/p(14-1))+p(17-1)*pow(p(16-1)/x8,p(17-1))*x(120-1)))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(113-1)-p(21-1)*pow(p(20-1)+x1,2)*x(116-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(114-1)-p(24-1)*pow(p(23-1)+x2,2)*x(117-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(115-1)-p(27-1)*pow(p(26-1)+x3,2)*x(118-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(116-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(119-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(117-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(119-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(120-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(118-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(120-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(117-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(119-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(120-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(121-1)),
		-(p(12-1)*x(122-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(127-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		(-2*p(13-1)*pow(P/p(14-1), p(15-1))*log(P/p(14-1))*((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))-p(17-1)*pow(p(16-1)/x8,p(17-1)))*x8+p(17-1)*pow(p(16-1)/x8,p(17-1))*x(8-1)))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),3)*x8)-p(18-1)*x(123-1)+(p(13-1)*(pow(P/p(14-1), p(15-1))*x8*log(P/p(14-1))+p(17-1)*pow(p(16-1)/x8,p(17-1))*x(128-1)))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(121-1)-p(21-1)*pow(p(20-1)+x1,2)*x(124-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(122-1)-p(24-1)*pow(p(23-1)+x2,2)*x(125-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(123-1)-p(27-1)*pow(p(26-1)+x3,2)*x(126-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(124-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(127-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(125-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(127-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(128-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(126-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(128-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(125-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(127-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(128-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(129-1)),
		-(p(12-1)*x(130-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(135-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		(p(13-1)*pow(p(17-1),2)*(1+pow(P/p(14-1), p(15-1))-pow(p(16-1)/x8,p(17-1)))*pow(p(16-1)/x8,p(17-1)-1)*x(8-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),3)*pow(x8,2))-p(18-1)*x(131-1)+(p(13-1)*p(17-1)*pow(p(16-1)/x8,-1+p(17-1))*(-((1+pow(P/p(14-1), p(15-1))+p(17-1)*(1+pow(P/p(14-1), p(15-1))-pow(p(16-1)/x8,p(17-1)))+pow(p(16-1)/x8,p(17-1)))*x8)+p(16-1)*(1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1)))*x(136-1)))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),3)*pow(x8,2)),
		(p(19-1)*p(20-1)*x(129-1)-p(21-1)*pow(p(20-1)+x1,2)*x(132-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(130-1)-p(24-1)*pow(p(23-1)+x2,2)*x(133-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(131-1)-p(27-1)*pow(p(26-1)+x3,2)*x(134-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(132-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(135-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(133-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(135-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(136-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(134-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(136-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(133-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(135-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(136-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(137-1)),
		-(p(12-1)*x(138-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(143-1))/(pow(1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1)),2)*x7),
		(-2*p(13-1)*pow(p(16-1)/x8,p(17-1))*log(p(16-1)/x8)*((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))-p(17-1)*pow(p(16-1)/x8,p(17-1)))*x8+p(17-1)*pow(p(16-1)/x8,p(17-1))*x(8-1)))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),3)*x8)-p(18-1)*x(139-1)-(p(13-1)*pow(p(16-1)/x8,p(17-1))*(x8-x8*log(p(16-1)/x8)+p(17-1)*x8*log(p(16-1)/x8)-(1+p(17-1)*log(p(16-1)/x8))*x(8-1)-p(17-1)*x(144-1)))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(137-1)-p(21-1)*pow(p(20-1)+x1,2)*x(140-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(138-1)-p(24-1)*pow(p(23-1)+x2,2)*x(141-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(139-1)-p(27-1)*pow(p(26-1)+x3,2)*x(142-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(140-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(143-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(141-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(143-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(144-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(142-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(144-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(141-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(143-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(144-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(145-1)),
		-(p(12-1)*x(146-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(151-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-x(3-1)-p(18-1)*x(147-1)+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(152-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(145-1)-p(21-1)*pow(p(20-1)+x1,2)*x(148-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(146-1)-p(24-1)*pow(p(23-1)+x2,2)*x(149-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(147-1)-p(27-1)*pow(p(26-1)+x3,2)*x(150-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(148-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(151-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(149-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(151-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(152-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(150-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(152-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(149-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(151-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(152-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(153-1)),
		-(p(12-1)*x(154-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(159-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(155-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(160-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(pow(x1,2)+p(20-1)*x1+p(19-1)*p(20-1)*x(153-1)-p(21-1)*pow(p(20-1)+x1,2)*x(156-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(154-1)-p(24-1)*pow(p(23-1)+x2,2)*x(157-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(155-1)-p(27-1)*pow(p(26-1)+x3,2)*x(158-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(156-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(159-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(157-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(159-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(160-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(158-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(160-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(157-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(159-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(160-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(161-1)),
		-(p(12-1)*x(162-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(167-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(163-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(168-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(-2*(p(19-1)*pow(x1,2)+p(19-1)*p(20-1)*x1-p(21-1)*pow(p(20-1)+x1,2)*x(4-1))+(p(20-1)+x1)*(-2*p(21-1)*(p(20-1)+x1)*x(4-1)+p(19-1)*(x1+p(20-1)*x(161-1))-p(21-1)*pow(p(20-1)+x1,2)*x(164-1)))/pow(p(20-1)+x1,3),
		(p(22-1)*p(23-1)*x(162-1)-p(24-1)*pow(p(23-1)+x2,2)*x(165-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(163-1)-p(27-1)*pow(p(26-1)+x3,2)*x(166-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(164-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(167-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(165-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(167-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(168-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(166-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(168-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(165-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(167-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(168-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(169-1)),
		-(p(12-1)*x(170-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(175-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(171-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(176-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(-(pow(p(20-1)+x1,2)*x(4-1))+p(19-1)*p(20-1)*x(169-1)-p(21-1)*pow(p(20-1)+x1,2)*x(172-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(170-1)-p(24-1)*pow(p(23-1)+x2,2)*x(173-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(171-1)-p(27-1)*pow(p(26-1)+x3,2)*x(174-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(172-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(175-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(173-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(175-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(176-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(174-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(176-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(173-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(175-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(176-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(177-1)),
		-(p(12-1)*x(178-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(183-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(179-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(184-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(177-1)-p(21-1)*pow(p(20-1)+x1,2)*x(180-1))/pow(p(20-1)+x1,2),
		(pow(x2,2)+p(23-1)*x2+p(22-1)*p(23-1)*x(178-1)-p(24-1)*pow(p(23-1)+x2,2)*x(181-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(179-1)-p(27-1)*pow(p(26-1)+x3,2)*x(182-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(180-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(183-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(181-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(183-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(184-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(182-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(184-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(181-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(183-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(184-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(185-1)),
		-(p(12-1)*x(186-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(191-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(187-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(192-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(185-1)-p(21-1)*pow(p(20-1)+x1,2)*x(188-1))/pow(p(20-1)+x1,2),
		(-2*(p(22-1)*pow(x2,2)+p(22-1)*p(23-1)*x2-p(24-1)*pow(p(23-1)+x2,2)*x(5-1))+(p(23-1)+x2)*(-2*p(24-1)*(p(23-1)+x2)*x(5-1)+p(22-1)*(x2+p(23-1)*x(186-1))-p(24-1)*pow(p(23-1)+x2,2)*x(189-1)))/pow(p(23-1)+x2,3),
		(p(25-1)*p(26-1)*x(187-1)-p(27-1)*pow(p(26-1)+x3,2)*x(190-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(188-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(191-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(189-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(191-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(192-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(190-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(192-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(189-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(191-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(192-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(193-1)),
		-(p(12-1)*x(194-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(199-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(195-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(200-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(193-1)-p(21-1)*pow(p(20-1)+x1,2)*x(196-1))/pow(p(20-1)+x1,2),
		(-(pow(p(23-1)+x2,2)*x(5-1))+p(22-1)*p(23-1)*x(194-1)-p(24-1)*pow(p(23-1)+x2,2)*x(197-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(195-1)-p(27-1)*pow(p(26-1)+x3,2)*x(198-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(196-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(199-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(197-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(199-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(200-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(198-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(200-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(197-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(199-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(200-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(201-1)),
		-(p(12-1)*x(202-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(207-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(203-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(208-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(201-1)-p(21-1)*pow(p(20-1)+x1,2)*x(204-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(202-1)-p(24-1)*pow(p(23-1)+x2,2)*x(205-1))/pow(p(23-1)+x2,2),
		(pow(x3,2)+p(26-1)*x(3-1)+p(25-1)*p(26-1)*x(203-1)-p(27-1)*pow(p(26-1)+x3,2)*x(206-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(204-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(207-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(205-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(207-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(208-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(206-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(208-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(205-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(207-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(208-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(209-1)),
		-(p(12-1)*x(210-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(215-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(211-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(216-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(209-1)-p(21-1)*pow(p(20-1)+x1,2)*x(212-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(210-1)-p(24-1)*pow(p(23-1)+x2,2)*x(213-1))/pow(p(23-1)+x2,2),
		(-2*(p(25-1)*pow(x3,2)+p(25-1)*p(26-1)*x(3-1)-p(27-1)*pow(p(26-1)+x3,2)*x(6-1))+(p(26-1)+x3)*(-2*p(27-1)*(p(26-1)+x3)*x(6-1)+p(25-1)*(x(3-1)+p(26-1)*x(211-1))-p(27-1)*pow(p(26-1)+x3,2)*x(214-1)))/pow(p(26-1)+x3,3),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(212-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(215-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(213-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(215-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(216-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(214-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(216-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(213-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(215-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(216-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(217-1)),
		-(p(12-1)*x(218-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(223-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(219-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(224-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(217-1)-p(21-1)*pow(p(20-1)+x1,2)*x(220-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(218-1)-p(24-1)*pow(p(23-1)+x2,2)*x(221-1))/pow(p(23-1)+x2,2),
		(-(pow(p(26-1)+x3,2)*x(6-1))+p(25-1)*p(26-1)*x(219-1)-p(27-1)*pow(p(26-1)+x3,2)*x(222-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(220-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(223-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(221-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(223-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(224-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(222-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(224-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(221-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(223-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(224-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(225-1)),
		-(p(12-1)*x(226-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(231-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(227-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(232-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(225-1)-p(21-1)*pow(p(20-1)+x1,2)*x(228-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(226-1)-p(24-1)*pow(p(23-1)+x2,2)*x(229-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(227-1)-p(27-1)*pow(p(26-1)+x3,2)*x(230-1))/pow(p(26-1)+x3,2),
		(p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(4-1)+(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*(x7-x(7-1))))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)+(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(228-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(231-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(229-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(231-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(232-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(230-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(232-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(229-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(231-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(232-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(233-1)),
		-(p(12-1)*x(234-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(239-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(235-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(240-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(233-1)-p(21-1)*pow(p(20-1)+x1,2)*x(236-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(234-1)-p(24-1)*pow(p(23-1)+x2,2)*x(237-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(235-1)-p(27-1)*pow(p(26-1)+x3,2)*x(238-1))/pow(p(26-1)+x3,2),
		(2*p(28-1)*p(30-1)*(p(30-1)+x7)*(-((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(4-1))-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*(x7-x(7-1))))/pow(p(30-1)*S+p(29-1)*(p(30-1)+x7),3)+(p(28-1)*p(30-1)*(p(30-1)*x4*x7+S*x4*x7+(S-x7)*(p(30-1)+x7)*x(4-1)-(p(30-1)+S)*x4*x(7-1)+p(29-1)*p(30-1)*S*x(236-1)+p(30-1)*pow(S,2)*x(236-1)-p(29-1)*p(30-1)*x7*x(236-1)+p(29-1)*S*x7*x(236-1)-p(30-1)*S*x7*x(236-1)-p(29-1)*pow(x7,2)*x(236-1)-p(29-1)*p(30-1)*x4*x(239-1)-p(29-1)*S*x4*x(239-1)-p(30-1)*S*x4*x(239-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(237-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(239-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(240-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(238-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(240-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(237-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(239-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(240-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(241-1)),
		-(p(12-1)*x(242-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(247-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(243-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(248-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(241-1)-p(21-1)*pow(p(20-1)+x1,2)*x(244-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(242-1)-p(24-1)*pow(p(23-1)+x2,2)*x(245-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(243-1)-p(27-1)*pow(p(26-1)+x3,2)*x(246-1))/pow(p(26-1)+x3,2),
		(-2*p(28-1)*p(30-1)*(p(29-1)+S)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(4-1)+(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*(x7-x(7-1))))/pow(p(30-1)*S+p(29-1)*(p(30-1)+x7),3)+(p(28-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(4-1)+(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*(x7-x(7-1))))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(28-1)*p(30-1)*(-(p(29-1)*x4*x7)-S*x4*x7-(p(29-1)+S)*(S-x7)*x(4-1)+(p(29-1)+S)*x4*x(7-1)-p(29-1)*p(30-1)*S*x(244-1)-p(30-1)*pow(S,2)*x(244-1)+p(29-1)*p(30-1)*x7*x(244-1)-p(29-1)*S*x7*x(244-1)+p(30-1)*S*x7*x(244-1)+p(29-1)*pow(x7,2)*x(244-1)+p(29-1)*p(30-1)*x4*x(247-1)+p(29-1)*S*x4*x(247-1)+p(30-1)*S*x4*x(247-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(245-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(247-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(248-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(246-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(248-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(245-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(247-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(248-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(249-1)),
		-(p(12-1)*x(250-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(255-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(251-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(256-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(249-1)-p(21-1)*pow(p(20-1)+x1,2)*x(252-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(250-1)-p(24-1)*pow(p(23-1)+x2,2)*x(253-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(251-1)-p(27-1)*pow(p(26-1)+x3,2)*x(254-1))/pow(p(26-1)+x3,2),
		(p(33-1)*(-((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(5-1))+x5*(p(32-1)*p(33-1)*(x7-x8)-(p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(7-1)+(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(8-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2)+(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(252-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(255-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(253-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(255-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(256-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(5-1)+x5*(p(32-1)*p(33-1)*(-x7+x8)+(p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(7-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(8-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2)+(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(254-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(256-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(253-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(255-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(256-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(257-1)),
		-(p(12-1)*x(258-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(263-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(259-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(264-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(257-1)-p(21-1)*pow(p(20-1)+x1,2)*x(260-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(258-1)-p(24-1)*pow(p(23-1)+x2,2)*x(261-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(259-1)-p(27-1)*pow(p(26-1)+x3,2)*x(262-1))/pow(p(26-1)+x3,2),
		(2*p(31-1)*p(33-1)*(p(33-1)+x8)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(5-1)+x5*(p(32-1)*p(33-1)*(-x7+x8)+(p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(7-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(8-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),3)+(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(260-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(263-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)+(p(31-1)*p(33-1)*(p(33-1)*x5*x7-p(33-1)*x5*x8-(x7-x8)*(p(33-1)+x8)*x(5-1)-x5*(p(33-1)+x8)*x(7-1)+p(33-1)*x5*x(8-1)+x5*x7*x(8-1)-p(32-1)*p(33-1)*x7*x(261-1)-p(33-1)*pow(x7,2)*x(261-1)+p(32-1)*p(33-1)*x8*x(261-1)-p(32-1)*x7*x8*x(261-1)+p(33-1)*x7*x8*x(261-1)+p(32-1)*pow(x8,2)*x(261-1)-p(32-1)*p(33-1)*x5*x(263-1)-p(32-1)*x5*x8*x(263-1)-p(33-1)*x5*x8*x(263-1)+p(32-1)*p(33-1)*x5*x(264-1)+p(32-1)*x5*x7*x(264-1)+p(33-1)*x5*x7*x(264-1)))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(2*p(31-1)*p(33-1)*(p(33-1)+x8)*(-((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(5-1))+x5*(p(32-1)*p(33-1)*(x7-x8)-(p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(7-1)+(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(8-1))))/pow((p(33-1)*x7+p(32-1)*(p(33-1)+x8)),3)+(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(262-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(264-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*(-(p(33-1)*x5*x7)+p(33-1)*x5*x8+(x7-x8)*(p(33-1)+x8)*x(5-1)+x5*(p(33-1)+x8)*x(7-1)-p(33-1)*x5*x(8-1)-x5*x7*x(8-1)+p(32-1)*p(33-1)*x7*x(261-1)+p(33-1)*pow(x7,2)*x(261-1)-p(32-1)*p(33-1)*x8*x(261-1)+p(32-1)*x7*x8*x(261-1)-p(33-1)*x7*x8*x(261-1)-p(32-1)*pow(x8,2)*x(261-1)+p(32-1)*p(33-1)*x5*x(263-1)+p(32-1)*x5*x8*x(263-1)+p(33-1)*x5*x8*x(263-1)-p(32-1)*p(33-1)*x5*x(264-1)-p(32-1)*x5*x7*x(264-1)-p(33-1)*x5*x7*x(264-1)))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(265-1)),
		-(p(12-1)*x(266-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(271-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(267-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(272-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(265-1)-p(21-1)*pow(p(20-1)+x1,2)*x(268-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(266-1)-p(24-1)*pow(p(23-1)+x2,2)*x(269-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(267-1)-p(27-1)*pow(p(26-1)+x3,2)*x(270-1))/pow(p(26-1)+x3,2),
		(2*p(31-1)*p(33-1)*(p(32-1)+x7)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(5-1)+x5*(p(32-1)*p(33-1)*(-x7+x8)+(p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(7-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(8-1))))/pow((p(33-1)*x7+p(32-1)*(p(33-1)+x8)),3)+(p(31-1)*(-((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(5-1))+x5*(p(32-1)*p(33-1)*(x7-x8)-(p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(7-1)+(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(8-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2)+(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(268-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(271-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*(-(p(32-1)*x5*x7)+p(32-1)*x5*x8+(p(32-1)+x7)*(x7-x8)*x(5-1)+x5*(p(32-1)+x8)*x(7-1)-p(32-1)*x5*x(8-1)-x5*x7*x(8-1)+p(32-1)*p(33-1)*x7*x(269-1)+p(33-1)*pow(x7,2)*x(269-1)-p(32-1)*p(33-1)*x8*x(269-1)+p(32-1)*x7*x8*x(269-1)-p(33-1)*x7*x8*x(269-1)-p(32-1)*pow(x8,2)*x(269-1)+p(32-1)*p(33-1)*x5*x(271-1)+p(32-1)*x5*x8*x(271-1)+p(33-1)*x5*x8*x(271-1)-p(32-1)*p(33-1)*x5*x(272-1)-p(32-1)*x5*x7*x(272-1)-p(33-1)*x5*x7*x(272-1)))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(31-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(5-1)+x5*(p(32-1)*p(33-1)*(-x7+x8)+(p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(7-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(8-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2)+(2*p(31-1)*p(33-1)*(p(32-1)+x7)*(-((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(5-1))+x5*(p(32-1)*p(33-1)*(x7-x8)-(p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(7-1)+(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(8-1))))/pow((p(33-1)*x7+p(32-1)*(p(33-1)+x8)),3)+(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(270-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(272-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*(-(p(32-1)*x5*x7)+p(32-1)*x5*x8+(p(32-1)+x7)*(x7-x8)*x(5-1)+x5*(p(32-1)+x8)*x(7-1)-p(32-1)*x5*x(8-1)-x5*x7*x(8-1)+p(32-1)*p(33-1)*x7*x(269-1)+p(33-1)*pow(x7,2)*x(269-1)-p(32-1)*p(33-1)*x8*x(269-1)+p(32-1)*x7*x8*x(269-1)-p(33-1)*x7*x8*x(269-1)-p(32-1)*pow(x8,2)*x(269-1)+p(32-1)*p(33-1)*x5*x(271-1)+p(32-1)*x5*x8*x(271-1)+p(33-1)*x5*x8*x(271-1)-p(32-1)*p(33-1)*x5*x(272-1)-p(32-1)*x5*x7*x(272-1)-p(33-1)*x5*x7*x(272-1)))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(273-1)),
		-(p(12-1)*x(274-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(279-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(275-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(280-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(273-1)-p(21-1)*pow(p(20-1)+x1,2)*x(276-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(274-1)-p(24-1)*pow(p(23-1)+x2,2)*x(277-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(275-1)-p(27-1)*pow(p(26-1)+x3,2)*x(278-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(276-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(279-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(277-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(279-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(280-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(6-1)+(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*(x8-x(8-1))))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(34-1)*p(36-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(278-1)-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*x(280-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(277-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(279-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(280-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(281-1)),
		-(p(12-1)*x(282-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(287-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(283-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(288-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(281-1)-p(21-1)*pow(p(20-1)+x1,2)*x(284-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(282-1)-p(24-1)*pow(p(23-1)+x2,2)*x(285-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(283-1)-p(27-1)*pow(p(26-1)+x3,2)*x(286-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(284-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(287-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(285-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(287-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(288-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(-2*p(34-1)*p(36-1)*(P+p(36-1))*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(6-1)+(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*(x8-x(8-1))))/pow((P*p(35-1)+p(36-1)*(p(35-1)+x8)),3)-(p(34-1)*p(36-1)*(-(P*x6*x8)-p(36-1)*x6*x8-(P+p(36-1))*(P-x8)*x(6-1)+(P+p(36-1))*x6*x(8-1)-pow(P,2)*p(35-1)*x(286-1)-P*p(35-1)*p(36-1)*x(286-1)+P*p(35-1)*x8*x(286-1)-P*p(36-1)*x8*x(286-1)+p(35-1)*p(36-1)*x8*x(286-1)+p(36-1)*pow(x8,2)*x(286-1)+P*p(35-1)*x6*x(288-1)+P*p(36-1)*x6*x(288-1)+p(35-1)*p(36-1)*x6*x(288-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(285-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(287-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(288-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		-(p(6-1)*x(289-1)),
		-(p(12-1)*x(290-1))+(p(11-1)*p(7-1)*pow(p(10-1)/x7,p(11-1))*x(295-1))/(pow((1+pow(P/p(8-1), p(9-1))+pow(p(10-1)/x7,p(11-1))),2)*x7),
		-(p(18-1)*x(291-1))+(p(13-1)*p(17-1)*pow(p(16-1)/x8,p(17-1))*x(296-1))/(pow((1+pow(P/p(14-1), p(15-1))+pow(p(16-1)/x8,p(17-1))),2)*x8),
		(p(19-1)*p(20-1)*x(289-1)-p(21-1)*pow(p(20-1)+x1,2)*x(292-1))/pow(p(20-1)+x1,2),
		(p(22-1)*p(23-1)*x(290-1)-p(24-1)*pow(p(23-1)+x2,2)*x(293-1))/pow(p(23-1)+x2,2),
		(p(25-1)*p(26-1)*x(291-1)-p(27-1)*pow(p(26-1)+x3,2)*x(294-1))/pow(p(26-1)+x3,2),
		(p(28-1)*p(30-1)*((S-x7)*(p(30-1)*S+p(29-1)*(p(30-1)+x7))*x(292-1)-(p(30-1)*S+p(29-1)*(p(30-1)+S))*x4*x(295-1)))/pow((p(30-1)*S+p(29-1)*(p(30-1)+x7)),2)-(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(293-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(295-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(296-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2),
		(2*p(34-1)*p(36-1)*(p(35-1)+x8)*(-((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(6-1))-(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*(x8-x(8-1))))/pow((P*p(35-1)+p(36-1)*(p(35-1)+x8)),3)+(p(34-1)*((P-x8)*(P*p(35-1)+p(36-1)*(p(35-1)+x8))*x(6-1)+(p(35-1)*p(36-1)+P*(p(35-1)+p(36-1)))*x6*(x8-x(8-1))))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)-(p(34-1)*p(36-1)*(-(P*x6*x8)-p(35-1)*x6*x8-(P-x8)*(p(35-1)+x8)*x(6-1)+(P+p(35-1))*x6*x(8-1)-pow(P,2)*p(35-1)*x(294-1)-P*p(35-1)*p(36-1)*x(294-1)+P*p(35-1)*x8*x(294-1)-P*p(36-1)*x8*x(294-1)+p(35-1)*p(36-1)*x8*x(294-1)+p(36-1)*pow(x8,2)*x(294-1)+P*p(35-1)*x6*x(296-1)+P*p(36-1)*x6*x(296-1)+p(35-1)*p(36-1)*x6*x(296-1)))/pow(P*p(35-1)+p(36-1)*(p(35-1)+x8),2)+(p(31-1)*p(33-1)*((x7-x8)*(p(33-1)*x7+p(32-1)*(p(33-1)+x8))*x(293-1)+x5*((p(33-1)*x8+p(32-1)*(p(33-1)+x8))*x(295-1)-(p(33-1)*x7+p(32-1)*(p(33-1)+x7))*x(296-1))))/pow(p(33-1)*x7+p(32-1)*(p(33-1)+x8),2);

		return result;
	}*/

	mat nonlinearOdes::toggle_switch(const double& t, const vec& x, const vec& u)
	{
		(void)t;

		double kL1 = u(0);
		double kL2 = u(1);
		double KT = u(2);
		double nT = u(3);
		double dL1 = u(4);
		double dL2 = u(5);

		double kT1 = u(6);
		double kT2 = u(7);
		double KL = u(8);
		double nL = u(9);
		double dT1 = u(10);
		double dT2 = u(11);

		mat result(4,1);
		result << 	kL1*pow(KT,nT)/(pow(KT,nT) + pow(x(3),nT)) - dL1*x(0),
					kL2*x(0) - dL2*x(1),
					kT1*pow(KL,nL)/(pow(KL,nL) + pow(x(1),nL)) - dT1*x(2),
					kT2*x(2) - dT2*x(3);

		return result;
	}
	mat nonlinearOdes::toggle_switch_config1(const double& t, const vec& x, const vec& u)
	{
		(void)t;

		double kL1 = u(0);
		double kL2 = u(1);
		double KT = u(2);
		double nT = 2;
		double dL1 = u(3);
		double dL2 = u(4);

		double kT1 = u(5);
		double kT2 = u(6);
		double KL = u(7);
		double nL = 2;
		double dT1 = u(8);
		double dT2 = u(9);


		mat result(4,1);
		result << 	kL1*pow(KT,nT)/(pow(KT,nT) + pow(x(3),nT)) - dL1*x(0),
					kL2*x(0) - dL2*x(1),
					kT1*pow(KL,nL)/(pow(KL,nL) + pow(x(1),nL)) - dT1*x(2),
					kT2*x(2) - dT2*x(3);
				nonlinearOdes::addCounter();
		return result;
	}
	mat nonlinearOdes::toggle_switch_config2(const double& t, const vec& x, const vec& u)
	{
		(void)t;

		double kL1 = u(0);
		double kL2 = u(1);
		double KTnT = u(2);
		double nT = 2;
		double dL1 = u(3);
		double dL2 = u(4);

		double kT1 = u(5);
		double kT2 = u(6);
		double KLnL = u(7);
		double nL = 2;
		double dT1 = u(8);
		double dT2 = u(9);

		mat result(4,1);
		result << 	kL1/(1 + pow(x(3),nT)/KTnT) - dL1*x(0),
					kL2*x(0) - dL2*x(1),
					kT1/(1 + pow(x(1),nL)/KLnL) - dT1*x(2),
					kT2*x(2) - dT2*x(3);
		nonlinearOdes::addCounter();
		return result;
	}
	mat nonlinearOdes::toggle_switch_config3(const double& t, const vec& x, const vec& u)
	{
		(void)t;

		double kL1 = u(0);
		double kL2 = u(1);
		double KTnT = u(2);
		double nT = u(3);
		double dL1 = u(4);
		double dL2 = u(5);

		double kT1 = u(6);
		double kT2 = u(7);
		double KLnL = u(8);
		double nL = u(9);
		double dT1 = u(10);
		double dT2 = u(11);


		mat result(4,1);
		result << 	kL1/(1 + pow(x(3),nT)/KTnT) - dL1*x(0),
					kL2*x(0) - dL2*x(1),
					kT1/(1 + pow(x(1),nL)/KLnL) - dT1*x(2),
					kT2*x(2) - dT2*x(3);
		nonlinearOdes::addCounter();
		return result;
	}
	mat nonlinearOdes::toggle_switch_config4(const double& t, const vec& x, const vec& u)
	{
		(void)t;

		double kL1 = u(0);
		double kL2 = u(1);
		double KT = 1;
		double nT = u(2);
		double dL1 = u(3);
		double dL2 = u(4);

		double kT1 = u(5);
		double kT2 = u(6);
		double KL = 1;
		double nL = u(7);
		double dT1 = u(8);
		double dT2 = u(9);


		mat result(4,1);
		result << 	kL1*pow(KT,nT)/(pow(KT,nT) + pow(x(3),nT)) - dL1*x(0),
					kL2*x(0) - dL2*x(1),
					kT1*pow(KL,nL)/(pow(KL,nL) + pow(x(1),nL)) - dT1*x(2),
					kT2*x(2) - dT2*x(3);
				nonlinearOdes::addCounter();
		return result;
	}
	mat nonlinearOdes::repressilator(const double& t, const vec& x, const vec& u)
	{
		(void) t;
		double a = u(0);
		double a0 = u(1);
		double n = u(2);
		double b = u(3);

		mat result(6,1);
		result << 	-x(0) + a/(1 + pow(x(5),n)) + a0,
					-b*(x(1) - x(0)),
					-x(2) + a/(1 + pow(x(1),n)) + a0,
					-b*(x(3) - x(2)),
					-x(4) + a/(1 + pow(x(3),n)) + a0,
					-b*(x(5) - x(4));
		nonlinearOdes::addCounter();
		return result;
	}

	mat nonlinearOdes::general_repressilator(const double& t, const vec& x, const vec& u)
	{
		(void)t;
		double a = u(0);
		double b = u(1);
		double n = u(2);
		double y = u(3);
		double d = u(4);

		double aB = u(5);
		double bB = u(6);
		double nB = u(7);
		double yB = u(8);
		double dB = u(9);

		double aC = u(10);
		double bC = u(11);
		double nC = u(12);
		double yC = u(13);
		double dC = u(14);

		mat result(6,1);
		result << 	-x(0) + a/(1 + pow(x(5),n)) + b,
					-y*x(1) + d*x(0),
					-x(2) + aB/(1 + pow(x(1),nB)) + bB,
					-yB*x(3) + dB*x(2),
					-x(4) + aC/(1 + pow(x(3),nC)) + bC,
					-yC*x(5) + dC*x(4);
					
		return result;
	}
}
