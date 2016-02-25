#include <mp_nonlinear_odes.h>
#include <complex>

using namespace mpfr;

namespace thesis{
	
	mpnonlinearOdes::mpnonlinearOdes(){
		mpOdeFuncMap["lotka_volterra"] = mpnonlinearOdes::lotka_volterra;
		//qLinFuncMap["lotka_volterra_linearization"] = mpnonlinearOdes::lotka_volterra_linearization;
		mpOdeFuncMap["pielou"] = mpnonlinearOdes::pielou;
		//qLinFuncMap["pielou_linearization"] = mpnonlinearOdes::pielou_linearization;
		mpOdeFuncMap["angiogenesis"] = mpnonlinearOdes::angiogenesis;
		mpOdeFuncMap["cancer"] = mpnonlinearOdes::cancer;
		//qLinFuncMap["angiogenesis_linearization"] = mpnonlinearOdes::angiogenesis_linearization;
		//mpOdeFuncMap["coral"] = mpnonlinearOdes::coral;
		mpOdeFuncMap["coral5"] = mpnonlinearOdes::coral5;
		//mpOdeFuncMap["coral_pw"] = mpnonlinearOdes::coral_pw;
		//mpOdeFuncMap["coral_two"] = mpnonlinearOdes::coral_two;
		//mpOdeFuncMap["coral_four"] = mpnonlinearOdes::coral_four;
		//qLinFuncMap["coral_linearization"] = mpnonlinearOdes::coral_linearization;
		mpOdeFuncMap["bistable_switch"] = mpnonlinearOdes::bistable_switch;
		mpOdeFuncMap["bistable_switch_two"] = mpnonlinearOdes::bistable_switch_two;
		//qLinFuncMap["bistable_switch_linearization"] = mpnonlinearOdes::bistable_switch_linearization;
		mpOdeFuncMap["eight_part"] = mpnonlinearOdes::eight_part;
		mpOdeFuncMap["eight_part_spc"] = mpnonlinearOdes::eight_part_spc;
		//qLinFuncMap["eight_part_linearization"] = mpnonlinearOdes::eight_part_linearization;
		mpOdeFuncMap["gen_switch"] = mpnonlinearOdes::gen_switch;
	}
	
	mp_mat mpnonlinearOdes::lotka_volterra(const mpreal& t, const mp_vec& x, const mp_vec& u)
	{
		mp_mat result(2,1); 
		result(0) = u(0)*x(0) - u(1)*x(0)*x(1);
		result(1) = u(1)*x(0)*x(1) - u(2)*x(1);
			
		return result;
	}
	
	mp_mat mpnonlinearOdes::pielou(const mpreal& t, const mp_vec& x, const mp_vec& u)
	{
		mpreal k = 250;
		
		mp_mat result(2,1); 
		result(0) = u(0)*(1-x(0)/k)*x(0) - u(1)*x(0)*x(1);
		result(1) = -u(2)*x(1) + u(3)*x(0)*x(1);
		
		return result;
	}
	
	mp_mat mpnonlinearOdes::angiogenesis(const mpreal& t, const mp_vec& x, const mp_vec& u)
	{
		mpreal m = u(0); // m>0
		mpreal n = 1; //u(1); // n = .25, = 1 
		mpreal a = u(2); // a>=0
		mpreal c = u(3); // c(t) = c = 0 => no treatment, c = const implies fixed treatment
		mpreal w = u(4); // w>=0
		mpreal g = u(5); // g>=0
		
		/*
		mpreal m = u(0); // m>0
		mpreal n = .25; // n = .25, = 1 
		mpreal a = u(1); // a>=0
		mpreal c = u(2); // c(t) = c = 0 => no treatment, c = const implies fixed treatment
		mpreal w = u(3); // w>=0
		mpreal g = u(4); // g>=0
		*/
		
		mpreal f;
		if(n == 0){
			f = -m*log(x(0)/x(1));
		}else{
			f = m*x(0)/n*(1-pow(x(0)/x(1),n));
		}
	
		mp_mat result(2,1); 
		result(0) = -a*c*x(1) + f;
		result(1) = -a*c*x(1) + w*x(0) - g*pow(x(0)*x(0),1/3)*x(1);
		
		return result;
	}
	
	mp_mat mpnonlinearOdes::cancer(const mpreal& t, const mp_vec& x, const mp_vec& u)
	{
		mpreal d   = u(0); // >=0
		mpreal tau = u(1); // >=0 
		mpreal g   = u(2); // >=0
		mpreal a   = u(3); // >=0, a <= g/2
		mpreal k   = u(4); // >=0		
		mpreal R   = 1.0;	
		mp_mat result(2,1); 
		result(0) = d*R - x(0)/tau - g*x(0)*x(0);
		result(1) = -(a*R + k*x(0))*x(1);
		
		return result;
	}
	
	mp_mat mpnonlinearOdes::coral5(const mpreal& t, const mp_vec& x, const mp_vec& u)
	{
		//input
		mpreal m = u(0);
		mpreal r = u(1);
		mpreal k1 = u(2);
		//mpreal p = u(0); 
		mpreal k2 = u(3);
		mpreal k3 = u(4);

		//Pre-set parameter
		mpreal a2 = 1.5;
		mpreal b2 = .0013;
		mpreal N = 1750;
		mpreal I0 = 2000;
		mpreal pi = atan(1)*4;
		
		mp_mat result(2,1); 
		result(0) = m*x(1-1) + r*x(1-1)*(1 - x(1-1)/N) - k1*x(2-1)*x(1-1);// - p*x(1-1);
		result(1) = k2*x(2-1)+ k3*a2*b2*pow(2,b2*sin(pi*t/12))*I0*log(2)/12*pi*cos(pi*t/12);
		
		return result;
	}
	
	mp_mat mpnonlinearOdes::bistable_switch(const mpreal& t, const mp_vec& x, const mp_vec& p)
	{
		mpreal alpha = p(0);
		mpreal u     = p(1);
		mpreal n     = p(2);

		mp_mat result(2,1); 
		result(0) = alpha/(1 + pow(u*x(2-1),n)) - x(1-1);
		result(1) = alpha/(1 + pow(x(1-1),n))- x(2-1);
			
		return result;
	}
	
	mp_mat mpnonlinearOdes::bistable_switch_two(const mpreal& t, const mp_vec& x, const mp_vec& p)
	{
		mpreal alpha = p(0);
		mpreal u     = 1;//3.2;
		mpreal n     = p(1);

		mp_mat result(2,1); 
		result(0) = alpha/(1 + pow(u*x(1),n)) - x(0);
		result(1) = alpha/(1 + pow(x(0),n))- x(1);
			
		return result;
	}

	mp_mat mpnonlinearOdes::eight_part(const mpreal& t, const mp_vec& x, const mp_vec& pp)
	{
		mpreal P = .05;
		mpreal S = 10;
		
		mp_vec p(36);
		if(pp.size() < 36){
			p << 1.0,1.0,2.0,1.0,2.0,1.0,1.0,1.0,2.0,1.0,2.0,1.0,1.0,1.0,2.0,1.0,2.0,1.0,0.1,1.0,0.1,0.1,1.0,0.1,0.1,1.0,0.1,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0;
			for(int i=0; i<pp.size(); i++){
				p(i) = pp(i);
			}
		}else{
			p << pp;
		}
		mp_mat result(8,1); 
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
	
	 mp_mat mpnonlinearOdes::eight_part_spc(const mpreal& t, const mp_vec& x, const mp_vec& pp)
	{
		mpreal P = 1;
		mpreal S = .1;
		
		mp_vec p(36);
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
		mp_mat result(8,1); 
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
	
	mp_mat mpnonlinearOdes::gen_switch(const mpreal& t, const mp_vec& x, const mp_vec& u)
	{
		
		mpreal kL1 = u(0); // 
		mpreal kL2 = u(1); // 
		mpreal KT  = 10.0;
		mpreal nT  = 2; // >= 2	
		mpreal dL1 = u(2); // ~5mins(2-8)
		mpreal dL2 = u(3); // ~30 mins(20-40)
		
		mpreal kT1 = u(4); // 
		mpreal kT2 = u(5); // 
		mpreal KL  = 10.0;
		mpreal nL  = 2; // >= 2	
		mpreal dT1 = u(6); //~5mins
		mpreal dT2 = u(7); //~30 mins
		
		/*
		mpreal kL1 = u(0);
		mpreal kL2 = u(1);
		mpreal KT = 10.0;
		mpreal nT = u(2);
		mpreal dL1 = u(3);
		mpreal dL2 = u(4);
		
		mpreal kT1 = u(5);
		mpreal kT2 = u(6);
		mpreal KL = 10.0;
		mpreal nL = u(7);
		mpreal dT1 = u(8);
		mpreal dT2 = u(9);
		*/
		/*
		mpreal kL1 = u(0);
		mpreal kL2 = u(1);
		mpreal KT = u(2);
		mpreal nT = u(3);
		mpreal dL1 = u(4);
		mpreal dL2 = u(5);
		
		mpreal kT1 = u(6);
		mpreal kT2 = u(7);
		mpreal KL = u(8);
		mpreal nL = u(9);
		mpreal dT1 = u(10);
		mpreal dT2 = u(11);
		*/
		mp_mat result(4,1); 
		result << 	kL1*pow(KT,nT)/(pow(KT,nT) + pow(x(3),nT)) - dL1*x(0),
					kL2*x(0) - dL2*x(1),
					kT1*pow(KL,nL)/(pow(KL,nL) + pow(x(1),nL)) - dT1*x(2),
					kT2*x(2) - dT2*x(3);
				
		return result;
	}
}