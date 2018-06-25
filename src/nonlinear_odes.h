#ifndef NONLINEAR_ODES_H
#define NONLINEAR_ODES_H

#include <map>
#include <dbg.h>
#include <misc.h>
#include <spline.h>

typedef mat  (*sys)(const double& t, const vec& x, const vec& u);
typedef mat (*qsys)(const double& t, const vec& x, const vec& u, vector<thesis::spline> xn);

namespace thesis{
	class nonlinearOdes{
		typedef mat (*odeFuncPtr)(const double& , const vec&, const vec&);
		typedef mat (*qLinFuncPtr)(const double& , const vec&, const vec&, std::vector<thesis::spline>);
		static int counter;
		
	public:
		map<string, odeFuncPtr> odeFuncMap;
		map<string, qLinFuncPtr> qLinFuncMap;

		nonlinearOdes();
		static void setCounter(int val){
			counter = val;
		}
		static int getCounter(){
			return counter;
		}
		static void addCounter(){
			counter++;
		}
		static mat lotka_volterra(const double& t, const vec& x, const vec& u);
		static mat lotka4(const double& t, const vec& x, const vec& u);
		static mat lorenz(const double& t, const vec& x, const vec& u);
		static mat general_lv(const double& t, const vec& x, const vec& u);
		static mat lotka_volterra_linearization(const double& t, const vec& x, const vec& u, std::vector<thesis::spline> xn);
		static mat pielou(const double& t, const vec& x, const vec& u);
		//static mat pielou_linearization(const double& t, const vec& x, const vec& u, const mat& xn, const vec& time);
		static mat angiogenesis(const double& t, const vec& x, const vec& u);
		//static mat angiogenesis_linearization(const double& t, const vec& x, const vec& u, const mat& xn, const vec& time);
		static mat cancer(const double& t, const vec& x, const vec& u);
		//static mat coral(const double& t, const vec& x, const vec& u);
		static mat coral5(const double& t, const vec& x, const vec& u);
		/*
		static mat coral_pw(const double& t, const vec& x, const vec& u);
		static mat coral_two(const double& t, const vec& x, const vec& u);
		static mat coral_four(const double& t, const vec& x, const vec& u);
		static mat coral_linearization(const double& t, const vec& x, const vec& u, const mat& xn, const vec& time);
		*/
		static mat bistable_switch(const double& t, const vec& x, const vec& u);
		static mat bistable_switch_two(const double& t, const vec& x, const vec& p);
		static mat bistable_switch_linearization(const double& t, const vec& x, const vec& u, std::vector<thesis::spline> xn);
		static mat eight_part(const double& t, const vec& x, const vec& u);
		static mat eight_part_spc(const double& t, const vec& x, const vec& u);
		/*
		static mat eight_part_linearization(const double& t, const vec& x, const vec& u, const mat& xn, const vec& time);
		*/
		static mat toggle_switch(const double& t, const vec& x, const vec& u);
		static mat toggle_switch_config1(const double& t, const vec& x, const vec& u);
		static mat toggle_switch_config2(const double& t, const vec& x, const vec& u);
		static mat toggle_switch_config3(const double& t, const vec& x, const vec& u);
		static mat toggle_switch_config4(const double& t, const vec& x, const vec& u);
		static mat repressilator(const double& t, const vec& x, const vec& u);
		static mat general_repressilator(const double& t, const vec& x, const vec& u);
	};
}
#endif
