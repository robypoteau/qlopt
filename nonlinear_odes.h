#ifndef NONLINEAR_ODES_H
#define NONLINEAR_ODES_H

#include <map>

#include "misc.h"
#include "spline.h"

namespace thesis{
	class nonlinearOdes{
		typedef mat (*odeFuncPtr)(const double& , const vec&, const vec&);
		typedef mat (*qLinFuncPtr)(const double& , const vec&, const vec&, const mat&, const vec&);

	public:
		map<string, odeFuncPtr> odeFuncMap;
		map<string, qLinFuncPtr> qLinFuncMap;
		
		nonlinearOdes();
		
		static mat lotka_volterra(const double& t, const vec& x, const vec& u);
		static mat lotka_volterra_linearization(const double& t, const vec& x, const vec& u, const mat& xn, const vec& time);
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
		static mat bistable_switch_linearization(const double& t, const vec& x, const vec& u, const mat& xn, const vec& time);
		static mat eight_part(const double& t, const vec& x, const vec& u);
		/* 
		static mat eight_part_spc(const double& t, const vec& x, const vec& u);
		static mat eight_part_linearization(const double& t, const vec& x, const vec& u, const mat& xn, const vec& time);
		*/
		static mat gen_switch(const double& t, const vec& x, const vec& u);
	};
}
#endif