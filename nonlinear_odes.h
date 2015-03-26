#ifndef NONLINEAR_ODES_H
#define NONLINEAR_ODES_H

#include <map>

#include "misc.h"
#include "spline.h"

namespace thesis{
	class nonlinearOdes{
		typedef mat (*odeFuncPtr)(const mpreal& , const vec&, const vec&);
		typedef mat (*qLinFuncPtr)(const mpreal& , const vec&, const vec&, const mat&, const vec&);

	public:
		map<string, odeFuncPtr> odeFuncMap;
		map<string, qLinFuncPtr> qLinFuncMap;
		
		nonlinearOdes();
		
		static mat lotka_volterra(const mpreal& t, const vec& x, const vec& u);
		static mat lotka_volterra_linearization(const mpreal& t, const vec& x, const vec& u, const mat& xn, const vec& time);
		static mat pielou(const mpreal& t, const vec& x, const vec& u);
		static mat pielou_linearization(const mpreal& t, const vec& x, const vec& u, const mat& xn, const vec& time);
		static mat angiogenesis(const mpreal& t, const vec& x, const vec& u);
		static mat angiogenesis_linearization(const mpreal& t, const vec& x, const vec& u, const mat& xn, const vec& time);
		static mat cancer(const mpreal& t, const vec& x, const vec& u);
		static mat coral(const mpreal& t, const vec& x, const vec& u);
		static mat coral5(const mpreal& t, const vec& x, const vec& u);
		static mat coral_pw(const mpreal& t, const vec& x, const vec& u);
		static mat coral_two(const mpreal& t, const vec& x, const vec& u);
		static mat coral_four(const mpreal& t, const vec& x, const vec& u);
		static mat coral_linearization(const mpreal& t, const vec& x, const vec& u, const mat& xn, const vec& time);
		static mat bistable_switch(const mpreal& t, const vec& x, const vec& u);
		static mat bistable_switch_two(const mpreal& t, const vec& x, const vec& p);
		static mat bistable_switch_linearization(const mpreal& t, const vec& x, const vec& u, const mat& xn, const vec& time);
		static mat eight_part(const mpreal& t, const vec& x, const vec& u);
		static mat eight_part_spc(const mpreal& t, const vec& x, const vec& u);
		static mat eight_part_linearization(const mpreal& t, const vec& x, const vec& u, const mat& xn, const vec& time);
		static mat jak_stat(const mpreal& t, const vec& x, const vec& u);
		static mat gen_switch(const mpreal& t, const vec& x, const vec& u);
	};
}
#endif