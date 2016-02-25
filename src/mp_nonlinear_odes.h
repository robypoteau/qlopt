#ifndef MP_NONLINEAR_ODES_H
#define MP_NONLINEAR_ODES_H

#include <map>
#include <misc.h>
#include <mp_spline.h>

namespace thesis{
	class mpnonlinearOdes{
		typedef mp_mat (*mpOdeFuncPtr)(const mpreal& , const mp_vec&, const mp_vec&);
		
	public:
		map<string, mpOdeFuncPtr> mpOdeFuncMap;
		
		mpnonlinearOdes();
		
		static mp_mat lotka_volterra(const mpreal& t, const mp_vec& x, const mp_vec& u);
		static mp_mat pielou(const mpreal& t, const mp_vec& x, const mp_vec& u);
		static mp_mat angiogenesis(const mpreal& t, const mp_vec& x, const mp_vec& u);
		static mp_mat cancer(const mpreal& t, const mp_vec& x, const mp_vec& u);
		static mp_mat coral5(const mpreal& t, const mp_vec& x, const mp_vec& u);
		static mp_mat bistable_switch(const mpreal& t, const mp_vec& x, const mp_vec& u); 
		static mp_mat bistable_switch_two(const mpreal& t, const mp_vec& x, const mp_vec& p);
		static mp_mat eight_part(const mpreal& t, const mp_vec& x, const mp_vec& u);
		static mp_mat eight_part_spc(const mpreal& t, const mp_vec& x, const mp_vec& u);
		static mp_mat gen_switch(const mpreal& t, const mp_vec& x, const mp_vec& u);
	};
}
#endif