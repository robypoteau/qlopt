#ifndef MP_SPLINE_H
#define MP_SPLINE_H

#include <misc.h>
#include <iostream>

using namespace mpfr;

namespace thesis{
	struct mpSizeDoNotMatchError: public exception
	{
	  const char* what() const throw()
	  {
		return "Eigen Matrix Sizes Do Not Match";
	  }
	};
	
	class mp_spline{
		private:
			mp_vec x, y, b, c, d;
			int n;
			
		public:
			mp_spline(){}
			mp_spline(mp_vec x, mp_vec y);
			void update(mp_vec x, mp_vec y);
			mpreal interpolate(mpreal ti);
			
			void setX(mp_vec x){ 
				if(x.size() != this->n){ throw mpSizeDoNotMatchError();}
				this->x = x;
			}
			void setY(mp_vec y){
				if(y.size() != this->n){ throw mpSizeDoNotMatchError();}
				this->y = y;
			}
	};
}
#endif