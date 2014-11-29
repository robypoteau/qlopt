#ifndef SPLINE_H
#define SPLINE_H

#include <unsupported/Eigen/MPRealSupport>
#include "misc.h"

namespace thesis{

	struct sizeDoNotMatchError: public exception
	{
	  const char* what() const throw()
	  {
		return "Eigen Matrix Sizes Do Not Match";
	  }
	};
	
	class spline{
		private:
			vec x, y, b, c, d;
			int n;
			//mpreal min(vec t);
			//int indexOf(vec t, mpreal ti);
			
		public:
			spline(){}
			spline(vec x, vec y);
			void update(vec x, vec y);
			mpreal interpolate(mpreal t);
			
			void setX(vec x){ 
				if(x.size() != this->n){ throw sizeDoNotMatchError();}
				this->x = x;
			}
			void setY(vec y){
				if(y.size() != this->n){ throw sizeDoNotMatchError();}
				this->y = y;
			}
			
			//vec getA(){ return y;}
			//vec getB(){ return b;}
			//vec getC(){ return c;}
			//vec getD(){ return d;}
			//int getN(){ return n;}
	};
	
	inline spline::spline(vec x, vec y){
		spline::update(x,y);
	}

	inline void spline::update(vec x, vec y){
		this->n = x.size();
		
		try{
			spline::setX(x);
			spline::setY(y);
		}
		catch(sizeDoNotMatchError& e){
			cout << e.what() << endl;
		}
		
		b = vec::Zero(n-1);
		d = b;
		c = vec::Zero(n);
		vec h(n-1), alpha(n-1), mu(n-1), l(n), z(n);
		//find h
		for(int i=0; i<(n-1); i++){
			h(i) = x(i+1) - x(i);
		}
		//find alphas
		for(int i=1; i<(n-1); i++){
			alpha(i) = 3/h(i)*(y(i+1) - y(i)) - 3/h(i-1)*(y(i) - y(i-1));
		}
		//init l, mu, z
		l(0) = l(n-1) = 1;
		mu(0) = 0;
		z(0) = z(n-1) = 0;
		c(n-1) = 0;
		
		//fill in l mu z
		for(int i=1; i<(n-1); i++){
			l(i) = 2*(x(i+1) - x(i-1)) - h(i-1)*mu(i-1);
			mu(i) = h(i)/l(i);
			z(i) = (alpha(i) - h(i-1)*z(i-1))/l(i);
		}
		// find c, b, d
		for(int i=(n-2); i>=0; i--){
			c(i) = z(i) - mu(i)*c(i+1);
			b(i) = (y(i+1) - y(i))/h(i) - h(i)*(c(i+1)+2*c(i))/3;
			d(i) = (c(i+1) - c(i))/(3*h(i));
		}
	}

	inline mpreal spline::interpolate(mpreal ti){
		int i = 0;
		
		for(int j=1; j<n-1 ; j++){ // turn into a function
			if(x(j) <= ti){
				i = j;
				break;
			}
		}
		
		mpreal diff = ti - x(i);
		return y(i) + diff*(b(i) + c(i)*diff + d(i)*diff*diff);
	}
}
#endif