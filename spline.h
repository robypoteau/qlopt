#ifndef SPLINE_H
#define SPLINE_H

#include "misc.h"
#include <iostream>

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
			
		public:
			spline(){}
			spline(vec x, vec y);
			void update(vec x, vec y);
			double interpolate(double ti);
			
			void setX(vec x){ 
				if(x.size() != this->n){ throw sizeDoNotMatchError();}
				this->x = x;
			}
			void setY(vec y){
				if(y.size() != this->n){ throw sizeDoNotMatchError();}
				this->y = y;
			}
			
			/*void printX(){
				cout << x.transpose() << endl;
			}
			
			void printA(){
				cout << y.transpose() << endl;
			}
			
			void printB(){
				cout << b.transpose() << endl;
			}
			
			void printC(){
				cout << c.transpose() << endl;
			}
			
			void printD(){
				cout << d.transpose() << endl;
			}
			
			void printAll(){
				printA();
				printB();
				printC();
				printD();
			}*/
	};
}
#endif