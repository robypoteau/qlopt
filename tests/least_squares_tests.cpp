#include <least_squares.h>
#define BOOST_TEST_MODULE LeastSquaresTest
#include <boost/test/included/unit_test.hpp>

using namespace thesis;

BOOST_AUTO_TEST_CASE(constructor_test)
{
	size_t n=7, p=3;
	//create and init
	lsquares ls;
	vec x(n);
	vec y(n);
	x << 1,2,3,4,5,6,7;
	y << 6,1,1,7,6,2,9;
	ls.init(n,p);

	try{
		lsquares ls2(n,p);
		ls2.update(x,y);
	}catch(exception& e){
		cout << "Second constructor failed" << '\n';
		BOOST_CHECK(false);
	}
}
