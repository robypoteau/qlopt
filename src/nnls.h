#ifndef NNLS_H
#define NNLS_H

#include <dbg.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>

#include <misc.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_linalg.h>

namespace thesis{
	class nnls {
		private:
			size_t num_obs;
			size_t order;
			gsl_matrix *X;
            mat A;
			gsl_vector *gslx, *gslc;

			void vecToGslVec(const vec& v, gsl_vector *gslv);
			mat gslMatToMat(gsl_matrix *X);
            //int argmax(const vec v, const int i);
			static void* init_numpy() {
		        import_array();
		        return NULL;
			}
            void generateX(gsl_matrix *M, gsl_vector *v);
			void generate_xi(double xi, gsl_vector *gslv);

		public:
			nnls(){}
			~nnls();
			nnls(const size_t n, const size_t p);
			void init(const size_t n, const size_t p);
			void update(const vec& x, const vec& y);
			double interpolate(double xi);
	};
}

#endif
