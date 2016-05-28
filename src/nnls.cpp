#include <nnls.h>

namespace thesis{
	nnls::nnls(const size_t n, const size_t p)
	{
		init(n, p);
		Py_Initialize();
		log_info("Python Mode Initializes.");
		init_numpy();
	}

	void nnls::init(const size_t n, const size_t p)
	{
		num_obs = n;
		order = p;

        X    = gsl_matrix_alloc(num_obs, order);
		gslx = gsl_vector_alloc(num_obs);
		gslc = gsl_vector_alloc(order);
	}

	nnls::~nnls()
	{
		gsl_matrix_free(X);
		gsl_vector_free(gslx);
		gsl_vector_free(gslc);
		log_info("Python Mode Finalizes.");
		Py_Finalize();
	}

	void nnls::update(const vec& x, const vec& y)
	{
		vecToGslVec(x, gslx);
        generateX(X, gslx);
        A = gslMatToMat(X);
		double *c_out;
		check((int) *(A.data()) == 1 && (int) *(A.data()+1) == 0, *(A.data()));
		//Py_Initialize();
		//log_info("Python Mode Initializes.");
		//init_numpy();

		PyObject *pModule = PyImport_ImportModule("scipy.optimize");
		check(pModule, "pModule is not initialized.");

		PyObject *pFunc, *pArgs, *pResult, *pArg1, *pArg2;
		PyArrayObject *pNpArray;

		npy_intp Adims[2]; Adims[0] = A.rows(); Adims[1] = A.cols();
		check(Adims[0] == num_obs && Adims[1] == order, "A Matrix dims are wrong.");
		npy_intp  bdim[1];  bdim[0] = y.size();
		check(bdim[0] == num_obs, "b vector dims are wrong.");

		if (pModule != NULL) {
            pFunc = PyObject_GetAttrString(pModule, (char*)"nnls");
			check(pFunc, "pFunc is not initialized.");

			if (pFunc && PyCallable_Check(pFunc)) {
				pArg1 = PyArray_SimpleNewFromData(2, Adims, NPY_DOUBLE, const_cast<double*> (A.data()));
				pArg2 = PyArray_SimpleNewFromData(1, bdim , NPY_DOUBLE, const_cast<double*> (y.data()));
				pArgs = PyTuple_Pack(2, pArg1, pArg2);
				check(pArgs, "pArgs is not initialized.");

				pResult = PyObject_CallObject(pFunc, pArgs);
				check(pResult, "pResult is not initialized.");

				pNpArray = PyArray_GETCONTIGUOUS((PyArrayObject*) PyTuple_GetItem(pResult, 0));
				check(pNpArray, "pResult is not initialized.");

				check(PyArray_DIM(pNpArray, 0) == order, "Incorrect Dim for gslc.");
				c_out = reinterpret_cast<double*>(PyArray_DATA(pNpArray));
				cout << "NNLS coefficient vector" << endl ;
				for(size_t i=0; i<order; i++){
					gsl_vector_set(gslc, i, *(c_out+i));
					cout << gsl_vector_get(gslc,i) << ", " ;
				}
				cout << endl ;
				//log_info(PyFloat_AS_DOUBLE(PyTuple_GetItem(pResult, 1)));

				Py_XDECREF(pNpArray);
				Py_XDECREF(pResult);
				Py_XDECREF(pArgs);
			}
			Py_DECREF(pFunc);
        }
		Py_DECREF(pModule);

	}

    /*int nnls::argmax(const vec v, const int i){
        int temp = i;
        for(int j=i+1; j<v.size(); j++){
            if(vec(j)>vec(temp)) temp = j;
        }
        return temp;
    }*/

	void nnls::vecToGslVec(const vec& v, gsl_vector* gslv)
	{
		for (size_t i=0; i<gslv->size; i++){
			gsl_vector_set(gslv, i, v(i));
		}
	}

	mat nnls::gslMatToMat(gsl_matrix *gslM){
		mat M(gslM->size1, gslM->size2);
		for (size_t i=0; i<gslM->size1; i++){
			for (size_t j=0; j<gslM->size2; j++){
				M(i,j) = gsl_matrix_get(gslM, i, j);
			}
		}
		return M;
	}

	void nnls::generateX(gsl_matrix *M, gsl_vector* gslv)
	{
		gsl_vector *temp = gsl_vector_alloc(order);
		for (size_t i=0; i<M->size1; i++){
			generate_xi(gsl_vector_get(gslv, i), temp);
			gsl_matrix_set_row(M, i, temp);
		}
	}

	void nnls::generate_xi(double xi, gsl_vector *gslv)
	{
		for (size_t i=0; i<order; i++){
			gsl_vector_set(gslv, i, pow(xi,i));
		}
	}

	double nnls::interpolate(double xi)
	{
		double y = gsl_vector_get(gslc,0);
		for(size_t i = 1; i<order; i++){
			y += gsl_vector_get(gslc,i)*pow(xi,i);
		}
		return y;
	}
}
