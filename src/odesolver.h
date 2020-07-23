#ifndef ODESOLVER_H
#define ODESOLVER_H

#include <misc.h>
#include <odeWrapper.h>
#include <chrono>
using namespace std::chrono;
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;


#include <iostream>
#include <cvodes/cvodes.h> // prototypes for CVODE fcts., consts.
//#include <cvodes/cvodes_dense.h>
//#include <cvodes/cvodes_lapack.h>
#include <cvodes/cvodes_bandpre.h>
//#include <cvodes/cvodes_spils.h> // access to CVSpils interface

#include <nvector/nvector_serial.h>  // access to serial N_Vector
#include <sunmatrix/sunmatrix_dense.h>

#include <sunlinsol/sunlinsol_dense.h>  //access to SPGMR SUNLinearSolver
#include <sunlinsol/sunlinsol_spgmr.h>  //access to SPGMR SUNLinearSolver
#include <sunlinsol/sunlinsol_spfgmr.h>  // use generic dense solver in precond
#include <sunlinsol/sunlinsol_spbcgs.h>  // use generic dense solver in precond
#include <sunlinsol/sunlinsol_sptfqmr.h>  // use generic dense solver in precond
#include <sunlinsol/sunlinsol_pcg.h>  // use generic dense solver in precond
//#include <sunlinsol/sunlinsol_lapackdense.h>  // use generic dense solver in precond

#include <sundials/sundials_dense.h>  // use generic dense solver in precond
#include <sundials/sundials_types.h>  // defs. of realtype, sunindextype
#include <sundials/sundials_math.h> // contains the macros ABS, SUNSQR, EXP
#include <sundials/sundials_nonlinearsolver.h>

#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

namespace thesis {

    typedef struct{
        OdeWrapper * pOdeWrapper;
        int N;
        realtype * p;
    } userData;

    // integrate_observer
    /*struct push_back_state_and_time
    {
        vector< vector<double> >& m_states;
        vector< double >& m_times;

        push_back_state_and_time(
            vector< vector<double> > &states,
            vector< double > &times )
        : m_states( states ) , m_times( times ) { }

        void operator()( const vector<double> &x , double t )
        {
            m_states.push_back( x );
            m_times.push_back( t );
        }
    };
    mat OdeIntWrapper(OdeWrapper sys, const vec& x, const vec& t, double tol);
    vec vectorToVec(vector<double> x);
    vector<double> vecToVector(const vec& x);*/
    mat OdeIntWrapper(OdeWrapper sys, const vec& x, const vec& t, double tol);
}

#endif
