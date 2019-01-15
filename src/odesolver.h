#ifndef ODESOLVER_H
#define ODESOLVER_H

#include <misc.h>
#include <odeWrapper.h>
#include <chrono>
using namespace std::chrono;
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

namespace thesis {
    // integrate_observer
    struct push_back_state_and_time
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
    mat OdeIntWrapper(OdeWrapper sys, const vec& x, const vec& t);
    vec vectorToVec(vector<double> x);
    vector<double> vecToVector(const vec& x);
}

#endif
