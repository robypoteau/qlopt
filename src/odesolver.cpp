#include <odesolver.h>

namespace thesis{

    mat OdeIntWrapper(OdeWrapper sys, const vec& x, const vec& t, double tol)
    {
        vector_type v = vecToVector(x);
        vector<vector<double>> x_vec;
        vector<double> times;

        //runge_kutta4<vector_type> stepper;
        //runge_kutta_dopri5< vector_type > stepper;
        //runge_kutta_fehlberg78<vector_type> stepper;

        bulirsch_stoer< vector_type > stepper( tol, tol, 1.0 , 1.0 );
        //bulirsch_stoer_dense_out< vector_type > stepper;

        //controlled_runge_kutta< runge_kutta_fehlberg78<vector_type> > stepper;
        //controlled_runge_kutta< runge_kutta_cash_karp54<vector_type> > stepper;
        //controlled_runge_kutta< runge_kutta_dopri5<vector_type> > stepper( 1E-4, 1E-4 );

        auto start = high_resolution_clock::now();
        size_t num_of_steps = integrate_const( stepper,
            [&sys](const vector_type &x , vector_type &dxdt , double t){sys.call(x,dxdt,t);},
            v , t(0), t(t.size()-1), t(1)-t(0),
           push_back_state_and_time(x_vec, times));

        // Get ending timepoint
        auto end = high_resolution_clock::now();
   		auto duration = duration_cast<microseconds>(end - start);
        cout << "OdeInt Time: " << duration.count()/1E6 << " s" << endl;
        cout << "Number of steps: " << num_of_steps << endl;
        size_t xlen = x_vec[0].size();
        size_t tlen = times.size();
        mat output(xlen, tlen);

        for( size_t i=0; i<tlen; i++ )
        {
           for(size_t j=0; j<xlen; j++ )
           {
               output(j, i) = x_vec[i][j];
           }
        }
        return output;
    }

    vec vectorToVec(vector<double> x)
    {
        Map<vec> v(x.data(),x.size());
        return v;
    }

    vector<double> vecToVector(const vec& x)
    {
        vector<double> v(x.data(), x.data()+x.size());
        return v;
    }
}
