#include <odesolver.h>

namespace thesis{

    mat OdeIntWrapper(OdeWrapper sys, const vec& x, const vec& t)
    {
        vector_type v = vecToVector(x);
        vector<vector<double>> x_vec;
        vector<double> times;

        //runge_kutta_fehlberg78<vector_type> stepper;
        //runge_kutta4<vector_type> stepper;
        runge_kutta_dopri5< vector_type > stepper;

		auto start = high_resolution_clock::now();
        size_t num_of_steps = integrate_adaptive( stepper,
            [&sys](const vector_type &x , vector_type &dxdt , double t){sys.call(x,dxdt,t);},
            v , t(0), t(t.size()-1), t(1)-t(0),
           push_back_state_and_time(x_vec, times));

        // Get ending timepoint
       	auto end = high_resolution_clock::now();
   		auto duration = duration_cast<microseconds>(end - start);

   		cout << "OdeInt Time: " << duration.count() << " ms" << endl;

        size_t xlen = x_vec[0].size();
        size_t tlen = times.size();
        mat output(xlen+1, tlen); //times and x values

        output.row(0) = vectorToVec(times);
        for( size_t i=0; i<tlen; i++ )
        {
           for(size_t j=0; j<xlen; j++ )
           {
               output(j+1, i) = x_vec[i][j];
           }
        }
        //cout << output.topRows(8).transpose() << endl; exit(0);
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
