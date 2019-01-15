#ifndef ODE_WRAPPER_H
#define ODE_WRAPPER_H

#include <misc.h>
#include <spline.h>

namespace thesis {
    class OdeWrapper
    {
    private:
        odefunction fun;
        vec control;
        vec u;
        vector<spline> Xn;

        vec qlinear(const double& t, const vec& x);
    	vec der(const vec& dx, const double& dt);
    	mat jac(double t, const vec& x, const double& h);
        vec vectorToVec(vector<double> x);
        vector<double> vecToVector(const vec& x);
        //vector_type vecToUVector(const vec& x);
        //vec uvectorToVec(const vector_type& x);
        //mat umatToMat(const matrix_type& m);
        //matrix_type matToUMat(const mat& m);

    public:
        OdeWrapper(odefunction of, vec input, vec parameters);
        void setControl(const vec& input);
        void setParameter(const vec& parameters);
        void setPreviousIteration(vector<spline> &prevIter);
        vec fhandle(double t, vec  x, vec u)
        {
            //TODO check that the control has been set.
            return fun(t, x, u, control);
        }
        void operator() (const vector<double> & x, vector<double> &dxdt, const double t){
            dxdt = vecToVector(qlinear(t, vectorToVec(x)));
        }

        void call( const vector_type &x , vector_type &dxdt , double t)
        {
            dxdt = vecToVector(qlinear(t, vectorToVec(x)));
        }

        //void jacobian( const vector_type & /* x */ , matrix_type &J , const double & t , vector_type &dfdt );
    };
}
#endif /* end of include guard: ODE_WRAPPER_H */
