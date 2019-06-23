#include <odeWrapper.h>
namespace thesis{
    OdeWrapper::OdeWrapper(odefunction of, vec input, vec parameters):
        fun(of), control(input), u(parameters){}

    /*vec OdeWrapper::fcvs(double t, vec  x)
    {
        return flinear(t, x);
        //return fun(t, x, u, control);
    }*/

    void OdeWrapper::setControl(const vec& input)
    {
        control = input;
    }

    void OdeWrapper::setParameter(const vec& parameters)
    {
        u = parameters;
    }

    void OdeWrapper::setPreviousIteration(std::vector<spline> &prevIter)
    {
        Xn = prevIter;
    }
    
    vec OdeWrapper::flinear(const double& t, const vec& x, const vec& u)
    {
        size_t m = u.size();
        size_t n = this->Xn.size();
        double step = sqrt(2.2E-16);

        vec xn1(n); // this is x_N-1
        vec dxn(n);

        for(size_t i=0; i<n; i++)
        {

            xn1(i) = this->Xn[i].interpolate(t);
            dxn(i) = x(i) - xn1(i); // x_N - x_{N-1}
        }

        mat dx(n,n);// for x derivative
        dx << mat::Identity(n,n)*step;

        //First few lines of linearization
        vec fx(n);
        fx = fhandle(t, xn1, u);
        vec ans(n+n*m);
        ans << vec::Zero(n+n*m);
        ans.segment(0, n) = fx + jac(t, xn1, step)*dxn;

        return ans;
    }

    vec OdeWrapper::qlinear(const double& t, const vec& x)
    {
        size_t m = u.size();
        size_t n = this->Xn.size();
        double step = sqrt(2.2E-16);

        vec xn1(n); // this is x_N-1
        vec dxn(n);

        for(size_t i=0; i<n; i++)
        {

            xn1(i) = this->Xn[i].interpolate(t);
            dxn(i) = x(i) - xn1(i); // x_N - x_{N-1}
        }

        mat dx(n,n);// for x derivative
        dx << mat::Identity(n,n)*step;

        //First few lines of linearization
        vec fx(n);
        fx = fhandle(t, xn1, u);
        vec ans(n+n*m);
        ans << vec::Zero(n+n*m);
        ans.segment(0, n) = fx + jac(t, xn1, step)*dxn;

        vec dfdx(n);
        vec dfdu(n);

        // the Un part of t the linearization
        mat dun(m, m);// for u derivative
        dun << mat::Identity(m,m)*step;

        size_t ind;
        for(size_t j=0; j<m; j++)
        {
            ind = (j+1)*n;
            dfdu = fhandle(t, xn1, u+dun.col(j)*u(j));
            ans.block(ind, 0, n, 1) = der(dfdu - fx, step*u(j)); // df/du
            for(size_t k=0; k<n; k++){
                dfdx = fhandle(t, xn1+dx.col(k)*xn1(k), u);
                ans.block(ind, 0, n, 1) += der(dfdx  - fx, step*xn1(k))*x(ind+k);//J*Un

                ans.block(ind, 0, n, 1) += der(
                      fhandle(t, xn1+dx.col(k)*xn1(k), u+dun.col(j)*u(j))
                    - fhandle(t, xn1+dx.col(k)*xn1(k), u-dun.col(j)*u(j))
                    - fhandle(t, xn1-dx.col(k)*xn1(k), u+dun.col(j)*u(j))
                    + fhandle(t, xn1-dx.col(k)*xn1(k), u-dun.col(j)*u(j)),
                    4*step*xn1(k)*step*u(j))*dxn(k); //phi_ij
            }
        }
        //cout << "ans: " << ans.transpose() << endl;
        return ans;
    }

    vec OdeWrapper::der(const vec& dx, const double& dt)
    {
        mat ans(dx.size(),1);
        ans << dx/dt;
        return ans;
    }

    mat OdeWrapper::jac(double t, const vec& x, const double& h)
    {
        int n = x.size();
        mat dx(n,n);
        dx << mat::Identity(n,n)*h;
        mat fprime(n,n);
        for(int j=0; j<n; j++)
        {
            fprime.col(j) = (-fhandle(t, x+2*dx.col(j)*x(j), u)
                + 8*fhandle(t, x+dx.col(j)*x(j), u)
                - 8*fhandle(t, x-dx.col(j)*x(j), u)
                + fhandle(t, x-2*dx.col(j)*x(j), u))/x(j);
        }

        return fprime/(12*h);
    }

    /*void OdeWrapper::jacobian( const vector_type &x  , matrix_type &J , const double & t , vector_type &dfdt ){
        size_t n = this->Xn.size();
        double step = sqrt(2.2E-16);
        vec xn1(n);

        for(size_t i=0; i<n; i++)
            xn1(i) = this->Xn[i].interpolate(t);

        J = matToUMat(jac(t, xn1, step));
        dfdt = vecToUVector(vec::Zero(n));
        //dfdt = vecToUVector(der(qlinear(t,uvectorToVec(x))-qlinear(t+step*t,uvectorToVec(x)),step*t));
    }*/

    vec OdeWrapper::vectorToVec(vector<double> x)
    {
        Map<vec> v(x.data(),x.size());
        return v;
    }

    vector<double> OdeWrapper::vecToVector(const vec& x)
    {
        vector<double> v(x.data(), x.data()+x.size());
        return v;
    }

    /*vector_type OdeWrapper::vecToUVector(const vec& x){
        vector_type v(x.size());
        for (unsigned i = 0; i < v.size (); ++ i)
            v(i) = x(i);
        return v;
    }
    vec OdeWrapper::uvectorToVec(const vector_type& x){
        vec v(x.size());
        for (unsigned i = 0; i < v.size (); ++ i)
            v(i) = x(i);
        return v;
    }
    mat OdeWrapper::umatToMat(const matrix_type& m){
        mat a(m.size1(), m.size2());
        for (unsigned i = 0; i < m.size1 (); ++ i)
            for (unsigned j = 0; j < m.size2 (); ++ j)
                a(i, j) = m(i,j);
        return a;
    }
    matrix_type OdeWrapper::matToUMat(const mat& m){
        matrix_type a(m.rows(), m.cols());
        for (unsigned i = 0; i < m.rows(); ++ i)
            for (unsigned j = 0; j < m.cols(); ++ j)
                a(i, j) = a(i,j);
        return a;
    }*/
}
