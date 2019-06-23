#include <regularization.h>

namespace thesis{

    double findAlpha1(const mat& A, const vec& P, const double Om, const double Nd){
        int m = A.cols();
        int ds = A.rows()/m - 1;
        vec sv = A.topRows(ds*m).jacobiSvd().singularValues();

        double divs = 24, upper = 6.0, lower = -6.0;
        double alpha = pow(10,lower), da = (upper - lower)/divs;

        mat I = mat::Identity(m,m);
        vec du(m);

        auto F = [A, sv, P, Om, I, Nd](double a){
            mat M = A;
            int m = A.cols();
            int ds = A.rows()/m - 1;
            vec temp = sv;
            temp.array() += a;
            temp = sv.array().square() / temp.array().square();
            M.bottomRows(m) = a*I;
            vec du = M.colPivHouseholderQr().solve(P);
            double top = Om;
            for(int i = 0; i < ds; i++)
            {
                top += du.transpose()*(M.middleRows(i*m, m)*du
                    - 2*P.segment(i*m, m));
            }
            return top/(Nd - temp.sum());
        };

        double FOld = F(alpha);
        double FNew;
        cout << "(" << pow(10,lower) << "," << FOld << ")"  << endl;

        for(double i=(lower+da); i<=upper; i+=da){

            FNew = F(pow(10,i));
            cout << "(" << pow(10,i) << "," << FNew << ")"  << endl;
            if(FNew < FOld){
                alpha = pow(10,i);
                FOld = FNew;
            }
        }

        //exit(0);
        return alpha;
    }

    double findAlpha2(const mat& A, const vec& P, const double Om, const double Nd){
        int m = A.cols();
        int ds = A.rows()/m - 1;
        vec sv = A.topRows(ds*m).jacobiSvd().singularValues();

        double divs = 24, upper = 6.0, lower = -6.0;
        double alpha = pow(10,lower), da = (upper - lower)/divs;
        vec du(m);
        mat I = mat::Identity(m,m);
        mat SumA = A.topRows(m);
        vec SumP = P.head(m);
        for(int i = 1; i < ds; i++)
        {
            SumA = SumA + A.middleRows(i*m, m);
            SumP = SumP + P.segment(i*m,m);
        }

        auto F = [sv, SumP, Om, I, Nd, SumA](double a){
            vec temp = sv;
            temp.array() += a;
            temp = sv.array().square() / temp.array().square();
            mat Aplus;
            Aplus = SumA + a*I;
            vec du = Aplus.colPivHouseholderQr().solve(SumP);
            double top = Om + du.transpose()*(SumA*du - 2*SumP);
            return top/(Nd - temp.sum());
        };

        double FOld = F(alpha);
        double FNew;
        cout << "(" << pow(10,lower) << "," << FOld << ")"  << endl;

        for(double i=(lower+da); i<=upper; i+=da){

            FNew = F(pow(10,i));
            cout << "(" << pow(10,i) << "," << FNew << ")"  << endl;
            if(FNew < FOld){
                alpha = pow(10,i);
                FOld = FNew;
            }
        }
        return alpha;
    }

    double findAlpha3(const mat& A, const vec& P, const double Om, const vec& ustar, const double Nd){
        int m = A.cols();
        int ds = A.rows()/m - 1;
        vec sv = A.topRows(ds*m).jacobiSvd().singularValues();

        double divs = 18, upper = 3.0, lower = -6.0;
        double alpha = pow(10,lower), da = (upper - lower)/divs;

        mat I;
        vec du(m);

        I = mat::Identity(m,m);
        mat SumA = A.topRows(m);
        vec SumP = P.head(m);
        for(int i = 1; i < ds; i++)
        {
            SumA = SumA + A.middleRows(i*m, m);
            SumP = SumP + P.segment(i*m,m);
        }

        auto F = [sv, SumP, Om, I, Nd, SumA, ustar](double a){
            vec temp = sv;
            temp.array() += a;
            temp = sv.array().square() / temp.array().square();
            mat Aplus;
            vec Pplus;
            Aplus = SumA + a*I;
            Pplus = SumP + a*ustar;
            vec du = Aplus.colPivHouseholderQr().solve(Pplus);
            double top = Om + du.transpose()*(Aplus*du - 2*Pplus);
            return top/((Nd - temp.sum()));
        };

        double FOld = F(alpha);
        double FNew;
        cout << "(" << pow(10,lower) << "," << FOld << ")"  << endl;

        for(double i=(lower+da); i<=upper; i+=da){

            FNew = F(pow(10,i));
            cout << "(" << pow(10,i) << "," << FNew << ")"  << endl;
            if(FNew < FOld){
                alpha = pow(10,i);
                FOld = FNew;
            }
        }
        return alpha;
    }

    double findGamma(const mat& A, const vec& P, const vec& uNot, const vec& u){
        mat M  = A;
        int m = M.cols();

        double divs = 240, upper = 4.0, lower = -4.0;
        double alpha = pow(10,lower), da = (upper - lower)/divs;

        mat I = mat::Identity(m,m);

        auto F = [M, P, I, uNot, u](double a){
            int m = M.cols();
            mat A = M;
            A.bottomRows(m) = a*I;

            vec top = uNot + A.colPivHouseholderQr().solve(P) - u ;
            return top.norm();
        };

        double FOld = F(alpha);
        double FNew;

        //cout << "(" << pow(10,lower) << "," << FOld << ")"  << endl;

        for(double i=(lower+da); i<=upper; i+=da){

            FNew = F(pow(10,i));
            //cout << "(" << pow(10,i) << "," << FNew << ")"  << endl;
            if((FNew - FOld) < 1E-6){
                alpha = pow(10,i);
                FOld = FNew;
            }
        }
        return alpha;
    }

    double findGamma2(const mat& A, const vec& P, const vec& uNot, const vec& uguess, const vec& u){
        mat M  = A;
        int m = M.cols();
        int ds = M.rows()/m - 1;

        double divs = 240, upper = 4.0, lower = -4.0;
        double alpha = pow(10,lower), da = (upper - lower)/divs;

        mat I = mat::Identity(m,m);
        mat SumA = A.topRows(m);
        vec SumP = P.head(m);
        for(int i = 1; i < ds; i++)
        {
            SumA = SumA + A.middleRows(i*m, m);
            SumP = SumP + P.segment(i*m,m);
        }

        auto F = [SumA, SumP, I, uNot, uguess, u](double a){
            mat Aplus;
            vec Pplus;
            Aplus = SumA + a*I;
            Pplus = SumP + a*uguess;

            vec top = uNot + Aplus.colPivHouseholderQr().solve(Pplus) - u ;
            return top.norm();
        };

        double FOld = F(alpha);
        double FNew;

        //cout << "(" << pow(10,lower) << "," << FOld << ")"  << endl;

        for(double i=(lower+da); i<=upper; i+=da){

            FNew = F(pow(10,i));
            //cout << "(" << pow(10,i) << "," << FNew << ")"  << endl;
            if((FNew - FOld) < 1E-6){
                alpha = pow(10,i);
                FOld = FNew;
            }
        }
        return alpha;
    }

    double findAlpha5(const mat& A, const vec& P, const vec& u0,
                OdeWrapper ow, const vector<mat>& msmt, const vector<vec>& input,
                vec y0, vec ts, const vector<vector<thesis::spline>>& spl_pairs){

        double divs = 9, upper = 3.0, lower = -6.0;
        double alpha = pow(10,lower), da = (upper - lower)/divs;

        auto F = [u0, msmt, input, P, ts, y0](double a, mat  A, OdeWrapper ow, vector<vector<thesis::spline>> spl_pairs){
            int m = A.cols(), n = msmt[0].rows();
            int ds = A.rows()/m - 1;
            int lt = ts.size();
            mat I = mat::Identity(m,m), bob, robert;
            vec temp;

            A.bottomRows(m) = a*I;
            vec du = A.colPivHouseholderQr().solve(P);
            double top = 0.0;
            /********************
                Construct top
            ********************/
            ow.setParameter(u0 + du);
            for(int i = 0; i < ds; i++)
            {
                ow.setControl(input[i]);
                ow.setPreviousIteration(spl_pairs[i]);
                robert = OdeIntWrapper(ow, y0, ts, 1E-6);
                bob = robert.bottomRows(n+n*m);
                temp = reshape(msmt[i] - bob.topRows(n), 1, n*lt).row(0).transpose();
                top += findO(ts, temp);

            }
            return top;
        };

        double FOld = F(alpha,A, ow, spl_pairs);
        double FNew;
        cout << "(" << pow(10,lower) << "," << FOld << ")"  << endl;

        for(double i=(lower+da); i<=upper; i+=da){

            FNew = F(pow(10,i),A, ow, spl_pairs);
            cout << "(" << pow(10,i) << "," << FNew << ")"  << endl;

            if(isnan(FOld)){
                alpha = pow(10,i);
                FOld = FNew;
            }
            else if(FNew < FOld){
                alpha = pow(10,i);
                FOld = FNew;
            }
        }

        //exit(0);
        return alpha;
	}

}
