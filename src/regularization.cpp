#include <regularization.h>

namespace thesis{

    double findAlpha(mat A, vec P, double Om, double Nd){
        int m = A.cols();
        int ds = A.rows()/m - 1;

        vec sv = A.topRows(ds*m).jacobiSvd().singularValues();

        double divs = 24, upper = 6.0, lower = -6.0;
        double alpha = pow(10,lower), da = (upper - lower)/divs;

        mat I;
        vec du(m);

        I = mat::Identity(m,m);
        A.bottomRows(m) = pow(10,lower)*I;
        du = A.colPivHouseholderQr().solve(P);

        auto F = [sv, P, Om, I, Nd](double a, mat  A){
            int m = A.cols();
            int ds = A.rows()/m - 1;
            vec temp = sv;
            temp.array() += a;
            temp = sv.array().square() / temp.array().square();
            A.bottomRows(m) = a*I;
            vec du = A.colPivHouseholderQr().solve(P);
            double top = Om;
            for(int i = 0; i < ds; i++)
            {
                top += du.transpose()*(A.middleRows(i*m, m)*du
                    - 2*P.segment(i*m, m));
            }
            return top/(Nd - temp.sum());
        };

        double FOld = F(alpha,A);
        double FNew;
        cout << "(" << pow(10,lower) << "," << FOld << ")"  << endl;

        for(double i=(lower+da); i<=upper; i+=da){

            FNew = F(pow(10,i),A);
            cout << "(" << pow(10,i) << "," << FNew << ")"  << endl;
            if(FNew < FOld){
                alpha = pow(10,i);
                FOld = FNew;
            }
        }

        //exit(0);
        return alpha;
    }

    double gcv(mat A, vec P, vec u0, OdeWrapper ow, vector<mat>& msmt, const vector<vec>& input, vec y0, vec ts, vector<vector<thesis::spline>> spl_pairs){
        int m = A.cols();
        int ds = A.rows()/m - 1;
        vec sv = A.topRows(ds*m).jacobiSvd().singularValues();

        double divs = 12, upper = 3.0, lower = -3.0;
        double alpha = pow(10,lower), da = (upper - lower)/divs;

        vec du(m);

        A.bottomRows(m) = pow(10,lower)*mat::Identity(m,m);
        du = A.colPivHouseholderQr().solve(P);

        auto F = [sv, u0, msmt, input, P, ts, y0](double a, mat  A, OdeWrapper ow, vector<vector<thesis::spline>> spl_pairs){
            int m = A.cols(), n = msmt[0].rows();
            int ds = A.rows()/m - 1;
            int lt = ts.size();
            mat I = mat::Identity(m,m), bob, robert;
            vec temp = sv, temp2;
            temp.array() += a;
            temp = sv.array().square() / temp.array().square();

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
                robert = OdeIntWrapper(ow, y0, ts);
                bob = robert.bottomRows(n+n*m);
                temp2 = reshape(msmt[i] - bob.topRows(n), 1, n*lt).row(0).transpose();
                top += findO(ts, temp2);

            }
            return top/(ds*lt - temp.sum());
        };

        double FOld = F(alpha,A, ow, spl_pairs);
        double FNew;
        cout << "(" << pow(10,lower) << "," << FOld << ")"  << endl;

        for(double i=(lower+da); i<=upper; i+=da){

            FNew = F(pow(10,i),A, ow, spl_pairs);
            cout << "(" << pow(10,i) << "," << FNew << ")"  << endl;
            if(FNew < FOld){
                alpha = pow(10,i);
                FOld = FNew;
            }
        }

        //exit(0);
        return alpha;
	}

}
