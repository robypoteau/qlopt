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
}
