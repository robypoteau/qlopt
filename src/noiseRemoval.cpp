#include <noiseRemoval.h>

namespace thesis {

    vec lsNoiseRemoval(const vec& y, double noise){
        size_t n = y.size();
        double alpha;

        mat I = mat::Identity(n,n);
        mat D = -I.topRows(n-1);
        for(int i=0; i<D.rows(); i++){
            D(i,i+1) = 1;
        }
        alpha = regs(y, noise);
        I = I + alpha*D.transpose()*D;
        vec x = I.colPivHouseholderQr().solve(y);
        return x;
    }

    double regs(const vec& y, double noise){
        size_t n = y.size();
        double tau = 1.0*noise;

        double divs = 240, upper = 6.0, lower = -6.0;
        double alpha = pow(10,lower), da = (upper - lower)/divs;

        mat I = mat::Identity(n,n);
        mat D = -I.topRows(n-1);
        for(int i=0; i<D.rows(); i++){
            D(i,i+1) = 1;
        }

        auto F = [I, D, y, tau](double alpha){
            mat A = I + alpha*D.transpose()*D;
            vec diff = A.colPivHouseholderQr().solve(y) - y;
            return diff.norm() - abs(tau)*y.norm();
        };

        double FOld = F(alpha);
        double FNew;

        for(double i=(lower+da); i<=upper; i+=da){

            FNew = F(pow(10,i));
            if(FNew < FOld){
                alpha = pow(10,i);
                FOld = FNew;
            }
        }
        cout << "Noise alpha = " << alpha << endl;
        return alpha;
    }
}
