#include <noiseRemoval.h>

namespace thesis {

    vec lsNoiseRemoval(const vec& y, double regParam){
        size_t n = y.size();

        mat I = mat::Identity(n,n);
        mat D = -regParam*I.topRows(n-1);
        for(int i=0; i<D.rows(); i++){
            D(i,i+1) = 1;
        }

        I = I + D.transpose()*D.transpose();
        vec x = I.fullPivHouseholderQr().solve(y);
        return x;
    }

}
