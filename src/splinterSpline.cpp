#include <splinterSpline.h>

namespace thesis{
    vec splinterSpline(vec x, vec y, double alpha)
    {
        // Create new DataTable to manage samples
        DataTable samples;

        // Sample the function
        DenseVector s(1);
        vec ans(x.size());

        for(int j = 0; j < x.size(); j++)
        {
            // Sample function at x
            s(0) = x(j);

            // Store sample
            samples.addSample(s,y(j));
        }
        // Build penalized B-spline (P-spline) that smooths the samples
        BSpline pspline = BSpline::Builder(samples)
                .degree(2)
                .smoothing(BSpline::Smoothing::PSPLINE)
                .alpha(alpha)
                .build();

        for(int j = 0; j < x.size(); j++)
        {
            s(0) = x(j);
            ans(j) = pspline.eval(s);
        }
        return ans;
    }

}
