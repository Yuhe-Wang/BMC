#pragma once

#include "zeusConfig.h"

namespace Zeus { 
    /**
    * \brief Namespace for the subroutines related to spline fitting.
    */
    class ZEUS_EXPORT Spline {

    public:

        /**
         \brief Subroutine to perform interpolation.

        The interpolating polynomial in the i-th interval, from
        x[i] to x[i+1], is pi[x]=a[i]+x*(b[i]+x*(c[i]+x*d[i])).

        Reference:
        M.J. Maron, 'Numerical Analysis: A Practical Approach',
        MacMillan Publ. Co., New York 1982.

        \param x Grid points in increasing order.
        \param y Corresponding function values.
        \param s1 Second derivative at x[0].
        \param sn Second derivative at x[n-1].
        \param a Array of coefficient a in the interpolating polynomial, one per grid point.
        \param b Array of coefficient b in the interpolating polynomial, one per grid point.
        \param c Array of coefficient c in the interpolating polynomial, one per grid point.
        \param d Array of coefficient d in the interpolating polynomial, one per grid point.
        \param n Number of grid points.
        */
        static void prepareSpline ( const double x[], const double y[], 
            double *a, double *b, double *c, double *d, 
            double s1, double sn, long n );
        /**
        * \brief Performs a spline interpolation using the coefficients previously prepared with prepareSpline()
        */
        static double interpolate(double s, int nbin, const double x[], const double a[], const double b[], const double c[], const double d[]) {
            int i = findBin(s,nbin,x);
            double s2 = s*s;
            return a[i] + s*b[i] + s2*c[i] + s*s2*d[i];
        };
        /**
        * \brief Returns the bin index in the array \a x[] in which the input variable \a s falls. Assumes x[] is in increasing order and uses a binary search
        */
        static int findBin(double s, int nbin, const double x[]) {
            int bin;
            if( s <= x[0] ) bin = 0;
            else if( s >= x[nbin-1] ) bin = nbin-1;
            else {
                int low = 0, high = nbin;
                while( high - low > 1 ) {
                    int av = (low + high)/2;
                    if( s < x[av] ) high = av; else low = av;
                }
                bin = high - 1;
            }
            return bin;
        };

    };

} // namespaces

