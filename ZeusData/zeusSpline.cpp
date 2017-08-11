#include "zeusSpline.h"
#include "zeusLogProviders.h"
#include "zeusExceptions.h"

namespace Zeus { 

    void Spline::prepareSpline ( const double x[], const double y[], 
            double *a, double *b, double *c, double *d, 
            double s1, double sn, long n ) {

        long n1, n2;
        long i, k;
        double r;
        double si, si1;
        double h, hi;

        if ( n < 4 ) {
            Zeus::LogProviders::FatalLog("Spline interpolation cannot be performed with %ld points\n",n);
            throw ( new Exceptions::NotEnoughPointsSplineException(n) );
        }

        n1 = n - 1;
        n2 = n - 2;

        /* auxiliary arrays h(=a) and delta(=d) */
        for ( i = 0 ; i < n1 ; i++ ) {
            if ( (x[i+1]-x[i]) < 1.0e-10 ) {
                Zeus::LogProviders::FatalLog("Spline x values not in increasing order\n");
                throw( new Exceptions::MonotonicSplineException() );
            }
            *(a+i) = x[i+1] - x[i];
            *(d+i) = (y[i+1]-y[i]) / *(a+i);
        }

        /* symmetric coefficient matrix (augmented) */
        for ( i = 0 ; i < n2 ; i++ ) {
            *(b+i) = 2.0 * (*(a+i) + *(a+i+1));
            k = n1-i-1;
            *(d+k) = 6.0 * (*(d+k) - *(d+k-1));
        }
        *(d+1) -= *(a)*s1;
        *(d+n1-1) -= *(a+n1-1)*sn;

        /* gauss solution of the tridiagonal system */
        for ( i = 1 ; i < n2 ; i++ ) {
            r = *(a+i) / *(b+i-1);
            *(b+i) -= r * *(a+i);
            *(d+i+1) -= r * *(d+i);
        }

        /* the sigma coefficients are stored in array d */
        *(d+n1-1) /= *(b+n2-1);
        for ( i = 1 ; i < n2 ; i++ ) {
            k = n1-i-1;
            *(d+k) = (*(d+k) - *(a+k)**(d+k+1)) / *(b+k-1);
        }
        *(d+n-1) = sn;

        /* spline coefficients */
        si1 = s1;
        for ( i = 0 ; i < n1 ; i++ ) {
            si = si1;
            si1 = *(d+i+1);
            h = *(a+i);
            hi = 1.0 / h;
            *(a+i) = hi/6.0 * (si * pow(x[i+1],3) - si1 * pow(x[i],3)) +
                hi * (y[i]*x[i+1] - y[i+1]*x[i]) + h/6.0 * (si1*x[i] - si*x[i+1]);
            *(b+i) = hi/2.0 * (si1 * pow(x[i],2) - si * pow(x[i+1],2)) +
                hi * (y[i+1] - y[i]) + h/6.0 * (si - si1);
            *(c+i) = hi/2.0 * (si*x[i+1] - si1*x[i]);
            *(d+i) = hi/6.0 * (si1 - si);
        }

    } // Spline

} // namespaces
