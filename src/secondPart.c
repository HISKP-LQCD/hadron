#include "zetaFunc.h"

#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dawson.h>

#include <complex.h>
#include <math.h>

/***************************************************
 * This program calculates the second part of the
 * zeta function:
 * 		delta_{L,0} *Y_00*gamma*PI^{3/2}
 * 		* [2*q^2* \int_0^Lamda exp{t*q^2}/sqrt{t} dt
 * 		   - 2*exp{Lamda*q^2}/sqrt{Lamda}
 ***************************************************/

double complex secondPart(const int l,
                          const double gamma,
                          const double Lambda,
                          const double qsq,
                          int *const rstatus,
                          const int verbose) {
  int s1 = 0, s2 = 0;
  double complex secondPartInt = 0 + 0 * I;
  if (l == 0) {
    // Arithmetic method.
    double inteCore = 0.0;
    if (qsq > 0) {
      // Dawson function for qsq > 0
      inteCore = 4 * sqrt(qsq) * exp(Lambda * qsq) * gsl_sf_dawson(sqrt(Lambda * qsq));
    } else if (fabs(qsq) < DBL_EPSILON) {
      inteCore = 0;
    } else if (qsq < 0) {
      // Error function for qsq < 0
      inteCore = -2 * sqrt(M_PI) * sqrt(-qsq) * erf(sqrt(-Lambda * qsq));
    }

    secondPartInt = spheHarm(0, 0, 0, 0, &s1) * gamma * pow(M_PI, 3.0 / 2.0) *
                    (inteCore - 2 * exp(Lambda * qsq) / sqrt(Lambda));
  }
  *rstatus = s1 + s2;
  if (verbose) {
    Rprintf("Second term = (%e, %e)\n", creal(secondPartInt), cimag(secondPartInt));
  }
  return secondPartInt;
}
