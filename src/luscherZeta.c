#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <gsl/gsl_errno.h>
#include "zetaFunc.h"
#include "luscherZeta.h"

void luscherZeta(double * res, const double qsq, const int l, const int m, const double gamma, const double lambda, 
                 double * const dvec, const double tol, const int verbose, int * rstatus) {

  // this switches off the default GSL error handler, now we have to do it by ourselves!
  gsl_error_handler_t * handler = gsl_set_error_handler_off();

  complex double sum = firstPart(tol, l, m, dvec, gamma, lambda, qsq, verbose, &rstatus[0]) + secondPart( l, gamma, lambda, qsq, &rstatus[1]) + thirdPart(tol, l, m, dvec, gamma, lambda, qsq, verbose, &rstatus[2]);

  res[0] = creal(sum);
  res[1] = cimag(sum);

  gsl_set_error_handler (handler);
  return;
}
