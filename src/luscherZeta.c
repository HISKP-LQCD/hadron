#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <gsl/gsl_errno.h>
#include "zetaFunc.h"
#include "luscherZeta.h"

void luscherZeta(double * res, const double qsq, const int l, const int m, const double gamma, const double lambda, int * dvec, const int N, int * rstatus) {

  gsl_error_handler_t * handler = gsl_set_error_handler_off();

  complex double sum = firstPart(N, l, m, dvec, gamma, lambda, qsq, &rstatus[0]) + secondPart( l, gamma, lambda, qsq, &rstatus[1]) + thirdPart(N, l, m, dvec, gamma, lambda, qsq, &rstatus[2]);

  res[0] = creal(sum);
  res[1] = cimag(sum);

  gsl_set_error_handler (handler);
  return;
}
