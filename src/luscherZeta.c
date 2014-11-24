#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "zetaFunc.h"
#include "luscherZeta.h"

void luscherZeta(double * res, const double qsq, const int l, const int m, const double gamma, const double lambda, int * dvec) {
  const int N = 3;

  complex double sum = firstPart(N, l, m, dvec, gamma, lambda, qsq) + secondPart( l, gamma, lambda, qsq) + thirdPart(N, l, m, dvec, gamma, lambda, qsq);

  res[0] = creal(sum);
  res[1] = cimag(sum);
  
  return;
}
