#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include <complex.h>
#include <math.h>

const double N = 16. * M_PI * M_PI;
const double fmGeV = 0.1973269631;

double w(const double x, const double lambda) {
  return (exp(-sqrt(1. + x * x) * lambda));
}

double complex Jb(const double complex x) {
  double complex sigma;

  if (cabs(x) < 1.e-10) {
    sigma = csqrt(1. - 4. / x);
    return ((sigma * clog((sigma - 1.) / (sigma + 1.)) + 2.) / N);
  }
  return (0. * I);
}

double complex Jb1(const double complex x) {
  double complex sigma;

  if (cabs(x) < 1.e-10) {
    sigma = csqrt(1. - 4. / x);
    return ((2. / x / sigma * clog((sigma - 1.) / (sigma + 1.)) - 1.) / N / x);
  }
  return (1. / 6. / N + 0. * I);
}

typedef struct {
  int n;
  int k;
  double lambda;
} pars;

double x1(double x, void *params_) {
  double complex K;
  double complex z = 2 * (1. + I * x);
  pars *params = (pars *)params_;

  if (params->n == 0) {
    K = Jb(z);
  } else {
    K = Jb1(z);
  }

  if (params->k % 2 == 0) {
    return (w(x, params->lambda) * pow(x, params->k) * cimag(K));
  }
  return (w(x, params->lambda) * pow(x, params->k) * creal(K));
}

void calc_R_int(double *R[], const double lambda) {
  double result, error;
  int i, j;
  gsl_function F;
  gsl_integration_workspace *w;
  pars params;
  const int intervals = 1000000;

  w = gsl_integration_workspace_alloc(intervals);
  F.function = &x1;
  F.params = (void *)&params;
  params.lambda = lambda;

  for (i = 0; i < 2; i++) {
    params.n = i;
    for (j = 0; j < 3; j++) {
      params.k = j;
      gsl_integration_qagi(&F, 0., 1.e-6, intervals, w, &result, &error);
      R[i][j] = N * result;
    }
  }

  gsl_integration_workspace_free(w);
}
