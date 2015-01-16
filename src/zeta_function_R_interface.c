#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "luscherZeta.h"
void pmode_free_arrays();

SEXP LuscherZetaArray(SEXP qsq_, SEXP n_, SEXP l_, SEXP m_, SEXP dvec_, SEXP gamma_, SEXP lambda_, SEXP tol_, SEXP verbose_) {

  PROTECT(qsq_ = AS_NUMERIC(qsq_));
  PROTECT(l_ = AS_INTEGER(l_));
  PROTECT(m_ = AS_INTEGER(m_));
  PROTECT(gamma_ = AS_NUMERIC(gamma_));
  PROTECT(lambda_ = AS_NUMERIC(lambda_));
  PROTECT(dvec_ = AS_NUMERIC(dvec_));
  PROTECT(n_ = AS_INTEGER(n_));
  PROTECT(tol_ = AS_NUMERIC(tol_));
  PROTECT(verbose_ = AS_INTEGER(verbose_));

  double * qsq = NUMERIC_POINTER(qsq_);
  const int l = INTEGER_POINTER(l_)[0];
  const int m = INTEGER_POINTER(m_)[0];
  const double gamma = NUMERIC_POINTER(gamma_)[0];
  const double lambda = NUMERIC_POINTER(lambda_)[0];
  const int n = INTEGER_POINTER(n_)[0];
  const double tol = NUMERIC_POINTER(tol_)[0];
  double * dvec = NUMERIC_POINTER(dvec_);
  const int verbose = INTEGER_POINTER(verbose_)[0];
  SEXP res;
  Rcomplex * resp;
  double ires[2];
  PROTECT(res = NEW_COMPLEX(n));
  resp = COMPLEX_POINTER(res);

  for(int i = 0; i < n; i ++) {
    int rstatus[3] = {0,0,0};

    luscherZeta(ires, qsq[i], l, m, gamma, lambda, dvec, tol, verbose, rstatus);

    resp[i].r = ires[0];
    resp[i].i = ires[1];

    if(rstatus[0]!=0) {
      warning("GSL error in evaluation of first contribution to Luescher Zeta function, status code %d for element %d\n", rstatus[0], i+1);
      resp[i].r = NA_REAL;
      resp[i].i = NA_REAL;
    }
    else if(rstatus[1]!=0) {
      warning("GSL error in evaluation of second contribution to Luescher Zeta function, status code %d for element %d\n", rstatus[1], i+1);
      resp[i].r = NA_REAL;
      resp[i].i = NA_REAL;
    }
    else if(rstatus[2]!=0) {
      warning("GSL error in evaluation of third contribution to Luescher Zeta function, status code %d for element %d\n", rstatus[2], i+1);
      resp[i].r = NA_REAL;
      resp[i].i = NA_REAL;
    }
  }

  UNPROTECT(10);
  pmode_free_arrays();
  return(res);

}
