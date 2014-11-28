#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "luscherZeta.h"

SEXP LuscherZeta(SEXP qsq_, SEXP l_, SEXP m_, SEXP dvec_, SEXP gamma_, SEXP lambda_, SEXP N_) {

  PROTECT(qsq_ = AS_NUMERIC(qsq_));
  PROTECT(l_ = AS_INTEGER(l_));
  PROTECT(m_ = AS_INTEGER(m_));
  PROTECT(gamma_ = AS_NUMERIC(gamma_));
  PROTECT(lambda_ = AS_NUMERIC(lambda_));
  PROTECT(dvec_ = AS_NUMERIC(dvec_));
  PROTECT(N_ = AS_INTEGER(N_));

  const double qsq = NUMERIC_POINTER(qsq_)[0];
  const int l = INTEGER_POINTER(l_)[0];
  const int m = INTEGER_POINTER(m_)[0];
  const double gamma = NUMERIC_POINTER(gamma_)[0];
  const double lambda = NUMERIC_POINTER(lambda_)[0];
  const int N = INTEGER_POINTER(N_)[0];
  double * dvec = NUMERIC_POINTER(dvec_);
  
  SEXP res;
  Rcomplex * resp;
  double ires[2] = {0.,0.};
  int rstatus[3];
  luscherZeta(ires, qsq, l, m, gamma, lambda, dvec, N, rstatus);

  if(rstatus[0]!=0) {
    error("GSL error in evaluation of first contribution to Luescher Zeta function, status code %d\n", rstatus[0]);
  }
  if(rstatus[1]!=0) {
    error("GSL error in evaluation of second contribution to Luescher Zeta function, status code %d\n", rstatus[1]);
  }
  if(rstatus[2]!=0) {
    error("GSL error in evaluation of third contribution to Luescher Zeta function, status code %d\n", rstatus[2]);
  }

  PROTECT(res = NEW_COMPLEX(1));
  resp = COMPLEX_POINTER(res);
  resp[0].r = ires[0];
  resp[0].i = ires[1];

  UNPROTECT(8);
  return(res);
}

SEXP LuscherZetaArray(SEXP qsq_, SEXP n_, SEXP l_, SEXP m_, SEXP dvec_, SEXP gamma_, SEXP lambda_, SEXP N_) {

  PROTECT(qsq_ = AS_NUMERIC(qsq_));
  PROTECT(l_ = AS_INTEGER(l_));
  PROTECT(m_ = AS_INTEGER(m_));
  PROTECT(gamma_ = AS_NUMERIC(gamma_));
  PROTECT(lambda_ = AS_NUMERIC(lambda_));
  PROTECT(dvec_ = AS_NUMERIC(dvec_));
  PROTECT(n_ = AS_INTEGER(n_));
  PROTECT(N_ = AS_INTEGER(N_));

  double * qsq = NUMERIC_POINTER(qsq_);
  const int l = INTEGER_POINTER(l_)[0];
  const int m = INTEGER_POINTER(m_)[0];
  const double gamma = NUMERIC_POINTER(gamma_)[0];
  const double lambda = NUMERIC_POINTER(lambda_)[0];
  const int n = INTEGER_POINTER(n_)[0];
  const int N = INTEGER_POINTER(N_)[0];
  double * dvec = NUMERIC_POINTER(dvec_);
  int rstatus[3];
  SEXP res;
  Rcomplex * resp;
  double ires[2];
  PROTECT(res = NEW_COMPLEX(n));
  resp = COMPLEX_POINTER(res);

  for(int i = 0; i < n; i ++) {
    luscherZeta(ires, qsq[i], l, m, gamma, lambda, dvec, N, rstatus);
    resp[i].r = ires[0];
    resp[i].i = ires[1];
    if(rstatus[0]!=0) {
      error("GSL error in evaluation of first contribution to Luescher Zeta function, status code %d\n", rstatus[0]);
      break;
    }
    if(rstatus[1]!=0) {
      error("GSL error in evaluation of second contribution to Luescher Zeta function, status code %d\n", rstatus[1]);
      break;
    }
    if(rstatus[2]!=0) {
      error("GSL error in evaluation of third contribution to Luescher Zeta function, status code %d\n", rstatus[2]);
      break;
    }

  }


  UNPROTECT(9);
  return(res);

}
