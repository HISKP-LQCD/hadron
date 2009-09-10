#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "alpha_s.h"

const double pi=3.1415926535897932384626433832;

void getmpssqpion(double * res, double * r0sqTwoBmu, double * par, const int N, 
		  const int npar, const int fitnnlo, const int fitkmf,
		  const int fitasq) {


  double xi, rln1, rln2, rln3;
  double asq = 1., fitk = 0.;
  int i, np = 3*N+6;

  if(fitasq < 0) {
    asq = 0;
    np = npar-1;
  }
  if(fitkmf) {
    fitk = 1.;
  }
  if(fitnnlo) {
    for(i = 0; i < N; i++) {
      rln1 = log(r0sqTwoBmu[i]/par[3*N+4]/par[3*N+4]);
      rln2 = log(r0sqTwoBmu[i]/par[3*N+5]/par[3*N+5]);
      rln3 = log(r0sqTwoBmu[i]/par[0]/par[0]);
      xi = r0sqTwoBmu[i]/(4.0*pi*par[2])/(4.0*pi*par[2]);
      res[i] = r0sqTwoBmu[i]*(1. + xi*rln3 + asq/par[3+fitasq]/par[3+fitasq]*par[npar-2] +
			      17./2.*xi*xi*(((-28.*rln1 -32.*rln2 + 9.*rln3 + 49.)/51.)*((-28.*rln1 -32.*rln2 + 9.*rln3 + 49.)/51.)
					    +fitk*8./17.*par[np]
					    )
			      );
    }
  }
  else {
    for(i = 0; i < N; i++) {
      rln3 = log(r0sqTwoBmu[i]/par[0]/par[0]);
      res[i] = r0sqTwoBmu[i]*(1.0 + asq/par[3+fitasq]/par[3+fitasq]*par[npar-2] + r0sqTwoBmu[i]*(rln3/(4.0*pi*par[2])/(4.0*pi*par[2]) ));
    }
  }
  return;
}

void getfpspion(double * res, double * r0sqTwoBmu, double * par, const int N, 
		const int npar, const int fitnnlo, const int fitkmf,
		const int fitasq) {

  double xi, rln1, rln2, rln3, rln4;
  double asq = 1., fitk = 0.;
  int i, np = 3*N+6;

  if(fitasq < 0) {
    asq = 0;
    np = npar-1;
  }
  if(fitkmf) {
    fitk = 1.;
  }
  if(fitnnlo) {
    for(i = 0; i < N; i++) {
      rln1 = log(r0sqTwoBmu[i]/par[3*N+4]/par[3*N+4]);
      rln2 = log(r0sqTwoBmu[i]/par[3*N+5]/par[3*N+5]);
      rln3 = log(r0sqTwoBmu[i]/par[0]/par[0]);
      rln4 = log(r0sqTwoBmu[i]/par[1]/par[1]);
      xi = r0sqTwoBmu[i]/(4.0*pi*par[2])/(4.0*pi*par[2]);
      res[i] = par[2]*(1. - 2.*xi*rln4 +  asq/par[3+fitasq]/par[3+fitasq]*par[npar-1] -
		       5*xi*xi*(((-14*rln1 - 16*rln2 - 6*rln3 + 6*rln4 + 23)/30.)*((-14*rln1 - 16*rln2 - 6*rln3 + 6*rln4 + 23)/30.)
				+ fitk*4./5.*par[np]
				)
		       );
    }
  }
  else {
    for(i = 0; i < N; i++) {
      rln3 = log(r0sqTwoBmu[i]/par[0]/par[0]);
      rln4 = log(r0sqTwoBmu[i]/par[1]/par[1]);
      res[i] = par[2]*(1.0 + asq/par[3+fitasq]/par[3+fitasq]*par[npar-1]  
		       - 2.0*r0sqTwoBmu[i]*(rln4)/(4.0*pi*par[2])/(4.0*pi*par[2]) );
    }
  }
  return;  
}


SEXP chisqrbody(SEXP data_, SEXP r0exp_, SEXP N_, SEXP par_, SEXP i_, SEXP ij_) {

  double * r0a, * dr0, * resp, * par, * mu, r0exp, chisum = 0.;
  int i, j, N, *ij;

  PROTECT(data_ = AS_VECTOR(data_));
  PROTECT(N_ = AS_INTEGER(N_));
  PROTECT(par_ = AS_NUMERIC(par_));
  PROTECT(r0exp_ = AS_NUMERIC(r0exp_));
  PROTECT(i_ = AS_INTEGER(i_));
  PROTECT(ij_ = AS_INTEGER(ij_));

  N = INTEGER_POINTER(N_)[0];
  par = NUMERIC_POINTER(par_);
  r0exp = NUMERIC_POINTER(r0exp_)[0];
  mu = NUMERIC_POINTER(VECTOR_ELT(data_, 0));
  r0a = NUMERIC_POINTER(VECTOR_ELT(data_, 12));
  dr0 = NUMERIC_POINTER(VECTOR_ELT(data_, 13));
  i = INTEGER_POINTER(i_)[0];
  ij = INTEGER_POINTER(ij_);

  for(j = 0; j < N; j++) {
    if(!ISNA(r0a[i])) {
      chisum += pow((par[3+i] + par[3+N+i]*pow(mu[j], r0exp) - r0a[i])/dr0[i], 2.);
    }
  }

  return ScalarReal(chisum);
}
