#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "alpha_s.h"

const double pi=3.1415926535897932384626433832;

static R_INLINE void getboth(double * res, double * r0sqTwoBmu, double * par, const int N, 
			     const int npar, const int dl, const int fitnnlo, const int fitkmf,
			     const int fitasq) {
  
  double xi, rln1, rln2, rln3, rln4;
  double asq = 1., fitk = 0.;
  int i, np = 2*N+8, np2 = 2*N+9;

  if(fitasq < 0) {
    asq = 0;
    np  = npar-1;
    np2 = np;
  }
  if(fitkmf) {
    fitk = 1.;
  }
  if(fitnnlo) {
    for(i = 0; i < dl; i++) {
      rln1 = log(r0sqTwoBmu[i]/par[2*N+6]/par[2*N+6]);
      rln2 = log(r0sqTwoBmu[i]/par[3*N+5]/par[3*N+5]);
      rln3 = log(r0sqTwoBmu[i]/par[0]/par[0]);
      rln4 = log(r0sqTwoBmu[i]/par[1]/par[1]);
      xi = r0sqTwoBmu[i]/(4.0*pi*par[2])/(4.0*pi*par[2]);
      res[i] = par[2]*(1. - 2.*xi*rln4 +  asq/par[3+fitasq]/par[3+fitasq]*par[npar-1] -
		       5*xi*xi*(((-14*rln1 - 16*rln2 - 6*rln3 + 6*rln4 + 23)/30.)*((-14*rln1 - 16*rln2 - 6*rln3 + 6*rln4 + 23)/30.)
				+ fitk*4./5.*par[np2]
				)
		       );
      res[i+dl] = sqrt(r0sqTwoBmu[i]*(1. + xi*rln3 + asq/par[3+fitasq]/par[3+fitasq]*par[npar-2] +
				      17./2.*xi*xi*(((-28.*rln1 -32.*rln2 + 9.*rln3 + 49.)/51.)*((-28.*rln1 -32.*rln2 + 9.*rln3 + 49.)/51.)
						    +fitk*8./17.*par[np]
						    )
				      )
		       );
    }
  }
  else {
    for(i = 0; i < dl; i++) {
      rln3 = log(r0sqTwoBmu[i]/par[0]/par[0]);
      rln4 = log(r0sqTwoBmu[i]/par[1]/par[1]);
      res[i] = par[2]*(1.0 + asq/par[3+fitasq]/par[3+fitasq]*par[npar-1]  
		       - 2.0*r0sqTwoBmu[i]*(rln4)/(4.0*pi*par[2])/(4.0*pi*par[2]) );
      res[i+dl] = sqrt(r0sqTwoBmu[i]*(1.0 + asq/par[3+fitasq]/par[3+fitasq]*par[npar-2] + r0sqTwoBmu[i]*(rln3/(4.0*pi*par[2])/(4.0*pi*par[2]) )));
    }
  }
  return;

}

void getmpssqpion(double * res, double * r0sqTwoBmu, double * par, const int N, 
		  const int npar, const int dl, const int fitnnlo, const int fitkmf,
		  const int fitasq) {


  double xi, rln1, rln2, rln3;
  double asq = 1., fitk = 0.;
  int i, np = 2*N+8;

  if(fitasq < 0) {
    asq = 0;
    np = npar-1;
  }
  if(fitkmf) {
    fitk = 1.;
  }
  if(fitnnlo) {
    for(i = 0; i < dl; i++) {
      rln1 = log(r0sqTwoBmu[i]/par[2*N+6]/par[2*N+6]);
      rln2 = log(r0sqTwoBmu[i]/par[2*N+7]/par[2*N+7]);
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
    for(i = 0; i < dl; i++) {
      rln3 = log(r0sqTwoBmu[i]/par[0]/par[0]);
      res[i] = r0sqTwoBmu[i]*(1.0 + asq/par[3+fitasq]/par[3+fitasq]*par[npar-2] + r0sqTwoBmu[i]*(rln3/(4.0*pi*par[2])/(4.0*pi*par[2]) ));
    }
  }
  return;
}

void getfpspion(double * res, double * r0sqTwoBmu, double * par, const int N, 
		const int npar, const int dl, const int fitnnlo, const int fitkmf,
		const int fitasq) {

  double xi, rln1, rln2, rln3, rln4;
  double asq = 1., fitk = 0.;
  int i, np = 2*N+9;

  if(fitasq < 0) {
    asq = 0;
    np = npar-1;
  }
  if(fitkmf) {
    fitk = 1.;
  }
  if(fitnnlo) {
    for(i = 0; i < dl; i++) {
      rln1 = log(r0sqTwoBmu[i]/par[2*N+6]/par[2*N+6]);
      rln2 = log(r0sqTwoBmu[i]/par[2*N+7]/par[2*N+7]);
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

SEXP getbothc(SEXP r0sqTwoBmu_, SEXP par_, SEXP N_, SEXP fitnnlo_, SEXP fitkmf_, SEXP fitasq_) {
  double *r0sqTwoBmu, *par, *res, fitasq;
  int N, fitnnlo, fitkmf, npar, dl;
  SEXP res_;
  
  PROTECT(r0sqTwoBmu_ = AS_NUMERIC(r0sqTwoBmu_));
  PROTECT(par_ = AS_NUMERIC(par_));
  PROTECT(fitasq_ = AS_NUMERIC(fitasq_));
  PROTECT(fitnnlo_ = AS_INTEGER(fitnnlo_));
  PROTECT(fitkmf_ = AS_INTEGER(fitkmf_));
  PROTECT(N_ = AS_INTEGER(N_));
  
  r0sqTwoBmu = NUMERIC_POINTER(r0sqTwoBmu_);
  par = NUMERIC_POINTER(par_);
  fitasq = NUMERIC_POINTER(fitasq_)[0];
  fitnnlo = INTEGER_POINTER(fitnnlo_)[0];
  fitkmf = INTEGER_POINTER(fitkmf_)[0];
  N = INTEGER_POINTER(N_)[0];
  npar = LENGTH(par_);
  dl = LENGTH(r0sqTwoBmu_);

  PROTECT(res_ = NEW_NUMERIC(2*dl));
  res = NUMERIC_POINTER(res_);
  getboth(res, r0sqTwoBmu, par, N, npar, dl, fitnnlo, fitkmf, (int)fitasq);  
  UNPROTECT(7);
  return(res_);
}


SEXP chisqrbody(SEXP data_, SEXP r0exp_, SEXP N_, SEXP par_, SEXP i_, SEXP ij_,
		SEXP fitnnlo_, SEXP fitkmf_, SEXP fitasqr_) {

  double * r0a, * dr0, * par, * mu, r0exp, chisum = 0.;
  int i, j, N, *ij, npar, fitnnlo, fitkmf, fitasqr;
  double *r0sqTwoBmu, *mpssq, *fpsV;

  PROTECT(data_ = AS_VECTOR(data_));
  PROTECT(N_ = AS_INTEGER(N_));
  PROTECT(par_ = AS_NUMERIC(par_));
  PROTECT(r0exp_ = AS_NUMERIC(r0exp_));
  PROTECT(i_ = AS_INTEGER(i_));
  PROTECT(ij_ = AS_INTEGER(ij_));
  PROTECT(fitnnlo_ = AS_INTEGER(fitnnlo_));
  PROTECT(fitkmf_ = AS_INTEGER(fitkmf_));
  PROTECT(fitasqr_ = AS_INTEGER(fitasqr_));

  N = INTEGER_POINTER(N_)[0];
  par = NUMERIC_POINTER(par_);
  r0exp = NUMERIC_POINTER(r0exp_)[0];
  mu = REAL(VECTOR_ELT(data_, 0));
  r0a = REAL(VECTOR_ELT(data_, 12));
  dr0 = REAL(VECTOR_ELT(data_, 13));
  i = INTEGER_POINTER(i_)[0];
  ij = INTEGER_POINTER(ij_);
  npar = LENGTH(par_);
  fitasqr = INTEGER_POINTER(fitasqr_)[0];
  fitnnlo = INTEGER_POINTER(fitnnlo_)[0];
  fitkmf = INTEGER_POINTER(fitkmf_)[0];

  r0sqTwoBmu = (double*) Calloc(N, double);
  mpssq = (double*) Calloc(N, double);
  fpsV = (double*) Calloc(N, double);

  for(j = 0; j < N; j++) {
    if(!ISNA(r0a[j]) && !ISNAN(r0a[j])) {
      printf("%e %e %e\n", r0a[j], dr0[j], mu[j]);
      chisum += pow((par[3+i] + par[3+N+i]*pow(mu[j], r0exp) - r0a[j])/dr0[j], 2.);
    }
    //      r0sqTwoBmu <- r0TwoB*data[[i]]$mu[ij]*par[4+i]/par[4+2*N+i]
    r0sqTwoBmu[j] = par[3]*mu[j]*par[3+i]/par[3+2+N+i];
  }
  //    mpssq <- getmpssq.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.a)
  //    fpsV <- getfps.pion(r0sqTwoBmu, par, N, fit.nnlo=fit.nnlo, fit.kmf=fit.kmf, fit.asq=fit.a)
  //getmpssqpion(mpssq, r0sqTwoBmu, par, N, npar, fitnnlo, fitkmf, fitasqr);
  //getfpspion(fpsV, r0sqTwoBmu, par, N, npar, fitnnlo, fitkmf, fitasqr);

  

  UNPROTECT(9);
  Free(r0sqTwoBmu);
  Free(mpssq);
  Free(fpsV);
  return ScalarReal(chisum);
}

