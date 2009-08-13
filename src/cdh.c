#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include "cdh.h"

double g1(double x) {
  double weights[20] = {6.,12.,8.,6.,24.,24.,0.,12.,30.,24.,24.,8.,24.,48.,0.,6.,48.,36.,24.,24.};
  double sex, res=0.;
  int i;

  for(i = 0; i < 20; i++) {
    sex = x*sqrt((double)(i+1.));
    res += 4.*weights[i]*gsl_sf_bessel_Kn(1, sex)/sex;
  }
  return(res);
}

int g1array(double * x, double * res, const int n) {
  double weights[20] = {6.,12.,8.,6.,24.,24.,0.,12.,30.,24.,24.,8.,24.,48.,0.,6.,48.,36.,24.,24.};
  double ex[20], sex[20];
  int i, j;

  for(i = 0; i < 20; i++) {
    ex[i] = sqrt((double)(i+1.));
  }
  for(j = 0; i < n; i++) {
    res[j] = 0.;
    for(i = 0; i < 20; i++) {
      sex[i] = x[j]*ex[i];
      res[j] += 4.*weights[i]*gsl_sf_bessel_Kn(1, sex[i])/sex[i];
    }
  }
  return(0);
}


void fscdh(double rev, double aLamb1, double aLamb2, double aLamb3, double aLamb4,
	   double aF0_, double a_fm, int * L, double * ampiV, double *afpiV, const int n,
	   double * mpiFV, double * fpiFV, const int printit, 
	   double * rtilde, const int incim6) {
  
  const double pi=3.1415926535897932384626433832;
  int i, j;
  double N, amrho_phys, z, Rfpi, Rmpi, lambda_pi;
  double gg[4];
  double mm[] = {6, 12, 8, 6, 24, 24, 0, 12, 30, 24, 24, 8, 24, 48, 0, 6, 48, 36, 24, 24};
  const int mm1 = 20;
  double *lb1, *lb2, *lb3, *lb4, *lpi, *M_P, *F_P, *xi_P, *mmB0, *mmB2;
  double *S4mpi, *S4fpi, *I4mpi, *I2mpi, *I2fpi, *I4fpi, *I6mpi;
  double aF0;

  aF0 = aF0_/sqrt(2.);

  gg[0] = 2. - pi/2.;
  gg[1] = pi/4. - 0.5;
  gg[2] = 0.5 - pi/8.;
  gg[3] = 3*pi/16. - 0.5;
  N = 16*pi*pi;
  amrho_phys = a_fm * 770.0/197.3;
  lb1 = (double*) malloc(n*sizeof(double));
  lb2 = (double*) malloc(n*sizeof(double));
  lb3 = (double*) malloc(n*sizeof(double));
  lb4 = (double*) malloc(n*sizeof(double));
  lpi = (double*) malloc(n*sizeof(double));
  for(i = 0; i < n; i++) {
    lb1[i] = log(aLamb1/ampiV[i]*aLamb1/ampiV[i]);
    lb2[i] = log(aLamb2/ampiV[i]*aLamb2/ampiV[i]);
    lb3[i] = log(aLamb3/ampiV[i]*aLamb3/ampiV[i]);
    lb4[i] = log(aLamb4/ampiV[i]*aLamb4/ampiV[i]);
    lpi[i] = log(ampiV[i]/amrho_phys*ampiV[i]/amrho_phys);
  }

  M_P = (double*)   malloc(n*sizeof(double));
  F_P = (double*)   malloc(n*sizeof(double));
  mmB0 = (double*)  malloc(n*sizeof(double));
  mmB2 = (double*)  malloc(n*sizeof(double));
  xi_P = (double*)  malloc(n*sizeof(double));
  S4mpi = (double*) malloc(n*sizeof(double));
  S4fpi = (double*) malloc(n*sizeof(double));
  I2mpi = (double*) malloc(n*sizeof(double));
  I4mpi = (double*) malloc(n*sizeof(double));
  I2fpi = (double*) malloc(n*sizeof(double));
  I4fpi = (double*) malloc(n*sizeof(double));
  I6mpi = (double*) malloc(n*sizeof(double));
  for(i = 0; i < n; i++) {
    M_P[i] = ampiV[i];
    F_P[i] = afpiV[i];
    xi_P[i] = (M_P[i]/(4*pi*aF0))*(M_P[i]/(4*pi*aF0)); 
    mmB0[i] = 0.;
    mmB2[i] = 0.;
  }
  for(j = 0; j < n; j++) {
    lambda_pi = ampiV[j]*L[j];
    for(i = 0; i < mm1; i++) {
      z = sqrt(1.+i)*lambda_pi;
      mmB0[j] += mm[i]*2.*gsl_sf_bessel_K1(z)/z;
      mmB2[j] += mm[i]*2.*gsl_sf_bessel_Kn(2, z)/z/z;
    }
  }
  for(i = 0; i < n; i++) {
    S4mpi[i] = (13./3.)*gg[0] * mmB0[i] - (1./3.)*(40.*gg[0] + 32.*gg[1] + 26.*gg[2])*mmB2[i];
    S4fpi[i] = (1./6.)*(8*gg[0] - 13.*gg[1])* mmB0[i] - 
      (1./3.)*(40.*gg[0] - 12.*gg[1] - 8.*gg[2] - 13.*gg[3])*mmB2[i];
  }
  for(i = 0; i < n; i++) {
    I2mpi[i] = -mmB0[i];
    I4mpi[i] = mmB0[i]*(-55./18. + 4.*lb1[i] + 8./3.*lb2[i] - 5./2.*lb3[i] -2.*lb4[i]) +
      mmB2[i]*(112./9. - (8./3.)*lb1[i] - (32./3.)*lb2[i]) + S4mpi[i];
    if(incim6) {
      I6mpi[i] = mmB0[i]*(10049./1296. - 13./72.*N + 20./9.*lb1[i] - 40./27.*lb2[i] 
			  - 3./4.*lb3[i] - 110./9.*lb4[i] - 5./2.*lb3[i]*lb3[i] 
			  - 5.*lb4[i]*lb4[i] + lb4[i]*(16*lb1[i] + 32./3.*lb2[i] - 11.*lb3[i])
			  + lpi[i]*(70./9.*lpi[i] + 12*lb1[i] + 32./9.*lb2[i] - lb3[i] 
				    + lb4[i] + 47./18.)
			  + 5*rtilde[0] + 4*rtilde[1] + 8*rtilde[2] + 8*rtilde[3] 
			  + 16*rtilde[4] + 16*rtilde[5]) +
	mmB2[i]*(3476./81. - 77./288.*N + 32./9.*lb1[i] + 464./27.*lb2[i] + 448./9.*lb4[i] 
		 - 32./3.*lb4[i]*(lb1[i]+4*lb2[i]) + lpi[i]*(100./9.*lpi[i] + 8./3.*lb1[i] 
							     + 176./9.*lb2[i] - 248./9.)
		 - 8*rtilde[2] - 56*rtilde[3] - 48*rtilde[4] + 16*rtilde[5]);
    }
    else {
      I6mpi[i] = 0.;
    }

    I2fpi[i] = -2*mmB0[i];
    I4fpi[i] = mmB0[i]*(-7./9. + 2*lb1[i] + (4./3.)*lb2[i] - 3*lb4[i]) + 
      mmB2[i]*(112./9. - (8./3.)*lb1[i] -(32./3.)*lb2[i]) + S4fpi[i];
  }
  for(i = 0; i < n; i++) {
/*     Rmpi = - (xi_P[i]/2.) * (ampiV[i]/M_P[i]) * (I2mpi[i] + xi_P[i] * I4mpi[i] + xi_P^2 * I6mpi[i]); */
    mpiFV[i] = ampiV[i] * (1. - rev *(xi_P[i]/2.) * (I2mpi[i] + xi_P[i] * I4mpi[i] + xi_P[i]*xi_P[i] * I6mpi[i]) );

/*     Rfpi = (xi_P[i])   * (afpiV[i]/F_P[i]) * (I2fpi[i] + xi_P[i] * I4fpi[i] + xi_P^2 * I6fpi[i]); */
    fpiFV[i] = afpiV[i]*(1. + rev * (xi_P[i])   * (I2fpi[i] + xi_P[i] * I4fpi[i]));
  }
  if(printit) {
    printf("Rmpi ");
    for(i = 0; i < n; i++) {
      printf("%f ", -(xi_P[i]/2.) * (I2mpi[i] + xi_P[i] * I4mpi[i] + xi_P[i]*xi_P[i] * I6mpi[i]));
    }
    printf("\nRfpi ");
    for(i = 0; i < n; i++) {
      printf("%f ", (xi_P[i])   * (I2fpi[i] + xi_P[i] * I4fpi[i]));
    }
    printf("\n");
  }

  free(lb1); free(lb2); free(lb3); free(lb4); free(lpi);
  free(M_P); free(F_P); free(xi_P); free(mmB0); free(mmB2);
  free(S4mpi); free(S4fpi); free(I4mpi); free(I2mpi);
  free(I2fpi); free(I4fpi); free(I6mpi);
  return;
}

void fscdhnew(double rev, double aLamb1, double aLamb2, double aLamb3, double aLamb4,
	      double aF0_, int * L, double * ampiV, double *afpiV, double *a2B0mu, 
	      const int n, double * mpiFV, double * fpiFV, const int printit) {
  
  const double pi=3.1415926535897932384626433832;
  int i, j;
  double N, amrho_phys, z, Rfpi, Rmpi, lambda_pi;
  double gg[4];
  double mm[] = {6, 12, 8, 6, 24, 24, 0, 12, 30, 24, 24, 8, 24, 48, 0, 6, 48, 36, 24, 24};
  const int mm1 = 20;
  double *lb1, *lb2, *lb3, *lb4, *lpi, *DeltaM, *DeltaF, *xi_P, *mmB0, *mmB1, *mmB2;
  double *S4mpi, *S4fpi, *I4mpi, *I2mpi, *I2fpi, *I4fpi;
  double aF0;

  aF0 = aF0_/sqrt(2.);

  gg[0] = 2. - pi/2.;
  gg[1] = pi/4. - 0.5;
  gg[2] = 0.5 - pi/8.;
  gg[3] = 3*pi/16. - 0.5;
  N = 16*pi*pi;
  lb1 = (double*) malloc(n*sizeof(double));
  lb2 = (double*) malloc(n*sizeof(double));
  lb3 = (double*) malloc(n*sizeof(double));
  lb4 = (double*) malloc(n*sizeof(double));

  for(i = 0; i < n; i++) {
    lb1[i] = log(aLamb1*aLamb1/a2B0mu[i]);
    lb2[i] = log(aLamb2*aLamb2/a2B0mu[i]);
    lb3[i] = log(aLamb3*aLamb3/a2B0mu[i]);
    lb4[i] = log(aLamb4*aLamb4/a2B0mu[i]);
  }

  DeltaM = (double*)   malloc(n*sizeof(double));
  DeltaF = (double*)   malloc(n*sizeof(double));
  mmB0 = (double*)  malloc(n*sizeof(double));
  mmB1 = (double*)  malloc(n*sizeof(double));
  mmB2 = (double*)  malloc(n*sizeof(double));
  xi_P = (double*)  malloc(n*sizeof(double));
  S4mpi = (double*) malloc(n*sizeof(double));
  S4fpi = (double*) malloc(n*sizeof(double));
  I2mpi = (double*) malloc(n*sizeof(double));
  I4mpi = (double*) malloc(n*sizeof(double));
  I2fpi = (double*) malloc(n*sizeof(double));
  I4fpi = (double*) malloc(n*sizeof(double));
  for(i = 0; i < n; i++) {
    xi_P[i] = a2B0mu[i]/(4*pi*aF0)/(4*pi*aF0); 
    mmB0[i] = 0.;
    mmB1[i] = 0.;
    mmB2[i] = 0.;
  }
  for(j = 0; j < n; j++) {
    for(i = 0; i < mm1; i++) {
      z = sqrt(1.+i)*sqrt(a2B0mu[j])*L[j];
      mmB0[j] += mm[i]*2.*gsl_sf_bessel_K0(z);
      mmB1[j] += mm[i]*2.*gsl_sf_bessel_K1(z)/z;
      mmB2[j] += mm[i]*2.*gsl_sf_bessel_Kn(2, z)/z/z;
    }
    DeltaM[j] = -lb3[j]/(2.*N);
    DeltaF[j] = 2.*lb4[j]/N;
  }

  for(i = 0; i < n; i++) {
    S4mpi[i] = (13./3.)*gg[0] * mmB1[i] - (1./3.)*(40.*gg[0] + 32.*gg[1] + 26.*gg[2])*mmB2[i];
    S4fpi[i] = (1./6.)*(8*gg[0] - 13.*gg[1])* mmB1[i] - 
      (1./3.)*(40.*gg[0] - 12.*gg[1] - 8.*gg[2] - 13.*gg[3])*mmB2[i];
  }
  for(i = 0; i < n; i++) {
    I2mpi[i] = -mmB1[i];
    I4mpi[i] = mmB1[i]*(-55./18. + 4.*lb1[i] + 8./3.*lb2[i] - 5./2.*lb3[i] -2.*lb4[i]) +
      mmB2[i]*(112./9. - (8./3.)*lb1[i] - (32./3.)*lb2[i]) + S4mpi[i] +
      N/2.*DeltaM[i]*mmB0[i] + N*DeltaF[i]*mmB1[i];
  }
  for(i = 0; i < n; i++) {
    I2fpi[i] = -2*mmB1[i];
    I4fpi[i] =  mmB1[i]*(-7./9. + 2.*lb1[i] + (4./3.)*lb2[i] - 3.*lb4[i]) + 
      mmB2[i]*(112./9. - (8./3.)*lb1[i] -(32./3.)*lb2[i]) + 
      S4fpi[i] + N*DeltaM[i]*mmB0[i] + 2*N*DeltaF[i]*mmB1[i];
  }
  for(i = 0; i < n; i++) {
    mpiFV[i] = ampiV[i] * (1 + rev* ( - (xi_P[i]/2) * (I2mpi[i] + xi_P[i] * I4mpi[i])));
    fpiFV[i] = afpiV[i] * (1 + rev* ((xi_P[i])   * (I2fpi[i] + xi_P[i] * I4fpi[i])));
  }

  if(printit) {
    printf("Rmpi ");
    for(i = 0; i < n; i++) {
      printf("%f ", -(xi_P[i]/2.) * (I2mpi[i] + xi_P[i] * I4mpi[i]));
    }
    printf("\nRfpi ");
    for(i = 0; i < n; i++) {
      printf("%f ", (xi_P[i])   * (I2fpi[i] + xi_P[i] * I4fpi[i]));
    }
    printf("\n");
  }

  free(lb1); free(lb2); free(lb3); free(lb4); free(lpi);
  free(DeltaM); free(DeltaF); free(xi_P); free(mmB0); free(mmB1); 
  free(S4mpi); free(S4fpi); free(I4mpi); free(I2mpi);
  free(I2fpi); free(I4fpi); free(mmB2);
  return;
}

SEXP cdh_c(SEXP rev, SEXP L1, SEXP L2, SEXP L3, SEXP L4, SEXP F0, SEXP a, SEXP L, 
	   SEXP mpi, SEXP fpi, SEXP printit, SEXP rtilde, SEXP incim6) {

  double *revp, *L1p, *L2p, *L3p, *L4p, *F0p, *ap, *mpip, *fpip, *resp, *rtildep;
  double *fpiFV, *mpiFV;
  int *Lp, *printitp, *incim6p;
  SEXP res;
  int N, i;

  PROTECT(rev = AS_NUMERIC(rev));
  PROTECT(L1 = AS_NUMERIC(L1));
  PROTECT(L2 = AS_NUMERIC(L2));
  PROTECT(L3 = AS_NUMERIC(L3));
  PROTECT(L4 = AS_NUMERIC(L4));
  PROTECT(F0 = AS_NUMERIC(F0));
  PROTECT(a = AS_NUMERIC(a));
  PROTECT(L = AS_INTEGER(L));
  PROTECT(mpi = AS_NUMERIC(mpi));
  PROTECT(fpi = AS_NUMERIC(fpi));
  PROTECT(printit = AS_INTEGER(printit));
  PROTECT(rtilde = AS_NUMERIC(rtilde));
  PROTECT(incim6 = AS_INTEGER(incim6));

  revp = NUMERIC_POINTER(rev);
  L1p = NUMERIC_POINTER(L1);
  L2p = NUMERIC_POINTER(L2);
  L3p = NUMERIC_POINTER(L3);
  L4p = NUMERIC_POINTER(L4);
  F0p = NUMERIC_POINTER(F0);
  ap = NUMERIC_POINTER(a);
  Lp = INTEGER_POINTER(L);
  mpip = NUMERIC_POINTER(mpi);
  fpip = NUMERIC_POINTER(fpi);
  printitp = INTEGER_POINTER(printit);
  rtildep = NUMERIC_POINTER(rtilde);
  incim6p = INTEGER_POINTER(incim6);

  N = LENGTH(mpi);

  PROTECT(res = NEW_NUMERIC(2*N));
  resp = NUMERIC_POINTER(res);
  mpiFV = (double*) malloc(N*sizeof(double));
  fpiFV = (double*) malloc(N*sizeof(double));

  fscdh(revp[0], L1p[0], L2p[0], L3p[0], L4p[0], F0p[0], ap[0], Lp, 
	mpip, fpip, N, mpiFV, fpiFV, printitp[0], rtildep, incim6p[0]);

  for(i = 0; i < N; i++) {
    resp[i] = mpiFV[i];
    resp[N+i] = fpiFV[i];
  }
  free(mpiFV);
  free(fpiFV);
  UNPROTECT(14);
  return(res);
}

SEXP cdhnew_c(SEXP rev, SEXP L1, SEXP L2, SEXP L3, SEXP L4, SEXP F0, SEXP a2B0mu,
	      SEXP L, SEXP mpi, SEXP fpi, SEXP printit) {

  double *revp, *L1p, *L2p, *L3p, *L4p, *F0p, *mpip, *fpip, *resp, *a2B0mup;
  double *fpiFV, *mpiFV;
  int *Lp, *printitp;
  SEXP res;
  int N, i;

  PROTECT(rev = AS_NUMERIC(rev));
  PROTECT(L1 = AS_NUMERIC(L1));
  PROTECT(L2 = AS_NUMERIC(L2));
  PROTECT(L3 = AS_NUMERIC(L3));
  PROTECT(L4 = AS_NUMERIC(L4));
  PROTECT(F0 = AS_NUMERIC(F0));
  PROTECT(L = AS_INTEGER(L));
  PROTECT(mpi = AS_NUMERIC(mpi));
  PROTECT(fpi = AS_NUMERIC(fpi));
  PROTECT(a2B0mu = AS_NUMERIC(a2B0mu));
  PROTECT(printit = AS_INTEGER(printit));

  revp = NUMERIC_POINTER(rev);
  L1p = NUMERIC_POINTER(L1);
  L2p = NUMERIC_POINTER(L2);
  L3p = NUMERIC_POINTER(L3);
  L4p = NUMERIC_POINTER(L4);
  F0p = NUMERIC_POINTER(F0);
  Lp = INTEGER_POINTER(L);
  mpip = NUMERIC_POINTER(mpi);
  fpip = NUMERIC_POINTER(fpi);
  a2B0mup = NUMERIC_POINTER(a2B0mu);
  printitp = INTEGER_POINTER(printit);

  N = LENGTH(mpi);

  PROTECT(res = NEW_NUMERIC(2*N));
  resp = NUMERIC_POINTER(res);
  mpiFV = (double*) malloc(N*sizeof(double));
  fpiFV = (double*) malloc(N*sizeof(double));

  fscdhnew(revp[0], L1p[0], L2p[0], L3p[0], L4p[0], F0p[0], Lp, mpip, fpip, a2B0mup, N, mpiFV, fpiFV, printitp[0]);

  for(i = 0; i < N; i++) {
    resp[i] = mpiFV[i];
    resp[N+i] = fpiFV[i];
  }
  free(mpiFV);
  free(fpiFV);
  UNPROTECT(12);
  return(res);
}
