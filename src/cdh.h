#pragma once

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

double g1(double x);
int g1array(double *x, double *res, const int n);
static R_INLINE void fscdh(double rev,
                           double aLamb1,
                           double aLamb2,
                           double aLamb3,
                           double aLamb4,
                           double *aF0,
                           double a_fm,
                           int *L,
                           double *ampiV,
                           double *afpiV,
                           const int n,
                           double *mpiFV,
                           double *fpiFV,
                           const int printit,
                           double *rtilde,
                           const int incim6);

static R_INLINE void fscdhnew(double rev,
                              double aLamb1,
                              double aLamb2,
                              double aLamb3,
                              double aLamb4,
                              double aF0,
                              int *L,
                              double *ampiV,
                              double *afpiV,
                              double *a2B0mu,
                              const int n,
                              double *mpiFV,
                              double *fpiFV,
                              const int printit);
