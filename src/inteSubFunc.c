#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <math.h>

/******************************************
 *
 * This program calculates the integral of
 * the second part:
 *  \int_0_Lamda exp{t*qSqur}/sqrt(t) dt
 *
 *  **************************************/

#define WSNUM 1000
#define EPSREL 1e-8

double integrandPart2(const double t, void *const params) {
  double alpha = *(double *)params;
  double f2 = exp(t * alpha) / sqrt(t);
  return f2;
}

/***************************************
 *
 * For qag---quadture adaptive genearl,
 * the key value means:
 * GSL_INTEG_GAUSS15  (key = 1)
 * GSL_INTEG_GAUSS21  (key = 2)
 * GSL_INTEG_GAUSS31  (key = 3)
 * GSL_INTEG_GAUSS41  (key = 4)
 * GSL_INTEG_GAUSS51  (key = 5)
 * GSL_INTEG_GAUSS61  (key = 6)
 *
 ***************************************/
double sndInteFunc(const double Lamda, const double qSqur, int *const rstatus) {
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(WSNUM);
  double result, error;
  double qsq = qSqur;
  gsl_function FuncPart2;
  FuncPart2.function = &integrandPart2;
  FuncPart2.params = &qsq;

  *rstatus =
      gsl_integration_qags(&FuncPart2, 0, Lamda, 0, EPSREL, WSNUM, w, &result, &error);
  // gsl_integration_qag(&FuncPart2, 0, Lamda, 0, EPSREL, WSNUM, 2, w, &result, &error);
  if (*rstatus != GSL_SUCCESS) {
    warning("gsl_integration_qags in sndInteFunc error code \"%s\"\n",
            gsl_strerror(*rstatus));
  }
  gsl_integration_workspace_free(w);

  return result;
}

// The integrand is
//(M_PI/t)^{3/2+l}*exp(t*qSqur-M_PI^2*|\hat{gamma}\vec{w}|^2/t)
//
// dVec[0]=paraArray[0], dVec[1]=paraArray[1], dVec[2]=paraArray[2],
// l=paraArray[3], qSqur=paraArray[4],
// nVec[0] = n1 = paraArray[5],
// nVec[1] = n2 = paraArray[6],
// nVec[2] = n3 = paraArray[7],
// gamma = paraArray[8]
double integrandPart3(const double t, void *const params) {
  double *paraArray = (double *)params;
  double f3;
  double dModSqur = paraArray[0] * paraArray[0] + paraArray[1] * paraArray[1] +
                    paraArray[2] * paraArray[2];

  if (dModSqur == 0) {
    f3 = pow(M_PI / t, 3.0 / 2.0 + (*(paraArray + 3))) *
         exp(t * (*(paraArray + 4)) -
             pow(M_PI, 2) *
                 (paraArray[5] * paraArray[5] + paraArray[6] * paraArray[6] +
                  paraArray[7] * paraArray[7]) /
                 t);
  } else {
    f3 = pow(M_PI / t, 3.0 / 2.0 + (*(paraArray + 3))) *
         exp(t * (*(paraArray + 4)) -
             pow(M_PI, 2) *
                 ((pow(paraArray[8], 2) - 1) *
                      pow(paraArray[0] * paraArray[5] + paraArray[1] * paraArray[6] +
                              paraArray[2] * paraArray[7],
                          2) /
                      (paraArray[0] * paraArray[0] + paraArray[1] * paraArray[1] +
                       paraArray[2] * paraArray[2]) +
                  paraArray[5] * paraArray[5] + paraArray[6] * paraArray[6] +
                  paraArray[7] * paraArray[7]) /
                 t);
  }
  return f3;
}

double trdInteFunc(const double Lamda,
                   double *const dVec,
                   const int l,
                   const double qSqur,
                   int *const nVec,
                   const double gamma,
                   int *const rstatus) {
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(WSNUM);
  double result, error;

  // dVec[0]=paraArray[0], dVec[1]=paraArray[1], dVec[2]=paraArray[2],
  // l=paraArray[3], qSqur=paraArray[4],
  // nVec[0] = n1 = paraArray[5],
  // nVec[1] = n2 = paraArray[6],
  // nVec[2] = n3 = paraArray[7],
  // gamma = paraArray[8]

  double paraArray[9] = {dVec[0],
                         dVec[1],
                         dVec[2],
                         (double)l,
                         qSqur,
                         (double)nVec[0],
                         (double)nVec[1],
                         (double)nVec[2],
                         gamma};
  gsl_function FuncPart3;
  FuncPart3.function = &integrandPart3;
  FuncPart3.params = paraArray;

  *rstatus =
      gsl_integration_qags(&FuncPart3, 0, Lamda, 0, EPSREL, WSNUM, w, &result, &error);
  if (*rstatus != GSL_SUCCESS) {
    warning("gsl_integration_qags in trdInteFunc error code \"%s\"\n",
            gsl_strerror(*rstatus));
  }
  gsl_integration_workspace_free(w);

  return result;
}
