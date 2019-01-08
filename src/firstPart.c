#include "zetaFunc.h"

#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_errno.h>

#include <complex.h>
#include <float.h>
#include <math.h>

double complex firstPart(const double Tolerance,
                         const int l,
                         const int m,
                         const double *dVec,
                         const double gamma,
                         const double Lamda,
                         const double qSqur,
                         const int verbose,
                         int *const rstatus) {
  double complex firstTerms = 0 + I * 0, pmodeSum = 0 + I * 0, firstPartSum = 0 + I * 0;
  int npmode[2] = {40, 72};
  int i_npmode = 1;

  double error = 1.0;
  int pmodeSqur;

  int n0, n1, n2;
  double dModSqur = dVec[0] * dVec[0] + dVec[1] * dVec[1] + dVec[2] * dVec[2];

  int *degnrtDOF = NULL;
  int *arrayPmode = NULL;

  int niter = 0;
  int genReturn;
  get_npmode(npmode, i_npmode);
  genReturn = gen_points_array(&degnrtDOF, &arrayPmode, npmode[0], npmode[1]);

  if (genReturn != 0) {
    Rprintf("Generated the points wrongly in firstPart!");
    *rstatus = 1000;
    return (firstPartSum);
  }

  if (verbose && 0) {
    for (int i = 0; i < npmode[0]; i++) {
      if (degnrtDOF[i] == 0)
        Rprintf("pmodeSqur=%d has no corresponding points.\n", i);
      else
        Rprintf("pmodeSqur=%d have %d degenaration DOF.\n", i, degnrtDOF[i]);
    }
  }

  pmodeSqur = 0;

  while (error > Tolerance) {
    if (pmodeSqur >= npmode[0]) {
      i_npmode++;
      get_npmode(npmode, i_npmode);
      genReturn = gen_points_array(&degnrtDOF, &arrayPmode, npmode[0], npmode[1]);
      if (verbose)
        REprintf("increased i_npmode to %d in firstPart %d %d\n",
                 i_npmode,
                 npmode[0],
                 npmode[1]);

      if (genReturn != 0) {
        Rprintf("Generated the points wrongly in firstPart!");
        *rstatus = 1000;
        return (firstPartSum);
      }

      if (i_npmode > 4) {
        if (verbose)
          REprintf(
              "NPmode and DimMax need to be larger than available in firstPart! "
              "Aborting...!\n");
        *rstatus = 2000;
        return (firstPartSum);
      }
    }

    pmodeSum = 0 + I * 0;

    // These pmodes has no contribution to the points sets.
    if (degnrtDOF[pmodeSqur] == 0) {
      pmodeSqur += 1;
      continue;
    }

    for (int i = 0; i < degnrtDOF[pmodeSqur]; i++) {
      n0 = arrayPmode[pmodeSqur * npmode[1] * 3 + i * 3 + 0];
      n1 = arrayPmode[pmodeSqur * npmode[1] * 3 + i * 3 + 1];
      n2 = arrayPmode[pmodeSqur * npmode[1] * 3 + i * 3 + 2];

      double r[3], rpar[3], rort[3];
      if (fabs(dModSqur) < DBL_EPSILON) {
        r[0] = n0 / gamma;
        r[1] = n1 / gamma;
        r[2] = n2 / gamma;
      } else {
        double nDotd = n0 * dVec[0] + n1 * dVec[1] + n2 * dVec[2];
        // we split the vector first into a parallel and orthogonal part w.r.t. dVec
        rpar[0] = nDotd / dModSqur * dVec[0];
        rpar[1] = nDotd / dModSqur * dVec[1];
        rpar[2] = nDotd / dModSqur * dVec[2];
        rort[0] = n0 - rpar[0];
        rort[1] = n1 - rpar[1];
        rort[2] = n2 - rpar[2];
        r[0] = (rpar[0] - 0.5 * dVec[0]) / gamma + rort[0];
        r[1] = (rpar[1] - 0.5 * dVec[1]) / gamma + rort[1];
        r[2] = (rpar[2] - 0.5 * dVec[2]) / gamma + rort[2];
      }
      // now we determine the spherical coordinates appropriate for
      // usage with GSL spherical harmonics
      double u, v, w, xy;
      xy = r[0] * r[0] + r[1] * r[1];
      u = sqrt(xy + r[2] * r[2]);
      v = atan2(sqrt(xy), r[2]);
      w = atan2(r[1], r[0]) * 180 / M_PI;
      r[0] = u;
      r[1] = cos(v);
      r[2] = w;

      firstTerms = exp(-Lamda * (pow(r[0], 2.0) - qSqur)) * pow(r[0], l) *
                   spheHarm(l, m, r[1], r[2], rstatus) / (pow(r[0], 2.0) - qSqur);

      if (*rstatus != 0) {
        REprintf("spheHarm produced error code \"%s\"\n", gsl_strerror(*rstatus));
        return (firstPartSum);
      }

      // Add every term within the same pmode into pmodeSum
      pmodeSum += firstTerms;

    }  // end of pmode loop

    firstPartSum += pmodeSum;
    // Both pmodeSum and firstPartSum are complex numbers,
    // cabs take the modulus of these variables.
    // only calculate new error if firstPartSum != 0.
    if (cabs(firstPartSum) > DBL_EPSILON)
      error = cabs(pmodeSum) / cabs(firstPartSum);

    if (verbose)
      Rprintf("first term: pmode %d error: %.16f result (%e, %e)\n",
              pmodeSqur,
              error,
              creal(firstPartSum),
              cimag(firstPartSum));

    // if the result is still zero after 4 iterations it is assumed to stay zero
    if (cabs(firstPartSum) < DBL_EPSILON && niter > 4)
      break;
    pmodeSqur += 1;
    ++niter;

  }  // end of while.
  if (verbose) {
    Rprintf("First term = (%e, %e)\n", creal(firstPartSum), cimag(firstPartSum));
  }

  return firstPartSum;
}
