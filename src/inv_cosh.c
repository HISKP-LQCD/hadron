#include <R.h>
#include <Rinternals.h>

#include <math.h>

SEXP invcosh(SEXP ratio, SEXP timeextent, SEXP t, SEXP eps, SEXP maxiter) {
  const double rat = asReal(ratio), epsilon = asReal(eps);
  const int dt0 = asInteger(timeextent) - 2 * asInteger(t), dt1 = dt0 + 2;
  const int n = asInteger(maxiter);
  double newmass = log(rat), mass = 0, r;
  int i;

  for (i = 0; i < n && fabs(mass - newmass) >= epsilon * mass; i++) {
    mass = newmass;
    r = (1 + exp(-mass * dt0)) / (1 + exp(-mass * dt1));
    newmass = log(rat * r);
  }

  return ScalarReal(newmass);
}
