#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <math.h>

SEXP alphas(SEXP mu_, SEXP nl_, SEXP lam0_, SEXP Nc_, SEXP Nf_) {
  const double pi = 3.1415926535897932384626433832;
  const double Z3 = 1.20206;
  double Cf, b0, b1, b2, b3, L2, LL2, als0, als1, als2, als3;
  double mu, lam0, *resp, Nc, Nf;
  int nl;
  SEXP res;

  PROTECT(mu_ = AS_NUMERIC(mu_));
  PROTECT(nl_ = AS_INTEGER(nl_));
  PROTECT(lam0_ = AS_NUMERIC(lam0_));
  PROTECT(Nc_ = AS_NUMERIC(Nc_));
  PROTECT(Nf_ = AS_NUMERIC(Nf_));

  mu = NUMERIC_POINTER(mu_)[0];
  nl = INTEGER_POINTER(nl_)[0];
  lam0 = NUMERIC_POINTER(lam0_)[0];
  Nc = NUMERIC_POINTER(Nc_)[0];
  Nf = NUMERIC_POINTER(Nf_)[0];

  Cf = (Nc * Nc - 1.) / 2. / Nc;

  b0 = 11. / 3. * Nc - 2. / 3. * Nf;
  b1 = 34. / 3. * Nc * Nc - 38. / 3. * Nf;
  b2 = 2857. / 54. * Nc * Nc * Nc + Cf * Cf * Nf - 205. / 18. * Cf * Nc * Nf -
       1415. / 54. * Nc * Nc * Nf + 11. / 9. * Cf * Nf * Nf + 79. / 54. * Nc * Nf * Nf;
  b3 = (150653. / 486. - 44. / 9. * Z3) * Nc * Nc * Nc * Nc +
       (-39143. / 162. + 68. / 3. * Z3) * Nc * Nc * Nc * Nf +
       (7073. / 486. - 328. / 9. * Z3) * Cf * Nc * Nc * Nf +
       (-2102. / 27. + 176. / 9. * Z3) * Cf * Cf * Nc * Nf + 23. * Cf * Cf * Cf * Nf +
       (3965. / 162. + 56. / 9. * Z3) * Nc * Nc * Nf * Nf +
       (338. / 27. - 176. / 9. * Z3) * Cf * Cf * Nf * Nf +
       (4288. / 243. + 112. / 9. * Z3) * Cf * Nc * Nf * Nf +
       53. / 243. * Nc * Nf * Nf * Nf + 154. / 243. * Cf * Nf * Nf * Nf +
       (-10. / 27. + 88. / 9. * Z3) * Nc * Nc * (Nc * Nc + 36.) +
       (32. / 27. - 104. / 9. * Z3) * Nc * (Nc * Nc + 6) * Nf +
       (-22. / 27. + 16. / 9. * Z3) * (Nc * Nc * Nc * Nc - 6. * Nc * Nc + 18.) /
           (Nc * Nc) * Nf * Nf;

  b1 = b1 / b0 / 4. / pi;
  b2 = b2 / b0 / 16. / pi / pi;
  b3 = b3 / b0 / 64. / pi / pi / pi;

  L2 = log(mu * mu / lam0 / lam0);
  LL2 = log(L2);

  als0 = 4. * pi / b0 / L2;
  als1 = als0 - als0 * als0 * b1 * LL2;
  als2 = als1 + als0 * als0 * als0 * (b1 * b1 * (LL2 * LL2 - LL2 - 1.) + b2);
  als3 = als2 + als0 * als0 * als0 * als0 *
                    (b1 * b1 * b1 *
                         (-LL2 * LL2 * LL2 + 5. / 2. * LL2 * LL2 + 2 * LL2 - 1. / 2.) -
                     3. * b1 * b2 * LL2 + b3 / 2.);
  PROTECT(res = NEW_NUMERIC(1));
  resp = NUMERIC_POINTER(res);

  if (nl == 0)
    resp[0] = (als0 / (4. * pi));
  else if (nl == 1)
    resp[0] = (als1 / (4. * pi));
  else if (nl == 2)
    resp[0] = (als2 / (4. * pi));
  else
    resp[0] = (als3 / (4. * pi));
  UNPROTECT(6);
  return (res);
}

double alphas_c(
    const double mu, const int nl, const double lam0, const int Nc, const int Nf) {
  const double pi = 3.1415926535897932384626433832;
  const double Z3 = 1.20206;
  double Cf, b0, b1, b2, b3, L2, LL2, als0, als1, als2, als3;

  Cf = (Nc * Nc - 1.) / 2. / Nc;

  b0 = 11. / 3. * Nc - 2. / 3. * Nf;
  b1 = 34. / 3. * Nc * Nc - 38. / 3. * Nf;
  b2 = 2857. / 54. * Nc * Nc * Nc + Cf * Cf * Nf - 205. / 18. * Cf * Nc * Nf -
       1415. / 54. * Nc * Nc * Nf + 11. / 9. * Cf * Nf * Nf + 79. / 54. * Nc * Nf * Nf;
  b3 = (150653. / 486. - 44. / 9. * Z3) * Nc * Nc * Nc * Nc +
       (-39143. / 162. + 68. / 3. * Z3) * Nc * Nc * Nc * Nf +
       (7073. / 486. - 328. / 9. * Z3) * Cf * Nc * Nc * Nf +
       (-2102. / 27. + 176. / 9. * Z3) * Cf * Cf * Nc * Nf + 23. * Cf * Cf * Cf * Nf +
       (3965. / 162. + 56. / 9. * Z3) * Nc * Nc * Nf * Nf +
       (338. / 27. - 176. / 9. * Z3) * Cf * Cf * Nf * Nf +
       (4288. / 243. + 112. / 9. * Z3) * Cf * Nc * Nf * Nf +
       53. / 243. * Nc * Nf * Nf * Nf + 154. / 243. * Cf * Nf * Nf * Nf +
       (-10. / 27. + 88. / 9. * Z3) * Nc * Nc * (Nc * Nc + 36.) +
       (32. / 27. - 104. / 9. * Z3) * Nc * (Nc * Nc + 6) * Nf +
       (-22. / 27. + 16. / 9. * Z3) * (Nc * Nc * Nc * Nc - 6. * Nc * Nc + 18.) /
           (Nc * Nc) * Nf * Nf;

  b1 = b1 / b0 / 4. / pi;
  b2 = b2 / b0 / 16. / pi / pi;
  b3 = b3 / b0 / 64. / pi / pi / pi;

  L2 = log(mu * mu / lam0 / lam0);
  LL2 = log(L2);

  als0 = 4. * pi / b0 / L2;
  als1 = als0 - als0 * als0 * b1 * LL2;
  als2 = als1 + als0 * als0 * als0 * (b1 * b1 * (LL2 * LL2 - LL2 - 1.) + b2);
  als3 = als2 + als0 * als0 * als0 * als0 *
                    (b1 * b1 * b1 *
                         (-LL2 * LL2 * LL2 + 5. / 2. * LL2 * LL2 + 2 * LL2 - 1. / 2.) -
                     3. * b1 * b2 * LL2 + b3 / 2.);

  if (nl == 0)
    return (als0 / (4. * pi));
  else if (nl == 1)
    return (als1 / (4. * pi));
  else if (nl == 2)
    return (als2 / (4. * pi));
  return (als3 / (4. * pi));
}
