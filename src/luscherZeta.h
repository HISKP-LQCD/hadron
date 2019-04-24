#pragma once

void luscherZeta(double *res,
                 const double qsq,
                 const int l,
                 const int m,
                 const double gamma,
                 const double lambda,
                 double *const dvec,
                 const double tol,
                 const int verbose,
                 int *rstatus);
