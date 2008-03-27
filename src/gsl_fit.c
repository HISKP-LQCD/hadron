#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>

/* #define _DEBUG */

struct data {
  int N;
  int tr;
  int Thalf;
  double * x;
  double * y;
  double * err;
};

int exp_S1_f(const gsl_vector * x, void *data, 
	      gsl_vector * f)
{
  int N = ((struct data *)data)->N;
  int tr = ((struct data *)data)->tr;
  double Time = 2.*((struct data *)data)->Thalf;
  double *t = ((struct data *)data)->x;
  double *y = ((struct data *)data)->y;
  double *err = ((struct data *) data)->err;
  double p1, p2, p3, p4, p5, p6, m;

  size_t i, j;
  double Y = 0.;
  double sign = 1.;
  double c = 0.;

  p1 = gsl_vector_get (x, 0);
  p2 = gsl_vector_get (x, 1);
  if(N > 2) {
    p3 = gsl_vector_get (x, 2);
    p4 = gsl_vector_get (x, 3);
  }
  if(N > 4) {
    p5 = gsl_vector_get (x, 4);
    p6 = gsl_vector_get (x, 5);
  }
  m  = gsl_vector_get (x, N);
  
  /* PP */
  for(j = 0; j < 3; j++) {
    if (j==0) c = p1*p1;
    else if (j==1) c = p1*p2;
    else c = p2*p2;
    for (i = 0; i < tr; i++) {
      Y = c*0.5*(exp(-m*(Time-t[i])) + sign*exp(-m*t[i]));
      gsl_vector_set (f, i+j*tr, (Y - y[i+j*tr])/err[i+j*tr]);
    }
  }
  

  if(N > 2) {
    /* PA, AP */
    sign = -1.;
    for(j = 0; j < 4; j++) {
      if (j==0) c = p1*p3;
      else if (j==1) c = p1*p4;
      else if (j==2) c = p2*p3;
      else c = p2*p4;
      for (i = 0; i < tr; i++) {
	Y = c*0.5*(exp(-m*(Time-t[i])) + sign*exp(-m*t[i]));
	gsl_vector_set (f, i+(3+j)*tr, (Y - y[i+(3+j)*tr])/err[i+(3+j)*tr]);
      }
    }
    /* AA */
    sign = 1.;
    for(j = 0; j < 3; j++) {
      if (j==0) c = p3*p3;
      else if (j==1) c = p3*p4;
      else c = p4*p4;
      for (i = 0; i < tr; i++) {
	Y = c*0.5*(exp(-m*(Time-t[i])) + sign*exp(-m*t[i]));
	gsl_vector_set (f, i+(7+j)*tr, (Y - y[i+(7+j)*tr])/err[i+(7+j)*tr]);
      }
    }
  }

  if(N > 4) {
    /* 44 */
    sign = 1.;
    for(j = 0; j < 3; j++) {
      if (j==0) c = p5*p5;
      else if (j==1) c = p5*p6;
      else c = p6*p6;
      for (i = 0; i < tr; i++) {
	Y = c*0.5*(exp(-m*(Time-t[i])) + sign*exp(-m*t[i]));
	gsl_vector_set (f, i+(10+j)*tr, (Y - y[i+(10+j)*tr])/err[i+(10+j)*tr]);
      }
    }
    /* P4 */
    sign = -1.;
    for(j = 0; j < 4; j++) {
      if (j==0) c = p1*p5;
      else if (j==1) c = p1*p6;
      else if (j==2) c = p2*p5;
      else c = p2*p6;
      for (i = 0; i < tr; i++) {
	Y = c*0.5*(exp(-m*(Time-t[i])) + sign*exp(-m*t[i]));
	gsl_vector_set (f, i+(13+j)*tr, (Y - y[i+(13+j)*tr])/err[i+(13+j)*tr]);
      }
    }
    /* 4A */
    sign = 1.;
    for(j = 0; j < 4; j++) {
      if (j==0) c = p3*p5;
      else if (j==1) c = p3*p6;
      else if (j==2) c = p4*p5;
      else c = p4*p6;
      for (i = 0; i < tr; i++) {
	Y = c*0.5*(exp(-m*(Time-t[i])) + sign*exp(-m*t[i]));
	gsl_vector_set (f, i+(17+j)*tr, (Y - y[i+(17+j)*tr])/err[i+(17+j)*tr]);
      }
    }
  }
  
  return GSL_SUCCESS;
}


int exp_S1_df (const gsl_vector * x, void *data, 
		  gsl_matrix * J)
{
  int tr = ((struct data *)data)->tr;
  int N = ((struct data *)data)->N;
  double Time = 2.*((struct data *)data)->Thalf;
  double *t = ((struct data *)data)->x;
  double *y = ((struct data *)data)->y;
  double *err = ((struct data *) data)->err;
  int npar = N+1;
  double p1, p2, p3, p4, p5, p6, m;
  p1 = gsl_vector_get (x, 0);
  p2 = gsl_vector_get (x, 1);
  m = gsl_vector_get (x, N);
  if(N>2) {
    p3 = gsl_vector_get (x, 2);
    p4 = gsl_vector_get (x, 3);
  }
  if(N>4) {
    p5 = gsl_vector_get (x, 4);
    p6 = gsl_vector_get (x, 5);
  }

  size_t i, j;
  double Y = 0., dY=0., sign = 1.;

  j = 0;
  for (i = 0; i < tr; i++) {
    Y = 0.5*(exp(-m*(Time-t[i])) + exp(-m*t[i]));
    dY = -0.5*((Time-t[i])*exp(-m*(Time-t[i])) + t[i]*exp(-m*t[i]));
    j = 0;
    gsl_matrix_set(J, i+j*tr, 0, 2.*p1*Y/err[i+j*tr]);
    gsl_matrix_set(J, i+j*tr, 1, 0.);
    gsl_matrix_set(J, i+j*tr, N, p1*p1*dY/err[i+j*tr]);
    j = 1;
    gsl_matrix_set(J, i+j*tr, 0, p2*Y/err[i+j*tr]);
    gsl_matrix_set(J, i+j*tr, 1, p1*Y/err[i+j*tr]);
    gsl_matrix_set(J, i+j*tr, N, p1*p2*dY/err[i+j*tr]);
    j = 2;
    gsl_matrix_set(J, i+j*tr, 0, 0.);
    gsl_matrix_set(J, i+j*tr, 1, 2.*p2*Y/err[i+j*tr]);
    gsl_matrix_set(J, i+j*tr, N, p1*p1*dY/err[i+j*tr]);
    
  }
  if (N > 2) {
    for (i = 0; i < tr; i++) {
      for(j = 0; j < 3; j++) {
	gsl_matrix_set(J, i+j*tr, 2, 0.);
	gsl_matrix_set(J, i+j*tr, 3, 0.);
      }
    }
    for (i = 0; i < tr; i++) {
      Y = 0.5*(exp(-m*(Time-t[i])) - exp(-m*t[i]));
      dY = -0.5*((Time-t[i])*exp(-m*(Time-t[i])) - t[i]*exp(-m*t[i]));
      j = 3;
      gsl_matrix_set(J, i+j*tr, 0, p3*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, p1*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 3, 0.);
      gsl_matrix_set(J, i+j*tr, N, p1*p3*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, p4*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, 0.);
      gsl_matrix_set(J, i+j*tr, 3, p1*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, N, p1*p4*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, p3*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 2, p2*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 3, 0.);
      gsl_matrix_set(J, i+j*tr, N, p2*p3*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, p4*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 2, 0.);
      gsl_matrix_set(J, i+j*tr, 3, p2*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, N, p2*p4*dY/err[i+j*tr]);
    }
    for (i = 0; i < tr; i++) {
      Y = 0.5*(exp(-m*(Time-t[i])) + exp(-m*t[i]));
      dY = -0.5*((Time-t[i])*exp(-m*(Time-t[i])) + t[i]*exp(-m*t[i]));
      j = 7;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, 2.*p3*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 3, 0.);
      gsl_matrix_set(J, i+j*tr, N, p3*p3*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, p4*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 3, p3*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, N, p3*p4*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, 0.);
      gsl_matrix_set(J, i+j*tr, 3, 2.*p4*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, N, p4*p4*dY/err[i+j*tr]);
    }
  }
  if (N > 4) {
    for (i = 0; i < tr; i++) {
      for(j = 0; j < 10; j++) {
	gsl_matrix_set(J, i+j*tr, 4, 0.);
	gsl_matrix_set(J, i+j*tr, 5, 0.);
      }
    }
    for (i = 0; i < tr; i++) {
      Y = 0.5*(exp(-m*(Time-t[i])) + exp(-m*t[i]));
      dY = -0.5*((Time-t[i])*exp(-m*(Time-t[i])) + t[i]*exp(-m*t[i]));
      j = 10;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, 0.);
      gsl_matrix_set(J, i+j*tr, 3, 0.);
      gsl_matrix_set(J, i+j*tr, 4, 2.*p5*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 5, 0.);
      gsl_matrix_set(J, i+j*tr, N, p5*p5*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, 0.);
      gsl_matrix_set(J, i+j*tr, 3, 0.);
      gsl_matrix_set(J, i+j*tr, 4, p6*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 5, p5*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, N, p5*p6*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, 0.);
      gsl_matrix_set(J, i+j*tr, 3, 0.);
      gsl_matrix_set(J, i+j*tr, 4, 0.);
      gsl_matrix_set(J, i+j*tr, 5, 2.*p6*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, N, p6*p6*dY/err[i+j*tr]);
    }
    for (i = 0; i < tr; i++) {
      Y = 0.5*(exp(-m*(Time-t[i])) - exp(-m*t[i]));
      dY = -0.5*((Time-t[i])*exp(-m*(Time-t[i])) - t[i]*exp(-m*t[i]));
      j = 13;
      gsl_matrix_set(J, i+j*tr, 0, p5*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, 0.);
      gsl_matrix_set(J, i+j*tr, 3, 0.);
      gsl_matrix_set(J, i+j*tr, 4, p1*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 5, 0.);
      gsl_matrix_set(J, i+j*tr, N, p1*p5*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, p6*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, 0.);
      gsl_matrix_set(J, i+j*tr, 3, 0.);
      gsl_matrix_set(J, i+j*tr, 4, 0.);
      gsl_matrix_set(J, i+j*tr, 5, p1*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, N, p1*p6*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, p5*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 2, 0.);
      gsl_matrix_set(J, i+j*tr, 3, 0.);
      gsl_matrix_set(J, i+j*tr, 4, p2*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 5, 0.);
      gsl_matrix_set(J, i+j*tr, N, p2*p5*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, p6*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 2, 0.);
      gsl_matrix_set(J, i+j*tr, 3, 0.);
      gsl_matrix_set(J, i+j*tr, 4, 0.);
      gsl_matrix_set(J, i+j*tr, 5, p2*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, N, p2*p6*dY/err[i+j*tr]);
    }
    for (i = 0; i < tr; i++) {
      Y = 0.5*(exp(-m*(Time-t[i])) + exp(-m*t[i]));
      dY = -0.5*((Time-t[i])*exp(-m*(Time-t[i])) + t[i]*exp(-m*t[i]));
      j = 17;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, p5*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 3, 0.);
      gsl_matrix_set(J, i+j*tr, 4, p3*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 5, 0.);
      gsl_matrix_set(J, i+j*tr, N, p3*p5*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, p6*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 3, 0.);
      gsl_matrix_set(J, i+j*tr, 4, 0.);
      gsl_matrix_set(J, i+j*tr, 5, p3*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, N, p3*p6*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, 0.);
      gsl_matrix_set(J, i+j*tr, 3, p5*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 4, p4*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 5, 0.);
      gsl_matrix_set(J, i+j*tr, N, p4*p5*dY/err[i+j*tr]);
      j++;
      gsl_matrix_set(J, i+j*tr, 0, 0.);
      gsl_matrix_set(J, i+j*tr, 1, 0.);
      gsl_matrix_set(J, i+j*tr, 2, 0.);
      gsl_matrix_set(J, i+j*tr, 3, p6*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 4, 0.);
      gsl_matrix_set(J, i+j*tr, 5, p4*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, N, p4*p6*dY/err[i+j*tr]);
    }
    
  }
  return GSL_SUCCESS;
}
     
     
int exp_S1_fdf (const gsl_vector * x, void *data,
               gsl_vector * f, gsl_matrix * J)
{
  exp_S1_f (x, data, f);
  exp_S1_df (x, data, J);
  
  return GSL_SUCCESS;
}

void Print_State_Mass_Fit_Helper_1(int iter, gsl_multifit_fdfsolver *solver)
{
  fprintf(stderr, "iter = %4d: p1 = %+9.6lf, p2 = %+9.6lf, m = %+9.6lf (chi_square = %13.6lf).\n", iter, gsl_vector_get(solver->x, 0), gsl_vector_get(solver->x, 1), gsl_vector_get(solver->x, 2), pow(gsl_blas_dnrm2(solver->f), 2.0));
}

SEXP multifit_S1(SEXP par, SEXP Thalf, SEXP x, SEXP y, SEXP err, SEXP tr, 
		 SEXP value, SEXP prec, SEXP N, SEXP max_iter)
{
  int npar, nx, ny, i, j, iter_max;
  double p1, p2, m;
  double *yp, *errp, *parp, *precp, *statep, *valuep;
  double chi_square, red_chi_square;
  int dof;
  int * xp, *Thalfp, *trp, *Np, *mip;
  SEXP state;
  const gsl_multifit_fdfsolver_type *solver_type =
    gsl_multifit_fdfsolver_lmsder;
  /*  Allocate the solver. */
  gsl_multifit_fdfsolver *solver;
  /*  Initialize the data structure. */
  struct data data_struct;
  gsl_multifit_function_fdf function_fdf;
  double * para_initial, c=1.;
  gsl_vector_view para_initial_;
  int status, iter=0, no_points=0;

  PROTECT(par = AS_NUMERIC(par));
  PROTECT(Thalf = AS_INTEGER(Thalf));
  PROTECT(x = AS_INTEGER(x));
  PROTECT(y = AS_NUMERIC(y));
  PROTECT(err = AS_NUMERIC(err));
  PROTECT(tr = AS_INTEGER(tr));
  PROTECT(value = AS_NUMERIC(value));
  PROTECT(prec = AS_NUMERIC(prec));
  PROTECT(N = AS_INTEGER(N));
  PROTECT(max_iter = AS_INTEGER(max_iter));

  xp = INTEGER_POINTER(x);
  Thalfp = INTEGER_POINTER(Thalf);
  trp = INTEGER_POINTER(tr);
  Np = INTEGER_POINTER(N);
  yp = NUMERIC_POINTER(y);
  errp = NUMERIC_POINTER(err);
  parp = NUMERIC_POINTER(par);
  precp = NUMERIC_POINTER(prec);
  valuep = NUMERIC_POINTER(value);
  mip = INTEGER_POINTER(max_iter);
  iter_max = mip[0];

  npar = LENGTH(par);
  nx = LENGTH(x);
  ny = LENGTH(y);

  assert(npar == Np[0]+1);
  PROTECT(state = NEW_NUMERIC(5));
  statep = NUMERIC_POINTER(state);
  
  if(Np[0] == 2) no_points = 3*trp[0];
  if(Np[0] == 4) no_points = 10*trp[0];
  if(Np[0] == 6) no_points = 21*trp[0];

  solver = gsl_multifit_fdfsolver_alloc(solver_type, ny, npar);

  data_struct.x = (double*) malloc(nx*sizeof(double));
  data_struct.y = (double*) malloc(ny*sizeof(double));
  data_struct.err = (double*) malloc(ny*sizeof(double));
  para_initial = (double*) malloc(npar*sizeof(double));
  for(i = 0; i < nx; i++) {
    data_struct.x[i] = (double)xp[i];
  }

  for(i = 0; i < ny; i++) {
    data_struct.y[i] = yp[i];
    data_struct.err[i] = errp[i];
  }

  data_struct.Thalf = Thalfp[0];
  data_struct.tr = trp[0];
  data_struct.N = Np[0];

  // The ansatz.
  function_fdf.f = &exp_S1_f;
  function_fdf.df = &exp_S1_df;
  function_fdf.fdf = &exp_S1_fdf;
  function_fdf.n = ny;
  function_fdf.p = npar;
  function_fdf.params = &data_struct;

  for(i = 0; i < npar; i++) {
    para_initial[i] = parp[i];
  }

  para_initial_ = gsl_vector_view_array(para_initial, npar);

  gsl_multifit_fdfsolver_set(solver, &function_fdf, &para_initial_.vector);

  // Perform the fit.
  // Print the initial state.
#ifdef _DEBUG
  Print_State_Mass_Fit_Helper_1(iter, solver);
#endif

  do {
    iter++;
    
    /*  Do a solver iteration. */
    status = gsl_multifit_fdfsolver_iterate(solver);
#ifdef _DEBUG
    fprintf(stderr, "status = %s.\n", gsl_strerror(status));
    Print_State_Mass_Fit_Helper_1(iter, solver);
#endif
    
    if(status) break;
    
    status = gsl_multifit_test_delta(solver->dx, solver->x,
				     precp[0], precp[1]);
  }
  while(status == GSL_CONTINUE && iter < iter_max);
#ifdef _DEBUG
  fprintf(stderr, "\n");
#endif
  
  // *****
  
  
  // Compute the covariance matrix.

#ifdef _DEBUG  
  gsl_matrix *covar = gsl_matrix_alloc(npar, npar);
  gsl_multifit_covar(solver->J, 0.0, covar);
#endif

  // Output.

  chi_square = pow(gsl_blas_dnrm2(solver->f), 2.0);
#ifdef _DEBUG
  fprintf(stderr, "chi_square = %13.6lf.\n", chi_square);
#endif
  dof = no_points - npar;
#ifdef _DEBUG
  fprintf(stderr, "dof = %d\n", dof);
#endif
  red_chi_square = chi_square / (double)dof;
#ifdef _DEBUG
  fprintf(stderr, "red_chi_square = %13.6lf.\n", red_chi_square);
  fprintf(stderr, "\n");
#endif
  p1 = gsl_vector_get(solver->x, 0);
  p2 = gsl_vector_get(solver->x, 1);
  m = gsl_vector_get(solver->x, npar-1);

  if(red_chi_square > 1.0)
    c = sqrt(red_chi_square);
#ifdef _DEBUG
  fprintf(stderr, "p1 = %+9.6lf +/- %9.6lf.\n",
      p1, c * sqrt(gsl_matrix_get(covar, 0, 0)));
  fprintf(stderr, "p2 = %+9.6lf +/- %9.6lf.\n",
      p2, c * sqrt(gsl_matrix_get(covar, 1, 1)));
  fprintf(stderr, "m = %+9.6lf +/- %9.6lf.\n",
      m, c * sqrt(gsl_matrix_get(covar, npar-1, npar-1)));
  fprintf(stderr, "\n");
  fprintf(stderr, "status = %s.\n", gsl_strerror(status));
  fprintf(stderr, "\n");
#endif

  for(i = 0; i < npar; i++) {
    parp[i] =  gsl_vector_get(solver->x, i);
  }

  valuep[0] = chi_square;
  statep[0] = chi_square;
  statep[1] = gsl_blas_dnrm2(solver->f);
  statep[2] = (double)iter;
  statep[3] = (double)dof;
  statep[4] = (double)status;

  gsl_multifit_fdfsolver_free(solver);
#ifdef _DEBUG
  gsl_matrix_free(covar);
#endif

  free(data_struct.x);
  free(data_struct.y);
  free(data_struct.err);
  free(para_initial);

  UNPROTECT(11);
  return(state);
}
