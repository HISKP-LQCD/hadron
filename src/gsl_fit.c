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
  int no_masses;
  int N;
  int tr;
  int Thalf;
  double * x;
  double * y;
  double * err;
};

int exp_f(const gsl_vector * x, void *data, 
	      gsl_vector * f)
{
  int no_masses = ((struct data *)data)->no_masses;
  int N = ((struct data *)data)->N;
  int tr = ((struct data *)data)->tr;
  double Time = 2.*((struct data *)data)->Thalf;
  double *t = ((struct data *)data)->x;
  double *y = ((struct data *)data)->y;
  double *err = ((struct data *) data)->err;
  double p[6][6], m[6];

  size_t i, j, k=0, id0, id1, kludge = 0;
  double Y = 0., sign = 1., c = 0.;

  for(i = 0; i < no_masses; i++) {
    for(j = 0; j < N; j++) {
      p[j][i] = gsl_vector_get (x, j + i*(N+1));
    }
    m[i]  = gsl_vector_get (x, N + i*(N+1));
  }

  /* PP */
  for(j = 0; j < 3; j++) {
    if(j == 0) {
      id0 = 0; id1 = 0;
    }
    if(j == 1) {
      id0 = 0; id1 = 1;
    }
    if(j == 2) {
      id0 = 1; id1 = 1;
    }
    for (i = 0; i < tr; i++) {
      Y = 0.;
      for(k = 0; k < no_masses; k++) {
	c = p[id0][k]*p[id1][k];
	Y += c*0.5*(exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]));
      } 
      if(kludge==1 && j == 1 && no_masses > 1) {
	Y += 1000000.*p[id0][1] + 1000000.*p[id1][1];
	printf("kludge\n");
      }
      gsl_vector_set (f, i+j*tr, (Y - y[i+j*tr])/err[i+j*tr]);
    }
  }
  
  if(N > 2) {
    /* PA, AP */
    sign = -1.;
    for(j = 0; j < 4; j++) {
      if(j == 0) {
	id0 = 0; id1 = 2;
      }
      if(j == 1) {
	id0 = 0; id1 = 3;
      }
      if(j == 2) {
	id0 = 1; id1 = 2;
      }
      if(j == 3) {
	id0 = 1; id1 = 3;
      }
      for (i = 0; i < tr; i++) {
	Y = 0.;
	for(k = 0; k < no_masses; k++) {
	  c = p[id0][k]*p[id1][k];	  
	  Y += c*0.5*(exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]));
	}
	gsl_vector_set (f, i+(3+j)*tr, (Y - y[i+(3+j)*tr])/err[i+(3+j)*tr]);
      }
    }
    /* AA */
    sign = 1.;
    for(j = 0; j < 3; j++) {
      if(j == 0) {
	id0 = 2; id1 = 2;
      }
      if(j == 1) {
	id0 = 2; id1 = 3;
      }
      if(j == 2) {
	id0 = 3; id1 = 3;
      }
      for (i = 0; i < tr; i++) {
	Y = 0.;
	for(k = 0; k < no_masses; k++) {
	  c = p[id0][k]*p[id1][k];	  
	  Y += c*0.5*(exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]));
	}
	gsl_vector_set (f, i+(7+j)*tr, (Y - y[i+(7+j)*tr])/err[i+(7+j)*tr]);
      }
    }
  }

  if(N > 4) {
    /* 44 */
    sign = 1.;
    for(j = 0; j < 3; j++) {
      if(j == 0) {
	id0 = 4; id1 = 4;
      }
      if(j == 1) {
	id0 = 4; id1 = 5;
      }
      if(j == 2) {
	id0 = 5; id1 = 5;
      }
      for (i = 0; i < tr; i++) {
	Y = 0.;
	for(k = 0; k < no_masses; k++) {
	  c = p[id0][k]*p[id1][k];	  
	  Y += c*0.5*(exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]));
	}
	gsl_vector_set (f, i+(10+j)*tr, (Y - y[i+(10+j)*tr])/err[i+(10+j)*tr]);
      }
    }
    /* P[3] */
    sign = -1.;
    for(j = 0; j < 4; j++) {
      if(j == 0) {
	id0 = 0; id1 = 4;
      }
      if(j == 1) {
	id0 = 0; id1 = 5;
      }
      if(j == 2) {
	id0 = 1; id1 = 4;
      }
      if(j == 3) {
	id0 = 1; id1 = 5;
      }
      for (i = 0; i < tr; i++) {
	Y = 0.;
	for(k = 0; k < no_masses; k++) {
	  c = p[id0][k]*p[id1][k];	  
	  Y += c*0.5*(exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]));
	}
	gsl_vector_set (f, i+(13+j)*tr, (Y - y[i+(13+j)*tr])/err[i+(13+j)*tr]);
      }
    }
    /* 4A */
    sign = 1.;
    for(j = 0; j < 4; j++) {
      if(j == 0) {
	id0 = 2; id1 = 4;
      }
      if(j == 1) {
	id0 = 2; id1 = 5;
      }
      if(j == 2) {
	id0 = 3; id1 = 4;
      }
      if(j == 3) {
	id0 = 3; id1 = 5;
      }
      for (i = 0; i < tr; i++) {
	Y = 0.;
	for(k = 0; k < no_masses; k++) {
	  c = p[id0][k]*p[id1][k];	  
	  Y += c*0.5*(exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]));
	}
	gsl_vector_set (f, i+(17+j)*tr, (Y - y[i+(17+j)*tr])/err[i+(17+j)*tr]);
      }
    }
  }
  
  return GSL_SUCCESS;
}


int exp_df (const gsl_vector * x, void *data, 
		  gsl_matrix * J)
{
  int no_masses = ((struct data *)data)->no_masses;
  int tr = ((struct data *)data)->tr;
  int N = ((struct data *)data)->N;
  double Time = 2.*((struct data *)data)->Thalf;
  double *t = ((struct data *)data)->x;
  double *y = ((struct data *)data)->y;
  double *err = ((struct data *) data)->err;
  int npar = N+1;
  double p[6][6], m[6];
  size_t i, j, k=0, kludge = 0;
  double Y = 0., dY=0., sign = 1.;

  for(i = 0; i < no_masses; i++) {
    for(j = 0; j < N; j++) {
      p[j][i] = gsl_vector_get (x, j + i*(N+1));
    }
    m[i]  = gsl_vector_get (x, N + i*(N+1));
  }

  j = 0;
  for (i = 0; i < tr; i++) {
    for(k = 0; k < no_masses; k++) {
      Y = 0.5*(exp(-m[k]*(Time-t[i])) + exp(-m[k]*t[i]));
      dY = -0.5*((Time-t[i])*exp(-m[k]*(Time-t[i])) + t[i]*exp(-m[k]*t[i]));
      j = 0;
      gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 2.*p[0][k]*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
      gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[0][k]*p[0][k]*dY/err[i+j*tr]);
      j = 1;
      gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), p[1][k]*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), p[0][k]*Y/err[i+j*tr]);
      if(k == 1 && kludge == 1) {
	printf("kludge df\n");
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), gsl_matrix_get(J, i+j*tr, 0 + k*(N+1)) + 1000000.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), gsl_matrix_get(J, i+j*tr, 1 + k*(N+1)) + 1000000.);
      }
      gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[0][k]*p[1][k]*dY/err[i+j*tr]);
      j = 2;
      gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
      gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 2.*p[1][k]*Y/err[i+j*tr]);
      gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[0][k]*p[0][k]*dY/err[i+j*tr]);
    }
  }
  if (N > 2) {
    for (i = 0; i < tr; i++) {
      for(j = 0; j < 3; j++) {
	for(k = 0; k < no_masses; k++) {
	  gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	  gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	}
      }
    }
    for (i = 0; i < tr; i++) {
      for (k = 0; k < no_masses; k++) {
	Y = 0.5*(exp(-m[k]*(Time-t[i])) - exp(-m[k]*t[i]));
	dY = -0.5*((Time-t[i])*exp(-m[k]*(Time-t[i])) - t[i]*exp(-m[k]*t[i]));
	j = 3;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), p[2][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), p[0][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[0][k]*p[2][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), p[3][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), p[0][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[0][k]*p[3][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), p[2][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), p[1][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[1][k]*p[2][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), p[3][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), p[1][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[1][k]*p[3][k]*dY/err[i+j*tr]);
      }
    }
    for (i = 0; i < tr; i++) {
      for (k = 0; k < no_masses; k++) {
	Y = 0.5*(exp(-m[k]*(Time-t[i])) + exp(-m[k]*t[i]));
	dY = -0.5*((Time-t[i])*exp(-m[k]*(Time-t[i])) + t[i]*exp(-m[k]*t[i]));
	j = 7;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 2.*p[2][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[2][k]*p[2][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), p[3][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), p[2][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[2][k]*p[3][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 2.*p[3][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[3][k]*p[3][k]*dY/err[i+j*tr]);
      }
    }
  }
  if (N > 4) {
    for (i = 0; i < tr; i++) {
      for(j = 0; j < 10; j++) {
	for(k = 0; k < no_masses; k++) {
	  gsl_matrix_set(J, i+j*tr, 4 + k*(N+1), 0.);
	  gsl_matrix_set(J, i+j*tr, 5 + k*(N+1), 0.);
	}
      }
    }
    for (i = 0; i < tr; i++) {
      for(k = 0; k < no_masses; k++) {
	Y = 0.5*(exp(-m[k]*(Time-t[i])) + exp(-m[k]*t[i]));
	dY = -0.5*((Time-t[i])*exp(-m[k]*(Time-t[i])) + t[i]*exp(-m[k]*t[i]));
	j = 10;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 4 + k*(N+1), 2.*p[4][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 5 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[4][k]*p[4][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 4 + k*(N+1), p[5][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 5 + k*(N+1), p[4][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[4][k]*p[5][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 4 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 5 + k*(N+1), 2.*p[5][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[5][k]*p[5][k]*dY/err[i+j*tr]);
      }
    }
    for (i = 0; i < tr; i++) {
      for(k = 0; k < no_masses; k++) {
	Y = 0.5*(exp(-m[k]*(Time-t[i])) - exp(-m[k]*t[i]));
	dY = -0.5*((Time-t[i])*exp(-m[k]*(Time-t[i])) - t[i]*exp(-m[k]*t[i]));
	j = 13;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), p[4][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 4 + k*(N+1), p[0][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 5 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[0][k]*p[4][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), p[5][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 4 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 5 + k*(N+1), p[0][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[0][k]*p[5][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), p[4][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 4 + k*(N+1), p[1][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 5 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[1][k]*p[4][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), p[5][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 4 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 5 + k*(N+1), p[1][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[1][k]*p[5][k]*dY/err[i+j*tr]);
      }
    }
    for (i = 0; i < tr; i++) {
      for(k = 0; k < no_masses; k++) {
	Y = 0.5*(exp(-m[k]*(Time-t[i])) + exp(-m[k]*t[i]));
	dY = -0.5*((Time-t[i])*exp(-m[k]*(Time-t[i])) + t[i]*exp(-m[k]*t[i]));
	j = 17;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), p[4][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 4 + k*(N+1), p[2][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 5 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[2][k]*p[4][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), p[5][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 4 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 5 + k*(N+1), p[2][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[2][k]*p[5][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), p[4][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 4 + k*(N+1), p[3][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 5 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[3][k]*p[4][k]*dY/err[i+j*tr]);
	j++;
	gsl_matrix_set(J, i+j*tr, 0 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 1 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 2 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 3 + k*(N+1), p[5][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, 4 + k*(N+1), 0.);
	gsl_matrix_set(J, i+j*tr, 5 + k*(N+1), p[3][k]*Y/err[i+j*tr]);
	gsl_matrix_set(J, i+j*tr, N + k*(N+1), p[3][k]*p[5][k]*dY/err[i+j*tr]);
      }
    }
  }
  return GSL_SUCCESS;
}
     
     
int exp_fdf (const gsl_vector * x, void *data,
               gsl_vector * f, gsl_matrix * J)
{
  exp_f (x, data, f);
  exp_df (x, data, J);
  
  return GSL_SUCCESS;
}

void Print_State_Mass_Fit_Helper_1(int iter, gsl_multifit_fdfsolver *solver)
{
  fprintf(stderr, "iter = %4d: p1 = %+9.6lf, p2 = %+9.6lf, m = %+9.6lf (chi_square = %13.6lf).\n", iter, gsl_vector_get(solver->x, 0), gsl_vector_get(solver->x, 1), gsl_vector_get(solver->x, 2), pow(gsl_blas_dnrm2(solver->f), 2.0));
}

SEXP multifit_cor(SEXP par, SEXP Thalf, SEXP x, SEXP y, SEXP err, SEXP tr, 
		  SEXP value, SEXP prec, SEXP N, SEXP max_iter, SEXP no_masses)
{
  int npar, nx, ny, i, j, iter_max;
  double p1, p2, m;
  double *yp, *errp, *parp, *precp, *statep, *valuep;
  double chi_square, red_chi_square;
  int dof;
  int * xp, *Thalfp, *trp, *Np, *mip, *nmp;
  SEXP state;
  const gsl_multifit_fdfsolver_type *solver_type =
    gsl_multifit_fdfsolver_lmsder;
  /*  Allocate the solver. */
  gsl_multifit_fdfsolver *solver;
  /*  Initialize the data structure. */
  struct data data_struct;
  gsl_multifit_function_fdf function_fdf;
  gsl_matrix *covar;
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
  PROTECT(no_masses = AS_INTEGER(no_masses));

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
  nmp = INTEGER_POINTER(no_masses);
  iter_max = mip[0];

  npar = LENGTH(par);
  nx = LENGTH(x);
  ny = LENGTH(y);

  assert(npar == nmp[0]*(Np[0]+1));
  PROTECT(state = NEW_NUMERIC(5));
  statep = NUMERIC_POINTER(state);

/*   PROTECT(gradient = allocMatrix(REALSXP, npar, npar)); */
  
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
  data_struct.no_masses = nmp[0];

  // The ansatz.
  function_fdf.f = &exp_f;
  function_fdf.df = &exp_df;
  function_fdf.fdf = &exp_fdf;
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
    
    if(status) {
      break;
    }

    status = gsl_multifit_test_delta(solver->dx, solver->x,
				     precp[0], precp[1]);

  }
  while(status == GSL_CONTINUE && iter < iter_max);
#ifdef _DEBUG
  fprintf(stderr, "\n");
#endif
  
  // *****
  
  
  // Compute the covariance matrix.

  covar = gsl_matrix_alloc(npar, npar);
  gsl_multifit_covar(solver->J, 0.0, covar);

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

  UNPROTECT(12);
  return(state);
}
