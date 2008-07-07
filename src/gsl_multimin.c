#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
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

double exp_fsq(const gsl_vector * x, void *data)
{
  int no_masses = ((struct data *)data)->no_masses;
  int N = ((struct data *)data)->N;
  int tr = ((struct data *)data)->tr;
  double Time = 2.*((struct data *)data)->Thalf;
  double *t = ((struct data *)data)->x;
  double *y = ((struct data *)data)->y;
  double *err = ((struct data *) data)->err;
  double p[6][6], m[6];
  double sum = 0.;

  size_t i, j, k=0, id0, id1;
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
      sum += ((Y - y[i+j*tr])/err[i+j*tr])*((Y - y[i+j*tr])/err[i+j*tr]);
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
	sum += ((Y - y[i+j*tr])/err[i+j*tr])*((Y - y[i+j*tr])/err[i+j*tr]);
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
	sum += ((Y - y[i+j*tr])/err[i+j*tr])*((Y - y[i+j*tr])/err[i+j*tr]);
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
	sum += ((Y - y[i+j*tr])/err[i+j*tr])*((Y - y[i+j*tr])/err[i+j*tr]);
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
	sum += ((Y - y[i+j*tr])/err[i+j*tr])*((Y - y[i+j*tr])/err[i+j*tr]);
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
	sum += ((Y - y[i+j*tr])/err[i+j*tr])*((Y - y[i+j*tr])/err[i+j*tr]);
      }
    }
  }
  return(sum);
}



void exp_dfsq (const gsl_vector * x, void *data, 
		  gsl_vector * g)
{
  int no_masses = ((struct data *)data)->no_masses;
  int tr = ((struct data *)data)->tr;
  int N = ((struct data *)data)->N;
  double Time = 2.*((struct data *)data)->Thalf;
  double *t = ((struct data *)data)->x;
  double *y = ((struct data *)data)->y;
  double *err = ((struct data *) data)->err;
  int npar = N+1, id0, id1;
  double p[6][6], m[6];
  size_t i, j, k=0, kludge = 0;
  double Y = 0., dY=0., sign = 1., c;
  double *res;

  res = (double*)malloc(no_masses*npar*sizeof(double));

  for(i = 0; i < no_masses; i++) {
    for(j = 0; j < N; j++) {
      p[j][i] = gsl_vector_get (x, j + i*(N+1));
    }
    m[i]  = gsl_vector_get (x, N + i*(N+1));
  }
  for(i = 0; i < npar; i++) {
    res[i] = 0.;
  }

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
      for(k = 0; k < no_masses; k++) {
	res[k*(N+1)+id0] += (Y- y[i+j*tr])*p[id1][k]*
	  (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	res[k*(N+1)+id1] += (Y- y[i+j*tr])*p[id0][k]*
	  (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	res[N + k*(N+1)] -= (Y- y[i+j*tr])*c/err[i+j*tr]*
	  ((Time-t[i])*exp(-m[k]*(Time-t[i])) + sign*t[i]*exp(-m[k]*t[i]))/err[i+j*tr];
      }
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
	for(k = 0; k < no_masses; k++) {
	  res[k*(N+1)+id0] += (Y- y[i+j*tr])*p[id1][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[k*(N+1)+id1] += (Y- y[i+j*tr])*p[id0][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[N + k*(N+1)] -= (Y- y[i+j*tr])*c/err[i+j*tr]*
	    ((Time-t[i])*exp(-m[k]*(Time-t[i])) + sign*t[i]*exp(-m[k]*t[i]))/err[i+j*tr];
	}
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
	for(k = 0; k < no_masses; k++) {
	  res[k*(N+1)+id0] += (Y- y[i+j*tr])*p[id1][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[k*(N+1)+id1] += (Y- y[i+j*tr])*p[id0][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[N + k*(N+1)] -= (Y- y[i+j*tr])*c/err[i+j*tr]*
	    ((Time-t[i])*exp(-m[k]*(Time-t[i])) + sign*t[i]*exp(-m[k]*t[i]))/err[i+j*tr];
	}
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
	for(k = 0; k < no_masses; k++) {
	  res[k*(N+1)+id0] += (Y- y[i+j*tr])*p[id1][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[k*(N+1)+id1] += (Y- y[i+j*tr])*p[id0][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[N + k*(N+1)] -= (Y- y[i+j*tr])*c/err[i+j*tr]*
	    ((Time-t[i])*exp(-m[k]*(Time-t[i])) + sign*t[i]*exp(-m[k]*t[i]))/err[i+j*tr];
	}
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
	for(k = 0; k < no_masses; k++) {
	  res[k*(N+1)+id0] += (Y- y[i+j*tr])*p[id1][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[k*(N+1)+id1] += (Y- y[i+j*tr])*p[id0][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[N + k*(N+1)] -= (Y- y[i+j*tr])*c/err[i+j*tr]*
	    ((Time-t[i])*exp(-m[k]*(Time-t[i])) + sign*t[i]*exp(-m[k]*t[i]))/err[i+j*tr];
	}
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
	for(k = 0; k < no_masses; k++) {
	  res[k*(N+1)+id0] += (Y- y[i+j*tr])*p[id1][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[k*(N+1)+id1] += (Y- y[i+j*tr])*p[id0][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[N + k*(N+1)] -= (Y- y[i+j*tr])*c/err[i+j*tr]*
	    ((Time-t[i])*exp(-m[k]*(Time-t[i])) + sign*t[i]*exp(-m[k]*t[i]))/err[i+j*tr];
	}
      }
    }
  }
  
  for(i = 0; i < npar; i++) {
    gsl_vector_set(g, i, res[i]);
  }
  free(res);
  return;
}

void exp_fdfsq (const gsl_vector * x, void *data, double *f, 
		gsl_vector * g)
{
  int no_masses = ((struct data *)data)->no_masses;
  int tr = ((struct data *)data)->tr;
  int N = ((struct data *)data)->N;
  double Time = 2.*((struct data *)data)->Thalf;
  double *t = ((struct data *)data)->x;
  double *y = ((struct data *)data)->y;
  double *err = ((struct data *) data)->err;
  int npar = N+1, id0, id1;
  double p[6][6], m[6];
  size_t i, j, k=0, kludge = 0;
  double Y = 0., dY=0., sign = 1., c, sum=0.;
  double *res;

  res = (double*)malloc(no_masses*npar*sizeof(double));

  for(i = 0; i < no_masses; i++) {
    for(j = 0; j < N; j++) {
      p[j][i] = gsl_vector_get (x, j + i*(N+1));
    }
    m[i]  = gsl_vector_get (x, N + i*(N+1));
  }
  for(i = 0; i < npar; i++) {
    res[i] = 0.;
  }

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
      sum += ((Y - y[i+j*tr])/err[i+j*tr])*((Y - y[i+j*tr])/err[i+j*tr]);
      for(k = 0; k < no_masses; k++) {
	res[k*(N+1)+id0] += (Y- y[i+j*tr])*p[id1][k]*
	  (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	res[k*(N+1)+id1] += (Y- y[i+j*tr])*p[id0][k]*
	  (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	res[N + k*(N+1)] -= (Y- y[i+j*tr])*c/err[i+j*tr]*
	  ((Time-t[i])*exp(-m[k]*(Time-t[i])) + sign*t[i]*exp(-m[k]*t[i]))/err[i+j*tr];
      }
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
	sum += ((Y - y[i+j*tr])/err[i+j*tr])*((Y - y[i+j*tr])/err[i+j*tr]);
	for(k = 0; k < no_masses; k++) {
	  res[k*(N+1)+id0] += (Y- y[i+j*tr])*p[id1][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[k*(N+1)+id1] += (Y- y[i+j*tr])*p[id0][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[N + k*(N+1)] -= (Y- y[i+j*tr])*c/err[i+j*tr]*
	    ((Time-t[i])*exp(-m[k]*(Time-t[i])) + sign*t[i]*exp(-m[k]*t[i]))/err[i+j*tr];
	}
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
	sum += ((Y - y[i+j*tr])/err[i+j*tr])*((Y - y[i+j*tr])/err[i+j*tr]);
	for(k = 0; k < no_masses; k++) {
	  res[k*(N+1)+id0] += (Y- y[i+j*tr])*p[id1][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[k*(N+1)+id1] += (Y- y[i+j*tr])*p[id0][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[N + k*(N+1)] -= (Y- y[i+j*tr])*c/err[i+j*tr]*
	    ((Time-t[i])*exp(-m[k]*(Time-t[i])) + sign*t[i]*exp(-m[k]*t[i]))/err[i+j*tr];
	}
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
	sum += ((Y - y[i+j*tr])/err[i+j*tr])*((Y - y[i+j*tr])/err[i+j*tr]);
	for(k = 0; k < no_masses; k++) {
	  res[k*(N+1)+id0] += (Y- y[i+j*tr])*p[id1][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[k*(N+1)+id1] += (Y- y[i+j*tr])*p[id0][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[N + k*(N+1)] -= (Y- y[i+j*tr])*c/err[i+j*tr]*
	    ((Time-t[i])*exp(-m[k]*(Time-t[i])) + sign*t[i]*exp(-m[k]*t[i]))/err[i+j*tr];
	}
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
	sum += ((Y - y[i+j*tr])/err[i+j*tr])*((Y - y[i+j*tr])/err[i+j*tr]);
	for(k = 0; k < no_masses; k++) {
	  res[k*(N+1)+id0] += (Y- y[i+j*tr])*p[id1][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[k*(N+1)+id1] += (Y- y[i+j*tr])*p[id0][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[N + k*(N+1)] -= (Y- y[i+j*tr])*c/err[i+j*tr]*
	    ((Time-t[i])*exp(-m[k]*(Time-t[i])) + sign*t[i]*exp(-m[k]*t[i]))/err[i+j*tr];
	}
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
	sum += ((Y - y[i+j*tr])/err[i+j*tr])*((Y - y[i+j*tr])/err[i+j*tr]);
	for(k = 0; k < no_masses; k++) {
	  res[k*(N+1)+id0] += (Y- y[i+j*tr])*p[id1][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[k*(N+1)+id1] += (Y- y[i+j*tr])*p[id0][k]*
	    (exp(-m[k]*(Time-t[i])) + sign*exp(-m[k]*t[i]))/err[i+j*tr]/err[i+j*tr];
	  res[N + k*(N+1)] -= (Y- y[i+j*tr])*c/err[i+j*tr]*
	    ((Time-t[i])*exp(-m[k]*(Time-t[i])) + sign*t[i]*exp(-m[k]*t[i]))/err[i+j*tr];
	}
      }
    }
  }
  
  for(i = 0; i < npar; i++) {
    gsl_vector_set(g, i, res[i]);
  }
  (*f) = sum;
  free(res);
  return;
}

/* void exp_fdfsq (const gsl_vector * x, void *data, double *f, gsl_vector * g) */
/* { */
/*   (*f) = exp_fsq (x, data); */
/*   exp_dfsq (x, data, g); */
  
/*   return; */
/* } */


void Print_State_Mass_Fit_Helper_2(int iter, gsl_multimin_fdfminimizer *minimizer)
{
    printf ("%5d %.5f %.5f %.5f %10.5f\n", iter,
	    gsl_vector_get (minimizer->x, 0), 
	    gsl_vector_get (minimizer->x, 1), 
	    gsl_vector_get (minimizer->x, 2), 
	    minimizer->f);

}

SEXP multimin_cor(SEXP par, SEXP Thalf, SEXP x, SEXP y, SEXP err, SEXP tr, 
		  SEXP prec, SEXP N, SEXP max_iter, SEXP no_masses)
{
  int npar, nx, ny, i, j, iter_max;
  double p1, p2, m;
  double *yp, *errp, *parp, *precp, *statep;
  double chi_square, red_chi_square;
  int dof;
  int * xp, *Thalfp, *trp, *Np, *mip, *nmp;
  SEXP state;
  const gsl_multimin_fdfminimizer_type *minimizer_type =
  gsl_multimin_fdfminimizer_vector_bfgs;
/*
  gsl_multimin_fdfminimizer_conjugate_fr;
  gsl_multimin_fdfminimizer_conjugate_pr;
  gsl_multimin_fdfminimizer_vector_bfgs2;
  gsl_multimin_fdfminimizer_vector_bfgs;
  gsl_multimin_fdfminimizer_steepest_descent;
  gsl_multimin_fminimizer_nmsimplex;
*/
  /*  Allocate the solver. */
  gsl_multimin_fdfminimizer *minimizer;
  /*  Initialize the data structure. */
  struct data data_struct;
  gsl_multimin_function_fdf function_fdf;
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
  mip = INTEGER_POINTER(max_iter);
  nmp = INTEGER_POINTER(no_masses);
  iter_max = mip[0];

  npar = LENGTH(par);
  nx = LENGTH(x);
  ny = LENGTH(y);

  assert(npar == nmp[0]*(Np[0]+1));
  PROTECT(state = NEW_NUMERIC(5+npar));
  statep = NUMERIC_POINTER(state);

  
  if(Np[0] == 2) no_points = 3*trp[0];
  if(Np[0] == 4) no_points = 10*trp[0];
  if(Np[0] == 6) no_points = 21*trp[0];

  minimizer = gsl_multimin_fdfminimizer_alloc(minimizer_type, npar);

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
  function_fdf.f = &exp_fsq;
  function_fdf.df = &exp_dfsq;
  function_fdf.fdf = &exp_fdfsq;
  function_fdf.n = npar;
  function_fdf.params = &data_struct;

  for(i = 0; i < npar; i++) {
    para_initial[i] = parp[i];
  }
  para_initial_ = gsl_vector_view_array(para_initial, npar);

  gsl_multimin_fdfminimizer_set(minimizer, &function_fdf, &para_initial_.vector, 0.001, 1.e-4);

  // Perform the fit.
  // Print the initial state.
#ifdef _DEBUG
  Print_State_Mass_Fit_Helper_2(iter, minimizer);
#endif

  do {
    iter++;
    
    /*  Do a solver iteration. */
    status = gsl_multimin_fdfminimizer_iterate(minimizer);
#ifdef _DEBUG
    fprintf(stderr, "status = %s.\n", gsl_strerror(status));
    Print_State_Mass_Fit_Helper_2(iter, minimizer);
#endif
    
    if(status) {
      break;
    }

    status = gsl_multimin_test_gradient(minimizer->gradient, precp[1]);
#ifdef _DEBUG
    if (status == GSL_SUCCESS)
      printf ("Minimum found at:\n");

    Print_State_Mass_Fit_Helper_2(iter, minimizer);
#endif

  }
  while(status == GSL_CONTINUE && iter < iter_max);
#ifdef _DEBUG
  fprintf(stderr, "\n");
#endif
  
  // *****
  dof = no_points - npar;
  chi_square = minimizer->f;
  red_chi_square = chi_square / (double)dof;

  for(i = 0; i < npar; i++) {
    statep[5+i] =  gsl_vector_get(minimizer->x, i);
  }

  statep[0] = red_chi_square;
  statep[1] = minimizer->f;
  statep[2] = (double)iter;
  statep[3] = (double)dof;
  statep[4] = (double)status;

  
  gsl_multimin_fdfminimizer_free(minimizer);

  free(data_struct.x);
  free(data_struct.y);
  free(data_struct.err);
  free(para_initial);

  UNPROTECT(11);
  return(state);
}
