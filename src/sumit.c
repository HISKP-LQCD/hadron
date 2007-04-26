#include<stdlib.h>
#include<math.h>
#include<R.h>

/*
sum(((y[ii]
- par[1]*par[1]*(cv1))/err[ii])^2)
*/

void sumit(double *y, double *par1, double *par2, 
	   double *v, double *err, int * n,
	   double *res) {

  int i;
  double temp, p=(*par1)*(*par2);
  (*res) = 0.;
  for(i = 0; i < (*n); i++) {
    temp = (y[i] - p * v[i])/err[i];
    (*res) += temp*temp;
  }
}
