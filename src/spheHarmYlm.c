#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>

//Calculating the spherical harmonics Y_lm(theta,phi).
//x=cos(theta);
//theta is the polar angle, phi is the azimuthal angle.
complex spheHarm(int l, int m, double x, double phi)
{
  gsl_sf_result result;
  double complex Ylm;

  int status = gsl_sf_legendre_sphPlm_e(l, m, x, &result);
  if(status == 0){
    Ylm = pow(-1,m)*result.val*cos(m*phi) + I*pow(-1,m)*result.val*sin(m*phi);
    
    //printf("Y_%d%d = %.12f\n", l, m, result.val);
    return Ylm;
  }
  else 
    {
      printf("Spherical Harmonics calculating wrong!");
      exit(1);
    }
}

//Calculating the azimutal angle with in the sperical coordinate
//system from the Descartes coordiante x and y.
double azimutalAngle(double x, double y)
{
  double phi = 0.;
  if(x==0 && y==0) {
    phi = 0;
  }
  else if(x==0 && y>0){
    phi = 90;
  }
  else if(x==0 && y<0){
    phi = 270;
  }
  else if(x>0 && y<0){
    phi = atan(y/x)*180/M_PI + 360;
  }
  else if(x>0 && y>=0){
    phi = atan(y/x)*180/M_PI;
  }
  else if (x<0)
    phi = atan(y/x)*180/M_PI + 180;
  return phi;
}


