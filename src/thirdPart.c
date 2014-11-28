#include <stdio.h>
#include <math.h>
#include <complex.h>

#include <gsl/gsl_math.h>
#include "zetaFunc.h"

double complex thirdPart(int N, int l, int m, int * dVec, double gamma, double Lamda, double qSqur, int * rstatus)
{
  //	double qSqur;
  //	int N;
  //	printf("Please input the N and qSqur:\n");
  //	scanf("%d %lf", &N, &qSqur);
  //l,m is the parameter of the spherical harmonics
  //	int l,m;
  //	printf("Please input the l and m:\n");
  //	scanf("%d %d", &l, &m);
  //d is the boost-vector,0 for cubic group.
  //	int dVec[3]={0,0,0};
  //Lorenz factor
  //	double gamma = 1 ;
  //Lamda is the partition point of integration range (0, Inf).
  //	double Lamda = 1;

  int n1,n2,n3;
  int dModSqur = dVec[0]*dVec[0]+dVec[1]*dVec[1]+dVec[2]*dVec[2];

  double complex thirdTerms = 0+I*0, thirdPartSum = 0+I*0;
  double cosPolarAngle=0,azAngle=0;
  double wVecMod=0;
  int s1=0, s2=0;

  for(n1=-N; n1<=N; n1++)
    for(n2=-N; n2<=N; n2++)
      for(n3=-N; n3<=N; n3++){
        double nSqur = n1*n1 + n2*n2 + n3*n3;
        if(nSqur <= N*N && !(n1==0 && n2==0 && n3==0)){
          //nVec needed by the integrand in the exponential.
          int nVec[3] = {n1, n2 ,n3};
          if( dModSqur == 0 ){
            wVecMod = sqrt( nSqur );
            cosPolarAngle =  n3/wVecMod ;
            azAngle = azimutalAngle(n1,n2) ;
            thirdTerms = gamma * pow(wVecMod, l)
              * spheHarm(l, m, cosPolarAngle, azAngle, &s1)
              * trdInteFunc(Lamda, dVec, l, qSqur, nVec, gamma, &s2);
          }
          else{
            double wDotd = n1*dVec[0]+n2*dVec[1]+n3*dVec[2];
            //wVecMod stands for |\hat{gamma} * \vec{w}|
            wVecMod = sqrt((pow(gamma,2)-1) * pow(wDotd,2) /dModSqur
                           + nSqur);
            cosPolarAngle=((gamma-1)*wDotd*dVec[2]/dModSqur + n3)/wVecMod;
            azAngle = azimutalAngle((gamma-1)*wDotd*dVec[0]/dModSqur + n1 ,
                                    (gamma-1)*wDotd*dVec[1]/dModSqur + n2);
            thirdTerms = gamma * (cos(M_PI*wDotd) - I*sin(M_PI*wDotd))
              * pow(wVecMod, l) *spheHarm(l, m, cosPolarAngle, azAngle, &s1)
              * trdInteFunc(Lamda, dVec, l, qSqur, nVec, gamma, &s2);
          }
          if(s1 != 0 || s2 != 0) {
            *rstatus = s1 + s2;
            return(thirdPartSum);
          }
          thirdPartSum += thirdTerms;
        }
      }

  *rstatus = s1 + s2;
  //printf("SumN=%d,qSqur=%.4f,  thirdPartSum=%.24lf %+.24lfI.\n", N, qSqur, creal(thirdPartSum),cimag(thirdPartSum));
  return thirdPartSum;
}

				








	
