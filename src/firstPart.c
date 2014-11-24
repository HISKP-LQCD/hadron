#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#include "zetaFunc.h"


complex double firstPart(int N, int l, int m, int * dVec, double gamma, double Lamda, double qSqur)
{
  //	double qSqur;
  //	int	N;
  //	printf("Please input the N and qSqur:\n");
  //	scanf("%d",&N);
  //	scanf("%lf",&qSqur);
  //l,m is the parameter of the spherical harmonics
  //	int l,m;
  //	printf("Please input the l and m:\n");
  //	scanf("%d %d", &l, &m);
  //d is the boost-vector,0 for cubic group.
  //	int dVec[3]={0,0,0};
  //Lorenz factor
  //	double gamma=1;
  //Lamda is the partition point of integration range (0, Inf).
  //	double Lamda=1;

	
  //n1,n2,n3 is the coordinate within the Descartes coordiante system
  int n1,n2,n3,nSqur;

  double complex firstTerms=0+I*0, firstPartSum=0+I*0;
  double cosPolarAngle,azAngle;
  double rVecMod;

  //	FILE * fp=NULL;
  //	fp=fopen("coordinate.txt","w");

  int	dModSqur=dVec[0]*dVec[0]+dVec[1]*dVec[1]+dVec[2]*dVec[2];
  for(n1=-N; n1<=N; n1++)
    for(n2=-N; n2<=N; n2++)
      for(n3=-N; n3<=N; n3++)
        {
          nSqur = n1*n1+n2*n2+n3*n3; 
          if( nSqur <= N*N ){

            //printf("Up to now no problem!\n");
            //fprintf(fp,"n1=%3d,n2=%3d,n3=%3d\n",n1,n2,n3);

            if( dModSqur == 0 ){
              rVecMod = sqrt( nSqur );
              if(n1==0 && n2==0 && n3==0){
                cosPolarAngle =0;
                azAngle =0;
              }
              else{
                cosPolarAngle =  n3/sqrt(nSqur) ;
                azAngle = azimutalAngle(n1,n2) ;
              }
            }
            else{
              double nDotd = n1*dVec[0]+n2*dVec[1]+n3*dVec[2];
              rVecMod = sqrt( (pow(nDotd,2)/dModSqur + dModSqur/4.0 - nDotd)/pow(gamma,2)
                              + nSqur + pow(nDotd,2)/dModSqur - 2*pow(nDotd,2)/dModSqur);
              cosPolarAngle = ((nDotd*dVec[2]/dModSqur-dVec[2]/2.0)/gamma + (n3-nDotd*dVec[2]/dModSqur))
                     / rVecMod;
              azAngle=azimutalAngle((nDotd*dVec[0]/dModSqur-dVec[0]/2.0)/gamma+(n1-nDotd*dVec[0]/dModSqur),
                                    (nDotd*dVec[1]/dModSqur-dVec[1]/2.0)/gamma + (n2-nDotd*dVec[1]/dModSqur));
            }

            if(fabs(cosPolarAngle) > 1) {
              // cosPolarAngle must not become larger than 1 
              // we check for this here and drop a warning if unexpectedly large
              if(fabs(1-cosPolarAngle) > DBL_EPSILON) fprintf(stderr, "Warning, cosPolarAngle > 1 by %e\n", 1-cosPolarAngle);
              cosPolarAngle /= fabs(cosPolarAngle);
            }

            firstTerms = exp(-Lamda*(pow(rVecMod,2)-qSqur)) * pow(rVecMod,l)
              * spheHarm(l, m, cosPolarAngle, azAngle)
              / (pow(rVecMod,2) - qSqur);
            //fprintf(fp,"n1=%d,n2=%d,n3=%d:  firstTerm=  %.24f %+.24fI\n",n1,n2,n3,creal(firstTerms), cimag(firstTerms));
          }
          firstPartSum += firstTerms;
        }
  //fclose(fp);
  //fp=NULL;
  //printf("SumN=%d,qSqur=%.4f,	 firstPartSum = %.24f %+.24fI.\n",N,qSqur,creal(firstPartSum), cimag(firstPartSum));
  return firstPartSum;
}
				






