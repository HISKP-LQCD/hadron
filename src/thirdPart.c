#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include "zetaFunc.h"

//NPmode should be bigger than the biggest pmodeSqur we may meet in the iteration,
//DimMax is the biggest degenation degree within the range of NPmode
//Usually we should set big enough value for these two parameters, 
//or the program will crash.
//For precision of 1e-8, (NPmode=40, DimMAX=72) is good enough.
//
//Other selections may be 
//(NPmode=70,DimMAX=96),
//(NPmode=100,DimMAX=120),
//(NPmode=145,DimMAX=168)
//These two parameters are set in the header file zetaFunc.h

double complex thirdPart(const double Tolerance, const int l, const int m, double * const dVec, 
                         const double gamma, const double Lamda, const double qSqur, const int verbose,
                         int * const rstatus)
{
  double complex thirdTerms = 0+I*0, pmodeSum=0+I*0, thirdPartSum = 0+I*0;
  double cosPolarAngle=0,azAngle=0;
  double wVecMod=0;

  double error = 1.0;
  int pmodeSqur;

  int n1,n2,n3;
  double dModSqur = dVec[0]*dVec[0]+dVec[1]*dVec[1]+dVec[2]*dVec[2];

  int s1=0, s2=0;

  int * degnrtDOF = NULL;
  int * arrayPmode= NULL;

	
  int niter = 0;
  int genReturn;
  genReturn = gen_points_array(&degnrtDOF, &arrayPmode, NPmode, DimMAX);
	
  if(genReturn != 0){
    printf("Generating the points wrong!");
    exit(-1);
  }
	
  //From the formula in the paper w!=0,so we start from pmodeSqur=1.
  pmodeSqur = 1;
	
  while(error > Tolerance){

    if(pmodeSqur > NPmode){
      printf("The tolerance requisition has exceeded the pmodeSqur upper limit set by NPmode!\nPlease increase the macro definition of NPmode and DimMax in the head file zetaFunc.h!\n");
      exit(-1);
    }

    pmodeSum = 0+I*0;

    //From the formula in the paper: w!=0
    if(pmodeSqur == 0)
      continue;
	
    //These pmodes has no contribution to the points sets.
    if(degnrtDOF[pmodeSqur] == 0){
      pmodeSqur += 1;
      continue;
    }

    for(int i=0; i<degnrtDOF[pmodeSqur]; i++){

      //n1,n2,n3 stands for the components of vector w.
      n1=arrayPmode[pmodeSqur*DimMAX*3 + i*3 + 0];
      n2=arrayPmode[pmodeSqur*DimMAX*3 + i*3 + 1];
      n3=arrayPmode[pmodeSqur*DimMAX*3 + i*3 + 2];

      if(verbose)
	printf("%3d %3d %3d\n", n1, n2, n3);

      //nVec needed by the integrand in the exponential
      int nVec[3] = {n1, n2 ,n3};
      double nSqur = n1*n1 + n2*n2 + n3*n3;

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
        
        
        if(fabs(cosPolarAngle) > 1) {
          // cosPolarAngle must not become larger than 1 
          // we check for this here and drop a warning if unexpectedly large
          if(fabs(1-cosPolarAngle) > DBL_EPSILON) fprintf(stderr, "Warning, cosPolarAngle > 1 by %e\n", 1-cosPolarAngle);
          cosPolarAngle /= fabs(cosPolarAngle);
        }
        
        
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

      //Add every term within the same pmode into pmodeSum
      pmodeSum += thirdTerms;
    }//end of pmode loop.
      
    thirdPartSum += pmodeSum;
    //Both pmodeSum and firstPartSum are complex numbers,
    //cabs take the mode of these variables.
    // only calculate new error if firstPartSum != 0.
    if (cabs(thirdPartSum) > DBL_EPSILON)
      error = cabs(pmodeSum) / cabs(thirdPartSum);
    
    if(verbose)
      printf("pmode=%d error: %.16f\n\n",pmodeSqur , error);
    
    // if the result is still zero after 4 iterations it is assumed to stay zero
    if (cabs(thirdPartSum) < DBL_EPSILON && niter > 4)
      break;
    pmodeSqur += 1;
    ++niter;
    
  }//end of while.
	
  *rstatus = s1 + s2;

  return thirdPartSum;
}



