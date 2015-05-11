#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>

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

double complex firstPart(const double Tolerance, const int l, const int m, const double * dVec, const double gamma, const double Lamda, const double qSqur, const int verbose, int * const rstatus)
{
  double complex firstTerms=0+I*0, pmodeSum=0+I*0, firstPartSum=0+I*0;
  double cosPolarAngle,azAngle;
  double rVecMod;
  
  double error = 1.0;
  int pmodeSqur;
  
  int n1,n2,n3;
  double dModSqur=dVec[0]*dVec[0]+dVec[1]*dVec[1]+dVec[2]*dVec[2];
  
  int * degnrtDOF = NULL;
  int * arrayPmode= NULL;
  
  int niter = 0;
  int genReturn;
  genReturn = gen_points_array(&degnrtDOF, &arrayPmode, NPmode, DimMAX);

  if(genReturn != 0){
    printf("Generating the points wrong!");
    exit(-1);
  }
  
  if(verbose){
    for(int i=0; i<NPmode; i++){
      if(degnrtDOF[i] == 0)
	printf("pmodeSqur=%d has no corresponding points.\n", i);
      else
	printf("pmodeSqur=%d have %d degenaration DOF.\n", i, degnrtDOF[i]);
    }
  }
  
  pmodeSqur = 0;
  
  while(error > Tolerance){
    
    if(pmodeSqur > NPmode){
      printf("The tolerance requisition has exceeded the pmodeSqur upper limit set by NPmode!\nPlease increase the macro definition of NPmode and DimMax in the head file zetaFunc.h!\n");
      exit(-1);
    }
    
    pmodeSum = 0+I*0;
    
    //These pmodes has no contribution to the points sets.
    if(degnrtDOF[pmodeSqur] == 0){
      pmodeSqur += 1;
      continue;
    }
    
    for(int i=0; i<degnrtDOF[pmodeSqur]; i++){

      n1=arrayPmode[pmodeSqur*DimMAX*3 + i*3 + 0];
      n2=arrayPmode[pmodeSqur*DimMAX*3 + i*3 + 1];
      n3=arrayPmode[pmodeSqur*DimMAX*3 + i*3 + 2];
      
      if(verbose)
	printf("%3d %3d %3d\n", n1, n2, n3);
      
      double nSqur = n1*n1+n2*n2+n3*n3; 
      
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
        rVecMod = sqrt( (pow(nDotd,2.0)/dModSqur + dModSqur/4.0 - nDotd)/pow(gamma,2.0)
                        + nSqur  - pow(nDotd,2.0)/dModSqur);
        cosPolarAngle = ((nDotd*dVec[2]/dModSqur-dVec[2]/2.0)/gamma + (n3-nDotd*dVec[2]/dModSqur))
          / rVecMod;
        azAngle=azimutalAngle((nDotd*dVec[0]/dModSqur-dVec[0]/2.0)/gamma+(n1-nDotd*dVec[0]/dModSqur),
                              (nDotd*dVec[1]/dModSqur-dVec[1]/2.0)/gamma + (n2-nDotd*dVec[1]/dModSqur));
      }

      if(fabs(cosPolarAngle) > 1) {
        // cosPolarAngle must not become larger than 1 
        // we check for this here and drop a warning if unexpectedly large
        if(fabs(1-fabs(cosPolarAngle)) > DBL_EPSILON*10) fprintf(stderr, "Warning, cosPolarAngle > 1 by %e\n", 1-fabs(cosPolarAngle));
        cosPolarAngle /= fabs(cosPolarAngle);
      }
      
      firstTerms = exp(-Lamda*(pow(rVecMod,2.0)-qSqur)) * pow(rVecMod,l)
        * spheHarm(l, m, cosPolarAngle, azAngle, rstatus)
        / (pow(rVecMod,2.0) - qSqur);

      if(*rstatus != 0) {
	fprintf(stderr, "spheHarm produced error code \"%s\"\n", gsl_strerror(*rstatus)); 
	return(firstPartSum);
      }

      //Add every term within the same pmode into pmodeSum
      pmodeSum += firstTerms;

    }//end of pmode loop
    
    firstPartSum += pmodeSum;
    //Both pmodeSum and firstPartSum are complex numbers,
    //cabs take the mode of these variables.
    // only calculate new error if firstPartSum != 0.
    if (cabs(firstPartSum) > DBL_EPSILON)
      error = cabs(pmodeSum) / cabs(firstPartSum);
    
    if(verbose)
      printf("pmode%d error: %.16f\n\n",pmodeSqur , error);
    
    // if the result is still zero after 4 iterations it is assumed to stay zero
    if (cabs(firstPartSum) < DBL_EPSILON && niter > 4)
      break;
    pmodeSqur += 1;
    ++niter;
    
  }//end of while.
  
  return firstPartSum;
}
				


