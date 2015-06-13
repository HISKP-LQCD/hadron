#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>

#include "zetaFunc.h"

double complex firstPart(const double Tolerance, const int l, const int m, const double * dVec, const double gamma, const double Lamda, const double qSqur, const int verbose, int * const rstatus)
{
  double complex firstTerms=0+I*0, pmodeSum=0+I*0, firstPartSum=0+I*0;
  double cosPolarAngle,azAngle;
  double rVecMod;
  int npmode[2] = {40, 72};
  int i_npmode=1;
  
  double error = 1.0;
  int pmodeSqur;
  
  int n1, n2, n3;
  double dModSqur=dVec[0]*dVec[0]+dVec[1]*dVec[1]+dVec[2]*dVec[2];
  
  int * degnrtDOF = NULL;
  int * arrayPmode= NULL;
  
  int niter = 0;
  int genReturn;
  get_npmode(npmode, i_npmode);
  genReturn = gen_points_array(&degnrtDOF, &arrayPmode, npmode[0], npmode[1]);

  if(genReturn != 0){
    printf("Generating the points wrong in firstPart!");
    exit(-1);
  }
  
  if(verbose){
    for(int i = 0; i < npmode[0]; i++){
      if(degnrtDOF[i] == 0)
	printf("pmodeSqur=%d has no corresponding points.\n", i);
      else
	printf("pmodeSqur=%d have %d degenaration DOF.\n", i, degnrtDOF[i]);
    }
  }
  
  pmodeSqur = 0;
  
  while(error > Tolerance){
    
    if(pmodeSqur > npmode[0]){
      i_npmode++;
      fprintf(stderr, "increased i_npmode to %d in thirdPart\n", i_npmode);
      get_npmode(npmode, i_npmode);
      genReturn = gen_points_array(&degnrtDOF, &arrayPmode, npmode[0], npmode[1]);

      if(genReturn != 0){
	printf("Generating the points wrong in firstPart!");
	exit(-1);
      }

      if(i_npmode > 4) {
	fprintf(stderr, "NPmode and DimMax need to be larger than available in firstPart! Aborting...!\n");
	exit(-1);
      }
    }
    
    pmodeSum = 0+I*0;
    
    //These pmodes has no contribution to the points sets.
    if(degnrtDOF[pmodeSqur] == 0){
      pmodeSqur += 1;
      continue;
    }
    
    for(int i = 0; i < degnrtDOF[pmodeSqur]; i++){

      n1 = arrayPmode[pmodeSqur*npmode[1]*3 + i*3 + 0];
      n2 = arrayPmode[pmodeSqur*npmode[1]*3 + i*3 + 1];
      n3 = arrayPmode[pmodeSqur*npmode[1]*3 + i*3 + 2];
      
      if(verbose)
	printf("%3d %3d %3d\n", n1, n2, n3);
      
      double nSqur = (double)n1*n1+(double)n2*n2+(double)n3*n3; 
      
      if( fabs(dModSqur) < DBL_EPSILON  ){
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
        if(fabs(1-fabs(cosPolarAngle)) > DBL_EPSILON*10) fprintf(stderr, "Warning, cosPolarAngle > 1 by %e in firstPart.c\n", 1-fabs(cosPolarAngle));
	if(fabs(1-fabs(cosPolarAngle)) > DBL_EPSILON*100) {
	  fprintf(stderr, "rVecMod: %e, dModSqur %e n: %d %d %d nSqur %e\n", rVecMod, dModSqur, n1, n2, n3, nSqur);
	  *rstatus = 13;
	  return(firstPartSum);
	}
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
      printf("pmode %d error: %.16f\n\n", pmodeSqur , error);
    
    // if the result is still zero after 4 iterations it is assumed to stay zero
    if (cabs(firstPartSum) < DBL_EPSILON && niter > 4)
      break;
    pmodeSqur += 1;
    ++niter;
    
  }//end of while.
  
  return firstPartSum;
}
				


