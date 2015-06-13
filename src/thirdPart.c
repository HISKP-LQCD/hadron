#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include "zetaFunc.h"

double complex thirdPart(const double Tolerance, const int l, const int m, double * const dVec, 
                         const double gamma, const double Lamda, const double qSqur, const int verbose,
                         int * const rstatus)
{
  double complex thirdTerms = 0+I*0, pmodeSum=0+I*0, thirdPartSum = 0+I*0;
  double cosPolarAngle=0,azAngle=0;
  double wVecMod=0;
  int npmode[2] = {40, 72};
  int i_npmode=1;

  double error = 1.0;
  int pmodeSqur;

  int n1,n2,n3;
  double dModSqur = dVec[0]*dVec[0]+dVec[1]*dVec[1]+dVec[2]*dVec[2];

  int s1=0, s2=0;

  int * degnrtDOF = NULL;
  int * arrayPmode= NULL;

	
  int niter = 0;
  int genReturn;
  get_npmode(npmode, i_npmode);
  genReturn = gen_points_array(&degnrtDOF, &arrayPmode, npmode[0], npmode[1]);
	
  if(genReturn != 0){
    printf("Generating the points wrong in thirdPart!");
    exit(-1);
  }
	
  //From the formula in the paper w!=0,so we start from pmodeSqur=1.
  pmodeSqur = 1;
	
  while(error > Tolerance){

    if(pmodeSqur > npmode[0]){
      i_npmode++;
      fprintf(stderr, "increased i_npmode to %ud in thirdPart\n", i_npmode);
      get_npmode(npmode, i_npmode);
      genReturn = gen_points_array(&degnrtDOF, &arrayPmode, npmode[0], npmode[1]);

      if(genReturn != 0){
	printf("Generating the points wrong in thirdPart!");
	exit(-1);
      }

      if(i_npmode > 4) {
	fprintf(stderr, "NPmode and DimMax need to be larger than available in thirdPart! Aborting...!\n");
	exit(-1);
      }
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

    for(int i = 0; i < degnrtDOF[pmodeSqur]; i++){

      //n1,n2,n3 stands for the components of vector w.
      n1 = arrayPmode[pmodeSqur*npmode[1]*3 + i*3 + 0];
      n2 = arrayPmode[pmodeSqur*npmode[1]*3 + i*3 + 1];
      n3 = arrayPmode[pmodeSqur*npmode[1]*3 + i*3 + 2];

      if(verbose)
	printf("%3d %3d %3d\n", n1, n2, n3);

      //nVec needed by the integrand in the exponential
      int nVec[3] = {n1, n2 ,n3};
      double nSqur = n1*n1 + n2*n2 + n3*n3;

      if( fabs(dModSqur) < DBL_EPSILON ){
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
          if(fabs(1-cosPolarAngle) > DBL_EPSILON*10) fprintf(stderr, "Warning, cosPolarAngle > 1 by %e in thirdPart\n", 1-cosPolarAngle);
	  if(fabs(1-fabs(cosPolarAngle)) > DBL_EPSILON*100) {
	    fprintf(stderr, "wVecMod: %e\n", wVecMod);
	    *rstatus = 13;
	    return(thirdPartSum);
	  }

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



