#include <stdio.h>
#include <complex.h>
#include <math.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include "zetaFunc.h"
/***************************************************
 * This program calculates the second part of the 
 * zeta function:
 * 		delta_{L,0} *Y_00*gamma*PI^{3/2}
 * 		* [2*q^2* \int_0^Lamda exp{t*q^2}/sqrt{t} dt
 * 		   - 2*exp{Lamda*q^2}/sqrt{Lamda}
 ***************************************************/

//int main(void)
double complex secondPart(int l, double gamma, double Lamda, double qSqur)
{
//	double Lamda = 1;
//	double gamma = 1;
//	double qSqur;
//	int l;
//	printf("Please input the qSqur and l:\n");
//	scanf("%lf %d",&qSqur ,&l);

	double complex secondPartInt = 0+0*I;
	if(l != 0){
		secondPartInt = 0.0;
	}
	else{
			//problem1: should delete creal and change it to a complex number.
			secondPartInt = spheHarm(0, 0, 0, 0) * gamma * pow(M_PI,3.0/2.0) 
				* ( 2 * qSqur * sndInteFunc(Lamda, qSqur)
				   -2 * exp(Lamda * qSqur)/sqrt(Lamda));
	}
	//printf("Lamda=%lf,qSqur = %.4f\nsecondPartInt = %.24f \n", Lamda, qSqur, secondPartInt);
	return secondPartInt;
}
















