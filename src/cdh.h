#ifndef _FSCDH_H
#define _FSCDH_H

double g1(double x);
int g1array(double * x, double * res, const int n);
void fscdh(double rev, double aLamb1, double aLamb2, double aLamb3, double aLamb4,
	   double aF0_, double a_fm, int * L, double * ampiV, double *afpiV, const int n,
	   double * mpiFV, double * fpiFV, const int printit);
#endif
