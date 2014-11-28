#ifndef ZETA_FUNC_H
#define ZETA_FUNC_H

//Main func:

int read_i(char *src,int *data);
int read_d(char *src,double *data);
int read_a(char *src,double *data);


//First part:

double azimutalAngle(const double x, const double y);

complex double spheHarm(const int l, const int m, const double x, const double phi, int * const rstatus);

double complex firstPart(int N, int l, int m, double * dVec, double gamma, double Lamda, double qSqur, int * rstatus);
//Second part:

double integrandPart2(const double t, void * const params);

double sndInteFunc(const double Lamda, const double qSqur, int * const rstatus);

double complex secondPart(const int l, const double gamma, const double Lamda, const double qSqur, int * const rstatus);
//Third part:

double integrandPart3(const double t, void * const params);

double trdInteFunc(const double Lamda, double * const dVec, const int l, const double qSqur, int * const nVec, const double gamma, int * const rstatus);

double complex thirdPart(const int N, const int l, const int m, double * const dVec, const double gamma, const double Lamda, const double qSqur, int * const rstatus);

#endif
