#ifndef ZETA_FUNC_H
#define ZETA_FUNC_H

//Main func:

int read_i(char *src,int *data);
int read_d(char *src,double *data);
int read_a(char *src,int *data);


//First part:

double azimutalAngle(double x, double y);

double spheHarm(int l, int m, double x, double phi, int * rstatus);

double complex firstPart(int N, int l, int m, int * dVec, double gamma, double Lamda, double qSqur, int * rstatus);
//Second part:

double integrandPart2(double t, void * params);

double sndInteFunc(double Lamda, double qSqur, int * rstatus);

double complex secondPart(int l, double gamma, double Lamda, double qSqur, int * rstatus);
//Third part:

double integrandPart3(double t, void * params);

double trdInteFunc(double Lamda, int * dVec, int l, double qSqur, int *nVec, double gamma, int * rstatus);

double complex thirdPart(int N, int l, int m, int * dVec, double gamma, double Lamda, double qSqur, int * rstatus);

#endif
