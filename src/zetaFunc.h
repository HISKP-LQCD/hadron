#ifndef ZETA_FUNC_H
#define ZETA_FUNC_H

//NPmode should be bigger than the biggest pmodeSqur we may meet in the iteration,
//DimMax is the biggest degenation degree within the range of NPmode
//Usually we should set big enough value for these two parameters or the program will crash.
//For precision of 1e-8, (NPmode=40, DimMAX=72) is good enough.
//Actually,by my test NPmode=40 can be selected up to precision 1e-11 for usual case.
//
//Other selections may be 
//(NPmode=70,DimMAX=96)
//NPmode=70 can be selected up to precision 1e-16 for usual case.
//
//(NPmode=100,DimMAX=120),
//(NPmode=145,DimMAX=168)
#define NPmode  70
#define DimMAX  96


//Main func:

int read_i(char *src, int *data);
int read_d(char *src, double *data);
int read_a(char *src, double *data);


int gen_points_array(int ** degnrtDOF, int ** arrayPmode, const int npmod, const int dimmax);
void pmode_free_arrays();

//First part:

double azimutalAngle(const double x, const double y);

complex double spheHarm(const int l, const int m, const double x, const double phi, int * const rstatus);

double complex firstPart(const double Tolerance, const int l, const int m, const double * dVec, const double gamma, const double Lamda, const double qSqur, const int verbose, int * const rstatus);
//Second part:

double integrandPart2(const double t, void * const params);

double sndInteFunc(const double Lamda, const double qSqur, int * const rstatus);

double complex secondPart(const int l, const double gamma, const double Lamda, const double qSqur, int * const rstatus);
//Third part:

double integrandPart3(const double t, void * const params);

double trdInteFunc(const double Lamda, double * const dVec, const int l, const double qSqur, int * const nVec, const double gamma, int * const rstatus);

double complex thirdPart(const double Tolerance, const int l, const int m, double * const dVec, const double gamma, const double Lamda, const double qSqur, const int verbose, int * const rstatus);

#endif



