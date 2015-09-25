#ifndef ZETA_FUNC_H
#define ZETA_FUNC_H


//Main func:

int read_i(char *src, int *data);
int read_d(char *src, double *data);
int read_a(char *src, double *data);


void get_npmode(int * pair, const int i);
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



