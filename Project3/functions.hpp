#include <math.h>
#include <fstream>
#include <iostream>
#define EPS 3.0e-14
#define MAXIT 10
using namespace std;
void gauss_laguerre(double *x, double *w, int n, double alf);
double gammln( double xx);
double int_func_montecarlo(double r1, double r2, double theta1, double theta2, double phi1, double phi2);
double int_func(double x1, double x2, double y1, double y2, double z1, double z2);
double int_func_polar_gaulag(double r1, double r2, double theta1, double theta2, double phi1, double phi2);
double int_func_nor1r2(double r1, double r2, double theta1, double theta2, double phi1, double phi2);
long long int * getParallelizationCoefficients(long long int N,int mynum, int totNum, int loopings);
bool fileExists(const char *fileName);
