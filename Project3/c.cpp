#include "lib.h"
#include <algorithm>
//#include <mpi.h>
#include <random>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <functional>
using namespace std;
double int_func(double x1, double x2, double y1, double y2, double z1, double z2){
  double exponential_part=exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2))); // Exponential part of the function
  double denom=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)); // Denominator part
  double returnval;
  //cout << denom << endl;
  if (fabs(denom)<= 1e-15){ //If denominator is zero, it is omitted
    returnval=0;
  }
  else{
    returnval= exponential_part/denom; //Final value of function
  }
  return returnval;
}
int main(int argc, char** argv){
  int N;
  if(argc>=2){
    N=atoi(argv[1]); //n needs to be stated
  }
  else{
    cout << "You need to state a number N" << endl;
    exit(1);
  }
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> RnG(0.0,1.0);
  double lambda=2; //length of the box
  double int_mc=0.; double variance=0.;
  double sum_sigma=0.; double sum_mc;
  double x[6];
  double fx;
  double jacobi_det=pow(2*lambda,6);

  for (int i=0;i<N;i++){
    for(int j=0;j<6;j++){
      x[j]=-lambda+2*lambda*RnG(gen); // makes x[j] a random number between -lambda and lamnbda
    }
    fx=int_func(x[0],x[1],x[2],x[3],x[4],x[5]);
    int_mc+=fx;
    sum_sigma+=fx*fx;
  }
  int_mc=int_mc/N;
  sum_sigma=sum_sigma/N;
  variance=sum_sigma-int_mc*int_mc;
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << " Monte carlo result= " << setw(10) << setprecision(8) << jacobi_det*int_mc;
  cout << " Sigma= " << setw(10) << setprecision(8) << jacobi_det*sqrt(variance/((double) N )) << endl;
}
