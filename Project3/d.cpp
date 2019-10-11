#include "lib.h"
#include "functions.hpp"
#include <algorithm>
#include <random>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <functional>
#define PI 3.14159265358979

using namespace std;

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
  double int_mc=0.0; double variance=0.0;
  double sum_sigma=0.0; double sum_mc=0.0;
  double x[6];
  double fx=0;
  double jacobi_det=4.0*pow(PI,4)/16.0;//pow(2*PI,2)*pow(PI,2)/16.0;'
  double rand1;
  double rand2;
  for (int i=0;i<N;i++){
    for(int j=0;j<2;j++){
      x[j]=-0.25*log(1.0-RnG(gen));
    }
    for(int j=2;j<4;j++){
      x[j]=PI*RnG(gen); // makes x[j] a random number between 0 and PI
    }
    for(int j=4;j<6;j++){
      x[j]=2*PI*RnG(gen); // makes x[j] a random number between 0 and 2PI
    }
    fx=int_func_montecarlo(x[0],x[1],x[2],x[3],x[4],x[5]);
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
