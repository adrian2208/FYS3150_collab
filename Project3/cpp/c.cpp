#include "lib.h"
#include "functions.hpp"
#include <random>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <time.h>
#include <fstream>
using namespace std;

int main(int argc, char** argv){
  ofstream outfile;
  outfile.open("../results/time_info.csv",ios::out | ios::app); //time-info file
  clock_t start, finish;
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
  std::uniform_real_distribution<double> RnG(0.0,1.0); //Random number generator
  double lambda=2; //length of the box
  double int_mc=0.; double variance=0.;
  double sum_sigma=0.;
  double x[6]; // In order to store 6 random numbers
  double fx;
  double jacobi_det=pow(2*lambda,6); //The factor in front, calculation can be done in the end
  start=clock();
  for (int i=0;i<N;i++){
    for(int j=0;j<6;j++){
      x[j]=-lambda+2*lambda*RnG(gen); // makes x[j] a random number between -lambda and lamnbda
    }
    fx=int_func(x[0],x[1],x[2],x[3],x[4],x[5]);
    int_mc+=fx;
    sum_sigma+=fx*fx;
  }
  int_mc=int_mc/N;//*jacobi_det to get real result
  sum_sigma=sum_sigma/N; //*jacobi-det² to get real result
  variance=sum_sigma-int_mc*int_mc; // This is E(X^2)-E(X)^2 *jacobi_det^2 to get real result
  finish=clock();
  double ellapsed_time=((finish-start)/(float)CLOCKS_PER_SEC);
  double correct_result=5*3.14159265359*3.14159265359/(16*16);
  double relative_error=fabs(int_mc*jacobi_det-correct_result)/correct_result;
  outfile<<"\nbadMonteCarlo,no,1,"<<N<<","<<ellapsed_time<<","<<int_mc*jacobi_det<<","<<relative_error<<","<<jacobi_det*sqrt(variance/((double) N ));
  outfile.close();
}
