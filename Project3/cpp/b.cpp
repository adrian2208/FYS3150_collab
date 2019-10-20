#include "lib.h"
#include "functions.hpp"
#include <time.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#define PI 3.14159265359
using namespace std;
void moveToLeft(double *arr,int N){
  /*Takes an array of length N+1 and moves each element one "to the left". Thus, the 1. element is overwritten, and the Last
  one occurs twice.*/
  for (int i=0;i<N;i++){
    arr[i]=arr[i+1];
  }
}
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
  double *theta = new double[N];
  double *phi = new double[N];
  double *r=new double[N+1];
  double *omega_theta=new double[N];
  double *omega_phi=new double[N];
  double *omega_r=new double[N+1];
  double int_gauss=0.0, add_var;
  gauleg(0,PI,theta,omega_theta,N);
  gauleg(0,2.0*PI,phi,omega_phi,N);
  gauss_laguerre(r,omega_r,N,0);
  moveToLeft(r,N);
  moveToLeft(omega_r,N);
  start=clock();
  for (int i=0;i<N;i++){
    for (int j = 0;j<N;j++){
      for (int k = 0;k<N;k++){
        for (int l = 0;l<N;l++){
           for (int m = 0;m<N;m++){
              for (int n = 0;n<N;n++){
                add_var=omega_r[i]*omega_r[j]*omega_theta[k]*omega_theta[l]*omega_phi[m]*omega_phi[n]*
                int_func_polar_gaulag(r[i],r[j],theta[k],theta[l],phi[m],phi[n]); //evaluation of function times weights
                int_gauss+=add_var;
     		      }
            }
        }
      }
    }
  }
  delete [] r; delete [] omega_r;delete [] phi; delete [] omega_phi;delete [] theta; delete [] omega_theta;
  finish=clock();
  double ellapsed_time=((finish-start)/(float)CLOCKS_PER_SEC);
  double correct_result=5*3.14159265359*3.14159265359/(16*16);
  double relative_error=fabs(int_gauss-correct_result)/correct_result;
  outfile<<"\ngaulag,no,1,"<<N<<","<<ellapsed_time<<","<<int_gauss<<","<<relative_error<<","<<0;
  outfile.close();
}
