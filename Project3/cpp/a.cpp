#include "lib.h"
#include "functions.hpp"
#include <math.h>
#include <iomanip>
#include <iostream>
#include <time.h>
#include <fstream>
using namespace std;
int main(int argc, char** argv){
  ofstream outfile;
  outfile.open("../results/time_info.csv",ios::out | ios::app); //time-info file
  //outfile<<"Type,Parrallelisation_flag,Amount_threads,N,Time_use,Result,relative_error,standard_deviation";
  int N;
  clock_t start, finish;
  if(argc>=2){
    N=atoi(argv[1]); //n needs to be stated
  }
  else{
    cout << "You need to state a number N" << endl;
    exit(1);
  }
  double *x = new double[N];
  double *w=new double[N];
  double lambda=2;
  double int_gauss=0.;
  double add_var;
  gauleg(-lambda,lambda,x,w,N); // Calls the lib.cpp function gauleg // THIS  LOOKS LIKE GULAG
  start=clock();
  for (int i=0;i<N;i++){
    for (int j = 0;j<N;j++){
      for (int k = 0;k<N;k++){
        for (int l = 0;l<N;l++){
           for (int m = 0;m<N;m++){
              for (int n = 0;n<N;n++){
                add_var=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_func(x[i],x[j],x[k],x[l],x[m],x[n]);
                /*This is clearly only allowed because all 6 dimensions have the same dimensionality */
                //cout <<"val of addvar " <<add_var<<endl;
                int_gauss+=add_var;
     		      }
            }
        }
      }
    }
  }
  finish=clock();
  delete [] x; delete [] w;
  double ellapsed_time=((finish-start)/(float)CLOCKS_PER_SEC);
  double correct_result=5*3.14159265359*3.14159265359/(16*16);
  double relative_error=fabs(int_gauss-correct_result)/correct_result;
  outfile<<"\ngauleg,no,1,"<<N<<","<<ellapsed_time<<","<<int_gauss<<","<<relative_error<<","<<0;
  outfile.close();
}
