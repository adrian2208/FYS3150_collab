#include "lib.h"
#include "functions.hpp"
#include <math.h>
#include <iomanip>
#include <iostream>

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
  double *x = new double[N];
  double *w=new double[N];
  double lambda=3;
  double int_gauss=0.;
  double add_var;
  gauleg(-lambda,lambda,x,w,N); // Calls the lib.cpp function gauleg // THIS  LOOKS LIKE GULAG
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
  cout <<" approach :"<< int_gauss << " correct: "<< 5*3.14159265359*3.14159265359/(16*16)<< endl;

}
