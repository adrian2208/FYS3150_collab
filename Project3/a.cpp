#include "lib.h"

#include <math.h>
#include <iomanip>
#include <iostream>

using namespace std;
int g;
double int_func(double x1, double x2, double y1, double y2, double z1, double z2){
  double exponential_part=exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2))); // Exponential part of the function
  double denom=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)); // Denominator part
  double returnval;
  //cout << denom << endl;
  if (fabs(denom)<= 1e-15){ //If denominator is zero, it is omitted
    returnval=0;
    g++;
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
  double *x = new double[N];
  double *w=new double[N];
  double *x2=new double[N];
  double *w2=new double[N];
  double lambda=3;
  double int_gauss=0.;
  double int_gauss2=0.;
  double add_var;
  double add_var2=0;
  gauleg(-lambda,lambda,x,w,N); // Calls the lib.cpp function gauleg // THIS  LOOKS LIKE GULAG
  gauleg(-lambda-1,lambda-1,x2,w2,N);
  for(int i=0; i<N;i++){
    cout<<w[i]<< " "<<x[i]<<endl;
  }
  cout <<endl;
  for(int i=0; i<N;i++){
    cout<<w2[i]<< " "<<x2[i]<<endl;
  }
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
  cout <<" approach 1:"<< int_gauss << " approach 2: "<< int_gauss2<< " correct: "<< 5*3.14159265359*3.14159265359/(16*16)<< endl;
  cout << g << endl;
}
