#include "vecop.hpp"
#include <cmath>
#include <algorithm>
using namespace std;
int main(){
  const double PI = atan(1.0)*4;
  int n=100;
  double **A=createNNMatrix(n);
  double rhomax=1;double rhomin=0;
  double h=(rhomax-rhomin)/n;
  double hh=h*h;
  double d=hh*2.0;
  double a=-1*hh;
  for(int i=0;i<n;i++){ //Fill diagonal with 2*hh
    A[i][i]=d;
  }
  for(int i=0;i<n-1;i++){ //Fill the other two rows with -1*hh
    A[i+1][i]=a; //a
    A[i][i+1]=a; // a
  }
  double solutions [n];
  double accurate_sol [n];
  for(int i=1;i<=n;i++){
    accurate_sol[i-1]=(d+2*a*cos(i*PI/(n+1.0)));
  }
  sort(accurate_sol,accurate_sol+n);
  jacobi_diag(A,n,10e-20);
  for(int i=0;i<n;i++){
    solutions[i]=A[i][i];
  }
  sort(solutions,solutions+n);
  for(int i=0; i<n;i++){
    cout << solutions[i] << " " << accurate_sol[i] << endl;
  }
}
