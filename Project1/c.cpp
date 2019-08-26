#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "vecop.hpp"
using namespace std;
int main(int argc, char** argv){
  if (argc<2){
    cout << "You need to state a number n" << endl;
    exit(1);
  }
  int n=atoi(argv[1]);
  ofstream outfile;
  outfile.open(createFileName("oppc_",n));
  double h= 1.0/(n+1);
  double hh=h*h;
  double *x=createArray(n+2);
  double *deriv=createArray(n+2); // The second derivative, values from 0 to n+1
  for(int i=0;i<=n+1;i++){
    x[i]=i*h;
  }
  fillArrayFunction(deriv,n+2,x,func);
  for (int i=0;i<n+2;i++){
    deriv[i]=deriv[i]*=hh;
  }

  double *solution=improvedSolve(deriv,n);
  outfile << "n: x: accurate_solution: approximate_solution: "<<"\n";
  for(int i=0;i<=n+1;i++){
    outfile <<n<<setprecision(16)<<" "<<x[i]<<" "<< 1-(1-exp(-10))*x[i]-exp(-10*x[i])<<" "<<solution[i]<<"\n";
  }
  outfile.close();
}
