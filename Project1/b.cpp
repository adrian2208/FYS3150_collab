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
  outfile.open(createFileName("oppb_",n));
  double h= 1.0/(n+1);
  double hh=h*h;
  double *x=createArray(n+2);
  double *a=createArray(n); //The 0th value is not used, but I add it for my own sake of overview
  double *b=createArray(n+1);
  double *c=createArray(n);
  double *deriv=createArray(n+2); // The second derivative, values from 0 to n+1
  fillArray(a,n,-1);
  fillArray(b,n+1,2);
  fillArray(c,n,-1);
  for(int i=0;i<=n+1;i++){
    x[i]=i*h;
  }
  fillArrayFunction(deriv,n+2,x,func);
  for (int i=0;i<n+2;i++){
    deriv[i]=deriv[i]*=hh;
  }
  /*
  for (int i=0; i<n+2;i++){
    cout << x[i]<<"  "<< a[i] << " "<<b[i] << " "<<c[i] << " "<<deriv[i] << " "<< endl;
  }*/
  double *solution=solve(a,b,c,deriv,n);
  outfile << "n: x: accurate_solution: approximate_solution: "<<"\n";
  for(int i=0;i<=n+1;i++){
    outfile <<n<<setprecision(16)<<" "<<x[i]<<" "<< 1-(1-exp(-10))*x[i]-exp(-10*x[i])<<" "<<solution[i]<<"\n";
  }
  outfile.close();
}
