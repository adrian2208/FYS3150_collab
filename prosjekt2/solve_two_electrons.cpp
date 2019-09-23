#include "vecop.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
using namespace std;
double potential(double r,double omega){
  return r*r*omega*omega+1/r;
}
int main(int argc, char** argv){
  int n;
  double omega;
  ofstream outfile;
  outfile.open("solutions_two_electrons.txt");
  const double PI = atan(1.0)*4;
  if(argc==2){
    n=atoi(argv[1]);
    omega=1.0;
  }
  else if(argc==3){
    n=atoi(argv[1]);
    omega=atof(argv[2]);
  }
  else{
    cout << "You need to state a number n and omega" << endl;
    exit(1);
  }
  double **A=createNNMatrix(n);
  double rhomax=11;double rhomin=0;
  double *rho=createArray(n);

  double h=(rhomax-rhomin)/(n+1);
  double hh=h*h;
  double d=2.0/hh;
  double a=-1/hh;
  for(int i=0;i<n-1;i++){ //Fill the other two rows with -1*hh
    A[i+1][i]=a; //a
    A[i][i+1]=a; // a
  }
  for(int i=0;i<n;i++){
    rho[i]=(i+1)*h;
  }
  for(int i=0;i<n;i++){ //Fill diagonal with 2*hh
    //A[i][i]=d+rho[i]*rho[i];
    A[i][i]=d+potential(rho[i],omega);
  }
  double solutions [n];
  double **R=createNNMatrix(n);
  jacobi_diag(A,R,n,10e-8);
  for(int i=0;i<n;i++){
    solutions[i]=A[i][i];
  }
  deleteNNMatrix(A,n);

  outfile << "n: " <<n<< endl;
  outfile << "rhomax: "<<rhomax<<endl;
  outfile << "rhomin: "<<rhomin<<endl;
  for(int i=0; i<n;i++){
    outfile <<solutions[i] <<" ";
  }
  outfile << endl;
  for(int i=0; i<n;i++){
    for (int j=0;j<n;j++){
      outfile << R[i][j]<< " ";
    }
    outfile << endl;
  }
  deleteNNMatrix(R,n);
  sort(solutions,solutions+n);
  for(int i=0; i<n;i++){
    cout << solutions[i] << endl;
  }
  outfile.close();
}
