#include "vecop.hpp"
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>
using namespace std;
int main(int argc, char** argv){
  int n;
  ofstream outfile;
  outfile.open("solutions_harmonic_oscillator.txt");
  const double PI = atan(1.0)*4;
  if(argc>=2){
    n=atoi(argv[1]);;
  }
  else{
    cout << "You need to state a number n" << endl;
    exit(1);
  }
  double **A=createNNMatrix(n);
  double rhomax=1;double rhomin=0;
  double h=(rhomax-rhomin)/n;
  double hh=h*h;
  double d=2.0/hh;
  double a=-1/hh;
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
    cout << solutions[i] << " " << accurate_sol[i] << endl;
  }
  outfile.close();
}
