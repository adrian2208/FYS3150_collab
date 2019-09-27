#include "vecop.hpp"
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include "time.h"
using namespace std;
int main(int argc, char** argv){
  clock_t start, finish; //Start and finish time
  int n;
  ofstream outfile;
  outfile.open("solutions_harmonic_oscillator.txt"); //Output file
  ofstream outfile_time;
  outfile_time.open("time_info.txt",ios::out | ios::app); //time-info file, shared with armadillo-one
  const double PI = atan(1.0)*4; //Because pi's not in there. Genious
  if(argc>=2){
    n=atoi(argv[1]); //n needs to be stated
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
  double accurate_sol [n]; //The correct values
  for(int i=1;i<=n;i++){
    accurate_sol[i-1]=(d+2*a*cos(i*PI/(n+1.0)));
  }
  sort(accurate_sol,accurate_sol+n); //Sorting for exposal
  double **R=createNNMatrix(n);
  start=clock();
  int amount=jacobi_diag2(A,R,n,1e-8); //
  finish=clock();
  double ellapsed_time=((finish-start)/(float)CLOCKS_PER_SEC);
  outfile_time <<"egen n: "<<n<< " amount_looping: "<<amount<<" ellapsed time: "<<ellapsed_time<<endl;
  outfile_time.close();
  for(int i=0;i<n;i++){
    solutions[i]=A[i][i]; //Gets the solution
  }
  deleteNNMatrix(A,n); //Done with A
  outfile << "n: " <<n<< endl;
  outfile << "rhomax: "<<rhomax<<endl;
  outfile << "rhomin: "<<rhomin<<endl;
  outfile << endl;
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
  deleteNNMatrix(R,n); //Done with R
  sort(solutions,solutions+n);
  for(int i=0; i<n;i++){
    cout << solutions[i]- accurate_sol[i]<< endl;
  }
  outfile.close();
}
