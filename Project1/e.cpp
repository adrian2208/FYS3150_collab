#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "vecop.hpp"
#include "time.h"
#include "lib.h"
#include "vecop.hpp"
using namespace std;
int main(int argc, char** argv){
  ofstream outfile;
  clock_t start, finish;
  start=clock();
  bool writeout=false;
  if (argc<2){
    cout << "You need to state a number n" << endl;
    cout << "If you want to make a file, add a second term"<<endl;
    exit(1);
  }
  if (argc>=3){
    writeout=true;
  }
  int n=atoi(argv[1]);
  if(writeout){
    outfile.open(createFileName("oppe_",n));
  }
  ofstream outfile_time;
  outfile_time.open("time_info.txt",ios::out | ios::app);
  double **A=createNNMatrix(n); //creates an empty n*n matrix
  double h= 1.0/(n+1);
  double hh=h*h;
  double *x=createArray(n);
  double *deriv=createArray(n); // The second derivative, values from 0 to n+1
  for(int i=1;i<=n;i++){
    x[i-1]=i*h;
  }

  fillArrayFunction(deriv,n,x,func);
  for (int i=0;i<n;i++){
    deriv[i]*=hh;
  }
  for(int i=0;i<n;i++){ //Fill diagonal with 2
    A[i][i]=2.0;
  }
  for(int i=0;i<n-1;i++){ //Fill the other two rows with -1
    A[i+1][i]=-1.0; //a and c
    A[i][i+1]=-1.0; // a and c
  }
  int          i,j, *indx;
  double       d, *col, **y;
  // allocate space in memory
  indx = new int[n];
  col  = new double[n];
  ludcmp(A,n,indx,&d);
  lubksb(A,n,indx,deriv);
  deleteNNMatrix(A,n);
  finish=clock();
  double ellapsed_time=((finish-start)/(float)CLOCKS_PER_SEC);
  outfile_time <<"e n: "<<n<<" ellapsed time: "<<ellapsed_time<<endl;
  outfile_time.close();
  if(writeout){
    outfile << "n: "<<n<<" ellapsed time: "<<ellapsed_time<<endl;
    outfile << "x: accurate_solution: approximate_solution relative_error:"<<endl;
    double accurate_sol;
    for(int i=0;i<=n+1;i++){
      accurate_sol=1-(1-exp(-10))*x[i]-exp(-10*x[i]);
      outfile <<setprecision(16)<<x[i]<<" "<<accurate_sol<<" "<<deriv[i]<<" "<<fabs((accurate_sol-deriv[i])/accurate_sol)<<"\n";
    }
    delete [] deriv;
    delete [] x;
    outfile.close();
  }
}
