#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <ctime>
#include "vecop.hpp"
using namespace std;
int main(int argc, char** argv){
  clock_t start, finish;
  start=clock();
  bool writeout=false;
  if (argc<2){
    cout << "You need to state a number n" << endl;
    cout << "If you want to create a data file, add a second term"<<endl;
    exit(1);
  }
  if (argc>=3){
    writeout=true;
  }
  int n=atoi(argv[1]);
  ofstream outfile;
  if(writeout){
    outfile.open(createFileName("oppc_",n));
  }
  ofstream outfile_time;
  outfile_time.open("time_info.txt",ios::out | ios::app);
  double h= 1.0/(n+1);
  double hh=h*h;
  double *x=createArray(n+2); //Length is given as n+2 because we have the two extra values for n=0 (being 0) and n=1 (being 0)
  double *deriv=createArray(n+2); // The second derivative, values from 0 to n+1
  for(int i=0;i<=n+1;i++){
    x[i]=i*h;
  }
  fillArrayFunction(deriv,n+2,x,func);
  for (int i=0;i<n+2;i++){
    deriv[i]*=hh;
  }

  double *solution=improvedSolve(deriv,n); //Solves using the improved algorithm
  finish=clock();
  double ellapsed_time=((finish-start)/(float)CLOCKS_PER_SEC);
  outfile_time <<"c n: "<<n<<" ellapsed time: "<<ellapsed_time<<endl;
  outfile_time.close();
  if(writeout){
    outfile << "n: "<<n<<" ellapsed time: "<<ellapsed_time<<endl;
    outfile << "x: accurate_solution: approximate_solution relative_error:"<<endl;
    double accurate_sol;
    for(int i=0;i<=n+1;i++){
      accurate_sol=1-(1-exp(-10))*x[i]-exp(-10*x[i]);
      outfile <<setprecision(16)<<x[i]<<" "<<accurate_sol<<" "<<solution[i]<<" "<<fabs((accurate_sol-solution[i])/accurate_sol)<<"\n";
    }
    delete [] solution;
    delete [] x;
    outfile.close();
  }
}
