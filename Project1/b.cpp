#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "vecop.hpp"
#include <ctime>
using namespace std;
int main(int argc, char** argv){
  clock_t start, finish;
  bool writeout=false; //If you want to write out a file containing the calculated solutions
  /*If argc=2, only an append is made to the time_info.txt file, if arc>=3, writeout is set to be true*/
  if (argc<2){
    cout << "You need to state a number n" << endl;
    cout << "If you want to make a file, add a second term"<<endl;
    exit(1);
  }
  if (argc>=3){
    writeout=true;
  }
  int n=atoi(argv[1]);
  ofstream outfile;
  if(writeout){
    outfile.open(createFileName("oppb_",n));
  }
  ofstream outfile_time;
  outfile_time.open("time_info.txt",ios::out | ios::app);
  double h= 1.0/(n+1);
  double hh=h*h; //efficiency
  double *x=createArray(n+2); //x values from 0 to 1 with steplength h
  double *a=createArray(n+1); //The 0th value and 1st value are not used, but I add it for my own sake of overview
  double *b=createArray(n+1); // from 0 to n, but the 0th value is ommitted
  double *c=createArray(n);
  double *deriv=createArray(n+2); // The second derivative, values from 0 to n+1
  fillArray(a,n+1,-1);
  fillArray(b,n+1,2);
  fillArray(c,n,-1);
  for(int i=0;i<=n+1;i++){
    x[i]=i*h;
  }
  fillArrayFunction(deriv,n+2,x,func); //fills the deriv array with the function defined in vecop.cpp
  for (int i=0;i<n+2;i++){
    deriv[i]*=hh;
  }

  double* solution=createArray(n+2);
  start=clock();
  solve(a,b,c,deriv,n,solution); //Solves, and does changes to the "solution"-array
  finish=clock();
  double ellapsed_time=((finish-start)/(float)CLOCKS_PER_SEC);
  outfile_time <<"b n: "<<n<<" ellapsed time: "<<ellapsed_time<<endl;
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
