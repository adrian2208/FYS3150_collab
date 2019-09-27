#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;
int main(int argc, char** argv){
  clock_t start, finish; //Start and end time
  int n;
  ofstream outfile_time;
  outfile_time.open("time_info.txt",ios::out | ios::app); //Appending to time-info
  const double PI = atan(1.0)*4;
  if(argc>=2){
    n=atoi(argv[1]);;
  }
  else{
    cout << "You need to state a number n" << endl;
    exit(1);
  }
  mat A(n,n); //initialize matrix
  double rhomax=1;double rhomin=0;
  double h=(rhomax-rhomin)/n;
  double hh=h*h;
  double d=2.0/hh;
  double a=-1/hh;
  for(int i=0;i<n;i++){ //Fill diagonal with 2*hh
    A(i,i)=d;
  }
  for(int i=0;i<n-1;i++){ //Fill the other two rows with -1*hh
    A(i+1,i)=a; //a
    A(i,i+1)=a; // a
  }
  vec eigenval;
  mat eigenvec;
  start=clock();
  eig_sym(eigenval,eigenvec,A);
  finish=clock();
  double ellapsed_time=((finish-start)/(float)CLOCKS_PER_SEC);
  for(int i=0;i<n;i++){
    cout << eigenval(i) << endl;
  }
  outfile_time <<"arma n: "<<n<<" ellapsed time: "<<ellapsed_time<<endl;

}
