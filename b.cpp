#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;
double* createArray(int n); //Creates an Array of length n, filled with zeros
void fillArray(double* arr,int n,double num); //fills an array (of length n) (by reference) with a given number (num)
void fillArrayFunction(double * arr, int n,double *x,double (*func)(double)); //same as fillArray, but doesnt take an umber, but a generic function
double func(double x);
double* solve(double *a, double *b, double*c, double *deriv, int n); //The actual solving algorithm, returns an array of length n+2 with solutions
int main(int argc, char** argv){
  if (argc<2){
    cout << "You need to state a number n" << endl;
    exit(1);
  }
  int n=atoi(argv[1]);
  ofstream outfile;
  outfile.open("oppb.txt");
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
  for (int i=0; i<n+2;i++){
    cout << x[i]<<"  "<< a[i] << " "<<b[i] << " "<<c[i] << " "<<deriv[i] << " "<< endl;
  }
  double *solution=solve(a,b,c,deriv,n);
  outfile << "n: x: accurate_solution: approximate_solution: "<<"\n";
  for(int i=0;i<=n+1;i++){
    outfile <<n<<setprecision(16)<<" "<<x[i]<<" "<< 1-(1-exp(-10))*x[i]-exp(-10*x[i])<<" "<<solution[i]<<"\n";
  }
}
double * createArray(int n){
  double *arr;
  arr=new double[n];
  for(int i=0;i<n;i++){
    arr[i]=0.0;
  }
  return arr;
}
void fillArray(double * arr, int n, double num){
  for (int i=0;i<n;i++){
    arr[i]=num;
  }
}
void fillArrayFunction(double * arr, int n,double *x,double (*func)(double)){
  for (int i=0;i<n;i++){
    arr[i]=(*func)(x[i]);
  }
}
double func(double x){
  return 100.0*exp(-10.0*x);
}
double* solve(double *a,double *b, double *c, double *deriv, int n){
  double* solutions;
  solutions=createArray(n+2); //First and last value are already set to zero
  for(int i=1;i<n;i++){
    b[i+1]=b[i+1]-a[i+1]*c[i]/b[i];
    deriv[i+1]=deriv[i+1]-a[i+1]/b[i]*deriv[i];
  }
  solutions[n]=deriv[n]/b[n];
  delete [] a; // A not needed anymore (set to 0 anyways)
  for(int i=n-1;i>=1;i--){
    solutions[i]=(deriv[i]-c[i]*solutions[i+1])/b[i];
  }
  delete [] c; delete [] b; //c and b not needed anymore
  return solutions;
}
