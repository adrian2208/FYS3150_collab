#include "vecop.hpp"
using namespace std;
string createFileName(string name,int n){
  return name.append(to_string(n)).append(".txt");
}
double** createNNMatrix(int n){
  double** A;
  A = new double*[n];
  for (int i = 0; i < n; i++){
    A[i] = new double[n];
  }
  return A;
}
void deleteNNMatrix(double** A,int n){
  for (int i = 0; i < n; i++){
    delete[] A[i];
  }
delete[] A;
}
// Creates an array of length n with all values set to zero
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
void solve(double *a,double *b, double *c, double *deriv, int n, double *solutions){
  double precalc;
  for(int i=1;i<n;i++){
    precalc=a[i+1]/b[i];
    b[i+1]=b[i+1]-precalc*c[i]; //b[i+1]=b[i+1]-a[i+1]*c[i]/b[i];
    deriv[i+1]=deriv[i+1]-precalc*deriv[i];
  }
  solutions[n]=deriv[n]/b[n];
  delete [] a; // A not needed anymore (set to 0 anyways)
  for(int i=n-1;i>=1;i--){
    solutions[i]=(deriv[i]-c[i]*solutions[i+1])/b[i];
  }
  delete [] c; delete [] b; delete [] deriv; //c and b and deriv are not needed anymore
}
/* Oppgave c start*/
void improvedSolve(double *deriv,double *b, int n, double *solutions){
  for(int i=2;i<=n;i++){
    deriv[i]=deriv[i]+deriv[i-1]/b[i-1];
  }
  solutions[n]=deriv[n]/b[n];
  for(int i=n-1;i>=1;i--){
    solutions[i]=(deriv[i]+solutions[i+1])/b[i];
  }
  delete [] b; delete [] deriv;
}
/* Oppgave c end*/
