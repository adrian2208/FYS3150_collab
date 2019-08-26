#include "vecop.hpp"
string createFileName(string name,int n){
  return name.append(to_string(n)).append(".txt");
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
/* Oppgave c start*/
double* improvedSolve(double *deriv, int n){
  double *solutions; double *b;
  solutions=createArray(n+2);
  b=createArray(n+1); //First element is zero, but we don't care about first element
  for(int i=1;i<=n;i++){
    b[i]=((double)(i+1))/((double) i);
  }
  for(int i=2;i<=n;i++){
    deriv[i]=deriv[i]+deriv[i-1]/b[i-1];
  }
  solutions[n]=deriv[n]/b[n];
  for(int i=n-1;i>=1;i--){
    solutions[i]=(deriv[i]+solutions[i+1])/b[i];
  }
  delete [] b;
  return solutions;
}
/* Oppgave c end*/
