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
/*PROJECT TWO*/
double ** matrixMult(int n, double **A, double **B){
  double **returnMatrix=createNNMatrix(n);
  for(int i=0 ; i < n ; i++) {
    for(int j=0 ; j < n ; j++) {
      for(int k=0 ; k < n ; k++) {
        returnMatrix[i][j]+=A[i][k]*B[k][j];
      }
    }
  }
  return returnMatrix;
}

double ** createJacobiMatrix(int n, int i, int j, double s,double c){
  double ** S= createNNMatrix(n);
  for (int k=0; k<n;k++){
    S[k][k]=1.0;
  }
  S[i][i]=c; S[j][j]=c; S[j][i]=-s; S[i][j]=s;
  return S;
}
void invertJacobiMatrix(double **S, int n, int i, int j){
  S[i][j]=-S[i][j];
  S[j][i]=-S[j][i];
}
double** similarity_transform(double**A,int n, int k, int l, double s, double c){
  double **S=createJacobiMatrix(n,k,l,s,c);
  double **SA=matrixMult(n,S,A);
  invertJacobiMatrix(S,n,k,l);
  A=matrixMult(n,SA,S);
  deleteNNMatrix(S,n); deleteNNMatrix(SA,n);
  return A;
}
void jacobi_diag(double** A, int n, double tol){
  double tau, tant, sint, cost;
  double rootElem;
  double* maxelem=findMax(A,n);
  int k= round(maxelem[1]); int l=round(maxelem[2]);
  while(maxelem[0]>tol){
    tau=(A[l][l]-A[k][k])/(2.0*A[k][l]);
    rootElem=sqrt(1+tau*tau);
    if (fabs(-tau+rootElem)<fabs(-tau-rootElem)){
      tant=-tau+rootElem;
    }
    else{
      tant=-tau-rootElem;
    }
    cost=1/sqrt(1+tant*tant);
    sint=tant*cost;
    A=similarity_transform(A,n,k,l,sint,cost);
    maxelem=findMax(A,n); k= round(maxelem[1]); l=round(maxelem[2]);
    printMatrix(A,n);
  }
}
void jacobi_diagNy(double **A, int n, double tol){
  double tau,tant, sint, cost;
  double* maxelem=findMax(A,n);
  int k=round(maxelem[1]); int l=round(maxelem[2]); double maxSquare=maxelem[0];
  while(maxSquare>tol){
    tau=(A[l][l]-A[k][k])/(2.0*A[k][l]);
    if ( tau > 0 ) {
      tant = 1.0/(tau + sqrt(1.0 + tau*tau));
    }
    else {
      tant = -1.0/( -tau + sqrt(1.0 + tau*tau));
    }
    cost=1/sqrt(1+tant*tant);
    sint=cost*tant;
    double a_kk,a_ll,a_ik,a_il;
    a_kk=A[k][k];
    a_ll=A[l][l];
    A[k][k]=cost*cost*a_kk-2.0*cost*sint*A[k][l]+sint*sint*a_ll;
    A[l][l]=sint*sint*a_kk+2.0*cost*sint*A[k][l]+cost*cost*a_ll;
    A[k][l]=0.0;
    A[l][k]=0.0;
    for(int i=0; i<n;i++){
      if(i != k && i != l){
        a_ik=A[i][k];
        a_il=A[i][l];
        A[i][k]=cost*a_ik-sint*a_il;
        A[k][i]=A[i][k];
        A[i][l]=cost*a_il+sint*a_ik;
        A[l][i]=A[i][l];
      }
    }
    maxelem=findMax(A,n);
    k=round(maxelem[1]);
    l=round(maxelem[2]);
    maxSquare=maxelem[0];
  }
}
double* findMax(double **A, int n){ // returns an array where arr[0] is the largest value, arr[1] the i element and arr[2] the j element
  double maxelem=0;
  double asquare;
  double* returnarray=createArray(3);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(j==i){
        continue;
      }
      asquare=A[i][j]*A[i][j];
      if(asquare>maxelem){
        maxelem=asquare;
        returnarray[0]=asquare; returnarray[1]=i; returnarray[2]=j;
      }
    }
  }
  return returnarray;
}
double ** copyMatrix(double**A,int n){
  double** B=createNNMatrix(n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      B[i][j]=A[i][j];
    }
  }
  return B;
}
void printMatrix(double **A, int n){
  std::cout << std::setprecision(5) << std::fixed;
  for (int i=0; i<n;i++){
    for(int j=0;j<n;j++){
      cout << A[i][j] << " ";
    }
    cout << endl;
  }
}
