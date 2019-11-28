#include "functions.h"
double** createNMatrix(int n,int m){
  double** A;
  A = new double*[n];
  for (int i = 0; i < n; i++){
    A[i] = new double[m];
  }
  for (int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      A[i][j]=i*m;
    }
  }
  return A;
}
