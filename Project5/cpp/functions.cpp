#include "functions.h"
double** createNMatrix(int n,int m){ //Create n*m integer matrix
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
double distance(double **pos){
  return distance(pos[0][0],pos[0][1],pos[0][2],pos[1][0],pos[1][1],pos[1][2]);
}
double distance(double x1, double y1, double z1, double x2, double y2, double z2){
  return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}
