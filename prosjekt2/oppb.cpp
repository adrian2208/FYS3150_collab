#include <iostream>
#include "vecop.hpp"
using namespace std;
void makeTestMatrix(double **A,int n){
  double elem[25]={2.0,0.0,1.0,2.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,4.0,3.0,3.0,2.0,1.0,3.0,4.0,3.0,1.0,1.0,3.0,3.0,3.0};
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      A[i][j]=elem[n*i+j];
    }
  }
}
void testJacobi(double **A, int n){
  jacobi_diagNy(A,n,1e-8);
  cout <<"no"<<endl;
  for(int i=0;i<n;i++){
    cout << A[i][i]<<endl;
  }
}
int main(){
  double **A=createNNMatrix(5);
  makeTestMatrix(A,5);
  printMatrix(A,5);
  double * b=findMax(A,5);
  cout << b[0]<<" " << b[1] <<" "<< b[2]<<endl;
  testJacobi(A,5);
}
