//#include "catch.hpp"
#include "vecop.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
using namespace std;
//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

double** makeTestMatrix(){
  double **A=createNNMatrix(5);
  double elem[25]={2.0,0.0,1.0,2.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,4.0,3.0,3.0,2.0,1.0,3.0,4.0,3.0,1.0,1.0,3.0,3.0,3.0};
  for(int i=0;i<5;i++){
    for(int j=0;j<5;j++){
      A[i][j]=elem[5*i+j];
    }
  }
  return A;
}
void testBiggestElement(){
  double **testMatrix=makeTestMatrix();
  testMatrix[1][3]=34.0;
  double *vecci=findMax(testMatrix,5);double maxelem=vecci[0];
  if(abs(maxelem-34.0*34.0)>1e-10){
    throw invalid_argument("Biggest Element finder is broken");
  }
  deleteNNMatrix(testMatrix,5);
}
void testEigenvalue(){
  double **testMatrix=makeTestMatrix();
  double correigenvalues [5]={0.1545,0.3691,0.8115,2.0324,10.6325}; // Solution from using MatLab
  double gottenEigenvalues[5];
  jacobi_diag(testMatrix,5,1e-8);
  for(int i=0;i<5;i++){
    gottenEigenvalues[i]=testMatrix[i][i];
  }
  sort(gottenEigenvalues,gottenEigenvalues+5);
  for(int i=0;i<5;i++){
    if(abs(gottenEigenvalues[i]-correigenvalues[i])>1e-4){
      cout<<"got: "<<gottenEigenvalues[i]<<" expected:"<< correigenvalues[i]<<endl;
      throw invalid_argument("Wrong eigenvalue");
      //throw invalid_argument(("got: %.4f expected: %.4f",correigenvalues[i],gottenEigenvalues[i]));
    }
  }
  deleteNNMatrix(testMatrix,5);
}
int main(){
  testBiggestElement();
  testEigenvalue();
}
