#include "catch.hpp"
#include "vecop.hpp"
#include <algorithm>
#include <cmath>
using namespace std;
double** makeTestMatrix(){ //Makes a predefined test matrix where we know eigenvalues and vectors.
  double **A=createNNMatrix(5);
  double elem[25]={2.0,0.0,1.0,2.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0,1.0,4.0,3.0,3.0,2.0,1.0,3.0,4.0,3.0,1.0,1.0,3.0,3.0,3.0};
  for(int i=0;i<5;i++){
    for(int j=0;j<5;j++){
      A[i][j]=elem[5*i+j];
    }
  }
  return A;
}
void testBiggestElement(){ //Tests if biggest element is recognized. I am well aware that I destroy the symmetry of the matrix.
  double **testMatrix=makeTestMatrix();
  testMatrix[1][3]=34.0;
  double *vecci=findMax(testMatrix,5);double maxelem=vecci[0];
  if(abs(maxelem-34.0*34.0)>1e-10){
    throw invalid_argument("Biggest Element finder is broken");
  }
  deleteNNMatrix(testMatrix,5);
}
TEST_CASE("Check if biggest element finder works"){
  double **testMatrix=makeTestMatrix();
  testMatrix[1][3]=34.0;
  double *vecci=findMax(testMatrix,5);double maxelem=vecci[0];
  REQUIRE((maxelem-34.0*34.0)<1e-10);
}
TEST_CASE("Test if correct eigenvalues are found, and wether eigenvectors are orthogonal"){
  double **testMatrix=makeTestMatrix();
  double **solvMatrix=createNNMatrix(5);
  double correigenvalues [5]={0.1545,0.3691,0.8115,2.0324,10.6325}; // Solution from using MatLab
  double gottenEigenvalues[5];
  jacobi_diag(testMatrix,solvMatrix,5,1e-8);
  for(int i=0;i<5;i++){
    gottenEigenvalues[i]=testMatrix[i][i]; //The elements on the diagonal are the eigenvalues
  }
  sort(gottenEigenvalues,gottenEigenvalues+5);
  for(int i=0;i<5;i++){
    REQUIRE((gottenEigenvalues[i]-correigenvalues[i])<1e-4); //The  correct eigenvalues and the eigenvalues I got need to be approximately identical
  }
  float sum=0.0;
  for(int i=0;i<5;i++){
    for(int j=i+1;j<5;j++){
      for (int k=0;k<5;k++){
        //cout <<solvMatrix[i][k] << " " <<  solvMatrix[j][k];
        sum+=solvMatrix[i][k]*solvMatrix[j][k];

      }
      REQUIRE(sum<1e-5); //The vector product needs to be very little, approximating zero
      sum=0;
    }
  }
}
