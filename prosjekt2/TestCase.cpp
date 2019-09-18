#include "catch.hpp"
#include "vecop.hpp"
#include <algorithm>
#include <cmath>
using namespace std;
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}

TEST_CASE( "Factorials are computed", "[factorial]" ) {
    REQUIRE( Factorial(1) == 1 );
    REQUIRE( Factorial(2) == 2 );
    REQUIRE( Factorial(3) == 6 );
    REQUIRE( Factorial(10) == 3628800 );
}
/*
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
TEST_CASE("Check if Jacobis algorithm is implemented correctly","[jacobi]"){
  //Set up a given 5*5 matrix, with known eigenvalues

  double **testMatrix=makeTestMatrix();
  double correigenvalues [5]={0.1545,0.3691,0.8115,2.0324,10.6325}; // Solution from using MatLab
  double gottenEigenvalues[5];
  jacobi_diag(testMatrix,5,1e-8);
  for(int i=0;i<5;i++){
    gottenEigenvalues[i]=testMatrix[i][i];
  }
  sort(gottenEigenvalues,gottenEigenvalues+5);
  for(int i=0;i<5;i++){
    REQUIRE(abs(gottenEigenvalues[i]-correigenvalues[i])<1e-4);
  }
}
*/
