#include "functions.hpp"
#include "lib.h"
#include <algorithm>
#include <random>
#include <cmath>
#include "catch.hpp"
using namespace std;
double example_functions(double x,double y){
  return x*y;
}
TEST_CASE("Test whether the correct numbers are given by the getParallelizationCoefficients-function for more than one loop"){
  int N=8; //a fictive 105 runs in total
  int n=3; //a fictive of 3 threads
  int loop=2; // a fictive of 2 loops
  long long * val=getParallelizationCoefficients(N,2,n,loop);
  int val_expected[n]={3,5,21};
  for(int i=0;i<(loop+1);i++){
    REQUIRE(val_expected[i]==val[i]);
  }
}
TEST_CASE("Test whether the correct numbers are given by the getParallelizationCoefficients-function for one loop"){
  int N=503; //a fictive 105 runs in total
  int n=3; //a fictive of 3 threads
  int loop=1; // a fictive of 2 loops
  long long * val=getParallelizationCoefficients(N,1,n,1);
  int val_expected[n]={168,168};

  for(int i=0;i<(loop+1);i++){
    REQUIRE(val_expected[i]==val[i]);
  }
}
TEST_CASE("Tests a simple integral,xydxdy from 0 to 5 for both integrands, using legrendre polynoms"){
  int N=100;
  double *x = new double[N];
  double *w=new double[N];
  double add_var=0;;
  gauleg(0,5,x,w,N);
  for (int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      add_var+=w[i]*w[j]*example_functions(x[i],x[j]);
    }
  }
  REQUIRE(fabs((add_var-625.0/4.0))<1e-3);
}
void moveToLeft(double *arr,int N){
  for (int i=0;i<N;i++){
    arr[i]=arr[i+1];
  }
}
TEST_CASE("Tests a simple integral,9x*y*e**(-x-3y)dxdy from 0 to infinity for both integrands, using laguerre polynoms"){
  int N=100;
  double *x = new double[N+1];
  double *w=new double[N+1];
  double add_var=0;
  moveToLeft(x,N);
  moveToLeft(w,N);
  gauss_laguerre(x,w,N,1);
  for (int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      add_var+=w[i]*w[j]*(9*exp(-2*x[i]));
    }
  }
  REQUIRE(fabs((add_var-1.0))<1e-3);
}
TEST_CASE("Tests a simple integral,xydxdy from 0 to 5 for both integrands, using simple monte carlo."){
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> RnG(0.0,1.0);
  int N=1e7;
  double x, y;
  double jacobi=25.0;
  double sum=0;
  for(int i=0;i<N;i++){
    x=5.0*RnG(gen);
    y=5.0*RnG(gen);
    sum+=example_functions(x,y);
  }
  sum=sum*jacobi/N;
  REQUIRE(fabs((sum-625.0/4.0))<1e-1);
}
