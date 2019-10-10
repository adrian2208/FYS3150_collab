#include "lib.h"
#include <algorithm>
#include <random>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <functional>
#define PI 3.14159265358979

using namespace std;
/*double int_func(double x1, double x2, double y1, double y2, double z1, double z2){
  double exponential_part=exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2))); // Exponential part of the function
  double denom=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)); // Denominator part
  double returnval;
  //cout << denom << endl;
  if (fabs(denom)<= 1e-15){ //If denominator is zero, it is omitted
    returnval=0;
  }
  else{
    returnval= exponential_part/denom; //Final value of function
  }
  return returnval;
}
double int_func(double r1, double r2, double theta1, double theta2, double phi1, double phi2){
  double cosb=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2); // the actor
  double top=r1*r1*r2*r2*sin(theta1)*sin(theta2)*exp(-4*(r1+r2));
  double denum;
  double root=r1*r1+r2*r2-2*r1*r2*cosb;
  if (root >1.0e-8){
    denum=sqrt(root);
  }
  else{
    return 0;
  }
  return top/denum;
}*/
double int_func_montecarlo(double r1, double r2, double theta1, double theta2, double phi1, double phi2){
  double cosb=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2); // the actor
  double top=r1*r1*r2*r2*sin(theta1)*sin(theta2);
  double root=r1*r1+r2*r2-2*r1*r2*cosb;
  double denum;
  if (root >1.0e-8){
    denum=sqrt(root);
  }
  else{
    return 0;
  }
  return top/denum;
}

int main(int argc, char** argv){
  int N;
  if(argc>=2){
    N=atoi(argv[1]); //n needs to be stated
  }
  else{
    cout << "You need to state a number N" << endl;
    exit(1);
  }
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> RnG(0.0,1.0);
  double int_mc=0.0; double variance=0.0;
  double sum_sigma=0.0; double sum_mc=0.0;
  double x[6];
  double fx=0;
  double jacobi_det=4.0*pow(PI,4)/16.0;//pow(2*PI,2)*pow(PI,2);'
  double rand1;
  double rand2;
  for (int i=0;i<N;i++){
    for(int j=0;j<2;j++){
      x[j]=-0.25*log(1.0-RnG(gen));
    }
    for(int j=2;j<4;j++){
      x[j]=PI*RnG(gen); // makes x[j] a random number between 0 and PI
    }
    for(int j=4;j<6;j++){
      x[j]=2*PI*RnG(gen); // makes x[j] a random number between 0 and 2PI
    }
    fx=int_func_montecarlo(x[0],x[1],x[2],x[3],x[4],x[5]);
    int_mc+=fx;
    sum_sigma+=fx*fx;
  }
  int_mc=int_mc/N;
  sum_sigma=sum_sigma/N;
  variance=sum_sigma-int_mc*int_mc;
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout << " Monte carlo result= " << setw(10) << setprecision(8) << jacobi_det*int_mc;
  cout << " Sigma= " << setw(10) << setprecision(8) << jacobi_det*sqrt(variance/((double) N )) << endl;
}
