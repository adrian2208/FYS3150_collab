#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
//#include "System.h"
#include "VMC.h"
#include "functions.h"

using namespace std;
int main(int argc, char** argv){
  int samplings=1000000;
  int skip=100000;
  double alpha=1.01;
  double omega=1;
  double dr=1;
  double beta;
  if (argc>=6){
    omega=atof(argv[1]);
    alpha=atof(argv[2]);
    beta=atof(argv[3]);
    dr=atof(argv[4]);
    samplings=atoi(argv[5]);
  }
  else{
    cout << "You need to give omega, alpha and beta, and dr, and the sampling size, as parameters"<<endl;
    exit(1);
  }
  double **pos=createNMatrix(2,3);pos[0][0]=1;pos[1][1]=1;pos[1][2]=1;pos[0][0]=0;pos[0][1]=0;pos[0][2]=0; //placement not based on anything
  Testsystem per=Testsystem(alpha,beta,omega); // alpha, beta (not relevant for system1) and omega
  VRMonteCarlo vrc=VRMonteCarlo(&per, dr,samplings,skip,0); //System, dr, amount of samplings, how many skips
  //VRMonteCarlo(System* system, double dr, int amount, int skip, int seed=0){
  double energy=0,energysquared=0,time=0;
  double sigma;
  for (int i=0;i<1;i++){
    vrc.sample(&energy,&energysquared,&time,pos);
    sigma=sqrt(energysquared-energy*energy);
    cout <<"Energy: "<<energy<<" Energy squared: "<<energysquared<< "sigma: "<<sigma<<endl;
    energy=energysquared=0;
  }
}
