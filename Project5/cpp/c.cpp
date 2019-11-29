#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
//#include "System.h"
#include "VMC.h"
#include "functions.h"
#include <fstream>
using namespace std;
int main(int argc, char** argv){
  ofstream outfile;
  outfile.open("../results/function1.csv",ios::out | ios::app); //Collected results (Numbers)
  double alpha_start=0.1;
  double alpha_end=2.0;
  double d_alpha=0.05;
  int n=(int)((alpha_end-alpha_start)/d_alpha)+1;
  double omega[3]={0.01,0.5,1.0};
  double * energy_list=new double[3*n];
  double * sigma_list=new double[3*n];
  double *distance=new double[3*n];
  int samplings=1e7;
  int skip=100000;
  double dr=1.0;
  double beta=0;
  double **pos=createNMatrix(2,3);pos[0][0]=1;pos[1][1]=1;pos[1][2]=1;pos[0][0]=0;pos[0][1]=0;pos[0][2]=0; //placement not based on anything
  System1 per=System1(0,0,0);
  VRMonteCarlo vrc=VRMonteCarlo(&per, 0,0,0,0);;
  //VRMonteCarlo(System* system, double dr, int amount, int skip, int seed=0){
  double energy=0,energysquared=0,time=0;
  double sigma;
  for (int j=0;j<3;j++){
    for (int i=0;i<n;i++){
      per=System1(alpha_start+i*d_alpha,beta,omega[j]); // alpha, beta (not relevant for system1) and omega
      vrc=VRMonteCarlo(&per, dr,samplings,skip,0); //System, dr, amount of samplings, how many skips
      vrc.sample(&energy,&energysquared,&time,pos);
      sigma=sqrt(energysquared-energy*energy);
      energy_list[j*n+i]=energy;
      sigma_list[j*n+i]=sigma;
      outfile<<omega[j] <<","<<alpha_start+i*d_alpha << ","<< energy << "," << sigma<<endl;
      cout <<"Energy: "<<energy<<" Energy squared: "<<energysquared<< "sigma: "<<sigma<<endl;
      energy=energysquared=0;
    }
  }
  outfile.close();
}
