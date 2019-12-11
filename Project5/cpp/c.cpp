#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
//#include "System.h"
#include "VMC.h"
#include "functions.h"
#include <fstream>
using namespace std;
void energy_equilibrium(){
  ofstream outfile;
  outfile.open("../results/energy_distribution.csv");
  double energy=0,energysquared=0,time=0;
  double alpha=0.88;
  int skip=2e5;
  int amount=5e7;
  double dr=1.3;
  double beta=0;
  double omega=1.0;
  double distance=0;
  double **pos=createNMatrix(2,3);pos[0][0]=1;pos[1][1]=1;pos[1][2]=1;pos[0][0]=0;pos[0][1]=0;pos[0][2]=0; //placement not based on anything
  System1 *per;//=System1(0,0,0);
  VRMonteCarlo *vrc;//=VRMonteCarlo(&per, 0,0,0,0);;
  per=new System1(alpha,beta,omega); // alpha, beta (not relevant for system1) and omega
  vrc=new VRMonteCarlo(per, dr,amount,skip,0); //System, dr, amount of samplings, how many skips
  double *energies=vrc->sample_detailed(&energy,&energysquared,&distance,&time,pos);
  outfile << "Index,Energy,omega,alpha"<<endl;
  for(int i=1000;i<amount;i+=1000){
    outfile << i << ","<<energies[i]<<","<<omega<<","<<alpha<<endl;
  }
  delete per;
  delete vrc;
}
int main(int argc, char** argv){
  energy_equilibrium();
  ofstream outfile;
  outfile.open("../results/function1.csv",ios::out | ios::app); //Collected results (Numbers)
  double alpha_start=0.2;
  double alpha_end=1.2;
  double d_alpha=0.01;
  int n=(int)((alpha_end-alpha_start)/d_alpha)+1;
  double omega[3]={0.01,0.5,1.0};//;,5.0};
  double * energy_list=new double[3*n];
  double * sigma_list=new double[3*n];
  double *distance_list=new double[3*n];
  int samplings=1e8;
  int skip=2e5;
  double dr=1.0;
  double beta=0;
  double **pos=createNMatrix(2,3);pos[0][0]=1;pos[1][1]=1;pos[1][2]=1;pos[0][0]=0;pos[0][1]=0;pos[0][2]=0; //placement not based on anything
  System1* per;//= new System1(0,0,0);
  VRMonteCarlo* vrc;///= new VRMonteCarlo(per, 0,0,0,0);
  //VRMonteCarlo(System* system, double dr, int amount, int skip, int seed=0){
  double energy=0,energysquared=0,time=0,distance=0,V=0;
  double sigma;
  for (int j=0;j<3;j++){
    for (int i=0;i<n;i++){
      per=new System1(alpha_start+i*d_alpha,beta,omega[j]); // alpha, beta (not relevant for system1) and omega
      vrc=new VRMonteCarlo(per, dr,samplings,skip,0); //System, dr, amount of samplings, how many skips
      vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
      sigma=sqrt(energysquared-energy*energy);
      energy_list[j*n+i]=energy;
      sigma_list[j*n+i]=sigma;
      distance_list[j*n+i]=distance;
      outfile<<omega[j] <<","<<alpha_start+i*d_alpha << ","<< energy << "," << sigma<<","<<distance<<endl;
      cout <<"Omega: "<<omega[j]<<" alpha: " <<alpha_start+i*d_alpha<< " Energy: "<<energy<<" Energy squared: "<<energysquared<< " sigma: "<<sigma<< "distance: "<<distance<<endl;
      energy=energysquared=distance=0;
      delete per;
      delete vrc;
    }
  }
  outfile.close();
}
