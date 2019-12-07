#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
//#include "System.h"
#include "VMC.h"
#include "functions.h"
#include <fstream>
#include <limits>
using namespace std;
void print(double alpha,double beta,double energy){
  cout << "alpha: "<< alpha <<endl<< "beta: "<< beta <<endl<< "energy: " << energy<<endl;
}
int main(int argc, char** argv){
  double omegas[3]={0.01,0.5,1.0};
  double alphas[3]={0.45,0.85,0.9};
  double betas[3]={0,0,0};
  double sigmas[3]={0,0,0};
  double distances[3]={0,0,0};
  double energies[3]={std::numeric_limits<float>::max(),std::numeric_limits<int>::max(),std::numeric_limits<int>::max()};
  int samplings=1e7;
  int skip=2e5;
  double dr=1.0;
  int counter=0;
  double alpha,beta;
  double dalpha,dbeta; // Change in alpha and beta
  double **pos=createNMatrix(2,3);pos[0][0]=1;pos[1][1]=1;pos[1][2]=1;pos[0][0]=0;pos[0][1]=0;pos[0][2]=0; //placement not based on anything
  System2 *per;//=new System2(0,0,0);
  VRMonteCarlo *vrc;//=new VRMonteCarlo(&per, 0,0,0,0);// double dr, int amount, int skip, int seed
  //VRMonteCarlo(System* system, double dr, int amount, int skip, int seed=0){
  double energy=0,energysquared=0,time=0,distance=0;
  double sigma=0;
  double direction=1;
  for (int j=0;j<3;j++){ //General idea to find the lowest beta
    alpha=alphas[j];
      for(beta=0.1;beta<2;beta+=0.1){
        per=new System2(alpha,beta,omegas[j]); // alpha, beta (not relevant for system1) and omega
        vrc=new VRMonteCarlo(per, dr,samplings,skip,0); //System, dr, amount of samplings, how many skips
        vrc->sample(&energy,&energysquared,&distance,&time,pos);
        if (energy<energies[j]){
          betas[j]=beta;
          energies[j]=energy;
          sigmas[j]=sqrt(energysquared-energy*energy);
          distances[j]=distance;
        }
        energy=energysquared=distance=0;
        delete per;
        delete vrc;
    }
    if(j==2){
      cout << "omega: " << omegas[j] << endl;
      print(alphas[j],betas[j],energies[j]);
    }
  }

  double energy_right,energy_left,energy_prev,energy_old,energy_new;
  for (int j=0;j<3;j++){
    if(j<2){
      continue;
    }
    dalpha=dbeta=1;
    do{
      energy_old=energy_prev=energies[j];
      dalpha*=0.1;
      per=new System2(alphas[j]+dalpha,betas[j],omegas[j]); // alpha, beta (not relevant for system1) and omega
      vrc=new VRMonteCarlo(per, dr,samplings,skip,0); //System, dr, amount of samplings, how many skips
      vrc->sample(&energy,&energysquared,&distance,&time,pos);
      energy_right=energy;
      energy=energysquared=distance=0;
      delete per; delete vrc;
      per=new System2(alphas[j]-dalpha,betas[j],omegas[j]); // alpha, beta (not relevant for system1) and omega
      vrc=new VRMonteCarlo(per, dr,samplings,skip,0); //System, dr, amount of samplings, how many skips
      vrc->sample(&energy,&energysquared,&distance,&time,pos);
      energy_left=energy;
      energy=energysquared=distance=0;
      delete per; delete vrc;
      energy_new=energy_right;
      if(energy_right>energy_left){
        direction=-1;
        energy_new=energy_left;
      }
      else{
        direction=1;
      }
      counter=1;
      while(energy_prev>energy_new){
        counter++;
        energy_prev=energy_new;
        per=new System2(alphas[j]+counter*direction*dalpha,betas[j],omegas[j]); // alpha, beta (not relevant for system1) and omega
        vrc=new VRMonteCarlo(per, dr,samplings,skip,0); //System, dr, amount of samplings, how many skips
        vrc->sample(&energy,&energysquared,&distance,&time,pos);
        energy_new=energy;
        energy=energysquared=distance=0;
        delete per; delete vrc;
      }
      alphas[j]=alphas[j]+(counter-1)*direction*dalpha;
      energies[j]=energy_prev;
      if(j==2){
        print(alphas[j],betas[j],energies[j]);
      }
      dbeta*=0.1;
      per=new System2(alphas[j],betas[j]+dbeta,omegas[j]); // alpha, beta (not relevant for system1) and omega
      vrc=new VRMonteCarlo(per, dr,samplings,skip,0); //System, dr, amount of samplings, how many skips
      vrc->sample(&energy,&energysquared,&distance,&time,pos);
      energy_right=energy;
      energy=energysquared=distance=0;
      delete per; delete vrc;
      per=new System2(alphas[j],betas[j]-dbeta,omegas[j]); // alpha, beta (not relevant for system1) and omega
      vrc=new VRMonteCarlo(per, dr,samplings,skip,0); //System, dr, amount of samplings, how many skips
      vrc->sample(&energy,&energysquared,&distance,&time,pos);
      energy_left=energy;
      energy=energysquared=distance=0;
      delete per; delete vrc;
      energy_new=energy_right;
      if(energy_right>energy_left){
        direction=-1;
        energy_new=energy_left;
      }
      else{
        direction=1;
      }
      counter=1;
      while(energy_prev>energy_new){
        counter++;
        energy_prev=energy_new;
        per=new System2(alphas[j],betas[j]+counter*direction*dbeta,omegas[j]); // alpha, beta (not relevant for system1) and omega
        vrc=new VRMonteCarlo(per, dr,samplings,skip,0); //System, dr, amount of samplings, how many skips
        vrc->sample(&energy,&energysquared,&distance,&time,pos);
        energy_new=energy;
        energy=energysquared=distance=0;
        delete per; delete vrc;
      }
      betas[j]=betas[j]+(counter-1)*direction*dbeta;
      energies[j]=energy_prev;
      if(j==2){
        print(alphas[j],betas[j],energies[j]);
      }
    }while(fabs(energies[j]-energy_old)>1e-4);
  }
}
