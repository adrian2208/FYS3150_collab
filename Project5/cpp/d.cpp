#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "VMC.h"
#include "functions.h"
#include <fstream>
#include <limits>
using namespace std;
/*Outdated, don't use*/
void print(double alpha,double beta,double energy){
  cout << "alpha: "<< alpha <<endl<< "beta: "<< beta <<endl<< "energy: " << energy<<endl;
}
int main(int argc, char** argv){
  ofstream outfile;
  outfile.open("../results/function2.csv",ios::out | ios::app);

  double omegas[3]={0.01,0.5,1.0};
  double alphas[3]={1.32,0.84,0.88};
  double betas[3]={0.03,0,0};
  double sigmas[3]={0,0,0};
  double distances[3]={0,0,0};
  double energies[3]={std::numeric_limits<float>::max(),std::numeric_limits<int>::max(),std::numeric_limits<int>::max()};
  int samplings=5e8;
  int skip=2e5;
  double dr=1.0;
  int counter=0;
  double alpha,beta;
  double dalpha,dbeta; // Change in alpha and beta
  double **pos=createNMatrix(2,3);pos[1][0]=1;pos[1][1]=1;pos[1][2]=1;pos[0][0]=0;pos[0][1]=0;pos[0][2]=0; //placement not based on anything
  System2 *per=new System2(1,1,1);//=new System2(0,0,0);
  VRMonteCarlo *vrc=new VRMonteCarlo(per, dr,samplings,skip,0);//=new VRMonteCarlo(&per, 0,0,0,0);// double dr, int amount, int skip, int seed
  //VRMonteCarlo(System* system, double dr, int amount, int skip, int seed=0){
  double energy=0,energysquared=0,time=0,distance=0,V=0;
  double sigma=0;
  double direction=1;
  for (int j=2;j<3;j++){ //General idea to find the lowest beta
    alpha=alphas[j];
      for(beta=0.1;beta<2;beta+=0.1){
        vrc->update(alpha,beta,omegas[j]);
        vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
        if (energy<energies[j]){
          betas[j]=beta;
          energies[j]=energy;
          sigmas[j]=sqrt(energysquared-energy*energy);
          distances[j]=distance;
        }
        energy=energysquared=distance=0;
        //delete per;
        //delete vrc;
    }
    if(j==2){
      cout << "omega: " << omegas[j] << endl;
      print(alphas[j],betas[j],energies[j]);
    }
  }

  double energy_right,energy_left,energy_prev,energy_old,energy_new;
  double distance_right,distance_left,distance_prev,distance_new;
  double sigma_right,sigma_left,sigma_prev,sigma_new;
  for (int j=2;j<3;j++){
    dalpha=dbeta=0.1;
    do{
      energy_old=energy_prev=energies[j];sigma_prev=sigmas[j];distance_prev=distances[j];
      vrc->update(alphas[j]+dalpha,betas[j],omegas[j]);
      vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
      energy_right=energy; sigma_right=sqrt((energy*energy-energysquared));distance_right=distance;
      energy=energysquared=distance=0;
      vrc->update(alphas[j]-dalpha,betas[j],omegas[j]); // alpha, beta (not relevant for system1) and omega
      vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
      energy_left=energy; sigma_left=sqrt((energy*energy-energysquared));distance_left=distance;
      energy=energysquared=distance=0;
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
        sigma_prev=sigma_new;
        vrc->update(alphas[j]+counter*direction*dalpha,betas[j],omegas[j]); // alpha, beta (not relevant for system1) and omega
        vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
        energy_new=energy;
        sigma_new=sqrt((energy*energy-energysquared));
        energy=energysquared=distance=0;
      }
      alphas[j]=alphas[j]+(counter-1)*direction*dalpha;
      energies[j]=energy_prev;
      if(j==2){
        print(alphas[j],betas[j],energies[j]);
      }
      if(fabs(betas[j]-0.01)<1e-3){
        dbeta=0.001;
      }
      vrc->update(alphas[j],betas[j]+dbeta,omegas[j]); // alpha, beta (not relevant for system1) and omega
      vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
      energy_right=energy;
      energy=energysquared=distance=0;
      vrc->update(alphas[j],betas[j]-dbeta,omegas[j]); // alpha, beta (not relevant for system1) and omega
      vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
      energy_left=energy;
      energy=energysquared=distance=0;
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
        vrc->update(alphas[j],betas[j]+counter*direction*dbeta,omegas[j]); // alpha, beta (not relevant for system1) and omega
        vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
        energy_new=energy;
        energy=energysquared=distance=0;
      }
      betas[j]=betas[j]+(counter-1)*direction*dbeta;
      energies[j]=energy_prev;
      if(j==2){
        print(alphas[j],betas[j],energies[j]);
      }
      if(energy_old==energies[j]){ // In case the energy hasn't been updated, let the loop rn again
        energy_old=0;
        dalpha*=0.1;
      }
    }while(fabs(energies[j]-energy_old)>sigmas[j]/sqrt(samplings));
  outfile << omegas[j]<<","<<alphas[j]<<","<<betas[j]<<energies[j]<<","<<sigmas[j]<<endl;
  }
}
