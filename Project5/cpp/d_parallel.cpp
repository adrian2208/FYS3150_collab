#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
//#include "System.h"
#include "VMC.h"
#include "functions.h"
#include <fstream>
#include <limits>
#include <mpi.h>
#include <iomanip>
using namespace std;
void print(double omega,double alpha,double beta,double energy, double sigma, double distance, ofstream * outfile){
  *outfile<<setprecision(8)  << omega<<","<<alpha <<","<< beta<<","<<energy<<","<<sigma<<","<<distance<<endl;
}
int main(int argc, char** argv){

  double omegas[4]={0.01,0.5,1.0,5.0};
  double alphas[4]={0.45,0.84,0.88,0.95};
  double betas[4]={0.01,0.14,0.31,0.80};
  double sigmas[4]={0,0,0,0};
  double distances[4]={0,0,0,0};
  double energies[4]={std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<int>::max(),std::numeric_limits<int>::max()};
  //double energies[4]={0.5,1.9776574,3.7081356,16.686586};
  int samplings=1e7;
  int skip=2e5;
  double dr=1.0;
  int counter=0;
  double alpha,beta;
  double dalpha,dbeta; // Change in alpha and beta
  double **pos=createNMatrix(2,3);pos[1][0]=1;pos[1][1]=1;pos[1][2]=1;pos[0][0]=0;pos[0][1]=0;pos[0][2]=0; //placement not based on anything
  System2 *per=new System2(1,1,1);//=new System2(0,0,0);
  int numprocs,my_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD,&my_rank);
  double time_start=MPI_Wtime();
  string filename="../results/function2_"+to_string(omegas[my_rank])+".csv";
  ofstream outfile;
  outfile.open(filename,ios::out | ios::app);
  //outfile << "omega,alpha,beta,energy,sigma,distance"<<endl;
  VRMonteCarlo *vrc=new VRMonteCarlo(per, dr,samplings,skip,my_rank);//=new VRMonteCarlo(&per, 0,0,0,0);// double dr, int amount, int skip, int seed
  //VRMonteCarlo(System* system, double dr, int amount, int skip, int seed=0){
  double energy=0,energysquared=0,time=0,distance=0,V=0;
  double sigma=0;
  double direction=1;
  int j=my_rank;
  bool bad_start=true;
    if(bad_start){
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
    }
    }
  print(omegas[j],alphas[j],betas[j],energies[j],sigmas[j] , distances[j], &outfile);
  double energy_right,energy_left,energy_prev,energy_old,energy_new;
  double distance_right,distance_left,distance_prev,distance_new;
  double sigma_right,sigma_left,sigma_prev,sigma_new;
  dalpha=0.1;
  dbeta=0.1;
    do{
      if(dalpha<0.005 && omegas[j]>0.1){
        break;
      }
      if(dbeta<0.005 && omegas[j]>0.1){
        break;
      }
      energy_old=energy_prev=energies[j];sigma_prev=sigmas[j];distance_prev=distances[j];
      vrc->update(alphas[j]+dalpha,betas[j],omegas[j]);
      vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
      energy_right=energy; sigma_right=sqrt((energysquared-energy*energy));distance_right=distance;
      energy=energysquared=distance=0;
      vrc->update(alphas[j]-dalpha,betas[j],omegas[j]); // alpha, beta (not relevant for system1) and omega
      vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
      energy_left=energy; sigma_left=sqrt((energysquared-energy*energy));distance_left=distance;
      energy=energysquared=distance=0;
      energy_new=energy_right; sigma_new=sigma_right; distance_new=distance_right;
      if(energy_right>energy_left){
        direction=-1;
        energy_new=energy_left;sigma_new=sigma_left; distance_new=distance_left;
      }
      else{
        direction=1;
      }
      counter=1;
      while(energy_prev>energy_new){
        counter++;
        energy_prev=energy_new;        sigma_prev=sigma_new; distance_prev=distance_new;
        vrc->update(alphas[j]+counter*direction*dalpha,betas[j],omegas[j]); // alpha, beta (not relevant for system1) and omega
        vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
        energy_new=energy;        sigma_new=sqrt((energysquared-energy*energy));   distance_new=distance;
        energy=energysquared=distance=0;
      }
      alphas[j]=alphas[j]+(counter-1)*direction*dalpha;
      energies[j]=energy_prev;  sigmas[j]=sigma_prev; distances[j]=distance_prev;
      if(energies[j]!=energy_old){
        print(omegas[j],alphas[j],betas[j],energies[j],sigmas[j] , distances[j], &outfile);
      }

      vrc->update(alphas[j],betas[j]+dbeta,omegas[j]); // alpha, beta (not relevant for system1) and omega
      vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
      energy_right=energy; sigma_right=sqrt((energysquared-energy*energy));distance_right=distance;
      energy=energysquared=distance=0;
      if(betas[j]-dbeta<1e-4){ //If the next beta value would be zero, the calculation must not be done, and the energy is worthless.
        //dbeta*=0.5;
        energy=std::numeric_limits<float>::max();
      }
      else{
        vrc->update(alphas[j],betas[j]-dbeta,omegas[j]); // alpha, beta (not relevant for system1) and omega
        vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
      }
      energy_left=energy; sigma_left=sqrt((energysquared-energy*energy));distance_left=distance;
      energy=energysquared=distance=0;
      energy_new=energy_right;
      if(energy_right>energy_left){
        direction=-1;
        energy_new=energy_left;sigma_new=sigma_left; distance_new=distance_left;
      }
      else{
        direction=1;
      }
      counter=1;
      while(energy_prev>energy_new){
        counter++;
        energy_prev=energy_new;        sigma_prev=sigma_new; distance_prev=distance_new;
        if(betas[j]+counter*direction*dbeta<1e-3){
          energy=std::numeric_limits<float>::max();
        }
        else{
          vrc->update(alphas[j],betas[j]+counter*direction*dbeta,omegas[j]); // alpha, beta (not relevant for system1) and omega
          vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
        }
        energy_new=energy;        sigma_new=sqrt((energysquared-energy*energy));   distance_new=distance;
        energy=energysquared=distance=0;
      }
      betas[j]=betas[j]+(counter-1)*direction*dbeta;

      energies[j]=energy_prev;  sigmas[j]=sigma_prev; distances[j]=distance_prev;
      if(energy_old==energies[j]){ // In case the energy hasn't been updated, let the loop rn again
        energy_old=0;
        dalpha*=0.1;
        dbeta*=0.1;
      }
      else{
        print(omegas[j],alphas[j],betas[j],energies[j],sigmas[j] , distances[j], &outfile);
      }
      if(betas[j]<0.15){
        dbeta=0.01;
      }
      if(betas[j]<0.015){
        dbeta=0.001;
      }
    if(alphas[j]>1 && dalpha<0.1 && omegas[j]<0.05){
      dalpha=0.1;
    }
    }while(fabs(energies[j]-energy_old)>sigmas[j]/sqrt(samplings));
    vrc->update(alphas[j],betas[j],omegas[j]); // alpha, beta (not relevant for system1) and omega
    vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
    print(omegas[j],alphas[j],betas[j],energy,sqrt(energysquared-energy*energy),distance, &outfile);
  delete per; delete vrc;
  outfile.close();
  cout << "Thread " << j << " finished" << endl;
  MPI_Finalize();
}
