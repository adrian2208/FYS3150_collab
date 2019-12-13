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
  /*Writes omega, alpha, beta, the energy, sigma and the distance for given parameters omega, alpha and beta, to file*/
  *outfile<<setprecision(8)  << omega<<","<<alpha <<","<< beta<<","<<energy<<","<<sigma<<","<<distance<<endl;
}
int main(int argc, char** argv){

  double omegas[4]={0.05,0.15,1.0/3.0,0.7}; //Original setup
  double alphas[4]={0.95,0.98,0.99,0.99}; //Either guesses from psi_T1, or "improved" results
  double betas[4]={0.09,0.13,0.18,0.25}; //Can be zero if guessed from psi_T2, otherwise: Best guess (or anything)
  double sigmas[4]={0,0,0,0}; //Zero at start
  double distances[4]={0,0,0,0}; //Zero at start
  //double energies[4]={std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<int>::max(),std::numeric_limits<int>::max()}; //Maximum at start so it's overwritten
  double energies[4]={0.28267408,0.7033741,1.3990892,2.7016245};
  int samplings=1e8; //Sampling size
  int skip=2e5;
  double dr=1.0;
  int counter=0;
  double alpha,beta;
  double dalpha,dbeta; // Change in alpha and beta
  double **pos=createNMatrix(2,3);pos[1][0]=1;pos[1][1]=1;pos[1][2]=1;pos[0][0]=0;pos[0][1]=0;pos[0][2]=0; //placement not based on anything
  System2 *per=new System2(1,1,1);//Initinal system. Numbers irrelvant here, will be updated anywas
  int numprocs,j; //J is "my_rank"
  MPI_Init(&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD,&j);
  string filename="../results/function2_"+to_string(omegas[j])+".csv";
  ofstream outfile;
  outfile.open(filename);//,ios::out | ios::app);
  outfile << "omega,alpha,beta,energy,sigma,distance"<<endl;
  VRMonteCarlo *vrc=new VRMonteCarlo(per, dr,samplings,skip,j);
  double energy=0,energysquared=0,time=0,distance=0,V=0,sigma=0; //
  double direction=1; //1 represents "right", -1 represents left
  bool bad_start=false; //this is important: If bad_start is true, an approximate value for beta between 0.1 and 2.0 will be found. If it's false, the beta from the list will be taken.
    if(bad_start){
    alpha=alphas[j];
      for(beta=0.1;beta<=1;beta+=0.1){
        vrc->update(alpha,beta,omegas[j]);
        vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
        if (energy<energies[j]){
          /*Find the variationally lowest beta.*/
          betas[j]=beta;
          energies[j]=energy;
          sigmas[j]=sqrt(energysquared-energy*energy);
          distances[j]=distance;
        }
        energy=energysquared=distance=0;
    }
    }
  print(omegas[j],alphas[j],betas[j],energies[j],sigmas[j] , distances[j], &outfile); //Write the start guesses to file.
  double energy_right,energy_left,energy_prev,energy_old,energy_new;
  double distance_right,distance_left,distance_prev,distance_new;
  double sigma_right,sigma_left,sigma_prev,sigma_new;
  dalpha=0.1; //Start value
  dbeta=0.1;
    do{
      if(dalpha<0.005){ //If alpha (and thereby also beta) is known by two digits precision, finish the loop
        break;
      }
      energy_old=energy_prev=energies[j];sigma_prev=sigmas[j];distance_prev=distances[j];
      vrc->update(alphas[j]+dalpha,betas[j],omegas[j]); //right side
      vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
      energy_right=energy; sigma_right=sqrt((energysquared-energy*energy));distance_right=distance;
      energy=energysquared=distance=0;
      vrc->update(alphas[j]-dalpha,betas[j],omegas[j]); //left side
      vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
      energy_left=energy; sigma_left=sqrt((energysquared-energy*energy));distance_left=distance;
      energy=energysquared=distance=0;
      energy_new=energy_right; sigma_new=sigma_right; distance_new=distance_right;
      if(energy_right>energy_left){ //Compare right side and left side
        direction=-1;
        energy_new=energy_left;sigma_new=sigma_left; distance_new=distance_left; //Take the "better side"
      }
      else{
        direction=1;
      }
      counter=1;
      while(energy_prev>energy_new){ //Require that the new found energy is lower than the previous
        counter++;
        energy_prev=energy_new;        sigma_prev=sigma_new; distance_prev=distance_new;
        vrc->update(alphas[j]+counter*direction*dalpha,betas[j],omegas[j]); // alpha, beta (not relevant for system1) and omega
        vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
        energy_new=energy;        sigma_new=sqrt((energysquared-energy*energy));   distance_new=distance;
        energy=energysquared=distance=0;
      }
      alphas[j]=alphas[j]+(counter-1)*direction*dalpha; //udpate variational parameter alpha
      energies[j]=energy_prev;  sigmas[j]=sigma_prev; distances[j]=distance_prev; //and energies, variance, distance.
      if(energies[j]!=energy_old){ //If a new, lower energy was found, write to file.
        print(omegas[j],alphas[j],betas[j],energies[j],sigmas[j] , distances[j], &outfile);
      }

      vrc->update(alphas[j],betas[j]+dbeta,omegas[j]); // right side
      vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
      energy_right=energy; sigma_right=sqrt((energysquared-energy*energy));distance_right=distance;
      energy=energysquared=distance=0;
      if(betas[j]-dbeta<1e-4){ //If the left side beta value would be zero, the calculation must not be done, and the energy is worthless. (assumed infinity)
        //dbeta*=0.5;
        energy=std::numeric_limits<float>::max();
      }
      else{
        vrc->update(alphas[j],betas[j]-dbeta,omegas[j]); //left side
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
        if(betas[j]+counter*direction*dbeta<1e-3){ //Beta can't be zero
          energy=std::numeric_limits<float>::max();
        }
        else{
          vrc->update(alphas[j],betas[j]+counter*direction*dbeta,omegas[j]); // alpha, beta (not relevant for system1) and omega
          vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
        }
        energy_new=energy;        sigma_new=sqrt((energysquared-energy*energy));   distance_new=distance;
        energy=energysquared=distance=0;
      }
      betas[j]=betas[j]+(counter-1)*direction*dbeta; //Update beta-parameter

      energies[j]=energy_prev;  sigmas[j]=sigma_prev; distances[j]=distance_prev;
      if(energy_old==energies[j]){ // In case the energy hasn't been updated at all, let the loop run again with decreased step size
        energy_old=0;
        dalpha*=0.1;
        dbeta*=0.1;
      }
      else{
        print(omegas[j],alphas[j],betas[j],energies[j],sigmas[j] , distances[j], &outfile); //otherwise, write to file
      }
    }while(fabs(energies[j]-energy_old)>sigmas[j]/sqrt(samplings)); //Rough estimate wether the "newly achieved" energy is within one standard deviation
    vrc->update(alphas[j],betas[j],omegas[j]); // calculate alpha and beta again
    vrc->sample(&energy,&energysquared,&V,&distance,&time,pos);
    print(omegas[j],alphas[j],betas[j],energy,sqrt(energysquared-energy*energy),distance, &outfile);
  delete per; delete vrc; //Free space
  outfile.close();
  cout << "Thread " << j << " finished" << endl;
  MPI_Finalize();
}
