/*
Compile and run as
 c++ -c virial_theorem.cpp
 c++ -o virial_theorem.exe virial_theorem.o VMC.o functions.o System.o
./virial_theorem.exe
*/
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "VMC.h"
#include "functions.h"
#include <fstream>
#include <limits>
#include <iomanip>
using namespace std;
int main(int argc, char** argv){
  int amount=8;
  double omegas[amount]={0.01,0.05,0.15,1.0/3.0,0.5,0.7,1.0,5.0}; //Omega-values
  double alphas_1[amount]={0.45,0.63,0.74,0.81,0.84,0.87,0.88,0.95}; //Alpha values from T1
  double alphas_2[amount]={0.89,0.95,0.98,0.99,0.99,0.99,1.0,1.00};//Alpha values from T2
  double betas[amount]={0.05,0.09,0.13,0.14,0.17,0.24,0.27,0.59}; //Beta values from T2
  double kintetic_energy=0, energy=0, energysquared=0, potential_energy=0, distance=0,time=0;
  int samplings=1e7;
  int skip=2e5;
  double dr=1.0;
  double kinetic_energies[amount*3]; //For the three functions run
  double potential_energies[amount*3];
  double **pos=createNMatrix(2,3);pos[0][0]=1;pos[1][1]=1;pos[1][2]=1;pos[0][0]=0;pos[0][1]=0;pos[0][2]=0; //placement not based on anything
  Testsystem* per0=new Testsystem(1.0,0,omegas[0]); // No coloumb repulsion
  VRMonteCarlo* vrc0=new VRMonteCarlo(per0,dr,samplings,skip,0); //coloumb repulsion, improved system
  System1 * per1=new System1(1,0,omegas[0]); //coloumb repulsion, bad system
  VRMonteCarlo* vrc1=new VRMonteCarlo(per1,dr,samplings,skip,0);
  System* per2=new System2(1.0,0,omegas[0]);//coloumb repulsion, improved system
  VRMonteCarlo* vrc2=new VRMonteCarlo(per2,dr,samplings,skip,0);
  for (int i=0;i<amount;i++){
    pos[0][0]=1;pos[1][1]=1;pos[1][2]=1;pos[0][0]=0;pos[0][1]=0;pos[0][2]=0;
    vrc0->update(1.0,0,omegas[i]);
    vrc1->update(alphas_1[i],0,omegas[i]);
    vrc2->update(alphas_2[i],betas[i],omegas[i]);
    vrc0->sample(&energy,&energysquared,&potential_energy,&distance,&time,pos);
    kintetic_energy=energy-potential_energy; //Calculate kinetic energy
    cout <<"Energy:"<<energy<< " Potential: "<<potential_energy << " Kinetic: " << kintetic_energy<<endl;
    kinetic_energies[i]=kintetic_energy;potential_energies[i]=potential_energy;
    kintetic_energy=potential_energy=energy=distance=energysquared=time=0; //Zero out
    vrc1->sample(&energy,&energysquared,&potential_energy,&distance,&time,pos);
    kintetic_energy=energy-potential_energy;
    kinetic_energies[amount+i]=kintetic_energy;potential_energies[amount+i]=potential_energy;
    kintetic_energy=potential_energy=energy=distance=energysquared=time=0;
    vrc2->sample(&energy,&energysquared,&potential_energy,&distance,&time,pos);
    kintetic_energy=energy-potential_energy;
    kinetic_energies[2*amount+i]=kintetic_energy;potential_energies[2*amount+i]=potential_energy;
    kintetic_energy=potential_energy=energy=distance=energysquared=time=0;
  }
  ofstream outfile;
  outfile.open("../results/virial.csv");
  outfile << "omega,kinetic_energy,potential_energy"<<endl;
  for(int i=0;i<3;i++){
    for(int j=0;j<amount;j++){
      outfile << omegas[j]<<","<<kinetic_energies[i*amount+j]<<","<<potential_energies[i*amount+j]<<endl;
    }
  }
  outfile.close();
}
