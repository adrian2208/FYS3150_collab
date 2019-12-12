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
  int amount=4;
    double omegas[amount]={0.01,0.5,1.0,5.0};
    double kintetic_energy=0, energy=0, energysquared=0, potential_energy=0, distance=0,time=0;
    int samplings=1e7;
    int skip=2e5;
    double dr=1.0;
    double **pos=createNMatrix(2,3);pos[0][0]=1;pos[1][1]=1;pos[1][2]=1;pos[0][0]=0;pos[0][1]=0;pos[0][2]=0; //placement not based on anything
    Testsystem* per0=new Testsystem(1.0,0,omegas[0]); // No coloumb repulsion
    VRMonteCarlo* vrc0=new VRMonteCarlo(per0,dr,samplings,skip,0); //coloumb repulsion, improved system
    for(int i=0;i<amount;i++){
      vrc0->update(1.0,0,omegas[i]);
      vrc0->sample(&energy,&energysquared,&potential_energy,&distance,&time,pos);
      cout <<"omega:"<<omegas[i]<< " Energy: "<< energy<< " Potential: "<< potential_energy<<endl;
      energy=potential_energy=energysquared=distance=0;
    }

}
