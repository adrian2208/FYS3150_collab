#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
//#include "System.h"
#include "VMC.h"
#include "functions.h"
#include "catch.hpp"
using namespace std;
TEST_CASE("Test that alpha=1.0 is the ideal solution when there's no repulsion"){
  double alpha=1.0, beta=0,omega=1,time=0,dr=1;int samplings=1e7,skip=2e5;
  double **pos=createNMatrix(2,3);pos[0][0]=1;pos[1][1]=1;pos[1][2]=1;pos[0][0]=0;pos[0][1]=0;pos[0][2]=0; //placement not based on anything
  Testsystem *per=new Testsystem(alpha,beta,omega); // alpha, beta (not relevant for system1) and omega
  VRMonteCarlo *vrc= new VRMonteCarlo(per, dr,samplings,skip,0); //System, dr, amount of samplings, how many skips
  //VRMonteCarlo(System* system, double dr, int amount, int skip, int seed=0){
  double energy_alpha=0,energysquared_alpha=0,distance=0,V=0;
  double sigma_alpha;
    vrc->sample(&energy_alpha,&energysquared_alpha,&V,&distance,&time,pos);
    sigma_alpha=sqrt(energysquared_alpha-energy_alpha*energy_alpha);
    REQUIRE(fabs(energy_alpha-3)<1e-5); //Energy should be 3
    REQUIRE(sigma_alpha<1e-3); // Sigma should be teeny-tiny
  double energy=0, energySquared=0,sigma=0,V=0;
  for (double alpha=0.1;alpha<2;alpha+=0.2){
    per=Testsystem(alpha,beta,omega);
    vrc=VRMonteCarlo(&per, dr,samplings);
    vrc.sample(&energy,&energySquared,&distance,&V,&time,pos);
    sigma=sqrt(energySquared-energy*energy);
    REQUIRE(energy>energy_alpha); //Energy should be more than 3
    REQUIRE(sigma>sigma_alpha); // Sigma should be larger than the accurate one
    energy=0,energySquared=0,sigma=0,V=0;
  }
}
TEST_Case("Test the virial Theorem for no-repulsion-systems")
