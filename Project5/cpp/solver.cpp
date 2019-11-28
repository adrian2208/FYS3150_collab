#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "System.h"
#include "VMC.h"
#include "functions.h"

using namespace std;
int main(int argc, char** argv){
  double **pos=createNMatrix(2,3);
  System1 per=System1(1.0,2.0,1.0);
  VRMonteCarlo vrc=VRMonteCarlo(&per, 0.5,10000,10,5);
  double energy=0,energysquared=0,time=0;
  for (int i=0;i<=100;i++){
    vrc.sample(&energy,&energysquared,&time,pos);
    cout << energy<<endl;
  }
}
