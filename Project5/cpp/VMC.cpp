#include "VMC.h"
using namespace std;
VRMonteCarlo::VRMonteCarlo(System* system, double dr, int amount, int skip, int seed=0){
      this->system=system;
      this->dr=dr;
      this->amount=amount;
      this->skip=skip;
      //mt_eng{std::random_device{}()}, prob_dist(0.0, 1.0) {}
      random_device rd;
      mt_eng= mt19937(rd()+seed); //Each thread gets a different seed, as the rank is included
      prob_dist= uniform_real_distribution<double>(0.0,1.0);
    }
void VRMonteCarlo::sample(double * energy, double * energysquared, double * time){
  cout << (*system).energy(rand(),rand(),rand())<<endl;
}
double VRMonteCarlo::rand(){
  return prob_dist(mt_eng);
}
