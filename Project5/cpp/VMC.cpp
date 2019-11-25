#include "VMC.h"
using namespace std;
VRMonteCarlo::VRMonteCarlo(System* system, double dr, int amount, int skip, int seed=0){
      this->system=system;
      this->dr=dr;
      this->amount=amount;
      this->skip=skip;
      std::random_device rd;
      std::mt19937_64 gen(rd()+seed); //Each thread gets a different seed, as the rank is included
      std::uniform_real_distribution<double> RnG(0.0,1.0);
    }
void VRMonteCarlo::sample(double * energy, double * energysquared, double * time){
  cout << RnG(gen)<<endl;
}
