#include "VMC.h"
VRMonteCarlo::VRMonteCarlo(System* system, double dr, int amount, int skip){
      this->system=system;
      this->dr=dr;
      this->amount=amount;
      this->skip=skip;
    }
