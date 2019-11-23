#ifndef VMC_H
#define VMC_H
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "System.h"
class VRMonteCarlo{
  private:
    System *system;
    int amount;
    int skip;
    double dr;
  public:
    VRMonteCarlo(System * system, double dr, int amount, int skip);
};
#endif
