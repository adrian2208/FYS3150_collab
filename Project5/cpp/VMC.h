#ifndef VMC_H
#define VMC_H
#include <random>
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
    std::mt19937 gen;
    std::uniform_real_distribution<double> RnG;
  public:
    VRMonteCarlo(System * system, double dr, int amount, int skip, int seed);
    void sample(double * energy, double * energysquared, double * time);
};
#endif
