#ifndef VMC_H
#define VMC_H
#include <random>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "System.h"
#include "functions.h"

#include <ctime>
class VRMonteCarlo{
  private:
    System *system;
    int amount;
    int skip;
    double dr;
    std::mt19937 mt_eng;
    std::uniform_real_distribution<double> prob_dist;
    double rand();
  public:
    VRMonteCarlo(System * system, double dr, int amount, int skip, int seed);
    void sample(double * energy, double * energysquared, double * time, double **pos);
};
#endif
