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
    int skip, skip_orig;
    double dr, dr_orig;
    std::mt19937 mt_eng;
    std::uniform_real_distribution<double> prob_dist;
    double rand();
  public:
    VRMonteCarlo(System * system, double dr, int amount, int skip=2e5, int seed=0);
    void sample(double * energy, double * energysquared,double *v,double *distance_, double * time, double **posold);
    double * sample_detailed(double * energy, double * energysquared,double *distance_, double * time, double **posold);
    void update(double alpha, double beta, double omega);
};
#endif
