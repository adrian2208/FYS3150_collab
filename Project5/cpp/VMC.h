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
    std::mt19937 mt_eng; //Random number generator
    std::uniform_real_distribution<double> prob_dist; //Distribution
    double rand();
  public:
    VRMonteCarlo(System * system, double dr, int amount, int skip=2e5, int seed=0); //Constructor
    void sample(double * energy, double * energysquared,double *v,double *distance_, double * time, double **posold); //Samplef
    double * sample_detailed(double * energy, double * energysquared,double *distance_, double * time, double **posold); //sample that retourns local energies
    void update(double alpha, double beta, double omega); //Updates the inner system.
};
#endif
