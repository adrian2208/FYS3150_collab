#include "VMC.h"
using namespace std;
/*
double** createNMatrix(int n,int m){
  double** A;
  A = new double*[n];
  for (int i = 0; i < n; i++){
    A[i] = new double[m];
  }
  for (int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      A[i][j]=i*m;
    }
  }
  return A;
}
*/
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
void VRMonteCarlo::sample(double * energy, double * energysquared, double * time, double **posold){
  int i,particle,dim;
  int accepted=0;
  double **posnew=createNMatrix(2,3);
  double wfold=(*system).functionCart(posold);
  double wfnew;
  double local_energy=(*system).energyCart(posold);
  for(i=0;i<amount+skip;i++){ // For each iteration
    /** Suggest new position**/
    for (particle=0;particle<2;particle++){ // For each particle step
      for(dim=0;dim<3;dim++){
        posnew[particle][dim]=posold[particle][dim]+dr*rand();
      }
    }
    wfnew=(*system).functionCart(posnew);
    if (wfnew*wfnew/(wfold*wfold) > rand()){
      local_energy=(*system).energyCart(posnew);
      for (particle=0;particle<2;particle++){ // For each particle step
        for(dim=0;dim<3;dim++){
          posold[particle][dim]=posnew[particle][dim];
        }
      }
      accepted++;
      if (i>skip){
        *energy+=local_energy;
        *energysquared+=local_energy*local_energy;
      }
    }
  }
  *energy=*energy/(float)(amount);
}
double VRMonteCarlo::rand(){ // Returns a random number between -0.5 and 0.5
  return prob_dist(mt_eng)-0.5;
}
