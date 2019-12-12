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

VRMonteCarlo::VRMonteCarlo(System* system, double dr, int amount, int skip, int seed){
      this->system=system;
      this->dr=dr; this->dr_orig=dr;
      this->amount=amount;
      this->skip=skip; this->skip_orig=skip;
      //mt_eng{std::random_device{}()}, prob_dist(0.0, 1.0) {}
      random_device rd;
      mt_eng= mt19937(rd()+seed); //Each thread gets a different seed, as the rank is included
      prob_dist= uniform_real_distribution<double>(0.0,1.0);
    }
void VRMonteCarlo::update(double alpha, double beta, double omega){
  skip=skip_orig;
  dr=dr_orig;
  system->update(alpha,beta,omega);
}
double * VRMonteCarlo::sample_detailed(double * energy, double * energysquared, double *distance_,double * time, double **posold){
  int original_skip=skip;
  int i,particle,dim;
  int accepted=0;
  double **posnew=createNMatrix(2,3);
  double wfold=(*system).functionCart(posold);
  double wfnew;
  double local_energy=system->energyCart(posold);
  double accepted_moves_last10000=0;
  bool equilibrium=false;
  double *e_avg=new double[amount];
  int check_size=10000;
  if (skip<check_size){
    skip=check_size;
  }
  for(i=0;i<amount+skip;i++){ // For each iteration
    /** Suggest new position**/
    for (particle=0;particle<2;particle++){ // For each particle step
      for(dim=0;dim<3;dim++){
        posnew[particle][dim]=posold[particle][dim]+dr*rand();
      }
    }
    wfnew=(*system).functionCart(posnew);
    if (wfnew*wfnew/(wfold*wfold) > rand()+0.5){
      accepted_moves_last10000++;
      if (i>skip){
        accepted++;
      }
      local_energy=(*system).energyCart(posnew);
      for (particle=0;particle<2;particle++){ // For each particle step
        for(dim=0;dim<3;dim++){
          posold[particle][dim]=posnew[particle][dim];
        }
      }
      wfold=wfnew;
    }
    if (i%check_size==0 && i!=0){ //each 10.000th time
      if (accepted_moves_last10000>0.6*check_size && i!=0){
        (this->dr)*=1.1;
        if(!equilibrium){
          skip+=check_size;
        }
      }
      else if (accepted_moves_last10000<0.4*check_size && i!=0){
        (this->dr)*=0.9;
        if(!equilibrium){
          skip+=check_size;
        }
      }
      else{
        equilibrium=true;
      }
      accepted_moves_last10000=0;
    }
    if (i>skip && equilibrium){
      *energy+=local_energy;
      *energysquared+=local_energy*local_energy;
      e_avg[i-skip]=*energy/(float)(i-skip);
    }
  }
  *energy=*energy/(float)(amount);
  *energysquared=*energysquared/(float)(amount);
  skip=original_skip;
  cout <<"accepted moves: " << accepted << endl;
  return e_avg;
}
void VRMonteCarlo::sample(double * energy, double * energysquared,double *v,double *distance_,double * time, double **posold){
  posold[1][0]=1;posold[1][1]=1;posold[1][2]=1;posold[0][0]=0;posold[0][1]=0;posold[0][2]=0;
  int original_skip=skip;
  int i,particle,dim;
  int accepted=0;
  double **posnew=createNMatrix(2,3);
  double wfold=(*system).functionCart(posold);
  double wfnew;
  double local_energy=system->energyCart(posold);
  double potential=system->potentialCart(posold);
  double accepted_moves_last10000=0;
  bool equilibrium=false;
  int check_size=10000;
  if (skip<check_size){
    skip=check_size;
  }
  for(i=0;i<amount+skip;i++){ // For each iteration
    /** Suggest new position**/
    for (particle=0;particle<2;particle++){ // For each particle step
      for(dim=0;dim<3;dim++){
        posnew[particle][dim]=posold[particle][dim]+dr*rand();
      }
    }
    wfnew=(*system).functionCart(posnew);
    if (wfnew*wfnew/(wfold*wfold) > rand()+0.5){
      accepted_moves_last10000++;
      if (i>skip){
        accepted++;
      }
      local_energy=(*system).energyCart(posnew);
      potential=(*system).potentialCart(posnew);
      for (particle=0;particle<2;particle++){ // For each particle step
        for(dim=0;dim<3;dim++){
          posold[particle][dim]=posnew[particle][dim];
        }
      }
      wfold=wfnew;
    }
    if (i%check_size==0 && i!=0){ //each 10.000th time
      if (accepted_moves_last10000>0.6*check_size && i!=0){
        (this->dr)*=1.1;
        if(!equilibrium){
          skip+=check_size;
        }
      }
      else if (accepted_moves_last10000<0.4*check_size && i!=0){
        (this->dr)*=0.9;
        if(!equilibrium){
          skip+=check_size;
        }
      }
      else{
        equilibrium=true;
      }
      accepted_moves_last10000=0;
    }
    if (i>skip && equilibrium){
      *distance_+=distance(posold);
      *energy=*energy+local_energy;
      *energysquared=*energysquared+local_energy*local_energy;
      *v=*v+potential;
    }
  }
  *energy=*energy/(float)(amount);
  *energysquared=*energysquared/(float)(amount);
  *distance_=*distance_/(float)(amount);
  *v=*v/(float)(amount);
  cout <<"accepted moves: " << accepted << "  Energy: " << *energy<<"Alpha: "<<system->alpha<<"Beta: "<<system->beta << endl;
  skip=original_skip;
}
double VRMonteCarlo::rand(){ // Returns a random number between -0.5 and 0.5
  return prob_dist(mt_eng)-0.5;
}
