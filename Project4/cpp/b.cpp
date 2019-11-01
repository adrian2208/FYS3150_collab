//more to add here
#include "vecop.hpp"
#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
using namespace std;
int getPeriodic(int i, int n){
  return (i+n)%n;
}

int findStartEnergy(int** A, int n){
  int tot_eng=0; // Total energy
  for (int i=0;i<n-1;i++){
    /*Calculate the energy of the neighbours under for each row*/
    for (int j=0;j<n-1;j++){
      tot_eng+=A[i][j]*A[i+1][j];
    }
  }
  for (int i=0;i<n-1;i++){
    /*Calculate the energy of the neighbours to the right for each row*/
    for (int j=0;j<n-1;j++){
      tot_eng+=A[i][j]*A[i][j+1];
    }
  }
  for(int i=0;i<n-1;i++){
    /*Calculate last row and last column, but not boundary conditions*/
    tot_eng+=A[i][n-1]*A[i+1][n-1];
    tot_eng+=A[n-1][i]*A[n-1][i+1];
  }
  for (int i=0;i<n;i++){
    /*Calculate the energy due to boundary conditions*/
    tot_eng+=A[i][n-1]*A[i][0];
    tot_eng+=A[n-1][i]*A[0][i];
  }
  return -tot_eng;
}
int findStartMagnetization(int** A, int n){
  int tot_mag=0;
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      tot_mag+=A[i][j];
    }
  }
  return tot_mag;
}
/*
void random_permute(int** A,int n,int* swap_positions,double* RnG(mt19937_64),mt19937_64* gen){
  swap_positions[0]=(int) (((*RnG)(*gen)) * n);
  swap_positions[1]=(int) (((*RnG)(*gen)) * n);
  A[swap_positions[0]][swap_positions[1]]*=-1;
}
*/
int main(int argc, char** argv){
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> RnG(0.0,1.0);
  ofstream single_results;
  ofstream all_results;
  single_results.open("../results/single_results.csv",ios::out | ios::app); //Collected results (Numbers)
  all_results.open("../results/all_results.csv");
  clock_t start, finish; //Start and end time
  int n;
  double temp;
  long long int amount;
  long long int accepted_configurations=0;
  bool all_results_write=false;
  int discard=0;
  //ofstream outfile;
  if(argc>=4){
    amount=atoi(argv[1]);
    n=atoi(argv[2]);
    temp=atof(argv[3]);
  }
  else{
    cout << "You need to state the amount of iterations, the matrix dimension and the temperature" << endl;
    exit(1);
  }
  int index;
  int when_dump=1000;
  double* energies=new double[amount/when_dump];
  double* absolute_magnetisations=new double[amount/when_dump];
  long long int* i_values=new long long int[amount/when_dump];
  long long int* accepted_configurations_arr=new long long int[amount/when_dump];
  if(argc>4){
    all_results_write=true;
  }
  double exponents[17];
  for(int i=-8;i<=8;i+=4){
    exponents[i+8]=exp(-i/temp);
  }
  //int ** A=setUpRandomMatrix(n);
  int ** A=setUpUpMatrix(n);
  int magnet=findStartMagnetization(A, n);
  int energy=findStartEnergy(A, n);
  int swap_i,swap_j;
  int deltaE,deltaM,newSpin;
  double result_energy=0,result_magnet=0,result_energySquared=0,result_magnetSquared=0,result_magnetAbs=0;

  for(int i=0;i<amount;i++){
    swap_i=(int)(n*RnG(gen));
    swap_j=(int)(n*RnG(gen));
    deltaE=2*A[swap_i][swap_j]*(A[getPeriodic(swap_i+1,n)][swap_j]+A[getPeriodic(swap_i-1,n)][swap_j]+A[swap_i][getPeriodic(swap_j-1,n)]+A[swap_i][getPeriodic(swap_j+1,n)]);
    newSpin=-A[swap_i][swap_j];
    deltaM=2*newSpin;
    if(exponents[deltaE+8]>=RnG(gen)){
      A[swap_i][swap_j]*=-1;
      magnet+=deltaM;
      energy+=deltaE;
      accepted_configurations++;
    }
    result_energy+=energy;
    result_energySquared+=energy*energy;
    result_magnetSquared+=magnet*magnet;
    result_magnetAbs+=fabs(magnet);
    result_magnet+=magnet;
    if(all_results_write && (i%when_dump==0)){ //When I want to write to file, and i is a multiple of "when_dump"
      index=i/when_dump;
      i_values[index]=i;
      accepted_configurations_arr[index]=accepted_configurations;
      if(index==0){
        energies[index]=result_energy;
        absolute_magnetisations[index]=result_magnet;

      }
      else{
        energies[index]=(result_energy)/(double)(i);
        absolute_magnetisations[index]=result_magnetAbs/(double)(i);
      }
    }
  }
  deleteNNMatrix_int(A,n);
  if(all_results_write){
    all_results <<"Temperature,matrix_size,index,energies,magnetisation,accepted_configurations";
    for(int i=0;i<amount/when_dump;i++){
      all_results <<"\n"<<temp<<","<<n<<","<<  i_values[i]<<","<<energies[i]<<","<<absolute_magnetisations[i]<<","<<accepted_configurations_arr[i];
    }
  }
  delete [] energies;
  all_results.close();
  result_energy=result_energy/amount;
  result_magnet=result_magnet/amount;
  result_magnetAbs=result_magnetAbs/amount;
  result_energySquared=result_energySquared/amount;
  result_magnetSquared=result_magnetSquared/amount;
  double stdev_energy=sqrt(result_energySquared-result_energy*result_energy);
  double stdev_magnet=sqrt(result_magnetSquared-result_magnet*result_magnet);
  double stdev_magnetAbs=sqrt(result_magnetSquared-result_magnetAbs*result_magnetAbs);
  single_results <<"\n"<<temp<<","<<n<<","<<amount<<","<<result_energy<<","<<stdev_energy<<",";
  single_results <<result_magnet<<","<<stdev_magnet<<","<<result_magnetAbs<<","<<stdev_magnetAbs;
  single_results.close();
  cout <<"Absolute magnet: " <<result_magnetAbs <<"Magnet: " <<result_magnet << "Energy: " << result_energy<<endl;
  cout << "result_magnetSquared: " <<result_magnetSquared<<" result_energySquared: "<<result_energySquared<<endl;
}
