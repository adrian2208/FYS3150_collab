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
void append_to_file(int total,int matrix_size,double temp, double *results, ofstream* ofile){ // Writes relevant data (per spin particle) to file
  int num_spins=matrix_size*matrix_size;
  double result_energy=results[0]/total;
  double result_magnet=results[2]/total;
  double result_magnetAbs=results[4]/total;
  double result_energySquared=results[1]/total;
  double result_magnetSquared=results[3]/total;
  double energy_variance = (result_energySquared- result_energy*result_energy)/num_spins;
  double Mvariance = (result_magnetSquared - result_magnet*result_magnet)/num_spins;
  double xi=Mvariance/temp;
  double cv=energy_variance/(temp*temp);
  *ofile << setiosflags(ios::showpoint | ios::uppercase);
  *ofile << setprecision(8) << total<<",";
  *ofile << setprecision(8) << matrix_size<<",";
  *ofile << setprecision(8) << temp<<",";
  *ofile << setprecision(8) << result_energy/num_spins<<",";
  *ofile << setprecision(8) << cv<<",";
  *ofile << setprecision(8) << result_magnet/num_spins<<",";
  *ofile << setprecision(8) << xi<<",";
  *ofile << setprecision(8) << result_magnetAbs/num_spins << endl;
  (*ofile).close();
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
  int warmUp=1000000; // How many runs are "ignored" before the system is in equilibrium
  if(argc>4){
    all_results_write=true;
    warmUp=0; //In case I want to write the intermediate results, I do not need the warmup time.
  }
  double exponents[17];
  for(int i=-8;i<=8;i+=4){
    exponents[i+8]=exp(-i/temp);
  }
  int ** A=setUpRandomMatrix(n);
  //int ** A=setUpUpMatrix(n);
  int magnet=findStartMagnetization(A, n);
  int energy=findStartEnergy(A, n);
  int swap_i,swap_j;
  int deltaE,deltaM,newSpin;
  double * results=new double[5]; //Sum of energies, sum of energies squared, magnetical moment, magnetical moment squared, absolute magnetical moment
  for(int i=0;i<5;i++){
    results[i]=0;
  }
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
    if (i>=warmUp){ // When the system is done equilbriating
      results[0]+=energy;
      results[1]+=energy*energy;
      results[3]+=magnet*magnet;
      results[4]+=fabs(magnet);
      results[2]+=magnet;
    }
    if(all_results_write && (i%when_dump==0)){ //When I want to write to file, and i is a multiple of "when_dump"
      index=i/when_dump;
      i_values[index]=i;
      accepted_configurations_arr[index]=accepted_configurations;
      if(index==0){
        energies[index]=results[0];
        absolute_magnetisations[index]=results[2];

      }
      else{
        energies[index]=results[0]/(double)(i);
        absolute_magnetisations[index]=results[4]/(double)(i);
      }
    }
  }
  deleteNNMatrix_int(A,n); // Free matrix space
  if(all_results_write){ //If everything is to be written to file
    all_results <<"Temperature,matrix_size,index,energies,magnetisation,accepted_configurations";
    for(int i=0;i<amount/when_dump;i++){
      all_results <<"\n"<<temp<<","<<n<<","<<  i_values[i]<<","<<energies[i]<<","<<absolute_magnetisations[i]<<","<<accepted_configurations_arr[i];
    }
  }
  delete [] energies;
  all_results.close();
  append_to_file(amount-warmUp,n,temp,results,&single_results);
}
