//more to add here
#include "vecop.hpp"
#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <algorithm> //to get "count" function
#include <string>
using namespace std;
void write_energies(int total,int matrix_size,double temp,int *energies){
  ofstream energyfile;
  string name="../results/individual_energies_"+std::to_string(temp)+".txt";
  energyfile.open(name);
  energyfile << setiosflags(ios::showpoint | ios::uppercase);
  energyfile << "Steps, Matrix_size,temp"<<endl;
  energyfile << setprecision(2)<< total<<","<<matrix_size<<","<<temp<<endl;
  int count;
  for(int i=-2*matrix_size*matrix_size;i<=2*matrix_size*matrix_size;i+=4){
    //As energy varies between -2*J*n_particles and 2*J*n_particles
    energyfile<<i <<"," <<std::count(energies,energies+total,i)<<endl;
  }
  energyfile.close();
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


int main(int argc, char** argv){
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> RnG(0.0,1.0);
  ofstream single_results;
  ofstream all_results;
  single_results.open("../results/single_results.csv",ios::out | ios::app); //Collected results (Numbers)

  string matrixtype="up"; //Random or up
  string filename;
  filename = "../results/all_results_"+matrixtype+".csv";
  all_results.open(filename);
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
  int warmUp=20000; // How many runs are "ignored" before the system is in equilibrium
  bool count_energy=true;
  if(argc>4){
    all_results_write=true;
    warmUp=0; //In case I want to write the intermediate results, I do not need the warmup time.
  }
  double exponents[17];
  for(int i=-8;i<=8;i+=4){
    exponents[i+8]=exp(-i/temp);
  }
  int **A;
  if(matrixtype.compare("random") == 0){
    A=setUpRandomMatrix(n);
  }
  else{
    A=setUpUpMatrix(n);
  }
  int magnet=absolute_magnetisations[0]=findStartMagnetization(A, n);
  int energy=energies[0]=findStartEnergy(A, n);
  int swap_i,swap_j;
  int deltaE,deltaM,newSpin;
  int * energy_count=new int[amount-warmUp]; // The occurence of the individual energies
  double * results=new double[5]; //Sum of energies, sum of energies squared, magnetical moment, magnetical moment squared, absolute magnetical moment
  for(int i=0;i<5;i++){
    results[i]=0;
  }
  double result_energy=0,result_magnet=0,result_energySquared=0,result_magnetSquared=0,result_magnetAbs=0;
  for(int i=0;i<amount;i++){
    for(int x=0;x<n;x++){
      for(int y=0;y<n;y++){
        swap_i=(int)((double)n*RnG(gen));
        swap_j=(int)((double)n*RnG(gen));
        deltaE=2*A[swap_i][swap_j]*(A[getPeriodic(swap_i+1,n)][swap_j]+A[getPeriodic(swap_i-1,n)][swap_j]+A[swap_i][getPeriodic(swap_j-1,n)]+A[swap_i][getPeriodic(swap_j+1,n)]);
        newSpin=-A[swap_i][swap_j];
        if(exponents[deltaE+8]>=RnG(gen)){
          deltaM=2*newSpin;
          A[swap_i][swap_j]*=-1;
          magnet+=deltaM;
          energy+=deltaE;
          accepted_configurations++;
        }
      }
    }
    if (i>=warmUp){ // When the system is done equilbriating
      if(count_energy){
        energy_count[i-warmUp]=energy;
      }
      results[0]+=energy;
      results[1]+=energy*energy;
      results[3]+=magnet*magnet;
      results[4]+=fabs(magnet);
      results[2]+=magnet;
    }
    if(all_results_write && (i%when_dump==0) && i!=0){ //When I want to write to file, and i is a multiple of "when_dump"
      index=i/when_dump;
      i_values[index]=i;
      accepted_configurations_arr[index]=accepted_configurations;
      energies[index]=results[0]/(double)(i);
      absolute_magnetisations[index]=results[4]/(double)(i);
    }
  }
  write_energies(amount-warmUp,n,temp,energy_count);
  delete [] energy_count;
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
