#include "vecop.hpp"
#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <mpi.h>
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
int main(int argc, char** argv){
  int amount; // One billion

  double t_start,t_end,dt;
  if(argc>=5){
    amount=atoi(argv[1]);
    t_start=atof(argv[2]);
    t_end=atof(argv[3]);
    dt=atof(argv[4]);
  }
  else{
    cout << "You need to state the amount of iterations, the start temperature, the end temperature and the steplength" << endl;
    exit(1);
  }
  int tot_temp=(int)(ceil((t_end-t_start)/dt)+1e-8); //Amount of temperature calculations
  int warmUp=1000000; // How many runs are "ignored" before the system is in equilibrium. One million seems about reasonabel
  int L[4]={20,20,20,20};
  double *temperatures=new double[tot_temp];

  int counter=0;double t_pos=t_start;
  while (t_pos<t_end-1e-10){
    temperatures[counter]=t_pos;
    cout << temperatures[counter] << endl;
    counter++;t_pos+=dt;
  }
  double time_start,time_end,total_time,temp;
  int numprocs,my_rank; // numprocs __needs__ to be 4, otherwise the program has to go through quite some changes...
  int magnet,energy,swap_i,swap_j,deltaE,deltaM,newSpin,accepted_configurations=0;
  double exponents[17];
  double * results=new double[5]; //Sum of energies, sum of energies squared, magnetical moment, magnetical moment squared, absolute magnetical moment
  for(int i=0;i<5;i++){
    results[i]=0;
  }
  double all_results_total[4*tot_temp][8]; // Array storing all results
  for(int i=0; i<4*tot_temp;i++){
    for(int j=0;j<4*tot_temp;j++){
      all_results_total[i][j]=0;
    }
  }
  double all_results[4*tot_temp][8]; // Array storing all results
  for(int i=0; i<4*tot_temp;i++){
    for(int j=0;j<4*tot_temp;j++){
      all_results[i][j]=0;
    }
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD,&my_rank);
  std::random_device rd;
  std::mt19937_64 gen(rd()+my_rank); //Each thread gets a different seed, as the rank is included
  std::uniform_real_distribution<double> RnG(0.0,1.0);
  int ** A=setUpRandomMatrix(L[my_rank]);
  magnet=findStartMagnetization(A, L[my_rank]);
  energy=findStartEnergy(A, L[my_rank]);
  int result_start=tot_temp*my_rank, actualpos;
  for (int tempcounter=0;tempcounter<tot_temp;tempcounter++){
    temp=temperatures[tempcounter];
    actualpos=result_start+tempcounter;
    for(int i=-8;i<=8;i+=4){
      exponents[i+8]=exp(-i/temp);
    }
    for(int i=0;i<amount;i++){
      swap_i=(int)(L[my_rank]*RnG(gen));
      swap_j=(int)(L[my_rank]*RnG(gen));
      deltaE=2*A[swap_i][swap_j]*(A[getPeriodic(swap_i+1,L[my_rank])][swap_j]+A[getPeriodic(swap_i-1,L[my_rank])][swap_j]+A[swap_i][getPeriodic(swap_j-1,L[my_rank])]+A[swap_i][getPeriodic(swap_j+1,L[my_rank])]);
      newSpin=-A[swap_i][swap_j];
      deltaM=2*newSpin;
      if(exponents[deltaE+8]>=RnG(gen)){
        A[swap_i][swap_j]*=-1;
        magnet+=deltaM;
        energy+=deltaE;
      }
      if (i>=warmUp){ // When the system is done equilbriating
        all_results[actualpos][3]+=energy;
        all_results[actualpos][4]+=energy*energy;
        all_results[actualpos][6]+=magnet*magnet;
        all_results[actualpos][7]+=fabs(magnet);
        all_results[actualpos][5]+=magnet;
      }

    }
    all_results[actualpos][0]=temp;
    all_results[actualpos][1]=L[my_rank];
    all_results[actualpos][2]=amount-warmUp;
    for(int l=3;l<8;l++){
      all_results[actualpos][l]/=(double)(amount-warmUp);
    }
  }
  for(int i=0;i<4*tot_temp;i++){
    for(int j=0;j<8;j++){
      MPI_Reduce(&all_results[i][j],&all_results_total[i][j],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    }
  }
  if (my_rank==0){
    ofstream outfile;
    outfile.open("../results/results_calculations.csv"); //time-info file
    outfile << "temperature,matrix_size,steps,energy_total,energySquared_total,magnetic_total,magneticSquared_total,magneticAbsolute_total\n";
    for(int i=0;i<4*tot_temp;i++){
      for(int j=0;j<7;j++){
        outfile << setprecision(8) << all_results_total[i][j] << ",";
      }
      outfile << setprecision(8) << all_results_total[i][7] << "\n";
    }
    outfile.close();
  }
  MPI_Finalize ();
}
