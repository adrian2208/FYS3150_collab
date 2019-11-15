#include "vecop.hpp"
#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <mpi.h>
using namespace std;
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
  int tot_temp=(int)((t_end-t_start)/dt+1e-8)+1; //Amount of temperature calculations
  int warmUp=20000; // How many runs are "ignored" before the system is in equilibrium. One million seems about reasonabel
  int L[4]={40,60,80,100};
  double *temperatures=new double[tot_temp];
  int counter=0;
  while (counter<tot_temp){
    /*Fill the temperature array with the temperatures from t_start to t_end with steplength dt*/
    temperatures[counter]=t_start+dt*counter; //
    counter++;
  }
  for(int i=0;i<=tot_temp;i++){
    cout << temperatures[i] << " "<<endl;
  }
  double time_start,time_end,total_time,temp,energy_variance,magnetic_variance;
  int numprocs,my_rank; // numprocs __needs__ to be 4, otherwise the program has to go through quite some changes...
  int magnet,energy,swap_i,swap_j,deltaE,deltaM,newSpin,accepted_configurations=0;
  double exponents[17];
  double * results=new double[5]; //Sum of energies, sum of energies squared, magnetical moment, magnetical moment squared, absolute magnetical moment
  for(int i=0;i<5;i++){
    results[i]=0;
  }
  double all_results_total[4*tot_temp][12]; // Array storing all results
  for(int i=0; i<4*tot_temp;i++){
    for(int j=0;j<12;j++){
      all_results_total[i][j]=0;
    }
  }

  double all_results[4*tot_temp][12]; // Array storing all results
  for(int i=0; i<4*tot_temp;i++){
    for(int j=0;j<12;j++){
      all_results[i][j]=0;
    }
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD,&my_rank);
  time_start=MPI_Wtime();
  std::random_device rd;
  std::mt19937_64 gen(rd()+my_rank); //Each thread gets a different seed, as the rank is included
  std::uniform_real_distribution<double> RnG(0.0,1.0);


  int actualpos;
  time_start=MPI_Wtime();
  for (int l=0;l<4;l++){
    for (int tempcounter=my_rank;tempcounter<tot_temp;tempcounter+=numprocs){
      int ** A=setUpUpMatrix(L[l]);
      magnet=findStartMagnetization(A, L[l]);
      energy=findStartEnergy(A, L[l]);
      temp=temperatures[tempcounter];
      actualpos=tempcounter+l*tot_temp;
      for(int i=-8;i<=8;i+=4){
        exponents[i+8]=exp(-i/temp);
      }
      for(int i=0;i<amount;i++){
        for(int x=0;x<L[l];x++){
          for(int y=0;y<L[l];y++){
            swap_i=(int)(L[l]*RnG(gen));
            swap_j=(int)(L[l]*RnG(gen));
            deltaE=2*A[swap_i][swap_j]*(A[getPeriodic(swap_i+1,L[l])][swap_j]+A[getPeriodic(swap_i-1,L[l])][swap_j]+A[swap_i][getPeriodic(swap_j-1,L[l])]+A[swap_i][getPeriodic(swap_j+1,L[l])]);
            newSpin=-A[swap_i][swap_j];
            deltaM=2*newSpin;
            if(exponents[deltaE+8]>=RnG(gen)){
              A[swap_i][swap_j]*=-1;
              magnet+=deltaM;
              energy+=deltaE;
            }
          }
        }
        if (i>=warmUp){ // When the system is done equilbriating
          all_results[actualpos][3]+=energy;
          all_results[actualpos][4]+=energy*energy;
          all_results[actualpos][6]+=magnet*magnet;
          all_results[actualpos][7]+=fabs(magnet);
          all_results[actualpos][5]+=magnet;
        }

      }

      for(int l=3;l<8;l++){
        all_results[actualpos][l]/=(double)(amount-warmUp);
      }
      cout << "Temp: " << temp << "Size: " << L[l] << endl;
      all_results[actualpos][0]=temp;
      all_results[actualpos][1]=L[l];
      all_results[actualpos][2]=amount-warmUp;
      energy_variance = (all_results[actualpos][4]- all_results[actualpos][3]*all_results[actualpos][3])/(L[l]*L[l]);
      magnetic_variance = (all_results[actualpos][6] - all_results[actualpos][7]*all_results[actualpos][7])/(L[l]*L[l]);

      all_results[actualpos][3]/=(double)(L[l]*L[l]); //Convert energy to per particle
      all_results[actualpos][7]/=(double)(L[l]*L[l]); //Convert magnetic to per particle
      all_results[actualpos][5]/=(double)(L[l]*L[l]); //Convert magnetic_abs to per particle
      all_results[actualpos][8]=energy_variance;
      all_results[actualpos][9]=magnetic_variance;
      all_results[actualpos][10]=energy_variance/(temp*temp); // Cv
      all_results[actualpos][11]=magnetic_variance/(temp); // Xi
    }
  }
  time_end=MPI_Wtime();
  total_time=time_end-time_start;
  for(int i=0;i<4*tot_temp;i++){
    for(int j=0;j<12;j++){
      MPI_Reduce(&all_results[i][j],&all_results_total[i][j],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    }
  }
  if (my_rank==0){ // Only one thread is allowed to do this 
    writeTime(total_time,tot_temp,"yes","none");
    ofstream outfile;
    cout << "Total time: " << total_time << " seconds"<<endl;
    outfile.open("../results/results_calculations.csv"); //time-info file
    outfile << "temperature,matrix_size,steps,energy_pP,energySquared_total,magnetic_pP,magneticSquared_total,magneticAbsolute_pP,E_var,M_abs_var,Cv,Xi\n";
    for(int i=0;i<4*tot_temp;i++){
      for(int j=0;j<11;j++){
        outfile << setprecision(8) << all_results_total[i][j] << ",";
      }
      outfile << setprecision(8) << all_results_total[i][11] << "\n";
    }
    outfile.close();
  }
  MPI_Finalize ();
}
