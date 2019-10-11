#include "lib.h"
#include "functions.hpp"
#include <algorithm>
#include <mpi.h>
#include <random>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <functional>
using namespace std;

int main(int argc, char** argv){
  int N;
  if(argc>=2){
    N=atoi(argv[1]); //n needs to be stated
  }
  else{
    cout << "You need to state a number N" << endl;
    exit(1);
  }
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> RnG(0.0,1.0);

  double lambda=2; //length of the box
  double int_mc=0.; double variance=0.;
  double sum_sigma=0.; double sum_mc;
  double x[6];
  double fx;
  double jacobi_det=pow(2*lambda,6);
  double int_mc_local,int_mc_total;
  double* local=new double[2];
  double* total=new double[2];
  long long int * val;
  int local_start,local_end,numprocs,my_rank;
  double time_start,time_end,total_time;
  double each;
  long long int counter;
  MPI_Init(&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD,&my_rank);
  each=((double)N/(double)numprocs);
  time_start=MPI_Wtime();
  val=getParallelizationCoefficients(N,my_rank,numprocs,1);
  counter=val[1];
  while(counter>0){
    counter--;
    for(int j=0;j<6;j++){
      x[j]=-lambda+2*lambda*RnG(gen); // makes x[j] a random number between -lambda and lamnbda
    }
    fx=int_func(x[0],x[1],x[2],x[3],x[4],x[5]);
    int_mc+=fx;
    sum_sigma+=fx*fx;
  }
  local[0]=int_mc;
  local[1]=sum_sigma;
  MPI_Reduce(local,total,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  //MPI_Reduce(&local[0],&total[0],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  //MPI_Reduce(&local[1],&total[1],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  time_end=MPI_Wtime();
  total_time=time_end-time_start;
  int_mc=total[0];
  sum_sigma=total[1];
  int_mc=int_mc/N;
  sum_sigma=sum_sigma/N;
  variance=sum_sigma-int_mc*int_mc;
  if (my_rank==0){
    cout << "Time = " <<  total_time << " on number of processors: "  << numprocs  << endl;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << " Monte carlo result= " << setw(10) << setprecision(8) << jacobi_det*int_mc;
    cout << " Sigma= " << setw(10) << setprecision(8) << jacobi_det*sqrt(variance/((double) N )) << endl;
  }
  MPI_Finalize ();
}
