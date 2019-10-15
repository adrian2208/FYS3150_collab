//#include "lib.h"
#include "functions.hpp"
//#include <algorithm>
#include <random>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <functional>
#include <mpi.h>
#define PI 3.14159265358979

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
  double int_mc=0.0; double variance=0.0;
  double sum_sigma=0.0; double sum_mc=0.0;
  double x[6];
  double fx=0;
  double jacobi_det=4.0*pow(PI,4)/16.0;//pow(2*PI,2)*pow(PI,2)/16.0;'

  double int_mc_local=0.0,int_mc_total;
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
  time_start=MPI_Wtime();
  val=getParallelizationCoefficients(N,my_rank,numprocs,1);
  counter=val[1];
  while(counter >0){
    counter--;
    for(int j=0;j<2;j++){
      x[j]=-0.25*log(1.0-RnG(gen));
    }
    for(int j=2;j<4;j++){
      x[j]=PI*RnG(gen); // makes x[j] a random number between 0 and PI
    }
    for(int j=4;j<6;j++){
      x[j]=2*PI*RnG(gen); // makes x[j] a random number between 0 and 2PI
    }
    fx=int_func_montecarlo(x[0],x[1],x[2],x[3],x[4],x[5]);
    int_mc_local+=fx;
    sum_sigma+=fx*fx;
  }
  local[0]=int_mc_local;
  local[1]=sum_sigma;
  MPI_Reduce(local,total,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  time_end=MPI_Wtime();
  delete [] val;delete [] local;
  total_time=time_end-time_start;
  int_mc=total[0];
  sum_sigma=total[1];
  delete [] total;
  int_mc=int_mc/N;
  sum_sigma=sum_sigma/N;
  variance=sum_sigma-int_mc*int_mc;
  if (my_rank==0){
    ofstream outfile;
    outfile.open("results/time_info.csv",ios::out | ios::app); //time-info file
    double correct_result=5*3.14159265359*3.14159265359/(16*16);
    double relative_error=fabs(int_mc*jacobi_det-correct_result)/correct_result;
    outfile<<"\ngoodMonteCarlo,no,"<<numprocs<<","<<N<<","<<total_time<<","<<int_mc*jacobi_det<<","<<relative_error<<","<<jacobi_det*sqrt(variance/((double) N ));
    cout << "Time = " <<  total_time << " on number of processors: "  << numprocs  << endl;
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << " Monte carlo result= " << setw(10) << setprecision(8) << jacobi_det*int_mc;
    cout << " Sigma= " << setw(10) << setprecision(8) << jacobi_det*sqrt(variance/((double) N )) << endl;
    outfile.close();
  }
  MPI_Finalize ();
}
