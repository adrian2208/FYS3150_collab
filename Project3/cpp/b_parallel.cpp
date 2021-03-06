#include "lib.h"
#include "functions.hpp"
#include <mpi.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#define PI 3.14159265359
using namespace std;

void moveToLeft(double *arr,int N){
  for (int i=0;i<N;i++){
    arr[i]=arr[i+1];
  }
}
int main(int argc, char** argv){
  int N;
  if(argc>=2){
    N=atoi(argv[1]); //n needs to be stated
  }
  else{
    cout << "You need to state a number N" << endl;
    exit(1);
  }
  /*The weights and meshpoints for the angles and r1 and r2*/
  double *theta = new double[N];
  double *phi = new double[N];
  double *r=new double[N+1];
  double *omega_theta=new double[N];
  double *omega_phi=new double[N];
  double *omega_r=new double[N+1];
  double local_sum=0,total_sum=0;
  long long int * val;
  long long int i,j,k,l,m,n,counter;//Future "entry points" for each thread. In theory, i-n don't need to be lli
  int numprocs,my_rank; //Necessary for parallelisation
  double time_start,time_end,total_time;
  gauleg(0,PI,theta,omega_theta,N);
  gauleg(0,2*PI,phi,omega_phi,N);
  gauss_laguerre(r,omega_r,N,0);
  moveToLeft(r,N); //So that the array goes from 0 to N, not 1 to N+1
  moveToLeft(omega_r,N);
  MPI_Init(&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD,&my_rank);
  time_start=MPI_Wtime();
  val=getParallelizationCoefficients(N,my_rank,numprocs,6);
  n=val[0]; m=val[1];l=val[2];k=val[3];j=val[4];i=val[5];counter=val[6];
  total_sum=0;
  while(counter>0){
    if(n>=N){
      n=0; m++;
      if(m>=N){
        m=0;l++;
      }
      if(l>=N){
        l=0;k++;
      }
      if(k>=N){
        k=0;j++;
      }
      if(j>=N){
        j=0;i++;
      }
    }
    local_sum+=omega_r[i]*omega_r[j]*omega_theta[k]*omega_theta[l]*omega_phi[m]*omega_phi[n]*
    int_func_polar_gaulag(r[i],r[j],theta[k],theta[l],phi[m],phi[n]);//actual function evaluation
    n++;
    counter--;
  }
  MPI_Reduce(&local_sum,&total_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); //Give result to 1. thread
  time_end=MPI_Wtime();
  delete [] r; delete [] omega_r;delete [] phi; delete [] omega_phi;delete [] theta; delete [] omega_theta;delete [] val;
  total_time=time_end-time_start;
  if (my_rank==0){
    ofstream outfile; //Only necessary for one thread
    outfile.open("../results/time_info.csv",ios::out | ios::app); //time-info file
    double correct_result=5*3.14159265359*3.14159265359/(16*16);
    double relative_error=fabs(total_sum-correct_result)/correct_result;
    outfile<<"\ngaulag,no,"<<numprocs<<","<<N<<","<<total_time<<","<<total_sum<<","<<relative_error<<","<<0;
    outfile.close();
  }
  MPI_Finalize ();
}
