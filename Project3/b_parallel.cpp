#include "lib.h"
#include "functions.hpp"
#include <algorithm>
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
  double *theta = new double[N];
  double *phi = new double[N];
  double *r=new double[N+1];
  double *omega_theta=new double[N];
  double *omega_phi=new double[N];
  double *omega_r=new double[N+1];
  double local_sum,total_sum, add_var;
  long long int * val;
  int i,j,k,l,m,n,counter;
  int local_n,numprocs,my_rank;
  double time_start,time_end,total_time;
  gauleg(0,PI,theta,omega_theta,N);
  gauleg(0,2*PI,phi,omega_phi,N);
  gauss_laguerre(r,omega_r,N,0);
  moveToLeft(r,N);
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
    int_func_polar_gaulag(r[i],r[j],theta[k],theta[l],phi[m],phi[n]);
    n++;
    counter--;
  }
  MPI_Reduce(&local_sum,&total_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  time_end=MPI_Wtime();
  total_time=time_end-time_start;
  if (my_rank==0){
    cout<<"total sum:" << total_sum << endl;
    cout << "Time = " <<  total_time << " on number of processors: "  << numprocs  << endl;
  }
  MPI_Finalize ();
}
