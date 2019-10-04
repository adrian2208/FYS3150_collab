#include "lib.h"
#include <mpi.h>
#include <math.h>
#include <iomanip>
#include <iostream>

using namespace std;
int g;
double int_func(double x1, double x2, double y1, double y2, double z1, double z2){
  double exponential_part=exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2))); // Exponential part of the function
  double denom=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)); // Denominator part
  double returnval;
  //cout << denom << endl;
  if (fabs(denom)<= 1e-15){ //If denominator is zero, it is omitted
    returnval=0;
    g++;
  }
  else{
    returnval= exponential_part/denom; //Final value of function
  }
  return returnval;
}
int main(int argc, char** argv){
  int N;
  N=atoi(argv[1]);
  double *x = new double[N];
  double *w=new double[N];
  double lambda=3;
  double add_var;
  gauleg(-lambda,lambda,x,w,N); // Calls the lib.cpp function gauleg // THIS  LOOKS LIKE GULAG
  int local_n,numprocs,my_rank;
  double time_start,time_end,total_time;
  double local_sum;
  double total_sum;
  int i;
  MPI_Init(&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD,&my_rank);
  time_start=MPI_Wtime();
  total_sum=0;
  i=my_rank;
    for (int j = 0;j<N;j++){
      for (int k = 0;k<N;k++){
        for (int l = 0;l<N;l++){
           for (int m = 0;m<N;m++){
              for (int n = 0;n<N;n++){
                local_sum+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_func(x[i],x[j],x[k],x[l],x[m],x[n]);
                /*This is clearly only allowed because all 6 dimensions have the same dimensionality */
                //cout <<"val of addvar " <<add_var<<endl;
     		      }
            }
        }
      }
    }
  MPI_Reduce(&local_sum,&total_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  time_end=MPI_Wtime();
  total_time=time_end-time_start;
  if (my_rank==0){
    cout<<"total sum:" << total_sum << endl;
    cout << "Time = " <<  total_time << " on number of processors: "  << numprocs  << endl;
    cout << g << endl;
  }
  MPI_Finalize ();
}
