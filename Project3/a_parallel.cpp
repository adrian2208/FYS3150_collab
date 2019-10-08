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
long long int * getParallelizationCoefficients(long long int N,int mynum, int totNum, int loopings){
  /* calculates the
  input:
  N - the number the user gave as input
  mynum - the process' ID
  totNum - the total number of processes
  loopings: How often the loop is run

  returns: An array of amount_loopings+1 with coefficients, where arr[0] is the innermost loop_start, arr[1] the outermost loop-end and so on. Last one is amount of loops.
  */
  long long int total_amount=1;
  long long int amount_per_thread,startamount,endamount;
  long long int *returnval= new long long int[N+1];
  total_amount=pow(N,loopings)+ 1e-9;
  cout << total_amount<<endl;
  amount_per_thread=total_amount/totNum;
  int division_rest=total_amount%totNum;
  if (division_rest>=(mynum+1)){
    amount_per_thread+=1;
  }
  startamount=amount_per_thread*mynum;
  endamount=amount_per_thread*(mynum+1);

  if (division_rest>=(mynum+1)){
    startamount+=mynum;
    endamount+=mynum+1;
  }
  else{
    startamount+=division_rest;
    endamount+=division_rest;
  }
  long long int start_temp=startamount, end_temp=endamount;
  long long int startrestend,endrest;
  for(int r=0;r<loopings;r++){ // loopings amount of times
    returnval[r]=start_temp%N;
    start_temp=start_temp/N;
  }
  returnval[loopings]=endamount-startamount;
  return returnval;
}
int main(int argc, char** argv){
  int N;
  N=atoi(argv[1]);
  double *x = new double[N];
  double *w=new double[N];
  double lambda=2;
  double add_var;
  gauleg(-lambda,lambda,x,w,N); // Calls the lib.cpp function gauleg // THIS  LOOKS LIKE GULAG
  int local_n,numprocs,my_rank;
  double time_start,time_end,total_time;
  double local_sum;
  double total_sum;
  long long int * val;
  int i,j,k,l,m,n,counter;
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
    local_sum+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_func(x[i],x[j],x[k],x[l],x[m],x[n]);
    n++;
    counter--;
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
