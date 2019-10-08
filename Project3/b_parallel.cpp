#include "lib.h"
#include <algorithm>
#include <mpi.h>

#include <math.h>
#include <iomanip>
#include <iostream>
#define EPS 3.0e-14
#define MAXIT 10
#define PI 3.14159265359
using namespace std;
double gammln( double xx){
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
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
void gauss_laguerre(double *x, double *w, int n, double alf)
{
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
	}
}
// end function gaulag
double int_func(double r1, double r2, double theta1, double theta2, double phi1, double phi2){
  double cosb=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2); // the actor
  double top=r1*r1*r2*r2*sin(theta1)*sin(theta2)*exp(-4*(r1+r2));
  double denum;
  double root=r1*r1+r2*r2-2*r1*r2*cosb;
  if (root >1e-8){
    denum=sqrt(root);
  }
  else{
    return 0;
  }
  return top/denum;
}
double int_func_gaulag(double r1, double r2, double theta1, double theta2, double phi1, double phi2){
  /*
  Same as int_func, but here we removed the r1 and r2 as well as one e^r1 and e^-r2
  in order to use laguerre-polynomials.
  */
  double cosb=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2); // the actor
  double top=sin(theta1)*sin(theta2)*exp(-3*(r1+r2));
  double denum;
  double root=r1*r1+r2*r2-2*r1*r2*cosb;
  if (root >1e-8){
    denum=sqrt(root);
  }
  else{
    return 0;
  }
  return top/denum;
}
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
  gauss_laguerre(r,omega_r,N,2);
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
    int_func_gaulag(r[i],r[j],theta[k],theta[l],phi[m],phi[n]);
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
