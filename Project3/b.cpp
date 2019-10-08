#include "lib.h"
#include <algorithm>
//#include <mpi.h>

#include <math.h>
#include <iomanip>
#include <iostream>
#define EPS 3.0e-14
#define MAXIT 10
#define PI 3.14159265359
using namespace std;
void gauss_laguerre(double *x, double *w, int n, double alf);
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
  if (root >1.0e-8){
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
  double top=r1*r1*r2*r2*sin(theta1)*sin(theta2)*exp(-3.0*(r1+r2));
  double denum;
  double root=r1*r1+r2*r2-2*r1*r2*cosb;
  if (root >1.0e-8){
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
  double int_gauss=0.0, add_var;
  gauleg(0,PI,theta,omega_theta,N);

  gauleg(0,2.0*PI,phi,omega_phi,N);
  gauss_laguerre(r,omega_r,N,0);
  for (int i=0;i<=N+1;i++){
    cout << omega_r[i]<< " ";
  }
  cout << endl;
  moveToLeft(r,N);
  moveToLeft(omega_r,N);
  for (int i=0;i<N;i++){
    for (int j = 0;j<N;j++){
      for (int k = 0;k<N;k++){
        for (int l = 0;l<N;l++){
           for (int m = 0;m<N;m++){
              for (int n = 0;n<N;n++){
                add_var=omega_r[i]*omega_r[j]*omega_theta[k]*omega_theta[l]*omega_phi[m]*omega_phi[n]*
                int_func_gaulag(r[i],r[j],theta[k],theta[l],phi[m],phi[n]);
                /*This is clearly only allowed because all 6 dimensions have the same dimensionality */
                //cout <<"val of addvar " <<add_var<<endl;
                int_gauss+=add_var;
     		      }
            }
        }
      }
    }
  }
  cout << "estimate: " << int_gauss << " correct: " << PI*PI*5/(16*16) << endl;
}
