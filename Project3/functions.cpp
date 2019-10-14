#include "functions.hpp"
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

  startamount=amount_per_thread*mynum;
  endamount=amount_per_thread*(mynum+1);
	if (division_rest>=(mynum+1)){
    amount_per_thread+=1;
  }
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
  //returnval[loopings]=endamount-startamount;
	returnval[loopings]=amount_per_thread;
  return returnval;
}
double int_func(double x1, double x2, double y1, double y2, double z1, double z2){
  double exponential_part=exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2))); // Exponential part of the function
  double denom=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)); // Denominator part
  double returnval;
  //cout << denom << endl;
  if (fabs(denom)<= 1e-15){ //If denominator is zero, it is omitted
    returnval=0;
  }
  else{
    returnval= exponential_part/denom; //Final value of function
  }
  return returnval;
}
double int_func_polar_gaulag(double r1, double r2, double theta1, double theta2, double phi1, double phi2){
  double cosb=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2); // the actor
  double top=r1*r1*r2*r2*sin(theta1)*sin(theta2)*exp(-3*(r1+r2));
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
double int_func_nor1r2(double r1, double r2, double theta1, double theta2, double phi1, double phi2){
  /*
  Same as int_func_polar_gaulag, but here we removed the r1 and r2 as well as one e^r1 and e^-r2
  in order to use laguerre-polynomials.
  */
  double cosb=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2); // the actor
  double top=sin(theta1)*sin(theta2)*exp(-3.0*(r1+r2));
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
double int_func_montecarlo(double r1, double r2, double theta1, double theta2, double phi1, double phi2){
  double cosb=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2); // the actor
  double top=r1*r1*r2*r2*sin(theta1)*sin(theta2);
  double root=r1*r1+r2*r2-2*r1*r2*cosb;
  double denum;
  if (root >1.0e-8){
    denum=sqrt(root);
  }
  else{
    return 0;
  }
  return top/denum;
}
