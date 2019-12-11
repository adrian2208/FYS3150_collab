#include "System.h"
System::System(double alpha1,double beta1, double omega1){
      alpha = alpha1;
      beta = beta1;
      omega = omega1;
    }
double System::potentialCart(double ** pos){
      return potentialCart(pos[0][0],pos[0][1],pos[0][2],pos[1][0],pos[1][1],pos[1][2]);
    }
double System::potentialCart(double x1, double y1, double z1, double x2, double y2, double z2){
      return potential(sqrt(x1*x1+y1*y1+z1*z1),sqrt(x2*x2+y2*y2+z2*z2),distance(x1,y1,z1,x2,y2,z2));
    }
double System::potential(double r1, double r2,double r12){
  
}
double System::functionCart(double ** pos){
      return functionCart(pos[0][0],pos[0][1],pos[0][2],pos[1][0],pos[1][1],pos[1][2]);
    }
double System::functionCart(double x1, double y1, double z1, double x2, double y2, double z2){
      return function(sqrt(x1*x1+y1*y1+z1*z1),sqrt(x2*x2+y2*y2+z2*z2),distance(x1,y1,z1,x2,y2,z2));
    }
double System::energyCart(double x1, double y1, double z1, double x2, double y2, double z2){
      return energy(sqrt(x1*x1+y1*y1+z1*z1),sqrt(x2*x2+y2*y2+z2*z2),distance(x1,y1,z1,x2,y2,z2));
    }
double System::energyCart(double ** pos){
  return energyCart(pos[0][0],pos[0][1],pos[0][2],pos[1][0],pos[1][1],pos[1][2]);
}
double System1::energy(double r1, double r2,double r12){
      return 0.5*omega*omega*(r1*r1+r2*r2)*(1-alpha*alpha)+3*alpha*omega+1.0/r12;
    }
double System1::function(double r1, double r2, double r12){
      return exp(-alpha*omega*(r1*r1+r2*r2)*0.5);
    }
double Testsystem::energy(double r1, double r2,double r12){

      return 0.5*omega*omega*(r1*r1+r2*r2)*(1-alpha*alpha)+3*alpha*omega;
}
double Testsystem::function(double r1, double r2, double r12){
      return exp(-alpha*omega*(r1*r1+r2*r2)*0.5);
}

double System2::energy(double r1, double r2,double r12){
      double EL1=system1.energy(r1,r2,r12);
      double enbetar12=(1.0+beta*r12);
      return EL1+1.0/(2.0*enbetar12*enbetar12)*(alpha*omega*r12-1/(2.0*enbetar12*enbetar12)-2/r12+2*beta/enbetar12);
    }
double System2::function(double r1, double r2, double r12){
      return (system1.function(r1,r2,r12))*exp(r12/(1.0+beta*r12));
    }
void System1::update(double alpha, double beta, double omega){
  this->alpha=alpha;
  this->beta=beta;
  this->omega=omega;
}
void Testsystem::update(double alpha, double beta, double omega){
  this->alpha=alpha;
  this->beta=beta;
  this->omega=omega;
}
void System2::update(double alpha, double beta, double omega){
  this->alpha=alpha;
  this->beta=beta;
  this->omega=omega;
  system1.update(alpha,beta,omega);
}
