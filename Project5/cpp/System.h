#ifndef SYSTEM_H
#define SYSTEM_H
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "functions.h"
using namespace std;
class System{
  protected:
    /*parameters*/
    double alpha;
    double beta;
    double omega;
    friend class VRMonteCarlo; //VRMonteCarlo has full access
  public:
    System(double alpha1,double beta1, double omega1);
    double functionCart(double **pos); //Function in Cartesian coordinates (input matrix)
    double functionCart(double x1, double y1, double z1, double x2, double y2, double z2);  //Function in Cartesian coordinates
    double energyCart(double x1, double y1, double z1, double x2, double y2, double z2);  //energy in Cartesian coordinates (matrix)
    double energyCart(double **pos); //energy in Cartesian coordinates
    virtual double energy(double r1, double r2, double r12)=0;  //energy in polar coordinates
    virtual double function(double r1,double r2,double r12)=0; //function in Cartesian coordinates
    virtual void update(double alpha, double beta, double omega)=0;  //update alpha, beta and omega
    double potentialCart(double x1, double y1, double z1, double x2, double y2, double z2); //potential in Cartesian coordinates
    double potentialCart(double **pos);  //potential in Cartesian coordinates (matrix)
    virtual double potential(double r1,double r2,double r12)=0;  //potential in polar coordinatews
};
class System1: public System{
  public:
    using System::System; //Same constructor
    double energy(double r1, double r2,double r12);;
    double function(double r1, double r2, double r12);
    void update(double alpha, double beta, double omega);
    double potential(double r1,double r2,double r12);
};
class System2: public System{
  public:
    using System::System;
    System1 system1=System1(alpha,beta,omega); //Has it's own system1
    double energy(double r1, double r2,double r12);
    double function(double r1, double r2, double r12);
    void update(double alpha, double beta, double omega);
    double potential(double r1,double r2,double r12);
};
class Testsystem: public System{
  using System::System;
  double energy(double r1, double r2,double r12);
  double function(double r1, double r2, double r12);
  void update(double alpha, double beta, double omega);
  double potential(double r1,double r2,double r12);
};
#endif
