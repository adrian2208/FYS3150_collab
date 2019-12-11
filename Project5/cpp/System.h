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
    double alpha;
    double beta;
    double omega;
  public:
    System(double alpha1,double beta1, double omega1);
    double functionCart(double **pos);
    double functionCart(double x1, double y1, double z1, double x2, double y2, double z2);
    double energyCart(double x1, double y1, double z1, double x2, double y2, double z2);
    double energyCart(double **pos);
    virtual double energy(double r1, double r2, double r12)=0;
    virtual double function(double r1,double r2,double r12)=0;
    virtual void update(double alpha, double beta, double omega)=0;
    double potentialCart(double x1, double y1, double z1, double x2, double y2, double z2);
    double potentialCart(double **pos);
    double potential(double r1,double r2,double r12);
};
class System1: public System{
  public:
    using System::System;
    double energy(double r1, double r2,double r12);;
    double function(double r1, double r2, double r12);
    void update(double alpha, double beta, double omega);
    double potential(double r1,double r2,double r12);
};
class System2: public System{
  public:
    using System::System;
    System1 system1=System1(alpha,beta,omega);
    //double energy(double x1, double y1, double z1, double x2, double y2, double z2);
    double energy(double r1, double r2,double r12);
    double function(double r1, double r2, double r12);
    void update(double alpha, double beta, double omega);
    double potential(double r1,double r2,double r12);
    //double energy(double ** pos);
    //double function(double ** pos);
};
#endif
class Testsystem: public System{
  using System::System;
  double energy(double r1, double r2,double r12);
  double function(double r1, double r2, double r12);
  void update(double alpha, double beta, double omega);
  double potential(double r1,double r2,double r12);
};
