#ifndef SYSTEM_H
#define SYSTEM_H
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
class System{
  protected:
    double alpha;
    double beta;
    double omega;
  public:
    System(double alpha,double beta, double omega);
    double function(double x1, double y1, double z1, double x2, double y2, double z2);
    double energy(double x1, double y1, double z1, double x2, double y2, double z2);
    virtual double energy(double r1, double r2, double r12)=0;
    virtual double function(double r1,double r2,double r12)=0;
};
class System1: public System{
  public:
    using System::System;
    double energy(double x1, double y1, double z1, double x2, double y2, double z2);
    double energy(double r1, double r2,double r12);
    double function(double r1, double r2, double r12);
};
class System2: public System{
  public:
    using System::System;
    System1 system1=System1(alpha,beta,omega);
    double energy(double x1, double y1, double z1, double x2, double y2, double z2);
    double energy(double r1, double r2,double r12);
    double function(double r1, double r2, double r12);
};
#endif
