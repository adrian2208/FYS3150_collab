#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
//#include "System.h"
#include "VMC.h"
#include "functions.h"
#include <fstream>
#include <limits>
#include <mpi>
using namespace std;
void print(double alpha,double beta,double energy, double sigma, double distance, double dalpha, double dbeta, ofstream * outfile){
  cout << "alpha: "<< alpha <<endl<< "beta: "<< beta <<endl<< "energy: " << energy<<endl;
}
