#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "System.h"
#include "VMC.h"
using namespace std;

int main(int argc, char** argv){
  System1 per=System1(1.0,2.0,1.0);
  for (int i=0;i<=100;i++){
    cout << per.energy(i+0.1,0.0,0.0,0.0,0.0,0.1)<<endl;
  }
}
