/*Needs to be compiled as
cpp -c unit_testing.cpp
cpp -c tests_main.cpp
cpp -o unit_testing_executable unit_testing.o vecop.o tests_main.o
./unit_testing_executable
*/
/*The sampling is explained better witth comments in e.cpp*/
#include <algorithm>
#include <random>
#include <cmath>
#include "catch.hpp"
#include "vecop.hpp"
using namespace std;
void analytical(double temp, double* results);
void simple_sampling(double temp, int n, int amount, double* results);
TEST_CASE("Test that the analytical results match with the simulated results for several temperatures"){
  double * sampling_results=new double[7];
  double * analytical_results=new double[7];
  double relative_error,largest_relative_error=0,largest_temp; int largest_i;
  cout << "Temp Energy Energy^2 abs(mag) mag mag^2 cv xi"<<endl;
  for(double temp=2.0;temp<3.0;temp+=0.1){
    for(int i=0;i<7;i++){
      sampling_results[i]=analytical_results[i]=0.0;
    }
    analytical(temp,analytical_results);
    simple_sampling(temp,2,1e6,sampling_results);
    cout << temp;
    for(int i=0;i<7;i++){
      cout << "   "<< analytical_results[i]/4;
      if (fabs(analytical_results[i])>1e-5){ //Most should be much larger
        relative_error=fabs((sampling_results[i]-analytical_results[i])/analytical_results[i]);
        if(relative_error>largest_relative_error){
          largest_relative_error=relative_error;
          largest_temp=temp;
          largest_i=i;
        }
        REQUIRE(relative_error<1e-2); //Relative error should be less than one percent
      }
      else{
        REQUIRE(fabs(sampling_results[i]-analytical_results[i])<1e-1); // absolute error if analytical results is close to zero (only magnetization)
      }
    }
    cout <<endl;
  }
  cout << "Largest relative error: " << largest_relative_error << " at temperature " << largest_temp << " in position " << largest_i<<endl;
}
void analytical(double temp, double* results){ // n is assumed to be zero
  int energies[16]; int magnetics[16]; double beta=1/temp; int energy; int magnet;
  int counter=0;
  int ** A=setUpUpMatrix(2);
  for(int i=-1;i<2;i+=2){
  	for(int j=-1;j<2;j+=2){
  		for(int k=-1;k<2;k+=2){
  			for(int l=-1;l<2;l+=2){
  				A[0][0]=i;
  				A[0][1]=j;
  				A[1][0]=k;
  				A[1][1]=l;
  				energy=findStartEnergy(A,2);
          energies[counter]=energy;
  				magnetics[counter]=findStartMagnetization(A,2);
          counter++;
        }
      }
    }
  }
  double Z=0; // Partition function
  for (int i=0;i<16;i++){
    Z+=exp(-beta*energies[i]); //Update for each energy
  }
  for (int i=0;i<16;i++){
    results[0]+=exp(-beta*energies[i])*energies[i];
    results[1]+=exp(-beta*energies[i])*energies[i]*energies[i];
    results[4]+=exp(-beta*energies[i])*magnetics[i]*magnetics[i];
    results[2]+=exp(-beta*energies[i])*fabs(magnetics[i]);
    results[3]+=exp(-beta*energies[i])*magnetics[i];
  }
  for(int i=0;i<5;i++){
    results[i]/=Z; //divide everything by partition function
  }
  results[5]=(results[1]-results[0]*results[0])/(temp*temp);
  //results[6]=(results[4]-results[2]*results[2])/temp;
  results[6]=(results[4]-results[3]*results[3])/temp;
}
void simple_sampling(double temp, int n, int amount, double* results){
  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_real_distribution<double> RnG(0.0,1.0);

  double exponents[17]; int deltaE; int swap_i; int swap_j; int newSpin; int deltaM; int magnet; int energy;
  int ** A=setUpUpMatrix(2);
  magnet=findStartMagnetization(A, 2);
  energy=findStartEnergy(A,2);
  int warmUp=20000;
  for(int i=-8;i<=8;i+=4){
    exponents[i+8]=exp(-i/temp);
  }
  for(int i=0;i<amount;i++){
    for(int x=0;x<n;x++){
      for(int y=0;y<n;y++){
        swap_i=(int)(n*RnG(gen));
        swap_j=(int)(n*RnG(gen));
        deltaE=2*A[swap_i][swap_j]*(A[getPeriodic(swap_i+1,n)][swap_j]+A[getPeriodic(swap_i-1,n)][swap_j]+A[swap_i][getPeriodic(swap_j-1,n)]+A[swap_i][getPeriodic(swap_j+1,n)]);
        newSpin=-A[swap_i][swap_j];
        deltaM=2*newSpin;
        if(exponents[deltaE+8]>=RnG(gen)){
          A[swap_i][swap_j]*=-1;
          magnet+=deltaM;
          energy+=deltaE;
        }
      }
    }
    if (i>=warmUp){ // When the system is done equilbriating
      results[0]+=energy;
      results[1]+=energy*energy;
      results[4]+=magnet*magnet;
      results[2]+=fabs(magnet);
      results[3]+=magnet;
    }
  }
  for(int l=0;l<5;l++){
  results[l]/=(double)(amount-warmUp);
  }
  results[5]=(results[1]-results[0]*results[0])/(temp*temp);
  //results[6]=(results[4]-results[2]*results[2])/temp;
  results[6]=(results[4]-results[3]*results[3])/temp;
}
