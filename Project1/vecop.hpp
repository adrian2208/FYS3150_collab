#include <iostream>
#include <cmath>
#include <string>
using namespace std;
string createFileName(string name,int n); // creates a string in the form *name*_*n*.txt
double* createArray(int n); //Creates an Array of length n, filled with zeros
void fillArray(double* arr,int n,double num); //fills an array (of length n) (by reference) with a given number (num)
void fillArrayFunction(double * arr, int n,double *x,double (*func)(double)); //same as fillArray, but doesnt take an umber, but a generic function
double func(double x); //The function we have to double integrate
double* solve(double *a, double *b, double*c, double *deriv, int n); //The actual solving algorithm for a general triangular matrix returns an array of length n+2 with solutions
/* Oppgave c start */
double* improvedSolve(double *deriv,int n); // The solving algorithm when a = c=-1 and b=2
/*Oppgave e start */
double** createNNMatrix(int n); //Creates a double two-pointer of size n*n (that is, an n*n matrix)
void deleteNNMatrix(double** A,int n); //Deletes a matrix created by createNNMatrix
