#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
using namespace std;
string createFileName(string name,double n); // creates a string in the form *name*_*n*.txt
double* createArray(int n); //Creates an Array of length n, filled with zeros
void fillArray(double* arr,int n,double num); //fills an array (of length n) (by reference) with a given number (num)
void fillArrayFunction(double * arr, int n,double *x,double (*func)(double)); //same as fillArray, but doesnt take an umber, but a generic function
double func(double x); //The function we have to double integrate
void solve(double *a, double *b, double*c, double *deriv, int n,double *solutions); //The actual solving algorithm for a general triangular matrix does changes to "solutions" (with length n+2)
/* Oppgave c start */
void improvedSolve(double *deriv,double *b, int n, double *solutions); // The solving algorithm when a = c=-1 and b=2
/*Oppgave e start */
double** createNNMatrix(int n); //Creates a double two-pointer of size n*n (that is, an n*n matrix)
void deleteNNMatrix(double** A,int n); //Deletes a matrix created by createNNMatrix
double ** createJacobiMatrix(int n, int i, int j, double theta); // Creates a Jacobi matrix
void invertJacobiMatrix(double **S,int n, int i, int j); //Inverts an existing Jacobi Matrix (that is, S is being overwritten)
double ** matrixMult(int n, double **A, double **B);
void jacobi_diag(double** A, double **R, int n, double tol);
double* findMax(double **A, int n); // returns an array where arr[0] is the largest value, arr[1] the i element and arr[2] the j element
void printMatrix(double **A, int n);
void jacobi_diag(double **A,double **R, int n, double tol);
int jacobi_diag2(double **A,double **R, int n, double tol);

double ** copyMatrix(double**A,int n);
