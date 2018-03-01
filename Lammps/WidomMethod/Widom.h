#include  <iostream>
#include  <fstream>
#include  <vector>
#include  <string>
#include <iterator>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include  <vector>

#define SQR(x) ((x)*(x))                  /* Defines the operation x^2*/
#define CUBE(x) ((x)*(x)*(x)) 				    /* Defines the operation x^3*/

using namespace std;

typedef vector<double> Vec;
typedef vector<Vec> Mat;
/*************** FUNCTIONS ***************/
Mat readIn2dData(const char* filename,int Skiprows);     /*Reads the input file*/
Mat EnergyShift();
double LJ(double r2, int i, int j);
double RAND();
double EnergyParticle(int i);
double TotalPe();
double ChemPot(double region[6], int type);
double Density(double region[6], int type);

//void PrintResults(); //Prints the Results

/*************** VARIABLES ***************/
//char Name[256]; //Input file name
//double data[rows];
//double correlation[rows];
//double results[3]; //contains the average,  autocorrelation time and the error.
double Lx,Ly,Lz;
double Beta;
double Rc;   //Cutoff Radius
int Ntypes; //Number of species
Mat Data; //Matrix with the input data
Mat Sigma, Epsilon, Eshift; //Matrices with the LJ parameters.
int n,m;  //Row and column number of the input data
