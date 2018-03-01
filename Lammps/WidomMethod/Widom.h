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
Vec average();
Mat readIn2dData(const char* filename,int Skiprows);     /*Reads the input file*/
double LJ(double r2);
double EnergyParticle(int i);
double TotalPe();
//void PrintResults(); //Prints the Results

/*************** VARIABLES ***************/
//char Name[256]; //Input file name
//double data[rows];
//double correlation[rows];
//double results[3]; //contains the average,  autocorrelation time and the error.
double Eshift;
double Lx,Ly,Lz;
double Sigma,Epsilon;
double Rc;   //Cutoff Radius
Mat Data; //Matrix with the input data
int n,m;  //Row and column number of the input data
