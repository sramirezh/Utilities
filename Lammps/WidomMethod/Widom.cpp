#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Widom.h"


int main(int argc, char *argv[])
{
/*All the scaling are respect to sigma11,Epsilon11.
The input parameter is the trajectory file.
*/
								double mu1Source,mu2Source,mu1Sink,mu2Sink;
								FILE *file;
								file=fopen("Output","w");
								srand48(time(NULL));
								clock_t start = clock();

								//Some input parameters
								double Box []= {-30,30,-10,10,-10,10};
								Beta=1.0;


								//Basic Calculations
								Lx=Box[1]-Box[0];
								Ly=Box[3]-Box[2];
								Lz=Box[5]-Box[4];

								printf("Reading the potential parameters...\n");
								Sigma=readIn2dData("Sigma.param",0);

								Epsilon=readIn2dData("Epsilon.param",0);

								Ntypes=n;

								printf("Reading the configuration file...\n");
								Data=readIn2dData(argv[1],2);

								cout<<"The input data has "<<n<<" rows and "<<m<<" columns\n";
								cout<<"Starting the calculations \n";

								//Assuming the Cutoff radius scales with the 00 interaction parameter
								Rc=pow(2.0,1.0/6.0)*Sigma[0][0];
								std::cout << "Rc " <<Rc<< '\n';
								Eshift=EnergyShift();
								//double Pe=TotalPe();
								//Computing the chemical potential

								double Source []= {-20,-10,-10,10,-10,10};
								double Sink [] ={10,20,-10,10,-10,10};

								mu1Source=ChemPot(Source,1);
								mu2Source=ChemPot(Source,2);
								mu1Sink=ChemPot(Sink,1);
								mu2Sink=ChemPot(Sink,2);
								
								std::cout << "Printing the chemical potentials to the file Output.dat" << '\n';
								fprintf(file," %lf %lf %lf %lf \n", mu1Source, mu2Source, mu1Sink, mu2Sink);
								fclose(file);
}

Mat EnergyShift()
{
	int i,j;

	Mat matrix;
	double E;
	for (i=0;i<Ntypes;i++)
	{
		Vec row;
		for (j=0;j<Ntypes;j++)
		{
		E=LJ(SQR(Rc),i,j);
		row.push_back(E);
		}
	matrix.push_back(row);

	}

	return matrix;
}

double ChemPot(double region[6],int type)
/* Passes the region where the atoms are going to be inserted and the atom specie (type).
Returns the chemical potential for that species*/
{
								int X=10000; //Trial insertions.
								int i;
								Vec Row(4,0);
								double lx,ly,lz,U,rho,mu;

								lx=region[1]-region[0];
								ly=region[3]-region[2];
								lz=region[5]-region[4];
								Data.push_back(Row); //adding a new entrance in the array containing all the particles
								Data[n][0]=type;

								double BF=0; //Boltzmann factor
								for (i=0; i<X; i++)
								{
																Data[n][1]=(lx*RAND()+region[0]);
																Data[n][2]=(ly*RAND()+region[2]);
																Data[n][3]=(lz*RAND()+region[4]);
																U=EnergyParticle(n);
																BF+=exp(-Beta*U);

								}
								BF=BF/X;
								rho=Density(region,type);
								mu=log(rho)-log(BF);

								return mu;
}


double Density(double region[6],int type)
{
								int i;
								int Npart=0;
								double Volume,rho;
								Volume=(region[1]-region[0])*(region[3]-region[2])*(region[5]-region[4]);
								for  (i=0; i<n; i++)
								{
									if(Data[i][1]<region[1] && Data[i][1]>region[0] && Data[i][0]==type)
									{
										Npart++;
									}
								}
								rho=Npart/Volume*CUBE(Sigma[0][0]);
								return(rho);

}

double TotalPe()
{
								double Pe=0;
								int j;
								for(j=0; j<n; j++)
								{
																Pe+=EnergyParticle(j);
								}
								Pe=Pe/n*0.5;
								cout <<"Potential Energy = "<< Pe << '\n';
								return Pe;
}

double EnergyParticle(int i)
/* Computes the energy of interaction of particle i with all the other particles
   takes into accound periodic boundary conditions */
{
								double r2;
								int j;
								double dr [3];
								double Ei=0;

								for(j=0; j<n; j++)
								{
																if (j!=i)
																{
																								dr[0]=Data[i][1]-Data[j][1];
																								dr[1]=Data[i][2]-Data[j][2];
																								dr[2]=Data[i][3]-Data[j][3];

																								//Apply boundary conditions

																								dr[0]-=Lx*rint(dr[0]/Lx);
																								dr[1]-=Ly*rint(dr[1]/Ly);
																								dr[2]-=Lz*rint(dr[2]/Lz);

																								r2=SQR(dr[0])+SQR(dr[1])+SQR(dr[2]);
																								//std::cout <<"r2= "<< r2 << '\n';
																								// calculate the energy
																								if(r2<SQR(Rc))
																								{
																									int type1=int(Data[i][0])-1;
																									int type2=int(Data[j][0])-1;
																									Ei+=LJ(r2,type1,type1)-Eshift[type1][type2];
																								}
																}
								}
								return Ei;
}


double LJ(double r2, int type1, int type2)
{
								double V;
								double s2r=SQR(Sigma[type1][type2])/r2;

								V=4.0*Epsilon[type1][type2]*(pow(s2r,6)-pow(s2r,3));

								return V;
}

double RAND()
{
								return drand48();
}


vector<vector<double> > readIn2dData(const char* filename,int Skiprows)
{
/* Function takes a char* filename argument and returns a
 * 2d dynamic array containing the data
 * Then after the matrix is created:



   table.size() is the number of rows

   table[i] is row i

   table[i].size() is the number of cols. in row i

   table[i][j] is the element at the j-th col. of row i

 */


								vector< vector<double> > table;
								fstream ifs;

								//open file
								ifs.open(filename);
								int i=0; //Skiprow counter
								while (true)
								{
																string line;
																double buf;
																getline(ifs, line);
																//cout<<line<<"\n";
																i++;

																stringstream ss(line, ios_base::out|ios_base::in|ios_base::binary);

																if (!ifs) // mainly catch EOF

																								break;
																if (line[0] == '#' || line.empty()||i <= Skiprows ) // catch empty lines or comment lines and skiprows
																								continue;

																vector<double> row;

																while (ss >> buf)
																								row.push_back(buf);

																table.push_back(row);

								}

								ifs.close();
								n=table.size();
								m=table[0].size();

								return table;
}

// void PrintResults()
// {
//  //Creating the output file with the average, error and autocorrelation time.
//  FILE *ptr_file;
//  ptr_file=fopen("Results.out","w");
//  int i;
//  for (i=0;i<m;i++)
//   {
//   fprintf(ptr_file, " %lf %lf %lf \n" ,Results[0][i],Results[1][i],Results[2][i]);
//   }
//
//  fclose(ptr_file);
// }
//
