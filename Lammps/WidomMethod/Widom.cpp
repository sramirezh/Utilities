#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Widom.h"


int main(int argc, char *argv[])
{
								double mu1;
								FILE *file;
								file=fopen("Output","w");
								srand48(time(NULL));
								clock_t start = clock();

//Some input parameters
								double Box []= {-30,30,-10,10,-10,10};
								Sigma=1;
								Epsilon=1;
								Beta=1;


								//Basic Calculations
								Lx=Box[1]-Box[0];
								Ly=Box[3]-Box[2];
								Lz=Box[5]-Box[4];


								Rc=pow(2.0,1.0/6.0)*Sigma;

								printf("Reading the configuration file...\n");
								Data=readIn2dData(argv[1],2);

								cout<<"The input data has "<<n<<" rows and "<<m<<" columns\n";
								cout<<"Starting the calculations \n";


								Eshift=LJ(SQR(Rc));
								//TotalPe();

								//Computing the chemical potential
								double Source []= {-20,-10,-10,10,-10,10};
								mu1=ChemPot(Source,1);
								std::cout <<  mu1<< '\n';
//while (SimulationTime<=tmax)
								fclose(file);
}


double ChemPot(double region[6],int type)
/* Passes the region where the atoms are going to be inserted and the atom specie (type).*/
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


// Row.push_back(type);
// Row.push_back(lx*RAND()+region[0]);
// Row.push_back(ly*RAND()+region[2]);
// Row.push_back(lz*RAND()+region[4]);
// Data.push_back(Row); //Creates a new entrance in the array containing all the particles
//
// U=EnergyParticle(n); //Energy due to the last particle inserted
// std::cout << U << '\n';
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
								rho=Npart/Volume*CUBE(Sigma);
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
																																Ei+=LJ(r2)-Eshift;
																								}
																}
								}
								return Ei;
}


double LJ(double r2)
{
								double V;
								double s2r=SQR(Sigma)/r2;

								V=4*Epsilon*(pow(s2r,6)-pow(s2r,3));

								return V;
}

double RAND()
{
								return drand48();
}



//
// void read_conf() {
//   int i,n;
//   FILE*f;
//   double x,y,z,vx,vy,vz;
//
//   f = fopen("initial_conf","r");
//
//   for(i=0; i<NumberOfParticles; i++) {
//     fscanf(f, "%lf %lf %lf %lf %lf %lf", &x, &y, &z, &vx, &vy, &vz);
//     Positions[i].x = x;
//     Positions[i].y = y;
//     Positions[i].z = z;
//     Velocities[i].x = vx;
//     Velocities[i].y = vy;
//     Velocities[i].z = vz;
//   }
//
//   fclose(f);
// }


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
// Vec average(){
// int i,j;
// Vec Average(m);
//
// for (j=0;j<m;j++)
// {
//  double Result=0;
//  for(i=0; i<n; i++)
//   {
//   Result=Result+InputData[i][j];
//   }
//  Average[j]=Result/n;
//
// }
//
// return(Average);
//
// }




// Calculate The Forces And Potential Energy
// void Force()
// {
//   int i,j;
//   double r2,Virial;
//   vector dr;
//
//   // set forces, potential energy and pressure to zero
//   for(i=0;i<NumberOfParticles;i++)
//   {
//     Forces[i].x=0.0;
//     Forces[i].y=0.0;
//     Forces[i].z=0.0;
//   }
//
//   EPotential=0.0;
//   Pressure=0.0;
//
//   // loop over all particle pairs
//   for(i=0;i<NumberOfParticles-1;i++)
//   {
//     for(j=i+1;j<NumberOfParticles;j++)
//     {
//  dr.x=Positions[i].x-Positions[j].x;
//  dr.y=Positions[i].y-Positions[j].y;
//  dr.z=Positions[i].z-Positions[j].z;
//
//       // apply boundary conditions
//  dr.x-=Length*rint(dr.x/Length);
//  dr.y-=Length*rint(dr.y/Length);
//  dr.z-=Length*rint(dr.z/Length);
//
//       r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
//
//       // check if the distance is within the cutoff radius
//       if(r2<SQR(Cutoff*Length))
//       {
//         EPotential+=(SQR(Sigma)*exp(-pow(r2,0.5)/Sigma))/r2; //Adimensional Potential Energy
//         Virial=(A*SQR(Sigma)*exp(-pow(r2,0.5)/Sigma))/r2*((pow(r2,0.5)/Sigma)+2.0);
//  Pressure+=Virial;
//
//         Forces[i].x+=Virial/r2*dr.x;
//         Forces[i].y+=Virial/r2*dr.y;
//         Forces[i].z+=Virial/r2*dr.z;
//
//         Forces[j].x-=Virial/r2*dr.x;
//         Forces[j].y-=Virial/r2*dr.y;
//         Forces[j].z-=Virial/r2*dr.z;
//       }
//     }
//   }
//   Pressure=(1.0*CUBE(Sigma)/(3.0*CUBE(Length)))*(2.0*(EKinetic/A)-Pressure);
//
//   //tail-correction
//   Pressure+=(2.0/3.0)*SQR(Rho)*pow(Sigma,5.0)*M_PI*A*exp(-Cutoff*Length/Sigma)*(Cutoff*Length+3.0*Sigma);
//   EPotential+=2.0*M_PI*NumberOfParticles*Rho*pow(Sigma,3.0)*exp(-Cutoff*Length/Sigma);
//
// }
