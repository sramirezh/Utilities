#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "MC_NVT2.h"


int main(void)
{

FILE *ptr_file;
ptr_file=fopen("Data","w");
Rho=0.5;
Delta=1.0/1.0;
A=1.0;
Cutoff=2.0;
Beta=1.0;
Sigma=1.0;
srand48(time(NULL));

clock_t start = clock(); 

	
Length=pow(NumberOfParticles/Rho,1.0/3.0)*Sigma;
Initialize();
EnergySystem();
fprintf(ptr_file, "0 %lf %lf %lf \n" ,TotalEnergy/NumberOfParticles, Pressure,(double)Naccepted/Nattempts);
//Performs the MC iterations
int i;
for (i=0;i<NtotIterations;i++){
	MC_iteration();
	EnergySystem();
	fprintf(ptr_file, "%d %lf %lf %lf \n" ,i+1,TotalEnergy/NumberOfParticles,Pressure,(double)Naccepted/Nattempts);
}

printf("Elapsed Time : %f \n", ((double)clock() - start) / CLOCKS_PER_SEC);  

fclose(ptr_file);

}


double RAND()
{	
	return drand48();

}

void Initialize() {
  	int i;
  	for(i=0; i<NumberOfParticles; i++)
	{
  		Positions[i].x = RAND()*Length;
  	  	Positions[i].y = RAND()*Length;
  	  	Positions[i].z = RAND()*Length;
  	}
}

void EnergyParticle(vector pos,int i,double *En,double *Vir)
{
double r2,Virij,Enij;
vector dr;
int j,k,l,m;
Enij=Virij=0.0;
for(j=0;j<NumberOfParticles;j++)
{
//if(i!=j){
	dr.x=pos.x-Positions[j].x;
	dr.y=pos.y-Positions[j].y;
	dr.z=pos.z-Positions[j].z;

	
	for (k=-1;k<2;k++){
		for (l=-1;l<2;l++){
			for (m=-1;m<2;m++){
				if( j==i && k==0 && l==0 && m==0){}
				else{
									
					// apply periodic boundary conditions
					dr.x-=Length*rint(dr.x/Length)+k*Length;
					dr.y-=Length*rint(dr.y/Length)+l*Length;
					dr.z-=Length*rint(dr.z/Length)+m*Length;
					r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
					// calculate the energy
					if(r2<SQR(Cutoff*Length))
					{
						Enij+=(A*SQR(Sigma)*exp(-pow(r2,0.5)/Sigma))/r2;
						Virij+=(A*SQR(Sigma)*exp(-pow(r2,0.5)/Sigma))/r2*((pow(r2,0.5)/Sigma)+2.0);
					}
				}
				
			}
		}
	}

}
*En=Enij;
*Vir=Virij;
//}
}

// calculates total system energy
void EnergySystem(void)
{
  double Eni,Viri;
  int i;

  TotalEnergy=0.0;
  TotalVirial=0.0;
  for(i=0;i<NumberOfParticles;i++)  //It doesnt calculate the last because it has already been taken into account with all the others.
  {
    EnergyParticle(Positions[i],i,&Eni,&Viri);

    TotalEnergy+=Eni;
    TotalVirial+=Viri;
  }
   //Avoiding Overcounting
   TotalVirial=TotalVirial*0.5;
   TotalEnergy=TotalEnergy*0.5;
   Pressure=Rho*CUBE(Sigma)+(1.0/(3.0*CUBE(Length)))*Beta*TotalVirial; 
  
  //tail-correction
  Pressure+=(2.0/3.0)*SQR(Rho*Sigma)*M_PI*A*exp(-Cutoff*Length/Sigma)*(Cutoff*Length+3.0*Sigma);	
  TotalEnergy+=2.0*M_PI*NumberOfParticles*Rho*A*pow(Sigma,3.0)*exp(-Cutoff*Length/Sigma);
}

void MC_move(int i){
{
  double EnergyNew,VirialNew,EnergyOld,VirialOld;
  vector NewPosition;
  Nattempts+=1.0;
  // calculate old energy
  EnergyParticle(Positions[i],i,&EnergyOld,&VirialOld);

  // give a random displacement
  NewPosition.x=Positions[i].x+(RAND()-0.5)*Delta*Length;
  NewPosition.y=Positions[i].y+(RAND()-0.5)*Delta*Length;
  NewPosition.z=Positions[i].z+(RAND()-0.5)*Delta*Length;

  // calculate new energy
  EnergyParticle(NewPosition,i,&EnergyNew,&VirialNew);
         
  if(RAND()<exp(-Beta*(EnergyNew-EnergyOld)))
  {
   Naccepted+=1.0;
   // update new position if movement is accepted
   Positions[i].x=NewPosition.x;
   Positions[i].y=NewPosition.y;
   Positions[i].z=NewPosition.z;
  }
}

}

void MC_iteration()
{
Naccepted=0.0;
Nattempts=0.0;
int i;
for (i=0;i<NumberOfParticles;i++)
{
	MC_move(i);
}
}


