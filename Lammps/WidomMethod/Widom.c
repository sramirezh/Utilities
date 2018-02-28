#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Widom.h"


int main(void)
{
FILE *file;
file=fopen("Output","w");
Cutoff=2.5;
T=1.0
srand48(time(NULL));
X=1000 #Number of trial insertions

clock_t start = clock();



//Performs the MD
printf("Starting the run...\n");
SimulationTime=0.0; //Initializing time

int i,j;

//fprintf(file, "t	U	K	E	T	P \n");

while (SimulationTime<=tmax)
{
	KineticEnergy();
	Force();
	fprintf(file, "%lf %lf %lf %lf %lf %lf \n" ,SimulationTime,EPotential/NumberOfParticles,EKinetic/NumberOfParticles,(EPotential+EKinetic)/NumberOfParticles,TemperatureK,Pressure);

	//Position Update.
	for(i=0;i<NumberOfParticles;i++)
	{

		Positions[i].x=Positions[i].x+Velocities[i].x*Deltat+0.5*(Forces[i].x)*Deltat*Deltat;
		Positions[i].y=Positions[i].y+Velocities[i].y*Deltat+0.5*(Forces[i].y)*Deltat*Deltat;
		Positions[i].z=Positions[i].z+Velocities[i].z*Deltat+0.5*(Forces[i].z)*Deltat*Deltat;
	}

	//Reassigning the forces
	for (j=0;j<NumberOfParticles;j++)
	{
		OldForces[j].x=Forces[j].x;
		OldForces[j].y=Forces[j].y;
		OldForces[j].z=Forces[j].z;
	}
	Force();
	//Velocity Update.
	for(i=0;i<NumberOfParticles;i++)
	{
		Velocities[i].x=Velocities[i].x+0.5*(Forces[i].x+OldForces[i].x)*Deltat;
		Velocities[i].y=Velocities[i].y+0.5*(Forces[i].y+OldForces[i].y)*Deltat;
		Velocities[i].z=Velocities[i].z+0.5*(Forces[i].z+OldForces[i].z)*Deltat;
	}

//Time Update.
SimulationTime+=Deltat;

}

fclose(file);


}

double RAND()
{
	return drand48();

}

void Initialize(int token) {
  	int i;
	double scale;

	if (token==1)
	{
	FILE *ptr_file;
	ptr_file=fopen("initial_conf","w");
		printf("Creating a new initial configuration... \n");
		vector vcm;
		vcm.x=vcm.y=vcm.z=0; //Initializing center of mass velocity vector

		//Setting the initial positions and velocities, generated from a random distribution.
	  	for(i=0; i<NumberOfParticles; i++)
		{

	  		Positions[i].x = RAND()*Length;
	  	  	Positions[i].y = RAND()*Length;
	  	  	Positions[i].z = RAND()*Length;
			Velocities[i].x=RAND()-0.5;
			Velocities[i].y=RAND()-0.5;
			Velocities[i].z=RAND()-0.5;
			vcm.x=vcm.x+Velocities[i].x;
			vcm.y=vcm.y+Velocities[i].y;
			vcm.z=vcm.z+Velocities[i].z;

	  	}
  		//Calculating the velocity of the center of mass.
		vcm.x=vcm.x/NumberOfParticles;
		vcm.y=vcm.y/NumberOfParticles;
		vcm.z=vcm.z/NumberOfParticles;

		//Correction due to the movement of the center of mass.
		for (i=0; i<NumberOfParticles; i++)
		{
			Velocities[i].x=Velocities[i].x-vcm.x;
			Velocities[i].y=Velocities[i].y-vcm.y;
			Velocities[i].z=Velocities[i].z-vcm.z;
		}

		//Scaling Velocities to set initial kinetic Energy
		KineticEnergy();
		scale=sqrt((NumberOfParticles)/(EKinetic));
		double sumx,sumy,sumz;
		sumx=sumy=sumz=0;
		for(i=0;i<NumberOfParticles;i++)
		{
			Velocities[i].x*=scale;
			Velocities[i].y*=scale;
			Velocities[i].z*=scale;
			sumx+=Velocities[i].x;
			sumy+=Velocities[i].y;
			sumz+=Velocities[i].z;
		}
		//Checking if the velocity definition and scaling is correct
		KineticEnergy();
		printf("Checking if the velocity definition and scaling is correct \n");
		printf("Kinetic Energy/ Equipartition: %lf \n",EKinetic*2.0/(3.0*NumberOfParticles));
		printf("Momentum in (x,y,z):(%lf,%lf,%lf) \n",sumx,sumy,sumz);

		//Saving the initial Configuration file.
		printf("Saving initial configuration file \n");
		for(i=0;i<NumberOfParticles;i++)
		{
			fprintf(ptr_file, "%lf %lf %lf %lf %lf %lf \n" ,Positions[i].x,Positions[i].y,Positions[i].z,Velocities[i].x,Velocities[i].y,Velocities[i].z);
		}
		fclose(ptr_file);
	}
	else if (token==0)
	{
	read_init_conf();
	printf("Reading initial configuration file... \n");
	}

}

void read_conf() {
  int i,n;
  FILE*f;
  double x,y,z,vx,vy,vz;

  f = fopen("initial_conf","r");

  for(i=0; i<NumberOfParticles; i++) {
    fscanf(f, "%lf %lf %lf %lf %lf %lf", &x, &y, &z, &vx, &vy, &vz);
    Positions[i].x = x;
    Positions[i].y = y;
    Positions[i].z = z;
    Velocities[i].x = vx;
    Velocities[i].y = vy;
    Velocities[i].z = vz;
  }

  fclose(f);
}

void KineticEnergy(){
int i;
EKinetic=0;
	for (i=0; i<NumberOfParticles; i++)
	{
		EKinetic+=(SQR(Velocities[i].x)+SQR(Velocities[i].y)+SQR(Velocities[i].z));
	}
	EKinetic=0.5*EKinetic/A; //Adimensional Kinetic Energy
	TemperatureK=2.0*EKinetic/(3.0*NumberOfParticles); //Adimensional Temperature;
}


// Calculate The Forces And Potential Energy
void Force()
{
  int i,j;
  double r2,Virial;
  vector dr;

  // set forces, potential energy and pressure to zero
  for(i=0;i<NumberOfParticles;i++)
  {
    Forces[i].x=0.0;
    Forces[i].y=0.0;
    Forces[i].z=0.0;
  }

  EPotential=0.0;
  Pressure=0.0;

  // loop over all particle pairs
  for(i=0;i<NumberOfParticles-1;i++)
  {
    for(j=i+1;j<NumberOfParticles;j++)
    {
	dr.x=Positions[i].x-Positions[j].x;
	dr.y=Positions[i].y-Positions[j].y;
	dr.z=Positions[i].z-Positions[j].z;

      // apply boundary conditions
	dr.x-=Length*rint(dr.x/Length);
	dr.y-=Length*rint(dr.y/Length);
	dr.z-=Length*rint(dr.z/Length);

      r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      // check if the distance is within the cutoff radius
      if(r2<SQR(Cutoff*Length))
      {
        EPotential+=(SQR(Sigma)*exp(-pow(r2,0.5)/Sigma))/r2; //Adimensional Potential Energy
        Virial=(A*SQR(Sigma)*exp(-pow(r2,0.5)/Sigma))/r2*((pow(r2,0.5)/Sigma)+2.0);
	Pressure+=Virial;

        Forces[i].x+=Virial/r2*dr.x;
        Forces[i].y+=Virial/r2*dr.y;
        Forces[i].z+=Virial/r2*dr.z;

        Forces[j].x-=Virial/r2*dr.x;
        Forces[j].y-=Virial/r2*dr.y;
        Forces[j].z-=Virial/r2*dr.z;
      }
    }
  }
  Pressure=(1.0*CUBE(Sigma)/(3.0*CUBE(Length)))*(2.0*(EKinetic/A)-Pressure);

  //tail-correction
  Pressure+=(2.0/3.0)*SQR(Rho)*pow(Sigma,5.0)*M_PI*A*exp(-Cutoff*Length/Sigma)*(Cutoff*Length+3.0*Sigma);
  EPotential+=2.0*M_PI*NumberOfParticles*Rho*pow(Sigma,3.0)*exp(-Cutoff*Length/Sigma);

}
