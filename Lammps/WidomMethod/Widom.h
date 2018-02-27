#define NumberOfParticles 50                                /* number of particle s*/
#define NtotIterations 4000                                 /* number of total MC iterations */
#define SQR(x) ((x)*(x))                                    /* Defines the operation x^2*/
#define CUBE(x) ((x)*(x)*(x)) 				    /* Defines the operation x^3*/


typedef struct {
  double x;
  double y;
  double z;
} vector;

/*************** FUNCTIONS ***************/

//void Initialize();                                      
double RAND();  //Generates a Random Number between [0,1]
void Initialize();  //Initializes the systems, particles start at random positions.                                              
void EnergyParticle(vector Positions,int i,double *En,double *Vir); // calculates the energy of particle i with all other particles
void EnergySystem(); // calculates total system energy

void MC_move(int i); // Attemps to displace the k-th particle.

void MC_iteration(); // Attemps to make N MC movements, where N=NumberOfParticles.
/*************** VARIABLES ***************/
vector Positions[NumberOfParticles];
double Sigma, A, Cutoff;
double Beta;
double Rho, Length;
double Eni, Viri;
double TotalEnergy, EnergyPerParticle;
double TotalVirial, Pressure;
double Etail, Ptail;
double AverageEnPerPar, AveragePressure;
double Delta;
int Naccepted, Nattempts;
double Acceptance;

