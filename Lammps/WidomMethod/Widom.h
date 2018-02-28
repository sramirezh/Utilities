#define NumberOfParticles 75                                /* number of particle s*/
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
void Initialize(int token);  	/*Initializes the systems, particles start at random positions, 
				with random velocities but the center of mass has zero velocity. token=1, to create new initial conf
				token=0 to read initial configuration file*/
void read_init_conf();          /*Reads initial configuration file*/
void KineticEnergy();           /*Computes the (Adimensional)kinetic energy, supposing the particles have equal mass=1 */
void Force();                   // Calculate The Forces And Potential Energy (Adimensional), don't take into account all the images as in the MC'
                              

/*************** VARIABLES ***************/
vector Positions[NumberOfParticles], Velocities[NumberOfParticles],Forces[NumberOfParticles];
vector NewPositions[NumberOfParticles], NewVelocities[NumberOfParticles], OldForces[NumberOfParticles];
double Sigma, A, Cutoff,Temperature, tmax; //,,,,Maximum time of simulation
double Beta;
double Rho, Length;
double Eni, Viri;
double EPotential,EKinetic;
double Pressure;
double Etail, Ptail;
double AverageEnPerPar, AveragePressure;
double Deltat;  //time step.
int Naccepted, Nattempts, NEW_INITIAL_POSITION; //,,A control variable to recalculate the initial configuration (1) or not (0)
double Acceptance;
double TemperatureK; //Temperature Calculated from the kinetic Energy.
double SimulationTime; //Simulation Time
