int dprimeIdx, dprime0Idx;
int num_out_constraint;

double fixed_theta, fixed_dprime;
double fixed_thetaM,fixed_thetaF; // Sex-specific analysis


double maxima_x[20];
double overallMOD= __DBL_MIN_10_EXP__ ;// replacing name : maximum_function_value = 0.0; 6/4/2009
double dprime0_MOD; //maximum_dprime0_value;
double theta0_MOD; //maximum_theta0_value;
double localmax_x[20];
double localMOD; // replacing name :localmax_value = 0.0; 6/4/2009

//double localMaxLR;

int total_dim = 0;

double initialProb2[3];
void *initialProbAddr2[3];
void *initialHetProbAddr[3];

double alpha[5][2] = { //{0.8, 1.0},  //This is for LOD not for HLOD
{0.04691, 0.118463443},
{0.230765, 0.239314335},
{0.5, 0.284444444},
{0.769235, 0.239314335},
{0.95309, 0.118463443}
};

SubLocusList savedLocusList;
SubLocusList traitLocusList;
SubLocusList markerLocusList;
int polynomialFlag;

dcuhre_state *s,init_state;
double *xl;   //xl[17] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0 };
double *xu;   //xu[17] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,1 };


/* st_DKMaxModel Moved to integrationGlobals.h 6/18/2009 */
