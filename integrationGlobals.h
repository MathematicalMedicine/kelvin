int dprimeIdx, dprime0Idx;
int num_out_constraint;

double fixed_theta, fixed_dprime;
double fixed_thetaM,fixed_thetaF; // Sex-specific analysis


double maxima_x[20];
double maximum_function_value = 0.0;
double localmax_x[20];
double localmax_value = 0.0;
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
double xl[17] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0 };
double xu[17] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,1 };

int print_point_flag = 0;
FILE *fphlod = NULL;

