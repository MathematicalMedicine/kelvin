int dprimeIdx, dprime0Idx;
int num_out_constraint;

double fixed_theta;
double fixed_dprime;
FILE *fpSeok = NULL;		// Not used anymore
FILE *fpSeok_theta = NULL;	// Not used anymore

double maxima_x[17];
double maximum_function_value = 0.0;
double localmax_x[17];
double localmax_value = 0.0;
int total_dim = 0;

int markerSetChanged;		/* flag for multipoint analysis */
int locusListChanged;		/* flag for multipoint analysis */
int prevFirstMarker;		/* first marker in the set for multipoint analysis */
int prevLastMarker;		/* last marker in the set for multipoint analysis */
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
double xl[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
double xu[15] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

int print_point_flag = 0;
FILE *fphlod = NULL;

