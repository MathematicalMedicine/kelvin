
/**********************************************************************
 * Kelvin - Linkage and Linkage Disequalibrium Analysis Program
 * Yungui Huang
 * Integration - Sang-Cheol Seok
 * Polynomial features - Hongling Wang
 * config.c and error logging modules - Alberto Maria Segre
 * Regex code - Nathan Burnette
 * RADSMM storage code - Martin Milder
 * 
 * Copyright 2007, Columbus Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/
#include "dkelvin.h"
#include "likelihood.h"
#include "pedlib/polynomial.h"
#include "saveResults.h"
#include "dcuhre.h"

extern Polynomial *constant1Poly;
extern char *likelihoodVersion, *locusVersion, *polynomialVersion;
struct swStopwatch *overallSW;
time_t startTime;
int currentVMK, maximumVMK;
char messageBuffer[MAXSWMSG];

#include <signal.h>		/* Signalled dumps, exits, whatever */
volatile sig_atomic_t signalSeen = 0;	/* If we wanted to be really gentle. */

/* Typing ^\ or issuing a kill -s SIGQUIT gets a dump of statistics.
   We used to set the signalSeen flag and watch for it in breaks in the
   code, but so long as we don't interfere with kelvin when we dump out
   statistics, that is unnecessary.

   P.S. - cygwin requires "stty quit ^C" first for this to work. */
void
quitSignalHandler (int signal)
{
  if (modelOptions.polynomial == TRUE)
    polyDynamicStatistics ("Signal received");
  else
    swDump (overallSW);
#ifdef DMTRACK
  swLogPeaks ("Signal");
#endif
}

#if defined (GPROF) || (GCOV)

/* We catch a SIGTERM to allow early exit() for profiling. */
void
termSignalHandler (int signal)
{
  fprintf (stderr, "Terminating early for gprof or gcov!\n");
  exit (EXIT_SUCCESS);
}
#endif

void
intSignalHandler (int signal)
{
  fprintf (stderr, "Terminating early via interrupt!\n");
  exit (EXIT_FAILURE);
}

pid_t childPID = 0;		/* For a child process producing timed statistics. */
/* Exit handler to clean-up after we hit any of our widely-distributed exit points. */
void
exit_kelvin ()
{
  swLogMsg ("Exiting");
  if (childPID != 0)
    kill (childPID, SIGKILL);	/* Sweep away any errant children */
}

char *programVersion = "V0.34.2";
char *dkelvinVersion = "$Id$";

void print_dryrun_stat (PedigreeSet * pSet, double pos);
void logStatistics(PedigreeSet *pSet, int posIdx);
void test_darray (double **);

#define checkpt() fprintf(stderr,"Checkpoint at line %d of file \"%s\"\n",__LINE__,__FILE__)

void compute_hlod_2p_dt (double x[], double *f);
void compute_hlod_mp_dt (double x[], double *f);
void compute_hlod_2p_qt (double x[], double *f);
void compute_hlod_mp_qt (double x[], double *f);

int kelvin_dcuhre_integrate (double *integral, double *abserr);

/* Variables became global from local */
PedigreeSet pedigreeSet;	/* Pedigrees. */
int loc1, loc2;
Locus *pLocus;
Locus *pLocus1;
Locus *pLocus2;
Trait *pTrait;
int traitLocus;
int totalLoci;
void *initialProbAddr[3];
char *tmpID;
LDLoci *pLDLoci = NULL;

int R_square_flag = FALSE;
double R_square = 0;
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

double alpha[5][2] = { {0.04691, 0.118463443},
{0.230765, 0.239314335},
{0.5, 0.284444444},
{0.769235, 0.239314335},
{0.95309, 0.118463443}
};

SubLocusList savedLocusList;
SubLocusList traitLocusList;
SubLocusList markerLocusList;
int polynomialFlag;

dcuhre_state *s;
double xl[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
double xu[15] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

int print_point_flag = 0;

/* Some default global values. */
char resultsprefix[KMAXFILENAMELEN + 1] = "./";
char markerfile[KMAXFILENAMELEN + 1] = "markers.dat";
char mapfile[KMAXFILENAMELEN + 1] = "mapfile.dat";
char pedfile[KMAXFILENAMELEN + 1] = "pedfile.dat";
char datafile[KMAXFILENAMELEN + 1] = "datafile.dat";
char loopsfile[KMAXFILENAMELEN + 1] = "loops.dat";
char ccfile[KMAXFILENAMELEN + 1] = "";	/* case control count file */
char outfile[KMAXFILENAMELEN + 1] = "lods.out";
char avghetfile[KMAXFILENAMELEN + 1] = "avghet.out";
char avghomofile[KMAXFILENAMELEN + 1] = "avghomo.out";
char pplfile[KMAXFILENAMELEN + 1] = "ppl.out";
char ldPPLfile[KMAXFILENAMELEN + 1] = "ldppl.out";
FILE *fpHet = NULL;		/* average HET LR file */

  //FILE *fpHomo = NULL;          /* average HOMO LR file */
FILE *fpPPL = NULL;		/* PPL output file */
int polynomialScale = 1;	/* Scale of static allocation and dynamic
				   growth in polynomial.c, 1-10 with 1 as
				   the default, and 10 the old standard. */
FILE *fphlod = NULL;

/* Model datastructures. modelOptions is defined in the pedigree library. */
ModelType modelType;
ModelRange modelRange;
ModelOptions modelOptions;

/* number of D primes 
 * if there are more than 2 alleles in the marker/trait, number of D primes
 * and D prime ranges are assumed to be the same to reduce complexity 
 * for initial phase of this project */
int num_of_d_prime;
double *d_prime;
int num_of_theta;

/* three dimensional array for the two point summary results *
 * first dimension is the D prime, for LE, D prime=0 with just one element
 * in this dimension 
 * second dimension is theta values 
 * third dimension is marker allele frequency, for LE, only one element in this dimension */
SUMMARY_STAT ***tp_result;	/*Leave these not to make errors */

/* two dimensional array per (dprime, theta) 
 * this will be used to calculate PPL
 */

/* storage for the NULL likelihood for the multipoint calculation under polynomial */
double ***likelihoodDT = NULL;
double *****likelihoodQT = NULL;	/*Leave these not to make errors */

//double **likelihoodDT = NULL;
//double **likelihoodQT = NULL;
double markerSetLikelihood;

/* for multipoint, we use genetic map positions on a chromosome */
double *map_position;
int num_of_map_position;

/* one dimensional array, indexing by map position 
 * for multipoint, we don't know how to incorporate LD in yet 
 * This map could be sex specific map or sex averaged map 
 * For two point, we don't have to distinguish sex specific/avearge 
 * as we use theta relative to marker during analysis and after analysis
 * (result) */
SUMMARY_STAT *mp_result;	/*Leave these not to make errors */
int numPositions;

XMission *nullMatrix;
XMission *altMatrix;
XMission *traitMatrix;
XMission *markerMatrix;

LambdaCell *pLambdaCell = NULL;
int prevNumDPrime = 0;
int loopMarkerFreqFlag = 0;
int total_count;

char *flexBuffer = NULL;
int flexBufferSize = 0;

/**********************************************************************
 * Usage:
 *    kelvin [-s][-c] config.dat
 *
 * The config.dat file gives information about the specific linkage
 * analysis run. All information about, e.g., which markers to use,
 * what outputs to calculate, and so on, are stored in this
 * configuration file.
 **********************************************************************/
int
main (int argc, char *argv[])
{
  int i, j, k;			/*index variables */
  char configfile[KMAXFILENAMELEN] = "";

  //int breakFlag = FALSE;

  // PedigreeSet pedigreeSet;    /* Pedigrees. */

  /* Local variable declaration */
  Pedigree *pPedigree;

  /*index variables */
  //double theta[2];            /* theta */
  //int paramIdx = -1;
  int liabIdx;
  int pedIdx;
  int posIdx;
  int prevTraitInd;
  double *prevPos, *currPos;	/* for MP */
  int locus;

  //int thresholdIdx = -1;
  //double threshold = 0;  
  double traitPos;		/* trait position for multipoint analysis */
  TraitLocus *pTraitLocus;
  int mkrFreqIdx;
  double mkrFreq;
  double *marker1Pos, *marker2Pos;
  double relativePos;
  int traitIndex = 0;
  int leftMarker = -1;

  /* PPL output */
  double ppl;
  double ldppl, ppld;


  /* other control variables */
  int status;			/* xmatrix status */
  double dist;			/*position distance between Marker 1 and 2 */

  /*Time vaiables */
  clock_t time0, time1, time2;

  //int numberOfCompute = 0;

  /*Polynomial variables, there are some more as global variables : Eventually put all of them in the same place! */
  Polynomial *initialProbPoly[3];
  Polynomial *initialProbPoly2[3];
  double initialProb[3];
  int exitDueToLoop = FALSE;
#ifdef _OPENMP
  char *envVar;
  int threadCount = 0;
#endif

  /* Variables for DCUHRE   added 1/2008 */
  double integral = 0.0, abserr = 0;
  int num_eval = 0;

  //int method_option=2;   * used to choose the integratio method
  // 1: kelvin_dcuhre with grid from the user
  // 2: kelvin_dcuhre with rule 13
  // */
  double low_theta_integral = 0.0, high_theta_integral = 0.0;
  double low_integral = 0.0, high_integral = 0.0;
  double low_ld_integral = 0.0;
  double volume_region = 1.0;

  /*Dcuhre rule points and weights for D' and theta or only theta */
  double dcuhre2[140][4] = { {0.0, 0.00234550385, 0.118463442, 0.0},
  {0.0, 0.0115382672, 0.239314336, 0.0},
  {0.0, 0.025, 0.284444444, 0.0},
  {0.0, 0.0384617328, 0.239314336, 0.0},
  {0.0, 0.0476544961, 0.118463442, 0.0},
  {0.0, 0.0711095347, 0.118463442, 0.0},
  {0.0, 0.153844405, 0.239314336, 0.0},
  {0.0, 0.275, 0.284444444, 0.0},
  {0.0, 0.396155595, 0.239314336, 0.0},
  {0.0, 0.478890465, 0.118463442, 0.0},
  {0.0000000000000, 0.0250000000000, 0.0084492309003, 0.0},
  {-0.2517129343453, 0.0250000000000, 0.0237714740190, 0.0},
  {0.2517129343453, 0.0250000000000, 0.0237714740190, 0.0},
  {-0.7013933644534, 0.0250000000000, 0.0294001617014, 0.0},
  {0.7013933644534, 0.0250000000000, 0.0294001617014, 0.0},
  {0.0000000000000, 0.0187071766414, 0.0237714740190, 0.0},
  {0.0000000000000, 0.0312928233586, 0.0237714740190, 0.0},
  {0.0000000000000, 0.0074651658887, 0.0294001617014, 0.0},
  {0.0000000000000, 0.0425348341113, 0.0294001617014, 0.0},
  {0.9590960631620, 0.0250000000000, 0.0066444364658, 0.0},
  {-0.9590960631620, 0.0250000000000, 0.0066444364658, 0.0},
  {0.0000000000000, 0.0489774015790, 0.0066444364658, 0.0},
  {0.0000000000000, 0.0010225984210, 0.0066444364658, 0.0},
  {0.9956010478552, 0.0250000000000, 0.0042536044255, 0.0},
  {-0.9956010478552, 0.0250000000000, 0.0042536044255, 0.0},
  {0.0000000000000, 0.0498900261964, 0.0042536044255, 0.0},
  {0.0000000000000, 0.0001099738036, 0.0042536044255, 0.0},
  {0.5000000000000, 0.0250000000000, 0.0000000000000, 0.0},
  {-0.5000000000000, 0.0250000000000, 0.0000000000000, 0.0},
  {0.0000000000000, 0.0375000000000, 0.0000000000000, 0.0},
  {0.0000000000000, 0.0125000000000, 0.0000000000000, 0.0},
  {0.1594544658298, 0.0289863616457, 0.0040664827466, 0.0},
  {-0.1594544658298, 0.0289863616457, 0.0040664827466, 0.0},
  {0.1594544658298, 0.0210136383543, 0.0040664827466, 0.0},
  {-0.1594544658298, 0.0210136383543, 0.0040664827466, 0.0},
  {0.3808991135940, 0.0345224778399, 0.0336223164632, 0.0},
  {-0.3808991135940, 0.0345224778399, 0.0336223164632, 0.0},
  {0.3808991135940, 0.0154775221601, 0.0336223164632, 0.0},
  {-0.3808991135940, 0.0154775221601, 0.0336223164632, 0.0},
  {0.6582769255267, 0.0414569231382, 0.0332008041365, 0.0},
  {-0.6582769255267, 0.0414569231382, 0.0332008041365, 0.0},
  {0.6582769255267, 0.0085430768618, 0.0332008041365, 0.0},
  {-0.6582769255267, 0.0085430768618, 0.0332008041365, 0.0},
  {0.8761473165029, 0.0469036829126, 0.0140936869250, 0.0},
  {-0.8761473165029, 0.0469036829126, 0.0140936869250, 0.0},
  {0.8761473165029, 0.0030963170874, 0.0140936869250, 0.0},
  {-0.8761473165029, 0.0030963170874, 0.0140936869250, 0.0},
  {0.9982431840532, 0.0499560796013, 0.0009770697703, 0.0},
  {-0.9982431840532, 0.0499560796013, 0.0009770697703, 0.0},
  {0.9982431840532, 0.0000439203987, 0.0009770697703, 0.0},
  {-0.9982431840532, 0.0000439203987, 0.0009770697703, 0.0},
  {0.9790222658168, 0.0412307108141, 0.0075319969436, 0.0},
  {-0.9790222658168, 0.0412307108141, 0.0075319969436, 0.0},
  {0.9790222658168, 0.0087692891859, 0.0075319969436, 0.0},
  {-0.9790222658168, 0.0087692891859, 0.0075319969436, 0.0},
  {0.6492284325645, 0.0494755566454, 0.0075319969436, 0.0},
  {-0.6492284325645, 0.0494755566454, 0.0075319969436, 0.0},
  {0.6492284325645, 0.0005244433546, 0.0075319969436, 0.0},
  {-0.6492284325645, 0.0005244433546, 0.0075319969436, 0.0},
  {0.8727421201131, 0.0339565366147, 0.0257718308672, 0.0},
  {-0.8727421201131, 0.0339565366147, 0.0257718308672, 0.0},
  {0.8727421201131, 0.0160434633853, 0.0257718308672, 0.0},
  {-0.8727421201131, 0.0160434633853, 0.0257718308672, 0.0},
  {0.3582614645881, 0.0468185530028, 0.0257718308672, 0.0},
  {-0.3582614645881, 0.0468185530028, 0.0257718308672, 0.0},
  {0.3582614645881, 0.0031814469972, 0.0257718308672, 0.0},
  {-0.3582614645881, 0.0031814469972, 0.0257718308672, 0.0},
  {0.5666666666667, 0.0301944444444, 0.0156250000000, 0.0},
  {-0.5666666666667, 0.0301944444444, 0.0156250000000, 0.0},
  {0.5666666666667, 0.0198055555556, 0.0156250000000, 0.0},
  {-0.5666666666667, 0.0198055555556, 0.0156250000000, 0.0},
  {0.2077777777778, 0.0391666666667, 0.0156250000000, 0.0},
  {-0.2077777777778, 0.0391666666667, 0.0156250000000, 0.0},
  {0.2077777777778, 0.0108333333333, 0.0156250000000, 0.0},
  {-0.2077777777778, 0.0108333333333, 0.0156250000000, 0.0},
  {0.0000000000000, 0.2750000000000, 0.0084492309003, 0.0},
  {-0.2517129343453, 0.2750000000000, 0.0237714740190, 0.0},
  {0.2517129343453, 0.2750000000000, 0.0237714740190, 0.0},
  {-0.7013933644534, 0.2750000000000, 0.0294001617014, 0.0},
  {0.7013933644534, 0.2750000000000, 0.0294001617014, 0.0},
  {0.0000000000000, 0.2183645897723, 0.0237714740190, 0.0},
  {0.0000000000000, 0.3316354102277, 0.0237714740190, 0.0},
  {0.0000000000000, 0.1171864929980, 0.0294001617014, 0.0},
  {0.0000000000000, 0.4328135070020, 0.0294001617014, 0.0},
  {0.9590960631620, 0.2750000000000, 0.0066444364658, 0.0},
  {-0.9590960631620, 0.2750000000000, 0.0066444364658, 0.0},
  {0.0000000000000, 0.4907966142114, 0.0066444364658, 0.0},
  {0.0000000000000, 0.0592033857886, 0.0066444364658, 0.0},
  {0.9956010478552, 0.2750000000000, 0.0042536044255, 0.0},
  {-0.9956010478552, 0.2750000000000, 0.0042536044255, 0.0},
  {0.0000000000000, 0.4990102357674, 0.0042536044255, 0.0},
  {0.0000000000000, 0.0509897642326, 0.0042536044255, 0.0},
  {0.5000000000000, 0.2750000000000, 0.0000000000000, 0.0},
  {-0.5000000000000, 0.2750000000000, 0.0000000000000, 0.0},
  {0.0000000000000, 0.3875000000000, 0.0000000000000, 0.0},
  {0.0000000000000, 0.1625000000000, 0.0000000000000, 0.0},
  {0.1594544658298, 0.3108772548117, 0.0040664827466, 0.0},
  {-0.1594544658298, 0.3108772548117, 0.0040664827466, 0.0},
  {0.1594544658298, 0.2391227451883, 0.0040664827466, 0.0},
  {-0.1594544658298, 0.2391227451883, 0.0040664827466, 0.0},
  {0.3808991135940, 0.3607023005587, 0.0336223164632, 0.0},
  {-0.3808991135940, 0.3607023005587, 0.0336223164632, 0.0},
  {0.3808991135940, 0.1892976994413, 0.0336223164632, 0.0},
  {-0.3808991135940, 0.1892976994413, 0.0336223164632, 0.0},
  {0.6582769255267, 0.4231123082435, 0.0332008041365, 0.0},
  {-0.6582769255267, 0.4231123082435, 0.0332008041365, 0.0},
  {0.6582769255267, 0.1268876917565, 0.0332008041365, 0.0},
  {-0.6582769255267, 0.1268876917565, 0.0332008041365, 0.0},
  {0.8761473165029, 0.4721331462132, 0.0140936869250, 0.0},
  {-0.8761473165029, 0.4721331462132, 0.0140936869250, 0.0},
  {0.8761473165029, 0.0778668537868, 0.0140936869250, 0.0},
  {-0.8761473165029, 0.0778668537868, 0.0140936869250, 0.0},
  {0.9982431840532, 0.4996047164120, 0.0009770697703, 0.0},
  {-0.9982431840532, 0.4996047164120, 0.0009770697703, 0.0},
  {0.9982431840532, 0.0503952835880, 0.0009770697703, 0.0},
  {-0.9982431840532, 0.0503952835880, 0.0009770697703, 0.0},
  {0.9790222658168, 0.4210763973270, 0.0075319969436, 0.0},
  {-0.9790222658168, 0.4210763973270, 0.0075319969436, 0.0},
  {0.9790222658168, 0.1289236026730, 0.0075319969436, 0.0},
  {-0.9790222658168, 0.1289236026730, 0.0075319969436, 0.0},
  {0.6492284325645, 0.4952800098088, 0.0075319969436, 0.0},
  {-0.6492284325645, 0.4952800098088, 0.0075319969436, 0.0},
  {0.6492284325645, 0.0547199901912, 0.0075319969436, 0.0},
  {-0.6492284325645, 0.0547199901912, 0.0075319969436, 0.0},
  {0.8727421201131, 0.3556088295323, 0.0257718308672, 0.0},
  {-0.8727421201131, 0.3556088295323, 0.0257718308672, 0.0},
  {0.8727421201131, 0.1943911704677, 0.0257718308672, 0.0},
  {-0.8727421201131, 0.1943911704677, 0.0257718308672, 0.0},
  {0.3582614645881, 0.4713669770255, 0.0257718308672, 0.0},
  {-0.3582614645881, 0.4713669770255, 0.0257718308672, 0.0},
  {0.3582614645881, 0.0786330229745, 0.0257718308672, 0.0},
  {-0.3582614645881, 0.0786330229745, 0.0257718308672, 0.0},
  {0.5666666666667, 0.3217500000000, 0.0156250000000, 0.0},
  {-0.5666666666667, 0.3217500000000, 0.0156250000000, 0.0},
  {0.5666666666667, 0.2282500000000, 0.0156250000000, 0.0},
  {-0.5666666666667, 0.2282500000000, 0.0156250000000, 0.0},
  {0.2077777777778, 0.4025000000000, 0.0156250000000, 0.0},
  {-0.2077777777778, 0.4025000000000, 0.0156250000000, 0.0},
  {0.2077777777778, 0.1475000000000, 0.0156250000000, 0.0},
  {-0.2077777777778, 0.1475000000000, 0.0156250000000, 0.0}
  };

/*********  end of local variable declaration   **********/

  overallSW = swCreate ("overall");	/* Overall performance stopwatch */
  startTime = time (NULL);

  /* Add an exit handler to deal with wayward children. */

  if (atexit (exit_kelvin)) {
    perror ("Could not register exit handler!");
    exit (EXIT_FAILURE);
  }

  /* Fork a child that loops sleeping several seconds and then signalling 
     us with SIGUSR1 to do an asynchronous dump of peak statistitics to stderr. */

  if ((maximumVMK = swGetMaximumVMK ()) != 0) {
    childPID = fork ();
    if (childPID == 0) {
      /* Code executed by child only! */
      pid_t parentPID = 0;
      /* Ignore QUIT signals, 'cause they're actually status requests for Mom. */
      signal(SIGQUIT, SIG_IGN);

      while (1) {
	sleep (30);
	parentPID = getppid ();
	if (parentPID == 1)
	  exit (EXIT_SUCCESS);
	currentVMK = swGetCurrentVMK (parentPID);
	fprintf (stderr, "%lus, %dKb (%2d%% of %2.1fGb)\n",
		 time (NULL) - startTime,
		 currentVMK, (currentVMK * 100) / maximumVMK,
		 maximumVMK / (1024.0 * 1024.0));
      }
    }
  }

  /* Setup signal handlers */
  struct sigaction quitAction, intAction;
  sigset_t quitBlockMask, intBlockMask;

#if defined (GPROF) || (GCOV)
  struct sigaction termAction;
  sigset_t termBlockMask;
#endif
  sigfillset (&quitBlockMask);
  quitAction.sa_handler = quitSignalHandler;
  quitAction.sa_mask = quitBlockMask;
  quitAction.sa_flags = 0;
  sigaction (SIGQUIT, &quitAction, NULL);

  sigfillset (&intBlockMask);
  intAction.sa_handler = intSignalHandler;
  intAction.sa_mask = intBlockMask;
  intAction.sa_flags = 0;
  sigaction (SIGINT, &intAction, NULL);

#if defined (GPROF) || (GCOV)
  sigfillset (&termBlockMask);
  termAction.sa_handler = termSignalHandler;
  termAction.sa_mask = termBlockMask;
  termAction.sa_flags = 0;
  sigaction (SIGTERM, &termAction, NULL);
#endif

  /* Annouce ourselves for performance tracking. */
  char currentWorkingDirectory[MAXSWMSG-32];
  sprintf (messageBuffer, "dkelvin %s built %s %s",
	   programVersion, __DATE__, __TIME__);
  swLogMsg (messageBuffer);
  swLogMsg (kelvinVersion);
  swLogMsg (likelihoodVersion);
  swLogMsg (locusVersion);
  swLogMsg (polynomialVersion);

#ifdef _OPENMP
  if ((envVar = getenv ("OMP_NUM_THREADS")) != NULL)
    threadCount = atoi (envVar);
  sprintf (messageBuffer, "OpenMP-enabled w/%d threads.", threadCount);
  swLogMsg (messageBuffer);
#endif

#ifdef DMTRACK
  swLogMsg
    ("Dynamic memory usage dumping is turned on, so performance will be poor!");
#endif
#ifdef GPROF
  sprintf (messageBuffer,
	   "GNU profiler (gprof) run, use \"kill -%d %d\" to finish early.",
	   SIGTERM, getpid ());
  swLogMsg (messageBuffer);
#endif
#ifdef GCOV
  sprintf (messageBuffer,
	   "GNU coverage analyzer (gcov) run, use \"kill -%d %d\" to finish early.",
	   SIGTERM, getpid ());
  swLogMsg (messageBuffer);
#endif
  fprintf (stderr,
	   "To force a dump of stats (at some risk), type CTRL-\\ or type \"kill -%d %d\".\n",
	   SIGQUIT, getpid ());
  swStart (overallSW);

  memset (&savedLocusList, 0, sizeof (savedLocusList));
  memset (&markerLocusList, 0, sizeof (markerLocusList));
  memset (&traitLocusList, 0, sizeof (traitLocusList));

  memset (&modelType, 0, sizeof (modelType));

#ifdef DEBUG
  //  mtrace();
#endif

  /* Initialize the logging system. */
  logInit ();
  /* logSet(LOGINPUTFILE, LOGDEBUG); */
  //logSet(LOGGENOELIM, LOGDEBUG);
  /* logSet(LOGPEELGRAPH, LOGFATAL); */
  //logSet (LOGLIKELIHOOD, LOGWARNING);
  //logSet(LOGLIKELIHOOD, LOGDEBUG); 
  //logSet(LOGPARENTALPAIR, LOGDEBUG); 
  //logSet(LOGSETRECODING, LOGDEBUG);

  /* Start by parsing command line arguments. Most essential: figure
   * out where the configuration file lives. */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case '?':
	/* Help */
	fprintf (stdout, "Usage:\n");
	fprintf (stdout, "  %s [-?] <configuration file>\nwhere:\n", argv[0]);
	fprintf (stdout, "      -? : this output;\n");
	fprintf (stdout,
		 "      <configuration file> : file containing run parameters.\n");
	exit (EXIT_FAILURE);
	break;
      }
    } else if (strlen (configfile) != 0) {
      /* Unexpected argument; we already have a configuration file! Punt. */
      KLOG (LOGDEFAULT, LOGFATAL,
	    "Unexpected command line argument '%s'; aborting.\n", argv[i]);
    } else if (strlen (argv[i]) >= KMAXFILENAMELEN) {
      /* Configuration file name too long! Punt. */
      KLOG (LOGDEFAULT, LOGFATAL,
	    "Configuration file name '%s' exceeds limit of %d; aborting.\n",
	    argv[i], KMAXFILENAMELEN);
    } else {
      /* Got a configuration file name. Copy it. */
      strncpy (configfile, argv[i], KMAXFILENAMELEN);
      getcwd (currentWorkingDirectory, sizeof(currentWorkingDirectory));
      sprintf(messageBuffer, "In %s w/%s", currentWorkingDirectory, configfile);
      swLogMsg(messageBuffer);
    }
    i++;
  }

  /* Check to see if the configuration file name was specified. */
  KASSERT ((strlen (configfile) > 0),
	   "No configuration file specified; aborting.\n");

  /* set the default unknown person ID */
  modelOptions.sUnknownPersonID = malloc (sizeof (char) * 2);
  strcpy (modelOptions.sUnknownPersonID, "0");

  /* set default values for PPL calculations */
  /* LRs are weighted heavier for theta less than the cutoff */
  modelOptions.thetaCutoff[0] = 0.05;
  modelOptions.thetaCutoff[1] = 0.05;
  /* weight ofr theta less than the cutoff */
  modelOptions.thetaWeight = 0.95;
  /* prior probability of linkage */
  modelOptions.prior = 0.02;
  /* prior probability of LD given close linkage */
  modelOptions.LDprior = 0.02;

  /* set default for QT */
  modelType.minOriginal = -999999999.00;
  modelType.maxOriginal = 999999999.00;
  modelType.minThreshold = -999999999.00;
  modelType.maxThreshold = 999999999.00;

  /* Parse the configuration file. */
  KASSERT (readConfigFile (configfile, &modelType, &modelRange, &modelOptions)
	   != ERROR, "Error in configuration file; aborting.\n");

  /* For now, reject all models we can't deal with. So use KASSERT to
   * check that we're looking at, e.g., twopoint dichotomous models
   * with direct eval (not polynomial --> add this to model?) and
   * give appropriate error message otherwise. */
  KASSERT (modelRange.nalleles == 2, "Only biallelic traits supported.\n");

  /* the difference between QT and CT is whether we use threshold or not. Under CT -  yes to
   * threshold, under QT - no threshold */
  if (modelRange.ntthresh > 0 && modelType.trait != DT) {
    modelType.trait = CT;
    KASSERT (modelType.minThreshold > -999999998
	     && modelType.maxThreshold < 999999998,
	     "Under QT threshold model, MIN and MAX of the QT threshold values need to be provided through keywords T_MIN and T_MAX.\n");
  }

  total_dim = 2;		// alpha gf
  total_dim += 3 * modelRange.nlclass;	//DD Dd dd
  if (modelType.type == TP) {
    total_dim += 1;		// theta;
    if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
      total_dim += 1;		// dprime
    }
  }

  if (modelType.trait != DT) {
    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      total_dim += 3 * modelRange.nlclass;	//SD_DD SD_Dd SD_dd
    }
    if (modelType.trait == CT) {
      total_dim++;
    }
  }

  fprintf (stderr, "total dim =%d\n", total_dim);

  if (modelType.trait == QT) {
    /* threshold value will not be used in any meaningful way, but we will use it for 
       the loop */
    modelRange.ntthresh = 1;
    modelType.minOriginal = 0;
    modelType.maxOriginal = 1;
    if (modelRange.tthresh == NULL) {
      modelRange.tthresh = (double **) malloc (sizeof (double *));
      for (i = 0; i < modelRange.nlclass; i++) {
	modelRange.tthresh[i] = malloc (sizeof (double));
      }
    }
  }

  if (modelType.trait != DT) {
    /* Setting ranges for each variables. Default is [0,1] */
    k = 1;
    for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	xl[k] = xl[k + 1] = xl[k + 2] = -3;
	xu[k] = xu[k + 1] = xu[k + 2] = 3;
      } else {
	xl[k] = xl[k + 1] = xl[k + 2] = 0.1;
	xu[k] = xu[k + 1] = xu[k + 2] = 30;
      }
      volume_region *= (xu[k] - xl[k]);
      volume_region *= (xu[k + 1] - xl[k + 1]);
      volume_region *= (xu[k + 2] - xl[k + 2]);
      k += 3;
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	xl[k] = xl[k + 1] = xl[k + 2] = 0.5;
	xu[k] = xu[k + 1] = xu[k + 2] = 3.0;
	volume_region *= (xu[k] - xl[k]);
	volume_region *= (xu[k + 1] - xl[k + 1]);
	volume_region *= (xu[k + 2] - xl[k + 2]);
	k += 3;
      }
      if (modelType.trait == CT) {
	xl[k] = 10;
	xu[k] = 30.0;
	volume_region *= (xu[k] - xl[k]);
	k++;
	//   fprintf(stderr, " in CT\n ");

      }
    }

    fprintf (stderr,
	     "The number of dimension for calculation of BR should be %d\n",
	     k);
  }

  fpHet = fopen (avghetfile, "w");
  KASSERT (fpHet != NULL,
	   "Error in opening file Theta result file for write.\n");
  //  fprintf (fpHet, "# Version %s\n", programVersion);

  if (print_point_flag)
    fphlod = fopen ("hlod.pts", "w");
  //  fprintf (fphlod, "# Version %s\n", programVersion);


  if (modelType.type == TP) {
    fpPPL = fopen (pplfile, "w");
  //  fprintf (fpPPL, "# Version %s\n", programVersion);
    KASSERT (fpPPL != NULL, "Error in opening file %s for write.\n", pplfile);
    fprintf (fpPPL, "%4s %15s %9s %6s ", "CHR", "MARKER", "cM", "PPL");
    if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
      fprintf (fpPPL, "%6s %6s ", "LD-PPL", "PPLD");
    }
    fprintf (fpPPL, " MOD \n");
    fflush (fpPPL);
  }
  if (modelOptions.polynomial == TRUE) {
    polynomialInitialization ();
    fprintf (stderr,
	     "!!!!!!!!!!!The Computation is done in polynomial mode!!!!!!!!!!!!!!!\n");
  } else {
    fprintf (stderr, "Polynomial is off!\n");
  }

  /* Read in the map file. */
  read_mapfile (mapfile);

  /* Initialize the locus list and read in the marker file. */
  memset (&originalLocusList, 0, sizeof (originalLocusList));
  /* read in what loci are in the pedigree file */
  read_datafile (datafile);


  /* The configuration has all the information about the disease trait if any */
  if (originalLocusList.numTraitLocus > 0) {
    /* we are not doing marker to marker analysis
     * Need to add the alleles into trait locus 
     * Assume the traitLoucs is 0 for now  - Need to fix this later */
    traitLocus = 0;
    pLocus = originalLocusList.ppLocusList[traitLocus];
    pTraitLocus = pLocus->pTraitLocus;
    add_allele (pLocus, "D", 0.5);
    add_allele (pLocus, "d", 0.5);
    /* fix number of trait variables at 1 for now */
    pTraitLocus->numTrait = 1;
    pTrait = add_trait (0, pTraitLocus, modelType.trait);
    pTrait->numLiabilityClass = modelRange.nlclass;
    if (modelType.trait == QT || modelType.trait == CT) {
      modelType.min = (modelType.minOriginal - modelType.mean) / modelType.sd;
      modelType.max = (modelType.maxOriginal - modelType.mean) / modelType.sd;
      pTrait->minFlag = modelType.minFlag;
      pTrait->maxFlag = modelType.maxFlag;
      pTrait->min = modelType.min;
      pTrait->max = modelType.max;
      pTrait->functionQT = modelType.distrib;
      if (modelType.distrib == QT_FUNCTION_T)
	pTrait->dfQT = modelType.constants[0];
      pTrait->sampleMean = modelType.mean;
      pTrait->sampleSD = modelType.sd;
      pTrait->unknownTraitValue =
	modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN];
      pTrait->lessCutoffFlag =
	modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED];
      pTrait->moreCutoffFlag =
	modelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED];
    }
  }

  /* read in marker allele frequencies */
  read_markerfile (markerfile, modelType.numMarkers);

  /* build allele set information */
  for (locus = 0; locus < originalLocusList.numLocus; locus++) {
    construct_original_allele_set_list (locus);
  }

  /* Initialize the pedigree set datastructure and read in the pedigrees. */
  memset (&pedigreeSet, 0, sizeof (PedigreeSet));
  read_pedfile (pedfile, &pedigreeSet);
  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
    if (pPedigree->currentLoopFlag) {
      fprintf(stderr, "Pedigree %s has at least one loop not broken yet\n", pPedigree->sPedigreeID);
      exitDueToLoop = TRUE;
    }
  }
  KASSERT (exitDueToLoop == FALSE, "Not all loops in pedigrees are broken.\n");

  /* read in case control file if provided */
  if (strlen (ccfile) > 0) {
    read_ccfile (ccfile, &pedigreeSet);
    modelType.ccFlag = 1;
  }
  flexBufferSize = 0;
  free (flexBuffer);
  fflush (stderr);
  fflush (stdout);

  /* allocate space for results */
  if (modelType.type == TP) {
    modelType.numMarkers = 1;
    totalLoci = 2;
    /* two point analysis */
    if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) {
      /* in order to simplify looping, even for LE, we add a fake LD parameter dprime=0, which
       * is LE */
      modelRange.ndprime = 1;
      modelRange.dprime = (double *) calloc (1, sizeof (double));
      modelRange.dprime[0] = 0;
      pLambdaCell = findLambdas (&modelRange, 2, 2);
      dprime0Idx = 0;
    }
  } else {
    /* we are doing multipoint analysis */
    totalLoci = modelType.numMarkers + originalLocusList.numTraitLocus;
    if (modelRange.tlmark == TRUE) {
      /* add marker positions to the list of positions we want to conduct analysis */
      for (i = 0; i < originalLocusList.numLocus; i++) {
	pLocus = originalLocusList.ppLocusList[i];
	if (pLocus->locusType == LOCUS_TYPE_TRAIT)
	  continue;
	addTraitLocus (&modelRange, pLocus->pMapUnit->mapPos[SEX_AVERAGED]);
      }
    }
  }

  /* allocate storage for keeping track of het locus in nuclear families */
  allocate_nucfam_het (&pedigreeSet, totalLoci);

  /* initialize some work space */
  initialize_parental_pair_workspace (&parentalPairSpace,
				      originalLocusList.numLocus);
  /* allocate transmission probability matrix */
  build_xmission_matrix (&nullMatrix, totalLoci);
  build_xmission_matrix (&altMatrix, totalLoci);
  build_xmission_matrix (&traitMatrix, 1);
  build_xmission_matrix (&markerMatrix, totalLoci - 1);
  xmissionMatrix = nullMatrix;
  tmpID = (char *) calloc (totalLoci, sizeof (char));

  /* initialize loci by doing genotype elimination, set recoding */
  initialize_loci (&pedigreeSet);

  if (modelOptions.dryRun != 0) {
    for (loc1 = 0; loc1 < originalLocusList.numLocus; loc1++) {
      fprintf (stderr, "Locus %d:\n", loc1);
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
	pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
	print_pedigree_locus_genotype_count (pPedigree, loc1);
      }

    }
  }
  if (modelOptions.polynomial == TRUE) {
    //      constant1Poly = constantExp (1);
    // constant0Poly = constantExp (0);
    for (k = 0; k < 3; k++) {
      initialProbPoly[k] = constant1Poly;
      initialProbPoly2[k] = constant1Poly;
      initialProbAddr[k] = initialProbPoly[k];
      initialProbAddr2[k] = initialProbPoly2[k];
      initialHetProbAddr[k] = NULL;
    }
  } else {
    for (k = 0; k < 3; k++) {
      initialProb[k] = 1.0;
      initialProb2[k] = 1.0;
      initialProbAddr[k] = &initialProb[k];
      initialProbAddr2[k] = &initialProb2[k];
      initialHetProbAddr[k] = NULL;
    }
  }


  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
    pPedigree->load_flag = 0;
  }

  /* only for multipoint - we don't handle LD under multipoint yet */
  /* DCUHRE do now use likelihoodDT or likelihoodQT to store null likelihoods */

  /* find out the max we need to allocate */
  /* after genotype lists have been built, we want to pre-allocate parental pair work space
   * it used to be done dynamically, but it's too costly 
   * parental pair is based on each nuclear family */
  stat_parental_pair_workspace (&pedigreeSet);

  /* after genotype elimination and tracking the max work space needed 
   * for constructing parental pair */
  allocate_parental_pair_workspace (&parentalPairSpace,
				    modelType.numMarkers + 1);

  /* conditional likelihood storage space for each individual */
  allocate_likelihood_space (&pedigreeSet, modelType.numMarkers + 1);

  time0 = clock ();
  time1 = clock ();


  if (modelType.trait == DT)
    fprintf (stderr, "Dichotomous Trait & ");
  else if (modelType.trait == QT)
    fprintf (stderr, "Quantitative Trait without threshold & ");
  else
    fprintf (stderr, "Quantitative Trait with threshold & ");

  fprintf (stderr, "%s\n",
	   (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) ? "LE" : "LD");
  fprintf (stderr, "Total number of markers in data: %d\n",
	   originalLocusList.numLocus - originalLocusList.numTraitLocus);
  fprintf (stderr, "Total number of trait locus in data: %d\n",
	   originalLocusList.numTraitLocus);
  fprintf (stderr, "Total number of families in pedigree file: %d\n",
	   pedigreeSet.numPedigree);
  fprintf (stderr, "Number of loci to use for analysis: %d\n",
	   modelType.numMarkers + originalLocusList.numTraitLocus);


  /* Initialize the connection to the infrastructure. */
  /* niceInit (); */

  /* assume the trait locus is the first one in the list */
  traitLocus = 0;
  pLocus = originalLocusList.ppLocusList[traitLocus];
  pTraitLocus = originalLocusList.ppLocusList[traitLocus]->pTraitLocus;
  pTrait = pTraitLocus->pTraits[traitLocus];
  if (modelType.type == TP) {
    /* Two point. */
    if (originalLocusList.pLDLoci == NULL) {
      originalLocusList.pLDLoci = (LDLoci *) malloc (sizeof (LDLoci));
      memset (originalLocusList.pLDLoci, 0, sizeof (LDLoci));
    }
    pLDLoci = &originalLocusList.pLDLoci[0];
    originalLocusList.numLDLoci = 1;

    if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) {
      /* fake some LD information to simplify looping */
      pLDLoci->numAllele1 = 2;
      pLDLoci->ppDPrime = (double **) malloc (sizeof (double *));
      pLDLoci->ppDPrime[0] = (double *) malloc (sizeof (double));
      pLDLoci->ppDValue = (double **) malloc (sizeof (double *));
      pLDLoci->ppDValue[0] = (double *) malloc (sizeof (double));
      pLDLoci->ppHaploFreq = (double **) malloc (sizeof (double *) * 2);
      pLDLoci->ppHaploFreq[0] = (double *) malloc (sizeof (double) * 2);
      pLDLoci->ppHaploFreq[1] = (double *) malloc (sizeof (double) * 2);

      /* initialize it */
      pLDLoci->ppDPrime[0][0] = 0;
    }

    locusList = &savedLocusList;
    savedLocusList.numLocus = 2;
    savedLocusList.pLocusIndex =
      (int *) malloc (sizeof (int) * savedLocusList.numLocus);
    for (i = 0; i < 3; i++) {
      savedLocusList.pPrevLocusDistance[i] =
	(double *) malloc (sizeof (double) * savedLocusList.numLocus);
      savedLocusList.pNextLocusDistance[i] =
	(double *) malloc (sizeof (double) * savedLocusList.numLocus);
      savedLocusList.pPrevLocusDistance[i][0] = -1;
      savedLocusList.pNextLocusDistance[i][1] = -1;
    }

    if (modelOptions.polynomial == TRUE) {
      /* populate the matrix */
      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
					 initialProbAddr2,	/* probability */
					 initialHetProbAddr, 0,	/* cell index */
					 -1,	/* last het locus */
					 -1,	/* last  pattern (P-1 or M-2) */
					 0);	/* current locus - start with 0 */
      fprintf (stderr,
	       "holdAllPolys from population of transmission matrix\n");
      holdAllPolys ();
    }
    //total_count = modelRange.npenet * modelRange.ngfreq * modelRange.nalpha;

    if (modelOptions.markerAnalysis == FALSE) {
      savedLocusList.traitLocusIndex = 0;
      savedLocusList.traitOrigLocus = 0;
    } else {
      savedLocusList.traitLocusIndex = -1;
      savedLocusList.traitOrigLocus = -1;
    }

    for (loc1 = 0; loc1 < originalLocusList.numLocus - 1; loc1++) {

      savedLocusList.pLocusIndex[0] = loc1;
      pLocus1 = originalLocusList.ppLocusList[loc1];
      if (modelOptions.markerAnalysis != FALSE
	  && pLocus1->locusType != LOCUS_TYPE_MARKER)
	continue;

      for (loc2 = loc1 + 1; loc2 < originalLocusList.numLocus; loc2++) {

	maximum_function_value = 0.0;

	pLocus2 = originalLocusList.ppLocusList[loc2];
	if (pLocus2->locusType != LOCUS_TYPE_MARKER)
	  continue;
	savedLocusList.pLocusIndex[1] = loc2;

	/* find out number of alleles this marker locus has *//* Check if this is okay with DCUHRE  ???????????? */
	if (modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM) {
	  /* get the LD parameters */
	  pLambdaCell =
	    findLambdas (&modelRange, pLocus1->numOriginalAllele,
			 pLocus2->numOriginalAllele);
	  reallocate_LD_loci (pLDLoci, pLocus1->numOriginalAllele,
			      pLocus2->numOriginalAllele);
	  pLDLoci->locus1 = loc1;
	  pLDLoci->locus2 = loc2;
	  pLDLoci->numAllele1 = pLocus1->numOriginalAllele;
	  pLDLoci->numAllele2 = pLocus2->numOriginalAllele;
	  if (pLocus1->numOriginalAllele == 2
	      && pLocus2->numOriginalAllele == 2)
	    R_square_flag = TRUE;
	  else
	    R_square_flag = FALSE;
	}

	loopMarkerFreqFlag = 0;
	if (modelRange.nafreq >= 2
	    && modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM
	    && pLocus2->numOriginalAllele == 2) {
	  loopMarkerFreqFlag = 1;
	} else if (modelRange.nafreq == 0) {
	  /* add a fake one to facilitate loops and other handlings */
	  addAlleleFreq (&modelRange, pLocus2->pAlleleFrequency[0]);
	} else {
	  modelRange.nafreq = 1;
	  modelRange.afreq[0] = pLocus2->pAlleleFrequency[0];
	}

	/* allocate/initialize result storage */
	// initialize_tp_result_storage ();

	/* we will force marker allele frequency loop to execute at least once */
	for (mkrFreqIdx = 0;
	     mkrFreqIdx == 0 || mkrFreqIdx < modelRange.nafreq;
	     mkrFreqIdx++) {
	  mkrFreq = pLocus2->pAlleleFrequency[0];
	  /* we should only loop over marker allele frequency under twopoint
	   * and when markers are SNPs (only have two alleles) */
	  if (loopMarkerFreqFlag) {
	    mkrFreq = modelRange.afreq[mkrFreqIdx];
	    /* update the locus */
	    pLocus2->pAlleleFrequency[0] = mkrFreq;
	    pLocus2->pAlleleFrequency[1] = 1 - mkrFreq;
	    if (modelOptions.polynomial == TRUE);
	    else
	      update_locus (&pedigreeSet, loc2);
	  }
	  fprintf (stderr, "mkrFreq =%d  model nafreq= %d \n",
		   mkrFreqIdx, modelRange.nafreq);


	  if (1 && modelOptions.markerAnalysis == FALSE) {

	    if (modelOptions.polynomial == TRUE);
	    else
	      update_locus (&pedigreeSet, loc1);
	  }

	  /* clear Dprime combination impossible flag */
	  memset (pLambdaCell->impossibleFlag, 0,
		  sizeof (int) * pLambdaCell->ndprime);
	  /* set up haplotype frequencies */
	  for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
	    if (isDPrime0
		(pLambdaCell->lambda[dprimeIdx], pLambdaCell->m,
		 pLambdaCell->n))
	      dprime0Idx = dprimeIdx;
	    status =
	      setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprimeIdx);
	    if (status < 0) {
	      pLambdaCell->impossibleFlag[dprimeIdx] = 1;
	    }
	  }

	  /* for each D prime and theta, print out average and maximizing model information - MOD */
	  fprintf (fpHet, "# %-d  %s %s \n", loc2, pLocus1->sName,
		   pLocus2->sName);
	  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
	    fprintf (fpHet, "Dprime ");
	  }
	  if (modelType.trait == DICHOTOMOUS) {
	    fprintf (fpHet, "%6s %6s %8s %8s %8s %6s %6s %5s %5s %5s \n",
		     "Theta", "COUNT", "BR", "ERR_EST", "MAX_HLOD", "ALPHA",
		     "DGF", "PEN_DD", "PEN_Dd", "PEN_dd");
	  } else {
	    fprintf (fpHet, "%6s %6s %8s %8s %8s %6s %6s %5s %5s %5s ",
		     "Theta", "COUNT", "BR", "ERR_EST", "MAX_HLOD", "ALPHA",
		     "DGF", "MEAN_DD", "MEAN_Dd", "MEAN_dd");
	    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	      fprintf (fpHet, " %5s %5s %5s ", "SD_DD", "SD_Dd", "SD_dd");
	    }
	    if (modelType.trait == CT) {
	      fprintf (fpHet, "  %5s", "t");
	    }
	    fprintf (fpHet, "\n");
	  }

	  low_theta_integral = 0.0;
	  high_theta_integral = 0.0;
	  low_integral = 0.0;
	  high_integral = 0.0;
	  low_ld_integral = 0.0;

	  for (i = 0; i < 140; i++) {
	    fixed_dprime = dcuhre2[i][0];
	    fixed_theta = dcuhre2[i][1];

	    integral = 0.0;
	    abserr = 0.0;
	    fprintf (stderr, "i=%d Dprime=%f theta=%f \n", i,
		     fixed_dprime, fixed_theta);
	    // fprintf(fpSeok_theta,"%f %f ",      fixed_dprime, fixed_theta);

	    if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
	      for (dprimeIdx = 0; dprimeIdx < 33; dprimeIdx++) {
		if (fabs
		    (pLambdaCell->lambda[dprimeIdx][0][0] -
		     fixed_dprime) < 0.0001) {
		  // fprintf(stderr,"dprimeIdx =%d with %15.13f which is matching wit fixed_dprime\n",dprimeIdx,pLambdaCell->lambda[dprimeIdx][0][0]);
		  break;
		}
	      }
	      if (dprimeIdx == 33) {
		fprintf (stderr, "dprimeIdx is %d\n", dprimeIdx);
		exit (0);
	      }
	    }

	    num_out_constraint = 0;
	    kelvin_dcuhre_integrate (&integral, &abserr);
	    for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	      integral *= 6;
	      abserr *= 6;
	    }
	    integral /= volume_region;
	    abserr /= volume_region;
	    dcuhre2[i][3] = integral;

	    if (modelType.trait == DICHOTOMOUS) {
	      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
		fprintf (fpHet, "%6.4f ", fixed_dprime);
	      }
	      fprintf (fpHet, "%6.4f %6d %8.4f %8.4f %8.4f %6.4f %6.4f  ",
		       fixed_theta, s->total_neval, integral, abserr,
		       log10 (localmax_value), localmax_x[1], localmax_x[0]);
	      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
		fprintf (fpHet, "%6.4f %6.4f %6.4f",
			 localmax_x[liabIdx * 3 + 2],
			 localmax_x[liabIdx * 3 + 3],
			 localmax_x[liabIdx * 3 + 4]);
	      }
	      fprintf (fpHet, "\n");
	      if (maximum_function_value < localmax_value) {
		maximum_function_value = localmax_value;
		maxima_x[0] = fixed_dprime;
		maxima_x[1] = fixed_theta;
		maxima_x[2] = localmax_x[0];
		maxima_x[3] = localmax_x[1];
		for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
		  maxima_x[liabIdx * 3 + 4] = localmax_x[liabIdx * 3 + 2];
		  maxima_x[liabIdx * 3 + 5] = localmax_x[liabIdx * 3 + 3];
		  maxima_x[liabIdx * 3 + 6] = localmax_x[liabIdx * 3 + 4];
		}
	      }
	    } else {		//QT
	      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
		fprintf (fpHet, "%6.4f ", fixed_dprime);
	      }
	      fprintf (fpHet, "%6.4f %6d %8.4f %8.4f %8.4f %6.4f %6.4f ",
		       fixed_theta, s->total_neval, integral, abserr,
		       log10 (localmax_value), localmax_x[1], localmax_x[0]);
	      j = 2;
	      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
		fprintf (fpHet, "%6.4f %6.4f %6.4f ", localmax_x[j],
			 localmax_x[j + 1], localmax_x[j + 2]);
		j += 3;
		if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
		  fprintf (fpHet, "%6.4f %6.4f %6.4f ", localmax_x[j],
			   localmax_x[j + 1], localmax_x[j + 2]);
		  j += 3;
		}
		if (modelType.trait == CT) {
		  fprintf (fpHet, "%6.4f", localmax_x[j++]);
		}
	      }
	      fprintf (fpHet, "  %d\n", num_out_constraint);

	      if (maximum_function_value < localmax_value) {
		maximum_function_value = localmax_value;
		maxima_x[0] = fixed_dprime;
		maxima_x[1] = fixed_theta;
		maxima_x[2] = localmax_x[0];	// gf
		maxima_x[3] = localmax_x[1];	// alpha
		j = 2;
		for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
		  maxima_x[j + 2] = localmax_x[j];	//mean_DD
		  maxima_x[j + 3] = localmax_x[j + 1];	//mean_Dd
		  maxima_x[j + 4] = localmax_x[j + 2];	//mean_dd
		  j += 3;
		  if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
		    maxima_x[j + 2] = localmax_x[j];	// SD_DD
		    maxima_x[j + 3] = localmax_x[j + 1];	// SD_Dd        
		    maxima_x[j + 4] = localmax_x[j + 2];	// SD_dd
		    j += 3;
		  }
		  if (modelType.trait == CT) {
		    maxima_x[j + 2] = localmax_x[j];	// t
		    j++;
		  }
		}
	      }
	    }			/* End of writing max */
	    fflush (fpHet);
	    fprintf (stderr, "tp result %f %f is %13.10f   \n",
		     fixed_theta, fixed_dprime, integral);

	    if (i < 5) {
	      low_theta_integral += integral * dcuhre2[i][2];
	    } else if (i < 10) {
	      high_theta_integral += integral * dcuhre2[i][2];
	    } else {
	      if (fixed_theta < modelOptions.thetaCutoff[0]) {
		low_integral += integral * dcuhre2[i][2];
	      } else {
		high_integral += integral * dcuhre2[i][2];
	      }
	    }

	    if ((modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) && (i == 9)) {
	      printf ("End of LE case\n");
	      i = 140;
	    }




	  }			/* end of for to calculate BR(theta, dprime) */


	  /*Calculate ppl, ppld and ldppl */
	  ppl =
	    modelOptions.thetaWeight * low_theta_integral + (1 -
							     modelOptions.
							     thetaWeight)
	    * high_theta_integral;
	  ppl = ppl / (ppl + (1 - modelOptions.prior) / modelOptions.prior);
	  fprintf (fpPPL, "%4d %15s %9.4f %8.6f ",
		   pLocus2->pMapUnit->chromosome, pLocus2->sName,
		   pLocus2->pMapUnit->mapPos[SEX_AVERAGED], ppl);
	  fprintf (stderr, "ppl is %f\n", ppl);

	  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
	    ldppl =
	      modelOptions.thetaWeight * low_integral + (1 -
							 modelOptions.
							 thetaWeight)
	      * high_integral;
	    ldppl =
	      ldppl / (ldppl + (1 - modelOptions.prior) / modelOptions.prior);
	    // this is temp in ppl.c
	    low_ld_integral =
	      low_integral * modelOptions.LDprior * modelOptions.thetaWeight;
	    ppld =
	      low_ld_integral / (low_ld_integral +
				 (1 -
				  modelOptions.LDprior) *
				 modelOptions.thetaWeight *
				 low_theta_integral + (1 -
						       modelOptions.
						       thetaWeight)
				 * high_theta_integral);
	    fprintf (fpPPL, "%6.4f %6.4f ", ldppl, ppld);
	    // fprintf (fpSeok, "ppl= %6.4f  ldppl= %6.4f  ppld= %6.4f\n", ppl,ldppl, ppld);     
	  }

	  fprintf (fpPPL, " %8.4f %6.4f %6.4f %6.4f %6.4f ",
		   log10 (maximum_function_value), maxima_x[0],
		   maxima_x[1], maxima_x[3], maxima_x[2]);
	  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	    fprintf (fpPPL, "%6.4f %6.4f %6.4f ",
		     maxima_x[liabIdx * 3 + 4],
		     maxima_x[liabIdx * 3 + 5], maxima_x[liabIdx * 3 + 6]);
	  }
	  fprintf (fpPPL, "\n");
	  fflush (fpPPL);
	  //fprintf(fpSeok,"low_integral =%f high integral =%f low thetat=%f high theta=%f low ld =%f\n",low_integral, high_integral,low_theta_integral,high_theta_integral,low_ld_integral);



	  /* only loop marker allele frequencies when doing LD */
	  if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
	    break;
	  /* we can only do SNPs when looping over marker allele frequency */
	  if (pLocus2->numOriginalAllele > 2)
	    break;
	}
	/* end of marker allele frequency looping */
	prevNumDPrime = pLambdaCell->ndprime;
	/* need to clear polynomial */

	if (modelOptions.polynomial == TRUE && modelType.ccFlag == 0) {
	  /* under case ctrl we don't clear up the polynomial */
	  pedigreeSetPolynomialClearance (&pedigreeSet);
	}


	if (modelOptions.markerAnalysis == ADJACENTMARKER)
	  loc2 = originalLocusList.numLocus;
      }				/* end of looping second locus - loc2 */
      /* if we are doing trait marker, then we are done */
      /* Used to read: modelOptions.markerToMarker != TRUE which
         is the same as markerAnalysis == FALSE as long as the old
         markerToMarker and adjacentMarker flags were truly
         orthogonal. Otherwise, it should be markerAnalysis !=
         ADJACENTMARKER. */
      if (modelOptions.markerAnalysis == FALSE)
	loc1 = originalLocusList.numLocus;
    }				/* end of looping first locus - loc1 */
    /* free two point result storage */
    //free_tp_result_storage (prevNumDPrime);
  }
  /* end of two point */
  else {
    /* multipoint */
    /* marker set locus list for each position */
    markerLocusList.maxNumLocus = modelType.numMarkers;
    markerLocusList.numLocus = modelType.numMarkers;
    markerLocusList.traitOrigLocus = -1;
    markerLocusList.traitLocusIndex = -1;
    markerLocusList.pLocusIndex =
      (int *) calloc (markerLocusList.maxNumLocus, sizeof (int));
    for (k = 0; k < 3; k++) {
      markerLocusList.pPrevLocusDistance[k] =
	(double *) calloc (markerLocusList.maxNumLocus, sizeof (double));
      markerLocusList.pNextLocusDistance[k] =
	(double *) calloc (markerLocusList.maxNumLocus, sizeof (double));
    }

    /* assuming we always have trait in the analysis - this may not be true 
     * need to add code to process marker to marker analysis under multipoin
     */
    savedLocusList.numLocus = modelType.numMarkers + 1;
    savedLocusList.maxNumLocus = modelType.numMarkers + 1;
    savedLocusList.pLocusIndex =
      (int *) calloc (savedLocusList.maxNumLocus, sizeof (int));
    for (k = 0; k < 3; k++) {
      savedLocusList.pPrevLocusDistance[k] =
	(double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
      savedLocusList.pNextLocusDistance[k] =
	(double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
    }

    /* calculate the trait likelihood independent of the trait position */
    traitLocusList.numLocus = 1;
    traitLocusList.maxNumLocus = 1;
    traitLocusList.traitLocusIndex = 0;
    traitLocusList.traitOrigLocus = traitLocus;
    traitLocusList.pLocusIndex =
      (int *) calloc (traitLocusList.maxNumLocus, sizeof (int));
    traitLocusList.pLocusIndex[0] = 0;
    for (k = 0; k < 3; k++) {
      traitLocusList.pPrevLocusDistance[k] =
	(double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
      traitLocusList.pNextLocusDistance[k] =
	(double *) calloc (savedLocusList.maxNumLocus, sizeof (double));

      traitLocusList.pPrevLocusDistance[k][0] = -1;
      traitLocusList.pNextLocusDistance[k][0] = -1;
    }
    /* populate the trait xmission matrix */
    locusList = &traitLocusList;
    xmissionMatrix = traitMatrix;
    status = populate_xmission_matrix (traitMatrix, 1, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1,	/* last he locus */
				       -1,	/* last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */

    if (modelOptions.polynomial == TRUE) {
      holdAllPolys ();
      fprintf (stderr,
	       "holdAllPolys from further population of transmission matrix\n");
    }

    /* for trait likelihood */
    fprintf (stderr, "MP start time: %f\n", (double) time0 / CLOCKS_PER_SEC);
    locusList = &traitLocusList;
    xmissionMatrix = traitMatrix;
    if (pTrait->type == DICHOTOMOUS) {

      /*call compute_likelihood with dummy numbers to build polynomials */
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	pTrait->penetrance[2][liabIdx][0][0] = 0.7;
	pTrait->penetrance[2][liabIdx][0][1] = 0.5;
	pTrait->penetrance[2][liabIdx][1][0] = 0.5;
	pTrait->penetrance[2][liabIdx][1][1] = 0.3;
	pTrait->penetrance[1][liabIdx][0][0] = 1 - 0.7;
	pTrait->penetrance[1][liabIdx][0][1] = 1 - 0.5;
	pTrait->penetrance[1][liabIdx][1][0] = 1 - 0.5;
	pTrait->penetrance[1][liabIdx][1][1] = 1 - 0.3;
      }

      if (modelOptions.polynomial == TRUE);
      else
	/* only need to update trait locus */
	update_penetrance (&pedigreeSet, traitLocus);

      pLocus->pAlleleFrequency[0] = 0.5;
      pLocus->pAlleleFrequency[1] = 1 - 0.5;

      if (modelOptions.polynomial == TRUE);
      else
	update_locus (&pedigreeSet, traitLocus);
      /* get the likelihood for the trait */
      KLOG (LOGLIKELIHOOD, LOGDEBUG, "Trait Likelihood\n");
      compute_likelihood (&pedigreeSet);	/* This builds polynomials with dummy numbers */

    } else {			// QT
      pLocus->pAlleleFrequency[0] = 0.5;
      pLocus->pAlleleFrequency[1] = 1 - 0.5;
      update_locus (&pedigreeSet, traitLocus);

      /*call compute_likelihood with dummy numbers to build polynomials */
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	pTrait->means[liabIdx][0][0] = 2.0;
	pTrait->means[liabIdx][0][1] = 1.0;
	pTrait->means[liabIdx][1][0] = 1.0;
	pTrait->means[liabIdx][1][1] = 0.0;
	pTrait->stddev[liabIdx][0][0] = 1.0;
	pTrait->stddev[liabIdx][0][1] = 1.0;
	pTrait->stddev[liabIdx][1][0] = 1.0;
	pTrait->stddev[liabIdx][1][1] = 1.0;

	/* threshold for QT */
	pTrait->cutoffValue[liabIdx] = 0.5;

      }				/* liability class Index */
      if (modelOptions.polynomial == TRUE);
      else
	update_penetrance (&pedigreeSet, traitLocus);
      KLOG (LOGLIKELIHOOD, LOGDEBUG, "Trait Likelihood\n");
      compute_likelihood (&pedigreeSet);

    }
    time2 = clock ();
    fprintf (stderr, "MP done trait: %f\n", (double) time2 / CLOCKS_PER_SEC);

    /* coply the polynomials built from above to traitLikelihoodPolynomials */
    if (modelOptions.polynomial == TRUE) {
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
	/* save the likelihood at trait */
	pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
	pPedigree->traitLikelihoodPolynomial =
	  pPedigree->likelihoodPolynomial;
	pPedigree->traitLikelihoodPolyList = pPedigree->likelihoodPolyList;
	pPedigree->likelihoodPolyList = NULL;
	pPedigree->likelihoodPolynomial = NULL;

	//fprintf(stderr,"Building traitPoly pedIdx =%d Null likelihood = %20.15f\n",pedIdx, pPedigree->likelihood);
	//fprintf(stderr,"pedIdx %d eType= %d\n", pedIdx, ((pPedigree->traitLikelihoodPolyList)->pList[0])->eType);
      }

    }



    /* get the trait locations we need to evaluate at */
    numPositions = modelRange.ntloc;
    mp_result = (SUMMARY_STAT *) calloc (numPositions, sizeof (SUMMARY_STAT));
    /* Need to output the results */
    if (modelType.trait == DICHOTOMOUS) {
      fprintf (fpHet,
	       "           pos       PPL         BR       error    num    markerList   maximum   alpha    gf     DD       Dd      dd\n");
    } else {
      fprintf (fpHet,
	       "           pos       PPL         BR       error    num    markerList   maximum   alpha    gf  meanDD  meanDd   meandd");
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	fprintf (fpHet, " %5s %5s %5s ", "SD_DD", "SD_Dd", "SD_dd");
      }
      if (modelType.trait == CT) {
	fprintf (fpHet, "  %5s", "t");
      }
      fprintf (fpHet, "\n");
    }
    fflush (fpHet);


    prevFirstMarker = -1;
    prevLastMarker = -1;
    prevTraitInd = -1;
    leftMarker = -1;
    for (posIdx = 0; posIdx < numPositions; posIdx++) {
      /* positions listed are sex average positions */
      traitPos = modelRange.tloc[posIdx];
      /* set the sex average position first 
       * the sex specific positions will be updated once markers are selected
       * as interpolation might be needed
       */
      pTraitLocus->mapPosition[0] = traitPos;
      pTraitLocus->mapPosition[1] = traitPos;
      pTraitLocus->mapPosition[2] = traitPos;
      /* initialize the locusList */
      locusList = &savedLocusList;
      memset (locusList->pLocusIndex, 0,
	      sizeof (int) * locusList->maxNumLocus);
      for (k = 0; k < 3; k++) {
	memset (&locusList->pPrevLocusDistance[k][0], 0,
		sizeof (double) * locusList->maxNumLocus);
	memset (&locusList->pNextLocusDistance[k][0], 0,
		sizeof (double) * locusList->maxNumLocus);
      }
      locusList->numLocus = 1;
      locusList->pLocusIndex[0] = traitLocus;
      for (k = 0; k < 3; k++) {
	locusList->pPrevLocusDistance[k][0] = -1;
	locusList->pNextLocusDistance[k][0] = -1;
      }
      /* select markers to be used for the multipoint analysis */
      add_markers_to_locuslist (locusList, modelType.numMarkers,
				&leftMarker, 0,
				originalLocusList.numLocus - 1, traitPos, 0);

      /* store the markers used */
      mp_result[posIdx].pMarkers =
	(int *) calloc (modelType.numMarkers, sizeof (int));
      k = 0;			/* marker index */
      for (i = 0; i < locusList->numLocus; i++) {
	j = locusList->pLocusIndex[i];
	if (originalLocusList.ppLocusList[j]->locusType == LOCUS_TYPE_MARKER) {
	  mp_result[posIdx].pMarkers[k] = j;
	  k++;
	} else {
	  mp_result[posIdx].trait = i;
	  traitIndex = i;
	}
      }
      locusList->traitLocusIndex = traitIndex;
      locusList->traitOrigLocus = traitLocus;
      markerSetChanged = FALSE;
      if (prevFirstMarker != mp_result[posIdx].pMarkers[0]
	  || prevLastMarker !=
	  mp_result[posIdx].pMarkers[modelType.numMarkers - 1]) {
	/* marker set has changed */
	markerSetChanged = TRUE;
	markerLocusList.pLocusIndex[0] = mp_result[posIdx].pMarkers[0];
	prevPos = get_map_position (markerLocusList.pLocusIndex[0]);
	/* set up marker set locus list */
	for (k = 1; k < modelType.numMarkers; k++) {
	  markerLocusList.pLocusIndex[k] = mp_result[posIdx].pMarkers[k];
	  currPos = get_map_position (markerLocusList.pLocusIndex[k]);
	  for (j = 0; j < 3; j++) {
	    markerLocusList.pPrevLocusDistance[j][k] =
	      markerLocusList.pNextLocusDistance[j][k - 1] =
	      cm_to_recombination_fraction (currPos[j] - prevPos[j],
					    map.mapFunction);
	  }
	  if (modelOptions.mapFlag == SA) {
	    for (j = 1; j <= 2; j++) {
	      markerLocusList.pPrevLocusDistance[j][k] =
		markerLocusList.pNextLocusDistance[j][k - 1] =
		markerLocusList.pPrevLocusDistance[0][k];
	    }
	  }
	  prevPos = currPos;
	}			/* end of loop over the markers to set up locus list */

	/* calculate likelihood for the marker set */
	locusList = &markerLocusList;
	xmissionMatrix = markerMatrix;
#if 1
	if (modelOptions.polynomial == TRUE) {
	  pedigreeSetPolynomialClearance (&pedigreeSet);
	}
#endif

	/* save the polynomial flag */
	polynomialFlag = modelOptions.polynomial;
#if 0
	if (polynomialFlag == TRUE) {
	  for (k = 0; k < 3; k++) {
	    initialProb[k] = 1.0;
	    initialProb2[k] = 1.0;
	    initialProbAddr[k] = &initialProb[k];
	    initialProbAddr2[k] = &initialProb2[k];
	    initialHetProbAddr[k] = NULL;
	  }
	}

	modelOptions.polynomial = FALSE;
#endif
	/* populate the matrix */
	status = populate_xmission_matrix (markerMatrix, markerLocusList.numLocus, initialProbAddr,	/* probability */
					   initialProbAddr2,	/* probability */
					   initialHetProbAddr, 0,	/* cell index */
					   -1,	/* last he locus */
					   -1,	/* last het pattern (P-1 or M-2) */
					   0);	/* current locus - start with 0 */

	if (modelOptions.polynomial == TRUE)
	  freePolys ();

	print_xmission_matrix (markerMatrix, markerLocusList.numLocus, 0, 0,
			       tmpID);

	KLOG (LOGLIKELIHOOD, LOGDEBUG, "Marker Likelihood\n");
	compute_likelihood (&pedigreeSet);
	modelOptions.polynomial = polynomialFlag;

	/* save the results for marker likelihood */
	for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
	  /* save the likelihood at null */
	  pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
	  pPedigree->markerLikelihood = pPedigree->likelihood;
	}
	pedigreeSet.markerLikelihood = pedigreeSet.likelihood;
	pedigreeSet.log10MarkerLikelihood = pedigreeSet.log10Likelihood;
      }
      /* end of marker set change */
      prevFirstMarker = mp_result[posIdx].pMarkers[0];
      prevLastMarker = mp_result[posIdx].pMarkers[modelType.numMarkers - 1];
      if (markerSetChanged || prevTraitInd != mp_result[posIdx].trait)
	locusListChanged = TRUE;
      else
	locusListChanged = FALSE;
      prevTraitInd = mp_result[posIdx].trait;

      locusList = &savedLocusList;
      xmissionMatrix = altMatrix;
      /* interpolate trait postion for sex specific analysis */
      if (modelOptions.mapFlag == SS) {
	if (traitIndex == 0) {
	  marker1Pos = get_map_position (locusList->pLocusIndex[1]);
	  /* trait is the first one in the list */
	  if (traitPos < ERROR_MARGIN && traitPos > -ERROR_MARGIN) {
	    /* trait is at 0 */
	    pTraitLocus->mapPosition[MAP_MALE] =
	      pTraitLocus->mapPosition[MAP_FEMALE] = 0;
	  } else {
	    /* get the relative position on the sex average map */
	    relativePos = traitPos / marker1Pos[0];
	    pTraitLocus->mapPosition[MAP_MALE] =
	      relativePos * marker1Pos[MAP_MALE];
	    pTraitLocus->mapPosition[MAP_FEMALE] =
	      relativePos * marker1Pos[MAP_FEMALE];
	  }
	  /* update the inter locus distance - sex averaged already done before */
	  for (k = 1; k < 3; k++) {
	    locusList->pNextLocusDistance[k][0] =
	      locusList->pPrevLocusDistance[k][1] =
	      cm_to_recombination_fraction (marker1Pos[k] -
					    pTraitLocus->
					    mapPosition[k], map.mapFunction);
	  }
	} else if (traitIndex == modelType.numMarkers) {
	  /* trait is the last one in the list */
	  marker1Pos =
	    get_map_position (locusList->
			      pLocusIndex[modelType.numMarkers - 2]);
	  marker2Pos =
	    get_map_position (locusList->
			      pLocusIndex[modelType.numMarkers - 1]);
	  /* get the relative position on the sex average map */
	  dist = marker2Pos[0] - marker1Pos[0];
	  if (dist > ERROR_MARGIN) {
	    relativePos = (traitPos - marker2Pos[0]) / dist;
	    pTraitLocus->mapPosition[MAP_MALE] =
	      relativePos * (marker2Pos[MAP_MALE] -
			     marker1Pos[MAP_MALE]) + marker2Pos[MAP_MALE];
	    pTraitLocus->mapPosition[MAP_FEMALE] =
	      relativePos * (marker2Pos[MAP_FEMALE] -
			     marker1Pos[MAP_FEMALE]) + marker2Pos[MAP_FEMALE];
	  } else {
	    pTraitLocus->mapPosition[MAP_MALE] =
	      traitPos - marker2Pos[0] + marker2Pos[MAP_MALE];
	    pTraitLocus->mapPosition[MAP_FEMALE] =
	      traitPos - marker2Pos[0] + marker2Pos[MAP_FEMALE];
	  }

	  /* update the inter locus distance - sex averaged already done before */
	  for (k = 1; k <= 2; k++) {
	    locusList->pNextLocusDistance[k][traitIndex - 1] =
	      locusList->pPrevLocusDistance[k][traitIndex] =
	      cm_to_recombination_fraction (pTraitLocus->
					    mapPosition[k] -
					    marker2Pos[k], map.mapFunction);
	  }
	} else {
	  /* trait is in between two markers */
	  marker1Pos =
	    get_map_position (locusList->pLocusIndex[traitIndex - 1]);
	  marker2Pos =
	    get_map_position (locusList->pLocusIndex[traitIndex + 1]);
	  /* get the relative position on the sex average map */
	  dist = marker2Pos[0] - marker1Pos[0];
	  if (dist > ERROR_MARGIN) {
	    relativePos = (traitPos - marker1Pos[0]) / dist;
	    pTraitLocus->mapPosition[MAP_MALE] =
	      relativePos * (marker2Pos[MAP_MALE] -
			     marker1Pos[MAP_MALE]) + marker1Pos[MAP_MALE];
	    pTraitLocus->mapPosition[MAP_FEMALE] =
	      relativePos * (marker2Pos[MAP_FEMALE] -
			     marker1Pos[MAP_FEMALE]) + marker1Pos[MAP_FEMALE];
	  } else {
	    pTraitLocus->mapPosition[MAP_MALE] = marker1Pos[MAP_MALE];
	    pTraitLocus->mapPosition[MAP_FEMALE] = marker1Pos[MAP_FEMALE];
	  }
	  /* update the inter locus distance - sex averaged already done before */
	  for (k = 1; k < 3; k++) {
	    locusList->pNextLocusDistance[k][traitIndex - 1] =
	      locusList->pPrevLocusDistance[k][traitIndex] =
	      cm_to_recombination_fraction (pTraitLocus->
					    mapPosition[k] -
					    marker1Pos[k], map.mapFunction);
	    locusList->pNextLocusDistance[k][traitIndex] =
	      locusList->pPrevLocusDistance[k][traitIndex + 1] =
	      cm_to_recombination_fraction (marker2Pos[k] -
					    pTraitLocus->
					    mapPosition[k], map.mapFunction);

	  }
	}

      }


      /*end of SS model */
      /* the locus list has been built, go on to the analysis 
       * multipoint DT */
      if (markerSetChanged || locusListChanged) {
	if (modelOptions.polynomial == TRUE) {
	  pedigreeSetPolynomialClearance (&pedigreeSet);
	  /* populate the matrix */
	  status = populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr,	/* probability */
					     initialProbAddr2,	/* probability */
					     initialHetProbAddr, 0,	/* cell index */
					     -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
					     0);	/* current locus - start with 0 */
	  print_xmission_matrix (altMatrix, savedLocusList.numLocus, 0, 0,
				 tmpID);
	  if (modelOptions.polynomial == TRUE)
	    freePolys ();
	}
      }

      if (modelOptions.polynomial != TRUE);
	/* populate the matrix */
	status = populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr,	/* probability */
					   initialProbAddr2,	/* probability */
					   initialHetProbAddr, 0,	/* cell index */
					   -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
					   0);	/* current locus - start with 0 */

      /* multipoint DT */
      integral = 0.0;
      abserr = 0.0;


      num_out_constraint = 0;
      num_eval = kelvin_dcuhre_integrate (&integral, &abserr);
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	integral *= 6;
	abserr *= 6;
      }
      integral /= volume_region;
      abserr /= volume_region;

      time2 = clock ();
      fprintf (stderr, "MP done ALT on pos %d: %f\n", posIdx,
	       (double) time2 / CLOCKS_PER_SEC);
      fprintf (stderr, "BR is %15.13f  and error is %15.13f \n", integral,
	       abserr);


      /* calculate imputed PPL and print the results */
      if (integral > 0.214)
	ppl =
	  (integral * integral) / (-5.77 + 54 * integral +
				   integral * integral);
      else
	ppl = 0;

      fprintf (fpHet, "\t %f  %6.4f %12.8f %12.8f %d  ", traitPos, ppl,
	       integral, abserr, num_eval);
      fprintf (stderr, "\t %f  %6.4f %12.8f %12.8f %d  ", traitPos, ppl,
	       integral, abserr, num_eval);
      /* print out markers used for this position */
      fprintf (fpHet, "(%d", mp_result[posIdx].pMarkers[0]);
      for (k = 1; k < modelType.numMarkers; k++) {
	fprintf (fpHet, ",%d", mp_result[posIdx].pMarkers[k]);
      }
      fprintf (fpHet, ") %f %f %f ", log10 (localmax_value), localmax_x[1],
	       localmax_x[0]);
      if (modelType.trait == DICHOTOMOUS) {
	for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	  fprintf (fpHet, " %f %f %f ", localmax_x[3 * liabIdx + 2],
		   localmax_x[3 * liabIdx + 3], localmax_x[3 * liabIdx + 4]);
	}
	fprintf (fpHet, "\n");
      } else {			//QT
	j = 2;
	for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	  fprintf (fpHet, "%6.4f %6.4f %6.4f ", localmax_x[j],
		   localmax_x[j + 1], localmax_x[j + 2]);
	  j += 3;
	  if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	    fprintf (fpHet, "%6.4f %6.4f %6.4f ", localmax_x[j],
		     localmax_x[j + 1], localmax_x[j + 2]);
	    j += 3;
	  }
	  if (modelType.trait == CT) {
	    fprintf (fpHet, "%6.4f", localmax_x[j++]);
	  }
	}
	fprintf (fpHet, "  %d\n", num_out_constraint);
      }				/* End of writing max */
      fflush (fpHet);

    }				/* end of walking down the chromosome */
  }				/* end of multipoint */


  time2 = clock ();

  fprintf (stderr, "Computation time:  %fs  %fs \n",
	   (double) (time1 - time0) / CLOCKS_PER_SEC,
	   (double) (time2 - time1) / CLOCKS_PER_SEC);

  if (modelOptions.polynomial == TRUE) {
//   dismantle();
  }

  set_removeGenotypeFlag (TRUE);

  if (modelType.type == TP) {
    if (tp_result != NULL) {
      /* two point */
      for (i = 0; i < pLambdaCell->ndprime; i++) {
	free (tp_result[i]);
      }
      free (tp_result);
    }
  } else {
    /* multipoint */
    for (posIdx = 0; posIdx < numPositions; posIdx++) {
      free (mp_result[posIdx].pMarkers);
    }
    free (mp_result);
  }

  /* free parental pair work space */
  free_parental_pair_workspace (&parentalPairSpace, modelType.numMarkers + 1);

  /* free transmission probability matrix */
  free (altMatrix);
  free (nullMatrix);

  if (modelType.type == TP)
    free_LD_loci (pLDLoci);

  if (modelOptions.polynomial == TRUE) {
    pedigreeSetPolynomialClearance (&pedigreeSet);
  }
  free_likelihood_storage (&pedigreeSet);
  free_likelihood_space (&pedigreeSet);
  free_pedigree_set (&pedigreeSet);
  free_sub_locus_list (&traitLocusList);
  free_sub_locus_list (&markerLocusList);
  free_sub_locus_list (&savedLocusList);
  free (modelOptions.sUnknownPersonID);
  final_cleanup ();


  fprintf (stderr, "Computation time:  %fs  %fs \n",
	   (double) (time1 - time0) / CLOCKS_PER_SEC,
	   (double) (time2 - time1) / CLOCKS_PER_SEC);

  /* Final dump and clean-up for performance. */
  swStop (overallSW);
  if (modelOptions.polynomial == TRUE)
    polyStatistics ("End of run");
  else
    swDump (overallSW);
#ifdef DMTRACK
  swLogPeaks ("End of run");
  swDumpHeldTotals ();
  swDumpSources ();
  //  swDumpCrossModuleChunks ();
#endif
  swLogMsg ("Finished run");

  /* close file pointers */
  if (modelType.type == TP) {
    fclose (fpPPL);
  }
  fclose (fpHet);
  //  fclose (fpHomo);

  return 0;
}

void
logStatistics(PedigreeSet *pSet, int posIdx) {
  int pedIdx, i;
  NuclearFamily *pNucFam;
  Pedigree *pPedigree;
  int l = 0, nf = 0, pg = 0, sg = 0;
  for (pedIdx = 0; pedIdx < pSet->numPedigree; pedIdx++) {
    pPedigree = pSet->ppPedigreeSet[pedIdx];
    l += pPedigree->numLoop;
    nf += pPedigree->numNuclearFamily;
    for (i = 0; i < pPedigree->numNuclearFamily; i++) {
      pNucFam = pPedigree->ppNuclearFamilyList[i];
      pg += pNucFam->totalNumPairGroups;
      sg += pNucFam->totalNumSimilarPairs;
    }
  }
  sprintf(messageBuffer, "At %d: p:%d, l:%d, nf:%d, pg:%d, sg:%d, n:%d",
	  posIdx, pSet->numPedigree, l, nf, pg, sg, nodeId);
  swLogMsg(messageBuffer);
}

void
print_dryrun_stat (PedigreeSet * pSet, double pos)
{
  int pedIdx;
  long subTotalPairGroups, subTotalSimilarPairs;
  long totalPairGroups, totalSimilarPairs;
  NuclearFamily *pNucFam;
  Pedigree *pPedigree;
  int i;

  totalPairGroups = 0;
  totalSimilarPairs = 0;
  for (pedIdx = 0; pedIdx < pSet->numPedigree; pedIdx++) {
    /* save the likelihood at null */
    pPedigree = pSet->ppPedigreeSet[pedIdx];
    fprintf (stderr, "Ped %s(%d) has %d loops, %d nuclear families.\n",
	     pPedigree->sPedigreeID, pedIdx,
	     pPedigree->numLoop, pPedigree->numNuclearFamily);
    subTotalPairGroups = 0;
    subTotalSimilarPairs = 0;
    for (i = 0; i < pPedigree->numNuclearFamily; i++) {
      pNucFam = pPedigree->ppNuclearFamilyList[i];
      fprintf (stderr,
	       "    Nuc %d w/ proband %s(%s) has %ld unique pp groups, %ld similar pp, total %ld.\n",
	       i, pNucFam->pProband->sID,
	       pNucFam->childProbandFlag ? "child" : "parent",
	       pNucFam->totalNumPairGroups, pNucFam->totalNumSimilarPairs,
	       pNucFam->totalNumPairGroups + pNucFam->totalNumSimilarPairs);
      subTotalPairGroups += pNucFam->totalNumPairGroups;
      subTotalSimilarPairs += pNucFam->totalNumSimilarPairs;
    }
    fprintf (stderr,
	     "    Ped has total %ld unique pp groups, %ld similar pp, total %ld.\n",
	     subTotalPairGroups, subTotalSimilarPairs,
	     subTotalPairGroups + subTotalSimilarPairs);
    totalPairGroups += subTotalPairGroups;
    totalSimilarPairs += subTotalSimilarPairs;
  }
  fprintf (stderr,
	   "POS %f has %ld unique pp groups, %ld similar pp, total %ld.\n",
	   pos, totalPairGroups, totalSimilarPairs,
	   totalPairGroups + totalSimilarPairs);
}


void
test_darray (double **tpl)
{
  int i, j;
  double *gene_tpl;

  for (i = 0; i < 6; i++) {
    gene_tpl = tpl[i];

    for (j = 0; j < 275; j++) {
      if (gene_tpl[j] > 1.0e40 || gene_tpl[j] < 1.0e-40) {
	printf ("gfId= %d penId%d  likelihood = %G\n", i, j, gene_tpl[j]);
	break;
      }
    }
  }


}


int
kelvin_dcuhre_integrate (double *integral, double *abserr)
{
  /* Local variables */
  //double a[15], b[15];
  int dim, return_val;
  dcuhre_state init_state;

  //extern /* Subroutine */ int ftest_();  

  localmax_value = 0.0;



  if (modelType.trait == DICHOTOMOUS) {

    dim = 1 + 3 * modelRange.nlclass;
    s = &init_state;
    initialize_state (s, xl, xu, dim);
    s->verbose = 1;
    s->nlclass = modelRange.nlclass;

    if (modelType.type == TP) {
      s->funsub = (U_fp) compute_hlod_2p_dt;
      s->mType = TP_DT;
    } else {
      s->funsub = (U_fp) compute_hlod_mp_dt;
      s->mType = MP_DT;
    }
  } else {			/*  QT or combined */

    dim = total_dim - 1;	// alpha
    if (modelType.type == TP) {
      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
	dim -= 2;		// theta and dprime
      } else {
	dim -= 1;		//theta
      }
    }

    s = &init_state;
    initialize_state (s, xl, xu, dim);
    s->verbose = 1;
    s->nlclass = modelRange.nlclass;

    if (modelType.type == TP) {
      s->funsub = (U_fp) compute_hlod_2p_qt;
      s->mType = TP_DT;
    } else {
      s->funsub = (U_fp) compute_hlod_mp_qt;
      s->mType = MP_DT;
    }
  }


  fprintf (stderr, "Starting DCUHRE with dim=%d\n", dim);
  return_val = dcuhre_ (s);
  if (return_val > 0 && return_val < 20) {
    fprintf (stderr, "Ending program with error! ifail =%d \n", s->ifail);
  }

  fprintf (stderr,
	   "Final result =%15.10f  with error =%15.10f and neval = %d\n",
	   s->result, s->error, s->total_neval);
  fprintf (stderr, "End of DCUHRE with ifail =%d\n", s->ifail);


  *integral = s->result;
  *abserr = s->error;
  return return_val;

}


void
compute_hlod_mp_qt (double x[], double *f)
{

  int k, j;
  int pedIdx, liabIdx = 0, status;
  double constraint = 0.0;
  double mean_DD;
  double mean_Dd;
  double mean_dd;
  double SD_DD = 0.0;
  double SD_Dd = 0.0;
  double SD_dd = 0.0;
  double gfreq;
  double alphaV;

//  double  theta;
  double threshold = 0.0;
  double log10_likelihood_null, log10_likelihood_alternative,
    log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;
  double log10Likelihood;

  /* for null likelihood calculation */
  double product_likelihood = 1;	/* product of the likelihoods for all the pedigrees */
  double sum_log_likelihood = 0;	/* sum of the log10(likelihood) for all the pedigrees */

  Pedigree *pPedigree;

  //  checkpt ();

  int origLocus = locusList->pLocusIndex[0];

  if (locusList->numLocus > 1)
    origLocus = locusList->pLocusIndex[1];

  gfreq = x[0];
  pLocus->pAlleleFrequency[0] = gfreq;
  pLocus->pAlleleFrequency[1] = 1 - gfreq;

  update_locus (&pedigreeSet, traitLocus);

  //  checkpt ();
  if (print_point_flag)
    fprintf (fphlod, "dp=%f th=%f gf=%f ", fixed_dprime, fixed_theta, gfreq);


  j = 1;			// j=0 for gfrequency

  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    mean_DD = x[j];
    mean_Dd =
      (x[j + 1] - xl[j + 1]) * (x[j] - xl[j]) / (xu[j] - xl[j]) + xl[j + 1];
    mean_dd =
      (x[j + 2] - xl[j + 2]) * (x[j + 1] - xl[j + 1]) / (xu[j + 1] -
							 xl[j + 1]) * (x[j] -
								       xl[j])
      / (xu[j] - xl[j]) + xl[j + 2];
    if (print_point_flag)
      fprintf (fphlod, "pe %f %f %f ", mean_DD, mean_Dd, mean_dd);
    j += 3;
    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      SD_DD = x[j];		//      modelRange.param[liabIdx][0][0][paramIdx];
      SD_Dd = x[j + 1];		//        modelRange.param[liabIdx][1][0][paramIdx];
      SD_dd = x[j + 2];		//        modelRange.param[liabIdx][2][0][paramIdx];
      if (print_point_flag)
	fprintf (fphlod, "sd %f %f %f ", SD_DD, SD_Dd, SD_dd);
      j += 3;
    }
    /* threshold for QT */
    if (modelType.trait == CT) {
      threshold = x[j];		// modelRange.tthresh[liabIdx][thresholdIdx];
      if (print_point_flag)
	fprintf (fphlod, "thd=%f ", threshold);
      j++;
    }

    /* check against the hard coded constraint */
    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      constraint =
	(1.0 - gfreq) * (1.0 - gfreq) * mean_dd * SD_dd + 2 * gfreq * (1.0 -
								       gfreq)
	* mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
      /* fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
         constraint, gfreq, mean_DD, SD_DD, mean_Dd, SD_DD, mean_dd, SD_dd);
       */
      if (constraint >= 3.0 || constraint <= -3.0) {
	// fprintf(stderr,"Constraint is %f \n", constraint);
	num_out_constraint++;
      }
    }

    pTrait->means[liabIdx][0][0] = mean_DD;
    pTrait->means[liabIdx][0][1] = mean_Dd;
    pTrait->means[liabIdx][1][0] = mean_Dd;
    pTrait->means[liabIdx][1][1] = mean_dd;
    pTrait->stddev[liabIdx][0][0] = SD_DD;
    pTrait->stddev[liabIdx][0][1] = SD_Dd;
    pTrait->stddev[liabIdx][1][0] = SD_Dd;
    pTrait->stddev[liabIdx][1][1] = SD_dd;

    /* threshold for QT */
    pTrait->cutoffValue[liabIdx] = threshold;

  }				/* liability class Index */
  //  checkpt ();
  if (modelOptions.polynomial == TRUE);
  else
    /* only need to update trait locus */
    update_penetrance (&pedigreeSet, traitLocus);
  //  update_penetrance (&pedigreeSet, traitLocus);

  if (print_point_flag)
    fprintf (fphlod, "\n");

  /*This is a temporary checking. */
  if (j != s->ndim) {
    fprintf (stderr, "j=%d  while dim for BR is %d\n", j, s->ndim);
    exit (0);
  }

  /* for trait likelihood */
  locusList = &traitLocusList;
  xmissionMatrix = traitMatrix;

  //fprintf(stderr,"Null likelihood computation is starting for %d pedigrees \n",pedigreeSet.numPedigree);
  /* compute the null likelihood with   */
  pedigreeSet.likelihood = 1;
  pedigreeSet.log10Likelihood = 0;


  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {

    /* save the likelihood at null */
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];

    //fprintf(stderr,"pedIdx =%d ", pedIdx); 
    if (modelOptions.polynomial == TRUE) {
      KASSERT (pPedigree->traitLikelihoodPolynomial != NULL, "Error in  \n");
      /* evaluate likelihood */
      // fprintf (stderr, "evaluaing poly %d pedigree\n", pedIdx);

      evaluatePoly (pPedigree->traitLikelihoodPolynomial,
		    pPedigree->traitLikelihoodPolyList,
		    &pPedigree->likelihood);

      // fprintf (stderr, " is done %f with %d pedigrees\n",
      //               pPedigree->likelihood, pedigreeSet.numPedigree);
    } else {
      initialize_multi_locus_genotype (pPedigree);
      status = compute_pedigree_likelihood (pPedigree);
    }

    /*pPedigree->likelihood is now computed and now check it */
    if (pPedigree->likelihood == 0.0) {
      KLOG (LOGLIKELIHOOD, LOGWARNING,
	    "Pedigree %s has likelihood of 0 or too small.\n",
	    pPedigree->sPedigreeID);
      fprintf (stderr, "Pedigree %s has likelihood of 0 or too small.\n",
	       pPedigree->sPedigreeID);
      product_likelihood = 0.0;
      sum_log_likelihood = -9999.99;
      break;
    } else if (pPedigree->likelihood < 0.0) {
      KASSERT (pPedigree->likelihood >= 0.0,
	       "Pedigree %s with NEGATIVE likelihood - This is CRAZY!!!.\n",
	       pPedigree->sPedigreeID);
      product_likelihood = 0.0;
      sum_log_likelihood = -9999.99;
      break;
    } else {
      if (pPedigree->pCount[origLocus] == 1) {
	product_likelihood *= pPedigree->likelihood;
	log10Likelihood = log10 (pPedigree->likelihood);
      } else {
	product_likelihood *=
	  pow (pPedigree->likelihood, pPedigree->pCount[origLocus]);
	log10Likelihood =
	  log10 (pPedigree->likelihood) * pPedigree->pCount[origLocus];
      }
      sum_log_likelihood += log10Likelihood;
    }
    pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
    // fprintf (stderr,
    //     "null likelihood pedIdx=%d is done %20.18f with product =%20.16f\n",
    //     pedIdx, pPedigree->likelihood, product_likelihood);
  }

  pedigreeSet.likelihood = product_likelihood;
  pedigreeSet.log10Likelihood = sum_log_likelihood;
  log10_likelihood_null = pedigreeSet.log10Likelihood;
  KLOG (LOGLIKELIHOOD, LOGDEBUG, "Sum of log Likelihood is: %e\n",
	sum_log_likelihood);

  //  checkpt ();

  /* This is for alternative likelihood */
  locusList = &savedLocusList;
  xmissionMatrix = altMatrix;
  if (modelOptions.polynomial == TRUE);
  else
    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */

  KLOG (LOGLIKELIHOOD, LOGDEBUG, "Likelihood\n");
  compute_likelihood (&pedigreeSet);

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;
  if (isnan (log10_likelihood_alternative))
    fprintf (stderr, "ALT likelihood is NAN.\n");
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
  } else {
    log10_likelihood_ratio =
      log10_likelihood_alternative - log10_likelihood_null -
      pedigreeSet.log10MarkerLikelihood;
  }
  /* check for overflow problem !!! */
  if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
    likelihood_ratio = DBL_MAX;
  } else {
    /* check for underflow problem too !!! */
    if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
      likelihood_ratio = 0;
    } else {
      likelihood_ratio = pow (10.0, log10_likelihood_ratio);
      //mp_result[posIdx].lr_total += likelihood_ratio;
    }
  }
  if (isnan (likelihood_ratio))
    fprintf (stderr, "LR for the pedigree set is NAN.\n");

  //checkpt ();
  /* caculating the HET */
  for (j = 0; j < 5; j++) {
    alphaV = alpha[j][0];
    alphaV2 = 1 - alphaV;
    if (alphaV2 < 0)
      alphaV2 = 0;

    log10HetLR = 0;
    for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
      homoLR =
	pPedigree->likelihood / (pedigreeSet.nullLikelihood[pedIdx] *
				 pPedigree->markerLikelihood);
      //fprintf(stderr,"j=%d pedIdx=%d  %20.18f %20.16f %20.16f %20.16f \n",j, pedIdx,pPedigree->likelihood,pedigreeSet.nullLikelihood[pedIdx] * pPedigree->markerLikelihood, homoLR ,log10HetLR);
      if (alphaV * homoLR + alphaV2 < 0)
	fprintf (stderr, "HET LR less than 0. Check!!!\n");
      log10HetLR += log10 (alphaV * homoLR + alphaV2);
    }
    if (log10HetLR >= DBL_MAX_10_EXP - 1) {
      hetLR = DBL_MAX;
    } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
      hetLR = 0;
    } else {
      hetLR = pow (10, log10HetLR);
    }
    alpha_integral += hetLR * alpha[j][1];

    if (print_point_flag)
      fprintf (fphlod, "al=%f Hlod=%f\n", alphaV, hetLR);
    /*Update local maximum as necessary */
    if (hetLR > localmax_value) {
      localmax_value = hetLR;
      localmax_x[0] = gfreq;
      localmax_x[1] = alphaV;
      k = 2;
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	localmax_x[k] = x[k - 1];
	localmax_x[k + 1] =
	  (x[k] - xl[k]) * (x[k - 1] - xl[k - 1]) / (xu[k - 1] - xl[k - 1]) +
	  xl[k];
	localmax_x[k + 2] =
	  (x[k + 1] - xl[k + 1]) * (x[k] - xl[k]) / (xu[k] -
						     xl[k]) * (x[k - 1] -
							       xl[k -
								  1]) /
	  (xu[k - 1] - xl[k - 1]) + xl[k + 1];
	k += 3;
	if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	  localmax_x[k] = x[k - 1];
	  localmax_x[k + 1] = x[k];
	  localmax_x[k + 2] = x[k + 1];
	  k += 3;
	}
	/* threshold for QT */
	if (modelType.trait == CT) {
	  localmax_x[k] = x[k - 1];
	  k++;
	}
      }
    }
  }
  //checkpt ();
  avg_hetLR = alpha_integral;

  //  if (constraint >= 3.0 || constraint <= -3.0) {
  //fprintf (stderr, "Constraint is %f with avg hetLR =%f\n", constraint,
  //     avg_hetLR);
  //}


  /* Jacobian */
  k = 1;
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    avg_hetLR *=
      (x[k] - xl[k]) / (xu[k] - xl[k]) * (x[k] - xl[k]) / (xu[k] -
							   xl[k]) * (x[k +
								       1] -
								     xl[k +
									1]) /
      (xu[k + 1] - xl[k + 1]);
    k += 3;
    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      k += 3;
    }
    if (modelType.trait == CT) {
      k++;
    }
  }
  /*This is a temporary checking. */
  if (k != s->ndim) {
    fprintf (stderr, "k=%d  while dim for BR is %d\n", k, s->ndim);
    exit (0);
  }

  f[0] = avg_hetLR;

}






void
compute_hlod_mp_dt (double x[], double *f)
{

  int j;
  int pedIdx, liabIdx, status;

  double pen_DD, pen_Dd, pen_dd, gfreq, alphaV;
  double log10_likelihood_null, log10_likelihood_alternative,
    log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;
  double log10Likelihood;

  /* for null likelihood calculation */
  double product_likelihood = 1;	/* product of the likelihoods for all the pedigrees */
  double sum_log_likelihood = 0;	/* sum of the log10(likelihood) for all the pedigrees */

  // double *cw_center, *cw_hwidth;  /*pointers for centers and hwidths of the currrent working subregion*/

  Pedigree *pPedigree;
  int origLocus = locusList->pLocusIndex[0];

  if (locusList->numLocus > 1)
    origLocus = locusList->pLocusIndex[1];

  //fprintf(stderr,"in compute hlod x %G %G %G %G\n", x[0],x[1],x[2],x[3]);
  gfreq = x[0];

  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    pen_DD = x[3 * liabIdx + 1];
    pen_Dd = x[3 * liabIdx + 2] * x[3 * liabIdx + 1];
    pen_dd = x[3 * liabIdx + 3] * x[3 * liabIdx + 1] * x[3 * liabIdx + 2];
    pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
    pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
    pTrait->penetrance[2][liabIdx][1][0] = pen_Dd;
    pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
    pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
    pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
    pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_Dd;
    pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;

    //fprintf(stderr,"pene   DD=%G Dd=%G dd=%G\n", pen_DD, pen_Dd, pen_dd);

  }



  if (modelOptions.polynomial == TRUE);
  else
    /* only need to update trait locus */
    update_penetrance (&pedigreeSet, traitLocus);

  pLocus->pAlleleFrequency[0] = gfreq;
  pLocus->pAlleleFrequency[1] = 1 - gfreq;


  if (modelOptions.polynomial == TRUE);
  else
    update_locus (&pedigreeSet, traitLocus);

  /* for trait likelihood */
  locusList = &traitLocusList;
  xmissionMatrix = traitMatrix;

  // fprintf(stderr, "Null likelihood computation is starting for %d pedigrees \n",pedigreeSet.numPedigree);

  /* compute the null likelihood with   */
  pedigreeSet.likelihood = 1;
  pedigreeSet.log10Likelihood = 0;


  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {

    /* save the likelihood at null */
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];


    if (modelOptions.polynomial == TRUE) {
      KASSERT (pPedigree->traitLikelihoodPolynomial != NULL, "Error in  \n");
      /* evaluate likelihood */
      //fprintf(stderr, "evaluaing poly ");
      // printAllVariables();
      //expTermPrinting(stderr, pPedigree->traitLikelihoodPolynomial, 1);
      //fprintf(stderr, "\n");

      evaluatePoly (pPedigree->traitLikelihoodPolynomial,
		    pPedigree->traitLikelihoodPolyList,
		    &pPedigree->likelihood);
      //fprintf(stderr, " is done %f with %d pedigrees\n",pPedigree->likelihood, pedigreeSet.numPedigree);

      if (isnan (pPedigree->likelihood)) {
	fprintf (stderr, " likelihood is nan\n");
	exit (1);
      }

    } else {
      initialize_multi_locus_genotype (pPedigree);
      status = compute_pedigree_likelihood (pPedigree);
    }



    /*pPedigree->likelihood is now computed and now check it */
    if (pPedigree->likelihood == 0.0) {
      KLOG (LOGLIKELIHOOD, LOGWARNING,
	    "Pedigree %s has likelihood of 0 or too small.\n",
	    pPedigree->sPedigreeID);
      fprintf (stderr, "Pedigree %s has likelihood of 0 or too small.\n",
	       pPedigree->sPedigreeID);
      product_likelihood = 0.0;
      sum_log_likelihood = -9999.99;
      break;
    } else if (pPedigree->likelihood < 0.0) {
      KASSERT (pPedigree->likelihood >= 0.0,
	       "Pedigree %s with NEGATIVE likelihood - This is CRAZY!!!.\n",
	       pPedigree->sPedigreeID);
      product_likelihood = 0.0;
      sum_log_likelihood = -9999.99;
      break;
    } else {
      if (pPedigree->pCount[origLocus] == 1) {
	product_likelihood *= pPedigree->likelihood;
	log10Likelihood = log10 (pPedigree->likelihood);
      } else {
	product_likelihood *=
	  pow (pPedigree->likelihood, pPedigree->pCount[origLocus]);
	log10Likelihood =
	  log10 (pPedigree->likelihood) * pPedigree->pCount[origLocus];
      }
      sum_log_likelihood += log10Likelihood;
    }
    pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
    //fprintf(stderr,"null likelihood pedIdx=%d is done %20.18f with product =%20.16f\n",pedIdx,pPedigree->likelihood,product_likelihood );
  }


  pedigreeSet.likelihood = product_likelihood;
  pedigreeSet.log10Likelihood = sum_log_likelihood;
  log10_likelihood_null = pedigreeSet.log10Likelihood;
  KLOG (LOGLIKELIHOOD, LOGDEBUG, "Sum of log Likelihood is: %e\n",
	sum_log_likelihood);

  //fprintf(stderr,"Null likelihood = %20.15f\n", pedigreeSet.likelihood);

  /* This is for alternative likelihood */
  locusList = &savedLocusList;
  xmissionMatrix = altMatrix;
  compute_likelihood (&pedigreeSet);

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
  } else {
    log10_likelihood_ratio =
      log10_likelihood_alternative - log10_likelihood_null;
  }
  /* check for overflow problem !!! */
  if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
    likelihood_ratio = DBL_MAX;
  } else {
    /* check for underflow problem too !!! */
    if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
      likelihood_ratio = 0;
    } else {
      likelihood_ratio = pow (10.0, log10_likelihood_ratio);
      //mp_result[posIdx].lr_total += likelihood_ratio;
    }
  }
  /* add the result to the right placeholder */
  //mp_result[posIdx].lr_count++;
  //fprintf(stderr," %f %f %f %f %20.15f  ", gfreq,pen_DD, pen_Dd,pen_dd, log10_likelihood_alternative);
  /* caculating the HET */
  for (j = 0; j < 5; j++) {
    alphaV = alpha[j][0];
    alphaV2 = 1 - alphaV;
    if (alphaV2 < 0)
      alphaV2 = 0;

    log10HetLR = 0;
    for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
      homoLR =
	pPedigree->likelihood / (pedigreeSet.nullLikelihood[pedIdx] *
				 pPedigree->markerLikelihood);
      //fprintf(stderr,"j=%d pedIdx=%d  %20.18f %20.16f %20.16f %20.16f \n",j, pedIdx,pPedigree->likelihood,pedigreeSet.nullLikelihood[pedIdx] * pPedigree->markerLikelihood, homoLR ,log10HetLR);
      if (alphaV * homoLR + alphaV2 < 0)
	fprintf (stderr, "HET LR less than 0. Check!!!\n");
      log10HetLR += log10 (alphaV * homoLR + alphaV2);
    }
    if (log10HetLR >= DBL_MAX_10_EXP - 1) {
      hetLR = DBL_MAX;
      //mp_result[posIdx].het_lr_total = DBL_MAX;
    } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
      hetLR = 0;
    } else {
      hetLR = pow (10, log10HetLR);
      //mp_result[posIdx].het_lr_total += hetLR;
    }

    //fprintf(stderr,"j=%d het LR=%15.10f ",j, hetLR);

    alpha_integral += hetLR * alpha[j][1];

    if (hetLR > localmax_value) {
      localmax_value = hetLR;
      localmax_x[0] = gfreq;
      localmax_x[1] = alphaV;

      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	localmax_x[3 * liabIdx + 2] = x[3 * liabIdx + 1];
	localmax_x[3 * liabIdx + 3] = x[3 * liabIdx + 2] * x[3 * liabIdx + 1];
	localmax_x[3 * liabIdx + 4] =
	  x[3 * liabIdx + 3] * x[3 * liabIdx + 1] * x[3 * liabIdx + 2];
      }
    }

  }				/* end of calculating HET LR */




  avg_hetLR = alpha_integral;
  //fprintf(stderr,"avg hetLR =%15.10f \n", avg_hetLR);

  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    avg_hetLR *= x[3 * liabIdx + 1] * x[3 * liabIdx + 1] * x[3 * liabIdx + 2];
  }

  //  fprintf(stderr,"   hetLR =%f   ", avg_hetLR);


  *f = avg_hetLR;

  // return avg_hetLR;

}

void
compute_hlod_2p_qt (double x[], double *f)
{

  int k, j;
  int pedIdx, liabIdx = 0, status;
  double constraint = 0.0;
  double mean_DD;
  double mean_Dd;
  double mean_dd;
  double SD_DD = 0.0;
  double SD_Dd = 0.0;
  double SD_dd = 0.0;
  double gfreq;
  double alphaV;
  double theta;
  double threshold = 0.0;
  double log10_likelihood_null, log10_likelihood_alternative,
    log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;

  Pedigree *pPedigree;

  gfreq = x[0];
  theta = fixed_theta;		//x[5];  

  if (1 && modelOptions.markerAnalysis == FALSE) {
    pLocus->pAlleleFrequency[0] = gfreq;
    pLocus->pAlleleFrequency[1] = 1 - gfreq;

    if (modelOptions.polynomial == TRUE);
    else
      update_locus (&pedigreeSet, loc1);

  }

  /* this should be MEAN + SD */
  if (print_point_flag)
    fprintf (fphlod, "dp=%f th=%f gf=%f ", fixed_dprime, fixed_theta, gfreq);
  j = 1;
  if (modelOptions.markerAnalysis == FALSE) {
    for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
      mean_DD = x[j];
      mean_Dd =
	(x[j + 1] - xl[j + 1]) * (x[j] - xl[j]) / (xu[j] - xl[j]) + xl[j + 1];
      mean_dd =
	(x[j + 2] - xl[j + 2]) * (x[j + 1] - xl[j + 1]) / (xu[j + 1] -
							   xl[j +
							      1]) * (x[j] -
								     xl[j]) /
	(xu[j] - xl[j]) + xl[j + 2];
      if (print_point_flag)
	fprintf (fphlod, "pe %f %f %f ", mean_DD, mean_Dd, mean_dd);
      j += 3;
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	SD_DD = x[j];		//      modelRange.param[liabIdx][0][0][paramIdx];
	SD_Dd = x[j + 1];	//        modelRange.param[liabIdx][1][0][paramIdx];
	SD_dd = x[j + 2];	//        modelRange.param[liabIdx][2][0][paramIdx];
	if (print_point_flag)
	  fprintf (fphlod, "sd %f %f %f ", SD_DD, SD_Dd, SD_dd);
	j += 3;
      }
      /* threshold for QT */
      if (modelType.trait == CT) {
	threshold = x[j];	// modelRange.tthresh[liabIdx][thresholdIdx];
	if (print_point_flag)
	  fprintf (fphlod, "thd=%f ", threshold);
	j++;
      }

      /* check against the hard coded constraint */
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	constraint =
	  (1.0 - gfreq) * (1.0 - gfreq) * mean_dd * SD_dd + 2 * gfreq * (1 -
									 gfreq)
	  * mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
	/* fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
	   constraint, gfreq, mean_DD, SD_DD, mean_Dd, SD_DD, mean_dd, SD_dd);
	 */
	if (constraint >= 3.0 || constraint <= -3.0) {
	  // fprintf(stderr,"Constraint is %f \n", constraint);
	  num_out_constraint++;
	}
      }

      pTrait->means[liabIdx][0][0] = mean_DD;
      pTrait->means[liabIdx][0][1] = mean_Dd;
      pTrait->means[liabIdx][1][0] = mean_Dd;
      pTrait->means[liabIdx][1][1] = mean_dd;
      pTrait->stddev[liabIdx][0][0] = SD_DD;
      pTrait->stddev[liabIdx][0][1] = SD_Dd;
      pTrait->stddev[liabIdx][1][0] = SD_Dd;
      pTrait->stddev[liabIdx][1][1] = SD_dd;

      /* threshold for QT */
      pTrait->cutoffValue[liabIdx] = threshold;

    }				/* liability class Index */

    if (modelOptions.polynomial == TRUE);
    else
      update_penetrance (&pedigreeSet, traitLocus);


  }				/* marker to marker analysis */
  if (print_point_flag)
    fprintf (fphlod, "\n");

  /*This is a temporary checking. */
  if (j != s->ndim) {
    fprintf (stderr, "j=%d  while dim for BR is %d\n", j, s->ndim);
    exit (0);
  }

  /* get the likelihood at 0.5 first and LD=0 */
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    set_null_dprime (pLDLoci);
    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprime0Idx]);
    copy_DValue (pLDLoci, pLambdaCell->DValue[dprime0Idx]);

    KASSERT (pLambdaCell->impossibleFlag[dprime0Idx] == 0,
	     "Haplotype frequency combination impossible at LE. Exiting!\n");
  }

  for (k = 0; k < 3; k++) {
    locusList->pNextLocusDistance[k][0] = 0.5;
    locusList->pPrevLocusDistance[k][1] = 0.5;
  }

  if (modelOptions.polynomial == TRUE);
  else
    /* populate the matrix */
    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */

  KLOG (LOGLIKELIHOOD, LOGDEBUG, "NULL Likelihood\n");
  compute_likelihood (&pedigreeSet);



  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    fprintf (stderr, "Theta 0.5 has likelihood 0\n");
    fprintf (stderr, "dgf=%f\n", gfreq);
    for (j = 1; j < s->ndim; j++) {
      fprintf (stderr, " %f", x[j]);
    }
    fprintf (stderr, "\n");

    exit (-1);
  }

  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
    /* save the likelihood at null */
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
    pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
  }

  log10_likelihood_null = pedigreeSet.log10Likelihood;
  //for (dprimeIdx = 0;dprimeIdx < pLambdaCell->ndprime;dprimeIdx++){
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);
//    if (pLambdaCell->impossibleFlag[dprimeIdx] != 0)
//        continue;
    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
    copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);
  }

  if (modelOptions.mapFlag == SA) {
    for (k = 0; k < 3; k++) {
      locusList->pNextLocusDistance[k][0] = theta;
      locusList->pPrevLocusDistance[k][1] = theta;
    }
  } else {
    printf ("mapflag sould be SA\n");
    exit (-1);
  }

  if (modelOptions.polynomial == TRUE);
  else
    /* populate the matrix */
    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */

  KLOG (LOGLIKELIHOOD, LOGDEBUG, "ALT Likelihood\n");
  compute_likelihood (&pedigreeSet);

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
  } else {
    log10_likelihood_ratio =
      log10_likelihood_alternative - log10_likelihood_null;
  }

  /* check for overflow problem !!! */
  if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
    likelihood_ratio = DBL_MAX;
  } else {
    /* check for underflow problem too !!! */
    if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
      likelihood_ratio = 0;
    } else {
      likelihood_ratio = pow (10.0, log10_likelihood_ratio);
    }
  }

  /* caculating the HET */
  for (j = 0; j < 5; j++) {
    alphaV = alpha[j][0];
    alphaV2 = 1 - alphaV;
    if (alphaV2 < 0)
      alphaV2 = 0;

    log10HetLR = 0;
    for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
      homoLR = pPedigree->likelihood / pedigreeSet.nullLikelihood[pedIdx];
      log10HetLR += log10 (alphaV * homoLR + alphaV2);
    }
    if (log10HetLR >= DBL_MAX_10_EXP - 1) {
      hetLR = DBL_MAX;
    } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
      hetLR = 0;
    } else {
      hetLR = pow (10, log10HetLR);
    }
    alpha_integral += hetLR * alpha[j][1];

    if (print_point_flag)
      fprintf (fphlod, "al=%f Hlod=%f\n", alphaV, hetLR);
    /*Update local maximum as necessary */
    if (hetLR > localmax_value) {
      localmax_value = hetLR;
      localmax_x[0] = gfreq;
      localmax_x[1] = alphaV;
      k = 2;
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	localmax_x[k] = x[k - 1];
	localmax_x[k + 1] =
	  (x[k] - xl[k]) * (x[k - 1] - xl[k - 1]) / (xu[k - 1] - xl[k - 1]) +
	  xl[k];
	localmax_x[k + 2] =
	  (x[k + 1] - xl[k + 1]) * (x[k] - xl[k]) / (xu[k] -
						     xl[k]) * (x[k - 1] -
							       xl[k -
								  1]) /
	  (xu[k - 1] - xl[k - 1]) + xl[k + 1];
	k += 3;
	if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	  localmax_x[k] = x[k - 1];
	  localmax_x[k + 1] = x[k];
	  localmax_x[k + 2] = x[k + 1];
	  k += 3;
	}
	/* threshold for QT */
	if (modelType.trait == CT) {
	  localmax_x[k] = x[k - 1];
	  k++;
	}
      }
    }
  }

  avg_hetLR = alpha_integral;

  //  if (constraint >= 3.0 || constraint <= -3.0) {
  //fprintf (stderr, "Constraint is %f with avg hetLR =%f\n", constraint,
  //   avg_hetLR);
//}


  /* Jacobian */
  k = 1;
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    avg_hetLR *=
      (x[k] - xl[k]) / (xu[k] - xl[k]) * (x[k] - xl[k]) / (xu[k] -
							   xl[k]) * (x[k +
								       1] -
								     xl[k +
									1]) /
      (xu[k + 1] - xl[k + 1]);
    k += 3;
    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      k += 3;
    }
    if (modelType.trait == CT) {
      k++;
    }
  }
  /*This is a temporary checking. */
  if (k != s->ndim) {
    printf ("k=%d  while dim for BR is %d\n", k, s->ndim);
    exit (0);
  }

  *f = avg_hetLR;


}


void
compute_hlod_2p_dt (double x[], double *f)
{
//double compute_hlod(PedigreeSet *pedigreeSet,double x[], int loc1, int loc2, Locus *pLocus, Trait *pTrait, int traitLocus, int totalLoci, double * initialProbAddr[3], Locus *pLocus1){

/*  Limit of this function

   modelOptions.type := TP  (Two points)
   modelOptions.trait := DT (DICHOTOMOUS);
   modelOptions.equilibrium :=LINKAGE_EQUILIBRIUM
   modelOptions.polynomial := TRUE
   modelOptions->markerAnalysis := FALSE;
   modelOptions.mapFlag := SA 
   modelRange.nafreq :=1
   
*/

  int k, j;
  int pedIdx, liabIdx = 0, status;

  double pen_DD, pen_Dd, pen_dd, gfreq, alphaV, theta;
  double log10_likelihood_null, log10_likelihood_alternative,
    log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, tmp, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;

  Pedigree *pPedigree;

  gfreq = x[0];
  theta = fixed_theta;		//x[5];  
  //printf("Calculating hetLR with gf=%f DD=%f Dd=%f dd=%f theta=%f\n", gfreq, pen_DD,pen_Dd, pen_dd, fixed_theta);
  if (1 && modelOptions.markerAnalysis == FALSE) {
    pLocus->pAlleleFrequency[0] = gfreq;
    pLocus->pAlleleFrequency[1] = 1 - gfreq;

    if (modelOptions.polynomial == TRUE);
    else
      update_locus (&pedigreeSet, loc1);

  }

  if (modelOptions.markerAnalysis == FALSE
      && pLocus1->locusType == LOCUS_TYPE_TRAIT) {
    for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
      pen_DD = x[3 * liabIdx + 1];
      pen_Dd = x[3 * liabIdx + 2] * x[3 * liabIdx + 1];
      pen_dd = x[3 * liabIdx + 3] * x[3 * liabIdx + 1] * x[3 * liabIdx + 2];
      pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
      pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
      pTrait->penetrance[2][liabIdx][1][0] = pen_Dd;
      pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
      pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
      pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
      pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_Dd;
      pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;
    }
  }
  if (modelOptions.polynomial == TRUE);
  else
    update_penetrance (&pedigreeSet, traitLocus);

  /* get the likelihood at 0.5 first and LD=0 */
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    set_null_dprime (pLDLoci);
    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprime0Idx]);
    copy_DValue (pLDLoci, pLambdaCell->DValue[dprime0Idx]);
    KASSERT (pLambdaCell->impossibleFlag[dprime0Idx] == 0,
	     "Haplotype frequency combination impossible at LE. Exiting!\n");
  }
  for (k = 0; k < 3; k++) {
    locusList->pNextLocusDistance[k][0] = 0.5;
    locusList->pPrevLocusDistance[k][1] = 0.5;
  }
  if (modelOptions.polynomial == TRUE);
  else
    /* populate the matrix */
    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */


  KLOG (LOGLIKELIHOOD, LOGDEBUG, "NULL Likelihood\n");
  compute_likelihood (&pedigreeSet);

  //printf("likelihood =%15.13f with theta 0.5 with %d pedigrees\n", pedigreeSet.likelihood, pedigreeSet.numPedigree);
  //scanf("%d ",&k);
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    fprintf (stderr, "Theta 0.5 has likelihood 0\n");
    fprintf (stderr, "dgf=%f\n", gfreq);
    exit (-1);
  }

  /* save the results for NULL */
  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
    /* save the likelihood at null */
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
    pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
  }

  log10_likelihood_null = pedigreeSet.log10Likelihood;

  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);
//      if (pLambdaCell->impossibleFlag[dprimeIdx] != 0)
//        continue;
    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
    copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);

    /* calculate R square if the marker is a SNP */
    if (R_square_flag == TRUE)
      R_square =
	calculate_R_square (pLocus1->pAlleleFrequency[0],
			    pLocus2->pAlleleFrequency[0],
			    pLDLoci->ppDValue[0][0]);
    else
      R_square = -1;
  }


  if (modelOptions.mapFlag == SA) {

    for (k = 0; k < 3; k++) {
      locusList->pNextLocusDistance[k][0] = theta;
      locusList->pPrevLocusDistance[k][1] = theta;

    }
  } else {
    printf ("mapflag sould be SA\n");
    exit (-1);
  }

  if (modelOptions.polynomial == TRUE);
  else
    /* populate the matrix */
    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */

  KLOG (LOGLIKELIHOOD, LOGDEBUG, "ALT Likelihood\n");
  compute_likelihood (&pedigreeSet);

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;

  //printf("likelihood =%15.13f with theta %f  %d pedigree\n", pedigreeSet.likelihood,fixed_theta, pedigreeSet.numPedigree);                                  
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
  } else {
    log10_likelihood_ratio =
      log10_likelihood_alternative - log10_likelihood_null;
  }


  /* check for overflow problem !!! */
  if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
    likelihood_ratio = DBL_MAX;
  } else {
    /* check for underflow problem too !!! */
    if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
      likelihood_ratio = 0;
    } else {
      likelihood_ratio = pow (10.0, log10_likelihood_ratio);
    }
  }
  /* caculating the HET */
  for (j = 0; j < 5; j++) {
    alphaV = alpha[j][0];
    //for (j = 0; j < modelRange.nalpha; j++){
    //alphaV = modelRange.alpha[j];         
    alphaV2 = 1 - alphaV;
    if (alphaV2 < 0)
      alphaV2 = 0;

    log10HetLR = 0;

    for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
      homoLR = pPedigree->likelihood / pedigreeSet.nullLikelihood[pedIdx];
      tmp = log10 (alphaV * homoLR + (1 - alphaV));
      log10HetLR += tmp * pPedigree->pCount[loc2];
    }

    if (log10HetLR >= __DBL_MAX_10_EXP__ - 1) {
      hetLR = __DBL_MAX__;
    } else if (log10HetLR <= __DBL_MIN_10_EXP__ + 1) {
      hetLR = 0;
    } else {
      hetLR = pow (10, log10HetLR);
    }

    alpha_integral += hetLR * alpha[j][1];
    //alpha_integral +=hetLR;

    if (hetLR > localmax_value) {
      localmax_value = hetLR;
      localmax_x[0] = gfreq;
      localmax_x[1] = alphaV;

      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	localmax_x[3 * liabIdx + 2] = x[3 * liabIdx + 1];
	localmax_x[3 * liabIdx + 3] = x[3 * liabIdx + 2] * x[3 * liabIdx + 1];
	localmax_x[3 * liabIdx + 4] =
	  x[3 * liabIdx + 3] * x[3 * liabIdx + 1] * x[3 * liabIdx + 2];
      }
    }
    // fprintf(fphlod,"%f %f %f %f %f %f %f\n", log10(hetLR*x[1]*x[1]*x[2]), gfreq, pen_DD,pen_Dd, pen_dd, alphaV,fixed_theta);
  }				//end of calculating the HET         


  avg_hetLR = alpha_integral;
  //avg_hetLR= alpha_integral/modelRange.nalpha;

  //printf("avg hetLR =%15.10f with gf=%f DD=%f Dd=%f dd=%f theta=%f\n", avg_hetLR, gfreq, pen_DD,pen_Dd, pen_dd, fixed_theta);

  // avg_hetLR *=x[1]*x[1]*x[2];
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    avg_hetLR *= x[3 * liabIdx + 1] * x[3 * liabIdx + 1] * x[3 * liabIdx + 2];
  }


  f[0] = avg_hetLR;

  //  return 0.0;

}
