char *programVersion = "V0.38.0";       ///< Overall kelvin version set upon release.
char *kelvinVersion = "$Id$";        ///< svn's version for kelvin.c

struct swStopwatch *overallSW;
char messageBuffer[MAXSWMSG];

/* Some default global values. */

/** Model datastructures. */
ModelType modelType;
ModelRange modelRange;
ModelOptions modelOptions;

/** Number of D primes. If there are more than 2 alleles in the
  marker/trait, number of D primes and D prime ranges are assumed to
  be the same to reduce complexity for initial phase of this
  project. */
int num_of_d_prime;
double *d_prime;
int num_of_theta;

/** Three dimensional array for the two point summary results.
  - first dimension is the D prime, for LE, D prime=0 with just one element
    in this dimension
  - second dimension is theta values 
  - third dimension is marker allele frequency, for LE, only one element in 
    this dimension */
SUMMARY_STAT ***tp_result;

/** Two dimensional array per (dprime, theta).
  This will be used to calculate PPL. */

/** Storage for the NULL likelihood for the multipoint calculation under polynomial. */
double markerSetLikelihood;

/** For multipoint, we use genetic map positions on a chromosome. */
double *map_position;
int num_of_map_position;

/** One dimensional array, indexing by map position.  For multipoint,
 we don't know how to incorporate LD yet.  This map could be sex
 specific map or sex averaged map. For two point, we don't have to
 distinguish sex specific/avearge as we use theta relative to marker
 during analysis and after analysis (result) */
SUMMARY_STAT *mp_result;
int numPositions;

/** Transmission matrices provide the pre-computed probability of
 inheritance of a a given combination of marker and trait alleles. */
XMission *nullMatrix;
XMission *altMatrix;
XMission *traitMatrix;
XMission *markerMatrix;

int markerSetChanged; /* Flag for multipoint analysis, did set of markers change? */
int locusListChanged; /* flag for multipoint analysis, did relative trait position or marker set change? */

int prevFirstMarker;		/* first marker in the set for multipoint analysis */
int prevLastMarker;		/* last marker in the set for multipoint analysis */

LambdaCell *pLambdaCell = NULL;
int loopMarkerFreqFlag = 0;
int total_count;

char *flexBuffer = NULL;
int flexBufferSize = 0;

double *****likelihoodQT = NULL;
double **likelihoodDT = NULL;   ///< This is now for homeLR


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

FILE *fpCond = NULL;    ///< Conditional LR for genetic counseling, global due to likelihood.c write!
