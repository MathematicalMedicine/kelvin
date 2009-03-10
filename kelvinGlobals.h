char *programVersion = "V0.37.1";       ///< Overall kelvin version set upon release.
char *kelvinVersion = "$Id$";        ///< svn's version for kelvin.c

struct swStopwatch *overallSW;
char messageBuffer[MAXSWMSG];

/* Some default global values. */

/* Storage and default names for files that are always opened (depending on analysis options) */
char markerfile[KMAXFILENAMELEN + 1] = "markers.dat";  /// Marker (frequency) file
char mapfile[KMAXFILENAMELEN + 1] = "mapfile.dat";     /// Map file
char pedfile[KMAXFILENAMELEN + 1] = "pedfile.dat";     /// Pedigree file
char datafile[KMAXFILENAMELEN + 1] = "datafile.dat";   /// Data (pedigree description) file
char avghetfile[KMAXFILENAMELEN + 1] = "br.out";       /// Bayes Ratio file
char pplfile[KMAXFILENAMELEN + 1] = "ppl.out";         /// PPL file
char condFile[KMAXFILENAMELEN + 1] = "condL.out";      /// Conditional LR file
char ldPPLfile[KMAXFILENAMELEN + 1] = "ldppl.out";     /// This appears to be unused

/* Storage for files that are only opened based on explicit directives */
char ccfile[KMAXFILENAMELEN + 1] = "";                 /// Case control count file
char modfile[KMAXFILENAMELEN + 1] = "";                /// MOD and maximizing model file
char maxmodelfile[KMAXFILENAMELEN + 1] = "";           /// verbose Maximizing Model file
char intermediatefile[KMAXFILENAMELEN + 1] = "";       /// Intermediate Result file
char dkelvinoutfile[KMAXFILENAMELEN + 1] = "";         /// DCHURE detail file
char resultsprefix[KMAXFILENAMELEN + 1] = "./"; ///< Path for SR directive result storage

FILE *fpHet = NULL;     ///< Average HET LR file (Bayes Ratio file) pointer
FILE *fpPPL = NULL;     ///< PPL output file pointer
FILE *fpCond = NULL;    ///< Conditional LR for genetic counseling
FILE *fpMOD = NULL;     // MOD and maximizing model information
FILE *fpTP = NULL;      ///< Ancillary Two-point output, used to go to stderr
FILE *fpIR = NULL;      ///< Intermediate results, used to go to stderr, normally dkelvin-only
FILE *fpDK = NULL;      // DCHURE detail file

int polynomialScale = 1;        ///< Scale of static allocation and dynamic growth in polynomial.c.

/** Model datastructures. modelOptions is defined in the pedigree library. */
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
