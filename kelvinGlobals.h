char *programVersion = "V0.37.0";       ///< Overall kelvin version set upon release.

/* Some default global values. */
char resultsprefix[KMAXFILENAMELEN + 1] = "./"; ///< Path for SR directive result storage
char markerfile[KMAXFILENAMELEN + 1] = "markers.dat";   ///< Default name (and storage) for marker file
char maxmodelfile[KMAXFILENAMELEN + 1] = "tp.out";   ///< Default name (and storage) for maximizing model file
char mapfile[KMAXFILENAMELEN + 1] = "mapfile.dat";      ///< Default name (and storage) for map file
char pedfile[KMAXFILENAMELEN + 1] = "pedfile.dat";      ///< Default name (and storage) for pedigree file
char datafile[KMAXFILENAMELEN + 1] = "datafile.dat";    ///< Default name (and storage) for marker data file
char ccfile[KMAXFILENAMELEN + 1] = "";  ///< Case control count file
char avghetfile[KMAXFILENAMELEN + 1] = "br.out";        ///< Default name (and storage) for Bayes Ratio file
char pplfile[KMAXFILENAMELEN + 1] = "ppl.out";  ///< Default name (and storage) for PPL file
char ldPPLfile[KMAXFILENAMELEN + 1] = "ldppl.out";
FILE *fpHet = NULL;     ///< Average HET LR file (Bayes Ratio file) pointer
FILE *fpPPL = NULL;     ///< PPL output file pointer
FILE *fpTP = NULL;      ///< Ancillary Two-point output, used to go to stderr
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
