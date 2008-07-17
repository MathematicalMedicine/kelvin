/**
@mainpage

  kelvin - Linkage and Linkage Disequilibrium Analysis Program.

  Implementation of the Elston-Stewart algorithm for linkage
  analysis. Currently supports two-point and multipoint analyses,
  dichotomous and quantitative traits, linkage equilibrium and
  disequilibrium, case control, and many other options.
  
  Copyright 2008, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  See http://hodgkin.ccri.net/software/kelvin/index.html for full
  documentation.

  AUTHORS

  - Yungui Huang - principle author
  - Hongling Wang - Polynomial features
  - Alberto Maria Segre - config and error logging modules
  - Nathan Burnette - Regex code
  - Bill Valentine-Cooper - clean-up, performance and tracking


@file kelvin.c

  Main driver and way too many loops.

  USAGE:
  <pre>
  kelvin <kelvin.conf>

  where <kelvin.conf> is a text file containing directives
  describing the locations of supporting files and specifying the
  nature of the analysis to be performed. See the provided documentation
  for details.
  </pre>
  LIMITATIONS

  Currently only handles biallelic traits.

  COMPILE-TIME CONDITIONALS

  There are numerous compilation conditions that affect kelvin. Many
  of them are diagnostic in nature.

  - _OPENMP - this conditional is not defined by the user, but rather
  set by the compiler when the -fopenmp flag is specified to enable
  OpenMP multithreading. Aspects of polynomial build and evaluation
  are currently setup to handle multithreading.  As a caveat --
  polynomial evaluation performance is only improved with
  multi-threading if the analysis has more unique pedigrees than
  threads. Performance will actually degrade when there are only a few
  unique pedigrees because term optimization will cause their
  polynomials to share most of their terms, and then per-pedigree
  evaluation threads will get into extensive synchronization conflicts
  traversing those polynomials. You can monitor this by watching the
  System CPU time as opposed to User CPU time. System CPU time is
  almost exclusively wasted time due to thread synchronization
  conflicts. Note that it is not, however, necessary to disable
  multi-threading by rebuilding kelvin for these situations. Simply
  set the OMP_NUM_THREADS environment variable to 1, and no locking
  conflicts will occur.

  - MEMSTATUS - handy when there is concern about memory capacity,
  but can really clutter-up the display. Runs in a child process and
  reports on the main process' memory usage every 30 seconds. This
  is only possible on systems with a supported pmap command.

  - MEMGRAPH - provides the same information as MEMSTATUS but in a
  format appropriate for graphing with a tool like gnuplot. The file
  is named kelvin_<pid>_memory.dat, where <pid> is the process id.
  These files will need to be cleaned-up periodically, so this 
  conditional is not on by default. This is only possible on systems
  with a supported pmap command.

  - SIMPLEPROGRESS - suppresses more complicated progress reports in
  favor of a single percentage progress and estimated remaining time.
  This simple progress indication cannot be linear, but is easier to
  understand. This conditional is on by default.

  - POLYSTATISTICS - enable or disable display of extensive polynomial
  statistics at 2Mp creation intervals as well as at key points in
  pedigree processing. This does not affect performance or processing.
  
  - DMTRACK - enable exhaustive memory management tracking. Keeps track
  of all allocations and frees by source code line, and can provide
  cross-references at any time by total utilization or retention. Can
  also report on cross-module operations, i.e. one module allocates and
  another module frees the same chunk of memory.

  - DMUSE
  - SOURCEDIGRAPH

*/
#include <signal.h>
#include <gsl/gsl_version.h>
#include "kelvin.h"
#include "likelihood.h"
#include "saveResults.h"
#include "trackProgress.h"

extern char *likelihoodVersion, *locusVersion, *polynomialVersion;
extern Polynomial *constant1Poly;

struct swStopwatch *overallSW;  ///< Performance timer used throughout code for overall statistics.
time_t startTime;
char messageBuffer[MAXSWMSG];   ///< Commonly-used message buffer sized to work with swLogMsg().
volatile sig_atomic_t statusRequestSignal = FALSE;      ///< Status update requested via signal

/**

  Handler for SIGQUIT.

  Typing ^\ or issuing a kill -s SIGQUIT gets a dump of statistics.

  P.S. - cygwin requires "stty quit ^C" first for this to work.

*/
void quitSignalHandler (int signal)
{
  statusRequestSignal = TRUE;
#ifdef POLYSTATISTICS
  if (modelOptions.polynomial == TRUE)
    polyDynamicStatistics ("Signal received");
  else
#endif
    swDump (overallSW);
#ifdef DMTRACK
  swLogPeaks ("Signal");
#endif
}

/**

  Handler for SIGUSR1.

  Handles signal is raised by our child process to get status update.
  Sets the signalSeen flag which is watched for in breaks in the
  code.

*/
void usr1SignalHandler (int signal)
{
  statusRequestSignal = TRUE;
}

#if defined (GPROF) || (GCOV)
/*

  Handler for SIGTERM.

  If we're profiling or doing coverage analysis, we catch a SIGTERM to 
  allow early exit(). An orderly exit like this (with EXIT_SUCCESS
  status, ensures that the profileing or coverage information completely
  written and fit for analysis.

*/
void termSignalHandler (int signal)
{
  fprintf (stderr, "Terminating early for gprof or gcov!\n");
  exit (EXIT_SUCCESS);
}
#endif

/// Handler for SIGINT
void intSignalHandler (int signal)
{
  fprintf (stderr, "Terminating early via interrupt!\n");
  exit (EXIT_FAILURE);
}

pid_t childPID = 0;     ///< For a child process producing timing (and memory?) statistics.
/**

  General-purpose exit handler

  Exit handler to clean-up after we hit any of our widely-distributed 
  exit points. Ensures that errant child processes are handled so
  we don't pester people with unnecessary messages.

*/
void exit_kelvin ()
{
  swLogMsg ("Exiting");
  if (childPID != 0)
    kill (childPID, SIGKILL);   /* Sweep away any errant children */
}

/**

  Global variables

*/
char *programVersion = "V0.35.0";       ///< Overall kelvin version set upon release.
char *kelvinVersion = "$Id$";        ///< svn's version for kelvin.c

/* Some default global values. */
char resultsprefix[KMAXFILENAMELEN + 1] = "./"; ///< Path for SR directive result storage
char markerfile[KMAXFILENAMELEN + 1] = "markers.dat";   ///< Default name (and storage) for marker file
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
//double ***likelihoodDT = NULL;   // This is now moved into each pedigree
double **likelihoodDT = NULL;   ///< This is now for homeLR
double *****likelihoodQT = NULL;        ///< This is now moved into each pedigree (really?)
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

/**

  Driver for all types of analyses.

  <pre>
  Usage:
     kelvin kelvin.conf
  </pre>
  The kelvin.conf file gives information about the specific linkage
  analysis run. All information about, e.g., which markers to use,
  what outputs to calculate, and so on, are stored in this
  configuration file.

*/
int main (int argc, char *argv[])
{
  int i, j;
  char configfile[KMAXFILENAMELEN] = "";
  int breakFlag = FALSE;
  double alphaV, alphaV2;
  int loc1, loc2;
  int cL[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 },    ///< Actual # calls to each instance of compute_likelihood
    eCL[9] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0};   ///< Est. final # calls to each instance of compute_likelihood

  PedigreeSet pedigreeSet;      /* Pedigrees. */

  /* Start GAW */
  Pedigree *pPedigree;
  double pen_DD, pen_Dd, pen_dd;
  double mean_DD, mean_Dd, mean_dd;
  double SD_DD, SD_Dd, SD_dd;
  double gfreq; /* disease gene frequency */
  double theta[2];      /* theta */
  int penIdx, liabIdx, gfreqInd, thetaInd;
  int paramIdx = -1;
  int dprimeIdx;

  //double likelihood_null, likelihood_alternative;
  double log10_likelihood_null, log10_likelihood_alternative;
  double likelihood_ratio;
  double log10_likelihood_ratio;
  Locus *pLocus;
  Locus *pLocus1, *pLocus2;
  Trait *pTrait;
  int pedIdx;

  //  int pedID;
  double homoLR, hetLR;
  double adjustedHetLR = 0;
  double log10HetLR;
  double max;
  double constraint;
  LDLoci *pLDLoci = NULL;
  double traitPos;      /* trait position for multipoint analysis */
  TraitLocus *pTraitLocus;
  int traitLocus;
  int leftMarker = -1;
  int posIdx;
  int k;
  double avgLR;
  double ppl;
  double ldppl, ppld;
  int markerSetChanged; /* flag for multipoint analysis */
  int locusListChanged; /* flag for multipoint analysis */
  int prevFirstMarker;  /* first marker in the set for multipoint analysis */
  int prevLastMarker;   /* last marker in the set for multipoint analysis */
  int prevTraitInd;
  double *prevPos, *currPos;    /* for MP */
  int locus;
  int thresholdIdx = -1;
  double threshold = 0;
  int R_square_flag = FALSE;
  double R_square = 0;
  int dprime0Idx = 0;
  int mkrFreqIdx;
  double mkrFreq;
  int maxThetaIdx = 0;
  int maxDPrimeIdx = 0;
  int maxDPrimeIdx_at_theta0 = 0;
  double max_at_theta0;
  double lr;
  int theta0Idx = 0;
  double max_at_dprime0;
  int maxTheta_at_dprime0 = -1;
  int totalLoci;
  int status;
  double *marker1Pos, *marker2Pos;
  double relativePos;
  int traitIndex = 0;
  double dist;
  int initialFlag = 0;
  double tmp;
  int polynomialFlag;

  SubLocusList savedLocusList;
  SubLocusList traitLocusList;
  SubLocusList markerLocusList;

  char **markerNameList = NULL;

  Polynomial *initialProbPoly[3];
  Polynomial *initialProbPoly2[3];
  double initialProb[3];
  void *initialProbAddr[3];
  double initialProb2[3];
  void *initialProbAddr2[3];
  void *initialHetProbAddr[3];
  char *tmpID;
  int exitDueToLoop = FALSE;

#ifdef _OPENMP
  char *envVar;
  int threadCount = 0;
#endif

  struct swStopwatch *combinedComputeSW,        ///< Combined likelihood compute stopwatch
   *combinedBuildSW;    ///< Combined likelihood polynomial build stopwatch

  overallSW = swCreate ("overall");
  combinedComputeSW = swCreate ("combinedComputeSW");
  combinedBuildSW = swCreate ("combinedBuildSW");
  startTime = time (NULL);

  /* Add an exit handler to deal with wayward children. */

  if (atexit (exit_kelvin)) {
    perror ("Could not register exit handler!");
    exit (EXIT_FAILURE);
  }

  /* Fork a child that loops sleeping several seconds and then signalling 
   * us with SIGUSR1 to do an asynchronous dump of peak statistitics to stderr. */

  int currentVMK, maximumVMK;   ///< Properly the property of the child process.
  if ((maximumVMK = swGetMaximumVMK ()) != 0) {
    childPID = fork ();
    if (childPID == 0) {
      /* Code executed by child only! */
      pid_t parentPID = 0;

      /* Ignore QUIT signals, 'cause they're actually status requests for Mom. */
      signal (SIGQUIT, SIG_IGN);
#ifdef MEMGRAPH
      FILE *graphFile;
      char graphFileName[64];
      sprintf (graphFileName, "kelvin_%d_memory.dat", getppid ());
      if ((graphFile = fopen (graphFileName, "w")) == NULL) {
        perror ("Cannot open memory graph file!");
        exit (EXIT_FAILURE);
      }
#endif
      int wakeCount = 0;
      while (1) {
        sleep (30);
        wakeCount++;
        /* See if we've been reparented due to Mom's demise. */
        parentPID = getppid ();
        if (parentPID == 1) {
#ifdef MEMGRAPH
          fclose (graphFile);
#endif
          exit (EXIT_SUCCESS);
        }
        if (wakeCount % 2)
          kill (getppid (), SIGUSR1);   // Send a status-updating signal to parent.
        currentVMK = swGetCurrentVMK (getppid ());
#ifdef MEMGRAPH
        fprintf (graphFile, "%lu, %d\n", time (NULL) - startTime, currentVMK);
        fflush (graphFile);
#endif
#ifdef MEMSTATUS
        fprintf (stdout, "%lus, %dKb (%.1f%% of %.1fGb)\n", time (NULL) - startTime, currentVMK,
                 currentVMK / (maximumVMK / 100.0), maximumVMK / (1024.0 * 1024.0));
#endif
      }
    }
  }

  /* Setup signal handlers */
  struct sigaction usr1Action, quitAction, intAction;
  sigset_t usr1BlockMask, quitBlockMask, intBlockMask;

#if defined (GPROF) || (GCOV)
  struct sigaction termAction;
  sigset_t termBlockMask;
#endif
  sigfillset (&usr1BlockMask);
  usr1Action.sa_handler = usr1SignalHandler;
  usr1Action.sa_mask = usr1BlockMask;
  usr1Action.sa_flags = 0;
  sigaction (SIGUSR1, &usr1Action, NULL);

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
  char currentWorkingDirectory[MAXSWMSG - 32];

  sprintf (messageBuffer, "kelvin %s built %s %s", programVersion, __DATE__, __TIME__);
  swLogMsg (messageBuffer);
  swLogMsg (kelvinVersion);
  swLogMsg (likelihoodVersion);
  swLogMsg (locusVersion);
  swLogMsg (polynomialVersion);
  sprintf (messageBuffer, "Compiler %s, GSL %s\n", __VERSION__, GSL_VERSION);
  swLogMsg (messageBuffer);

#ifdef _OPENMP
  if ((envVar = getenv ("OMP_NUM_THREADS")) != NULL)
    threadCount = atoi (envVar);
  sprintf (messageBuffer, "OpenMP-enabled w/%d threads.", threadCount);
  swLogMsg (messageBuffer);
#endif
#ifdef FAKEEVALUATE
  swLogMsg ("Polynomial evaluation is being SKIPPED FOR TESTING, results will be wrong!");
#endif
#ifdef DMUSE
  swLogMsg ("Experimental static memory handling enabled!");
#endif
#ifdef DMTRACK
  swLogMsg ("Dynamic memory usage dumping is turned on, so performance will be poor!");
#endif
#ifdef GPROF
  sprintf (messageBuffer, "GNU profiler (gprof) run, use \"kill -%d %d\" to finish early.", SIGTERM, getpid ());
  swLogMsg (messageBuffer);
#endif
#ifdef GCOV
  sprintf (messageBuffer, "GNU coverage analyzer (gcov) run, use \"kill -%d %d\" to finish early.", SIGTERM, getpid ());
  swLogMsg (messageBuffer);
#endif
  fprintf (stdout, "To check status (at some risk), type CTRL-\\ or type \"kill -%d %d\".\n", SIGQUIT, getpid ());
  swStart (overallSW);

  memset (&savedLocusList, 0, sizeof (savedLocusList));
  memset (&markerLocusList, 0, sizeof (markerLocusList));
  memset (&traitLocusList, 0, sizeof (traitLocusList));

  memset (&modelType, 0, sizeof (modelType));

  // Initialize the logging system.
  logInit ();

  /* Start by parsing command line arguments. Most essential: figure
   * out where the configuration file lives. */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-')
      switch (argv[i][1]) {
      case '?':
        /* Help */
        fprintf (stdout, "Usage:\n");
        fprintf (stdout, "  %s [-?] <configuration file>\nwhere:\n", argv[0]);
        fprintf (stdout, "      -? : this output;\n");
        fprintf (stdout, "      <configuration file> : file containing run parameters.\n");
        exit (EXIT_FAILURE);
        break;
    } else if (strlen (configfile) != 0) {
      /* Unexpected argument; we already have a configuration file! Punt. */
      KLOG (LOGDEFAULT, LOGFATAL, "Unexpected command line argument '%s'; aborting.\n", argv[i]);
    } else if (strlen (argv[i]) >= KMAXFILENAMELEN) {
      /* Configuration file name too long! Punt. */
      KLOG (LOGDEFAULT, LOGFATAL, "Config file name '%s' too long (>=%d); aborting.\n", argv[i], KMAXFILENAMELEN);
    } else {
      /* Got a configuration file name. Copy it. */
      strncpy (configfile, argv[i], KMAXFILENAMELEN);
      getcwd (currentWorkingDirectory, sizeof (currentWorkingDirectory));
      sprintf (messageBuffer, "In %s w/%s", currentWorkingDirectory, configfile);
      swLogMsg (messageBuffer);
    }
    i++;
  }
  /* Check to see if the configuration file name was specified. */
  KASSERT ((strlen (configfile) > 0), "No configuration file specified; aborting.\n");

  // IMHO, ALL OF THE FOLLOWING UP TO OPENING FILES BELONGS IN CONFIG.C

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

  /* For now, reject all models we can't deal with. */
  KASSERT (modelRange.nalleles == 2, "Only biallelic traits supported.\n");

  /* the difference between QT and CT is whether we use threshold or not. Under CT -  yes to
   * threshold, under QT - no threshold */
  if (modelRange.ntthresh > 0 && modelType.trait != DT) {
    modelType.trait = CT;
    KASSERT (modelType.minThreshold > -999999998 &&
             modelType.maxThreshold < 999999998,
             "Under QT threshold model, MIN and MAX of the QT threshold values need to be provided through keywords T_MIN and T_MAX.\n");
  }
  if (modelType.trait == QT) {
    /* threshold value will not be used in any meaningful way, but we will use it for 
     * the loop */
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

  /* Open output files */
  fpHet = fopen (avghetfile, "w");
  KASSERT (fpHet != NULL, "Error in opening file %s for write.\n", avghetfile);
  fprintf (fpHet, "# Version %s\n", programVersion);
  if (modelType.type == TP) {
    fpTP = fopen ("tp.out", "w");
    KASSERT (fpTP != NULL, "Error in opening file %s for write.\n", "tp.out");
    fpPPL = fopen (pplfile, "w");
    fprintf (fpPPL, "# Version %s\n", programVersion);
    KASSERT (fpPPL != NULL, "Error in opening file %s for write.\n", pplfile);
    fprintf (fpPPL, "Chr Marker Position PPL");
    if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
      fprintf (fpPPL, "LD-PPL PPLD");
    }
    fprintf (fpPPL, "\n");
    fflush (fpPPL);
  }
  if (modelOptions.polynomial == TRUE) {
    polynomialInitialization ();
    swLogMsg ("Computation is done in polynomial mode");
  } else {
    swLogMsg ("Computation is done in non-polynomial (direct evaluation) mode");
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
      pTrait->unknownTraitValue = modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN];
      pTrait->lessCutoffFlag = modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED];
      pTrait->moreCutoffFlag = modelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED];
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
    if (pPedigree->currentLoopFlag)
      exitDueToLoop = TRUE;
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

  /* Estimate number of calls to each (appropriate) instance of compute_likelihood for
   * use in progress reporting, and display model information at this point since markers have
   * already been added to locus list */
  swLogMsg (estimateIterations (modelType, modelOptions, modelRange, eCL));

  /* allocate storage for keeping track of het locus in nuclear families */
  allocate_nucfam_het (&pedigreeSet, totalLoci);

  /* initialize some work space */
  initialize_parental_pair_workspace (&parentalPairSpace, originalLocusList.numLocus);
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
    pPedigree->load_flag = 0;   /* Initially 0 and changes to 1 when marker or 
                                 * alternative likelihood values are retrieved */
    pPedigree->polynomialFunction = NULL;
    pPedigree->polynomialFunctionHandle = NULL;
    pPedigree->polynomialFunctionName = NULL;
  }

  /* only for multipoint - we don't handle LD under multipoint yet */
  if (modelType.type == MP) {
    /* allocate space to save temporary results */
    markerNameList = (char **) calloc (sizeof (char *), modelType.numMarkers);
    if (modelType.trait == DT) {

      /* likelihoodDT is for homoLR */
      likelihoodDT = (double **) calloc (sizeof (double *), modelRange.ngfreq);
      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
        /* second dimension is penetrance */
        likelihoodDT[gfreqInd] = (double *) calloc (sizeof (double), modelRange.npenet);
      }
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        /* first dimension is gene freq */
        pPedigree->traitLikelihoodDT = (double **) calloc (sizeof (double *), modelRange.ngfreq);
        pPedigree->alternativeLikelihoodDT = (double **) calloc (sizeof (double *), modelRange.ngfreq);
        for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
          /* second dimension is penetrance */
          pPedigree->traitLikelihoodDT[gfreqInd] = (double *) calloc (sizeof (double), modelRange.npenet);
          pPedigree->alternativeLikelihoodDT[gfreqInd] = (double *) calloc (sizeof (double), modelRange.npenet);
        }
      }
    } else {    /* QT */

      /* first dimension is pedigree */
      likelihoodQT = (double *****) calloc (sizeof (double ****), pedigreeSet.numPedigree + 1);
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree + 1; pedIdx++) {
        /* second dimension is gene freq */
        likelihoodQT[pedIdx] = (double ****) calloc (sizeof (double ***), modelRange.ngfreq);
        for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {

          /* third dimension is mean */
          likelihoodQT[pedIdx][gfreqInd] = (double ***) calloc (sizeof (double **), modelRange.npenet);
          for (penIdx = 0; penIdx < modelRange.npenet; penIdx++) {
            /* fourth dimension is SD */
            likelihoodQT[pedIdx][gfreqInd][penIdx] = (double **) calloc (sizeof (double *), modelRange.nparam);
            for (paramIdx = 0; paramIdx < modelRange.nparam; paramIdx++) {
              /* 5th dimension is threshold */
              likelihoodQT[pedIdx][gfreqInd][penIdx][paramIdx] = (double *) calloc (sizeof (double), modelRange.ntthresh);
            }
          }
        }
      }
    }
  }

  /* find out the max we need to allocate */
  /* after genotype lists have been built, we want to pre-allocate parental pair work space
   * it used to be done dynamically, but it's too costly 
   * parental pair is based on each nuclear family */
  stat_parental_pair_workspace (&pedigreeSet);

  /* after genotype elimination and tracking the max work space needed 
   * for constructing parental pair */
  allocate_parental_pair_workspace (&parentalPairSpace, modelType.numMarkers + 1);

  /* conditional likelihood storage space for each individual */
  allocate_likelihood_space (&pedigreeSet, modelType.numMarkers + 1);

  /* Assume the trait locus is the first one in the list */
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
    savedLocusList.pLocusIndex = (int *) malloc (sizeof (int) * savedLocusList.numLocus);
    for (i = 0; i < 3; i++) {
      savedLocusList.pPrevLocusDistance[i] = (double *) malloc (sizeof (double) * savedLocusList.numLocus);
      savedLocusList.pNextLocusDistance[i] = (double *) malloc (sizeof (double) * savedLocusList.numLocus);
      savedLocusList.pPrevLocusDistance[i][0] = -1;
      savedLocusList.pNextLocusDistance[i][1] = -1;
    }

    if (modelOptions.polynomial == TRUE) {
      status =
        populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0,
                                  -1, -1, 0);
      holdAllPolys ();
    }

    total_count = modelRange.npenet * modelRange.ngfreq * modelRange.nalpha;

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
      if (modelOptions.markerAnalysis != FALSE && pLocus1->locusType != LOCUS_TYPE_MARKER)
        continue;

      for (loc2 = loc1 + 1; loc2 < originalLocusList.numLocus; loc2++) {
        pLocus2 = originalLocusList.ppLocusList[loc2];
        if (pLocus2->locusType != LOCUS_TYPE_MARKER)
          continue;
        savedLocusList.pLocusIndex[1] = loc2;

#ifndef SIMPLEPROGRESS
        fprintf (stdout, "Starting w/loci %s and %s (%d of %d pairs)\n", pLocus1->sName, pLocus2->sName,
                 loc2, originalLocusList.numLocus - 1);
#endif

        /* find out number of alleles this marker locus has */
        if (modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM) {
          /* get the LD parameters */
          pLambdaCell = findLambdas (&modelRange, pLocus1->numOriginalAllele, pLocus2->numOriginalAllele);
          reallocate_LD_loci (pLDLoci, pLocus1->numOriginalAllele, pLocus2->numOriginalAllele);
          pLDLoci->locus1 = loc1;
          pLDLoci->locus2 = loc2;
          pLDLoci->numAllele1 = pLocus1->numOriginalAllele;
          pLDLoci->numAllele2 = pLocus2->numOriginalAllele;
          if (pLocus1->numOriginalAllele == 2 && pLocus2->numOriginalAllele == 2)
            R_square_flag = TRUE;
          else
            R_square_flag = FALSE;
        }

        loopMarkerFreqFlag = 0;
        if (modelRange.nafreq >= 2 && modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM && pLocus2->numOriginalAllele == 2) {
          loopMarkerFreqFlag = 1;
        } else if (modelRange.nafreq == 0) {
          /* add a fake one to facilitate loops and other handlings */
          addAlleleFreq (&modelRange, pLocus2->pAlleleFrequency[0]);
        } else {
          modelRange.nafreq = 1;
          modelRange.afreq[0] = pLocus2->pAlleleFrequency[0];
        }

        /* allocate/initialize result storage */
        initialize_tp_result_storage ();

        /* we will force marker allele frequency loop to execute at least once */
        for (mkrFreqIdx = 0; mkrFreqIdx == 0 || mkrFreqIdx < modelRange.nafreq; mkrFreqIdx++) {
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
          /* Loop over the penetrances, genefrequencies, thetas and call
           * the likelihood calculation, storing each value obtained to
           * disk. */
          for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
            gfreq = modelRange.gfreq[gfreqInd];
            // WHAT ON EARTH IS THIS ALL ABOUT? &&&
            if (1 && modelOptions.markerAnalysis == FALSE) {
              pLocus->pAlleleFrequency[0] = gfreq;
              pLocus->pAlleleFrequency[1] = 1 - gfreq;
              if (modelOptions.polynomial == TRUE);
              else
                update_locus (&pedigreeSet, loc1);
            }

            /* clear Dprime combination impossible flag */
            memset (pLambdaCell->impossibleFlag, 0, sizeof (int) * pLambdaCell->ndprime);
            /* set up haplotype frequencies */
            for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
              if (isDPrime0 (pLambdaCell->lambda[dprimeIdx], pLambdaCell->m, pLambdaCell->n))
                dprime0Idx = dprimeIdx;
              status = setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprimeIdx);
              if (status < 0) {
                pLambdaCell->impossibleFlag[dprimeIdx] = 1;
              }
            }

            if (modelType.trait == DICHOTOMOUS) {

              for (penIdx = 0; penIdx < modelRange.npenet; penIdx++) {
                if (modelOptions.markerAnalysis == FALSE && pLocus1->locusType == LOCUS_TYPE_TRAIT) {
                  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                    pen_DD = modelRange.penet[liabIdx][0][penIdx];
                    pen_Dd = modelRange.penet[liabIdx][1][penIdx];
                    pen_dd = modelRange.penet[liabIdx][2][penIdx];
                    pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
                    pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
                    pTrait->penetrance[2][liabIdx][1][0] = pen_Dd;
                    pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
                    pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
                    pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
                    pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_Dd;
                    pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;
                  }
                  if (modelOptions.polynomial == TRUE);
                  else
                    update_penetrance (&pedigreeSet, traitLocus);
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
                  status =
                    populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2,
                                              initialHetProbAddr, 0, -1, -1, 0);

                /* If we're not on the first iteration, it's not a polynomial build, so
                 * show progress at 1 minute intervals. Have a care to avoid division by zero. */
                strcpy (partialPolynomialFunctionName, "cL0_P%s");
                if (gfreqInd != 0 || penIdx != 0) {
                  swStart (combinedComputeSW);
                  compute_likelihood (&pedigreeSet);
                  cL[0]++;
                  swStop (combinedComputeSW);
                  if (statusRequestSignal) {
                    statusRequestSignal = FALSE;
                    if (cL[0] > 1) {    // The first time thru we have no basis for estimation
                      fprintf (stdout, "%s %d%% complete (~%ld min left)\r",
                               "Calculations", (cL[0] + cL[1]) * 100 / (eCL[0] + eCL[1]),
                               ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                                (eCL[0] + eCL[1]) / (cL[0] + cL[1]) -
                                (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
                      fflush (stdout);
                    }
                  }
                } else { // This _is_ the first iteration
		  swStart (combinedBuildSW);
		  compute_likelihood (&pedigreeSet);
		  cL[0]++;
		  swStop (combinedBuildSW);
		  fprintf (stdout, "%s %d%% complete\r", "Calculations", (cL[0] + cL[1]) * 100 / (eCL[0] + eCL[1]));
		  fflush (stdout);
		}
                if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
                  fprintf (stderr, "Theta 0.5 has likelihood 0\n");
                  fprintf (stderr, "dgf=%f\n", gfreq);
                  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                    pen_DD = modelRange.penet[liabIdx][0][penIdx];
                    pen_Dd = modelRange.penet[liabIdx][1][penIdx];
                    pen_dd = modelRange.penet[liabIdx][2][penIdx];
                    fprintf (stderr, "Liab %d penentrance %f %f %f\n", liabIdx + 1, pen_DD, pen_Dd, pen_dd);
                  }
                  exit (EXIT_FAILURE);
                }
                /* save the results for NULL */
                for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                  /* save the likelihood at null */
                  pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                  pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
                }

                log10_likelihood_null = pedigreeSet.log10Likelihood;
                for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
                  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
                    copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);
                    if (pLambdaCell->impossibleFlag[dprimeIdx] != 0)
                      continue;
                    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
                    copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);
                    /* calculate R square if the marker is a SNP */
                    if (R_square_flag == TRUE)
                      R_square =
                        calculate_R_square (pLocus1->
                                            pAlleleFrequency[0], pLocus2->pAlleleFrequency[0], pLDLoci->ppDValue[0][0]);
                    else
                      R_square = -1;
                  }

                  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++) {
                    if (modelOptions.mapFlag == SA) {
                      theta[0] = modelRange.theta[0][thetaInd];
                      theta[1] = modelRange.theta[1][thetaInd];
                      for (k = 0; k < 3; k++) {
                        locusList->pNextLocusDistance[k][0] = theta[0];
                        locusList->pPrevLocusDistance[k][1] = theta[0];
                      }
                    } else {
                      locusList->pNextLocusDistance[MAP_MALE][0] =
                        locusList->pPrevLocusDistance[MAP_MALE][1] = modelRange.theta[0][thetaInd];
                      locusList->pNextLocusDistance[MAP_FEMALE][0] =
                        locusList->pPrevLocusDistance[MAP_FEMALE][1] = modelRange.theta[1][thetaInd];
                    }

                    if (modelOptions.polynomial == TRUE);
                    else
                      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,
                                                         initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);

                    strcpy (partialPolynomialFunctionName, "cL1_P%s");
                    swStart (combinedComputeSW);
                    compute_likelihood (&pedigreeSet);
                    cL[1]++;
                    swStop (combinedComputeSW);
                    if (statusRequestSignal) {
                      statusRequestSignal = FALSE;
                      if (cL[1] > 1) {  // The first time thru we have no basis for estimation
                        fprintf (stdout, "%s %d%% complete (~%ld min left)\r",
                                 "Calculations", (cL[0] + cL[1]) * 100 / (eCL[0] + eCL[1]),
                                 ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                                  (eCL[0] + eCL[1]) / (cL[0] + cL[1]) -
                                  (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
                        fflush (stdout);
                      }
                    }

                    log10_likelihood_alternative = pedigreeSet.log10Likelihood;
                    if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
                      log10_likelihood_ratio = 0;
                    } else {
                      log10_likelihood_ratio = log10_likelihood_alternative - log10_likelihood_null;
                    }
                    /* check for overflow problem !!! */
                    if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
                      likelihood_ratio = DBL_MAX;
                      tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_total = DBL_MAX;
                    } else
                      /* check for underflow problem too !!! */
                    if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
                      likelihood_ratio = 0;
                    } else {
                      likelihood_ratio = pow (10.0, log10_likelihood_ratio);
                      tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_total += likelihood_ratio;
                    }
                    tp_result[dprimeIdx][thetaInd]
                      [mkrFreqIdx].lr_count++;

                    /* caculating the HET */
                    for (j = 0; j < modelRange.nalpha; j++) {
                      alphaV = modelRange.alpha[j];
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
                      if (log10HetLR >= DBL_MAX_10_EXP - 1) {
                        hetLR = DBL_MAX;
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total = DBL_MAX;
                      } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
                        hetLR = 0;
                      } else {
                        hetLR = pow (10, log10HetLR);
                        if (modelType.ccFlag)
                          /* scale it to prevent overflow */
                          tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total += hetLR / total_count;
                        else
                          tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total += hetLR;
                      }
                      if (tp_result[dprimeIdx][thetaInd]
                          [mkrFreqIdx].max_penIdx < 0 || hetLR > tp_result[dprimeIdx][thetaInd]
                          [mkrFreqIdx].max_lr) {
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_lr = hetLR;
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_alpha = alphaV;
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_gfreq = gfreq;
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_penIdx = penIdx;
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].R_square = R_square;
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_mf = mkrFreq;
                      }
                    }   /* end of calculating HET LR */
                  }     /* end of theta loop */
                }       /* end of D prime loop */
                if (modelOptions.markerAnalysis != FALSE) {
                  /* marker to marker analysis, marker allele frequency is fixed */
                  gfreqInd = modelRange.ngfreq;
                  break;
                }
                if (modelOptions.markerAnalysis != FALSE) {
                  /* marker to marker analysis, penetrance stays at 1 */
                  break;
                }
              } /* end of penetrance loop */
            } /* end of TP */
            else
              /* should be QT or COMBINED - twopoint */
            {
              /* this should be MEAN + SD */
              for (paramIdx = 0; (paramIdx == 0 && modelType.distrib == QT_FUNCTION_CHI_SQUARE)
                   || (modelType.distrib != QT_FUNCTION_CHI_SQUARE && paramIdx < modelRange.nparam); paramIdx++) {
                for (penIdx = 0; penIdx < modelRange.npenet; penIdx++) {
                  breakFlag = FALSE;
                  for (thresholdIdx = 0; thresholdIdx < modelRange.ntthresh; thresholdIdx++) {
                    if (modelOptions.markerAnalysis == FALSE) {
                      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                        mean_DD = modelRange.penet[liabIdx][0][penIdx];
                        mean_Dd = modelRange.penet[liabIdx][1][penIdx];
                        mean_dd = modelRange.penet[liabIdx][2][penIdx];
                        SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
                        SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
                        SD_dd = modelRange.param[liabIdx][2][0][paramIdx];
                        /* threshold for QT */
                        threshold = modelRange.tthresh[liabIdx][thresholdIdx];
                        /* check against the hard coded constraint */
                        if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
                          constraint =
                            (1 - gfreq) * (1 - gfreq) * mean_dd * SD_dd + 2 * gfreq * (1 - gfreq) * mean_Dd * SD_Dd +
                            gfreq * gfreq * mean_DD * SD_DD;
                          if (constraint >= 3.0 || constraint <= -3.0) {
                            breakFlag = TRUE;
                            break;
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

                      } /* liability class Index */
                      if (breakFlag == TRUE)
                        continue;
                      if (modelOptions.polynomial == TRUE);
                      else
                        update_penetrance (&pedigreeSet, traitLocus);
                    }
                    /* marker to marker analysis */
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
                      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,
                                                         initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);
                    KLOG (LOGLIKELIHOOD, LOGDEBUG, "NULL Likelihood\n");
                    strcpy (partialPolynomialFunctionName, "cL2_P%s");
                    compute_likelihood (&pedigreeSet);
                    cL[2]++;

                    if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
                      fprintf (stderr, "Theta 0.5 has likelihood 0\n");
                      fprintf (stderr, "dgf=%f\n", gfreq);
                      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                        pen_DD = modelRange.penet[liabIdx][0][penIdx];
                        pen_Dd = modelRange.penet[liabIdx][1][penIdx];
                        pen_dd = modelRange.penet[liabIdx][2][penIdx];
                        fprintf (stderr, "Liab %d penentrance %f %f %f\n", liabIdx + 1, pen_DD, pen_Dd, pen_dd);
                      }
                      exit (EXIT_FAILURE);
                    }
                    for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                      /* save the likelihood at null */
                      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                      pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
                    }
                    log10_likelihood_null = pedigreeSet.log10Likelihood;
                    for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
                      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
                        copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);
                        if (pLambdaCell->impossibleFlag[dprimeIdx] != 0)
                          continue;
                        copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
                        copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);
                      }
                      for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++) {

                        if (modelOptions.mapFlag == SA) {
                          theta[0] = modelRange.theta[0][thetaInd];
                          theta[1] = modelRange.theta[1][thetaInd];
                          for (k = 0; k < 3; k++) {
                            locusList->pNextLocusDistance[k][0] = theta[0];
                            locusList->pPrevLocusDistance[k][1] = theta[0];
                          }
                        } else {
                          locusList->pNextLocusDistance[MAP_MALE][0] =
                            locusList->pPrevLocusDistance[MAP_MALE][1] = modelRange.theta[0][thetaInd];
                          locusList->pNextLocusDistance[MAP_FEMALE][0] =
                            locusList->pPrevLocusDistance[MAP_FEMALE][1] = modelRange.theta[1][thetaInd];
                        }

                        if (modelOptions.polynomial == TRUE);
                        else
                          status =
                            populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2,
                                                      initialHetProbAddr, 0, -1, -1, 0);

                        KLOG (LOGLIKELIHOOD, LOGDEBUG, "ALT Likelihood\n");
                        strcpy (partialPolynomialFunctionName, "cL3_P%s");
                        compute_likelihood (&pedigreeSet);
                        cL[3]++;
                        log10_likelihood_alternative = pedigreeSet.log10Likelihood;
                        if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
                          log10_likelihood_ratio = 0;
                        } else {
                          log10_likelihood_ratio = log10_likelihood_alternative - log10_likelihood_null;
                        }
                        /* check for overflow problem !!! */
                        if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
                          likelihood_ratio = DBL_MAX;
                          tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_total = DBL_MAX;
                        } else
                          /* check for underflow problem too !!! */
                        if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
                          likelihood_ratio = 0;
                        } else {
                          likelihood_ratio = pow (10.0, log10_likelihood_ratio);
                          tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_total += likelihood_ratio;
                        }
                        tp_result[dprimeIdx][thetaInd]
                          [mkrFreqIdx].lr_count++;
                        /* caculating the HET */
                        for (j = 0; j < modelRange.nalpha; j++) {
                          alphaV = modelRange.alpha[j];
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
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total = DBL_MAX;
                          } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
                            hetLR = 0;
                          } else {
                            adjustedHetLR = hetLR = pow (10, log10HetLR);
                            /* for threshold parameter, we need to make sure the weighting is even */
                            if (1 || modelType.distrib == QT_FUNCTION_CHI_SQUARE) {
                              if (modelRange.ntthresh == 1) {
                                adjustedHetLR *= 2 * (modelType.maxThreshold - modelType.minThreshold);
                              } else if (thresholdIdx == modelRange.ntthresh - 1) {
                                adjustedHetLR *=
                                  (2 * modelType.maxThreshold - threshold - modelRange.tthresh[0][thresholdIdx - 1]);
                              } else if (thresholdIdx == 0) {
                                adjustedHetLR *=
                                  (threshold + modelRange.tthresh[0][thresholdIdx + 1] - 2 * modelType.minThreshold);
                              } else {
                                adjustedHetLR *=
                                  (modelRange.tthresh[0][thresholdIdx + 1] - modelRange.tthresh[0][thresholdIdx - 1]);
                              }
                            }
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total += adjustedHetLR;
                          }
                          if (tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_penIdx < 0
                              || hetLR > tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_lr) {
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_lr = hetLR;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_alpha = alphaV;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_gfreq = gfreq;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_penIdx = penIdx;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_paramIdx = paramIdx;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_thresholdIdx = thresholdIdx;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].R_square = R_square;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_mf = mkrFreq;
                          }
                        }
                      } /* end of theta */
                    }   /* end of D prime */
                    if (modelOptions.markerAnalysis != FALSE)
                      break;
                  }     /* end of threshold loop */
                  if (modelOptions.markerAnalysis != FALSE)
                    break;
                }       /* end of penetrance loop */
                if (modelOptions.markerAnalysis != FALSE)
                  break;
              } /* end of parameter loop */
              if (modelOptions.markerAnalysis != FALSE)
                break;
            }   /* end of QT */
          }     /* end of gene freq */
          /* only loop marker allele frequencies when doing LD */
          if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
            break;
          /* we can only do SNPs when looping over marker allele frequency */
          if (pLocus2->numOriginalAllele > 2)
            break;
        }       /* end of marker allele frequency looping */

        /* calculate the average BR */
        get_average_LR (tp_result);

        /* for each D prime and theta, print out average and maximizing model information - MOD */
        fprintf (fpHet, "# %-d  %s %s \n", loc2, pLocus1->sName, pLocus2->sName);
        fprintf (fpHet, "Chr ");
        if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
          for (i = 0; i < pLocus1->numOriginalAllele - 1; i++)
            for (j = 0; j < pLocus2->numOriginalAllele - 1; j++)
              fprintf (fpHet, "D%1d%1d ", i + 1, j + 1);
        fprintf (fpHet, "Theta(M,F) BayesRatio MOD R2 Alpha DGF MF ");
	for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
	  if (modelType.trait == DT)
	    fprintf (fpHet, "LC%dPV(DD,Dd,dd) ", liabIdx);
	  else
	    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	      fprintf (fpHet, "LC%dPV(DDMean,DdMean,ddMean,DDSD,DdSD,ddSD,Thresh) ", liabIdx);
	    else
	      fprintf (fpHet, "LC%dPV(DDDF,DdDF,ddDF,Thresh) ", liabIdx);
	fprintf (fpHet, "\n");
        for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
          for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++) {
            if (tp_result[dprimeIdx][thetaInd]
                [modelRange.nafreq].lr_count == 0)
              continue;
            theta[0] = modelRange.theta[0][thetaInd];
            theta[1] = modelRange.theta[1][thetaInd];
            max = log10 (tp_result[dprimeIdx][thetaInd]
                         [modelRange.nafreq].max_lr);
            gfreq = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_gfreq;
            alphaV = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_alpha;
            penIdx = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_penIdx;
            paramIdx = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_paramIdx;
            thresholdIdx = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_thresholdIdx;
            R_square = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].R_square;
            fprintf (fpHet, "%d ", pLocus2->pMapUnit->chromosome);
            if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
              for (i = 0; i < pLocus1->numOriginalAllele - 1; i++)
                for (j = 0; j < pLocus2->numOriginalAllele - 1; j++) {
                  fprintf (fpHet, "%.2f ", pLambdaCell->lambda[dprimeIdx][i][j]);
                }
            }
            fprintf (fpHet, "(%.4f,%.4f) %.6e %.4f %.4f %.2f %.4f %.4f",
                     theta[0], theta[1],
                     tp_result[dprimeIdx][thetaInd][modelRange.nafreq].het_lr_avg, max,
                     tp_result[dprimeIdx][thetaInd][modelRange.nafreq].R_square, alphaV, gfreq,
                     tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_mf);
            for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
              pen_DD = modelRange.penet[liabIdx][0][penIdx];
              pen_Dd = modelRange.penet[liabIdx][1][penIdx];
              pen_dd = modelRange.penet[liabIdx][2][penIdx];
	      fprintf (fpHet, " (%.3f,%.3f,%.3f", pen_DD, pen_Dd, pen_dd);
              if (modelType.trait != DT && modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
                SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
                SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
                SD_dd = modelRange.param[liabIdx][2][0][paramIdx];
                fprintf (fpHet, ",%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dd);
              }
              if (modelType.trait != DT) {
                threshold = modelRange.tthresh[liabIdx][thresholdIdx];
                fprintf (fpHet, ",%.3f)", threshold);
              } else
		fprintf (fpHet, ")");
            }
            fprintf (fpHet, "\n");
          }     /* theta loop */
        }       /* dprime loop */
        fprintf (fpTP, "# %-d  %s %s Max Het LR\n", loc2, pLocus2->sName, pLocus1->sName);
        initialFlag = 1;
        max = -99999;
        max_at_theta0 = -99999;
        max_at_dprime0 = -99999;
        for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
          //dprime = pLambdaCell->lambda[dprimeIdx][0][0];
          for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++) {
            theta[0] = modelRange.theta[0][thetaInd];
            theta[1] = modelRange.theta[1][thetaInd];
            lr = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_lr;
            if (initialFlag || lr > max) {
              /* overall max */
              max = lr;
              maxDPrimeIdx = dprimeIdx;
              maxThetaIdx = thetaInd;
            }
            if (initialFlag || (-ERROR_MARGIN <= theta[0] && theta[0] <= ERROR_MARGIN && -ERROR_MARGIN <= theta[1]
                                && theta[1] <= ERROR_MARGIN)) {
              /* find the max for models with theta equal to 0 */
              theta0Idx = thetaInd;
              if (lr > max_at_theta0) {
                max_at_theta0 = lr;
                maxDPrimeIdx_at_theta0 = dprimeIdx;
              }
            }
            if (dprimeIdx == dprime0Idx) {
              if (initialFlag || maxTheta_at_dprime0 < 0 || lr > max_at_dprime0) {
                max_at_dprime0 = lr;
                maxTheta_at_dprime0 = thetaInd;
              }
            }
            initialFlag = 0;
          }
          initialFlag = 0;
        }
        fprintf (fpTP, "Chr     Marker   Position   MOD   DPrime Theta R2 ALPHA DGF MF PEN_DD PEN_Dd PEN_dd\n");
        /* overall maximizing model - MOD */
        fprintf (fpTP, "# Overall MOD maximizing model:\n");
        theta[0] = modelRange.theta[0][maxThetaIdx];
        theta[1] = modelRange.theta[1][maxThetaIdx];
        gfreq = tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].max_gfreq;
        mkrFreq = tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].max_mf;
        alphaV = tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].max_alpha;
        penIdx = tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].max_penIdx;
        R_square = tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].R_square;
        paramIdx = tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].max_paramIdx;
        thresholdIdx = tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].max_thresholdIdx;
        fprintf (fpTP,
                 "%4d %15s %8.4f %8.4f %5.2f (%6.4f, %6.4f) %5.3f %4.2f %6.4f %6.4f",
                 pLocus2->pMapUnit->chromosome, pLocus2->sName,
                 pLocus2->pMapUnit->mapPos[SEX_AVERAGED], log10 (max),
                 pLambdaCell->lambda[maxDPrimeIdx][0][0], theta[0], theta[1], R_square, alphaV, gfreq, mkrFreq);
        for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
          pen_DD = modelRange.penet[liabIdx][0][penIdx];
          pen_Dd = modelRange.penet[liabIdx][1][penIdx];
          pen_dd = modelRange.penet[liabIdx][2][penIdx];
          fprintf (fpTP, " %5.3f %5.3f %5.3f ", pen_DD, pen_Dd, pen_dd);
          if (modelType.trait != DT && modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
            SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
            SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
            SD_dd = modelRange.param[liabIdx][2][0][paramIdx];
            fprintf (fpTP, " %5.3f %5.3f %5.3f ", SD_DD, SD_Dd, SD_dd);
          }
          if (modelType.trait != DT) {
            threshold = modelRange.tthresh[liabIdx][thresholdIdx];
            fprintf (fpTP, " %5.3f ", threshold);
          }
        }
        fprintf (fpTP, "\n");
        fflush (fpTP);

        /* maximizing model at theta equal to 0 - MOD */
        fprintf (fpTP, "# MOD maximizing model for theta=0:\n");
        gfreq = tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.nafreq].max_gfreq;
        mkrFreq = tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.nafreq].max_mf;
        alphaV = tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.nafreq].max_alpha;
        penIdx = tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.nafreq].max_penIdx;
        R_square = tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.nafreq].R_square;
        paramIdx = tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.nafreq].max_paramIdx;
        thresholdIdx = tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.nafreq].max_thresholdIdx;
        fprintf (fpTP,
                 "%4d %15s %8.4f %8.4f %5.2f %6.4f %5.3f %4.2f %6.4f %6.4f",
                 pLocus2->pMapUnit->chromosome, pLocus2->sName,
                 pLocus2->pMapUnit->mapPos[SEX_AVERAGED],
                 log10 (max_at_theta0),
                 pLambdaCell->lambda[maxDPrimeIdx_at_theta0][0][0], 0.0, R_square, alphaV, gfreq, mkrFreq);
        for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
          pen_DD = modelRange.penet[liabIdx][0][penIdx];
          pen_Dd = modelRange.penet[liabIdx][1][penIdx];
          pen_dd = modelRange.penet[liabIdx][2][penIdx];
          fprintf (fpTP, " %5.3f %5.3f %5.3f ", pen_DD, pen_Dd, pen_dd);
          if (modelType.trait != DT && modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
            SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
            SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
            SD_dd = modelRange.param[liabIdx][2][0][paramIdx];
            fprintf (fpTP, " %5.3f %5.3f %5.3f ", SD_DD, SD_Dd, SD_dd);
          }
          if (modelType.trait != DT) {
            threshold = modelRange.tthresh[liabIdx][thresholdIdx];
            fprintf (fpTP, " %5.3f ", threshold);
          }
        }
        fprintf (fpTP, "\n");
        fflush (fpTP);

        /* maximizing model at d prime equal to 0 - MOD */
        fprintf (fpTP, "# MOD maximizing model for dprime=0:\n");
        gfreq = tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].max_gfreq;
        mkrFreq = tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].max_mf;
        alphaV = tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].max_alpha;
        penIdx = tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].max_penIdx;
        R_square = tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].R_square;
        paramIdx = tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].max_paramIdx;
        thresholdIdx = tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].max_thresholdIdx;
        fprintf (fpTP,
                 "%4d %15s %8.4f %8.4f %5.2f %6.4f %5.3f %4.2f %6.4f %6.4f ",
                 pLocus2->pMapUnit->chromosome, pLocus2->sName,
                 pLocus2->pMapUnit->mapPos[SEX_AVERAGED], log10 (max_at_dprime0), 0.0, 0.0, R_square, alphaV, gfreq, mkrFreq);
        for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
          pen_DD = modelRange.penet[liabIdx][0][penIdx];
          pen_Dd = modelRange.penet[liabIdx][1][penIdx];
          pen_dd = modelRange.penet[liabIdx][2][penIdx];
          fprintf (fpTP, " %5.3f %5.3f %5.3f ", pen_DD, pen_Dd, pen_dd);
          if (modelType.trait != DT && modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
            SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
            SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
            SD_dd = modelRange.param[liabIdx][2][0][paramIdx];
            fprintf (fpTP, " %5.3f %5.3f %5.3f ", SD_DD, SD_Dd, SD_dd);
          }
          if (modelType.trait != DT) {
            threshold = modelRange.tthresh[liabIdx][thresholdIdx];
            fprintf (fpTP, " %5.3f ", threshold);
          }
        }
        fprintf (fpTP, "\n");
        fflush (fpTP);

        /* find the overall maximizing theta and dprime - LR
         * with the other parameter integrated out */
        max = -9999.99;
        max_at_dprime0 = -9999.99;
        max_at_theta0 = -9999.99;
        for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
          //              dprime = pLambdaCell->lambda[dprimeIdx][0][0];
          for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++) {
            theta[0] = modelRange.theta[0][thetaInd];
            theta[1] = modelRange.theta[1][thetaInd];
            lr = tp_result[dprimeIdx][thetaInd][0].het_lr_avg;
            if (lr > max) {
              max = lr;
              maxThetaIdx = thetaInd;
              maxDPrimeIdx = dprimeIdx;
            }
            if (-ERROR_MARGIN <= theta[0]
                && theta[0] <= ERROR_MARGIN && -ERROR_MARGIN <= theta[1]
                && theta[1] <= ERROR_MARGIN) {
              if (lr > max_at_theta0) {
                max_at_theta0 = lr;
                maxDPrimeIdx_at_theta0 = dprimeIdx;
              }
            }
            if (dprime0Idx == dprimeIdx) {
              if (lr > max_at_dprime0) {
                max_at_dprime0 = lr;
                maxTheta_at_dprime0 = thetaInd;
              }
            }
          }
        }
        /* overall maximizing model - LR */
        fprintf (fpTP, "# Overall LR maximizing model:\n");
        theta[0] = modelRange.theta[0][maxThetaIdx];
        theta[1] = modelRange.theta[1][maxThetaIdx];
        //gfreq = tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].max_gfreq;
        fprintf (fpTP,
                 "%4d %15s %8.4f %8.4f %5.2f (%6.4f %6.4f)\n",
                 pLocus2->pMapUnit->chromosome, pLocus2->sName,
                 pLocus2->pMapUnit->mapPos[SEX_AVERAGED], log10 (max),
                 pLambdaCell->lambda[maxDPrimeIdx][0][0], theta[0], theta[1]);
        fflush (fpTP);

        /* maximizing model at theta equal to 0 - LR */
        fprintf (fpTP, "# LR maximizing model for theta (0, 0):\n");
        fprintf (fpTP,
                 "%4d %15s %8.4f %8.4f %5.2f %6.4f\n",
                 pLocus2->pMapUnit->chromosome, pLocus2->sName,
                 pLocus2->pMapUnit->mapPos[SEX_AVERAGED],
                 log10 (max_at_theta0), pLambdaCell->lambda[maxDPrimeIdx_at_theta0][0][0], 0.0);
        fflush (fpTP);

        /* maximizing model at d prime equal to 0 - LR */
        fprintf (fpTP, "# LR maximizing model for dprime=0:\n");
        fprintf (fpTP,
                 "%4d %15s %8.4f %8.4f %5.2f %6.4f\n",
                 pLocus2->pMapUnit->chromosome, pLocus2->sName,
                 pLocus2->pMapUnit->mapPos[SEX_AVERAGED], log10 (max_at_dprime0), 0.0, 0.0);
        fflush (fpTP);

        /* output PPL now */
        /* chromosome, marker name, position, PPL */
        ppl = calculate_PPL (tp_result[dprime0Idx]);
        fprintf (fpPPL, "%d %s %.4f %.*f ",
                 pLocus2->pMapUnit->chromosome, pLocus2->sName, pLocus2->pMapUnit->mapPos[SEX_AVERAGED],
                 ppl >= .025 ? 2 : 3, ppl >= .025 ? rint (ppl * 100.) / 100. : rint (ppl * 1000.) / 1000.);
        fflush (fpPPL);
        /* output LD-PPL now if needed */
        if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
          /* calculate the LD LR average first */
          get_average_LD_LR (tp_result);
          /* calculate the LD-PPL - posterior probability of linkage allowing for LD */
          ldppl = calculate_PPL (tp_result[pLambdaCell->ndprime]);
          /* now calculate the PPLD - posterior probability of LD given linkage */
          ppld = calculate_PPLD (tp_result);
          fprintf (fpPPL, "%.*f %.*f ",
                   ldppl >= .025 ? 2 : 3, ldppl >= .025 ? rint (ldppl * 100.) / 100. : rint (ldppl * 1000.) / 1000.,
                   ppld >= .025 ? 2 : 3, ppld >= .025 ? rint (ppld * 100.) / 100. : rint (ppld * 1000.) / 1000.);
        }
        fprintf (fpPPL, "\n");
        fflush (fpPPL);

        prevNumDPrime = pLambdaCell->ndprime;
        /* need to clear polynomial */

        if (modelOptions.polynomial == TRUE && modelType.ccFlag == 0) {
          /* under case ctrl we don't clear up the polynomial */
          pedigreeSetPolynomialClearance (&pedigreeSet);
        }


        if (modelOptions.markerAnalysis == ADJACENTMARKER)
          loc2 = originalLocusList.numLocus;

#ifndef SIMPLEPROGRESS
        fprintf (stdout, "\n");
#endif
      } /* end of looping second locus - loc2 */
      /* if we are doing trait marker, then we are done */
      /* Used to read: modelOptions.markerToMarker != TRUE which
       * is the same as markerAnalysis == FALSE as long as the old
       * markerToMarker and adjacentMarker flags were truly
       * orthogonal. Otherwise, it should be markerAnalysis !=
       * ADJACENTMARKER. */
      if (modelOptions.markerAnalysis == FALSE)
        loc1 = originalLocusList.numLocus;
    }   /* end of looping first locus - loc1 */
    /* free two point result storage */
    free_tp_result_storage (prevNumDPrime);
  } /* end of two point */
  else {        /* multipoint */

    /* marker set locus list for each position */
    markerLocusList.maxNumLocus = modelType.numMarkers;
    markerLocusList.numLocus = modelType.numMarkers;
    markerLocusList.traitOrigLocus = -1;
    markerLocusList.traitLocusIndex = -1;
    markerLocusList.pLocusIndex = (int *) calloc (markerLocusList.maxNumLocus, sizeof (int));
    for (k = 0; k < 3; k++) {
      markerLocusList.pPrevLocusDistance[k] = (double *) calloc (markerLocusList.maxNumLocus, sizeof (double));
      markerLocusList.pNextLocusDistance[k] = (double *) calloc (markerLocusList.maxNumLocus, sizeof (double));
    }

    /* assuming we always have trait in the analysis - this may not be true 
     * need to add code to process marker to marker analysis under multipoin
     */
    savedLocusList.numLocus = modelType.numMarkers + 1;
    savedLocusList.maxNumLocus = modelType.numMarkers + 1;
    savedLocusList.pLocusIndex = (int *) calloc (savedLocusList.maxNumLocus, sizeof (int));
    for (k = 0; k < 3; k++) {
      savedLocusList.pPrevLocusDistance[k] = (double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
      savedLocusList.pNextLocusDistance[k] = (double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
    }

    /* Allocate storage to calculate the trait likelihood independent of the trait position */
    traitLocusList.numLocus = 1;
    traitLocusList.maxNumLocus = 1;
    traitLocusList.traitLocusIndex = 0;
    traitLocusList.traitOrigLocus = traitLocus;
    traitLocusList.pLocusIndex = (int *) calloc (traitLocusList.maxNumLocus, sizeof (int));
    traitLocusList.pLocusIndex[0] = 0;
    for (k = 0; k < 3; k++) {
      traitLocusList.pPrevLocusDistance[k] = (double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
      traitLocusList.pNextLocusDistance[k] = (double *) calloc (savedLocusList.maxNumLocus, sizeof (double));

      traitLocusList.pPrevLocusDistance[k][0] = -1;
      traitLocusList.pNextLocusDistance[k][0] = -1;
    }
    /* populate the trait xmission matrix */
    locusList = &traitLocusList;
    xmissionMatrix = traitMatrix;
    status = populate_xmission_matrix (traitMatrix, 1, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);
    if (modelOptions.polynomial == TRUE)
      holdAllPolys ();

    /* For trait likelihood */
#ifndef SIMPLEPROGRESS
    fprintf (stdout, "Determining trait likelihood...\n");
#else
    fprintf (stdout, "Calculations 0%% complete\r");
    fflush (stdout);
#endif

    locusList = &traitLocusList;
    xmissionMatrix = traitMatrix;
    if (pTrait->type == DICHOTOMOUS) {
      /* load all saved trait likelihood */
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        if (modelOptions.saveResults == TRUE) {
          pPedigree->load_flag = restoreTrait (modelOptions.sexLinked, pPedigree->sPedigreeID, pPedigree->traitLikelihoodDT);
        } else
          pPedigree->load_flag = 0;
      }

      for (penIdx = 0; (penIdx == 0) || (modelOptions.dryRun == 0 && penIdx < modelRange.npenet); penIdx++) {
        for (liabIdx = 0; (liabIdx == 0) || (modelOptions.dryRun == 0 && liabIdx < modelRange.nlclass); liabIdx++) {
          pen_DD = modelRange.penet[liabIdx][0][penIdx];
          pen_Dd = modelRange.penet[liabIdx][1][penIdx];
          pen_dd = modelRange.penet[liabIdx][2][penIdx];
          pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
          pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
          pTrait->penetrance[2][liabIdx][1][0] = pen_Dd;
          pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
          pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
          pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
          pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_Dd;
          pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;
        }

        if (modelOptions.polynomial == TRUE);
        else
          /* only need to update trait locus */
          update_penetrance (&pedigreeSet, traitLocus);

        /* Iterate over gene frequencies, but only one time thru if doing a dry-run. */
        for (gfreqInd = 0; (gfreqInd == 0) || (modelOptions.dryRun == 0 && gfreqInd < modelRange.ngfreq); gfreqInd++) {

          /* updated trait locus allele frequencies */
          gfreq = modelRange.gfreq[gfreqInd];
          pLocus->pAlleleFrequency[0] = gfreq;
          pLocus->pAlleleFrequency[1] = 1 - gfreq;

          if (modelOptions.polynomial == TRUE);
          else
            update_locus (&pedigreeSet, traitLocus);

          /* Compute the likelihood for the trait */
          sprintf (partialPolynomialFunctionName, "T_P%%sSL%d", modelOptions.sexLinked);
          compute_likelihood (&pedigreeSet);
          cL[4]++;
#ifndef SIMPLEPROGRESS
          if (cL[4] % MAX (1, eCL[4] / 5) == 1) {
            fprintf (stdout, "Trait likelihood evaluations %d%% complete\r", cL[4] * 100 / eCL[4]);
            fflush (stdout);
          }
#endif
          if (modelOptions.dryRun != 0)
            continue;

          if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
            fprintf (stderr, "Trait has likelihood 0\n");
            fprintf (stderr, "dgf=%f\n", gfreq);
            for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
              pen_DD = modelRange.penet[liabIdx][0][penIdx];
              pen_Dd = modelRange.penet[liabIdx][1][penIdx];
              pen_dd = modelRange.penet[liabIdx][2][penIdx];
              fprintf (stderr, "Liab %d penentrance %f %f %f\n", liabIdx + 1, pen_DD, pen_Dd, pen_dd);
            }
            exit (EXIT_FAILURE);
          }
          /* save the results for NULL */
          for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
            /* save the likelihood at null */
            pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
            if (pPedigree->load_flag == 0) {    /*update only for the pedigrees which were add for this run */
              pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
              pPedigree->traitLikelihoodDT[gfreqInd][penIdx] = pPedigree->likelihood;
            }
          }

          log10_likelihood_null = pedigreeSet.log10Likelihood;
          likelihoodDT[gfreqInd][penIdx] = log10_likelihood_null;
        }       /* gfreq */
      } /* pen */

      /* save all  trait likelihood which were created in this run */
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        /* save the likelihood at null */
        pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        if ((modelOptions.saveResults == TRUE) && (pPedigree->load_flag == 0)) {        /*save only for the pedigrees which were add for this run */
          pPedigree->load_flag = saveTrait (modelOptions.sexLinked, pPedigree->sPedigreeID, pPedigree->traitLikelihoodDT);
        } else {
          pPedigree->load_flag = 0;
        }
      }
    } else
      /* multipoint QT or COMBINED */
    {
      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
        gfreq = modelRange.gfreq[gfreqInd];
        pLocus->pAlleleFrequency[0] = gfreq;
        pLocus->pAlleleFrequency[1] = 1 - gfreq;

        update_locus (&pedigreeSet, traitLocus);
        /* this should be MEAN + SD */
        for (paramIdx = 0; paramIdx < modelRange.nparam; paramIdx++) {
          for (penIdx = 0; penIdx < modelRange.npenet; penIdx++) {
            breakFlag = FALSE;
            for (thresholdIdx = 0; thresholdIdx < modelRange.ntthresh; thresholdIdx++) {
              for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                mean_DD = modelRange.penet[liabIdx][0][penIdx];
                mean_Dd = modelRange.penet[liabIdx][1][penIdx];
                mean_dd = modelRange.penet[liabIdx][2][penIdx];
                SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
                SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
                SD_dd = modelRange.param[liabIdx][2][0][paramIdx];
                threshold = modelRange.tthresh[liabIdx][thresholdIdx];

                /* check against the hard coded constraint */
                if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
                  constraint =
                    (1 - gfreq) * (1 -
                                   gfreq) * mean_dd *
                    SD_dd + 2 * gfreq * (1 - gfreq) * mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
                  /*      fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
                   * constraint, gfreq, mean_DD, SD_DD, 
                   * mean_Dd, SD_DD, 
                   * mean_dd, SD_dd);
                   */
                  if (constraint >= 3.0 || constraint <= -3.0) {
                    breakFlag = TRUE;
                    break;
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

              } /* liability class Index */
              if (breakFlag == TRUE)
                continue;
              if (modelOptions.polynomial == TRUE);
              else
                update_penetrance (&pedigreeSet, traitLocus);
              sprintf (partialPolynomialFunctionName, "T_P%%sSL%d", modelOptions.sexLinked);
              compute_likelihood (&pedigreeSet);
              cL[5]++;
#ifndef SIMPLEPROGRESS
              if (cL[5] % MAX (1, eCL[5] / 5) == 1) {
                fprintf (stdout, "Trait likelihood evaluations %d%% complete\r", cL[5] * 100 / eCL[5]);
                fflush (stdout);
              }
#endif
              if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
                fprintf (stderr, "Trait has likelihood 0\n");
                fprintf (stderr, "dgf=%f\n", gfreq);
                for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                  pen_DD = modelRange.penet[liabIdx][0][penIdx];
                  pen_Dd = modelRange.penet[liabIdx][1][penIdx];
                  pen_dd = modelRange.penet[liabIdx][2][penIdx];
                  fprintf (stderr, "Liab %d penentrance %f %f %f\n", liabIdx + 1, pen_DD, pen_Dd, pen_dd);
                }
                exit (EXIT_FAILURE);
              }

              for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                /* save the likelihood at null */
                pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
                likelihoodQT[pedIdx][gfreqInd][penIdx]
                  [paramIdx][thresholdIdx] = pPedigree->likelihood;
              }

              log10_likelihood_null = pedigreeSet.log10Likelihood;
              if (isnan (log10_likelihood_null))
                fprintf (stderr, "Trait likelihood is NAN.\n");
              likelihoodQT[pedigreeSet.numPedigree][gfreqInd][penIdx]
                [paramIdx][thresholdIdx] = log10_likelihood_null;
            }   /* thresholdIdx */
          }     /* penIdx */
        }       /* paramIdx */
      } /* gfreq */
    }   /* end of QT */

#ifndef SIMPLEPROGRESS
    fprintf (stderr, "Trait likelihood evaluations 100%% complete\n");
#endif

    if (modelOptions.polynomial == TRUE) {
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        /* save the likelihood at trait */
        pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        pPedigree->traitLikelihoodPolynomial = pPedigree->likelihoodPolynomial;
        pPedigree->traitLikelihoodPolyList = pPedigree->likelihoodPolyList;
      }
    }

    /* print out some statistics under dry run */
    if (modelOptions.dryRun != 0) {
      print_dryrun_stat (&pedigreeSet, -1);
    }

    /* get the trait locations we need to evaluate at */
    numPositions = modelRange.ntloc;
    mp_result = (SUMMARY_STAT *) calloc (numPositions, sizeof (SUMMARY_STAT));
    /* Need to output the results */
    fprintf (fpHet, "Chr Position PPL BayesRatio MOD Alpha DGF ");
    for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
      if (modelType.trait == DT)
        fprintf (fpHet, "LC%dPV(DD,Dd,dd) ", liabIdx);
      else
	if (modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	  fprintf (fpHet, "LC%dPV(DDMean,DdMean,ddMean,DDSD,DdSD,ddSD,Thresh) ", liabIdx);
	else
	  fprintf (fpHet, "LC%dPV(DDDF,DdDF,ddDF,Thresh) ", liabIdx);
    fprintf (fpHet, "MarkerList(0");
    for (k = 1; k < modelType.numMarkers; k++)
      fprintf (fpHet, ",%d", k);
    fprintf (fpHet, ")\n");

    prevFirstMarker = -1;
    prevLastMarker = -1;
    prevTraitInd = -1;
    leftMarker = -1;

    /* Iterate over all positions in the analysis. */
    for (posIdx = 0; posIdx < numPositions; posIdx++) {
      /* positions listed are sex average positions */
      traitPos = modelRange.tloc[posIdx];
      /* Set the sex-averaged position first. The sex-specific positions will be updated 
       * once markers are selected since interpolation might be needed. */
      pTraitLocus->mapPosition[0] = traitPos;
      pTraitLocus->mapPosition[1] = traitPos;
      pTraitLocus->mapPosition[2] = traitPos;
      /* initialize the locusList */
      locusList = &savedLocusList;
      memset (locusList->pLocusIndex, 0, sizeof (int) * locusList->maxNumLocus);
      for (k = 0; k < 3; k++) {
        memset (&locusList->pPrevLocusDistance[k][0], 0, sizeof (double) * locusList->maxNumLocus);
        memset (&locusList->pNextLocusDistance[k][0], 0, sizeof (double) * locusList->maxNumLocus);
      }
      locusList->numLocus = 1;
      locusList->pLocusIndex[0] = traitLocus;
      for (k = 0; k < 3; k++) {
        locusList->pPrevLocusDistance[k][0] = -1;
        locusList->pNextLocusDistance[k][0] = -1;
      }
      /* select markers to be used for the multipoint analysis */
      add_markers_to_locuslist (locusList, modelType.numMarkers, &leftMarker, 0, originalLocusList.numLocus - 1, traitPos, 0);
      /* store the markers used */
      mp_result[posIdx].pMarkers = (int *) calloc (modelType.numMarkers, sizeof (int));
      k = 0;    /* marker index */
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

#ifndef SIMPLEPROGRESS
      /* Say where we're at with the trait locus and markers. */
      fprintf (stdout, "Starting w/trait locus at %.2f (%d/%d positions) with", traitPos, posIdx + 1, numPositions);
#endif

      markerSetChanged = FALSE;
      if (prevFirstMarker != mp_result[posIdx].pMarkers[0] ||
          prevLastMarker != mp_result[posIdx].pMarkers[modelType.numMarkers - 1]) {
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
              cm_to_recombination_fraction (currPos[j] - prevPos[j], map.mapFunction);
          }
          if (modelOptions.mapFlag == SA) {
            for (j = 1; j <= 2; j++) {
              markerLocusList.pPrevLocusDistance[j][k] =
                markerLocusList.pNextLocusDistance[j][k - 1] = markerLocusList.pPrevLocusDistance[0][k];
            }
          }
          prevPos = currPos;
        }       /* end of loop over the markers to set up locus list */

#ifndef SIMPLEPROGRESS
        fprintf (stdout, " new markers");
        for (k = 0; k < modelType.numMarkers; k++)
          fprintf (stdout, " %d(%.2f)", markerLocusList.pLocusIndex[k], *get_map_position (markerLocusList.pLocusIndex[k]));
        fprintf (stdout, "\n");

        /* Calculate likelihood for the marker set */
        fprintf (stdout, "Determining marker set likelihood...\n");
#endif

        locusList = &markerLocusList;
        xmissionMatrix = markerMatrix;
        if (modelOptions.polynomial == TRUE) {
          pedigreeSetPolynomialClearance (&pedigreeSet);
        }

        /* save the polynomial flag */
        polynomialFlag = modelOptions.polynomial;
        status = populate_xmission_matrix (markerMatrix, markerLocusList.numLocus, initialProbAddr, initialProbAddr2,
                                           initialHetProbAddr, 0, -1, -1, 0);
        if (modelOptions.polynomial == TRUE)
          freePolys ();

        print_xmission_matrix (markerMatrix, markerLocusList.numLocus, 0, 0, tmpID);
        for (k = 0; k < modelType.numMarkers; k++) {
          markerNameList[k] = (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[k]])->sName;
        }
        for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
          /* save the marker likelihood   */
          pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
          if (modelOptions.saveResults == TRUE) {
            pPedigree->load_flag =
              restoreMarker (pPedigree->sPedigreeID,
                             (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
                             modelType.numMarkers, markerNameList, &(pPedigree->markerLikelihood));
          } else {
            pPedigree->load_flag = 0;
          }
        }
        sprintf (partialPolynomialFunctionName, "ML_P%%sC%dFM%dof%d",
                 (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
                 mp_result[posIdx].pMarkers[0], modelType.numMarkers);
        compute_likelihood (&pedigreeSet);
        cL[6]++;
#ifndef SIMPLEPROGRESS
        fprintf (stdout, "Marker set likelihood evaluations %d%% complete...\n",
                 MAX (cL[6] * 100 / eCL[6], (posIdx + 1) * 100 / numPositions));
#endif

        modelOptions.polynomial = polynomialFlag;

        /* print out some statistics under dry run */
        if (modelOptions.dryRun != 0) {
          print_dryrun_stat (&pedigreeSet, -1);
        } else {
          /* save the results for marker likelihood */
          for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
            /* save the likelihood at null */
            pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];

            //      fprintf(stderr, "pedIdx=%d  markerpediLikehood %G\n", pedIdx, pPedigree->likelihood);
            if (modelOptions.saveResults == TRUE) {
              if (pPedigree->load_flag == 0) {  /*save only for the pedigrees which were add for this run */
                pPedigree->markerLikelihood = pPedigree->likelihood;
                pPedigree->load_flag =
                  saveMarker (pPedigree->sPedigreeID,
                              (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
                              modelType.numMarkers, markerNameList, &(pPedigree->markerLikelihood));
              }
            } else {
              pPedigree->markerLikelihood = pPedigree->likelihood;
            }
            pPedigree->load_flag = 0;
          }
          pedigreeSet.log10MarkerLikelihood = pedigreeSet.log10Likelihood;
        }
      } /* end of marker set change */
      else
#ifndef SIMPLEPROGRESS
        fprintf (stdout, " same markers\n");
#else
        ;
#endif
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
            pTraitLocus->mapPosition[MAP_MALE] = pTraitLocus->mapPosition[MAP_FEMALE] = 0;
          } else {
            /* get the relative position on the sex average map */
            relativePos = traitPos / marker1Pos[0];
            pTraitLocus->mapPosition[MAP_MALE] = relativePos * marker1Pos[MAP_MALE];
            pTraitLocus->mapPosition[MAP_FEMALE] = relativePos * marker1Pos[MAP_FEMALE];
          }
          /* update the inter locus distance - sex averaged already done before */
          for (k = 1; k < 3; k++) {
            locusList->pNextLocusDistance[k][0] =
              locusList->pPrevLocusDistance[k][1] =
              cm_to_recombination_fraction (marker1Pos[k] - pTraitLocus->mapPosition[k], map.mapFunction);
          }
        } else if (traitIndex == modelType.numMarkers) {
          /* trait is the last one in the list */
          marker1Pos = get_map_position (locusList->pLocusIndex[modelType.numMarkers - 2]);
          marker2Pos = get_map_position (locusList->pLocusIndex[modelType.numMarkers - 1]);
          /* get the relative position on the sex average map */
          dist = marker2Pos[0] - marker1Pos[0];
          if (dist > ERROR_MARGIN) {
            relativePos = (traitPos - marker2Pos[0]) / dist;
            pTraitLocus->mapPosition[MAP_MALE] =
              relativePos * (marker2Pos[MAP_MALE] - marker1Pos[MAP_MALE]) + marker2Pos[MAP_MALE];
            pTraitLocus->mapPosition[MAP_FEMALE] =
              relativePos * (marker2Pos[MAP_FEMALE] - marker1Pos[MAP_FEMALE]) + marker2Pos[MAP_FEMALE];
          } else {
            pTraitLocus->mapPosition[MAP_MALE] = traitPos - marker2Pos[0] + marker2Pos[MAP_MALE];
            pTraitLocus->mapPosition[MAP_FEMALE] = traitPos - marker2Pos[0] + marker2Pos[MAP_FEMALE];
          }

          /* update the inter locus distance - sex averaged already done before */
          for (k = 1; k <= 2; k++) {
            locusList->pNextLocusDistance[k][traitIndex - 1] =
              locusList->pPrevLocusDistance[k][traitIndex] =
              cm_to_recombination_fraction (pTraitLocus->mapPosition[k] - marker2Pos[k], map.mapFunction);
          }

        } else {
          /* trait is in between two markers */
          marker1Pos = get_map_position (locusList->pLocusIndex[traitIndex - 1]);
          marker2Pos = get_map_position (locusList->pLocusIndex[traitIndex + 1]);
          /* get the relative position on the sex average map */
          dist = marker2Pos[0] - marker1Pos[0];
          if (dist > ERROR_MARGIN) {
            relativePos = (traitPos - marker1Pos[0]) / dist;
            pTraitLocus->mapPosition[MAP_MALE] =
              relativePos * (marker2Pos[MAP_MALE] - marker1Pos[MAP_MALE]) + marker1Pos[MAP_MALE];
            pTraitLocus->mapPosition[MAP_FEMALE] =
              relativePos * (marker2Pos[MAP_FEMALE] - marker1Pos[MAP_FEMALE]) + marker1Pos[MAP_FEMALE];
          } else {
            pTraitLocus->mapPosition[MAP_MALE] = marker1Pos[MAP_MALE];
            pTraitLocus->mapPosition[MAP_FEMALE] = marker1Pos[MAP_FEMALE];
          }
          /* update the inter locus distance - sex averaged already done before */
          for (k = 1; k < 3; k++) {
            locusList->pNextLocusDistance[k][traitIndex - 1] =
              locusList->pPrevLocusDistance[k][traitIndex] =
              cm_to_recombination_fraction (pTraitLocus->mapPosition[k] - marker1Pos[k], map.mapFunction);
            locusList->pNextLocusDistance[k][traitIndex] =
              locusList->pPrevLocusDistance[k][traitIndex + 1] =
              cm_to_recombination_fraction (marker2Pos[k] - pTraitLocus->mapPosition[k], map.mapFunction);
          }
        }
      }

      /* the locus list has been built, go on to the analysis 
       * multipoint DT */
      if (markerSetChanged || locusListChanged) {
        if (modelOptions.polynomial == TRUE) {
          pedigreeSetPolynomialClearance (&pedigreeSet);
          status =
            populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr,
                                      0, -1, -1, 0);
          print_xmission_matrix (altMatrix, savedLocusList.numLocus, 0, 0, tmpID);
          if (modelOptions.polynomial == TRUE)
            freePolys ();
        }
      }

      if (modelOptions.polynomial != TRUE)
        status =
          populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);

      /* For alternative */
#ifndef SIMPLEPROGRESS
      fprintf (stdout, "Determining combined likelihood...\n");
#endif

      if (pTrait->type == DICHOTOMOUS) {
        for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
          pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
          /* load stored alternative likelihood if they were already stored */
          if (modelOptions.saveResults == TRUE)
            pPedigree->load_flag =
              restoreAlternative (pPedigree->sPedigreeID,
                                  (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
                                  traitPos, pPedigree->alternativeLikelihoodDT);
          else
            pPedigree->load_flag = 0;
        }

        for (penIdx = 0; (penIdx == 0) || (modelOptions.dryRun == 0 && penIdx < modelRange.npenet); penIdx++) {
          for (liabIdx = 0; (liabIdx == 0) || (modelOptions.dryRun == 0 && liabIdx < modelRange.nlclass); liabIdx++) {
            pen_DD = modelRange.penet[liabIdx][0][penIdx];
            pen_Dd = modelRange.penet[liabIdx][1][penIdx];
            pen_dd = modelRange.penet[liabIdx][2][penIdx];
            pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
            pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
            pTrait->penetrance[2][liabIdx][1][0] = pen_Dd;
            pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
            pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
            pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
            pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_Dd;
            pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;
          }

          if (modelOptions.polynomial != TRUE)
            update_penetrance (&pedigreeSet, traitLocus);       // Only need to update trait locus

          /* Iterate over gene frequencies -- just one loop for dry-runs. */
          for (gfreqInd = 0; (gfreqInd == 0) || (modelOptions.dryRun == 0 && gfreqInd < modelRange.ngfreq); gfreqInd++) {

            /* Updated trait locus allele frequencies */
            gfreq = modelRange.gfreq[gfreqInd];
            pLocus->pAlleleFrequency[0] = gfreq;
            pLocus->pAlleleFrequency[1] = 1 - gfreq;

            if (modelOptions.polynomial != TRUE)
              update_locus (&pedigreeSet, traitLocus);

            /* If we're not on the first iteration, it's not a polynomial build, so
             * show progress at 1 minute intervals. Have a care to avoid division by zero. */
            sprintf (partialPolynomialFunctionName, "CL_P%%sC%dFM%dof%d",
                     (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
                     mp_result[posIdx].pMarkers[0], modelType.numMarkers);
            if (gfreqInd != 0 || penIdx != 0) {
              swStart (combinedComputeSW);
              compute_likelihood (&pedigreeSet);
              cL[7]++;
              swStop (combinedComputeSW);
              if (statusRequestSignal) {
                statusRequestSignal = FALSE;
                if (cL[7] > 1) {        // The first time thru we have no basis for estimation
#ifndef SIMPLEPROGRESS
                  fprintf (stdout, "%s %d%% complete (~%ld min left)\r",
                           "Combined likelihood evaluations", cL[7] * 100 / eCL[7],
                           ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                            eCL[7] / cL[7] - (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
#else
                  fprintf (stdout, "%s %d%% complete (~%ld min left)\r",
                           "Calculations", (cL[6] + cL[7]) * 100 / (eCL[6] + eCL[7]),
                           ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                            (eCL[6] + eCL[7]) / (cL[6] + cL[7]) -
                            (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
#endif
                  fflush (stdout);
                }
              }
            } else {     // This _is_ the first iteration
              swStart (combinedBuildSW);
              compute_likelihood (&pedigreeSet);
              cL[7]++;
              swStop (combinedBuildSW);
#ifndef SIMPLEPROGRESS
              fprintf (stdout, "%s %d%% complete\r", "Combined likelihood evaluations", cL[7] * 100 / eCL[7]);
#else
              fprintf (stdout, "%s %d%% complete\r", "Calculations", (cL[6] + cL[7]) * 100 / (eCL[6] + eCL[7]));
#endif
              fflush (stdout);
            }
            /* print out some statistics under dry run */
            if (modelOptions.dryRun != 0) {
              print_dryrun_stat (&pedigreeSet, traitPos);
            } else {

              log10_likelihood_alternative = pedigreeSet.log10Likelihood;
              if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99)
                log10_likelihood_ratio = 0;
              else
                log10_likelihood_ratio =
                  log10_likelihood_alternative - likelihoodDT[gfreqInd][penIdx] - pedigreeSet.log10MarkerLikelihood;
              /* check for overflow problem !!! */
              if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
                likelihood_ratio = DBL_MAX;
                mp_result[posIdx].lr_total += DBL_MAX;
              } else
                /* check for underflow problem too !!! */
              if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
                likelihood_ratio = 0;
              } else {
                likelihood_ratio = pow (10.0, log10_likelihood_ratio);
                mp_result[posIdx].lr_total += likelihood_ratio;
              }
              /* add the result to the right placeholder */
              mp_result[posIdx].lr_count++;
              for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                if (pPedigree->load_flag == 0) {
                  pPedigree->alternativeLikelihoodDT[gfreqInd]
                    [penIdx] = pPedigree->likelihood;
                }
              }
              /* caculating the Het */
              for (j = 0; j < modelRange.nalpha; j++) {
                alphaV = modelRange.alpha[j];
                alphaV2 = 1 - alphaV;
                if (alphaV2 < 0)
                  alphaV2 = 0;
                log10HetLR = 0;
                for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                  pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                  homoLR = pPedigree->alternativeLikelihoodDT[gfreqInd]
                    [penIdx] / (pPedigree->traitLikelihoodDT[gfreqInd][penIdx] * pPedigree->markerLikelihood);
                  /*              if (homoLR > 1.0e40 || homoLR < 1.0e-40) {
                   * fprintf(stderr, "homoLR %G, alt %G, trait %G, mrk %G\n",
                   * homoLR, pPedigree->alternativeLikelihoodDT[gfreqInd][penIdx],
                   * pPedigree->traitLikelihoodDT[gfreqInd][penIdx],
                   * pPedigree->markerLikelihood);
                   * } */
                  if (alphaV * homoLR + alphaV2 < 0)
                    fprintf (stderr, "HET LR less than 0. Check!!!\n");
                  log10HetLR += log10 (alphaV * homoLR + alphaV2);
                  // if (log10HetLR > 10 || log10HetLR < -40) {
                  /*if(gfreqInd ==0 && j==0){
                   * fprintf(stderr, "gf=%d pen=%d log10HetLR %G, homoLR %G, alt %G, trait %G, mrk %G\n",
                   * gfreqInd, penIdx,log10HetLR,
                   * homoLR, pPedigree->alternativeLikelihoodDT[gfreqInd][penIdx],
                   * pPedigree->traitLikelihoodDT[gfreqInd][penIdx],
                   * pPedigree->markerLikelihood);
                   * //  exit(0);
                   * } */
                }
                if (log10HetLR >= DBL_MAX_10_EXP - 1) {
                  hetLR = DBL_MAX;
                  mp_result[posIdx].het_lr_total = DBL_MAX;
                } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
                  hetLR = 0;
                } else {
                  hetLR = pow (10, log10HetLR);
                  mp_result[posIdx].het_lr_total += hetLR;
                }
                if (mp_result[posIdx].max_penIdx < 0 || hetLR > mp_result[posIdx].max_lr) {
                  mp_result[posIdx].max_lr = hetLR;
                  mp_result[posIdx].max_alpha = alphaV;
                  mp_result[posIdx].max_gfreq = gfreq;
                  mp_result[posIdx].max_penIdx = penIdx;
                }
              } /* end of calculating HET LR */
            }
          }     /* end of genFreq loop */
        }


        /* end of penetrance loop */
        /* save the alternative likelihood */
        for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
          pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
          if ((modelOptions.saveResults == TRUE) && (pPedigree->load_flag == 0)) {      /*save only for the pedigrees which were add for this run */
            pPedigree->load_flag =
              saveAlternative (pPedigree->sPedigreeID,
                               (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome, traitPos,
                               pPedigree->alternativeLikelihoodDT);
          }
          pPedigree->load_flag = 0;
        }

#ifndef SIMPLEPROGRESS
        fprintf (stdout, "%s %d%% complete (~%ld min left)\n",
                 "Combined likelihood evaluations", cL[7] * 100 / eCL[7],
                 (combinedComputeSW->swAccumWallTime * eCL[7] / cL[7] - combinedComputeSW->swAccumWallTime) / 60);
#else
        fprintf (stdout, "%s %d%% complete (~%ld min left)\r",
                 "Calculations", (cL[6] + cL[7]) * 100 / (eCL[6] + eCL[7]),
                 (combinedComputeSW->swAccumWallTime * (eCL[6] + eCL[7]) / (cL[6] + cL[7]) -
                  combinedComputeSW->swAccumWallTime) / 60);
#endif

      } /* end of TP */
      else
        /* multipoint QT or COMBINED */
      {
        for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
          gfreq = modelRange.gfreq[gfreqInd];
          pLocus->pAlleleFrequency[0] = gfreq;
          pLocus->pAlleleFrequency[1] = 1 - gfreq;

          update_locus (&pedigreeSet, traitLocus);
          /* this should be MEAN + SD */
          for (paramIdx = 0; paramIdx < modelRange.nparam; paramIdx++) {
            for (penIdx = 0; penIdx < modelRange.npenet; penIdx++) {
              breakFlag = FALSE;
              for (thresholdIdx = 0; thresholdIdx < modelRange.ntthresh; thresholdIdx++) {
                for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                  mean_DD = modelRange.penet[liabIdx][0][penIdx];
                  mean_Dd = modelRange.penet[liabIdx][1][penIdx];
                  mean_dd = modelRange.penet[liabIdx][2][penIdx];
                  SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
                  SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
                  SD_dd = modelRange.param[liabIdx][2][0][paramIdx];
                  threshold = modelRange.tthresh[liabIdx][thresholdIdx];

                  if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
                    /* check against the hard coded constraint */
                    constraint =
                      (1 - gfreq) * (1 -
                                     gfreq) * mean_dd *
                      SD_dd + 2 * gfreq * (1 - gfreq) * mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
                    /*      fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
                     * constraint, gfreq, mean_DD, SD_DD, 
                     * mean_Dd, SD_DD, 
                     * mean_dd, SD_dd);
                     */
                    if (constraint >= 3.0 || constraint <= -3.0) {
                      breakFlag = TRUE;
                      break;
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

                }       /* liability class Index */
                if (breakFlag == TRUE)
                  continue;
                if (modelOptions.polynomial == TRUE);
                else
                  update_penetrance (&pedigreeSet, traitLocus);
                /* ready for the alternative hypothesis */
                locusList = &savedLocusList;
                xmissionMatrix = altMatrix;
                if (modelOptions.polynomial == TRUE);
                else
                  status =
                    populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2,
                                              initialHetProbAddr, 0, -1, -1, 0);

                /* If we're not on the first iteration, it's not a polynomial build, so
                 * show progress at 1 minute intervals. Have a care to avoid division by zero. */
                sprintf (partialPolynomialFunctionName, "CL_P%%sC%dFM%dof%d",
                         (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
                         mp_result[posIdx].pMarkers[0], modelType.numMarkers);
                if (gfreqInd != 0 || paramIdx != 0 || penIdx != 0) {
                  swStart (combinedComputeSW);
                  compute_likelihood (&pedigreeSet);
                  cL[8]++;
                  swStop (combinedComputeSW);
                  if (statusRequestSignal) {
                    statusRequestSignal = FALSE;
                    if (cL[8] > 1) {    // The first time thru we have no basis for estimation
#ifndef SIMPLEPROGRESS
                      fprintf (stdout, "%s %d%% complete (~%ld min left)\r",
                               "Combined likelihood evaluations", cL[8] * 100 / eCL[8],
                               ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                                eCL[8] / cL[8] - (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
#else
                      fprintf (stdout, "%s %d%% complete (~%ld min left)\r",
                               "Calculations", (cL[6] + cL[8]) * 100 / (eCL[6] + eCL[8]),
                               ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                                (eCL[6] + eCL[8]) / (cL[6] + cL[8]) -
                                (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
#endif
                      fflush (stdout);
                    }
                  }
                } else {  // This _is_ the first iteration
                  swStart (combinedBuildSW);
                  compute_likelihood (&pedigreeSet);
                  cL[8]++;
                  swStop (combinedBuildSW);
#ifndef SIMPLEPROGRESS
                  fprintf (stdout, "%s %d%% complete\r", "Combined likelihood evaluations", cL[8] * 100 / eCL[8]);
#else
                  fprintf (stdout, "%s %d%% complete\r", "Calculations", (cL[6] + cL[8]) * 100 / (eCL[6] + eCL[8]));
#endif
                  fflush (stdout);
                }
                log10_likelihood_alternative = pedigreeSet.log10Likelihood;
                if (isnan (log10_likelihood_alternative))
                  fprintf (stderr, "ALT likelihood is NAN.\n");
                if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
                  log10_likelihood_ratio = 0;
                } else {
                  log10_likelihood_ratio = log10_likelihood_alternative - likelihoodQT[pedigreeSet.numPedigree][gfreqInd]
                    [penIdx][paramIdx][thresholdIdx] - pedigreeSet.log10MarkerLikelihood;
                }
                /* check for overflow problem !!! */
                if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
                  likelihood_ratio = DBL_MAX;
                  mp_result[posIdx].lr_total += DBL_MAX;
                } else
                  /* check for underflow problem too !!! */
                if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
                  likelihood_ratio = 0;
                } else {
                  likelihood_ratio = pow (10.0, log10_likelihood_ratio);
                  mp_result[posIdx].lr_total += likelihood_ratio;
                }
                /* add the result to the right placeholder */
                mp_result[posIdx].lr_count++;

                if (isnan (likelihood_ratio))
                  fprintf (stderr, "LR for the pedigree set is NAN.\n");
                /* caculating the HET */
                for (j = 0; j < modelRange.nalpha; j++) {
                  alphaV = modelRange.alpha[j];
                  alphaV2 = 1 - alphaV;
                  if (alphaV2 < 0)
                    alphaV2 = 0;
                  log10HetLR = 0;
                  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                    homoLR =
                      pPedigree->likelihood / (likelihoodQT[pedIdx][gfreqInd][penIdx][paramIdx][thresholdIdx] *
                                               pPedigree->markerLikelihood);
                    log10HetLR += log10 (alphaV * homoLR + alphaV2);
                  }
                  if (log10HetLR >= DBL_MAX_10_EXP - 1) {
                    hetLR = DBL_MAX;
                    mp_result[posIdx].het_lr_total = DBL_MAX;
                  } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
                    hetLR = 0;
                  } else {
                    adjustedHetLR = hetLR = pow (10, log10HetLR);
                    /* for threshold parameter, we need to make sure the weighting is even */
                    if (1 || modelType.distrib == QT_FUNCTION_CHI_SQUARE) {
                      if (modelRange.ntthresh == 1) {
                        adjustedHetLR *= 2 * (modelType.maxThreshold - modelType.minThreshold);
                      } else if (thresholdIdx == modelRange.ntthresh - 1) {
                        adjustedHetLR *= (2 * modelType.maxThreshold - threshold - modelRange.tthresh[0][thresholdIdx - 1]);
                      } else if (thresholdIdx == 0) {
                        adjustedHetLR *= (threshold + modelRange.tthresh[0][thresholdIdx + 1] - 2 * modelType.minThreshold);
                      } else
                        adjustedHetLR *= modelRange.tthresh[0][thresholdIdx + 1] - modelRange.tthresh[0][thresholdIdx - 1];
                    }
                    mp_result[posIdx].het_lr_total += adjustedHetLR;
                  }
                  if (mp_result[posIdx].max_penIdx < 0 || hetLR > mp_result[posIdx].max_lr) {
                    mp_result[posIdx].max_lr = hetLR;
                    mp_result[posIdx].max_alpha = alphaV;
                    mp_result[posIdx].max_gfreq = gfreq;
                    mp_result[posIdx].max_penIdx = penIdx;
                    mp_result[posIdx].max_paramIdx = paramIdx;
                    mp_result[posIdx].max_thresholdIdx = thresholdIdx;
                  }
                }
              } /* end of threshold loop */
            }   /* end of penetrance loop */
          }     /* end of parameter loop */
        }       /* end of gene freq */

#ifndef SIMPLEPROGRESS
        fprintf (stdout, "%s %d%% complete (~%ld min left)\n",
                 "Combined likelihood evaluations", cL[8] * 100 / eCL[8],
                 (combinedComputeSW->swAccumWallTime * eCL[8] / cL[8] - combinedComputeSW->swAccumWallTime) / 60);
#else
        fprintf (stdout, "%s %d%% complete (~%ld min left)\r",
                 "Calculations", (cL[6] + cL[8]) * 100 / (eCL[6] + eCL[8]),
                 (combinedComputeSW->swAccumWallTime * (eCL[6] + eCL[8]) / (cL[6] + cL[8]) -
                  combinedComputeSW->swAccumWallTime) / 60);
#endif

      } /* end of QT */

      /* print out average and log10(max) and maximizing parameters */
      //      if (modelType.trait == DT || modelType.distrib != QT_FUNCTION_CHI_SQUARE)
      if (modelType.trait == DT)
        avgLR = mp_result[posIdx].het_lr_total / (modelRange.nalpha * mp_result[posIdx].lr_count);
      else
        /* under QT CHISQ, threshold parameter has been evenly weighted */
        avgLR =
          mp_result[posIdx].het_lr_total / (modelRange.nalpha * (mp_result[posIdx].lr_count / modelRange.ntthresh) * 2 *
                                            (modelType.maxThreshold - modelType.minThreshold));

      if (avgLR > 0.214)
        ppl = (avgLR * avgLR) / (-5.77 + 54 * avgLR + avgLR * avgLR);
      else
        ppl = 0;
      max = mp_result[posIdx].max_lr;
      gfreq = mp_result[posIdx].max_gfreq;
      alphaV = mp_result[posIdx].max_alpha;
      penIdx = mp_result[posIdx].max_penIdx;
      paramIdx = mp_result[posIdx].max_paramIdx;
      thresholdIdx = mp_result[posIdx].max_thresholdIdx;
      fprintf (fpHet, "%d %f %.*f %.6e %.6f %f %f",
               (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
               traitPos, ppl >= .025 ? 2 : 3, ppl >= .025 ? rint (ppl * 100.) / 100. : rint (ppl * 1000.) / 1000.,
               avgLR, log10 (max), alphaV, gfreq);
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
        pen_DD = modelRange.penet[liabIdx][0][penIdx];
        pen_Dd = modelRange.penet[liabIdx][1][penIdx];
        pen_dd = modelRange.penet[liabIdx][2][penIdx];
	fprintf (fpHet, " (%.3f,%.3f,%.3f", pen_DD, pen_Dd, pen_dd);
        if (modelType.trait != DT && modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
          SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
          SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
          SD_dd = modelRange.param[liabIdx][2][0][paramIdx];
          fprintf (fpHet, ",%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dd);
        }
        if (modelType.trait != DT) {
          threshold = modelRange.tthresh[liabIdx][thresholdIdx];
          fprintf (fpHet, ",%.3f)", threshold);
        } else
	  fprintf (fpHet, ")");
      }
      /* print out markers used for this position */
      fprintf (fpHet, " (%d", mp_result[posIdx].pMarkers[0]);
      for (k = 1; k < modelType.numMarkers; k++) {
        fprintf (fpHet, ",%d", mp_result[posIdx].pMarkers[k]);
      }
      fprintf (fpHet, ")\n");
      fflush (fpHet);
    }   /* end of walking down the chromosome */
  }     /* end of multipoint */

  /* only for multipoint - deallocate memory  */
  if (modelType.type == MP) {
    /* allocate space to save temporary results */
    if (modelType.trait == DT) {

      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
        free (likelihoodDT[gfreqInd]);
      }
      free (likelihoodDT);
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];

        for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
          free (pPedigree->traitLikelihoodDT[gfreqInd]);
          free (pPedigree->alternativeLikelihoodDT[gfreqInd]);
        }
        free (pPedigree->traitLikelihoodDT);
        free (pPedigree->alternativeLikelihoodDT);
      }
      free (markerNameList);
    }
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

#ifdef SOURCEDIGRAPH
  if (modelOptions.polynomial == TRUE)
    dumpSourceParenting ();
#endif

  /* Final dump and clean-up for performance. */
  swStop (overallSW);
  swDump (overallSW);
#ifdef POLYSTATISTICS
  if (modelOptions.polynomial == TRUE)
    polyStatistics ("End of run");
#endif
#ifdef DMUSE
  fprintf (stderr, "Missed/Used %d/%d 24s, %d/%d 48s, %d/%d 100s\n",
           missed24s, used24s, missed48s, used48s, missed100s, used100s);
#endif
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
    fclose (fpTP);
  }
  fclose (fpHet);
  return 0;
}
