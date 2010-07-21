#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> // For index on some platforms
#include <errno.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "config.h"
#include "../utils/utils.h"
#include "../utils/pageManagement.h"
#include "../pedlib/pedlib.h"

#ifdef STUDYDB
  #include "../database/StudyDB.h"
#endif

#define BUFFSIZE 256

#define fault(...) \
do { \
  int length; \
  char message[MAXLOGMSG + 1], *pMessage = message; \
  pMessage += length = snprintf (message, MAXLOGMSG, "ERROR, "); \
  snprintf (pMessage, MAXLOGMSG - length,  __VA_ARGS__); \
  swLogMsg (stderr, message); \
  fault++; \
} while(0)

/* Structure for the configuration parser dispatch table.
 * 'key' is the configuration directive
 * 'parse' is an int-returning function that takes (char**, int, void*)
 * 'hint' is a pointer that will be passed as the last argument to 'parse'
 */
typedef struct {
  char *key;
  int (*parse) (char **, int, void *);
  void *hint;
} st_dispatch;

/* There are some directives that we use not only as keys in the
 * dispatch table, but we also compare against inside the parsing
 * routines, either to facilitate consolidating largely-identical
 * routines, or because they are also valid arguments to the 
 * 'Constrain' directive. To prevent the literal strings used in 
 * the comparisons being changed in one place, but not in other
 * places, we defined symbols to respresnt those strings here.
 */
#define QT_STR            "QT"
#define QTT_STR           "QTT"
#define THETA_STR         "Theta"
#define PENETRANCE_STR    "Penetrance"
#define MEAN_STR          "Mean"
#define STANDARDDEV_STR   "StandardDev"
#define DEGOFFREEDOM_STR  "DegreesOfFreedom"
#define THRESHOLD_STR     "Threshold"

/* Structure for returning results from expand_vals() in terms of the
 * contents of each comma separated substring, instead of a list of doubles.
 */
typedef struct {
  int type;
  union {
    double val;
    struct {
      double start, end, incr;
    } range;
    int symbol;
  } vun;
} st_valuelist;

/* Legal values for the 'type' field of st_valuelist. */
#define VL_VALUE          0
#define VL_RANGE          1
#define VL_RANGE_SYMBEND  2
#define VL_SYMBOL         4

/* Legal values for the vun.symbol field of st_valuelist, if type is VL_SYMBOL */
#define VL_SYM_MARKER      0

typedef struct {
  int traitPositions,
    constraints,
    maxclass,
    penetrance,
    mean,
    standardDev,
    degOfFreedom;
} st_observed;

/* Some trait parameters (Mean, StandardDev, DegreesOfFreedom) can be specified
 * with or without a trait genotype. We need to keep track of which, if any, were
 * specified with trait genotypes.
 */

#define PARM_BARE 1
#define PARM_GENO 2

/** @defgroup vettedGlobals Vetted globals that will stay globals
    @{
  Globals that have been reviewed and found to be appropriately
  global. They're gathered here in order to consolidate documentation.
  Eventually, when we've eliminated all unnecessary globals, the "Globals"
  tab in doxygen will be useful again, and this group can go away. Right
  now it's just too cluttered. */

ModelOptions *modelOptions; ///< All configuration options for the working model
ModelRange *modelRange; ///< All expanded parameter ranges for the working model
ModelType *modelType; ///< All typing information for the working model
#ifdef STUDYDB
struct StudyDB studyDB; ///< All parameters relating to database access for storing studies
#endif
/*@}*/

/* Module globals */

static ModelOptions staticModelOptions; // Statically-allocated version for John's dispatch table
static ModelRange staticModelRange; // Statically-allocated version for John's dispatch table
static ModelType staticModelType; // Statically-allocated version for John's dispatch table
static char buff[BUFFSIZE] = "";   /* These two are global to provide context to getNextTokgroup */
static char *buffptr = NULL;
static st_observed observed;       /* track non-obvious directives */
static char *conffilename=NULL;
static int lineno = 0;             /* so failure messages are more useful to users */
static int fault = 0;              /* For tracking number of configuration faults */

/* prototypes for non-public routines */
void initializeDefaults ();
int expandVals (char **toks, int numtoks, double **vals_h, st_valuelist **vlist_h);
int lookupDispatch (char *key, st_dispatch *table);
int compareDispatch (const void *a, const void *b);
int getNextTokgroup (FILE *fp, char ***tokgroup_h, int *tokgroupsize);
int tokenizeLine (char *line, char ***tokgroup_h, int *tokgroupsize);
int singleDigit (char *str);
void dumpmodelOptions (ModelOptions *mo);
void bail (char *fmt, char *arg);

/* functions for use in the dispatch table */
int set_optionfile (char **toks, int numtoks, void *filename);
int set_flag (char **toks, int numtoks, void *flag);
int clear_flag (char **toks, int numtoks, void *flag);
int set_int (char **toks, int numtoks, void *field);

int set_traitPositions (char **toks, int numtoks, void *unused);
int set_alleleFreq (char **toks, int numtoks, void *unused);
int set_geneFreq (char **toks, int numtoks, void *unused);
int set_dprime (char **toks, int numtoks, void *unused);
int set_theta (char **toks, int numtoks, void *unused);
int set_alpha (char **toks, int numtoks, void *unused);
int set_penetrance (char **toks, int numtoks, void *unused);
int set_constraint (char **toks, int numtoks, void *unused);
int set_multipoint (char **toks, int numtoks, void *unused);
int set_markerAnalysis (char **toks, int numtoks, void *unused);
int set_mapFlag (char **toks, int numtoks, void *unused);
int set_disequilibrium (char **toks, int numtoks, void *unused);
int set_quantitative (char **toks, int numtoks, void *unused);
int set_qt_mean (char **toks, int numtoks, void *unused);
int set_qt_degfreedom (char **toks, int numtoks, void *unused);
int set_qt_standarddev (char **toks, int numtoks, void *unused);
int set_qt_threshold (char **toks, int numtoks, void *unused);
int set_qt_truncation (char **toks, int numtoks, void *unused);
int set_affectionStatus (char **toks, int numtoks, void *unused);
int set_study_parameters (char **toks, int numtoks, void *unused);
int set_resultsprefix (char **toks, int numtoks, void *unused);

int skip_analysis (char **toks, int numtoks, void *unused);
int noop (char **toks, int numtoks, void *unused);


st_dispatch dispatchTable[] = { {"FrequencyFile", set_optionfile, &staticModelOptions.markerfile},
				{"MapFile", set_optionfile, &staticModelOptions.mapfile},
				{"PedigreeFile", set_optionfile, &staticModelOptions.pedfile},
				{"LocusFile", set_optionfile, &staticModelOptions.datafile},
				{"BayesRatioFile", set_optionfile, &staticModelOptions.avghetfile},
				{"PPLFile", set_optionfile, &staticModelOptions.pplfile},
				{"CountFile", set_optionfile, &staticModelOptions.ccfile},
				{"MODFile", set_optionfile, &staticModelOptions.modfile},
				{"SurfaceFile", set_optionfile, &staticModelOptions.intermediatefile},
				{"NIDetailFile", set_optionfile, &staticModelOptions.dkelvinoutfile},

				{"NonPolynomial", clear_flag, &staticModelOptions.polynomial},
				{"Imprinting", set_flag, &staticModelOptions.imprintingFlag},
				{"SexLinked", set_flag, &staticModelOptions.sexLinked},
#ifndef DISTRIBUTION
				{"FixedModels", clear_flag, &staticModelOptions.integration},
#endif
				{"DryRun", set_flag, &staticModelOptions.dryRun},
				{"ExtraMODs", set_flag, &staticModelOptions.extraMODs},
				{"ForceBRFile", set_flag, &staticModelOptions.forceAvghetFile},

				{"PolynomialScale", set_int, &staticModelOptions.polynomialScale},
				{"LiabilityClasses", set_int, &staticModelRange.nlclass},
				{"DiseaseAlleles", set_int, &staticModelRange.nalleles},
				{"MaxIterations", set_int, &staticModelOptions.maxIterations},

				{"TraitPositions", set_traitPositions, NULL},
				{"MarkerAlleleFrequency", set_alleleFreq, NULL},
				{"DiseaseGeneFrequency", set_geneFreq, NULL},
				{"DPrime", set_dprime, NULL},
				{THETA_STR, set_theta, NULL},
				{"Alpha", set_alpha, NULL},
				{PENETRANCE_STR, set_penetrance, NULL},
				{"Constraint", set_constraint, NULL},
				{"Multipoint", set_multipoint, NULL},
				{"MarkerToMarker", set_markerAnalysis, NULL},
				{"SexSpecific", set_mapFlag, NULL},
				{"LD", set_disequilibrium, NULL},
				{QT_STR, set_quantitative, NULL},
				{QTT_STR, set_quantitative, NULL},
				{MEAN_STR, set_qt_mean, NULL},
				{DEGOFFREEDOM_STR, set_qt_degfreedom, NULL},
				{STANDARDDEV_STR, set_qt_standarddev, NULL},
				{THRESHOLD_STR, set_qt_threshold, NULL},
				{"Truncate", set_qt_truncation, NULL},
				{"PhenoCodes", set_affectionStatus, NULL},
				{"SurfacesPath", set_resultsprefix, NULL},
				// {"condfile", set_condrun, &staticModelOptions.condFile},
				{"ProgressDelaySeconds", set_int, &swProgressDelaySeconds},
				{"ProgressLevel", set_int, &swProgressLevel},
				{"Study", set_study_parameters, NULL},

				{"SkipEstimation", noop, NULL},
				{"SkipPedCount", noop, NULL},
				{"SkipAnalysis", skip_analysis, NULL}
};


void initializeDefaults ()
{
  /* Initialize the the global configuration structures to default values */
  memset (&staticModelOptions, 0, sizeof (ModelOptions));
  memset (&staticModelRange, 0, sizeof (ModelRange));
  memset (&staticModelType, 0, sizeof (ModelType));

  strcpy (staticModelOptions.markerfile, DEFAULTMARKERFILENAME);
  strcpy (staticModelOptions.mapfile, DEFAULTMAPFILENAME);
  strcpy (staticModelOptions.pedfile, DEFAULTPEDFILENAME);
  strcpy (staticModelOptions.datafile, DEFAULTDATAFILENAME);
  strcpy (staticModelOptions.avghetfile, DEFAULTAVGHETFILENAME);
  strcpy (staticModelOptions.resultsprefix, DEFAULTRESULTSPREFIX);

  MALCHOKE(staticModelOptions.sUnknownPersonID, sizeof (char) * 2,void *);
  strcpy (staticModelOptions.sUnknownPersonID, "0");

  staticModelOptions.equilibrium = LINKAGE_EQUILIBRIUM;
  staticModelOptions.markerAnalysis = FALSE;
  staticModelOptions.saveResults = FALSE;
  staticModelOptions.polynomial = TRUE;
  staticModelOptions.integration = TRUE;
  staticModelOptions.maxIterations = -1;
  staticModelOptions.imprintingFlag = FALSE;
  staticModelOptions.mapFlag = SA;
  staticModelOptions.sexLinked = FALSE;
  staticModelOptions.dryRun = FALSE;
  staticModelOptions.forceAvghetFile = FALSE;
  staticModelOptions.polynomialScale = 0;
  staticModelOptions.extraMODs = FALSE;
  staticModelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = -DBL_MAX;
  staticModelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] = -DBL_MAX;
  staticModelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED] = -DBL_MAX;

  /* Set default values for PPL calculations. Maybe should be configuration parameters, esp. if
   we want to call them defaults, which implies the possibility of non-default settings. */
  staticModelOptions.thetaCutoff[0] = 0.05;   /* LRs when theta < cutoff are weighted heavier */
  staticModelOptions.thetaCutoff[1] = 0.05;
  staticModelOptions.thetaWeight = 0.95;      /* Weight for LRs when theta < cutoff */
  staticModelOptions.prior = 0.02;            /* Prior probability of linkage */
  staticModelOptions.LDprior = 0.02;          /* Prior probability of LD given close linkage */
  
  staticModelRange.nalleles = 2;
  staticModelRange.nlclass = 1;
  staticModelRange.npardim = 0;
  staticModelRange.nlambdas = 0;
  staticModelRange.maxnlambdas = 0;
  staticModelRange.tlocRangeStart = -1;
  staticModelRange.tlocRangeIncr = -1;
  staticModelRange.tlmark = FALSE;
  staticModelRange.atypicalQtTrait = FALSE;
  staticModelRange.microsats = FALSE;

  staticModelType.type = TP;
  staticModelType.trait = DT;
  staticModelType.distrib = -1;
  /* set default for QT */
  staticModelType.mean = -DBL_MAX;
  staticModelType.sd = -DBL_MAX;
  staticModelType.minOriginal = -999999999.00;
  staticModelType.maxOriginal = 999999999.00;

  /* Sort the dispatch table so the binary search works */
  qsort (dispatchTable, sizeof (dispatchTable) / sizeof (st_dispatch), sizeof (st_dispatch),
	 compareDispatch);
  memset (&observed, 0, sizeof (st_observed));
  return;
}


void readConfigFile (char *config)
{
  int numtoks, tokgroupsize=0, va;
  char **toks = NULL;
  FILE *conffp;

  conffilename = config;
  if ((conffp = fopen (config, "r")) == NULL)
    ERROR ("Unable to open '%s', %s\n", config, strerror (errno));

  while ((numtoks = getNextTokgroup (conffp, &toks, &tokgroupsize)) > 0) {
    for (va = 0; va < numtoks; va++)
      DIAG (INPUTFILE, 1, {fprintf (stderr, "token %d: %s\n", va, toks[va]);});
    if ((va = lookupDispatch (toks[0], dispatchTable)) >= 0) {
      DIAG (INPUTFILE, 1, {fprintf (stderr, "directive '%s' matches at index %d\n", toks[0], va);});
      (*dispatchTable[va].parse) (toks, numtoks, dispatchTable[va].hint);
    } else
      ERROR ("Directive '%s' on line %d is %s", toks[0],
	     lineno, (va == -1) ? "unknown" : "not unique");
  }

  if (tokgroupsize > 0)
    free (toks);
  return;
}


/* Here we have to marshall command line arguments into a format that looks like
 * we read them out of a file. Directives on the command line must be prefixed
 * with '--'; for each directive, we concatenate the directive and all the arguments
 * up to the next directive into a single buffer, and run that buffer through 
 * permuteLine() and tokenizeLine(), and get back the same kind of list of tokens
 * that getNextTokGroup() returns.
 */
void parseCommandLine (int argc, char *argv[])
{
  int curidx, bufflen=0, numtoks, tokgroupsize=0, va;
  char **toks=NULL;
  
  if (argc == 0)
    return;

  if (strncmp (argv[0], "--", 2) != 0)
    ERROR ("Expected directive on command line, found '%s'", argv[0]);
  bufflen = strlen (argv[0]) - 2;
  strcpy (buff, argv[0]+2);
  
  curidx = 1;
  while (curidx < argc) {
    if (strncmp (argv[curidx], "--", 2) == 0) {
      DIAG (INPUTFILE, 1, { fprintf (stderr, "Buffer is '%s'\n", buff);});
      permuteLine (buff, BUFFSIZE);
      numtoks = tokenizeLine (buff, &toks, &tokgroupsize);
      if ((va = lookupDispatch (toks[0], dispatchTable)) >= 0) {
	DIAG (INPUTFILE, 1, { fprintf (stderr, "Directive '%s' matches at index %d\n", toks[0], va);});
	(*dispatchTable[va].parse) (toks, numtoks, dispatchTable[va].hint);
      } else
	ERROR ("Directive '%s' on command line is %s", toks[0],
		(va == -1) ? "unknown" : "not unique");
      bufflen = strlen (argv[curidx]) - 2;
      strcpy (buff, argv[curidx]+2);

    } else {
      buff[bufflen++] = ' ';
      buff[bufflen] = '\0';
      bufflen += strlen (argv[curidx]);
      strcat (buff, argv[curidx]);
    }
    curidx++;
  }
    
  DIAG (INPUTFILE, 1, { fprintf (stderr, "Buffer is '%s'\n", buff);});
  permuteLine (buff, BUFFSIZE);
  numtoks = tokenizeLine (buff, &toks, &tokgroupsize);
  if ((va = lookupDispatch (toks[0], dispatchTable)) >= 0) {
    DIAG (INPUTFILE, 1, { fprintf (stderr, "Directive '%s' matches at index %d\n", toks[0], va);});
    (*dispatchTable[va].parse) (toks, numtoks, dispatchTable[va].hint);
  } else
    ERROR ("Directive '%s' on command line is %s\n", toks[0],
	    (va == -1) ? "unknown" : "not unique");
  
  if (tokgroupsize > 0)
    free (toks);
  return;
}


void validateConfig ()
{
  /* Check that the configuration directives are both compatible and sufficient. We may
   * check the same combinations from multiple directions. For the most part, we don't
   * imply anything anymore (ex. implying a sex specific map by configuring male theta
   * values). The user must specify what they want, we won't guess. One major exception
   * is MarkerToMarker, which silently turns on FixedModels, if it's not on already.
   */

  if (staticModelOptions.polynomialScale && ! staticModelOptions.polynomial)
    fault ("PolynomialScale is incompatible with NonPolynomial\n");
    
  if (staticModelOptions.markerAnalysis) {
    /* MarkerToMarker is a special case. It only supports TP, LD, fixed grid thetas
     * and D-primes. Since only markers are considered, we disallow any directives
     * related to the trait. LD implies no sex-specific, and TP means no Multipoint,
     * so those directives are out, too. Once we're done, we return immediately, so
     * we don't have to worry about MarkerToMarker vs. trait-to-marker later on.
     */
    if (staticModelOptions.imprintingFlag)
      fault ("Trait directives (Imprinting) are incompatible with MarkerToMarker\n");
    if (staticModelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] != -DBL_MAX)
      fault ("Trait directives (PhenoCodes) are incompatible with MarkerToMarker\n");
    if (staticModelRange.nalleles != 2)    
      fault ("Trait directives (DiseaseAlleles) are incompatible with MarkerToMarker\n");
    if (staticModelRange.nlclass != 1)    
      fault ("Trait directives (LiabilityClasses) are incompatible with MarkerToMarker\n");
    if (staticModelType.trait == QT) 
      fault ("Trait directives (%s) are incompatible with MarkerToMarker\n", QT_STR);
    if (staticModelType.trait == CT) 
      fault ("Trait directives (%s) are incompatible with MarkerToMarker\n", QTT_STR);
    if (staticModelRange.ntthresh > 0)
      fault ("Trait directives (%s) are incompatible with MarkerToMarker\n", THRESHOLD_STR);
    if (staticModelType.minOriginal != -999999999.00 || staticModelType.maxOriginal != 999999999.00)
      fault ("Trait directives (Truncate) are incompatible with MarkerToMarker\n");
    if (staticModelType.type == MP)
      fault ("Multipoint is incompatible with MarkerToMarker\n");
    if (observed.traitPositions)
      fault ("Multipoint directives (TraitPositions) are incompatible with MarkerToMarker\n");
    if (staticModelRange.nafreq > 0)
      fault ("MarkerAlleleFrequency is incompatible with MarkerToMarker\n");
    if (staticModelRange.ngfreq > 0)
      fault ("Trait directives (DiseaseGeneFrequency) are incompatible with MarkerToMarker\n");
    if (staticModelRange.nalpha > 0)
      fault ("Trait directives (Alpha) are incompatible with MarkerToMarker\n");
    if (observed.penetrance)
      fault ("Trait directives (%s) are incompatible with MarkerToMarker\n", PENETRANCE_STR);
    if (observed.mean)
      fault ("Trait directives (%s) are incompatible with MarkerToMarker\n", MEAN_STR);
    if (observed.standardDev)
      fault ("Trait directives (%s) are incompatible with MarkerToMarker\n", STANDARDDEV_STR);
    if (observed.degOfFreedom)
      fault ("Trait directives (%s) are incompatible with MarkerToMarker\n", DEGOFFREEDOM_STR);
    if (observed.constraints)
      fault ("Trait directives (Constraint) are incompatible with MarkerToMarker\n");
    if (staticModelOptions.avghetfile[0] != '\0' && ! staticModelOptions.forceAvghetFile)
      INFO ("MarkerToMarker will write no output to BayesRatioFile");
    if (staticModelOptions.dkelvinoutfile[0] != '\0')
      INFO ("MarkerToMarker will write no output to NIDetailFile");

    if (! staticModelOptions.integration) {
      if (staticModelRange.ndprime == 0 && staticModelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM)
	fault ("FixedModels and LD require DPrime\n");
      if (staticModelRange.ndprime > 0 && staticModelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
	fault ("FixedModels and DPrime requires LD\n");
      if (staticModelRange.thetacnt == NULL)
	fault ("MarkerToMarker and FixedModels require %s\n", THETA_STR);
    } else {
      if (staticModelRange.ndprime > 0)
	fault ("MarkerToMarker and DPrime require FixedModels\n");
      if (staticModelRange.thetacnt != NULL)
	fault ("MarkerToMarker and %s require FixedModels\n", THETA_STR);
    }
    if (fault)
      ERROR ("%d configuration errors detected", fault);
    return;
  } 
  /* Everything hereafter is trait-to-marker */
  /* First, try to rule out the simplest invalid combinations of options */

  /* We only handle bi-allelic traits for now */
  if (staticModelRange.nalleles != 2)
    fault ("DiseaseAlleles must be set to 2; polyallelic traits are not supported\n");
  
  /* set_affectionStatus() guarantees that 0, 1 or all of these will be set */
  if (staticModelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] != -DBL_MAX) {
    if (staticModelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] != -DBL_MAX) {
      if (staticModelType.trait == QT)
	fault ("PhenoCodes with 3 arguments is incompatible with QT\n");
    } else {
      if (staticModelType.trait != QT)
	fault ("PhenoCodes with 1 argument requires QT\n");
    }
  }
  
  if (staticModelType.type == MP) {
    /* Multipoint */
    if (staticModelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM)
      fault ("LD is incompatible with Multipoint\n");
    if (staticModelOptions.extraMODs)
      fault ("ExtraMODs is incompatible with Multipoint\n");
    if (staticModelRange.nafreq > 0)
      fault ("MarkerAlleleFrquency is incompatible with Multipoint\n");
    if (staticModelOptions.pplfile[0] != '\0')
      WARNING ("Multipoint analysis will write no output to PPLFile '%s', but will write PPLs in BayesRatioFile '%s'", staticModelOptions.pplfile, staticModelOptions.avghetfile);
    if (! observed.traitPositions) 
      fault ("Multipoint requires TraitPositions\n");
  } else {
    /* Two point */
    if (observed.traitPositions) 
      fault ("TraitPositions requires Multipoint\n");
  }
  
  if (staticModelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM && staticModelOptions.mapFlag == SS) 
    fault ("SexSpecific is not supported with LD\n");

  if (observed.mean && (staticModelType.trait == DT || staticModelType.distrib != QT_FUNCTION_T))
    fault ("%s requires %s Normal or %s Normal\n", MEAN_STR, QT_STR, QTT_STR);
  if (observed.standardDev && (staticModelType.trait == DT || staticModelType.distrib != QT_FUNCTION_T))
    fault ("%s requires %s Normal or %s Normal\n", STANDARDDEV_STR, QT_STR, QTT_STR);
#ifdef DISTRIBUTION
  if (observed.degOfFreedom)
    fault ("ChiSq distribution is still under development and unavailable at this time\n");
#endif
  if (observed.degOfFreedom && (staticModelType.trait == DT || staticModelType.distrib != QT_FUNCTION_CHI_SQUARE))
    fault ("%s requires %s ChiSq or %s ChiSq\n", STANDARDDEV_STR, QT_STR, QTT_STR);
  
  if (staticModelOptions.integration) {
    /* Dynamic sampling, so disallow all fixed model directives */
    if (observed.penetrance)
      fault ("%s requires FixedModels\n", PENETRANCE_STR);
    if (observed.mean & PARM_GENO)
      fault ("%s with trait genotypes requires FixedModels\n", MEAN_STR);
    if (observed.degOfFreedom & PARM_GENO)
      fault ("%s with trait genotypes requires FixedModels\n", DEGOFFREEDOM_STR);
    if (observed.constraints)
      fault ("Constraint requires FixedModels\n");
    if (staticModelRange.thetacnt != NULL)
      fault ("%s requires FixedModels\n", THETA_STR);
    if (staticModelRange.ndprime > 0)
      fault ("DPrime requires FixedModels\n");
    if (staticModelRange.ngfreq > 0)
      fault ("DiseaseGeneFrequency requires FixedModels\n");
    if (staticModelRange.nafreq > 0)
      fault ("MarkerAlleleFrequency requires FixedModels\n");
    if (staticModelRange.nalpha > 0)
      fault ("Alpha requires FixedModels\n");
    
    if (staticModelType.trait == CT && staticModelRange.ntthresh > 0) {
      if (staticModelRange.ntthresh != 2)
	fault ("QTT allows exactly two Threshold values (min and max)\n");
    }

    if (fault)
      ERROR ("%d configuration errors detected", fault);
    return;
  }
  
  /* So much for the low-hanging fruit... */
  
  if (staticModelOptions.dkelvinoutfile[0] != '\0')
    WARNING ("Analysis with FixedModels will write no output to NIDetailFile '%s'", staticModelOptions.dkelvinoutfile);

  if (staticModelOptions.mapFlag == SS) {
    if (staticModelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM)
      fault ("SexSpecific is not supported with LD\n");
  }

  if (staticModelType.type == MP) {
    /* Multipoint */
    if (staticModelRange.thetacnt != NULL)
      fault ("%s is incompatible with Multipoint\n", THETA_STR);
    if (staticModelRange.ndprime > 0)
      fault ("DPrime is incompatible with Multipoint\n");
  } else {
    /* Two point */
    if (staticModelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM && staticModelRange.ndprime == 0)
      fault ("FixedModels with LD requires DPrime\n");
    if (staticModelOptions.equilibrium == LINKAGE_EQUILIBRIUM && staticModelRange.ndprime > 0)
      fault ("FixedModels with DPrime requires LD\n");
    if (staticModelRange.thetacnt == NULL)
      fault ("FixedModels without Multipoint requires %s\n", THETA_STR);
  }

  if (staticModelRange.ngfreq == 0)
    fault ("FixedModels requires DiseaseGeneFrequency\n");
  if (! staticModelRange.alpha)
    fault ("FixedModels requires Alpha\n");

  if (staticModelType.trait == DT) {
    if (! observed.penetrance)
      fault ("Dichotomous trait requires %s\n", PENETRANCE_STR);
    if (observed.mean)
      fault ("%s requires %s Normal or %s Normal\n", MEAN_STR, QT_STR, QTT_STR);
    if (observed.standardDev)
      fault ("%s requires %s Normal or %s Normal\n", STANDARDDEV_STR, QT_STR, QTT_STR);
    if (observed.degOfFreedom)
      fault ("%s requires %s ChiSq or %s ChiSq\n", DEGOFFREEDOM_STR, QT_STR, QTT_STR);
    if (staticModelRange.ntthresh > 0)
      fault ("%s requires %s\n", THRESHOLD_STR, QTT_STR);
  } else {
    if (observed.penetrance)
      fault ("%s is incompatible with %s\n", PENETRANCE_STR, staticModelType.trait == QT ? QT_STR : QTT_STR);
    if (staticModelType.distrib == QT_FUNCTION_T) {
      if (! observed.mean) {
	fault ("%s Normal requires %s\n", staticModelType.trait == QT ? QT_STR : QTT_STR, MEAN_STR);
      } else if (! observed.mean & PARM_GENO) {
      	fault ("FixedModels and %s Normal requires %s with trait genotypes\n", staticModelType.trait == QT ? QT_STR : QTT_STR, MEAN_STR);
      }
      if (! observed.standardDev) 
	fault ("%s Normal requires %s\n", staticModelType.trait == QT ? QT_STR : QTT_STR, STANDARDDEV_STR);
    }
    if (staticModelType.distrib == QT_FUNCTION_CHI_SQUARE) {
      if (! observed.degOfFreedom) {
	fault ("%s ChiSq requires %s\n", staticModelType.trait == QT ? QT_STR : QTT_STR, DEGOFFREEDOM_STR);
      } else if (! observed.degOfFreedom & PARM_GENO) {
      	fault ("FixedModels and %s ChiSq requires %s with trait genotypes\n", staticModelType.trait == QT ? QT_STR : QTT_STR, DEGOFFREEDOM_STR);
      }
    }
    
    if (staticModelType.trait == CT) {
      if (staticModelRange.ntthresh == 0)
	fault ("%s requires %s\n", QTT_STR, THRESHOLD_STR);
    } else {
      if (staticModelRange.ntthresh > 0)
	fault ("%s requires %s\n", THRESHOLD_STR, QTT_STR);
    }
  }
  /* Look for faults before we call checkImprintingPenets(), which doesn't do well
   * if the user hasn't specified any penetrance values at all.
   */
  if (fault)
    ERROR ("%d configuration errors detected", fault);
  
  if ((checkImprintingPenets (&staticModelRange, staticModelOptions.imprintingFlag) < 0)) {
    if (staticModelOptions.imprintingFlag) {
      if (staticModelType.trait == DT) {
	fault ("Imprinting requires Penetrance values for the dD trait genotype\n");
      } else if (staticModelType.distrib == QT_FUNCTION_T) {
	fault ("Imprinting requires Mean values for the dD trait genotype\n"); 
      } else { /* QT_FUNCTION_CHI_SQUARE */	
	fault ("Imprinting requires DegreesOfFreedom values for the dD trait genotype\n");
      }
    } else 
      if (staticModelType.trait == DT) {
	fault ("Penetrance values for the dD trait genotype require Imprinting\n");
      } else if (staticModelType.distrib == QT_FUNCTION_T) {
	fault ("Mean values for the dD trait genotype require Imprinting\n");
      } else { /* QT_FUNCTION_CHI_SQUARE */
	fault ("DegreesOfFreedom values for the dD trait genotype require Imprinting\n");
      }
  }
  if (observed.maxclass != 0 && staticModelRange.nlclass < observed.maxclass)
    fault ("A Constraint references a liability class %d that is not specified with LiabilityClass\n", observed.maxclass);
  
  if (fault)
    ERROR ("%d configuration errors detected", fault);

  return;
}


void fillConfigDefaults (ModelRange *modelRange, ModelOptions *modelOptions, ModelType *modelType)
{
  int i;
  /*  double integrationLDDPrimeValues[33] =
    {-0.9982431840532, -0.9956010478552, -0.9790222658168, -0.9590960631620, -0.8761473165029,
     -0.8727421201131, -0.7013933644534, -0.6582769255267, -0.6492284325645, -0.5666666666667,
     -0.5000000000000, -0.3808991135940, -0.3582614645881, -0.2517129343453, -0.2077777777778,
     -0.1594544658298, 0.0000000000000, 0.1594544658298, 0.2077777777778, 0.2517129343453,
     0.3582614645881, 0.3808991135940, 0.5000000000000, 0.5666666666667, 0.6492284325645,
     0.6582769255267, 0.7013933644534, 0.8727421201131, 0.8761473165029, 0.9590960631620,
     0.9790222658168, 0.9956010478552, 0.9982431840532}; old */
  double integrationLDDPrimeValues[67] = 
    {-0.9991215920266,-0.9978005239276,-0.9895111329084,-0.9795480315810,-0.9380736582515,
     -0.9363710600566,-0.8506966822267,-0.8291384627634,-0.8246142162823,-0.7833333333333,
     -0.7500000000000,-0.6904495567970,-0.6791307322941,-0.6258564671727,-0.6038888888889,
     -0.5797272329149,-0.5000000000000,-0.4202727670851,-0.3961111111111,-0.3741435328273,
     -0.3208692677059,-0.3095504432030,-0.2500000000000,-0.2166666666667,-0.1753857837177,
     -0.1708615372366,-0.1493033177733,-0.0636289399434,-0.0619263417485,-0.0204519684190,
     -0.0104888670916,-0.0021994760724,-0.0008784079734,0.0000000000, 0.0008784079734,
     0.0021994760724,0.0104888670916,0.0204519684190,0.0619263417485,0.0636289399434,
     0.1493033177733,0.1708615372366,0.1753857837177,0.2166666666667,0.2500000000000,
     0.3095504432030,0.3208692677059,0.3741435328273,0.3961111111111,0.4202727670851,
     0.5000000000000,0.5797272329149,0.6038888888889,0.6258564671727,0.6791307322941,
     0.6904495567970,0.7500000000000,0.7833333333333,0.8246142162823,0.8291384627634,
     0.8506966822267,0.9363710600566,0.9380736582515,0.9795480315810,0.9895111329084,
     0.9978005239276,0.9991215920266};







  /* Fill in default values for fields that could have been configured, if they weren't */

  if (staticModelType.trait == DT) {
    if (staticModelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] == -DBL_MAX) {
      staticModelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = AFFECTION_STATUS_UNKNOWN;
      staticModelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] = AFFECTION_STATUS_UNAFFECTED;
      staticModelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED] = AFFECTION_STATUS_AFFECTED;
    }
  } else {
    /* QT or CT */
    if (staticModelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] == -DBL_MAX)
      staticModelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = -99.99;
    if (staticModelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] == -DBL_MAX) {
      staticModelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] = -88.88;
      staticModelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED] = 88.88;
    }
  }

  if (staticModelOptions.polynomial && ! staticModelOptions.polynomialScale)
    staticModelOptions.polynomialScale = 1;
  if ((staticModelType.type == TP) && (staticModelOptions.pplfile[0] == '\0'))
    strcpy (staticModelOptions.pplfile, DEFAULTPPLFILENAME);

  /* MarkerToMaker: validateConfig should have already weeded out patently
   * incompatible options. Here, force LD, FixedModels and fill in default
   * Theta and DPrime values, if needed.
   */
  if (staticModelOptions.markerAnalysis) {
    staticModelOptions.integration = FALSE;
    if (staticModelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM && staticModelRange.ndprime == 0)
      /* Default range of DPrimes is -1 to 1 in steps of 0.02 */
      for (i = -50; i <= 50; i++)
	addDPrime (&staticModelRange, 0.02 * i);
    if (staticModelRange.thetacnt == NULL || staticModelRange.thetacnt[SEXAV] == 0)
      /* Default range of Thetas if 0 to 0.5 in steps of 0.01 */
      for (i = 0; i < 50; i++)
	addTheta (&staticModelRange, THETA_AVG, 0.01 * i);
  }

  /* Fix up the DPrimes if LD is turned on*/
  if (staticModelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM) {
    if (staticModelOptions.integration == TRUE) {
      /* If integration (that is, dynamic grid) is turned on, then there shouldn't
       * be any DPrimes at all, and we need to insert some. These are MAGIC DPrimes,
       * statically declared at the top of this function.
       */
      for (i = 0; i < 67; i++)
	addDPrime (&staticModelRange, integrationLDDPrimeValues[i]);
      
    } else {
      /* If not integration (that is, fixed models) make sure the user
       * didn't omit 0 from the range of DPrimes. Silently add one if needed.
       */
      for (i=0; i<staticModelRange.ndprime; i++)
	if (fabs(staticModelRange.dprime[i]) <= ERROR_MARGIN) break;
      if (i == staticModelRange.ndprime)
	addDPrime (&staticModelRange, (double) 0.0);
    }
  }
  if (staticModelOptions.integration == TRUE) {
    if (staticModelType.distrib == QT_FUNCTION_T) {
      if (! observed.mean) {
	addPenetrance (&staticModelRange, 0, -3.0);
	addPenetrance (&staticModelRange, 0, 3.0);
      }
      if (! observed.standardDev) {
	addParameter (&staticModelRange, 0, 0.7);
	addParameter (&staticModelRange, 0, 1.0);
      }
      if (staticModelType.trait == CT && staticModelRange.ntthresh == 0) {
	addTraitThreshold (&staticModelRange, 0);
	addTraitThreshold (&staticModelRange, 3);
      }
    }
    if (staticModelType.distrib == QT_FUNCTION_CHI_SQUARE) {
      if (! observed.degOfFreedom) {
	addPenetrance (&staticModelRange, 0, 0.05);
	addPenetrance (&staticModelRange, 0, 30);
      }
      if (staticModelType.trait == CT && staticModelRange.ntthresh == 0) {
	addTraitThreshold (&staticModelRange, 0.05);
	addTraitThreshold (&staticModelRange, 30);
      }
    }
    if (staticModelType.trait != DT) 
      duplicatePenets (&staticModelRange, staticModelOptions.imprintingFlag);
    
  } else { /* no integration (fixed grid) */
    /* For 2-point and fixed models, make sure there's a Theta of 0.5 */
    if (staticModelType.type == TP) {
      for (i = 0; i < staticModelRange.thetacnt[SEXAV]; i++)
	if (fabs (0.05 - staticModelRange.theta[SEXAV][i]) <= ERROR_MARGIN) break;
      if (i == staticModelRange.thetacnt[SEXAV])
	addTheta (&staticModelRange, THETA_AVG, 0.5);
      
      /* If SexSpecific is on, copy the sex-average (male) thetas for females */
      if (staticModelOptions.mapFlag == SS) {
	for (i = 0; i < staticModelRange.thetacnt[SEXAV]; i++)
	  addTheta (&staticModelRange, THETA_FEMALE, staticModelRange.theta[SEXAV][i]);
      }
    }
    /* Sync param values for hetrozygous genotypes in non-imprinting runs */
    if ((staticModelType.trait != DT) && (staticModelOptions.imprintingFlag != TRUE))
      addConstraint (PARAMC, PEN_dD, 0, 1, EQ, PEN_Dd, 0, 1, FALSE);
  }
  MALCHOKE (staticModelRange.lclassLabels, sizeof (int) * staticModelRange.nlclass, int *);
  for (i = 0; i < staticModelRange.nlclass; i++)
    staticModelRange.lclassLabels[i] = i+1;

  /* Copy our statically-allocated structures over to their global page-allocated
   * counterparts so we can protect them from monkeying.
   */
  memcpy(modelOptions, &staticModelOptions, sizeof(ModelOptions));
  memcpy(modelRange, &staticModelRange, sizeof(ModelRange));
  memcpy(modelType, &staticModelType, sizeof(ModelType));
  return;
}


void finishConfig (ModelRange *modelRange, ModelType *modelType)
{
  /* Sort the values in the final model. Sorted values better support
   * the application of constraints. */
  sortRange (modelRange);

  /* Once sorted, removing duplicates is easy. */
  uniqRange (modelRange);

#if FALSE
  /* Show the unexpanded model. At level 0, all elements are sorted
   * and unique, but we may have:
   *  gfreq: ok
   *  thetas: nonuniform lengths of male/female values
   *  penet: nonuniform lengths of values by allele, lclass=0
   *  param: nonuniform lengths of values by dimension, lclass=0, allele=0 */
  showRange (modelRange, modelType, 0);
  /* Show the constraints. */
  showConstraints ();
#endif

  /* Expand the model, honoring constraints. */
  expandRange (modelRange);
#if FALSE
  /* Show the partially expanded model. At level 1, following
   * expandRange(), we will have refined the model specification while
   * enforcing all specified constraints that do not involve liability
   * classes. */
  showRange (modelRange, modelType, 1);
#endif

  /* Expand the liability classes, but only if necessary and always
   * honoring inter-class constraints. */
  if (staticModelRange.nlclass > 1) {
    expandClassThreshold (modelRange);
    expandClassPenet (modelRange);
  }
  
  //#if FALSE
  /* At level 2, all constraints (including those between classes) are
   * honored, but penet[][][], param[][][][] are not yet fully
   * "factored". */
  //showRange (modelRange, modelType, 2);
  //#endif

  /* Tidy up */
  cleanupRange ();

  return;
}


int set_optionfile (char **toks, int numtoks, void *filename)
{
  if (numtoks < 2)
    bail ("missing filename argument to directive '%s'\n", toks[0]);
  if (numtoks > 2)
    bail ("extra arguments to directive '%s'\n", toks[0]);
  strcpy ((char *) filename, toks[1]);
  return (0);
}


int set_flag (char **toks, int numtoks, void *flag)
{
  if (numtoks > 1)
    bail ("extra arguments to directive '%s'\n", toks[0]);
  *((int *) flag) = TRUE;
  return (0);
}


int clear_flag (char **toks, int numtoks, void *flag)
{
  if (numtoks > 1)
    bail ("extra arguments to directive '%s'\n", toks[0]);
  *((int *) flag) = FALSE;
  return (0);
}


int set_int (char **toks, int numtoks, void *field)
{
  int value;
  char *ptr = NULL;

  if (numtoks < 2)
    bail ("missing integer argument to directive '%s'\n", toks[0]);
  if (numtoks > 2)
    bail ("extra arguments to directive '%s'\n", toks[0]);
  value = (int) strtol (toks[1], &ptr, 10);
  if ((toks[1] == ptr) || (*ptr != '\0'))
    bail ("directive '%s' requires an integer argument\n", toks[0]);
  *((int *) field) = value;
  return (0);
}


int set_traitPositions (char **toks, int numtoks, void *unused)
{
  int numvals, va, vb;
  st_valuelist *vlist;
  double val;

  if (numtoks < 2)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if ((numvals = expandVals (&toks[1], numtoks-1, NULL, &vlist)) <= 0)
    bail ("illegal argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++) {
    if (vlist[va].type == VL_VALUE) {
      addTraitLocus (&staticModelRange, vlist[va].vun.val);
    } else if (vlist[va].type == VL_RANGE) {
      vb = 0;
      while ((val = vlist[va].vun.range.start + (vb++ * vlist[va].vun.range.incr)) <=
	     vlist[va].vun.range.end)
	addTraitLocus (&staticModelRange, val);
    } else if (vlist[va].type == VL_RANGE_SYMBEND) {
      staticModelRange.tlocRangeStart = vlist[va].vun.range.start;
      staticModelRange.tlocRangeIncr = vlist[va].vun.range.incr;
    } else if ((vlist[va].type == VL_SYMBOL) && (vlist[va].vun.symbol == VL_SYM_MARKER))
      staticModelRange.tlmark = TRUE;
  }
  free (vlist);
  observed.traitPositions = 1;
  return (0);
}


int set_alleleFreq (char **toks, int numtoks, void *unused)
{
  int numvals, va=0;
  double *vals;

  if (numtoks < 2)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if ((numvals = expandVals (&toks[1], numtoks-1, &vals, NULL)) <= 0)
    bail ("illegal argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++)
    addAlleleFreq (&staticModelRange, vals[va]);
  free (vals);
  return (0);
}


int set_geneFreq (char **toks, int numtoks, void *unused)
{
  int numvals, va=0;
  double *vals;

  if (numtoks < 2)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if ((numvals = expandVals (&toks[1], numtoks-1, &vals, NULL)) <= 0)
    bail ("illegal argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++)
    addGeneFreq (&staticModelRange, vals[va]);
  free (vals);
  return (0);
}


int set_dprime (char **toks, int numtoks, void *unused)
{
  int numvals, va=0;
  double *vals;

  if (numtoks < 2)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if ((numvals = expandVals (&toks[1], numtoks-1, &vals, NULL)) <= 0)
    bail ("illegal argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++)
    addDPrime (&staticModelRange, vals[va]);
  free (vals);
  return (0);
}


int set_theta (char **toks, int numtoks, void *unused)
{
  int numvals, va=0, type=0;
  double *vals;

  if (numtoks < 2)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if ((numvals = expandVals (&toks[1], numtoks-1, &vals, NULL)) <= 0)
    bail ("illegal argument to directive '%s'\n", toks[0]);
  type = THETA_AVG;
  for (va = 0; va < numvals; va++)
    addTheta (&staticModelRange, type, vals[va]);
  free (vals);
  return (0);
}


int set_alpha (char **toks, int numtoks, void *unused)
{
  int numvals, va=0;
  double *vals;

  if (numtoks < 2)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if ((numvals = expandVals (&toks[1], numtoks-1, &vals, NULL)) <= 0)
    bail ("illegal argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++)
    addAlpha (&staticModelRange, vals[va]);
  free (vals);
  return (0);
}


int set_penetrance (char **toks, int numtoks, void *unused)
{
  int numvals, va=0, geno;
  double *vals;

  if (numtoks < 3)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if ((geno = lookup_modelparam (toks[1])) == -1)
    bail ("illegal model parameter argument to directive '%s'\n", toks[0]);
  if ((numvals = expandVals (&toks[2], numtoks-2, &vals, NULL)) <= 0)
    bail ("illegal argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++)
    addPenetrance (&staticModelRange, geno-PEN_DD, vals[va]);
  observed.penetrance |= PARM_BARE;
  free (vals);
  return (0);
}


int set_constraint (char **toks, int numtoks, void *unused)
{
  int len, first=2, type=-1;
  int oper=0, geno1=0, geno2=0, class1, class2, disjunct=0, maxclass;

  //printf ("set_constraint:");
  //for (oper = 1; oper < numtoks; oper++)
  //  printf (" %s", toks[oper]);
  //printf (", first %d, numtoks %d\n", first, numtoks);
  
  /* We'll just preemptively set this here, since we either return success or die */
  observed.constraints = 1;
  if (numtoks < 2)
    bail ("missing argument to directive '%s'\n", toks[0]);
  len = strlen (toks[first]);
  if ((strncasecmp (toks[1], PENETRANCE_STR, len) == 0) || 
      (strncasecmp (toks[1], MEAN_STR, len) == 0) ||
      (strncasecmp (toks[1], DEGOFFREEDOM_STR, len) == 0)) {
    while (1) {
      class1 = class2 = 0;
      if ((numtoks >= first + 3) &&
	  ((geno1 = lookup_modelparam (toks[first])) != -1) &&
	  ((oper = lookup_comparator (toks[first+1])) != -1) &&
	  ((geno2 = lookup_modelparam (toks[first+2])) != -1)) {
	if ((type != -1) && (type != SIMPLE))
	  bail ("illegal combination of %s and Simple constraints\n", contype_strs[type]);
	type = SIMPLE;
	first += 3;
	//printf ("  SIMPLE: %s %s %s, first %d\n", mp_strs[geno1], op_strs[oper], mp_strs[geno2],
	//first);
	
      } else if ((numtoks >= first + 5) &&
		 ((geno1 = lookup_modelparam (toks[first])) != -1) &&
		 ((class1 = singleDigit (toks[first+1])) > 0) &&
		 ((oper = lookup_comparator (toks[first+2])) != -1) &&
		 ((geno2 = lookup_modelparam (toks[first+3])) != -1) &&
		 ((class2 = singleDigit (toks[first+4])) > 0)) {
	
	if ((type != -1) && (type != CLASSC))
	  bail ("illegal combination of %s and Liability Class constraints\n", contype_strs[type]);
	type = CLASSC;
	first += 5;
	//printf ("  CLASSC: %s %d %s %s %d, first %d\n", mp_strs[geno1], class1, op_strs[oper],
	//mp_strs[geno2], class2, first);
      } else 
	bail ("illegal argument to directive '%s'\n", toks[0]);

      addConstraint (type, geno1, class1, 0, oper, geno2, class2, 0, disjunct);
      maxclass = (class1 > class2) ? class1 : class2;
      observed.maxclass = (observed.maxclass > maxclass) ? observed.maxclass : maxclass;
      disjunct = 1;
      if (numtoks <= first)
	return (0);
      if (strcmp (toks[first++], ",") != 0)
	bail ("illegal argument to directive '%s'\n", toks[0]);
    }

  } else if (strncasecmp (toks[1], STANDARDDEV_STR, len) == 0) {
    while (1) {
      class1 = class2 = 0;
      if ((numtoks >= first + 3) &&
	  ((geno1 = lookup_modelparam (toks[first])) != -1) &&
	  ((oper = lookup_comparator (toks[first+1])) != -1) &&
	  ((geno2 = lookup_modelparam (toks[first+2])) != -1)) {
	if ((type != -1) && (type != PARAMC))
	  bail ("illegal combination of %s and Parameter constraints\n", contype_strs[type]);
	type = PARAMC;
	first += 3;
	//printf ("  PARAMC: %s %s %s, first %d\n", mp_strs[geno1], op_strs[oper], mp_strs[geno2],
	//first);
	
      } else if ((numtoks >= first + 5) &&
		 ((geno1 = lookup_modelparam (toks[first])) != -1) &&
		 ((class1 = singleDigit (toks[first+1])) > 0) &&
		 ((oper = lookup_comparator (toks[first+2])) != -1) &&
		 ((geno2 = lookup_modelparam (toks[first+3])) != -1) &&
		 ((class2 = singleDigit (toks[first+4])) > 0)) {
	
	if ((type != -1) && (type != PARAMCLASSC))
	  bail ("illegal combination of %s and Liability Class constraints\n", contype_strs[type]);
	type = PARAMCLASSC;
	first += 5;
	//printf ("  PARAMCLASSC: %s %d %s %s %d, first %d\n", mp_strs[geno1], class1,
	//op_strs[oper], mp_strs[geno2], class2, first);
      } else 
	bail ("illegal argument to directive '%s'\n", toks[0]);

      addConstraint (type, geno1, class1, 1, oper, geno2, class2, 1, disjunct);
      maxclass = (class1 > class2) ? class1 : class2;
      observed.maxclass = (observed.maxclass > maxclass) ? observed.maxclass : maxclass;
      disjunct = 1;
      if (numtoks <= first)
	return (0);
      if (strcmp (toks[first++], ",") != 0)
	bail ("illegal argument to directive '%s'\n", toks[0]);
    }

  } else if (strncasecmp (toks[1], THRESHOLD_STR, len) == 0) {
    class1 = class2 = 0;
    while (1) {
      if (! ((numtoks >= first + 3) &&
	     ((class1 = singleDigit (toks[first])) > 0) &&
	     ((oper = lookup_comparator (toks[first+1])) != -1) &&
	     ((class2 = singleDigit (toks[first+2])) > 0)))
	bail ("illegal arguments to directive '%s'\n", toks[0]);
      addConstraint (SIMPLE, THRESHOLD, class1, 0, oper, THRESHOLD, class2, 0, disjunct);
      maxclass = (class1 > class2) ? class1 : class2;
      observed.maxclass = (observed.maxclass > maxclass) ? observed.maxclass : maxclass;
      disjunct = 1;
      first += 3;
      if (numtoks <= first) {
	observed.constraints = 1;
	return (0);
      }
      if (strcmp (toks[first++], ",") != 0)
	bail ("illegal argument to directive '%s'\n", toks[0]);
    }
  } else 
    bail ("illegal argument to directive '%s'\n", toks[0]);
  
  /* we'll never get this far, but just to shut GCC up... */
  return (0);
}


int set_multipoint (char **toks, int numtoks, void *unused)
{
  int value;
  char *ptr = NULL;

  if (numtoks < 2)
    bail ("missing integer argument to directive '%s'\n", toks[0]);
  if (numtoks > 2)
    bail ("extra arguments to directive '%s'\n", toks[0]);
  staticModelType.type = MP;
  value = (int) strtol (toks[1], &ptr, 10);
  if ((toks[1] == ptr) || (*ptr != '\0'))
    bail ("directive '%s' requires an integer argument\n", toks[0]);
  staticModelType.numMarkers = value;
  return (0);
}


int set_markerAnalysis (char **toks, int numtoks, void *unused)
{
  if (numtoks < 2)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if (numtoks > 2)
    bail ("extra arguments to directive '%s'\n", toks[0]);
  if (strncasecmp (toks[1], "All", strlen (toks[1])) == 0)
    staticModelOptions.markerAnalysis = MM;
  else if (strncasecmp (toks[1], "Adjacent", strlen (toks[1])) == 0)
    staticModelOptions.markerAnalysis = AM;
  else
    bail ("unknown argument to directive '%s'\n", toks[1]);
  return (0);
}


int set_mapFlag (char **toks, int numtoks, void *unused)
{
  if (numtoks > 1)
    bail ("extra arguments to directive '%s'\n", toks[0]);
  staticModelOptions.mapFlag = SS;
  return (0);
}


int set_disequilibrium (char **toks, int numtoks, void *unused)
{
  if (numtoks > 1)
    bail ("extra arguments to directive '%s'\n", toks[0]);
  staticModelOptions.equilibrium = LINKAGE_DISEQUILIBRIUM;
  return (0);
}


int set_quantitative (char **toks, int numtoks, void *unused)
{
  int numvals;
  double *vals=NULL;
  
  if (numtoks < 2)
    bail ("missing arguments to directive '%s'\n", toks[0]);

  if (strcasecmp (toks[0], QT_STR) == 0) {
    staticModelType.trait = QT;
  } else if (strcasecmp (toks[0], QTT_STR) == 0) {
    staticModelType.trait = CT;
  } else {
    bail ("set_quantitative called with bad directive '%s'\n", toks[0]);
  }

  if (strcasecmp (toks[1], "normal") == 0) {
    staticModelType.distrib = QT_FUNCTION_T;
    /* I think this is degrees of freedom; anyway, YH sez: fix it at 30 */
    REALCHOKE(staticModelType.constants, 1 * sizeof (int), void *);
    staticModelType.constants[0] = 30;
    staticModelRange.npardim = 1;
    if (numtoks == 2)
      return (0);
    if ((numtoks != 3) && ((numvals = expandVals (&toks[2], numtoks-2, &vals, NULL)) != 2))
      bail ("illegal arguments to directive '%s'\n", toks[0]);
    staticModelType.mean = vals[0];
    staticModelType.sd = vals[1];
    free (vals);
  } else if (strcasecmp (toks[1], "chisq") == 0) {
    if (numtoks > 2) 
      bail ("illegal arguments to directive '%s'\n", toks[0]);
    staticModelType.distrib = QT_FUNCTION_CHI_SQUARE;
    staticModelType.mean = 0;
    staticModelType.sd = 1;
    /* A non-empty range for this parameter triggers a loop elsewhere */
    staticModelRange.npardim = 1;
    addParameter (&staticModelRange, 0, 1.0);
  } else
    bail ("illegal arguments to directive '%s'\n", toks[0]);
  return (0);
}


/* qt_set_mean: two forms are acceptable: no trait genotype and exactly
 * two values (the form for 'advanced' users under dynamic sampling),
 * or a trait genotype and at least one value (the form for fixed models).
 */
int set_qt_mean (char **toks, int numtoks, void *unused)
{
  int numvals, va=0, geno;
  double *vals;

  if (numtoks < 3)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if ((geno = lookup_modelparam (toks[1])) == -1) {
    if ((numvals = expandVals (&toks[1], numtoks-1, &vals, NULL)) == 0)
      bail ("illegal model parameter argument to directive '%s'\n", toks[0]);
    if (numvals != 2)
      bail ("illegal argument to directive '%s'\n", toks[0]);
    addPenetrance (&staticModelRange, 0, vals[0]);
    addPenetrance (&staticModelRange, 0, vals[1]);
    observed.mean |= PARM_BARE;

  } else {
    if ((numvals = expandVals (&toks[2], numtoks-2, &vals, NULL)) <= 0)
      bail ("illegal argument to directive '%s'\n", toks[0]);
    for (va = 0; va < numvals; va++)
      addPenetrance (&staticModelRange, geno-PEN_DD, vals[va]);
    observed.mean |= PARM_GENO;
  }
  free (vals);
  return (0);
}


/* Basically the same as set_qt_mean, above. */
int set_qt_degfreedom (char **toks, int numtoks, void *unused)
{
  int numvals, va=0, geno;
  double *vals;

  if (numtoks < 3)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if ((geno = lookup_modelparam (toks[1])) == -1) {
    if ((numvals = expandVals (&toks[1], numtoks-1, &vals, NULL)) == 0)
      bail ("illegal model parameter argument to directive '%s'\n", toks[0]);
    if (numvals != 2)
      bail ("illegal argument to directive '%s'\n", toks[0]);
    addPenetrance (&staticModelRange, 0, vals[0]);
    addPenetrance (&staticModelRange, 0, vals[1]);
    observed.degOfFreedom |= PARM_BARE;

  } else {
    if ((numvals = expandVals (&toks[2], numtoks-2, &vals, NULL)) <= 0)
      bail ("illegal argument to directive '%s'\n", toks[0]);
    for (va = 0; va < numvals; va++)
      addPenetrance (&staticModelRange, geno-PEN_DD, vals[va]);
    observed.degOfFreedom |= PARM_GENO;
  }
  free (vals);
  return (0);
}


int set_qt_standarddev (char **toks, int numtoks, void *unused)
{
  int numvals, va=0;
  double *vals;

  if (numtoks < 2)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if ((numvals = expandVals (&toks[1], numtoks-1, &vals, NULL)) <= 0)
    bail ("illegal argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++)
    addParameter (&staticModelRange, 0, vals[va]);
  observed.standardDev |= PARM_BARE;
  free (vals);
  return (0);
}


int set_qt_threshold (char **toks, int numtoks, void *unused)
{
  int numvals, va;
  double *vals=NULL;

  if (numtoks < 2)
    bail ("missing arguments to directive '%s'\n", toks[0]);
  if ((numvals = expandVals (&toks[1], numtoks-1, &vals, NULL)) <= 0)
    bail ("illegal argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++)
    addTraitThreshold (&staticModelRange, vals[va]);
  free (vals);
  return (0);
}


int set_qt_truncation (char **toks, int numtoks, void *unused)
{
  int va=1;
  double val;
  char *ca;

  if (numtoks < 2)
    bail ("missing arguments to directive '%s'\n", toks[0]);
  while (1) {
    if (va + 2 > numtoks)
      bail ("illegal arguments to directive '%s'\n", toks[0]);
    val = strtod (toks[va+1], &ca);
    if ((ca == NULL) || (*ca != '\0'))
      bail ("illegal arguments to directive '%s'\n", toks[0]);
    if (strcasecmp (toks[va], "left") == 0) {
      staticModelType.minOriginal = val;
    } else if (strcasecmp (toks[va], "right") == 0) {
      staticModelType.maxOriginal = val;
    } else {
      bail ("illegal arguments to directive '%s'\n", toks[0]);
    }
    if ((va += 2) >= numtoks)
      break;
    if (strcmp (toks[va++], ",") != 0)
      bail ("illegal arguments to directive '%s'\n", toks[0]);
  }
  return (0);
}

#ifdef STUDYDB
int set_study_parameters (char **toks, int numtoks, void *unused)
{
  char *ptr = NULL;
  // Want studyId, role, hostname, dbName, username, password
  if (numtoks != 7)
    bail ("inappropriate number of arguments to directive '%s'\n", toks[0]);
  studyDB.inStudyId = (int) strtol (toks[1], &ptr, 10);
  if ((toks[1] == ptr) || (*ptr != '\0'))
    bail ("directive '%s' requires an initial integer argument\n", toks[0]);
  strncpy (studyDB.role, toks[2], sizeof (studyDB.role));
  strncpy (studyDB.hostname, toks[3], sizeof (studyDB.hostname));
  strncpy (studyDB.dBName, toks[4], sizeof (studyDB.dBName));
  strncpy (studyDB.username, toks[5], sizeof (studyDB.username));
  strncpy (studyDB.password, toks[6], sizeof (studyDB.password));
  return (0);
}
#else
int set_study_parameters (char **toks, int numtoks, void *unused)
{
  // Completely ignore it
  return (0);
}
#endif

int set_affectionStatus (char **toks, int numtoks, void *unused)
{
  int numvals;
  double *vals=NULL;

  if (numtoks < 2)
    bail ("missing arguments to directive '%s'\n", toks[0]);
  if ((numvals = expandVals (&toks[1], numtoks-1, &vals, NULL)) == 1) {
    /* This is legal for QT analyses, just to set the 'undefined' pheno code */
    staticModelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = vals[0];
  } else if (numvals == 3) {
    /* Everything else (DT, CT) requires three pheno codes */
    staticModelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = vals[0];
    staticModelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] = vals[1];
    staticModelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED] = vals[2];
  } else
    bail ("illegal arguments to directive '%s'\n", toks[0]);
  
  free (vals);
  return (0);
}


int set_resultsprefix (char **toks, int numtoks, void *unused)
{
  int len;

  if (numtoks < 2)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if (numtoks > 2)
    bail ("extra arguments to directive '%s'\n", toks[0]);
  if ((len = strlen (toks[1])) > PATH_MAX - 1)
    bail ("argument to directive '%s' is too long\n", toks[0]);
  strcpy (staticModelOptions.resultsprefix, toks[1]);
  if (staticModelOptions.resultsprefix[len-1] != '/')
    strcat (staticModelOptions.resultsprefix, "/");
  return (0);
}


int skip_analysis (char **toks, int numtoks, void *unused)
{
  WARNING ("Skipping analysis due to directive");
  return (0);
}


int noop (char **toks, int numtoks, void *unused)
{
  return (0);
}


int expandVals (char **toks, int numtoks, double **vals_h, st_valuelist **vlist_h)
{
  int numvals=0, listsize=10, tokidx=0, va;
  char *ca = NULL, *cb = NULL;
  double start, end=-1, incr, val, *vals=NULL;
  st_valuelist *vlist=NULL;

  /* Sanity check */
  if (numtoks == 0)
    return (0);

  /* Exactly one return value pointer must be non-NULL */
  if (((vlist_h == NULL) && (vals_h == NULL)) || ((vlist_h != NULL) && (vals_h != NULL)))
    return (-1);
  
  if (vlist_h != NULL)
    MALCHOKE (vlist, sizeof (st_valuelist) * listsize, void *);
  else
    MALCHOKE (vals, sizeof (double) * listsize, void *);

  while (1) {
    if (numvals >= listsize) {
      if (((vlist != NULL) &&
	   ((vlist = realloc (vlist, sizeof (st_valuelist) * (listsize += 10))) == NULL)) || 
	  ((vals = realloc (vals, sizeof (double) * (listsize += 10))) == NULL))
	ERROR ("Complex config parsing realloc failed");
    }
    
    /* Skip the first character of the token to avoid the leading '-' of a negative number */
    if ((cb = index (toks[tokidx]+1, '-')) == NULL) {
      /* This should be a single value or symbol */
      if (vlist != NULL) {
	/* returning a list of st_valuelist */
	if (strcasecmp (toks[tokidx], "Marker") == 0) {
	  vlist[numvals].type = VL_SYMBOL;
	  vlist[numvals++].vun.symbol = VL_SYM_MARKER;
	} else {
	  vlist[numvals].type = VL_VALUE;
	  vlist[numvals++].vun.val = strtod (toks[tokidx], &cb);
	  if ((cb == toks[tokidx]) || (*cb != '\0'))
	    break;
	}
      } else {
	/* returning a list of doubles */
	vals[numvals++] = strtod (toks[tokidx], &cb);
	if ((cb == toks[tokidx]) || (*cb != '\0'))
	  break;
      }

    } else {
      /* This should be a range, either 'i-j:k' or 'i-end:k' */
      *cb = '\0';
      // printf ("  range format: first substring '%s', remainder, '%s'\n", toks[tokidx], cb+1);
      start = strtod (toks[tokidx], &cb);
      if ((cb == toks[tokidx]) || (*cb != '\0'))
	break;
      // printf ("  start value: %f\n", start);
      ca = cb + 1;
      if ((*ca == '\0') || ((cb = index (ca, ':')) == NULL))
	break;
      // printf ("                next substring '%s', remainder, '%s'\n", ca, cb+1);
      *cb = '\0';
      if (vlist != NULL) {
	if (strcasecmp (ca, "End") == 0) {
	  vlist[numvals].type = VL_RANGE_SYMBEND;
	} else {
	  vlist[numvals].type = VL_RANGE;
	  end = strtod (ca, &cb);
	  if ((cb == ca) || (*cb != '\0'))
	    break;
	}
      } else {
	end = strtod (ca, &cb);
	if ((cb == ca) || (*cb != '\0'))
	  break;
	// printf ("  end value: %f\n", end);
      }
      ca = cb + 1;
      if (*ca == '\0')
	break;
      // printf ("                last substring '%s'\n", ca);
      incr = strtod (ca, &cb);
      if (*cb != '\0')
	break;
      // printf ("  incr value: %f\n", incr);
      
      if (vlist != NULL) {
	if (vlist[numvals].type == VL_RANGE) {
	  if (start > end)
	    break;
	  vlist[numvals].vun.range.end = end;
	}
	if (incr <= 0)
	  break;
	vlist[numvals].vun.range.start = start;
	vlist[numvals++].vun.range.incr = incr;
      } else {
	if ((start > end) || (incr <= 0))
	  break;
	va = 0;
	while ((val = start + (va++ * incr)) <= (end + ERROR_MARGIN)) {
	  if ((numvals >= listsize) && 
	      ((vals = realloc (vals, sizeof (double) * (listsize += 10))) == NULL))
	    ERROR ("Complex config parsing realloc failed");
	  vals[numvals++] = val;
	}
      }
    }
    if (++tokidx >= numtoks) {
      if (vlist_h != NULL)
	*vlist_h = vlist;
      else 
	*vals_h = vals;
      return (numvals);
    }
    if (strcmp (toks[tokidx++], ",") != 0)
      break;
  }

  if (vlist != NULL)
    free (vlist);
  else 
    free (vals);
  return (-1);
}


/* A binary search, so the dispatch table must be sorted. Given a search key of
 * length n, a table entry is a match if the search key and the table key are
 * identical for the first n characters, and no other table keys are also identical
 * for the first n characters; or, if the table key is also length n, and the search
 * key is identical to the table key.
 */
int lookupDispatch (char *key, st_dispatch *table)
{
  int tablen, keylen, hi, lo, mid, res;

  tablen = sizeof (dispatchTable) / sizeof (st_dispatch);
  keylen = strlen (key);

  lo = 0;
  hi = tablen - 1;
  mid = (hi + lo) / 2;
  
  while (hi >= lo) {
    if ((res = strncasecmp (key, table[mid].key, keylen)) == 0)
      break;
    else if (res < 0)
      hi = mid - 1;
    else
      lo = mid + 1;
    mid = (hi + lo) / 2;
  }
  if (res != 0)
    return (-1);

  lo = hi = mid;
  while ((lo > 0) && (strncasecmp (key, table[lo-1].key, keylen) == 0))
    lo--;
  while ((hi < tablen -1) && (strncasecmp (key, table[hi+1].key, keylen) == 0))
    hi++;
  if (lo == hi)
    return (mid);
  for (mid = lo; mid <= hi; mid++)
    if (strlen (table[mid].key) == keylen)
      return (mid);
  return (-2);
}


/* We use this two ways: when we initially qsort the dispatch table, and when we
 * use bsearch to find an entry in the table. 
 */
int compareDispatch (const void *a, const void *b)
{
  int len;

  len = strlen (((st_dispatch *)a)->key);
  return (strncasecmp (((st_dispatch *)a)->key, ((st_dispatch *)b)->key, len));
}


/* Reads input from fp, using permuteLine to collapse blank and comment
 * lines, and splits input into substrings at newlines or semicolons.
 * The substrings are tokenized and returned as a list of tokens. Maintains
 * context in between calls, so the tokens from multiple semicolon-
 * separated substrings will be returned on successive calls. 
 */
int getNextTokgroup (FILE *fp, char ***tokgroup_h, int *tokgroupsize)
{
  int numtoks;
  char *semicolon;
  
  while (1) {
    if (buffptr == NULL || strlen (buffptr) == 0) {
      if (fgets (buff, BUFFSIZE, fp) == 0) {
	lineno = -1;
	return (0);
      }
      lineno++;
      permuteLine (buffptr = buff, BUFFSIZE);
    }
    if ((semicolon = index (buffptr, ';')) != NULL)
      *(semicolon++) = '\0';
    numtoks = tokenizeLine (buffptr, tokgroup_h, tokgroupsize);
    buffptr = semicolon;
    if (numtoks == 0)
      continue;
    return (numtoks);
  }
}


/* Splits a character string into a list of tokens. The line should already
 * have been passed through permuteLine(), so extra whitespace has already been
 * trimmed.
 */
int tokenizeLine (char *line, char ***tokgroup_h, int *tokgroupsize)
{
  int numtoks=0;
  char **tokgroup, *ca, *cb = NULL;
  
  if ((tokgroup = *tokgroup_h) == NULL) {
    *tokgroupsize = 10;
    MALCHOKE (tokgroup, sizeof (char *) * *tokgroupsize, void *);
  }
  
  ca = strtok_r (line, " ", &cb);
  while (ca != NULL) {
    if (numtoks + 1 >= *tokgroupsize) {
      *tokgroupsize += 10;
      REALCHOKE (tokgroup, sizeof (char *) * *tokgroupsize, void *);
    }
    tokgroup[numtoks++] = ca;
    ca = strtok_r (NULL, " ", &cb);
  }

  *tokgroup_h = tokgroup;
  return (numtoks);
}


int singleDigit (char *str)
{
  int digit;
  char *ca;

  if ((*str != '\0') && ((digit = (int) strtol (str, &ca, 10)) >= 1) &&
      (*ca == '\0') && (digit <= 9))
    return (digit);
  return (-1);
}

/* Dumps a subset of the fields in staticModelOptions.
 */
void dumpmodelOptions (ModelOptions *mo)
{
  printf ("%18s : %s\n", "markerfile", mo->markerfile);
  printf ("%18s : %s\n", "mapfile", mo->mapfile);
  printf ("%18s : %s\n", "pedfile", mo->pedfile);
  printf ("%18s : %s\n", "datafile", mo->datafile);
  printf ("%18s : %s\n", "avghetfile", mo->avghetfile);
  printf ("%18s : %s\n", "pplfile", mo->pplfile);
  //  printf ("%18s : %s\n", "condFile", mo->condFile);
  printf ("%18s : %s\n", "ccfile", mo->ccfile);
  printf ("%18s : %s\n", "modfile", mo->modfile);
  printf ("%18s : %s\n", "maxmodelfile", mo->maxmodelfile);
  printf ("%18s : %s\n", "intermediatefile", mo->intermediatefile);
  printf ("%18s : %s\n", "dkelvinoutfile", mo->dkelvinoutfile);
}


void bail (char *fmt, char *arg)
{
  char newfmt[BUFFSIZE];

  if (lineno > 0) {
    strcpy (newfmt, "'%s' line %d: ");
    strcat (newfmt, fmt);
    ERROR (newfmt, conffilename, lineno, arg);
  } else {
    strcpy (newfmt, "on command line: ");
    strcat (newfmt, fmt);
    ERROR (newfmt, arg);
  }
}
