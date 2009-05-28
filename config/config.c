#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include "config.h"
#include "../utils/utils.h"
#include "../pedlib/pedlib.h"

#define BUFFSIZE 256

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
#define MALETHETA_STR     "MaleTheta"
#define FEMALETHETA_STR   "FemaleTheta"
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
  int sexSpecificThetas,
    traitLoci,
    constraints,
    penetrance,
    mean,
    standarddev,
    degfreedom;
} st_observed;

ModelOptions modelOptions;
ModelRange modelRange;
ModelType modelType;


/* Globals */
char buff[BUFFSIZE] = "";   /* These two are global to provide context to getNextTokgroup */
char *buffptr = NULL;
st_observed observed;       /* track non-obvious directives */
char *conffilename=NULL;
int lineno = 0;             /* so failure messages are more useful to users */

/* prototypes for non-public routines */
void initializeDefaults ();
int expandVals (char **toks, int numtoks, double **vals_h, st_valuelist **vlist_h);
int lookupDispatch (char *key, st_dispatch *table);
int compareDispatch (const void *a, const void *b);
int getNextTokgroup (FILE *fp, char ***tokgroup_h, int *tokgroupsize);
int tokenizeLine (char *line, char ***tokgroup_h, int *tokgroupsize);
int singleDigit (char *str);
void dumpModelOptions (ModelOptions *mo);
void bail (char *fmt, char *arg);

/* functions for use in the dispatch table */
int set_optionfile (char **toks, int numtoks, void *filename);
int set_flag (char **toks, int numtoks, void *flag);
int clear_flag (char **toks, int numtoks, void *flag);
int set_int (char **toks, int numtoks, void *field);

int set_traitLoci (char **toks, int numtoks, void *unused);
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
int set_qt_standarddev (char **toks, int numtoks, void *unused);
int set_qt_degfreedom (char **toks, int numtoks, void *unused);
int set_qt_threshold (char **toks, int numtoks, void *unused);
int set_qt_truncation (char **toks, int numtoks, void *unused);
int set_affectionStatus (char **toks, int numtoks, void *unused);
int set_resultsprefix (char **toks, int numtoks, void *unused);
int set_logLevel (char **toks, int numtoks, void *unused);


st_dispatch dispatchTable[] = { {"FrequencyFile", set_optionfile, &modelOptions.markerfile},
				{"MapFile", set_optionfile, &modelOptions.mapfile},
				{"PedigreeFile", set_optionfile, &modelOptions.pedfile},
				{"LocusFile", set_optionfile, &modelOptions.datafile},
				{"BayesRatioFile", set_optionfile, &modelOptions.avghetfile},
				{"PPLFile", set_optionfile, &modelOptions.pplfile},
				{"CountFile", set_optionfile, &modelOptions.ccfile},
				{"MODFile", set_optionfile, &modelOptions.modfile},
				{"NIDetailFile", set_optionfile, &modelOptions.dkelvinoutfile},

				{"NonPolynomial", clear_flag, &modelOptions.polynomial},
				{"Imprinting", set_flag, &modelOptions.imprintingFlag},
				{"SexLinked", set_flag, &modelOptions.sexLinked},
				{"FixedModels", clear_flag, &modelOptions.integration},
				{"DryRun", set_flag, &modelOptions.dryRun},
				{"ExtraMODs", set_flag, &modelOptions.extraMODs},

				{"PolynomialScale", set_int, &modelOptions.polynomialScale},
				{"LiabilityClasses", set_int, &modelRange.nlclass},
				{"DiseaseAlleles", set_int, &modelRange.nalleles},

				{"TraitLoci", set_traitLoci, NULL},
				{"DiseaseGeneFrequency", set_geneFreq, NULL},
				{"DPrime", set_dprime, NULL},
				{THETA_STR, set_theta, NULL},
				{MALETHETA_STR, set_theta, NULL},
				{FEMALETHETA_STR, set_theta, NULL},
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
				{STANDARDDEV_STR, set_qt_standarddev, NULL},
				{DEGOFFREEDOM_STR, set_qt_degfreedom, NULL},
				{THRESHOLD_STR, set_qt_threshold, NULL},
				{"Truncate", set_qt_truncation, NULL},
				{"PhenoCodes", set_affectionStatus, NULL},
				{"SurfacesPath", set_resultsprefix, NULL},
				/*{"condfile", set_condrun, &modelOptions.condFile},*/
				{"Log", set_logLevel, NULL}
};


#if 0
main (int argc, char *argv[])
{
  logInit ();

  if (argc == 1)
    logMsg (LOGDEFAULT, LOGFATAL, "%s: no configuration file specified\n", argv[0]);
  initializeDefaults ();
  my_readConfigFile (argv[1]);
  if (argc > 2) {
    parseCommandLine (argc-2, &argv[2]);
  }
  dumpModelOptions (&modelOptions);
  validateConfig ();
  finishConfig ();
}
#endif 


void initializeDefaults ()
{
  /* Initialize the the global configuration structures to default values */
  memset (&modelOptions, 0, sizeof (ModelOptions));
  memset (&modelRange, 0, sizeof (ModelRange));
  memset (&modelType, 0, sizeof (ModelType));

  strcpy (modelOptions.markerfile, DEFAULTMARKERFILENAME);
  strcpy (modelOptions.mapfile, DEFAULTMAPFILENAME);
  strcpy (modelOptions.pedfile, DEFAULTPEDFILENAME);
  strcpy (modelOptions.datafile, DEFAULTDATAFILENAME);
  strcpy (modelOptions.avghetfile, DEFAULTAVGHETFILENAME);
  strcpy (modelOptions.condFile, DEFAULTCONDFILENAME);
  strcpy (modelOptions.resultsprefix, DEFAULTRESULTSPREFIX);

  modelOptions.sUnknownPersonID = malloc (sizeof (char) * 2);
  strcpy (modelOptions.sUnknownPersonID, "0");

  modelOptions.equilibrium = LINKAGE_EQUILIBRIUM;
  modelOptions.markerAnalysis = FALSE;
  modelOptions.saveResults = FALSE;
  modelOptions.polynomial = TRUE;
  modelOptions.integration = TRUE;
  modelOptions.imprintingFlag = FALSE;
  modelOptions.mapFlag = SA;
  modelOptions.sexLinked = FALSE;
  modelOptions.dryRun = FALSE;
  modelOptions.polynomialScale = 0;
  modelOptions.extraMODs = FALSE;
  modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = -DBL_MAX;
  modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] = -DBL_MAX;
  modelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED] = -DBL_MAX;

  /* Set default values for PPL calculations */
  modelOptions.thetaCutoff[0] = 0.05;   /* LRs when theta < cutoff are weighted heavier */
  modelOptions.thetaCutoff[1] = 0.05;
  modelOptions.thetaWeight = 0.95;      /* Weight for LRs when theta < cutoff */
  modelOptions.prior = 0.02;            /* Prior probability of linkage */
  modelOptions.LDprior = 0.02;          /* Prior probability of LD given close linkage */
  
  modelRange.nalleles = 2;
  modelRange.nlclass = 1;
  modelRange.npardim = 0;
  modelRange.nlambdas = 0;
  modelRange.maxnlambdas = 0;
  modelRange.tlocRangeStart = -1;
  modelRange.tlocRangeIncr = -1;
  modelRange.tlmark = FALSE;

  modelType.type = TP;
  modelType.trait = DT;
  modelType.distrib = -1;
  /* set default for QT */
  modelType.minOriginal = -999999999.00;
  modelType.maxOriginal = 999999999.00;

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
    logMsg (LOGDEFAULT, LOGFATAL, "open '%s' failed, %s\n", config, strerror (errno));

  while ((numtoks = getNextTokgroup (conffp, &toks, &tokgroupsize)) > 0) {
    for (va = 0; va < numtoks; va++) {
      printf ("tok %d: %s\n", va, toks[va]);
    }
    if ((va = lookupDispatch (toks[0], dispatchTable)) >= 0) {
      printf ("directive '%s' matches at index %d\n", toks[0], va);
      (*dispatchTable[va].parse) (toks, numtoks, dispatchTable[va].hint);
    } else {
      logMsg (LOGDEFAULT, LOGFATAL, "directive '%s' on line %d is %s\n", toks[0],
	      lineno, (va == -1) ? "unknown" : "not unique");
    }
    //printf ("\n");
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
  char buff[BUFFSIZE], **toks=NULL;
  
  if (argc == 0)
    return;

  if (strncmp (argv[0], "--", 2) != 0)
    logMsg (LOGDEFAULT, LOGFATAL, "expected directive on command line, found '%s'\n", argv[0]);
  bufflen = strlen (argv[0] - 2);
  strcpy (buff, argv[0]+2);
  
  curidx = 1;
  while (curidx < argc) {
    if (strncmp (argv[curidx], "--", 2) == 0) {
      permuteLine (buff, BUFFSIZE);
      numtoks = tokenizeLine (buff, &toks, &tokgroupsize);
      if ((va = lookupDispatch (toks[0], dispatchTable)) >= 0) {
	//printf ("directive '%s' matches at index %d\n", toks[0], va);
	(*dispatchTable[va].parse) (toks, numtoks, dispatchTable[va].hint);
      } else
	logMsg (LOGDEFAULT, LOGFATAL, "directive '%s' on command line is %s\n", toks[0],
		(va == -1) ? "unknown" : "not unique");
      bufflen = strlen (argv[0] - 2);
      strcpy (buff, argv[0]+2);

    } else {
      buff[bufflen++] = ' ';
      buff[bufflen++] = '\0';
      bufflen += strlen (argv[curidx]);
      strcat (buff, argv[curidx]);
    }
    curidx++;
  }
    
  permuteLine (buff, BUFFSIZE);
  numtoks = tokenizeLine (buff, &toks, &tokgroupsize);
  if ((va = lookupDispatch (toks[0], dispatchTable)) >= 0) {
    //printf ("directive '%s' matches at index %d\n", toks[0], va);
    (*dispatchTable[va].parse) (toks, numtoks, dispatchTable[va].hint);
  } else
    logMsg (LOGDEFAULT, LOGFATAL, "directive '%s' on command line is %s\n", toks[0],
	    (va == -1) ? "unknown" : "not unique");
  
  if (tokgroupsize > 0)
    free (toks);
  return;
}


void validateConfig ()
{
  /* Check what was explicitly specified against what was implied. Implying model options
   * by setting parameters (ex. implying a sex specific map by configuring male theta
   * values) is no longer acceptable. The user must specify what they want, we won't guess.
   */

  /* Theta-related checks */
  if (modelRange.theta && modelOptions.integration)
    logMsg (LOGDEFAULT, LOGFATAL, "Don't specify Theta values without FixedModels\n");
  if (modelRange.theta && modelType.type == MP)
    logMsg (LOGDEFAULT, LOGFATAL, "Don't specify Theta values for Multipoint analyses\n");
  if (observed.sexSpecificThetas && modelOptions.mapFlag != SS)
    logMsg (LOGDEFAULT, LOGFATAL, "Don't specify sex-specific Thetas without SexSpecific\n");

  /* DPrime-related checks */
  if (modelRange.dprime && modelOptions.integration)
    logMsg (LOGDEFAULT, LOGFATAL, "Don't specify DPrime values without FixedModels\n");
  if (modelRange.dprime && modelType.type == MP)
    logMsg (LOGDEFAULT, LOGFATAL, "Don't specify DPrime values for Multipoint analyses\n");
  if (modelRange.dprime && modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
    logMsg (LOGDEFAULT, LOGFATAL, "Don't specify DPrime values without LD\n");

  /* Other fixed-model-related checks */
  if (modelRange.alpha && modelOptions.integration)
    logMsg (LOGDEFAULT, LOGFATAL, "Don't specify Alpha values without FixedModels\n");
  if (modelRange.gfreq && modelOptions.integration)
    logMsg (LOGDEFAULT, LOGFATAL, "Don't specify DiseaseGeneFrequency values without FixedModels\n");
  if (observed.constraints && modelOptions.integration)
    logMsg (LOGDEFAULT, LOGFATAL, "Don't specify Constraints values without FixedModels\n");

  /* FIXME: Should prolly make sure only two values (hi and lo limits) are specified here */
  if ((modelRange.penet && modelOptions.integration) &&
      !(modelType.distrib == QT_FUNCTION_CHI_SQUARE && observed.degfreedom))
    logMsg (LOGDEFAULT, LOGFATAL, "Don't specify Penetrance values without FixedModels\n");
  if (modelRange.tthresh && modelOptions.integration && modelType.trait != CT)
    logMsg (LOGDEFAULT, LOGFATAL, "Don't specify Threshold values without FixedModels\n");
  /* FIXME: And make sure that two threshold values are specified for dkelvin+QTT */

  /* Things that don't work with Multipoint */
  if (modelOptions.markerAnalysis && modelType.type == MP)
    logMsg (LOGDEFAULT, LOGFATAL, "MarkerToMarker incompatible with Multipoint\n");
  if ((modelOptions.pplfile[0] != '\0') && (modelType.type == MP))
    logMsg (LOGDEFAULT, LOGFATAL, "PPLFile is incompatible with Multipoint\n");
  if (modelOptions.extraMODs && (modelType.type == MP))
    logMsg (LOGDEFAULT, LOGFATAL, "ExtraMods is incompatible with Multipoint\n");
  if (modelType.type == MP && modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM)
    logMsg (LOGDEFAULT, LOGFATAL, "LD is incompatible with Multipoint\n");

  /* Miscellaneous */
  if (observed.traitLoci && modelType.type != MP)
    logMsg (LOGDEFAULT, LOGFATAL, "Don't specify TraitLoci without Multipoint\n");
  if (modelOptions.polynomialScale && ! modelOptions.polynomial)
    logMsg (LOGDEFAULT, LOGFATAL, "PolynomialScale is incompatible with NonPolynomial\n");
  if ((modelOptions.dkelvinoutfile[0] != '\0') && ! modelOptions.integration)
    logMsg (LOGDEFAULT, LOGFATAL, "NIDetailFile is incompatible with FixedModels\n");
  if ((! modelOptions.integration) &&
      (checkImprintingPenets (&modelRange, modelOptions.imprintingFlag) < 0)) {
    if (modelOptions.imprintingFlag) 
      logMsg (LOGDEFAULT, LOGFATAL, "Imprinting requires Penetrance values for the dD trait genotype\n");
    else 
      logMsg (LOGDEFAULT, LOGFATAL, "Don't specify Penetrance values for the dD trait genotype without Imprinting\n");
  }
  return;
}


void finishConfig ()
{
  int i;
  double integrationLDDPrimeValues[33] =
    {-0.9982431840532, -0.9956010478552, -0.9790222658168, -0.9590960631620, -0.8761473165029,
     -0.8727421201131, -0.7013933644534, -0.6582769255267, -0.6492284325645, -0.5666666666667,
     -0.5000000000000, -0.3808991135940, -0.3582614645881, -0.2517129343453, -0.2077777777778,
     -0.1594544658298, 0.0000000000000, 0.1594544658298, 0.2077777777778, 0.2517129343453,
     0.3582614645881, 0.3808991135940, 0.5000000000000, 0.5666666666667, 0.6492284325645,
     0.6582769255267, 0.7013933644534, 0.8727421201131, 0.8761473165029, 0.9590960631620,
     0.9790222658168, 0.9956010478552, 0.9982431840532};
  
  /* Fill in default values for fields that could have been configured, if they weren't */

  if (modelType.trait == DT) {
    if (modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] == -DBL_MAX) {
      modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = AFFECTION_STATUS_UNKNOWN;
      modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] = AFFECTION_STATUS_UNAFFECTED;
      modelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED] = AFFECTION_STATUS_AFFECTED;
    }
  } else {
    /* QT or CT */
    if (modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] == -DBL_MAX)
      modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = -99.99;
    if (modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] == -DBL_MAX) {
      modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] = -88.88;
      modelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED] = 88.88;
    }
  }

  if (modelOptions.polynomial && ! modelOptions.polynomialScale)
    modelOptions.polynomialScale = 1;
  if ((modelType.type == TP) && (modelOptions.pplfile[0] == '\0'))
    strcpy (modelOptions.pplfile, DEFAULTPPLFILENAME);

  /* Fix up the DPrimes if LD is turned on*/
  if (modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM) {
    /* FIXME: this bit is just for testing */
    strcpy (modelOptions.maxmodelfile, "max.out");
    if (modelOptions.integration == TRUE) {
      /* If integration (that is, dynamic grid) is turned on, then there shouldn't
       * be any DPrimes at all, and we need to insert some. These are MAGIC DPrimes,
       * statically declared at the top of this function.
       */
      for (i = 0; i < 33; i++)
	addDPrime (&modelRange, integrationLDDPrimeValues[i]);

    } else {
      /* If not integration (that is, fixed models) make sure the user
       * didn't omit 0 from the range of DPrimes. Silently add one if needed.
       */
      for (i=0; i<modelRange.ndprime; i++)
        if (fabs(modelRange.dprime[i]) <= ERROR_MARGIN) break;
      if (i == modelRange.ndprime)
        addDPrime (&modelRange, (double) 0.0);
    }
  }

  /* Sync param values for hetrozygous genotypes in non-imprinting runs */
  if ((modelType.trait != DT) && (modelOptions.imprintingFlag != TRUE))
    addConstraint (PARAMC, PEN_dD, 0, 1, EQ, PEN_Dd, 0, 1, FALSE);

  /* Sort the values in the final model. Sorted values better support
   * the application of constraints. */
  sortRange (&modelRange);

  /* Once sorted, removing duplicates is easy. */
  uniqRange (&modelRange);

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
  expandRange (&modelRange);
#if FALSE
  /* Show the partially expanded model. At level 1, following
   * expandRange(), we will have refined the model specification while
   * enforcing all specified constraints that do not involve liability
   * classes. */
  showRange (&modelRange, &modelType, 1);
#endif

  /* Expand the liability classes, but only if necessary and always
   * honoring inter-class constraints. */
  if (modelRange.nlclass > 1)
    expandClass (&modelRange);
  
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


int set_traitLoci (char **toks, int numtoks, void *unused)
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
      addTraitLocus (&modelRange, vlist[va].vun.val);
    } else if (vlist[va].type == VL_RANGE) {
      vb = 0;
      while ((val = vlist[va].vun.range.start + (vb++ * vlist[va].vun.range.incr)) <=
	     vlist[va].vun.range.end)
	addTraitLocus (&modelRange, val);
    } else if (vlist[va].type == VL_RANGE_SYMBEND) {
      modelRange.tlocRangeStart = vlist[va].vun.range.start;
      modelRange.tlocRangeIncr = vlist[va].vun.range.incr;
    } else if ((vlist[va].type == VL_SYMBOL) && (vlist[va].vun.symbol == VL_SYM_MARKER))
      modelRange.tlmark = TRUE;
  }
  free (vlist);
  observed.traitLoci = 1;
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
    addGeneFreq (&modelRange, vals[va]);
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
    addDPrime (&modelRange, vals[va]);
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
  if (strncasecmp (toks[0], THETA_STR, strlen (toks[0])) == 0)
    type = THETA_AVG;
  else if (strncasecmp (toks[0], MALETHETA_STR, strlen (toks[0])) == 0) {
    type = THETA_MALE;
    observed.sexSpecificThetas = 1;
  } else if (strncasecmp (toks[0], FEMALETHETA_STR, strlen (toks[0])) == 0) {
    type = THETA_FEMALE;
    observed.sexSpecificThetas = 1;
  } else
    KLOG (LOGDEFAULT, LOGFATAL, "set_theta called with unexpected directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++)
    addTheta (&modelRange, type, vals[va]);
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
    addAlpha (&modelRange, vals[va]);
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
    addPenetrance (&modelRange, geno-PEN_DD, vals[va]);
  observed.penetrance = 1;
  free (vals);
  return (0);
}


int set_constraint (char **toks, int numtoks, void *unused)
{
  int len, first=2, type=-1;
  int oper, geno1=0, geno2=0, class1, class2, disjunct=0;

  printf ("set_constraint:");
  for (oper = 1; oper < numtoks; oper++)
    printf (" %s", toks[oper]);
  printf (", first %d, numtoks %d\n", first, numtoks);
  
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
	printf ("  SIMPLE: %s %s %s, first %d\n", mp_strs[geno1], op_strs[oper], mp_strs[geno2],
		first);
	
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
	printf ("  CLASSC: %s %d %s %s %d, first %d\n", mp_strs[geno1], class1, op_strs[oper],
		mp_strs[geno2], class2, first);
      } else 
	bail ("illegal argument to directive '%s'\n", toks[0]);

      addConstraint (type, geno1, class1, 0, oper, geno2, class2, 0, disjunct);
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
	printf ("  PARAMC: %s %s %s, first %d\n", mp_strs[geno1], op_strs[oper], mp_strs[geno2],
		first);
	
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
	printf ("  PARAMCLASSC: %s %d %s %s %d, first %d\n", mp_strs[geno1], class1, op_strs[oper],
		mp_strs[geno2], class2, first);
      } else 
	bail ("illegal argument to directive '%s'\n", toks[0]);

      addConstraint (type, geno1, class1, 1, oper, geno2, class2, 1, disjunct);
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
  modelType.type = MP;
  value = (int) strtol (toks[1], &ptr, 10);
  if ((toks[1] == ptr) || (*ptr != '\0'))
    bail ("directive '%s' requires an integer argument\n", toks[0]);
  modelType.numMarkers = value;
  return (0);
}


int set_markerAnalysis (char **toks, int numtoks, void *unused)
{
  if (numtoks < 2)
    bail ("missing argument to directive '%s'\n", toks[0]);
  if (numtoks > 2)
    bail ("extra arguments to directive '%s'\n", toks[0]);
  if (strncasecmp (toks[1], "All", strlen (toks[1])) == 0)
    modelOptions.markerAnalysis = MM;
  else if (strncasecmp (toks[1], "Adjacent", strlen (toks[1])) == 0)
    modelOptions.markerAnalysis = AM;
  else
    bail ("unknown argument to directive '%s'\n", toks[1]);
  return (0);
}


int set_mapFlag (char **toks, int numtoks, void *unused)
{
  if (numtoks > 1)
    bail ("extra arguments to directive '%s'\n", toks[0]);
  modelOptions.mapFlag = SS;
  return (0);
}


int set_disequilibrium (char **toks, int numtoks, void *unused)
{
  if (numtoks > 1)
    bail ("extra arguments to directive '%s'\n", toks[0]);
  modelOptions.equilibrium = LINKAGE_DISEQUILIBRIUM;
  return (0);
}


int set_quantitative (char **toks, int numtoks, void *unused)
{
  int numvals;
  double *vals=NULL;
  
  if (numtoks < 2)
    bail ("missing arguments to directive '%s'\n", toks[0]);

  if (strcasecmp (toks[0], QT_STR) == 0) {
    modelType.trait = QT;
  } else if (strcasecmp (toks[0], QTT_STR) == 0) {
    modelType.trait = CT;
  } else {
    bail ("set_quantitative called with bad directive '%s'\n", toks[0]);
  }

  if (strcasecmp (toks[1], "normal") == 0) {
    if ((numtoks < 3) || ((numvals = expandVals (&toks[2], numtoks-2, &vals, NULL)) != 2))
      bail ("illegal arguments to directive '%s'\n", toks[0]);
    modelType.distrib = QT_FUNCTION_T;
    /* I think this is degrees of freedom; anyway, YH sez: fix it at 30 */
    modelType.constants = realloc (modelType.constants, 1 * sizeof (int));
    modelType.constants[0] = 30;
    modelType.mean = vals[0];
    modelType.sd = vals[1];
    modelRange.npardim = 1;
    free (vals);
  } else if (strcasecmp (toks[1], "chisq") == 0) {
    if (numtoks > 2) 
      bail ("illegal arguments to directive '%s'\n", toks[0]);
    modelType.distrib = QT_FUNCTION_CHI_SQUARE;
    modelType.mean = 0;
    modelType.sd = 1;
    /* A non-empty range for this parameter triggers a loop elsewhere */
    modelRange.npardim = 1;
    addParameter (&modelRange, 0, 1.0);
  } else
    bail ("illegal arguments to directive '%s'\n", toks[0]);
  return (0);
}


int set_qt_mean (char **toks, int numtoks, void *unused)
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
    addPenetrance (&modelRange, geno-PEN_DD, vals[va]);
  observed.mean = 1;
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
    addParameter (&modelRange, 0, vals[va]);
  observed.standarddev = 1;
  free (vals);
  return (0);
}


int set_qt_degfreedom (char **toks, int numtoks, void *unused)
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
    addPenetrance (&modelRange, geno-PEN_DD, vals[va]);
  observed.degfreedom = 1;
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
    addTraitThreshold (&modelRange, vals[va]);
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
      modelType.minOriginal = val;
    } else if (strcasecmp (toks[va], "right") == 0) {
      modelType.maxOriginal = val;
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


int set_affectionStatus (char **toks, int numtoks, void *unused)
{
  int numvals;
  double *vals=NULL;

  if (numtoks < 2)
    bail ("missing arguments to directive '%s'\n", toks[0]);
  if ((numvals = expandVals (&toks[1], numtoks-1, &vals, NULL)) == 1) {
    /* This is legal for QT analyses, just to set the 'undefined' pheno code */
    modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = vals[0];
  } else if (numvals == 3) {
    /* Everything else (DT, CT) requires three pheno codes */
    modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = vals[0];
    modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] = vals[1];
    modelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED] = vals[2];
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
  if ((len = strlen (toks[1])) > KMAXFILENAMELEN - 2)
    bail ("argument to directive '%s' is too long\n", toks[0]);
  strcpy (modelOptions.resultsprefix, toks[1]);
  if (modelOptions.resultsprefix[len-1] != '/')
    strcat (modelOptions.resultsprefix, "/");
  return (0);
}


int set_logLevel (char **toks, int numtoks, void *filename)
{
  int logType=0, logLevel=0;

  if (numtoks < 3)
    bail ("missing argument(s) to directive '%s'\n", toks[0]);
  if (numtoks > 3)
    bail ("extra arguments to directive '%s'\n", toks[0]);

  if (!strcasecmp (toks[1], "pedfile"))
    logType = LOGPEDFILE;
  else if (!strcasecmp (toks[1], "inputfile"))
    logType = LOGINPUTFILE;
  else if (!strcasecmp (toks[1], "genoelim"))
    logType = LOGGENOELIM;
  else if (!strcasecmp (toks[1], "parentalpair"))
    logType = LOGPARENTALPAIR;
  else if (!strcasecmp (toks[1], "peelgraph"))
    logType = LOGPEELGRAPH;
  else if (!strcasecmp (toks[1], "likelihood"))
    logType = LOGLIKELIHOOD;
  else if (!strcasecmp (toks[1], "setrecoding"))
    logType = LOGSETRECODING;
  else if (!strcasecmp (toks[1], "memory"))
    logType = LOGMEMORY;
  else if (!strcasecmp (toks[1], "integration"))
    logType = LOGINTEGRATION;
  else if (!strcasecmp (toks[1], "default"))
    logType = LOGDEFAULT;
  else
    bail ("unknown log facility '%s'", toks[1]);
  
  if (!strcasecmp (toks[2], "fatal"))
    logLevel = LOGFATAL;
  else if (!strcasecmp (toks[2], "error"))
    logLevel = LOGERROR;
  else if (!strcasecmp (toks[2], "warning"))
    logLevel = LOGWARNING;
  else if (!strcasecmp (toks[2], "advise"))
    logLevel = LOGADVISE;
  else if (!strcasecmp (toks[2], "debug"))
    logLevel = LOGDEBUG;
  else
    bail ("unknown log severity '%s'", toks[2]);

  logSet (logType, logLevel);
  return (0);
}


int expandVals (char **toks, int numtoks, double **vals_h, st_valuelist **vlist_h)
{
  int numvals=0, listsize=10, tokidx=0, va;
  char *ca, *cb;
  double start, end=-1, incr, val, *vals=NULL;
  st_valuelist *vlist=NULL;

  /* Sanity check */
  if (numtoks == 0)
    return (0);

  /* Exactly one return value pointer must be non-NULL */
  if (((vlist_h == NULL) && (vals_h == NULL)) || ((vlist_h != NULL) && (vals_h != NULL)))
    return (-1);
  
  if (vlist_h != NULL) {
    if ((vlist = malloc (sizeof (st_valuelist) * listsize)) == NULL)
      logMsg (LOGDEFAULT, LOGFATAL, "malloc failed\n");
  } else {
    if ((vals = malloc (sizeof (double) * listsize)) == NULL)
      logMsg (LOGDEFAULT, LOGFATAL, "malloc failed\n");
  }

  // printf ("starting\n");
  while (1) {
    if (numvals >= listsize) {
      if (((vlist != NULL) &&
	   ((vlist = realloc (vlist, sizeof (st_valuelist) * (listsize += 10))) == NULL)) || 
	  ((vals = realloc (vals, sizeof (double) * (listsize += 10))) == NULL))
	logMsg (LOGDEFAULT, LOGFATAL, "malloc failed\n");
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
	while ((val = start + (va++ * incr)) <= end) {
	  if ((numvals >= listsize) && 
	      ((vals = realloc (vals, sizeof (double) * (listsize += 10))) == NULL))
	    logMsg (LOGDEFAULT, LOGFATAL, "realloc failed\n");
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
  
  while ((buffptr == NULL) || (strlen (buffptr) == 0)) {
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
  return (numtoks);
}


/* Splits a character string into a list of tokens. The line should already
 * have been passed through permuteLine(), so extra whitespace has already been
 * trimmed.
 */
int tokenizeLine (char *line, char ***tokgroup_h, int *tokgroupsize)
{
  int numtoks=0;
  char **tokgroup, *ca, *cb;
  
  if ((tokgroup = *tokgroup_h) == NULL) {
    *tokgroupsize = 10;
    if ((tokgroup = malloc (sizeof (char *) * *tokgroupsize)) == NULL) {
      logMsg (LOGDEFAULT, LOGFATAL, "malloc failed\n");
      exit (-1);
    }
  }
  
  ca = strtok_r (line, " ", &cb);
  while (ca != NULL) {
    if (numtoks + 1 >= *tokgroupsize) {
      *tokgroupsize += 10;
      if ((tokgroup = realloc (tokgroup, sizeof (char *) * *tokgroupsize)) == NULL)
	logMsg (LOGDEFAULT, LOGFATAL, "malloc failed\n");
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

/* Dumps a subset of the fields in modelOptions.
 */
void dumpModelOptions (ModelOptions *mo)
{
  printf ("%18s : %s\n", "markerfile", mo->markerfile);
  printf ("%18s : %s\n", "mapfile", mo->mapfile);
  printf ("%18s : %s\n", "pedfile", mo->pedfile);
  printf ("%18s : %s\n", "datafile", mo->datafile);
  printf ("%18s : %s\n", "avghetfile", mo->avghetfile);
  printf ("%18s : %s\n", "pplfile", mo->pplfile);
  printf ("%18s : %s\n", "condFile", mo->condFile);
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
    logMsg (LOGDEFAULT, LOGFATAL, newfmt, conffilename, lineno, arg);
  } else {
    strcpy (newfmt, "on command line: ");
    strcat (newfmt, fmt);
    logMsg (LOGDEFAULT, LOGFATAL, newfmt, arg);
  }
}
