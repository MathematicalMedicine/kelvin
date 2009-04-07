#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
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


/* Globals */
char buff[BUFFSIZE] = "";   /* These two are global to provide context to getNextTokgroup */
char *buffptr = NULL;

/* These shouldn't be here unless we're using the local main() */
ModelOptions modelOptions;
ModelRange modelRange;
ModelType modelType;

/* prototypes for non-public routines */
void initializeDefaults ();
int expand_vals (char **toks, int numtoks, double **vals_h, st_valuelist **vlist_h);
int lookupDispatch (char *key, st_dispatch *table);
int compareDispatch (const void *a, const void *b);
int getNextTokgroup (FILE *fp, char ***tokgroup_h, int *tokgroupsize);
int tokenizeLine (char *line, char ***tokgroup_h, int *tokgroupsize);
int single_digit (char *str);
void dumpModelOptions (ModelOptions *mo);

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

				{"NonPolynomial", clear_flag, &modelOptions.integration},
				{"Imprinting", set_flag, &modelOptions.imprintingFlag},
				{"SexLinked", set_flag, &modelOptions.sexLinked},
				{"FixedModels", clear_flag, &modelOptions.integration},
				{"DryRun", set_flag, &modelOptions.dryRun},
				{"ExtraMODs", clear_flag, &modelOptions.extraMODs},

				{"PolynomialScale", set_int, &modelOptions.polynomialScale},
				{"LiabilityClasses", set_int, &modelRange.nlclass},
				{"DiseaseAlleles", set_int, &modelRange.nalleles},

				{"TraitLoci", set_traitLoci, NULL},
				{"DiseaseGeneFrequency", set_geneFreq, NULL},
				{"DPrime", set_dprime, NULL},
				{"Theta", set_theta, NULL},
				{"MaleTheta", set_theta, NULL},
				{"FemaleTheta", set_theta, NULL},
				{"Alpha", set_alpha, NULL},
				{"Penetrance", set_penetrance, NULL},
				{"Constraint", set_constraint, NULL},
				{"Multipoint", set_multipoint, NULL},
				{"MarkerToMarker", set_markerAnalysis, NULL},
				{"SexSpecific", set_mapFlag, NULL},
				{"LD", set_disequilibrium, NULL},
				{"PhenoCodes", set_affectionStatus, NULL},
				{"SurfacesPath", set_resultsprefix, NULL},
				/*{"condfile", set_condrun, &modelOptions.condFile},*/
				{"Log", set_logLevel, NULL}
};


main (int argc, char *argv[])
{
  logInit ();

  if (argc == 1)
    logMsg (LOGDEFAULT, LOGFATAL, "%s: no configuration file specified\n", argv[0]);
  initializeDefaults ();
  my_readConfigFile (argv[1]);
  if (argc > 2) {
    my_parseCommandLine (argc-2, &argv[2]);
  }
  dumpModelOptions (&modelOptions);
}


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
  strcpy (modelOptions.pplfile, DEFAULTPPLFILENAME);
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
  modelOptions.polynomialScale = 1;
  modelOptions.extraMODs = FALSE;
  modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = AFFECTION_STATUS_UNKNOWN;
  modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] = AFFECTION_STATUS_UNAFFECTED;
  modelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED] = AFFECTION_STATUS_AFFECTED;

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
  /* set default for QT */
  modelType.minOriginal = -999999999.00;
  modelType.maxOriginal = 999999999.00;
  modelType.minThreshold = -999999999.00;
  modelType.maxThreshold = 999999999.00;

  return;
}


void my_readConfigFile (char *config)
{
  int numtoks, tokgroupsize=0, va;
  char **toks = NULL;
  FILE *conffp;

  if ((conffp = fopen (config, "r")) == NULL)
    logMsg (LOGDEFAULT, LOGFATAL, "open '%s' failed, %s\n", config, strerror (errno));

  /* Sort the dispatch table so the binary search works */
  qsort (dispatchTable, sizeof (dispatchTable) / sizeof (st_dispatch), sizeof (st_dispatch),
	 compareDispatch);

  for (va = 0; va < sizeof (dispatchTable) / sizeof (st_dispatch); va++)
    printf ("idx %2d: %s\n", va, dispatchTable[va].key);
  
  while (numtoks = getNextTokgroup (conffp, &toks, &tokgroupsize)) {
    for (va = 0; va < numtoks; va++) {
      printf ("tok %d: %s\n", va, toks[va]);
    }
    if ((va = lookupDispatch (toks[0], dispatchTable)) >= 0) {
      printf ("directive '%s' matches at index %d\n", toks[0], va);
      (*dispatchTable[va].parse) (toks, numtoks, dispatchTable[va].hint);
    } else {
      printf ("directive '%s' is %s\n", toks[0], (va == -1) ? "unknown" : "not unique");
    }
    printf ("\n");
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
	printf ("directive '%s' matches at index %d\n", toks[0], va);
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
    printf ("directive '%s' matches at index %d\n", toks[0], va);
    (*dispatchTable[va].parse) (toks, numtoks, dispatchTable[va].hint);
  } else
    logMsg (LOGDEFAULT, LOGFATAL, "directive '%s' on command line is %s\n", toks[0],
	    (va == -1) ? "unknown" : "not unique");
  
  if (tokgroupsize > 0)
    free (toks);
  return;
}


int set_optionfile (char **toks, int numtoks, void *filename)
{
  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing filename argument to directive '%s'\n", toks[0]);
  if (numtoks > 2)
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  strcpy ((char *) filename, toks[1]);
  return (0);
}


int set_flag (char **toks, int numtoks, void *flag)
{
  if (numtoks > 1)
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  *((int *) flag) = TRUE;
  return (0);
}


int clear_flag (char **toks, int numtoks, void *flag)
{
  if (numtoks > 1)
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  *((int *) flag) = FALSE;
  return (0);
}


int set_int (char **toks, int numtoks, void *field)
{
  int value;
  char *ptr = NULL;

  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing integer argument to directive '%s'\n", toks[0]);
  if (numtoks > 2)
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  value = (int) strtol (toks[1], &ptr, 10);
  if (toks[0] == ptr)
    logMsg (LOGDEFAULT, LOGFATAL, "argument '%s' to directive '%s' is not an integer\n",
	    toks[1], toks[0]);
  *((int *) field) = value;
  return (0);
}


int set_markerAnalysis (char **toks, int numtoks, void *unused)
{
  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing argument to directive '%s'\n", toks[0]);
  if (numtoks > 2)
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  if (strncasecmp (toks[1], "All", strlen (toks[1])) == 0)
    modelOptions.markerAnalysis = MM;
  else if (strncasecmp (toks[1], "Adjacent", strlen (toks[1])) == 0)
    modelOptions.markerAnalysis = AM;
  else
    logMsg (LOGDEFAULT, LOGFATAL, "set_markerAnalysis called with bad argument '%s'\n", toks[1]);
  return (0);
}


int set_mapFlag (char **toks, int numtoks, void *unused)
{
  if (numtoks > 1)
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  modelOptions.mapFlag = SS;
  return (0);
}


int set_multipoint (char **toks, int numtoks, void *unused)
{
  int value;
  char *ptr = NULL;

  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing integer argument to directive '%s'\n", toks[0]);
  if (numtoks > 2)
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  modelType.type = MP;
  value = (int) strtol (toks[1], &ptr, 10);
  if (toks[0] == ptr)
    logMsg (LOGDEFAULT, LOGFATAL, "argument '%s' to directive '%s' is not an integer\n",
	    toks[1], toks[0]);
  modelType.numMarkers = value;
  return (0);
}


int set_disequilibrium (char **toks, int numtoks, void *unused)
{
  if (numtoks > 1)
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  modelOptions.equilibrium = LINKAGE_DISEQUILIBRIUM;
  return (0);
}


int set_affectionStatus (char **toks, int numtoks, void *unused)
{
  int numvals;
  double *vals=NULL;

  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing list argument to directive '%s'\n", toks[0]);
  if ((numvals = expand_vals (&toks[1], numtoks-1, &vals, NULL)) != 3)
    logMsg (LOGDEFAULT, LOGFATAL, "bad arguments to directive '%s'\n", toks[0]);
  modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = vals[0];
  modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] = vals[1];
  modelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED] = vals[2];
  free (vals);

  return (0);
}


int set_traitLoci (char **toks, int numtoks, void *unused)
{
  int numvals, va=0;
  st_valuelist *vlist;
  double val;

  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing argument to directive '%s'\n", toks[0]);
  if ((numvals = expand_vals (&toks[1], numtoks-1, NULL, &vlist)) <= 0)
    logMsg (LOGDEFAULT, LOGFATAL, "illegal argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++) {
    if (vlist[va].type == VL_VALUE) {
      addTraitLocus (&modelRange, vlist[va].vun.val);
    } else if (vlist[va].type == VL_RANGE) {
      while ((val = vlist[va].vun.range.start + (va++ * vlist[va].vun.range.incr)) <=
	     vlist[va].vun.range.end)
	addTraitLocus (&modelRange, val);
    } else if (vlist[va].type == VL_RANGE_SYMBEND) {
      modelRange.tlocRangeStart = vlist[va].vun.range.start;
      modelRange.tlocRangeIncr = vlist[va].vun.range.incr;
    } else if ((vlist[va].type == VL_SYMBOL) && (vlist[va].vun.symbol == VL_SYM_MARKER))
      modelRange.tlmark = TRUE;
  }
  free (vlist);

  return (0);
}


int set_geneFreq (char **toks, int numtoks, void *unused)
{
  int numvals, va=0;
  double *vals;

  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing list argument to directive '%s'\n", toks[0]);
  if ((numvals = expand_vals (&toks[1], numtoks-1, &vals, NULL)) <= 0)
    logMsg (LOGDEFAULT, LOGFATAL, "illegal argument to directive '%s'\n", toks[0]);
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
    logMsg (LOGDEFAULT, LOGFATAL, "missing list argument to directive '%s'\n", toks[0]);
  if ((numvals = expand_vals (&toks[1], numtoks-1, &vals, NULL)) <= 0)
    logMsg (LOGDEFAULT, LOGFATAL, "illegal argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++)
    addDPrime (&modelRange, vals[va]);
  free (vals);

  return (0);
}


int set_theta (char **toks, int numtoks, void *unused)
{
  int numvals, va=0, type;
  double *vals;

  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing list argument to directive '%s'\n", toks[0]);
  if ((numvals = expand_vals (&toks[1], numtoks-1, &vals, NULL)) <= 0)
    logMsg (LOGDEFAULT, LOGFATAL, "illegal argument to directive '%s'\n", toks[0]);
  /* This is evil; I shouldn't presume to know what the legal values for toks[0] are */
  if (strcasecmp (toks[0], "Theta") == 0)
    type = THETA_AVG;
  else if (strcasecmp (toks[0], "MaleTheta") == 0)
    type = THETA_MALE;
  else if (strcasecmp (toks[0], "MaleTheta") == 0)
    type = THETA_FEMALE;
  else
    logMsg (LOGDEFAULT, LOGFATAL, "set_theta called with unexpected directive '%s'\n", toks[0]);
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
    logMsg (LOGDEFAULT, LOGFATAL, "missing list argument to directive '%s'\n", toks[0]);
  if ((numvals = expand_vals (&toks[1], numtoks-1, &vals, NULL)) <= 0)
    logMsg (LOGDEFAULT, LOGFATAL, "illegal argument to directive '%s'\n", toks[0]);
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
    logMsg (LOGDEFAULT, LOGFATAL, "missing argument to directive '%s'\n", toks[0]);
  if ((geno = lookup_modelparam (toks[1])) == -1)
    logMsg (LOGDEFAULT, LOGFATAL, "illegal argument to directive '%s'\n", toks[0]);
  if ((numvals = expand_vals (&toks[2], numtoks-2, &vals, NULL)) <= 0)
    logMsg (LOGDEFAULT, LOGFATAL, "illegal argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++)
    addPenetrance (&modelRange, geno-PEN_DD, vals[va]);
  free (vals);

  return (0);
}


int set_constraint (char **toks, int numtoks, void *unused)
{
  int first=1, type=-1, oper, geno1, geno2, class1, class2, param1, param2, disjunct=0;
  char *ca, *cb;

  printf ("set_constraint:");
  for (oper = 1; oper < numtoks; oper++)
    printf (" %s", toks[oper]);
  printf (", first %d, numtoks %d\n", first, numtoks);

  while (1) {
    class1 = class2 = param1 = param2 = 0;
    if ((numtoks >= first + 3) &&
	((geno1 = lookup_modelparam (toks[first])) != -1) &&
	((oper = lookup_comparator (toks[first+1])) != -1) &&
	((geno2 = lookup_modelparam (toks[first+2])) != -1)) {
      if ((type != -1) && (type != SIMPLE))
	logMsg (LOGDEFAULT, LOGFATAL, "illegal combination of constraints\n");
      type = SIMPLE;
      first += 3;
      printf ("  SIMPLE: %s %s %s, first %d\n", mp_strs[geno1], op_strs[oper], mp_strs[geno2],
	      first);
      
    } else if ((numtoks >= first + 5) &&
	       ((geno1 = lookup_modelparam (toks[first])) != -1) &&
	       ((class1 = single_digit (toks[first+1])) > 0) &&
	       ((oper = lookup_comparator (toks[first+2])) != -1) &&
	       ((geno2 = lookup_modelparam (toks[first+3])) != -1) &&
	       ((class2 = single_digit (toks[first+4])) > 0)) {

      if ((type != -1) && (type != CLASSC))
	logMsg (LOGDEFAULT, LOGFATAL, "illegal combination of constraints\n");
      type = CLASSC;
      first += 5;
      printf ("  CLASSC: %s %d %s %s %d, first %d\n", mp_strs[geno1], class1, op_strs[oper],
	      mp_strs[geno2], class2, first);

    } else if ((numtoks >= first + 5) &&
	       (toks[first][0] == 'P') && ((param1 = single_digit (toks[first]+1)) > 0) &&
	       ((geno1 = lookup_modelparam (toks[first+1])) != -1) &&
	       ((oper = lookup_comparator (toks[first+2])) != -1) &&
	       (toks[first+3][0] == 'P') && ((param1 = single_digit (toks[first+3]+1)) > 0) &&
	       ((geno2 = lookup_modelparam (toks[first+4])) != -1)) {

      if ((type != -1) && (type != PARAMC))
	logMsg (LOGDEFAULT, LOGFATAL, "illegal combination of constraints\n");
      type = PARAMC;
      first += 5;
      printf ("  PARAMC: P%d %s %s P%d %s, first %d\n", param1, mp_strs[geno1], op_strs[oper],
	      param2, mp_strs[geno2], first);

    } else if ((numtoks >= first + 7) &&
	       (toks[first][0] == 'P') && ((param1 = single_digit (toks[first]+1)) > 0) &&
	       ((geno1 = lookup_modelparam (toks[first+1])) != -1) &&
	       ((class1 = single_digit (toks[first+2])) > 0) &&
	       ((oper = lookup_comparator (toks[first+3])) != -1) &&
	       (toks[first+4][0] == 'P') && ((param1 = single_digit (toks[first+4]+1)) > 0) &&
	       ((geno2 = lookup_modelparam (toks[first+5])) != -1) &&
	       ((class2 = single_digit (toks[first+6])) > 0)) {

      if ((type != -1) && (type != PARAMCLASSC))
	logMsg (LOGDEFAULT, LOGFATAL, "illegal combination of constraints\n");
      type = PARAMCLASSC;
      first += 7;
      printf ("  PARAMCLASSC: P%d %s %d %s P%d %s %d, first %d\n", param1, mp_strs[geno1],
	      class1, op_strs[oper], param2, mp_strs[geno2], class2, first);

    } else
      logMsg (LOGDEFAULT, LOGFATAL, "illegal argument to directive '%s'\n", toks[0]);

    addConstraint (type, geno1, class1, param1, oper, geno2, class2, param2, disjunct);
    disjunct = 1;
    if (numtoks <= first)
      return (0);
    if (strcmp (toks[first++], ",") != 0)
      logMsg (LOGDEFAULT, LOGFATAL, "illegal argument to directive '%s'\n", toks[0]);
  }
}


int set_resultsprefix (char **toks, int numtoks, void *unused)
{
  int len;

  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing argument to directive '%s'\n", toks[0]);
  if (numtoks > 2)
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  if ((len = strlen (toks[1])) > KMAXFILENAMELEN - 2)
    logMsg (LOGDEFAULT, LOGFATAL, "argument to directive '%s' is too long\n", toks[0]);
  strcpy (modelOptions.resultsprefix, toks[1]);
  if (modelOptions.resultsprefix[len-1] != '/')
    strcat (modelOptions.resultsprefix, "/");
  return (0);
}


int set_logLevel (char **toks, int numtoks, void *filename)
{
  int logType, logLevel;

  if (numtoks < 3)
    logMsg (LOGDEFAULT, LOGFATAL, "missing argument(s) to directive '%s'\n", toks[0]);
  if (numtoks > 3)
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);

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
    logMsg (LOGDEFAULT, LOGFATAL, "unknown log facility '%s'", toks[2]);
  
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
    logMsg (LOGDEFAULT, LOGFATAL, "unknown log severity '%s'", toks[2]);

  logSet (logType, logLevel);
  return (0);
}


int expand_vals (char **toks, int numtoks, double **vals_h, st_valuelist **vlist_h)
{
  int numvals=0, listsize=10, tokidx=0, va;
  char *ca, *cb, *cc, sub[BUFFSIZE];
  double start, end, incr, val, *vals=NULL;
  st_valuelist *vlist=NULL;

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
	  vlist[numvals++].vun.val = strtod (sub, &cb);
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
    if ((res = strncasecmp (key, table[mid].key, keylen)) == 0) {
      if ((strlen (table[mid].key) != keylen) &&
	  (((mid > 0) && (strncasecmp (key, table[mid-1].key, keylen) == 0)) || 
	   ((mid < tablen - 1) && (strncasecmp (key, table[mid+1].key, keylen) == 0))))
	return (-2);
      return (mid);
    } else if (res < 0) {
      hi = mid - 1;
    } else {
      lo = mid + 1;
    }
    mid = (hi + lo) / 2;
  }
  return (-1);
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
    if (fgets (buff, BUFFSIZE, fp) == 0)
      return (0);
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


int single_digit (char *str)
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
