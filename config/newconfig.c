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


/* Legal values for the 'type' field of st_valuelist. */
#define VL_VALUE          0
#define VL_RANGE          1
#define VL_RANGE_SYMBEND  2
#define VL_SYMBOL         4

/* Legal values for the vun.symbol field of st_valuelist, if type is VL_SYMBOL */
#define VL_SYM_MARKER      0


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


/* Globals */
int exactTokenCount = 1;    /* controls whether token parsers will tolerate extra tokens */
char buff[BUFFSIZE] = "";   /* These two are global to provide context to getNextTokgroup */
char *buffptr = NULL;

/* These shouldn't be here unless we're using the local main() */
ModelOptions modelOptions;
ModelRange modelRange;
ModelType modelType;

/* prototypes for non-public routines */
void initializeDefaults ();
int expand_vals (char *str, double **vals, st_valuelist **vlist);
int lookupDispatch (char *key, st_dispatch *table);
int compareDispatch (const void *a, const void *b);
int getNextTokgroup (FILE *fp, char ***tokgroup_h, int *tokgroupsize);
void dumpModelOptions (ModelOptions *mo);

/* functions for use in the dispatch table */
int set_optionfile (char **toks, int numtoks, void *filename);
int set_flag (char **toks, int numtoks, void *flag);
int clear_flag (char **toks, int numtoks, void *flag);
int set_int (char **toks, int numtoks, void *field);
int set_markerAnalysis (char **toks, int numtoks, void *unused);
int set_mapFlag (char **toks, int numtoks, void *unused);
int set_disequilibrium (char **toks, int numtoks, void *unused);
int set_affectionStatus (char **toks, int numtoks, void *unused);
int set_traitLoci (char **toks, int numtoks, void *unused);
int set_geneFreq (char **toks, int numtoks, void *unused);
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

				{"MarkerToMarker", set_markerAnalysis, NULL},
				{"SexSpecific", set_mapFlag, NULL},
				{"LD", set_disequilibrium, NULL},
				{"PhenoCodes", set_affectionStatus, NULL},
				{"SurfacesPath", set_resultsprefix, NULL},
				{"TraitLoci", set_traitLoci, NULL},
				{"DiseaseGeneFrequency", set_geneFreq, NULL},
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
    exactTokenCount = 0;
    my_parseCommandLine (argc-2, &argv[2]);
  }
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
  dumpModelOptions (&modelOptions);
}


void my_parseCommandLine (int argc, char *argv[])
{
  int curidx=0, used, va;

  while (curidx < argc) {
    /* in case we don't get as far as lookupDispatch in the following conditional */
    va == -1;
    /* FIXME: Do we want to force leading dashes? */
    if ((strncmp (argv[curidx], "--", 2) != 0) || 
	(memmove (argv[curidx], argv[curidx]+2, strlen (argv[curidx]) - 2) == NULL) ||
	((va = lookupDispatch (argv[curidx], dispatchTable)) < 0)) {
      if (va == -1)
	logMsg (LOGDEFAULT, LOGFATAL, 
		"found '%s' on command line when expecting directive\n", argv[curidx]);
      else
	logMsg (LOGDEFAULT, LOGFATAL, 
		"directive '%s' on command line in not unique\n", argv[curidx]);
    }
    printf ("directive '%s' matches at index %d\n", argv[curidx], va);
    curidx += (*dispatchTable[va].parse) (&argv[curidx], argc - curidx, dispatchTable[va].hint);
  }
  return;
}


int set_optionfile (char **toks, int numtoks, void *filename)
{
  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing filename argument to directive '%s'\n", toks[0]);
  if ((numtoks > 2) && (exactTokenCount))
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  strcpy ((char *) filename, toks[1]);
  return (2);
}


int set_flag (char **toks, int numtoks, void *flag)
{
  if ((numtoks > 1) && (exactTokenCount))
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  *((int *) flag) = TRUE;
  return (1);
}


int clear_flag (char **toks, int numtoks, void *flag)
{
  if ((numtoks > 1) && (exactTokenCount))
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  *((int *) flag) = FALSE;
  return (1);
}


int set_int (char **toks, int numtoks, void *field)
{
  int value;
  char *ptr = NULL;

  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing integer argument to directive '%s'\n", toks[0]);
  if ((numtoks > 2) && (exactTokenCount))
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  value = (int) strtol (toks[1], &ptr, 10);
  if (toks[0] == ptr)
    logMsg (LOGDEFAULT, LOGFATAL, "argument '%s' to directive '%s' is not an integer\n",
	    toks[1], toks[0]);
  *((int *) field) = value;
  return (2);
}


int set_markerAnalysis (char **toks, int numtoks, void *unused)
{
  if ((numtoks > 2) && (exactTokenCount))
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  if (strncasecmp (toks[1], "All", strlen (toks[1])) == 0)
    modelOptions.markerAnalysis = MM;
  else if (strncasecmp (toks[1], "Adjacent", strlen (toks[1])) == 0)
    modelOptions.markerAnalysis = AM;
  else
    logMsg (LOGDEFAULT, LOGFATAL, "set_markerAnalysis called with bad token '%s'\n", toks[1]);
  return (1);
}


int set_mapFlag (char **toks, int numtoks, void *unused)
{
  if ((numtoks > 1) && (exactTokenCount))
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  modelOptions.mapFlag = SS;
  return (1);
}


int set_disequilibrium (char **toks, int numtoks, void *unused)
{
  if ((numtoks > 1) && (exactTokenCount))
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  modelOptions.equilibrium = LINKAGE_DISEQUILIBRIUM;
  return (1);
}


int set_affectionStatus (char **toks, int numtoks, void *unused)
{
  int numvals;
  double *vals=NULL;

  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing list argument to directive '%s'\n", toks[0]);
  if ((numtoks > 2) && (exactTokenCount))
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  if ((numvals = expand_vals (toks[1], &vals, NULL)) != 3)
    logMsg (LOGDEFAULT, LOGFATAL, "directive '%s' requires three values\n", toks[0]);
  modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN] = vals[0];
  modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED] = vals[1];
  modelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED] = vals[2];
  free (vals);

  return (2);
}


int set_traitLoci (char **toks, int numtoks, void *unused)
{
  int numvals, va=0;
  st_valuelist *vlist;
  double val;

  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing list argument to directive '%s'\n", toks[0]);
  if ((numtoks > 2) && (exactTokenCount))
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  if ((numvals = expand_vals (toks[1], NULL, &vlist)) <= 0)
    logMsg (LOGDEFAULT, LOGFATAL, "illegal or missing argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++) {
    if (vlist[va].type == VL_VALUE) {
      addTratiLoci (&modelRange, vlist[va].vun.val);
    } else if (vlist[va].type == VL_RANGE) {
      while ((val = vlist[va].vun.range.start + (va++ * vlist[va].vun.range.incr)) <=
	     vlist[va].vun.range.end)
	addTratiLoci (&modelRange, val);
    } else if (vlist[va].type == VL_RANGE_SYMBEND) {
      modelRange.tlocRangeStart = vlist[va].vun.range.start;
      modelRange.tlocRangeIncr = vlist[va].vun.range.incr;
    } else if ((vlist[va].type == VL_SYMBOL) && (vlist[va].vun.symbol == VL_SYM_MARKER))
      modelRange.tlmark = TRUE;
  }
  free (vlist);

  return (2);
}


int set_geneFreq (char **toks, int numtoks, void *unused)
{
  int numvals, va=0;
  double *vals;

  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing list argument to directive '%s'\n", toks[0]);
  if ((numtoks > 2) && (exactTokenCount))
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  if ((numvals = expand_vals (toks[1], NULL, &vlist)) <= 0)
    logMsg (LOGDEFAULT, LOGFATAL, "illegal or missing argument to directive '%s'\n", toks[0]);
  for (va = 0; va < numvals; va++)
    addGeneFreq (&modelRange, val);
  free (vals);

  return (2);
}


int set_resultsprefix (char **toks, int numtoks, void *unused)
{
  int len;

  if (numtoks < 2)
    logMsg (LOGDEFAULT, LOGFATAL, "missing argument to directive '%s'\n", toks[0]);
  if ((numtoks > 2) && (exactTokenCount))
    logMsg (LOGDEFAULT, LOGFATAL, "extra arguments to directive '%s'\n", toks[0]);
  if ((len = strlen (toks[1])) > KMAXFILENAMELEN - 2)
    logMsg (LOGDEFAULT, LOGFATAL, "argument to directive '%s' is too long\n", toks[0]);
  strcpy (modelOptions.resultsprefix, toks[1]);
  if (modelOptions.resultsprefix[len-1] != '/')
    strcat (modelOptions.resultsprefix, "/");
  return (2);
}


int set_logLevel (char **toks, int numtoks, void *filename)
{
  int logType, logLevel;

  if (numtoks < 3)
    logMsg (LOGDEFAULT, LOGFATAL, "missing argument(s) to directive '%s'\n", toks[0]);
  if ((numtoks > 3) && (exactTokenCount))
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
  return (3);
}


int expand_vals (char *str, double **vals_h, st_valuelist **vlist_h)
{
  int numvals=0, listsize=10, va;
  char *ca, *cb, *cc, sub[BUFFSIZE];
  double start, end, incr, val, *vals=NULL;
  st_valuelist *vlist=NULL;
  
  if (((vlist_h == NULL) && (vals_h == NULL)) || ((vlist_h != NULL) && (vals_h != NULL)))
    return (-1);
  
  if (vlist_h != NULL) {
    if ((vlist = malloc (sizeof (st_valuelist) * listsize)) == NULL)
      logMsg (LOGDEFAULT, LOGFATAL, "malloc failed\n");
  } else if (vals_h != NULL) {
    if ((vals = malloc (sizeof (double) * listsize)) == NULL)
      logMsg (LOGDEFAULT, LOGFATAL, "malloc failed\n");
  } else
    return (-1);
  
  ca = str;
  // printf ("starting\n");
  while (1) {
    if (*ca == '\0') {
      // printf ("ca is empty\n");
      if (vlist_h != NULL)
	*vlist_h = vlist;
      else 
	*vals_h = vals;
      return (numvals);
    }

    // printf ("ca is '%s'\n", ca);
    /* Copy the next comma-separated substring, or the entire remaining string if no comma */
    memset (sub, 0, BUFFSIZE);
    if ((cb = index (ca, ',')) != NULL) {
      strncpy (sub, ca, cb - ca);
      ca = cb + 1;
    } else {
      strcpy (sub, ca);
      ca = sub + strlen (sub);
    }
    // printf (" sub is '%s'\n", sub);
    
    if (numvals >= listsize) {
      if (((vlist != NULL) &&
	   ((vlist = realloc (vlist, sizeof (st_valuelist) * (listsize += 10))) == NULL)) || 
	  ((vals = realloc (vals, sizeof (double) * (listsize += 10))) == NULL))
	logMsg (LOGDEFAULT, LOGFATAL, "malloc failed\n");
    }
    
    /* Skip the first character of sub to avoid the leading '-' of a negative number */
    if ((cc = index (sub+1, '-')) == NULL) {
      /* This should be a single value or symbol */
      // printf ("  single value: %f\n", vals[numvals-1]);
      
      if (vlist != NULL) {
	/* returning a list of st_valuelist */
	if (strcasecmp (sub, "Marker") == 0) {
	  vlist[numvals].type = VL_SYMBOL;
	  vlist[numvals++].vun.symbol = VL_SYM_MARKER;
	} else {
	  vlist[numvals].type = VL_VALUE;
	  vlist[numvals++].vun.val = strtod (sub, &cc);
	  if ((cc == sub) || (*cc != '\0'))
	    break;
	}
      } else {
	/* returning a list of doubles */
	vals[numvals++] = strtod (sub, &cc);
	if ((cc == sub) || (*cc != '\0'))
	  break;
      }
      
    } else {
      /* This should be a range of values in the form 'start-end:increment' */
      *cc = '\0';
      // printf ("  range format: first substring '%s', remainder, '%s'\n", sub, cc+1);
      start = strtod (sub, &cc);
      if ((cc == sub) || (*cc != '\0'))
	break;
      // printf ("  start value: %f\n", start);
      cb = cc + 1;
      if ((*cb == '\0') || ((cc = index (cb, ':')) == NULL))
	break;
      // printf ("                next substring '%s', remainder, '%s'\n", cb, cc+1);
      *cc = '\0';
      if (vlist != NULL) {
	if (strcasecmp (cb, "End") == 0) {
	  vlist[numvals].type = VL_RANGE_SYMBEND;
	} else {
	  vlist[numvals].type = VL_RANGE;
	  end = strtod (cb, &cc);
	  if ((cc == cb) || (*cc != '\0'))
	    break;
	}
      } else {
	end = strtod (cb, &cc);
	if ((cc == cb) || (*cc != '\0'))
	  break;
	// printf ("  end value: %f\n", end);
      }
      cb = cc + 1;
      if (*cb == '\0')
	break;
      // printf ("                last substring '%s'\n", cb);
      incr = strtod (cb, &cc);
      if (*cc != '\0')
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


int getNextTokgroup (FILE *fp, char ***tokgroup_h, int *tokgroupsize)
{
  int numtoks=0, size;
  char **tokgroup, *semicolon, *ca, *cb;
  
  if ((tokgroup = *tokgroup_h) == NULL) {
    *tokgroupsize = 10;
    if ((tokgroup = malloc (sizeof (char *) * *tokgroupsize)) == NULL) {
      logMsg (LOGDEFAULT, LOGFATAL, "malloc failed\n");
      exit (-1);
    }
  }

  while ((buffptr == NULL) || (strlen (buffptr) == 0)) {
    if (fgets (buff, BUFFSIZE, fp) == 0) {
      free (tokgroup);
      *tokgroupsize = 0;
      return (0);
    }
    permuteLine (buffptr = buff);
  }
  if ((semicolon = index (buffptr, ';')) != NULL)
    *(semicolon++) = '\0';

  if ((ca = strtok_r (buffptr, " ", &cb)) == NULL) {
    tokgroup[0] = buffptr;
    numtoks = 1;
  } else {
    while (ca != NULL) {
      if (numtoks + 1 >= *tokgroupsize) {
	*tokgroupsize += 10;
	if ((tokgroup = realloc (tokgroup, sizeof (char *) * *tokgroupsize)) == NULL) {
	  logMsg (LOGDEFAULT, LOGFATAL, "malloc failed\n");
	  exit (-1);
	}
      }
      tokgroup[numtoks++] = ca;
      ca = strtok_r (NULL, " ", &cb);
    }
  }
  
  *tokgroup_h = tokgroup;
  buffptr = semicolon;
  return (numtoks);
}


#define STARTOFLINE  0
#define INSTRING     1
#define INWHITESPACE 2
#define INSEPARATOR  3

/* Permutes a line of input. Leading and trailing whitespace is deleted;
 * whitespace that brackets 'separator' characters (comma, hyphen, semicolon or
 * colon) is deleted; all other whitespace is reduced to a single space;
 * comment characters (pound sign) and all subsequent characters up to end of
 * line are deleted, unless the pound sign is 'escacped' with a backslash.
 */
void permuteLine (char *line)
{
  int state, va, vb;

  va = vb = 0;
  state = STARTOFLINE;

  while (1) {
    //printf ("va %d is '%c', vb %d is '%c' state is %d -> ", va, line[va], vb, line[vb], state);
    if (index (" \t", line[vb]) != NULL) {
      if (state == INSTRING)
	line[va++] = ' ';
      vb++;
      if (state != STARTOFLINE)
	state = INWHITESPACE;
      
    } else if (index ("-:;,", line[vb]) != NULL) {
      if (state == INWHITESPACE)
	line[va-1] = line[vb++];
      else 
	line[va++] = line[vb++];
      state = INSEPARATOR;

    } else if (line[vb] == '#') {
      if ((vb == 0) || (line[vb-1] != '\\')) {
	if (state == INWHITESPACE)
	  va--;
	line[va] = '\0';
	break;
      }
      line[va++] = line[vb++];
      state = INSTRING;

    } else if ((index ("\n\r", line[vb]) != NULL) || (line[vb] == '\0')) {
      if (state == INWHITESPACE)
	va--;
      line[va] = '\0';
      break;

    } else {
      line[va++] = line[vb++];
      state = INSTRING;
    }
    //printf ("va %d is '%c', vb %d is '%c' state is %d\n", va, line[va], vb, line[vb], state);
  }
  //printf ("all done: va %d, vb %d, line = '%s'\n", va, vb, line);
  return;
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
