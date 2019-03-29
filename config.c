/**********************************************************************
 * Multiprocessor Linkage Analysis
 * Alberto Maria Segre
 * Regex code originally by Nathan Burnette
 * 
 * Copyright 2006, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software 
 * for non-profit educational purposes only.
 **********************************************************************/
#include "kelvin.h"
#include <sys/types.h>		/* C regexps */
#include <regex.h>		/* C regexps */

/**********************************************************************
 * Process configuration file and set up model datastructures
 * (modelType, range, modelOptions) describing the experimental
 * trial.
 *
 * TODO: check regexps to extract more, rather than just recognize, input
 * TODO: need to validate inputs: e.g., negative numbers allowed for DD when QT
 * TODO: handling of multiallelic diseases; parsing e.g., 00, 21, 11 etc
 * TODO: handling of LC for AF, and extend AF to more than 2 marker alleles (non SNP)
 * TODO: check distribs with nparam>1 (not sure it works)
 * TODO: liability class should be the slowest moving index so that we
 *       can censor at run-time if ped has no such LC?
 **********************************************************************/

/**********************************************************************
 * Regular expressions used to parse configuration file
 * lines. 
 **********************************************************************/
/* DIRECTIVE matches: GF, DD, Dd, dd, Th, Tm, Tf, 00..99, AF, AL, MM, AM, TL, TM, TT, LD */
#define DIRECTIVE "^[[:space:]]*[[:digit:]AdDLMGT][[:digit:]dDfFhLMmT]"
/* PARAMETER matches: P[0-9] (note implicit 1 digit limit on distrib params index) */
#define PARAMETER "^[[:space:]]*P([[:digit:]])"
/* Match 1 or 3 floating point numbers with optional trailing semicolon. */
#define DOUBLE1 "[[:space:]]*([-]?[[:digit:]]*.?[[:digit:]]+)[[:space:]]*[;]?"
#define DOUBLE3 "[[:space:]]*([-]?[[:digit:]]*.?[[:digit:]]+)[[:space:]]+([-]?[[:digit:]]*.?[[:digit:]]+)[[:space:]]+([-]?[[:digit:]]*.?[[:digit:]]+)[[:space:]]*[;]?"
/* Match a constraint of the form, e.g., DD < Dd */
#define CONSTRAINT2 "[[:space:]]*([[:digit:]dDT][[:digit:]dDfmT])[[:space:]]*([>]|[>][=]|[=][=]|[!][=])[[:space:]]*([[:digit:]dDT][[:digit:]dDfmT])[[:space:]]*[;]?"
/* Match a constraint with liability classes, e.g., DD 1 < DD 2 */
#define CLASSCONSTRAINT2 "[[:space:]]*([[:digit:]dDT][[:digit:]dDfmT])[[:space:]]*([[:digit:]])[[:space:]]*([>]|[>][=]|[=][=]|[!][=])[[:space:]]*([[:digit:]dDT][[:digit:]dDfmT])[[:space:]]*([[:digit:]])[[:space:]]*[;]?"
/* Match a constraint of the form, e.g., Px DD < Px Dd */
#define PCONSTRAINT2 "[[:space:]]*P([[:digit:]])[[:space:]]+([[:digit:]dDT][[:digit:]dDfm])[[:space:]]*([>]|[>][=]|[=][=]|[!][=])[[:space:]]*P([[:digit:]])[[:space:]]+([[:digit:]dDT][[:digit:]dDfm])[[:space:]]*[;]?"
/* Match a constraint of the form, e.g., Px DD 1 < Px Dd 2 */
#define PCLASSCONSTRAINT2 "[[:space:]]*P([[:digit:]])[[:space:]]+([[:digit:]dDT][[:digit:]dDfm])[[:space:]]*([[:digit:]])[[:space:]]*([>]|[>][=]|[=][=]|[!][=])[[:space:]]*P([[:digit:]])[[:space:]]+([[:digit:]dDT][[:digit:]dDfm])[[:space:]]*([[:digit:]])[[:space:]]*[;]?"

/* Strings used to parse configuration file lines (see getType). */
#define DIRECTIVES "GFThTmTfALTLTTLDAFDDDddd"
#define GF 0			/* gene frequency */
#define Th 1			/* theta */
#define Tm 2			/* male theta */
#define Tf 3			/* female theta */
#define AL 4			/* alpha */
#define TL 5			/* trait loci */
#define TT 6			/* trait threshold */
#define LD 7			/* d prime */
#define AF 8
#define DD 9			/* Must be last for numeric directives to work */
#define Dd 10
#define dd 11
/* Strings used to parse configuration file constraints (see
 * getOperator). Matches:
 *   ==, !=, >, >=
 *
 * TODO: Should probably allow a single = as an alias for == to allow
 * for non-programmers?
 *
 * TODO: Add <, <= ?? */
#define OPERATORS "==!=> >="
#define EQ 0
#define NE 1
#define GT 2
#define GE 3

/* Used in parsing configuration file constraints. */
#define SEXML 0
#define SEXFM 1
#define SEXAV 0

/**********************************************************************
 * Structure used for linked lists of constraints. A constraint can be
 * of several different forms:
 *   a1 op a2 			DD > Dd			SIMPLE
 *   a1 c1 op a2 c2		Dd 1 != dd 2		CLASSC
 *   p1 a1 op p2 a2    		P2 DD > P1 DD		PARAMC
 *   p1 a1 c1 op p2 a2 c2       P1 dd 1 >= P2 dd 2	PARAMCLASSC
 * but they are all stored in a uniform constraint structure.
 * where arg is, e.g., Tm or DD and op is one of the operators defined
 * below. If more than one constraint appears on a line, then we take
 * that to be an implicit OR, and all but the last constraint on such
 * a line will be marked alt=TRUE.
 ***********************************************************************/
#define SIMPLE 0
#define CLASSC 1
#define PARAMC 2
#define PARAMCLASSC 3
typedef struct constraint
{
  int type;			/* SIMPLE, CLASSC, PARAMC, PARAMCLASSC */
  int a1;
  int c1;
  int p1;
  int a2;
  int c2;
  int p2;
  int op;
  int alt;			/* Is this one of a disjunctive set? */
}
Constraint;

/**********************************************************************
 * Some global variables used only within config.c -- avoids having to
 * carry these around in the model structure, since they are only
 * useful while setting up the penetrances, etc.
 **********************************************************************/
int maxgfreq;			/* Max number of gfreqs in model (for dynamic alloc) */
int maxalpha;			/* Max number of alphas in model (for dynamic alloc) */
int maxtloc;			/* Max number of trait loci in model (for dynamic alloc) */
int maxtthresh;			/* Max number of trait thresholds in model (for dynamic alloc) */
int maxafreq;			/* Max number of afreqs in model (for dynamic alloc) */
int maxdprime;			/* Max number of dprimes in model (for dynamic alloc) */
int *penetcnt;			/* Array of number of penetrances */
int *penetmax;			/* Array of max penetrances */
int *paramcnt;			/* Number of QT/CT parameters */
int *parammax;			/* Max QT/CT parameters */
int *thetacnt;			/* Array of number of thetas */
int *thetamax;			/* Array of max thetas */
Constraint *constraints[4];	/* Array of constraints by type */
int constmax[4] = { 0, 0, 0, 0 };	/* Max constraints in array (for dynamic alloc) */
int constcnt[4] = { 0, 0, 0, 0 };	/* Current number of constraints in array */

/* Chunk size used in reallocating arrays of gene frequencies,
 * penetrances, thetas, and constraints. */
#define CHUNKSIZE 64

/**********************************************************************
 * Some internal prototypes.
 **********************************************************************/
int getType (char *line);
int getOperator (char *line);
int getInteger (char *line);
void addRange (ModelRange * range, int type, double lo, double hi,
	       double incr);
void addTheta (ModelRange * range, int type, double val);
void addPenetrance (ModelRange * range, int type, double val);
void addGeneFreq (ModelRange * range, double val);
void addAlpha (ModelRange * range, double val);
/* void addTraitLocus (ModelRange *range, double val); */
void addTraitThreshold (ModelRange * range, double val);
void addDPrime (ModelRange * range, double val);
void addAlleleFreq (ModelRange * range, double val);
void addConstraint (int type, int a1, int c1, int p1,
		    int op, int a2, int c2, int p2, int disjunct);
void addParameter (ModelRange * range, int dim, double val);
int checkThetas (ModelRange * range, int i);
int checkPenets (ModelRange * range, int i);
int checkClassPenets (ModelRange * range, int i);
int checkParams (ModelRange * range, int i);
int checkClassParams (ModelRange * range, int i);
int checkClassThreshold (ModelRange * range, int i);
void sortRange (ModelRange * range);
inline void swap (double *array, int i, int j);
void quicksort (double *array, int lo, int hi);
void uniqRange (ModelRange * range);
inline int uniquify (double *array, int len);
void expandRange (ModelRange * range, ModelType * type);
void expandClass (ModelRange * range, ModelType * type);
LambdaCell *findLambdas (ModelRange * range, int m, int n);
void showRange (ModelRange * range, ModelType * type, int level);
void showConstraints ();

/**********************************************************************
 * Read configuration file. Returns ERROR if there is a problem
 * opening the file.
 *
 * Parsing of the configuration file is done mostly with scanf's, but
 * we do use C regular expressions to parse the arguments where
 * necessary. This gives us the flexibility we need to read in
 * arbitrary patterns of arguments from a single line where necessary.
 **********************************************************************/
int
readConfigFile (char *file, ModelType * modelType,
		ModelRange * modelRange, ModelOptions * modelOptions)
{
  FILE *fp;			/* Filepointer. */
  int i = 0;			/* Number of lines read. */
  char *start;
  char line[KMAXLINELEN + 1];	/* Current line. */
  int dir1, dir2, op;		/* For parsing constraints. */
  int a1, a2, c1, c2, p1;

  regex_t *buffer0;
  regex_t *buffer1;
  regex_t *buffer2;
  regex_t *buffer3;
  regex_t *buffer4;
  regex_t *buffer5;
  regex_t *buffer6;
  regex_t *buffer7;
  regmatch_t match[8];		/* Store extracted values. */

  /* Set up the default model values. */
  modelType->type = TP;
  modelType->trait = DT;
  modelOptions->equilibrium = LINKAGE_EQUILIBRIUM;
  modelOptions->markerAnalysis = FALSE;
  modelOptions->polynomial = FALSE;
  modelRange->nalleles = 2;
  modelRange->nlclass = 1;
  modelRange->npardim = 0;
  modelRange->nlambdas = 0;
  modelRange->maxnlambdas = 0;
  modelRange->tlmark = FALSE;

  /* Allocate storage for pattern spaces. */
  buffer0 = malloc (KMAXLINELEN + 1);
  buffer1 = malloc (KMAXLINELEN + 1);
  buffer2 = malloc (KMAXLINELEN + 1);
  buffer3 = malloc (KMAXLINELEN + 1);
  buffer4 = malloc (KMAXLINELEN + 1);
  buffer5 = malloc (KMAXLINELEN + 1);
  buffer6 = malloc (KMAXLINELEN + 1);
  buffer7 = malloc (KMAXLINELEN + 1);

  /* Open the configuration file for reading and punt if you can't. */
  KASSERT ((fp = fopen (file, "r")),
	   "Can't open configuration file %s; aborting.\n", file);

  /* Set up the regular expression pattern buffer we'll use to parse
   * the lines in the file. The flags REG_EXTENDED and REG_ICASE mean
   * we'll use Posix extended regular expressions and we'll ignore
   * case. */
  KASSERT ((regcomp (buffer0, DIRECTIVE, (REG_EXTENDED | REG_NEWLINE)) == 0),
	   "Internal error in regular expression; aborting.\n");
  KASSERT ((regcomp (buffer1, DOUBLE1, (REG_EXTENDED | REG_NEWLINE)) == 0),
	   "Internal error in regular expression; aborting.\n");
  KASSERT ((regcomp (buffer2, DOUBLE3, (REG_EXTENDED | REG_NEWLINE)) == 0),
	   "Internal error in regular expression; aborting.\n");
  KASSERT ((regcomp (buffer3, CONSTRAINT2, (REG_EXTENDED | REG_NEWLINE)) ==
	    0), "Internal error in regular expression; aborting.\n");
  KASSERT ((regcomp (buffer4, CLASSCONSTRAINT2, (REG_EXTENDED | REG_NEWLINE))
	    == 0), "Internal error in regular expression; aborting.\n");
  KASSERT ((regcomp (buffer5, PARAMETER, (REG_EXTENDED | REG_NEWLINE)) == 0),
	   "Internal error in regular expression; aborting.\n");
  KASSERT ((regcomp (buffer6, PCONSTRAINT2, (REG_EXTENDED | REG_NEWLINE)) ==
	    0), "Internal error in regular expression; aborting.\n");
  KASSERT ((regcomp (buffer7, PCLASSCONSTRAINT2, (REG_EXTENDED | REG_NEWLINE))
	    == 0), "Internal error in regular expression; aborting.\n");

  /* Start scanning each line of configuration file input. We'll parse
   * each line by looking for the easy cases first, then using C
   * regexps to parse the harder configuration directives. */
  while (fgets (line, KMAXLINELEN, fp))
    {
      /* Before we try to parse it, check to see if the line may be
       * too long, and give up if it is. */
      KASSERT ((strlen (line) < KMAXLINELEN),
	       "Line %d in configuration file %s exceeds %d characters; aborting.\n",
	       i, file, KMAXLINELEN);

      /* Flush lines starting with a comment character or consisting
       * of only a newline. */
      if ((strlen (line) == 1) || (strncmp (line, "#", 1) == 0))
	continue;

      /* TODO: flushing lines containing only whitespace characters
       * will make everything slightly more efficient. */

      /* Type of analysis; 2 point or multipoint. */
      if (strncmp (line, "TP", 2) == 0)
	{
	  modelType->type = TP;	/* 2 point (default) */
	  KLOG (LOGINPUTFILE, LOGDEBUG, "Configuring for 2 point analysis\n");
	  continue;
	}
      if (sscanf (line, "SS %d", &(modelType->numMarkers)) == 1)
	{
	  modelType->type = MP;	/* Multipoint */
	  modelOptions->mapFlag = SS;
	  KLOG (LOGINPUTFILE, LOGDEBUG,
		"Configuring for sex-specific multipoint analysis\n");
	  continue;
	}
      if (sscanf (line, "SA %d", &(modelType->numMarkers)) == 1)
	{
	  modelType->type = MP;	/* Multipoint */
	  modelOptions->mapFlag = SA;
	  KLOG (LOGINPUTFILE, LOGDEBUG,
		"Configuring for sex-averaged multipoint analysis\n");
	  continue;
	}
      if (strncmp (line, "TM", 2) == 0)
	{
	  modelRange->tlmark = TRUE;
	  KLOG (LOGINPUTFILE, LOGDEBUG,
		"Configuring for on-marker trait loci\n");
	  continue;
	}
      if (strncmp (line, "MM", 2) == 0)
	{
	  modelOptions->markerAnalysis = MM;
	  KLOG (LOGINPUTFILE, LOGDEBUG,
		"Configuring for marker to marker analysis\n");
	  continue;
	}
      if (strncmp (line, "AM", 2) == 0)
	{
	  modelOptions->markerAnalysis = AM;
	  KLOG (LOGINPUTFILE, LOGDEBUG,
		"Configuring for adjacent marker analysis\n");
	  continue;
	}
      if (strncmp (line, "PE", 2) == 0)
	{
	  modelOptions->polynomial = TRUE;	/* Polynomial evaluation */
	  KLOG (LOGINPUTFILE, LOGDEBUG,
		"Configuring for polynomial evaluation\n");
	  continue;
	}
      if (sscanf (line, "LC %d", &modelRange->nlclass) == 1)
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG,
		"Configuring for %d liability classes\n",
		modelRange->nlclass);
	  continue;
	}

      if (sscanf (line, "DA %d", &modelRange->nalleles) == 1)	/* Disease alleles */
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG,
		"Configuring for %d disease alleles\n", modelRange->nalleles);
	  continue;
	}

      if (sscanf (line, "AS %lg %lg %lg",	/* Affection status values */
		  &(modelOptions->affectionStatus[AFFECTION_STATUS_UNKNOWN]),
		  &(modelOptions->
		    affectionStatus[AFFECTION_STATUS_UNAFFECTED]),
		  &(modelOptions->
		    affectionStatus[AFFECTION_STATUS_AFFECTED])) == 3)
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG,
		"Resetting affection status values (%g, %g, %g)\n",
		modelOptions->affectionStatus[AFFECTION_STATUS_UNKNOWN],
		modelOptions->affectionStatus[AFFECTION_STATUS_UNAFFECTED],
		modelOptions->affectionStatus[AFFECTION_STATUS_AFFECTED]);
	  continue;
	}

      /* Dichotomous trait directive. */
      if (strncmp (line, "DT", 2) == 0)
	{
	  modelType->trait = DT;	/* Dichotomous trait */
	  /* Establish the default affected, unaffected, and unknown
	   * values for DT. These can be overridden elsewhere. */
	  modelOptions->affectionStatus[AFFECTION_STATUS_UNKNOWN] =
	    AFFECTION_STATUS_UNKNOWN;
	  modelOptions->affectionStatus[AFFECTION_STATUS_UNAFFECTED] =
	    AFFECTION_STATUS_UNAFFECTED;
	  modelOptions->affectionStatus[AFFECTION_STATUS_AFFECTED] =
	    AFFECTION_STATUS_AFFECTED;
	  KLOG (LOGINPUTFILE, LOGDEBUG,
		"Configuring for dichotomous traits\n");
	  continue;
	}

      /* Quantitative trait directives; each different distribution
       * may require a different pattern of parameters. */
      if (sscanf (line, "QT normal %lf %lf", &modelType->mean, &modelType->sd)
	  == 2)
	{
	  modelType->trait = QT;	/* Quantitative trait */
	  modelType->distrib = NORMAL_DISTRIBUTION;
	  /* Establish the default affected, unaffected, and unknown
	   * values for QT/CT. These can be overridden elsewhere. */
	  modelOptions->affectionStatus[AFFECTION_STATUS_UNKNOWN] = NAN;
	  modelOptions->affectionStatus[AFFECTION_STATUS_UNAFFECTED] =
	    -9999.99;
	  modelOptions->affectionStatus[AFFECTION_STATUS_AFFECTED] = 9999.99;
	  /* The normal distribution has a two distributional
	   * parameters, mean (specified as the penetrance) and std
	   * dev, specified as the first additional parameter P1. */
	  modelRange->npardim = 1;
	  KLOG (LOGINPUTFILE, LOGDEBUG,
		"Configuring for quantitative traits (normal distribution)\n");
	  continue;
	}
      if (sscanf
	  (line, "QT T %d %lf %lf", &a1, &modelType->mean,
	   &modelType->sd) == 3)
	{
	  modelType->trait = QT;	/* Quantitative trait */
	  modelType->distrib = T_DISTRIBUTION;
	  /* The T distribution has single distribution constant, the
	   * degrees of freedom, which is fixed at definition. We'll
	   * use integer a1 temporarily so as to set up the
	   * appropriate number of constants in modelType. */
	  modelType->constants =
	    realloc (modelType->constants, 1 * sizeof (int));
	  modelType->constants[0] = a1;
	  /* Establish the default affected, unaffected, and unknown
	   * values for QT/CT. These can be overridden elsewhere. */
	  modelOptions->affectionStatus[AFFECTION_STATUS_UNKNOWN] = NAN;
	  modelOptions->affectionStatus[AFFECTION_STATUS_UNAFFECTED] =
	    -9999.99;
	  modelOptions->affectionStatus[AFFECTION_STATUS_AFFECTED] = 9999.99;

	  /* The T distribution, like the normal distribution, also
	   * has two distributional parameters, the mean (specified as
	   * the penetrance) and the std dev, specified as the first
	   * additional parameter P1. */
	  modelRange->npardim = 1;
	  KLOG (LOGINPUTFILE, LOGDEBUG,
		"Configuring for quantitative traits (T distribution)\n");
	  continue;
	}

      /* Directives that take a single string argument. */
      if (sscanf (line, "PD %s", pedfile) == 1)	/* Pedigree file */
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG, "Configure pedigree file %s\n",
		pedfile);
	  continue;
	}
      if (sscanf (line, "DF %s", datafile) == 1)	/* Data file - a list of loci */
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG, "Configure data file %s\n", datafile);
	  continue;
	}
      if (sscanf (line, "MK %s", markerfile) == 1)	/* Marker file - marker allele frequencies */
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG, "Configure marker file %s\n",
		markerfile);
	  continue;
	}
      if (sscanf (line, "MP %s", mapfile) == 1)	/* Map file */
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG, "Configure map file %s\n", mapfile);
	  continue;
	}
#if FALSE
      if (sscanf (line, "LP %s", loopfile) == 1)
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG, "Configure loop file %s\n", loopfile);
	  continue;
	}
#endif
      if (sscanf (line, "OF %s", outfile) == 1)	/* Output file */
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG, "Configure output file %s\n",
		outfile);
	  continue;
	}
      if (sscanf (line, "HE %s", avghetfile) == 1)	/* Average hetergeneity LR file */
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG, "Configure output file %s\n",
		avghetfile);
	  continue;
	}
      if (sscanf (line, "HO %s", avghomofile) == 1)	/* Average homogeneity LR file */
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG, "Configure output file %s\n",
		avghomofile);
	  continue;
	}
      if (sscanf (line, "PF %s", pplfile) == 1)	/* PPL output file */
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG, "Configure PPL output file %s\n",
		pplfile);
	  continue;
	}

      if (sscanf (line, "LF %s", ldPPLfile) == 1)	/* LD-PPL output file */
	{
	  KLOG (LOGINPUTFILE, LOGDEBUG, "Configure LD-PPL output file %s\n",
		ldPPLfile);
	  continue;
	}

      if (sscanf (line, "UP %d", &dir1) == 1)
	{
	  /* Allocate space for the identifier based on how many
	   * digits in the scanned integer. Remember to leave a space
	   * for the terminator. */
	  modelOptions->sUnknownPersonID =
	    realloc (modelOptions->sUnknownPersonID,
		     (dir1 / 10) + 2 * sizeof (char));
	  snprintf (modelOptions->sUnknownPersonID, (dir1 / 10) + 2, "%d",
		    dir1);
	  KLOG (LOGINPUTFILE, LOGDEBUG, "Unknown person identifier is %s\n",
		modelOptions->sUnknownPersonID);
	  continue;
	}

      /* Complex directives that require regular expression matching
       * to handle the arguments. Most of these take at least one
       * floating point argument, and sometimes two or three.
       *
       * We'll handle each type of directive separately. But, first,
       * "chomp" the newline to remove the newline and ensure null
       * termination. */
      line[strlen (line) - 1] = '\0';

      /* First, look for simple constraints on thetas, gene
       * frequencies, or penetrances. */
      if (regexec (buffer3, line, 4, &match[0], 0) != REG_NOMATCH)
	{
	  /* Matching constraint line (CONSTRAINT2):
	   *   DD < Dd 
	   * Find out what type of constraint it is by looking it up
	   * in DIRECTIVES. */
	  KASSERT (((dir1 = getType (line)) != ERROR),
		   "Bad constraint directive '%s'; aborting.\n", line);
	  KASSERT (((op = getOperator (line + match[2].rm_so)) != ERROR),
		   "Bad constraint operator '%s'; aborting.\n", line);
	  KASSERT (((dir2 = getType (line + match[3].rm_so)) != ERROR),
		   "Bad constraint directive '%s'; aborting.\n", line);
	  /* Constraints not supported for AL or TL */
	  KASSERT (((dir1 != AL) && (dir2 != AL)),
		   "Bad constraint directive '%s'; aborting.\n", line);
	  KASSERT (((dir1 != TL) && (dir2 != TL)),
		   "Bad constraint directive '%s'; aborting.\n", line);

	  /* Add the constraint. */
	  addConstraint (SIMPLE, dir1, 0, 0, op, dir2, 0, 0, FALSE);

	  /* Next, loop through any semicolon-separated disjuncts. */
	  start = line;
	  while ((start = start + match[0].rm_eo) < line + strlen (line))
	    {
	      /* Add the disjunction to the first constraint. */
	      if (regexec (buffer3, start, 4, &match[0], REG_NOTBOL) !=
		  REG_NOMATCH)
		{
		  KASSERT (((dir1 =
			     getType (start + match[1].rm_so)) != ERROR),
			   "Bad constraint directive '%s'; aborting.\n",
			   start);
		  KASSERT (((op =
			     getOperator (start + match[2].rm_so)) != ERROR),
			   "Bad constraint directive '%s'; aborting.\n",
			   start);
		  KASSERT (((dir2 =
			     getType (start + match[3].rm_so)) != ERROR),
			   "Bad constraint directive '%s'; aborting.\n",
			   start);
		  addConstraint (SIMPLE, dir1, 0, 0, op, dir2, 0, 0, TRUE);
		}
	      else
		break;
	    }
	  continue;
	}
      else if (regexec (buffer4, line, 6, &match[0], 0) != REG_NOMATCH)
	{
	  /* Next, look for class constraints on penetrances, gene
	   * frequencies, thetas, or thresholds (CLASSCONSTRAINT2):
	   *   Px DD 1 < Px Dd 2 
	   * Find out what type of constraint it is by looking it up
	   * in DIRECTIVES. */
	  KASSERT (((dir1 = getType (line)) != ERROR),
		   "Bad constraint directive '%s'; aborting.\n", line);
	  KASSERT (((op = getOperator (line + match[3].rm_so)) != ERROR),
		   "Bad constraint operator '%s'; aborting.\n", line);
	  KASSERT (((dir2 = getType (line + match[4].rm_so)) != ERROR),
		   "Bad constraint directive '%s'; aborting.\n", line);
	  /* Constraints not supported for AL or TL */
	  KASSERT (((dir1 != AL) && (dir2 != AL)),
		   "Bad constraint directive '%s'; aborting.\n", line);
	  KASSERT (((dir1 != TL) && (dir2 != TL)),
		   "Bad constraint directive '%s'; aborting.\n", line);

	  /* Add the constraint. */
	  addConstraint (CLASSC,
			 dir1, getInteger (line + match[2].rm_so), 0,
			 op,
			 dir2, getInteger (line + match[5].rm_so), 0, FALSE);

	  /* Next, loop through any semicolon-separated disjuncts. */
	  start = line;
	  while ((start = start + match[0].rm_eo) < line + strlen (line))
	    {
	      /* Add the disjunction to the first constraint. */
	      if (regexec (buffer4, start, 6, &match[0], REG_NOTBOL) !=
		  REG_NOMATCH)
		{
		  KASSERT (((dir1 =
			     getType (start + match[1].rm_so)) != ERROR),
			   "Bad constraint directive '%s'; aborting.\n",
			   start);
		  KASSERT (((op =
			     getOperator (start + match[3].rm_so)) != ERROR),
			   "Bad constraint directive '%s'; aborting.\n",
			   start);
		  KASSERT (((dir2 =
			     getType (start + match[4].rm_so)) != ERROR),
			   "Bad constraint directive '%s'; aborting.\n",
			   start);
		  addConstraint (CLASSC, dir1,
				 getInteger (start + match[2].rm_so), -1, op,
				 dir2, getInteger (start + match[5].rm_so),
				 -1, TRUE);
		}
	      else
		break;
	    }
	  continue;
	}
      else if (regexec (buffer6, line, 6, &match[0], 0) != REG_NOMATCH)
	{
	  /* Now look for constraints on parameters (PCONSTRAINT2):
	   *   Px DD < Px Dd */
	  dir1 = getInteger (line + match[1].rm_so);
	  a1 = getType (line + match[2].rm_so);
	  op = getOperator (line + match[3].rm_so);
	  dir2 = getInteger (line + match[4].rm_so);
	  a2 = getType (line + match[5].rm_so);

	  /* SHOULD CHECK CONSTRAINT FORMAT! */

	  /* Add the constraint */
	  addConstraint (PARAMC, a1, 0, dir1, op, a2, 0, dir2, FALSE);

	  /* Next, loop through any semicolon-separated disjuncts. */
	  start = line;
	  while ((start = start + match[0].rm_eo) < line + strlen (line))
	    {
	      /* Add the disjunction to the first constraint. */
	      if (regexec (buffer6, start, 6, &match[0], REG_NOTBOL) !=
		  REG_NOMATCH)
		{
		  /* Matching parameter constraint line. */
		  dir1 = getInteger (start + match[1].rm_so);
		  a1 = getType (start + match[2].rm_so);
		  op = getOperator (start + match[3].rm_so);
		  dir2 = getInteger (start + match[4].rm_so);
		  a2 = getType (start + match[5].rm_so);

		  /* Add the constraint */
		  addConstraint (PARAMC, a1, 0, dir1, op, a2, 0, dir2, TRUE);
		}
	      else
		break;
	    }
	  continue;
	}
      else if (regexec (buffer7, line, 8, &match[0], 0) != REG_NOMATCH)
	{
	  /* Finally, look for cross-class constraints on
	   * parameters: (PCLASSCONSTRAINT2):
	   *   Px DD 1 < Px Dd 2 */
	  dir1 = getInteger (line + match[1].rm_so);
	  a1 = getType (line + match[2].rm_so);
	  c1 = getInteger (line + match[3].rm_so);
	  op = getOperator (line + match[4].rm_so);
	  dir2 = getInteger (line + match[5].rm_so);
	  a2 = getType (line + match[6].rm_so);
	  c2 = getInteger (line + match[7].rm_so);

	  /* SHOULD CHECK CONSTRAINT FORMAT! */

	  /* Add the constraint */
	  addConstraint (PARAMCLASSC, a1, c1, dir1, op, a2, c2, dir2, FALSE);

	  /* Next, loop through any semicolon-separated disjuncts. */
	  start = line;
	  while ((start = start + match[0].rm_eo) < line + strlen (line))
	    {
	      /* Add the disjunction to the first constraint. */
	      if (regexec (buffer6, start, 6, &match[0], REG_NOTBOL) !=
		  REG_NOMATCH)
		{
		  /* Matching parameter constraint line. */
		  dir1 = getInteger (start + match[1].rm_so);
		  a1 = getType (start + match[2].rm_so);
		  c1 = getInteger (start + match[3].rm_so);
		  op = getOperator (start + match[4].rm_so);
		  dir2 = getInteger (start + match[5].rm_so);
		  a2 = getType (start + match[6].rm_so);
		  c2 = getInteger (start + match[7].rm_so);

		  /* Add the constraint */
		  addConstraint (PARAMCLASSC, dir1, a1, c1, op,
				 dir2, a2, c2, TRUE);
		}
	      else
		break;
	    }
	  continue;
	}
      else if (regexec (buffer0, line, 1, &match[0], 0) != REG_NOMATCH)
	{
	  /* Look for value specifications (DIRECTIVE). 
	   *   xx <number or numbers>
	   * Find out what type of directive it is by looking it up in
	   * DIRECTIVES. */
	  dir1 = getType (line);

	  /* Next, loop through the semicolon-separated arguments. */
	  start = line + match[0].rm_eo;
	  while (TRUE)
	    {
	      if (regexec (buffer2, start, 4, &(match[1]), REG_NOTBOL) !=
		  REG_NOMATCH)
		{
		  /* Matches three numbers. */
		  addRange (modelRange, dir1,
			    strtod (start + match[2].rm_so, NULL),
			    strtod (start + match[3].rm_so, NULL),
			    strtod (start + match[4].rm_so, NULL));
		}
	      else if (regexec (buffer1, start, 2, &(match[1]), REG_NOTBOL) !=
		       REG_NOMATCH)
		{
		  /* Matches one number. */
		  if (dir1 == Th || dir1 == Tm || dir1 == Tf)
		    addTheta (modelRange, dir1,
			      strtod (start + match[2].rm_so, NULL));
		  else if (dir1 == GF)
		    addGeneFreq (modelRange,
				 strtod (start + match[2].rm_so, NULL));
		  else if (dir1 == AL)
		    addAlpha (modelRange,
			      strtod (start + match[2].rm_so, NULL));
		  else if (dir1 == TL)
		    addTraitLocus (modelRange,
				   strtod (start + match[2].rm_so, NULL));
		  else if (dir1 == TT)
		    addTraitThreshold (modelRange,
				       strtod (start + match[2].rm_so, NULL));
		  else if (dir1 == LD)
		    addDPrime (modelRange,
			       strtod (start + match[2].rm_so, NULL));
		  else if (dir1 == AF)
		    addAlleleFreq (modelRange,
				   strtod (start + match[2].rm_so, NULL));
		  else
		    /* All parsing is done only in pre-expansion mode,
		     * hence class is always moot. Also, subtract base
		     * value of DD from dir1 to get 0-offset value.*/
		    addPenetrance (modelRange, dir1 - DD,
				   strtod (start + match[2].rm_so, NULL));
		}
	      else
		{
		  /* Doesn't match anything. This is a badly formatted
		   * line which should be flagged. */
		  KASSERT (FALSE,
			   "Ill-formed line in configuration file: '%s'\n",
			   line);
		}

	      /* Update the start pointer and check for termination. */
	      start = start + match[1].rm_eo;
	      if (start >= line + strlen (line))
		break;
	    }
	  continue;
	}
      else if (regexec (buffer5, line, 1, &match[0], 0) != REG_NOMATCH)
	{
	  /* Matching a parameter specification.
	   *   P[0-9] <number or numbers>
	   * Find out which parameter it is by parsing the parameter
	   * directive (the +1 ensures we skip over the leading P in,
	   * e.g., P1).  Remember, the parameter number read in will
	   * be 1-indexed, and you need to make it 0-indexed first. */
	  p1 = strtod (line + match[0].rm_so + 1, NULL) - 1;
	  /* Make sure that the parameter number is within the
	   * expected number of parameters for this distribution. */
	  KASSERT ((p1 < modelRange->npardim), "Illegal parameter P%d.\n",
		   p1 + 1);

	  /* Next, loop through the semicolon-separated arguments. */
	  start = line + match[0].rm_eo;
	  while (TRUE)
	    {
	      if (regexec (buffer2, start, 4, &(match[1]), REG_NOTBOL) !=
		  REG_NOMATCH)
		{
		  /* Matches three numbers. We'll fake out addRange by
		   * giving it a negative type to distinguish it from
		   * the other types of values we handle. So, for
		   * example, parameter 2 would be -3. */
		  addRange (modelRange, -(p1 + 1),
			    strtod (start + match[2].rm_so, NULL),
			    strtod (start + match[3].rm_so, NULL),
			    strtod (start + match[4].rm_so, NULL));
		}
	      else if (regexec (buffer1, start, 2, &(match[1]), REG_NOTBOL) !=
		       REG_NOMATCH)
		{
		  /* Matches one number. */
		  addParameter (modelRange, p1,
				strtod (start + match[2].rm_so, NULL));
		}
	      else
		{
		  /* Doesn't match anything. This is a badly formatted
		   * line which should be flagged. */
		  KASSERT (FALSE,
			   "Ill-formed line in configuration file: '%s'\n",
			   line);
		}

	      /* Update the start pointer and check for termination. */
	      start = start + match[1].rm_eo;
	      if (start >= line + strlen (line))
		break;
	    }
	  continue;
	}
      else
	KASSERT (FALSE, "Ill-formed line in configuration file: '%s'\n",
		 line);
    }

  /* Done. Release the pattern spaces and free the buffers. */
  regfree (buffer0);
  regfree (buffer1);
  regfree (buffer2);
  regfree (buffer3);
  regfree (buffer4);
  regfree (buffer5);
  regfree (buffer6);
  regfree (buffer7);
  free (buffer0);
  free (buffer1);
  free (buffer2);
  free (buffer3);
  free (buffer4);
  free (buffer5);
  free (buffer6);
  free (buffer7);

  /* Clean up after yourself. Here, we set/reset whatever model
   * options are implict by the configuration parameters given. The
   * only example that comes to mind is linkage disequilibrium, but
   * there may later be others. */
  if (modelRange->dprime && modelOptions->equilibrium == LINKAGE_EQUILIBRIUM)
    {
      /* If dprime exists, we must have set it explicitly, and we must
       * be doing LD. */
      modelOptions->equilibrium = LINKAGE_DISEQUILIBRIUM;	/* Linkage disequilibrium */
      KLOG (LOGINPUTFILE, LOGDEBUG,
	    "Configuring for linkage disequilibrium\n");
    }

  /* Now check the integrity of the parameters you've read. Here is
   * where you check for things like, e.g., no parameters specified
   * for QT/CT, or parameters specified for DT. */
  KASSERT ((modelType->trait != QT || modelRange->param),
	   "Failure to provide distribution parameters for quantitative trait.\n");
  KASSERT ((modelType->trait != CT || modelRange->param),
	   "Failure to provide distribution parameters for combined trait.\n");
  KASSERT ((modelType->trait != DT || !modelRange->param),
	   "Attempt to provide distribution parameters for dichotomous trait.\n");
  KASSERT ((modelType->type == TP || !modelRange->dprime),
	   "Linkage disequilibrium only supported for two point analysis.\n");
  KASSERT ((modelType->type == TP || !modelRange->theta),
	   "Theta specification only for two point analysis.\n");
  KASSERT ((modelType->type == TP || !modelRange->afreq),
	   "Marker allele frequencies only supported for two point analysis.\n");
  KASSERT ((modelType->type != TP || !modelRange->tloc),
	   "Trait loci specification only for multipoint analysis.\n");
  KASSERT ((modelType->type != TP || !modelRange->tlmark),
	   "On-marker trait locus specification only for multipoint analysis.\n");

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
  expandRange (modelRange, modelType);
#if FALSE
  /* Show the partially expanded model. At level 1, following
   * expandRange(), we will have refined the model specification while
   * enforcing all specified constraints that do not involve liability
   * classes. */
  showRange (modelRange, modelType, 1);
#endif

  /* Expand the liability classes, but only if necessary and always
   * honoring inter-class constraints. */
  if (modelRange->nlclass > 1)
    expandClass (modelRange, modelType);

#if FALSE
  /* At level 2, all constraints (including those between classes) are
   * honored, but penet[][][], param[][][][] are not yet fully
   * "factored". */
  showRange (modelRange, modelType, 2);
#endif

  /* Done. Free the constraints; you're done with them. */
  for (i = 0; i < 4; i++)
    if (constraints[i])
      free (constraints[i]);

  /* TODO: fix return expression to yield number of LODS; also for QT
   * and CT models. */
  /* fprintf (stderr, "%d * %d * %d\n", model->ngfreq, model->npenet, model->ntheta); */
  return (TRUE);
}

/**********************************************************************
 * Find the index corresponding to the two-letter code of the type of
 * record we are adding by searching through DIRECTIVES. Returns an
 * integer corresponding to the entry in DIRECTIVES or ERROR if no match
 * is found.
 **********************************************************************/
int
getType (char *line)
{
  char types[25] = DIRECTIVES;
  char *ptr = types;

  while (ptr < types + strlen (types))
    {
      if (strncmp (ptr, line, 2) == 0)
	break;
      ptr += 2;
    }
  /* If no match is found, return ERROR. */
  if (ptr - types >= strlen (types))
    return (ERROR);
  else
    return ((ptr - types) / 2);
}

/**********************************************************************
 * Find the index corresponding to the two-letter relation by
 * searching through OPERATORS. Returns an integer corresponding to
 * the entry in OPERATORS or ERROR if no match is found.
 **********************************************************************/
int
getOperator (char *line)
{
  char types[13] = OPERATORS;
  char *ptr = types;

  while (ptr < types + strlen (types))
    {
      if (strncmp (ptr, line, 2) == 0)
	break;
      ptr += 2;
    }
  /* If no match is found, return ERROR. */
  if (ptr - types >= strlen (types))
    return (ERROR);
  else
    return ((ptr - types) / 2);
}

/**********************************************************************
 * Return the integer found at the beginning of line.
 **********************************************************************/
int
getInteger (char *line)
{
  return ((int) strtol (line, NULL, 10));
}

/**********************************************************************
 * Add one more element to theta vector.  May need to allocate or
 * reallocate memory in order to allow additional room for more
 * values. Maintains values in sorted order for ease of interpretation
 * during debugging.
 **********************************************************************/
void
addTheta (ModelRange * range, int type, double val)
{
  int i;
  double **dtmp;
  int *itmp;

  /* Validate value. */
  KASSERT ((val >= 0 && val <= 0.5), "Bad theta value %g; aborting.\n", val);

  /* Initialize the structure if first access. */
  if (!range->theta)
    {
      if (type == Th)
	{
	  range->ngender = 1;
	  range->theta = malloc (range->ngender * sizeof (double *));
	  range->theta[SEXAV] = malloc (CHUNKSIZE * sizeof (double));
	  thetacnt = malloc (range->ngender * sizeof (int *));
	  thetamax = malloc (range->ngender * sizeof (int *));
	  thetacnt[SEXAV] = thetamax[SEXAV] = 0;
	}
      else
	{
	  range->ngender = 2;
	  range->theta = malloc (range->ngender * sizeof (double *));
	  range->theta[SEXML] = malloc (CHUNKSIZE * sizeof (double));
	  range->theta[SEXFM] = malloc (CHUNKSIZE * sizeof (double));
	  thetacnt = malloc (range->ngender * sizeof (int *));
	  thetamax = malloc (range->ngender * sizeof (int *));
	  thetacnt[SEXML] = thetamax[SEXML] = 0;
	  thetacnt[SEXFM] = thetamax[SEXFM] = 0;
	}
    }

  /* Second, if you are giving sex-specific thetas but have thus far
   * only handled sex-averaged thetas, you need to split and clone the
   * existing sex-averaged thetas. Relies on SEXAV = SEXML indexing. */
  if (type != Th && range->ngender == 1)
    {
      range->ngender = 2;
      dtmp = malloc (range->ngender * sizeof (double *));
      dtmp[0] = range->theta[0];
      free (range->theta);
      range->theta = dtmp;
      itmp = malloc (range->ngender * sizeof (int *));
      itmp[0] = thetacnt[0];
      free (thetacnt);
      thetacnt = itmp;
      itmp = malloc (range->ngender * sizeof (int *));
      itmp[0] = thetamax[0];
      free (thetamax);
      thetacnt = itmp;
      /* Now clone the existing sex-averaged values. */
      range->theta[SEXFM] = malloc (thetamax[SEXAV] * sizeof (double));
      for (i = 0; i < thetacnt[SEXAV]; i++)
	range->theta[SEXFM][i] = range->theta[SEXAV][i];
      thetacnt[SEXFM] = thetacnt[SEXAV];
      thetamax[SEXFM] = thetamax[SEXAV];
    }

  /* Finally, add the thetas as specified. */
  if (type == Th && range->ngender == 1)
    {
      /* Sex-averaged thetas. Only need to update one array. First,
       * enlarge the array if there isn't enough room. */
      if (thetacnt[SEXAV] == thetamax[SEXAV])
	{
	  range->theta[SEXAV] = realloc (range->theta[SEXAV],
					 (thetamax[SEXAV] +
					  CHUNKSIZE) * sizeof (double));
	  thetamax[SEXAV] = thetamax[SEXAV] + CHUNKSIZE;
	}
      /* Add the element. */
      range->theta[SEXAV][thetacnt[SEXAV]] = val;
      thetacnt[SEXAV]++;
    }
  else
    {
      /* Sex-specific thetas. Update one or both arrays separately. */
      if (type == Th || type == Tf)
	{
	  /* Enlarge array if necessary. */
	  if (thetacnt[SEXFM] == thetamax[SEXFM])
	    {
	      range->theta[SEXFM] = realloc (range->theta[SEXFM],
					     (thetamax[SEXFM] +
					      CHUNKSIZE) * sizeof (double));
	      thetamax[SEXFM] = thetamax[SEXFM] + CHUNKSIZE;
	    }
	  /* Add the element. */
	  range->theta[SEXFM][thetacnt[SEXFM]] = val;
	  thetacnt[SEXFM]++;
	}
      /* Next, update the other array. */
      if (type == Th || type == Tm)
	{
	  /* Enlarge array if necessary. */
	  if (thetacnt[SEXML] == thetamax[SEXML])
	    {
	      range->theta[SEXML] = realloc (range->theta[SEXML],
					     (thetamax[SEXML] +
					      CHUNKSIZE) * sizeof (double));
	      thetamax[SEXML] = thetamax[SEXML] + CHUNKSIZE;
	    }
	  /* Add the element. */
	  range->theta[SEXML][thetacnt[SEXML]] = val;
	  thetacnt[SEXML]++;
	}
    }
}

/**********************************************************************
 * Add one more element to appropriate penetrance vector.  May need to
 * allocate or reallocate memory in order to allow additional room for
 * more values. 
 *
 * Recall that penetrances are initially specified as if there are no
 * liability classes. Later, if liability classes are in use, we
 * "expand" the penetrance array to its full three dimensions. So
 * addPenetrance() is only ever used in "preexpansion" mode, where the
 * first dimension of the array (corresponds to liability class) is
 * always 1.
 *
 * The penetrance array: range->penet[nlclass][nallele][nparam]
 **********************************************************************/
void
addPenetrance (ModelRange * range, int type, double val)
{
  int i, j;

  /* Validate value. This is problematic because for DT we know values
   * should be between 0 and 1, but for QT/CT values given here
   * represent the mean of a distribution, and could be just about
   * anything. */
  /* KASSERT ((val >= 0 && val <= 1.0), "Bad penetrance value %g; aborting.\n", val); */

  /* Initialize the structure if first access. */
  if (!range->penet)
    {
      /* Initialize: remember, if you are a pre-expansion penet[][][]
       * array, the first dimension (liability class) will always have
       * a dimension of size 1. */
      range->penet = malloc (sizeof (double **));
      i = NPENET (range->nalleles);
      range->penet[0] = malloc (i * sizeof (double *));
      for (j = 0; j < i; j++)
	range->penet[0][j] = malloc (CHUNKSIZE * sizeof (double));

      /* Remember, liability class is not a true "independent" index
       * in the sense that we are storing combinations across
       * liability classes; hence we only need NPENET() individual
       * values for the counters (indeed, we only need one of each,
       * but that would complicate things too much I suspect). */
      penetcnt = malloc (i * sizeof (int));
      penetmax = malloc (i * sizeof (int));
      for (j = 0; j < i; j++)
	{
	  penetcnt[j] = 0;
	  penetmax[j] = CHUNKSIZE;
	}
    }

  /* Next, add the penetrance to the appropriate location in the
   * penet[][][] array, which may entail allocating more space. Recall
   * type is an integer, where DD=00=0; Dd=01=1; and dd=10=2. We
   * get the type by subtracting DD, the "base type" from the
   * integer count at invocation. */
  if (penetcnt[type] == penetmax[type])
    {
      range->penet[0][type] = realloc (range->penet[0][type],
				       (penetmax[type] +
					CHUNKSIZE) * sizeof (double));
      penetmax[type] = penetmax[type] + CHUNKSIZE;
    }
  /* Add the element. */
  range->penet[0][type][penetcnt[type]] = val;
  penetcnt[type]++;
}

/**********************************************************************
 * Add one more element to gene frequency vector.  May need to
 * allocate or reallocate memory in order to allow additional room for
 * more values.
 **********************************************************************/
void
addGeneFreq (ModelRange * range, double val)
{
  /* Validate value. */
  KASSERT ((val >= 0
	    && val <= 1.0), "Bad gene frequency value %g; aborting.\n", val);

  /* Initialize the structure if first access. */
  if (!range->gfreq)
    range->ngfreq = maxgfreq = 0;
  /* Enlarge array if necessary. */
  if (range->ngfreq == maxgfreq)
    {
      range->gfreq = realloc (range->gfreq,
			      (maxgfreq + CHUNKSIZE) * sizeof (double));
      maxgfreq = maxgfreq + CHUNKSIZE;
    }
  /* Add the element. */
  range->gfreq[range->ngfreq] = val;
  range->ngfreq++;
}

/**********************************************************************
 * Add one more element to alpha vector.  May need to allocate or
 * reallocate memory in order to allow additional room for more
 * values.
 **********************************************************************/
void
addAlpha (ModelRange * range, double val)
{
  /* Validate value. */
  KASSERT ((val >= 0 && val <= 1.0), "Bad alpha value %g; aborting.\n", val);

  /* Initialize the structure if first access. */
  if (!range->alpha)
    range->nalpha = maxalpha = 0;
  /* Enlarge array if necessary. */
  if (range->nalpha == maxalpha)
    {
      range->alpha = realloc (range->alpha,
			      (maxalpha + CHUNKSIZE) * sizeof (double));
      maxalpha = maxalpha + CHUNKSIZE;
    }
  /* Add the element. */
  range->alpha[range->nalpha] = val;
  range->nalpha++;
}

/**********************************************************************
 * Add one more element to trait loci vector.  May need to allocate or
 * reallocate memory in order to allow additional room for more
 * values.
 *
 * The only ugly thing here is that since the TM directive requires
 * calling addTraitLocus() again and again (once for each eventual
 * marker location) we should maintain this array sorted and unique as
 * we go.
 **********************************************************************/
void
addTraitLocus (ModelRange * range, double val)
{
  int i = 0, j;
  /* Validate value. But trait loci may be negative! So no reasonable
   * validation is really possible. */
  /* KASSERT ((val >=0), "Bad trait locus %g; aborting.\n", val); */

  /* Initialize the structure if first access. */
  if (!range->tloc)
    range->ntloc = maxtloc = 0;
  /* Enlarge array if necessary. */
  if (range->ntloc == maxtloc)
    {
      range->tloc = realloc (range->tloc,
			     (maxtloc + CHUNKSIZE) * sizeof (double));
      maxtloc = maxtloc + CHUNKSIZE;
    }

  /* Add the element. First, cue up to where the new element belongs. */
  while (i < range->ntloc && range->tloc[i] < val)
    i++;

  /* Second, if the element is already there, just quit. */
  if (i < range->ntloc && range->tloc[i] == val)
    return;

  /* Third, make room for the new element. */
  for (j = range->ntloc; j > i; j--)
    range->tloc[j] = range->tloc[j - 1];

  /* Fourth, add the element. */
  range->tloc[i] = val;
  range->ntloc++;
}

/**********************************************************************
 * Add one more element to trait threshold vector.  May need to
 * allocate or reallocate memory in order to allow additional room for
 * more values. Note that these are "raw" values; we'll have to expand
 * them by liability and check for any threshold constraints later.
 *
 * Recall that trait thresholds are initially specified as if there
 * are no liability classes. Later, if liability classes are in use,
 * we "expand" the threshold array to its full two dimensions. So
 * addTraitThreshold() is only ever used in "preexpansion" mode, where
 * the first dimension of the array (corresponds to liability class)
 * is always 1.
 **********************************************************************/
void
addTraitThreshold (ModelRange * range, double val)
{
  /* Validate value. The trait threshold should be between the lowest
   * and highest means. But since trait thresholds are subject to
   * liability classes, we won't be able to impose this constraint
   * until later. */
  /* KASSERT ((val >=0), "Bad trait threshold %g; aborting.\n", val); */

  /* Initialize the structure if first access. */
  if (!range->tthresh)
    {
      /* Initialize: remember, if you are a pre-expansion tthresh[][]
       * array, the first dimension (liability class) will always have
       * a dimension of size 1. */
      range->tthresh = malloc (sizeof (double *));
      range->tthresh[0] = malloc (CHUNKSIZE * sizeof (double));
      maxtthresh = CHUNKSIZE;
      range->ntthresh = 0;
    }

  /* Enlarge array if necessary. */
  if (range->ntthresh == maxtthresh)
    {
      range->tthresh[0] = realloc (range->tthresh[0],
				   (maxtthresh +
				    CHUNKSIZE) * sizeof (double));
      maxtthresh = maxtthresh + CHUNKSIZE;
    }
  /* Add the element. */
  range->tthresh[0][range->ntthresh] = val;
  range->ntthresh++;
}

/**********************************************************************
 * Add one more element to dprime vector.  May need to allocate or
 * reallocate memory in order to allow additional room for more
 * values. 
 *
 * Recall that dprime values must be between -1 and 1, where 0
 * corresponds to LE.
 **********************************************************************/
void
addDPrime (ModelRange * range, double val)
{
  /* Validate value. DPrime values must be between -1 and 1. */
  KASSERT ((val >= -1 && val <= 1), "Bad D prime value %g; aborting.\n", val);

  /* Initialize the structure if first access. */
  if (!range->dprime)
    range->ndprime = maxdprime = 0;
  /* Enlarge array if necessary. */
  if (range->ndprime == maxdprime)
    {
      range->dprime = realloc (range->dprime,
			       (maxdprime + CHUNKSIZE) * sizeof (double));
      maxdprime = maxdprime + CHUNKSIZE;
    }
  /* Add the element. */
  range->dprime[range->ndprime] = val;
  range->ndprime++;
}

/**********************************************************************
 * Add one more element to allele frequency vector.  May need to
 * allocate or reallocate memory in order to allow additional room for
 * more values.
 **********************************************************************/
void
addAlleleFreq (ModelRange * range, double val)
{
  /* Validate value. */
  KASSERT ((val >= 0
	    && val <= 1.0), "Bad allele frequency value %g; aborting.\n",
	   val);

  /* Initialize the structure if first access. */
  if (!range->afreq)
    range->nafreq = maxafreq = 0;
  /* Enlarge array if necessary. */
  if (range->nafreq == maxafreq)
    {
      range->afreq = realloc (range->afreq,
			      (maxafreq + CHUNKSIZE) * sizeof (double));
      maxafreq = maxafreq + CHUNKSIZE;
    }
  /* Add the element. */
  range->afreq[range->nafreq] = val;
  range->nafreq++;
}

/**********************************************************************
 * Add one more element to appropriate parameter vector.  May need to
 * allocate or reallocate memory in order to allow additional room for
 * more values.
 *
 * Recall that parameters are only specified in the case of QT or CT
 * traits, and that the number of parameters (npardim) is determined
 * by the distribution type. Moreover, parameters are specified only
 * once and applied uniformly to all allele combinations as if there
 * are no liability classes. Later, if liability classes are in use,
 * we "expand" the parameter array to its full three dimensions. So
 * addParameter() is only ever used in "preexpansion" mode, where the
 * first dimension of the array (corresponds to liability class) is
 * always 0, and the alleles are always 0, too.
 *
 * The parameter array: range->param[nlclass][nallele][pardim][nparam] 
 *
 * Note that paramcnt and parammax are just simple arrays since as we
 * read the config file, parameter values are shared across all allele
 * combinations, and we don't handle liability classes until later:
 * hence only the parameter dimension will be variable.
 **********************************************************************/
void
addParameter (ModelRange * range, int dim, double val)
{
  int i;

  /* First, if this is the first access to the structure, you must
   * initialize it. */
  if (!range->param)
    {
      /* Initialize: remember, if you are a pre-expansion
       * param[][][][] array, the first dimension (liability class)
       * will always have a dimension of size 1, and the second
       * dimension (allele) will also always have a dimension of size
       * 1, since we "share" values, at least before applying
       * constraints, across all allele combinations. */
      range->param = malloc (sizeof (double ***));
      range->param[0] = malloc (sizeof (double **));
      range->param[0][0] = malloc (range->npardim * sizeof (double *));
      for (i = 0; i < range->npardim; i++)
	range->param[0][0][i] = malloc (CHUNKSIZE * sizeof (double));

      /* The only "real" dimension of param[][][][] reflected in
       * paramcnt and parammax is the dim dimension, since that's the
       * only one that will vary at input time. Recall that all
       * parameters are shared across all alleles and liability
       * classes.  */
      paramcnt = malloc (range->npardim * sizeof (int));
      parammax = malloc (range->npardim * sizeof (int));
      for (i = 0; i < range->npardim; i++)
	{
	  paramcnt[i] = 0;
	  parammax[i] = CHUNKSIZE;
	}
    }

  /* Next, add the parameter to the appropriate location in the
   * param[][][][] array, which may entail allocating more space. */
  if (paramcnt[dim] == parammax[dim])
    {
      range->param[0][0][dim] = realloc (range->param[0][0][dim],
					 (parammax[dim] +
					  CHUNKSIZE) * sizeof (double));
      parammax[dim] = parammax[dim] + CHUNKSIZE;
    }
  /* Add the element. */
  range->param[0][0][dim][paramcnt[dim]] = val;
  paramcnt[dim]++;
}

/**************q********************************************************
 * Add one more range to genefrequencies, penetrances, parameter,
 * theta, allelefrequencies, or LD. Parameters are handled by passing
 * a negative type; so -3 correponds to the third parameter.
 **********************************************************************/
void
addRange (ModelRange * range, int type, double lo, double hi, double incr)
{
  double val = lo;

  /* First check to make sure increment moves in the right direction. */
  KASSERT (((hi >= lo && incr > 0) || (hi <= lo && incr < 0)),
	   "Bad range specification from %g to %g by %g; aborting.\n", lo, hi,
	   incr);

  /* Now step through the increment and add individual values as
   * needed. Note that we need to be careful about roundoff issues,
   * like 0, 4, 1 stopping because you get to 4.00001; we'll use
   * ERROR_MARGIN from locus.h as our safety margin, and clean up
   * values that go astray before adding them. */
  while (val <= hi + ERROR_MARGIN)
    {
      /* Clean things up before adding the value. */
      val = MIN (hi, val);

      /* Add value as needed. */
      if (type == Th || type == Tm || type == Tf)
	addTheta (range, type, val);
      else if (type == GF)
	addGeneFreq (range, val);
      else if (type < 0)
	addParameter (range, -type - 1, val);
      else if (type == AL)
	addAlpha (range, val);
      else if (type == TL)
	addTraitLocus (range, val);
      else if (type == TT)
	addTraitThreshold (range, val);
      else if (type == LD)
	addDPrime (range, val);
      else if (type == AF)
	addAlleleFreq (range, val);
      else
	/* addRange() is only used in pre-expansion mode, hence class
	 * is always ERROR. Also, subtract base value of DD from type
	 * to get 0-offset value.*/
	addPenetrance (range, type - DD, val);

      /* Go on to next value. */
      val = val + incr;
    }
}

/**********************************************************************
 * Add a constraint on penetrances or thetas. We're just going to keep
 * these in an array, which we scan when checking penetrance or theta
 * values. Clunky, but OK since it will only be done at setup. A
 * better way would perhaps be to keep theta and penetrance
 * constraints separately.
 **********************************************************************/
void
addConstraint (int type, int a1, int c1, int p1,
	       int op, int a2, int c2, int p2, int disjunct)
{
  char types[25] = DIRECTIVES;

  /* Check for meaningless constraints. TODO: do more of this! */
  KASSERT ((((a1 == Tm && a2 == Tf) || (a1 == Tf && a2 == Tm))
	    || (a1 == TT && a2 == TT) || (a1 >= DD && a2 >= DD && a1 <= dd
					  && a2 <= dd)),
	   "Meaningless constraint %c%c %s %c%c %s; aborting.\n",
	   types[a1 * 2], types[a1 * 2 + 1],
	   ((op == NE) ? "!=" : (op == GE) ? ">=" : ">"), types[a2 * 2],
	   types[a2 * 2 + 1], (disjunct == TRUE) ? "*" : "");

  /* Allocate more space if necessary. */
  if (constmax[type] == constcnt[type])
    {
      constraints[type] = (Constraint *) realloc (constraints[type],
						  (constmax[type] +
						   CHUNKSIZE) *
						  sizeof (Constraint));
      constmax[type] = constmax[type] + CHUNKSIZE;
    }

  /* Mark previous constraint as having a disjunct, if applicable. */
  if (disjunct && constcnt[type] > 0)
    constraints[type][constcnt[type] - 1].alt = TRUE;
  /* Set up current constraint. */
  constraints[type][constcnt[type]].a1 = a1;
  constraints[type][constcnt[type]].op = op;
  constraints[type][constcnt[type]].a2 = a2;
  constraints[type][constcnt[type]].alt = FALSE;

  if (type == CLASSC || type == PARAMCLASSC)
    {
      constraints[type][constcnt[type]].c1 = c1;
      constraints[type][constcnt[type]].c2 = c2;
    }
  if (type == PARAMC || type == PARAMCLASSC)
    {
      constraints[type][constcnt[type]].p1 = p1;
      constraints[type][constcnt[type]].p2 = p2;
    }
  /* Increment constraint count. */
  constcnt[type]++;
}

/**********************************************************************
 * Return TRUE if the ith thetas satisfy all theta constraints, else
 * FALSE. The general idea is to scan through each constraint, exiting
 * whenever you encounter an unsatisfied one (that does not have a
 * disjunct). If there is a disjunct, then step through those until
 * you find one that satisfies the constraint; if you get to the end
 * of the line without satisfying the disjunct, you fail.
 **********************************************************************/
int
checkThetas (ModelRange * range, int i)
{
  int j = 0;
#if FALSE
  char types[15] = DIRECTIVES;
#endif

  /* Scan through each constraint. */
  while (j < constcnt[SIMPLE])
    {
      /* Relies on Tm being the lowest numbered penetrance. */
      if ((constraints[SIMPLE][j].a1 == Tm && constraints[SIMPLE][j].a2 == Tf)
	  || (constraints[SIMPLE][j].a1 == Tf
	      && constraints[SIMPLE][j].a2 == Tm))
	{
	  if ((constraints[SIMPLE][j].op == EQ &&
	       range->theta[constraints[SIMPLE][j].a1 - Tm][i] !=
	       range->theta[constraints[SIMPLE][j].a2 - Tm][i]) ||
	      (constraints[SIMPLE][j].op == NE &&
	       range->theta[constraints[SIMPLE][j].a1 - Tm][i] !=
	       range->theta[constraints[SIMPLE][j].a2 - Tm][i]) ||
	      (constraints[SIMPLE][j].op == GT &&
	       range->theta[constraints[SIMPLE][j].a1 - Tm][i] >
	       range->theta[constraints[SIMPLE][j].a2 - Tm][i]) ||
	      (constraints[SIMPLE][j].op == GE &&
	       range->theta[constraints[SIMPLE][j].a1 - Tm][i] >=
	       range->theta[constraints[SIMPLE][j].a2 - Tm][i]))
	    {
	      /* Satisfied; skip other disjuncts. */
	      while (constraints[SIMPLE][j].alt == TRUE)
		j++;
	    }
	  else if (constraints[SIMPLE][j].alt == FALSE)
	    {
#if FALSE
	      fprintf (stderr, "%c%c %s %c%c %s => %g %g %g\n",
		       types[(constraints[SIMPLE][j].a1) * 2],
		       types[(constraints[SIMPLE][j].a1) * 2 + 1],
		       ((constraints[SIMPLE][j].op ==
			 EQ) ? "==" : (constraints[SIMPLE][j].op ==
				       NE) ? "!=" : (constraints[SIMPLE][j].
						     op == GE) ? ">=" : ">"),
		       types[(constraints[SIMPLE][j].a2) * 2],
		       types[(constraints[SIMPLE][j].a2) * 2 + 1],
		       (constraints[SIMPLE][j].alt == TRUE) ? "*" : "",
		       range->penet[0][i], range->penet[1][i],
		       range->penet[2][i]);
#endif
	      return (FALSE);
	    }
	}
      j++;
    }
  return (TRUE);
}

/**********************************************************************
 * Return TRUE if the ith penetrances satisfy all penetrance
 * constraints, else FALSE. The general idea is to scan through each
 * constraint, exiting whenever you encounter an unsatisfied one (that
 * does not have a disjunct). If there is a disjunct, then step
 * through those until you find one that satisfies the constraint; if
 * you get to the end of the line without satisfying the disjunct, you
 * fail.
 **********************************************************************/
int
checkPenets (ModelRange * range, int i)
{
  int j = 0;
#if FALSE
  char types[15] = DIRECTIVES;
#endif

  /* Scan through each constraint. */
  while (j < constcnt[SIMPLE])
    {
      if (constraints[SIMPLE][j].a1 == TT || constraints[SIMPLE][j].a2 == TT)
	{
	  /* Skip threshold constraints; they are irrelevant for penetrances. */
	  j++;
	  continue;
	}
      else if (constraints[SIMPLE][j].a1 >= DD
	       && constraints[SIMPLE][j].a2 >= DD)
	{
	  /* Relies on DD being the lowest numbered penetrance. */
	  if ((constraints[SIMPLE][j].op == EQ &&
	       range->penet[0][constraints[SIMPLE][j].a1 - DD][i] ==
	       range->penet[0][constraints[SIMPLE][j].a2 - DD][i]) ||
	      (constraints[SIMPLE][j].op == NE &&
	       range->penet[0][constraints[SIMPLE][j].a1 - DD][i] !=
	       range->penet[0][constraints[SIMPLE][j].a2 - DD][i]) ||
	      (constraints[SIMPLE][j].op == GT &&
	       range->penet[0][constraints[SIMPLE][j].a1 - DD][i] >
	       range->penet[0][constraints[SIMPLE][j].a2 - DD][i]) ||
	      (constraints[SIMPLE][j].op == GE &&
	       range->penet[0][constraints[SIMPLE][j].a1 - DD][i] >=
	       range->penet[0][constraints[SIMPLE][j].a2 - DD][i]))
	    {
	      /* Satisfied; skip other disjuncts. */
	      while (constraints[SIMPLE][j].alt == TRUE)
		j++;
	    }
	  else if (constraints[SIMPLE][j].alt == FALSE)
	    {
#if FALSE
	      fprintf (stderr, "%c%c %s %c%c %s => %g %g %g\n",
		       types[(constraints[SIMPLE][j].a1) * 2],
		       types[(constraints[SIMPLE][j].a1) * 2 + 1],
		       ((constraints[SIMPLE][j].op ==
			 EQ) ? "==" : (constraints[SIMPLE][j].op ==
				       NE) ? "!=" : (constraints[SIMPLE][j].
						     op == GE) ? ">=" : ">"),
		       types[(constraints[SIMPLE][j].a2) * 2],
		       types[(constraints[SIMPLE][j].a2) * 2 + 1],
		       (constraints[SIMPLE][j].alt == TRUE) ? "*" : "",
		       range->penet[0][0][i], range->penet[0][1][i],
		       range->penet[0][2][i]);
#endif
	      return (FALSE);
	    }
	}
      j++;
    }
  return (TRUE);
}

/**********************************************************************
 * Return TRUE if the ith penetrances satisfy all inter-liability
 * class penetrance constraints, else FALSE. The general idea is to
 * scan through each constraint, exiting whenever you encounter an
 * unsatisfied one (that does not have a disjunct). If there is a
 * disjunct, then step through those until you find one that satisfies
 * the constraint; if you get to the end of the line without
 * satisfying the disjunct, you fail.
 *
 * This is somewhat ugly: while it does reduce the number of overall
 * legal combinations, some extra LODs will still be written to disk
 * in the case that a pedigree does not contain an individual from
 * each liability class. This case will have to be detected and dealt
 * with dynamically?
 **********************************************************************/
int
checkClassPenets (ModelRange * range, int i)
{
  int j = 0;
#if FALSE
  char types[15] = DIRECTIVES;
#endif

  /* Scan through each constraint: inter-class constraints are on
   * penetrances or trait thresholds, not, e.g., gene frequencies or
   * thetas. */
  while (j < constcnt[CLASSC])
    {
      if (constraints[CLASSC][j].a1 == TT || constraints[CLASSC][j].a2 == TT)
	{
	  /* Skip threshold constraints; they are irrelevant for penetrances. */
	  j++;
	  continue;
	}
      else if ((constraints[CLASSC][j].op == EQ &&
		range->penet[constraints[CLASSC][j].c1 -
			     1][constraints[CLASSC][j].a1 - DD][i] ==
		range->penet[constraints[CLASSC][j].c2 -
			     1][constraints[CLASSC][j].a2 - DD][i])
	       || (constraints[CLASSC][j].op == NE
		   && range->penet[constraints[CLASSC][j].c1 -
				   1][constraints[CLASSC][j].a1 - DD][i] !=
		   range->penet[constraints[CLASSC][j].c2 -
				1][constraints[CLASSC][j].a2 - DD][i])
	       || (constraints[CLASSC][j].op == GT
		   && range->penet[constraints[CLASSC][j].c1 -
				   1][constraints[CLASSC][j].a1 - DD][i] >
		   range->penet[constraints[CLASSC][j].c2 -
				1][constraints[CLASSC][j].a2 - DD][i])
	       || (constraints[CLASSC][j].op == GE
		   && range->penet[constraints[CLASSC][j].c1 -
				   1][constraints[CLASSC][j].a1 - DD][i] >=
		   range->penet[constraints[CLASSC][j].c2 -
				1][constraints[CLASSC][j].a2 - DD][i]))
	{
	  /* Satisfied; skip other disjuncts. */
	  while (constraints[CLASSC][j].alt == TRUE)
	    j++;
	}
      else if (constraints[CLASSC][j].alt == FALSE)
	{
#if FALSE
	  fprintf (stderr, "!%c%c %d %s %c%c %d %s => %g %s %g\n",
		   types[(constraints[CLASSC][j].a1) * 2],
		   types[(constraints[CLASSC][j].a1) * 2 + 1],
		   constraints[CLASSC][i].c1,
		   ((constraints[CLASSC][j].op ==
		     EQ) ? "==" : (constraints[CLASSC][j].op ==
				   NE) ? "!=" : (constraints[CLASSC][j].op ==
						 GE) ? ">=" : ">"),
		   types[(constraints[CLASSC][j].a2) * 2],
		   types[(constraints[CLASSC][j].a2) * 2 + 1],
		   constraints[CLASSC][i].c2,
		   (constraints[CLASSC][j].alt == TRUE) ? "*" : "",
		   (range->
		    penet[constraints[CLASSC][j].c1 -
			  1][constraints[CLASSC][j].a1 - DD][i]),
		   ((constraints[CLASSC][j].op ==
		     EQ) ? "==" : (constraints[CLASSC][j].op ==
				   NE) ? "!=" : (constraints[CLASSC][j].op ==
						 GE) ? ">=" : ">"),
		   (range->
		    penet[constraints[CLASSC][j].c2 -
			  1][constraints[CLASSC][j].a2 - DD][i]));
#endif
	  return (FALSE);
	}
      j++;
    }
  return (TRUE);
}

/**********************************************************************
 * Return TRUE if the ith parameters satisfy all parameter
 * constraints, else FALSE. The general idea is to scan through each
 * constraint, exiting whenever you encounter an unsatisfied one (that
 * does not have a disjunct). If there is a disjunct, then step
 * through those until you find one that satisfies the constraint; if
 * you get to the end of the line without satisfying the disjunct, you
 * fail.
 **********************************************************************/
int
checkParams (ModelRange * range, int i)
{
  int j = 0;
#if FALSE
  char types[15] = DIRECTIVES;
#endif

  /* Scan through each constraint. */
  while (j < constcnt[PARAMC])
    {
#if FALSE
      fprintf (stderr, "%d: P%d %c%c %s P%d %c%c %s => ", i,
	       constraints[PARAMC][j].a1,
	       types[(constraints[PARAMC][j].a1) * 2],
	       types[(constraints[PARAMC][j].a1) * 2 + 1],
	       ((constraints[PARAMC][j].op ==
		 EQ) ? "==" : (constraints[PARAMC][j].op ==
			       NE) ? "!=" : (constraints[PARAMC][j].op ==
					     GE) ? ">=" : ">"),
	       constraints[PARAMC][j].a2,
	       types[(constraints[PARAMC][j].a2) * 2],
	       types[(constraints[PARAMC][j].a2) * 2 + 1],
	       (constraints[PARAMC][j].alt == TRUE) ? "*" : "");
#endif
      /* Relies on DD being the lowest numbered penetrance. */
      if (constraints[PARAMC][j].a1 >= DD && constraints[PARAMC][j].a2 >= DD)
	{
	  if ((constraints[PARAMC][j].op == EQ &&
	       range->param[0][constraints[PARAMC][j].a1 -
			       DD][constraints[PARAMC][j].p1 - 1][i] ==
	       range->param[0][constraints[PARAMC][j].a2 -
			       DD][constraints[PARAMC][j].p2 - 1][i])
	      || (constraints[PARAMC][j].op == NE
		  && range->param[0][constraints[PARAMC][j].a1 -
				     DD][constraints[PARAMC][j].p1 - 1][i] !=
		  range->param[0][constraints[PARAMC][j].a2 -
				  DD][constraints[PARAMC][j].p2 - 1][i])
	      || (constraints[PARAMC][j].op == GT
		  && range->param[0][constraints[PARAMC][j].a1 -
				     DD][constraints[PARAMC][j].p1 - 1][i] >
		  range->param[0][constraints[PARAMC][j].a2 -
				  DD][constraints[PARAMC][j].p2 - 1][i])
	      || (constraints[PARAMC][j].op == GE
		  && range->param[0][constraints[PARAMC][j].a1 -
				     DD][constraints[PARAMC][j].p1 - 1][i] >=
		  range->param[0][constraints[PARAMC][j].a2 -
				  DD][constraints[PARAMC][j].p2 - 1][i]))
	    {
	      /* Satisfied; skip other disjuncts. */
#if FALSE
	      fprintf (stderr, "satisfied\n");
#endif
	      while (constraints[PARAMC][j].alt == TRUE)
		j++;
	    }
	  else if (constraints[PARAMC][j].alt == FALSE)
	    {
#if FALSE
	      fprintf (stderr, "fail\n");
#endif
	      return (FALSE);
	    }
	}
      j++;
    }
  return (TRUE);
}

/**********************************************************************
 * Return TRUE if the ith parameters satisfy all inter-liability class
 * parameter constraints, else FALSE. The general idea is to scan
 * through each constraint, exiting whenever you encounter an
 * unsatisfied one (that does not have a disjunct). If there is a
 * disjunct, then step through those until you find one that satisfies
 * the constraint; if you get to the end of the line without
 * satisfying the disjunct, you fail.
 *
 * This is somewhat ugly: while it does reduce the number of overall
 * legal combinations, some extra LODs will still be written to disk
 * in the case that a pedigree does not contain an individual from
 * each liability class. This case will have to be detected and dealt
 * with dynamically?
 **********************************************************************/
int
checkClassParams (ModelRange * range, int i)
{
  int j = 0;
#if FLASE
  char types[15] = DIRECTIVES;
#endif

  /* Scan through each constraint. */
  while (j < constcnt[PARAMCLASSC])
    {
      /* Relies on DD being the lowest numbered penetrance. */
      if (constraints[PARAMCLASSC][j].a1 >= DD
	  && constraints[PARAMCLASSC][j].a2 >= DD)
	{
	  if ((constraints[PARAMCLASSC][j].op == EQ &&
	       range->param[constraints[PARAMCLASSC][j].c1 -
			    1][constraints[PARAMCLASSC][j].a1 -
			       DD][constraints[PARAMCLASSC][j].p1 - 1][i] ==
	       range->param[constraints[PARAMCLASSC][j].c2 -
			    1][constraints[PARAMCLASSC][j].a2 -
			       DD][constraints[PARAMCLASSC][j].p2 - 1][i])
	      || (constraints[PARAMCLASSC][j].op == NE
		  && range->param[constraints[PARAMCLASSC][j].c1 -
				  1][constraints[PARAMCLASSC][j].a1 -
				     DD][constraints[PARAMCLASSC][j].p1 -
					 1][i] !=
		  range->param[constraints[PARAMCLASSC][j].c2 -
			       1][constraints[PARAMCLASSC][j].a2 -
				  DD][constraints[PARAMCLASSC][j].p2 - 1][i])
	      || (constraints[PARAMCLASSC][j].op == GT
		  && range->param[constraints[PARAMCLASSC][j].c1 -
				  1][constraints[PARAMCLASSC][j].a1 -
				     DD][constraints[PARAMCLASSC][j].p1 -
					 1][i] >
		  range->param[constraints[PARAMCLASSC][j].c2 -
			       1][constraints[PARAMCLASSC][j].a2 -
				  DD][constraints[PARAMCLASSC][j].p2 - 1][i])
	      || (constraints[PARAMCLASSC][j].op == GE
		  && range->param[constraints[PARAMCLASSC][j].c1 -
				  1][constraints[PARAMCLASSC][j].a1 -
				     DD][constraints[PARAMCLASSC][j].p1 -
					 1][i] >=
		  range->param[constraints[PARAMCLASSC][j].c2 -
			       1][constraints[PARAMCLASSC][j].a2 -
				  DD][constraints[PARAMCLASSC][j].p2 - 1][i]))
	    {
	      /* Satisfied; skip other disjuncts. */
#if FALSE
	      fprintf (stderr, "satisfied\n");
#endif
	      while (constraints[PARAMCLASSC][j].alt == TRUE)
		j++;
	    }
	  else if (constraints[PARAMCLASSC][j].alt == FALSE)
	    {
#if FALSE
	      fprintf (stderr, "fail\n");
#endif
	      return (FALSE);
	    }
	}
      j++;
    }
  return (TRUE);
}

/**********************************************************************
 * Return TRUE if the ith threshold satisfy all inter-liability class
 * threshold constraints, else FALSE. The general idea is to scan
 * through each constraint, exiting whenever you encounter an
 * unsatisfied one (that does not have a disjunct). If there is a
 * disjunct, then step through those until you find one that satisfies
 * the constraint; if you get to the end of the line without
 * satisfying the disjunct, you fail.
 *
 * This is somewhat ugly: while it does reduce the number of overall
 * legal combinations, some extra LODs will still be written to disk
 * in the case that a pedigree does not contain an individual from
 * each liability class. This case will have to be detected and dealt
 * with dynamically?
 *
 * Note: there is no corresponding checkThreshold () function, because
 * the only possible constraints you can impose on trait thresholds
 * are inter liability class constraints.
 *
 * Note: We should also be performing "external" (that is, implicit
 * constraints) validation here, since we want these thresholds to be
 * within the min and max values of the mean for a given class.
 **********************************************************************/
int
checkClassThreshold (ModelRange * range, int i)
{
  int j = 0;
#if FALSE
  char types[15] = DIRECTIVES;
#endif

  /* Scan through each constraint: inter-class constraints are on
   * penetrances or trait thresholds, not, e.g., gene frequencies or
   * thetas. */
  while (j < constcnt[CLASSC])
    {
      if (constraints[CLASSC][j].a1 != TT || constraints[CLASSC][j].a2 != TT)
	{
	  /* Skip non threshold constraints. */
	  j++;
	  continue;
	}
      else if ((constraints[CLASSC][j].op == EQ &&
		range->tthresh[constraints[CLASSC][j].c1 - 1][i] ==
		range->tthresh[constraints[CLASSC][j].c2 - 1][i]) ||
	       (constraints[CLASSC][j].op == NE &&
		range->tthresh[constraints[CLASSC][j].c1 - 1][i] !=
		range->tthresh[constraints[CLASSC][j].c2 - 1][i]) ||
	       (constraints[CLASSC][j].op == GT &&
		range->tthresh[constraints[CLASSC][j].c1 - 1][i] >
		range->tthresh[constraints[CLASSC][j].c2 - 1][i]) ||
	       (constraints[CLASSC][j].op == GE &&
		range->tthresh[constraints[CLASSC][j].c1 - 1][i] >=
		range->tthresh[constraints[CLASSC][j].c2 - 1][i]))
	{
	  /* Satisfied; skip other disjuncts. */
	  while (constraints[CLASSC][j].alt == TRUE)
	    j++;
	}
      else if (constraints[CLASSC][j].alt == FALSE)
	{
#if FALSE
	  fprintf (stderr, "!%c%c %d %s %c%c %d %s => %g %s %g\n",
		   types[(constraints[CLASSC][j].a1) * 2],
		   types[(constraints[CLASSC][j].a1) * 2 + 1],
		   constraints[CLASSC][i].c1,
		   ((constraints[CLASSC][j].op ==
		     EQ) ? "==" : (constraints[CLASSC][j].op ==
				   NE) ? "!=" : (constraints[CLASSC][j].op ==
						 GE) ? ">=" : ">"),
		   types[(constraints[CLASSC][j].a2) * 2],
		   types[(constraints[CLASSC][j].a2) * 2 + 1],
		   constraints[CLASSC][i].c2,
		   (constraints[CLASSC][j].alt == TRUE) ? "*" : "",
		   ((constraints[CLASSC][j].type !=
		     TT) ? (range->penet[constraints[CLASSC][j].c1 -
					 1][constraints[CLASSC][j].a1 -
					    DD][i]) : (range->
						       tthresh[constraints
							       [CLASSC][j].
							       a1][i])) ((constraints[CLASSC][j].op == EQ) ? "==" : (constraints[CLASSC][j].op == NE) ? "!=" : (constraints[CLASSC][j].op == GE) ? ">=" : ">"), ((constraints[CLASSC][j].type != TT) ? (range->penet[constraints[CLASSC][j].c2 - 1][constraints[CLASSC][j].a2 - DD][i]) : (range->tthresh[constraints[CLASSC][j].a2][i])));
#endif
	  return (FALSE);
	}
      j++;
    }
  return (TRUE);
}

/**********************************************************************
 * Sort the model values so that we can apply constraints more
 * easily. Note: sortRange() is only applied prior to class expansion,
 * so we needn't worry, ever, about non-zero liability classes in
 * penet[][][]. Similarly, we needn't worry about liability classes or
 * allele counts in param[][][].
 *
 * Note that range->tloc is handled differently, because it may be
 * called back (in the event of a TM parameter) and so sorting and
 * uniquifying must be done on the fly for this array only.
 **********************************************************************/
void
sortRange (ModelRange * range)
{
  int i;
  quicksort (range->gfreq, 0, range->ngfreq);
  if (range->penet[0])
    for (i = 0; i < NPENET (range->nalleles); i++)
      quicksort (range->penet[0][i], 0, penetcnt[i]);
  if (range->param)
    for (i = 0; i < range->npardim; i++)
      quicksort (range->param[0][0][i], 0, paramcnt[i]);
  if (range->theta)
    for (i = 0; i < range->ngender; i++)
      quicksort (range->theta[i], 0, thetacnt[i]);
  /* if (range->tloc)
     quicksort (range->tloc, 0, range->ntloc); */
  if (range->tthresh)
    quicksort (range->tthresh[0], 0, range->ntthresh);
  if (range->alpha)
    quicksort (range->alpha, 0, range->nalpha);
}

/**********************************************************************
 * Quicksort routines, used to sort values (both doubles and floats)
 * in the model specification.
 **********************************************************************/
#define RANDINT(n) (random() % (n))

/* Swap two doubles in array[] */
inline void
swap (double *array, int i, int j)
{
  double temp = array[i];
  array[i] = array[j];
  array[j] = temp;
}

/* Recursive quicksort for array of doubles. */
void
quicksort (double *array, int lo, int hi)
{
  int mid, i;
  double pivot;

  if (lo >= hi - 1)
    return;

  i = RANDINT ((hi - lo)) + lo;
  swap (array, lo, i);
  pivot = array[lo];

  mid = lo;
  for (i = lo + 1; i < hi; i++)
    if (array[i] < pivot)
      {
	mid++;
	swap (array, mid, i);
      }
  swap (array, lo, mid);

  quicksort (array, lo, mid);
  quicksort (array, mid + 1, hi);
}

/**********************************************************************
 * Remove duplicate values from the model for efficiency's sake.  Both
 * double and float versions are given.
 **********************************************************************/
/* Remove exact duplicates in an array. */
inline int
uniquify (double *array, int len)
{
  int i = 0, j;
  while (i < len - 1)
    {
      if (array[i] == array[i + 1])
	{
	  for (j = i + 1; j < len - 1; j++)
	    array[j] = array[j + 1];
	  len--;
	}
      else
	i++;
    }
  return (len);
}

/* Uniquify the model's arrays. This is only done prior to expanding
 * the ranges, so we needn't worry, ever, about non-zero liability
 * classes in penet[][][]. Similarly, we needn't worry about liability
 * classes or allele counts in param[][][]. 
 *
 * Note that range->tloc is handled differently, because it may be
 * called back (in the event of a TM parameter) and so sorting and
 * uniquifying must be done on the fly for this array only. */
void
uniqRange (ModelRange * range)
{
  int i;
  range->ngfreq = uniquify (range->gfreq, range->ngfreq);
  for (i = 0; i < NPENET (range->nalleles); i++)
    penetcnt[i] = uniquify (range->penet[0][i], penetcnt[i]);
  if (range->param)
    for (i = 0; i < range->npardim; i++)
      paramcnt[i] = uniquify (range->param[0][0][i], paramcnt[i]);
  if (range->theta)
    for (i = 0; i < range->ngender; i++)
      thetacnt[i] = uniquify (range->theta[i], thetacnt[i]);
  /* if (range->tloc)
     range->ntloc = uniquify (range->tloc, range->ntloc); */
  if (range->tthresh)
    range->ntthresh = uniquify (range->tthresh[0], range->ntthresh);
  if (range->alpha)
    range->nalpha = uniquify (range->alpha, range->nalpha);
}

/**********************************************************************
 * Fully expand the penetrance (by DD, Dd, dd) and theta values (by
 * Tm, Tf) in the model while honoring all constraints. Using the
 * fully expanded version allows us to accurately determine the size
 * of the space a priori, rather than using the constraints to censor
 * values as we go. The expansion helps us avoid allocating excess
 * space on disk.
 *
 * If you are doing QT or CT, you'll need to expand the parameters as
 * well. Expanded parameters will also increase the size of the
 * penetrance arrays.
 *
 * Note that liability class constraints are handled later.  So you
 * won't need to expand trait threshold values for CT, since these
 * only "factor" by liability class.
 **********************************************************************/
void
expandRange (ModelRange * range, ModelType * type)
{
  int i, j, k, l, m;
  double ****tmp4;
  double ***tmp3;
  double **tmp2;

  /* Start with thetas, but only if DT. QT/CT use trait loci, which
   * also may need expansion, but we'll need to wait until later (when
   * the markers are all read in) to do so.
   *
   * Since thetas are not subject to liability class constraints, we
   * can finalize the thetas here by fully "factoring" male and female
   * thetas, if available, recording the number of combinations in
   * range->ntheta, and freeing thetacnt and thetamax.
   *
   * Keep the current theta array values, and create a new theta
   * array. We'll allocate enough space for all combinations, even
   * though, in the end, constraint application may reduce this
   * number. The actualy number of combinations will be stored in
   * range->nthetas. */
  if ((tmp2 = range->theta))
    {
      range->theta = malloc (range->ngender * sizeof (double *));
      if (range->ngender == 1)
	{
	  /* Easy case is that you have only one gender anyway. Since
	   * there is only one dimension, there can be no constraints
	   * worth enforcing. */
	  range->theta[0] = tmp2[0];
	  range->ntheta = thetacnt[SEXAV];
	}
      else
	{
	  /* If you have 2 genders, you need to factor their respective
	   * values while checking constraints. */
	  for (i = 0; i < range->ngender; i++)
	    range->theta[i] =
	      malloc ((thetacnt[SEXML] * thetacnt[SEXFM]) * sizeof (double));
	  range->ntheta = 0;
	  for (i = 0; i < thetacnt[SEXML]; i++)
	    for (j = 0; j < thetacnt[SEXFM]; j++)
	      {
		range->theta[SEXML][range->ntheta] = tmp2[SEXML][i];
		range->theta[SEXFM][range->ntheta] = tmp2[SEXFM][j];

		if (checkThetas (range, range->ntheta))
		  range->ntheta++;
	      }
	  /* Free old copies of range->theta[i]. */
	  for (i = 0; i < range->ngender; i++)
	    free (tmp2[i]);
	}
      /* Free old copy of range->theta. */
      free (tmp2);
      free (thetacnt);
      free (thetamax);
    }

  /* Next, we expand the penetrances. The goal here is just to get the
   * combinations right, while ignoring the liability classes. It's a
   * stickier problem than the thetas, mostly because these values
   * will have to be expanded yet again when dealing with liability
   * classes, but also because the code must work for multiallelic
   * diseases (where there are more than 3 allele combinations, and,
   * therefore, more than 3 rows in the penetrance array).
   *
   * Stash pointers to penet, penet[0], and penetcnt away
   * temporarily. We'll need these to free things appropriately. */
  tmp3 = range->penet;
  tmp2 = range->penet[0];

  /* Set up the penetrance array, ignoring the first dimension for
   * now. We know there will be PROD_i(penetcnt[i])^#allcombo values
   * in the expanded array, unless the constraints rule some out. 
   *
   * We'll keep track of the resulting number of combinations in
   * range->npenet. */
  range->npenet = 0;
  i = 1;
  for (j = 0; j < NPENET (range->nalleles); j++)
    i = i * penetcnt[j];
  range->penet = malloc (sizeof (double **));
  range->penet[0] = malloc (NPENET (range->nalleles) * sizeof (double *));
  for (j = 0; j < NPENET (range->nalleles); j++)
    range->penet[0][j] = malloc (i * sizeof (double));

  /* OK, now populate the array. */
  for (k = 0; k < i; k++)
    {
      l = 1;
      for (m = 0; m < NPENET (range->nalleles); m++)
	{
	  range->penet[0][m][range->npenet] =
	    tmp2[m][((int) (k / l)) % penetcnt[m]];
	  l = l * penetcnt[m];
	}
      /* Check the constraints. */
      if (checkPenets (range, range->npenet))
	range->npenet++;
    }

  /* OK, we're done with penetcnt and penetmax, since we're using
   * range->npenet to keep track of the number of combinations. We
   * can also ditch the stashed values for penet[][][]. */
  for (i = 0; i < NPENET (range->nalleles); i++)
    free (tmp2[i]);
  free (tmp2);
  free (tmp3);
  free (penetmax);
  free (penetcnt);

  /* Finally, expand the parameters for QT or CT. Recall param[][][][]
   * at this point only has 2 real indexes: dimension and the number
   * of values. The goal here is to factor these into 3 real indexes,
   * including allele combinations, while respecting the
   * constraints. We'll worry about liability classes later.
   *
   * Note that the resulting param[0][][][] array will be of uniform
   * dimension, that being ((#value)^#dim)^#allcombos). */
  if (range->param)
    {
      /* Stash the original values somewhere, to be freed
       * appropriately. */
      tmp4 = range->param;
      tmp3 = range->param[0];
      tmp2 = range->param[0][0];

      /* Since the resulting parameter array will be uniform, we can
       * allocate, in advance, the appropriate amount of memory. This
       * will be large! Hopefully, the constraints will limit how many
       * elements are actually in use, which will be stored in
       * range->nparam. */

      /* Let i be the max number of entries you might see for each
       * dimension, then raise it to the power of the dimensions. */
      i = pow (paramcnt[0], NPENET (range->nalleles));
      i = pow (i, range->npardim);

      /* Since the parameter vector is to be replicated across all of
       * the dimensions, they will all initially have the same number
       * of parameter values, although this will surely change as the
       * constraints are applied. Go ahead and allocate what you might
       * need. */
      range->param = malloc (sizeof (double ***));
      range->param[0] =
	malloc (NPENET (range->nalleles) * sizeof (double **));
      for (j = 0; j < NPENET (range->nalleles); j++)
	{
	  range->param[0][j] = malloc (range->npardim * sizeof (double *));
	  /* Corrected by Yungui: used to read:
	   *  for (k = 0; k < paramcnt[0]; k++) */
	  for (k = 0; k < range->npardim; k++)
	    range->param[0][j][k] = malloc (i * sizeof (double));
	}

      /* Time to populate the range->param array. Recall the first
       * index will always be 0 since we're not yet dealing with
       * liability classes. We'll use range->nparam to store the
       * number of valid combinations. */
      range->nparam = 0;
      for (j = 0; j < i; j++)
	{
	  for (k = 0; k < NPENET (range->nalleles); k++)
	    for (l = 0; l < range->npardim; l++)
	      range->param[0][k][l][range->nparam] =
		tmp2[l][((int)
			 (j /
			  (pow (paramcnt[l], (l + k * range->npardim))))) %
			paramcnt[l]];
	  /* Check the constraints. */
	  if (checkParams (range, range->nparam))
	    range->nparam++;
	}

      /* Done. Free copies of original parameters you'd stashed away,
       * including the paramcnt and parammax arrays; you won't be
       * needing them again, since the (uniform) number of entries
       * stored is given by range->nparam. */
      free (tmp2);
      free (tmp3);
      free (tmp4);
      free (paramcnt);
      free (parammax);
    }
}

/**********************************************************************
 * Fully expand the threshold, penetrance and parameter values by
 * liability class while honoring any inter-class constraints. This
 * only gets called if we are using liability classes (i.e.,
 * modelRange->nlclass is greater than 1).
 **********************************************************************/
void
expandClass (ModelRange * range, ModelType * type)
{
  int i, j, k, l, m, n, o;
  double *tmp1;
  double **tmp2;
  double ***tmp3;
  double ****tmp4;

  /* Threshold expansion by liability class. */
  if (range->tthresh)
    {
      /* Stash pointers so you can free properly. Recall
       * range->tthresh[0] is the only dimension originally allocated;
       * we'll use the values stored there during expansion, producing
       * a new multidimensional range->tthresh array. */
      tmp2 = range->tthresh;
      tmp1 = range->tthresh[0];
      i = range->ntthresh;

      /* Barring constraints, how many values might you have? */
      j = pow (i, range->nlclass);

      /* Allocate a new tthresh array structure. */
      range->tthresh = malloc (range->nlclass * sizeof (double *));
      for (k = 0; k < range->nlclass; k++)
	range->tthresh[k] = malloc (j * sizeof (double));

      /* OK, now populate the array. */
      range->ntthresh = 0;
      for (k = 0; k < j; k++)
	{
	  l = 1;
	  for (m = 0; m < range->nlclass; m++)
	    {
	      range->tthresh[m][range->ntthresh] = tmp1[((int) (k / l)) % i];
	      l = l * i;
	    }
	  /* Check the class constraints. */
	  if (checkClassThreshold (range, range->ntthresh))
	    range->ntthresh++;
	}
      /* Done. Free up the old copy of the array. */
      free (tmp1);
      free (tmp2);
    }

  /* Penetrance expansion by liability class. Here, we want to factor
   * the existing penetrances over multiple liability classes.
   *
   * Stash pointers to the current penet array and its size. Recall
   * range->penet[0] is the only dimension originally allocated; we'll
   * use the values stored there during expansion, producing a new
   * multidimensional range->penet array. Be sure to free everything
   * properly when done. */
  tmp3 = range->penet;
  tmp2 = range->penet[0];
  i = range->npenet;

  /* Barring constraints, how many values might you have? */
  j = pow (range->npenet, range->nlclass);

  /* Allocate a new penet array structure. */
  range->penet = malloc (range->nlclass * sizeof (double **));
  for (k = 0; k < range->nlclass; k++)
    {
      range->penet[k] = malloc (NPENET (range->nalleles) * sizeof (double *));
      for (l = 0; l < NPENET (range->nalleles); l++)
	range->penet[k][l] = malloc (j * sizeof (double));
    }

  /* OK, now populate the array. */
  range->npenet = 0;
  for (k = 0; k < j; k++)
    {
      l = 1;
      for (m = 0; m < range->nlclass; m++)
	{
	  for (n = 0; n < NPENET (range->nalleles); n++)
	    range->penet[m][n][range->npenet] =
	      tmp2[n][(((int) (k / l)) % i)];
	  l = l * i;
	}
      /* Check the class constraints. */
      if (checkClassPenets (range, range->npenet))
	range->npenet++;
    }
  /* Done. Free up the old copy of the range->penet array. */
  for (i = 0; i < NPENET (range->nalleles); i++)
    free (tmp2[i]);
  free (tmp2);
  free (tmp3);

  /* Ready to work on the param array, if necessary.
   *
   * Again, stash a copy of the current param array and its size. */
  if (range->param)
    {
      tmp4 = range->param;
      tmp3 = range->param[0];
      i = range->nparam;

      /* Barring constraints, how many values might you have? */
      j = pow (range->nparam, range->nlclass);

      /* Allocate a new param array structure. */
      range->param = malloc (range->nlclass * sizeof (double ***));
      for (k = 0; k < range->nlclass; k++)
	{
	  range->param[k] =
	    malloc (NPENET (range->nalleles) * sizeof (double **));
	  for (l = 0; l < NPENET (range->nalleles); l++)
	    {
	      range->param[k][l] =
		malloc (range->npardim * sizeof (double *));
	      for (m = 0; m < range->npardim; m++)
		range->param[k][l][m] = malloc (j * sizeof (double));
	    }
	}

      /* OK, now populate the array. */
      range->nparam = 0;
      for (k = 0; k < j; k++)
	{
	  l = 1;
	  for (m = 0; m < range->nlclass; m++)
	    {
	      for (n = 0; n < NPENET (range->nalleles); n++)
		for (o = 0; o < range->npardim; o++)
		  range->param[m][n][o][range->nparam] =
		    tmp3[n][o][((int) (k / l)) % i];
	      l = l * i;
	    }

	  /* Check the class constraints. */
	  if (checkClassParams (range, range->nparam))
	    range->nparam++;
	}

      /* Done. Free up the old copy of the array. */
      for (i = 0; i < NPENET (range->nalleles); i++)
	{
	  for (j = 0; j < range->npardim; j++)
	    free (tmp3[i][j]);
	  free (tmp3[i]);
	}
      free (tmp4);

#if FALSE
      /* This block "factors" the QT parameters (stdev) with the
       * penetrance array (mean). There is little reason to do this,
       * because there are no constraints imposed at this factoring
       * step, so you may as well just go ahead and loop over the
       * parameters and penetrances together. Also, this code won't
       * work for QT model distributions that have more than one
       * parameter. Probably should just delete it.
       *
       * If you do want it, though, you'll also need to uncomment
       * variable tmp3 defined above.
       *
       * OK, you're almost done. Just need to factor the parameter and
       * penetrance arrays together. Since there can be no constraints
       * enforced here, these just factor together brute force. 
       *
       * There will be range->npenet * range->nparam combinations. */
      tmp2 = range->penet;
      tmp3 = range->param;
      i = range->npenet;

      /* Barring constraints, how many values might you have? */
      j = range->npenet * range->nparam;

      /* Allocate new penet and param array structures. */
      range->penet = malloc (range->nlclass * sizeof (double **));
      for (k = 0; k < range->nlclass; k++)
	{
	  range->penet[k] =
	    malloc (NPENET (range->nalleles) * sizeof (double *));
	  for (l = 0; l < NPENET (range->nalleles); l++)
	    range->penet[k][l] = malloc (j * sizeof (double));
	}
      range->param = malloc (range->nlclass * sizeof (double ***));
      for (k = 0; k < range->nlclass; k++)
	{
	  range->param[k] =
	    malloc (NPENET (range->nalleles) * sizeof (double **));
	  for (l = 0; l < NPENET (range->nalleles); l++)
	    {
	      range->param[k][l] =
		malloc (range->npardim * sizeof (double *));
	      for (m = 0; m < range->npardim; m++)
		range->param[k][l][m] = malloc (j * sizeof (double));
	    }
	}

      /* Populate the arrays. */
      range->npenet = 0;
      for (k = 0; k < j; k++)
	{
	  l = 1;
	  for (m = 0; m < range->nlclass; m++)
	    for (n = 0; n < NPENET (range->nalleles); n++)
	      {
		range->penet[m][n][range->npenet] =
		  tmp2[m][n][(((int) (k / l)) % i)];
		for (o = 0; o < range->npardim; o++)
		  range->param[m][n][o][range->npenet] =
		    tmp3[m][n][o][(((int) (k / l)) % i)];
	      }
	  range->npenet++;
	}
      /* Make sure nparam returns the same number. */
      range->nparam = range->npenet;
#if FALSE
      /* Free up the old copies. */
      for (i = 0; i < nlclass; i++)
	{
	  for (j = 0; j < NPENET (range->nalleles); j++)
	    {
	      free (tmp2[i][j]);
	      for (k = 0; k < range->npardim; k++)
		free (tmp3[i][j][k]);
	      free (tmp3[i][j]);
	    }
	  free (tmp2[i]);
	  free (tmp3[i]);
	}
      free (tmp2);
      free (tmp3);
#endif
#endif
    }
}

/**********************************************************************
 * Expand the dprime array into the appropriate lambda array. Unlike
 * the other range expansion functions, this function is not called
 * upfront, but rather as needed when operating under linkage
 * disequilibrium. The two extra parameters are the number of alleles
 * for the two loci under consideration. 
 *
 * Since n and m may be repeated over the course of an analysis, we
 * cache the lambda arrays produced from the dprime values for each
 * value of n and m in modelRange->lambdas, an array of structures of
 * type lambdaCell. That way, we can retrieve the appropriate array if
 * its already been generated, otherwise, we build the array from the
 * dprimes, cache it, and return a pointer to it.
 *
 * TODO: since the array is symmetric regardless of n and m, we could
 * save some storage by reconfiguring the existing array if the mxn
 * version (but not the nxm version) already exists. This shouldn't
 * happen that often, as long as we assume that the first variable
 * corresponds to the disease or trait and the second corresponds to
 * the marker; in this situation, since we are usually not dealing
 * with multiallelic diseases, m=2 and n>=2.
 **********************************************************************/
LambdaCell *
findLambdas (ModelRange * range, int m, int n)
{
  int i = 0, j, k, l = pow (range->ndprime, ((m - 1) * (n - 1)));

  /* First, see if the array of lambdas already exists. */
  while (i < range->nlambdas)
    {
      if (range->lambdas[i].m == m && range->lambdas[i].n == n)
	/* return (range->lambdas[i].lambda); */
	return (&range->lambdas[i]);
      i++;
    }

  /* OK, no matching cached array found. Check to make sure there's
   * room for a new one. */
  if (range->nlambdas == range->maxnlambdas)
    {
      range->lambdas =
	realloc (range->lambdas,
		 (range->maxnlambdas + CHUNKSIZE) * sizeof (LambdaCell));
      range->maxnlambdas = range->maxnlambdas + CHUNKSIZE;
    }
  /* Create the new entry. */
  range->lambdas[range->nlambdas].m = m;
  range->lambdas[range->nlambdas].n = n;
  range->lambdas[range->nlambdas].ndprime = l;
  range->lambdas[range->nlambdas].lambda =
    (double ***) malloc (l * sizeof (double **));
  for (i = 0; i < l; i++)
    {
      range->lambdas[range->nlambdas].lambda[i] =
	(double **) malloc ((m - 1) * sizeof (double *));
      for (j = 0; j < (m - 1); j++)
	range->lambdas[range->nlambdas].lambda[i][j] =
	  (double *) malloc ((n - 1) * sizeof (double));
    }
  range->nlambdas++;

  /* Now populate the values in the appropriate array. */
  for (i = 0; i < l; i++)
    for (j = 0; j < (m - 1); j++)
      for (k = 0; k < (n - 1); k++)
	range->lambdas[range->nlambdas - 1].lambda[i][j][k] =
	  range->dprime[((int) (i / pow (range->ndprime, k))) %
			range->ndprime];

  /* Return a pointer to the appropriate 3 dimensional lambda array. */
  /* return (range->lambdas[range->nlambdas-1].lambda); */
  return (&range->lambdas[range->nlambdas - 1]);
}

/**********************************************************************
 * Dump the model and constraints for debugging purposes.  The level
 * argument indicates if we are pre-expansion (level = 0),
 * post-expansion (level = 1), or post-class-expansion (level =
 * 2). This is ugly, but is only used for debugging so it needn't be
 * excessively pretty!
 **********************************************************************/
void
showRange (ModelRange * range, ModelType * type, int level)
{
  int i, j, k, l;
  printf
    ("======================================================================\n");
  printf ("LEVEL %d MODEL\n", level);
  printf
    ("======================================================================\n");
  printf ("%d GF=", range->ngfreq);
  for (i = 0; i < range->ngfreq; i++)
    printf ("%3.2g ", range->gfreq[i]);

  if (type->trait == DT)
    {
      if (level > 0)
	printf ("\n%d Thetas", range->ntheta);
      for (i = 0; i < range->ngender; i++)
	{
	  if (level == 0)
	    printf ("\n%d Theta[%d]=",
		    (level == 0 ? thetacnt[i] : range->ntheta), i);
	  else
	    printf ("\n  Theta[%d]=", i);
	  for (j = 0; j < (level == 0 ? thetacnt[i] : range->ntheta); j++)
	    printf ("%3.2g ", range->theta[i][j]);
	}
    }
  else
    {
      if (range->tloc)
	{
	  printf ("\n%d Trait Loci", range->ntloc);
	  printf ("\n  Trait Loci=");
	  for (i = 0; i < range->ntloc; i++)
	    printf ("%3.2g ", range->tloc[i]);
	}
      if (range->tthresh)
	{
	  printf ("\n%d Trait Thresholds", range->ntthresh);
	  for (i = 0; i < (level == 2 ? range->nlclass : 1); i++)
	    {
	      printf ("\n  Trait Threshold[%d]=", i);
	      for (j = 0; j < range->ntthresh; j++)
		printf ("%3.2g ", range->tthresh[i][j]);
	    }
	}
    }

  printf ("\n%d Penetrances", (level == 0 ? penetcnt[0] : range->npenet));
  for (i = 0; i < (level == 2 ? range->nlclass : 1); i++)
    for (j = 0; j < NPENET (range->nalleles); j++)
      {
	printf ("\n  Penet[%d][%d]=", i, j);
	for (k = 0; k < (level == 0 ? penetcnt[j] : range->npenet); k++)
	  printf ("%3.2g ", range->penet[i][j][k]);
      }

  if (range->param)
    {
      printf ("\n%d Parameters", ((level == 2 ? range->nlclass : 1) *
				  (level ==
				   0 ? paramcnt[0] : range->nparam)));
      /* for (i = 0; i < (level==2?range->nlclass:1); i++) */
      for (i = 0; i < (level == 2 ? range->nlclass : 1); i++)
	for (j = 0; j < (level == 0 ? 1 : NPENET (range->nalleles)); j++)
	  for (k = 0; k < (level > 0 ? range->npardim : 1); k++)
	    {
	      printf ("\n  Param[%d][%d][%d]=", i, j, k);
	      for (l = 0; l < (level == 0 ? paramcnt[k] : range->nparam); l++)
		printf ("%3.2g ", range->param[i][j][k][l]);
	    }
    }

  if (range->afreq)
    {
      printf ("\n%d AF=", range->nafreq);
      for (i = 0; i < range->nafreq; i++)
	printf ("%3.2g ", range->afreq[i]);
    }

  if (range->dprime)
    {
      printf ("\n%d DPrime", range->ndprime);
      printf ("\n  DPrime=");
      for (i = 0; i < range->ndprime; i++)
	printf ("%3.2g ", range->dprime[i]);
#if FALSE
      /* Just for kicks, show a few off. Recall findLambdas() is
       * called on the fly, so usually these would not be precomputed
       * at all. */
      findLambdas (range, 2, 4);
      findLambdas (range, 2, 3);
      findLambdas (range, 2, 4);
#endif
    }
  if (range->alpha)
    {
      printf ("\n%d Alphas", range->nalpha);
      printf ("\n  Alpha=");
      for (i = 0; i < range->nalpha; i++)
	printf ("%3.2g ", range->alpha[i]);
    }
  printf ("\n");
}

void
showConstraints ()
{
  int i, j;
  char types[25] = DIRECTIVES;

  printf
    ("======================================================================\n");
  printf ("CONSTRAINTS:\n");
  printf
    ("======================================================================\n");

  for (i = 0; i < 4; i++)
    {
      if (constcnt[i] == 0)
	continue;
      printf ("%s constraints:\n",
	      ((i ==
		3) ? "PCC" : ((i == 2) ? "PC" : ((i == 1) ? "CC" : "C"))));
      for (j = 0; j < constcnt[i]; j++)
	{
	  if (i == SIMPLE)
	    printf (" %c%c %s %c%c %s\n",
		    types[constraints[i][j].a1 * 2],
		    types[constraints[i][j].a1 * 2 + 1],
		    ((constraints[i][j].op ==
		      EQ) ? "==" : (constraints[i][j].op ==
				    NE) ? "!=" : (constraints[i][j].op ==
						  GE) ? ">=" : ">"),
		    types[constraints[i][j].a2 * 2],
		    types[constraints[i][j].a2 * 2 + 1],
		    (constraints[i][j].alt == TRUE) ? "|" : "");
	  else if (i == CLASSC)
	    printf (" %c%c %d %s %c%c %d %s\n",
		    types[constraints[i][j].a1 * 2],
		    types[constraints[i][j].a1 * 2 + 1],
		    constraints[i][j].c1,
		    ((constraints[i][j].op ==
		      EQ) ? "==" : (constraints[i][j].op ==
				    NE) ? "!=" : (constraints[i][j].op ==
						  GE) ? ">=" : ">"),
		    types[constraints[i][j].a2 * 2],
		    types[constraints[i][j].a2 * 2 + 1], constraints[i][j].c2,
		    (constraints[i][j].alt == TRUE) ? "|" : "");
	  else if (i == PARAMC)
	    printf (" P%d %c%c %s P%d %c%c %s\n",
		    constraints[i][j].p1,
		    types[constraints[i][j].a1 * 2],
		    types[constraints[i][j].a1 * 2 + 1],
		    ((constraints[i][j].op ==
		      EQ) ? "==" : (constraints[i][j].op ==
				    NE) ? "!=" : (constraints[i][j].op ==
						  GE) ? ">=" : ">"),
		    constraints[i][j].p2, types[constraints[i][j].a2 * 2],
		    types[constraints[i][j].a2 * 2 + 1],
		    (constraints[i][j].alt == TRUE) ? "|" : "");
	  else if (i == PARAMCLASSC)
	    printf (" P%d %c%c %d %s P%d %c%c %d %s\n",
		    constraints[i][j].p1,
		    types[constraints[i][j].a1 * 2],
		    types[constraints[i][j].a1 * 2 + 1],
		    constraints[i][j].c1,
		    ((constraints[i][j].op ==
		      EQ) ? "==" : (constraints[i][j].op ==
				    NE) ? "!=" : (constraints[i][j].op ==
						  GE) ? ">=" : ">"),
		    constraints[i][j].p2, types[constraints[i][j].a2 * 2],
		    types[constraints[i][j].a2 * 2 + 1], constraints[i][j].c2,
		    (constraints[i][j].alt == TRUE) ? "|" : "");
	}
    }
}
