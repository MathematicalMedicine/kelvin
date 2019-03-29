/**********************************************************************
 * Multiprocessor Linkage Analysis
 * Alberto Maria Segre
 * Regex code Nathan Burnette
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
 * TODO: need to read lambda parameters for diseq
 * TODO: handling of multiallelic diseases; parsing e.g., 00, 21, 11 etc
 * TODO: check distribs with nparam>1 (not sure it works)
 * TODO: liability class should be the slowest moving index so that we
 *       can censor at run-time if ped has no such LC
 **********************************************************************/

/**********************************************************************
 * Regular expressions used to parse configuration file
 * lines. 
 **********************************************************************/
/* DIRECTIVE matches: GF, DD, Dd, dd, Th, Tm, Tf, 00..99 */
#define DIRECTIVE "^[[:space:]]*[[:digit:]dDGT][[:digit:]dDfFhm]"
/* PARAMETER matches: P[0-9] (note implicit 1 digit limit on distrib params index) */
#define PARAMETER "^[[:space:]]*P([[:digit:]])"
/* Match 1 or 3 floating point numbers with optional trailing semicolon. */
#define DOUBLE1 "[[:space:]]+([-]?[[:digit:]]*.[[:digit:]]+)[;]?"
#define DOUBLE3 "[[:space:]]+([-]?[[:digit:]]*.[[:digit:]]+)[[:space:]]+([-]?[[:digit:]]*.[[:digit:]]+)[[:space:]]+([-]?[[:digit:]]*.[[:digit:]]+)[;]?"
/* Match a constraint of the form, e.g., DD < Dd */
#define CONSTRAINT2 "[[:space:]]*([[:digit:]dDT][[:digit:]dDfm])[[:space:]]*([>]|[>][=]|[=][=]|[!][=])[[:space:]]*([[:digit:]dDT][[:digit:]dDfm])[;]?"
/* Match a constraint with liability classes, e.g., DD 1 < DD 2 */
#define CLASSCONSTRAINT2 "[[:space:]]*([[:digit:]dDT][[:digit:]dDfm])[[:space:]]*([[:digit:]])[[:space:]]*([>]|[>][=]|[=][=]|[!][=])[[:space:]]*([[:digit:]dDT][[:digit:]dDfm])[[:space:]]*([[:digit:]])[;]?"
/* Match a constraint of the form, e.g., Px DD < Px Dd */
#define PCONSTRAINT2 "[[:space:]]*P([[:digit:]])[[:space:]]+([[:digit:]dDT][[:digit:]dDfm])[[:space:]]*([>]|[>][=]|[=][=]|[!][=])[[:space:]]*P([[:digit:]])[[:space:]]+([[:digit:]dDT][[:digit:]dDfm])[;]?"
/* Match a constraint of the form, e.g., Px DD 1 < Px Dd 2 */
#define PCLASSCONSTRAINT2 "[[:space:]]*P([[:digit:]])[[:space:]]+([[:digit:]dDT][[:digit:]dDfm])[[:space:]]*([[:digit:]])[[:space:]]*([>]|[>][=]|[=][=]|[!][=])[[:space:]]*P([[:digit:]])[[:space:]]+([[:digit:]dDT][[:digit:]dDfm])[[:space:]]*([[:digit:]])[;]?"

/* Strings used to parse configuration file lines (see getType). */
#define DIRECTIVES "GFThTmTfDDDddd"
#define GF 0
#define Th 1
#define Tm 2
#define Tf 3
#define DD 4
#define Dd 5
#define dd 6
/* Strings used to parse configuration file constraints (see
 * getOperator). Matches:
 *   ==, !=, >, >=
 *
 * TODO: Should probably allow a single = as alias for === to allow
 * for non-programmers. 
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
 * but they are all stored in a uniform constriant structure.
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
  int type;	/* SIMPLE, CLASSC, PARAMC, PARAMCLASSC */
  int a1;
  int c1;
  int p1;
  int a2;
  int c2;
  int p2;
  int op;
  int alt;	/* Is this one of a disjunctive set? */
} 
Constraint;

/**********************************************************************
 * Some global variables used only within config.c -- avoids having to
 * carry these around in the model structure, since they are only
 * useful while setting up the penetrances, etc.
 **********************************************************************/
int maxgfreq;	/* Max number of gfreqs in model (for dynamic alloc) */
int *penetcnt;	/* Array of number of penetrances */
int *penetmax; 	/* Array of max penetrances */
int *paramcnt;  /* Number of QT/CT parameters */
int *parammax;  /* Max QT/CT parameters */
int *thetacnt;  /* Array of number of thetas */
int *thetamax;	/* Array of max thetas */
Constraint *constraints[4];   /* Array of constraints by type */
int constmax[4] = {0,0,0,0};  /* Max constraints in array (for dynamic alloc) */
int constcnt[4] = {0,0,0,0}; /* Current number of constraints in array */

/* Chunk size used in reallocating arrays of gene frequencies,
 * penetrances, thetas, and constraints. */
#define CHUNKSIZE 64

/**********************************************************************
 * Some internal prototypes.
 **********************************************************************/
int getType (char *line);
int getOperator (char *line);
int getInteger (char *line);
void addRange (ModelRange *range, int type, double lo, double hi, double incr);
void addTheta (ModelRange *range, int type, double val);
void addPenetrance (ModelRange *range, int type, double val);
void addGeneFreq (ModelRange *range, double val);
void addConstraint (int type, int a1, int c1, int p1, 
		    int op, int a2, int c2, int p2, int disjunct);
void addParameter (ModelRange *range, int dim, double val);
int checkThetas (ModelRange *range, int i);
int checkPenets (ModelRange *range, int i);
int checkClassPenets (ModelRange *range, int i);
int checkParams (ModelRange *range, int i);
int checkClassParams (ModelRange *range, int i);
void sortRange (ModelRange *range);
inline void swap (double *array, int i, int j);
void quicksort (double *array, int lo, int hi);
void uniqRange (ModelRange *range);
inline int uniquify (double *array, int len);
void expandRange (ModelRange *range);
void expandClass (ModelRange *range);
void showRange (ModelRange *range, int level);
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
readConfigFile (char *file, ModelType *modelType, 
		ModelRange *modelRange, ModelOptions *modelOptions)
{
  FILE *fp;				/* Filepointer. */
  int i = 0;				/* Number of lines read. */
  char *start ;
  char line[KMAXLINELEN + 1];		/* Current line. */
  int dir1, dir2, op;			/* For parsing constraints. */
  int a1, a2, c1, c2, p1;

  regex_t *buffer0;
  regex_t *buffer1;
  regex_t *buffer2;
  regex_t *buffer3;
  regex_t *buffer4;
  regex_t *buffer5;
  regex_t *buffer6;
  regex_t *buffer7;
  regmatch_t match[8];			/* Store extracted values. */

  /* Set up the default model values. */
  modelType->type = TP;
  modelType->trait = DT;
  modelOptions->equilibrium = LE;
  modelRange->nalleles = 2;
  modelRange->nlclass = 1;
  modelRange->npardim = 0;

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
  KASSERT((fp = fopen (file, "r")),
	   "Can't open configuration file %s; aborting.\n", file);

  /* Set up the regular expression pattern buffer we'll use to parse
   * the lines in the file. The flags REG_EXTENDED and REG_ICASE mean
   * we'll use Posix extended regular expressions and we'll ignore
   * case. */
  KASSERT((regcomp(buffer0, DIRECTIVE, (REG_EXTENDED | REG_NEWLINE)) == 0),
	   "Internal error in regular expression; aborting.\n");
  KASSERT((regcomp(buffer1, DOUBLE1, (REG_EXTENDED | REG_NEWLINE)) == 0),
	   "Internal error in regular expression; aborting.\n");
  KASSERT((regcomp(buffer2, DOUBLE3, (REG_EXTENDED | REG_NEWLINE)) == 0),
	   "Internal error in regular expression; aborting.\n");
  KASSERT((regcomp(buffer3, CONSTRAINT2, (REG_EXTENDED | REG_NEWLINE)) == 0),
	   "Internal error in regular expression; aborting.\n");
  KASSERT((regcomp(buffer4, CLASSCONSTRAINT2, (REG_EXTENDED | REG_NEWLINE)) == 0),
	   "Internal error in regular expression; aborting.\n");
  KASSERT((regcomp(buffer5, PARAMETER, (REG_EXTENDED | REG_NEWLINE)) == 0),
	   "Internal error in regular expression; aborting.\n");
  KASSERT((regcomp(buffer6, PCONSTRAINT2, (REG_EXTENDED | REG_NEWLINE)) == 0),
	   "Internal error in regular expression; aborting.\n");
  KASSERT((regcomp(buffer7, PCLASSCONSTRAINT2, (REG_EXTENDED | REG_NEWLINE)) == 0),
	   "Internal error in regular expression; aborting.\n");

  /* Start scanning each line of configuration file input. We'll parse
   * each line by looking for the easy cases first, then using C
   * regexps to parse the harder configuration directives. */
  while (fgets (line, KMAXLINELEN, fp))
    {
      /* Before we try to parse it, check to see if the line may be
       * too long, and give up if it is. */
      KASSERT((strlen (line) < KMAXLINELEN),
	       "Line %d in configuration file %s exceeds %d characters; aborting.\n", i,
	       file, KMAXLINELEN);

      /* Flush lines starting with a comment character or consisting
       * of only a newline. */
      if ((strlen(line) == 1) || (strncmp (line, "#", 1) == 0))
	continue;

      /* TODO: flush lines containing only whitespace characters. */

      /* Directives that take no arguments. */
      if (strncmp (line, "TP", 2) == 0)
	{
	  modelType->type = TP;		/* 2 point */
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configuring for 2 point analysis\n"); */
	  continue;
	}
      if (strncmp (line, "MP", 2) == 0)
	{
	  modelType->type = MP;		/* Multipoint */
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configuring for multipoint analysis\n"); */
	  continue;
	}

      if (strncmp (line, "LE", 2) == 0)
	{
	  modelOptions->equilibrium = LE;		/* Linkage equilibrium */
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configuring for linkage equilibrium\n"); */
	  continue;
	}
      if (strncmp (line, "LD", 2) == 0)
	{
	  modelOptions->equilibrium = LD;		/* Linkage disequilibrium */
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configuring for linkage disequilibrium\n"); */
	  continue;
	}

      if (sscanf (line, "LC %d", &modelRange->nlclass) == 1)
	{
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configuring for %d liability classes\n", modelRange->nlclass); */
	  continue;
	}

      if (strncmp (line, "DT", 2) == 0)
	{
	  modelType->trait = DT;		/* Dichotomous trait */
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configuring for dichotomous trait\n"); */
	  continue;
	}

      /* Directives that take a single integer argument. */
      if (sscanf (line, "DA %d", &modelRange->nalleles) == 1)	/* Disease alleles */
	{
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configuring for %d disease alleles\n", model->nalleles); */
	  continue;
	}
      /* TODO: add clauses (or generalize the following two clauses)
       * for other-than-normal distributions. */
      if (strncmp (line, "QT normal", 9) == 0)
	{
	  modelType->trait = QT;		/* Quantitative trait */
	  modelType->distrib = NORMAL_DISTRIBUTION;
	  modelRange->npardim = 1;
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configuring for quantitative trait\n"); */
	  continue;
	}
#if FALSE
      if (sscanf (line, "CT normal %lg", &modelType->thresh) == 1)
	{
	  modelType->trait = CT;		/* Combined trait */
	  modelType->distrib = NORMAL_DISTRIBUTION;
	  modelRange->npardim = 1;
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configuring for combined trait\n"); */
	  continue;
	}
#endif
      /* Directives that take a single string argument. */
      if (sscanf (line, "PD %s", pedfile) == 1)		/* Pedigree file */
	{
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configure pedigree file %s\n", pedfile); */
	  continue;
	}
      if (sscanf (line, "MK %s", markerfile) == 1)	/* Marker file */
	{
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configure marker file %s\n", markerfile); */
	  continue;
	}
      if (sscanf (line, "MP %s", mapfile) == 1)		/* Map file */
	{
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configure map file %s\n", mapfile); */
	  continue;
	}
#if FALSE
      if (sscanf (line, "LP %s", loopfile) == 1)	
	{
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configure loop file %s\n", loopfile); */
	  continue;
	}
#endif
      if (sscanf (line, "OF %s", outfile) == 1)		/* Output file */
	{
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Configure output file %s\n", outfile); */
	  continue;
	}
      if (sscanf (line, "UP %d", &dir1) == 1)
	{
	  /* Allocate space for the identifier based on how many
	   * digits in the scanned integer. Remember to leave a space
	   * for the terminator. */
	  modelOptions->sUnknownPersonID = realloc (modelOptions->sUnknownPersonID,
						   (dir1/10)+2*sizeof (char));
	  snprintf (modelOptions->sUnknownPersonID, (dir1/10)+2, "%d", dir1);
	  /* KLOG(LOGDEFAULT, LOGDEBUG, "Unknown person identifier is %s\n", modelOptions->sUnknownPersonID); */
	  continue;
	}

      /* Complex directives that require regular expression matching
       * to handle the arguments. Most of these takes at least one
       * floating point argument, and sometimes two or three.
       *
       * We'll handle each type of directive separately. But, first,
       * "chomp" the newline to remove the newline and ensure null
       * termination. */
      line[strlen(line) - 1] = '\0';

      /* First, look for simple constraints on thetas, gene
       * frequencies, or penetrances. */
      if (regexec (buffer3, line, 4, &match[0], 0) != REG_NOMATCH)
	{
	  /* Matching constraint line. Find out what type of
	   * constraint it is by looking it up in DIRECTIVES. */
	  KASSERT(((dir1 = getType (line)) != ERROR),
		  "Bad constraint directive '%s'; aborting.\n", line);
	  KASSERT(((op = getOperator (line+match[2].rm_so)) != ERROR),
		  "Bad constraint operator '%s'; aborting.\n", line);
	  KASSERT(((dir2 = getType (line+match[3].rm_so)) != ERROR),
		  "Bad constraint directive '%s'; aborting.\n", line);
 
	  /* Add the constraint. */
	  addConstraint (SIMPLE, dir1, 0, 0, op, dir2, 0, 0, FALSE);

	  /* Next, loop through any semicolon-separated disjuncts. */
	  start = line;
	  while ((start = start + match[0].rm_eo) < line+strlen(line))
	    {
	      /* Add the disjunction to the first constraint. */
	      if (regexec (buffer3, start, 4, &match[0], REG_NOTBOL) != REG_NOMATCH)
		{
		  KASSERT(((dir1 = getType (start+match[1].rm_so)) != ERROR),
			  "Bad constraint directive '%s'; aborting.\n", start);
		  KASSERT(((op = getOperator (start+match[2].rm_so)) != ERROR),
			  "Bad constraint directive '%s'; aborting.\n", start);
		  KASSERT(((dir2 = getType (start+match[3].rm_so)) != ERROR),
			  "Bad constraint directive '%s'; aborting.\n", start);
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
	   * frequencies, or thetas. Find out what type of constraint
	   * it is by looking it up in DIRECTIVES. */
	  KASSERT(((dir1 = getType (line)) != ERROR),
		  "Bad constraint directive '%s'; aborting.\n", line);
	  KASSERT(((op = getOperator (line+match[3].rm_so)) != ERROR),
		  "Bad constraint operator '%s'; aborting.\n", line);
	  KASSERT(((dir2 = getType (line+match[4].rm_so)) != ERROR),
		  "Bad constraint directive '%s'; aborting.\n", line);
 
	  /* Add the constraint. */
	  addConstraint (CLASSC, 
			 dir1, getInteger(line+match[2].rm_so), 0,
			 op, 
			 dir2, getInteger(line+match[5].rm_so), 0,
			 FALSE);
	  
	  /* Next, loop through any semicolon-separated disjuncts. */
	  start = line;
	  while ((start = start + match[0].rm_eo) < line+strlen(line))
	    {
	      /* Add the disjunction to the first constraint. */
	      if (regexec (buffer4, start, 6, &match[0], REG_NOTBOL) != REG_NOMATCH)
		{
		  KASSERT(((dir1 = getType (start+match[1].rm_so)) != ERROR),
			  "Bad constraint directive '%s'; aborting.\n", start);
		  KASSERT(((op = getOperator (start+match[3].rm_so)) != ERROR),
			  "Bad constraint directive '%s'; aborting.\n", start);
		  KASSERT(((dir2 = getType (start+match[4].rm_so)) != ERROR),
			  "Bad constraint directive '%s'; aborting.\n", start);
		  addConstraint (CLASSC, 
				 dir1, getInteger(start+match[2].rm_so), -1, 
				 op, 
				 dir2, getInteger(start+match[5].rm_so), -1, 
				 TRUE);
		}
	      else
		break;
	    }
	  continue;
	}
      else if (regexec (buffer6, line, 6, &match[0], 0) != REG_NOMATCH)
	{
	  /* Now look for constraints on parameters. */
	  dir1 = getInteger (line+match[1].rm_so);
	  a1 = getType (line+match[2].rm_so);
	  op = getOperator (line+match[3].rm_so);
	  dir2 = getInteger (line+match[4].rm_so);
	  a2 = getType (line+match[5].rm_so);

	  /* SHOULD CHECK CONSTRAINT FORMAT! */

	  /* Add the constraint */
	  addConstraint (PARAMC, a1, 0, dir1, op, a2, 0, dir2, FALSE);

	  /* Next, loop through any semicolon-separated disjuncts. */
	  start = line;
	  while ((start = start + match[0].rm_eo) < line+strlen(line))
	    {
	      /* Add the disjunction to the first constraint. */
	      if (regexec (buffer6, start, 6, &match[0], REG_NOTBOL) != REG_NOMATCH)
		{
		  /* Matching parameter constraint line. */
		  dir1 = getInteger (start+match[1].rm_so);
		  a1 = getType (start+match[2].rm_so);
		  op = getOperator (start+match[3].rm_so);
		  dir2 = getInteger (start+match[4].rm_so);
		  a2 = getType (start+match[5].rm_so);
		  
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
	   * parameters. */
	  dir1 = getInteger (line+match[1].rm_so);
	  a1 = getType (line+match[2].rm_so);
	  c1 = getInteger(line+match[3].rm_so);
	  op = getOperator (line+match[4].rm_so);
	  dir2 = getInteger (line+match[5].rm_so);
	  a2 = getType (line+match[6].rm_so);
	  c2 = getInteger(line+match[7].rm_so);

	  /* SHOULD CHECK CONSTRAINT FORMAT! */

	  /* Add the constraint */
	  addConstraint (PARAMCLASSC, a1, c1, dir1, op, a2, c2, dir2, FALSE);

	  /* Next, loop through any semicolon-separated disjuncts. */
	  start = line;
	  while ((start = start + match[0].rm_eo) < line+strlen(line))
	    {
	      /* Add the disjunction to the first constraint. */
	      if (regexec (buffer6, start, 6, &match[0], REG_NOTBOL) != REG_NOMATCH)
		{
		  /* Matching parameter constraint line. */
		  dir1 = getInteger (start+match[1].rm_so);
		  a1 = getType (start+match[2].rm_so);
		  c1 = getInteger(start+match[3].rm_so);
		  op = getOperator (start+match[4].rm_so);
		  dir2 = getInteger (start+match[5].rm_so);
		  a2 = getType (start+match[6].rm_so);
		  c2 = getInteger(start+match[7].rm_so);
		  
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
	  /* Done with constraints, now look for value
	   * specifications. Find out what type of directive it is by
	   * looking it up in DIRECTIVES. */
	  dir1 = getType (line);

	  /* Next, loop through the semicolon-separated arguments. */
	  start = line + match[0].rm_eo;
	  while (TRUE)
	    {
	      if (regexec (buffer2, start, 4, &(match[1]), REG_NOTBOL) != REG_NOMATCH)
		{
		  /* Matches three numbers. */
		  addRange (modelRange, dir1,
				strtod(start+match[2].rm_so, NULL),
				strtod(start+match[3].rm_so, NULL),
				strtod(start+match[4].rm_so, NULL));
		}
	      else if (regexec (buffer1, start, 2, &(match[1]), REG_NOTBOL) != REG_NOMATCH)
		{
		  /* Matches one number. */
		  if (dir1 == Th || dir1 == Tm || dir1 == Tf)
		    addTheta (modelRange, dir1, strtod(start+match[2].rm_so, NULL));
		  else if (dir1 == GF)
		    addGeneFreq (modelRange, strtod(start+match[2].rm_so, NULL));
		  else
		    /* All parsing is done only in pre-expansion mode,
		     * hence class is always moot. Also, subtract base
		     * value of DD from dir1 to get 0-offset value.*/
		    addPenetrance (modelRange, dir1-DD, strtod(start+match[2].rm_so, NULL));
		}
	      else
		{
		  /* Doesn't match anything. This is a badly formatted
		   * line which should be flagged. */
		  KASSERT(FALSE, "Ill-formed line in configuration file: '%s'\n", line);
		}

	      /* Update the start pointer and check for termination. */
	      start = start + match[1].rm_eo;
	      if (start >= line+strlen(line))
		break;
	    }
	  continue;
	}
      else if (regexec (buffer5, line, 1, &match[0], 0) != REG_NOMATCH)
	{
	  /* Matching a parameter specification. Find out which
	   * parameter it is by parsing the parameter directive, then
	   * make sure that the parameter number is within the
	   * expected number of parameters for this distribution. 
	   *
	   * Remember, the parameter number here will be 1-indexed. */
	  p1 = strtod(line+match[0].rm_so, NULL);
	  KASSERT((p1<=modelRange->npardim), "Illegal parameter %d.\n", p1);

	  /* Next, loop through the semicolon-separated arguments. */
	  start = line + match[0].rm_eo;
	  while (TRUE)
	    {
	      if (regexec (buffer2, start, 4, &(match[1]), REG_NOTBOL) != REG_NOMATCH)
		{
		  /* Matches three numbers. We'll fake out addRange by
		   * giving it a negative type to distinguish it from
		   * the other types of values we handle. So, for
		   * example, parameter 2 would be -3. */
		  addRange (modelRange, -(p1+1),
				strtod(start+match[2].rm_so, NULL),
				strtod(start+match[3].rm_so, NULL),
				strtod(start+match[4].rm_so, NULL));
		}
	      else if (regexec (buffer1, start, 2, &(match[1]), REG_NOTBOL) != REG_NOMATCH)
		{
		  /* Matches one number. */
		  addParameter (modelRange, p1, strtod(start+match[2].rm_so, NULL));
		}
	      else
		{
		  /* Doesn't match anything. This is a badly formatted
		   * line which should be flagged. */
		  KASSERT(FALSE, "Ill-formed line in configuration file: '%s'\n", line);
		}

	      /* Update the start pointer and check for termination. */
	      start = start + match[1].rm_eo;
	      if (start >= line+strlen(line))
		break;
	    }
	  continue;
	}
      else
	KASSERT(FALSE, "Ill-formed line in configuration file: '%s'\n", line);
    }

  /* Done. Release the pattern spaces. */
  regfree (buffer0);
  regfree (buffer1);
  regfree (buffer2);
  regfree (buffer3);
  regfree (buffer4);
  regfree (buffer5);
  regfree (buffer6);
  regfree (buffer7);

  /* Check the integrity of what you've read. Here is where you check
   * for things like, e.g., no parameters specified for QT/CT, or
   * parameters specified for DT. */
  KASSERT((modelType->trait != QT || modelRange->param),
	  "Failure to provide distribution parameters for quantitative trait.\n");
  KASSERT((modelType->trait != CT || modelRange->param),
	  "Failure to provide distribution parameters for combined trait.\n");
  KASSERT((modelType->trait != DT || !modelRange->param),
	  "Attempt to provide distribution parameters for dichotomous trait.\n");

  /* Sort the values in the final model. Sorted values better support
   * the application of constraints. */
  sortRange (modelRange);

  /* Once sorted, removing duplicates is easy. */
  uniqRange (modelRange);

#if TRUE
  /* Show the unexpanded model. At level 0, all elements are sorted
   * and unique, but we may have:
   *  gfreq: ok
   *  thetas: nonuniform lengths of male/female values
   *  penet: nonuniform lengths of values by allele, lclass=0
   *  param: nonuniform lengths of values by dimension, lclass=0, allele=0 */
  showRange (modelRange, 0); 
  /* Show the constraints. */
  showConstraints ();
#endif  

  /* Expand the model, honoring constraints. */
  expandRange (modelRange);
#if TRUE
  /* Show the partially expanded model. At level 1, following
   * expandRange(), we will have refined the model specification while
   * enforcing all specified constraints that do not involve liability
   * classes:
   *  gfreq: ok
   *  thetas: ok (male/female combinations fully "factored")
   *  penet:
   *  param: */
  showRange (modelRange, 1);
#endif  

  /* Expand the liability classes, but only if necessary and always
   * honoring inter-class constraints. */
  if (modelRange->nlclass > 1)
    expandClass (modelRange);

#if TRUE
  /* At level 2, all constraints (including those between classes) are
   * honored, but penet[][][], param[][][][] are not yet fully
   * "factored". */
  showRange (modelRange, 2);
#endif  

  /* Done. Free the constraints; you're done with them. */
#if FALSE
  for (i = 0; i < 4; i++)
    if (constraints[i])
      free (constraints[i]);
  free(constraints);
  free(constcnt);
  free(constmax);
#endif

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
  char types[15] = DIRECTIVES;
  char *ptr = types;

  while (ptr < types + strlen(types))
    {
      if (strncmp (ptr, line, 2) == 0)
	break;
      ptr += 2;
    }
  /* If no match is found, return ERROR. */
  if (ptr - types >= strlen(types))
    return (ERROR);
  else
    return ((ptr-types)/2);
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

  while (ptr < types + strlen(types))
    {
      if (strncmp (ptr, line, 2) == 0)
	break;
      ptr += 2;
    }
  /* If no match is found, return ERROR. */
  if (ptr - types >= strlen(types))
    return (ERROR);
  else
    return ((ptr-types)/2);
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
addTheta (ModelRange *range, int type, double val)
{
  int i;
  double **dtmp;
  int *itmp;
  
  /* First, if this is the first access to the structure, you must
   * initialize it. */
  if (!range->theta)
    {
      if (type == Th)
	{
	  range->ngender = 1;
	  range->theta = malloc (range->ngender * sizeof (double *));
	  range->theta[SEXAV] = malloc (CHUNKSIZE * sizeof(double));
	  thetacnt = malloc (range->ngender * sizeof (int *));
	  thetamax = malloc (range->ngender * sizeof (int *));
	  thetacnt[SEXAV] = thetamax[SEXAV] = 0;
	}
      else
	{
	  range->ngender = 2;
	  range->theta = malloc (range->ngender * sizeof (double *));
	  range->theta[SEXML] = malloc (CHUNKSIZE * sizeof(double));
	  range->theta[SEXFM] = malloc (CHUNKSIZE * sizeof(double));
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
					 (thetamax[SEXAV] + CHUNKSIZE) * sizeof(double));
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
					     (thetamax[SEXFM] + CHUNKSIZE) * sizeof(double));
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
					     (thetamax[SEXML] + CHUNKSIZE) * sizeof(double));
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
addPenetrance (ModelRange *range, int type, double val)
{
  int i, j;

  /* First, if this is the first access to the structure, you must
   * initialize it. */
  if (!range->penet)
    {
      /* Initialize: remember, if you are a pre-expansion penet[][][]
       * array, the first dimension (liability class) will always have
       * a dimension of size 1. */
      range->penet = malloc (sizeof (double **));
      i = NPENET(range->nalleles);
      range->penet[0] = malloc (i * sizeof (double *));
      for (j = 0; j < i; j++)
	range->penet[0][j] = malloc (CHUNKSIZE * sizeof(double));

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
   * type is an integer, where DD=00=0; Dd=01=1; and dd=10=2. We get
   * the type by subtracting DD, the "base type" from the integer
   * count at invocation. */
  if (penetcnt[type] == penetmax[type])
    {
      range->penet[0][type] = realloc (range->penet[0][type], 
				 (penetmax[type] + CHUNKSIZE) * sizeof(double));
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
addGeneFreq (ModelRange *range, double val)
{
  /* Initialize the structure if first access. */
  if (!range->gfreq)
    range->ngfreq = maxgfreq = 0;
  /* Enlarge array if necessary. */
  if (range->ngfreq == maxgfreq)
    {
      range->gfreq = realloc (range->gfreq, 
			      (maxgfreq + CHUNKSIZE) * sizeof(double));
      maxgfreq = maxgfreq + CHUNKSIZE;
    }
  /* Add the element. */
  range->gfreq[range->ngfreq] = val;
  range->ngfreq++;
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
addParameter (ModelRange *range, int dim, double val)
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
	range->param[0][0][i] = malloc (CHUNKSIZE * sizeof(double));

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
					(parammax[dim] + CHUNKSIZE) * sizeof(double));
      parammax[dim] = parammax[dim] + CHUNKSIZE;
    }
  /* Add the element. */
  range->param[0][0][dim][paramcnt[dim]] = val;
  paramcnt[dim]++;
}

/**************q********************************************************
 * Add one more range to genefrequencies, penetrances, parameter, or
 * theta. Parameters are handled by passing a negative type; so -3
 * correponds to the third parameter.
 **********************************************************************/
void
addRange (ModelRange *range, int type, double lo, double hi, double incr)
{
  double val = lo;
  char types[15] = DIRECTIVES;
  while (val <= hi)
    {
      KASSERT(((val >= 0) &&
		(((type == Th || type == Tm || type == Tf) && val <= 0.5) ||
		 (type != Th && type != Tm && type != Tf && val <= 1.0))), 
	       "Bad configuration parameter %c%c=%g; aborting.\n",
	       types[type*2], types[type*2+1], val);

      /* OK, add it. */
      if (type == Th || type == Tm || type == Tf)
	addTheta (range, type, val);
      else if (type == GF)
	addGeneFreq (range, val);
      else if (type < 0)
	addParameter (range, -type-1, val);
      else
	/* addRange() is only used in pre-expansion mode, hence class
	 * is always ERROR. Also, subtract base value of DD from type
	 * to get 0-offset value.*/
	addPenetrance (range, type-DD, val);

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
  char types[15] = DIRECTIVES;

  /* Check for meaningless constraints. TODO: do more of this! */
  KASSERT((((a1==Tm && a2==Tf) || (a1==Tf && a2==Tm)) ||
	    (a1>=DD && a2>=DD && a1<=dd && a2<=dd)), 
	   "Meaningless constraint %c%c %s %c%c %s; aborting.\n",
	    types[a1*2],
	    types[a1*2+1],
	    ((op==NE)?"!=":(op==GE)?">=":">"),
	    types[a2*2],
	    types[a2*2+1],
	    (disjunct==TRUE)?"*":"");

  /* Allocate more space if necessary. */
  if (constmax[type] == constcnt[type])
    {
      constraints[type] = (Constraint *) realloc(constraints[type], 
						 (constmax[type] + CHUNKSIZE) * sizeof(Constraint));
      constmax[type] = constmax[type]+CHUNKSIZE;
    }

  /* Mark previous constraint as having a disjunct, if applicable. */
  if (disjunct && constcnt[type] > 0)
    constraints[type][constcnt[type]-1].alt = TRUE;
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
checkThetas (ModelRange *range, int i)
{
  int j = 0;
#if FALSE
  char types[15]=DIRECTIVES;
#endif

  /* Scan through each constraint. */
  while (j < constcnt[SIMPLE])
    {
      /* Relies on Tm being the lowest numbered penetrance. */
      if ((constraints[SIMPLE][j].a1 == Tm && constraints[SIMPLE][j].a2 == Tf) ||
	  (constraints[SIMPLE][j].a1 == Tf && constraints[SIMPLE][j].a2 == Tm))
	{
	  if ((constraints[SIMPLE][j].op == EQ && 
	       range->theta[constraints[SIMPLE][j].a1-Tm][i] !=
	       range->theta[constraints[SIMPLE][j].a2-Tm][i]) ||
	      (constraints[SIMPLE][j].op == NE && 
	       range->theta[constraints[SIMPLE][j].a1-Tm][i] !=
	       range->theta[constraints[SIMPLE][j].a2-Tm][i]) ||
	      (constraints[SIMPLE][j].op == GT &&
	       range->theta[constraints[SIMPLE][j].a1-Tm][i] >
	       range->theta[constraints[SIMPLE][j].a2-Tm][i]) ||
	      (constraints[SIMPLE][j].op == GE &&
	       range->theta[constraints[SIMPLE][j].a1-Tm][i] >=
	       range->theta[constraints[SIMPLE][j].a2-Tm][i]))
	    {
	      /* Satisfied; skip other disjuncts. */
	      while (constraints[SIMPLE][j].alt == TRUE)
		j++;
	    }
	  else if (constraints[SIMPLE][j].alt == FALSE)
	    {
#if FALSE
	      fprintf (stderr, "%c%c %s %c%c %s => %g %g %g\n", 
		       types[(constraints[SIMPLE][j].a1)*2],
		       types[(constraints[SIMPLE][j].a1)*2+1],
		       ((constraints[SIMPLE][j].op==EQ)?"==":(constraints[SIMPLE][j].op==NE)?"!=":(constraints[SIMPLE][j].op==GE)?">=":">"),
		       types[(constraints[SIMPLE][j].a2)*2],
		       types[(constraints[SIMPLE][j].a2)*2+1],
		       (constraints[SIMPLE][j].alt==TRUE)?"*":"",
		       range->penet[0][i],
		       range->penet[1][i],
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
checkPenets (ModelRange *range, int i)
{
  int j = 0;
#if FALSE
  char types[15]=DIRECTIVES;
#endif

  /* Scan through each constraint. */
  while (j < constcnt[SIMPLE])
    {
      /* Relies on DD being the lowest numbered penetrance. */
      if (constraints[SIMPLE][j].a1 >= DD && constraints[SIMPLE][j].a2 >= DD)
	{
	  if ((constraints[SIMPLE][j].op == EQ && 
	       range->penet[0][constraints[SIMPLE][j].a1-DD][i] ==
	       range->penet[0][constraints[SIMPLE][j].a2-DD][i]) ||
	      (constraints[SIMPLE][j].op == NE && 
	       range->penet[0][constraints[SIMPLE][j].a1-DD][i] !=
	       range->penet[0][constraints[SIMPLE][j].a2-DD][i]) ||
	      (constraints[SIMPLE][j].op == GT &&
	       range->penet[0][constraints[SIMPLE][j].a1-DD][i] >
	       range->penet[0][constraints[SIMPLE][j].a2-DD][i]) ||
	      (constraints[SIMPLE][j].op == GE &&
	       range->penet[0][constraints[SIMPLE][j].a1-DD][i] >=
	       range->penet[0][constraints[SIMPLE][j].a2-DD][i]))
	    {
	      /* Satisfied; skip other disjuncts. */
	      while (constraints[SIMPLE][j].alt == TRUE)
		j++;
	    }
	  else if (constraints[SIMPLE][j].alt == FALSE)
	    {
#if FALSE
	      fprintf (stderr, "%c%c %s %c%c %s => %g %g %g\n", 
		       types[(constraints[SIMPLE][j].a1)*2],
		       types[(constraints[SIMPLE][j].a1)*2+1],
		       ((constraints[SIMPLE][j].op==EQ)?"==":(constraints[SIMPLE][j].op==NE)?"!=":(constraints[SIMPLE][j].op==GE)?">=":">"),
		       types[(constraints[SIMPLE][j].a2)*2],
		       types[(constraints[SIMPLE][j].a2)*2+1],
		       (constraints[SIMPLE][j].alt==TRUE)?"*":"",
		       range->penet[0][0][i],
		       range->penet[0][1][i],
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
 * with dynamically.
 **********************************************************************/
int
checkClassPenets (ModelRange *range, int i)
{
  int j = 0;
#if FALSE
  char types[15]=DIRECTIVES;
#endif

  /* Scan through each constraint: inter-class constraints are alaways
   * on penetrances, not, e.g., gene frequencies or thetas. */
  while (j < constcnt[CLASSC])
    {
      /* Now check the constraint. */
      if ((constraints[CLASSC][j].op == EQ && 
	   range->penet[constraints[CLASSC][j].c1-1][constraints[CLASSC][j].a1-DD][i] ==
	   range->penet[constraints[CLASSC][j].c2-1][constraints[CLASSC][j].a2-DD][i]) ||
	  (constraints[CLASSC][j].op == NE && 
	   range->penet[constraints[CLASSC][j].c1-1][constraints[CLASSC][j].a1-DD][i] !=
	   range->penet[constraints[CLASSC][j].c2-1][constraints[CLASSC][j].a2-DD][i]) ||
	  (constraints[CLASSC][j].op == GT &&
	   range->penet[constraints[CLASSC][j].c1-1][constraints[CLASSC][j].a1-DD][i] >
	   range->penet[constraints[CLASSC][j].c2-1][constraints[CLASSC][j].a2-DD][i]) ||
	  (constraints[CLASSC][j].op == GE &&
	   range->penet[constraints[CLASSC][j].c1-1][constraints[CLASSC][j].a1-DD][i] >=
	   range->penet[constraints[CLASSC][j].c2-1][constraints[CLASSC][j].a2-DD][i]))
	{
	  /* Satisfied; skip other disjuncts. */
	  while (constraints[CLASSC][j].alt == TRUE)
	    j++;
	}
      else if (constraints[CLASSC][j].alt == FALSE)
	{
#if FALSE
	  fprintf (stderr, "!%c%c %d %s %c%c %d %s => %g %s %g\n", 
		   types[(constraints[CLASSC][j].a1)*2],
		   types[(constraints[CLASSC][j].a1)*2+1],
		   constraints[CLASSC][i].c1,
		   ((constraints[CLASSC][j].op==EQ)?"==":(constraints[CLASSC][j].op==NE)?"!=":(constraints[CLASSC][j].op==GE)?">=":">"),
		   types[(constraints[CLASSC][j].a2)*2],
		   types[(constraints[CLASSC][j].a2)*2+1],
		   constraints[CLASSC][i].c2,
		   (constraints[CLASSC][j].alt==TRUE)?"*":"",
		   range->penet[constraints[CLASSC][j].c1-1][constraints[CLASSC][j].a1-DD][i],
		   ((constraints[CLASSC][j].op==EQ)?"==":(constraints[CLASSC][j].op==NE)?"!=":(constraints[CLASSC][j].op==GE)?">=":">"),
		   range->penet[constraints[CLASSC][j].c2-1][constraints[CLASSC][j].a2-DD][i]);
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
checkParams (ModelRange *range, int i)
{
  int j = 0;
#if FLASE
  char types[15]=DIRECTIVES;
#endif

  /* Scan through each constraint. */
  while (j < constcnt[PARAMC])
    {
#if FALSE
      fprintf (stderr, "%d: P%d %c%c %s P%d %c%c %s => ", i,
	       constraints[PARAMC][j].a1,
	       types[(constraints[PARAMC][j].a1)*2],
	       types[(constraints[PARAMC][j].a1)*2+1],
	       ((constraints[PARAMC][j].op==EQ)?"==":(constraints[PARAMC][j].op==NE)?"!=":(constraints[PARAMC][j].op==GE)?">=":">"),
	       constraints[PARAMC][j].a2,
	       types[(constraints[PARAMC][j].a2)*2],
	       types[(constraints[PARAMC][j].a2)*2+1],
	       (constraints[PARAMC][j].alt==TRUE)?"*":"");
#endif
      /* Relies on DD being the lowest numbered penetrance. */
      if (constraints[PARAMC][j].a1 >= DD && constraints[PARAMC][j].a2 >= DD)
	{
	  if ((constraints[PARAMC][j].op == EQ && 
	       range->param[0][constraints[PARAMC][j].a1-DD][constraints[PARAMC][j].p1-1][i] ==
	       range->param[0][constraints[PARAMC][j].a2-DD][constraints[PARAMC][j].p2-1][i]) ||
	      (constraints[PARAMC][j].op == NE && 
	       range->param[0][constraints[PARAMC][j].a1-DD][constraints[PARAMC][j].p1-1][i] !=
	       range->param[0][constraints[PARAMC][j].a2-DD][constraints[PARAMC][j].p2-1][i]) ||
	      (constraints[PARAMC][j].op == GT &&
	       range->param[0][constraints[PARAMC][j].a1-DD][constraints[PARAMC][j].p1-1][i] >
	       range->param[0][constraints[PARAMC][j].a2-DD][constraints[PARAMC][j].p2-1][i]) ||
	      (constraints[PARAMC][j].op == GE &&
	       range->param[0][constraints[PARAMC][j].a1-DD][constraints[PARAMC][j].p1-1][i] >=
	       range->param[0][constraints[PARAMC][j].a2-DD][constraints[PARAMC][j].p2-1][i]))
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
 * with dynamically.
 **********************************************************************/
int
checkClassParams (ModelRange *range, int i)
{
  int j = 0;
#if FLASE
  char types[15]=DIRECTIVES;
#endif

  /* Scan through each constraint. */
  while (j < constcnt[PARAMCLASSC])
    {
      /* Relies on DD being the lowest numbered penetrance. */
      if (constraints[PARAMCLASSC][j].a1 >= DD && constraints[PARAMCLASSC][j].a2 >= DD)
	{
	  if ((constraints[PARAMCLASSC][j].op == EQ && 
	       range->param[constraints[PARAMCLASSC][j].c1-1][constraints[PARAMCLASSC][j].a1-DD][constraints[PARAMCLASSC][j].p1-1][i] ==
	       range->param[constraints[PARAMCLASSC][j].c2-1][constraints[PARAMCLASSC][j].a2-DD][constraints[PARAMCLASSC][j].p2-1][i]) ||
	      (constraints[PARAMCLASSC][j].op == NE && 
	       range->param[constraints[PARAMCLASSC][j].c1-1][constraints[PARAMCLASSC][j].a1-DD][constraints[PARAMCLASSC][j].p1-1][i] !=
	       range->param[constraints[PARAMCLASSC][j].c2-1][constraints[PARAMCLASSC][j].a2-DD][constraints[PARAMCLASSC][j].p2-1][i]) ||
	      (constraints[PARAMCLASSC][j].op == GT &&
	       range->param[constraints[PARAMCLASSC][j].c1-1][constraints[PARAMCLASSC][j].a1-DD][constraints[PARAMCLASSC][j].p1-1][i] >
	       range->param[constraints[PARAMCLASSC][j].c2-1][constraints[PARAMCLASSC][j].a2-DD][constraints[PARAMCLASSC][j].p2-1][i]) ||
	      (constraints[PARAMCLASSC][j].op == GE &&
	       range->param[constraints[PARAMCLASSC][j].c1-1][constraints[PARAMCLASSC][j].a1-DD][constraints[PARAMCLASSC][j].p1-1][i] >=
	       range->param[constraints[PARAMCLASSC][j].c2-1][constraints[PARAMCLASSC][j].a2-DD][constraints[PARAMCLASSC][j].p2-1][i]))
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
 * Sort the model values so that we can apply constraints more
 * easily. Note: sortRange() is only applied prior to class expansion,
 * so we needn't worry, ever, about non-zero liability classes in
 * penet[][][]. Similarly, we needn't worry about liability classes or
 * allele counts in param[][][].
 **********************************************************************/
void
sortRange (ModelRange *range)
{
  int i;
  quicksort (range->gfreq, 0, range->ngfreq);
  if (range->penet[0])
    for (i = 0; i < NPENET(range->nalleles); i++)
      quicksort (range->penet[0][i], 0, penetcnt[i]);
  if (range->param)
    for (i = 0; i < range->npardim; i++)
      quicksort (range->param[0][0][i], 0, paramcnt[i]);
  for (i = 0; i < range->ngender; i++)
    quicksort (range->theta[i], 0, thetacnt[i]);
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
  while (i < len-1)
    {
      if (array[i] == array[i+1])
	{
	  for (j = i+1; j<len-1; j++)
	    array[j] = array[j+1];
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
 * classes or allele counts in param[][][]. */
void
uniqRange (ModelRange *range)
{
  int i;
  range->ngfreq = uniquify (range->gfreq, range->ngfreq);
  for (i = 0; i < NPENET(range->nalleles); i++)
    penetcnt[i] = uniquify (range->penet[0][i], penetcnt[i]);
  if (range->param)
    for (i = 0; i < range->npardim; i++)
      paramcnt[i] = uniquify (range->param[0][0][i], paramcnt[i]);
  for (i = 0; i < range->ngender; i++)
    thetacnt[i] = uniquify (range->theta[i], thetacnt[i]);
}

/**********************************************************************
 * Fully expand the penetrance and theta values in the model while
 * honoring all constraints. Using the fully expanded version allows
 * us to accurately determine the size of the space a priori, rather
 * than using the constraints to censor values as we go. The expansion
 * helps us avoid allocating excess space on disk.
 *
 * If you are doing QT or CT, you'll need to expand the parameters as
 * well. Expanded parameters will also increase the size of the
 * penetrance arrays.
 *
 * Note that liability class constraints are handled later.
 **********************************************************************/
void 
expandRange (ModelRange *range)
{
  int i, j, k, l, m;
  double **tmp;

  /* Start with thetas. Since these are not subject to liability class
   * constraints, we can finalize the thetas here by fully "factoring"
   * male and female thetas, if available, recording the number of
   * combinations in range->ntheta, and freeing thetacnt and thetamax. */

  /* Keep the current theta array values, and create a new theta
   * array. We'll allocate enough space for all combinations, even
   * though, in the end, constraint application may reduce this
   * number. The actualy number of combinations will be stored in
   * range->nthetas. */
  tmp = range->theta;
  range->theta = malloc (range->ngender * sizeof (double *));
  if (range->ngender == 1)
    {
      /* Easy case is that you have only one gender anyway. Since
       * there is only one dimension, there can be no constraints
       * worth enforcing. */
      range->theta[0] = tmp[0];
      range->ntheta = thetacnt[SEXAV];
    }
  else
    {
      /* If you have 2 genders, you need to factor their respective
       * values while checking constraints. */
      for (i = 0; i < range->ngender; i++)
	range->theta[i] = malloc ((thetacnt[SEXML]*thetacnt[SEXFM]) * sizeof (double));
      range->ntheta = 0;
      for (i = 0; i < thetacnt[SEXML]; i++)
	for (j = 0; j < thetacnt[SEXFM]; j++)
	  {
	    range->theta[SEXML][range->ntheta] = tmp[SEXML][i];
	    range->theta[SEXFM][range->ntheta] = tmp[SEXFM][j];
	    
	    if (checkThetas (range, range->ntheta))
	      range->ntheta++;
	  }
      /* Free old theta array, which you had stashed in tmp. */
      for (i = 1; i < range->ngender; i++)
	free (tmp[i]);
      free (tmp);
    }
#if FALSE
  /* SEGFAULTS WHEN YOU DO THESE!! */
  /* In either case, you're done with thetacnt and thetamax. */
  free(thetacnt);
  free(thetamax);
#endif
  /* Next, we expand the penetrances. The goal here is just to get the
   * combinations right, while ignoring the liability classes. It's a
   * stickier problem than the thetas, mostly because these values
   * will have to be expanded yet again when dealing with liability
   * classes, but also because the code must work for multiallelic
   * diseases (where there are more than 3 allele combinations, and,
   * therefore, more than 3 rows in the penetrance array).
   *
   * Stash penet[0] and penetcnt away temporarily. */
  tmp = range->penet[0];
#if FALSE
  free(range->penet); 
#endif
  /* Set up the penetrance array, ignoring the first dimension for
   * now. We know there will be PROD_i(penetcnt[i])^#allcombo values
   * in the expanded array, unless the constraints rule some out. 
   *
   * We'll keep track of the resulting number of combinations in
   * range->npenet. */
  range->npenet = 0;
  i = 1;
  for (j = 0; j < NPENET(range->nalleles); j++) 
    i = i*penetcnt[j];
  range->penet = malloc (sizeof (double **));
  range->penet[0] = malloc (NPENET(range->nalleles) * sizeof (double *));
  for (j = 0; j < NPENET(range->nalleles); j++)
    range->penet[0][j] = malloc (i * sizeof (double));

  /* OK, now populate the array. */
  for (k = 0; k < i; k++)
    {
      l = 1;
      for (m = 0; m < NPENET(range->nalleles); m++)
	{
	  range->penet[0][m][range->npenet] = 
	    tmp[m][((int) (k/l))%penetcnt[m]];
	  l = l*penetcnt[m];
	}
      /* Check the constraints. */
      if (checkPenets (range, range->npenet))
	range->npenet++;
    }

  /* OK, we're done with penetcnt and penetmax, since we're using
   * range->npenet to keep track of the numeber of combinations. We
   * can also ditch the stashed values for penet[][][]. */
#if FALSE
  for (i = 0; i < NPENET(range->nalleles); i++)
    free(tmp[i]);
  free (tmp);
  free(penetmax);
  free(penetcnt);
#endif

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
      /* Stash the original values somewhere. */
      tmp = range->param[0][0];
#if FALSE
      free(range->param[0]); 
      free(range->param);
#endif
      /* Since the resulting parameter array will be uniform, we can
       * allocate, in advance, the appropriate amount of memory. This
       * will be large! Hopefully, the constraints will limit how many
       * elements are actually in use, which will be stored in
       * range->nparam. */

      /* Let i be the max number of entries you might see for each
       * dimension, then raise it to the power of the dimensions. */
      i = pow(paramcnt[0], NPENET(range->nalleles));
      i = pow(i, range->npardim);
      
      /* Since the parameter vector is to be replicated across all of
       * the dimensions, they will all initially have the same number
       * of parameter values, although this will surely change as the
       * constraints are applied. Go ahead and allocate what you might
       * need. */
      range->param = malloc (sizeof (double ***));
      range->param[0] = malloc (NPENET(range->nalleles) * sizeof (double **));
      for (j = 0; j < NPENET(range->nalleles); j++)
	{
	  range->param[0][j] = malloc (range->npardim * sizeof (double *));
	  for (k = 0; k < paramcnt[0]; k++)
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
	    for (l =0; l < range->npardim; l++)
	      range->param[0][k][l][range->nparam] =
		tmp[l][((int)(j/(pow(paramcnt[l],(l+k*range->npardim)))))%paramcnt[l]]; 
	  /* Check the constraints. */
	  if (checkParams (range, range->nparam))
	    range->nparam++;
	}

      /* Done. Free copies of original parameters you'd stashed away,
       * including the paramcnt and parammax arrays; you won't be
       * needing them again, since the (uniform) number of entries
       * stored is given by range->nparam. */
#if FALSE
      free (paramcnt);
      free (parammax);
      for (i = 0; i < range->npardim; i++)
	free (tmp[i]);
      free (tmp);
#endif
    }
}

/**********************************************************************
 * Fully expand the penetrance and parameter values by liability class
 * while honoring any inter-class constraints. This only gets called
 * if we are using liability classes (i.e., modelRange->nlclass is
 * greater than 1). 
 **********************************************************************/
void 
expandClass (ModelRange *range)
{
  int i, j, k, l, m, n, o;
  double **tmp1;
  double ***tmp2;
#if FALSE
  double ****tmp3;
#endif  

  /* Penetrance expansion by liability class. Here, we want to factor
   * the existing penetrances over multiple liability classes.
   *
   * Stash a copy of the current penet array and its size.  */
  tmp1 = range->penet[0];
#if FALSE
  free (range->penet);
#endif
  i = range->npenet;
  
  /* Barring constraints, how many values might you have? */
  j = pow(range->npenet, range->nlclass);
  
  /* Allocate a new penet array structure. */
  range->penet = malloc (range->nlclass * sizeof (double **));
  for (k = 0; k < range->nlclass; k++)
    {
      range->penet[k] = malloc (NPENET(range->nalleles) * sizeof (double *));
      for (l = 0; l < NPENET(range->nalleles); l++)
	range->penet[k][l] = malloc (j * sizeof (double));
    }

  /* OK, now populate the array. */
  range->npenet=0;
  for (k = 0; k < j; k++)
    {
      l = 1;
      for (m = 0; m < range->nlclass; m++)
	{
	  for (n = 0; n < NPENET(range->nalleles); n++)
	    range->penet[m][n][range->npenet] = tmp1[n][(((int) (k/l))%i)];
	  l = l*i;
	}
      
      /* Check the class constraints. */
      if (checkClassPenets (range, range->npenet))
	range->npenet++;
    }
  /* Done. Free up the old copy of the array. */
#if FALSE
  for (i = 0; i < NPENET(range->nalleles); i++)
    free(tmp1[i]);
  free (tmp1);
#endif

  /* Ready to work on the param array, if necessary.
   *
   * Again, stash a copy of the current param array and its size. */
  if (range->param)
    {
      tmp2 = range->param[0];
      i = range->nparam;
      
      /* Barring constraints, how many values might you have? */
      j = pow(range->nparam, range->nlclass);
      
      /* Free the old param array structure and allocate a new one. */
#if FALSE
      free (range->param);
#endif
      range->param = malloc (range->nlclass * sizeof (double ***));
      for (k = 0; k < range->nlclass; k++)
	{
	  range->param[k] = malloc (NPENET(range->nalleles) * sizeof (double **));
	  for (l = 0; l < NPENET(range->nalleles); l++)
	    {
	      range->param[k][l] = malloc (range->npardim * sizeof (double *));
	      for (m = 0; m < range->npardim; m++)
		range->param[k][l][m] = malloc (j * sizeof (double));
	    }
	}

      /* OK, now populate the array. */
      range->nparam=0;
      for (k = 0; k < j; k++)
	{
	  l = 1;
	  for (m = 0; m < range->nlclass; m++)
	    {
	      for (n = 0; n < NPENET(range->nalleles); n++)
		for (o = 0; o < range->npardim; o++)
		  range->param[m][n][o][range->nparam] =
		    tmp2[n][o][((int) (k/l))%i];
	      l = l*i;
	    }
	  
	  /* Check the class constraints. */
	  if (checkClassParams (range, range->nparam))
	    range->nparam++;
	}
      
      /* Done. Free up the old copy of the array. */
#if FALSE
      for (i = 0; i < NPENET(range->nalleles); i++)
	{
	  for (j = 0; j < range->npardim; j++)
	    free(tmp2[i][j]);
	  free (tmp2[i]);
	}
      free (tmp2);
#endif

      /* YUNGUI: UNCOMMENT THIS BLOCK IF YOU WANT TO "FACTOR" THE QT
       * PARAMETERS (stdev) WITH THE PENETRANCE ARRAY (mean). LEAVE IT
       * OUT IF YOU ARE HAPPY TO LOOP OVER BOTH STRUCTURES; SINCE
       * THERE ARE NO CONSTRAINTS IMPOSED AT THIS FACOTRING STEP THERE
       * IS LITTLE REASON TO DO IT. YOU'LL ALSO NEED TO UNCOMENT
       * VARIABLE tmp3 DEFINITION ABOVE. */
#if FALSE
      /* OK, you're almost done. Just need to factor the parameter and
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
	  range->penet[k] = malloc (NPENET(range->nalleles) * sizeof (double *));
	  for (l = 0; l < NPENET(range->nalleles); l++)
	    range->penet[k][l] = malloc (j * sizeof (double));
	}
      range->param = malloc (range->nlclass * sizeof (double ***));
      for (k = 0; k < range->nlclass; k++)
	{
	  range->param[k] = malloc (NPENET(range->nalleles) * sizeof (double **));
	  for (l = 0; l < NPENET(range->nalleles); l++)
	    {
	      range->param[k][l] = malloc (range->npardim * sizeof (double *));
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
	    for (n = 0; n < NPENET(range->nalleles); n++)
	      {
		range->penet[m][n][range->npenet] = tmp2[m][n][(((int) (k/l))%i)];
		for (o = 0; o < range->npardim; o++)
		  range->param[m][n][o][range->npenet] =
		    tmp3[m][n][o][(((int) (k/l))%i)];
	      }
	  range->npenet++;
	}
      /* Make sure nparam returns the same number. */
      range->nparam = range->npenet;
#if FALSE
      /* Free up the old copies. */
      for (i = 0; i < nlclass; i++)
	{
	  for (j = 0; j < NPENET(range->nalleles); j++)
	    {
	      free(tmp2[i][j]);
	      for (k = 0; k < range->npardim; k++)
		free(tmp3[i][j][k]);
	      free(tmp3[i][j]);
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
 * Dump the model and constraints for debugging purposes.  The level
 * argument indicates if we are pre-expansion (level = 0),
 * post-expansion (level = 1), or post-class-expansion (level =
 * 2). This is ugly, but is only for debugging.
 **********************************************************************/
void
showRange (ModelRange *range, int level)
{
  int i, j, k, l;
  printf ("======================================================================\n");
  printf ("LEVEL %d MODEL\n", level);
  printf ("======================================================================\n");
  printf ("%d GF=", range->ngfreq);
  for (i = 0; i < range->ngfreq; i++)
    printf ("%3.2g ", range->gfreq[i]);
 
  if (level > 0)
    printf ("\n%d Thetas", range->ntheta);
  for (i = 0; i < range->ngender; i++)
    {
      if (level == 0)
	printf ("\n%d Theta[%d]=", (level==0?thetacnt[i]:range->ntheta),i);
      else
	printf ("\n  Theta[%d]=", i);
      for (j = 0; j < (level==0?thetacnt[i]:range->ntheta); j++)
	printf ("%3.2g ", range->theta[i][j]);
    }

  printf ("\n%d Penetrances", (level==0?penetcnt[0]:range->npenet));
  for (i = 0; i < (level==2?range->nlclass:1); i++)
    for (j = 0; j < NPENET(range->nalleles); j++)
      {
	printf ("\n  Penet[%d][%d]=", i, j);
	for (k = 0; k < (level==0?penetcnt[j]:range->npenet); k++)
	  printf ("%3.2g ", range->penet[i][j][k]);
      }

  if (range->param)
    {
      printf ("\n%d Parameters", ((level==2?range->nlclass:1)*
				  (level==0?paramcnt[0]:range->nparam)));
      /* for (i = 0; i < (level==2?range->nlclass:1); i++) */
      for (i = 0; i < (level==2?range->nlclass:1); i++)
	for (j = 0; j < (level==0?1:NPENET(range->nalleles)); j++)
	  for (k = 0; k < (level>0?range->npardim:1); k++)
	    {
	      printf ("\n  Param[%d][%d][%d]=", i, j, k);
	      for (l = 0; l < (level==0?paramcnt[k]:range->nparam); l++)
		printf ("%3.2g ", range->param[i][j][k][l]);
	    }
    }
  printf ("\n");
}
void
showConstraints ()
{
  int i, j;
  char types[15] = DIRECTIVES;

  printf ("======================================================================\n");
  printf ("CONSTRAINTS:\n");
  printf ("======================================================================\n");

  for (i = 0; i < 4; i++)
    {
      if (constcnt[i] == 0)
	continue;
      printf ("%s constraints:\n", ((i==3)?"PCC":((i==2)?"PC":((i==1)?"CC":"C"))));
      for (j = 0; j < constcnt[i]; j++)
	{
	  if (i == SIMPLE)
	    printf (" %c%c %s %c%c %s\n", 
		    types[constraints[i][j].a1*2],
		    types[constraints[i][j].a1*2+1],
		    ((constraints[i][j].op==EQ)?"==":(constraints[i][j].op==NE)?"!=":(constraints[i][j].op==GE)?">=":">"),
		    types[constraints[i][j].a2*2],
		    types[constraints[i][j].a2*2+1],
		    (constraints[i][j].alt==TRUE)?"|":"");
	  else if (i == CLASSC)
	    printf (" %c%c %d %s %c%c %d %s\n", 
		    types[constraints[i][j].a1*2],
		    types[constraints[i][j].a1*2+1],
		    constraints[i][j].c1,
		    ((constraints[i][j].op==EQ)?"==":(constraints[i][j].op==NE)?"!=":(constraints[i][j].op==GE)?">=":">"),
		    types[constraints[i][j].a2*2],
		    types[constraints[i][j].a2*2+1],
		    constraints[i][j].c2,
		    (constraints[i][j].alt==TRUE)?"|":"");
	  else if (i == PARAMC)
	    printf (" P%d %c%c %s P%d %c%c %s\n", 
		    constraints[i][j].p1,
		    types[constraints[i][j].a1*2],
		    types[constraints[i][j].a1*2+1],
		    ((constraints[i][j].op==EQ)?"==":(constraints[i][j].op==NE)?"!=":(constraints[i][j].op==GE)?">=":">"),
		    constraints[i][j].p2,
		    types[constraints[i][j].a2*2],
		    types[constraints[i][j].a2*2+1],
		    (constraints[i][j].alt==TRUE)?"|":"");
	  else if (i == PARAMCLASSC)
	    printf (" P%d %c%c %d %s P%d %c%c %d %s\n", 
		    constraints[i][j].p1,
		    types[constraints[i][j].a1*2],
		    types[constraints[i][j].a1*2+1],
		    constraints[i][j].c1,
		    ((constraints[i][j].op==EQ)?"==":(constraints[i][j].op==NE)?"!=":(constraints[i][j].op==GE)?">=":">"),
		    constraints[i][j].p2,
		    types[constraints[i][j].a2*2],
		    types[constraints[i][j].a2*2+1],
		    constraints[i][j].c2,
		    (constraints[i][j].alt==TRUE)?"|":"");
	}
    }
}
