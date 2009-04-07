#include <stdlib.h>
#include <string.h>
#include "model_range.h"
#include "../utils/utils.h"

#define CHUNKSIZE 64

/* Legal values for Constraint.type */
#define SIMPLE 0
#define CLASSC 1
#define PARAMC 2
#define PARAMCLASSC 3

/* Legal values for Constraint.op; these need correspond to the elements in op_strs */
#define EQ 0
#define NE 1
#define GT 2
#define GE 3

/* used by lookup_comparator(); these need to correspond to the EQ, NE, etc. symbols */
char *op_strs[] = { "==", "!=", ">", ">=", NULL };

/* used by lookup_modelparam(); these need to correspond to the PEN_DD, etc. symbols
 * defined in model_range.h
 */
char *mp_strs[] = { "", "", "", "", "", "", "", "", "", "DD", "Dd", "dD", "dd", NULL };


/* Global variables for tracking actual size of the various dynamically allocated
 * lists inside ModelRange. Only interesting during configuration parsing, so no
 * point cluttering up the ModelRange structure.
 */
static int maxgfreq=0;
static int maxtloc=0;
static int maxalpha=0;
static int maxdprime=0;
static int maxtheta[2]={0,0};       /* increase dimension if there's more than 2 genders (!) */
int *penetcnt;                      /* Array of number of penetrances */
int *penetmax;                      /* Array of max penetrances */
Constraint *constraints[4];         /* Array of constraints by type */
int constmax[4] = { 0, 0, 0, 0 };   /* Max constraints in array (for dynamic alloc) */
int constcnt[4] = { 0, 0, 0, 0 };   /* Current number of constraints in array */


/* We make an effort here to keep the list of trait loci sorted as we
 * insert elements, so we can guarantee uniqueness at the same time.
 */
void addTraitLocus (ModelRange * range, double val)
{
  int i = 0, j;
  
  /* Validate value. But trait loci may be negative! So no reasonable
   * validation is really possible. */
  /* KASSERT ((val >=0), "Bad trait locus %g; aborting.\n", val); */

  /* Initialize the structure if first access. */
  if (!range->tloc)
    range->ntloc = maxtloc = 0;
  /* Enlarge array if necessary. */
  if (range->ntloc == maxtloc) {
    range->tloc = realloc (range->tloc, (maxtloc + CHUNKSIZE) * sizeof (double));
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


void addGeneFreq (ModelRange * range, double val)
{
  /* Validate value. */
  KASSERT ((val >= 0
	    && val <= 1.0), "Bad gene frequency value %g; aborting.\n", val);

  /* Initialize the structure if first access. */
  if (!range->gfreq)
    range->ngfreq = maxgfreq = 0;
  /* Enlarge array if necessary. */
  if (range->ngfreq == maxgfreq) {
    range->gfreq = realloc (range->gfreq, (maxgfreq + CHUNKSIZE) * sizeof (double));
    maxgfreq = maxgfreq + CHUNKSIZE;
  }
  /* Add the element. */
  range->gfreq[range->ngfreq] = val;
  range->ngfreq++;
}


void addAlpha (ModelRange * range, double val)
{
  /* Validate value. */
  KASSERT ((val >= 0 && val <= 1.0), "Bad alpha value %g; aborting.\n", val);

  /* Initialize the structure if first access. */
  if (!range->alpha)
    range->nalpha = maxalpha = 0;
  /* Enlarge array if necessary. */
  if (range->nalpha == maxalpha) {
    range->alpha = realloc (range->alpha,
			    (maxalpha + CHUNKSIZE) * sizeof (double));
    maxalpha = maxalpha + CHUNKSIZE;
  }
  /* Add the element. */
  range->alpha[range->nalpha] = val;
  range->nalpha++;
}


void addDPrime (ModelRange * range, double val)
{
  /* Validate value. DPrime values must be between -1 and 1. */
  KASSERT ((val >= -1 && val <= 1), "Bad D prime value %g; aborting.\n", val);

  /* Initialize the structure if first access. */
  if (!range->dprime)
    range->ndprime = maxdprime = 0;
  /* Enlarge array if necessary. */
  if (range->ndprime == maxdprime) {
    range->dprime = realloc (range->dprime,
			     (maxdprime + CHUNKSIZE) * sizeof (double));
    maxdprime = maxdprime + CHUNKSIZE;
  }
  /* Add the element. */
  range->dprime[range->ndprime] = val;
  range->ndprime++;
}


void addTheta (ModelRange * range, int type, double val)
{
  /* Validate value. */
  KASSERT ((val >= 0 && val <= 0.5), "Bad theta value %g; aborting.\n", val);

  /* Initialize the structure if first access. */
  if (!range->theta) {
    range->ngender = 2;
    range->theta = malloc (range->ngender * sizeof (double *));
    range->theta[SEXML] = malloc (CHUNKSIZE * sizeof (double));
    range->theta[SEXFM] = malloc (CHUNKSIZE * sizeof (double));
    range->thetacnt = malloc (range->ngender * sizeof (int *));
    range->thetacnt[SEXML] = maxtheta[SEXML] = 0;
    range->thetacnt[SEXFM] = maxtheta[SEXFM] = 0;
  }

  /* Finally, add the thetas as specified. */
  /* Sex-specific thetas. Update one or both arrays separately. */
  if (type == THETA_AVG || type == THETA_FEMALE) {
    /* Enlarge array if necessary. */
    if (range->thetacnt[SEXFM] == maxtheta[SEXFM]) {
      range->theta[SEXFM] = realloc (range->theta[SEXFM],
				     (maxtheta[SEXFM] +
				      CHUNKSIZE) * sizeof (double));
      maxtheta[SEXFM] = maxtheta[SEXFM] + CHUNKSIZE;
    }
    /* Add the element. */
    range->theta[SEXFM][range->thetacnt[SEXFM]] = val;
    range->thetacnt[SEXFM]++;
  }
  /* Next, update the other array. */
  if (type == THETA_AVG || type == THETA_MALE) {
    /* Enlarge array if necessary. */
    if (range->thetacnt[SEXML] == maxtheta[SEXML]) {
      range->theta[SEXML] = realloc (range->theta[SEXML],
				     (maxtheta[SEXML] +
				      CHUNKSIZE) * sizeof (double));
      maxtheta[SEXML] = maxtheta[SEXML] + CHUNKSIZE;
    }
    /* Add the element. */
    range->theta[SEXML][range->thetacnt[SEXML]] = val;
    range->thetacnt[SEXML]++;
  }
}


/**********************************************************************
 * Recall that penetrances are initially specified as if there are no
 * liability classes. Later, if liability classes are in use, we
 * "expand" the penetrance array to its full three dimensions. So
 * addPenetrance() is only ever used in "preexpansion" mode, where the
 * first dimension of the array (corresponds to liability class) is
 * always 1.
 *
 * The penetrance array: range->penet[nlclass][nallele][nparam]
 **********************************************************************/
void addPenetrance (ModelRange * range, int type, double val)
{
  int i, j;

  /* Validate value. This is problematic because for DT we know values
   * should be between 0 and 1, but for QT/CT values given here
   * represent the mean of a distribution, and could be just about
   * anything. */
  /* KASSERT ((val >= 0 && val <= 1.0), "Bad penetrance value %g; aborting.\n", val); */

  /* Initialize the structure if first access. */
  if (!range->penet) {
    /* Initialize: remember, if you are a pre-expansion penet[][][]
     * array, the first dimension (liability class) will always have
     * a dimension of size 1. */
    range->penet = malloc (sizeof (double **));
    i = NPENET (range->nalleles);
    range->penet[0] = malloc (i * sizeof (double *));
    range->penetLimits = malloc (i * sizeof (double *)); /* Space for limits */
    for (j = 0; j < i; j++) {
      range->penet[0][j] = malloc (CHUNKSIZE * sizeof (double));
      range->penetLimits[j] = malloc (2 * sizeof (double)); /* A min and a max for each */
      range->penetLimits[j][0] = 999999999.00;
      range->penetLimits[j][1] = -999999999.00;
    }
    /* Remember, liability class is not a true "independent" index
     * in the sense that we are storing combinations across
     * liability classes; hence we only need NPENET() individual
     * values for the counters (indeed, we only need one of each,
     * but that would complicate things too much I suspect). */
    penetcnt = malloc (i * sizeof (int));
    penetmax = malloc (i * sizeof (int));
    for (j = 0; j < i; j++) {
      penetcnt[j] = 0;
      penetmax[j] = CHUNKSIZE;
    }
  }

  /* Next, add the penetrance to the appropriate location in the
   * penet[][][] array, which may entail allocating more space. Recall
   * type is an integer, where DD=00=0; Dd=01=1; dD=10=2; and dd=11=3. We
   * get the type by subtracting DD, the "base type" from the
   * integer count at invocation. */
  if (penetcnt[type] == penetmax[type]) {
    range->penet[0][type] = realloc (range->penet[0][type],
				     (penetmax[type] +
				      CHUNKSIZE) * sizeof (double));
    penetmax[type] = penetmax[type] + CHUNKSIZE;
  }
  /* Add the element. */
  range->penet[0][type][penetcnt[type]] = val;

  /* See if it's a raw maximum or minimum for this type (allele) */
  if (range->penetLimits[type][0] > val)
    range->penetLimits[type][0] = val;
  if (range->penetLimits[type][1] < val)
    range->penetLimits[type][1] = val;

  penetcnt[type]++;
}


/**********************************************************************
 * Add a constraint on penetrances or thetas. We're just going to keep
 * these in an array, which we scan when checking penetrance or theta
 * values. Clunky, but OK since it will only be done at setup. A
 * better way would perhaps be to keep theta and penetrance
 * constraints separately.
 **********************************************************************/
void addConstraint (int type, int a1, int c1, int p1,
		    int op, int a2, int c2, int p2, int disjunct)
{
  /* Check for meaningless constraints. TODO: do more of this! */
  KASSERT ((((a1 == THETA_MALE && a2 == THETA_FEMALE) ||
	     (a1 == THETA_FEMALE && a2 == THETA_MALE)) ||
	    (a1 == THRESHOLD && a2 == THRESHOLD) ||
	    (a1 >= PEN_DD && a2 >= PEN_DD && a1 <= PEN_dd && a2 <= PEN_dd)),

	   "Meaningless constraint %s %s %s %s; aborting.\n",
	   mp_strs[a1], op_strs[op], mp_strs[a2], (disjunct == TRUE) ? "*" : "");

  /* Allocate more space if necessary. */
  if (constmax[type] == constcnt[type]) {
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

  if (type == CLASSC || type == PARAMCLASSC) {
    constraints[type][constcnt[type]].c1 = c1;
    constraints[type][constcnt[type]].c2 = c2;
  }
  if (type == PARAMC || type == PARAMCLASSC) {
    constraints[type][constcnt[type]].p1 = p1;
    constraints[type][constcnt[type]].p2 = p2;
  }
  /* Increment constraint count. */
  constcnt[type]++;
}


/* Search the (short) list of legal comparators for the given string, return the
 * index if found, -1 otherwise. */
int lookup_comparator (char *str)
{
  int va=0;

  while (op_strs[va] != NULL) {
    if (strcmp (op_strs[va], str) == 0)
      return (va);
    va++;
  }
  return (-1);
}


/* Search the list of legal model parameters, for the given string, return the index
 * if found, -1 otherwise. Note that we depend here (as in soooo many other places) 
 * on DD, Dd, dD and dd being contiguous and in that order.
 */
int lookup_modelparam (char *str)
{
  int va=0;

  while (mp_strs[va] != NULL) {
    if ((va >= PEN_DD) && (va <= PEN_dd)) {
      if (strcmp (mp_strs[va], str) == 0)
	return (va);
    } else if (strcasecmp (mp_strs[va], str) == 0)
      return (va);
    va++;
  }
  return (-1);
}
