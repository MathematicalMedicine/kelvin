/* Copyright (C) 2009, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model_range.h"
#include "../utils/utils.h"

#define CHUNKSIZE 64

/* String representations of the legal values for Constraint.type. Only
 * for prettifying logging output.
 */
char *contype_strs[] = { "Simple", "Liability Class", "Trait Parameter", "Class/Parameter" };

/* String representations of the legal values for Constraint.op. Used by
 * lookup_comparator()
 */
char *op_strs[] = { "==", "!=", ">", ">=", NULL };

/* String representations of the legal values for Constraint.a1 and a2. Used by
 * used by lookup_modelparam()
 */
char *mp_strs[] = { "", "Theta", "", "", "Threshold", "DD", "Dd", "dD", "dd", NULL };

/* Global variables for tracking actual size of the various dynamically allocated
 * lists inside ModelRange. Only interesting during configuration parsing, so no
 * point cluttering up the ModelRange structure.
 */
static int maxgfreq=0;
static int maxafreq=0;
static int maxtloc=0;
static int maxalpha=0;
static int maxdprime=0;
static int maxtheta[2]={0,0};       /* increase dimension if there's more than 2 genders (!) */
static int maxtthresh=0;
static int *penetcnt=NULL;          /* Array of number of penetrances */
static int *penetmax=NULL;          /* Array of max penetrances */
static int *paramcnt=NULL;          /* Number of QT/CT parameters */
static int *parammax=NULL;          /* Max QT/CT parameters */

Constraint *constraints[4];         /* Array of constraints by type */
static int constmax[4] = { 0, 0, 0, 0 };   /* Max constraints in array (for dynamic alloc) */
static int constcnt[4] = { 0, 0, 0, 0 };   /* Current number of constraints in array */

/* Protoypes for private routines */
inline void swap (double *array, int i, int j);
void quicksort (double *array, int lo, int hi);
inline int uniquify (double *array, int len);


/* We make an effort here to keep the list of trait loci sorted as we
 * insert elements, so we can guarantee uniqueness at the same time.
 */
void addTraitLocus (ModelRange * range, double val)
{
  int i = 0, j;
  
  /* Validate value. But trait loci may be negative! So no reasonable
   * validation is really possible. */
  /* ASSERT ((val >=0), "Bad trait locus %g", val); */

  /* Initialize the structure if first access. */
  if (!range->tloc)
    range->ntloc = maxtloc = 0;
  /* Enlarge array if necessary. */
  if (range->ntloc == maxtloc) {
    REALCHOKE(range->tloc, (maxtloc + CHUNKSIZE) * sizeof (double), void *);
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
  ASSERT ((val > 0 && val < 1.0), "Bad gene frequency value %g", val);

  /* Initialize the structure if first access. */
  if (!range->gfreq)
    range->ngfreq = maxgfreq = 0;
  /* Enlarge array if necessary. */
  if (range->ngfreq == maxgfreq) {
    REALCHOKE(range->gfreq, (maxgfreq + CHUNKSIZE) * sizeof (double), void *);
    maxgfreq = maxgfreq + CHUNKSIZE;
  }
  /* Add the element. */
  range->gfreq[range->ngfreq] = val;
  range->ngfreq++;
}


void addAlleleFreq (ModelRange * range, double val)
{
  /* Validate value. */
  ASSERT ((val >= 0 && val <= 1.0), "Bad allele frequency value %g", val);

  /* Initialize the structure if first access. */
  if (!range->afreq)
    range->nafreq = maxafreq = 0;
  /* Enlarge array if necessary. */
  if (range->nafreq == maxafreq) {
    REALCHOKE(range->afreq, (maxafreq + CHUNKSIZE) * sizeof (double), void *);
    maxafreq = maxafreq + CHUNKSIZE;
  }
  /* Add the element. */
  range->afreq[range->nafreq] = val;
  range->nafreq++;
}


void addAlpha (ModelRange * range, double val)
{
  /* Validate value. */
  ASSERT ((val >= 0 && val <= 1.0), "Bad alpha value %g", val);

  /* Initialize the structure if first access. */
  if (!range->alpha)
    range->nalpha = maxalpha = 0;
  /* Enlarge array if necessary. */
  if (range->nalpha == maxalpha) {
    REALCHOKE(range->alpha, (maxalpha + CHUNKSIZE) * sizeof (double), void *);
    maxalpha = maxalpha + CHUNKSIZE;
  }
  /* Add the element. */
  range->alpha[range->nalpha] = val;
  range->nalpha++;
}


void addDPrime (ModelRange * range, double val)
{
  /* Validate value. DPrime values must be between -1 and 1. */
  ASSERT ((val >= -1 && val <= 1), "Bad D prime value %g", val);

  /* Initialize the structure if first access. */
  if (!range->dprime)
    range->ndprime = maxdprime = 0;
  /* Enlarge array if necessary. */
  if (range->ndprime == maxdprime) {
    REALCHOKE(range->dprime, (maxdprime + CHUNKSIZE) * sizeof (double), void *);
    maxdprime = maxdprime + CHUNKSIZE;
  }
  /* Add the element. */
  range->dprime[range->ndprime] = val;
  range->ndprime++;
}


void addTheta (ModelRange * range, int type, double val)
{
  /* Validate value. */
  ASSERT ((val >= 0 && val <= 0.5), "Bad theta value %g", val);

  /* Initialize the structure if first access. */
  if (!range->theta) {
    range->ngender = 2;
    MALCHOKE(range->theta, range->ngender * sizeof (double *), void *);
    range->theta[SEXML] = range->theta[SEXFM] = NULL;
    MALCHOKE(range->thetacnt, range->ngender * sizeof (int *), void *);
    range->thetacnt[SEXML] = maxtheta[SEXML] = 0;
    range->thetacnt[SEXFM] = maxtheta[SEXFM] = 0;
  }

  /* Put sex-average thetas and male thetas in the same place */
  if (type == THETA_AVG || type == THETA_MALE) {
    /* Enlarge array if necessary. */
    if (range->thetacnt[SEXML] == maxtheta[SEXML]) {
      REALCHOKE(range->theta[SEXML], (maxtheta[SEXML] + CHUNKSIZE) * sizeof (double), void *);
      maxtheta[SEXML] = maxtheta[SEXML] + CHUNKSIZE;
    }
    /* Add the element. */
    range->theta[SEXML][range->thetacnt[SEXML]] = val;
    range->thetacnt[SEXML]++;
  }

  /* Female thetas go someplace else */
  if (type == THETA_FEMALE) {
    /* Enlarge array if necessary. */
    if (range->thetacnt[SEXFM] == maxtheta[SEXFM]) {
      REALCHOKE(range->theta[SEXFM], (maxtheta[SEXFM] + CHUNKSIZE) * sizeof (double), void *);
      maxtheta[SEXFM] = maxtheta[SEXFM] + CHUNKSIZE;
    }
    /* Add the element. */
    range->theta[SEXFM][range->thetacnt[SEXFM]] = val;
    range->thetacnt[SEXFM]++;
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
  /* ASSERT ((val >= 0 && val <= 1.0), "Bad penetrance value %g", val); */

  /* Initialize the structure if first access. */
  if (!range->penet) {
    /* Initialize: remember, if you are a pre-expansion penet[][][]
     * array, the first dimension (liability class) will always have
     * a dimension of size 1. */
    MALCHOKE(range->penet, sizeof (double **), void *);
    i = NPENET (range->nalleles);
    MALCHOKE(range->penet[0], i * sizeof (double *), void *);
    MALCHOKE(range->penetLimits, i * sizeof (double *), void *);
    for (j = 0; j < i; j++) {
      MALCHOKE(range->penet[0][j], CHUNKSIZE * sizeof (double), void *);
      MALCHOKE(range->penetLimits[j], 2 * sizeof (double), void *);
      range->penetLimits[j][0] = 999999999.00;
      range->penetLimits[j][1] = -999999999.00;
    }
    /* Remember, liability class is not a true "independent" index
     * in the sense that we are storing combinations across
     * liability classes; hence we only need NPENET() individual
     * values for the counters (indeed, we only need one of each,
     * but that would complicate things too much I suspect). */
    MALCHOKE(penetcnt, i * sizeof (int), void *);
    MALCHOKE(penetmax, i * sizeof (int), void *);
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
    REALCHOKE(range->penet[0][type], (penetmax[type] + CHUNKSIZE) * sizeof (double), void *);
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
  ASSERT ((((a1 == THETA_MALE && a2 == THETA_FEMALE) ||
	     (a1 == THETA_FEMALE && a2 == THETA_MALE)) ||
	    (a1 == THRESHOLD && a2 == THRESHOLD) ||
	    (a1 >= PEN_DD && a2 >= PEN_DD && a1 <= PEN_dd && a2 <= PEN_dd)),
	   "Meaningless constraint %s %s %s %s",
	   mp_strs[a1], op_strs[op], mp_strs[a2], (disjunct == TRUE) ? "*" : "");

  /* Allocate more space if necessary. */
  if (constmax[type] == constcnt[type]) {
    REALCHOKE(constraints[type], (constmax[type] + CHUNKSIZE) * sizeof (Constraint), Constraint *);
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
void addParameter (ModelRange * range, int dim, double val)
{
  int i;

  /* First, if this is the first access to the structure, you must
   * initialize it. */
  if (!range->param) {
    /* Initialize: remember, if you are a pre-expansion
     * param[][][][] array, the first dimension (liability class)
     * will always have a dimension of size 1, and the second
     * dimension (allele) will also always have a dimension of size
     * 1, since we "share" values, at least before applying
     * constraints, across all allele combinations. */
    MALCHOKE(range->param, sizeof (double ***), void *);
    MALCHOKE(range->param[0], sizeof (double **), void *);
    MALCHOKE(range->param[0][0], range->npardim * sizeof (double *), void *);

    /* Initialize the paramLimits min/max values, too */
    range->paramLimits[0] = 999999999.00;
    range->paramLimits[1] = -999999999.00;

    for (i = 0; i < range->npardim; i++)
      MALCHOKE(range->param[0][0][i], CHUNKSIZE * sizeof (double), void *);

    /* The only "real" dimension of param[][][][] reflected in
     * paramcnt and parammax is the dim dimension, since that's the
     * only one that will vary at input time. Recall that all
     * parameters are shared across all alleles and liability
     * classes.  */
    MALCHOKE(paramcnt, range->npardim * sizeof (int), void *);
    MALCHOKE(parammax, range->npardim * sizeof (int), void *);
    for (i = 0; i < range->npardim; i++) {
      paramcnt[i] = 0;
      parammax[i] = CHUNKSIZE;
    }
  }

  /* Next, add the parameter to the appropriate location in the
   * param[][][][] array, which may entail allocating more space. */
  if (paramcnt[dim] == parammax[dim]) {
    REALCHOKE(range->param[0][0][dim], (parammax[dim] + CHUNKSIZE) * sizeof (double), void *);
    parammax[dim] = parammax[dim] + CHUNKSIZE;
  }
  /* Add the element. */
  range->param[0][0][dim][paramcnt[dim]] = val;
  paramcnt[dim]++;

  /* See if it's a raw maximum or minimum */
  if (range->paramLimits[0] > val)
    range->paramLimits[0] = val;
  if (range->paramLimits[1] < val)
    range->paramLimits[1] = val;
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
void addTraitThreshold (ModelRange * range, double val)
{
  /* Validate value. The trait threshold should be between the lowest
   * and highest means. But since trait thresholds are subject to
   * liability classes, we won't be able to impose this constraint
   * until later. */
  /* ASSERT ((val >=0), "Bad trait threshold %g", val); */

  /* Initialize the structure if first access. */
  if (!range->tthresh) {
    /* Initialize: remember, if you are a pre-expansion tthresh[][]
     * array, the first dimension (liability class) will always have
     * a dimension of size 1. */
    MALCHOKE(range->tthresh, sizeof (double *), void *);
    range->tthresh[0] = NULL;
    range->ntthresh = 0;
  }

  /* Enlarge array if necessary. */
  if (range->ntthresh == maxtthresh) {
    REALCHOKE(range->tthresh[0], (maxtthresh + CHUNKSIZE) * sizeof (double), void *);
    maxtthresh = maxtthresh + CHUNKSIZE;
  }
  /* Add the element. */
  range->tthresh[0][range->ntthresh] = val;
  range->ntthresh++;
}


/* We check on penetrance values for PEN_dD here, to avoid having to export
 * penetcnt, and instead force configuration parsing to call out to this
 * routine. Virtuous? Meh.
 */
int checkImprintingPenets (ModelRange *range, int imprinting)
{
  int i;

  if (imprinting) {
    /* If Imprinting is on, penetrance values for PEN_dD better have been specified */
    if (penetcnt[PEN_dD-PEN_DD] == 0)
      return (-1);
    return (0);

  } else {
    /* If imprinting is off, penetrance values for PEN_dD should not have been specified */
    if (penetcnt[PEN_dD-PEN_DD] != 0)
      return (-1);

    /* Otherwise, copy penetrances for PEN_Dd to PEN_dD */
    for (i=0; i<penetcnt[PEN_Dd-PEN_DD]; i++) 
      addPenetrance (range, PEN_dD-PEN_DD, range->penet[0][PEN_Dd-PEN_DD][i]);
    addConstraint (SIMPLE, PEN_dD, 0, 0, EQ, PEN_Dd, 0, 0, FALSE);
    
    /*
    MALCHOKE(range->penet[0][PEN_dD-PEN_DD], (penetmax[PEN_Dd-PEN_DD]) * sizeof (double), void *);
    penetcnt[PEN_dD-PEN_DD] = penetcnt[PEN_Dd-PEN_DD];
    range->penetLimits[PEN_dD-PEN_DD][0] = range->penetLimits[PEN_Dd-PEN_DD][0];
    range->penetLimits[PEN_dD-PEN_DD][1] = range->penetLimits[PEN_Dd-PEN_DD][1];
    for (i=0; i<penetcnt[PEN_Dd-PEN_DD]; i++) 
      range->penet[0][PEN_dD-PEN_DD][i] = range->penet[0][PEN_Dd-PEN_DD][i];
    */

    return (0);
  }
}


/* Under dynamic sampling and QT/CT, there will be one min/max pair of
 * "penetrance" values (mean or degrees of freedom), either specified by
 * the user or filled in as defaults. These are stored as the penetrance
 * for the DD trait genotype. Those values need to be copied to all the
 * other trait genotypes.
 */
void duplicatePenets (ModelRange *range, int imprinting)
{
  addPenetrance (range, PEN_Dd-PEN_DD, range->penet[0][PEN_DD-PEN_DD][0]);
  addPenetrance (range, PEN_Dd-PEN_DD, range->penet[0][PEN_DD-PEN_DD][1]);
  if (imprinting) {
    addPenetrance (range, PEN_dD-PEN_DD, range->penet[0][PEN_DD-PEN_DD][0]);
    addPenetrance (range, PEN_dD-PEN_DD, range->penet[0][PEN_DD-PEN_DD][1]);
  }
  addPenetrance (range, PEN_dd-PEN_DD, range->penet[0][PEN_DD-PEN_DD][0]);
  addPenetrance (range, PEN_dd-PEN_DD, range->penet[0][PEN_DD-PEN_DD][1]);
  return;
}


/**********************************************************************
 * Return TRUE if the ith thetas satisfy all theta constraints, else
 * FALSE. The general idea is to scan through each constraint, exiting
 * whenever you encounter an unsatisfied one (that does not have a
 * disjunct). If there is a disjunct, then step through those until
 * you find one that satisfies the constraint; if you get to the end
 * of the line without satisfying the disjunct, you fail.
 **********************************************************************/
int checkThetas (ModelRange * range, int i)
{
  int j = 0;
  Constraint *conptr = constraints[SIMPLE];

  /* Scan through each constraint. */
  while (j < constcnt[SIMPLE]) {
    /* Relies on THETA_MALE being the lowest numbered penetrance. */
    if ((conptr[j].a1 == THETA_MALE && conptr[j].a2 == THETA_FEMALE) ||
	(conptr[j].a1 == THETA_FEMALE && conptr[j].a2 == THETA_MALE)) {
      if ((conptr[j].op == EQ &&
	   range->theta[conptr[j].a1 - THETA_MALE][i] ==
	   range->theta[conptr[j].a2 - THETA_MALE][i]) ||
	  (conptr[j].op == NE &&
	   range->theta[conptr[j].a1 - THETA_MALE][i] !=
	   range->theta[conptr[j].a2 - THETA_MALE][i]) ||
	  (conptr[j].op == GT &&
	   range->theta[conptr[j].a1 - THETA_MALE][i] >
	   range->theta[conptr[j].a2 - THETA_MALE][i]) ||
	  (conptr[j].op == GE &&
	   range->theta[conptr[j].a1 - THETA_MALE][i] >=
	   range->theta[conptr[j].a2 - THETA_MALE][i])) {
	/* Satisfied; skip other disjuncts. */
	while (conptr[j].alt == TRUE)
	  j++;
      } else if (conptr[j].alt == FALSE) {
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
int checkPenets (ModelRange * range, int i)
{
  int j = 0;
  Constraint *conptr = constraints[SIMPLE];

  /* Scan through each constraint. */
  while (j < constcnt[SIMPLE]) {
    if ((conptr[j].a1 < PEN_DD) || (conptr[j].a1 > PEN_dd) ||
	(conptr[j].a2 < PEN_DD) || (conptr[j].a2 > PEN_dd)) {
      /* Skip non-penetrance constraints */
      j++;
      continue;
    } 
    /* Relies on PEN_DD being the lowest numbered penetrance. */
    if ((conptr[j].op == EQ &&
	 range->penet[0][conptr[j].a1 - PEN_DD][i] ==
	 range->penet[0][conptr[j].a2 - PEN_DD][i]) ||
	(conptr[j].op == NE &&
	 range->penet[0][conptr[j].a1 - PEN_DD][i] !=
	 range->penet[0][conptr[j].a2 - PEN_DD][i]) ||
	(conptr[j].op == GT &&
	 range->penet[0][conptr[j].a1 - PEN_DD][i] >
	 range->penet[0][conptr[j].a2 - PEN_DD][i]) ||
	(conptr[j].op == GE &&
	 range->penet[0][conptr[j].a1 - PEN_DD][i] >=
	 range->penet[0][conptr[j].a2 - PEN_DD][i])) {
      /* Satisfied; skip other disjuncts. */
      while (conptr[j].alt == TRUE)
	j++;
    } else if (conptr[j].alt == FALSE) {
      return (FALSE);
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
int checkClassPenets (ModelRange * range, int i)
{
  int j = 0;
  Constraint *conptr = constraints[CLASSC];

  /* Scan through each constraint: inter-class constraints are on
   * penetrances or trait thresholds, not, e.g., gene frequencies or
   * thetas. */
  while (j < constcnt[CLASSC]) {
    if ((conptr[j].a1 < PEN_DD) || (conptr[j].a1 > PEN_dd) ||
	(conptr[j].a2 < PEN_DD) || (conptr[j].a2 > PEN_dd)) {
      /* Skip non-penetrace constraints */
      j++;
      continue;
    } 
    if ((conptr[j].op == EQ &&
	 range->penet[conptr[j].c1-1][conptr[j].a1 - PEN_DD][i] ==
	 range->penet[conptr[j].c2-1][conptr[j].a2 - PEN_DD][i]) ||
	(conptr[j].op == NE &&
	 range->penet[conptr[j].c1-1][conptr[j].a1 - PEN_DD][i] !=
	 range->penet[conptr[j].c2-1][conptr[j].a2 - PEN_DD][i]) ||
	(conptr[j].op == GT &&
	 range->penet[conptr[j].c1-1][conptr[j].a1 - PEN_DD][i] >
	 range->penet[conptr[j].c2-1][conptr[j].a2 - PEN_DD][i]) ||
	(conptr[j].op == GE &&
	 range->penet[conptr[j].c1-1][conptr[j].a1 - PEN_DD][i] >=
	 range->penet[conptr[j].c2-1][conptr[j].a2 - PEN_DD][i])) {
      /* Satisfied; skip other disjuncts. */
      while (conptr[j].alt == TRUE)
	j++;
    } else if (conptr[j].alt == FALSE) {
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
int checkParams (ModelRange * range, int i)
{
  int j = 0;
  Constraint *conptr = constraints[PARAMC];

  /* Scan through each constraint. */
  while (j < constcnt[PARAMC]) {
    if ((conptr[j].a1 < PEN_DD) || (conptr[j].a1 > PEN_dd) ||
	(conptr[j].a2 < PEN_DD) || (conptr[j].a2 > PEN_dd)) {
      /* Skip non-penetrace constraints */
      j++;
      continue;
    } 
    if ((conptr[j].op == EQ &&
	 range->param[0][conptr[j].a1 - PEN_DD][conptr[j].p1 - 1][i] ==
	 range->param[0][conptr[j].a2 - PEN_DD][conptr[j].p2 - 1][i]) ||
	(conptr[j].op == NE &&
	 range->param[0][conptr[j].a1 - PEN_DD][conptr[j].p1 - 1][i] !=
	 range->param[0][conptr[j].a2 - PEN_DD][conptr[j].p2 - 1][i]) ||
	(conptr[j].op == GT && 
	 range->param[0][conptr[j].a1 - PEN_DD][conptr[j].p1 - 1][i] >
	 range->param[0][conptr[j].a2 - PEN_DD][conptr[j].p2 - 1][i]) ||
	(conptr[j].op == GE &&
	 range->param[0][conptr[j].a1 - PEN_DD][conptr[j].p1 - 1][i] >=
	 range->param[0][conptr[j].a2 - PEN_DD][conptr[j].p2 - 1][i])) {
      /* Satisfied; skip other disjuncts. */
      while (conptr[j].alt == TRUE)
	j++;
    } else if (conptr[j].alt == FALSE) {
      return (FALSE);
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
int checkClassParams (ModelRange * range, int i)
{
  int j = 0;
  Constraint *conptr = constraints[PARAMCLASSC];

  /* Scan through each constraint. */
  while (j < constcnt[PARAMCLASSC]) {
    if ((conptr[j].a1 < PEN_DD) || (conptr[j].a1 > PEN_dd) ||
	(conptr[j].a2 < PEN_DD) || (conptr[j].a2 > PEN_dd)) {
      /* Skip non-penetrace constraints */
      j++;
      continue;
    } 
    
    /* Relies on DD being the lowest numbered penetrance. */
    if ((conptr[j].op == EQ &&
	 range->param[conptr[j].c1 - 1][conptr[j].a1 - PEN_DD][conptr[j].p1 - 1][i] ==
	 range->param[conptr[j].c2 - 1][conptr[j].a2 - PEN_DD][conptr[j].p2 - 1][i]) ||
	(conptr[j].op == NE &&
	 range->param[conptr[j].c1 - 1][conptr[j].a1 - PEN_DD][conptr[j].p1 - 1][i] !=
	 range->param[conptr[j].c2 - 1][conptr[j].a2 - PEN_DD][conptr[j].p2 - 1][i]) || 
	(conptr[j].op == GT &&
	 range->param[conptr[j].c1 - 1][conptr[j].a1 - PEN_DD][conptr[j].p1 - 1][i] >
	 range->param[conptr[j].c2 - 1][conptr[j].a2 - PEN_DD][conptr[j].p2 - 1][i]) ||
	(conptr[j].op == GE &&
	 range->param[conptr[j].c1 - 1][conptr[j].a1 - PEN_DD][conptr[j].p1 - 1][i] >=
	 range->param[conptr[j].c2 - 1][conptr[j].a2 - PEN_DD][conptr[j].p2 - 1][i])) {
      /* Satisfied; skip other disjuncts. */
      while (conptr[j].alt == TRUE)
	j++;
    } else if (conptr[j].alt == FALSE) {
      return (FALSE);
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
int checkClassThreshold (ModelRange * range, int i)
{
  int j = 0;

  /* Scan through each constraint: inter-class constraints are on
   * penetrances or trait thresholds, not, e.g., gene frequencies or
   * thetas. */
  while (j < constcnt[CLASSC]) {
    if (constraints[CLASSC][j].a1 != THRESHOLD || constraints[CLASSC][j].a2 != THRESHOLD) {
      /* Skip non threshold constraints. */
      j++;
      continue;
    } else if ((constraints[CLASSC][j].op == EQ &&
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
		range->tthresh[constraints[CLASSC][j].c2 - 1][i])) {
      /* Satisfied; skip other disjuncts. */
      while (constraints[CLASSC][j].alt == TRUE)
	j++;
    } else if (constraints[CLASSC][j].alt == FALSE) {
      return (FALSE);
    }
    j++;
  }
  return (TRUE);
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
void expandRange (ModelRange *range)
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
  if ((tmp2 = range->theta)) {
    MALCHOKE(range->theta, range->ngender * sizeof (double *), void *);
    if (range->thetacnt[SEXFM] == 0) {
      /* Easy case is that you have only one gender anyway. Since
       * there is only one dimension, there can be no constraints
       * worth enforcing. Since addTheta() doesn't fill in the female
       * theta list when sex-averaged thetas are added, we duplicate 
       * the sex-averaged/male theta list real fast. */
      
      range->theta[SEXML] = tmp2[SEXML];
      MALCHOKE(range->theta[SEXFM], range->thetacnt[SEXML] * sizeof (double), void *);
      for (i = 0; i < range->thetacnt[SEXML]; i++)
	range->theta[SEXFM][i] = range->theta[SEXML][i];
      range->thetacnt[SEXFM] = range->thetacnt[SEXML];
      range->ntheta = range->thetacnt[SEXML];

    } else {
      /* sex specific */
      /* If you have 2 genders, you need to factor their respective
       * values while checking constraints. */
      for (i = 0; i < range->ngender; i++)
	MALCHOKE(range->theta[i], (range->thetacnt[SEXML] * range->thetacnt[SEXFM]) * sizeof (double), void *);
      range->ntheta = 0;
      for (i = 0; i < range->thetacnt[SEXML]; i++)
	for (j = 0; j < range->thetacnt[SEXFM]; j++) {
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
  if ((tmp3 = range->penet) != NULL) {
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
    MALCHOKE(range->penet, sizeof (double **), void *);
    MALCHOKE(range->penet[0], NPENET (range->nalleles) * sizeof (double *), void *);
    for (j = 0; j < NPENET (range->nalleles); j++)
      MALCHOKE(range->penet[0][j], i * sizeof (double), void *);
    
    /* OK, now populate the array. */
    for (k = 0; k < i; k++) {
      l = 1;
      for (m = 0; m < NPENET (range->nalleles); m++) {
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
  }

  /* Finally, expand the parameters for QT or CT. Recall param[][][][]
   * at this point only has 2 real indexes: dimension and the number
   * of values. The goal here is to factor these into 3 real indexes,
   * including allele combinations, while respecting the
   * constraints. We'll worry about liability classes later.
   *
   * Note that the resulting param[0][][][] array will be of uniform
   * dimension, that being ((#value)^#dim)^#allcombos). */
  if (range->param) {
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
    MALCHOKE(range->param, sizeof (double ***), void *);
    MALCHOKE(range->param[0], NPENET (range->nalleles) * sizeof (double **), void *);
    for (j = 0; j < NPENET (range->nalleles); j++) {
      MALCHOKE(range->param[0][j], range->npardim * sizeof (double *), void *);
      /* Corrected by Yungui: used to read:
       *  for (k = 0; k < paramcnt[0]; k++) */
      for (k = 0; k < range->npardim; k++)
	MALCHOKE(range->param[0][j][k], i * sizeof (double), void *);
    }

    /* Time to populate the range->param array. Recall the first
     * index will always be 0 since we're not yet dealing with
     * liability classes. We'll use range->nparam to store the
     * number of valid combinations. */
    range->nparam = 0;
    for (j = 0; j < i; j++) {
      for (k = 0; k < NPENET (range->nalleles); k++)
	for (l = 0; l < range->npardim; l++)
	  range->param[0][k][l][range->nparam] =
	    tmp2[l][((int) (j / (pow (paramcnt[l], (l + k * range->npardim))))) % paramcnt[l]];
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
  }
}


/**********************************************************************
 * Fully expand the threshold values by liability class while honoring
 * any inter-class constraints. This only gets called if we are using
 * liability classes (i.e., modelRange->nlclass is greater than 1).
 **********************************************************************/
void expandClassThreshold (ModelRange * range)
{
  int i, j, k, l, m;
  double *tmp1;
  double **tmp2;

  /* Threshold expansion by liability class. */
  if (range->tthresh) {
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
    MALCHOKE(range->tthresh, range->nlclass * sizeof (double *), void *);
    for (k = 0; k < range->nlclass; k++)
      MALCHOKE(range->tthresh[k], j * sizeof (double), void *);

    /* OK, now populate the array. */
    range->ntthresh = 0;
    for (k = 0; k < j; k++) {
      l = 1;
      for (m = 0; m < range->nlclass; m++) {
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
}


/**********************************************************************
 * Fully expand the penetrance and parameter values by liability class
 * while honoring any inter-class constraints. This only gets called
 * if we are using liability classes (i.e., modelRange->nlclass is
 * greater than 1).
 **********************************************************************/
void expandClassPenet (ModelRange * range)
{
  int i, j, k, l, m, n, o;
  double **tmp2;
  double ***tmp3;
  double ****tmp4;

  /* Penetrance expansion by liability class. Here, we want to factor
   * the existing penetrances over multiple liability classes.
   *
   * Stash pointers to the current penet array and its size. Recall
   * range->penet[0] is the only dimension originally allocated; we'll
   * use the values stored there during expansion, producing a new
   * multidimensional range->penet array. Be sure to free everything
   * properly when done. */
  if ((tmp3 = range->penet) != NULL) {
    tmp2 = range->penet[0];
    i = range->npenet;
    
    /* Barring constraints, how many values might you have? */
    j = pow (range->npenet, range->nlclass);
    
    /* Allocate a new penet array structure. */
    MALCHOKE(range->penet, range->nlclass * sizeof (double **), void *);
    for (k = 0; k < range->nlclass; k++) {
      MALCHOKE(range->penet[k], NPENET (range->nalleles) * sizeof (double *), void *);
      for (l = 0; l < NPENET (range->nalleles); l++)
	MALCHOKE(range->penet[k][l], j * sizeof (double), void *);
    }
    
    /* OK, now populate the array. */
    range->npenet = 0;
    for (k = 0; k < j; k++) {
      l = 1;
      for (m = 0; m < range->nlclass; m++) {
	for (n = 0; n < NPENET (range->nalleles); n++)
	  range->penet[m][n][range->npenet] = tmp2[n][(((int) (k / l)) % i)];
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
  }

  /* Ready to work on the param array, if necessary.
   *
   * Again, stash a copy of the current param array and its size. */
  if (range->param) {
    tmp4 = range->param;
    tmp3 = range->param[0];
    i = range->nparam;

    /* Barring constraints, how many values might you have? */
    j = pow (range->nparam, range->nlclass);

    /* Allocate a new param array structure. */
    MALCHOKE(range->param, range->nlclass * sizeof (double ***), void *);
    for (k = 0; k < range->nlclass; k++) {
      MALCHOKE(range->param[k], NPENET (range->nalleles) * sizeof (double **), void *);
      for (l = 0; l < NPENET (range->nalleles); l++) {
	MALCHOKE(range->param[k][l], range->npardim * sizeof (double *), void *);
	for (m = 0; m < range->npardim; m++)
	  MALCHOKE(range->param[k][l][m], j * sizeof (double), void *);
      }
    }

    /* OK, now populate the array. */
    range->nparam = 0;
    for (k = 0; k < j; k++) {
      l = 1;
      for (m = 0; m < range->nlclass; m++) {
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
    for (i = 0; i < NPENET (range->nalleles); i++) {
      for (j = 0; j < range->npardim; j++)
	free (tmp3[i][j]);
      free (tmp3[i]);
    }
    free (tmp4);

  }
}


/**

  Either find or build and return a lambda array for the specified
  disease/marker allele count m and marker allele count n.

  Expand the provided dprime array into an appropriate generic lambda
  array. Unlike the other range expansion functions, this function is
  not called upfront, but rather as needed when operating under
  linkage disequilibrium. The parameters m and n are the number of
  alleles for the two loci under consideration.

  Since the lambda array is independent of actual allele frequencies,
  and only depends upon the number of alleles for each locus, it can
  frequently be reused, so we &&&
  preserve the lambda arrays produced from the dprime values for each
  value of n and m for reuse in modelRange->lambdas, a vector of
  modelRange->nlambdas structures of type lambdaCell distinguished by
  the values of m and n. That way, we can retrieve the appropriate
  array if its already been generated, otherwise, we build the array
  from the dprimes, cache it, and return a pointer to it.
 
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

  /* First, see if an m X n array of lambdas already exists. */
  while (i < range->nlambdas) {
    if (range->lambdas[i].m == m && range->lambdas[i].n == n) {
      /* return (range->lambdas[i].lambda); */
      return (&range->lambdas[i]);
    }
    i++;
  }

  /* OK, no matching m X n array found. Check to make sure there's
   * room for a new one. */
  if (range->nlambdas == range->maxnlambdas) {
    REALCHOKE(range->lambdas, (range->maxnlambdas + CHUNKSIZE) * sizeof (LambdaCell), void *);
    range->maxnlambdas = range->maxnlambdas + CHUNKSIZE;
  }
  /* Create the new entry. */
  range->lambdas[range->nlambdas].m = m;
  range->lambdas[range->nlambdas].n = n;
  range->lambdas[range->nlambdas].ndprime = l;
  MALCHOKE(range->lambdas[range->nlambdas].lambda, l * sizeof (double **), double ***);
  CALCHOKE(range->lambdas[range->nlambdas].impossibleFlag, (size_t) 1, l * sizeof (int), int *);
  MALCHOKE(range->lambdas[range->nlambdas].haploFreq, l * sizeof (double **), double ***);
  MALCHOKE(range->lambdas[range->nlambdas].DValue, l * sizeof (double **), double ***);

  for (i = 0; i < l; i++) {
    MALCHOKE(range->lambdas[range->nlambdas].lambda[i], (m - 1) * sizeof (double *), double **);

    for (j = 0; j < (m - 1); j++)
      MALCHOKE(range->lambdas[range->nlambdas].lambda[i][j], (n - 1) * sizeof (double), double *);

    MALCHOKE(range->lambdas[range->nlambdas].haploFreq[i], m * sizeof (double *),double **);
    MALCHOKE(range->lambdas[range->nlambdas].DValue[i], m * sizeof (double *),double **);

    for (j = 0; j < m; j++) {
      MALCHOKE(range->lambdas[range->nlambdas].haploFreq[i][j], n * sizeof (double), double *);
      MALCHOKE(range->lambdas[range->nlambdas].DValue[i][j], n * sizeof (double), double *);
    }

  }
  range->nlambdas++;

  /* Now populate the values in the appropriate array. */
  for (i = 0; i < l; i++)
    for (j = 0; j < (m - 1); j++)
      for (k = 0; k < (n - 1); k++) {
	range->lambdas[range->nlambdas - 1].lambda[i][j][k] =
	  range->dprime[((int) (i / pow (range->ndprime, k))) %
			range->ndprime];
	//	fprintf (stderr, "lambda[i=%d][j=%d][k=%d] is %g\n", i, j, k,
	//		 range->lambdas[range->nlambdas - 1].lambda[i][j][k]);
      }
  /* Return a pointer to the appropriate 3 dimensional lambda array. */
  /* return (range->lambdas[range->nlambdas-1].lambda); */
  return (&range->lambdas[range->nlambdas - 1]);
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
void sortRange (ModelRange * range)
{
  int i;

  if (range->gfreq)
    quicksort (range->gfreq, 0, range->ngfreq);
  if ((range->penet) && (range->penet[0]))
    for (i = 0; i < NPENET (range->nalleles); i++)
      quicksort (range->penet[0][i], 0, penetcnt[i]);
  if (range->param)
    for (i = 0; i < range->npardim; i++)
      quicksort (range->param[0][0][i], 0, paramcnt[i]);
  if (range->theta)
    for (i = 0; i < range->ngender; i++)
      quicksort (range->theta[i], 0, range->thetacnt[i]);
  if (range->tthresh)
    quicksort (range->tthresh[0], 0, range->ntthresh);
  if (range->alpha)
    quicksort (range->alpha, 0, range->nalpha);
  if (range->afreq)
    quicksort (range->afreq, 0, range->nafreq);
  if (range->dprime)
    quicksort (range->dprime, 0, range->ndprime);
}


/* Uniquify the model's arrays. This is only done prior to expanding
 * the ranges, so we needn't worry, ever, about non-zero liability
 * classes in penet[][][]. Similarly, we needn't worry about liability
 * classes or allele counts in param[][][]. 
 *
 * Note that range->tloc is handled differently, because it may be
 * called back (in the event of a TM parameter) and so sorting and
 * uniquifying must be done on the fly for this array only. */
void uniqRange (ModelRange * range)
{
  int i;

  if (range->gfreq)
    range->ngfreq = uniquify (range->gfreq, range->ngfreq);
  if (range->penet)
    for (i = 0; i < NPENET (range->nalleles); i++)
      penetcnt[i] = uniquify (range->penet[0][i], penetcnt[i]);
  if (range->param)
    for (i = 0; i < range->npardim; i++)
      paramcnt[i] = uniquify (range->param[0][0][i], paramcnt[i]);
  if (range->theta)
    for (i = 0; i < range->ngender; i++)
      range->thetacnt[i] = uniquify (range->theta[i], range->thetacnt[i]);
  if (range->tthresh)
    range->ntthresh = uniquify (range->tthresh[0], range->ntthresh);
  if (range->alpha)
    range->nalpha = uniquify (range->alpha, range->nalpha);
}


/**********************************************************************
 * Dump the model and constraints for debugging purposes.  The level
 * argument indicates if we are pre-expansion (level = 0),
 * post-expansion (level = 1), or post-class-expansion (level =
 * 2). This is ugly, but is only used for debugging so it needn't be
 * excessively pretty!
 **********************************************************************/
void showRange (ModelRange * range, ModelType * type, int level)
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

  if (type->trait == DT) {
    if (level > 0)
      printf ("\n%d Thetas", range->ntheta);
    for (i = 0; i < range->ngender; i++) {
      if (level == 0)
	printf ("\n%d Theta[%d]=", (level == 0 ? range->thetacnt[i] : range->ntheta), i);
      else
	printf ("\n  Theta[%d]=", i);
      for (j = 0; j < (level == 0 ? range->thetacnt[i] : range->ntheta); j++)
	printf ("%3.2g ", range->theta[i][j]);
    }
  } else {
    if (range->tloc) {
      printf ("\n%d Trait Loci", range->ntloc);
      printf ("\n  Trait Loci=");
      for (i = 0; i < range->ntloc; i++)
	printf ("%3.2g ", range->tloc[i]);
    }
    if (range->tthresh) {
      printf ("\n%d Trait Thresholds", range->ntthresh);
      for (i = 0; i < (level == 2 ? range->nlclass : 1); i++) {
	printf ("\n  Trait Threshold[%d]=", i);
	for (j = 0; j < range->ntthresh; j++)
	  printf ("%3.2g ", range->tthresh[i][j]);
      }
    }
  }
  
  printf ("\n%d Penetrances", (level == 0 ? penetcnt[0] : range->npenet));
  for (i = 0; i < (level == 2 ? range->nlclass : 1); i++)
    for (j = 0; j < NPENET (range->nalleles); j++) {
      printf ("\n  Penet[%d][%d]=", i, j);
      for (k = 0; k < (level == 0 ? penetcnt[j] : range->npenet); k++)
	printf ("%3.2g ", range->penet[i][j][k]);
    }
  
  if (range->param) {
    printf ("\n%d Parameters", ((level == 2 ? range->nlclass : 1) *
				(level == 0 ? paramcnt[0] : range->nparam)));
    /* for (i = 0; i < (level==2?range->nlclass:1); i++) */
    for (i = 0; i < (level == 2 ? range->nlclass : 1); i++)
      for (j = 0; j < (level == 0 ? 1 : NPENET (range->nalleles)); j++)
	for (k = 0; k < (level > 0 ? range->npardim : 1); k++) {
	  printf ("\n  Param[%d][%d][%d]=", i, j, k);
	  for (l = 0; l < (level == 0 ? paramcnt[k] : range->nparam); l++)
	    printf ("%3.2g ", range->param[i][j][k][l]);
	}
  }
  
  if (range->afreq) {
    printf ("\n%d AF=", range->nafreq);
    for (i = 0; i < range->nafreq; i++)
      printf ("%3.2g ", range->afreq[i]);
  }
  
  if (range->dprime) {
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
  if (range->alpha) {
    printf ("\n%d Alphas", range->nalpha);
    printf ("\n  Alpha=");
    for (i = 0; i < range->nalpha; i++)
      printf ("%3.2g ", range->alpha[i]);
  }
  printf ("\n");
}


void showConstraints ()
{
  int i, j;

  printf 
    ("======================================================================\n");
  printf ("CONSTRAINTS:\n");
  printf
    ("======================================================================\n");

  for (i = 0; i < 4; i++) {
    if (constcnt[i] == 0)
      continue;
    printf ("%s constraints:\n", contype_strs[i]);
    for (j = 0; j < constcnt[i]; j++) {
      if (i == SIMPLE)
	printf (" %s %s %s %s\n", mp_strs[constraints[i][j].a1],
		op_strs[constraints[i][j].op], mp_strs[constraints[i][j].a2],
		(constraints[i][j].alt == TRUE) ? "|" : "");
      
      else if (i == CLASSC)
	printf (" %s %d %s %s %d %s\n", mp_strs[constraints[i][j].a1],
		constraints[i][j].c1, op_strs[constraints[i][j].op],
		mp_strs[constraints[i][j].a2], constraints[i][j].c2,
		(constraints[i][j].alt == TRUE) ? "|" : "");
      
      else if (i == PARAMC)
	printf (" P%d %s %s P%d %s %s\n", constraints[i][j].p1,
		mp_strs[constraints[i][j].a1], op_strs[constraints[i][j].op],
		constraints[i][j].p2, mp_strs[constraints[i][j].a2],
		(constraints[i][j].alt == TRUE) ? "|" : "");
      
      else if (i == PARAMCLASSC)
	printf (" P%d %s %d %s P%d %s %d %s\n", constraints[i][j].p1,
		mp_strs[constraints[i][j].a1], constraints[i][j].c1,
		op_strs[constraints[i][j].op], constraints[i][j].p2,
		mp_strs[constraints[i][j].a2], constraints[i][j].c2,
		(constraints[i][j].alt == TRUE) ? "|" : "");
    }
  }
}


/* Free any transient dynamic memory allocated during configuration */
void cleanupRange ()
{
  int i;

  if (penetmax)
    free (penetmax);
  if (penetcnt)
    free (penetcnt);
  if (paramcnt)
    free (paramcnt);
  if (parammax)
    free (parammax);
  for (i = 0; i < 4; i++)
    if (constraints[i])
      free (constraints[i]);
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


/**********************************************************************
 * Quicksort routines, used to sort values (both doubles and floats)
 * in the model specification.
 **********************************************************************/
#define RANDINT(n) (random() % (n))

/* Swap two doubles in array[] */
inline void swap (double *array, int i, int j)
{
  double temp = array[i];

  array[i] = array[j];
  array[j] = temp;
}

/* Recursive quicksort for array of doubles. */
void quicksort (double *array, int lo, int hi)
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
    if (array[i] < pivot) {
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
inline int uniquify (double *array, int len)
{
  int i = 0, j;

  while (i < len - 1) {
    if (array[i] == array[i + 1]) {
      for (j = i + 1; j < len - 1; j++)
	array[j] = array[j + 1];
      len--;
    } else
      i++;
  }
  return (len);
}
