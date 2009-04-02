#include <stdlib.h>
#include "model_range.h"
#include "../utils/utils.h"

#define CHUNKSIZE 64

/* Global variables for tracking actual size of the various dynamically allocated
 * lists inside ModelRange. Only interesting during configuration parsing, so no
 * point cluttering up the ModelRange structure.
 */
static int maxgfreq=0;
static int maxtloc=0;
static int maxalpha=0;
static int maxdprime=0;
static int maxtheta[2]={0,0};  /* increase this dimension if there's more than 2 genders (!) */

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
