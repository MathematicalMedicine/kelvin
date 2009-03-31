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
