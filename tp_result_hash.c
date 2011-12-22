/**
@file tp_result_hash.c

 Maintain a hash of tp_result structures primarily for dynamic grid, but used by
 fixed grid as well to enable use of uniform output routines.

 Written by Bill Valentine-Cooper.

 Copyright 2010, Nationwide Children's Research Institute.  All rights
 reserved.  Permission is hereby given to use this software for
 non-profit educational purposes only.

 See the include file tp_result_hash.h for documentation.

 Includes a test driver. Build with:

  gcc -DMAIN -g -o tp_result_hash tp_result_hash.c -I utils/ -lutils

*/
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include "utils/hashtab.h"
#include "utils/utils.h"

#include "tp_result_hash.h"

// Hash table to find them individually
htab *tp_result_hash = NULL;
#define TPRHASHSIZE 15
#define TPRKEYSIZE 12

// Structure to hold key and next pointers as well as actual data
struct tp_result_hash_entry {
  SUMMARY_STAT *tp_result;
  char keyString[32];
  void *next;
};

// List of tp_result structures so we can qsort them and provide them sequentially
SUMMARY_STAT **tp_result_list = NULL;
int tp_result_list_count = 0; // How many there actually are
int tp_result_list_limit = 0; // Amount of space allocated

// Dump the tp_result_hash.
void dump_tp_result_hash ()
{
  struct tp_result_hash_entry *newTPRhe;
  int i;

  if (tp_result_hash == NULL)
    tp_result_hash = hcreate (TPRHASHSIZE);

  if (hfirst (tp_result_hash))
    do {
      ub1 *tp_result_hash_key;

      newTPRhe = (struct tp_result_hash_entry *) hstuff (tp_result_hash);
      tp_result_hash_key = hkey (tp_result_hash);
      printf ("key ");
      for (i = 0; i < 4; i++) {
	printf ("%02x ", tp_result_hash_key[i]);
      }
      printf ("for indices %d, %d and %d\n",
	      newTPRhe->tp_result->dprimeIdx, newTPRhe->tp_result->thetaIdx, newTPRhe->tp_result->mkrFreqIdx);
    }
    while (hnext (tp_result_hash));
}

// Return a pointer to the tp_result structure for the given indices, or NULL if not found.
SUMMARY_STAT *get_tp_result (int dprimeIdx, int thetaIdx, int mkrFreqIdx)
{
  struct tp_result_hash_entry targetTPRhe, *oldTPRhe;

  if (tp_result_hash == NULL)
    tp_result_hash = hcreate (TPRHASHSIZE);

  sprintf (targetTPRhe.keyString, "%d,%d,%d", dprimeIdx, thetaIdx, mkrFreqIdx);
  if (hfind (tp_result_hash, (ub1 *)targetTPRhe.keyString, strlen(targetTPRhe.keyString)) == FALSE) {
    fprintf (stderr, "Warning - non-existant tp_result referenced with indices (%d, %d and %d), returning NULL\n",
	     dprimeIdx, thetaIdx, mkrFreqIdx);
    return NULL;
  } else {
    oldTPRhe = (struct tp_result_hash_entry *) hstuff (tp_result_hash);

    while (strcmp (oldTPRhe->keyString, targetTPRhe.keyString)) {
      if (oldTPRhe->next != NULL) {
	oldTPRhe = oldTPRhe->next;
      } else {
	fprintf (stderr, "Warning - non-existant tp_result referenced with indices (%d, %d and %d), returning NULL\n",
	    dprimeIdx, thetaIdx, mkrFreqIdx);
	return NULL;
      }
    }
    return oldTPRhe->tp_result;
  }
}

// Sort the list of tp_result structures by the three indices. Or whatever.
int
compareTPRsByIndices (const void *left, const void *right)
{
  struct SUMMARY_STAT *leftTPR, *rightTPR;
  leftTPR = (struct SUMMARY_STAT *) (* (struct SUMMARY_STAT **) left);
  rightTPR = (struct SUMMARY_STAT *) (* (struct SUMMARY_STAT **) right);

  if (leftTPR->dprimeIdx == rightTPR->dprimeIdx)
    if (leftTPR->thetaIdx == rightTPR->thetaIdx)
      if (leftTPR->mkrFreqIdx == rightTPR->mkrFreqIdx)
	return 0;
      else
	return leftTPR->mkrFreqIdx < rightTPR->mkrFreqIdx ? -1 : 1  ;
    else
      return leftTPR->thetaIdx < rightTPR->thetaIdx ? -1 : 1;
  else
    return leftTPR->dprimeIdx < rightTPR->dprimeIdx ? -1 : 1;
}	

void
sort_tp_result_by_indices ()
{
  qsort (tp_result_list, tp_result_list_count,
	 sizeof (struct SUMMARY_STAT *), compareTPRsByIndices);
}

// Return a pointer to the next (offset) tp_result structure, or NULL if there are no more.
SUMMARY_STAT *get_next_tp_result (int *offset)
{
  if (*offset >= tp_result_list_count)
    return NULL;
  else
    return tp_result_list[(*offset)++];
}

// Create a new tp_result structure for the given indices and return a pointer to it.
SUMMARY_STAT *new_tp_result (int dprimeIdx, int thetaIdx, int mkrFreqIdx)
{
  struct tp_result_hash_entry *newTPRhe, *oldTPRhe;

  if (tp_result_hash == NULL)
    tp_result_hash = hcreate (TPRHASHSIZE);

  // Make sure we have enough room for an addition
  if (tp_result_list_count >= tp_result_list_limit) {
    tp_result_list_limit += 1024;
    if (tp_result_list != NULL) {
      REALCHOKE(tp_result_list, tp_result_list_limit * sizeof(struct SUMMARY_STAT *), struct SUMMARY_STAT **);
    } else {
      MALCHOKE(tp_result_list, tp_result_list_limit * sizeof (struct SUMMARY_STAT *), struct SUMMARY_STAT **);
    }
  }

  MALCHOKE(newTPRhe, sizeof (struct tp_result_hash_entry), struct tp_result_hash_entry *);
  CALCHOKE(newTPRhe->tp_result, sizeof (SUMMARY_STAT), 1, SUMMARY_STAT *);
  newTPRhe->next = NULL;

  // Combine the 3 integer indexes of tp_result into a single sequence of bytes and generate a key.
  newTPRhe->tp_result->dprimeIdx = dprimeIdx; newTPRhe->tp_result->thetaIdx = thetaIdx;
  newTPRhe->tp_result->mkrFreqIdx = mkrFreqIdx;

  /* See if we can add the new tp_result entry */
  sprintf (newTPRhe->keyString, "%d,%d,%d", dprimeIdx, thetaIdx, mkrFreqIdx);
  if (hadd (tp_result_hash, (ub1 *) newTPRhe->keyString, strlen(newTPRhe->keyString), newTPRhe) == FALSE) {
    // Thump! Collision, go for next!
    oldTPRhe = (struct tp_result_hash_entry *) hstuff (tp_result_hash);
    while (strcmp (oldTPRhe->keyString, newTPRhe->keyString)) {
      if (oldTPRhe->next == NULL) {
	oldTPRhe->next = newTPRhe;
	tp_result_list[tp_result_list_count++] = newTPRhe->tp_result;
	return newTPRhe->tp_result;
      }
    }
    // Reuse! Whine about it, free what we've allocated and return the old one
    fprintf (stderr, "Warning - new_tp_result called previously with same indices (%d, %d and %d), returning original\n",
	    dprimeIdx, thetaIdx, mkrFreqIdx);
    free (newTPRhe->tp_result);
    free (newTPRhe);
    return oldTPRhe->tp_result;
  }
  // Add worked!
  tp_result_list[tp_result_list_count++] = newTPRhe->tp_result;
  return newTPRhe->tp_result;
}

// The whole point to this routine is to make life simpler for "fkelvin" static allocation replacement.
// Since we can't assume this was called before new_tp_result, we have to use new_tp_result to do the work.
void
new_bulk_tp_result (int dprimeDim, int thetaDim, int mkrFreqDim)
{
  int i, j, k;
  for (i=0; i<dprimeDim; i++)
    for (j=0; j<thetaDim; j++)
      for (k=0; k<mkrFreqDim; k++)
	new_tp_result (i, j, k);
}

#ifdef MAIN

int main (int argc, char *argv[]) {
  SUMMARY_STAT *my_tp_result;

  my_tp_result = new_tp_result (1, 6, 3);
  my_tp_result->ppl = 16.3;

  my_tp_result = new_tp_result (3, 2, 3);
  my_tp_result->ppl = 32.3;

  my_tp_result = new_tp_result (32768, 2, 3);
  my_tp_result->ppl = 327682.3;

  my_tp_result = new_tp_result (32777, 2, 3);
  my_tp_result->ppl = 327772.3;
  
  my_tp_result = new_tp_result (3, 2, 4);
  my_tp_result->ppl = 32.4;
  
  my_tp_result = new_tp_result (3, 2, 12345678);
  my_tp_result->ppl = 32.12345678;
  
  my_tp_result = new_tp_result (2, 2, 6);
  my_tp_result->ppl = 22.6;

  my_tp_result = new_tp_result (1, 6, 3);
  my_tp_result->ppl = 16.4;

  my_tp_result = new_tp_result (1923, 677, 3812);
  my_tp_result->ppl = 1923677.3812;

  my_tp_result = new_tp_result (923, 677, 3812);
  my_tp_result->ppl = 923677.3812;

  my_tp_result = new_tp_result (1, 6, 3);
  my_tp_result->ppl = 16.5;

  dump_tp_result_hash ();

  if ((my_tp_result = get_tp_result (3, 2, 3)) == NULL) {
    printf ("Couldn't find  3, 2, 3\n");
  } else
    printf ("3, 2, 3 has %g\n", my_tp_result->ppl);

  if ((my_tp_result = get_tp_result (2, 2, 6)) == NULL) {
    printf ("Couldn't find  2, 2, 6\n");
  } else
    printf ("2, 2, 6 has %g\n", my_tp_result->ppl);

  if ((my_tp_result = get_tp_result (2, 2, 77)) == NULL) {
    printf ("Couldn't find  2, 2, 77\n");
  } else
    printf ("2, 2, 77 has %g\n", my_tp_result->ppl);

  /* The following approach allows easy conversion from an indexed structure, but is less
     efficient if the pointer returned can be reused. Probably doesn't matter with tp_result
     since it's only updated when we have results that took a long time to produce! */

  printf ("Before assignment, ppl is %g\n", (get_tp_result(2, 2, 6))->ppl);
  (get_tp_result(2, 2, 6))->ppl = 18.5;
  printf ("Trying again, ppl is %g\n", (get_tp_result(2, 2, 6))->ppl);

  sort_tp_result_by_indices();

  int offset = 0;
  while ((my_tp_result = get_next_tp_result (&offset)) != NULL) {
    printf ("%d, %d, %d -> %g\n", my_tp_result->dprimeIdx, my_tp_result->thetaIdx, my_tp_result->mkrFreqIdx, my_tp_result->ppl);
  }

  new_bulk_tp_result (10, 20, 30);

  (get_tp_result(1, 19, 29))->ppl = 273.6;
  printf ("Just set bulk item to %g\n", (get_tp_result(1, 19, 29))->ppl);

  (get_tp_result(1, 19, 33))->ppl = 1972.83;

}
#endif
