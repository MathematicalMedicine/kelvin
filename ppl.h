#ifndef __ppl_h__
#define __ppl_h__

#include "summary_result.h"

/* The elements of the various LD-related PPL statistics. */
typedef struct {
  double ld_small_theta,
    ld_big_theta,
    ld_unlinked,
    le_small_theta,
    le_big_theta,
    le_unlinked;
} LDVals;

/* routines to caclulate the LD statistics */
int get_LDVals (SUMMARY_STAT ***result, LDVals *ldvals);
double calc_ppl_allowing_ld (LDVals *ldvals, double prior);
double calc_ppld_given_linkage (LDVals *ldvals, double prior);
double calc_ppld_allowing_l (LDVals *ldvals, double prior);
/* using the average at each theta to calculate PPL - posterior probability of linkage */
double calculate_PPL (SUMMARY_STAT ** result);

/* it shouldn't be residing here, will find appropriate place for this function later */
double calculate_R_square (double p1, double q1, double d);
void free_likelihood_storage ();


#endif
