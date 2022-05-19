#ifndef __ppl_h__
#define __ppl_h__
/* Copyright (C) 2007, 2008, 2010, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include "kelvinGlobals.h"
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
void free_likelihood_storage (PedigreeSet *);


#endif
