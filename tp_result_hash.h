/**
@file tp_result_hash.c

 Maintain a hash of tp_result structures primarily for dynamic grid, but used by
 fixed grid as well to enable use of uniform output routines.

 Written by Bill Valentine-Cooper.

 Copyright (C) 2009, 2010, 2022 Mathematical Medicine LLC
 This program is free software: you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the Free
 Software Foundation, either version 3 of the License, or (at your option)
 any later version.
 You should have received a copy of the GNU General Public License along
 with this program. If not, see <https://www.gnu.org/licenses/>.

 Provides:

   For instantiation:

     SUMMARY_STAT *new_tp_result(int dprimeIdx, int thetaIdx, int mkrFreqIdx) - to allocate
     and return an individual tp_result structure in the hash. Prints a warning to stderr and
     returns a pointer to the original structure when the same set of indices is used more 
     than once.

   or

     bulk_new_tp_result(int dprimeDim, int thetaDim, int mkrFreqDim) - to allocate a fixed
     grid of tp_result structures in the hash.

   For access:

     SUMMARY_STAT *get_tp_result(int dprimeIdx, int thetaIdx, int mkrFreqIdx) - to retrieve
     a pointer to an existing tp_result structure in the hash. If you're confident
     that it exists, you can embed the function call in an assignment or reference
     as in:
             (get_tp_result(2, 5, 5))->ppl = 0.0234;
        or:
             printf ("PPL for %d, %d, %d is %g\n", i, j, k, get_tp_result(i, j, k));

     Prints a warning to stderr and returns a NULL when a non-existant tp_result is requested.

   For iteration:

     SUMMARY_STAT *get_next_tp_result(int *n) - to retrieve the nth tp_result structure in
     order of creation and increment n, or call:

     sort_tp_result_by_indices() - to change the order of entries returned to ascending
     numerically by indices, dprimeIdx, thetaIdx and finally mkrFreqIdx.

     returns a NULL when there are no more entries.

 See examples of use in the sample MAIN at the bottom of the source file.     

 Uses the hashtab routines in util, just like sw.c. Currently only
 provides next entries in the order in which they were created.

 Currently not thread-safe.

*/

#include "summary_result.h"

SUMMARY_STAT *new_tp_result(int dprimeIdx, int thetaIdx, int mkrFreqIdx);
void bulk_new_tp_result(int dprimeDim, int thetaDim, int mkrFreqDim);
SUMMARY_STAT *get_tp_result(int dprimeIdx, int thetaIdx, int mkrFreqIdx);
SUMMARY_STAT *get_next_tp_result(int *n);
void sort_tp_result_by_indices();
