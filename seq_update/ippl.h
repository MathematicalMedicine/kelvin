#ifndef _IPPL_H

/* ippl.h 
 *
 * John Burian - john.burian@nationwidechildrens.org
 *
 * Copyright 2010, The Research Institute at Nationwide Children's Hospital
 * All rights reserved. Permission is granted to use this software for
 * non-profit educational purposes only.
 *
 *
 * Routines for handling PPLs to be used for linkage-LD updates (called 'iPPLs' since
 * the multipoint PPLs, generally calculated at regular intervals, are linerally interpolated
 * to the cM positions of the markers from the two-point LD analysis; hence interpolated
 * PPLs, or iPPLs).
 */

typedef struct {
  char chr[8];
  int numpos;
  double *pos,
    *ppl;
} st_pplchr;

void read_ppls (char *filename);
void insert_ppls (char *chr, double pos, double ppl);
double get_ippl (char *chr, double pos);
void free_ppls ();

#define _IPPL_H
#endif
