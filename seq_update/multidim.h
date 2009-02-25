#ifndef _MULTIDIM_H

/* ippl.h 
 *
 * John Burian - john.burian@nationwidechildrens.org
 *
 * Copyright 2008, The Research Institute at Nationwide Children's Hospital
 * All rights reserved. Permission is granted to use this software for
 * non-profit educational purposes only.
 *
 *
 * Routines for managing multidimensional arrays with arbitrary index values.
 */

typedef struct {
  int arrsize,
    numelems,
    lastidx,
    dimsize;
  double *arr;
} st_dim;

typedef struct {
  short numdims;
  int totalelems;
  st_dim *dims;
} st_multidim;

int multi_insert (st_multidim *md, double *vals, int num);
int multi_find (st_multidim *md, double *vals, int num);
int insert (st_dim *dim, double val);
int find (st_dim *dim, double val, int *found);
void multi_free (st_multidim *md);
void multi_dump (st_multidim *md);

#define _MULTIDIM_H
#endif
