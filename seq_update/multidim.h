#ifndef _MULTIDIM_H

/* multidim.h 
 *
 * John Burian - john.burian@nationwidechildrens.org
 *
 * Copyright (C) 2009, 2010, 2022 Mathematical Medicine LLC
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
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
