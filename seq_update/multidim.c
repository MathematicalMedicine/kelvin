/* Copyright (C) 2009, 2010, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <multidim.h>

int multi_insert (st_multidim *md, double *vals, int num)
{
  int va, idx, offset=0;
  st_dim *tmp;

  if (md->numdims != num) {
    if ((tmp = realloc (md->dims, sizeof (st_dim) * num)) == NULL)
      return (-1);
    md->dims = tmp;
    for (va = md->numdims; va < num; va++)
      memset (&(md->dims[va]), 0, sizeof (st_dim));
    md->numdims = num;
  }
  
  for (va = md->numdims - 1; va >= 0; va--) {
    if ((idx = insert (&(md->dims[va]), vals[va])) == -1)
      return (-1);
    if (va == md->numdims - 1)
      md->dims[va].dimsize = 1;
    else
      md->dims[va].dimsize = md->dims[va+1].numelems * md->dims[va+1].dimsize;
    offset += md->dims[va].dimsize * idx;
  }
  md->totalelems = md->dims[0].dimsize * md->dims[0].numelems;
  return (offset);
}


int multi_find (st_multidim *md, double *vals, int num)
{
  int va, idx, found, offset=0;

  if (md->numdims != num)
    return (-1);
  for (va = md->numdims - 1; va >= 0; va--) {
    idx = find (&(md->dims[va]), vals[va], &found);
    if (! found)
      return (-1);
    offset += idx * md->dims[va].dimsize;
  }
  return (offset);
}


int insert (st_dim *dim, double val)
{
  double *tmp;
  int idx, found;

  idx = find (dim, val, &found);
  if (found)
    return (idx);

  if (dim->arrsize < dim->numelems + 1) {
    if ((tmp = realloc (dim->arr, sizeof (double) * (dim->arrsize + 10))) == NULL)
      return (-1);
    dim->arr = tmp;
    dim->arrsize += 10;
  }
  if (idx < dim->numelems)
    memmove (dim->arr+idx+1, dim->arr+idx, sizeof (double) * (dim->numelems - idx));
  dim->arr[idx] = val;
  dim->numelems++;
  return (idx);
}


int find (st_dim *dim, double val, int *found)
{
  int va;

  *found = 0;
  if (dim->numelems == 0)
    return (dim->lastidx = 0);
  if (val == dim->arr[dim->lastidx]) {
    *found = 1;
    return (dim->lastidx);
  }

  if (val > dim->arr[dim->lastidx]) {
    for (va = dim->lastidx + 1; va < dim->numelems; va++) {
      if (val < dim->arr[va]) {
	return (dim->lastidx = va);
      } else if (val == dim->arr[va]) {
	*found = 1;
	return (dim->lastidx = va);
      }
    }
    return (dim->lastidx = dim->numelems);
    
  } else {   /*  val < dim->arr[dim->lastidx]  */
    for (va = 0; va <= dim->lastidx; va++) {
      if (val < dim->arr[va]) {
	return (dim->lastidx = va);
      } else if (val == dim->arr[va]) {
	*found = 1;
	return (dim->lastidx = va);
      }
    }
  }
  /* Should be impossible to reach this point */
  return (-1);
}


void multi_free (st_multidim *md)
{
  int va;

  for (va = 0; va < md->numdims; va++)
    free (md->dims[va].arr);
  free (md->dims);
  return;
}


void multi_dump (st_multidim *md)
{
  int va, vb;

  for (va = 0; va < md->numdims; va++) {
    printf ("dim %d:", va);
    for (vb = 0; vb < md->dims[va].arrsize; vb++) {
      if (vb < md->dims[va].numelems) {
	printf (" %5.2f", md->dims[va].arr[vb]);
      } else {
	printf (" .....");
      }
    }
    printf (" (lastidx %d, dimsize %d)\n", md->dims[va].lastidx, md->dims[va].dimsize);
  }
  printf ("totalelems %d\n", md->totalelems);
  return;
}
