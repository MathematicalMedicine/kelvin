/* Copyright (C) 2011, 2022 Mathematical Medicine LLC
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
#include <errno.h>
#include "positionlist.h"

st_poschr *poschrs = NULL;
int numposchr = 0;

int curchr = 0, curpos = 0;

void insert_pos (char *chr, double pos, int id)
{
  int va;
  st_poschr *ptr = NULL;

  for (va = 0; va < numposchr; va++) {
    if (strcmp (poschrs[va].chr, chr) == 0) {
      ptr = &poschrs[va];
      break;
    }
  }
  if (ptr == NULL) {
    if ((poschrs = (st_poschr *) realloc (poschrs, sizeof (st_poschr) * (numposchr + 1))) == NULL) {
      fprintf (stderr, "realloc poschrs failed, %s\n", strerror (errno));
      exit (-1);
    }
    ptr = &poschrs[numposchr];
    numposchr++;
    memset (ptr, 0, sizeof (st_poschr));
    strcpy (ptr->chr, chr);
  }

  for (va = 0; va < ptr->numpos; va++) {
    if (pos == ptr->pos[va]) {
      va = -1;
      break;
    }
    if (pos < ptr->pos[va])
      break;
  }
  if (va != -1) {
    if ((ptr->pos = (double *) realloc (ptr->pos, sizeof (double) * (ptr->numpos + 1))) == NULL) {
      fprintf (stderr, "realloc positions failed, %s\n", strerror (errno));
      exit (-1);
    }
    if (va < ptr->numpos)
      memmove (&ptr->pos[va+1], &ptr->pos[va], (ptr->numpos - va) * sizeof (double));
    ptr->pos[va] = pos;
    ptr->numpos++;
  }

  for (va = 0; va < ptr->numlimit; va++) {
    if (id == ptr->limits[va].id) {
      if (ptr->limits[va].min > pos)
	ptr->limits[va].min = pos;
      if (ptr->limits[va].max < pos)
	ptr->limits[va].max = pos;
      return;
    }
  }
  if ((ptr->limits = (st_poslimit *) realloc (ptr->limits, sizeof (st_poslimit) * (ptr->numlimit + 1))) == NULL) {
    fprintf (stderr, "realloc limits failed, %s\n", strerror (errno));
    exit (-1);
  }
  ptr->limits[ptr->numlimit].id = id;
  ptr->limits[ptr->numlimit].min = ptr->limits[ptr->numlimit].max = pos;
  ptr->numlimit++;
  return;
}


int get_next_pos (char *chr, double *pos)
{
  if (numposchr == 0)
    return (-1);

  if (strlen (chr) == 0)
    curchr = curpos = 0;

  if (curchr < numposchr && curpos < poschrs[curchr].numpos) {
    strcpy (chr, poschrs[curchr].chr);
    *pos = poschrs[curchr].pos[curpos++];
    return (0);
  } else if (curchr + 1 < numposchr) {
    curchr++;
    curpos = 0;
    strcpy (chr, poschrs[curchr].chr);
    *pos = poschrs[curchr].pos[curpos++];
    return (0);
  } else {
    return (-1);
  }
}


void free_poslist ()
{
  int va;

  for (va = 0; va < numposchr; va++) {
    free (poschrs[va].pos);
    free (poschrs[va].limits);
  }
  free (poschrs);
}
