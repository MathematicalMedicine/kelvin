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
#include "iplist.h"



void iplist_insert (st_iplist *list, char *chr, double pos, double val)
{
  int va;
  st_ipchr *ptr = NULL;

  if (list == NULL) {
    fprintf (stderr, "iplist pointer is NULL\n");
    exit (-1);
  }
  for (va = 0; va < list->numipchr; va++) {
    if (strcmp (list->ipchr[va].chr, chr) == 0) {
      ptr = &list->ipchr[va];
      break;
    }
  }
  if (ptr == NULL) {
    if ((list->ipchr = 
	 (st_ipchr *) realloc (list->ipchr, sizeof (st_ipchr) * (list->numipchr + 1))) == NULL) {
      fprintf (stderr, "realloc ipchr failed, %s\n", strerror (errno));
      exit (-1);
    }
    ptr = &list->ipchr[list->numipchr];
    list->numipchr++;
    memset (ptr, 0, sizeof (st_ipchr));
    strcpy (ptr->chr, chr);
  }

  if (((ptr->pos = (double *) realloc (ptr->pos, sizeof (double) * (ptr->numpos+1))) == NULL) || 
      ((ptr->val = (double *) realloc (ptr->val, sizeof (double) * (ptr->numpos+1))) == NULL)) {
    fprintf (stderr, "realloc failed, %s\n", strerror (errno));
    exit (-1);
  }

  if (ptr->numpos == 0) {
    va = 0;
  } else if (pos > ptr->pos[ptr->numpos-1]) {
    va = ptr->numpos;
  } else {
    for (va = ptr->numpos - 1; va >= 0; va--) {
      if (pos == ptr->pos[va]) {
	fprintf (stderr, "duplicate chr '%s', position '%f' in input\n", chr, pos);
	exit (-1);
      }
      if (pos < ptr->pos[va])
	break;
    }
    memmove (&ptr->pos[va+1], &ptr->pos[va], (ptr->numpos - va) * sizeof (double));
    memmove (&ptr->val[va+1], &ptr->val[va], (ptr->numpos - va) * sizeof (double));
  }

  ptr->pos[va] = pos;
  ptr->val[va] = val;
  ptr->numpos++;
  return;
}

int iplist_lookup (st_iplist *list, char *chr, double pos, double *val)
{
  int va;
  st_ipchr *ptr = NULL;

  if (list == NULL) {
    fprintf (stderr, "iplist pointer is NULL\n");
    exit (-1);
  }
  for (va = 0; va < list->numipchr; va++) {
    if (strcmp (list->ipchr[va].chr, chr) == 0) {
      ptr = &list->ipchr[va];
      break;
    }
  }
  if (ptr == NULL)
    return (-1);
  for (va = 0; va < ptr->numpos; va++) {
    if (ptr->pos[va] == pos) {
      *val = ptr->val[va];
      return (0);
    }
    if (ptr->pos[va] > pos)
      return (-1);
  }
  return (-1);
}

double iplist_interpolate (st_iplist *list, char *chr, double pos)
{
  int va;
  double ival;
  st_ipchr *ptr = NULL;

  if (list == NULL) {
    fprintf (stderr, "iplist pointer is NULL\n");
    exit (-1);
  }
  for (va = 0; va < list->numipchr; va++) {
    if (strcmp (list->ipchr[va].chr, chr) == 0) {
      ptr = &list->ipchr[va];
      break;
    }
  }
  if (ptr == NULL) {
    fprintf (stderr, "iplist doesn't contain data for chromosome '%s'\n", chr);
    exit (-1);
  }
  if (pos <= ptr->pos[0])
    return (ptr->val[0]);
  if (pos >= ptr->pos[ptr->numpos-1])
    return (ptr->val[ptr->numpos-1]);
  for (va = 1; va < ptr->numpos; va++) {
    if ((ptr->pos[va-1] <= pos) && (pos <= ptr->pos[va])) {
      
      ival = ptr->val[va-1] + (ptr->val[va] - ptr->val[va-1]) *
	((pos - ptr->pos[va-1]) / (ptr->pos[va] - ptr->pos[va-1]));
      return (ival);
    }
  }
  // if all else fails
  fprintf (stderr, "iplist_interpolate had no value to return; this should be impossible\n");
  exit (-1);
}


void iplist_free (st_iplist *list)
{
  int va;

  if (list == NULL)
    return;
  for (va = 0; va < list->numipchr; va++) {
    if (list->ipchr[va].numpos > 0) {
      free (list->ipchr[va].pos);
      free (list->ipchr[va].val);
    }
  }
  if (list->numipchr > 0)
    free (list->ipchr);

  return;
}
