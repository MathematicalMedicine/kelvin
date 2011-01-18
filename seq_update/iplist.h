#ifndef _IPLIST_H

/* iplist.h 
 *
 * John Burian - john.burian@nationwidechildrens.org
 *
 * Copyright 2011, The Research Institute at Nationwide Children's Hospital
 * All rights reserved. Permission is granted to use this software for
 * non-profit educational purposes only.
 *
 *
 * Routines for creating and managing a list of floating point values, organized
 * by chromosome and sorted by position. The values in this list can be sampled
 * at arbitrary positions, and the value returned will be linearly interpolated
 * based on the values at the actual positions that bracket the sample position.
 * Samples that fall outside the range of positions in the last will return the
 * value for the first or last actual position, as appropriate.
 */

typedef struct {
  char chr[8];
  int numpos;
  double *pos,
    *val;
} st_ipchr;

typedef struct {
  int numipchr;
  st_ipchr *ipchr;
} st_iplist;

void iplist_insert (st_iplist *list, char *chr, double pos, double val);
int iplist_lookup (st_iplist *list, char *chr, double pos, double *val);
double iplist_interpolate (st_iplist *list, char *chr, double pos);
void iplist_free (st_iplist *list);

#define _IPLIST_H
#endif
