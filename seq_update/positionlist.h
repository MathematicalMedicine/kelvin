#ifndef _POSITIONLIST_H

/* positionlist.h 
 *
 * John Burian - john.burian@nationwidechildrens.org
 *
 * Copyright 2011, The Research Institute at Nationwide Children's Hospital
 * All rights reserved. Permission is granted to use this software for
 * non-profit educational purposes only.
 *
 *
 * Routines for creating and querying a sorted list of cM positions.
 */

typedef struct {
  int id;
  double min,
    max;
} st_poslimit;

typedef struct {
  char chr[8];
  int numpos,
    numlimit;
  double *pos;
  st_poslimit *limits;
} st_poschr;


void insert_pos (char *chr, double pos, int id);
int get_next_pos (char *chr, double *pos);
void free_poslist ();

#define _POSITIONLIST_H
#endif
