#ifndef _POSITIONLIST_H

/* positionlist.h 
 *
 * John Burian - john.burian@nationwidechildrens.org
 *
 * Copyright (C) 2011, 2022 Mathematical Medicine LLC
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
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
