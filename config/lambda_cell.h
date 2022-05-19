#ifndef __LAMBDA_CELL_H__
#define __LAMBDA_CELL_H__
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

/* LD cell type (used in modelRange below). */
typedef struct lambdaCell
{
  int m;			/* Trait/disease. */
  int n;			/* Marker. */
  int ndprime;			/* number of D prime combinations */
  double ***lambda;		/* Lambda array. */
  int *impossibleFlag;		/* whether combinatoin of D's is possible */
  double ***haploFreq;		/* haplotype frequency */
  double ***DValue;		/* D value */
}
LambdaCell;

#endif
