
/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#ifndef __LAMBDA_CELL_H__
#define __LAMBDA_CELL_H__

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
