/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#ifndef __MODEL_RANGE_H__
#define __MODEL_RANGE_H__

#include "lambda_cell.h"

/* Information about the analysis parameters. */
typedef struct ModelRange
{
  double *gfreq;		/* Gene frequencies */
  int ngfreq;			/* Number of gene frequencies */

  double ***penet;		/* 3D array of penetrance values */
  int nlclass;			/* Default is 1 liability class */
  int nalleles;			/* Default is 2 trait alleles D,d */
  int npenet;			/* Number of penetrance records */
  double **penetLimits;         /* 2D array of raw penetrance limits by allele and min/max */

  double ****param;		/* 4D array of extra parameters for QT and CT */
  int npardim;			/* Number of parameters for given distribution */
  int nparam;			/* Number of parameter records */

  double **theta;		/* Array of theta arrays */
  int ngender;			/* Default is sex-averaged */
  int ntheta;			/* Number of theta records */
  int *thetacnt;		/* number of theta for male and female */

  /* Currently can only handle SNP marker frequencies. */
  double *afreq;		/* Allele frequencies */
  int nafreq;			/* Number of allele frequencies */

  double *alpha;		/* Array of alpha values */
  int nalpha;			/* Number of alpha values */

  double *tloc;			/* Array of trait locations (multipoint only) */
  int ntloc;			/* Number of trait locations */
  int tlocRangeStart;           /* Starting position i for ranges specified as 'i-end:j' */
  int tlocRangeIncr;            /* Increment j for ranges specified as 'i-end:j' */
  int tlmark;			/* Include on-marker trait loci automatically. */

  double **tthresh;		/* Array of trait thresholds (combined traits only) */
  int ntthresh;			/* Number of trait thresholds */

  /* Currently can only handle two point analysis with LD between
   * trait locus and marker, or between marker and marker. */
  double *dprime;		/* Array of dprime values for LD. */
  int ndprime;			/* Number of dprime values for LD. */

  /* For linkage disequilibrium: these are the dprime combinations
   * needed for n and m alleles at the two loci in
   * disequilibrium. Only used for two point. */
  struct lambdaCell *lambdas;	/* Cached arrays of lambdas for disequilibrium. */
  int nlambdas;			/* Number of cached arrays. */
  int maxnlambdas;		/* Maximum number of cached arrays. */
} ModelRange;

/* Until the old config parser is dead, these defines have to line up with the
 * Th, Tf, and Tm symbols in config.c, so addTheta doesn't get confused.
 */
#define THETA_AVG    1   /* same as Th */
#define THETA_MALE   2   /* same as Tm */
#define THETA_FEMALE 3   /* same as Tf */

/* Used in parsing configuration file constraints. */
#define SEXAV 0
#define SEXML 0
#define SEXFM 1

void addTraitLocus (ModelRange * range, double val);
void addGeneFreq (ModelRange * range, double val);
void addAlpha (ModelRange * range, double val);
void addDPrime (ModelRange * range, double val);
void addTheta (ModelRange * range, int type, double val);

#endif
