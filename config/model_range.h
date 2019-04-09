/**********************************************************************
 * Copyright 2010, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#ifndef __MODEL_RANGE_H__
#define __MODEL_RANGE_H__

#include "model_type.h"
#include "lambda_cell.h"

/* Information about the analysis parameters. */
typedef struct ModelRange
{
  double *gfreq;		/* Gene frequencies */
  int ngfreq;			/* Number of gene frequencies */

  double ***penet;		/* 3D array of penetrance values */
  int nlclass;			/* Default is 1 liability class */
  int *lclassLabels;            /* numeric labels for liability classes */
  int *lclassEmpty;             /* flag for empty liability classes */
  int nalleles;			/* Default is 2 trait alleles D,d */
  int npenet;			/* Number of penetrance records */
  double **penetLimits;         /* 2D array of raw penetrance limits by allele and min/max */

  double ****param;		/* 4D array of extra parameters for QT and CT */
  double paramLimits[2];        /* array of raw parameter limits by min/max */
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
  double tlocRangeStart;        /* Starting position i for ranges specified as 'i-end:j' */
  double tlocRangeIncr;         /* Increment j for ranges specified as 'i-end:j' */
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
  int atypicalQtTrait;          /* QT Trait is (or not) outside typical range */
  int microsats;                /* Microsatellites detected */
} ModelRange;


/**********************************************************************
 * Structure used for lists of constraints. A constraint can be
 * of several different forms:
 *   a1 op a2 			DD > Dd			SIMPLE
 *   a1 c1 op a2 c2		Dd 1 != dd 2		CLASSC
 *   p1 a1 op p2 a2    		P2 DD > P1 DD		PARAMC
 *   p1 a1 c1 op p2 a2 c2       P1 dd 1 >= P2 dd 2	PARAMCLASSC
 * but they are all stored in a uniform constraint structure.
 * If more than one constraint appears on a line, then we take
 * that to be an implicit OR, and all but the last constraint on such
 * a line will be marked alt=TRUE.
 ***********************************************************************/
typedef struct constraint
{
  int type;        /* one of SIMPLE, CLASSC, PARAMC, PARAMCLASSC */
  int a1;          /* one of DD, Dd, etc. */
  int c1;          /* liability class number */
  int p1;          /* QT parameter numer */
  int a2;          /* one of DD, Dd, etc. */
  int c2;          /* liability class number */
  int p2;          /* QT parameter numer */
  int op;          /* comparator */
  int alt;         /* Is this one of a disjunctive set? */
}
Constraint;

/* Legal values for Contraint.a1 and Constraint.a2. These have to
 * line up with the strings in mp_strs[] (defined in model_range.c)
 */
#define THETA_AVG     1   /* same as Th */
#define THETA_MALE    2   /* same as Tm */
#define THETA_FEMALE  3   /* same as Tf */
#define THRESHOLD     4   /* same as TT */
#define PEN_DD        5   /* same as DD */
#define PEN_Dd        6   /* same as Dd */
#define PEN_dD        7   /* same as dD */
#define PEN_dd        8   /* same as dd */

/* Legal values for Constraint.type. These have to line up with the
 * strings in contype_strs[] (defined in model_range.c)
*/
#define SIMPLE 0
#define CLASSC 1
#define PARAMC 2
#define PARAMCLASSC 3

/* Legal values for Constraint.op; these have to line up with the
 * strings in op_strs[] (defined in model_range.c)
 */
#define EQ 0
#define NE 1
#define GT 2
#define GE 3

#define NPENET(x) ((x)*(x))

/* Indices into the thetacnt and theta members of ModelRange. NOTE
 * explicitly NOT interchangeable with THETA_AVG, etc.
 */
#define SEXAV 0
#define SEXML 0
#define SEXFM 1

/* Default min/max values for QT/QTT parameters, under dynamic sampling only */
/* For Normal (T) distribution */
#define QT_INTEG_MIN_MEAN -3.0
#define QT_INTEG_MAX_MEAN  3.0
/* "Classic" StDev range was 0.7 - 1.0 */
#define QT_INTEG_MIN_STDEV 0.5   // changed from 0 on April 9, 2019 to avoid possible numerical problems
#define QT_INTEG_MAX_STDEV 2.5   // changed from 4 on April 9, 2019 
#define QT_INTEG_NORM_MIN_THRESH  0.0
#define QT_INTEG_NORM_MAX_THRESH  3.0
/* fro Chi-squared distribution */
#define QT_INTEG_MIN_DOF   0.05
#define QT_INTEG_MAX_DOF  30.0
#define QT_INTEG_CHI_MIN_THRESH   0.05
#define QT_INTEG_CHI_MAX_THRESH  30.0


extern char *op_strs[];
extern char *mp_strs[];
extern char *contype_strs[];

void addTraitLocus (ModelRange * range, double val);
void addGeneFreq (ModelRange * range, double val);
void addAlleleFreq (ModelRange * range, double val);
void addAlpha (ModelRange * range, double val);
void addDPrime (ModelRange * range, double val);
void addTheta (ModelRange * range, int type, double val);
void addPenetrance (ModelRange * range, int type, double val);
void addConstraint (int type, int a1, int c1, int p1, int op,
		    int a2, int c2, int p2, int disjunct);
void addParameter (ModelRange * range, int dim, double val);
void addTraitThreshold (ModelRange * range, double val);
int checkImprintingPenets (ModelRange *range, int imprinting);
void duplicatePenets (ModelRange *range, int imprinting);
int checkThetas (ModelRange * range, int i);
int checkPenets (ModelRange * range, int i);
int checkClassPenets (ModelRange * range, int i);
int checkParams (ModelRange * range, int i);
int checkClassParams (ModelRange * range, int i);
int checkClassThreshold (ModelRange * range, int i);
void expandRange (ModelRange *range);
void expandClassThreshold (ModelRange * range);
void expandClassPenet (ModelRange * range);
LambdaCell *findLambdas (ModelRange * range, int m, int n);
void sortRange (ModelRange * range);
void uniqRange (ModelRange * range);
void showRange (ModelRange * range, ModelType * type, int level);
void showConstraints ();
void cleanupRange ();

int lookup_comparator (char *str);
int lookup_modelparam (char *str);

#endif
