
/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#include <stdio.h>		/* Needed for printf */
#include <stdlib.h>		/* EXIT_SUCCESS, EXIT_ERROR */
#include <unistd.h>		/* Needed for getpid */
#include <string.h>
#include <time.h>		/* Needed for sclock */
#include <signal.h>		/* Needed for signal */
#include <fcntl.h>		/* Needed for open */
#include <limits.h>		/* USHRT_MAX */
#include <math.h>		/* for calculating logs */
#include <float.h>
//#include <bits/nan.h>         /* NAN (not included by math.h?) */ 
#if FALSE
#include "niceaux.h"		/* NICE utilities. */
#include "nicecom.h"
#include "niceapi.h"
#endif
#include "utils.h"		/* Kelvin utilities. */
#include "pedlib.h"		/* Pedigree library. */
#include "trackProgress.h"
#ifdef RADSMM
#include "RADSMM.h"		/* RADSMM library. */
#endif

#include "sw.h"
#ifdef DMTRACK
#warning "Dynamic memory usage dumping is turned on, so performance will be poor!"
#define malloc(X) swMalloc((X), __FILE__, __LINE__)
#define calloc(X,Y) swCalloc((X),(Y), __FILE__, __LINE__)
#define realloc(X,Y) swRealloc((X),(Y), __FILE__, __LINE__)
#define free(X) swFree((X), __FILE__, __LINE__)
#endif

/**********************************************************************
 * Some internal defines.
 **********************************************************************/

/* Maximum number of characters in an KELVIN filename. */
#define KMAXFILENAMELEN 64

/* Maximum number of characters in a KELVIN input line. */
#define KMAXLINELEN 1024

/**********************************************************************
 * Nagging parameters.
 **********************************************************************/

/* Sleep interval when your master tells you to back off. */
#define NAGSLEEP	60	/* In seconds. */

/* Nagging messages. */
#define NQIDIDLE    	1	/* I need a problem to work on. */
#define NQIDSOLN	2	/* I have some solution values. */
#define NQIDABORT	3	/* Master wants me to abort. */

/**********************************************************************
 * Checkpointing parameters.
 **********************************************************************/
#ifdef DEBUG
#define CKPTINT    	600	/* In seconds. */
#else
#define CKPTINT		3600	/* In seconds. */
#endif

/**********************************************************************
 * Progress meters.
 **********************************************************************/
extern unsigned long int nLodsTotal;
extern unsigned long int nLodsSlave;
extern unsigned long int nLodsRoot;
extern struct swStopwatch *overallSW;

/* Performance knobs. */
extern int polynomialScale;

/**********************************************************************
 * Input files. These are not saved in the model structure because
 * they are really only ever needed by the root process.
 ***********************************************************************/
extern char mapfile[KMAXFILENAMELEN + 1];
extern char markerfile[KMAXFILENAMELEN + 1];
extern char maxmodelfile[KMAXFILENAMELEN + 1];
extern char resultsprefix[KMAXFILENAMELEN + 1];
extern char pedfile[KMAXFILENAMELEN + 1];
extern char datafile[KMAXFILENAMELEN + 1];

#if FALSE
extern char loopsfile[KMAXFILENAMELEN + 1];
#endif
extern char ccfile[KMAXFILENAMELEN + 1];
extern char avghetfile[KMAXFILENAMELEN + 1];
extern char condFile[KMAXFILENAMELEN + 1];
extern char avghomofile[KMAXFILENAMELEN + 1];
extern char pplfile[KMAXFILENAMELEN + 1];
extern char intermediatefile[KMAXFILENAMELEN + 1];
extern char modfile[KMAXFILENAMELEN + 1];
extern char dkelvinoutfile[KMAXFILENAMELEN + 1];

/**********************************************************************
 * Model structures (ModelType, ModelRange, ModelOptions). These
 * defines are tied to the pedlib defines wherever applicable. Note:
 * the ModelOptions structure is part of the pedigree library, while
 * the other two are defined here.
 ***********************************************************************/
#define TP TWOPOINT		/* 2 point */
#define MP MULTIPOINT		/* multipoint */
//#define LE LINKAGE_EQUILIBRIUM  /* linkage equilibrium */
//#define LD LINKAGE_DISEQUILIBRIUM       /* linkage disequilibrium */
#define SA SEX_AVERAGED
#define SS SEX_SPECIFIC
#define DT DICHOTOMOUS		/* dichotomous trait */
#define QT QUANTITATIVE		/* quantitative trait */
#define CT COMBINED		/* combined dichotomous/quantitative trait */
#define MM MARKERTOMARKER	/* marker to marker analysis type */
#define AM ADJACENTMARKER
#define ND NORMAL_DISTRIBUTION	/* normal distribution */
#define TD T_DISTRIBUTION	/* t distribution */
#define NPENET(x) ((x)*(x))

/* The elements of the various LD-related PPL statistics. */
typedef struct {
  double ld_small_theta,
    ld_big_theta,
    ld_unlinked,
    le_small_theta,
    le_big_theta,
    le_unlinked;
} LDVals;

/**********************************************************************
 * Function prototypes.
 **********************************************************************/
int readConfigFile (char *file);
void addTraitLocus (ModelRange * range, double val);
LambdaCell *findLambdas (ModelRange * range, int m, int n);
void addAlleleFreq (ModelRange * range, double val);

/* summary statistics */
typedef struct SUMMARY_STAT
{
  /* for calculating average */
  double lr_total;		/* sum of log(LR) */
  int lr_count;			/* number of models for the (ld, theta) pair */
  double het_lr_avg;		/* average of HET log(LR) */
  double het_lr_total;		/* with heterogeneity - alpha, parameter */

  /* for max MOD */
  double max_lr;		/* max het lr - MOD */
  double max_alpha;		/* alpha value that maximizes lr */
  double max_gfreq;		/* gfreq that maximizes lr */
  int max_penIdx;		/* penetrance that maximizes lr */
  double max_mf;		/* marker allele frequency that maximizes lr 
				 * only applies in bi-allelic case */
  double R_square;		/* only applies in bi-allelic, i.e. SNPs */
  double max_paramIdx;		/* parameter - SD index */
  double max_thresholdIdx;	/* threshold index */

  /* for max BR */
  double max_br_lr;		/* max for BR - with trait parameters integrated out */
  double max_br_mf;		/* marker allele freq that gives the max_br_lr */
  double max_br_theta;		/* theta that gives the max BR for the given D prime */
  double max_br_dprime;		/* dprime gives the overal max BR */

  /* for MP marker list */
  int *pMarkers;
  double ppl;			/* imputed MP ppl based on average likelihood ratio */
  int trait;			/* trait locus index in the locus list for MP */
} SUMMARY_STAT;


void free_likelihood_storage ();

void compute_hlod_2p_dt (double x[], double *f);
void compute_hlod_mp_dt (double x[], double *f);
void compute_hlod_2p_qt (double x[], double *f);
void compute_hlod_mp_qt (double x[], double *f);
int kelvin_dcuhre_integrate (double *integral, double *abserr, double);

/* allocate two point analysis result space */
int initialize_tp_result_storage ();
int free_tp_result_storage ();
double calculate_R_square (double p1, double q1, double d);

/* using the average at each theta to calculate PPL - posterior probability of linkage */
double calculate_PPL (SUMMARY_STAT ** result);

/* routines to caclulate the LD statistics */
int get_LDVals (SUMMARY_STAT ***result, LDVals *ldvals);
double calc_ppl_allowing_ld (LDVals *ldvals, double prior);
double calc_ppld_given_linkage (LDVals *ldvals, double prior);
double calc_ppld_allowing_l (LDVals *ldvals, double prior);
#define KROUND(dbl) dbl >= 0.025 ? rint (dbl * 100.0) / 100.0 : rint (dbl * 10000.0) / 10000.0

/* integrate out marker allele frequencies and get the max MOD */
int get_average_LR (SUMMARY_STAT *** result);

/* Model datastructures. modelOptions is defined in the pedigree library. */
extern ModelType modelType;
extern ModelRange modelRange;
extern ModelOptions modelOptions;

/* storage for the NULL likelihood for the multipoint calculation under polynomial */
extern double **likelihoodDT;
extern double *****likelihoodQT;

/* three dimensional array for the two point summary results *
 * first dimension is the D prime, for LE, D prime=0 with just one element
 * in this dimension 
 * second dimension is theta values 
 * third dimension is marker allele frequency, for LE, only one element in this dimension */
extern SUMMARY_STAT ***tp_result;
extern LambdaCell *pLambdaCell;
