#ifndef __summary_result__
#define __summary_result__

#include "pedlib/pedlib.h"
#include "iterationMain.h" 

#define SCALE_RESERVE 50

/* summary statistics */
typedef struct SUMMARY_STAT
{
  /* for calculating average */
  //  double lr_total;		/* sum of log(LR) */
  int scale;                    /* scale of the rescaled LR - to handle overflow */
  int lr_count;			/* number of models for the (ld, theta) pair */
  double het_lr_avg;		/* average of HET log(LR) */
  double het_lr_total;		/* with heterogeneity - alpha, parameter */
  double het_lr_avg_orig;	/* average of HET log(LR) */
  double het_lr_total_orig;	/* with heterogeneity - alpha, parameter */
  double het_lr_avg_orig2;	/* average of HET log(LR) */
  int scale_orig;               /* scale per (D', theta) pair */
  int scale_orig2;              /* scale per D' across theta */

  /* for max MOD */
  //  double max_lr;		/* max het lr - MOD */
  double max_log10_lr;          /* using max_log10_lr is to get around overflow */
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

/* three dimensional array for the two point summary results *
 * first dimension is the D prime, for LE, D prime=0 with just one element
 * in this dimension 
 * second dimension is theta values 
 * third dimension is marker allele frequency, for LE, only one element in this dimension */
extern SUMMARY_STAT ***tp_result;

/** One dimensional array, indexing by map position.  For multipoint,
 we don't know how to incorporate LD yet.  This map could be sex
 specific map or sex averaged map. For two point, we don't have to
 distinguish sex specific/avearge as we use theta relative to marker
 during analysis and after analysis (result) */
extern SUMMARY_STAT *mp_result;

/* storage for the NULL likelihood for the multipoint calculation under polynomial */
extern double *****likelihoodQT;
extern double **likelihoodDT;

int initialize_tp_result_storage ();
int free_tp_result_storage ();
void initialize_max_scale();
int record_tp_result(int callStatus, PedigreeSet *pedigreeSet, ParamStruct *param, int loc2);
int record_mp_result(int callStatus, PedigreeSet *pedigreeSet, ParamStruct *param, int posIdx);
void rescale_tp_result(int maxscale);
void rescale_tp_result_dprime0(int dprime0Idx);
int get_average_LR (SUMMARY_STAT *** result);


#endif
