#ifndef __MODEL_TYPE_H__
#define __MODEL_TYPE_H__
/* Copyright (C) 2014, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/* Analysis type in a different perspective: 2point or multipoint */
#define TWOPOINT        0   /* 2 point */
#define TP TWOPOINT
#define MULTIPOINT      1   /* multipoint */
#define MP MULTIPOINT

/* Trait value type: dichotomous, quantitative, or combined
 *
 * The main difference between DT and QT are in the way how conditional
 * penetrances are determined. For DT, penetrances are given, while for
 * QT, they are calculated using means and variances.
 *
 * Combination of dichotomous and quantitative trait 
 * some individuals in the pedigree has known quantitative trait values
 * individuals with unknown trait values are handled in a way similar to 
 * under dichotomouse trait (using some cutoff values to calculate conditional
 * penetrance given genotypes with CDF)
 */
#define DICHOTOMOUS     0   /* dichotomous trait */
#define DT DICHOTOMOUS
#define QUANTITATIVE    1   /* quantitative trait */
#define QT QUANTITATIVE
#define COMBINED        2   /* combined dichotomous/quantitative trait */
#define CT COMBINED

/* QT underlying distribution function */
#define QT_FUNCTION_NORMAL		0
#define QT_FUNCTION_T			1
#define QT_FUNCTION_CHI_SQUARE          2

/* Information about the type of analysis. The array of constants is
 * only used for some types of QT/CT models (for example, the degrees
 * of freedom of the T distribution are stored in constants[0]). */

typedef struct ModelType
{
  int type;			/* TP, MP: default is 2 point */
  int numMarkers;		/* number of markers to use for MP */
  int trait;			/* DT, QT, CT: default is dichotomous */
  int distrib;			/* Distribution type (QT/CT only) */
  double mean;			/* mean (QT/CT only) */
  double sd;			/* standard deviation (QT/CT only) */
  double minOriginal;		/* QT - minimum value of QT - left truncate */
  double min;			/* QT - min after standardization */
  double maxOriginal;		/* QT - maximum value of QT - right truncate */
  double max;			/* QT - max after standardization */
  int minFlag;			/* flag for the existence of lower bound */
  int maxFlag;			/* flag for the existence of upper bound */
  /* remove threshold adjustment code 
     YH 04/14/2009 */
#if 0
  double minThreshold;		/* minimum threshold value - in standardized unit */
  double maxThreshold;		/* maximum threshold value - in standardized unit */
#endif 
  int *constants;		/* Array of distribution constants (certain QT/CT distributions only) */
  int ccFlag;			/* Case Ctrl flag */
} ModelType;


double deNormalizeMean (ModelType *mt, double mean);
double deNormalizeStdev (ModelType *mt, double stdev);

#endif
