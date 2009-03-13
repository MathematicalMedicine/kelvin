/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#ifndef __MODEL_TYPE_H__
#define __MODEL_TYPE_H__

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
  double minThreshold;		/* minimum threshold value - in standardized unit */
  double maxThreshold;		/* maximum threshold value - in standardized unit */
  int *constants;		/* Array of distribution constants (certain QT/CT distributions only) */
  int ccFlag;			/* Case Ctrl flag */
} ModelType;

#endif
