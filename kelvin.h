/**********************************************************************
 * Multiprocessor Linkage Analysis
 * Alberto Maria Segre, Yungui Huang
 * RADSMM storage code Martin Milder
 * 
 * Copyright 2006, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software 
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
#include "niceaux.h"		/* NICE utilities. */
#include "nicecom.h"
#include "niceapi.h"
#include "utils.h"		/* Kelvin utilities. */
#include "pedlib.h"		/* Pedigree library. */
#include "RADSMM.h"		/* RADSMM library. */

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

/**********************************************************************
 * Input files. These are not saved in the model structure because
 * they are really only ever needed by the root process.
 ***********************************************************************/
extern char mapfile[KMAXFILENAMELEN + 1];
extern char markerfile[KMAXFILENAMELEN + 1];
extern char pedfile[KMAXFILENAMELEN + 1];
#if FALSE
extern char loopsfile[KMAXFILENAMELEN + 1];
#endif
extern char outfile[KMAXFILENAMELEN + 1];

/**********************************************************************
 * Model structures (ModelType, ModelRange, ModelOptions). These
 * defines are tied to the pedlib defines wherever applicable. Note:
 * the ModelOptions structure is part of the pedigree library, while
 * the other two are defined here.
 ***********************************************************************/
#define TP TWOPOINT		/* 2 point */
#define MP MULTIPOINT		/* multipoint */
#define LE LINKAGE_EQUILIBRIUM	/* linkage equilibrium */
#define LD LINKAGE_DISEQUILIBRIUM	/* linkage disequilibrium */
#define DT DICHOTOMOUS		/* dichotomous trait */
#define QT QUANTITATIVE		/* quantitative trait */
#define CT COMBINED		/* combined dichotomous/quantitative trait */
#define ND NORMAL_DISTRIBUTION	/* normal distribution */
#define TD T_DISTRIBUTION	/* t distribution */
#define NPENET(x) (((x)*((x)+1))/2)

/* Information about the type of analysis. */
typedef struct modelType
{
  int type;			/* TP, MP: default is 2 point */
  int trait;			/* DT, QT, CT: default is dichotomous */
  int distrib;			/* Distribution type (quantitative or combined only) */
  double thresh;		/* Threshold (combined only) */
}
ModelType;

/* Information about the analysis parameters. */
typedef struct modelRange
{
  double *gfreq;		/* Gene frequencies */
  int ngfreq;			/* Number of gene frequencies */

  double ***penet;		/* 3D array of penetrance values */
  int nlclass;			/* Default is 1 liability class */
  int nalleles;			/* Default is 2 trait alleles D,d */
  int npenet;			/* Number of penetrance records */

  double ****param;		/* 4D array of extra parameters for QT and CT */
  int npardim;			/* Number of parameters for given distribution */
  int nparam;			/* Number of parameter records */

  double **theta;		/* Array of theta arrays */
  int ngender;			/* Default is sex-averaged */
  int ntheta;			/* Number of theta records */

  double **lambda;		/* Array of lambdas for disequilibrium */
  int nlamdim;			/* Lambda dimensions = (m-1)(n-1) for m trait alleles and n markers */
  int nlambda;			/* Number of lambda records */
}
ModelRange;

/**********************************************************************
 * Function prototypes.
 **********************************************************************/
int readConfigFile (char *file, ModelType * modelType,
		    ModelRange * modelRange, ModelOptions * modelOptions);
