/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#ifndef __MODEL_OPTIONS_H__
#define __MODEL_OPTIONS_H__

#include <stdio.h>

/* Maximum number of characters in an KELVIN filename. */
#define KMAXFILENAMELEN 256 

#define DEFAULTMARKERFILENAME "markers.dat"
#define DEFAULTMAPFILENAME "mapfile.dat"
#define DEFAULTPEDFILENAME "pedfile.dat"
#define DEFAULTDATAFILENAME "datafile.dat"
#define DEFAULTAVGHETFILENAME "br.out"
#define DEFAULTPPLFILENAME "ppl.out"
#define DEFAULTCONDFILENAME "condL.out"
#define DEFAULTLDPPLFILENAME "ldppl.out"
#define DEFAULTRESULTSPREFIX "./"

/* Analysis mode (LE or LD) */
#define LINKAGE_EQUILIBRIUM       0
#define LINKAGE_DISEQUILIBRIUM    1

/* Marker to marker analysis types. FALSE means no analysis. */
#define MARKERTOMARKER            1   /* marker to marker analysis */
#define MM MARKERTOMARKER
#define ADJACENTMARKER            2   /* marker to adjacent marker only analysis */
#define AM ADJACENTMARKER

/* Specifiy how the map should be used */
#define SEX_AVERAGED              0   /* analysis uses sex-averaged map positions */
#define SA SEX_AVERAGED
#define SEX_SPECIFIC              1   /* analysis uses sex-specific map positions */
#define SS SEX_SPECIFIC

/* Affection status for a dichotomous trait for a person */
#define AFFECTION_STATUS_UNKNOWN        0
#define AFFECTION_STATUS_UNAFFECTED     1
#define AFFECTION_STATUS_AFFECTED       2

typedef struct ModelOptions
{
  /* analysis type - Linkage Equalibrium (LE) or Linkage Disequalibrium (LD) */
  int equilibrium;

  int polynomial;		/* TRUE if using polynomial to represent likelihood */
  int polynomialScale;		/* Optional argument to PE directive to set initial memory size */

  int integration;              /* TRUE if using dkelvin integration method */
  int maxIterations;            /* limit the number of per-BR dKelvin iterations */

  /* whether the chromosome is X chromosome
   * For some segments mainly at the telomeres, they can cross over with Y chromosome
   * so they behave like an autosomal segment, we call it pseudo autosomal regions
   * For X chromosome, male (XY) 's X chromosome will be doubled as 
   * homozygotes, so to handle both males and females the same way 
   */
  int sexLinked;

  /* whether to use sex-averaged or sex-specific map: only for
   * QT/CT. */
  int mapFlag;

  int imprintingFlag;		/* 1 - imprinting pen(1|2) may not be the same as pen(2|1) */

  /* affection status: DT defaults are {0,1,2} while QT/CT defaults
   * are {NAN, -INFINITY, INFINITY}. values are cast to double in both
   * cases. */
  double affectionStatus[3];

  /* unknown person ID */
  char *sUnknownPersonID;

  /* set recoding (super alleles) or not - we always use set recoding to speed up */
  /* int alleleSetRecodingFlag; *//* Unused */

  /* save results flag; either TRUE or FALSE */
  int saveResults;

  /* marker analysis flag; either MARKERTOMARKER or ADJACENTMARKER */
  int markerAnalysis;

  /* for calculating PPLs */
  double thetaCutoff[2];	/* using step function for the weight of each point */
  double thetaWeight;		/* weight for theta less than the thetaCutoff */
  double prior;			/* prior probability of linkage */
  double LDprior;		/* prior probability of LD (D' not 0) given linkage and theta
				 *  is within thetaCutoff */

  int dryRun;                   /* 1 - dry run to get statistics for complexity */
  int forceAvghetFile;          /* 1 - open a BR file, regardless of other directives */
  int conditionalRun;           /* 1 - print out proband's conditional LR */
  int loopCondRun;              /* 1 - print out loop breaker's conditional LR */
  int extraMODs;                /* 1 - put Theta==0 and D'==0 max models in MOD file */
  char loopBreaker[16];         /* loop breaker's individual ID */

  /* Storage and default names for files that are always opened (depending on analysis options) */
  char markerfile[KMAXFILENAMELEN + 1];  /// Marker (frequency) file
  char mapfile[KMAXFILENAMELEN + 1];     /// Map file
  char pedfile[KMAXFILENAMELEN + 1];     /// Pedigree file
  char datafile[KMAXFILENAMELEN + 1];   /// Data (pedigree description) file
  char avghetfile[KMAXFILENAMELEN + 1];       /// Bayes Ratio file
  char pplfile[KMAXFILENAMELEN + 1];         /// PPL file
  char condFile[KMAXFILENAMELEN + 1];      /// Conditional LR file
  char ldPPLfile[KMAXFILENAMELEN + 1];     /// This appears to be unused

  /* Storage for files that are only opened based on explicit directives */
  char ccfile[KMAXFILENAMELEN + 1];                 /// Case control count file
  char modfile[KMAXFILENAMELEN + 1];                /// MOD and maximizing model file
  char maxmodelfile[KMAXFILENAMELEN + 1];           /// verbose Max Model file, obsolete?
  char intermediatefile[KMAXFILENAMELEN + 1];       /// Intermediate Result file
  char dkelvinoutfile[KMAXFILENAMELEN + 1];         /// DCHURE detail file
  char resultsprefix[KMAXFILENAMELEN + 1]; ///< Path for SR directive result storage
  
} ModelOptions;

#endif

