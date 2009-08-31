#ifndef __MODEL_OPTIONS_H__
#define __MODEL_OPTIONS_H__

/**
@file model_options.h

  model_options structure definition - contains all analysis
  options parsed from the kelvin configuration file.

  Currently only one instance of model_options is created; a global
  called modelOptions. This is the sole repository of configuration
  information, not just trait model or any other type of model. 
  Includes definition of default values for directives.

  Copyright &copy; 2009, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

*/

#include <stdio.h>

/// Maximum number of characters in an KELVIN file name.
#define KMAXFILENAMELEN 256 

/// File name used for marker data if none specified in the configuration file.
#define DEFAULTMARKERFILENAME "markers.dat"
/// File name used for map data if none specified in the configuration file.
#define DEFAULTMAPFILENAME "mapfile.dat"
/// File name used for pedigree data if none specified in the configuration file.
#define DEFAULTPEDFILENAME "pedfile.dat"
/// File name used for marker frequency data if none specified in the configuration file.
#define DEFAULTDATAFILENAME "datafile.dat"
/// File name used for output of Bayes Ratio results if none specified in the configuration file.
#define DEFAULTAVGHETFILENAME "br.out"
/// File name used for output of PPL results if none specified in the configuration file.
#define DEFAULTPPLFILENAME "ppl.out"
/// Path name used for locating saved results if none specified in the configuration file.
#define DEFAULTRESULTSPREFIX "./"

/// model_options.equilibrium value indicating linkage equilibrium analysis mode
#define LINKAGE_EQUILIBRIUM       0
/// model_options.equilibrium value indicating linkage disequilibrium analysis mode
#define LINKAGE_DISEQUILIBRIUM    1

/// model_options.markerAnalysis value indicating 2pt marker-to-marker for all possible pairs of markers
#define MARKERTOMARKER            1
#define MM MARKERTOMARKER
/// model_options.markerAnalysis value indicating 2pt marker-to-marker for adjacent markers only
#define ADJACENTMARKER            2
#define AM ADJACENTMARKER

/// model_options.mapFlag value indicating analysis uses sex-averaged map positions.
#define SEX_AVERAGED              0
#define SA SEX_AVERAGED
/// model_options.mapFlag value indicating analysis uses sex-specific map positions.
#define SEX_SPECIFIC              1
#define SS SEX_SPECIFIC

/// Phenotype code (affection status) unknown for a dichotomous trait.
#define AFFECTION_STATUS_UNKNOWN        0
/// Phenotype code (affection status) for a dichotomous trait known to be unaffected.
#define AFFECTION_STATUS_UNAFFECTED     1
/// Phenotype code (affection status) for a dichotomous trait known to be affected.
#define AFFECTION_STATUS_AFFECTED       2

typedef struct ModelOptions
{
  int equilibrium;  //< Analysis type - LINKAGE_EQUALIBRIUM or LINKAGE_DISEQUILIBRIUM
  int polynomial;  //< Flag to indicate if polynomial representation is being used.
  int polynomialScale;	//< Optional scaling factor for polynomial construction storage
  int integration;      //< Flag to indicate if dynamic grid (integration) is being performed.
  int maxIterations;    //< Optional upper limit the number of per-BR dynamic grid iterations
  int sexLinked;   //< Flag to indicate if trait is sex-linked

  /* &&& whether the chromosome is X chromosome
   * For some segments mainly at the telomeres, they can cross over with Y chromosome
   * so they behave like an autosomal segment, we call it pseudo autosomal regions
   * For X chromosome, male (XY) 's X chromosome will be doubled as 
   * homozygotes, so to handle both males and females the same way 
   */

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

