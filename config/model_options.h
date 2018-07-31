#ifndef __MODEL_OPTIONS_H__
#define __MODEL_OPTIONS_H__

/**
@file model_options.h

  model_options structure definition - contains all analysis
  options parsed from the kelvin configuration file, and a few
  that maybe should be in the configuration file but aren't yet.

  Currently only one instance of model_options is created; a
  sanctioned global called modelOptions. This is the sole repository
  of configuration information, not just trait model or any other type
  of model.  Includes definition of default values for directives.

  Copyright &copy; 2010, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

*/

#include <stdio.h>
#include <limits.h>

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
/// File name used for output of MOD results if none specified in the configuration file.
#define DEFAULTMODFILENAME "mod.out"
/// Path name used for locating saved results if none specified in the configuration file.
#define DEFAULTRESULTSPREFIX "./"

/// model_options.equilibrium value indicating linkage equilibrium analysis mode.
#define LINKAGE_EQUILIBRIUM       0
/// model_options.equilibrium value indicating linkage disequilibrium analysis mode.
#define LINKAGE_DISEQUILIBRIUM    1

enum MarkerAnalysis {
  TRAITTOMARKER = 0, ///< 2pt analysis is trait-to-marker. Should stay zero so tests for FALSE will work.
  MARKERTOMARKER = 1, ///< 2pt analysis is marker-to-marker for all possible pairings of markers.
  MM = 1, ///< Synonym for MARKERTOMARKER.
  ADJACENTMARKER = 2, ///< 2pt analysis is marker-to-marker for adjacent markers  only.
  AM = 2 ///< Synonym for ADJACENTMARKER.
};

enum MapFlag {
  SEX_AVERAGED = 0, ///< model_options.mapFlag value indicating analysis uses sex-averaged map positions.
  SA = 0, //< Synonym for SEX_AVERAGED.
  SEX_SPECIFIC = 1, ///< model_options.mapFlag value indicating analysis uses sex-specific map positions.
  SS = 1 //< Synonym for SEX_SPECIFIC.
};

/// Phenotype code (affection status) unknown for a dichotomous trait.
#define AFFECTION_STATUS_UNKNOWN        0
/// Phenotype code (affection status) for a dichotomous trait known to be unaffected.
#define AFFECTION_STATUS_UNAFFECTED     1
/// Phenotype code (affection status) for a dichotomous trait known to be affected.
#define AFFECTION_STATUS_AFFECTED       2

/* For experimental dynamic-grid QT-normal analysis, used to control if means
 * or standard deviations can vary (parameter is sampled, and may vary between
 * trait genotypes), are the same (parameter is sampled, but does not vary 
 * betgween trait genotypes) or is fixed (is not sampled, but is set to a 
 * fixed value for all models).
 */
#define QT_MODE_VARY    1    // Sampled and varies
#define QT_MODE_SAME    2    // Sampled but does not vary
#define QT_MODE_FIXED   3    // Fixed

typedef struct ModelOptions
{
  int equilibrium;  ///< Analysis type - LINKAGE_EQUILIBRIUM or LINKAGE_DISEQUILIBRIUM
  int polynomial;  ///< Flag to indicate if polynomial representation is being used.
  int polynomialScale;	///< Optional scaling factor for polynomial construction storage
  int integration;      ///< Flag to indicate if dynamic grid (integration) is being performed.
  int maxIterations;    ///< Optional upper limit the number of per-BR dynamic grid iterations
  int sexLinked;   /**< Flag to indicate that analysis focuses on non-pseudoautosomal region of
		      the X chromosome. Should not be specified for X chromosome analyses that
		   focus on the pseudoautosomal regions because they behave like autosomes. */
  enum MapFlag mapFlag; ///< Flag to indicate whether to use the sex-averaged or sex-specific map
  int imprintingFlag;	/**< Flag to indicate if imprinting is to be considered. Requires additional
			   penetrance values (a set for each phase, e.g. dD and Dd) in the 
			   configuration file. */
  double affectionStatus[3]; /**< The three special values in the affection status column of the
				pedigree file that will indicate unknown, unaffected and affected
				individuals. There are different defaults for dichotomous and 
				quantitative trait models. For dichotomous traits, there can be
				some confusion because the symbols used for the	values 
				(AFFECTION_STATUS_UNKNOWN, etc) do double-duty as the indices
				of this vector, i.e. affectionStatus[AFFECTION_STATUS_UNKNOWN] =
				AFFECTION_STATUS_UNKNOWN. For quantitative traits, the default
				values for unknown, unaffected and affected are -99.99, -88.88
				and 88.88. It is unfortunate that these are "magic" values, but
				changing this is non-trivial. */
  char *sUnknownPersonID; ///< Unknown person ID, defaults to zero.
  int saveResults; /**< Flag to indicate whether or not we should save the results of likelihood
		      calculations in a file for reuse in subsequent runs, in lieu of 
		      recalculation. */
  enum MarkerAnalysis markerAnalysis; ///< Flag to indicate which loci are involved in a 2pt analysis

  /* For calculating PPLs. Not in the configuration file yet. */
  double thetaCutoff[2];	///< Theta weighting cutoffs for sex-average or male [0] and sex-specific female [1].
  double thetaWeight;		///< Weight for theta when it's less than the appropriate thetaCutoff.
  double prior;			///< Prior probability of linkage.
  double LDprior;		///< Prior probability of LD (D' not 0) given linkage and theta within thetaCutoff.

  /* For experimental LOD maximization algorithm */
  double modThreshold;          ///< Threshold over which LOD maximization kicks in

  /* Flag for experimental alternative QT mode, which may become our default in the future (but not now). */
  int alternativeQTFlag;        ///< Flag to indicate if experimental alternative QT algorithm is to be used.

  /* For experimental QT analysis integrating over multiple SD and/or means. */
  int qtMeanMode;
  int qtStandardDevMode;
  int qtThresholdMode;
  
  int dryRun;                   ///< Flag indicating dry run to get statistics for complexity.
  int forceAvghetFile;          ///< Flag to force open a BR file, regardless of other directives.
  int conditionalRun;           ///< Flag indicating to print out proband's conditional LR.
  int loopCondRun;              ///< Flag indicating to print out loop breaker's conditional LR.
  char loopBreaker[16];         ///< When loopCondRun set, use this ID to identify the loop breaker */
  int extraMODs;                ///< Flag indicating to put Theta==0 and D'==0 max models in MOD file.
  int dropEmptyClasses;         ///< Flag to indicate that empty liability classes should be dropped
  int physicalMap;              ///< Set when Map is read to indicate physical positions are available

  /* Storage and default names for files that are always opened (depending on analysis options) */
  char markerfile[PATH_MAX];  ///< Marker (frequency) file
  char mapfile[PATH_MAX];     ///< Map file
  char pedfile[PATH_MAX];     ///< Pedigree file
  char datafile[PATH_MAX];   ///< Data (pedigree description) file
  char avghetfile[PATH_MAX];       ///< Bayes Ratio file
  char pplfile[PATH_MAX];         ///< PPL file
  char condFile[PATH_MAX];      ///< Conditional LR file


  /* Storage for files that are only opened based on explicit directives */
  char ccfile[PATH_MAX];                 ///< Case control count file
  char modfile[PATH_MAX];                ///< MOD and maximizing model file
  char maxmodelfile[PATH_MAX];           ///< verbose Max Model file, obsolete?
  char intermediatefile[PATH_MAX];       ///< Intermediate Result file
  char dkelvinoutfile[PATH_MAX];         ///< DCHURE detail file
  char resultsprefix[PATH_MAX]; ///< Path for SR directive result storage
  
} ModelOptions;

#endif

