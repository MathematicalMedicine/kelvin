#ifndef __PROGRAM_OPTION_H__
#define __PROGRAM_OPTION_H__

typedef struct ModelOptions
{
  /* analysis type - Linkage Equalibrium (LE) or Linkage Disequalibrium (LD) */
  int equilibrium;

  int polynomial;		/* TRUE if using polynomial to represent likelihood */

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

  /* affection status: DT defaults are {0,1,2} while QT/CT defaults
   * are {NAN, -INFINITY, INFINITY}. values are cast to double in both
   * cases. */
  double affectionStatus[3];
#if FALSE
  /* Replaced with affectionStatus[]: 1/3/07 ams */
  /* unknown trait value */
  double unknownTraitValue;
#endif

  /* unknown person ID */
  char *sUnknownPersonID;

  /* set recoding (super alleles) or not - we always use set recoding to speed up */
  /* int alleleSetRecodingFlag; *//* Unused */

  /* marker analysis flag; either MARKERTOMARKER or ADJACENTMARKER */
  int markerAnalysis;

  /* for calculating PPLs */
  double thetaCutoff[2];	/* using step function for the weight of each point */
  double thetaWeight;		/* weight for theta less than the thetaCutoff */
  double prior;			/* prior probability of linkage */
  double LDprior;		/* prior probability of LD (D' not 0) given linkage and theta
				 *  is within thetaCutoff */

} ModelOptions;

extern ModelOptions modelOptions;

#endif /* __PROGRAM_OPTION_H__ */
