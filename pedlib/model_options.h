#ifndef __PROGRAM_OPTION_H__
#define __PROGRAM_OPTION_H__

typedef struct ModelOptions {
  /* analysis type - Linkage Equalibrium (LE) or Linkage Disequalibrium (LD) */
  int equilibrium; /* WAS analysisOption */

  int polynomial; /* missing? TRUE if using polynomial, I guess */

  /* initially we only handle linkmap ? */
  /* int analysisOption; */ /* Unused */

  /* whether the chromosome is X chromosome 
   * For X chromosome, male (XY) 's X chromosome will be doubled as 
   * homozygotes, so to handle both males and females the same way 
   *
   * was (wrongly) called secChromosome in likelihood.c */
  int sexLinked;

  /* unknown trait value */
  double unknownTraitValue;

  /* unknown person ID */
  char *sUnknownPersonID;

  /* set recoding (super alleles) or not */
  /* int alleleSetRecodingFlag; */ /* Unused */

} ModelOptions;

extern ModelOptions modelOptions;

#endif /* __PROGRAM_OPTION_H__ */
