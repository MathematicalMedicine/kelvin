#ifndef __PROGRAM_OPTION_H__
#define __PROGRAM_OPTION_H__

/* Analysis type */
#define ANALYSIS_LINKAGE_EQUILIBRIUM            0
#define ANALYSIS_LINKAGE_DISEQUILIBRIUM         1

typedef struct ProgramOption{
  /* analysis type - Linkage Equalibrium (LE) or Linkage Disequalibrium (LD) */
  int analysisType;

  /* initially we only handle linkmap ? */
  int analysisOption;

  /* whether the chromosome is X chromosome 
   * For X chromosome, male (XY) 's X chromosome will be doubled as 
   * homozygotes, so to handle both males and females the same way */
  int sexLinked;

  /* unknown trait value */
  double unknownTraitValue;

  /* unknown person ID */
  char *sUnknownPersonID;

  /* set recoding (super alleles) or not */
  int alleleSetRecodingFlag;

}ProgramOption;

extern ProgramOption programOption;

#endif /* __PROGRAM_OPTION_H__ */
