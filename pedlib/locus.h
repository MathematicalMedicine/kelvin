
/**********************************************************************
 * Copyright 2010, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#ifndef __LOCUS_H__
#define __LOCUS_H__

#include "../config/lambda_cell.h"
#include "../utils/polynomial.h"

#define ERROR_MARGIN            1.0E-9

#define LOCUS_TYPE_MARKER	0
#define LOCUS_TYPE_TRAIT	1

/* these affect and limit the size of pentrance matrix */

/* When we do categorial traits, there are potentially more than 3 possible values */
#define MAX_NUM_AFFECTION_STATUS 	3
#define MAX_NUM_LIABILITY_CLASS		50

/* usually only 2 alleles for a trait locus */
#define MAX_NUM_TRAIT_ALLELE		5

/* max number of trait locus genotypes - unphased */
#define MAX_NUM_TRAIT_GENOTYPE		(MAX_NUM_TRAIT_ALLELE * (MAX_NUM_TRAIT_ALLELE + 1) /2 )

/* max number of trait variables - usually 1, or bi-variate ...
 * can't see it go up to 3 or more */
#define MAX_NUM_TRAIT			2

/* covariance between trait variables */
#define MAX_NUM_TRAIT_COVARIANCE	(MAX_NUM_TRAIT* (MAX_NUM_TRAIT-1) /2 )


#define DEF_LOCUS_MALLOC_INCREMENT	10
#define MAX_LOCUS_NAME_LEN		30

#define LD_E				0.000001

/* map function */
#define MAP_FUNCTION_KOSAMBI		0
#define MAP_FUNCTION_HALDANE 		1

enum MapPosIndex {
  MAP_POS_SEX_AVERAGE,
  MAP_POS_MALE,
  MAP_POS_FEMALE,
  SIZE_MAP_POS
};

#define DIRECTION_LEFT                  0
#define DIRECTION_RIGHT                 1

/* PI */
#ifndef PI
#define PI				3.1415926535898e0
#endif

/* Genotype structure to store a pair of alleles represents a genotype
 * A Person has a list of possible genotypes for each locus - marker or disease
 * this Genotype is phased!!!
 * Original Vitesse structure is: v_v_glist & GLIST
 * For each locus, the possible genotypes are linked together
 * through a linklist. 
 * */
typedef struct Genotype
{
  /* paternal and maternal allele number. 
   * an integer should be big enough to handle 
   * any locus which can have huge number of alleles */
  short allele[2];
  /* bit is set for the corresponding paternalAllele
   * example: allele=3, then pAlleleBits[0]=8 
   * since an interger can only handle 32 alleles represented by bit 
   * mask, so potentially you need an array of pbits, so pointer is 
   * used. The size of the array of pbits is determined by 
   * len_set_recode_vector in the locus structure 
   * the value should be determined by the maximum number of alleles
   * for that locus plus some offset for super alleles 
   * */
  unsigned int *pAlleleBits[2];

  /* Link to next possible genotype for the same marker (Note: NOT for 
   * next marker genotype) 
   * through this pointer, current active genotypes are linked together 
   */
  struct Genotype *pNext;

  /* this pointer links the master list of genotypes together 
   * when doing loops, with fixed loop breaker genotype, a subset of active genotypes 
   * are linked through pNext, but the original master list is reflected in pSavedNext */
  struct Genotype *pSavedNext;

  /* This is used during likelihood calculation 
   * under one parental pair, only certain subset of the original 
   * genotype lists of the child's is valid. this subset is 
   * linked through pShadowNext */
  struct Genotype *pShadowNext;

  /* this is for likelihood calculation 
   * 1 - matches paternal
   * 2 - matches maternal
   * 3 - matches both
   */
  short inheritance[2];

  /* index on the list of the genotypes 
   * The value is populated when the initial set recoding and genotype
   * elimination are done, but before likelihood calculations 
   * The value is used when storing and retrieving multi locus 
   * conditional genotype likelihood */
  short position;

  /* the pointer and position of the genotype that is only different than 
   * this one by phase */
  short dualPosition;
  struct Genotype *pDualGenotype;

  /* this field can have different flags for different things, such
   * as homozygote bit, save bit, delete bit, ....
   * */
  int flagField;

  union
  {
    /* penetrance factor base on this person's phenotype and this genotype 
     * for marker locus, this is always 1 */
    double penetrance;
    Polynomial *penetrancePolynomial;
  } penslot;

  union
  {
    /* for a homozygous genotype, the weight is p*p
     * for a heterozygous unphased genotype, the weight is 2pq
     * but in this implementation, genotype is phased, so the weight is pq
     * where p and q are the allele frequencies */
    double weight;
    Polynomial *weightPolynomial;
  } wtslot;

} Genotype;

/* structure for allele set - super allele */
typedef struct AlleleSet
{
  /* identifier for the allele set
   * for singleton allele set, the allele ID matches the allele no. 
   * for non-singleton allele set, it is higher than the max 
   * number of alleles for that locus */
  int alleleID;
  /* number of alleles in the set */
  int numAllele;
  /* bit vector for the alleles representing through bit mask
   * allele 3 is represented by setting bit 3 of alleleBits[0]
   * this dynamically allocated space allows more than 32 alleles
   * represented by bit mask */
  unsigned int *pAlleleBits;
  /* max freq of the constituent allele in the set */
  double maxFreq;
  /* the sum of freqs of all the constituent allees in the set */
  double sumFreq;
  Polynomial *sumFreqPolynomial;
  /* an array of the alleles */
  int *pAlleles;

  /* next allele set 
   * usually paternal alleles are linked together
   * and maternal alleles are linked together during 
   * allele set recoding and genotype elimination */
  struct AlleleSet *pNext;
} AlleleSet;

/* each locus can be for marker, for trait and the trait can be 
 * of different value type - binary, quantitative etc.
 * trait can also have liability class associated with it
 * This struct assumes all these loci are on the same chromosome
 * */
typedef struct LocusList
{
  /* number of loci in this list */
  int numLocus;
  /* number of trait loci in this list */
  int numTraitLocus;
  /* maximum number of integers to represent alleles in bit mask for
   * any locus */
  int alleleSetLen;
  /* an array of pointers to the loci in the list */
  struct Locus **ppLocusList;

  /* LD frequencies between two loci */
  struct LDLoci *pLDLoci;
  int numLDLoci;
  int maxNumLDLoci;

  /* internal counter for memory allocateion */
  int maxNumLocus;

} LocusList;

/* A subset of the above LocusList */
typedef struct SubLocusList
{
  int numLocus;
  /* trait locus index in this locus list. If no trait locus, this index should be -1 */
  int traitLocusIndex;
  /* trait locus index in the original locus list. If no trait locus, this index should be -1 */
  int traitOrigLocus;
  /* an array of index of the original LocusList index 
   * for example: 2, 1, 4 means
   * original locus 2, 1 and 4 in such given order */
  int *pLocusIndex;
  /* distance is expressed in recombination fraction */
  double *pPrevLocusDistance[3];
  double *pNextLocusDistance[3];
  struct polynomial *pPrevLocusDistancePolynomial[3];
  struct polynomial *pNextLocusDistancePolynomial[3];

  /* internal counter for memory allocation */
  int maxNumLocus;
} SubLocusList;


/* Locus structure for each locus - marker or trait */
typedef struct Locus
{
  /* type of locus - trait or marker */
  short locusType;

  /* number of alleles for this locus including super alleles */
  short numAllele;

  /* number of original alleles for this locus without super alleles */
  short numOriginalAllele;
  /* number of integers required to set allele bits for this locus
   * if we need to do set recoding 
   * if < 32 alleles, just 1 */
  short alleleSetLen;

  /* locus name */
  char sName[MAX_LOCUS_NAME_LEN];

  /* assumed frequencies of these alleles 
   * pointer to an array of doubles of the actual frequencies */
  double *pAlleleFrequency;
  struct polynomial *pAlleleFrequencyPolynomial;
  /* actual count of these alleles in the pedigree data provided */
  short *pAlleleCount;
  /* original allele names */
  char **ppAlleleNames;

  /* allele set */
  AlleleSet **ppAlleleSetList;
  short numAlleleSet;
  /* internal memeory allocation counter */
  short maxNumAlleleSet;

  /* code vector ???? */
  //int **ppCodeVector;
  //int numCodeSize;

  /* if this is a disease locus, the following apply, otherwise ignore */
  struct TraitLocus *pTraitLocus;

  /* if this is a marker locus, then we have marker map information */
  struct MapUnit *pMapUnit;

  /* mark this locus is in LD with another locus */
  int LDLocus;

} Locus;

/* Trait structure is only for disease locus */
typedef struct TraitLocus
{
  /* map position in terms of cM - applies to trait only in multipoint analysis */
  double mapPosition[SIZE_MAP_POS];
  /* number of trait variables for this trait locus 
   * i.e. all the traits are assumed to associated with single 
   * disease gene (locus) */
  int numTrait;
  struct Trait *pTraits[MAX_NUM_TRAIT];
  //double covariance[MAX_NUM_TRAIT_COVARIANCE][MAX_NUM_TRAIT_GENOTYPE];
  double covariance[MAX_NUM_TRAIT][MAX_NUM_TRAIT][MAX_NUM_TRAIT_ALLELE]
    [MAX_NUM_TRAIT_ALLELE];
} TraitLocus;

typedef struct Trait
{
  /* trait locus value type - Binary or Quantitative(Continuous) or combined */
  int type;
  /* number of liability class - default should be 1 i.e. everyone 
   * has the same risk - only applicable for affection status traits */
  int numLiabilityClass;
  /* 4-dimension penetrance matrix for affection status trait
   * 1st - affection code
   * 2nd - liability class (usually only 1, so the index is 0 )
   * 3rd - allele 1
   * 4th - allele 2
   * */
  double penetrance[MAX_NUM_AFFECTION_STATUS][MAX_NUM_LIABILITY_CLASS]
    [MAX_NUM_TRAIT_ALLELE][MAX_NUM_TRAIT_ALLELE];
  /* 3-dimension matrix for quantitative trait
   * means & standard deviation
   * */
  double means[MAX_NUM_LIABILITY_CLASS][MAX_NUM_TRAIT_ALLELE]
    [MAX_NUM_TRAIT_ALLELE];
  /* standard deviation */
  double stddev[MAX_NUM_LIABILITY_CLASS][MAX_NUM_TRAIT_ALLELE]
    [MAX_NUM_TRAIT_ALLELE];
  /* unknown trait value used in the input pedigree file */
  double unknownTraitValue;
  /* flag indicating using the combined trait model
   * LESS indicates using the CDF(c), the area less than c,
   * MORE indicates using 1-CDF(c), the area more than c
   * For example, if lessCutoffFlag is -88.88, and an individual's trait
   * value is -88.88, the we need to determine the conditional penetrance
   * value given genotype using CDF(c)
   * NOTE: cutoff is another parameter which will be given under each locus */
  double lessCutoffFlag;
  double moreCutoffFlag;
  /* the cut off value */
  double cutoffValue[MAX_NUM_LIABILITY_CLASS];
  /* QT function flag : 0-normal, 1-t-distribution */
  int functionQT;
  /* degree of freedom for t-distribution */
  double dfQT;
  double sampleMean;
  double sampleSD;
  double min;
  double max;
  double minFlag;
  double maxFlag;
} Trait;

/* structure to keep marker map information */
typedef struct MapUnit
{
  /* index of this MapUnit in the map list */
  int mapIndex;
  /* chromosome */
  int chromosome;
  /* marker name */
  char sName[MAX_LOCUS_NAME_LEN];
  /* D number - ex. D1S243 - chromosome 1, segment 243?? */
  char sDNumber[MAX_LOCUS_NAME_LEN];
  /* male - 1
   * female - 2
   * sex average map position of this marker - 0 */
  double mapPos[3];
  /* base pair location if known  
   * the longest chromosome may have around 250 million base pairs, so a regular int can hold it
   * Be careful, if you do genome wide base pair location, then you ought to use "long int"
   */
  int basePairLocation;
  /* distance to previous marker in cM */
  //  double prevDistance;
  /* distance to next marker in cM */
  // double nextDistance;
} MapUnit;

typedef struct Map
{
  /* map function indication */
  int mapFunction;

  /* number of markers in this map */
  int count;
  /* pointer to an array of pointers to the actual maker info */
  MapUnit **ppMapUnitList;

  /* internal counter for memory allocation of the ppMapUnitList */
  int maxUnit;
} Map;

/* structure to keep the LD haplotype frequencies between two loci */
typedef struct LDLoci
{
  /* the locus number should be the original locus number */
  int locus1;
  int locus2;
  int numAllele1;		/* number of alleles for locus 1 */
  int numAllele2;		/* number of alleles for locus 2 */
  /* D primes between alleles of the two loci 
   * size = (m-1) * (n-1) if the two loci have m and n alleles respectively */
  double **ppDPrime;
  double **ppDValue;
  /* haplotype frequency between alleles of the two loci 
   * size = m * n */
  double **ppHaploFreq;
} LDLoci;

/* conditional likelihood of observing the corresponding nuclear family
 * conidtional on the multi locus genotype of this person */
typedef struct ConditionalLikelihood
{
  union
  {
    double likelihood;
    struct polynomial *likelihoodPolynomial;
  } lkslot;
  /* possibility of observing this multi locus genotype 
   * with parents - transmission probability * penetrance
   * without parents - genotype possibility * penetrance 
   * */
  union
  {
    double weight;
    struct polynomial *weightPolynomial;
  } wtslot;
  union
  {
    /* to save for likelihood calculations with only phase differences */
    double tmpLikelihood;
    struct polynomial *tmpLikelihoodPolynomial;
  } tmpslot;
  short touchedFlag;
  short tmpTouched;
} ConditionalLikelihood;

/* global function prototypes */
int read_mapfile (char *sMapfileName, int sexspecific);
int read_datafile (char *sMarkerfileName);
int read_markerfile (char *sMarkerfileName, int requiredMarkerCount);
int add_allele (Locus * pLocus, char *sAlleleName, double freq);
Trait *add_trait (int trait, TraitLocus * pTraitLocus, int traitType);
int create_baseline_trait_genotypes (int locus, Pedigree * pPedigree);
int create_baseline_marker_genotypes (int locus, Pedigree * pPedigree);
void print_person_locus_genotype_list (Person * pPerson, int locus);
void print_pedigree_locus_genotype_list (Pedigree * pPedigree, int locus);
void print_pedigree_locus_genotype_count (Pedigree * pPedigree, int locus);
int remove_genotype (Genotype ** pHead, Genotype * pGenotype, int *pCount);
LDLoci *find_LD_loci (int locus1, int locus2);
int allocate_multi_locus_genotype_storage (Pedigree * pPedigree,
					   int numLocus);
int count_multi_locus_genotype_storage (Pedigree * pPedigree);
int free_multi_locus_genotype_storage (Pedigree * pPedigree);
int initialize_multi_locus_genotype (Pedigree * pPedigree);
int set_genotype_weight (Pedigree * pPedigree, int locus);
int set_genotype_position (Pedigree * pPedigree, int locus);
int initialize_loci (PedigreeSet * pPedigreeSet);
int update_locus (PedigreeSet * pPedigreeSet, int locus);
int update_penetrance (PedigreeSet * pPedigreeSet, int locus);
double cm_to_recombination_fraction (double distance, int mapFunctionFlag);

/* int setup_LD_haplotype_freq (LDLoci * pLDLoci); */
int setup_LD_haplotype_freq (LDLoci * pLDLoci, LambdaCell * pCell,
			     int dprimeIdx);
int isDPrime0 (double **ppDrime, int m, int n);
double *get_map_position (int locus);
int add_analysis_locus (SubLocusList * pLocusList, int locus,
			int directionFlag, int mapFlag);
int add_markers_to_locuslist (SubLocusList * pLocusList, int numMarkers,
			      int *pLeftMarker, int start, int end,
			      double traitPosition, int mapFlag);
void free_sub_locus_list (SubLocusList * pLocusList);
void final_cleanup ();
void free_pedigree_set (PedigreeSet * pPedigreeSet);

int find_allele (int locus, char *sAlleleName);
int set_allele_bit (int allele, unsigned int *pBits);
int construct_original_allele_set_list (int locus);
int free_LD_loci (LDLoci * pLocus);
int reallocate_LD_loci (LDLoci * pLocus, int m, int n);
int set_null_dprime (LDLoci * pLocus);
int copy_dprime (LDLoci * pLocus, double **pSrc);
int copy_haploFreq (LDLoci * pLocus, double **pSrc);
int copy_DValue (LDLoci * pLocus, double **pSrc);
int find_locus (LocusList * pLocusList, char *sName);
int interpolate_physical_location (double pos);

void populate_pedigree_saved_genotype_link (int locus, Pedigree * pPed);
void populate_saved_genotype_link (PedigreeSet * pSet);
void restore_pedigree_genotype_link_from_saved (Pedigree * pPed);
void set_removeGenotypeFlag (int flag);

/* global variable */
extern Map map;

/* a master locus list with both marker and trait loci in */
extern LocusList originalLocusList;

/* The analysis locus list - this will be different than the originalLocusList in MP */
extern SubLocusList *analysisLocusList;
extern SubLocusList savedLocusList;

/* for MP null hypothesis */
extern SubLocusList nullLocusList;

extern LDLoci *pLDLoci;

#endif /* __LOCUS_H__ */
