/* avoid structures re-define because of multiple inclusion of the files */
#ifndef __PEDIGREE_H__
#define __PEDIGREE_H__

#ifndef NO_POLYNOMAIL
#include "polynomial.h"
#endif
/* Most of the structures defined in this file were taken from the original
 * vitesse codeset v_prog.h and trimed down to the bare neccessities
 * */

/* maximum line length in pedigree file */
#define MAX_LINE_LEN            8192
/* maximum length for pedigree label, individual IDs */
#define MAX_PED_LABEL_LEN       128
/* delimiter to separate fields in the pedfile */
#define PED_DELIM               " \t\n"
/* default pedigree malloc increment */
#define DEF_PED_MALLOC_INCREMENT 10

/* A person's sex */
#define MALE 		1
#define FEMALE		2

/* direction of the peeling */
#define PEDIGREE_UP			0
#define PEDIGREE_DOWN			1

/* parents */
//typedef enum { DAD=0, MOM=1} PARENTS;

#define DADID	"DAD"
#define MOMID	"MOM"

#define PATERNAL 	DAD
#define MATERNAL 	MOM

/* pedigree list data structure 
 * This structure groups all the pedigrees together by arrays of pointers 
 * */
typedef struct PedigreeSet {
  /* Number of pedigrees */
  int numPedigree;
  /* points to an array of pointers to pedigrees */
  struct Pedigree **ppPedigreeSet;
  /* if any of the pedigree has loop */
  int loopFlag;

  /* max number of founders in any pedigree of the pedigree set 
   * not used right now, will we ever??? */
  int maxNumFounder;

  /* if fields such as mom or dad matches this, it means that person's
   * parent(mom or dad) is unknown */
  char sUnknownID[MAX_PED_LABEL_LEN];

  /* likelihood holder */
  double *nullLikelihood;
  double likelihood;
  double log10Likelihood;
  /* assuming alpha from 0 to 1 with step size 0.05 */
  double hetLR[21];
  double log10HetLR[21];

#ifndef NO_POLYNOMIAL
  struct polynomial *likelihoodPolynomial;
#endif

  /* for internal memory allocation tracking purpose 
   * This is the number of pedigrees we have allocated space for the list */
  int maxNumPedigree;

}PedigreeSet;
/* Pedigree structure. Each pedigree is represented by this structure
 * It has information for each pedigree such as:
 *   number of persons 
 *   number of nuclear families
 *   number of founders etc. 
 * */
typedef struct Pedigree {
  /* Links to the pedigree set this pedigree is in */
  PedigreeSet *pPedigreeSet;
  /* Pedigree Index for the pedigree list in the pedigree set */
  int pedigreeIndex;

  /* Pedigree ID */
  char sPedigreeID[MAX_PED_LABEL_LEN];
  /* original pedigree lable 
   * if we are working on a post-makeped file, makeped may have changed
   * (renumbered) individuals in the pre-makeped file
   * need to keep track of the original pedigree lable which reflects
   * the pre-makeped pedigree lable for report, diagnostic purpose
   * */
  char sOriginalID[MAX_PED_LABEL_LEN];
  
  /* Number of individual in this pedigree */
  int numPerson;
  /* All individuals in this pedigree 
   * Points to an array of pointers to individuals */
  struct Person **ppPersonList;

  /* Number of founders who don't have parents in the pedigree */
  int numFounder;
  /* Founder list */
  struct Person **ppFounderList;
  /* internal memory allocation counter */
  int maxNumFounder;

  /* Number of nuclear families within this pedigree 
   * Nuclear families could overlap with connected individuals */
  int numNuclearFamily;
  /*  Points to an array of pointers to nuclear families */
  struct NuclearFamily **ppNuclearFamilyList;
  /* Number of nuclear families with founder parents */
  int numFounderNuclearFamily; 
  /* Founder nuclear family list */
  struct NuclearFamily **ppFounderNuclearFamilyList;

  /* internal tracking for memory allocation */
  int maxNumFounderNuclearFamily;

  /* 1 - if there is loop in this pedigree   0 - no loop 
   * loop is present if there are at least two paths between any pair
   * of individuals */
  int loopFlag;
  /* Number of loops in this pedigree */
  int numLoop;
  /* Number of loop breakers in this pedigree 
   * One same person can be in two, three... loops 
   * but it should only be counted as one loop breaker */
  int numLoopBreaker;

  /* Peeling will go to this selected person */
  struct Person *pPeelingProband;
  /* Peeling direction for this pedigree: PEDIGREE_UP or PEDIGREE_DOWN */
  int peelingDirection;
  /* peeling nuclear family - the one we start and end with */
  struct NuclearFamily *pPeelingNuclearFamily;

  /* likelihood holder */
  double likelihood;
#ifndef NO_POLYNOMIAL
  struct polynomial *likelihoodPolynomial;
  struct polyList   *likelihoodPolyList;
#endif

  /* Internal counters for memory allocation */
  int maxNumPerson;
  int maxNuclearFamily;
  
}Pedigree;

/* PERSON is a structure of individual information and links related 
 * individuals together in a pedigree 
 * Each person in a pedigree is representeted in this structure 
 * Matching Vitesse structure is: PERSON (v_v_person)
 * */
typedef struct Person {
  /* links to the pedigree this person belongs to */
  struct Pedigree *pPedigree;
  /* person index in the person array for this pedigree */
  int personIndex;

  /* individual ID - memory will be allocated as the person is populated */
  char sID[MAX_PED_LABEL_LEN];
  /* father and mother individual ID 
   * memory will be allocated dynamically as this person is populated */
 // char sDadID[MAX_PED_LABEL_LEN];
 // char sMomID[MAX_PED_LABEL_LEN];
  char sParentID[2][MAX_PED_LABEL_LEN];
  /* Links to parents */ 
//  struct Person *pDad, *pMom;
  struct Person *pParents[2];
  /* sex of this individual MALE or FEMALE - macro defined */
  int sex;
  /* affection status of this individual if the trait is dichotomous
   * AFFECTION_STATUS_UNKNOWN
   * AFFECTION_STATUS_UNAFFECTED
   * AFFECTION_STATUS_AFFECTED
   * */
  double **ppTraitValue; 
  /* mark whether the trait value is known or not 
   * at each locus and each trait - note for trait locus, there might 
   * be more than 1 trait associated with the locus */
  int **ppTraitKnown;
  
  /* Some analysis are done by grouping individuals based on some risk
   * factors etc. Number of liability class is really case dependent
   * Usually less than 5 for a moderate size of pedigree files */
  int **ppLiabilityClass; 
 
  /* Links to first offspring, other offsprings can be retrieved from
   * this offspring's maternal siblings if current person is the mom or 
   * this offspring's paternal if current person is the dad */ 
  struct Person *pFirstChild;
  char sFirstChildID[MAX_PED_LABEL_LEN];
  /* Links to next parternal and maternal siblings 
   * A full sibling will appear in both paternal and maternal sibiling links
   * A half sibling will only show up in one
   * */ 
  //struct Person *pNextPaternalSib, *pNextMaternalSib;
  struct Person *pNextSib[2];
  //char sNextPaternalSibID[MAX_PED_LABEL_LEN];
  //char sNextMaternalSibID[MAX_PED_LABEL_LEN];
  char sNextSibID[2][MAX_PED_LABEL_LEN];

  /* spouse list (possible multiple marriages) 
   * Jeff has advised that no need to keep track of spouses
   * multiple spouses come into play during peeling through
   * nuclear families (connectors)
   * */
  struct Person **ppSpouseList;
  /* number of spouse - usually its 0 or 1 
   * > 1 if multiple marraiges */
  int numSpouse;
  /* internal meter for memory allocation */
  int maxNumSpouse;
  
  /* 1 - indicates this person is a founder (without parents ) */
  int founderFlag;
  
  /* If this person is marked as proband, the likelihood peelings will 
   * go to this person? 
   * 1 - peeling proband 
   * 2 - loop breaker */
  int proband; 
  /* A flag whether this person is doubled: 1-doubled 0:regular*/
  int loopBreaker;
  /* ???? original individual ID - this gives information whether this 
   * individual is doubled by examing other individuals in the pedigree
   * and see whether they share the same original individual IDs 
   * */
  char sOriginalID[MAX_PED_LABEL_LEN];
  /* for a loop breaker, it points to the original person
   * for others, this field is NULL */
  struct Person *pOriginalPerson; 
  

  /* phenotype pairs (paternal & maternal )for each locus */
  int *pPhenotypeList[2];
  /* whether the above phenotype is phased or not */
  int *pPhasedFlag;
  /* whether the person is typed at each locus */
  int *pTypedFlag;

  /* phased genotype list 
   * points to an array of pointers to phased genotype pairs
   * for each locus
   * genotypes at each locus for the same person are linked together */
  struct Genotype **ppGenotypeList;

  /* tally of number of genotypes for each locus */
  int *pNumGenotype;
 
  /* saved copy of the phased genotype list 
   * this is used on probands during likelihood calculation */
  struct Genotype **ppSavedGenotypeList;
  int *pSavedNumGenotype;
  
  /* This is used during likelihood calculation
   * under one parental pair, only certain subset of the original
   * genotype lists of the child's is valid. this subset is
   * linked through pShadowNext */
  struct Genotype **ppShadowGenotypeList;
  int *pShadowGenotypeListLen;
  
  /* These structures are used to keep track what alleles are seen transmitted
   * and what not in bit mask
   * one for maternal, one for paternal 
   * point to an array of array
   * inner array is for allele set - as more than 32 alleles is possible 
   * outer array is for each locus 
   * The contents are only useful for an untyped individual */
  int **ppTransmittedAlleles[2];
  /* this is the complement set of the above set with the whole
   * set being all the possible alleles for this person (not for the locus), 
   * so it's a smaller set */
  int **ppNonTransmittedAlleles[2];
  /* whether there is a need to do set recoding for this person */
  int recodeFlag; 
  
  /* number of nuclear family this person belongs to */
  int numNuclearFamily;
  /* an array of pointers to the nuclear families this person belongs to */
  struct NuclearFamily **ppNuclearFamilyList;
  /* internal counter to keep track the size of array allocated to store
   * the points to nuclear families */
  int maxNumNuclearFamily;

  /* multi locus conditional genotype likelihood 
   * It's a flatted multi dimensional array
   * For example 3-locus, for this person locus1 has 5 possible genotypes,
   * locus 2 has 3 possible genotypes and locus 3 has 2 possible genotypes,
   * then this flattened array should be allocated with 5x3x2 elements to
   * represent each possible multi locus genotype 
   * the value stored in each element is the conditional likelihood of
   * the multi locus genotype 
   * The array is allocated when the set recoding and genotype eliminations
   * are done, but before likelihood calculations start */
  struct ConditionalLikelihood *pLikelihood;
  int numConditionals;

  /* current haplotype */
  struct Genotype **ppHaplotype;
}Person;


/* Data structure for nuclear families - two parents and their children 
 * Nuclear families are weaved together through connected individuals
 * bi-directionally (parent in one family connected (pUpFamilies) to child 
 * (same person) in the connected family. And reversely through pDownFamilies)
 * This structure is used for peeling during likelihood calculation 
 * Original Vitesse structure is: v_v_nuc_fam and NUC_FAM */
typedef struct NuclearFamily {
  /* links back to the pedigree this nuclear family belongs to */
  struct Pedigree *pPedigree;
  /* index of this nuclear family */
  int nuclearFamilyIndex;
  
  /* parents of this nuclear family */
  //Person *pDad;
  //Person *pMom;
  Person *pParents[2];
  /* number of children in this nuclear family */
  int numChildren;
  /* A list of pointers to all the children in this family */
  Person **ppChildrenList;
  /* Internal tracking for space allocation */
  int maxNumChildren; 

  /* individual we are peeling to 
   * how to determine who is proband for each nuclear family? 
   * parents or children?
   * not setting it yet 
   * */
  Person *pProband;
  /* family ID - not setting these yet - do we need them? */
  int familyID;
  //int momID;
  //int dadID;
  int parentID[2];
  /* point to a list of children IDs - not setting them yet */
  int *childrenID;
  int probandID;
  int probandSpouseID;

  /* When peeling is done with this family, flag is set to 1, otherwise 0 ?? */
  int doneFlag;

  /* Connected families through parents */
  struct NuclearFamilyConnector *pUpConnectors;
  /* Connected families through children */
  struct NuclearFamilyConnector *pDownConnectors;

  /* for likelihood calculation */
  /* parental pair for each locus 
   * linked list for each locus */
  struct ParentalPair **ppParentalPair;
  int *pNumParentalPair;
  /* haplotype pairs */
  struct HaplotypePair *pHaplotypePairs;
  int numHaplotypePairs;

  double likelihood;
#ifndef NO_POLYNOMIAL
  Polynomial *likelihoodPolynomial;
#endif

}NuclearFamily;

typedef struct NuclearFamilyConnector {
  /* Through this person, two nuclear families are connected */
  Person *pConnectedPerson;
  /* Connected nuclear family */
  NuclearFamily *pConnectedNuclearFamily;
  /* A family can have multiple individuals connected to other nuclear
   * families, thus multiple connectors */
  struct NuclearFamilyConnector *pNextConnector;
}NuclearFamilyConnector;


int print_nuclear_family(FILE *fp, Pedigree *pPed);

#endif /* __PEDIGREE_H__ */
