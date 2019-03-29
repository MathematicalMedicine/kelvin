#ifndef __PARENTALPAIR_H__
#define __PARENTALPAIR_H__

#include "polynomial.h"
/* each nuclear family can have many ParentalPair depends on the 
 * possible genotypes of the father and mother 
 * these parental pairs are phased single locus . */
typedef struct ParentalPair{
  /* a pair of genotypes of the two parents */
  struct Genotype *pGenotype[2];
  /* each child's compatible genotype list 
   * points to an array of pointers, each element represents each child
   * points to an array of pointers of genotypes - they represent
   * the genotype array */
  struct Genotype ***pppChildGenoList;
  int *pChildGenoLen;
  /* this array maintains the inheritance pattern of the child genotype 
   * for the corresponding genotype
   * 1 - matches paternal only 
   * 2 - matches maternal only
   * 3 - matches both paternal & maternal
   *
   * first dimension is the child index
   * second dimension is the genotype index 
   */
  int **ppChildInheritance[2];


  /* link to next parental pair for this nuclear family on this locus */
  struct ParentalPair *pNext;
} ParentalPair;

#if 0
/* multi-locus parental pair */
typedef struct HaplotypePair{
  /* parental pairs for each locus */
  struct ParentalPair **ppParentalPairs;
  /* likelihood of the nuclear family conditional on this haplotype pair */
  double likelihood;

#ifndef NO_POLYNOMAIL
  Polynomial *likelihoodPolynomial;
#endif

}HaplotypePair;
#endif

typedef struct ParentalPairSpace{
  /* actual number of parental pairs in this work space for each locus */
  int *pNumParentalPair;
  /* in each ParentalPair, there are information on number of children, their genotype list length */
  ParentalPair **ppParentalPair;
  /* this is used to remember the haplotype (indices of parental pair at each locus ) */
  int *pParentalPairInd;
  int numHaplotypePair;
  double likelihood;
#ifndef NO_POLYNOMAIL
  Polynomial *likelihoodPolynomial;
#endif

  /* internal counter for each locus*/
  int maxNumParentalPair;
  int maxNumChildren;
  int maxNumChildGenotype;
}ParentalPairSpace;


extern ParentalPairSpace parentalPairSpace;

int construct_parental_pair(NuclearFamily *pNucFam, Person *pProband, 
		            int locus);
int stat_parental_pair_workspace (PedigreeSet *);
int initialize_parental_pair_workspace(ParentalPairSpace *pSpace, int numLocus);
int allocate_parental_pair_workspace(ParentalPairSpace *pSpace, int numLocus);


#endif
