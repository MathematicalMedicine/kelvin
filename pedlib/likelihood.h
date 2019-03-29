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

  /* link to next parental pair for this nuclear family on this locus */
  struct ParentalPair *pNext;
} ParentalPair;

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


int construct_parental_pair(NuclearFamily *pNucFam, Person *pProband, 
		            int locus);
#endif
