#ifndef __PARENTALPAIR_H__
#define __PARENTALPAIR_H__
/* Copyright (C) 2006, 2010, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/* each nuclear family can have many ParentalPair depends on the 
 * possible genotypes of the father and mother 
 * these parental pairs are phased single locus . */
typedef struct ParentalPair
{
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

  /* the following two helped to remember phase differences between the previous parental pair if they 
   * are the same except the phase */
  short phase[2];

  /* link to next parental pair for this nuclear family on this locus */
  struct ParentalPair *pNext;
} ParentalPair;

typedef struct ParentalPairSpace
{
  /* actual number of parental pairs in this work space for each locus */
  int *pNumParentalPair;
  /* in each ParentalPair, there are information on number of children, their genotype list length */
  ParentalPair **ppParentalPair;
  /* this is used to remember the haplotype (indices of parental pair at each locus ) */
  int *pParentalPairInd;
  /* remembers the child genotype location */
  int *pChildGenoInd;
  short *phase[2];
  int numHaplotypePair;

  /* internal counter for each locus */
  int maxNumParentalPair;
  int maxNumChildren;
  int maxNumChildGenotype;
} ParentalPairSpace;

typedef struct XMission
{
  /* transmission probability */
  union
  {
    Polynomial *probPoly[3];
    double prob[3];
  } slot;
} XMission;

extern ParentalPairSpace parentalPairSpace;
extern XMission *xmissionMatrix;
extern double *half_pow;
extern char partialPolynomialFunctionName[];


int build_likelihood_polynomial (Pedigree * pPedigree);
int compute_likelihood (char *fileName, int lineNo, PedigreeSet * pPedigreeList);
int construct_parental_pair (NuclearFamily * pNucFam, Person * pProband,
			     int locus);
int stat_parental_pair_workspace (PedigreeSet *);
int initialize_parental_pair_workspace (ParentalPairSpace * pSpace,
					int numLocus);
int allocate_parental_pair_workspace (ParentalPairSpace * pSpace,
				      int numLocus);
int free_parental_pair_workspace (ParentalPairSpace * pSpace, int numLocus);
int allocate_likelihood_space (PedigreeSet * pPedigreeList, int numLocus);
int count_likelihood_space (PedigreeSet * pPedigreeList);
void free_likelihood_space (PedigreeSet * pPedigreeList);
void allocate_nucfam_het (PedigreeSet * pPedigreeList, int numLocus);
int build_xmission_matrix (XMission ** ppMatrix, int totalLoci);

int populate_xmission_matrix (char *fileName, int lineNo, XMission * pMatrix, int totalLoci,
			      void *prob[3], void *prob2[3],
			      void *hetProb[3],
			      int cellIndex,
			      int lastHetLoc, int prevPattern, int loc);

void
print_xmission_matrix (XMission * pMatrix, int totalLoci, int loc,
		       int cellIndex, char *pID);
void print_xmission_matrix_differences (XMission *, XMission *, int, int, int, char *);

#endif
