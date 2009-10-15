
/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

/* This file contains functions to facilitate likelihood calculations
 * using parental pair algorithm 
 * That is base on each multi-locus genotype pairs of the parents, 
 * calculate the likelihood of observing the children's phenotypes 
 * so the result of the likelihood is conditional likelihood based
 * on parents' genotypes */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "pedlib.h"
#include "pedigree.h"
#include "locus.h"
#include "../utils/utils.h"		/* for logging */
#include "likelihood.h"
#include "genotype_elimination.h"

ParentalPairSpace parentalPairSpace;

int shadow_genotype_elimnation (Genotype * pGenotype1, Genotype * pGenotype2,
				NuclearFamily * pNucFam, int locus);
void fill_parental_pair (int locus, int *numPair, NuclearFamily * pNucFam,
			 Genotype * pDad, Genotype * pMom,
			 int dadAdjust, int momAdjust);

/* construct the parental pair list for one locus */

/* locus is the index in the sub locusList */
int
construct_parental_pair (NuclearFamily * pNucFam, Person * pProband,
			 int locus)
{
  Genotype *pGenotype[2];
  Genotype *pFirstGenotype[2];
  int i;
  int status;
  Person *pParents[2];
  int origLocus = locusList->pLocusIndex[locus];
  int numPair;
  int adjust[2];
  int ourIndex[2];
  Genotype *initialGeno[2];
  int head, spouse;
  ParentalPair *pPair;

  /* initialize first */
  parentalPairSpace.pNumParentalPair[locus] = 0;
  numPair = 0;

  head = pNucFam->head;
  spouse = pNucFam->spouse;
  for (i = DAD; i <= MOM; i++) {
    pParents[i] = pNucFam->pParents[i];
    if (pParents[i]->loopBreaker >= 1 && pParents[i]->pParents[DAD] == NULL) {
      pGenotype[i] = pParents[i]->pOriginalPerson->ppGenotypeList[origLocus];
    } else
      pGenotype[i] = pParents[i]->ppGenotypeList[origLocus];
    pFirstGenotype[i] = pGenotype[i];
  }
  /* now find the genotype pairs */
  KLOG (LOGPARENTALPAIR, LOGDEBUG,
	"Nuc Fam %d (Dad: %s - Mom: %s) Proband: %s origLocus %d, #ofChildren: %d\n",
	pNucFam->nuclearFamilyIndex, pParents[DAD]->sID, pParents[MOM]->sID,
	pProband->sID, origLocus, pNucFam->numChildren);
  while (pGenotype[head]) {
    initialGeno[head] = pGenotype[head];
    pGenotype[spouse] = pFirstGenotype[spouse];
    while (pGenotype[spouse]) {
      /* do genotype elimination conditional on the pair 
       * to get compatible children genotype lists */
      status =
	shadow_genotype_elimnation (pGenotype[DAD], pGenotype[MOM],
				    pNucFam, locus);
      /* if any child is not compatible with the pair, 
       * then reject this parental pair and go on to next parental pair */
      if (status != -1) {
	adjust[head] = 0;
	ourIndex[head] = 0;
	initialGeno[spouse] = pGenotype[spouse];
	while ((ourIndex[head] < 2)
	       && (ourIndex[head] == 0
		   || (pGenotype[head]->pDualGenotype != NULL
		       && pGenotype[head]->pNext ==
		       pGenotype[head]->pDualGenotype))) {
	  if (ourIndex[head] != 0) {
	    pGenotype[head] = pGenotype[head]->pDualGenotype;
	    adjust[head] = 1;
	  }
	  ourIndex[spouse] = 0;
	  adjust[spouse] = 0;
	  pGenotype[spouse] = initialGeno[spouse];
	  while ((ourIndex[spouse] < 2)
		 && (ourIndex[spouse] == 0
		     || (pGenotype[spouse]->pDualGenotype != NULL
			 && pGenotype[spouse]->pDualGenotype ==
			 pGenotype[spouse]->pNext))) {
	    if (ourIndex[spouse] != 0
		&& pGenotype[spouse]->pDualGenotype != NULL) {
	      pGenotype[spouse] = pGenotype[spouse]->pDualGenotype;
	      if (ourIndex[spouse] != 0)
		adjust[spouse] = 1;
	    }
	    pPair = &parentalPairSpace.ppParentalPair[locus][numPair];
	    /* 0 - no change in phase  1 - flip of the original phase */
	    pPair->phase[head] = adjust[head];
	    pPair->phase[spouse] = adjust[spouse];
	    fill_parental_pair (locus, &numPair, pNucFam,
				pGenotype[DAD], pGenotype[MOM],
				adjust[DAD], adjust[MOM]);
	    ourIndex[spouse]++;
	  }
	  ourIndex[head]++;
	}

      } /* valid pair found */
      else {
	if (pGenotype[spouse]->pDualGenotype != NULL
	    && pGenotype[spouse]->pDualGenotype == pGenotype[spouse]->pNext)
	  pGenotype[spouse] = pGenotype[spouse]->pNext;
      }
      /* reset dad genotype pointer */
      pGenotype[head] = initialGeno[head];

      /* move on to next genotype */
      pGenotype[spouse] = pGenotype[spouse]->pNext;
    }
    if (pGenotype[head]->pDualGenotype != NULL
	&& pGenotype[head]->pDualGenotype == pGenotype[head]->pNext)
      pGenotype[head] = pGenotype[head]->pNext->pNext;
    else
      pGenotype[head] = pGenotype[head]->pNext;

  }

  return 0;
}

void
fill_parental_pair (int locus, int *numPair, NuclearFamily * pNucFam,
		    Genotype * pDad, Genotype * pMom,
		    int dadAdjust, int momAdjust)
{
  ParentalPair *pPair;
  Person *pChild;
  int genoLen;
  int i, j;
  Genotype *pChildGeno;
  int origLocus = locusList->pLocusIndex[locus];

  pPair = &parentalPairSpace.ppParentalPair[locus][*numPair];
  pPair->pGenotype[DAD] = pDad;
  pPair->pGenotype[MOM] = pMom;

  /* increase the count of the parental pairs */
  (*numPair)++;
  parentalPairSpace.pNumParentalPair[locus]++;
  /* copy all children's valid genotype list into the parental pair's
   * children genotype array */
  KLOG (LOGPARENTALPAIR, LOGDEBUG,
	"Parental Pair(%d) %s(%d, %d) && %s(%d, %d)\n",
	parentalPairSpace.pNumParentalPair[locus],
	pNucFam->pParents[DAD]->sID,
	pDad->allele[DAD], pDad->allele[MOM],
	pNucFam->pParents[MOM]->sID, pMom->allele[DAD], pMom->allele[MOM]);
  for (i = 0; i < pNucFam->numChildren; i++) {
    pChild = pNucFam->ppChildrenList[i];
    genoLen = pChild->pShadowGenotypeListLen[origLocus];
    pPair->pChildGenoLen[i] = genoLen;
    pChildGeno = pChild->ppShadowGenotypeList[origLocus];
    KLOG (LOGPARENTALPAIR, LOGDEBUG,
	  "Child %s #ofGenotypes: %d\n", pChild->sID, genoLen);
    for (j = 0; j < genoLen; j++) {
      pPair->pppChildGenoList[i][j] = pChildGeno;
      if (dadAdjust == 0)
	pPair->ppChildInheritance[DAD][i][j] = pChildGeno->inheritance[DAD];
      else
	pPair->ppChildInheritance[DAD][i][j] =
	  pChildGeno->inheritance[DAD] ^ 3;

      if (momAdjust == 0)
	pPair->ppChildInheritance[MOM][i][j] = pChildGeno->inheritance[MOM];
      else
	pPair->ppChildInheritance[MOM][i][j] =
	  pChildGeno->inheritance[MOM] ^ 3;


      KLOG (LOGPARENTALPAIR, LOGDEBUG, "  (%d, %d) inheritance (%d, %d) \n",
	    pChildGeno->allele[DAD], pChildGeno->allele[MOM], 
	    pPair->ppChildInheritance[DAD][i][j], pPair->ppChildInheritance[MOM][i][j]);
      pChildGeno = pChildGeno->pShadowNext;
    }
  }


}

/* locus is the index in the locusList */
int
shadow_genotype_elimnation (Genotype * pGenotype1, Genotype * pGenotype2,
			    NuclearFamily * pNucFam, int locus)
{
  Person *pChild;
  int numChildren = pNucFam->numChildren;
  int i;
  Genotype *pGenotype;
  int origLocus = locusList->pLocusIndex[locus];

  /* check each child's genotype list against the parental pair */
  for (i = 0; i < numChildren; i++) {
    pChild = pNucFam->ppChildrenList[i];
    pGenotype = pChild->ppGenotypeList[origLocus];
    pChild->ppShadowGenotypeList[origLocus] = NULL;
    pChild->pShadowGenotypeListLen[origLocus] = 0;
    while (pGenotype) {
      if (is_parent_child_genotype_compatible
	  (origLocus, DAD, pChild->sex, pGenotype1, pGenotype)
	  && is_parent_child_genotype_compatible (origLocus, MOM,
						  pChild->sex, pGenotype2,
						  pGenotype)) {
	pGenotype->pShadowNext = pChild->ppShadowGenotypeList[origLocus];
	pChild->ppShadowGenotypeList[origLocus] = pGenotype;
	pChild->pShadowGenotypeListLen[origLocus] += 1;
      }
      pGenotype = pGenotype->pNext;
    }				/* end of looping one child's genotype list */
    /* if we didn't find any compatible genotype for the current child, 
     * then return failure */
    if (pChild->ppShadowGenotypeList[origLocus] == NULL)
      return -1;
  }				/* done looping all children */

  /* we are here because every child has at least one compatible genotype */


  return 0;
}

/* allocate the parental pair list for one locus */

/* locus is the index in the originalLocusList */
int
stat_parental_pair_workspace (PedigreeSet * pPedigreeList)
{
  int maxNumParentalPair;
  int numChildren;
  int i;
  ParentalPairSpace *pParentalPairSpace;
  NuclearFamily *pNucFam;
  Pedigree *pPedigree;
  int ped;			/* pedigree index */
  int fam;			/* nuclear family index */
  int locus;
  int numGenotype[2];

  pParentalPairSpace = &parentalPairSpace;

  for (ped = 0; ped < pPedigreeList->numPedigree; ped++) {
    pPedigree = pPedigreeList->ppPedigreeSet[ped];
    for (fam = 0; fam < pPedigree->numNuclearFamily; fam++) {
      pNucFam = pPedigree->ppNuclearFamilyList[fam];
      for (locus = 0; locus < originalLocusList.numLocus; locus++) {
	/* be generous to pre-allocating work space */
	for (i = DAD; i <= MOM; i++) {
	  if (pNucFam->pParents[i]->loopBreaker >= 1
	      && pNucFam->pParents[i]->pParents[DAD] == NULL) {
	    numGenotype[i] =
	      pNucFam->pParents[i]->pOriginalPerson->pSavedNumGenotype[locus];
	  } else
	    numGenotype[i] = pNucFam->pParents[i]->pSavedNumGenotype[locus];
	}

	maxNumParentalPair = numGenotype[DAD] * numGenotype[MOM];

	if (maxNumParentalPair > parentalPairSpace.maxNumParentalPair) {
	  parentalPairSpace.maxNumParentalPair = maxNumParentalPair;
	}

	numChildren = pNucFam->numChildren;
	if (numChildren > parentalPairSpace.maxNumChildren) {
	  parentalPairSpace.maxNumChildren = pNucFam->numChildren;
	}

	for (i = 0; i < numChildren; i++) {
	  if (pNucFam->ppChildrenList[i]->pSavedNumGenotype[locus] >
	      parentalPairSpace.maxNumChildGenotype) {
	    pParentalPairSpace->maxNumChildGenotype =
	      pNucFam->ppChildrenList[i]->pSavedNumGenotype[locus];
	  }
	}
      }				/* move to next locus */
    }				/* move to next nuclear family */
  }				/* move to next pedigree */

  return 0;
}

/* Initialize this workspace for tracking what's the maximum number of parental pairs 
 * for each locus */
int
initialize_parental_pair_workspace (ParentalPairSpace * pSpace, int numLocus)
{
  CALCHOKE(pSpace->pNumParentalPair, sizeof (int), (size_t) numLocus, int *);
  CALCHOKE(pSpace->pParentalPairInd, sizeof (int), (size_t) numLocus, int *);
  CALCHOKE(pSpace->pChildGenoInd, sizeof (int), (size_t) numLocus, int *);
  CALCHOKE(pSpace->phase[DAD], sizeof (short), (size_t) numLocus, short *);
  CALCHOKE(pSpace->phase[MOM], sizeof (short), (size_t) numLocus, short *);
  CALCHOKE(pSpace->ppParentalPair, sizeof (ParentalPair *), (size_t) numLocus, ParentalPair **);
  pSpace->maxNumParentalPair = 0;
  pSpace->maxNumChildren = 0;
  pSpace->maxNumChildGenotype = 0;

  return 0;
}


/* once we have stat all nuclear families for all loci
 * we can allocate the max 
 * numLocus - the number of loci analyzing at one time
 */
int
allocate_parental_pair_workspace (ParentalPairSpace * pSpace, int numLocus)
{
  int locus;
  int i, j;

  for (locus = 0; locus < numLocus; locus++) {
    /* allocate parental pair space first */
    CALCHOKE(pSpace->ppParentalPair[locus], sizeof (ParentalPair), (size_t) pSpace->maxNumParentalPair, ParentalPair *);
    /* allocate space for children genotype list */
    for (i = 0; i < pSpace->maxNumParentalPair; i++) {
      CALCHOKE(pSpace->ppParentalPair[locus][i].pChildGenoLen, sizeof (int), (size_t) pSpace->maxNumChildren, int *);
      /* first dimension is the index of children */
      CALCHOKE(pSpace->ppParentalPair[locus][i].pppChildGenoList, sizeof (Genotype **), (size_t) pSpace->maxNumChildren, Genotype ***);
      CALCHOKE(pSpace->ppParentalPair[locus][i].ppChildInheritance[DAD], sizeof (int *), (size_t) pSpace->maxNumChildren, int **);
      CALCHOKE(pSpace->ppParentalPair[locus][i].ppChildInheritance[MOM], sizeof (int *), (size_t) pSpace->maxNumChildren, int **);
      /* for each child allocate enough space for genotypes */
      for (j = 0; j < pSpace->maxNumChildren; j++) {
	/* next dimension is the index of the different genotypes */
	CALCHOKE(pSpace->ppParentalPair[locus][i].pppChildGenoList[j], sizeof (Genotype *), (size_t) pSpace->maxNumChildGenotype, Genotype **);
	CALCHOKE(pSpace->ppParentalPair[locus][i].ppChildInheritance[DAD][j], sizeof (int), (size_t) pSpace->maxNumChildGenotype, int *);
	CALCHOKE(pSpace->ppParentalPair[locus][i].ppChildInheritance[MOM][j], sizeof (int), (size_t) pSpace->maxNumChildGenotype, int *);
      }
    }
  }


  return 0;
}

int
free_parental_pair_workspace (ParentalPairSpace * pSpace, int numLocus)
{
  int locus;
  int i, j;

  for (locus = 0; locus < numLocus; locus++) {
    /* free space for children genotype list */
    for (i = 0; i < pSpace->maxNumParentalPair; i++) {
      /* for each child allocate enough space for genotypes */
      for (j = 0; j < pSpace->maxNumChildren; j++) {
	/* next dimension is the index of the different genotypes */
	free (pSpace->ppParentalPair[locus][i].pppChildGenoList[j]);
	free (pSpace->ppParentalPair[locus][i].ppChildInheritance[DAD][j]);
	free (pSpace->ppParentalPair[locus][i].ppChildInheritance[MOM][j]);
      }
      free (pSpace->ppParentalPair[locus][i].pChildGenoLen);

      /* first dimension is the index of children */
      free (pSpace->ppParentalPair[locus][i].pppChildGenoList);
      free (pSpace->ppParentalPair[locus][i].ppChildInheritance[DAD]);
      free (pSpace->ppParentalPair[locus][i].ppChildInheritance[MOM]);
    }
    /* allocate parental pair space first */
    free (pSpace->ppParentalPair[locus]);
  }

  free (pSpace->phase[DAD]);
  free (pSpace->phase[MOM]);
  free (pSpace->pChildGenoInd);
  free (pSpace->pNumParentalPair);
  free (pSpace->pParentalPairInd);
  free (pSpace->ppParentalPair);

  return 0;
}
