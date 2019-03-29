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
#include "utils.h"		/* for logging */
#include "tools.h"
#include "likelihood.h"
#include "genotype_elimination.h"

ParentalPairSpace parentalPairSpace;

int shadow_genotype_elimnation (Genotype * pGenotype1, Genotype * pGenotype2,
				NuclearFamily * pNucFam, int locus);

/* construct the parental pair list for one locus */
/* locus is the index in the sub locusList */
 int
construct_parental_pair (NuclearFamily * pNucFam, Person * pProband,
			 int locus)
{
  Genotype *pGenotype[2];
  ParentalPair *pPair;
  int numChildren;
  Genotype *pChildGeno;
  int i, j;
  int genoLen;
  int status;
  Person *pChild;
  Person *pParents[2];
  int origLocus = locusList->pLocusIndex[locus];
  int numPair;

  /* initialize first */
  parentalPairSpace.pNumParentalPair[locus] = 0;
  numPair = 0;

  for (i = DAD; i <= MOM; i++)
    {
      pParents[i] = pNucFam->pParents[i];
      pGenotype[i] = pParents[i]->ppGenotypeList[origLocus];
    }
  /* now find the genotype pairs */
  numChildren = pNucFam->numChildren;
  KLOG (LOGPARENTALPAIR, LOGDEBUG,
	  "Nuc Fam %d (Dad: %s - Mom: %s) Proband: %s origLocus %d, #ofChildren: %d\n",
	  pNucFam->nuclearFamilyIndex, pParents[DAD]->sID, pParents[MOM]->sID,
	  pProband->sID, origLocus, pNucFam->numChildren);
  while (pGenotype[DAD])
    {
      pGenotype[MOM] = pParents[MOM]->ppGenotypeList[origLocus];
      while (pGenotype[MOM])
	{
	  /* do genotype elimination conditional on the pair 
	   * to get compatible children genotype lists */
	  status =
	    shadow_genotype_elimnation (pGenotype[DAD], pGenotype[MOM],
					pNucFam, locus);
	  /* if any child is not compatible with the pair, 
	   * then reject this parental pair and go on to next parental pair */
	  if (status != -1)
	    {
	      /* this parental pair is a valid one, construct it */
	      pPair = &parentalPairSpace.ppParentalPair[locus][numPair];
	      pPair->pGenotype[DAD] = pGenotype[DAD];
	      pPair->pGenotype[MOM] = pGenotype[MOM];
	      /* increase the count of the parental pairs */
	      numPair++;
	      parentalPairSpace.pNumParentalPair[locus]++;
	      /* copy all children's valid genotype list into the parental pair's
	       * children genotype array */
	      KLOG (LOGPARENTALPAIR, LOGDEBUG,
		      "Parental Pair(%d) %s(%d, %d) && %s(%d, %d)\n",
		      parentalPairSpace.pNumParentalPair[locus],
		      pNucFam->pParents[DAD]->sID,
		      pGenotype[DAD]->allele[DAD],
		      pGenotype[DAD]->allele[MOM],
		      pNucFam->pParents[MOM]->sID,
		      pGenotype[MOM]->allele[DAD],
		      pGenotype[MOM]->allele[MOM]);
	      for (i = 0; i < numChildren; i++)
		{
		  pChild = pNucFam->ppChildrenList[i];
		  genoLen = pChild->pShadowGenotypeListLen[origLocus];
		  pPair->pChildGenoLen[i] = genoLen;
		  pChildGeno = pChild->ppShadowGenotypeList[origLocus];
		  KLOG (LOGPARENTALPAIR, LOGDEBUG,
			  "Child %s #ofGenotypes: %d\n", pChild->sID,
			  genoLen);
		  for (j = 0; j < genoLen; j++)
		    {
		      pPair->pppChildGenoList[i][j] = pChildGeno;
		      pPair->ppChildInheritance[DAD][i][j] = pChildGeno->inheritance[DAD];
		      pPair->ppChildInheritance[MOM][i][j] = pChildGeno->inheritance[MOM];
		      KLOG (LOGPARENTALPAIR, LOGDEBUG, "  (%d, %d)\n",
			      pChildGeno->allele[DAD],
			      pChildGeno->allele[MOM]);
		      pChildGeno = pChildGeno->pShadowNext;
		    }
		}

	    }
	  /* move on to next genotype */
	  pGenotype[MOM] = pGenotype[MOM]->pNext;
	}
      pGenotype[DAD] = pGenotype[DAD]->pNext;
    }

  return 0;
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
  for (i = 0; i < numChildren; i++)
    {
      pChild = pNucFam->ppChildrenList[i];
      pGenotype = pChild->ppGenotypeList[origLocus];
      pChild->ppShadowGenotypeList[origLocus] = NULL;
      pChild->pShadowGenotypeListLen[origLocus] = 0;
      while (pGenotype)
	{
	  if (is_parent_child_genotype_compatible
	      (origLocus, DAD, pGenotype1, pGenotype)
	      && is_parent_child_genotype_compatible (origLocus, MOM, pGenotype2,
						      pGenotype))
	    {
	      pGenotype->pShadowNext = pChild->ppShadowGenotypeList[origLocus];
	      pChild->ppShadowGenotypeList[origLocus] = pGenotype;
	      pChild->pShadowGenotypeListLen[origLocus] += 1;
	    }
	  pGenotype = pGenotype->pNext;
	}			/* end of looping one child's genotype list */
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
  int ped;  /* pedigree index */
  int fam;  /* nuclear family index */
  int locus;

  pParentalPairSpace = &parentalPairSpace;

  for(ped=0; ped < pPedigreeList->numPedigree; ped++)
    {
      pPedigree = pPedigreeList->ppPedigreeSet[ped];
      for(fam=0; fam < pPedigree->numNuclearFamily; fam++)
	{
	  pNucFam = pPedigree->ppNuclearFamilyList[fam];
	  for(locus=0; locus < originalLocusList.numLocus; locus++)
	    {
	      /* be generous to pre-allocating work space */
	      maxNumParentalPair = pNucFam->pParents[DAD]->pSavedNumGenotype[locus] *
		pNucFam->pParents[MOM]->pSavedNumGenotype[locus];
	      
	      if(maxNumParentalPair > parentalPairSpace.maxNumParentalPair)
		{
		  parentalPairSpace.maxNumParentalPair = maxNumParentalPair;
		}
	      
	      numChildren = pNucFam->numChildren;
	      if(numChildren > parentalPairSpace.maxNumChildren)
		{
		  parentalPairSpace.maxNumChildren = pNucFam->numChildren;
		}
	      
	      for(i=0; i < numChildren; i++)
		{
		  if(pNucFam->ppChildrenList[i]->pSavedNumGenotype[locus] > 
		     parentalPairSpace.maxNumChildGenotype)
		    {
		      pParentalPairSpace->maxNumChildGenotype = 
			pNucFam->ppChildrenList[i]->pSavedNumGenotype[locus];
		    }
		}
	    } /* move to next locus */
	} /* move to next nuclear family */
    } /* move to next pedigree */
  
  return 0;
}

/* Initialize this workspace for tracking what's the maximum number of parental pairs 
 * for each locus */
int
initialize_parental_pair_workspace(ParentalPairSpace *pSpace, int numLocus)
{
  pSpace->pNumParentalPair = (int *)calloc(sizeof(int), numLocus);
  pSpace->pParentalPairInd = (int *)calloc(sizeof(int), numLocus);
  pSpace->ppParentalPair = (ParentalPair **) calloc(sizeof(ParentalPair *), numLocus);
  pSpace->maxNumParentalPair = 0;
  pSpace->maxNumChildren = 0;
  pSpace->maxNumChildGenotype = 0;

  return 0;
}


/* once we have stat all nuclear families for all loci
 * we can allocate the max */
int
allocate_parental_pair_workspace(ParentalPairSpace *pSpace, int numLocus)
{
  int locus;
  int i, j;

  for(locus=0; locus < numLocus; locus++)
    {
      /* allocate parental pair space first */
      pSpace->ppParentalPair[locus] = (ParentalPair *)calloc(sizeof(ParentalPair), 
							     pSpace->maxNumParentalPair);
      /* allocate space for children genotype list */
      for(i = 0; i < pSpace->maxNumParentalPair; i++)
	{
	  pSpace->ppParentalPair[locus][i].pChildGenoLen = (int *) calloc(sizeof(int), 
									  pSpace->maxNumChildren);
	  /* first dimension is the index of children */
	  pSpace->ppParentalPair[locus][i].pppChildGenoList = (Genotype ***) calloc(sizeof(Genotype **), 
										 pSpace->maxNumChildren);
	  pSpace->ppParentalPair[locus][i].ppChildInheritance[DAD] = (int **) calloc(sizeof(int *), 
										pSpace->maxNumChildren);
	  pSpace->ppParentalPair[locus][i].ppChildInheritance[MOM] = (int **) calloc(sizeof(int *), 
										pSpace->maxNumChildren);
	  /* for each child allocate enough space for genotypes */
	  for(j=0; j < pSpace->maxNumChildren; j++)
	    {
	      /* next dimension is the index of the different genotypes */
	      pSpace->ppParentalPair[locus][i].pppChildGenoList[j] = (Genotype **) calloc(sizeof(Genotype *), 
										       pSpace->maxNumChildGenotype);
	      pSpace->ppParentalPair[locus][i].ppChildInheritance[DAD][j] = (int *) calloc(sizeof(int),
										      pSpace->maxNumChildGenotype);
	      pSpace->ppParentalPair[locus][i].ppChildInheritance[MOM][j] = (int *) calloc(sizeof(int),
										      pSpace->maxNumChildGenotype);
	    }
	}
   }


  return 0;
}

int
free_parental_pair_workspace(ParentalPairSpace *pSpace, int numLocus)
{
  int locus;
  int i, j;

  for(locus=0; locus < numLocus; locus++)
    {
      /* free space for children genotype list */
      for(i = 0; i < pSpace->maxNumParentalPair; i++)
	{
	  /* for each child allocate enough space for genotypes */
	  for(j=0; j < pSpace->maxNumChildren; j++)
	    {
	      /* next dimension is the index of the different genotypes */
	      free(pSpace->ppParentalPair[locus][i].pppChildGenoList[j]);
	      free(pSpace->ppParentalPair[locus][i].ppChildInheritance[DAD][j]);
	      free(pSpace->ppParentalPair[locus][i].ppChildInheritance[MOM][j]);
	    }
	  free(pSpace->ppParentalPair[locus][i].pChildGenoLen);

	  /* first dimension is the index of children */
	  free(pSpace->ppParentalPair[locus][i].pppChildGenoList);
	  free(pSpace->ppParentalPair[locus][i].ppChildInheritance[DAD]);
	  free(pSpace->ppParentalPair[locus][i].ppChildInheritance[MOM]);
	}
      /* allocate parental pair space first */
      free(pSpace->ppParentalPair[locus]);
   }

  free(pSpace->pNumParentalPair);
  free(pSpace->pParentalPairInd);
  free(pSpace->ppParentalPair);

  return 0;
}
