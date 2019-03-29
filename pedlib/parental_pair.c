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

int shadow_genotype_elimnation (Genotype * pGenotype1, Genotype * pGenotype2,
				NuclearFamily * pNucFam, int locus);
/* construct the parental pair list for one locus */
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

  for (i = DAD; i <= MOM; i++)
    {
      pParents[i] = pNucFam->pParents[i];
      pGenotype[i] = pParents[i]->ppGenotypeList[locus];
    }
  /* now find the genotype pairs */
  numChildren = pNucFam->numChildren;
  logMsg (LOGPARENTALPAIR, LOGDEBUG,
	  "Nuc Fam %d (Dad: %s - Mom: %s) Proband: %s locus %d, #ofChildren: %d\n",
	  pNucFam->nuclearFamilyIndex, pParents[DAD]->sID, pParents[MOM]->sID,
	  pProband->sID, locus, pNucFam->numChildren);
  while (pGenotype[DAD])
    {
      pGenotype[MOM] = pParents[MOM]->ppGenotypeList[locus];
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
	      pPair =
		(ParentalPair *) MALLOC ("pPair", sizeof (ParentalPair));
	      memset (pPair, 0, sizeof (ParentalPair));
	      pPair->pppChildGenoList =
		(Genotype ***) MALLOC ("pPair->pppChildGenoList",
				       sizeof (Genotype **) * numChildren);
	      pPair->pChildGenoLen =
		(int *) MALLOC ("pPair->pChildGenoLen",
				sizeof (int) * numChildren);
	      pPair->pGenotype[DAD] = pGenotype[DAD];
	      pPair->pGenotype[MOM] = pGenotype[MOM];
	      /* add the pair to the list under this nuclear family */
	      pPair->pNext = pNucFam->ppParentalPair[locus];
	      pNucFam->ppParentalPair[locus] = pPair;
	      /* increase the count of the parental pairs */
	      pNucFam->pNumParentalPair[locus]++;
	      /* copy all children's valid genotype list into the parental pair's
	       * children genotype array */
	      logMsg (LOGPARENTALPAIR, LOGDEBUG,
		      "Parental Pair(%d) %s(%d, %d) && %s(%d, %d)\n",
		      pNucFam->pNumParentalPair[locus],
		      pNucFam->pParents[DAD]->sID,
		      pGenotype[DAD]->allele[DAD],
		      pGenotype[DAD]->allele[MOM],
		      pNucFam->pParents[MOM]->sID,
		      pGenotype[MOM]->allele[DAD],
		      pGenotype[MOM]->allele[MOM]);
	      for (i = 0; i < numChildren; i++)
		{
		  pChild = pNucFam->ppChildrenList[i];
		  genoLen = pChild->pShadowGenotypeListLen[locus];
		  pPair->pppChildGenoList[i] =
		    (Genotype **) MALLOC ("pPair->pppChildGenoList[i]",
					  sizeof (Genotype *) * genoLen);
		  pPair->pChildGenoLen[i] = genoLen;
		  pChildGeno = pChild->ppShadowGenotypeList[locus];
		  logMsg (LOGPARENTALPAIR, LOGDEBUG,
			  "Child %s #ofGenotypes: %d\n", pChild->sID,
			  genoLen);
		  for (j = 0; j < genoLen; j++)
		    {
		      pPair->pppChildGenoList[i][j] = pChildGeno;
		      logMsg (LOGPARENTALPAIR, LOGDEBUG, "  (%d, %d)\n",
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

int
shadow_genotype_elimnation (Genotype * pGenotype1, Genotype * pGenotype2,
			    NuclearFamily * pNucFam, int locus)
{
  Person *pChild;
  int numChildren = pNucFam->numChildren;
  int i;
  Genotype *pGenotype;

  /* check each child's genotype list against the parental pair */
  for (i = 0; i < numChildren; i++)
    {
      pChild = pNucFam->ppChildrenList[i];
      pGenotype = pChild->ppGenotypeList[locus];
      pChild->ppShadowGenotypeList[locus] = NULL;
      pChild->pShadowGenotypeListLen[locus] = 0;
      while (pGenotype)
	{
	  if (is_parent_child_genotype_compatible
	      (locus, DAD, pGenotype1, pGenotype)
	      && is_parent_child_genotype_compatible (locus, MOM, pGenotype2,
						      pGenotype))
	    {
	      pGenotype->pShadowNext = pChild->ppShadowGenotypeList[locus];
	      pChild->ppShadowGenotypeList[locus] = pGenotype;
	      pChild->pShadowGenotypeListLen[locus] += 1;
	    }
	  pGenotype = pGenotype->pNext;
	}			/* end of looping one child's genotype list */
      /* if we didn't find any compatible genotype for the current child, 
       * then return failure */
      if (pChild->ppShadowGenotypeList[locus] == NULL)
	return -1;
    }				/* done looping all children */

  /* we are here because every child has at least one compatible genotype */


  return 0;
}


void *
MALLOC (char *description, size_t size)
{
  void *ptr = NULL;
  ptr = malloc (size);
  if (ptr == NULL)
    {
      fprintf (stderr, "Memory allocation failed.\n");
    }
  KASSERT (ptr != NULL, "Can't allocate memory for %s of size %d.\n",
	   description, size);
  logMsg (LOGMEMORY, LOGDEBUG, "Allocate %p of size %d - %s \n",
	  ptr, size, description);
  return ptr;
}

void
FREE (char *description, void *ptr)
{
  logMsg (LOGMEMORY, LOGDEBUG, "Free %p - %s \n", ptr, description);
  free (ptr);
}

void *
REALLOC (char *description, void *ptr, size_t size)
{
  void *oldPtr = ptr;
  void *newPtr = NULL;

  newPtr = realloc (ptr, size);
  if (newPtr == NULL)
    {
      fprintf (stderr, "Memory reallocation failed.\n");
    }
  KASSERT (newPtr != NULL, "Can't re-allocate memory for %s of size %d.\n",
	   description, size);
  logMsg (LOGMEMORY, LOGDEBUG, "Reallocate %p(from %p) of size %d - %s \n",
	  newPtr, oldPtr, size, description);
  return newPtr;
}
