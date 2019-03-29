/* This file contains functions to conduct genotype elimination based on
 * Mendel's law - a child's paternal allele has to come from father
 * and maternal allele has to come from mother. That means either the 
 * paternal or maternal allele of the father has to match the child's 
 * paternal allele, similarly for the mother. 
 * */
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
#include "genotype_elimination.h"
#include "likelihood.h"

/* some internal functions */
int nuclear_family_genotype_elimination (int locus, NuclearFamily * pNucFam);
int parent_children_genotype_elimination (int locus, NuclearFamily * pNucFam,
					  int parent);
int parent_parent_children_genotype_elimination (int locus,
						 NuclearFamily * pNucFam,
						 int parent1, int parent2);
int child_parents_genotype_elimination (int locus, NuclearFamily * pNucFam,
					int child);
int child_children_parents_genotype_elimination (int locus,
						 NuclearFamily * pNucFam,
						 int child);


/* This function go through the nuclear families in the pedigree
 * to eliminate any genotype that is against Mendel's law 
 * ex. if a child is 3/4, then neither parents can be 1/2
 * as a parent must contribute to one allele of the child
 * so the parent must have at least 3 or 4 in the genotype 
 * constraints go from parent to child and vice versa 
 * genotype elimination will go through the pedigree in 
 * as many times as needed so there is no more elimination to
 * anyone's genotype */
int
pedigree_genotype_elimination (int locus, Pedigree * pPedigree)
{
  NuclearFamily *pNucFam;
  int i;
  int doneFlag;
  int ret = 0;
  int passNo;

  doneFlag = FALSE;
  passNo = 1;
  while (doneFlag != TRUE)
    {
      ret = 0;
      /* go through the entire list of nuclear family */
      for (i = 0; i < pPedigree->numNuclearFamily; i++)
	{
	  pNucFam = pPedigree->ppNuclearFamilyList[i];
	  /* return will tell us whether there was any genotype got eliminated */
	  ret += nuclear_family_genotype_elimination (locus, pNucFam);
	}
      /* there was no change in the genotype list in any person */
      if (ret == 0)
	doneFlag = TRUE;
      /* for debug purpose, print out pedigree genotype list */
      KLOG (LOGGENOELIM, LOGDEBUG, "Pass %d:\n", passNo);
      //print_pedigree_locus_genotype_list(pPedigree, locus);
      passNo++;
    }


  return 0;
}

/* genotype elimination is done through the unit of nuclear family */
int
nuclear_family_genotype_elimination (int locus, NuclearFamily * pNucFam)
{
  int ret = 0;
  int child;

  /* check each parent against children first */
  ret += parent_children_genotype_elimination (locus, pNucFam, DAD);
  ret += parent_children_genotype_elimination (locus, pNucFam, MOM);
  /* check each child against both parents */
  for (child = 0; child < pNucFam->numChildren; child++)
    {
      ret += child_parents_genotype_elimination (locus, pNucFam, child);
    }
  /* check one parent against the other parent and all children */
  ret +=
    parent_parent_children_genotype_elimination (locus, pNucFam, DAD, MOM);
  ret +=
    parent_parent_children_genotype_elimination (locus, pNucFam, MOM, DAD);

  /* check each child against the other children and both parents */
  for (child = 0; child < pNucFam->numChildren; child++)
    {
      ret +=
	child_children_parents_genotype_elimination (locus, pNucFam, child);
    }

  return ret;
}

/* check this parent's genotype against available children's genotypes 
 * as long as there is at least one genotype of each child compatible
 * with parent's genotype, it will be kept, otherwise eliminated
 * */
int
parent_children_genotype_elimination (int locus, NuclearFamily * pNucFam,
				      int parent)
{
  Person *pParent;
  int i;
  Person *pChild;
  Genotype *pGenotype;
  Genotype *pChildGenotype;
  Genotype *pElimGenotype;
  int keepFlag;
  int elimFlag;
  int ret = 0;
  Locus *pLocus = originalLocusList.ppLocusList[locus];

  /* it could be mom or dad */
  pParent = pNucFam->pParents[parent];

  /* go through all genotypes in this parent's genotype list */
  pGenotype = pParent->ppGenotypeList[locus];
  while (pGenotype != NULL)
    {
      /* keep flag is set to true initially and for each child, as long as we 
       * find one compatible genotype in each child's genotype list
       * we can't eliminate the genotype (yet) 
       * as soon as we find one incompatible child, then we will need to 
       * eliminate this parent's genotype */
      keepFlag = TRUE;
      /* go through each child */
      i = 0;
      for (i = 0; (i < pNucFam->numChildren) && keepFlag == TRUE; i++)
	{
	  pChild = pNucFam->ppChildrenList[i];
	  /* go through this child's genotype list */
	  pChildGenotype = pChild->ppGenotypeList[locus];
	  elimFlag = TRUE;
	  while (pChildGenotype != NULL && elimFlag == TRUE)
	    {
	      /* let's see whether the according allele in this child, i.e., 
	       * paternal or materal allele depends which parent we are working on,
	       * matches either of this parent's allele */
	      if (is_parent_child_genotype_compatible
		  (locus, parent, pChild->sex, pGenotype, pChildGenotype) == TRUE)
		/* found a compatible one, move on to next child */
		elimFlag = FALSE;
	      /* move on to next genotype */
	      pChildGenotype = pChildGenotype->pNext;
	    }			/* end of child's genotype loop */
	  /* if the elimination flag is still set to TRUE after going through the entire list
	   * of child's genotype, then we haven't found any
	   * compatible genotype in this child. Eliminate parent's genotype  */
	  if (elimFlag == TRUE)
	    keepFlag = FALSE;
	}			/* end of children loop */
      pElimGenotype = pGenotype;
      /* move to the next genotype for this parent */
      pGenotype = pGenotype->pNext;
      /* if keepFlag is not true, then this parent's genotype is not compatible */
      if (keepFlag != TRUE)
	{
	  remove_genotype (&pParent->ppGenotypeList[locus], pElimGenotype,
			   &pParent->pNumGenotype[locus]);
	  /* increase the counter of elimination */
	  ret++;
	}
    }				/* end of looping of parent genotypes */

  KASSERT (pParent->ppGenotypeList[locus] != NULL,
	   "Pedigree %s Person %s is not compatible at locus %s.\n",
	   pParent->pPedigree->sPedigreeID, pParent->sID, pLocus->sName);

  return ret;
}

/* check a child's genotype against both parents' 
 * if no genotype has been eliminated for this child, return 0
 * otherwise return number of genotypes that have been eliminated 
 *
 * as long as the child's genotype is compatible with at least one 
 * genotype for both parents, there will be no elimination of this genotype
 *
 * */
int
child_parents_genotype_elimination (int locus, NuclearFamily * pNucFam,
				    int child)
{
  int ret = 0;
  Person *pChild = pNucFam->ppChildrenList[child];
  Genotype *pChildGenotype;
  Genotype *pParent1Genotype;
  Genotype *pParent2Genotype;
  Genotype *pElimGenotype;
  Person *pParent1, *pParent2;
  int parent1, parent2;
  int elimFlag;
  int keepFlag;
  int compatibleFlag;
  Locus *pLocus = originalLocusList.ppLocusList[locus];

  /* loop through all genotypes for this child */
  pChildGenotype = pChild->ppGenotypeList[locus];
  while (pChildGenotype != NULL)
    {
      keepFlag = TRUE;
      /* now find a compatible set from parents */
      elimFlag = FALSE;
      parent1 = DAD;
      pParent1 = pNucFam->pParents[parent1];
      parent2 = MOM;
      pParent2 = pNucFam->pParents[parent2];
      pParent1Genotype = pParent1->ppGenotypeList[locus];
      compatibleFlag = FALSE;
      while (pParent1Genotype != NULL && compatibleFlag == FALSE &&
	     elimFlag == FALSE)
	{
	  /* for each parent, go through all the genotypes and make sure 
	   * the child's according allele appears in at least one of the
	   * genotype for this parent */
	  pParent2Genotype = pParent2->ppGenotypeList[locus];
	  while (pParent2Genotype != NULL && compatibleFlag == FALSE &&
		 elimFlag == FALSE)
	    {
	      if (is_parent_child_genotype_compatible (locus, parent1, pChild->sex,
						       pParent1Genotype,
						       pChildGenotype) == TRUE
		  && is_parent_child_genotype_compatible (locus, parent2, pChild->sex,
							  pParent2Genotype,
							  pChildGenotype) ==
		  TRUE)
		{
		  /* found compatibility with both parents */
		  compatibleFlag = TRUE;
		  /* we have exhausted all children, if elimination flag is still true
		   * then need to move on next parental pair, otherwise we
		   * all children are compatible with the current parental pair */
		}		/* end of if compatible parental pair */
	      /* move on to this parent's next genotype */
	      pParent2Genotype = pParent2Genotype->pNext;
	    }			/* end of parent2's genotype list looping */
	  pParent1Genotype = pParent1Genotype->pNext;
	}			/* end of parent1's genotype list looping */
      /* we have exhausted all parental pair 
       * if compatibleFlag is still FALSE, we need to eliminate this child's 
       * genotype */
      pElimGenotype = pChildGenotype;
      /* move on to this child's next genotype */
      pChildGenotype = pChildGenotype->pNext;
      if (compatibleFlag == FALSE)
	{
	  remove_genotype (&pChild->ppGenotypeList[locus], pElimGenotype,
			   &pChild->pNumGenotype[locus]);
	  /* increase the counter of elimination */
	  ret++;
	}
    }				/* end of looping of child's genotype list */

  KASSERT (pChild->ppGenotypeList[locus] != NULL,
	   "Pedigree %s Person %s is not compatible at locus %s.\n", 
	   pChild->pPedigree->sPedigreeID, 
	   pChild->sID, 
	   pLocus->sName);
  return ret;
}

/* check a parent's genotype using information in the other parent and also
 * children. if together with at least one genotype in the other parent is
 * compatible with at least one genotype in each child, then no need for
 * elimination */
int
parent_parent_children_genotype_elimination (int locus,
					     NuclearFamily * pNucFam,
					     int parent1, int parent2)
{
  Genotype *pParent1Geno;
  Genotype *pParent2Geno;
  Genotype *pChildGeno;
  Genotype *pElimGeno;
  int ret = 0;			/* number of elimination has been performed */
  int keepFlag;
  //int elimFlag;
  int child;
  int parent2CompatibleFlag;
  int compatibleFlag = TRUE;
  int skipFlag;
  Person *pParent1;
  Person *pChild;
  Locus *pLocus = originalLocusList.ppLocusList[locus];

  pParent1 = pNucFam->pParents[parent1];
  pParent1Geno = pParent1->ppGenotypeList[locus];

  /* check each genotype of parent1's for elimination */
  while (pParent1Geno != NULL)
    {
      /* we assume we keep the genotype */
      keepFlag = TRUE;
      pParent2Geno = pNucFam->pParents[parent2]->ppGenotypeList[locus];
      /* go through parent2's genotype list to construct parent pair
       * genotypes */
      parent2CompatibleFlag = FALSE;
      while (pParent2Geno != NULL && keepFlag == TRUE &&
	     parent2CompatibleFlag == FALSE)
	{
	  /* check the current parent genotype pair, see whether we can 
	   * find at least one compatible genotype in each child */
	  skipFlag = FALSE;
	  for (child = 0; child < pNucFam->numChildren && skipFlag == FALSE;
	       child++)
	    {
	      pChild = pNucFam->ppChildrenList[child];
	      pChildGeno =
		pNucFam->ppChildrenList[child]->ppGenotypeList[locus];
	      compatibleFlag = FALSE;
	      while (pChildGeno != NULL && compatibleFlag == FALSE)
		{
		  if (is_parent_child_genotype_compatible
		      (locus, parent1, pChild->sex, pParent1Geno, pChildGeno) == TRUE
		      && is_parent_child_genotype_compatible (locus, parent2, pChild->sex,
							      pParent2Geno,
							      pChildGeno) ==
		      TRUE)
		    {
		      /* current parent pair genotypes compatible with each child on 
		       * at least one genotype, move on */
		      compatibleFlag = TRUE;
		    }

		  pChildGeno = pChildGeno->pNext;
		}
	      /* if compatibleFlag is FALSE, then no genotype for current child
	       * is compatible with the parents' genotype pair, move on next 
	       * parent pair, skipping the rest of the children */
	      if (compatibleFlag == FALSE)
		{
		  skipFlag = TRUE;
		}
	    }			/* end of looping each child */
	  if (compatibleFlag == TRUE)
	    {
	      /* no elimination need, skip the rest of the other parent's genotype */
	      parent2CompatibleFlag = TRUE;
	    }
	  pParent2Geno = pParent2Geno->pNext;
	}			/* end of looping spouse genotype list */
      pElimGeno = pParent1Geno;
      pParent1Geno = pParent1Geno->pNext;
      /* if parent2Compatible flag is still set to FALSE, then eliminate parent1's
         current gentoype */
      if (parent2CompatibleFlag == FALSE)
	{
	  remove_genotype (&pNucFam->pParents[parent1]->ppGenotypeList[locus],
			   pElimGeno,
			   &pNucFam->pParents[parent1]->pNumGenotype[locus]);
	  ret++;
	}

    }				/* end of looping parent's genotype list */
  KASSERT (pParent1->ppGenotypeList[locus] != NULL,
	   "Pedigree %s Person %s is not compatible at locus %s.\n", 
	   pParent1->pPedigree->sPedigreeID, pParent1->sID, pLocus->sName);

  return ret;
}

/* check a child's genotype against all other children and both parents' 
 * if no genotype has been eliminated for this child, return 0
 * otherwise return number of genotypes that have been eliminated 
 *
 * as long as the child's genotype is compatible with at least one 
 * set of genotypes for both parents and the rest of the children, 
 * there will be no elimination of this genotype
 *
 * */
int
child_children_parents_genotype_elimination (int locus,
					     NuclearFamily * pNucFam,
					     int child)
{
  int ret = 0;
  Person *pChild = pNucFam->ppChildrenList[child];
  Person *pOtherChild;
  Genotype *pChildGenotype;
  Genotype *pOtherChildGenotype;
  Genotype *pParent1Genotype;
  Genotype *pParent2Genotype;
  Genotype *pElimGenotype;
  Person *pParent1, *pParent2;
  int parent1, parent2;
  int elimFlag;
  int keepFlag;
  int compatibleFlag;
  int skipFlag;
  int i;
  Locus *pLocus = originalLocusList.ppLocusList[locus];

  /* loop through all genotypes for this child */
  pChildGenotype = pChild->ppGenotypeList[locus];
  while (pChildGenotype != NULL)
    {
      keepFlag = TRUE;
      /* now find a compatible set from parents */
      elimFlag = FALSE;
      parent1 = DAD;
      pParent1 = pNucFam->pParents[parent1];
      parent2 = MOM;
      pParent2 = pNucFam->pParents[parent2];
      pParent1Genotype = pParent1->ppGenotypeList[locus];
      compatibleFlag = FALSE;
      while (pParent1Genotype != NULL && compatibleFlag == FALSE &&
	     elimFlag == FALSE)
	{
	  /* for each parent, go through all the genotypes and make sure 
	   * the child's according allele appears in at least one of the
	   * genotype for this parent */
	  pParent2Genotype = pParent2->ppGenotypeList[locus];
	  while (pParent2Genotype != NULL && compatibleFlag == FALSE &&
		 elimFlag == FALSE)
	    {
	      if (is_parent_child_genotype_compatible (locus, parent1, pChild->sex,
						       pParent1Genotype,
						       pChildGenotype) == TRUE
		  && is_parent_child_genotype_compatible (locus, parent2, pChild->sex,
							  pParent2Genotype,
							  pChildGenotype) ==
		  TRUE)
		{
		  /* found compatibility with both parents 
		   * now go through each child 
		   * if the current parent pair doesn't match at least one genotype
		   * of any child, then set elimination flag to TRUE */
		  skipFlag = FALSE;
		  for (i = 0; (i < pNucFam->numChildren) && skipFlag == FALSE;
		       i++)
		    {
		      if (i == child)
			continue;
		      /* assume we need to eliminate the genotype, 
		       * if we can't find a compatible genotype in each child 
		       * as soon as we found a child with no genotype compatible
		       * with the current parental pair, we should skip the 
		       * rest of the children and move on to next parental pair
		       * */
		      pOtherChild = pNucFam->ppChildrenList[i];
		      elimFlag = TRUE;
		      pOtherChildGenotype =pOtherChild->ppGenotypeList[locus];
		      while (pOtherChildGenotype != NULL && elimFlag == TRUE)
			{
			  if (is_parent_child_genotype_compatible
			      (locus, parent1, pOtherChild->sex, pParent1Genotype,
			       pOtherChildGenotype) == TRUE
			      && is_parent_child_genotype_compatible (locus,
								      parent1, pOtherChild->sex,
								      pParent1Genotype,
								      pOtherChildGenotype)
			      == TRUE)
			    {
			      /* found a compatible genotype in this child, move on 
			       * to next child */
			      elimFlag = FALSE;
			    }
			  pOtherChildGenotype = pOtherChildGenotype->pNext;
			}	/* end of looping other child's genotype list */
		      /* we have exhausted current child's genotype list 
		       * if elimination flag is true, we have no luck finding a 
		       * matching genotype, skip the rest of the children
		       * and move on to next paternal pair */
		      if (elimFlag == TRUE)
			skipFlag = TRUE;
		      else
			/* we have made sure every child is compatible with the
			 * current parental pair, done */
			compatibleFlag = TRUE;
		    }		/* end of looping all children */
		  /* we have exhausted all children, if elimination flag is still true
		   * then need to move on next parental pair, otherwise we
		   * all children are compatible with the current parental pair */
		}		/* end of if compatible parental pair */
	      /* move on to this parent's next genotype */
	      pParent2Genotype = pParent2Genotype->pNext;
	    }			/* end of parent2's genotype list looping */
	  pParent1Genotype = pParent1Genotype->pNext;
	}			/* end of parent1's genotype list looping */
      /* we have exhausted all parental pair 
       * if compatibleFlag is still FALSE, we need to eliminate this child's 
       * genotype */
      pElimGenotype = pChildGenotype;
      /* move on to this child's next genotype */
      pChildGenotype = pChildGenotype->pNext;
      if (elimFlag == TRUE)
	{
	  remove_genotype (&pChild->ppGenotypeList[locus], pElimGenotype,
			   &pChild->pNumGenotype[locus]);
	  /* increase the counter of elimination */
	  ret++;
	}
    }				/* end of looping of child's genotype list */

  KASSERT (pChild->ppGenotypeList[locus] != NULL,
	   "Pedigree %s Person %s is not compatible at locus %s.\n", 
	   pChild->pPedigree->sPedigreeID,
	   pChild->sID, pLocus->sName);

  return ret;
}

/* return TRUE if the parent's genotype is compatible with 
 * the child's. 
 * For example if the parent is DAD, then the paternal allele 
 * of the child's genotype has to match either of DAD's alleles
 *
 * we could use allele number to check when we are only considering 
 * singleton allele, 
 * if we are using super allele, we should really make sure the allele
 * sets for the parents are subsets of the children's
 *
 * well, we should always take the route of checking allele set
 * */
inline int
is_parent_child_genotype_compatible (int locus, int parent, int childSex,
				     Genotype * pParentGeno,
				     Genotype * pChildGeno)
{
  /* number of integers used to represent the bit masks for 
   * the alleles in this locus */
  int alleleSetLen = originalLocusList.alleleSetLen;
  int i;
  unsigned int *pParentAlleleSet1 = pParentGeno->pAlleleBits[DAD];
  unsigned int *pParentAlleleSet2 = pParentGeno->pAlleleBits[MOM];
  unsigned int *pChildAlleleSet = pChildGeno->pAlleleBits[parent];
  int compatible[2];

  /* instead of using allele number, we use allele bit masks */
  /*
     if(pChildGeno->allele[parent] == pParentGeno->allele[DAD] ||
     pChildGeno->allele[parent] == pParentGeno->allele[MOM] )
     return TRUE;
     else
     return FALSE;
   */
  /* exact match */
#if 0
  for (i = 0; i < alleleSetLen; i++)
    {
      if (*pChildAlleleSet != *pParentAlleleSet1 &&
	  *pChildAlleleSet != *pParentAlleleSet2)
	{
	  /* incompatibility found, return FALSE */
	  return FALSE;
	}
    }
#endif

  compatible[DAD] = 1;
  compatible[MOM] = 2;
  if(modelOptions.sexLinked && parent == DAD && childSex+1 == MALE)
    {
      /* male child doesn't inherit X chromosome from dad, so no need to check */
      pChildGeno->inheritance[parent] = compatible[DAD] + compatible[MOM];
      return TRUE;
    }

#if 1
  /* make sure the parent's allele set is a subset of the child's 
   * if the child's allele is not compatible with either maternal 
   * and paternal alleles, then return error */
  for (i = 0; i < alleleSetLen; i++)
    {
      if ((*pChildAlleleSet & *pParentAlleleSet1) != *pParentAlleleSet1)
	{
	  /* this child's genotype is not compatible with parent's paternal allele */
	  compatible[DAD] = 0;
	}
      if ((*pChildAlleleSet & *pParentAlleleSet2) != *pParentAlleleSet2)
	{
	  /* this child's genotype is not compatible with parent's maternal allele */
	  compatible[MOM] = 0;
	}
      if (compatible[DAD] == 0 && compatible[MOM] == 0)
	{
	  /* incompatibility found, return FALSE */
	  return FALSE;
	}
    }
#endif

  pChildGeno->inheritance[parent] = compatible[DAD] + compatible[MOM];

  /* no mismatch found, return TRUE */
  return TRUE;
}

inline int
is_parent_child_allele_compatible (int alleleSetLen,
				   int *pParentAlleleSet,
				   int parent, int childSex, Genotype * pChildGenotype)
{
  int i;
  unsigned int *pChildAlleleSet = pChildGenotype->pAlleleBits[parent];

  if(modelOptions.sexLinked && parent == DAD && childSex+1 == MALE)
    {
      /* male child doesn't inherit X chromosome from dad, so no need to check */
      return TRUE;
    }
  for (i = 0; i < alleleSetLen; i++)
    {
      if ((*pChildAlleleSet & *pParentAlleleSet) != *pParentAlleleSet)
	{
	  return FALSE;
	}
    }

  return TRUE;
}
