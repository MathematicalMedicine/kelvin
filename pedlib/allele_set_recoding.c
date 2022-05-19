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
#include "../utils/polynomial.h"

int initialize_transmitted_alleles (int locus, Pedigree * pPedigree);
int identify_transmitted_alleles (int locus, Pedigree * pPedigree);
Person *identify_child_transmitted_alleles (Person * pPerson,
					    Person * pChild, int locus);
int count_alleles (unsigned int *pAlleleSet, int alleleSetLen);

//int set_allele_bit (int allele, unsigned int *pBits);
int recode_genotype (int locus, Pedigree * pPedigree);
int print_pedigree_allele_sets (int locus, Pedigree * pPedigree);
void print_person_allele_set (Person * pPerson, int locus, int alleleSetLen);
int print_allele_set (unsigned int *pAllele, int len);
int print_locus_allele_set (Locus * pLocus, int alleleSetLen);

int add_allele_set (int locus, unsigned int *pAlleleBits);
int find_allele_set (Locus * pLocus, unsigned int *pAlleleBits,
		     int alleleSetLen);
int is_subset (unsigned int *pAlleleBits1, unsigned int *pAlleleBits2,
	       int alleleSetLen);

/* allele set recoding keeps track of what alleles have been transmitted
 * along founder lines and what not. By combining alleles not seen 
 * transmitted in the given pedigree, it will simplify likelihood 
 * calculation - less genotype combinations for sure if there are 
 * more than 2 alleles not seen transmitted for any person in the pedigree
 *
 * allele set recoding only makes sense when the locus has more than
 * two possible alleles, i.e. the maximum we can achieve in this situation
 * by clumping untransmitted alleles together in a super allele set is 1 
 * Note: we DO allele set recoding for SNPs as well, since some unknown individuals
 *       can be recoded as "3|3" in certain situations
 */
int
allele_set_recoding (int locus, Pedigree * pPedigree)
{
  int numPerson;
  int *pDonePerson;

#if 0
  /* it only makes sense to do locus with at least 3 alleles */
  if (originalLocusList.ppLocusList[locus]->numAllele <= 2)
    return 0;
#endif

  if (originalLocusList.ppLocusList[locus]->locusType == LOCUS_TYPE_TRAIT)
    return 0;


  /* get number of people for this pedigree */
  numPerson = pPedigree->numPerson;
  pDonePerson = pPedigree->pPedigreeSet->pDonePerson;
  memset (pDonePerson, 0, numPerson * sizeof (int));

  initialize_transmitted_alleles (locus, pPedigree);

  /* first identify transmitted and non-transmitted alleles */
  identify_transmitted_alleles (locus, pPedigree);

  if (pPedigree->loopFlag) {
    /* since loop breaker duplicates are tracked separately, make sure it got 
     * back to the original */
    identify_transmitted_alleles (locus, pPedigree);
  }

  /* Print out some debug messages if the debug level is turned on */
  DIAG (ALLELE_SET_RECODING, 1, {print_pedigree_allele_sets (locus, pPedigree);});

  /* now we actually recode the genotypes based on the 
   * non-transmitted allele sets */
  recode_genotype (locus, pPedigree);

  return 0;
}

int
initialize_transmitted_alleles (int locus, Pedigree * pPedigree)
{
  Person *pPerson;
  int numPerson = pPedigree->numPerson;
  int i, j, k, m;
  Locus *pLocus;

  /* number of alleles */
  int numAllele;

  /* allele set length for this locus */
  int alleleSetLen;
  Genotype *pGenotype;
  int nonTransmitted[2];
  int alleleCount;
  unsigned long int mask;

  /* get number of people for this pedigree */
  numPerson = pPedigree->numPerson;

  /* get the locus structure */
  pLocus = originalLocusList.ppLocusList[locus];

  /* get number of alleles for this locus */
  numAllele = pLocus->numAllele;

  /* how many integers needed to represent one allele in bit mask 
   * on this locus */
  alleleSetLen = originalLocusList.alleleSetLen;

  /* and initialize the transmitted and non-transmitted allele sets */
  /* loop through each individual to set the initial values of the transmitted 
   * and non-transmitted allele sets 
   * basically for the transmitted allele*/
  for (i = 0; i < numPerson; i++) {
    pPerson = pPedigree->ppPersonList[i];
    /* pass the loop breaker duplicates - the tracking will be done 
     * at the original person 
     */
    if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] == NULL)
      continue;

    /* 
     * There shouldn't be any loops in the pedigree as we should have 
     * broken the loops already 
     */
    for (k = 0; k < alleleSetLen; k++) {
      pPerson->pTransmittedAlleles[MOM][k] = 0;
      pPerson->pTransmittedAlleles[DAD][k] = 0;
      pPerson->pNonTransmittedAlleles[MOM][k] = 0;
      pPerson->pNonTransmittedAlleles[DAD][k] = 0;
    }


    /* go through each genotype in the person's list for this locus
     * add the possible alleles to basic non-transmitted allele list
     * the logic is that we will assume all these alleles are not 
     * actually seen transmitted, until we found the allele transmitted
     * by the line of descent at which time we will unmask the allele bit. */
    nonTransmitted[DAD] = FALSE;
    nonTransmitted[MOM] = FALSE;
    pGenotype = pPerson->ppGenotypeList[locus];
    while (pGenotype) {
      /* set the possible alleles in the non-transmitted list */
      for (k = 0; k < alleleSetLen; k++) {
	pPerson->pNonTransmittedAlleles[DAD][k] |=
	  pGenotype->pAlleleBits[DAD][k];
	pPerson->pNonTransmittedAlleles[MOM][k] |=
	  pGenotype->pAlleleBits[MOM][k];
      }
      /* count number of allele bits set in this genotype 
       * if multiple are set, then it's a set recoded genotype */
      alleleCount = count_alleles (pGenotype->pAlleleBits[DAD], alleleSetLen);
      if (alleleCount > 0)
	nonTransmitted[DAD] = TRUE;
      alleleCount = count_alleles (pGenotype->pAlleleBits[MOM], alleleSetLen);
      if (alleleCount > 0)
	nonTransmitted[MOM] = TRUE;

      pGenotype = pGenotype->pNext;
    }				/* end of loop of genotypes */
    /* if there was a genotype with set recoded alleles, we will set 
     * the possible allele set to contain all alleles for this person 
     * why ???? */
    for (j = DAD; j <= MOM; j++) {
      if (nonTransmitted[j] != TRUE)
	continue;
      for (m = 1; m <= numAllele; m++) {
	mask = 1 << ((m - 1) % INT_BITS);
	pPerson->pNonTransmittedAlleles[j][(m - 1) / INT_BITS] |= mask;
      }
    }
  }				/* end of loop of person */

  return 0;
}

/* This function pulls transmitted and non-transmitted alleles together 
 * for each untyped individual through the line of descent
 * This is a recursive process from parent to children until we exhaust 
 * the pedigree or until we hit a typed individual, through which 
 * transmitted alleles will be reflected to the parents
 * Maternal and paternal alleles are kept track of separately 
 * */
int
identify_transmitted_alleles (int locus, Pedigree * pPedigree)
{

  Person *pPerson, *pChild;
  int numPerson = pPedigree->numPerson;
  int i;

  /* allele set length for this locus */
  int alleleSetLen;
  int alleleCount;
  int *pDonePerson;

  /* how many integers needed to represent one allele in bit mask 
   * on this locus */
  alleleSetLen = originalLocusList.alleleSetLen;

  /* initialize to mark every individual as NOT DONE YET */
  pDonePerson = pPedigree->pPedigreeSet->pDonePerson;
  for (i = 0; i < numPerson; i++) {
    pDonePerson[i] = FALSE;
  }

  for (i = 0; i < numPerson; i++) {
    pPerson = pPedigree->ppPersonList[i];
    pChild = pPerson->pFirstChild;
    if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] == NULL) {
      pPerson = pPerson->pOriginalPerson;
    }
    /* only do untyped person */
    if (pPerson->pTypedFlag[locus] == TRUE)
      continue;
    if (pDonePerson[i] == TRUE)
      continue;
    pDonePerson[i] = TRUE;
    /* we will determine what alleles are transmitted by the line of 
     * descent */
    while (pChild != NULL) {
      /* this function will go through this child's siblings and own and 
       * their descendents */
      pChild = identify_child_transmitted_alleles (pPerson, pChild, locus);
    }

    /* by looking at how many untransmitted alleles, we can determine
     * whether set recode for this person is necessary */
    pPerson->recodeFlag = TRUE;
    /* if this person is a loop breaker, we have to do set recoding */
    if (pPerson->loopBreaker != FALSE) {
      if ((alleleCount =
	   count_alleles (pPerson->pNonTransmittedAlleles[DAD],
			  alleleSetLen)) <= 1)
	pPerson->recodeFlag = FALSE;
      if ((alleleCount =
	   count_alleles (pPerson->pNonTransmittedAlleles[MOM],
			  alleleSetLen)) <= 1)
	pPerson->recodeFlag = FALSE;
    }

  }

  return 0;
}

Person *
identify_child_transmitted_alleles (Person * pFounder,
				    Person * pChild, int locus)
{
  /* number of alleles */
  /* allele set length for this locus */
  int alleleSetLen;
  Genotype *pGenotype;
  int j, k;
  Person *pChildChild;
  int sex = pFounder->sex;

  /* get allele set length */
  alleleSetLen = originalLocusList.alleleSetLen;
  if (pChild == NULL)
    return NULL;

  /* if we are doing sex linked chromosome, and the parent is a male, and 
   * the child is a male, then we know the dad can't possibly pass a 
   * X chromosome to the son */
  /* if this child is not typed, then continue the line of descent through
   * this child */
  if (pChild->pTypedFlag[locus] != TRUE) {
    pChildChild = pChild->pFirstChild;
    while (pChildChild != NULL) {
      pChildChild = identify_child_transmitted_alleles (pChild,
							pChildChild, locus);
    }
    if (!(modelOptions->sexLinked == TRUE &&
	  pFounder->sex + 1 == MALE && pChild->sex + 1 == MALE)) {
      /* we are done with this child and this child's descendents, 
       * now update parent's */
      for (j = DAD; j <= MOM; j++) {
	for (k = 0; k < alleleSetLen; k++) {
	  pFounder->pNonTransmittedAlleles[j][k] &=
	    pChild->pNonTransmittedAlleles[sex][k];
	  pFounder->pTransmittedAlleles[j][k] |=
	    pChild->pTransmittedAlleles[sex][k];
	}
      }

    }
  } else {
    /* child is typed - pass on the child's alleles to the parents' 
     * transmitted allele list 
     * for unphased person, there are two genotypes 
     * can't be more than 2, right???? */
    pGenotype = pChild->ppGenotypeList[locus];
    while (pGenotype != NULL) {
      /* we don't know whether the paternal or maternal allele of
       * the parent has been transmitted down */
      for (j = DAD; j <= MOM; j++) {
	for (k = 0; k < alleleSetLen; k++) {
	  pFounder->pNonTransmittedAlleles[j][k] &=
	    ~pGenotype->pAlleleBits[sex][k];
	  pFounder->pTransmittedAlleles[j][k] |= pGenotype->pAlleleBits[j][k];
	}
      }
      pGenotype = pGenotype->pNext;
    }				/* end of genotype loops */
  }				/* end of child is typed */


  /* mark this person as done */
  if (pChild->loopBreaker == 0)
    pChild->pPedigree->pPedigreeSet->pDonePerson[pChild->personIndex] = TRUE;

  return pChild->pNextSib[pFounder->sex];

}

/* for each possible genotype in the list for the locus, set the allele bits 
 * allele >= 1 */
int
set_allele_bit (int allele, unsigned int *pBits)
{
  unsigned int mask;

  /* set up a mask representing this allele, which may be >= 32 */
  mask = 1;
  /* find out the right bit */
  mask = 1 << ((allele - 1) % INT_BITS);
  /* set the right integer in the set */
  pBits[(allele - 1) / INT_BITS] |= mask;

  return 0;
}

/* count how many bits are set in a allele set. 
 * If multiple alleles are set, it indicates it's a super allele set
 * usually it is just singleton with one bit set only 
 * Empty set is not expected, so this function should return >=1 
 * remember allele bit mask representation may use more than 1 integer when
 * the allele number is more than 32
 * so an array of integers are needed
 * */
int
count_alleles (unsigned int *pAlleleSet, int alleleSetLen)
{
  int i, j;
  int count = 0;
  unsigned long int mask;

  for (i = 0; i < alleleSetLen; i++) {
    mask = 1;
    for (j = 0; j < INT_BITS; j++) {
      if ((pAlleleSet[i] & mask) != 0)
	count++;
      mask = mask << 1;
    }
  }
  return count;
}

/* print transmitted and non-transmitted alleles for each person 
 * for debug purpose */
int
print_pedigree_allele_sets (int locus, Pedigree * pPedigree)
{
  int i, j, k;
  NuclearFamily *pNucFam;
  Person *pPerson;
  char *sParent;
  int alleleSetLen;

  /* how many integers needed to represent one allele in bit mask 
   * on this locus */
  alleleSetLen = originalLocusList.alleleSetLen;

  /* going through nuclear families seems to be an easier way 
   * to see the allele sets and how they changed along the 
   * line of descent */
  for (i = 0; i < pPedigree->numNuclearFamily; i++) {
    pNucFam = pPedigree->ppNuclearFamilyList[i];
    fprintf (stderr, "Nuclear Family %d: \n", i);

    /* Print out parents first */
    for (k = DAD; k <= MOM; k++) {
      if (k == DAD)
	sParent = DADID;
      else
	sParent = MOMID;
      pPerson = pNucFam->pParents[k];
      fprintf (stderr, "  %s %s\n", sParent, pPerson->sID);
      if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] == NULL)
	pPerson = pPerson->pOriginalPerson;
      print_person_allele_set (pPerson, locus, alleleSetLen);
    }

    /* Print out children next */
    for (j = 0; j < pNucFam->numChildren; j++) {
      pPerson = pNucFam->ppChildrenList[j];
      fprintf (stderr, "    Child %s\n", pPerson->sID);
      print_person_allele_set (pPerson, locus, alleleSetLen);
    }
  }
  return 0;
}

/* This function go through the existing genotype list for untyped 
 * individuals and construct a super allele made of all non-transmitted
 * alleles. This super allele is added to the allele list for 
 * the locus, the allele frequency for this super allele will be 
 * the sum of all the non-transmitted aleleles 
 * The new allele number will be greater than the origianl number of
 * alleles */
int
recode_genotype (int locus, Pedigree * pPedigree)
{
  Person *pPerson;
  int numPerson = pPedigree->numPerson;
  int i, j;
  Locus *pLocus;

  /* number of alleles */
  int numAllele;

  /* allele set length for this locus */
  int alleleSetLen;
  Genotype *pGenotype;
  Genotype *pGenotype2, *pPrevGenotype;
  unsigned int *pAlleleBits;
  int allele;

  /* get the locus structure */
  pLocus = originalLocusList.ppLocusList[locus];

  /* how many integers needed to represent one allele in bit mask 
   * on this locus */
  alleleSetLen = originalLocusList.alleleSetLen;

  /* go through untyped individuals */
  for (i = 0; i < numPerson; i++) {
    pPerson = pPedigree->ppPersonList[i];
    if (pPerson->pTypedFlag[locus] == TRUE)
      continue;
    if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] == NULL)
      continue;

    /* go through both paternal and maternal allele sets */
    for (j = DAD; j <= MOM; j++) {

//        fprintf (stderr, "recode_genotype i=%d j=%d\n", i, j);

      /* see whether there are more than two non-transmitted alleles */
      pAlleleBits = pPerson->pNonTransmittedAlleles[j];
      numAllele = count_alleles (pAlleleBits, alleleSetLen);

//        fprintf (stderr, "numAllele=%d\n", numAllele);

      if (numAllele <= 1)
	continue;
      /* check whether we need to create a new allele, there is a chance
       * we have already created a super allele with exactly the same
       * sets of alleles */

      allele = find_allele_set (pLocus, pAlleleBits, alleleSetLen);
      if (allele == -1) {
	//allele = add_allele_set(pAlleleBits, alleleSetLen);
	allele = add_allele_set (locus, pAlleleBits);
      }
      /* now go through the genotype list of this person to recode
       * using the super allele */
      pGenotype = pPerson->ppGenotypeList[locus];

//        fprintf (stderr, "Locus=%d\n", locus);

      while (pGenotype != NULL) {
	if (is_subset (pGenotype->pAlleleBits[j], pAlleleBits, alleleSetLen)) {
	  /* need to recode this person */
	  pGenotype->allele[j] = allele;
	  memcpy (pGenotype->pAlleleBits[j], pAlleleBits,
		  sizeof (unsigned int) * alleleSetLen);
	}
	pGenotype = pGenotype->pNext;
      }				/* end of genotype loop */

    }				/* end of paternal and maternal loop */
    /* now go through genotype list again to remove duplicates after 
     * both paternal and maternal alleles have been recoded */
    pGenotype = pPerson->ppGenotypeList[locus];
    while (pGenotype != NULL) {
      /* check the rest of the genotypes against this genotype and
       * remove any duplicates of this genotype */
      pGenotype2 = pGenotype->pNext;
      pPrevGenotype = pGenotype;
      while (pGenotype2 != NULL) {
	if (pGenotype2->allele[DAD] == pGenotype->allele[DAD] &&
	    pGenotype2->allele[MOM] == pGenotype->allele[MOM]) {
	  pPrevGenotype->pNext = pGenotype2->pNext;
	  free (pGenotype2->pAlleleBits[DAD]);
	  free (pGenotype2->pAlleleBits[MOM]);
	  free (pGenotype2);
	  pPerson->pNumGenotype[locus]--;
	} else
	  pPrevGenotype = pGenotype2;
	pGenotype2 = pPrevGenotype->pNext;
      }
      pGenotype = pGenotype->pNext;
    }

  }				/* end of person loop */

  /* for debug purpose, print out locus allele set information */
  DIAG (ALLELE_SET_RECODING, 1, {
      fprintf (stderr, "Genotype after set recoding:\n");

      print_locus_allele_set (pLocus, alleleSetLen);

      /* Print out person genotype */
      for (i = 0; i < numPerson; i++) {
	pPerson = pPedigree->ppPersonList[i];
	print_person_locus_genotype_list (pPerson, locus);
      }
    });

  return 0;
}

void
print_person_allele_set (Person * pPerson, int locus, int alleleSetLen)
{
  Genotype *pGenotype;

  if (pPerson->pTypedFlag[locus] == TRUE) {
    fprintf (stderr, "    Typed:  ");
    pGenotype = pPerson->ppGenotypeList[locus];
    while (pGenotype != NULL) {
      fprintf (stderr, "(%d,%d) ", pGenotype->allele[DAD], pGenotype->allele[MOM]);
      pGenotype = pGenotype->pNext;
    }
    fprintf (stderr, "\n");
  } else {
    /* Paternal transmitted alleles */
    fprintf (stderr, "    Paternal Transmitted Alleles: \n");
    print_allele_set (pPerson->pTransmittedAlleles[DAD], alleleSetLen);
    /* Maternal transmitted alleles */
    fprintf (stderr, "    Maternal Transmitted Alleles: \n");
    print_allele_set (pPerson->pTransmittedAlleles[MOM], alleleSetLen);
    /* Paternal non-transmitted alleles */
    fprintf (stderr, "    Paternal Non-Transmitted Alleles: \n");
    print_allele_set (pPerson->pNonTransmittedAlleles[DAD], alleleSetLen);
    /* Maternal non-transmitted alleles */
    fprintf (stderr, "    Maternal Non-Transmitted Alleles: \n");
    print_allele_set (pPerson->pNonTransmittedAlleles[MOM], alleleSetLen);
  }
}

int
print_allele_set (unsigned int *pAlleleSet, int len)
{
  int i, j;
  int allele;
  unsigned long int mask;
  int base = 0;

  fprintf (stderr, "        ");
  for (i = 0; i < len; i++) {
    mask = 1;
    for (j = 0; j < INT_BITS; j++) {
      if ((mask & pAlleleSet[i]) > 0) {
	allele = base + j + 1;
	/* To avoid printing out file name and line number info... */
	fprintf (stderr, "%d ", allele);
      }
      mask = mask << 1;
    }
    base += INT_BITS;
  }
  fprintf (stderr, "\n");

  return 0;
}


int
print_locus_allele_set (Locus * pLocus, int alleleSetLen)
{
  int i;

  fprintf (stderr, "Locus %s\n", pLocus->sName);
  for (i = 0; i < pLocus->numAlleleSet; i++) {
    fprintf (stderr, "  AlleleSet %d(%f): ", i + 1, pLocus->ppAlleleSetList[i]->sumFreq);
    print_allele_set (pLocus->ppAlleleSetList[i]->pAlleleBits, alleleSetLen);
  }
  return 0;
}

/* create the initial singleton allele set list based on
 * the original allelels i.e. no super allele */
int
construct_original_allele_set_list (int locus)
{
  Locus *pLocus;
  int numAlleleSet;
  int maxNumAlleleSet;
  int i;
  AlleleSet *pAlleleSet;
  int alleleSetLen = originalLocusList.alleleSetLen;

  pLocus = originalLocusList.ppLocusList[locus];
  numAlleleSet = pLocus->numOriginalAllele;
  pLocus->numAlleleSet = numAlleleSet;
  maxNumAlleleSet =
    (numAlleleSet / DEF_LOCUS_MALLOC_INCREMENT +
     1) * DEF_LOCUS_MALLOC_INCREMENT;
  pLocus->maxNumAlleleSet = maxNumAlleleSet;
  /* allocate the space for the array of pointers to actual allele set */
  CALCHOKE(pLocus->ppAlleleSetList, (size_t) 1, sizeof (AlleleSet *) * maxNumAlleleSet, AlleleSet **);
  for (i = 0; i < numAlleleSet; i++) {
    /* allocate space for the actual allele set */
    CALCHOKE(pAlleleSet, (size_t) 1, sizeof (AlleleSet), AlleleSet *);
    pLocus->ppAlleleSetList[i] = pAlleleSet;
    pAlleleSet->alleleID = i + 1;
    pAlleleSet->numAllele = 1;
    MALCHOKE(pAlleleSet->pAlleles, sizeof (int), int *);
    pAlleleSet->pAlleles[0] = pAlleleSet->alleleID;
    CALCHOKE(pAlleleSet->pAlleleBits, (size_t) 1, sizeof (unsigned int) * alleleSetLen, unsigned int *);
    set_allele_bit (pAlleleSet->alleleID, pAlleleSet->pAlleleBits);

    pAlleleSet->maxFreq = pLocus->pAlleleFrequency[i];

    if (modelOptions->polynomial == TRUE) {
      char vName[100];

      sprintf (vName, "freq_l%d_i%d", locus, i);

      if (pLocus->locusType == LOCUS_TYPE_MARKER)
	pAlleleSet->sumFreqPolynomial =
	  constantExp (pLocus->pAlleleFrequency[i]);
      else
	pAlleleSet->sumFreqPolynomial =
	  variableExp (&pLocus->pAlleleFrequency[i], NULL, 'D', vName);
      //fprintf(stderr, "variable: %s value: %f \n",pLocus->pAlleleFrequency[i]);
    } else
      pAlleleSet->sumFreq = pAlleleSet->maxFreq;

  }

  return 0;
}

/* add new allele set to the allele set list */
int
add_allele_set (int locus, unsigned int *pAlleleBits)
{
  Locus *pLocus;
  int numAlleleSet;
  int maxNumAlleleSet;
  AlleleSet *pAlleleSet;
  int alleleSetLen = originalLocusList.alleleSetLen;
  double freq;
  int allele;
  unsigned long int mask;
  int i, k, j;

  pLocus = originalLocusList.ppLocusList[locus];
  numAlleleSet = pLocus->numAlleleSet;
  pLocus->numAlleleSet++;
  maxNumAlleleSet = pLocus->maxNumAlleleSet;
  if (numAlleleSet >= maxNumAlleleSet) {
    maxNumAlleleSet += DEF_LOCUS_MALLOC_INCREMENT;

    REALCHOKE(pLocus->ppAlleleSetList, sizeof (AlleleSet *) * maxNumAlleleSet, AlleleSet **);
    pLocus->maxNumAlleleSet = maxNumAlleleSet;

  }
  /* allocate space for the actual allele set */
  CALCHOKE(pAlleleSet, (size_t) 1,  sizeof (AlleleSet), AlleleSet *);
  pLocus->ppAlleleSetList[numAlleleSet] = pAlleleSet;


  if (modelOptions->polynomial == TRUE) {
    pAlleleSet->sumFreqPolynomial = constantExp (0);
  }

  pAlleleSet->alleleID = numAlleleSet + 1;
  CALCHOKE(pAlleleSet->pAlleleBits, (size_t) 1, sizeof (unsigned int) * alleleSetLen, unsigned int *);

  pAlleleSet->numAllele = count_alleles (pAlleleBits, alleleSetLen);
  MALCHOKE(pAlleleSet->pAlleles, sizeof (int) * pAlleleSet->numAllele, int *);
  i = 0;
  /* go through each allele in the set */
  for (k = 0; k < alleleSetLen; k++) {
    pAlleleSet->pAlleleBits[k] = pAlleleBits[k];
    mask = 1;
    for (j = 0; j < INT_BITS; j++) {

//        fprintf (stderr, "K=%d of %d j=%d of %d condition=%d\n", k,
//                 alleleSetLen, j, INT_BITS,
//                 ((mask & pAlleleBits[k]) == mask));


      if ((mask & pAlleleBits[k]) == mask) {
	/* the bit is set */
	allele = k * INT_BITS + j + 1;
	pAlleleSet->pAlleles[i] = allele;
	i++;
	freq = pLocus->pAlleleFrequency[allele - 1];

//            fprintf (stderr,
//                     "freq=%f pAlleleSet->maxFreq=%f allele=%d i=%d\n",
//                     freq, pAlleleSet->maxFreq, allele, i);

	if (pAlleleSet->maxFreq < freq)
	  pAlleleSet->maxFreq = freq;

	if (modelOptions->polynomial == TRUE) {
	  char vName[100];

	  sprintf (vName, "f[%d][%d]", locus, allele - 1);

//                fprintf (stderr, "STarting printing polynomials\n");
//                expPrinting (pAlleleSet->sumFreqPolynomial);
//                fprintf (stderr, "\n");
//                expPrinting (variableExp
//                             (&pLocus->pAlleleFrequency[allele - 1],
//                              vName));
//                fprintf (stderr, "\n");
	  if (pLocus->locusType == LOCUS_TYPE_MARKER)
	    pAlleleSet->sumFreqPolynomial =
	      plusExp (2, 1.0, pAlleleSet->sumFreqPolynomial, 1.0,
		       constantExp (pLocus->pAlleleFrequency[allele - 1]), 1);
	  else
	    pAlleleSet->sumFreqPolynomial =
	      plusExp (2, 1.0, pAlleleSet->sumFreqPolynomial, 1.0,
		       variableExp (&pLocus->
				    pAlleleFrequency[allele - 1],
				    NULL, 'D', vName), 1);

//                fprintf (stderr, "variable: %s value: %f \n", vName,
//                         pLocus->pAlleleFrequency[allele - 1]);

	} else
	  pAlleleSet->sumFreq += freq;


      }
      mask = mask << 1;
    }
  }

  return numAlleleSet + 1;
}

int
find_allele_set (Locus * pLocus, unsigned int *pAlleleBits, int alleleSetLen)
{
  int i;
  int j;
  AlleleSet *pAlleleSet;
  int foundFlag = TRUE;


  for (i = 0; i < pLocus->numAlleleSet; i++) {
    pAlleleSet = pLocus->ppAlleleSetList[i];
    foundFlag = TRUE;
    for (j = 0; j < alleleSetLen; j++) {
      if (pAlleleBits[j] != pAlleleSet->pAlleleBits[j]) {
	/* not a match, move on */
	foundFlag = FALSE;
	continue;

      }
    }
    if (foundFlag == TRUE) {
      return i + 1;
    }
  }
  return -1;
}

/* return TRUE if allele set 1 is a subset of allele set 2 */
int
is_subset (unsigned int *pAlleleBits1, unsigned int *pAlleleBits2,
	   int alleleSetLen)
{
  int i;

  for (i = 0; i < alleleSetLen; i++) {
    if ((pAlleleBits1[i] & pAlleleBits2[i]) != pAlleleBits1[i])
      return FALSE;
  }
  return TRUE;
}
