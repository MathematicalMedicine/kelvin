/**********************************************************************
 * Copyright 2007, Columbus Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

/* This file contains functions to  compute likelihood for all the pedigrees,
 * peeling procedure etc. 
 * */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>

#include "pedlib.h"
#include "locus.h"
#include "utils.h"		/* for logging */
#include "tools.h"
#include "likelihood.h"
#include "genotype_elimination.h"
#include "polynomial.h"

Genotype **pTempGenoVector;

/* transmission probability matrix */
XMission *xmissionMatrix = NULL;

/* temporary likelihood results for related parental pairs */
typedef struct PPairElement
{
  /* index to the conditional likelihood array of the proband's */
  u_int32_t likelihoodIndex;
  /* number of the same likelihood */
  short count;
  /* temporary likelihood result */
  union
  {
    double likelihood;
#ifndef NO_POLYNOMIAL
    /* likelihood polynomial under polynomial mode */
    Polynomial *likelihoodPolynomial;
#endif
  } slot;
} PPairElement;

PPairElement **ppairMatrix = NULL;
int ppairMatrixRowSize;
int ppairMatrixNumLocus;
int *bitMask = NULL;

int peel_graph (NuclearFamily * pNucFam, Person * pProband,
		int peelingDirection);
int compute_nuclear_family_likelihood (NuclearFamily * pNucFam,
				       Person * pProband,
				       int peelingDirection);
int
loop_parental_pair (NuclearFamily * pNucFam,
		    Person * pPerson,
		    int locus, ParentalPairSpace * pHaplo, 
		    int multiLocusIndex[2],
		    void *dWeight[2]);
int
loop_child_multi_locus_genotype (Person * pPChild, Person * pProband,
				 ParentalPairSpace * pHaplo, int child,
				 int locus, int multiLocusIndex,
				 void *pSum, int xmissionIndex[2]);
void
get_haplotype_freq (int locus, int parent, void *freqPtr,
		    ParentalPairSpace * pHaplo);

int
loop_child_proband_genotype (NuclearFamily * pNucFam, Person * pProband,
			     int peelingDirection, int locus,
			     int multiLocusIndex);

void
loop_phases (NuclearFamily * pNucFam, ParentalPairSpace * pHaplo,
	     Person * pProband, int locus, int multiLocusIndex[2],
	     int multiLocusPhase[2], void *dWeight[2]);
int
calculate_likelihood (NuclearFamily * pNucFam, ParentalPairSpace * pHaplo,
		      Person * pProband, int locus, int multiLocusIndex[2],
		      int multiLocusPhase[2], void *dWeight[2]);
		      
inline void clear_ppairMatrix (PPairElement ** ppMatrix);
inline void initialize_proband_tmpLikelihood (Person * pPerson);
void populate_pedigree_loopbreaker_genotype_vector(Pedigree *pPed);
void populate_loopbreaker_genotype_vector(Person *pLoopBreaker, int locus);
int set_next_loopbreaker_genotype_vector(Pedigree *pPed, int initialFlag);


/* before likelihood calculation, pre-allocate space to store conditional likelihoods 
 * numLocus - number of loci we analyze at a time
 */
int
allocate_likelihood_space (PedigreeSet * pPedigreeList, int numLocus)
{
  Pedigree *pPedigree;
  int i;
  int size;

  pTempGenoVector = (Genotype **) malloc(sizeof(Genotype *) * numLocus);

  for (i = 0; i < pPedigreeList->numPedigree; i++)
    {
      pPedigree = pPedigreeList->ppPedigreeSet[i];
      /* allocate conditional likelihood storage for each person/pedigree */
      allocate_multi_locus_genotype_storage (pPedigree, numLocus);
    }

  /* allocate storage for temporarily stored likelihood for similar parental pairs (only with
   * phase differences 
   * either likelihood itself will be stored there or a pointer to likelihood polynomial will be */
  size = pow (2, numLocus);
  ppairMatrixNumLocus = numLocus;
  ppairMatrixRowSize = size * sizeof (PPairElement);
  ppairMatrix = (PPairElement **) calloc (size, sizeof (PPairElement *));
  for (i = 0; i < size; i++)
    {
      ppairMatrix[i] = (PPairElement *) calloc (size, sizeof (PPairElement));
    }
  if (bitMask != NULL)
    free (bitMask);
  bitMask = malloc (sizeof (int) * (numLocus + 1));
  for (i = 0; i <= numLocus; i++)
    {
      bitMask[i] = pow (2, i) - 1;
    }

  return 0;
}

/* free the storage space for conditionals */
void
free_likelihood_space (PedigreeSet * pPedigreeList)
{
  Pedigree *pPedigree;
  int i;

  for (i = 0; i < pPedigreeList->numPedigree; i++)
    {
      pPedigree = pPedigreeList->ppPedigreeSet[i];
      free_multi_locus_genotype_storage (pPedigree);
    }

  /* free storage for temporary likelihood for similar parental pairs */
  for (i = 0; i < pow(2, locusList->numLocus); i++)
    {
      free (ppairMatrix[i]);
    }
  free (ppairMatrix);
  ppairMatrix = NULL;
  free (bitMask);
  bitMask = NULL;
  free(pTempGenoVector);
  pTempGenoVector = NULL;
}

/* the main API to compute likelihood for all pedigrees in a data set */
int
compute_likelihood (PedigreeSet * pPedigreeList)
{
  Pedigree *pPedigree;
  int i;
  int status;			/* return status of function calls */
  double product_likelihood = 1;	/* product of the likelihoods for all the pedigrees */
  double sum_log_likelihood = 0;	/* sum of the log10(likelihood) for all the pedigrees */
  double log10Likelihood;
  int origLocus = locusList->pLocusIndex[1];	/* locus index in the original locus list 
						 * this is used to find out the pedigree counts
						 * mainly for case control analyses 
						 */


  /* initialization */
  sum_log_likelihood = 0;
  product_likelihood = 1;
  pPedigreeList->likelihood = 1;
  pPedigreeList->log10Likelihood = 0;

  /* loop over pedigrees in the data set */
  for (i = 0; i < pPedigreeList->numPedigree; i++)
    {
      pPedigree = pPedigreeList->ppPedigreeSet[i];

#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	{
	  if (pPedigree->likelihoodPolynomial == NULL)
	    {
	      /* only build likelihood polynomial once, if the ptr is not NULL, it means
	       * the polynomial has been constructed */
//            fprintf(stderr,"The polynomial building for this pedigree should be only once\n");
	      /* initialize likelihood space for each pedigree */
	      initialize_multi_locus_genotype (pPedigree);
//                fprintf(stderr,"Start polynomial building\n");
	      /* put a stamp in the polynomial list to mark the beginning of likelihood build
	       * for this pedigree */
	      makePolynomialStamp2 ();
	      status = compute_pedigree_likelihood (pPedigree);

//                expPrinting(pPedigree->likelihoodPolynomial);
//                fprintf(stderr,"\n");
	      pPedigree->likelihoodPolyList = buildPolyList ();
	      polyListSorting (pPedigree->likelihoodPolynomial,
			       pPedigree->likelihoodPolyList);
	      /* clean up polynomials that are not used in the final pedigree likelihood */
	      partialPolynomialClearance2 ();
	    }
	  /* evaluate likelihood */
	  pPedigree->likelihood =
	    evaluatePoly (pPedigree->likelihoodPolynomial,
			  pPedigree->likelihoodPolyList);

	}
      else
	{
	  initialize_multi_locus_genotype (pPedigree);
	  status = compute_pedigree_likelihood (pPedigree);
	}
#else
      initialize_multi_locus_genotype (pPedigree);
      status = compute_pedigree_likelihood (pPedigree);
#endif

      if (pPedigree->likelihood == 0.0)
	{
	  KLOG (LOGLIKELIHOOD, LOGWARNING,
		"Pedigree %s has likelihood of 0 or too small.\n",
		pPedigree->sPedigreeID);
	  fprintf (stderr,
		   "Pedigree %s has likelihood of 0 or too small.\n",
		   pPedigree->sPedigreeID);
	  product_likelihood = 0.0;
	  sum_log_likelihood = -9999.99;
	  break;
	}
      else if (pPedigree->likelihood < 0.0)
	{
	  KASSERT (pPedigree->likelihood >= 0.0,
		   "Pedigree %s with NEGATIVE likelihood - This is CRAZY!!!.\n",
		   pPedigree->sPedigreeID);
	  product_likelihood = 0.0;
	  sum_log_likelihood = -9999.99;
	  break;
	}
      else
	{
	  if (pPedigree->pCount[origLocus] == 1)
	    {
	      product_likelihood *= pPedigree->likelihood;
	      log10Likelihood = log10 (pPedigree->likelihood);
	    }
	  else
	    {
	      product_likelihood *=
		pow (pPedigree->likelihood, pPedigree->pCount[origLocus]);
	      log10Likelihood =
		log10 (pPedigree->likelihood) * pPedigree->pCount[origLocus];
	    }
	  /*
	     if(log10Likelihood <= __DBL_MIN_10_EXP__ + 1)
	     fprintf(stderr, "Pedigree %s has likelihood that's too small.\n",
	     pPedigree->sPedigreeID);
	   */
	  sum_log_likelihood += log10Likelihood;
	}
    }

  pPedigreeList->likelihood = product_likelihood;
  pPedigreeList->log10Likelihood = sum_log_likelihood;
  KLOG (LOGLIKELIHOOD, LOGDEBUG, "Sum of log Likelihood is: %e\n",
	sum_log_likelihood);

  return 0;
}

/* release polynomial for all pedigrees */
void
pedigreeSetPolynomialClearance (PedigreeSet * pPedigreeList)
{

#ifndef NO_POLYNOMIAL
  Pedigree *pPedigree;
  int i;

  if (modelOptions.polynomial == TRUE)
    {
      for (i = 0; i < pPedigreeList->numPedigree; i++)
	{
	  pPedigree = pPedigreeList->ppPedigreeSet[i];
	  if (pPedigree->likelihoodPolynomial != NULL)
	    {
	      pPedigree->likelihoodPolynomial = NULL;
	      free (pPedigree->likelihoodPolyList->pList);
	      free (pPedigree->likelihoodPolyList);
	      pPedigree->likelihoodPolyList = NULL;
	    }
	}
    }
#endif
}

/* compute likelihood for a given pedigree */
int
compute_pedigree_likelihood (Pedigree * pPedigree)
{
  int i;
  NuclearFamily *pNucFam;	/* nuclear families within the pedigree */
  int status;			/* function return status */
  Person *pProband;		/* peeling proband */
  double likelihood;
  double tmpLikelihood;
#ifndef NO_POLYNOMIAL
  Polynomial *pLikelihoodPolynomial = NULL;
  Polynomial *pTmpLikelihoodPolynomial = NULL;
#endif
  int ret = 0;

  fprintf(stderr, "PEDIGREE: %s (%d/%d)\n", 
       pPedigree->sPedigreeID, pPedigree->pedigreeIndex+1,
       pPedigree->pPedigreeSet->numPedigree);
  if(pPedigree->loopFlag == TRUE)
    {
      populate_pedigree_loopbreaker_genotype_vector(pPedigree);
      ret = set_next_loopbreaker_genotype_vector(pPedigree, TRUE);
    }

#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      pLikelihoodPolynomial = constantExp (0);
    }
  else
    likelihood = 0;
#else
  likelihood = 0;
#endif

  while(ret == 0)
    {
      initialize_multi_locus_genotype (pPedigree);
      /* initialize all the nuclear families before peeling starts
       * in many cases, multiple likelihoods are computed for the same family
       * with different parameters, we need to clean up before (or after) 
       * each calculation */
      for (i = 0; i < pPedigree->numNuclearFamily; i++)
	{
	  pNucFam = pPedigree->ppNuclearFamilyList[i];
	  pNucFam->doneFlag = 0;
	}
      
      /* peeling starts from the peeling proband and eventually will come back to 
       * the same proband 
       * this process will obtain the conditional likelihoods for the proband 
       */
      status = peel_graph (pPedigree->pPeelingNuclearFamily,
			   pPedigree->pPeelingProband,
			   pPedigree->peelingDirection);
      
      /* done peeling, need to add up the conditional likelihood for the leading peeling proband */
      pProband = pPedigree->pPeelingProband;
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	{
	  pTmpLikelihoodPolynomial = constantExp (0);
	}
      else
	tmpLikelihood = 0;
#else
      tmpLikelihood = 0;
#endif
      
      /* loop over all conditional likelihoods */
      for (i = 0; i < pProband->numConditionals; i++)
	{
	  /* Get the joint likelihood = marginal * p(G) 
	   * when the proband is a founder, the weight will be the multilocus genotype probabilities
	   * under LE or haplotype frequencies under LD
	   * otherwise the weight should be 1
	   *
	   */
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      /* build likelihood polynomial */
	      pLikelihoodPolynomial =
		plusExp (2,
			 1.0,
			 pLikelihoodPolynomial,
			 1.0,
			 timesExp (2,
				   pProband->pLikelihood[i].lkslot.
				   likelihoodPolynomial, 1,
				   pProband->pLikelihood[i].wtslot.
				   weightPolynomial, 1, 0), 1);
	    }
	  else
	    {
	      tmpLikelihood += pProband->pLikelihood[i].lkslot.likelihood *
		pProband->pLikelihood[i].wtslot.weight;
	      likelihood += pProband->pLikelihood[i].lkslot.likelihood *
		pProband->pLikelihood[i].wtslot.weight;
	    }
#else
	  tmpLikelihood += pProband->pLikelihood[i].lkslot.likelihood *
	    pProband->pLikelihood[i].wtslot.weight;
	  likelihood += pProband->pLikelihood[i].lkslot.likelihood *
	    pProband->pLikelihood[i].wtslot.weight;
#endif
	}				/* end of looping over all conditionals */
      
      if(pPedigree->loopFlag == TRUE)
	{
	  restore_pedigree_genotype_link_from_saved(pPedigree);
	  ret = set_next_loopbreaker_genotype_vector(pPedigree, FALSE);
	  KLOG (LOGLIKELIHOOD, LOGDEBUG, "Log Likelihood for this fixed looped pedigree %s is: %e\n",
		pPedigree->sPedigreeID, log10 (tmpLikelihood));
	}
      else
	{
	  ret = -1;
	}
    }
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    /* save the polynomial to the pedigree structure */
    pPedigree->likelihoodPolynomial = pLikelihoodPolynomial;
  else
    {
      /* save the likelihood in the pedigree structure */
      pPedigree->likelihood = likelihood;
      KLOG (LOGLIKELIHOOD, LOGDEBUG, "log Likelihood for pedigree %s is: %e\n",
	    pPedigree->sPedigreeID, log10 (likelihood));
    }
#else
  pPedigree->likelihood = likelihood;
  KLOG (LOGLIKELIHOOD, LOGDEBUG, "log Likelihood for pedigree %s is: %e\n",
	pPedigree->sPedigreeID, log10 (likelihood));
#endif

  return 0;
}

/* recursive procedure to go through all nuclear families within one pedigree 
 * pNucFam -- input nuclear family, the top layer is the peeling nuclear family which 
 *            is the nuclear family contains the peeling proband 
 * pProband -- connector 
 * peelingDireciton -- UP/DOWN -- currently it is not used at all 
 */
int
peel_graph (NuclearFamily * pNucFam, Person * pProband, int peelingDirection)
{
  NuclearFamilyConnector *pConnector;	/* connector structure which represents connector 
					   from one nuc to the other within a pedigree */
  NuclearFamily *pNucFam2;	/* another nuc family input nuc family is connected to */
  Person *pConnectPerson;	/* connector individual */
  int i;

  /* if this nuclear family has been processed or in the middle of process, skip it */
  if (pNucFam->doneFlag == TRUE)
    return 0;

  /* mark this nuclear family as done to avoid potential endless recurisve calls */
  pNucFam->doneFlag = TRUE;

  /* go up through the connectors if any */
  pConnector = pNucFam->pUpConnectors;
  while (pConnector)
    {
      pNucFam2 = pConnector->pConnectedNuclearFamily;
      pConnectPerson = pConnector->pConnectedPerson;
      if (pConnectPerson == pNucFam2->pParents[DAD] ||
	  pConnectPerson == pNucFam2->pParents[MOM])
	{
	  /* these two families are connected through multiple marraige */
	  peel_graph (pConnector->pConnectedNuclearFamily,
		      pConnector->pConnectedPerson, PEDIGREE_UP);
	}
      else
	{
	  peel_graph (pConnector->pConnectedNuclearFamily,
		      pConnector->pConnectedPerson, PEDIGREE_DOWN);
	}
      pConnector = pConnector->pNextConnector;
    }

  /* go down through the down connectors if any */
  pConnector = pNucFam->pDownConnectors;
  while (pConnector)
    {
      /* peel up to the proband */
      peel_graph (pConnector->pConnectedNuclearFamily,
		  pConnector->pConnectedPerson, PEDIGREE_UP);

      pConnector = pConnector->pNextConnector;
    }

  /* we are done with the up or down linked nuclear families or we are a leave
   * or top. ready for likelihood computation for this nuclear family 
   * save the original genotype list first 
   * during likelihood calculation, we limit the child proband's genotype to one at a time
   * once done, the original list will be copied back
   * loop breaker's duplicate can't be a proband as the duplicate is the one doesn't 
   * have parents, so it can't be a connector and can't be a proband 
   */
  memcpy (&pProband->ppProbandGenotypeList[0],
	  &pProband->ppGenotypeList[0],
	  sizeof (Genotype *) * originalLocusList.numLocus);
  memcpy (&pProband->pProbandNumGenotype[0],
	  &pProband->pNumGenotype[0],
	  sizeof (int) * originalLocusList.numLocus);

  KLOG (LOGPARENTALPAIR, LOGDEBUG, "\t Proband (%s) haplotype: \n",
	pProband->sID);
  if (pProband->ppHaplotype == NULL)
    {
      /* allocate space for storing proband's haplotype if not already done */
      pProband->ppHaplotype = MALLOC ("pProband->ppHaplotype",
				      sizeof (Genotype *) * sizeof (int) *
				      ppairMatrixNumLocus);
    }

  if (pProband->touchedFlag == TRUE)
    initialize_proband_tmpLikelihood (pProband);

  if (pNucFam->pParents[DAD] != pProband
      && pNucFam->pParents[MOM] != pProband)
    {
      /* proband is a child */
      pNucFam->childProbandFlag = TRUE;
    }
  else
    {
      /* proband is a parent */
      pNucFam->childProbandFlag = FALSE;
    }
  if (pNucFam->childProbandFlag == TRUE)
    {
      /* proband is a child
       * construct all possible multilocus genotype for the proband
       * for each of them, calculate the conditional likelihood of the nuclear family
       */
      loop_child_proband_genotype (pNucFam, pProband, peelingDirection, 0, 0);
    }
  else				/* A parent is the proband */
    {
      compute_nuclear_family_likelihood (pNucFam, pProband, peelingDirection);
      /* copy the temporary results back */
      for (i = 0; i < pProband->numConditionals; i++)
	{
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      if (pProband->touchedFlag == FALSE)
		{
		  pProband->pLikelihood[i].lkslot.likelihoodPolynomial =
		    constantExp (1.0);
		}
	      pProband->pLikelihood[i].lkslot.likelihoodPolynomial = 
		timesExp (2, 
			  pProband->pLikelihood[i].lkslot.likelihoodPolynomial, 
			  1, 
			  pProband->pLikelihood[i].tmpslot.tmpLikelihoodPolynomial, 
			  1, 
			  0);	//Dec 24
	    }
	  else			/* PE is not enabled */
	    {
	      if (pProband->touchedFlag == FALSE)
		pProband->pLikelihood[i].lkslot.likelihood = 1;
	      pProband->pLikelihood[i].lkslot.likelihood *=
		pProband->pLikelihood[i].tmpslot.tmpLikelihood;
	    }
#else
	  if (pProband->touchedFlag == FALSE)
	    pProband->pLikelihood[i].lkslot.likelihood = 1;
	  pProband->pLikelihood[i].lkslot.likelihood *=
	    pProband->pLikelihood[i].tmpslot.tmpLikelihood;
#endif
	}
    }

#if DEBUG
      /* debug purpose only */
      for (i = 0; i < pProband->numConditionals; i++)
	{
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      KLOG (LOGLIKELIHOOD, LOGDEBUG,
		    "Proband %s Conditional Likelihood (%d) = %e. Weight = %e\n",
		    pProband->sID, i,
		    evaluateValue(pProband->pLikelihood[i].lkslot.likelihoodPolynomial),
		    evaluateValue(pProband->pLikelihood[i].wtslot.weightPolynomial));
	    }
	  else
	    {
	      KLOG (LOGLIKELIHOOD, LOGDEBUG,
		    "Proband %s Conditional Likelihood (%d) = %e. Weight = %e\n",
		    pProband->sID, i,
		    pProband->pLikelihood[i].lkslot.likelihood,
		    pProband->pLikelihood[i].wtslot.weight);
	    }
#else
	  KLOG (LOGLIKELIHOOD, DEBUG,
		"Proband %s Conditional Likelihood (%d) = %e. Weight = %e \n",
		pProband->sID, i, pProband->pLikelihood[i].lkslot.likelihood,
		pProband->pLikelihood[i].wtslot.weight);
#endif
	}
#endif

  /* mark the proband as been touched - we have done some likelihood calculation on this person */
  pProband->touchedFlag = TRUE;

  /* copy back the genotypes for the proband */
  memcpy (&pProband->ppGenotypeList[0],
	  &pProband->ppProbandGenotypeList[0],
	  sizeof (Genotype *) * originalLocusList.numLocus);
  memcpy (&pProband->pNumGenotype[0],
	  &pProband->pProbandNumGenotype[0],
	  sizeof (int) * originalLocusList.numLocus);


  return 0;
}

/* This function loops over child proband's multilocus genotypes recursively and 
 * calls compute_nuclear_family_likelihood for each multilocus genotype 
 * pNucFam - the nuclear family we are working on 
 * locus - current working locus 
 * multiLocusIndex - index to the flattened array of conditional likelihood for the proband
 */
int
loop_child_proband_genotype (NuclearFamily * pNucFam,
			     Person * pProband,
			     int peelingDirection,
			     int locus, int multiLocusIndex)
{
  int origLocus = locusList->pLocusIndex[locus];	/* locus index in the original locus list */
  Genotype *pGenotype;
  int position;			/* genotype position */
  Genotype *pNextGenotype;
  int numGenotype;		/* number of possible genotypes for this proband at this locus */
  int multiLocusIndex2;
  int traitLocus;

  double penetrance = 1;
#ifndef NO_POLYNOMIAL
  Polynomial *penetrancePolynomial = NULL;
#endif

  /* we loop over the genotypes of the proband to condition the 
   * likelihood calculation on it */
  numGenotype = pProband->pSavedNumGenotype[origLocus];
  /* calculate the flattened conditional likelihood array index */
  multiLocusIndex2 = multiLocusIndex * numGenotype;
  pGenotype = pProband->ppProbandGenotypeList[origLocus];
  while (pGenotype != NULL)
    {
      /* record this locus's genotype in the haplotype structure - it's really just phased
       * multilocus genotype (not single chromosome haplotype)
       */
      pProband->ppHaplotype[locus] = pGenotype;
      KLOG (LOGLIKELIHOOD, LOGDEBUG, "\t proband (%s) %d|%d \n",
	    pProband->sID,
	    pGenotype->allele[DAD], pGenotype->allele[MOM]);
      /* temporarilly set the next pointer to NULL so to restrict
       * the genotype on the proband to current genotype only */
      pProband->ppGenotypeList[origLocus] = pGenotype;
      pNextGenotype = pGenotype->pNext;
      pGenotype->pNext = NULL;
      pProband->pNumGenotype[origLocus] = 1;
      position = pGenotype->position;
      /* calculate the flattened conditional likelihood array index */
      multiLocusIndex = multiLocusIndex2 + position;

      if (locus < locusList->numLocus - 1)
	{
	  loop_child_proband_genotype (pNucFam, pProband,
				       peelingDirection, locus + 1,
				       multiLocusIndex);
	}
      else
	{
	  /* we have got the entire multi-locus genotype for the proband */
	  compute_nuclear_family_likelihood (pNucFam, pProband,
					     peelingDirection);
	  /* store the likelihood in the corresponding flattened array */
	  if (pProband->touchedFlag == FALSE)
	    {
	      /* if trait locus exists, we need to retrieve the penetrance */
	      traitLocus = locusList->traitLocusIndex;
	      if (traitLocus >= 0)
		{
#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      penetrancePolynomial =
			pProband->ppHaplotype[traitLocus]->penslot.
			penetrancePolynomial;
		    }
		  else
		    {
		      penetrance =
			pProband->ppHaplotype[traitLocus]->penslot.penetrance;
		    }
#else
		  penetrance =
		    pProband->ppHaplotype[traitLocus]->penslot.penetrance;
#endif
		}
	      else		/* no trait locus */
		{
#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      penetrancePolynomial = constantExp (1);
		    }
		  else
		    {
		      penetrance = 1;
		    }
#else
		  penetrance = 1;
#endif
		}		/* no trait locus */

#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  pProband->pLikelihood[multiLocusIndex].lkslot.likelihoodPolynomial = timesExp (2, penetrancePolynomial, 1, pNucFam->likelihoodPolynomial, 1, 0);	//Dec 24

		  /* for a child, the weight should be 1 */
		  pProband->pLikelihood[multiLocusIndex].wtslot.
		    weightPolynomial = constantExp (1);
		}
	      else
		{
		  /* need to update the penetrance factors */
		  pProband->pLikelihood[multiLocusIndex].lkslot.likelihood =
		    penetrance * pNucFam->likelihood;
		  /* for a child, the weight should be 1 */
		  pProband->pLikelihood[multiLocusIndex].wtslot.weight = 1;
		}
#else
	      /* need to update the penetrance factors */
	      pProband->pLikelihood[multiLocusIndex].lkslot.likelihood =
		penetrance * pNucFam->likelihood;
	      /* for a child, the weight should be 1 */
	      pProband->pLikelihood[multiLocusIndex].wtslot.weight = 1;
#endif
	    }			/* first time updating the likelihood for this phased multilocus genotype */
	  else			/* NOT first time updating the likelihood for this phased multilocus genotypes */
	    {
	      /* no need to consider penetrance anymore */
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  pProband->pLikelihood[multiLocusIndex].lkslot.
		    likelihoodPolynomial =
		    timesExp (2,
			      pProband->pLikelihood[multiLocusIndex].lkslot.
			      likelihoodPolynomial, 1,
			      pNucFam->likelihoodPolynomial, 1, 1);
		  //fprintf(stderr,"Likelihood for this entire multi-locus genotype %f %f\n",
		  //evaluateValue(pNucFam->likelihoodPolynomial),
		  //evaluateValue(pProband->pLikelihood[multiLocusIndex].likelihoodPolynomial));
		  KLOG (LOGLIKELIHOOD, LOGDEBUG,
			"Proband %s Conditional Likelihood (%d) = %e.\n",
			pProband->sID, multiLocusIndex,
			evaluateValue(pProband->pLikelihood[multiLocusIndex].lkslot.likelihoodPolynomial));
		}
	      else
		{
		  pProband->pLikelihood[multiLocusIndex].lkslot.likelihood *=
		    pNucFam->likelihood;
		  /*  fprintf (stderr,
		     "Likelihood for this entire multi-locus genotype %f %f\n",
		     pNucFam->likelihood,
		     pProband->pLikelihood[multiLocusIndex].likelihood);
		   */
		  KLOG (LOGLIKELIHOOD, LOGDEBUG,
			"Proband %s Conditional Likelihood (%d) = %e.\n",
			pProband->sID, multiLocusIndex,
			pProband->pLikelihood[multiLocusIndex].lkslot.likelihood);
		}
#else
	      pProband->pLikelihood[multiLocusIndex].lkslot.likelihood *=
		pNucFam->likelihood;
	      KLOG (LOGLIKELIHOOD, LOGDEBUG,
		    "Proband %s Conditional Likelihood (%d) = %e.\n",
		    pProband->sID, multiLocusIndex,
		    pProband->pLikelihood[multiLocusIndex].lkslot.likelihood);
#endif
	    }
	}			/* end of processing complete child multilocus genotype  */
      /*when we are done, need to restore the genotype link list pointer */
      pGenotype->pNext = pNextGenotype;
      pGenotype = pNextGenotype;
    }				/* loop over all possible genotypes */

  return 0;
}

/* for now, we only use parental pair algorithm 
 * This function goes through all possible parental pairs, for each parental pair, it 
 * calculates likelihood
 */
int
compute_nuclear_family_likelihood (NuclearFamily * pNucFam,
				   Person * pProband, int peelingDirection)
{
  int locus;			/* locus index to construct parental pairs */
  double weight[2] = { 1, 1 };	/* weight for the two parents */
#ifndef NO_POLYNOMIAL
  /* need to define some terms for the polynomail operations */
  Polynomial *weightPolynomial[2];
#endif
  int numHaplotypePair = 1;	/* number of multilocus genotypes */
  int numChild;			/* number of children in this nuclear family */
  int i, j;
  int multiLocusIndex[2] = {0, 0};

  /* initialize the weight for each parent */
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      weightPolynomial[0] = constantExp (1);
      weightPolynomial[1] = constantExp (1);
      pNucFam->likelihoodPolynomial = constantExp (0);
    }
#endif
  pNucFam->likelihood = 0;

  /* the following is to help to set the order which parent's genotype to flip first */
  if (pProband == pNucFam->pParents[MOM])
    {
      /* MOM is the proband */
      pNucFam->head = MOM;
      pNucFam->spouse = DAD;
    }
  else
    {
      /* either DAD is the proband or the child is the proband */
      pNucFam->head = DAD;
      pNucFam->spouse = MOM;
    }

  numChild = pNucFam->numChildren;

  /* first construct the parental pair for this nuclear family locus
   * by locus */
  numHaplotypePair = 1;
  for (locus = 0; locus < locusList->numLocus; locus++)
    {
      /* construct parental pair locus by locus */
      construct_parental_pair (pNucFam, pProband, locus);
      /* calculate number of possible multilocus genotypes */
      numHaplotypePair *= parentalPairSpace.pNumParentalPair[locus];
    }

  for (i = DAD; i <= MOM; i++)
    {
      pNucFam->numHetLocus[i] = 0;
      pNucFam->firstHetLocus[i] = -1;
      for (j = 0; j < locusList->numLocus; j++)
	pNucFam->hetFlag[i][j] = 0;
    }

  /* now we can construct haplotypes and get likelihood computed */
  KLOG (LOGPARENTALPAIR, LOGDEBUG, "Haplotype for nuclear family No. %d:\n",
	pNucFam->nuclearFamilyIndex);
  KLOG (LOGPARENTALPAIR, LOGDEBUG, "\t\t\t DAD(%s)\t\t\t MOM(%s)\n",
	pNucFam->pParents[DAD]->sID, pNucFam->pParents[MOM]->sID);
  /* recursively call loop_parental_pair to get a complete multlocus genotype */
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      loop_parental_pair (pNucFam, pProband, 0, &parentalPairSpace,
			  multiLocusIndex, (void *) weightPolynomial);
//      fprintf(stderr,"Conditional likelihood for nuclear family %d is: %e\n",
//            pNucFam->nuclearFamilyIndex, evaluateValue(sumPolynomial));
    }
  else
    {
      loop_parental_pair (pNucFam, pProband, 0, &parentalPairSpace,
			  multiLocusIndex, (void *) weight);
      /*
         KLOG (LOGLIKELIHOOD, LOGDEBUG,
         "Conditional likelihood for nuclear family %d is: %e\n",
         pNucFam->nuclearFamilyIndex, sum);
       */
    }
#else
  loop_parental_pair (pNucFam, pProband, 0, &parentalPairSpace,
		      multiLocusIndex, (void *) weight);
  /*
     KLOG (LOGLIKELIHOOD, LOGDEBUG,
     "Conditional likelihood for nuclear family %d is: %e\n",
     pNucFam->nuclearFamilyIndex, sum);
   */
#endif

  return 0;
}

/* using the parental pair algorithm to go through all possible 
 * parental pairs and calculate the likelihood of each parental pair 
 * nested looping is for the multi-locus */
int
loop_parental_pair (NuclearFamily * pNucFam,
		    Person * pProband,
		    int locus, ParentalPairSpace * pHaplo, 
		    int multiLocusIndex[2], void *dWeight[2])
{
  ParentalPair *pPair;		/* parental pair for current locus */
  int origLocus;		/* locus index in the original locus list */
  int i, j, k;
  int multiLocusIndex2[2];
  Person *pParent[2];		/* two parents */
  int numGenotype[2];		/* number of genotypes at this locus */
  int numPair;			/* index to the list of parental pairs */
  double newWeight[2];		/* genotype weight */
#ifndef NO_POLYNOMIAL
  Polynomial *newWeightPolynomial[2];
#endif
  int head;			/* proband if no child is a proband, otherwise DAD  */
  int spouse;			/* spouse of the head */
  int likelihoodIndex;		/* likelihood index for the proband */
  int multiLocusPhase[2] = { 0, 0 };	/* index to the related parental pair matrix */

  head = pNucFam->head;
  spouse = pNucFam->spouse;

  origLocus = locusList->pLocusIndex[locus];
  for (i = DAD; i <= MOM; i++)
    {
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	newWeightPolynomial[i] = (Polynomial *) dWeight[i];
      else
	newWeight[i] = *((double *) dWeight + i);
#else
      newWeight[i] = *((double *) dWeight + i);
#endif
      pParent[i] = pNucFam->pParents[i];
      /* find the max number of possible genotypes for this parent */
      if(pParent[i]->loopBreaker >= 1 && pParent[i]->pParents[DAD] == NULL)
	numGenotype[i] = pParent[i]->pOriginalPerson->pSavedNumGenotype[origLocus];
      else
	numGenotype[i] = pParent[i]->pSavedNumGenotype[origLocus];
      /* calculate this parent's conditional likelihood's offset */
      multiLocusIndex2[i] = multiLocusIndex[i] * numGenotype[i];
    }
  /* parental pair index for this locus */
  numPair = -1;
  while ((numPair + 1) < pHaplo->pNumParentalPair[locus])
    {
      /* find related parental pairs */
      numPair++;
      /* get the parental pair */
      pPair = &pHaplo->ppParentalPair[locus][numPair];
      multiLocusIndex2[DAD] = multiLocusIndex[DAD] * numGenotype[DAD] + pPair->pGenotype[DAD]->position;
      multiLocusIndex2[MOM] = multiLocusIndex[MOM] * numGenotype[MOM] + pPair->pGenotype[MOM]->position;
      pHaplo->pParentalPairInd[locus] = numPair;
      /* set the het flag */
      for (i = DAD; i <= MOM; i++)
	{
	  if (locus == 0)
	    {
	      /* keep track of how many heterogenous loci there are at and before this locus */
	      pNucFam->tmpNumHet[i][0] = 0;
	    }
	  else
	    {
	      pNucFam->tmpNumHet[i][locus] = pNucFam->tmpNumHet[i][locus - 1];
	    }
	  if (pNucFam->firstHetLocus[i] >= locus)
	    {
	      /* initialize  */
	      pNucFam->firstHetLocus[i] = -1;
	    }
	  if (isHet (pPair->pGenotype[i]))
	    {
	      pNucFam->hetFlag[i][locus] = 1;
	      if (pNucFam->firstHetLocus[i] == -1)
		pNucFam->firstHetLocus[i] = locus;
	      pNucFam->tmpNumHet[i][locus]++;
	    }
	  else
	    pNucFam->hetFlag[i][locus] = 0;
	}
      /* record the start of the related parental pairs for this locus */
      pNucFam->relatedPPairStart[locus] = numPair;
      pNucFam->numRelatedPPair[locus] = 1;
      /* find the related parental pairs that have same pair of genotypes only differ in phases */
      while ((numPair + 1) < pHaplo->pNumParentalPair[locus] &&
	     (pHaplo->ppParentalPair[locus][numPair + 1].phase[DAD] > 0 ||
	      pHaplo->ppParentalPair[locus][numPair + 1].phase[MOM] > 0))
	{
	  numPair++;
	  pNucFam->numRelatedPPair[locus]++;
	}
      if(locus == 0)
	{
	  pNucFam->totalRelatedPPair[locus] = pNucFam->numRelatedPPair[0];
	}
      else
	{
	  pNucFam->totalRelatedPPair[locus] = pNucFam->totalRelatedPPair[locus-1] * 
	    pNucFam->numRelatedPPair[locus];
	}

      for (i = DAD; i <= MOM; i++)
	{
	  if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM
	      && pParent[i]->pParents[DAD] == NULL)
	    {
	      if(pParent[i]->loopBreaker >= 1)
		continue;
	      /* 
	       * For founders:
	       *   under LE, we just multiply the founder weights, 
	       *   under LD, we need to use haplotype freq - which we will do in loop_phases
	       * For non-founders, we don't even care, they should remain as 1 as initialized
	       */
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  newWeightPolynomial[i] =
		    timesExp (2, (Polynomial *) dWeight[i], 1,
			      pPair->pGenotype[i]->wtslot.weightPolynomial, 1,
			      0);
		}
	      else
		{
		  newWeight[i] =
		    *((double *) dWeight + i) * pPair->pGenotype[i]->wtslot.weight;
		}
#else
	      newWeight[i] =
		*((double *) dWeight + i) * pPair->pGenotype[i]->wtslot.weight;
#endif
	    }
	}			/* looping dad and mom genotypes */


      if (locus < locusList->numLocus - 1)
	{
	  /* recursively calling this function to get a complete multilocus genotype */
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {

	      loop_parental_pair (pNucFam, pProband, locus + 1, pHaplo,
				  multiLocusIndex2, (void *) newWeightPolynomial);
	    }
	  else
	    {
	      loop_parental_pair (pNucFam, pProband, locus + 1, pHaplo,
				  multiLocusIndex2, (void *) newWeight);
	    }
#else
	  loop_parental_pair (pNucFam, pProband, locus + 1, pHaplo,
			      multiLocusIndex2, (void *) newWeight);
#endif
	}
      else			/* got a complete set of parental pairs */
	{

	  if(pNucFam->totalRelatedPPair[locus] == 1)
	    {
	      multiLocusPhase[DAD] = 0;
	      multiLocusPhase[MOM] = 0;
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		calculate_likelihood(pNucFam, pHaplo, pProband, locus, 
				     multiLocusIndex2, multiLocusPhase, (void *)newWeightPolynomial);
	      else
		calculate_likelihood(pNucFam, pHaplo, pProband, locus, 
				     multiLocusIndex2, multiLocusPhase, (void *)newWeight);
#else
	      calculate_likelihood(pNucFam, pHaplo, pProband, locus, 
				   multiLocusIndex2, multiLocusPhase, (void *)newWeight);
#endif

	      if(pNucFam->childProbandFlag == TRUE)
		{
#ifndef NO_POLYNOMIAL
		  if(modelOptions.polynomial == TRUE)
		    {
		      pNucFam->likelihoodPolynomial =
			plusExp (2,
				 1.0,
				 pNucFam->likelihoodPolynomial,
				 1.0,
				 ppairMatrix[0][0].slot.likelihoodPolynomial, 1); 
		    }
		  else
		    pNucFam->likelihood += ppairMatrix[0][0].slot.likelihood;
#else
		  pNucFam->likelihood += ppairMatrix[0][0].slot.likelihood;
#endif
		}
	      else /* one of the parent is the proband */
		{
#ifndef NO_POLYNOMIAL
		  if(modelOptions.polynomial == TRUE)
		    {
		      pProband->pLikelihood[multiLocusIndex2[head]].tmpslot.
			tmpLikelihoodPolynomial =
			plusExp (2,
				 1.0,
				 pProband->pLikelihood[multiLocusIndex2[head]].tmpslot.
				 tmpLikelihoodPolynomial, 
				 1.0,
				 ppairMatrix[0][0].slot.likelihoodPolynomial, 1); 
		    }
		  else
		    pProband->pLikelihood[multiLocusIndex2[head]].tmpslot.
		      tmpLikelihood += ppairMatrix[0][0].slot.likelihood;
#else
		  pProband->pLikelihood[multiLocusIndex2[head]].tmpslot.
		    tmpLikelihood += ppairMatrix[0][0].slot.likelihood;
#endif

		}
	    }
	  else if(pNucFam->totalRelatedPPair[locus] > 0)
	    {
	  /* initialize the phase matrix */
	  clear_ppairMatrix (ppairMatrix);

	  /* set some information */
	  pNucFam->numHetLocus[head] = pNucFam->tmpNumHet[head][locus];
	  pNucFam->numHetLocus[spouse] = pNucFam->tmpNumHet[spouse][locus];
	  /* set bit mask - all bits set to 1 - number of bits == number of het loci */
	  pNucFam->hetLocusBits[head] = bitMask[pNucFam->numHetLocus[head]];
	  pNucFam->hetLocusBits[spouse] =
	    bitMask[pNucFam->numHetLocus[spouse]];

	  multiLocusIndex2[DAD] = 0;
	  multiLocusIndex2[MOM] = 0;
	  multiLocusPhase[DAD] = 0;
	  multiLocusPhase[MOM] = 0;
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    loop_phases (pNucFam, pHaplo, pProband, 0, multiLocusIndex2,
			 multiLocusPhase, (void *) newWeightPolynomial);
	  else
	    loop_phases (pNucFam, pHaplo, pProband, 0, multiLocusIndex2,
			 multiLocusPhase, (void *) newWeight);
#else
	  loop_phases (pNucFam, pHaplo, pProband, 0, multiLocusIndex2,
		       multiLocusPhase, (void *) newWeight);
#endif

	  /* post processing of results of similar parental pairs and store them 
	   * in the proband's likelihood space */
	  if (pNucFam->childProbandFlag == TRUE)
	    {
	      /* child is the proband, sum likelihood across all rows and columns  
	       * save it at the tmpLikelihood
	       */
	      for (j = 0; j <= bitMask[pNucFam->numHetLocus[head]]; j++)
		{
		  for (k = 0; k <= bitMask[pNucFam->numHetLocus[spouse]]; k++)
		    {
		      if (ppairMatrix[j][k].count > 1)
			{
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    {
			      pNucFam->likelihoodPolynomial =
				plusExp (2,
					 1.0,
					 pNucFam->likelihoodPolynomial,
					 1.0,
					 timesExp (2,
						   ppairMatrix[j][k].slot.
						   likelihoodPolynomial, 1,
						   constantExp (ppairMatrix[j]
								[k].count), 1,
						   0), 
					 1);
			    }
			  else
			    pNucFam->likelihood +=
			      ppairMatrix[j][k].slot.likelihood *
			      ppairMatrix[j][k].count;
#else
			  pNucFam->likelihood +=
			    ppairMatrix[j][k].slot.likelihood *
			    ppairMatrix[j][k].count;
#endif
			}
		      else if (ppairMatrix[j][k].count > 0)	/* count == 1 */
			{
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    {
			      pNucFam->likelihoodPolynomial =
				plusExp (2,
					 1.0,
					 pNucFam->likelihoodPolynomial,
					 1.0,
					 ppairMatrix[j][k].slot.
					 likelihoodPolynomial, 1);
			    }
			  else
			    pNucFam->likelihood +=
			      ppairMatrix[j][k].slot.likelihood;
#else
			  pNucFam->likelihood +=
			    ppairMatrix[j][k].slot.likelihood;
#endif
			}
		    }
		}
	    }
	  else
	    {			/* one of the parent is the proband */
	      for (j = 0; j <= bitMask[pNucFam->numHetLocus[head]]; j++)
		{
		  /* get the proband's conditional likelihood index for this row
		   * sum up across this row and save them to the proband's likelihood storage
		   */
		  likelihoodIndex = ppairMatrix[j][0].likelihoodIndex;
		  for (k = 0; k <= bitMask[pNucFam->numHetLocus[spouse]]; k++)
		    {
		      if (ppairMatrix[j][k].count > 1)
			{
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    {
			      pProband->pLikelihood[likelihoodIndex].tmpslot.
				tmpLikelihoodPolynomial =
				plusExp (2, 1.0,
					 pProband->
					 pLikelihood[likelihoodIndex].tmpslot.
					 tmpLikelihoodPolynomial, 1.0,
					 timesExp (2,
						   ppairMatrix[j][k].slot.
						   likelihoodPolynomial, 1,
						   constantExp (ppairMatrix[j]
								[k].count), 1,
						   0), 
					 1);
			    }
			  else
			    pProband->pLikelihood[likelihoodIndex].tmpslot.
			      tmpLikelihood +=
			      ppairMatrix[j][k].slot.likelihood *
			      ppairMatrix[j][k].count;
#else
			  pProband->pLikelihood[likelihoodIndex].tmpslot.
			    tmpLikelihood +=
			    ppairMatrix[j][k].slot.likelihood *
			    ppairMatrix[j][k].count;
#endif
			}	/* count > 1 */
		      else if (ppairMatrix[j][k].count > 0)
			{
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    {
			      pProband->pLikelihood[likelihoodIndex].tmpslot.
				tmpLikelihoodPolynomial =
				plusExp (2, 1.0,
					 pProband->
					 pLikelihood[likelihoodIndex].tmpslot.
					 tmpLikelihoodPolynomial, 1.0,
					 ppairMatrix[j][k].slot.
					 likelihoodPolynomial, 1);
			    }
			  else
			    pProband->pLikelihood[likelihoodIndex].tmpslot.
			      tmpLikelihood +=
			      ppairMatrix[j][k].slot.likelihood;
#else
			  pProband->pLikelihood[likelihoodIndex].tmpslot.
			    tmpLikelihood +=
			    ppairMatrix[j][k].slot.likelihood;
#endif
			}	/* count > 0 */
		    }		/* loop through column */
		}		/* loop through row */
	    }			/* one parent is proband */
	    }
	}			/* end of processing one parental pair */
    }

  return 0;
}



/* Compute likelihood for parental multilocus genotypes pairs that are only different in phases 
 * Input: parental pairs for this nuclear family
 *        parental pairs with different phases are next to each other
 *        weight - genotype weights for the parental genotypes 
 *                 LE:
 *                   founder - phase doesn't matter
 *                   non-founder - phase will come to play
 *                 LD: 
 *                   phase matters for both founders and non-founders
 *        penetrance - is not an input as it can be different for different phases if we consider
 *                     imprinting effect
 *                     penetrance only applies to disease locus if it's included in the list 
*/
void
loop_phases (NuclearFamily * pNucFam, ParentalPairSpace * pHaplo,
	     Person * pProband, int locus, int multiLocusIndex[2],
	     int multiLocusPhase[2], void *dWeight[2])
{
  /* conditional likelihood index for each parent */
  int multiLocusIndex2[2];
  /* index to related parental pair matrix */
  int multiLocusPhase2[2];
  /* index of the flip of current pattern in related parental pair matrix */
  int multiLocusPhaseFlip[2];
  Person *pParent[2];
  int numGenotype[2];
  int i;
  int numPair;
  int end;
  int origLocus = locusList->pLocusIndex[locus];	/* locus index in the original locus list */
  ParentalPair *pPair;
  /* whether a fresh calculation is needed or we could find existing pattern result */
  int calculateFlag;
  int phase[2];
  int proband;
  int spouse;
  int likelihoodIndex;
  /*
  double newWeight[2];
  double sum;
  int child;
  int xmissionIndex[2];
  double childProduct;
  double penetrance[2];
  int traitLocus;
  int genoIndex;
#ifndef NO_POlYNOMIAL
  Polynomial *newWeightPolynomial[2];
  Polynomial *penetrancePolynomial[2];
  Polynomial *sumPolynomial;
  Polynomial *childProductPolynomial = NULL;
#endif
  */

  for (i = DAD; i <= MOM; i++)
    {
      pParent[i] = pNucFam->pParents[i];
      /* find the max number of possible genotypes for this parent */
      if(pParent[i]->loopBreaker >= 1 && pParent[i]->pParents[DAD] == NULL)
	numGenotype[i] = pParent[i]->pOriginalPerson->pSavedNumGenotype[origLocus];
      else
	numGenotype[i] = pParent[i]->pSavedNumGenotype[origLocus];
      /* likelihood storage index */
      multiLocusIndex2[i] = multiLocusIndex[i] * numGenotype[i];
      if (pNucFam->hetFlag[i][locus] == 1)
	{
	  /* phase combination index */
	  multiLocusPhase2[i] = multiLocusPhase[i] * 2;
	}
      else
	{
	  multiLocusPhase2[i] = multiLocusPhase[i];
	}
    }

  proband = pNucFam->head;
  spouse = pNucFam->spouse;
  /* get the start of the related pair */
  numPair = pNucFam->relatedPPairStart[locus] - 1;
  end = pNucFam->numRelatedPPair[locus] + pNucFam->relatedPPairStart[locus];
  while ((numPair += 1) < end)
    {
      pPair = &pHaplo->ppParentalPair[locus][numPair];
      pHaplo->pParentalPairInd[locus] = numPair;
      KLOG (LOGLIKELIHOOD, LOGDEBUG,
	    "(%s) %2d->\t %2d|%-2d --X-- %2d|%-2d  (%s)\n",
	    pNucFam->pParents[DAD]->sID,
	    origLocus,
	    pPair->pGenotype[0]->allele[DAD],
	    pPair->pGenotype[0]->allele[MOM],
	    pPair->pGenotype[1]->allele[DAD],
	    pPair->pGenotype[1]->allele[MOM], pNucFam->pParents[MOM]->sID);
      for (i = DAD; i <= MOM; i++)
	{
	  pHaplo->phase[i][locus] = pPair->phase[i];
	  multiLocusPhase[i] = multiLocusPhase2[i] + pPair->phase[i];
	  multiLocusIndex[i] =
	    multiLocusIndex2[i] + pPair->pGenotype[i]->position;
	}
      if (locus < locusList->numLocus - 1)
	{
	  loop_phases (pNucFam, pHaplo, pProband, locus + 1, multiLocusIndex,
		       multiLocusPhase, dWeight);
	}
      else			/* got a complete multilocus parental pair */
	{
	  calculateFlag = 1;
	  ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
	    likelihoodIndex = multiLocusIndex[proband];
	  if (pNucFam->childProbandFlag == TRUE)
	    {
	      phase[DAD] = multiLocusPhase[DAD];
	      phase[MOM] = multiLocusPhase[MOM];
	      /* child is the proband */
	      for (i = DAD; i <= MOM; i++)
		{
		  if (pNucFam->firstHetLocus[i] >= 0 &&
		      pHaplo->phase[i][pNucFam->firstHetLocus[i]] != 0)
		    {
		      /* first het locus has a reverse phase than the origninal 
		       * potentially we could benefit from directly using the likelihood calculated 
		       * with different phase
		       */
		      /* if this parent is a founder or a proband, then flip==original */
		      if (pNucFam->pParents[i]->pParents[DAD] == NULL)
			{
			  multiLocusPhaseFlip[i] =
			    ~multiLocusPhase[i] & pNucFam->hetLocusBits[i];
			  phase[i] = multiLocusPhaseFlip[i];
			  if (ppairMatrix[phase[proband]][phase[spouse]].
			      count > 0)
			    calculateFlag = 0;
			}
		    }
		}
	      if (calculateFlag == 0)
		{
		  ppairMatrix[phase[proband]][phase[spouse]].count++;
#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		    KLOG (LOGLIKELIHOOD, LOGDEBUG,
			  "\t\t likelihood (%d) = %e\n",
			  ppairMatrix[phase[proband]][phase[spouse]].
			  likelihoodIndex,
			  evaluateValue(ppairMatrix[phase[proband]][phase[spouse]].slot.
					likelihoodPolynomial));
		    }
		  else
		    KLOG (LOGLIKELIHOOD, LOGDEBUG,
			  "\t\t likelihood (%d) = %e\n",
			  ppairMatrix[phase[proband]][phase[spouse]].
			  likelihoodIndex,
			  ppairMatrix[phase[proband]][phase[spouse]].slot.
			  likelihood);
#else
		  KLOG (LOGLIKELIHOOD, LOGDEBUG,
			"\t\t likelihood (%d) = %e\n",
			ppairMatrix[phase[proband]][phase[spouse]].
			likelihoodIndex,
			ppairMatrix[phase[proband]][phase[spouse]].slot.
			likelihood);
#endif
		}
	    }			/* a child is the proband */
	  else			/* proband is a parent */
	    {
	      if (pNucFam->firstHetLocus[spouse] >= 0 &&
		  pHaplo->phase[spouse][pNucFam->firstHetLocus[spouse]] != 0)
		{
		  if (pNucFam->pParents[spouse]->pParents[DAD] == NULL)
		    {
		      /* non-proband parent is a founder */
		      /* find the reverse pattern */
		      multiLocusPhaseFlip[spouse] =
			~multiLocusPhase[spouse] & pNucFam->
			hetLocusBits[spouse];
		      if (ppairMatrix[multiLocusPhase[proband]]
			  [multiLocusPhaseFlip[spouse]].count > 0)
			{

			  /* increase the count on the original  */
			  ppairMatrix[multiLocusPhase[proband]]
			    [multiLocusPhaseFlip[spouse]].count++;
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == FALSE)
			    KLOG (LOGLIKELIHOOD, LOGDEBUG,
				  "\t\t likelihood (%d) = %e\n",
				  ppairMatrix[multiLocusPhase[proband]]
				  [multiLocusPhaseFlip[spouse]].
				  likelihoodIndex,
				  ppairMatrix[multiLocusPhase[proband]]
				  [multiLocusPhaseFlip[spouse]].slot.
				  likelihood);
#else
			  KLOG (LOGLIKELIHOOD, LOGDEBUG,
				"\t\t likelihood (%d) = %e\n",
				ppairMatrix[multiLocusPhase[proband]]
				[multiLocusPhaseFlip[spouse]].likelihoodIndex,
				ppairMatrix[multiLocusPhase[proband]]
				[multiLocusPhaseFlip[spouse]].slot.
				likelihood);
#endif
			  calculateFlag = 0;
			}
		    }
		}
	      else if (pNucFam->firstHetLocus[proband] >= 0 &&
		       pHaplo->phase[proband][pNucFam->
					      firstHetLocus[proband]] != 0)
		{
		  /* find the reverse pattern */
		  multiLocusPhaseFlip[proband] =
		    ~multiLocusPhase[proband] & pNucFam->
		    hetLocusBits[proband];
		  if (ppairMatrix[multiLocusPhaseFlip[proband]]
		      [multiLocusPhase[spouse]].count > 0)
		    {
		      /* make sure we have calculated for this pattern before */
		      ppairMatrix[multiLocusPhase[proband]]
			[multiLocusPhase[spouse]].count = 1;
		      likelihoodIndex =
			ppairMatrix[multiLocusPhaseFlip[proband]]
			[multiLocusPhase[spouse]].likelihoodIndex;
#ifndef NO_POLYNOMIAL
		      if (modelOptions.polynomial == TRUE)
			{
			  ppairMatrix[multiLocusPhase[proband]]
			    [multiLocusPhase[spouse]].slot.
			    likelihoodPolynomial =
			    ppairMatrix[multiLocusPhaseFlip[proband]]
			    [multiLocusPhase[spouse]].slot.
			    likelihoodPolynomial;
			  pProband->pLikelihood[multiLocusIndex[proband]].
			    wtslot.weightPolynomial =
			    pProband->pLikelihood[likelihoodIndex].wtslot.
			    weightPolynomial;
			}
		      else
			{
			  ppairMatrix[multiLocusPhase[proband]]
			    [multiLocusPhase[spouse]].slot.likelihood =
			    ppairMatrix[multiLocusPhaseFlip[proband]]
			    [multiLocusPhase[spouse]].slot.likelihood;
			  pProband->pLikelihood[multiLocusIndex[proband]].
			    wtslot.weight =
			    pProband->pLikelihood[likelihoodIndex].wtslot.
			    weight;
			  KLOG (LOGLIKELIHOOD, LOGDEBUG,
				"\t\t likelihood (%d) = %e\n",
				ppairMatrix[multiLocusPhase[proband]]
				[multiLocusPhase[spouse]].likelihoodIndex,
				ppairMatrix[multiLocusPhase[proband]]
				[multiLocusPhase[spouse]].slot.likelihood);
			}
#else
		      ppairMatrix[multiLocusPhase[proband]]
			[multiLocusPhase[spouse]].slot.likelihood =
			ppairMatrix[multiLocusPhaseFlip[proband]]
			[multiLocusPhase[spouse]].slot.likelihood;
		      pProband->pLikelihood[multiLocusIndex[proband]].wtslot.
			weight =
			pProband->pLikelihood[likelihoodIndex].wtslot.weight;
		      KLOG (LOGLIKELIHOOD, LOGDEBUG,
			    "\t\t likelihood (%d) = %e\n",
			    ppairMatrix[multiLocusPhase[proband]]
			    [multiLocusPhase[spouse]].likelihoodIndex,
			    ppairMatrix[multiLocusPhase[proband]]
			    [multiLocusPhase[spouse]].likelihood);
#endif
		      calculateFlag = 0;
		    }
		}		/* end of finding patterns */
	    }			/* end of proband is a parent */

	  if (calculateFlag == 1)
	    {
#if 1
	      calculate_likelihood(pNucFam, pHaplo, pProband, locus, 
				   multiLocusIndex, multiLocusPhase, dWeight);
#endif 

#if 0
	      for (i = DAD; i <= MOM; i++)
		{
		  if (pParent[i]->touchedFlag == TRUE)
		    {
		      /* we have worked on this parent before */
		      if (pParent[i] == pProband)
			{
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    newWeightPolynomial[i] = constantExp (1);
			  else
			    newWeight[i] = 1.0;
#else
			  newWeight[i] = 1.0;
#endif
			}
		      else
			{
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    newWeightPolynomial[i] = timesExp (2, pParent[i]->pLikelihood[multiLocusIndex[i]].lkslot.likelihoodPolynomial, 1, pParent[i]->pLikelihood[multiLocusIndex[i]].wtslot.weightPolynomial, 1, 0);	//Dec 24
			  else
			    newWeight[i] =
			      pParent[i]->pLikelihood[multiLocusIndex[i]].
			      lkslot.likelihood *
			      pParent[i]->pLikelihood[multiLocusIndex[i]].
			      wtslot.weight;
#else
			  newWeight[i] =
			    pParent[i]->pLikelihood[multiLocusIndex[i]].
			    lkslot.likelihood *
			    pParent[i]->pLikelihood[multiLocusIndex[i]].
			    wtslot.weight;
#endif

			}
		    }
		  else		/* first time we work on this parent */
		    {
		      if (pNucFam->pParents[i]->pParents[DAD] == NULL)
			{
			  /* this parent is a founder. the weight is just the multiplication of 
			     genotype weights at each locus under LE. The weight should have been passed 
			     in as an input */
			  if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
			    {
#ifndef NO_POLYNOMIAL

			      if (modelOptions.polynomial == TRUE)
				newWeightPolynomial[i] =
				  (Polynomial *) dWeight[i];
			      else
				newWeight[i] = *((double *) dWeight + i);
#else
			      newWeight[i] = *((double *) dWeight + i);
#endif
			    }
			  else	/* founder under LD */
			    {
#ifndef NO_POlYNOMIAL
			      if (modelOptions.polynomial == TRUE)
				{
				  get_haplotype_freq (locus, i,
						      &newWeightPolynomial[i],
						      pHaplo);
				}
			      else
				{
				  get_haplotype_freq (locus, i, &newWeight[i],
						      pHaplo);
				}
#else
			      get_haplotype_freq (locus, i, &newWeight[i],
						  pHaplo);
#endif
			    }	/* end of founder and LD */
			}	/* founder */
		      else	/* first time on this nonfounder parent */
			{
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    newWeightPolynomial[i] = constantExp (1.0);
			  else
			    newWeight[i] = 1.0;
#else
			  newWeight[i] = 1.0;
#endif
			}
		    }		/* end of first time on this parent */
		  if (pParent[i]->touchedFlag != TRUE
		      && pParent[i] == pProband)
		    {
#ifndef NO_POLYNOMIAL
		      if (modelOptions.polynomial == TRUE)
			{
			  pParent[i]->pLikelihood[multiLocusIndex[i]].wtslot.
			    weightPolynomial = newWeightPolynomial[i];
			  newWeightPolynomial[i] = constantExp (1.0);
			}
		      else
			{
			  pParent[i]->pLikelihood[multiLocusIndex[i]].wtslot.
			    weight = newWeight[i];
			  newWeight[i] = 1.0;
			}
#else
		      pParent[i]->pLikelihood[multiLocusIndex[i]].wtslot.
			weight = newWeight[i];
		      newWeight[i] = 1.0;
#endif
		    }

		  /* need to multiply the penetrance if disease locus is in and 
		   * we haven't calculated any likelihood on this parent before 
		   */
		  traitLocus = locusList->traitLocusIndex;
		  if (traitLocus >= 0 && pParent[i]->touchedFlag != TRUE)
		    {
		      genoIndex = pHaplo->pParentalPairInd[traitLocus];
#ifndef NO_POLYNOMIAL
		      if (modelOptions.polynomial == TRUE)
			penetrancePolynomial[i] =
			  pHaplo->ppParentalPair[traitLocus][genoIndex].
			  pGenotype[i]->penslot.penetrancePolynomial;
		      else
			penetrance[i] =
			  pHaplo->ppParentalPair[traitLocus][genoIndex].
			  pGenotype[i]->penslot.penetrance;
#else
		      penetrance[i] =
			pHaplo->ppParentalPair[traitLocus][genoIndex].
			pGenotype[i]->penslot.penetrance;
#endif
		    }
		  else
		    {
#ifndef NO_POLYNOMIAL
		      if (modelOptions.polynomial == TRUE)
			penetrancePolynomial[i] = constantExp (1);
		      else
			penetrance[i] = 1.0;
#else
		      penetrance[i] = 1.0;
#endif

		    }
		}		/* loop over each parent */

	      /* now work on the children conditional on this parental pair */
	      childProduct = 1;
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		childProductPolynomial = constantExp (1);
#endif
	      for (child = 0; child < pNucFam->numChildren; child++)
		{
		  xmissionIndex[DAD] = 0;
		  xmissionIndex[MOM] = 0;

#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      sumPolynomial = constantExp (0.0);
		      loop_child_multi_locus_genotype (pNucFam->
						       ppChildrenList[child],
						       pProband, pHaplo,
						       child, 0, 0,
						       &sumPolynomial,
						       xmissionIndex);
		      childProductPolynomial =
			timesExp (2,
				  childProductPolynomial,
				  1, sumPolynomial, 1, 0);
		    }
		  else
		    {
		      sum = 0;
		      loop_child_multi_locus_genotype (pNucFam->
						       ppChildrenList[child],
						       pProband, pHaplo,
						       child, 0, 0,
						       &sum, xmissionIndex);

		      childProduct *= sum;
		    }
#else
		  sum = 0;
		  loop_child_multi_locus_genotype (pNucFam->
						   ppChildrenList[child],
						   pProband, pHaplo,
						   child, 0, 0,
						   &sum, xmissionIndex);

		  childProduct *= sum;
#endif
		}		/* looping over all children */
	      ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
		likelihoodIndex = multiLocusIndex[proband];
	      ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
		count = 1;
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  ppairMatrix[multiLocusPhase[proband]]
		    [multiLocusPhase[spouse]].slot.likelihoodPolynomial =
		    timesExp (5,
			      newWeightPolynomial[proband],
			      1,
			      newWeightPolynomial[spouse],
			      1,
			      penetrancePolynomial[proband],
			      1,
			      penetrancePolynomial[spouse],
			      1, childProductPolynomial, 1, 0);
		  KLOG (LOGLIKELIHOOD, LOGDEBUG,
			"\t\t likelihood (%d) = %e\n",
			ppairMatrix[multiLocusPhase[proband]][multiLocusPhase
							      [spouse]].
			likelihoodIndex,
			evaluateValue(ppairMatrix[multiLocusPhase[proband]]
				      [multiLocusPhase[spouse]].slot.
				      likelihoodPolynomial));
		}
	      else
		{
		  /* save it */
		  ppairMatrix[multiLocusPhase[proband]][multiLocusPhase
							[spouse]].slot.
		    likelihood =
		    newWeight[proband] * newWeight[spouse] *
		    penetrance[proband] * penetrance[spouse] * childProduct;
		  KLOG (LOGLIKELIHOOD, LOGDEBUG,
			"\t\t likelihood (%d) = %e\n",
			ppairMatrix[multiLocusPhase[proband]][multiLocusPhase
							      [spouse]].
			likelihoodIndex,
			ppairMatrix[multiLocusPhase[proband]][multiLocusPhase
							      [spouse]].slot.
			likelihood);
		}
#else
	      /* save it */
	      ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
		slot.likelihood = newWeight[proband] * newWeight[spouse] *
		penetrance[proband] * penetrance[spouse] * childProduct;
	      KLOG (LOGLIKELIHOOD, LOGDEBUG, "\t\t likelihood (%d) = %e\n",
		    ppairMatrix[multiLocusPhase[proband]][multiLocusPhase
							  [spouse]].
		    likelihoodIndex,
		    ppairMatrix[multiLocusPhase[proband]][multiLocusPhase
							  [spouse]].slot.
		    likelihood);
#endif
#endif /* if 0 */

	    }			/* fresh likelihood calculation */

	}			/* end of processing a complete multilocus parental pair */

    }				/* move on to next pair on this locus */

}				/* end of loop_phases() */

int
calculate_likelihood (NuclearFamily * pNucFam, ParentalPairSpace * pHaplo,
		      Person * pProband, int locus, int multiLocusIndex[2],
		      int multiLocusPhase[2], void *dWeight[2])
		      
{
  int i;
  int proband;
  int spouse;
  Person *pParent[2];
  double newWeight[2];
  double sum;
  int child;
  int xmissionIndex[2];
  double childProduct;
  double penetrance[2];
  int traitLocus;
  int genoIndex;
#ifndef NO_POlYNOMIAL
  Polynomial *newWeightPolynomial[2];
  Polynomial *penetrancePolynomial[2];
  Polynomial *sumPolynomial;
  Polynomial *childProductPolynomial = NULL;
#endif

  pParent[DAD] = pNucFam->pParents[DAD];
  pParent[MOM] = pNucFam->pParents[MOM];
  proband = pNucFam->head;
  spouse = pNucFam->spouse;
  for (i = DAD; i <= MOM; i++)
    {
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	newWeightPolynomial[i] = constantExp (1);
      else
	newWeight[i] = 1.0;
#else
      newWeight[i] = 1.0;
#endif
      if (pParent[i]->touchedFlag == TRUE)
	{
	  /* we have worked on this parent before */
	  if (pParent[i] != pProband)
	    {
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		newWeightPolynomial[i] = timesExp (2, pParent[i]->pLikelihood[multiLocusIndex[i]].lkslot.likelihoodPolynomial, 1, pParent[i]->pLikelihood[multiLocusIndex[i]].wtslot.weightPolynomial, 1, 0);	//Dec 24
	      else
		newWeight[i] =
		  pParent[i]->pLikelihood[multiLocusIndex[i]].lkslot.
		  likelihood *
		  pParent[i]->pLikelihood[multiLocusIndex[i]].wtslot.weight;
#else
	      newWeight[i] =
		pParent[i]->pLikelihood[multiLocusIndex[i]].lkslot.
		likelihood *
		pParent[i]->pLikelihood[multiLocusIndex[i]].wtslot.weight;
#endif

	    }
	}
      else			/* first time we work on this parent */
	{
	  if (pNucFam->pParents[i]->pParents[DAD] == NULL)
	    {
	      /* this parent is a founder. the weight is just the multiplication of 
	         genotype weights at each locus under LE. The weight should have been passed 
	         in as an input */
	      if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
		{
#ifndef NO_POLYNOMIAL

		  if (modelOptions.polynomial == TRUE)
		    newWeightPolynomial[i] = (Polynomial *) dWeight[i];
		  else
		    newWeight[i] = *((double *) dWeight + i);
#else
		  newWeight[i] = *((double *) dWeight + i);
#endif
		}
	      else if(pParent[i]->loopBreaker == 0)		/* founder under LD */
		{
#ifndef NO_POlYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      get_haplotype_freq (locus, i, &newWeightPolynomial[i],
					  pHaplo);
		    }
		  else
		    {
		      get_haplotype_freq (locus, i, &newWeight[i], pHaplo);
		    }
#else
		  get_haplotype_freq (locus, i, &newWeight[i], pHaplo);
#endif
		}		/* end of founder and LD */
	    }			/* founder */
	}			/* end of first time on this parent */
      if (pParent[i]->touchedFlag != TRUE && pParent[i] == pProband)
	{
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      pParent[i]->pLikelihood[multiLocusIndex[i]].wtslot.
		weightPolynomial = newWeightPolynomial[i];
	      newWeightPolynomial[i] = constantExp (1.0);
	    }
	  else
	    {
	      pParent[i]->pLikelihood[multiLocusIndex[i]].wtslot.weight =
		newWeight[i];
	      newWeight[i] = 1.0;
	    }
#else
	  pParent[i]->pLikelihood[multiLocusIndex[i]].wtslot.weight =
	    newWeight[i];
	  newWeight[i] = 1.0;
#endif
	}

      /* need to multiply the penetrance if disease locus is in and 
       * we haven't calculated any likelihood on this parent before 
       */
      traitLocus = locusList->traitLocusIndex;
      if (traitLocus >= 0 && pParent[i]->touchedFlag != TRUE && 
	  (pParent[i]->loopBreaker == 0 || pParent[i]->pParents[DAD] != NULL))
	{
	  genoIndex = pHaplo->pParentalPairInd[traitLocus];
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    penetrancePolynomial[i] =
	      pHaplo->ppParentalPair[traitLocus][genoIndex].
	      pGenotype[i]->penslot.penetrancePolynomial;
	  else
	    penetrance[i] =
	      pHaplo->ppParentalPair[traitLocus][genoIndex].
	      pGenotype[i]->penslot.penetrance;
#else
	  penetrance[i] =
	    pHaplo->ppParentalPair[traitLocus][genoIndex].
	    pGenotype[i]->penslot.penetrance;
#endif
	}
      else
	{
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    penetrancePolynomial[i] = constantExp (1);
	  else
	    penetrance[i] = 1.0;
#else
	  penetrance[i] = 1.0;
#endif

	}
    }				/* loop over each parent */

  /* now work on the children conditional on this parental pair */
  childProduct = 1;
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    childProductPolynomial = constantExp (1);
#endif
  for (child = 0; child < pNucFam->numChildren; child++)
    {
      xmissionIndex[DAD] = 0;
      xmissionIndex[MOM] = 0;

#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	{
	  sumPolynomial = constantExp (0.0);
	  loop_child_multi_locus_genotype (pNucFam->
					   ppChildrenList[child],
					   pProband, pHaplo,
					   child, 0, 0,
					   &sumPolynomial, xmissionIndex);
	  childProductPolynomial =
	    timesExp (2, childProductPolynomial, 1, sumPolynomial, 1, 0);
	}
      else
	{
	  sum = 0;
	  loop_child_multi_locus_genotype (pNucFam->
					   ppChildrenList[child],
					   pProband, pHaplo,
					   child, 0, 0, &sum, xmissionIndex);

	  childProduct *= sum;
	}
#else
      sum = 0;
      loop_child_multi_locus_genotype (pNucFam->
				       ppChildrenList[child],
				       pProband, pHaplo,
				       child, 0, 0, &sum, xmissionIndex);

      childProduct *= sum;
#endif
    }				/* looping over all children */

  /* results processing */
  ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
    likelihoodIndex = multiLocusIndex[proband];
  ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].count = 1;
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      ppairMatrix[multiLocusPhase[proband]]
	[multiLocusPhase[spouse]].slot.likelihoodPolynomial =
	timesExp (5,
		  newWeightPolynomial[proband],
		  1,
		  newWeightPolynomial[spouse],
		  1,
		  penetrancePolynomial[proband],
		  1,
		  penetrancePolynomial[spouse],
		  1, childProductPolynomial, 1, 0);
      KLOG (LOGLIKELIHOOD, LOGDEBUG, "\t\t likelihood (%d) = %e\n",
	    ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
	    likelihoodIndex,
	    evaluateValue(ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
			  slot.likelihoodPolynomial));
    }
  else
    {
      /* save it */
      ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
	slot.likelihood = newWeight[proband] * newWeight[spouse] *
	penetrance[proband] * penetrance[spouse] * childProduct;
      KLOG(LOGLIKELIHOOD, LOGDEBUG, "Parents: DAD(%s) weight %e   MOM(%s) weight %e \n",
	   pParent[DAD]->sID, newWeight[proband], pParent[MOM]->sID, newWeight[spouse]);
      KLOG(LOGLIKELIHOOD, LOGDEBUG, "Parents: DAD(%s) pen %e   MOM(%s) pen %e \n",
	   pParent[DAD]->sID, penetrance[proband], pParent[MOM]->sID, penetrance[spouse]);
      KLOG (LOGLIKELIHOOD, LOGDEBUG, "\t\t likelihood (%d) = %e\n",
	    ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
	    likelihoodIndex,
	    ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
	    slot.likelihood);
    }
#else
  /* save it */
  ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
    slot.likelihood = newWeight[proband] * newWeight[spouse] *
    penetrance[proband] * penetrance[spouse] * childProduct;
  KLOG (LOGLIKELIHOOD, LOGDEBUG, "\t\t likelihood (%d) = %e\n",
	ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
	likelihoodIndex,
	ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].slot.
	likelihood);
#endif

  return 0;
}


/*
 * Get haplotype frequency for founders under LD analysis */
void
get_haplotype_freq (int locus, int parent, void *freqPtr,
		    ParentalPairSpace * pHaplo)
{
  int origLocus1, origLocus2;	/* locus indices in the original locus list for the two loci in LD */
  double freq[2] = { 0, 0 };	/* variable to store the calculated frequency */
#ifndef NO_POLYNOMIAL
  Polynomial *freqPolynomial[2] = { NULL, NULL };
  char vName[100];
#endif
  ParentalPair *pPair1;		/* parental pair for one locus */
  ParentalPair *pPair2;		/* parental pair for the other locus */
  Locus *pLocus1;		/* first locus */
  Locus *pLocus2;		/* the other locus */
  int i, k, l;
  LDLoci *pLDLoci;		/* structure contains LD parameter values */
  int alleleID1, alleleID2;	/* allele IDs */
  AlleleSet *pAlleleSet1, *pAlleleSet2;	/* allele sets */
  int allele1, allele2;		/* */


#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      freqPolynomial[0] = constantExp (0);
      freqPolynomial[1] = constantExp (0);
    }
#endif


  /* locus index in the original locus list for the first locus */
  origLocus1 = locusList->pLocusIndex[locus - 1];
  /* locus index in the original locus list for the second locus */
  origLocus2 = locusList->pLocusIndex[locus];
  pLocus1 = originalLocusList.ppLocusList[origLocus1];
  pLocus2 = originalLocusList.ppLocusList[origLocus2];
  /* find the parameter values for these two loci */
  pLDLoci = find_LD_loci (origLocus1, origLocus2);
  KASSERT (pLDLoci != NULL,
	   "Can't find LD parameter between loci %d,%d.\n",
	   origLocus1, origLocus2);
  /* now find the corresponding haplotype frequency : 2 haplotypes 
   * paternal haplotype & maternal haplotype */
  pPair1 =
    &pHaplo->ppParentalPair[locus - 1][pHaplo->pParentalPairInd[locus - 1]];
  pPair2 = &pHaplo->ppParentalPair[locus][pHaplo->pParentalPairInd[locus]];
  for (i = DAD; i <= MOM; i++)
    {
      /* allele ID in the first locus */
      alleleID1 = pPair1->pGenotype[parent]->allele[i];
      /* allele ID in the second locus */
      alleleID2 = pPair2->pGenotype[parent]->allele[i];
      pAlleleSet1 = pLocus1->ppAlleleSetList[alleleID1 - 1];
      pAlleleSet2 = pLocus2->ppAlleleSetList[alleleID2 - 1];
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	{
	  freqPolynomial[i] = constantExp (0);
	}
      else
	freq[i] = 0;
#else
      freq[i] = 0;
#endif

      for (k = 0; k < pAlleleSet1->numAllele; k++)
	{
	  for (l = 0; l < pAlleleSet2->numAllele; l++)
	    {
	      allele1 = pAlleleSet1->pAlleles[k];
	      allele2 = pAlleleSet2->pAlleles[l];
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  sprintf (vName, "ppHaploFreq[%d][%d]",
			   allele1 - 1, allele2 - 1);
		  freqPolynomial[i] =
		    plusExp (2, 1.0, freqPolynomial[i], 1.0,
			     variableExp (&pLDLoci->
					  ppHaploFreq[allele1 - 1][allele2 -
								   1], NULL,
					  'D', vName), 1);
		}
	      else
		freq[i] += pLDLoci->ppHaploFreq[allele1 - 1][allele2 - 1];
#else
	      freq[i] += pLDLoci->ppHaploFreq[allele1 - 1][allele2 - 1];
#endif
	    }
	}

#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	{

	  *(Polynomial **) freqPtr =
	    timesExp (2, freqPolynomial[0], 1, freqPolynomial[1], 1, 0);
	}
      else
	{
	  *(double *) freqPtr = freq[0] * freq[1];
	}
#else
      *(double *) freqValue = freq[0] * freq[1];
#endif
    }				/* end of loop of parents */
}

/* loop over a child's list of genotypes that are compatible with the parental pair
 * retrieve the transmission probability saved in the transmission matrix 
 * sum the likelihood for each genotype configuration
 */
int
loop_child_multi_locus_genotype (Person * pChild,
				 Person * pProband,
				 ParentalPairSpace * pHaplo,
				 int child,
				 int locus,
				 int multiLocusIndex,
				 void *pSum, int xmissionIndex[2])
{
  int i;
  Genotype *pGeno;		/* this child's genotype at current locus */
  int genoIndex;		/* index of the geno in the list */
  /* parental pair we are working on now */
  ParentalPair *pParentalPair =
    &pHaplo->ppParentalPair[locus][pHaplo->pParentalPairInd[locus]];
  int parent;			/* parent - DAD or MOM */
  int origLocus = locusList->pLocusIndex[locus];
  double newProb = 1;
#ifndef NO_POLYNOMIAL
  Polynomial *newProbPolynomial = NULL;
#endif
  int newChromosome[2];
  int numGenotype;
  int multiLocusIndex2;
  int newXmissionIndex[2];
  ParentalPair *pTraitParentalPair;

  /* number of possible genotypes at this locus for this child */
  numGenotype = pChild->pSavedNumGenotype[origLocus];
  /* child's conditional likelihood offset for the multilocus genotype this function is building */
  multiLocusIndex2 = multiLocusIndex * numGenotype;
  /* build the index to xmission matrix for paternal inheritance and maternal inheritance */
  xmissionIndex[DAD] *= 3;
  xmissionIndex[MOM] *= 3;

  /* loop through all of this child's compatible genotypes at this locus */
  for (i = 0; i < pParentalPair->pChildGenoLen[child]; i++)
    {
      pGeno = pParentalPair->pppChildGenoList[child][i];
      /* record the index to the genotype list for this child */
      pHaplo->pChildGenoInd[locus] = i;
      KLOG (LOGLIKELIHOOD, LOGDEBUG,
	    "\t child (%s) locus %4d -> %4d|%-4d \n",
	    pChild->sID, locusList->pLocusIndex[locus], pGeno->allele[DAD],
	    pGeno->allele[MOM]);
      /* record this child's conditional likelihood index */
      multiLocusIndex = multiLocusIndex2 + pGeno->position;

      if (locus < locusList->numLocus)
	{
	  /* check the transmission probability */
	  for (parent = DAD; parent <= MOM; parent++)
	    {
	      newChromosome[parent] =
		pParentalPair->ppChildInheritance[parent][child][i];
	      /* xmissionIndex has already been multiplied by 3 before the loop */
	      newXmissionIndex[parent] =
		xmissionIndex[parent] + newChromosome[parent] - 1;
	    }			/* looping paternal and maternal chromosomes */
	  if (locus < locusList->numLocus - 1)
	    {
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
//                fprintf(stderr,"AAlocus=%d newPen=%f pSum=%f newProb=%f i=%d\n",
//                        locus,evaluateValue(newPenetrancePolynomial),
//                evaluateValue(*(Polynomial **)pSum),evaluateValue(newProbPolynomial),i);      
		  loop_child_multi_locus_genotype (pChild, pProband, pHaplo,
						   child, locus + 1,
						   multiLocusIndex,
						   pSum, newXmissionIndex);
//                fprintf(stderr,"BBlocus=%d newPen=%f pSum=%f newProb=%f\n",
//                        locus,evaluateValue(newPenetrancePolynomial),
//                evaluateValue(*(Polynomial **)pSum),evaluateValue(newProbPolynomial));                   
		}
	      else
		{
//                   fprintf(stderr,"AAlocus=%d newPen=%f, pSum=%f newProb=%f i=%d\n",
//                                   locus, newPenetrance, *((double *)pSum),newProb,i);
		  loop_child_multi_locus_genotype (pChild, pProband, pHaplo,
						   child, locus + 1,
						   multiLocusIndex,
						   pSum, newXmissionIndex);
//                   fprintf(stderr,"BBlocus=%d newPen=%f, pSum=%f newProb=%f\n",
//                                   locus,newPenetrance, *((double *)pSum),newProb);
		}
#else
	      loop_child_multi_locus_genotype (pChild, pProband, pHaplo,
					       child, locus + 1,
					       multiLocusIndex,
					       pSum, newXmissionIndex);
#endif
	    }
	  else
	    {
	      /* get the transmission probability from the matrix */

#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  newProbPolynomial =
		    timesExp (2,
			      xmissionMatrix[newXmissionIndex[DAD]].slot.
			      probPoly[DAD + 1], 1,
			      xmissionMatrix[newXmissionIndex[MOM]].slot.
			      probPoly[MOM + 1], 1, 0);
		  KLOG (LOGLIKELIHOOD, LOGDEBUG,
			"\t xmission prob: %f = %f * %f\n", 
			evaluateValue(newProbPolynomial),
			evaluateValue(xmissionMatrix[newXmissionIndex[DAD]].slot.probPoly[DAD +
											  1]),
			evaluateValue(xmissionMatrix[newXmissionIndex[MOM]].slot.probPoly[MOM +
									    1]));
		}
	      else
		{
		  newProb =
		    xmissionMatrix[newXmissionIndex[DAD]].slot.prob[DAD + 1] *
		    xmissionMatrix[newXmissionIndex[MOM]].slot.prob[MOM + 1];
		  KLOG (LOGLIKELIHOOD, LOGDEBUG,
			"\t xmission prob: %f = %f * %f\n", newProb,
			xmissionMatrix[newXmissionIndex[DAD]].slot.prob[DAD +
									1],
			xmissionMatrix[newXmissionIndex[MOM]].slot.prob[MOM +
									1]);
//                 fprintf(stderr,"newProb=%f newPenetrance=%f\n",newProb,newPenetrance);
		}
#else
	      newProb =
		xmissionMatrix[newXmissionIndex[DAD]].slot.prob[DAD +
								1] *
		xmissionMatrix[newXmissionIndex[MOM]].slot.prob[MOM + 1];
	      KLOG (LOGLIKELIHOOD, LOGDEBUG,
		    "\t xmission prob: %f = %f * %f\n", newProb,
		    xmissionMatrix[newXmissionIndex[DAD]].slot.prob[DAD + 1],
		    xmissionMatrix[newXmissionIndex[MOM]].slot.prob[MOM + 1]);
#endif

	      /* we have completed one multilocus genotype for this child */
	      /* check whether we already have some information about this kid
	       * we should have if this kid is a connector to another nuclear
	       * family we have processed before */
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  if (pChild != pProband)
		    {
		      /* the child is not a proband */
		      if (pChild->touchedFlag == 1)
			{
			  /* some likelihood calculation has been calculated for this child */
			  *(Polynomial **) pSum = plusExp (2, 1.0, *(Polynomial **) pSum, 1.0, timesExp (2, newProbPolynomial, 1, pChild->pLikelihood[multiLocusIndex].lkslot.likelihoodPolynomial, 1, 0),	//end of timesExp
							   1);	//end of plusExp
			  KLOG(LOGLIKELIHOOD, LOGDEBUG, 
			       "\t use already calculated child prob %e \n", 
			       evaluateValue(pChild->pLikelihood[multiLocusIndex].lkslot.likelihoodPolynomial));
			}
		      else if (locusList->traitLocusIndex >= 0)
			/* first time working on this child's current multilocus genotype 
			 * and we need to consider penetrance 
			 */
			{
			  genoIndex =
			    pHaplo->pChildGenoInd[locusList->traitLocusIndex];
			  pTraitParentalPair =
			    &pHaplo->ppParentalPair[locusList->
						    traitLocusIndex][pHaplo->
								     pParentalPairInd
								     [locusList->
								      traitLocusIndex]];
			  *(Polynomial **) pSum = plusExp (2, 1.0, *(Polynomial **) pSum, 1.0, timesExp (2, newProbPolynomial, 1, pTraitParentalPair->pppChildGenoList[child][genoIndex]->penslot.penetrancePolynomial, 1, 0),	//end of timesExp
							   1);	//end of plusExp
			}
		      else
			{
			  /* no trait locus and new to this child */
			  *(Polynomial **) pSum = plusExp (2,
							   1.0,
							   *(Polynomial **)
							   pSum, 1.0,
							   newProbPolynomial,
							   1);
			}
		    }
		  else		/* this child is not proband */
		    {
		      *(Polynomial **) pSum = plusExp (2,
						       1.0,
						       *(Polynomial **) pSum,
						       1.0,
						       newProbPolynomial, 1);
		    }
		  KLOG(LOGLIKELIHOOD, LOGDEBUG, 
		       "\t child sum %e \n", 
		       evaluateValue(*(Polynomial **)pSum));
		}
	      else		/* PE is not turned on */
		{
		  if (pChild != pProband)
		    {
		      /* the child is not a proband */
		      if (pChild->touchedFlag == 1)
			{
			  /* some likelihood calculation has been done for this child */
			  *(double *) pSum += newProb *
			    pChild->pLikelihood[multiLocusIndex].lkslot.
			    likelihood;
			  KLOG(LOGLIKELIHOOD, LOGDEBUG, 
			       "\t use already calculated child prob %e \n", 
			       pChild->pLikelihood[multiLocusIndex].lkslot.likelihood);
			}
		      else if (locusList->traitLocusIndex >= 0)
			/* first time working on this child's current multilocus genotype 
			 * and we need to consider penetrance 
			 */
			{
			  genoIndex =
			    pHaplo->pChildGenoInd[locusList->traitLocusIndex];
			  pTraitParentalPair =
			    &pHaplo->ppParentalPair[locusList->
						    traitLocusIndex][pHaplo->
								     pParentalPairInd
								     [locusList->
								      traitLocusIndex]];
			  *(double *) pSum +=
			    newProb *
			    pTraitParentalPair->
			    pppChildGenoList[child][genoIndex]->penslot.
			    penetrance;
			  KLOG(LOGLIKELIHOOD, LOGDEBUG, "child penetrance %e\n", 
			       pTraitParentalPair->
			       pppChildGenoList[child][genoIndex]->penslot.
			       penetrance);
			       
			}
		      else
			{
			  *(double *) pSum += newProb;
			}
		    }
		  else		/* this child is proband */
		    {
		      /* penetrance if applicable will be figured into later */
		      *(double *) pSum += newProb;
		    }
		  KLOG(LOGLIKELIHOOD, LOGDEBUG, 
		       "\t child sum %e \n", 
		       *(double *)pSum);
		}
#else /* PE is not compiled in */
	      if (pChild != pProband)
		{
		  /* the child is not a proband */
		  if (pChild->touchedFlag == 1)
		    {
		      /* some likelihood calculation has been done for this child */
		      *(double *) pSum += newProb *
			pChild->pLikelihood[multiLocusIndex].lkslot.
			likelihood;
		      KLOG(LOGLIKELIHOOD, LOGDEBUG, 
			   "\t use already calculated child prob %e \n", 
			   pChild->pLikelihood[multiLocusIndex].lkslot.likelihood);
		    }
		  else if (locusList->traitLocusIndex >= 0)
		    /* first time working on this child's current multilocus genotype 
		     * and we need to consider penetrance 
		     */
		    {
		      genoIndex =
			pHaplo->pChildGenoInd[locusList->traitLocusIndex];
		      pTraitParentalPair =
			&pHaplo->ppParentalPair[locusList->
						traitLocusIndex][pHaplo->
								 pParentalPairInd
								 [locusList->
								  traitLocusIndex]];
		      *(double *) pSum +=
			newProb *
			pTraitParentalPair->
			pppChildGenoList[child][genoIndex]->penslot.
			penetrance;
		    }
		  else
		    {
		      *(double *) pSum += newProb;
		    }
		}
	      else		/* this child is proband */
		{
		  /* penetrance if applicable will be figured into later */
		  *(double *) pSum += newProb;
		}
		  KLOG(LOGLIKELIHOOD, LOGDEBUG, 
		       "\t child sum %e \n", 
		       *(double *)pSum);
#endif

	    }			/* end of processing one complete multilocus genotype */
	}
    }				/* end of looping the genotype list for this locus */
  return 0;
}


/* A recursive call to build transmission probability matrix 
 * pMatrix - pass in matrix pointer - This should have been pre-allocated
 * totalLoci - total number of loci
 * loc - current locus 
 * prob - 
 */

int
populate_xmission_matrix (XMission * pMatrix, int totalLoci,
			  void *prob[3], void *prob2[3],
			  void *hetProb[3],
			  int cellIndex,
			  int lastHetLoc, int prevPattern, int loc)
{
  int pattern;
#ifndef NO_POLYNOMIAL
  Polynomial *newProbPoly[3];
  Polynomial *newProbPoly2[3];
  Polynomial *newHetProbPoly[3];
#endif
  double newProb[3];
  double *newProbPtr[3] = { &newProb[0], &newProb[1], &newProb[2] };
  double newProb2[3];
  double *newProbPtr2[3] = { &newProb2[0], &newProb2[1], &newProb2[2] };
  double *newHetProbPtr[3] = { hetProb[0], hetProb[1], hetProb[2] };
  int newCellIndex;
  int newLastHetLoc;
  int i;
  char vName1[100];

  /* at each locus, the inheritance could be paternal only (1), maternal only (2), and
   * both (3) which indicates the parent is homozygous at that locus */
  for (pattern = 1; pattern <= 3; pattern++)
    {
      /* sex averaged or sex specific map */
      for (i = 0; i < 3; i++)
	{
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      newProbPoly[i] = (Polynomial *) prob[i];
	      newProbPoly2[i] = (Polynomial *) prob2[i];
	      newHetProbPoly[i] = (Polynomial *) hetProb[i];
	    }
	  else
	    {
	      newProb[i] = *((double *) prob[i]);
	      newProb2[i] = *((double *) prob2[i]);
	      newHetProbPtr[i] = hetProb[i];
	    }
#else
	  newProb[i] = *((double *) prob[i]);
	  newProb2[i] = *((double *) prob2[i]);
	  newHetProbPtr[i] = hetProb[i];
#endif
	}
      newCellIndex = cellIndex * 3 + pattern - 1;
      newLastHetLoc = lastHetLoc;
      if (pattern != 3)
	{
	  /* parent is not homozygous */
	  if (lastHetLoc != -1)
	    {
	      if (prevPattern != 3)
		{
		  /* previous locus for the parent is het and 
		     current locus pattern is either paternal or maternal */
		  if (prevPattern == pattern)
		    {
		      /* no recombination */
		      for (i = 0; i < 3; i++)
			{
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    {
			      sprintf (vName1, "theta%d_%d", i, loc);
			      newProbPoly[i] = timesExp (2, newProbPoly[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, variableExp (&locusList->pPrevLocusDistance[i][loc], NULL, 'D', vName1), 0), 1, 0);	//Dec 24
			    }
			  else
			    {
			      newProb[i] *=
				(1 - locusList->pPrevLocusDistance[i][loc]);
			    }
#else
			  newProb[i] *=
			    (1 - locusList->pPrevLocusDistance[i][loc]);
#endif
			}
		    }
		  else
		    {
		      /* recombination */
		      for (i = 0; i < 3; i++)
#ifndef NO_POLYNOMIAL
			if (modelOptions.polynomial == TRUE)
			  {
			    sprintf (vName1, "theta%d_%d", i, loc);
			    newProbPoly[i] = timesExp (2, newProbPoly[i], 1, variableExp (&locusList->pPrevLocusDistance[i][loc], NULL, 'D', vName1), 1, 0);	//Dec 24

			  }
			else
			  {
			    newProb[i] *=
			      locusList->pPrevLocusDistance[i][loc];
			  }
#else
			newProb[i] *= locusList->pPrevLocusDistance[i][loc];
#endif
		    }
		}
	      else
		{
		  /* previous locus at parent is homo and current locus is het */
		  for (i = 0; i < 3; i++)
		    {
		      if (pattern == 1)
			{
			  /* paternal inheritance for this locus
			     either no recombination from previous paternal strand 
			     or recombination from previous maternal strand */
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    {
			      sprintf (vName1, "theta%d_%d", i, loc);
			      newProbPoly[i] = plusExp (2,
							1.0,
							timesExp (2,
								  (Polynomial
								   *) prob[i],
								  1,
								  plusExp (2,
									   1.0,
									   constantExp
									   (1.0),
									   -1.0,
									   variableExp
									   (&locusList->
									    pPrevLocusDistance
									    [i]
									    [loc],
									    NULL,
									    'D',
									    vName1),
									   0),
								  1, 0), 1.0,
							timesExp (2,
								  (Polynomial
								   *)
								  prob2[i], 1,
								  variableExp
								  (&locusList->
								   pPrevLocusDistance
								   [i][loc],
								   NULL, 'D',
								   vName1), 1,
								  0), 0);
			    }
			  else
			    newProb[i] = *((double *) prob[i]) *
			      (1 - locusList->pPrevLocusDistance[i][loc]) +
			      *((double *) prob2[i]) *
			      locusList->pPrevLocusDistance[i][loc];
#else
			  newProb[i] = *((double *) prob[i]) *
			    (1 - locusList->pPrevLocusDistance[i][loc]) +
			    *((double *) prob2[i]) *
			    locusList->pPrevLocusDistance[i][loc];
#endif
			}
		      else
			{
			  /* has to be maternal */
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    {
			      sprintf (vName1, "theta%d_%d", i, loc);
			      newProbPoly[i] = plusExp (2,
							1.0,
							timesExp (2,
								  (Polynomial
								   *)
								  prob2[i], 1,
								  plusExp (2,
									   1.0,
									   constantExp
									   (1.0),
									   -1.0,
									   variableExp
									   (&locusList->
									    pPrevLocusDistance
									    [i]
									    [loc],
									    NULL,
									    'D',
									    vName1),
									   0),
								  1, 0), 1.0,
							timesExp (2,
								  (Polynomial
								   *) prob[i],
								  1,
								  variableExp
								  (&locusList->
								   pPrevLocusDistance
								   [i][loc],
								   NULL, 'D',
								   vName1), 1,
								  0), 0);

			    }
			  else
			    newProb[i] = *((double *) prob2[i]) *
			      (1 - locusList->pPrevLocusDistance[i][loc]) +
			      *((double *) prob[i]) *
			      locusList->pPrevLocusDistance[i][loc];

#else
			  newProb[i] = *((double *) prob2[i]) *
			    (1 - locusList->pPrevLocusDistance[i][loc]) +
			    *((double *) prob[i]) *
			    locusList->pPrevLocusDistance[i][loc];
#endif
			}
		    }

		}		/* end of prevPattern is homo and current pattern is het */
	    }			/* end of prevHetLoc != -1 */
	  else
	    {
	      /* we don't have any het locus yet, this locus is the first het */
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  for (i = 0; i < 3; i++)
		    newProbPoly[i] = constantExp (0.5);
		}
	      else
		{
		  for (i = 0; i < 3; i++)
		    newProb[i] = 0.5;
		}
#else
	      for (i = 0; i < 3; i++)
		newProb[i] = 0.5;
#endif
	    }
	  newLastHetLoc = loc;


#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      for (i = 0; i < 3; i++)
		newHetProbPoly[i] = newProbPoly[i];
	    }
	  else
	    {
	      for (i = 0; i < 3; i++)
		newHetProbPtr[i] = newProbPtr[i];
	    }
#else
	  for (i = 0; i < 3; i++)
	    newHetProbPtr[i] = newProbPtr[i];
#endif


	}			/* end of current pattern is not homo */
      else
	{
	  /* current pattern is homo */
	  if (lastHetLoc == -1)
	    /* nothing needs to be done for this locus */
	    ;
	  else
	    {
	      if (loc == totalLoci - 1)
		{
		  /* this is the last locus and it's homo, take the previous het locus */
		  for (i = 0; i < 3; i++)
		    {
#ifndef NO_POLYNOMIAL
		      if (modelOptions.polynomial == TRUE)
			{
			  newProbPoly[i] = (Polynomial *) hetProb[i];
			}
		      else
			newProb[i] = *(double *) hetProb[i];
#else
		      newProb[i] = *(double *) hetProb[i];
#endif
		    }
		}
	      else
		{
		  if (prevPattern == 3)
		    {
		      /* previous locus pattern is homo */
		      for (i = 0; i < 3; i++)
			{

#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    {
			      sprintf (vName1, "theta%d_%d", i, loc);
			      newProbPoly[i] = plusExp (2,
							1.0,
							timesExp (2,
								  (Polynomial
								   *) prob[i],
								  1,
								  plusExp (2,
									   1.0,
									   constantExp
									   (1.0),
									   -1.0,
									   variableExp
									   (&locusList->
									    pPrevLocusDistance
									    [i]
									    [loc],
									    NULL,
									    'D',
									    vName1),
									   0),
								  1, 0), 1.0,
							timesExp (2,
								  (Polynomial
								   *)
								  prob2[i], 1,
								  variableExp
								  (&locusList->
								   pPrevLocusDistance
								   [i][loc],
								   NULL, 'D',
								   vName1), 1,
								  0), 0);
			      newProbPoly2[i] =
				plusExp (2, 1.0,
					 timesExp (2, (Polynomial *) prob2[i],
						   1, plusExp (2, 1.0,
							       constantExp
							       (1.0), -1.0,
							       variableExp
							       (&locusList->
								pPrevLocusDistance
								[i][loc],
								NULL, 'D',
								vName1), 0),
						   1, 0), 1.0, timesExp (2,
									 (Polynomial
									  *)
									 prob
									 [i],
									 1,
									 variableExp
									 (&locusList->
									  pPrevLocusDistance
									  [i]
									  [loc],
									  NULL,
									  'D',
									  vName1),
									 1,
									 0),
					 0);

			    }
			  else
			    {
			      newProb[i] = *(double *) prob[i] *
				(1 - locusList->pPrevLocusDistance[i][loc]) +
				*((double *) prob2[i]) *
				locusList->pPrevLocusDistance[i][loc];

			      newProb2[i] = *(double *) prob2[i] *
				(1 - locusList->pPrevLocusDistance[i][loc]) +
				*((double *) prob[i]) *
				locusList->pPrevLocusDistance[i][loc];
			    }
#else
			  newProb[i] = *(double *) prob[i] *
			    (1 - locusList->pPrevLocusDistance[i][loc]) +
			    *((double *) prob2[i]) *
			    locusList->pPrevLocusDistance[i][loc];

			  newProb2[i] = *(double *) prob2[i] *
			    (1 - locusList->pPrevLocusDistance[i][loc]) +
			    *((double *) prob[i]) *
			    locusList->pPrevLocusDistance[i][loc];
#endif
			}
		    }
		  else
		    {
		      for (i = 0; i < 3; i++)
			{
			  if (prevPattern == 1)
			    {
#ifndef NO_POLYNOMIAL
			      if (modelOptions.polynomial == TRUE)
				{
				  sprintf (vName1, "theta%d_%d", i, loc);
				  newProbPoly[i] =
				    timesExp (2, (Polynomial *) prob[i], 1,
					      plusExp (2, 1.0,
						       constantExp (1.0),
						       -1.0,
						       variableExp
						       (&locusList->
							pPrevLocusDistance[i]
							[loc], NULL, 'D',
							vName1), 0), 1, 0);
				  newProbPoly2[i] =
				    timesExp (2, (Polynomial *) prob[i], 1,
					      variableExp (&locusList->
							   pPrevLocusDistance
							   [i][loc], NULL,
							   'D', vName1), 1,
					      0);
				}
			      else
				{
				  newProb[i] = *(double *) prob[i] *
				    (1 -
				     locusList->pPrevLocusDistance[i][loc]);
				  newProb2[i] =
				    *(double *) prob[i] *
				    locusList->pPrevLocusDistance[i][loc];
				}
#else
			      newProb[i] = *(double *) prob[i] *
				(1 - locusList->pPrevLocusDistance[i][loc]);

			      newProb2[i] = *(double *) prob[i] *
				locusList->pPrevLocusDistance[i][loc];
#endif
			    }
			  else
			    {
#ifndef NO_POLYNOMIAL
			      if (modelOptions.polynomial == TRUE)
				{
				  sprintf (vName1, "theta%d_%d", i, loc);
				  newProbPoly2[i] =
				    timesExp (2, (Polynomial *) prob[i], 1,
					      plusExp (2, 1.0,
						       constantExp (1.0),
						       -1.0,
						       variableExp
						       (&locusList->
							pPrevLocusDistance[i]
							[loc], NULL, 'D',
							vName1), 0), 1, 0);
				  newProbPoly[i] =
				    timesExp (2, (Polynomial *) prob[i], 1,
					      variableExp (&locusList->
							   pPrevLocusDistance
							   [i][loc], NULL,
							   'D', vName1), 1,
					      0);
				}
			      else
				{
				  newProb2[i] = *(double *) prob[i] *
				    (1 -
				     locusList->pPrevLocusDistance[i][loc]);
				  newProb[i] =
				    *(double *) prob[i] *
				    locusList->pPrevLocusDistance[i][loc];
				}
#else
			      newProb2[i] = *(double *) prob[i] *
				(1 - locusList->pPrevLocusDistance[i][loc]);

			      newProb[i] = *(double *) prob[i] *
				locusList->pPrevLocusDistance[i][loc];
#endif
			    }

			}
		    }
		}
	    }
	}

      if (loc == totalLoci - 1)
	{
	  /* we have a complete set of multilocus inheritance pattern */

	  for (i = 0; i < 3; i++)
	    {
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  pMatrix[newCellIndex].slot.probPoly[i] = newProbPoly[i];
		}
	      else
		pMatrix[newCellIndex].slot.prob[i] = newProb[i];
#else
	      pMatrix[newCellIndex].slot.prob[i] = newProb[i];
#endif
	    }

	}
      else
	{
	  /* move on to next locus */
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      populate_xmission_matrix (pMatrix, totalLoci,
					(void *) newProbPoly,
					(void *) newProbPoly2,
					(void *) newHetProbPoly, newCellIndex,
					newLastHetLoc, pattern, loc + 1);
	    }
	  else
	    populate_xmission_matrix (pMatrix, totalLoci,
				      (void *) newProbPtr,
				      (void *) newProbPtr2,
				      (void *) newHetProbPtr, newCellIndex,
				      newLastHetLoc, pattern, loc + 1);
#else
	  populate_xmission_matrix (pMatrix, totalLoci,
				    newProbPtr, (void *) newProbPtr2,
				    (void *) newHetProbPtr,
				    newCellIndex,
				    newLastHetLoc, pattern, loc + 1);
#endif
	}
    }
  return 0;
}

int
build_xmission_matrix (XMission ** ppMatrix, int totalLoci)
{
  int size;
  //  int i;

  *ppMatrix = NULL;
  size = pow (3, totalLoci);
  /* minimal two loci */
  if (size < 9)
    size = 9;
  *ppMatrix = (XMission *) calloc (size, sizeof (XMission));
  if (*ppMatrix == NULL)
    /* memory allocation failed */
    return -1;

  return 0;
}

void
print_xmission_matrix (XMission * pMatrix, int totalLoci, int loc,
		       int cellIndex, char *pID)
{
  int pattern;
  int newCellIndex;
  int i;

  for (pattern = 1; pattern <= 3; pattern++)
    {
      newCellIndex = cellIndex * 3 + pattern - 1;
      if (pattern == 1)
	{
	  pID[loc] = 'P';
	}
      else if (pattern == 2)
	{
	  pID[loc] = 'M';
	}
      else
	{
	  pID[loc] = 'B';
	}
      if (loc != totalLoci - 1)
	{
	  /* not complete multi-locus haplotype yet */
	  print_xmission_matrix (pMatrix, totalLoci, loc + 1, newCellIndex,
				 pID);
	}
      else
	{
	  /* print the xmission probability out */
	  for (i = 0; i <= loc; i++)
	    {
	      fprintf (stderr, "%c", pID[i]);
	    }
	  fprintf (stderr, ":%f\n", pMatrix[newCellIndex].slot.prob[0]);
	}
    }
}

/* allocate storage for keeping track of het locus in nuclear families 
 * numLocus - number of loci analyzing at a time
 */
void
allocate_nucfam_het (PedigreeSet * pPedigreeList, int numLocus)
{
  int ped;
  Pedigree *pPedigree;
  int fam;
  NuclearFamily *pNucFam;

  for (ped = 0; ped < pPedigreeList->numPedigree; ped++)
    {
      pPedigree = pPedigreeList->ppPedigreeSet[ped];
      for (fam = 0; fam < pPedigree->numNuclearFamily; fam++)
	{
	  pNucFam = pPedigree->ppNuclearFamilyList[fam];
	  if (pNucFam->hetFlag[DAD] == NULL)
	    {
	      pNucFam->hetFlag[DAD] = (int *) calloc (sizeof (int), numLocus);
	      pNucFam->hetFlag[MOM] = (int *) calloc (sizeof (int), numLocus);
	      pNucFam->tmpNumHet[DAD] =
		(int *) calloc (sizeof (int), numLocus);
	      pNucFam->tmpNumHet[MOM] =
		(int *) calloc (sizeof (int), numLocus);
	      pNucFam->relatedPPairStart =
		(int *) calloc (sizeof (int), numLocus);
	      pNucFam->numRelatedPPair =
		(int *) calloc (sizeof (int), numLocus);
	      pNucFam->totalRelatedPPair =
		(int *) calloc (sizeof (int), numLocus);
	    }
	}
    }

}

inline void
clear_ppairMatrix (PPairElement ** ppMatrix)
{
  int i;

  for (i = 0; i <= bitMask[ppairMatrixNumLocus]; i++)
    {
      memset (ppMatrix[i], 0, ppairMatrixRowSize);
    }
}

inline void
initialize_proband_tmpLikelihood (Person * pPerson)
{
  int i;
  int size = pPerson->numConditionals;

  for (i = 0; i < size; i++)
    {
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	pPerson->pLikelihood[i].tmpslot.tmpLikelihoodPolynomial =
	  constantExp (0);
      else
	pPerson->pLikelihood[i].tmpslot.tmpLikelihood = 0;
#else
      pPerson->pLikelihood[i].tmpslot.tmpLikelihood = 0;
#endif
    }
}

void populate_pedigree_loopbreaker_genotype_vector(Pedigree *pPed)
{
  int numLoopBreaker = pPed->numLoopBreaker;
  int i; 
  Person *pLoopBreaker;
  
  for(i=0; i < numLoopBreaker; i++)
    {
      pLoopBreaker = pPed->loopBreakerList[i];
      pLoopBreaker->loopBreakerStruct->numGenotype = 0;
      populate_loopbreaker_genotype_vector(pLoopBreaker, 0);
      pLoopBreaker->loopBreakerStruct->genotypeIndex = 0;
    }
}

void populate_loopbreaker_genotype_vector(Person *pLoopBreaker, int locus)
{
  int origLocus = locusList->pLocusIndex[locus];	/* locus index in the original locus list */
  Genotype *pGenotype;
  int index;
  LoopBreaker *pLoopBreakerStruct;

  pGenotype = pLoopBreaker->ppSavedGenotypeList[origLocus];
  while(pGenotype != NULL)
    {
      pTempGenoVector[locus] = pGenotype;
      if(locus < locusList->numLocus -1)
	{
	  populate_loopbreaker_genotype_vector(pLoopBreaker, locus+1);
	}
      else
	{
	  /* one complete multilocus genotype */
	  pLoopBreakerStruct = pLoopBreaker->loopBreakerStruct;
	  index = pLoopBreakerStruct->numGenotype;
	  memcpy(pLoopBreakerStruct->genotype[index], pTempGenoVector, 
		 sizeof(Genotype *) * locusList->numLocus);
	  pLoopBreakerStruct->numGenotype++;
	  
	}
      pGenotype = pGenotype->pSavedNext;
    }
}

/* this function is not needed - so not finished */
void sync_loopbreaker_duplicates(Pedigree *pPed)
{
  int i;
  Person *pPerson;

  for(i=0; i < pPed->numPerson; i++)
    {
      pPerson = pPed->ppPersonList[i];
      if(pPerson->loopBreaker >=1 && pPerson->pParents[DAD] == NULL)
	{
	}
    }
}

/* initialFlag - TRUE - first vector (index all 0) */
int set_next_loopbreaker_genotype_vector(Pedigree *pPed, int initialFlag)
{
  int numLoopBreaker = pPed->numLoopBreaker;
  int i; 
  Person *pLoopBreaker;
  int found;
  LoopBreaker *loopStruct;
  int index;
  int origLocus;
  int locus;

  /* find the next genotype vector for at least one of the loop breaker */
  KLOG(LOGLIKELIHOOD, LOGDEBUG, "Set next loop breaker genotype\n");
  found = FALSE;
  if(initialFlag != TRUE)
    {
      for(i=0; i < numLoopBreaker; i++)
	{
	  pLoopBreaker = pPed->loopBreakerList[i];
	  loopStruct = pLoopBreaker->loopBreakerStruct;
	  /* increase index */
	  loopStruct->genotypeIndex++;
	  if(loopStruct->genotypeIndex >= loopStruct->numGenotype)
	    {
	      loopStruct->genotypeIndex=0;
	    }
	  else
	    {
	      found = TRUE;
	      break;
	    }
	}
      
    }
  else
    {
      found = TRUE;
#if 0
      for(i=0; i < numLoopBreaker; i++)
	{
	  pLoopBreaker = pPed->loopBreakerList[i];
	  loopStruct = pLoopBreaker->loopBreakerStruct;
	  loopStruct->genotypeIndex = loopStruct->numGenotype -1;
	}
#endif
    }

  if(found == FALSE)
    return -1;

  /* set the genotype list with the selected genotype vector */
  for(i=0; i < numLoopBreaker; i++)
    {
      pLoopBreaker = pPed->loopBreakerList[i];
      loopStruct = pLoopBreaker->loopBreakerStruct;
      index = loopStruct->genotypeIndex;
      KLOG(LOGLIKELIHOOD, LOGDEBUG, 
	   "Fix pedigree %s loop breaker %s to the genotype below (%d/%d):\n", 
	   pPed->sPedigreeID, pLoopBreaker->sID, 
	   loopStruct->genotypeIndex+1, loopStruct->numGenotype);
      for(locus = 0; locus < locusList->numLocus; locus++)
	{
	  origLocus = locusList->pLocusIndex[locus];
	  pLoopBreaker->ppGenotypeList[origLocus] = loopStruct->genotype[index][locus];
	  loopStruct->genotype[index][locus]->pNext = NULL;
	  pLoopBreaker->pNumGenotype[origLocus] = 1;
	  KLOG(LOGLIKELIHOOD, LOGDEBUG, "\t %d-> %d|%d \n", 
	       locus, loopStruct->genotype[index][locus]->allele[DAD],
	       loopStruct->genotype[index][locus]->allele[MOM]);
	}
    }

  /* as the loop breakers are set with fixed genotype, redo genotype elimination (without
   * actually remove any genotype - only the links will get updated */
  for(i=0; i < locusList->numLocus; i++)
    {
      origLocus = locusList->pLocusIndex[i];
      pedigree_genotype_elimination(origLocus, pPed);
    }

  return 0;
}
