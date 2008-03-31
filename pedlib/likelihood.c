
/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.
 * All rights reserved.
 * Permission is hereby given to use this software
 * for non-profit educational purposes only.
 **********************************************************************/

char *likelihoodVersion = "0.34.1";

/*
 * This file contains functions to  compute likelihood for all the pedigrees,
 * peeling procedure etc.
 */
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

/* this is for working out loop breaker's multilocus genotypes */
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

/* this matrix is used to keep the counts of likelihood under different parental pair phases */
PPairElement **ppairMatrix = NULL;

/* maximum is 2^n * sizeof(PPairElement) */
int ppairMatrixRowSize;

/* this is number of loci we analyze at a time */
int ppairMatrixNumLocus;

/* this is used to flip parental pair phases */
int *bitMask = NULL;

#ifndef NO_POLYNOMIAL
Polynomial *constant0Poly;
Polynomial *constant1Poly;
#endif

/* these are work space for likelihood calculation 
 * reluctant to use global variables, but they do save stack and heap space, time to allocate
 */
Person *pChild;
Person *pProband;
ParentalPairSpace *pHaplo;
int child;
void *childSum;
Genotype *pGenotype;
int traitGenoIndex;
ParentalPair *pTraitParentalPair;
int parent;
double newProb = 1;

#ifndef NO_POLYNOMIAL
Polynomial *newProbPolynomial = NULL;
#endif
int newChromosome[2];
int numLocus;
NuclearFamily *pNucFam;

/* 0 - don't keep result, 1 - keep result, 2 - use result */
int calcFlag;

/* the following will facilitate parenatl pattern flip */
typedef struct ChildElement
{
  /* xmission index of the paternal and maternal haplotypes */
  int xmissionIndex[2];
  /* the factor could be 1, or penetrance or existing likelihood for this child */
  union
  {
    double factor;
    Polynomial *factorPolynomial;
  } fslot;
} ChildElement;

ChildElement *likelihoodChildElements;
int maxChildElements;

int *likelihoodChildCount;
int maxChildren;

int multCount;

/* function prototypes */
void recalculate_child_likelihood (int flipMask[2], void *childProduct);
int peel_graph (NuclearFamily * pNucFam, Person * pProband,
		int peelingDirection);
int compute_nuclear_family_likelihood (int peelingDirection);
int loop_parental_pair (int locus, int multiLocusIndex[2], void *dWeight[2]);
int loop_child_multi_locus_genotype (int locus, int multiLocusIndex,
				     int xmissionIndex[2]);
void get_haplotype_freq (int locus, int parent, void *freqPtr);

int
loop_child_proband_genotype (int peelingDirection, int locus,
			     int multiLocusIndex);

void
loop_phases (int locus, int multiLocusIndex[2], int multiLocusPhase[2],
	     int flipMask[2], void *dWeight[2]);
int calculate_likelihood (int multiLocusIndex[2], int multiLocusPhase[2],
			  void *dWeight[2], void *childProductPtr);

inline void clear_ppairMatrix (PPairElement ** ppMatrix);
inline void initialize_proband_tmpLikelihood (Person * pPerson);
void populate_pedigree_loopbreaker_genotype_vector (Pedigree * pPed);
void populate_loopbreaker_genotype_vector (Person * pLoopBreaker, int locus);
int set_next_loopbreaker_genotype_vector (Pedigree * pPed, int initialFlag);


/*
 * before likelihood calculation, pre-allocate space to store conditional
 * likelihoods 
 * numLocus - number of loci we analyze at a time
 */
int
allocate_likelihood_space (PedigreeSet * pPedigreeList, int numLocus1)
{
  Pedigree *pPedigree;
  int i;
  int size;

  numLocus = numLocus1;
  /* this is for loop breaker multilocus genotypes */
  pTempGenoVector = (Genotype **) malloc (sizeof (Genotype *) * numLocus);

  for (i = 0; i < pPedigreeList->numPedigree; i++) {
    pPedigree = pPedigreeList->ppPedigreeSet[i];
    /*
     * allocate conditional likelihood storage for each
     * person/pedigree
     */
    allocate_multi_locus_genotype_storage (pPedigree, numLocus);
  }

  /*
   * allocate storage for temporarily stored likelihood for similar
   * parental pairs (only with phase differences either likelihood
   * itself will be stored there or a pointer to likelihood polynomial
   * will be
   */
  size = pow (2, numLocus);
  ppairMatrixNumLocus = numLocus;
  ppairMatrixRowSize = size * sizeof (PPairElement);
  ppairMatrix = (PPairElement **) calloc (size, sizeof (PPairElement *));
  for (i = 0; i < size; i++) {
    ppairMatrix[i] = (PPairElement *) calloc (size, sizeof (PPairElement));
  }
  if (bitMask != NULL)
    free (bitMask);
  bitMask = malloc (sizeof (int) * (numLocus + 1));
  for (i = 0; i <= numLocus; i++) {
    bitMask[i] = pow (2, i) - 1;
  }

  /* likelihood work space */
  pHaplo = &parentalPairSpace;

  /* pre allocate child likelihood elements */
  maxChildElements = 1024;
  likelihoodChildElements =
    (ChildElement *) calloc (sizeof (ChildElement), maxChildElements);

  maxChildren = 20;
  likelihoodChildCount = (int *) calloc (sizeof (int), maxChildren);

  return 0;
}

/* free the storage space for conditionals */
void
free_likelihood_space (PedigreeSet * pPedigreeList)
{
  Pedigree *pPedigree;
  int i;

  for (i = 0; i < pPedigreeList->numPedigree; i++) {
    pPedigree = pPedigreeList->ppPedigreeSet[i];
    free_multi_locus_genotype_storage (pPedigree);
  }

  /* free storage for temporary likelihood for similar parental pairs */
  for (i = 0; i < pow (2, locusList->numLocus); i++) {
    free (ppairMatrix[i]);
  }
  free (ppairMatrix);
  ppairMatrix = NULL;
  free (bitMask);
  bitMask = NULL;
  free (pTempGenoVector);
  pTempGenoVector = NULL;
  free (likelihoodChildElements);
  free (likelihoodChildCount);
  likelihoodChildElements = NULL;
  likelihoodChildCount = NULL;
}

/* the main API to compute likelihood for all pedigrees in a data set */
int
compute_likelihood (PedigreeSet * pPedigreeList)
{
  Pedigree *pPedigree;
  int i;
  int status;			/* return status of function calls */
  double product_likelihood = 1;	/* product of the likelihoods
					 * for all the pedigrees */
  double sum_log_likelihood = 0;	/* sum of the
					 * log10(likelihood) for all
					 * the pedigrees */
  double log10Likelihood;
  int origLocus = 0;		/* locus index in the original locus list
				 * this is used to find out the pedigree
				 * counts mainly for case control analyses */
  clock_t time2;

  if (locusList->numLocus > 1)
    origLocus = locusList->pLocusIndex[1];
  numLocus = locusList->numLocus;

  /* initialization */
  sum_log_likelihood = 0;
  product_likelihood = 1;
  pPedigreeList->likelihood = 1;
  pPedigreeList->log10Likelihood = 0;

  /* loop over pedigrees in the data set */
  for (i = 0; i < pPedigreeList->numPedigree; i++) {
    pPedigree = pPedigreeList->ppPedigreeSet[i];

    if (pPedigree->load_flag == 0) {

#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE) {
	if (pPedigree->likelihoodPolynomial == NULL) {
	  /*
	   * only build likelihood polynomial once, if
	   * the ptr is not NULL, it means the
	   * polynomial has been constructed
	   */
	  //fprintf(stderr, "The polynomial building for this pedigree should be only once\n");
	  /*
	   * initialize likelihood space for each
	   * pedigree
	   */
	  initialize_multi_locus_genotype (pPedigree);
	  //fprintf(stderr, "Start polynomial building\n");
	  status = compute_pedigree_likelihood (pPedigree);
	  //fprintf(stderr, "holdPoly for the pedigree polynomial\n");
	  //expTermPrinting(stderr, pPedigree->likelihoodPolynomial, 1);
	  //fprintf(stderr, "\n");
	  holdPoly (pPedigree->likelihoodPolynomial);
	  //fprintf(stderr, "freeKeptPolys after likelihood build and hold for pedigree\n");
	  freeKeptPolys ();
	  //              printAllPolynomials();
	  pPedigree->likelihoodPolyList = buildPolyList ();
	  polyListSorting (pPedigree->likelihoodPolynomial,
			   pPedigree->likelihoodPolyList);
	  if (i == pPedigreeList->numPedigree - 1) {
	    time2 = clock ();
	    fprintf (stderr, "Finished polynomial building: %f\n",
		     (double) time2 / CLOCKS_PER_SEC);
	  }
	}
	/* evaluate likelihood */
	evaluatePoly (pPedigree->likelihoodPolynomial,
		      pPedigree->likelihoodPolyList, &pPedigree->likelihood);
      } else {
	initialize_multi_locus_genotype (pPedigree);
	status = compute_pedigree_likelihood (pPedigree);
      }
#else
      initialize_multi_locus_genotype (pPedigree);
      status = compute_pedigree_likelihood (pPedigree);
#endif

      if (modelOptions.dryRun == 0) {
	if (pPedigree->likelihood == 0.0) {
	  KLOG (LOGLIKELIHOOD, LOGWARNING,
		"Pedigree %s has likelihood of 0 or too small.\n",
		pPedigree->sPedigreeID);
	  fprintf (stderr,
		   "Pedigree %s has likelihood of 0 or too small.\n",
		   pPedigree->sPedigreeID);
	  product_likelihood = 0.0;
	  sum_log_likelihood = -9999.99;
	  break;
	} else if (pPedigree->likelihood < 0.0) {
	  KASSERT (pPedigree->likelihood >= 0.0,
		   "Pedigree %s with NEGATIVE likelihood - This is CRAZY!!!.\n",
		   pPedigree->sPedigreeID);
	  product_likelihood = 0.0;
	  sum_log_likelihood = -9999.99;
	  break;
	} else {
	  if (pPedigree->pCount[origLocus] == 1) {
	    product_likelihood *= pPedigree->likelihood;
	    log10Likelihood = log10 (pPedigree->likelihood);
	  } else {
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

  if (modelOptions.polynomial == TRUE) {
    for (i = 0; i < pPedigreeList->numPedigree; i++) {
      pPedigree = pPedigreeList->ppPedigreeSet[i];
      if (pPedigree->likelihoodPolynomial != NULL) {
	unHoldPoly (pPedigree->likelihoodPolynomial);
	pPedigree->likelihoodPolynomial = NULL;
	free (pPedigree->likelihoodPolyList->pList);
	free (pPedigree->likelihoodPolyList);
	pPedigree->likelihoodPolyList = NULL;
      }
    }
    freeKeptPolys ();		/* Because holds overlapped keeps. */
    //    fprintf (stderr, "Post-pedigree free minimum:\n");
    //    polyStatistics ();
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
  double likelihood = 0;
  double tmpLikelihood = 0;
  ConditionalLikelihood *pConditional;

#ifndef NO_POLYNOMIAL
  Polynomial *pLikelihoodPolynomial = NULL;
#endif
  int ret = 0;

#if 1
  if (modelOptions.dryRun == 0) {
    fprintf (stderr, "PEDIGREE: %s (%d/%d)\n",
	     pPedigree->sPedigreeID, pPedigree->pedigreeIndex + 1,
	     pPedigree->pPedigreeSet->numPedigree);
  }
#endif

  for (i = 0; i < pPedigree->numNuclearFamily; i++) {
    pNucFam = pPedigree->ppNuclearFamilyList[i];
    pNucFam->totalNumPairGroups = 0;
    pNucFam->totalNumSimilarPairs = 0;
  }

  if (pPedigree->loopFlag == TRUE) {
    populate_pedigree_loopbreaker_genotype_vector (pPedigree);
    ret = -2;
    while (ret == -2) {
      ret = set_next_loopbreaker_genotype_vector (pPedigree, TRUE);
      if (ret == -2)
	restore_pedigree_genotype_link_from_saved (pPedigree);
    }
  }
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE) {
    pLikelihoodPolynomial = constant0Poly;
  } else
    likelihood = 0;
#else
  likelihood = 0;
#endif

  while (ret == 0) {
    initialize_multi_locus_genotype (pPedigree);
    /*
     * initialize all the nuclear families before peeling starts
     * in many cases, multiple likelihoods are computed for the
     * same family with different parameters, we need to clean up
     * before (or after) each calculation
     */
    for (i = 0; i < pPedigree->numNuclearFamily; i++) {
      pNucFam = pPedigree->ppNuclearFamilyList[i];
      pNucFam->doneFlag = 0;
      pNucFam->numPairGroups = 0;
      pNucFam->numSimilarPairs = 0;
    }

    /*
     * peeling starts from the peeling proband and eventually
     * will come back to the same proband this process will
     * obtain the conditional likelihoods for the proband
     */
    status = peel_graph (pPedigree->pPeelingNuclearFamily,
			 pPedigree->pPeelingProband,
			 pPedigree->peelingDirection);

    /*
     * done peeling, need to add up the conditional likelihood
     * for the leading peeling proband
     */
    pProband = pPedigree->pPeelingProband;
#ifndef NO_POLYNOMIAL
    if (modelOptions.polynomial != TRUE)
      tmpLikelihood = 0;
#else
    tmpLikelihood = 0;
#endif

    /* loop over all conditional likelihoods */
    for (i = 0; i < pProband->numConditionals; i++) {
      pConditional = &pProband->pLikelihood[i];
      if (pConditional->touchedFlag == 0)
	continue;
      /*
       * Get the joint likelihood = marginal * p(G) when
       * the proband is a founder, the weight will be the
       * multilocus genotype probabilities under LE or
       * haplotype frequencies under LD otherwise the
       * weight should be 1
       *
       */
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE) {
	/* build likelihood polynomial */
	pLikelihoodPolynomial =
	  plusExp (2,
		   1.0,
		   pLikelihoodPolynomial,
		   1.0,
		   timesExp (2,
			     pConditional->lkslot.
			     likelihoodPolynomial, 1,
			     pConditional->wtslot.weightPolynomial, 1, 0), 1);
      } else {
	tmpLikelihood += pConditional->lkslot.likelihood *
	  pConditional->wtslot.weight;
      }
#else
      tmpLikelihood += pConditional->lkslot.likelihood *
	pConditional->wtslot.weight;
#endif
    }				/* end of looping over all conditionals */


#ifndef NO_POLYNOMIAL
    if (modelOptions.polynomial != TRUE)
      likelihood += tmpLikelihood;
#else
    likelihood += tmpLikelihood;
#endif

    if (pPedigree->loopFlag == TRUE) {
      if (modelOptions.polynomial == TRUE)
	keepPoly (pLikelihoodPolynomial);
      else
	KLOG (LOGLIKELIHOOD, LOGDEBUG,
	      "Log Likelihood for this fixed looped pedigree %s is: %e\n",
	      pPedigree->sPedigreeID, log10 (tmpLikelihood));
      ret = -2;
      while (ret == -2) {
	restore_pedigree_genotype_link_from_saved (pPedigree);
	ret = set_next_loopbreaker_genotype_vector (pPedigree, FALSE);
      }
    } else {
      ret = -1;
    }

    for (i = 0; i < pPedigree->numNuclearFamily; i++) {
      pNucFam = pPedigree->ppNuclearFamilyList[i];
      pNucFam->totalNumPairGroups += pNucFam->numPairGroups;
      pNucFam->totalNumSimilarPairs += pNucFam->numSimilarPairs;
    }
  }
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    /* save the polynomial to the pedigree structure */
    pPedigree->likelihoodPolynomial = pLikelihoodPolynomial;
  else {
    /* save the likelihood in the pedigree structure */
    pPedigree->likelihood = likelihood;
    KLOG (LOGLIKELIHOOD, LOGDEBUG,
	  "log Likelihood for pedigree %s is: %e\n", pPedigree->sPedigreeID,
	  log10 (likelihood));
  }
#else
  pPedigree->likelihood = likelihood;
  KLOG (LOGLIKELIHOOD, LOGDEBUG, "log Likelihood for pedigree %s is: %e\n",
	pPedigree->sPedigreeID, log10 (likelihood));
#endif

  return 0;
}

/*
 * recursive procedure to go through all nuclear families within one pedigree
 * pNucFam -- input nuclear family, the top layer is the peeling nuclear
 * family which is the nuclear family contains the peeling proband pProband
 * -- connector peelingDireciton -- UP/DOWN -- currently it is not used at
 * all
 */
int
peel_graph (NuclearFamily * pNucFam1, Person * pProband1,
	    int peelingDirection)
{
  NuclearFamilyConnector *pConnector;	/* connector structure which
					 * represents connector from one nuc
					 * to the other within a pedigree */
  NuclearFamily *pNucFam2;	/* another nuc family input nuc family is
				 * connected to */
  Person *pConnectPerson;	/* connector individual */
  int i;
  ConditionalLikelihood *pConditional;

  /*
   * if this nuclear family has been processed or in the middle of
   * process, skip it
   */
  if (pNucFam1->doneFlag == TRUE)
    return 0;

  /*
   * mark this nuclear family as done to avoid potential endless
   * recurisve calls
   */
  pNucFam1->doneFlag = TRUE;

  /* go up through the connectors if any */
  pConnector = pNucFam1->pUpConnectors;
  while (pConnector) {
    pNucFam2 = pConnector->pConnectedNuclearFamily;
    pConnectPerson = pConnector->pConnectedPerson;
    if (pConnectPerson == pNucFam2->pParents[DAD] ||
	pConnectPerson == pNucFam2->pParents[MOM]) {
      /*
       * these two families are connected through multiple
       * marraige
       */
      peel_graph (pConnector->pConnectedNuclearFamily,
		  pConnector->pConnectedPerson, PEDIGREE_UP);
    } else {
      peel_graph (pConnector->pConnectedNuclearFamily,
		  pConnector->pConnectedPerson, PEDIGREE_DOWN);
    }
    pConnector = pConnector->pNextConnector;
  }

  /* go down through the down connectors if any */
  pConnector = pNucFam1->pDownConnectors;
  while (pConnector) {
    /* peel up to the proband */
    peel_graph (pConnector->pConnectedNuclearFamily,
		pConnector->pConnectedPerson, PEDIGREE_UP);

    pConnector = pConnector->pNextConnector;
  }

  /*
   * we are done with the up or down linked nuclear families or we are
   * a leave or top. ready for likelihood computation for this nuclear
   * family save the original genotype list first during likelihood
   * calculation, we limit the child proband's genotype to one at a
   * time once done, the original list will be copied back loop
   * breaker's duplicate can't be a proband as the duplicate is the one
   * doesn't have parents, so it can't be a connector and can't be a
   * proband
   */
  pProband = pProband1;
  pNucFam = pNucFam1;
  pNucFam->pProband = pProband;

  memcpy (&pProband->ppProbandGenotypeList[0],
	  &pProband->ppGenotypeList[0],
	  sizeof (Genotype *) * originalLocusList.numLocus);
  memcpy (&pProband->pProbandNumGenotype[0],
	  &pProband->pNumGenotype[0],
	  sizeof (int) * originalLocusList.numLocus);

  KLOG (LOGPARENTALPAIR, LOGDEBUG, "\t Proband (%s) haplotype: \n",
	pProband->sID);
  if (pProband->ppHaplotype == NULL) {
    /*
     * allocate space for storing proband's haplotype if not
     * already done
     */
    pProband->ppHaplotype = MALLOC ("pProband->ppHaplotype",
				    sizeof (Genotype *) * sizeof (int) *
				    ppairMatrixNumLocus);
  }
  if (pProband->touchedFlag == TRUE)
    initialize_proband_tmpLikelihood (pProband);

  if (pNucFam->pParents[DAD] != pProband
      && pNucFam->pParents[MOM] != pProband) {
    /* proband is a child */
    pNucFam->childProbandFlag = TRUE;
  } else {
    /* proband is a parent */
    pNucFam->childProbandFlag = FALSE;
  }
  if (pNucFam->childProbandFlag == TRUE) {
    /*
     * proband is a child construct all possible multilocus
     * genotype for the proband for each of them, calculate the
     * conditional likelihood of the nuclear family
     */
    loop_child_proband_genotype (peelingDirection, 0, 0);
  } else {			/* A parent is the proband */
    compute_nuclear_family_likelihood (peelingDirection);
    if (modelOptions.dryRun == 0) {
      /* copy the touched temporary results back */
      for (i = 0; i < pProband->numTmpLikelihood; i++) {
	pConditional =
	  &pProband->pLikelihood[pProband->pTmpLikelihoodIndex[i]];
	//      if (pConditional->tmpTouched == FALSE)
	//        continue;
	pConditional->touchedFlag = TRUE;
#ifndef NO_POLYNOMIAL
	if (modelOptions.polynomial == TRUE) {
	  if (pProband->touchedFlag == FALSE) {
	    pConditional->lkslot.likelihoodPolynomial = constant1Poly;
	  }
	  pConditional->lkslot.likelihoodPolynomial =
	    timesExp (2,
		      pConditional->lkslot.likelihoodPolynomial, 1,
		      pConditional->tmpslot.tmpLikelihoodPolynomial, 1, 0);
	  pConditional->tmpslot.tmpLikelihoodPolynomial = constant0Poly;
#if 0
	  KLOG (LOGLIKELIHOOD, LOGDEBUG,
		"Proband %s Conditional Likelihood (%d) = %e. Weight = %e\n",
		pProband->sID, i,
		evaluateValue (pConditional->lkslot.
			       likelihoodPolynomial),
		evaluateValue (pConditional->wtslot.weightPolynomial));
#endif
	} else {		/* PE is not enabled */
	  if (pProband->touchedFlag == FALSE)
	    pConditional->lkslot.likelihood = 1;
	  pConditional->lkslot.likelihood *=
	    pConditional->tmpslot.tmpLikelihood;
	  pConditional->tmpslot.tmpLikelihood = 0;
	  KLOG (LOGLIKELIHOOD, LOGDEBUG,
		"Proband %s Conditional Likelihood (%d) = %e. Weight = %e\n",
		pProband->sID, i,
		pConditional->lkslot.likelihood, pConditional->wtslot.weight);
	}
#else
	if (pProband->touchedFlag == FALSE)
	  pConditional->lkslot.likelihood = 1;
	pConditional->lkslot.likelihood *=
	  pConditional->tmpslot.tmpLikelihood;
	pConditional->tmpslot.tmpLikelihood = 0;
	KLOG (LOGLIKELIHOOD, DEBUG,
	      "Proband %s Conditional Likelihood (%d) = %e. Weight = %e \n",
	      pProband->sID, pProband->pTmpLikelihoodIndex[i],
	      pConditional->lkslot.likelihood, pConditional->wtslot.weight);
#endif
	pConditional->tmpTouched = FALSE;
      }
      pProband->numTmpLikelihood = 0;
    }
  }

  KLOG (LOGLIKELIHOOD, LOGDEBUG,
	"Nuclear Family %d with parents %s x %s.\n",
	pNucFam->nuclearFamilyIndex, pNucFam->pParents[DAD]->sID,
	pNucFam->pParents[MOM]->sID);

  /*
   * mark the proband as been touched - we have done some likelihood
   * calculation on this person
   */
  pProband->touchedFlag = TRUE;

  /* copy back the genotypes for the proband */
  memcpy (&pProband->ppGenotypeList[0],
	  &pProband->ppProbandGenotypeList[0],
	  sizeof (Genotype *) * originalLocusList.numLocus);
  memcpy (&pProband->pNumGenotype[0],
	  &pProband->pProbandNumGenotype[0],
	  sizeof (int) * originalLocusList.numLocus);

#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE) {
    //fprintf(stderr, "keepPoly for the likelihood and weight polynomials\n");
    for (i = 0; i < pProband->numConditionals; i++) {
      if (pProband->touchedFlag == TRUE) {
	keepPoly (pProband->pLikelihood[i].lkslot.likelihoodPolynomial);
	keepPoly (pProband->pLikelihood[i].wtslot.weightPolynomial);
      }
    }
    //fprintf(stderr, "freePolys after keeping the likelihood and weight polynomials\n");
    freePolys ();
  }
#endif

  return 0;
}

/*
 * This function loops over child proband's multilocus genotypes recursively
 * and calls compute_nuclear_family_likelihood for each multilocus genotype
 * pNucFam - the nuclear family we are working on locus - current working
 * locus multiLocusIndex - index to the flattened array of conditional
 * likelihood for the proband
 */
int
loop_child_proband_genotype (int peelingDirection,
			     int locus, int multiLocusIndex)
{
  int origLocus = locusList->pLocusIndex[locus];	/* locus index in the
							 * original locus list */
  Genotype *pGenotype;
  int position;			/* genotype position */
  Genotype *pNextGenotype;
  int numGenotype;		/* number of possible genotypes for this
				 * proband at this locus */
  int multiLocusIndex2;
  int traitLocus;
  ConditionalLikelihood *pConditional;
  double penetrance = 1;

#ifndef NO_POLYNOMIAL
  Polynomial *penetrancePolynomial = NULL;
#endif

  /*
   * we loop over the genotypes of the proband to condition the
   * likelihood calculation on it
   */
  numGenotype = pProband->pSavedNumGenotype[origLocus];
  /* calculate the flattened conditional likelihood array index */
  multiLocusIndex2 = multiLocusIndex * numGenotype;
  pGenotype = pProband->ppProbandGenotypeList[origLocus];
  while (pGenotype != NULL) {
    /*
     * record this locus's genotype in the haplotype structure -
     * it's really just phased multilocus genotype (not single
     * chromosome haplotype)
     */
    pProband->ppHaplotype[locus] = pGenotype;
    KLOG (LOGLIKELIHOOD, LOGDEBUG, "\t proband (%s) %d|%d \n",
	  pProband->sID, pGenotype->allele[DAD], pGenotype->allele[MOM]);
    /*
     * temporarilly set the next pointer to NULL so to restrict
     * the genotype on the proband to current genotype only
     */
    pProband->ppGenotypeList[origLocus] = pGenotype;
    pNextGenotype = pGenotype->pNext;
    pGenotype->pNext = NULL;
    pProband->pNumGenotype[origLocus] = 1;
    position = pGenotype->position;
    /* calculate the flattened conditional likelihood array index */
    multiLocusIndex = multiLocusIndex2 + position;

    if (locus < locusList->numLocus - 1) {
      loop_child_proband_genotype (peelingDirection, locus + 1,
				   multiLocusIndex);
    } else {
      /*
       * we have got the entire multi-locus genotype for
       * the proband
       */
      compute_nuclear_family_likelihood (peelingDirection);
      pConditional = &pProband->pLikelihood[multiLocusIndex];
      pConditional->touchedFlag = TRUE;
      if (modelOptions.dryRun == 0) {
	/*
	 * store the likelihood in the corresponding
	 * flattened array
	 */
	if (pProband->touchedFlag == FALSE) {
	  /*
	   * if trait locus exists, we need to retrieve
	   * the penetrance
	   */
	  traitLocus = locusList->traitLocusIndex;
	  if (traitLocus >= 0) {
#ifndef NO_POLYNOMIAL
	    if (modelOptions.polynomial == TRUE) {
	      penetrancePolynomial =
		pProband->ppHaplotype[traitLocus]->penslot.
		penetrancePolynomial;
	    } else {
	      penetrance =
		pProband->ppHaplotype[traitLocus]->penslot.penetrance;
	    }
#else
	    penetrance =
	      pProband->ppHaplotype[traitLocus]->penslot.penetrance;
#endif
	  } else {		/* no trait locus */
#ifndef NO_POLYNOMIAL
	    if (modelOptions.polynomial == TRUE) {
	      penetrancePolynomial = constant1Poly;
	    } else {
	      penetrance = 1;
	    }
#else
	    penetrance = 1;
#endif
	  }			/* end of trait locus processing */

#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE) {
	    pConditional->lkslot.likelihoodPolynomial =
	      timesExp (2, penetrancePolynomial, 1,
			pNucFam->likelihoodPolynomial, 1, 0);

	    /*
	     * for a child, the weight should be
	     * 1
	     */
	    pConditional->wtslot.weightPolynomial = constant1Poly;
	  } else {
	    /*
	     * need to update the penetrance
	     * factors
	     */
	    pConditional->lkslot.likelihood =
	      penetrance * pNucFam->likelihood;
	    /*
	     * for a child, the weight should be
	     * 1
	     */
	    pConditional->wtslot.weight = 1;
	  }
#else
	  /* need to update the penetrance factors */
	  pConditional->lkslot.likelihood = penetrance * pNucFam->likelihood;
	  /* for a child, the weight should be 1 */
	  pConditional->wtslot.weight = 1;
#endif
	}
	/*
	 * first time updating the likelihood for this phased
	 * multilocus genotype
	 */
	else {			/* NOT first time updating the likelihood for
				 * this phased multilocus genotypes */
	  /* no need to consider penetrance anymore */
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE) {
	    pConditional->lkslot.
	      likelihoodPolynomial =
	      timesExp (2,
			pConditional->lkslot.likelihoodPolynomial, 1,
			pNucFam->likelihoodPolynomial, 1, 1);
	    //fprintf(stderr, "Likelihood for this entire multi-locus genotype %f %f\n",
	    //evaluateValue(pNucFam->likelihoodPolynomial),
	    //evaluateValue(pProband->pLikelihood[multiLocusIndex].likelihoodPolynomial));
#if 0
	    KLOG (LOGLIKELIHOOD, LOGDEBUG,
		  "Proband %s Conditional Likelihood (%d) = %e.\n",
		  pProband->sID, multiLocusIndex,
		  evaluateValue (pConditional->lkslot.likelihoodPolynomial));
#endif
	  } else {
	    pConditional->lkslot.likelihood *= pNucFam->likelihood;
	    KLOG (LOGLIKELIHOOD, LOGDEBUG,
		  "Proband %s Conditional Likelihood (%d) = %e.\n",
		  pProband->sID, multiLocusIndex,
		  pConditional->lkslot.likelihood);
	  }
#else
	  pConditional->lkslot.likelihood *= pNucFam->likelihood;
	  KLOG (LOGLIKELIHOOD, LOGDEBUG,
		"Proband %s Conditional Likelihood (%d) = %e.\n",
		pProband->sID, multiLocusIndex,
		pConditional->lkslot.likelihood);
#endif
	}
      }
    }				/* end of processing complete child
				 * multilocus genotype  */

    /*
     * when we are done, need to restore the genotype link list
     * pointer
     */
    pGenotype->pNext = pNextGenotype;
    pGenotype = pNextGenotype;
  }				/* loop over all possible genotypes */

  return 0;
}

/*
 * for now, we only use parental pair algorithm This function goes through
 * all possible parental pairs, for each parental pair, it calculates
 * likelihood
 */
int
compute_nuclear_family_likelihood (int peelingDirection)
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
  int multiLocusIndex[2] = { 0, 0 };

  /* initialize the weight for each parent */
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE) {
    weightPolynomial[0] = constant1Poly;
    weightPolynomial[1] = constant1Poly;
    pNucFam->likelihoodPolynomial = constant0Poly;
  }
#endif
  pNucFam->likelihood = 0;

  /*
   * the following is to help to set the order which parent's genotype
   * to flip first
   */
  if (pProband == pNucFam->pParents[MOM]) {
    /* MOM is the proband */
    pNucFam->head = MOM;
    pNucFam->spouse = DAD;
  } else {
    /* either DAD is the proband or the child is the proband */
    pNucFam->head = DAD;
    pNucFam->spouse = MOM;
  }

  numChild = pNucFam->numChildren;

  /*
   * first construct the parental pair for this nuclear family locus by
   * locus
   */
  numHaplotypePair = 1;
  for (locus = 0; locus < locusList->numLocus; locus++) {
    /* construct parental pair locus by locus */
    construct_parental_pair (pNucFam, pProband, locus);
    /* calculate number of possible multilocus genotypes */
    numHaplotypePair *= parentalPairSpace.pNumParentalPair[locus];
  }

  for (i = DAD; i <= MOM; i++) {
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
  /*
   * recursively call loop_parental_pair to get a complete multlocus
   * genotype
   */
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE) {
    loop_parental_pair (0, multiLocusIndex, (void *) weightPolynomial);
    //fprintf(stderr, "Conditional likelihood for nuclear family %d is: %e\n",
    //pNucFam->nuclearFamilyIndex, evaluateValue(sumPolynomial));
  } else {
    loop_parental_pair (0, multiLocusIndex, (void *) weight);
    /*
       KLOG (LOGLIKELIHOOD, LOGDEBUG,
       "Conditional likelihood for nuclear family %d is: %e\n",
       pNucFam->nuclearFamilyIndex, sum);
     */
  }
#else
  loop_parental_pair (0, multiLocusIndex, (void *) weight);
  /*
     KLOG (LOGLIKELIHOOD, LOGDEBUG,
     "Conditional likelihood for nuclear family %d is: %e\n",
     pNucFam->nuclearFamilyIndex, sum);
   */
#endif

  return 0;
}

/*
 * using the parental pair algorithm to go through all possible parental
 * pairs and calculate the likelihood of each parental pair nested looping is
 * for the multi-locus
 */
int
loop_parental_pair (int locus, int multiLocusIndex[2], void *dWeight[2])
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
  int head;			/* proband if no child is a proband,
				 * otherwise DAD  */
  int spouse;			/* spouse of the head */
  int likelihoodIndex;		/* likelihood index for the proband */
  int multiLocusPhase[2] = { 0, 0 };	/* index to the related
					 * parental pair matrix */
  int flipMask[2] = { 0, 0 };
  ConditionalLikelihood *pConditional;
  PPairElement *pElement;


  head = pNucFam->head;
  spouse = pNucFam->spouse;

  origLocus = locusList->pLocusIndex[locus];
  for (i = DAD; i <= MOM; i++) {
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
    if (pParent[i]->loopBreaker >= 1 && pParent[i]->pParents[DAD] == NULL)
      numGenotype[i] =
	pParent[i]->pOriginalPerson->pSavedNumGenotype[origLocus];
    else
      numGenotype[i] = pParent[i]->pSavedNumGenotype[origLocus];
    /* calculate this parent's conditional likelihood's offset */
    multiLocusIndex2[i] = multiLocusIndex[i] * numGenotype[i];
  }
  /* parental pair index for this locus */
  numPair = -1;
  while ((numPair + 1) < pHaplo->pNumParentalPair[locus]) {
    /* find related parental pairs */
    numPair++;
    /* get the parental pair */
    pPair = &pHaplo->ppParentalPair[locus][numPair];
    multiLocusIndex2[DAD] =
      multiLocusIndex[DAD] * numGenotype[DAD] +
      pPair->pGenotype[DAD]->position;
    multiLocusIndex2[MOM] =
      multiLocusIndex[MOM] * numGenotype[MOM] +
      pPair->pGenotype[MOM]->position;
    pHaplo->pParentalPairInd[locus] = numPair;
    /* set the het flag */
    for (i = DAD; i <= MOM; i++) {
      if (locus == 0) {
	/*
	 * keep track of how many heterogenous loci
	 * there are at and before this locus
	 */
	pNucFam->tmpNumHet[i][0] = 0;
      } else {
	pNucFam->tmpNumHet[i][locus] = pNucFam->tmpNumHet[i][locus - 1];
      }
      if (pNucFam->firstHetLocus[i] >= locus) {
	/* initialize  */
	pNucFam->firstHetLocus[i] = -1;
      }
      if (isHet (pPair->pGenotype[i])) {
	pNucFam->hetFlag[i][locus] = 1;
	if (pNucFam->firstHetLocus[i] == -1)
	  pNucFam->firstHetLocus[i] = locus;
	pNucFam->tmpNumHet[i][locus]++;
      } else
	pNucFam->hetFlag[i][locus] = 0;
    }
    /*
     * record the start of the related parental pairs for this
     * locus
     */
    pNucFam->relatedPPairStart[locus] = numPair;
    pNucFam->numRelatedPPair[locus] = 1;
    /*
     * find the related parental pairs that have same pair of
     * genotypes only differ in phases
     */
    while ((numPair + 1) < pHaplo->pNumParentalPair[locus] &&
	   (pHaplo->ppParentalPair[locus][numPair + 1].phase[DAD] > 0 ||
	    pHaplo->ppParentalPair[locus][numPair + 1].phase[MOM] > 0)) {
      numPair++;
      pNucFam->numRelatedPPair[locus]++;
    }
    if (locus == 0) {
      pNucFam->totalRelatedPPair[locus] = pNucFam->numRelatedPPair[0];
    } else {
      pNucFam->totalRelatedPPair[locus] =
	pNucFam->totalRelatedPPair[locus -
				   1] * pNucFam->numRelatedPPair[locus];
    }

    for (i = DAD; i <= MOM; i++) {
      if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM
	  && pParent[i]->pParents[DAD] == NULL) {
	if (pParent[i]->loopBreaker >= 1)
	  continue;
	/*
	 * For founders: under LE, we just multiply
	 * the founder weights, under LD, we need to
	 * use haplotype freq - which we will do in
	 * loop_phases For non-founders, we don't
	 * even care, they should remain as 1 as
	 * initialized
	 */
#ifndef NO_POLYNOMIAL
	if (modelOptions.polynomial == TRUE) {
	  newWeightPolynomial[i] =
	    timesExp (2, (Polynomial *) dWeight[i], 1,
		      pPair->pGenotype[i]->wtslot.weightPolynomial, 1, 0);
	} else {
	  newWeight[i] =
	    *((double *) dWeight + i) * pPair->pGenotype[i]->wtslot.weight;
	}
#else
	newWeight[i] =
	  *((double *) dWeight + i) * pPair->pGenotype[i]->wtslot.weight;
#endif
      }
    }				/* looping dad and mom genotypes */


    if (locus < locusList->numLocus - 1) {
      /*
       * recursively calling this function to get a
       * complete multilocus genotype
       */
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE) {

	loop_parental_pair (locus + 1, multiLocusIndex2,
			    (void *) newWeightPolynomial);
      } else {
	loop_parental_pair (locus + 1, multiLocusIndex2, (void *) newWeight);
      }
#else
      loop_parental_pair (locus + 1, multiLocusIndex2, (void *) newWeight);
#endif
    } else {			/* got a complete set of parental pairs */
      pNucFam->numPairGroups++;
      pNucFam->numSimilarPairs += pNucFam->totalRelatedPPair[locus] - 1;
      if (modelOptions.dryRun != 0) {
	/* no actual calculations under dry run, continue to next set of ppair */
	continue;
      }

      if (pNucFam->totalRelatedPPair[locus] == 1) {
	multiLocusPhase[DAD] = 0;
	multiLocusPhase[MOM] = 0;
	calcFlag = 0;
#ifndef NO_POLYNOMIAL
	if (modelOptions.polynomial == TRUE)
	  calculate_likelihood (multiLocusIndex2, multiLocusPhase,
				(void *) newWeightPolynomial, NULL);
	else
	  calculate_likelihood (multiLocusIndex2, multiLocusPhase,
				(void *) newWeight, NULL);
#else
	calculate_likelihood (multiLocusIndex2, multiLocusPhase,
			      (void *) newWeight, NULL);
#endif

	/* when there is only one parental pair, the likelihood is saved at the first cell,
	 * as the phase index passed in was 0,0
	 */
	pElement = &ppairMatrix[0][0];
	if (pNucFam->childProbandFlag == TRUE) {
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE) {
	    pNucFam->likelihoodPolynomial =
	      plusExp (2,
		       1.0,
		       pNucFam->likelihoodPolynomial,
		       1.0, pElement->slot.likelihoodPolynomial, 1);
	  } else
	    pNucFam->likelihood += pElement->slot.likelihood;
#else
	  pNucFam->likelihood += pElement->slot.likelihood;
#endif
	} else {		/* one of the parent is the proband */
	  likelihoodIndex = multiLocusIndex2[head];
	  pConditional = &pProband->pLikelihood[likelihoodIndex];
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE) {
	    pConditional->tmpslot.tmpLikelihoodPolynomial =
	      plusExp (2,
		       1.0,
		       pConditional->tmpslot.
		       tmpLikelihoodPolynomial, 1.0,
		       pElement->slot.likelihoodPolynomial, 1);
	  } else
	    pConditional->tmpslot.tmpLikelihood += pElement->slot.likelihood;
#else
	  pConditional->tmpslot.tmpLikelihood += pElement->slot.likelihood;
#endif
	  if (pConditional->tmpTouched == FALSE) {
	    pConditional->tmpTouched = TRUE;
	    pProband->pTmpLikelihoodIndex[pProband->
					  numTmpLikelihood] = likelihoodIndex;
	    pProband->numTmpLikelihood++;
	  }
	}
	/* reset count */
	pElement->count = 0;
      } else if (pNucFam->totalRelatedPPair[locus] > 0) {
	/*
	 * initialize the phase matrix - we will
	 * count on count being reset after each use
	 * - this is to save time
	 */
	//clear_ppairMatrix(ppairMatrix);

	/* set some information */
	pNucFam->numHetLocus[head] = pNucFam->tmpNumHet[head][locus];
	pNucFam->numHetLocus[spouse] = pNucFam->tmpNumHet[spouse][locus];
	/*
	 * set bit mask - all bits set to 1 - number
	 * of bits == number of het loci
	 */
	pNucFam->hetLocusBits[head] = bitMask[pNucFam->numHetLocus[head]];
	pNucFam->hetLocusBits[spouse] = bitMask[pNucFam->numHetLocus[spouse]];

	multiLocusIndex2[DAD] = 0;
	multiLocusIndex2[MOM] = 0;
	multiLocusPhase[DAD] = 0;
	multiLocusPhase[MOM] = 0;
#ifndef NO_POLYNOMIAL
	if (modelOptions.polynomial == TRUE)
	  loop_phases (0, multiLocusIndex2, multiLocusPhase, flipMask,
		       (void *) newWeightPolynomial);
	else
	  loop_phases (0, multiLocusIndex2, multiLocusPhase, flipMask,
		       (void *) newWeight);
#else
	loop_phases (0, multiLocusIndex2, multiLocusPhase, flipMask,
		     (void *) newWeight);
#endif

	/*
	 * post processing of results of similar
	 * parental pairs and store them in the
	 * proband's likelihood space
	 */
	if (pNucFam->childProbandFlag == TRUE) {
	  /*
	   * child is the proband, sum
	   * likelihood across all rows and
	   * columns  save in the nuclear
	   * family
	   */
	  for (j = 0; j <= bitMask[pNucFam->numHetLocus[head]]; j++) {
	    for (k = 0; k <= bitMask[pNucFam->numHetLocus[spouse]]; k++) {
	      pElement = &ppairMatrix[j][k];
	      if (pElement->count > 1) {
#ifndef NO_POLYNOMIAL
		if (modelOptions.polynomial == TRUE) {
		  pNucFam->likelihoodPolynomial =
		    plusExp (2,
			     1.0,
			     pNucFam->likelihoodPolynomial,
			     1.0,
			     timesExp (2,
				       pElement->slot.
				       likelihoodPolynomial,
				       1,
				       constantExp
				       (pElement->count), 1, 0), 1);
		} else
		  pNucFam->likelihood +=
		    pElement->slot.likelihood * pElement->count;
#else
		pNucFam->likelihood +=
		  pElement->slot.likelihood * pElement->count;
#endif
	      } else if (pElement->count > 0) {	/* count == 1 */
#ifndef NO_POLYNOMIAL
		if (modelOptions.polynomial == TRUE) {
		  pNucFam->likelihoodPolynomial =
		    plusExp (2,
			     1.0,
			     pNucFam->likelihoodPolynomial,
			     1.0, pElement->slot.likelihoodPolynomial, 1);
		} else
		  pNucFam->likelihood += pElement->slot.likelihood;
#else
		pNucFam->likelihood += pElement->slot.likelihood;
#endif
	      }
	      /* end of count = 1 */
	      /* reset count to 0 */
	      pElement->count = 0;
	    }
	  }
	} else {		/* one of the parent is the proband */
	  for (j = 0; j <= bitMask[pNucFam->numHetLocus[head]]; j++) {
	    /*
	     * get the proband's
	     * conditional likelihood
	     * index for this row sum up
	     * across this row and save
	     * them to the proband's
	     * likelihood storage
	     */
	    likelihoodIndex = ppairMatrix[j][0].likelihoodIndex;
	    pConditional = &pProband->pLikelihood[likelihoodIndex];
	    for (k = 0; k <= bitMask[pNucFam->numHetLocus[spouse]]; k++) {
	      pElement = &ppairMatrix[j][k];
	      if (pElement->count == 0)
		continue;
	      if (pConditional->tmpTouched == FALSE) {
		pConditional->tmpTouched = TRUE;
		pProband->pTmpLikelihoodIndex[pProband->numTmpLikelihood]
		  = likelihoodIndex;
		pProband->numTmpLikelihood++;
	      }
	      if (pElement->count > 1) {
#ifndef NO_POLYNOMIAL
		if (modelOptions.polynomial == TRUE) {
		  pConditional->tmpslot.
		    tmpLikelihoodPolynomial =
		    plusExp (2, 1.0,
			     pConditional->tmpslot.
			     tmpLikelihoodPolynomial, 1.0,
			     timesExp (2,
				       pElement->slot.
				       likelihoodPolynomial,
				       1,
				       constantExp (pElement->
						    count), 1, 0), 1);
		} else
		  pConditional->tmpslot.tmpLikelihood +=
		    pElement->slot.likelihood * pElement->count;
#else
		pConditional->tmpslot.tmpLikelihood +=
		  pElement->slot.likelihood * pElement->count;
#endif
	      }
	      /* count > 1 */
	      else {		/* count==1 */
#ifndef NO_POLYNOMIAL
		if (modelOptions.polynomial == TRUE) {
		  pConditional->tmpslot.
		    tmpLikelihoodPolynomial =
		    plusExp (2, 1.0,
			     pConditional->tmpslot.
			     tmpLikelihoodPolynomial, 1.0,
			     pElement->slot.likelihoodPolynomial, 1);
		} else
		  pConditional->tmpslot.tmpLikelihood +=
		    pElement->slot.likelihood;
#else
		pConditional->tmpslot.tmpLikelihood +=
		  pElement->slot.likelihood;
#endif
	      }			/* count > 0 */
	      /* reset the count */
	      pElement->count = 0;
	    }			/* loop through column */
	  }			/* loop through row */
	}			/* one parent is proband */
      }
    }				/* end of processing one parental pair */
  }

  return 0;
}



/*
 * Compute likelihood for parental multilocus genotypes pairs that are only
 * different in phases Input: parental pairs for this nuclear family parental
 * pairs with different phases are next to each other weight - genotype
 * weights for the parental genotypes LE: founder - phase doesn't matter
 * non-founder - phase will come to play LD: phase matters for both founders
 * and non-founders penetrance - is not an input as it can be different for
 * different phases if we consider imprinting effect penetrance only applies
 * to disease locus if it's included in the list
 */
void
loop_phases (int locus, int multiLocusIndex[2], int multiLocusPhase[2],
	     int flipMask[2], void *dWeight[2])
{
  /* conditional likelihood index for each parent */
  int multiLocusIndex2[2];

  /* index to related parental pair matrix */
  int multiLocusPhase2[2];

  /*
   * index of the flip of current pattern in related parental pair
   * matrix
   */
  int multiLocusPhaseFlip[2];
  Person *pParent[2];
  int numGenotype[2];
  int i;
  int numPair;
  int end;
  int origLocus = locusList->pLocusIndex[locus];	/* locus index in the
							 * original locus list */
  ParentalPair *pPair;

  /*
   * whether a fresh calculation is needed or we could find existing
   * pattern result
   */
  int calculateFlag;
  int phase[2];
  int proband;
  int spouse;
  int likelihoodIndex;
  int newFlipMask[2];
  double childProduct;
  Polynomial *childProductPoly;
  void *childProductPtr;

  if (modelOptions.polynomial == TRUE)
    childProductPtr = &childProductPoly;
  else
    childProductPtr = &childProduct;

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

  for (i = DAD; i <= MOM; i++) {
    pParent[i] = pNucFam->pParents[i];
    /* find the max number of possible genotypes for this parent */
    if (pParent[i]->loopBreaker >= 1 && pParent[i]->pParents[DAD] == NULL)
      numGenotype[i] =
	pParent[i]->pOriginalPerson->pSavedNumGenotype[origLocus];
    else
      numGenotype[i] = pParent[i]->pSavedNumGenotype[origLocus];
    /* likelihood storage index */
    multiLocusIndex[i] *= numGenotype[i];
    if (pNucFam->hetFlag[i][locus] == 1) {
      /* phase combination index */
      multiLocusPhase[i] <<= 1;
    }
  }

  flipMask[DAD] <<= 2;
  flipMask[MOM] <<= 2;

  proband = pNucFam->head;
  spouse = pNucFam->spouse;
  /* get the start of the related pair */
  numPair = pNucFam->relatedPPairStart[locus] - 1;
  end = pNucFam->numRelatedPPair[locus] + pNucFam->relatedPPairStart[locus];
  while ((numPair += 1) < end) {
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
    for (i = DAD; i <= MOM; i++) {
      pHaplo->phase[i][locus] = pPair->phase[i];
      multiLocusPhase2[i] = multiLocusPhase[i] | pPair->phase[i];
      if (pPair->phase[i] == 0)
	newFlipMask[i] = flipMask[i];
      else
	newFlipMask[i] = flipMask[i] | 3;
      multiLocusIndex2[i] =
	multiLocusIndex[i] + pPair->pGenotype[i]->position;
    }
    if (locus < locusList->numLocus - 1) {
      loop_phases (locus + 1, multiLocusIndex2, multiLocusPhase2, newFlipMask,
		   dWeight);
    } else {			/* got a complete multilocus parental pair */
      calculateFlag = 1;
      ppairMatrix[multiLocusPhase2[proband]][multiLocusPhase2[spouse]].
	likelihoodIndex = multiLocusIndex2[proband];
      if (pNucFam->childProbandFlag == TRUE) {	/* a child is the proband */
	phase[DAD] = multiLocusPhase2[DAD];
	phase[MOM] = multiLocusPhase2[MOM];
	/* child is the proband */
	for (i = DAD; i <= MOM; i++) {
	  if (pNucFam->firstHetLocus[i] >= 0 &&
	      pHaplo->phase[i][pNucFam->firstHetLocus[i]] != 0) {
	    /*
	     * first het locus has a
	     * reverse phase than the
	     * origninal potentially we
	     * could benefit from
	     * directly using the
	     * likelihood calculated with
	     * different phase
	     */
	    /*
	     * if this parent is a
	     * founder or a proband, then
	     * flip==original
	     */
	    if (pNucFam->pParents[i]->pParents[DAD] == NULL) {
	      multiLocusPhaseFlip[i] =
		multiLocusPhase2[i] ^ pNucFam->hetLocusBits[i];
	      phase[i] = multiLocusPhaseFlip[i];
	      if (ppairMatrix[phase[proband]][phase[spouse]].count > 0)
		calculateFlag = 0;
	    }
	  }
	}
	if (calculateFlag == 0) {
	  ppairMatrix[phase[proband]][phase[spouse]].count++;
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE) {
#if 0
	    KLOG (LOGLIKELIHOOD, LOGDEBUG,
		  "\t\t likelihood (%d) = %e\n",
		  ppairMatrix[phase[proband]][phase[spouse]].
		  likelihoodIndex, evaluateValue (ppairMatrix[phase[proband]]
						  [phase[spouse]].slot.
						  likelihoodPolynomial));
#endif
	  } else
	    KLOG (LOGLIKELIHOOD, LOGDEBUG,
		  "\t\t likelihood (%d) = %e\n",
		  ppairMatrix[phase[proband]][phase[spouse]].
		  likelihoodIndex,
		  ppairMatrix[phase[proband]][phase[spouse]].slot.likelihood);
#else
	  KLOG (LOGLIKELIHOOD, LOGDEBUG,
		"\t\t likelihood (%d) = %e\n",
		ppairMatrix[phase[proband]][phase[spouse]].
		likelihoodIndex,
		ppairMatrix[phase[proband]][phase[spouse]].slot.likelihood);
#endif
	}
      } else {			/* proband is a parent */
	if (pNucFam->firstHetLocus[spouse] >= 0 &&
	    pHaplo->phase[spouse][pNucFam->firstHetLocus[spouse]] != 0) {
	  if (pNucFam->pParents[spouse]->pParents[DAD] == NULL) {
	    /*
	     * non-proband parent is a
	     * founder
	     */
	    /* find the reverse pattern */
	    multiLocusPhaseFlip[spouse] =
	      multiLocusPhase2[spouse] ^ pNucFam->hetLocusBits[spouse];
	    if (ppairMatrix[multiLocusPhase2[proband]]
		[multiLocusPhaseFlip[spouse]].count > 0) {

	      /*
	       * increase the count
	       * on the original
	       */
	      ppairMatrix[multiLocusPhase2[proband]]
		[multiLocusPhaseFlip[spouse]].count++;
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == FALSE)
		KLOG (LOGLIKELIHOOD, LOGDEBUG,
		      "\t\t likelihood (%d) = %e\n",
		      ppairMatrix[multiLocusPhase2[proband]]
		      [multiLocusPhaseFlip[spouse]].
		      likelihoodIndex, ppairMatrix[multiLocusPhase2[proband]]
		      [multiLocusPhaseFlip[spouse]].slot.likelihood);
#else
	      KLOG (LOGLIKELIHOOD, LOGDEBUG,
		    "\t\t likelihood (%d) = %e\n",
		    ppairMatrix[multiLocusPhase2[proband]]
		    [multiLocusPhaseFlip[spouse]].likelihoodIndex,
		    ppairMatrix[multiLocusPhase2[proband]]
		    [multiLocusPhaseFlip[spouse]].slot.likelihood);
#endif
	      calculateFlag = 0;
	    }
	  }
	} else if (pNucFam->firstHetLocus[proband] >= 0 &&
		   pHaplo->phase[proband][pNucFam->
					  firstHetLocus[proband]] != 0) {
	  /* find the reverse pattern */
	  multiLocusPhaseFlip[proband] =
	    multiLocusPhase2[proband] ^ pNucFam->hetLocusBits[proband];
	  if (ppairMatrix[multiLocusPhaseFlip[proband]]
	      [multiLocusPhase2[spouse]].count > 0) {
	    /*
	     * make sure we have
	     * calculated for this
	     * pattern before
	     */
	    ppairMatrix[multiLocusPhase2[proband]]
	      [multiLocusPhase2[spouse]].count = 1;
	    likelihoodIndex = ppairMatrix[multiLocusPhaseFlip[proband]]
	      [multiLocusPhase2[spouse]].likelihoodIndex;
#ifndef NO_POLYNOMIAL
	    if (modelOptions.polynomial == TRUE) {
	      ppairMatrix[multiLocusPhase2[proband]]
		[multiLocusPhase2[spouse]].slot.
		likelihoodPolynomial =
		ppairMatrix[multiLocusPhaseFlip[proband]]
		[multiLocusPhase2[spouse]].slot.likelihoodPolynomial;
	      pProband->pLikelihood[multiLocusIndex2[proband]].
		wtslot.weightPolynomial =
		pProband->pLikelihood[likelihoodIndex].wtslot.
		weightPolynomial;
	    } else {
	      ppairMatrix[multiLocusPhase2[proband]]
		[multiLocusPhase2[spouse]].slot.likelihood =
		ppairMatrix[multiLocusPhaseFlip[proband]]
		[multiLocusPhase2[spouse]].slot.likelihood;
	      pProband->pLikelihood[multiLocusIndex2[proband]].
		wtslot.weight =
		pProband->pLikelihood[likelihoodIndex].wtslot.weight;
	      KLOG (LOGLIKELIHOOD, LOGDEBUG,
		    "\t\t likelihood (%d) = %e\n",
		    ppairMatrix[multiLocusPhase2[proband]]
		    [multiLocusPhase2[spouse]].likelihoodIndex,
		    ppairMatrix[multiLocusPhase2[proband]]
		    [multiLocusPhase2[spouse]].slot.likelihood);
	    }
#else
	    ppairMatrix[multiLocusPhase2[proband]]
	      [multiLocusPhase2[spouse]].slot.likelihood =
	      ppairMatrix[multiLocusPhaseFlip[proband]]
	      [multiLocusPhase2[spouse]].slot.likelihood;
	    pProband->pLikelihood[multiLocusIndex2[proband]].wtslot.
	      weight = pProband->pLikelihood[likelihoodIndex].wtslot.weight;
	    KLOG (LOGLIKELIHOOD, LOGDEBUG,
		  "\t\t likelihood (%d) = %e\n",
		  ppairMatrix[multiLocusPhase2[proband]]
		  [multiLocusPhase2[spouse]].likelihoodIndex,
		  ppairMatrix[multiLocusPhase2[proband]]
		  [multiLocusPhase2[spouse]].likelihood);
#endif
	    calculateFlag = 0;
	  }
	}			/* end of finding patterns */
      }				/* end of proband is a parent */

      if (calculateFlag == 1) {
	/* 0 - don't keep result, 1 - keep result, 2 - use result */
	if (multiLocusPhase2[DAD] == 0 && multiLocusPhase2[MOM] == 0) {
	  memset (likelihoodChildCount, 0, sizeof (int) * maxChildren);
	  calcFlag = 1;
	} else {
	  calcFlag = 2;
	  recalculate_child_likelihood (newFlipMask, childProductPtr);
	}
	calculate_likelihood (multiLocusIndex2, multiLocusPhase2,
			      dWeight, childProductPtr);
      }				/* fresh likelihood calculation */
    }				/* end of processing a complete multilocus
				 * parental pair */

  }				/* move on to next pair on this locus */

}				/* end of loop_phases() */


/* calculating childProduct base on previous pattern */
void
recalculate_child_likelihood (int flipMask[2], void *childProduct)
{
  int i, j;
  double childSum = 0;
  Polynomial *childSumPoly = NULL;
  ChildElement *pElement;
  int xmissionIndex[2];

  multCount = 0;
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    *(Polynomial **) childProduct = constant1Poly;
  else
    *(double *) childProduct = 1;
#endif

  for (i = 0; i < pNucFam->numChildren; i++) {
#ifndef NO_POLYNOMIAL
    if (modelOptions.polynomial == TRUE)
      childSumPoly = constant0Poly;
    else
      childSum = 0;
#endif
    for (j = 0; j < likelihoodChildCount[i]; j++) {
      pElement = &likelihoodChildElements[multCount + j];
      for (parent = DAD; parent <= MOM; parent++) {
	xmissionIndex[parent] =
	  pElement->xmissionIndex[parent] ^ flipMask[parent];
      }

#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE) {
	childSumPoly =
	  plusExp (2,
		   1.0, childSumPoly,
		   1.0, timesExp (3,
				  xmissionMatrix[xmissionIndex[DAD]].slot.
				  probPoly[1], 1,
				  xmissionMatrix[xmissionIndex[MOM]].slot.
				  probPoly[2], 1,
				  pElement->fslot.factorPolynomial, 1, 0), 1);
      } else {
	childSum += xmissionMatrix[xmissionIndex[DAD]].slot.prob[1] *
	  xmissionMatrix[xmissionIndex[MOM]].slot.prob[2] *
	  pElement->fslot.factor;
      }
#endif
    }
#ifndef NO_POLYNOMIAL
    if (modelOptions.polynomial == TRUE)
      *(Polynomial **) childProduct =
	timesExp (2, *(Polynomial **) childProduct, 1, childSumPoly, 1, 1);
    else
      *(double *) childProduct *= childSum;
#endif
    multCount += likelihoodChildCount[i];

  }
}

int
calculate_likelihood (int multiLocusIndex[2], int multiLocusPhase[2],
		      void *dWeight[2], void *childProductPtr)
{
  int i;
  int proband;
  int spouse;
  Person *pParent[2];
  double newWeight[2];
  double childProduct = 1;
  double penetrance[2];
  int traitLocus;
  int genoIndex;
  double sum;

#ifndef NO_POLYNOMIAL
  Polynomial *newWeightPolynomial[2];
  Polynomial *penetrancePolynomial[2];
  Polynomial *childProductPolynomial = NULL;
  Polynomial *sumPolynomial = NULL;
#endif
  ConditionalLikelihood *pConditional;

  pParent[DAD] = pNucFam->pParents[DAD];
  pParent[MOM] = pNucFam->pParents[MOM];
  proband = pNucFam->head;
  spouse = pNucFam->spouse;
  int xmissionIndex[2] = { 0, 0 };

#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE) {
    childSum = &sumPolynomial;
    if (calcFlag == 2)
      childProductPolynomial = *(Polynomial **) childProductPtr;
  } else {
    childSum = &sum;
    if (calcFlag == 2)
      childProduct = *(double *) childProductPtr;
  }
#else
  childSum = &sum;
#endif

  for (i = DAD; i <= MOM; i++) {
#ifndef NO_POLYNOMIAL
    if (modelOptions.polynomial == TRUE)
      newWeightPolynomial[i] = constant1Poly;
    else
      newWeight[i] = 1.0;
#else
    newWeight[i] = 1.0;
#endif
    pConditional = &pParent[i]->pLikelihood[multiLocusIndex[i]];
    if (pParent[i]->touchedFlag == TRUE) {
      /* we have worked on this parent before */
      if (pParent[i] != pProband) {
#ifndef NO_POLYNOMIAL
	if (modelOptions.polynomial == TRUE)
	  newWeightPolynomial[i] =
	    timesExp (2,
		      pConditional->lkslot.likelihoodPolynomial, 1,
		      pConditional->wtslot.weightPolynomial, 1, 0);
	//Dec 24
	else
	  newWeight[i] =
	    pConditional->lkslot.likelihood * pConditional->wtslot.weight;
#else
	newWeight[i] =
	  pConditional->lkslot.likelihood * pConditional->wtslot.weight;
#endif

      }
    } else {			/* first time we work on this parent */
      if (pNucFam->pParents[i]->pParents[DAD] == NULL) {
	/*
	 * this parent is a founder. the weight is
	 * just the multiplication of genotype
	 * weights at each locus under LE. The weight
	 * should have been passed in as an input
	 */
	if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) {
#ifndef NO_POLYNOMIAL

	  if (modelOptions.polynomial == TRUE)
	    newWeightPolynomial[i] = (Polynomial *) dWeight[i];
	  else
	    newWeight[i] = *((double *) dWeight + i);
#else
	  newWeight[i] = *((double *) dWeight + i);
#endif
	} else if (pParent[i]->loopBreaker == 0) {	/* founder under LD */
#ifndef NO_POlYNOMIAL
	  if (modelOptions.polynomial == TRUE) {
	    get_haplotype_freq (numLocus - 1, i, &newWeightPolynomial[i]);
	  } else {
	    get_haplotype_freq (numLocus - 1, i, &newWeight[i]);
	  }
#else
	  get_haplotype_freq (numLocus - 1, i, &newWeight[i]);
#endif
	}			/* end of founder and LD */
      }				/* founder */
    }				/* end of first time on this parent */
    if (pParent[i]->touchedFlag != TRUE && pParent[i] == pProband) {
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE) {
	pConditional->wtslot.weightPolynomial = newWeightPolynomial[i];
	newWeightPolynomial[i] = constant1Poly;
      } else {
	pConditional->wtslot.weight = newWeight[i];
	newWeight[i] = 1.0;
      }
#else
      pConditional->weight = newWeight[i];
      newWeight[i] = 1.0;
#endif
    }
    /*
     * need to multiply the penetrance if disease locus is in and
     * we haven't calculated any likelihood on this parent before
     */
    traitLocus = locusList->traitLocusIndex;
    if (traitLocus >= 0 && pParent[i]->touchedFlag != TRUE &&
	(pParent[i]->loopBreaker == 0 || pParent[i]->pParents[DAD] != NULL)) {
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
    } else {
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	penetrancePolynomial[i] = constant1Poly;
      else
	penetrance[i] = 1.0;
#else
      penetrance[i] = 1.0;
#endif

    }
  }				/* loop over each parent */

  /* when calcFlag ==2, we just use existing results */
  if (calcFlag != 2) {
    /* now work on the children conditional on this parental pair */
    childProduct = 1;
    multCount = 0;
#ifndef NO_POLYNOMIAL
    if (modelOptions.polynomial == TRUE)
      childProductPolynomial = constant1Poly;
#endif
    for (child = 0; child < pNucFam->numChildren; child++) {
      pChild = pNucFam->ppChildrenList[child];

      xmissionIndex[DAD] = 0;
      xmissionIndex[MOM] = 0;

#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE) {
	sumPolynomial = constant0Poly;
	loop_child_multi_locus_genotype (0, 0, xmissionIndex);
	childProductPolynomial =
	  timesExp (2, childProductPolynomial, 1,
		    (Polynomial *) sumPolynomial, 1, 0);
      } else {
	sum = 0;
	loop_child_multi_locus_genotype (0, 0, xmissionIndex);
	childProduct *= sum;
      }
#else
      sum = 0;
      loop_child_multi_locus_genotype (0, 0, xmissionIndex);
      childProduct *= sum;
#endif
    }				/* looping over all children */
  }

  /* results processing */
  ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
    likelihoodIndex = multiLocusIndex[proband];
  ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].count = 1;
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE) {
    ppairMatrix[multiLocusPhase[proband]]
      [multiLocusPhase[spouse]].slot.likelihoodPolynomial =
      timesExp (5,
		childProductPolynomial, 1,
		newWeightPolynomial[proband], 1,
		newWeightPolynomial[spouse], 1,
		penetrancePolynomial[proband], 1,
		penetrancePolynomial[spouse], 1, 1);
#if 0
    KLOG (LOGLIKELIHOOD, LOGDEBUG, "\t\t likelihood (%d) = %e\n",
	  ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
	  likelihoodIndex,
	  evaluateValue (ppairMatrix[multiLocusPhase[proband]]
			 [multiLocusPhase[spouse]].slot.
			 likelihoodPolynomial));
#endif
  } else {
    /* save it */
    ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].
      slot.likelihood = newWeight[proband] * newWeight[spouse] *
      penetrance[proband] * penetrance[spouse] * childProduct;
    KLOG (LOGLIKELIHOOD, LOGDEBUG,
	  "Parents: DAD(%s) weight %e   MOM(%s) weight %e \n",
	  pParent[DAD]->sID, newWeight[proband], pParent[MOM]->sID,
	  newWeight[spouse]);
    KLOG (LOGLIKELIHOOD, LOGDEBUG,
	  "Parents: DAD(%s) pen %e   MOM(%s) pen %e \n", pParent[DAD]->sID,
	  penetrance[proband], pParent[MOM]->sID, penetrance[spouse]);
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
get_haplotype_freq (int locus, int parent, void *freqPtr)
{
  int origLocus1, origLocus2;	/* locus indices in the
				 * original locus list for
				 * the two loci in LD */
  double freq[2] = { 0, 0 };	/* variable to store the calculated
				 * frequency */
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
  if (modelOptions.polynomial == TRUE) {
    freqPolynomial[0] = constant0Poly;
    freqPolynomial[1] = constant0Poly;
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
  /*
   * now find the corresponding haplotype frequency : 2 haplotypes
   * paternal haplotype & maternal haplotype
   */
  pPair1 =
    &pHaplo->ppParentalPair[locus - 1][pHaplo->pParentalPairInd[locus - 1]];
  pPair2 = &pHaplo->ppParentalPair[locus][pHaplo->pParentalPairInd[locus]];
  for (i = DAD; i <= MOM; i++) {
    /* allele ID in the first locus */
    alleleID1 = pPair1->pGenotype[parent]->allele[i];
    /* allele ID in the second locus */
    alleleID2 = pPair2->pGenotype[parent]->allele[i];
    pAlleleSet1 = pLocus1->ppAlleleSetList[alleleID1 - 1];
    pAlleleSet2 = pLocus2->ppAlleleSetList[alleleID2 - 1];
#ifndef NO_POLYNOMIAL
    if (modelOptions.polynomial == TRUE) {
      freqPolynomial[i] = constant0Poly;
    } else
      freq[i] = 0;
#else
    freq[i] = 0;
#endif

    for (k = 0; k < pAlleleSet1->numAllele; k++) {
      for (l = 0; l < pAlleleSet2->numAllele; l++) {
	allele1 = pAlleleSet1->pAlleles[k];
	allele2 = pAlleleSet2->pAlleles[l];
#ifndef NO_POLYNOMIAL
	if (modelOptions.polynomial == TRUE) {
	  sprintf (vName, "ppHaploFreq[%d][%d]", allele1 - 1, allele2 - 1);
	  freqPolynomial[i] =
	    plusExp (2, 1.0, freqPolynomial[i], 1.0,
		     variableExp (&pLDLoci->
				  ppHaploFreq[allele1 - 1][allele2 -
							   1], NULL,
				  'D', vName), 1);
	} else
	  freq[i] += pLDLoci->ppHaploFreq[allele1 - 1][allele2 - 1];
#else
	freq[i] += pLDLoci->ppHaploFreq[allele1 - 1][allele2 - 1];
#endif
      }
    }

#ifndef NO_POLYNOMIAL
    if (modelOptions.polynomial == TRUE) {

      *(Polynomial **) freqPtr =
	timesExp (2, freqPolynomial[0], 1, freqPolynomial[1], 1, 0);
    } else {
      *(double *) freqPtr = freq[0] * freq[1];
    }
#else
    *(double *) freqValue = freq[0] * freq[1];
#endif
  }				/* end of loop of parents */
}

/*
 * loop over a child's list of genotypes that are compatible with the
 * parental pair retrieve the transmission probability saved in the
 * transmission matrix sum the likelihood for each genotype configuration
 */
int
loop_child_multi_locus_genotype (int locus, int multiLocusIndex,
				 int xmissionIndex[2])
{
  int i;
  int newMultiLocusIndex;
  int newXmissionIndex[2];
  ParentalPair *pParentalPair;

  /* number of possible genotypes at this locus for this child */
  /* child's conditional likelihood offset for the multilocus genotype this function is building */
  multiLocusIndex *= pChild->pSavedNumGenotype[locusList->pLocusIndex[locus]];
  /* build the index to xmission matrix for paternal inheritance and maternal inheritance */
  xmissionIndex[DAD] <<= 2;
  xmissionIndex[MOM] <<= 2;

  /* loop through all of this child's compatible genotypes at this locus */
  pParentalPair =
    &pHaplo->ppParentalPair[locus][pHaplo->pParentalPairInd[locus]];
  for (i = 0; i < pParentalPair->pChildGenoLen[child]; i++) {
    pGenotype = pParentalPair->pppChildGenoList[child][i];
    /* record the index to the genotype list for this child */
    pHaplo->pChildGenoInd[locus] = i;
    KLOG (LOGLIKELIHOOD, LOGDEBUG,
	  "\t child %s locus %4d -> %4d|%-4d \n",
	  pChild->sID, locusList->pLocusIndex[locus],
	  pGenotype->allele[DAD], pGenotype->allele[MOM]);
    /* record this child's conditional likelihood index */
    newMultiLocusIndex = multiLocusIndex + pGenotype->position;

    /* check the transmission probability */
    for (parent = DAD; parent <= MOM; parent++) {
      newChromosome[parent] =
	pParentalPair->ppChildInheritance[parent][child][i];
      /* xmissionIndex has already been multiplied by 4 before the loop */
      newXmissionIndex[parent] =
	xmissionIndex[parent] | newChromosome[parent];
    }				/* looping paternal and maternal chromosomes */
    if (locus < locusList->numLocus - 1) {
      loop_child_multi_locus_genotype (locus + 1, newMultiLocusIndex,
				       newXmissionIndex);
    } else {

      /* get the transmission probability from the matrix */
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE) {
	newProbPolynomial =
	  timesExp (2,
		    xmissionMatrix[newXmissionIndex[DAD]].slot.probPoly[1], 1,
		    xmissionMatrix[newXmissionIndex[MOM]].slot.probPoly[2], 1,
		    0);
#if 0
	KLOG (LOGLIKELIHOOD, LOGDEBUG,
	      "\t xmission prob: %f = %f * %f\n",
	      evaluateValue (newProbPolynomial),
	      evaluateValue (xmissionMatrix[newXmissionIndex[DAD]].
			     slot.probPoly[1]),
	      evaluateValue (xmissionMatrix[newXmissionIndex[MOM]].
			     slot.probPoly[2]));
#endif
      } else {
	newProb =
	  xmissionMatrix[newXmissionIndex[DAD]].slot.prob[1] *
	  xmissionMatrix[newXmissionIndex[MOM]].slot.prob[2];
	KLOG (LOGLIKELIHOOD, LOGDEBUG,
	      "\t xmission prob: %f = %f * %f\n", newProb,
	      xmissionMatrix[newXmissionIndex[DAD]].slot.prob[1],
	      xmissionMatrix[newXmissionIndex[MOM]].slot.prob[2]);
	//fprintf(stderr, "newProb=%f newPenetrance=%f\n", newProb, newPenetrance);
      }
#else
      newProb =
	xmissionMatrix[newXmissionIndex[DAD]].slot.prob[1] *
	xmissionMatrix[newXmissionIndex[MOM]].slot.prob[2];
      KLOG (LOGLIKELIHOOD, LOGDEBUG,
	    "\t xmission prob: %f = %f * %f\n", newProb,
	    xmissionMatrix[newXmissionIndex[DAD]].slot.prob[1],
	    xmissionMatrix[newXmissionIndex[MOM]].slot.prob[2]);
#endif

      /* we have completed one multilocus genotype for this child */
      /*
       * check whether we already have some information about this
       * kid we should have if this kid is a connector to another
       * nuclear family we have processed before
       */
      if (calcFlag == 1 && multCount >= maxChildElements) {
	/* resizing likelihoodchildElements array */
	maxChildElements += 1024;
	likelihoodChildElements =
	  (ChildElement *) realloc (likelihoodChildElements,
				    sizeof (ChildElement) * maxChildElements);

      }
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE) {
	if (pChild != pProband) {
	  /* the child is not a proband */
	  if (pChild->touchedFlag == 1) {
	    /*
	     * some likelihood calculation has
	     * been calculated for this child
	     */
	    *(Polynomial **) childSum = plusExp (2, 1.0, *(Polynomial **) childSum, 1.0, timesExp (2, newProbPolynomial, 1, pChild->pLikelihood[newMultiLocusIndex].lkslot.likelihoodPolynomial, 1, 0),	//end of timesExp
						 1);
	    if (calcFlag == 1) {
	      likelihoodChildElements[multCount].fslot.factorPolynomial =
		pChild->pLikelihood[newMultiLocusIndex].lkslot.
		likelihoodPolynomial;
	    }
	    //end of plusExp
#if 0
	    KLOG (LOGLIKELIHOOD, LOGDEBUG,
		  "\t use already calculated child prob %e \n",
		  evaluateValue (pChild->
				 pLikelihood[newMultiLocusIndex].
				 lkslot.likelihoodPolynomial));
#endif
	  } else if (locusList->traitLocusIndex >= 0)
	    /*
	     * first time working on this child's
	     * current multilocus genotype and we
	     * need to consider penetrance
	     */
	  {
	    traitGenoIndex =
	      pHaplo->pChildGenoInd[locusList->traitLocusIndex];
	    pTraitParentalPair =
	      &pHaplo->ppParentalPair[locusList->
				      traitLocusIndex][pHaplo->
						       pParentalPairInd
						       [locusList->
							traitLocusIndex]];
	    *(Polynomial **) childSum =
	      plusExp (2,
		       1.0, *(Polynomial **) childSum,
		       1.0, timesExp (2,
				      newProbPolynomial, 1,
				      pTraitParentalPair->
				      pppChildGenoList[child]
				      [traitGenoIndex]->penslot.penetrancePolynomial, 1, 0),	//end of timesExp
		       1);

	    if (calcFlag == 1) {
	      likelihoodChildElements[multCount].fslot.factorPolynomial =
		pTraitParentalPair->pppChildGenoList[child]
		[traitGenoIndex]->penslot.penetrancePolynomial;
	    }
	  } else {
	    /*
	     * no trait locus and new to this
	     * child
	     */
	    *(Polynomial **) childSum = plusExp (2,
						 1.0,
						 *(Polynomial
						   **) childSum,
						 1.0, newProbPolynomial, 1);
	    if (calcFlag == 1) {
	      likelihoodChildElements[multCount].fslot.factorPolynomial =
		constant1Poly;
	    }
	  }
	} else {		/* this child is proband */
	  *(Polynomial **) childSum = plusExp (2, 1.0, *(Polynomial **)
					       childSum, 1.0,
					       newProbPolynomial, 1);
	  if (calcFlag == 1) {
	    likelihoodChildElements[multCount].fslot.factorPolynomial =
	      constant1Poly;
	  }
	}
#if 0
	KLOG (LOGLIKELIHOOD, LOGDEBUG,
	      "\t child sum %e \n",
	      evaluateValue (*(Polynomial **) childSum));
#endif
      } else {			/* PE is not turned on */
	if (pChild != pProband) {
	  /* the child is not a proband */
	  if (pChild->touchedFlag == 1) {
	    /*
	     * some likelihood calculation has
	     * been done for this child
	     */
	    *(double *) childSum += newProb *
	      pChild->pLikelihood[newMultiLocusIndex].lkslot.likelihood;
	    KLOG (LOGLIKELIHOOD, LOGDEBUG,
		  "\t use already calculated child prob %e \n",
		  pChild->pLikelihood[newMultiLocusIndex].lkslot.likelihood);
	    if (calcFlag == 1) {
	      likelihoodChildElements[multCount].fslot.factor =
		pChild->pLikelihood[newMultiLocusIndex].lkslot.likelihood;
	    }
	  } else if (locusList->traitLocusIndex >= 0)
	    /*
	     * first time working on this child's
	     * current multilocus genotype and we
	     * need to consider penetrance
	     */
	  {
	    traitGenoIndex =
	      pHaplo->pChildGenoInd[locusList->traitLocusIndex];
	    pTraitParentalPair =
	      &pHaplo->ppParentalPair[locusList->
				      traitLocusIndex][pHaplo->
						       pParentalPairInd
						       [locusList->
							traitLocusIndex]];
	    *(double *) childSum +=
	      newProb *
	      pTraitParentalPair->
	      pppChildGenoList[child][traitGenoIndex]->penslot.penetrance;
	    if (calcFlag == 1) {
	      likelihoodChildElements[multCount].fslot.factor =
		pTraitParentalPair->pppChildGenoList[child][traitGenoIndex]->
		penslot.penetrance;
	    }
	    KLOG (LOGLIKELIHOOD, LOGDEBUG,
		  "child penetrance %e\n",
		  pTraitParentalPair->
		  pppChildGenoList[child][traitGenoIndex]->
		  penslot.penetrance);

	  } else {
	    *(double *) childSum += newProb;
	    if (calcFlag == 1) {
	      likelihoodChildElements[multCount].fslot.factor = 1;
	    }
	  }
	} else {		/* this child is proband */
	  /*
	   * penetrance if applicable will be figured
	   * into later
	   */
	  *(double *) childSum += newProb;
	  if (calcFlag == 1) {
	    likelihoodChildElements[multCount].fslot.factor = 1;
	  }
	}
	KLOG (LOGLIKELIHOOD, LOGDEBUG,
	      "\t child sum %e \n", *(double *) childSum);
      }
#else /* PE is not compiled in */
      if (pChild != pProband) {
	/* the child is not a proband */
	if (pChild->touchedFlag == 1) {
	  /*
	   * some likelihood calculation has been done
	   * for this child
	   */
	  *(double *) childSum += newProb *
	    pChild->pLikelihood[newMultiLocusIndex].lkslot.likelihood;
	  KLOG (LOGLIKELIHOOD, LOGDEBUG,
		"\t use already calculated child prob %e \n",
		pChild->pLikelihood[newMmultiLocusIndex].lkslot.likelihood);
	} else if (locusList->traitLocusIndex >= 0)
	  /*
	   * first time working on this child's current
	   * multilocus genotype and we need to
	   * consider penetrance
	   */
	{
	  traitGenoIndex = pHaplo->pChildGenoInd[locusList->traitLocusIndex];
	  pTraitParentalPair =
	    &pHaplo->ppParentalPair[locusList->
				    traitLocusIndex][pHaplo->
						     pParentalPairInd
						     [locusList->
						      traitLocusIndex]];
	  *(double *) childSum +=
	    newProb *
	    pTraitParentalPair->
	    pppChildGenoList[child][traitGenoIndex]->penslot.penetrance;
	} else {
	  *(double *) childSum += newProb;
	}
      } else {			/* this child is proband */
	/*
	 * penetrance if applicable will be figured into
	 * later
	 */
	*(double *) childSum += newProb;
      }
      KLOG (LOGLIKELIHOOD, LOGDEBUG,
	    "\t child sum %e \n", *(double *) childSum);
#endif


      if (calcFlag == 1) {
	likelihoodChildElements[multCount].xmissionIndex[DAD] =
	  newXmissionIndex[DAD];
	likelihoodChildElements[multCount].xmissionIndex[MOM] =
	  newXmissionIndex[MOM];
	likelihoodChildCount[child]++;
	multCount++;
      }
    }				/* end of processing one complete multilocus genotype */

  }				/* loop possible genotypes at this locus */
  return 0;
}


/*
 * A recursive call to build transmission probability matrix pMatrix - pass
 * in matrix pointer - This should have been pre-allocated totalLoci - total
 * number of loci loc - current locus prob -
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
   * both (3) which indicates the parent is homozygous at that locus 
   * added 0 for easy handle of pattern flip, 0 is equivalent of 3 */
  for (pattern = 0; pattern <= 3; pattern++) {
    /* sex averaged or sex specific map */
    for (i = 0; i < 3; i++) {

#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE) {
	newProbPoly[i] = (Polynomial *) prob[i];
	newProbPoly2[i] = (Polynomial *) prob2[i];
	newHetProbPoly[i] = (Polynomial *) hetProb[i];
      } else {
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
    newCellIndex = cellIndex * 4 + pattern;
    newLastHetLoc = lastHetLoc;
    if (pattern != 3 && pattern != 0) {
      /* parent is not homozygous */
      if (lastHetLoc != -1) {
	if (prevPattern != 3 && prevPattern != 0) {
	  /* previous locus for the parent is het and 
	     current locus pattern is either paternal or maternal */
	  if (prevPattern == pattern) {
	    /* no recombination */
	    for (i = 0; i < 3; i++) {
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE) {
		if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		  newProbPoly[i] = newProbPoly[0];
		} else {
		  sprintf (vName1, "theta%d_%d", i, loc);
		  newProbPoly[i] =
		    timesExp (2, newProbPoly[i], 1,
			      plusExp (2, 1.0,
				       constantExp
				       (1.0), -1.0,
				       variableExp
				       (&locusList->
					pPrevLocusDistance
					[i][loc], NULL,
					'D', vName1), 0), 1, 0);
		}
	      } else {
		newProb[i] *= (1 - locusList->pPrevLocusDistance[i][loc]);
	      }
#else
	      newProb[i] *= (1 - locusList->pPrevLocusDistance[i][loc]);
#endif
	    }
	  } else {
	    /* recombination */
	    for (i = 0; i < 3; i++)
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE) {
		if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		  newProbPoly[i] = newProbPoly[0];
		} else {
		  sprintf (vName1, "theta%d_%d", i, loc);
		  newProbPoly[i] =
		    timesExp (2, newProbPoly[i], 1,
			      variableExp (&locusList->
					   pPrevLocusDistance
					   [i][loc], NULL,
					   'D', vName1), 1, 0);

		}
	      } else {
		if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		  newProb[i] = newProb[0];
		} else {
		  newProb[i] *= locusList->pPrevLocusDistance[i][loc];
		}
	      }
#else
	      newProb[i] *= locusList->pPrevLocusDistance[i][loc];
#endif
	  }
	} else {
	  /* previous locus at parent is homo and current locus is het */
	  for (i = 0; i < 3; i++) {
	    if (pattern == 1) {
	      /* paternal inheritance for this locus
	         either no recombination from previous paternal strand 
	         or recombination from previous maternal strand */
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE) {
		if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		  newProbPoly[i] = newProbPoly[0];
		} else {
		  sprintf (vName1, "theta%d_%d", i, loc);
		  newProbPoly[i] =
		    plusExp (2, 1.0, timesExp (2, (Polynomial *)
					       prob[i], 1,
					       plusExp (2, 1.0,
							constantExp
							(1.0),
							-1.0,
							variableExp
							(&locusList->
							 pPrevLocusDistance[i]
							 [loc],
							 NULL,
							 'D',
							 vName1),
							0), 1,
					       0), 1.0,
			     timesExp (2, (Polynomial *)
				       prob2[i], 1,
				       variableExp
				       (&locusList->
					pPrevLocusDistance
					[i][loc], NULL,
					'D', vName1), 1, 0), 0);
		}
	      } else {
		if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		  newProb[i] = newProb[0];
		} else {
		  newProb[i] = *((double *) prob[i]) *
		    (1 -
		     locusList->
		     pPrevLocusDistance[i][loc]) +
		    *((double *) prob2[i]) *
		    locusList->pPrevLocusDistance[i][loc];
		}
	      }
#else
	      if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		newProb[i] = newProb[0];
	      } else {
		newProb[i] = *((double *) prob[i]) *
		  (1 -
		   locusList->
		   pPrevLocusDistance[i][loc]) +
		  *((double *) prob2[i]) *
		  locusList->pPrevLocusDistance[i][loc];
	      }
#endif
	    } else {
	      /* has to be maternal */
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE) {
		if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		  newProbPoly[i] = newProbPoly[0];
		} else {
		  sprintf (vName1, "theta%d_%d", i, loc);
		  newProbPoly[i] =
		    plusExp (2, 1.0, timesExp (2, (Polynomial *)
					       prob2[i], 1,
					       plusExp (2, 1.0,
							constantExp
							(1.0),
							-1.0,
							variableExp
							(&locusList->
							 pPrevLocusDistance[i]
							 [loc],
							 NULL,
							 'D',
							 vName1),
							0), 1,
					       0), 1.0,
			     timesExp (2, (Polynomial *)
				       prob[i], 1,
				       variableExp
				       (&locusList->
					pPrevLocusDistance
					[i][loc], NULL,
					'D', vName1), 1, 0), 0);
		}

	      } else if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		newProb[i] = newProb[0];
	      } else {
		newProb[i] = *((double *) prob2[i]) *
		  (1 -
		   locusList->
		   pPrevLocusDistance[i][loc]) +
		  *((double *) prob[i]) *
		  locusList->pPrevLocusDistance[i][loc];
	      }

#else
	      if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		newProb[i] = newProb[0];
	      } else {
		newProb[i] = *((double *) prob2[i]) *
		  (1 -
		   locusList->
		   pPrevLocusDistance[i][loc]) +
		  *((double *) prob[i]) *
		  locusList->pPrevLocusDistance[i][loc];
	      }
#endif
	    }
	  }

	}			/* end of prevPattern is homo and current pattern is het */
      } /* end of prevHetLoc != -1 */
      else {
	/* we don't have any het locus yet, this locus is the first het */
#ifndef NO_POLYNOMIAL
	if (modelOptions.polynomial == TRUE) {
	  for (i = 0; i < 3; i++)
	    newProbPoly[i] = constantExp (0.5);
	} else {
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
      if (modelOptions.polynomial == TRUE) {
	for (i = 0; i < 3; i++)
	  newHetProbPoly[i] = newProbPoly[i];
      } else {
	for (i = 0; i < 3; i++)
	  newHetProbPtr[i] = newProbPtr[i];
      }
#else
      for (i = 0; i < 3; i++)
	newHetProbPtr[i] = newProbPtr[i];
#endif


    } /* end of current pattern is not homo */
    else {			/* current pattern is homo */
      if (lastHetLoc == -1)
	/* nothing needs to be done for this locus */
	;
      else {
	if (loc == totalLoci - 1) {
	  /* this is the last locus and it's homo, take the previous het locus */
	  for (i = 0; i < 3; i++) {
#ifndef NO_POLYNOMIAL
	    if (modelOptions.polynomial == TRUE) {
	      newProbPoly[i] = (Polynomial *) hetProb[i];
	    } else
	      newProb[i] = *(double *) hetProb[i];
#else
	    newProb[i] = *(double *) hetProb[i];
#endif
	  }
	} else {
	  if (prevPattern == 3 || prevPattern == 0) {	/* previous locus pattern is homo */
	    for (i = 0; i < 3; i++) {

#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE) {
		if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		  newProbPoly[i] = newProbPoly[0];
                  newProbPoly2[i] = newProbPoly2[0];
		} else {
		  sprintf (vName1, "theta%d_%d", i, loc);
		  newProbPoly[i] =
		    plusExp (2, 1.0, timesExp (2, (Polynomial *)
					       prob[i], 1,
					       plusExp (2, 1.0,
							constantExp
							(1.0),
							-1.0,
							variableExp
							(&locusList->
							 pPrevLocusDistance[i]
							 [loc],
							 NULL,
							 'D',
							 vName1),
							0), 1,
					       0), 1.0,
			     timesExp (2, (Polynomial *)
				       prob2[i], 1,
				       variableExp
				       (&locusList->
					pPrevLocusDistance
					[i][loc], NULL,
					'D', vName1), 1, 0), 0);
		  newProbPoly2[i] =
		    plusExp (2, 1.0, timesExp (2, (Polynomial *)
					       prob2[i], 1,
					       plusExp (2, 1.0,
							constantExp
							(1.0),
							-1.0,
							variableExp
							(&locusList->
							 pPrevLocusDistance[i]
							 [loc],
							 NULL,
							 'D',
							 vName1),
							0), 1,
					       0), 1.0,
			     timesExp (2, (Polynomial *)
				       prob[i], 1,
				       variableExp
				       (&locusList->
					pPrevLocusDistance
					[i][loc], NULL,
					'D', vName1), 1, 0), 0);
		}

	      } else {
		if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		  newProb[i] = newProb[0];
		  newProb2[i] = newProb2[0];
		} else {
		  newProb[i] = *(double *) prob[i] *
		    (1 -
		     locusList->
		     pPrevLocusDistance[i][loc]) +
		    *((double *) prob2[i]) *
		    locusList->pPrevLocusDistance[i][loc];

		  newProb2[i] = *(double *) prob2[i] *
		    (1 -
		     locusList->
		     pPrevLocusDistance[i][loc]) +
		    *((double *) prob[i]) *
		    locusList->pPrevLocusDistance[i][loc];
		}
	      }
#else
	      if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		newProb[i] = newProb[0];
		newProb2[i] = newProb2[0];
	      } else {
		newProb[i] = *(double *) prob[i] *
		  (1 -
		   locusList->
		   pPrevLocusDistance[i][loc]) +
		  *((double *) prob2[i]) *
		  locusList->pPrevLocusDistance[i][loc];

		newProb2[i] = *(double *) prob2[i] *
		  (1 -
		   locusList->
		   pPrevLocusDistance[i][loc]) +
		  *((double *) prob[i]) *
		  locusList->pPrevLocusDistance[i][loc];
	      }
#endif
	    }
	  } else {		/* prev pattern is het */
	    for (i = 0; i < 3; i++) {
	      if (prevPattern == 1) {
#ifndef NO_POLYNOMIAL
		if (modelOptions.polynomial == TRUE) {
		  if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		    newProbPoly[i] = newProbPoly[0];
		    newProbPoly2[i] = newProbPoly2[0];
		  } else {
		    sprintf (vName1, "theta%d_%d", i, loc);
		    newProbPoly[i] = timesExp (2, (Polynomial *)
					       prob[i], 1,
					       plusExp (2, 1.0,
							constantExp
							(1.0), -1.0,
							variableExp
							(&locusList->
							 pPrevLocusDistance
							 [i][loc],
							 NULL, 'D',
							 vName1), 0), 1, 0);
		    newProbPoly2[i] = timesExp (2, (Polynomial *)
						prob[i], 1,
						variableExp
						(&locusList->
						 pPrevLocusDistance[i]
						 [loc], NULL, 'D',
						 vName1), 1, 0);
		  }
		} else {
		  if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		    newProb[i] = newProb[0];
		    newProb2[i] = newProb2[0];
		  } else {
		    newProb[i] = *(double *) prob[i] *
		      (1 - locusList->pPrevLocusDistance[i][loc]);
		    newProb2[i] =
		      *(double *) prob[i] *
		      locusList->pPrevLocusDistance[i][loc];
		  }
		}
#else
		if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		  newProb[i] = newProb[0];
		  newProb2[i] = newProb2[0];
		} else {
		  newProb[i] = *(double *) prob[i] *
		    (1 - locusList->pPrevLocusDistance[i][loc]);

		  newProb2[i] = *(double *) prob[i] *
		    locusList->pPrevLocusDistance[i][loc];
		}
#endif
	      } else {
#ifndef NO_POLYNOMIAL
		if (modelOptions.polynomial == TRUE) {
		  if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		    newProbPoly[i] = newProbPoly[0];
		    newProbPoly2[i] = newProbPoly2[0];
		  } else {
		    sprintf (vName1, "theta%d_%d", i, loc);
		    newProbPoly2[i] = timesExp (2, (Polynomial *)
						prob[i], 1,
						plusExp (2, 1.0,
							 constantExp
							 (1.0), -1.0,
							 variableExp
							 (&locusList->
							  pPrevLocusDistance
							  [i][loc],
							  NULL, 'D',
							  vName1), 0), 1, 0);
		    newProbPoly[i] = timesExp (2, (Polynomial *)
					       prob[i], 1,
					       variableExp
					       (&locusList->
						pPrevLocusDistance[i]
						[loc], NULL, 'D',
						vName1), 1, 0);
		  }
		} else {
		  if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		    newProb[i] = newProb[0];
		    newProb2[i] = newProb2[0];
		  } else {
		    newProb2[i] =
		      *(double *) prob[i] * (1 -
					     locusList->
					     pPrevLocusDistance[i][loc]);
		    newProb[i] =
		      *(double *) prob[i] *
		      locusList->pPrevLocusDistance[i][loc];
		  }
		}
#else
		if (i > 0 && modelOptions.mapFlag == SEX_AVERAGED) {
		  newProb[i] = newProb[0];
		  newProb2[i] = newProb2[0];
		} else {
		  newProb2[i] = *(double *) prob[i] *
		    (1 - locusList->pPrevLocusDistance[i][loc]);

		  newProb[i] = *(double *) prob[i] *
		    locusList->pPrevLocusDistance[i][loc];
		}
#endif
	      }

	    }
	  }
	}
      }
    }

    if (loc == totalLoci - 1) {
      /* we have a complete set of multilocus inheritance pattern */

      for (i = 0; i < 3; i++) {
#ifndef NO_POLYNOMIAL
	if (modelOptions.polynomial == TRUE) {
	  pMatrix[newCellIndex].slot.probPoly[i] = newProbPoly[i];
	  holdPoly(newProbPoly[i]);
	} else
	  pMatrix[newCellIndex].slot.prob[i] = newProb[i];
#else
	pMatrix[newCellIndex].slot.prob[i] = newProb[i];
#endif
      }

    } else {
      /* move on to next locus */
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE) {
	populate_xmission_matrix (pMatrix, totalLoci,
				  (void *) newProbPoly,
				  (void *) newProbPoly2,
				  (void *) newHetProbPoly,
				  newCellIndex, newLastHetLoc,
				  pattern, loc + 1);
      } else
	populate_xmission_matrix (pMatrix, totalLoci,
				  (void *) newProbPtr,
				  (void *) newProbPtr2,
				  (void *) newHetProbPtr,
				  newCellIndex, newLastHetLoc,
				  pattern, loc + 1);
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

  //int           i;

  *ppMatrix = NULL;
  size = pow (4, totalLoci);
  /* minimal two loci */
  //if (size < 9)
  //size = 9;
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

  for (pattern = 0; pattern <= 2; pattern++) {
    newCellIndex = cellIndex * 4 + pattern;
    if (pattern == 1) {
      pID[loc] = 'P';
    } else if (pattern == 2) {
      pID[loc] = 'M';
    } else {
      pID[loc] = 'B';
    }
    if (loc != totalLoci - 1) {
      /* not complete multi-locus haplotype yet */
      print_xmission_matrix (pMatrix, totalLoci, loc + 1, newCellIndex, pID);
    } else {
      /* print the xmission probability out */
      for (i = 0; i <= loc; i++) {
	fprintf (stderr, "%c", pID[i]);
      }
      fprintf(stderr, ": ");
      /* print out sex averaged xmission probability */
      if(modelOptions.polynomial == TRUE)
	{
	  expTermPrinting(stderr, pMatrix[newCellIndex].slot.probPoly[0], 16);
	  fprintf(stderr, "\n");
	}
      else
        fprintf (stderr, ":%f\n", pMatrix[newCellIndex].slot.prob[0]);
    }
  }
}

/*
 * allocate storage for keeping track of het locus in nuclear families
 * numLocus - number of loci analyzing at a time
 */
void
allocate_nucfam_het (PedigreeSet * pPedigreeList, int numLocus)
{
  int ped;
  Pedigree *pPedigree;
  int fam;
  NuclearFamily *pNucFam;

  for (ped = 0; ped < pPedigreeList->numPedigree; ped++) {
    pPedigree = pPedigreeList->ppPedigreeSet[ped];
    for (fam = 0; fam < pPedigree->numNuclearFamily; fam++) {
      pNucFam = pPedigree->ppNuclearFamilyList[fam];
      if (pNucFam->hetFlag[DAD] == NULL) {
	pNucFam->hetFlag[DAD] = (int *) calloc (sizeof (int), numLocus);
	pNucFam->hetFlag[MOM] = (int *) calloc (sizeof (int), numLocus);
	pNucFam->tmpNumHet[DAD] = (int *) calloc (sizeof (int), numLocus);
	pNucFam->tmpNumHet[MOM] = (int *) calloc (sizeof (int), numLocus);
	pNucFam->relatedPPairStart = (int *) calloc (sizeof (int), numLocus);
	pNucFam->numRelatedPPair = (int *) calloc (sizeof (int), numLocus);
	pNucFam->totalRelatedPPair = (int *) calloc (sizeof (int), numLocus);
      }
    }
  }

}

inline void
clear_ppairMatrix (PPairElement ** ppMatrix)
{
  int i;

  for (i = 0; i <= bitMask[ppairMatrixNumLocus]; i++) {
    memset (ppMatrix[i], 0, ppairMatrixRowSize);
  }
}

inline void
initialize_proband_tmpLikelihood (Person * pPerson)
{
  int i;

  //  int             size = pPerson->numConditionals;
  ConditionalLikelihood *pConditional;

  for (i = 0; i < pPerson->numTmpLikelihood; i++) {
    pConditional = &pPerson->pLikelihood[pPerson->pTmpLikelihoodIndex[i]];
    pConditional->tmpTouched = FALSE;
#ifndef NO_POLYNOMIAL
    if (modelOptions.polynomial == TRUE)
      pConditional->tmpslot.tmpLikelihoodPolynomial = constant0Poly;
    else
      pConditional->tmpslot.tmpLikelihood = 0;
#else
    pConditional->tmpslot.tmpLikelihood = 0;
#endif
  }
  pPerson->numTmpLikelihood = 0;
}

void
populate_pedigree_loopbreaker_genotype_vector (Pedigree * pPed)
{
  int numLoopBreaker = pPed->numLoopBreaker;
  int i;
  Person *pLoopBreaker;

  for (i = 0; i < numLoopBreaker; i++) {
    pLoopBreaker = pPed->loopBreakerList[i];
    pLoopBreaker->loopBreakerStruct->numGenotype = 0;
    populate_loopbreaker_genotype_vector (pLoopBreaker, 0);
    pLoopBreaker->loopBreakerStruct->genotypeIndex = 0;
  }
}

void
populate_loopbreaker_genotype_vector (Person * pLoopBreaker, int locus)
{
  Genotype *pGenotype;
  int index;
  LoopBreaker *pLoopBreakerStruct;

  pGenotype =
    pLoopBreaker->ppSavedGenotypeList[locusList->pLocusIndex[locus]];
  while (pGenotype != NULL) {
    pTempGenoVector[locus] = pGenotype;
    if (locus < locusList->numLocus - 1) {
      populate_loopbreaker_genotype_vector (pLoopBreaker, locus + 1);
    } else {
      /* one complete multilocus genotype */
      pLoopBreakerStruct = pLoopBreaker->loopBreakerStruct;
      index = pLoopBreakerStruct->numGenotype;
      memcpy (pLoopBreakerStruct->genotype[index], pTempGenoVector,
	      sizeof (Genotype *) * locusList->numLocus);
      pLoopBreakerStruct->numGenotype++;

    }
    pGenotype = pGenotype->pSavedNext;
  }
}

/* this function is not needed - so not finished */
void
sync_loopbreaker_duplicates (Pedigree * pPed)
{
  int i;
  Person *pPerson;

  for (i = 0; i < pPed->numPerson; i++) {
    pPerson = pPed->ppPersonList[i];
    if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] == NULL) {
    }
  }
}

/* initialFlag - TRUE - first vector (index all 0) */
int
set_next_loopbreaker_genotype_vector (Pedigree * pPed, int initialFlag)
{
  int numLoopBreaker = pPed->numLoopBreaker;
  int i;
  Person *pLoopBreaker;
  int found;
  LoopBreaker *loopStruct;
  int index;
  int origLocus;
  int locus;
  int ret;

  /* find the next genotype vector for at least one of the loop breaker */
  KLOG (LOGLIKELIHOOD, LOGDEBUG, "Set next loop breaker genotype\n");
  found = FALSE;
  if (initialFlag != TRUE) {
    for (i = 0; i < numLoopBreaker; i++) {
      pLoopBreaker = pPed->loopBreakerList[i];
      loopStruct = pLoopBreaker->loopBreakerStruct;
      /* increase index */
      loopStruct->genotypeIndex++;
      if (loopStruct->genotypeIndex >= loopStruct->numGenotype) {
	loopStruct->genotypeIndex = 0;
      } else {
	found = TRUE;
	break;
      }
    }

  } else {
    found = TRUE;
#if 0
    for (i = 0; i < numLoopBreaker; i++) {
      pLoopBreaker = pPed->loopBreakerList[i];
      loopStruct = pLoopBreaker->loopBreakerStruct;
      loopStruct->genotypeIndex = loopStruct->numGenotype - 1;
    }
#endif
  }

  if (found == FALSE)
    return -1;

  /* set the genotype list with the selected genotype vector */
  for (i = 0; i < numLoopBreaker; i++) {
    pLoopBreaker = pPed->loopBreakerList[i];
    loopStruct = pLoopBreaker->loopBreakerStruct;
    index = loopStruct->genotypeIndex;
    KLOG (LOGLIKELIHOOD, LOGDEBUG,
	  "Fix pedigree %s loop breaker %s to the genotype below (%d/%d):\n",
	  pPed->sPedigreeID, pLoopBreaker->sID,
	  loopStruct->genotypeIndex + 1, loopStruct->numGenotype);
    for (locus = 0; locus < locusList->numLocus; locus++) {
      origLocus = locusList->pLocusIndex[locus];
      pLoopBreaker->ppGenotypeList[origLocus] =
	loopStruct->genotype[index][locus];
      loopStruct->genotype[index][locus]->pNext = NULL;
      pLoopBreaker->pNumGenotype[origLocus] = 1;
      KLOG (LOGLIKELIHOOD, LOGDEBUG, "\t %d-> %d|%d \n",
	    locus, loopStruct->genotype[index][locus]->allele[DAD],
	    loopStruct->genotype[index][locus]->allele[MOM]);
    }
  }

  /*
   * as the loop breakers are set with fixed genotype, redo genotype
   * elimination (without actually remove any genotype - only the links
   * will get updated
   */
  for (i = 0; i < locusList->numLocus; i++) {
    origLocus = locusList->pLocusIndex[i];
    ret = pedigree_genotype_elimination (origLocus, pPed);
    if (ret < 0)
      return -2;
  }

  return 0;
}
