/**********************************************************************
 * Copyright 2010, Nationwide Children's Research Institute.
 * All rights reserved.
 * Permission is hereby given to use this software
 * for non-profit educational purposes only.
 **********************************************************************/

char *likelihoodVersion = "$Id$";

/*
 * This file contains functions to  compute likelihood for all the pedigrees,
 * peeling procedure etc.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <sys/types.h>

#include "pedlib.h"
#include "locus.h"
#include "../utils/utils.h"     /* for logging */
#include "../utils/sw.h"
#include "likelihood.h"
#include "genotype_elimination.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include <dlfcn.h>

#ifdef STUDYDB
#include "../database/databaseSupport.h"
#endif

//extern FILE *fpCond;
extern struct polynomial **variableList;

char partialPolynomialFunctionName[MAX_PFN_LEN + 1];

/* this is for working out loop breaker's multilocus genotypes */
Genotype **pTempGenoVector;

/* transmission probability matrix */
XMission *xmissionMatrix = NULL;

/* temporary likelihood results for related parental pairs */
typedef struct PPairElement
{
  /* index to the conditional likelihood array of the proband's */
  int likelihoodIndex;
  /* number of the same likelihood */
  short count;
  /* temporary likelihood result */
  union
  {
    double likelihood;
    /* likelihood polynomial under polynomial mode */
    Polynomial *likelihoodPolynomial;
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

Polynomial *constant0Poly;
Polynomial *constant1Poly;

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

Polynomial *newProbPolynomial = NULL;
int newChromosome[2];
int numLocus;
NuclearFamily *pNucFam;

/**
  Genotypes for parental pairs are said to be "related" when they are
  identical except for a phase change. When this is the case, the same
  likelihood calculation can be used for them all, and a trick is used to
  identify the different xmission matrix entries that need to be
  referenced to account for the phase difference.

  The calcFlag is used to keep track of whether we are reusing old
  calculation results from a related parental pair genotype (calcFlag = 2),
  producing  results we need to keep for reuse on related parental pair 
  genotypes (calcFlag = 1), or producing one-time only results that don't 
  need to be kept because there are no related parental pair genotypes
  (calcFlag = 0).

*/
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

typedef struct probandCondL
{
  char *pPedigreeID;
  char *pProbandID;
  char *pAllele1;
  char *pAllele2;
  int allele1;
  int allele2;
  double trait;
  double condL;
  float prob;
} probandCondL;
probandCondL *pCondSet = NULL;
int numCond;

/* function prototypes */
void recalculate_child_likelihood (int[], void *);
int peel_graph (NuclearFamily *, Person *, int);
int compute_nuclear_family_likelihood (int);
int loop_parental_pair (int, int[], void *[]);
int loop_child_multi_locus_genotype (int, int, int[]);
void get_haplotype_freq (int, int, void *);

int loop_child_proband_genotype (int, int, int);

void loop_phases (int, int[], int[], int[], void *[]);
int calculate_likelihood (int[], int[], void *[], void *);

inline void clear_ppairMatrix (PPairElement **);
inline void initialize_proband_tmpLikelihood (Person *);
void populate_pedigree_loopbreaker_genotype_vector (Pedigree *);
void populate_loopbreaker_genotype_vector (Person *, int);
int set_next_loopbreaker_genotype_vector (Pedigree *, int);

#include "../kelvinGlobals.h"

/*
 * before likelihood calculation, pre-allocate space to store conditional
 * likelihoods 
 * numLocus - number of loci we analyze at a time
 */
int allocate_likelihood_space (PedigreeSet * pPedigreeList, int numLocus1)
{
  Pedigree *pPedigree;
  int i;
  int size;

  numLocus = numLocus1;
  /* this is for loop breaker multilocus genotypes */
  MALCHOKE (pTempGenoVector, sizeof (Genotype *) * numLocus, Genotype **);

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
  CALCHOKE (ppairMatrix, (size_t) size, sizeof (PPairElement *), PPairElement **);
  for (i = 0; i < size; i++) {
    CALCHOKE (ppairMatrix[i], (size_t) size, sizeof (PPairElement), PPairElement *);
  }
  if (bitMask != NULL)
    free (bitMask);
  MALCHOKE (bitMask, sizeof (int) * (numLocus + 1), void *);
  for (i = 0; i <= numLocus; i++) {
    bitMask[i] = pow (2, i) - 1;
  }

  /* likelihood work space */
  pHaplo = &parentalPairSpace;

  /* pre allocate child likelihood elements */
  maxChildElements = 1024;
  CALCHOKE (likelihoodChildElements, sizeof (ChildElement), (size_t) maxChildElements, ChildElement *);

  maxChildren = 20;
  CALCHOKE (likelihoodChildCount, sizeof (int), (size_t) maxChildren, int *);

  if (pCondSet == NULL) {
    MALCHOKE (pCondSet, sizeof (probandCondL) * 4, probandCondL *);
    numCond = 4;
  }

  return 0;
}

/* free the storage space for conditionals */
void free_likelihood_space (PedigreeSet * pPedigreeList)
{
  Pedigree *pPedigree;
  int i;

  for (i = 0; i < pPedigreeList->numPedigree; i++) {
    pPedigree = pPedigreeList->ppPedigreeSet[i];
    free_multi_locus_genotype_storage (pPedigree);
  }

  /* free storage for temporary likelihood for similar parental pairs */
  for (i = 0; i < pow (2, analysisLocusList->numLocus); i++) {
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

int build_likelihood_polynomial (PedigreeSet * pPedigreeList)
{
  Pedigree *pPedigree;
  int i;
  int status;

  char polynomialFunctionName[MAX_PFN_LEN + 1];

  if (modelOptions->polynomial != TRUE)
    return EXIT_FAILURE;

  /* First build (or restore) all of the pedigrees */
  for (i = 0; i < pPedigreeList->numPedigree; i++) {
    pPedigree = pPedigreeList->ppPedigreeSet[i];
    if (pPedigree->load_flag == 0) {    // Skip everything, we're using saved results
      if (pPedigree->likelihoodPolynomial == NULL) {    // There's no polynomial, so come up with one
#ifdef POLYCHECK_DL
        pPedigree->cLikelihoodPolynomial = NULL;        // Get rid of the old tree comparison poly
#endif
        // Construct polynomialFunctionName for attempted load and maybe compilation
        sprintf (polynomialFunctionName, partialPolynomialFunctionName, pPedigree->sPedigreeID);
#ifdef POLYUSE_DL
        if ((pPedigree->likelihoodPolynomial = restoreExternalPoly (polynomialFunctionName)) == NULL) {
#endif
          // Failed to load, construct it
          initialize_multi_locus_genotype (pPedigree);
          status = compute_pedigree_likelihood (pPedigree);
#ifdef POLYCODE_DL
          // Used to skip compilation of simple polys here, but there are simple ones that are tough builds.
          pPedigree->likelihoodPolyList = buildPolyList ();
          polyListSorting (pPedigree->likelihoodPolynomial, pPedigree->likelihoodPolyList);
          codePoly (pPedigree->likelihoodPolynomial, pPedigree->likelihoodPolyList, polynomialFunctionName);
#ifdef POLYUSE_DL
          // Try to load the DL and ditch the tree polynomial
#ifdef POLYCHECK_DL
          // Squirrel-away the tree polynomial for later check
          pPedigree->cLikelihoodPolynomial = pPedigree->likelihoodPolynomial;
          holdPoly (pPedigree->cLikelihoodPolynomial);
          pPedigree->cLikelihoodPolyList = pPedigree->likelihoodPolyList;
#endif
#ifdef POLYCOMP_DL
          if ((pPedigree->likelihoodPolynomial = restoreExternalPoly (polynomialFunctionName)) == NULL)
            FATAL ("Couldn't load compiled likelihood polynomial that was just created!");
#else
#ifdef FAKEEVALUATE
          pPedigree->likelihoodPolynomial = constantExp (.05);
#endif
#endif
#endif
#endif
          // Notice we are normally holding only the external (compiled) poly!
          holdPoly (pPedigree->likelihoodPolynomial);
          freeKeptPolys ();
#ifdef POLYSTATISTICS
          if (i == pPedigreeList->numPedigree - 1)
            polyDynamicStatistics ("Post-build");
#endif
#ifdef POLYUSE_DL
        }
#endif
	// We still need to build a list even if there's only the external for the DL.
	pPedigree->likelihoodPolyList = buildPolyList ();
	polyListSorting (pPedigree->likelihoodPolynomial, pPedigree->likelihoodPolyList);
      }
    }
  }
  return EXIT_SUCCESS;
}


/* the main API to compute likelihood for all pedigrees in a data set */
int compute_likelihood (PedigreeSet * pPedigreeList)
{
  Pedigree *pPedigree;
  int i;
  int status;   /* return status of function calls */
  double product_likelihood = 1;        /* product of the likelihoods
                                         * for all the pedigrees */
  double sum_log_likelihood = 0;        /* sum of the
                                         * log10(likelihood) for all
                                         * the pedigrees */
  double log10Likelihood;
  int origLocus = 0;    /* locus index in the original locus list
                         * this is used to find out the pedigree
                         * counts mainly for case control analyses */
  int ret = 0;

  if (analysisLocusList->numLocus > 1)
    origLocus = analysisLocusList->pLocusIndex[1];
  numLocus = analysisLocusList->numLocus;

  /* Initialization */
  sum_log_likelihood = 0;
  product_likelihood = 1;
  pPedigreeList->likelihood = 1;
  pPedigreeList->log10Likelihood = 0;

#ifdef STUDYDB
  /* Here we can intercept all calls to compute_likelihood and suborn them to our purposes!

     If we're a client, we want to call GetLOD with all the appropriate likelhood variables.
     If we get a good LOD, then we just pass it on back up to the caller. If we didn't then
     the request is in and we pass back a dummy LOD and set a flag indicating that the
     calculation is to be ignored.

     If we're a server, we want to sign-in with full range of pedigree/positions that we can
     cover with the current set of markers, and start asking for work. When we're out of work
     we can return to the caller until we hit a new set of markers. A properly-configured
     config file can make this fast (only TM).

     Twopoint variables

     theta -- is both what is being varied by the driver and used by likelihood calculation.
     It is the recomination frequency between the trait position and the marker. Nothing
     fancy is required - we just store it as a client, and retrieve it as a server.

     Multipoint variables:

     trait position -- is what is being varied by the driver, and when we're a client, it
     has already been converted into inter-position theta values that are used for the
     likelihood calculation. Multipoint theta values specify distances between the trait 
     and two nearest markers from the set of relative trait/marker positions. E.g. when 
     for 3MP (2 markers  and 1 trait), we can have T-theta1-M-theta2-M, M-theta1-T-theta2-M, 
     and M-theta1-M-theta2-T. There are theta values between other markers, but they're
     constants.  The problem here lies in the LOD server situation. Execution will proceed
     to this chunk of code whereupon we will have a set of relative trait/marker
     positions for which we can calculate likelihoods, and a single (1) trait position and
     corresponding pair of thetas. We will ask the database for work within the span of
     positions covered by the current set of relative trait/marker positions, and get
     work requests in terms of the trait position only. This trait position will have to be
     converted to the two theta values that the likelihood calculation expects.

     penetrance vector -- 

     Where <liability-class> varies from 0-2, 
     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][<liability-class>][0][0] is pen_DD
     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][<liability-class>][0][1] is pen_Dd
     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][<liability-class>][1][0] is pen_dD
     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][<liability-class>][1][1] is pen_dd
     pTrait->penetrance[AFFECTION_STATUS_UNAFFECTED][<liability-class>][0][0] is 1 - pen_DD
     pTrait->penetrance[AFFECTION_STATUS_UNAFFECTED][<liability-class>][0][1] is 1 - pen_Dd
     pTrait->penetrance[AFFECTION_STATUS_UNAFFECTED][<liability-class>][1][0] is 1 - pen_dD
     pTrait->penetrance[AFFECTION_STATUS_UNAFFECTED][<liability-class>][1][1] is 1 - pen_dd

  */
  for (i = 0; i < pPedigreeList->numPedigree; i++) {
    pPedigree = pPedigreeList->ppPedigreeSet[i];
    GetDLOD (GetPedPosId (pPedigree->sPedigreeID, 40, KROUND(modelRange->tloc[dk_curModel.posIdx])),
	     dk_curModel.dgf, 
	     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][0][0][0],
	     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][0][0][1], 
	     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][0][1][0], 
	     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][0][1][1],
	     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][1][0][0],
	     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][1][0][1],
	     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][1][1][0],
	     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][1][1][1],
	     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][2][0][0],
	     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][2][0][1],
	     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][2][1][0],
	     pTrait->penetrance[AFFECTION_STATUS_AFFECTED][2][1][1],
	     1);
  }
#endif

  if (modelOptions->polynomial == TRUE) {

    /* Make sure they exist. Need to check until the day we have all builds separate, and
     * then we pull this to improve evaulation performance. */
    for (i = 0; i < pPedigreeList->numPedigree; i++) {
      pPedigree = pPedigreeList->ppPedigreeSet[i];
      if (pPedigree->load_flag == 0) {  // Skip everything if we're using saved results
        if (pPedigree->likelihoodPolynomial == NULL) {  // There's no polynomial, so come up with one
          build_likelihood_polynomial (pPedigreeList);
          break;        // Don't need to keep checking if we've done the build.
        }
      }
    }

    /* Now evaluate them all */
#ifdef _OPENMP
#pragma omp parallel for private(pPedigree)
#endif
    for (i = 0; i < pPedigreeList->numPedigree; i++) {
      pPedigree = pPedigreeList->ppPedigreeSet[i];
      if (pPedigree->load_flag == 0) {  // Skip everything if we're using saved results
#ifdef FAKEEVALUATE
        pPedigree->likelihood = .05;
#else
#ifdef MANYSMALLEVALUATE
        pPedigree->likelihood = evaluateValue (pPedigree->likelihoodPolynomial);
#else
        evaluatePoly (pPedigree->likelihoodPolynomial, pPedigree->likelihoodPolyList, &pPedigree->likelihood);
#endif
#ifdef POLYCHECK_DL
        // Check the result only if we actually just built a polynomial to check
        if (pPedigree->cLikelihoodPolynomial != NULL) {
          evaluatePoly (pPedigree->cLikelihoodPolynomial, pPedigree->cLikelihoodPolyList, &pPedigree->cLikelihood);
          if (fabs (pPedigree->likelihood - pPedigree->cLikelihood) > 1E-15)
            WARNING ("Discrepency between eV of %.*g and v of %.*g for %d\n", DBL_DIG, pPedigree->likelihood, DBL_DIG, pPedigree->cLikelihood, pPedigree->likelihoodPolynomial->id);
        }
#endif
#endif
      }
    }
  }
  /* Now (optional non-poly) and incorporate results */
  for (i = 0; i < pPedigreeList->numPedigree; i++) {
    pPedigree = pPedigreeList->ppPedigreeSet[i];
    if (pPedigree->load_flag == 0) {
      if (modelOptions->polynomial == FALSE) {
        initialize_multi_locus_genotype (pPedigree);
        status = compute_pedigree_likelihood (pPedigree);
      }
      if (modelOptions->dryRun == 0) {
        if (pPedigree->likelihood == 0.0) {
	  if (!modelRange->atypicalQtTrait)
	    WARNING ("Pedigree %s has likelihood %G thats too small.", pPedigree->sPedigreeID, pPedigree->likelihood);
          if (modelOptions->polynomial)
            DIAG (LIKELIHOOD, 1, {fprintf (stderr, "Polynomial terms (depth of 1):\n"); expTermPrinting (stderr, pPedigree->likelihoodPolynomial, 1); fprintf (stderr, "\n");});
          ret = -1;
          product_likelihood = 0.0;
          sum_log_likelihood = -9999.99;
          //break;
        } else if (pPedigree->likelihood < 0.0) {
          ASSERT (pPedigree->likelihood >= 0.0, "Pedigree %s has a NEGATIVE likelihood", pPedigree->sPedigreeID);
          product_likelihood = 0.0;
          sum_log_likelihood = -9999.99;
          ret = -2;
          break;
        } else {
          if (pPedigree->pCount[origLocus] == 1) {
            product_likelihood *= pPedigree->likelihood;
            log10Likelihood = log10 (pPedigree->likelihood);
          } else {
            product_likelihood *= pow (pPedigree->likelihood, pPedigree->pCount[origLocus]);
            log10Likelihood = log10 (pPedigree->likelihood) * pPedigree->pCount[origLocus];
          }
          sum_log_likelihood += log10Likelihood;
        }
      }
    }
  }

  if (ret < 0) {
    product_likelihood = 0.0;
    sum_log_likelihood = -9999.99;
    pPedigreeList->likelihood = product_likelihood;
    pPedigreeList->log10Likelihood = sum_log_likelihood;
  } else {
    pPedigreeList->likelihood = product_likelihood;
    pPedigreeList->log10Likelihood = sum_log_likelihood;
  }
  DIAG (LIKELIHOOD, 1, {fprintf (stderr, "Sum of log Likelihood is: %e\n", sum_log_likelihood);});
  return ret;
}

/* release polynomial for all pedigrees */
void pedigreeSetPolynomialClearance (PedigreeSet * pPedigreeList)
{

  Pedigree *pPedigree;
  int i;

  if (modelOptions->polynomial == TRUE) {
    for (i = 0; i < pPedigreeList->numPedigree; i++) {
      pPedigree = pPedigreeList->ppPedigreeSet[i];
      if (pPedigree->likelihoodPolynomial != NULL) {
        unHoldPoly (pPedigree->likelihoodPolynomial);
        if (pPedigree->likelihoodPolynomial->eType == T_EXTERNAL)
          releaseExternalPoly (pPedigree->likelihoodPolynomial);
        pPedigree->likelihoodPolynomial = NULL;
        free (pPedigree->likelihoodPolyList->pList);
        free (pPedigree->likelihoodPolyList);
        pPedigree->likelihoodPolyList = NULL;
      }
    }
    freeKeptPolys ();   /* Because holds overlapped keeps. */
#ifdef POLYSTATISTICS
    polyStatistics ("Post-pedigree free");
#endif
  }
}

/* compute likelihood for a given pedigree */
int compute_pedigree_likelihood (Pedigree * pPedigree)
{
  int i;
  NuclearFamily *pMyNucFam;     /* nuclear families within the pedigree */
  int status;   /* function return status */
  Person *pMyProband;   /* peeling proband */
  double likelihood = 0;
  double tmpLikelihood = 0;
  ConditionalLikelihood *pConditional;
  Polynomial *pLikelihoodPolynomial = NULL;
  int ret = 0;
  int origLocus = analysisLocusList->pLocusIndex[0];
  Genotype *pMyGenotype = NULL;
  //  int genoIdx = 0;
  Locus *pLocus = originalLocusList.ppLocusList[origLocus];
  int allele1, allele2;
  double sumFreq1, sumFreq2, freq1, freq2;
  AlleleSet *alleleSet1, *alleleSet2;
  double condL, condL2, sumCondL;
  int k, l, j, condIdx, idx;
  Person *pLoopBreaker;
  LoopBreaker *loopStruct;

  condIdx = 0;
  sumCondL = 0;
  if (modelOptions->dryRun == 0 && modelOptions->polynomial == TRUE)
    DETAIL (0, "Building polynomial w/pedigree: %s (%d/%d)", pPedigree->sPedigreeID, pPedigree->pedigreeIndex + 1, pPedigree->pPedigreeSet->numPedigree);

  for (i = 0; i < pPedigree->numNuclearFamily; i++) {
    pMyNucFam = pPedigree->ppNuclearFamilyList[i];
    pMyNucFam->totalNumPairGroups = 0;
    pMyNucFam->totalNumSimilarPairs = 0;
  }

  if (pPedigree->loopFlag) {
    populate_pedigree_loopbreaker_genotype_vector (pPedigree);
    ret = -2;
    while (ret == -2) {
      ret = set_next_loopbreaker_genotype_vector (pPedigree, TRUE);
      if (ret == -2)
        restore_pedigree_genotype_link_from_saved (pPedigree);
    }
  }
  if (modelOptions->polynomial == TRUE) {
    pLikelihoodPolynomial = constant0Poly;
  } else
    likelihood = 0;

  while (ret == 0) {
    initialize_multi_locus_genotype (pPedigree);
    /*
     * initialize all the nuclear families before peeling starts
     * in many cases, multiple likelihoods are computed for the
     * same family with different parameters, we need to clean up
     * before (or after) each calculation
     */
    for (i = 0; i < pPedigree->numNuclearFamily; i++) {
      pMyNucFam = pPedigree->ppNuclearFamilyList[i];
      pMyNucFam->doneFlag = 0;
      pMyNucFam->numPairGroups = 0;
      pMyNucFam->numSimilarPairs = 0;
    }

    /*
     * peeling starts from the peeling proband and eventually
     * will come back to the same proband this process will
     * obtain the conditional likelihoods for the proband
     */
    status = peel_graph (pPedigree->pPeelingNuclearFamily, pPedigree->pPeelingProband, pPedigree->peelingDirection);

    /*
     * done peeling, need to add up the conditional likelihood
     * for the leading peeling proband
     */
    pMyProband = pPedigree->pPeelingProband;
    if (modelOptions->polynomial != TRUE)
      tmpLikelihood = 0;

    pMyGenotype = pMyProband->ppProbandGenotypeList[origLocus];
    /* loop over all conditional likelihoods */
    for (i = 0; i < pMyProband->numConditionals; i++) {
      pConditional = &pMyProband->pLikelihood[i];
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
      if (modelOptions->polynomial == TRUE) {
        /* build likelihood polynomial */
        pLikelihoodPolynomial = plusExp (2, 1.0, pLikelihoodPolynomial, 1.0, timesExp (2, pConditional->lkslot.likelihoodPolynomial, 1, pConditional->wtslot.weightPolynomial, 1, 1), 1);
      } else {
        condL = pConditional->lkslot.likelihood * pConditional->wtslot.weight;
        tmpLikelihood += condL;
        sumCondL += condL;
        if (modelOptions->conditionalRun == 1) {
          /* find out the genotype index for the first locus - assuming that's the locus
           * we are interested in */
          //      genoIdx = i / pMyProband->pSavedNumGenotype[origLocus]; 
          alleleSet1 = pLocus->ppAlleleSetList[pMyGenotype->allele[DAD] - 1];
          alleleSet2 = pLocus->ppAlleleSetList[pMyGenotype->allele[MOM] - 1];
          sumFreq1 = alleleSet1->sumFreq;
          sumFreq2 = alleleSet2->sumFreq;
          for (k = 0; k < alleleSet1->numAllele; k++) {
            allele1 = alleleSet1->pAlleles[k];
            freq1 = pLocus->pAlleleFrequency[allele1 - 1];
            for (j = 0; j < alleleSet2->numAllele; j++) {
              allele2 = alleleSet2->pAlleles[j];
              freq2 = pLocus->pAlleleFrequency[allele2 - 1];
              if (condIdx >= numCond) {
                REALCHOKE (pCondSet, sizeof (probandCondL) * (numCond + 4), probandCondL *);
                numCond += 4;
              }
              condL2 = condL * (freq1 / sumFreq1) * (freq2 / sumFreq2);
              /* if heterzygous, find the matching one if there */
              for (idx = 0; idx < condIdx; idx++) {
                if ((pCondSet[idx].allele1 == allele2 && pCondSet[idx].allele2 == allele1) || (pCondSet[idx].allele1 == allele1 && pCondSet[idx].allele2 == allele2)) {
                  /* found the matching one */
                  condL2 += pCondSet[idx].condL;
                  /* remove it */
                  for (; idx < condIdx - 1; idx++) {
                    memcpy (&pCondSet[idx], &pCondSet[idx + 1], sizeof (probandCondL));
                  }
                  condIdx--;
                  break;
                }
              }

              for (idx = 0; idx < condIdx; idx++) {
                if (pCondSet[idx].condL <= condL2)
                  break;
              }
              for (l = condIdx; l > idx; l--) {
                memcpy (&pCondSet[l], &pCondSet[l - 1], sizeof (probandCondL));
              }
              pCondSet[idx].pPedigreeID = pPedigree->sPedigreeID;
              pCondSet[idx].pProbandID = pMyProband->sID;
              if (modelOptions->markerAnalysis == FALSE)
                pCondSet[idx].trait = pMyProband->ppTraitValue[0][0];
              pCondSet[idx].pAllele1 = pLocus->ppAlleleNames[allele1 - 1];
              pCondSet[idx].pAllele2 = pLocus->ppAlleleNames[allele2 - 1];
              pCondSet[idx].allele1 = allele1;
              pCondSet[idx].allele2 = allele2;
              pCondSet[idx].condL = condL2;
              condIdx++;
            }
          }
          // fflush(fpCond);
          pMyGenotype = pMyGenotype->pNext;
        }       /* end of conditional Run flag is on */
      } /* non polynomial mode */
    }   /* end of looping over all conditionals */

    if (modelOptions->polynomial != TRUE)
      likelihood += tmpLikelihood;

    if (pPedigree->loopFlag) {
      if (modelOptions->polynomial == TRUE)
        keepPoly (pLikelihoodPolynomial);
      else {
        DIAG (LIKELIHOOD, 1, {fprintf (stderr, "Log Likelihood for this fixed looped pedigree %s is: %e\n", pPedigree->sPedigreeID, log10 (tmpLikelihood));});
        condL = tmpLikelihood;
        if (modelOptions->loopCondRun == 1) {
          /* find the loop breaker we want */
          for (i = 0; i < pPedigree->numLoopBreaker; i++) {
            pLoopBreaker = pPedigree->loopBreakerList[i];
            if (strcmp (pLoopBreaker->sID, modelOptions->loopBreaker) != 0)
              continue;
            loopStruct = pLoopBreaker->loopBreakerStruct;
            pMyGenotype = loopStruct->genotype[loopStruct->genotypeIndex][0];
            alleleSet1 = pLocus->ppAlleleSetList[pMyGenotype->allele[DAD] - 1];
            alleleSet2 = pLocus->ppAlleleSetList[pMyGenotype->allele[MOM] - 1];
            sumFreq1 = alleleSet1->sumFreq;
            sumFreq2 = alleleSet2->sumFreq;
            for (k = 0; k < alleleSet1->numAllele; k++) {
              allele1 = alleleSet1->pAlleles[k];
              freq1 = pLocus->pAlleleFrequency[allele1 - 1];
              for (j = 0; j < alleleSet2->numAllele; j++) {
                allele2 = alleleSet2->pAlleles[j];
                freq2 = pLocus->pAlleleFrequency[allele2 - 1];
                if (condIdx >= numCond) {
                  pCondSet = (probandCondL *) realloc (pCondSet, sizeof (probandCondL) * (numCond + 4));
                  numCond += 4;
                }
                condL2 = condL * (freq1 / sumFreq1) * (freq2 / sumFreq2);
                /* find the matching one if there */
                for (idx = 0; idx < condIdx; idx++) {
                  if ((pCondSet[idx].allele1 == allele2 && pCondSet[idx].allele2 == allele1) || (pCondSet[idx].allele1 == allele1 && pCondSet[idx].allele2 == allele2)) {
                    /* found the matching one */
                    condL2 += pCondSet[idx].condL;
                    /* remove it */
                    for (; idx < condIdx - 1; idx++) {
                      memcpy (&pCondSet[idx], &pCondSet[idx + 1], sizeof (probandCondL));
                    }
                    condIdx--;
                    break;
                  }
                }

                /* find the right place to put it */
                for (idx = 0; idx < condIdx; idx++) {
                  if (pCondSet[idx].condL <= condL2)
                    break;
                }
                for (l = condIdx; l > idx; l--) {
                  memcpy (&pCondSet[l], &pCondSet[l - 1], sizeof (probandCondL));
                }
                pCondSet[idx].pPedigreeID = pPedigree->sPedigreeID;
                pCondSet[idx].pProbandID = pLoopBreaker->sID;
                pCondSet[idx].trait = pLoopBreaker->ppTraitValue[0][0];
                pCondSet[idx].pAllele1 = pLocus->ppAlleleNames[allele1 - 1];
                pCondSet[idx].pAllele2 = pLocus->ppAlleleNames[allele2 - 1];
                pCondSet[idx].allele1 = allele1;
                pCondSet[idx].allele2 = allele2;
                pCondSet[idx].condL = condL2;
                condIdx++;
              }
            }

          }
        }
      }

      ret = -2;
      while (ret == -2) {
        restore_pedigree_genotype_link_from_saved (pPedigree);
        ret = set_next_loopbreaker_genotype_vector (pPedigree, FALSE);
      }
    } else {
      ret = -1;
    }

    for (i = 0; i < pPedigree->numNuclearFamily; i++) {
      pMyNucFam = pPedigree->ppNuclearFamilyList[i];
      pMyNucFam->totalNumPairGroups += pMyNucFam->numPairGroups;
      pMyNucFam->totalNumSimilarPairs += pMyNucFam->numSimilarPairs;
    }
  }
  if (modelOptions->polynomial == TRUE)
    /* save the polynomial to the pedigree structure */
    pPedigree->likelihoodPolynomial = pLikelihoodPolynomial;
  else {
    /* save the likelihood in the pedigree structure */
    pPedigree->likelihood = likelihood;
    DIAG (LIKELIHOOD, 1, {fprintf (stderr, "log Likelihood for pedigree %s is: %e\n", pPedigree->sPedigreeID, log10 (likelihood));});
    //    if ((modelOptions->loopCondRun == 1 || modelOptions->conditionalRun == 1)
    //        && modelOptions->polynomial != TRUE) {
      //      for (k = 0; k < condIdx; k++) {
      //        fprintf (fpCond, "%s %s %d %s %s %e %5.1f%%\n", pCondSet[k].pPedigreeID, pCondSet[k].pProbandID, (int) pCondSet[k].trait, pCondSet[k].pAllele1, pCondSet[k].pAllele2, pCondSet[k].condL, pCondSet[k].condL / sumCondL * 100);
      //      }
    //    }
  }

  return 0;
}

/*
 * recursive procedure to go through all nuclear families within one pedigree
 * pNucFam -- input nuclear family, the top layer is the peeling nuclear
 * family which is the nuclear family contains the peeling proband pProband
 * -- connector peelingDireciton -- UP/DOWN -- currently it is not used at
 * all
 */
int peel_graph (NuclearFamily * pNucFam1, Person * pProband1, int peelingDirection)
{
  NuclearFamilyConnector *pConnector;   /* connector structure which
                                         * represents connector from one nuc
                                         * to the other within a pedigree */
  NuclearFamily *pNucFam2;      /* another nuc family input nuc family is
                                 * connected to */
  Person *pConnectPerson;       /* connector individual */
  int i;
  ConditionalLikelihood *pConditional;

  /*
   * if this nuclear family has been processed or in the middle of
   * process, skip it
   */
  if (pNucFam1->doneFlag == TRUE)
    return 0;

  swLogProgress(3, 0, "Start peeling nuclear family %d with parents %s x %s",
		pNucFam1->nuclearFamilyIndex, pNucFam1->pParents[DAD]->sID, pNucFam1->pParents[MOM]->sID);

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
    if (pConnectPerson == pNucFam2->pParents[DAD] || pConnectPerson == pNucFam2->pParents[MOM]) {
      /*
       * these two families are connected through multiple
       * marraige
       */
      peel_graph (pConnector->pConnectedNuclearFamily, pConnector->pConnectedPerson, PEDIGREE_UP);
    } else {
      peel_graph (pConnector->pConnectedNuclearFamily, pConnector->pConnectedPerson, PEDIGREE_DOWN);
    }
    pConnector = pConnector->pNextConnector;
  }

  /* go down through the down connectors if any */
  pConnector = pNucFam1->pDownConnectors;
  while (pConnector) {
    /* peel up to the proband */
    peel_graph (pConnector->pConnectedNuclearFamily, pConnector->pConnectedPerson, PEDIGREE_UP);

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

  memcpy (&pProband->ppProbandGenotypeList[0], &pProband->ppGenotypeList[0], sizeof (Genotype *) * originalLocusList.numLocus);
  memcpy (&pProband->pProbandNumGenotype[0], &pProband->pNumGenotype[0], sizeof (int) * originalLocusList.numLocus);

  DIAG (LIKELIHOOD, 1, {fprintf (stderr, "\t Proband (%s) haplotype: \n", pProband->sID); });
  if (pProband->ppHaplotype == NULL) {
    /*
     * allocate space for storing proband's haplotype if not
     * already done
     */
    MALCHOKE (pProband->ppHaplotype, sizeof (Genotype *) * sizeof (int) * ppairMatrixNumLocus, void *);
  }
  if (pProband->touchedFlag == TRUE)
    initialize_proband_tmpLikelihood (pProband);

  if (pNucFam->pParents[DAD] != pProband && pNucFam->pParents[MOM] != pProband) {
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
  } else {      /* A parent is the proband */
    compute_nuclear_family_likelihood (peelingDirection);
    if (modelOptions->dryRun == 0) {
      /* copy the touched temporary results back */
      for (i = 0; i < pProband->numTmpLikelihood; i++) {
        pConditional = &pProband->pLikelihood[pProband->pTmpLikelihoodIndex[i]];
        //      if (pConditional->tmpTouched == FALSE)
        //        continue;
        pConditional->touchedFlag = TRUE;
        if (modelOptions->polynomial == TRUE) {
          if (pProband->touchedFlag == FALSE) {
            pConditional->lkslot.likelihoodPolynomial = constant1Poly;
          }
          pConditional->lkslot.likelihoodPolynomial = timesExp (2, pConditional->lkslot.likelihoodPolynomial, 1, pConditional->tmpslot.tmpLikelihoodPolynomial, 1, 0);
          pConditional->tmpslot.tmpLikelihoodPolynomial = constant0Poly;
#if 0
          DIAG (LIKELIHOOD, 1, {fprintf (stderr, "Proband %s Conditional Likelihood (%d) = %e. Weight = %e\n",
					 pProband->sID, i, evaluateValue (pConditional->lkslot.likelihoodPolynomial), evaluateValue (pConditional->wtslot.weightPolynomial));});
#endif
        } else {        /* PE is not enabled */
          if (pProband->touchedFlag == FALSE)
            pConditional->lkslot.likelihood = 1;
          pConditional->lkslot.likelihood *= pConditional->tmpslot.tmpLikelihood;
          pConditional->tmpslot.tmpLikelihood = 0;
          DIAG (LIKELIHOOD, 1, {fprintf (stderr, "Proband %s Conditional Likelihood (%d) = %e. Weight = %e\n", pProband->sID, i, pConditional->lkslot.likelihood, pConditional->wtslot.weight);});
        }
        pConditional->tmpTouched = FALSE;
      }
      pProband->numTmpLikelihood = 0;
    }
  }

  DIAG (LIKELIHOOD, 1, {fprintf (stderr, "Nuclear Family %d with parents %s x %s.\n", pNucFam->nuclearFamilyIndex, pNucFam->pParents[DAD]->sID, pNucFam->pParents[MOM]->sID);});

  /*
   * mark the proband as been touched - we have done some likelihood
   * calculation on this person
   */
  pProband->touchedFlag = TRUE;

  /* copy back the genotypes for the proband */
  memcpy (&pProband->ppGenotypeList[0], &pProband->ppProbandGenotypeList[0], sizeof (Genotype *) * originalLocusList.numLocus);
  memcpy (&pProband->pNumGenotype[0], &pProband->pProbandNumGenotype[0], sizeof (int) * originalLocusList.numLocus);

  if (modelOptions->polynomial == TRUE) {
    for (i = 0; i < pProband->numConditionals; i++) {
      if (pProband->touchedFlag == TRUE) {
        keepPoly (pProband->pLikelihood[i].lkslot.likelihoodPolynomial);
        keepPoly (pProband->pLikelihood[i].wtslot.weightPolynomial);
      }
    }
    freePolys ();
  }

  swLogProgress(3, 0, "Finished peeling nuclear family %d with parents %s x %s",
		pNucFam->nuclearFamilyIndex, pNucFam->pParents[DAD]->sID, pNucFam->pParents[MOM]->sID);

  return 0;
}

/*
 * This function loops over child proband's multilocus genotypes recursively
 * and calls compute_nuclear_family_likelihood for each multilocus genotype
 * pNucFam - the nuclear family we are working on locus - current working
 * locus multiLocusIndex - index to the flattened array of conditional
 * likelihood for the proband
 */
int loop_child_proband_genotype (int peelingDirection, int locus, int multiLocusIndex)
{
  int origLocus = analysisLocusList->pLocusIndex[locus];        /* locus index in the
                                                         * original locus list */
  Genotype *pMyGenotype;
  int position; /* genotype position */
  Genotype *pNextGenotype;
  int numGenotype;      /* number of possible genotypes for this
                         * proband at this locus */
  int multiLocusIndex2;
  int traitLocus;
  ConditionalLikelihood *pConditional;
  double penetrance = 1;

  Polynomial *penetrancePolynomial = NULL;

  /*
   * we loop over the genotypes of the proband to condition the
   * likelihood calculation on it
   */
  numGenotype = pProband->pSavedNumGenotype[origLocus];
  /* calculate the flattened conditional likelihood array index */
  multiLocusIndex2 = multiLocusIndex * numGenotype;
  pMyGenotype = pProband->ppProbandGenotypeList[origLocus];
  while (pMyGenotype != NULL) {
    /*
     * record this locus's genotype in the haplotype structure -
     * it's really just phased multilocus genotype (not single
     * chromosome haplotype)
     */
    pProband->ppHaplotype[locus] = pMyGenotype;
    DIAG (LIKELIHOOD, 1, {fprintf (stderr, "\t proband (%s) %d|%d \n", pProband->sID, pMyGenotype->allele[DAD], pMyGenotype->allele[MOM]);});
    /*
     * temporarilly set the next pointer to NULL so to restrict
     * the genotype on the proband to current genotype only
     */
    pProband->ppGenotypeList[origLocus] = pMyGenotype;
    pNextGenotype = pMyGenotype->pNext;
    pMyGenotype->pNext = NULL;
    pProband->pNumGenotype[origLocus] = 1;
    position = pMyGenotype->position;
    /* calculate the flattened conditional likelihood array index */
    multiLocusIndex = multiLocusIndex2 + position;

    if (locus < analysisLocusList->numLocus - 1) {
      loop_child_proband_genotype (peelingDirection, locus + 1, multiLocusIndex);
    } else {
      /*
       * we have got the entire multi-locus genotype for
       * the proband
       */
      compute_nuclear_family_likelihood (peelingDirection);
      pConditional = &pProband->pLikelihood[multiLocusIndex];
      pConditional->touchedFlag = TRUE;
      if (modelOptions->dryRun == 0) {
        /*
         * store the likelihood in the corresponding
         * flattened array
         */
        if (pProband->touchedFlag == FALSE) {
          /*
           * if trait locus exists, we need to retrieve
           * the penetrance
           */
          traitLocus = analysisLocusList->traitLocusIndex;
          if (traitLocus >= 0) {
            if (modelOptions->polynomial == TRUE) {
              penetrancePolynomial = pProband->ppHaplotype[traitLocus]->penslot.penetrancePolynomial;
            } else {
              penetrance = pProband->ppHaplotype[traitLocus]->penslot.penetrance;
            }
          } else {      /* no trait locus */
            if (modelOptions->polynomial == TRUE) {
              penetrancePolynomial = constant1Poly;
            } else {
              penetrance = 1;
            }
          }     /* end of trait locus processing */

          if (modelOptions->polynomial == TRUE) {
            pConditional->lkslot.likelihoodPolynomial = timesExp (2, penetrancePolynomial, 1, pNucFam->likelihoodPolynomial, 1, 0);

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
            pConditional->lkslot.likelihood = penetrance * pNucFam->likelihood;
            pConditional->wtslot.weight = 1;    // For a child, the weight should be 1
          }
        }
        /*
         * first time updating the likelihood for this phased
         * multilocus genotype
         */
        else {  /* NOT first time updating the likelihood for
                 * this phased multilocus genotypes */
          /* no need to consider penetrance anymore */
          if (modelOptions->polynomial == TRUE) {
            pConditional->lkslot.likelihoodPolynomial = timesExp (2, pConditional->lkslot.likelihoodPolynomial, 1, pNucFam->likelihoodPolynomial, 1, 1);
#if 0
            DIAG (LIKELIHOOD, 2, {
                  fprintf (stderr, "Likelihood for this entire multi-locus genotype %f %f\n",
                      evaluateValue (pNucFam->likelihoodPolynomial),
                      evaluateValue (pProband->pLikelihood[multiLocusIndex].likelihoodPolynomial));
                  fprintf (stderr, "Proband %s Conditional Likelihood (%d) = %e.\n", pProband->sID, multiLocusIndex, evaluateValue (pConditional->lkslot.likelihoodPolynomial));
                }
            );
#endif
          } else {
            pConditional->lkslot.likelihood *= pNucFam->likelihood;
            DIAG (LIKELIHOOD, 1, {fprintf (stderr, "Proband %s Conditional Likelihood (%d) = %e.\n", pProband->sID, multiLocusIndex, pConditional->lkslot.likelihood);});
          }
        }
      }
    }   /* end of processing complete child
         * multilocus genotype  */

    /*
     * when we are done, need to restore the genotype link list
     * pointer
     */
    pMyGenotype->pNext = pNextGenotype;
    pMyGenotype = pNextGenotype;
  }     /* loop over all possible genotypes */

  return 0;
}

/*
 * for now, we only use parental pair algorithm This function goes through
 * all possible parental pairs, for each parental pair, it calculates
 * likelihood
 */
int compute_nuclear_family_likelihood (int peelingDirection)
{
  int locus;    /* locus index to construct parental pairs */
  double weight[2] = { 1, 1 };  /* weight for the two parents */
  /* need to define some terms for the polynomail operations */
  Polynomial *weightPolynomial[2];
  int numHaplotypePair = 1;     /* number of multilocus genotypes */
  int numChild; /* number of children in this nuclear family */
  int i, j;
  int multiLocusIndex[2] = { 0, 0 };

  /* initialize the weight for each parent */
  if (modelOptions->polynomial == TRUE) {
    weightPolynomial[0] = constant1Poly;
    weightPolynomial[1] = constant1Poly;
    pNucFam->likelihoodPolynomial = constant0Poly;
  }
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
  for (locus = 0; locus < analysisLocusList->numLocus; locus++) {
    /* construct parental pair locus by locus */
    construct_parental_pair (pNucFam, pProband, locus);
    /* calculate number of possible multilocus genotypes */
    numHaplotypePair *= parentalPairSpace.pNumParentalPair[locus];
  }

  for (i = DAD; i <= MOM; i++) {
    pNucFam->numHetLocus[i] = 0;
    pNucFam->firstHetLocus[i] = -1;
    for (j = 0; j < analysisLocusList->numLocus; j++)
      pNucFam->hetFlag[i][j] = 0;
  }

  /* now we can construct haplotypes and get likelihood computed */
  DIAG (LIKELIHOOD, 1, {
      fprintf (stderr, "Haplotype for nuclear family No. %d:\n", pNucFam->nuclearFamilyIndex);
      fprintf (stderr, "\t\t\t DAD(%s)\t\t\t MOM(%s)\n", pNucFam->pParents[DAD]->sID, pNucFam->pParents[MOM]->sID);
      }
  );
  /*
   * recursively call loop_parental_pair to get a complete multlocus
   * genotype
   */
  if (modelOptions->polynomial == TRUE)
    loop_parental_pair (0, multiLocusIndex, (void *) weightPolynomial);
  else
    loop_parental_pair (0, multiLocusIndex, (void *) weight);

  return 0;
}

/*
 * using the parental pair algorithm to go through all possible parental
 * pairs and calculate the likelihood of each parental pair nested looping is
 * for the multi-locus
 */
int loop_parental_pair (int locus, int multiLocusIndex[2], void *dWeight[2])
{
  ParentalPair *pPair;  /* parental pair for current locus */
  int origLocus;        /* locus index in the original locus list */
  int i, j, k;
  int multiLocusIndex2[2];
  Person *pParent[2];   /* two parents */
  int numGenotype[2];   /* number of genotypes at this locus */
  int numPair;  /* index to the list of parental pairs */
  double newWeight[2];  /* genotype weight */

  Polynomial *newWeightPolynomial[2];
  int head;     /* proband if no child is a proband,
                 * otherwise DAD  */
  int spouse;   /* spouse of the head */
  int likelihoodIndex;  /* likelihood index for the proband */
  int multiLocusPhase[2] = { 0, 0 };    /* index to the related
                                         * parental pair matrix */
  int flipMask[2] = { 0, 0 };
  ConditionalLikelihood *pConditional;
  PPairElement *pElement;


  head = pNucFam->head;
  spouse = pNucFam->spouse;

  origLocus = analysisLocusList->pLocusIndex[locus];
  for (i = DAD; i <= MOM; i++) {
    if (modelOptions->polynomial == TRUE)
      newWeightPolynomial[i] = (Polynomial *) dWeight[i];
    else
      newWeight[i] = *((double *) dWeight + i);
    pParent[i] = pNucFam->pParents[i];
    /* find the max number of possible genotypes for this parent */
    if (pParent[i]->loopBreaker >= 1 && pParent[i]->pParents[DAD] == NULL)
      numGenotype[i] = pParent[i]->pOriginalPerson->pSavedNumGenotype[origLocus];
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
    multiLocusIndex2[DAD] = multiLocusIndex[DAD] * numGenotype[DAD] + pPair->pGenotype[DAD]->position;
    multiLocusIndex2[MOM] = multiLocusIndex[MOM] * numGenotype[MOM] + pPair->pGenotype[MOM]->position;
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
    while ((numPair + 1) < pHaplo->pNumParentalPair[locus] && (pHaplo->ppParentalPair[locus][numPair + 1].phase[DAD] > 0 || pHaplo->ppParentalPair[locus][numPair + 1].phase[MOM] > 0)) {
      numPair++;
      pNucFam->numRelatedPPair[locus]++;
    }
    if (locus == 0) {
      pNucFam->totalRelatedPPair[locus] = pNucFam->numRelatedPPair[0];
    } else {
      pNucFam->totalRelatedPPair[locus] = pNucFam->totalRelatedPPair[locus - 1] * pNucFam->numRelatedPPair[locus];
    }

    for (i = DAD; i <= MOM; i++) {
      if (modelOptions->equilibrium == LINKAGE_EQUILIBRIUM && pParent[i]->pParents[DAD] == NULL) {
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
        if (modelOptions->polynomial == TRUE) {
          newWeightPolynomial[i] = timesExp (2, (Polynomial *) dWeight[i], 1, pPair->pGenotype[i]->wtslot.weightPolynomial, 1, 0);
        } else {
          newWeight[i] = *((double *) dWeight + i) * pPair->pGenotype[i]->wtslot.weight;
        }
      }
    }   /* looping dad and mom genotypes */


    if (locus < analysisLocusList->numLocus - 1) {
      /*
       * recursively calling this function to get a
       * complete multilocus genotype
       */
      if (modelOptions->polynomial == TRUE) {

        loop_parental_pair (locus + 1, multiLocusIndex2, (void *) newWeightPolynomial);
      } else {
        loop_parental_pair (locus + 1, multiLocusIndex2, (void *) newWeight);
      }
    } else {    /* got a complete set of parental pairs */
      pNucFam->numPairGroups++;
      pNucFam->numSimilarPairs += pNucFam->totalRelatedPPair[locus] - 1;
      if (modelOptions->dryRun != 0) {
        /* no actual calculations under dry run, continue to next set of ppair */
        continue;
      }

      if (pNucFam->totalRelatedPPair[locus] == 1) {
        multiLocusPhase[DAD] = 0;
        multiLocusPhase[MOM] = 0;
        calcFlag = 0;
        if (modelOptions->polynomial == TRUE)
          calculate_likelihood (multiLocusIndex2, multiLocusPhase, (void *) newWeightPolynomial, NULL);
        else
          calculate_likelihood (multiLocusIndex2, multiLocusPhase, (void *) newWeight, NULL);

        /* when there is only one parental pair, the likelihood is saved at the first cell,
         * as the phase index passed in was 0,0
         */
        pElement = &ppairMatrix[0][0];
        if (pNucFam->childProbandFlag == TRUE) {
          if (modelOptions->polynomial == TRUE) {
            pNucFam->likelihoodPolynomial = plusExp (2, 1.0, pNucFam->likelihoodPolynomial, 1.0, pElement->slot.likelihoodPolynomial, 1);
          } else
            pNucFam->likelihood += pElement->slot.likelihood;
        } else {        /* one of the parent is the proband */
          likelihoodIndex = multiLocusIndex2[head];
          pConditional = &pProband->pLikelihood[likelihoodIndex];
          if (modelOptions->polynomial == TRUE) {
            pConditional->tmpslot.tmpLikelihoodPolynomial = plusExp (2, 1.0, pConditional->tmpslot.tmpLikelihoodPolynomial, 1.0, pElement->slot.likelihoodPolynomial, 1);
          } else
            pConditional->tmpslot.tmpLikelihood += pElement->slot.likelihood;
          if (pConditional->tmpTouched == FALSE) {
            pConditional->tmpTouched = TRUE;
            pProband->pTmpLikelihoodIndex[pProband->numTmpLikelihood] = likelihoodIndex;
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
        if (modelOptions->polynomial == TRUE)
          loop_phases (0, multiLocusIndex2, multiLocusPhase, flipMask, (void *) newWeightPolynomial);
        else
          loop_phases (0, multiLocusIndex2, multiLocusPhase, flipMask, (void *) newWeight);

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
                if (modelOptions->polynomial == TRUE) {
                  pNucFam->likelihoodPolynomial = plusExp (2, 1.0, pNucFam->likelihoodPolynomial, 1.0, timesExp (2, pElement->slot.likelihoodPolynomial, 1, constantExp (pElement->count), 1, 1), 1);
                } else
                  pNucFam->likelihood += pElement->slot.likelihood * pElement->count;
              } else if (pElement->count > 0) { /* count == 1 */
                if (modelOptions->polynomial == TRUE) {
                  pNucFam->likelihoodPolynomial = plusExp (2, 1.0, pNucFam->likelihoodPolynomial, 1.0, pElement->slot.likelihoodPolynomial, 1);
                } else
                  pNucFam->likelihood += pElement->slot.likelihood;
              }
              /* end of count = 1 */
              /* reset count to 0 */
              pElement->count = 0;
            }
          }
        } else {        /* one of the parent is the proband */
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
                if (modelOptions->polynomial == TRUE) {
                  pConditional->tmpslot.tmpLikelihoodPolynomial = plusExp (2, 1.0, pConditional->tmpslot.tmpLikelihoodPolynomial, 1.0, timesExp (2, pElement->slot.likelihoodPolynomial, 1, constantExp (pElement->count), 1, 0), 1);
                } else
                  pConditional->tmpslot.tmpLikelihood += pElement->slot.likelihood * pElement->count;
              }
              /* count > 1 */
              else {    /* count==1 */
                if (modelOptions->polynomial == TRUE) {
                  pConditional->tmpslot.tmpLikelihoodPolynomial = plusExp (2, 1.0, pConditional->tmpslot.tmpLikelihoodPolynomial, 1.0, pElement->slot.likelihoodPolynomial, 1);
                } else
                  pConditional->tmpslot.tmpLikelihood += pElement->slot.likelihood;
              } /* count > 0 */
              /* reset the count */
              pElement->count = 0;
            }   /* loop through column */
          }     /* loop through row */
        }       /* one parent is proband */
      }
    }   /* end of processing one parental pair */
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
void loop_phases (int locus, int multiLocusIndex[2], int multiLocusPhase[2], int flipMask[2], void *dWeight[2])
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
  int origLocus = analysisLocusList->pLocusIndex[locus];        /* locus index in the
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

  if (modelOptions->polynomial == TRUE)
    childProductPtr = &childProductPoly;
  else
    childProductPtr = &childProduct;

  /*
   * double newWeight[2];
   * double sum;
   * int child;
   * int xmissionIndex[2];
   * double childProduct;
   * double penetrance[2];
   * int traitLocus;
   * int genoIndex;
   * #ifndef NO_POlYNOMIAL
   * Polynomial *newWeightPolynomial[2];
   * Polynomial *penetrancePolynomial[2];
   * Polynomial *sumPolynomial;
   * Polynomial *childProductPolynomial = NULL;
   * #endif
   */

  for (i = DAD; i <= MOM; i++) {
    pParent[i] = pNucFam->pParents[i];
    /* find the max number of possible genotypes for this parent */
    if (pParent[i]->loopBreaker >= 1 && pParent[i]->pParents[DAD] == NULL)
      numGenotype[i] = pParent[i]->pOriginalPerson->pSavedNumGenotype[origLocus];
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
    DIAG (LIKELIHOOD, 1, {fprintf (stderr, "(%s) %2d->\t %2d|%-2d --X-- %2d|%-2d  (%s)\n",
              pNucFam->pParents[DAD]->sID, origLocus, pPair->pGenotype[0]->allele[DAD], pPair->pGenotype[0]->allele[MOM], pPair->pGenotype[1]->allele[DAD], pPair->pGenotype[1]->allele[MOM], pNucFam->pParents[MOM]->sID);
        }
    );
    for (i = DAD; i <= MOM; i++) {
      pHaplo->phase[i][locus] = pPair->phase[i];
      multiLocusPhase2[i] = multiLocusPhase[i] | pPair->phase[i];
      if (pPair->phase[i] == 0)
        newFlipMask[i] = flipMask[i];
      else
        newFlipMask[i] = flipMask[i] | 3;
      multiLocusIndex2[i] = multiLocusIndex[i] + pPair->pGenotype[i]->position;
    }
    if (locus < analysisLocusList->numLocus - 1) {
      loop_phases (locus + 1, multiLocusIndex2, multiLocusPhase2, newFlipMask, dWeight);
    } else {    /* got a complete multilocus parental pair */
      calculateFlag = 1;
      ppairMatrix[multiLocusPhase2[proband]][multiLocusPhase2[spouse]].likelihoodIndex = multiLocusIndex2[proband];
      if (pNucFam->childProbandFlag == TRUE) {  /* a child is the proband */
        phase[DAD] = multiLocusPhase2[DAD];
        phase[MOM] = multiLocusPhase2[MOM];
        /* child is the proband */
        for (i = DAD; i <= MOM; i++) {
          if (modelOptions->imprintingFlag != TRUE && pNucFam->firstHetLocus[i] >= 0 && pHaplo->phase[i][pNucFam->firstHetLocus[i]] != 0) {
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
              multiLocusPhaseFlip[i] = multiLocusPhase2[i] ^ pNucFam->hetLocusBits[i];
              phase[i] = multiLocusPhaseFlip[i];
              if (ppairMatrix[phase[proband]][phase[spouse]].count > 0)
                calculateFlag = 0;
            }
          }
        }
        if (calculateFlag == 0) {
          ppairMatrix[phase[proband]][phase[spouse]].count++;
          if (modelOptions->polynomial == TRUE) {
#if 0
            DIAG (LIKELIHOOD, 1, {
                  fprintf (stderr, "\t\t likelihood (%d) = %e\n", ppairMatrix[phase[proband]][phase[spouse]].likelihoodIndex, evaluateValue (ppairMatrix[phase[proband]]
                          [phase[spouse]].slot.likelihoodPolynomial));
                }
            );
#endif
          } else
            DIAG (LIKELIHOOD, 1, {
                fprintf (stderr, "\t\t likelihood (%d) = %e\n", ppairMatrix[phase[proband]][phase[spouse]].likelihoodIndex, ppairMatrix[phase[proband]][phase[spouse]].slot.likelihood);
                }
          );
        }
      } else {  /* proband is a parent */
        if (modelOptions->imprintingFlag != TRUE && pNucFam->firstHetLocus[spouse] >= 0 && pHaplo->phase[spouse][pNucFam->firstHetLocus[spouse]] != 0) {
          if (pNucFam->pParents[spouse]->pParents[DAD] == NULL) {
            /*
             * non-proband parent is a
             * founder
             */
            /* find the reverse pattern */
            multiLocusPhaseFlip[spouse] = multiLocusPhase2[spouse] ^ pNucFam->hetLocusBits[spouse];
            if (ppairMatrix[multiLocusPhase2[proband]]
                [multiLocusPhaseFlip[spouse]].count > 0) {

              /*
               * increase the count
               * on the original
               */
              ppairMatrix[multiLocusPhase2[proband]]
                  [multiLocusPhaseFlip[spouse]].count++;
              if (modelOptions->polynomial == FALSE)
                DIAG (LIKELIHOOD, 1, {
                    fprintf (stderr, "\t\t likelihood (%d) = %e\n", ppairMatrix[multiLocusPhase2[proband]]
                        [multiLocusPhaseFlip[spouse]].likelihoodIndex, ppairMatrix[multiLocusPhase2[proband]]
                        [multiLocusPhaseFlip[spouse]].slot.likelihood);
                    }
              );
              calculateFlag = 0;
            }
          }
        } else if (modelOptions->imprintingFlag != TRUE && pNucFam->firstHetLocus[proband] >= 0 && pHaplo->phase[proband][pNucFam->firstHetLocus[proband]] != 0) {
          /* find the reverse pattern */
          multiLocusPhaseFlip[proband] = multiLocusPhase2[proband] ^ pNucFam->hetLocusBits[proband];
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
            if (modelOptions->polynomial == TRUE) {
              ppairMatrix[multiLocusPhase2[proband]]
                  [multiLocusPhase2[spouse]].slot.likelihoodPolynomial = ppairMatrix[multiLocusPhaseFlip[proband]]
                  [multiLocusPhase2[spouse]].slot.likelihoodPolynomial;
              pProband->pLikelihood[multiLocusIndex2[proband]].wtslot.weightPolynomial = pProband->pLikelihood[likelihoodIndex].wtslot.weightPolynomial;
            } else {
              ppairMatrix[multiLocusPhase2[proband]]
                  [multiLocusPhase2[spouse]].slot.likelihood = ppairMatrix[multiLocusPhaseFlip[proband]]
                  [multiLocusPhase2[spouse]].slot.likelihood;
              pProband->pLikelihood[multiLocusIndex2[proband]].wtslot.weight = pProband->pLikelihood[likelihoodIndex].wtslot.weight;
              DIAG (LIKELIHOOD, 1, {
                    fprintf (stderr, "\t\t likelihood (%d) = %e\n", ppairMatrix[multiLocusPhase2[proband]]
                        [multiLocusPhase2[spouse]].likelihoodIndex, ppairMatrix[multiLocusPhase2[proband]]
                        [multiLocusPhase2[spouse]].slot.likelihood);
                  }
              );
            }
            calculateFlag = 0;
          }
        }       /* end of finding patterns */
      } /* end of proband is a parent */

      if (calculateFlag == 1) {
        /* 0 - don't keep result, 1 - keep result, 2 - use result */
        if (modelOptions->imprintingFlag == TRUE) {
          calcFlag = 0;
        } else {
          if (multiLocusPhase2[DAD] == 0 && multiLocusPhase2[MOM] == 0) {
            memset (likelihoodChildCount, 0, sizeof (int) * maxChildren);
            calcFlag = 1;
          } else {
            calcFlag = 2;
            recalculate_child_likelihood (newFlipMask, childProductPtr);
          }
        }
        calculate_likelihood (multiLocusIndex2, multiLocusPhase2, dWeight, childProductPtr);
      } /* fresh likelihood calculation */
    }   /* end of processing a complete multilocus
         * parental pair */

  }     /* move on to next pair on this locus */

}       /* end of loop_phases() */


/* calculating childProduct base on previous pattern */
void recalculate_child_likelihood (int flipMask[2], void *childProduct)
{
  int i, j;
  double myChildSum = 0;
  Polynomial *childSumPoly = NULL;
  ChildElement *pElement;
  int xmissionIndex[2];

  multCount = 0;
  if (modelOptions->polynomial == TRUE)
    *(Polynomial **) childProduct = constant1Poly;
  else
    *(double *) childProduct = 1;

  for (i = 0; i < pNucFam->numChildren; i++) {
    if (modelOptions->polynomial == TRUE)
      childSumPoly = constant0Poly;
    else
      myChildSum = 0;
    for (j = 0; j < likelihoodChildCount[i]; j++) {
      pElement = &likelihoodChildElements[multCount + j];
      for (parent = DAD; parent <= MOM; parent++) {
        xmissionIndex[parent] = pElement->xmissionIndex[parent] ^ flipMask[parent];
      }

      if (modelOptions->polynomial == TRUE) {
        childSumPoly = plusExp (2, 1.0, childSumPoly, 1.0, timesExp (3, xmissionMatrix[xmissionIndex[DAD]].slot.probPoly[1], 1, xmissionMatrix[xmissionIndex[MOM]].slot.probPoly[2], 1, pElement->fslot.factorPolynomial, 1, 0), 1);
      } else {
        myChildSum += xmissionMatrix[xmissionIndex[DAD]].slot.prob[1] * xmissionMatrix[xmissionIndex[MOM]].slot.prob[2] * pElement->fslot.factor;
      }
    }
    if (modelOptions->polynomial == TRUE)
      *(Polynomial **) childProduct = timesExp (2, *(Polynomial **) childProduct, 1, childSumPoly, 1, 1);
    else
      *(double *) childProduct *= myChildSum;
    multCount += likelihoodChildCount[i];

  }
}

/**
   Basic likelihood calculation for a nuclear family conditional on one parental pair.


   Globals Referenced:

   - modelOptions
     - .polynomial
   - calcFlag

   - pParent

   Globals Changed:

   - childSum

   - ppairMatrix, Output, indexed by phases of proband and spouse.
     Changed members:
     - .likelihoodIndex gets multiLocusIndex of the proband.
     - .count gets set to 1.
     In polynomial mode, changed members:
     - .slot.likelihoodPolynomial gets polynomial product of childPolynomial and weights and penetrances of parents.
     In non-polynomial mode, changed members:
     - .slot.likelihood gets product of childProduct and weights and penetrances of parents.

*/

int calculate_likelihood (int multiLocusIndex[2],       ///< Input, index into pParent[i]->pLikelihood
    int multiLocusPhase[2],     ///< Input, proband and spouse phases used as index into ppairMatrix
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

  Polynomial *newWeightPolynomial[2] = { NULL, NULL };
  Polynomial *penetrancePolynomial[2] = { NULL, NULL };
  Polynomial *childProductPolynomial = NULL;
  Polynomial *sumPolynomial = NULL;
  ConditionalLikelihood *pConditional;

  pParent[DAD] = pNucFam->pParents[DAD];
  pParent[MOM] = pNucFam->pParents[MOM];
  proband = pNucFam->head;
  spouse = pNucFam->spouse;
  int xmissionIndex[2] = { 0, 0 };

  if (modelOptions->polynomial == TRUE) {
    childSum = &sumPolynomial;
    if (calcFlag == 2)  // If we're actually using the result
      childProductPolynomial = *(Polynomial **) childProductPtr;
  } else {
    childSum = &sum;
    if (calcFlag == 2)
      childProduct = *(double *) childProductPtr;
  }

  for (i = DAD; i <= MOM; i++) {
    if (modelOptions->polynomial == TRUE)
      newWeightPolynomial[i] = constant1Poly;   // Non-destructive assignment.
    else
      newWeight[i] = 1.0;
    pConditional = &pParent[i]->pLikelihood[multiLocusIndex[i]];
    if (pParent[i]->touchedFlag == TRUE) {
      /* we have worked on this parent before */
      if (pParent[i] != pProband) {
        if (modelOptions->polynomial == TRUE)
          newWeightPolynomial[i] = timesExp (2, pConditional->lkslot.likelihoodPolynomial, 1, pConditional->wtslot.weightPolynomial, 1, 1 /* Resetting, so freeing */ );
        else
          newWeight[i] = pConditional->lkslot.likelihood * pConditional->wtslot.weight;

      }
    } else {    /* first time we work on this parent */
      if (pNucFam->pParents[i]->pParents[DAD] == NULL) {
        /*
         * this parent is a founder. the weight is
         * just the multiplication of genotype
         * weights at each locus under LE. The weight
         * should have been passed in as an input
         */
        if (modelOptions->equilibrium == LINKAGE_EQUILIBRIUM) {

          if (modelOptions->polynomial == TRUE)
            newWeightPolynomial[i] = (Polynomial *) dWeight[i]; // Never overwritten unless constant1Poly
          else
            newWeight[i] = *((double *) dWeight + i);
        } else if (pParent[i]->loopBreaker == 0) {      /* founder under LD */
          if (modelOptions->polynomial == TRUE) {
            get_haplotype_freq (numLocus - 1, i, &newWeightPolynomial[i]);
          } else {
            get_haplotype_freq (numLocus - 1, i, &newWeight[i]);
          }
        }       /* end of founder and LD */
      } /* founder */
    }   /* end of first time on this parent */
    if (pParent[i]->touchedFlag != TRUE && pParent[i] == pProband) {
      if (modelOptions->polynomial == TRUE) {
        pConditional->wtslot.weightPolynomial = newWeightPolynomial[i];
        newWeightPolynomial[i] = constant1Poly; // Assigning, but stored first.
      } else {
        pConditional->wtslot.weight = newWeight[i];
        newWeight[i] = 1.0;
      }
    }
    /*
     * need to multiply the penetrance if disease locus is in and
     * we haven't calculated any likelihood on this parent before
     */
    traitLocus = analysisLocusList->traitLocusIndex;
    if (traitLocus >= 0 && pParent[i]->touchedFlag != TRUE && (pParent[i]->loopBreaker == 0 || pParent[i]->pParents[DAD] != NULL)) {
      genoIndex = pHaplo->pParentalPairInd[traitLocus];
      if (modelOptions->polynomial == TRUE)
        penetrancePolynomial[i] = pHaplo->ppParentalPair[traitLocus][genoIndex].pGenotype[i]->penslot.penetrancePolynomial;
      else
        penetrance[i] = pHaplo->ppParentalPair[traitLocus][genoIndex].pGenotype[i]->penslot.penetrance;
    } else {
      if (modelOptions->polynomial == TRUE)
        penetrancePolynomial[i] = constant1Poly;        // Non-destructive assignment
      else
        penetrance[i] = 1.0;

    }
  }     /* loop over each parent */

  /* when calcFlag ==2, we just use existing results */
  if (calcFlag != 2) {
    /* now work on the children conditional on this parental pair */
    childProduct = 1;
    multCount = 0;
    if (modelOptions->polynomial == TRUE)
      childProductPolynomial = constant1Poly;   // Non-destructive assignment
    for (child = 0; child < pNucFam->numChildren; child++) {
      pChild = pNucFam->ppChildrenList[child];

      xmissionIndex[DAD] = 0;
      xmissionIndex[MOM] = 0;

      if (modelOptions->polynomial == TRUE) {
        sumPolynomial = constant0Poly;
        loop_child_multi_locus_genotype (0, 0, xmissionIndex);
        childProductPolynomial = timesExp (2, childProductPolynomial, 1, (Polynomial *) sumPolynomial, 1, 1);
      } else {
        sum = 0;
        loop_child_multi_locus_genotype (0, 0, xmissionIndex);
        childProduct *= sum;
      }
    }   /* looping over all children */
  }

  /* results processing */
  ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].likelihoodIndex = multiLocusIndex[proband];
  ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].count = 1;
  if (modelOptions->polynomial == TRUE) {
    ppairMatrix[multiLocusPhase[proband]]
        [multiLocusPhase[spouse]].slot.likelihoodPolynomial =
        timesExp (5, childProductPolynomial, 1, newWeightPolynomial[proband], 1, newWeightPolynomial[spouse], 1, penetrancePolynomial[proband], 1, penetrancePolynomial[spouse], 1, 0 /* End of call, discarding */ );
#if 0
    DIAG (LIKELIHOOD, 1, {
        fprintf (stderr, "\t\t likelihood (%d) = %e\n",}
        ); ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].likelihoodIndex, evaluateValue (ppairMatrix[multiLocusPhase[proband]]
            [multiLocusPhase[spouse]].slot.likelihoodPolynomial));
#endif
    //    discardPoly (newWeightPolynomial[proband]); // This blows-up
    //    discardPoly (newWeightPolynomial[spouse]); // This blows-up
    //    discardPoly (penetrancePolynomial[proband]);
    //    discardPoly (penetrancePolynomial[spouse]);
  } else {
    /* save it */
    ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].slot.likelihood = newWeight[proband] * newWeight[spouse] * penetrance[proband] * penetrance[spouse] * childProduct;
    DIAG (LIKELIHOOD, 1, {
          fprintf (stderr, "Parents: DAD(%s) weight %e   MOM(%s) weight %e \n", pParent[DAD]->sID, newWeight[proband], pParent[MOM]->sID, newWeight[spouse]);
        }
    );
    DIAG (LIKELIHOOD, 1, {
          fprintf (stderr, "Parents: DAD(%s) pen %e   MOM(%s) pen %e \n", pParent[DAD]->sID, penetrance[proband], pParent[MOM]->sID, penetrance[spouse]);
        }
    );
    DIAG (LIKELIHOOD, 1, {
          fprintf (stderr, "\t\t likelihood (%d) = %e\n", ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].likelihoodIndex, ppairMatrix[multiLocusPhase[proband]][multiLocusPhase[spouse]].slot.likelihood);
        }
    );
  }

  return 0;
}


/*
 * Get haplotype frequency for founders under LD analysis */
void get_haplotype_freq (int locus, int myParent, void *freqPtr)
{
  int origLocus1, origLocus2;   /* locus indices in the
                                 * original locus list for
                                 * the two loci in LD */
  double freq[2] = { 0, 0 };    /* variable to store the calculated
                                 * frequency */
  Polynomial *freqPolynomial[2] = { NULL, NULL };
  char vName[100];
  ParentalPair *pPair1; /* parental pair for one locus */
  ParentalPair *pPair2; /* parental pair for the other locus */
  Locus *pLocus1;       /* first locus */
  Locus *pLocus2;       /* the other locus */
  int i, k, l;
  LDLoci *pLDLociLocal; /* structure contains LD parameter values */
  int alleleID1, alleleID2;     /* allele IDs */
  AlleleSet *pAlleleSet1, *pAlleleSet2; /* allele sets */
  int allele1, allele2; /* */


  if (modelOptions->polynomial == TRUE) {
    freqPolynomial[0] = constant0Poly;
    freqPolynomial[1] = constant0Poly;
  }


  /* locus index in the original locus list for the first locus */
  origLocus1 = analysisLocusList->pLocusIndex[locus - 1];
  /* locus index in the original locus list for the second locus */
  origLocus2 = analysisLocusList->pLocusIndex[locus];
  pLocus1 = originalLocusList.ppLocusList[origLocus1];
  pLocus2 = originalLocusList.ppLocusList[origLocus2];
  /* find the parameter values for these two loci */
  pLDLociLocal = find_LD_loci (origLocus1, origLocus2);
  ASSERT (pLDLociLocal != NULL, "Can't find LD parameter between loci %d,%d.\n", origLocus1, origLocus2);
  /*
   * now find the corresponding haplotype frequency : 2 haplotypes
   * paternal haplotype & maternal haplotype
   */
  pPair1 = &pHaplo->ppParentalPair[locus - 1][pHaplo->pParentalPairInd[locus - 1]];
  pPair2 = &pHaplo->ppParentalPair[locus][pHaplo->pParentalPairInd[locus]];
  for (i = DAD; i <= MOM; i++) {
    /* allele ID in the first locus */
    alleleID1 = pPair1->pGenotype[myParent]->allele[i];
    /* allele ID in the second locus */
    alleleID2 = pPair2->pGenotype[myParent]->allele[i];
    pAlleleSet1 = pLocus1->ppAlleleSetList[alleleID1 - 1];
    pAlleleSet2 = pLocus2->ppAlleleSetList[alleleID2 - 1];
    if (modelOptions->polynomial == TRUE) {
      freqPolynomial[i] = constant0Poly;
    } else
      freq[i] = 0;

    for (k = 0; k < pAlleleSet1->numAllele; k++) {
      for (l = 0; l < pAlleleSet2->numAllele; l++) {
        allele1 = pAlleleSet1->pAlleles[k];
        allele2 = pAlleleSet2->pAlleles[l];
        if (modelOptions->polynomial == TRUE) {
          sprintf (vName, "ppHaploFreq_lA%d_rA%d", allele1 - 1, allele2 - 1);
          freqPolynomial[i] = plusExp (2, 1.0, freqPolynomial[i], 1.0, variableExp (&pLDLociLocal->ppHaploFreq[allele1 - 1][allele2 - 1], NULL, 'D', vName), 1);
        } else
          freq[i] += pLDLociLocal->ppHaploFreq[allele1 - 1][allele2 - 1];
      }
    }

    /* for xchr and father, there is only one haplotype */
    if ((modelOptions->sexLinked != 0) && myParent == DAD) {

      if (modelOptions->polynomial == TRUE) {
        if (freqPolynomial[i + 1] != NULL)
          discardPoly (freqPolynomial[i + 1]);
        freqPolynomial[i + 1] = constant1Poly;
      } else {
        freq[i + 1] = 1.0;
      }
      break;
    }

  }     /* end of loop of parents */

  if (modelOptions->polynomial == TRUE) {

    *(Polynomial **) freqPtr = timesExp (2, freqPolynomial[0], 1, freqPolynomial[1], 1, 0);
  } else {
    *(double *) freqPtr = freq[0] * freq[1];
  }

}

/*
 * loop over a child's list of genotypes that are compatible with the
 * parental pair retrieve the transmission probability saved in the
 * transmission matrix sum the likelihood for each genotype configuration
 */
int loop_child_multi_locus_genotype (int locus, int multiLocusIndex, int xmissionIndex[2])
{
  double newProb = 1;
  int i;
  int newMultiLocusIndex;
  int newXmissionIndex[2];
  ParentalPair *pParentalPair;

  /* number of possible genotypes at this locus for this child */
  /* child's conditional likelihood offset for the multilocus genotype this function is building */
  multiLocusIndex *= pChild->pSavedNumGenotype[analysisLocusList->pLocusIndex[locus]];
  /* build the index to xmission matrix for paternal inheritance and maternal inheritance */
  xmissionIndex[DAD] <<= 2;
  xmissionIndex[MOM] <<= 2;

  /* loop through all of this child's compatible genotypes at this locus */
  pParentalPair = &pHaplo->ppParentalPair[locus][pHaplo->pParentalPairInd[locus]];
  for (i = 0; i < pParentalPair->pChildGenoLen[child]; i++) {
    pGenotype = pParentalPair->pppChildGenoList[child][i];
    /* record the index to the genotype list for this child */
    pHaplo->pChildGenoInd[locus] = i;
    DIAG (LIKELIHOOD, 1, {
          fprintf (stderr, "\t child %s locus %4d -> %4d|%-4d \n", pChild->sID, analysisLocusList->pLocusIndex[locus], pGenotype->allele[DAD], pGenotype->allele[MOM]);
        }
    );
    /* record this child's conditional likelihood index */
    newMultiLocusIndex = multiLocusIndex + pGenotype->position;

    /* check the transmission probability */
    for (parent = DAD; parent <= MOM; parent++) {
      newChromosome[parent] = pParentalPair->ppChildInheritance[parent][child][i];
      /* xmissionIndex has already been multiplied by 4 before the loop */
      newXmissionIndex[parent] = xmissionIndex[parent] | newChromosome[parent];
    }   /* looping paternal and maternal chromosomes */
    if (locus < analysisLocusList->numLocus - 1) {
      loop_child_multi_locus_genotype (locus + 1, newMultiLocusIndex, newXmissionIndex);
    } else {

      /* get the transmission probability from the matrix */
      if (modelOptions->polynomial == TRUE) {
        if ((modelOptions->sexLinked != 0) && pChild->sex + 1 == MALE) {
          newProbPolynomial = xmissionMatrix[newXmissionIndex[MOM]].slot.probPoly[2];
        } else {
          newProbPolynomial = timesExp (2, xmissionMatrix[newXmissionIndex[DAD]].slot.probPoly[1], 1, xmissionMatrix[newXmissionIndex[MOM]].slot.probPoly[2], 1, 0);
        }
#if 0
        DIAG (LIKELIHOOD, 1, {
            fprintf (stderr, "\t xmission prob: %f = %f * %f\n", evaluateValue (newProbPolynomial),}
            ); evaluateValue (xmissionMatrix[newXmissionIndex[DAD]].slot.probPoly[1]), evaluateValue (xmissionMatrix[newXmissionIndex[MOM]].slot.probPoly[2]));
#endif
      } else {
        if ((modelOptions->sexLinked != 0) && pChild->sex + 1 == MALE) {
          newProb = xmissionMatrix[newXmissionIndex[MOM]].slot.prob[2];
          DIAG (LIKELIHOOD, 1, {
                fprintf (stderr, "\t xmission prob: %f = %f\n", newProb, xmissionMatrix[newXmissionIndex[MOM]].slot.prob[2]);
              }
          );
        } else {
          newProb = xmissionMatrix[newXmissionIndex[DAD]].slot.prob[1] * xmissionMatrix[newXmissionIndex[MOM]].slot.prob[2];
          DIAG (LIKELIHOOD, 1, {
                fprintf (stderr, "\t xmission prob: %f = %f * %f\n", newProb, xmissionMatrix[newXmissionIndex[DAD]].slot.prob[1], xmissionMatrix[newXmissionIndex[MOM]].slot.prob[2]);
              }
          );
        }
      }

      /* we have completed one multilocus genotype for this child */
      /*
       * check whether we already have some information about this
       * kid we should have if this kid is a connector to another
       * nuclear family we have processed before
       */
      if (calcFlag == 1 && multCount >= maxChildElements) {
        /* resizing likelihoodchildElements array */
        maxChildElements += 1024;
        REALCHOKE (likelihoodChildElements, sizeof (ChildElement) * maxChildElements, ChildElement *);
      }
      if (modelOptions->polynomial == TRUE) {
        if (pChild != pProband) {
          /* the child is not a proband */
          if (pChild->touchedFlag == 1) {
            /*
             * some likelihood calculation has
             * been calculated for this child
             */
            *(Polynomial **) childSum = plusExp (2, 1.0, *(Polynomial **) childSum, 1.0, timesExp (2, newProbPolynomial, 1, pChild->pLikelihood[newMultiLocusIndex].lkslot.likelihoodPolynomial, 1, 1), //end of timesExp
                1);
            if (calcFlag == 1) {
              likelihoodChildElements[multCount].fslot.factorPolynomial = pChild->pLikelihood[newMultiLocusIndex].lkslot.likelihoodPolynomial;
            }
            //end of plusExp
#if 0
            DIAG (LIKELIHOOD, 1, {
                  fprintf (stderr, "\t use already calculated child prob %e \n", evaluateValue (pChild->pLikelihood[newMultiLocusIndex].lkslot.likelihoodPolynomial));
                }
            );
#endif
          } else if (analysisLocusList->traitLocusIndex >= 0)
            /*
             * first time working on this child's
             * current multilocus genotype and we
             * need to consider penetrance
             */
          {
            traitGenoIndex = pHaplo->pChildGenoInd[analysisLocusList->traitLocusIndex];
            pTraitParentalPair = &pHaplo->ppParentalPair[analysisLocusList->traitLocusIndex][pHaplo->pParentalPairInd[analysisLocusList->traitLocusIndex]];
            *(Polynomial **) childSum = plusExp (2, 1.0, *(Polynomial **) childSum, 1.0, timesExp (2, newProbPolynomial, 1, pTraitParentalPair->pppChildGenoList[child]
                    [traitGenoIndex]->penslot.penetrancePolynomial, 1, 1),      //end of timesExp
                1);

            if (calcFlag == 1) {
              likelihoodChildElements[multCount].fslot.factorPolynomial = pTraitParentalPair->pppChildGenoList[child]
                  [traitGenoIndex]->penslot.penetrancePolynomial;
            }
          } else {
            /*
             * no trait locus and new to this
             * child
             */
            *(Polynomial **) childSum = plusExp (2, 1.0, *(Polynomial **) childSum, 1.0, newProbPolynomial, 1);
            if (calcFlag == 1) {
              likelihoodChildElements[multCount].fslot.factorPolynomial = constant1Poly;
            }
          }
        } else {        /* this child is proband */
          *(Polynomial **) childSum = plusExp (2, 1.0, *(Polynomial **) childSum, 1.0, newProbPolynomial, 1);
          if (calcFlag == 1) {
            likelihoodChildElements[multCount].fslot.factorPolynomial = constant1Poly;
          }
        }
#if 0
        DIAG (LIKELIHOOD, 1, {
              fprintf (stderr, "\t child sum %e \n", evaluateValue (*(Polynomial **) childSum));
            }
        );
#endif
      } else {  /* PE is not turned on */
        if (pChild != pProband) {
          /* the child is not a proband */
          if (pChild->touchedFlag == 1) {
            /*
             * some likelihood calculation has
             * been done for this child
             */
            *(double *) childSum += newProb * pChild->pLikelihood[newMultiLocusIndex].lkslot.likelihood;
            DIAG (LIKELIHOOD, 1, {
                  fprintf (stderr, "\t use already calculated child prob %e \n", pChild->pLikelihood[newMultiLocusIndex].lkslot.likelihood);
                });
            if (calcFlag == 1) {
              likelihoodChildElements[multCount].fslot.factor = pChild->pLikelihood[newMultiLocusIndex].lkslot.likelihood;
            }
          } else if (analysisLocusList->traitLocusIndex >= 0)
            /*
             * first time working on this child's
             * current multilocus genotype and we
             * need to consider penetrance
             */
          {
            traitGenoIndex = pHaplo->pChildGenoInd[analysisLocusList->traitLocusIndex];
            pTraitParentalPair = &pHaplo->ppParentalPair[analysisLocusList->traitLocusIndex][pHaplo->pParentalPairInd[analysisLocusList->traitLocusIndex]];
            *(double *) childSum += newProb * pTraitParentalPair->pppChildGenoList[child][traitGenoIndex]->penslot.penetrance;
            if (calcFlag == 1) {
              likelihoodChildElements[multCount].fslot.factor = pTraitParentalPair->pppChildGenoList[child][traitGenoIndex]->penslot.penetrance;
            }
            DIAG (LIKELIHOOD, 1, {
                  fprintf (stderr, "child penetrance %e\n", pTraitParentalPair->pppChildGenoList[child][traitGenoIndex]->penslot.penetrance);
                }
            );

          } else {
            *(double *) childSum += newProb;
            if (calcFlag == 1) {
              likelihoodChildElements[multCount].fslot.factor = 1;
            }
          }
        } else {        /* this child is proband */
          /*
           * penetrance if applicable will be figured
           * into later
           */
          *(double *) childSum += newProb;
          if (calcFlag == 1) {
            likelihoodChildElements[multCount].fslot.factor = 1;
          }
        }
        DIAG (LIKELIHOOD, 1, {
              fprintf (stderr, "\t child sum %e \n", *(double *) childSum);
            });
      }


      if (calcFlag == 1) {
        likelihoodChildElements[multCount].xmissionIndex[DAD] = newXmissionIndex[DAD];
        likelihoodChildElements[multCount].xmissionIndex[MOM] = newXmissionIndex[MOM];
        likelihoodChildCount[child]++;
        multCount++;
      }
    }   /* end of processing one complete multilocus genotype */

  }     /* loop possible genotypes at this locus */
  return 0;
}


/*
 * A recursive call to build transmission probability matrix pMatrix - pass
 * in matrix pointer - This should have been pre-allocated totalLoci - total
 * number of loci loc - current locus prob -
 */

int do_populate_xmission_matrix (XMission * pMatrix, int totalLoci, void *prob[3], void *prob2[3], void *hetProb[3], int cellIndex, int lastHetLoc, int prevPattern, int loc)
{
  int pattern;

  Polynomial *newProbPoly[3];
  Polynomial *newProbPoly2[3];
  Polynomial *newHetProbPoly[3];
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

      if (modelOptions->polynomial == TRUE) {
        newProbPoly[i] = (Polynomial *) prob[i];
        newProbPoly2[i] = (Polynomial *) prob2[i];
        newHetProbPoly[i] = (Polynomial *) hetProb[i];
      } else {
        newProb[i] = *((double *) prob[i]);
        newProb2[i] = *((double *) prob2[i]);
        newHetProbPtr[i] = hetProb[i];
      }
    }
    newCellIndex = cellIndex * 4 + pattern;
    newLastHetLoc = lastHetLoc;
    if (pattern != 3 && pattern != 0) {
      /* parent is not homozygous */
      if (lastHetLoc != -1) {
        if (prevPattern != 3 && prevPattern != 0) {
          /* previous locus for the parent is het and 
           * current locus pattern is either paternal or maternal */
          if (prevPattern == pattern) {
            /* no recombination */
            for (i = 0; i < 3; i++) {
              if (modelOptions->polynomial == TRUE) {
                if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                  newProbPoly[i] = newProbPoly[0];
                } else {
                  if (analysisLocusList->traitLocusIndex < 0 || /* no trait locus in the list */
                      /* trait locus is not current locus or previous locus */
                      (analysisLocusList->traitLocusIndex != loc && analysisLocusList->traitLocusIndex != loc - 1)) {
                    /* theta is constant between marker loci 
                     * prob * (1-th)
                     */
                    newProbPoly[i] = timesExp (2, newProbPoly[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, constantExp (analysisLocusList->pPrevLocusDistance[i][loc]), 1), 1, 1);

                  } else {
                    /* prob * (1-th) */
                    sprintf (vName1, "theta_i%d_l%d", i, loc);
                    newProbPoly[i] = timesExp (2, newProbPoly[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, variableExp (&analysisLocusList->pPrevLocusDistance[i][loc], NULL, 'D', vName1), 1), 1, 1);
                  }
                }
              } else {
                newProb[i] *= (1 - analysisLocusList->pPrevLocusDistance[i][loc]);
              }
            }
          } else {
            /* recombination */
            for (i = 0; i < 3; i++)
              if (modelOptions->polynomial == TRUE) {
                if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                  newProbPoly[i] = newProbPoly[0];
                } else {
                  if (analysisLocusList->traitLocusIndex < 0 || /* no trait locus in the list */
                      /* trait locus is not current locus or previous locus */
                      (analysisLocusList->traitLocusIndex != loc && analysisLocusList->traitLocusIndex != loc - 1)) {
                    /* theta is constant between marker loci */
                    newProbPoly[i] = timesExp (2, newProbPoly[i], 1, constantExp (analysisLocusList->pPrevLocusDistance[i][loc]), 1, 0);
                  } else {
                    sprintf (vName1, "theta_i%d_l%d", i, loc);
                    newProbPoly[i] = timesExp (2, newProbPoly[i], 1, variableExp (&analysisLocusList->pPrevLocusDistance[i][loc], NULL, 'D', vName1), 1, 0);
                  }
                }
              } else {
                if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                  newProb[i] = newProb[0];
                } else {
                  newProb[i] *= analysisLocusList->pPrevLocusDistance[i][loc];
                }
              }
          }
        } else {
          /* previous locus at parent is homo and current locus is het */
          for (i = 0; i < 3; i++) {
            if (pattern == 1) {
              /* paternal inheritance for this locus
               * either no recombination from previous paternal strand 
               * or recombination from previous maternal strand */
              if (modelOptions->polynomial == TRUE) {
                if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                  newProbPoly[i] = newProbPoly[0];
                } else {
                  if (analysisLocusList->traitLocusIndex < 0 || /* no trait locus in the list */
                      /* trait locus is not current locus or previous locus */
                      (analysisLocusList->traitLocusIndex != loc && analysisLocusList->traitLocusIndex != loc - 1)) {
                    newProbPoly[i] =
                        plusExp (2,
                        1.0, timesExp (2,
                            (Polynomial *) prob[i], 1,
                            plusExp (2, 1.0, constantExp (1.0), -1.0, constantExp (analysisLocusList->pPrevLocusDistance[i][loc]), 1), 1, 1), 1.0, timesExp (2, (Polynomial *) prob2[i], 1, constantExp (analysisLocusList->pPrevLocusDistance[i][loc]), 1, 1), 0);

                  } else {
                    sprintf (vName1, "theta_i%d_l%d", i, loc);
                    newProbPoly[i] = plusExp (2, 1.0, timesExp (2, (Polynomial *)
                            prob[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, variableExp (&analysisLocusList->pPrevLocusDistance[i]
                                    [loc], NULL, 'D', vName1), 1), 1, 1), 1.0, timesExp (2, (Polynomial *)
                            prob2[i], 1, variableExp (&analysisLocusList->pPrevLocusDistance[i][loc], NULL, 'D', vName1), 1, 1), 0);
                  }
                }
              } else {
                if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                  newProb[i] = newProb[0];
                } else {
                  newProb[i] = *((double *) prob[i]) * (1 - analysisLocusList->pPrevLocusDistance[i][loc]) + *((double *) prob2[i]) * analysisLocusList->pPrevLocusDistance[i][loc];
                }
              }
            } else {
              /* has to be maternal */
              if (modelOptions->polynomial == TRUE) {
                if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                  newProbPoly[i] = newProbPoly[0];
                } else {
                  if (analysisLocusList->traitLocusIndex < 0 || /* no trait locus in the list */
                      /* trait locus is not current locus or previous locus */
                      (analysisLocusList->traitLocusIndex != loc && analysisLocusList->traitLocusIndex != loc - 1)) {
                    newProbPoly[i] =
                        plusExp (2,
                        1.0, timesExp (2,
                            (Polynomial *) prob2[i], 1,
                            plusExp (2, 1.0, constantExp (1.0), -1.0, constantExp (analysisLocusList->pPrevLocusDistance[i][loc]), 1), 1, 1), 1.0, timesExp (2, (Polynomial *) prob[i], 1, constantExp (analysisLocusList->pPrevLocusDistance[i][loc]), 1, 1), 0);
                  } else {
                    sprintf (vName1, "theta_i%d_l%d", i, loc);
                    newProbPoly[i] = plusExp (2, 1.0, timesExp (2, (Polynomial *)
                            prob2[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, variableExp (&analysisLocusList->pPrevLocusDistance[i]
                                    [loc], NULL, 'D', vName1), 1), 1, 1), 1.0, timesExp (2, (Polynomial *)
                            prob[i], 1, variableExp (&analysisLocusList->pPrevLocusDistance[i][loc], NULL, 'D', vName1), 1, 1), 0);
                  }
                }

              } else if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                newProb[i] = newProb[0];
              } else {
                newProb[i] = *((double *) prob2[i]) * (1 - analysisLocusList->pPrevLocusDistance[i][loc]) + *((double *) prob[i]) * analysisLocusList->pPrevLocusDistance[i][loc];
              }

            }
          }

        }       /* end of prevPattern is homo and current pattern is het */
      } /* end of prevHetLoc != -1 */
      else {
        /* we don't have any het locus yet, this locus is the first het */
        if (modelOptions->polynomial == TRUE) {
          for (i = 0; i < 3; i++)
            newProbPoly[i] = constantExp (0.5);
        } else {
          for (i = 0; i < 3; i++)
            newProb[i] = 0.5;
        }
      }
      newLastHetLoc = loc;


      if (modelOptions->polynomial == TRUE) {
        for (i = 0; i < 3; i++)
          newHetProbPoly[i] = newProbPoly[i];
      } else {
        for (i = 0; i < 3; i++)
          newHetProbPtr[i] = newProbPtr[i];
      }


    } /* end of current pattern is not homo */
    else {      /* current pattern is homo */
      if (lastHetLoc == -1)
        /* nothing needs to be done for this locus */
        ;
      else {
        if (loc == totalLoci - 1) {
          /* this is the last locus and it's homo, take the previous het locus */
          for (i = 0; i < 3; i++) {
            if (modelOptions->polynomial == TRUE) {
              newProbPoly[i] = (Polynomial *) hetProb[i];
            } else
              newProb[i] = *(double *) hetProb[i];
          }
        } else {
          if (prevPattern == 3 || prevPattern == 0) {   /* previous locus pattern is homo */
            for (i = 0; i < 3; i++) {

              if (modelOptions->polynomial == TRUE) {
                if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                  newProbPoly[i] = newProbPoly[0];
                  newProbPoly2[i] = newProbPoly2[0];
                } else {
                  if (analysisLocusList->traitLocusIndex < 0 || /* no trait locus in the list */
                      /* trait locus is not current locus or previous locus */
                      (analysisLocusList->traitLocusIndex != loc && analysisLocusList->traitLocusIndex != loc - 1)) {
                    /* theta is constant between marker loci */
                    newProbPoly[i] = plusExp (2, 1.0, timesExp (2, (Polynomial *)
                            prob[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, constantExp (analysisLocusList->pPrevLocusDistance[i]
                                    [loc]), 1), 1, 1), 1.0, timesExp (2, (Polynomial *)
                            prob2[i], 1, constantExp (analysisLocusList->pPrevLocusDistance[i][loc]), 1, 1), 0);
                    /* prevProb2 * (1-th1) + prevProb * th1 */
                    newProbPoly2[i] = plusExp (2, 1.0, timesExp (2, (Polynomial *)
                            prob2[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, constantExp (analysisLocusList->pPrevLocusDistance[i]
                                    [loc]), 1), 1, 1), 1.0, timesExp (2, (Polynomial *)
                            prob[i], 1, constantExp (analysisLocusList->pPrevLocusDistance[i][loc]), 1, 1), 0);

                  } else {      /* dealing with trait locus */

                    sprintf (vName1, "theta_i%d_l%d", i, loc);
                    /* prevProb * (1-th1) + prevProb2 * th1 */
                    newProbPoly[i] = plusExp (2, 1.0, timesExp (2, (Polynomial *)
                            prob[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, variableExp (&analysisLocusList->pPrevLocusDistance[i]
                                    [loc], NULL, 'D', vName1), 1), 1, 1), 1.0, timesExp (2, (Polynomial *)
                            prob2[i], 1, variableExp (&analysisLocusList->pPrevLocusDistance[i][loc], NULL, 'D', vName1), 1, 1), 0);
                    /* prevProb2 * (1-th1) + prevProb * th1 */
                    newProbPoly2[i] = plusExp (2, 1.0, timesExp (2, (Polynomial *)
                            prob2[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, variableExp (&analysisLocusList->pPrevLocusDistance[i]
                                    [loc], NULL, 'D', vName1), 1), 1, 1), 1.0, timesExp (2, (Polynomial *)
                            prob[i], 1, variableExp (&analysisLocusList->pPrevLocusDistance[i][loc], NULL, 'D', vName1), 1, 1), 0);
                  }
                }

              } else {
                if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                  newProb[i] = newProb[0];
                  newProb2[i] = newProb2[0];
                } else {
                  newProb[i] = *(double *) prob[i] * (1 - analysisLocusList->pPrevLocusDistance[i][loc]) + *((double *) prob2[i]) * analysisLocusList->pPrevLocusDistance[i][loc];

                  newProb2[i] = *(double *) prob2[i] * (1 - analysisLocusList->pPrevLocusDistance[i][loc]) + *((double *) prob[i]) * analysisLocusList->pPrevLocusDistance[i][loc];
                }
              }
            }
          } else {      /* prev pattern is het */
            for (i = 0; i < 3; i++) {
              if (prevPattern == 1) {
                if (modelOptions->polynomial == TRUE) {
                  if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                    newProbPoly[i] = newProbPoly[0];
                    newProbPoly2[i] = newProbPoly2[0];
                  } else {
                    if (analysisLocusList->traitLocusIndex < 0 ||       /* no trait locus in the list */
                        /* trait locus is not current locus or previous locus */
                        (analysisLocusList->traitLocusIndex != loc && analysisLocusList->traitLocusIndex != loc - 1)) {
                      /* theta is constant between marker loci */
                      newProbPoly[i] = timesExp (2, (Polynomial *)
                          prob[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, constantExp (analysisLocusList->pPrevLocusDistance[i][loc]), 1), 1, 0);
                      newProbPoly2[i] = timesExp (2, (Polynomial *)
                          prob[i], 1, constantExp (analysisLocusList->pPrevLocusDistance[i]
                              [loc]), 1, 0);

                    } else {    /* dealing with trait locus */

                      sprintf (vName1, "theta_i%d_l%d", i, loc);
                      newProbPoly[i] = timesExp (2, (Polynomial *)
                          prob[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, variableExp (&analysisLocusList->pPrevLocusDistance[i][loc], NULL, 'D', vName1), 1), 1, 0);
                      newProbPoly2[i] = timesExp (2, (Polynomial *)
                          prob[i], 1, variableExp (&analysisLocusList->pPrevLocusDistance[i]
                              [loc], NULL, 'D', vName1), 1, 0);
                    }
                  }
                } else {
                  if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                    newProb[i] = newProb[0];
                    newProb2[i] = newProb2[0];
                  } else {
                    newProb[i] = *(double *) prob[i] * (1 - analysisLocusList->pPrevLocusDistance[i][loc]);
                    newProb2[i] = *(double *) prob[i] * analysisLocusList->pPrevLocusDistance[i][loc];
                  }
                }
              } else {
                if (modelOptions->polynomial == TRUE) {
                  if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                    newProbPoly[i] = newProbPoly[0];
                    newProbPoly2[i] = newProbPoly2[0];
                  } else {
                    if (analysisLocusList->traitLocusIndex < 0 ||       /* no trait locus in the list */
                        /* trait locus is not current locus or previous locus */
                        (analysisLocusList->traitLocusIndex != loc && analysisLocusList->traitLocusIndex != loc - 1)) {
                      /* theta is constant between marker loci */
                      newProbPoly2[i] = timesExp (2, (Polynomial *)
                          prob[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, constantExp (analysisLocusList->pPrevLocusDistance[i][loc]), 1), 1, 0);
                      newProbPoly[i] = timesExp (2, (Polynomial *)
                          prob[i], 1, constantExp (analysisLocusList->pPrevLocusDistance[i]
                              [loc]), 1, 0);

                    } else {    /* dealing with trait locus */

                      sprintf (vName1, "theta_i%d_l%d", i, loc);
                      newProbPoly2[i] = timesExp (2, (Polynomial *)
                          prob[i], 1, plusExp (2, 1.0, constantExp (1.0), -1.0, variableExp (&analysisLocusList->pPrevLocusDistance[i][loc], NULL, 'D', vName1), 1), 1, 0);
                      newProbPoly[i] = timesExp (2, (Polynomial *)
                          prob[i], 1, variableExp (&analysisLocusList->pPrevLocusDistance[i]
                              [loc], NULL, 'D', vName1), 1, 0);
                    }
                  }
                } else {
                  if (i > 0 && modelOptions->mapFlag == SEX_AVERAGED) {
                    newProb[i] = newProb[0];
                    newProb2[i] = newProb2[0];
                  } else {
                    newProb2[i] = *(double *) prob[i] * (1 - analysisLocusList->pPrevLocusDistance[i][loc]);
                    newProb[i] = *(double *) prob[i] * analysisLocusList->pPrevLocusDistance[i][loc];
                  }
                }
              }

            }
          }
        }
      }
    }

    if (loc == totalLoci - 1) {
      /* we have a complete set of multilocus inheritance pattern */

      for (i = 0; i < 3; i++) {
        if (modelOptions->polynomial == TRUE) {
          if (pMatrix[newCellIndex].slot.probPoly[i] != NULL)
            unHoldPoly (pMatrix[newCellIndex].slot.probPoly[i]);
          pMatrix[newCellIndex].slot.probPoly[i] = newProbPoly[i];
          holdPoly (newProbPoly[i]);
        } else
          pMatrix[newCellIndex].slot.prob[i] = newProb[i];
      }

    } else {
      /* move on to next locus */
      if (modelOptions->polynomial == TRUE) {
        do_populate_xmission_matrix (pMatrix, totalLoci, (void *) newProbPoly, (void *) newProbPoly2, (void *) newHetProbPoly, newCellIndex, newLastHetLoc, pattern, loc + 1);
      } else
        do_populate_xmission_matrix (pMatrix, totalLoci, (void *) newProbPtr, (void *) newProbPtr2, (void *) newHetProbPtr, newCellIndex, newLastHetLoc, pattern, loc + 1);
    }
  }
  return 0;
}

int populate_xmission_matrix (XMission * pMatrix, int totalLoci, void *prob[3], void *prob2[3], void *hetProb[3], int cellIndex, int lastHetLoc, int prevPattern, int loc)
{

  do_populate_xmission_matrix (pMatrix, totalLoci, prob, prob2, hetProb, cellIndex, lastHetLoc, prevPattern, loc);

  return 0;
}

int build_xmission_matrix (XMission ** ppMatrix, int totalLoci)
{
  int size;

  //int           i;

  *ppMatrix = NULL;
  size = pow (4, totalLoci);
  /* minimal two loci */
  //if (size < 9)
  //size = 9;
  CALCHOKE (*ppMatrix, (size_t) size, sizeof (XMission), XMission *);
  if (*ppMatrix == NULL)
    /* memory allocation failed */
    return -1;

  return 0;
}

void print_xmission_matrix (XMission * pMatrix, int totalLoci, int loc, int cellIndex, char *pID)
{
  int pattern;
  int newCellIndex;
  int i;

  for (pattern = 0; pattern <= 2; pattern++) {
    newCellIndex = cellIndex * 4 + pattern;
    if (pattern == 1)
      pID[loc] = 'P';
    else if (pattern == 2)
      pID[loc] = 'M';
    else
      pID[loc] = 'B';

    if (loc != totalLoci - 1) {
      /* Not complete multi-locus haplotype yet */
      print_xmission_matrix (pMatrix, totalLoci, loc + 1, newCellIndex, pID);
    } else {
      /* Print the xmission probability */
      for (i = 0; i <= loc; i++)
        fprintf (stderr, "%c", pID[i]);
      fprintf (stderr, ": ");
      /* Print sex averaged xmission probability */
      if (modelOptions->polynomial == TRUE) {
        expTermPrinting (stderr, pMatrix[newCellIndex].slot.probPoly[0], 4);
        fprintf (stderr, "\n");
      } else
        fprintf (stderr, "%f\n", pMatrix[newCellIndex].slot.prob[0]);
    }
  }
}

/*
 * allocate storage for keeping track of het locus in nuclear families
 * numLocus - number of loci analyzing at a time
 */
void allocate_nucfam_het (PedigreeSet * pPedigreeList, int myNumLocus)
{
  int ped;
  Pedigree *pPedigree;
  int fam;
  NuclearFamily *pMyNucFam;

  for (ped = 0; ped < pPedigreeList->numPedigree; ped++) {
    pPedigree = pPedigreeList->ppPedigreeSet[ped];
    for (fam = 0; fam < pPedigree->numNuclearFamily; fam++) {
      pMyNucFam = pPedigree->ppNuclearFamilyList[fam];
      if (pMyNucFam->hetFlag[DAD] == NULL) {
        CALCHOKE (pMyNucFam->hetFlag[DAD], sizeof (int), (size_t) myNumLocus, int *);
        CALCHOKE (pMyNucFam->hetFlag[MOM], sizeof (int), (size_t) myNumLocus, int *);
        CALCHOKE (pMyNucFam->tmpNumHet[DAD], sizeof (int), (size_t) myNumLocus, int *);
        CALCHOKE (pMyNucFam->tmpNumHet[MOM], sizeof (int), (size_t) myNumLocus, int *);
        CALCHOKE (pMyNucFam->relatedPPairStart, sizeof (int), (size_t) myNumLocus, int *);
        CALCHOKE (pMyNucFam->numRelatedPPair, sizeof (int), (size_t) myNumLocus, int *);
        CALCHOKE (pMyNucFam->totalRelatedPPair, sizeof (int), (size_t) myNumLocus, int *);
      }
    }
  }

}

inline void clear_ppairMatrix (PPairElement ** ppMatrix)
{
  int i;

  for (i = 0; i <= bitMask[ppairMatrixNumLocus]; i++) {
    memset (ppMatrix[i], 0, ppairMatrixRowSize);
  }
}

inline void initialize_proband_tmpLikelihood (Person * pPerson)
{
  int i;

  //  int             size = pPerson->numConditionals;
  ConditionalLikelihood *pConditional;

  for (i = 0; i < pPerson->numTmpLikelihood; i++) {
    pConditional = &pPerson->pLikelihood[pPerson->pTmpLikelihoodIndex[i]];
    pConditional->tmpTouched = FALSE;
    if (modelOptions->polynomial == TRUE)
      pConditional->tmpslot.tmpLikelihoodPolynomial = constant0Poly;
    else
      pConditional->tmpslot.tmpLikelihood = 0;
  }
  pPerson->numTmpLikelihood = 0;
}

void populate_pedigree_loopbreaker_genotype_vector (Pedigree * pPed)
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

void populate_loopbreaker_genotype_vector (Person * pLoopBreaker, int locus)
{
  Genotype *pMyGenotype;
  int i;
  LoopBreaker *pLoopBreakerStruct;

  pMyGenotype = pLoopBreaker->ppSavedGenotypeList[analysisLocusList->pLocusIndex[locus]];
  while (pMyGenotype != NULL) {
    pTempGenoVector[locus] = pMyGenotype;
    if (locus < analysisLocusList->numLocus - 1) {
      populate_loopbreaker_genotype_vector (pLoopBreaker, locus + 1);
    } else {
      /* one complete multilocus genotype */
      pLoopBreakerStruct = pLoopBreaker->loopBreakerStruct;
      i = pLoopBreakerStruct->numGenotype;
      memcpy (pLoopBreakerStruct->genotype[i], pTempGenoVector, sizeof (Genotype *) * analysisLocusList->numLocus);
      pLoopBreakerStruct->numGenotype++;

    }
    pMyGenotype = pMyGenotype->pSavedNext;
  }
}

/* this function is not needed - so not finished */
void sync_loopbreaker_duplicates (Pedigree * pPed)
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
int set_next_loopbreaker_genotype_vector (Pedigree * pPed, int initialFlag)
{
  int numLoopBreaker = pPed->numLoopBreaker;
  int i, j;
  Person *pLoopBreaker;
  int found;
  LoopBreaker *loopStruct;
  int origLocus;
  int locus;
  int ret;

  /* find the next genotype vector for at least one of the loop breaker */
  DIAG (LIKELIHOOD, 1, {
        fprintf (stderr, "Set next loop breaker genotype\n");
      });
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
    j = loopStruct->genotypeIndex;
    DIAG (LIKELIHOOD, 1, {
          fprintf (stderr, "Fix pedigree %s loop breaker %s to the genotype below (%d/%d):\n", pPed->sPedigreeID, pLoopBreaker->sID, loopStruct->genotypeIndex + 1, loopStruct->numGenotype);
        }
    );
    for (locus = 0; locus < analysisLocusList->numLocus; locus++) {
      origLocus = analysisLocusList->pLocusIndex[locus];
      pLoopBreaker->ppGenotypeList[origLocus] = loopStruct->genotype[j][locus];
      loopStruct->genotype[j][locus]->pNext = NULL;
      pLoopBreaker->pNumGenotype[origLocus] = 1;
      DIAG (LIKELIHOOD, 1, {
            fprintf (stderr, "\t %d-> %d|%d \n", locus, loopStruct->genotype[j][locus]->allele[DAD], loopStruct->genotype[j][locus]->allele[MOM]);
          }
      );
    }
  }

  /*
   * as the loop breakers are set with fixed genotype, redo genotype
   * elimination (without actually remove any genotype - only the links
   * will get updated
   */
  for (i = 0; i < analysisLocusList->numLocus; i++) {
    origLocus = analysisLocusList->pLocusIndex[i];
    ret = pedigree_genotype_elimination (origLocus, pPed);
    if (ret < 0)
      return -2;
  }

  return 0;
}
