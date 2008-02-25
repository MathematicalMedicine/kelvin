/**********************************************************************
 * Copyright 2007, Columbus Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

/* This file contains functions to  compute likelihood for all the pedigrees
 * peeling procedure etc. 
 * */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "pedlib.h"
#include "locus.h"
#include "utils.h"		/* for logging */
#include "tools.h"
#include "likelihood.h"
#include "genotype_elimination.h"
#include "polynomial.h"

char *likelihoodVersion = "0.0.29";

//double *half_pow = NULL;
XMission *xmissionMatrix = NULL;

int peel_graph (NuclearFamily * pNucFam, Person * pProband,
		int peelingDirection);
int compute_nuclear_family_likelihood (NuclearFamily * pNucFam,
				       Person * pProband,
				       int peelingDirection);
int
loop_parental_pair (NuclearFamily * pNucFam,
		    Person * pPerson,
		    int locus,
		    ParentalPairSpace * pHaplo,
		    int multiLocusIndex[2],
		    void *dWeight[2], void *dPenetrance[2], void *pParentSum);
int
loop_child_multi_locus_genotype (Person * pPChild, Person * pProband,
				 ParentalPairSpace * pHaplo, int child,
				 int locus, int multiLocusIndex,
				 void *pPenetrance,
				 void *pSum, void *pProb,
				 int prevInformativeLocus[2],
				 int chromosome[2]);
int loop_proband_genotype (NuclearFamily * pNucFam, Person * pProband,
			   int peelingDirection, int locus,
			   int multiLocusIndex,
			   void *penetrance, void *weight);

int
allocate_likelihood_space (PedigreeSet * pPedigreeList, int numLocus)
{
  Pedigree *pPedigree;
  int i;

  for (i = 0; i < pPedigreeList->numPedigree; i++)
    {
      pPedigree = pPedigreeList->ppPedigreeSet[i];
      allocate_multi_locus_genotype_storage (pPedigree, numLocus);
    }

  return 0;
}

int
count_likelihood_space (PedigreeSet * pPedigreeList)
{
  Pedigree *pPedigree;
  int i;

  for (i = 0; i < pPedigreeList->numPedigree; i++)
    {
      pPedigree = pPedigreeList->ppPedigreeSet[i];
      count_multi_locus_genotype_storage (pPedigree);
    }

  return 0;
}

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
}

int
compute_likelihood (PedigreeSet * pPedigreeList)
{
  Pedigree *pPedigree;
  int i;
  int status;
  double product_likelihood = 1;
  double sum_log_likelihood = 0;
  double log10Likelihood;
  int origLocus = locusList->pLocusIndex[1];
  clock_t time2;

  sum_log_likelihood = 0;
  product_likelihood = 1;
  pPedigreeList->likelihood = 1;
  pPedigreeList->log10Likelihood = 0;

  for (i = 0; i < pPedigreeList->numPedigree; i++)
    {
      pPedigree = pPedigreeList->ppPedigreeSet[i];

#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	{
//         fprintf(stderr,"modelOptions.polynomial == TRUE\n");
	  if (pPedigree->likelihoodPolynomial == NULL)
	    {
//                fprintf(stderr,"The polynomial building for this pedigree should be only once\n");
//            free_multi_locus_genotype_storage (pPedigree);
//            allocate_multi_locus_genotype_storage (pPedigree);
	      initialize_multi_locus_genotype (pPedigree);
//                fprintf(stderr,"Start polynomial building\n");
	      makePolynomialStamp2 ();
	      status = compute_pedigree_likelihood (pPedigree);

//                expPrinting(pPedigree->likelihoodPolynomial);
//                fprintf(stderr,"\n");
	      polyStatistics(stderr);
	      pPedigree->likelihoodPolyList = buildPolyList ();
	      polyListSorting (pPedigree->likelihoodPolynomial,
			       pPedigree->likelihoodPolyList);
	      partialPolynomialClearance2 ();
	      if(i == pPedigreeList->numPedigree -1)
		{
		  time2 = clock();
		  fprintf(stderr, "Finished polynomial building: %f\n", 
			  (double)time2/CLOCKS_PER_SEC);
		}
	    }
	  pPedigree->likelihood =
	    evaluatePoly (pPedigree->likelihoodPolynomial,
			  pPedigree->likelihoodPolyList);

	}
      else
	{
	  //      free_multi_locus_genotype_storage (pPedigree);
	  //      allocate_multi_locus_genotype_storage (pPedigree);
	  initialize_multi_locus_genotype (pPedigree);
	  status = compute_pedigree_likelihood (pPedigree);
	  //      free_multi_locus_genotype_storage(pPedigree);
	}
#else
      //      free_multi_locus_genotype_storage (pPedigree);
      //      allocate_multi_locus_genotype_storage (pPedigree);
      initialize_multi_locus_genotype (pPedigree);
      status = compute_pedigree_likelihood (pPedigree);
      //      free_multi_locus_genotype_storage(pPedigree);
#endif

      /* log likelihood is additive 
         sum_log_likelihood += log10(pPedigree->likelihood);
       */
      if (pPedigree->likelihood == 0.0)
	{
	  KLOG (LOGLIKELIHOOD, LOGWARNING,
		"Pedigree %s with likelihood 0 or too small.\n",
		pPedigree->sPedigreeID);
	  fprintf (stderr,
		   "Pedigree %s has likelihood that's too small or 0.\n",
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


int
compute_pedigree_likelihood (Pedigree * pPedigree)
{
  int i;
  NuclearFamily *pNucFam;
  int status;
  Person *pProband;
  double likelihood;
#ifndef NO_POLYNOMIAL
  Polynomial *pLikelihoodPolynomial = NULL;
#endif

  fprintf(stderr, "PEDIGREE: %s (%d/%d)\n", 
       pPedigree->sPedigreeID, pPedigree->pedigreeIndex+1,
       pPedigree->pPedigreeSet->numPedigree);

  /* first, we need to set the genotype weight
   * this weight is the probability observing the genotype 
   * note: this will not be used in LD cases, as haplotype frequencies 
   * should be used in that case */
  //if (modelOptions.analysisType == LINKAGE_EQUILIBRIUM)
  //  set_genotype_weight (pPedigree);

  /* initialize all the nuclear families before peeling starts
   * in many cases, multiple likelihoods are computed for the same family
   * with different parameters, we need to clean up before (or after) 
   * each calculation */
  for (i = 0; i < pPedigree->numNuclearFamily; i++)
    {
      pNucFam = pPedigree->ppNuclearFamilyList[i];
      pNucFam->doneFlag = 0;
    }

  /* peeling starts from the peeling proband and eventually will end 
   * at proband */
  status = peel_graph (pPedigree->pPeelingNuclearFamily,
		       pPedigree->pPeelingProband,
		       pPedigree->peelingDirection);

  /* done peeling, need to add up the conditional likelihood for the proband */
  pProband = pPedigree->pPeelingProband;
  likelihood = 0;

#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    pLikelihoodPolynomial = constantExp (0);
#endif

  for (i = 0; i < pProband->numConditionals; i++)
    {
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	{
	  /* build likelihood polynomial */

	  pLikelihoodPolynomial =
	    plusExp (2, 1.0, pLikelihoodPolynomial,
		     1.0, timesExp (2,
				    pProband->pLikelihood[i].
				    likelihoodPolynomial, 1,
				    pProband->pLikelihood[i].weightPolynomial,
				    1, 0), 1);
	}
      else
	{
	  likelihood += pProband->pLikelihood[i].likelihood *
	    pProband->pLikelihood[i].weight;
	}
#else
      likelihood += pProband->pLikelihood[i].likelihood *
	pProband->pLikelihood[i].weight;
#endif
    }
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    /* save the polynomial to the pedigree structure */
    pPedigree->likelihoodPolynomial = pLikelihoodPolynomial;
  else
    pPedigree->likelihood = likelihood;
#else
  pPedigree->likelihood = likelihood;
  KLOG (LOGLIKELIHOOD, LOGDEBUG, "log Likelihood for pedigree %d is: %e\n",
	pPedigree->pedigreeIndex + 1, log10 (likelihood));
#endif


  return 0;
}

/* recursive procedure to go through all nuclear families for one pedigree */
int
peel_graph (NuclearFamily * pNucFam, Person * pProband, int peelingDirection)
{
  NuclearFamilyConnector *pConnector;
  NuclearFamily *pNucFam2;
  Person *pConnectPerson;
  double weight = 1;
  double penetrance = 1;
#ifndef NO_POLYNOMIAL
  Polynomial *weightPolynomial;
  Polynomial *penetrancePolynomial;
#endif

  if (pNucFam->doneFlag == TRUE)
    return 0;

  /* mark this nuclear family as done */
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
      peel_graph (pConnector->pConnectedNuclearFamily, pConnector->pConnectedPerson, PEDIGREE_UP);	/* peeling up to the proband */

      pConnector = pConnector->pNextConnector;
    }

  /* we are done with the up or down linked nuclear families or we are a leave
   * or top. ready for likelihood computation for this nuclear family 
   * save the original genotype list first */
  memcpy (&pProband->ppSavedGenotypeList[0],
	  &pProband->ppGenotypeList[0],
	  sizeof (Genotype *) * originalLocusList.numLocus);
  memcpy (&pProband->pSavedNumGenotype[0],
	  &pProband->pNumGenotype[0],
	  sizeof (int) * originalLocusList.numLocus);

  KLOG (LOGPARENTALPAIR, LOGDEBUG, "\t Proband (%s) haplotype: \n",
	pProband->sID);
  if (pProband->ppHaplotype == NULL)
    {
      pProband->ppHaplotype = MALLOC ("pProband->ppHaplotype",
				      sizeof (Genotype *) * sizeof (int) *
				      originalLocusList.numLocus);
    }

#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      /* build the constant polynomial */
      penetrancePolynomial = constantExp (1);
      weightPolynomial = constantExp (1);
      loop_proband_genotype (pNucFam, pProband, peelingDirection, 0, 0,
			     penetrancePolynomial, weightPolynomial);
    }
  else
    loop_proband_genotype (pNucFam, pProband, peelingDirection, 0, 0,
			   &penetrance, &weight);
#else

  loop_proband_genotype (pNucFam, pProband, peelingDirection, 0, 0,
			 &penetrance, &weight);
#endif

  /* copy back the genotypes for the proband */
  memcpy (&pProband->ppGenotypeList[0],
	  &pProband->ppSavedGenotypeList[0],
	  sizeof (Genotype *) * originalLocusList.numLocus);
  memcpy (&pProband->pNumGenotype[0],
	  &pProband->pSavedNumGenotype[0],
	  sizeof (int) * originalLocusList.numLocus);


  return 0;
}

/* This function only works on one nuclear family with given proband 
 * and go though each locus of the proband */
int
loop_proband_genotype (NuclearFamily * pNucFam, Person * pProband,
		       int peelingDirection, int locus, int multiLocusIndex,
		       void *penetrance, void *weight)
{
  int origLocus = locusList->pLocusIndex[locus];
  Locus *pLocus =
    originalLocusList.ppLocusList[locusList->pLocusIndex[locus]];
  int position;			/* genotype position */
  Genotype *pNextGenotype;
  Genotype *pGenotype;
  int numGenotype;
  int multiLocusIndex2;
  double newPenetrance;
  double newWeight;

#ifndef NO_POLYNOMIAL
  Polynomial *newPenetrancePolynomial = NULL;
  Polynomial *newWeightPolynomial = NULL;
#endif

  int origLocus1, origLocus2;
  Locus *pLocusTemp;
  LDLoci *pLDLoci;
  AlleleSet *pAlleleSet1, *pAlleleSet2;
  int allele1, allele2;
  int alleleID1, alleleID2;
  int j, k, l;
  double freq = 0;
#ifndef NO_POLYNOMIAL
  Polynomial *freqPolynomial = NULL;
#endif

  /* we loop over the genotypes of the proband to condition the 
   * likelihood calculation on it */
  numGenotype = pProband->pSavedNumGenotype[origLocus];
  multiLocusIndex2 = multiLocusIndex * numGenotype;
  pGenotype = pProband->ppSavedGenotypeList[origLocus];
  while (pGenotype != NULL)
    {
      pProband->ppHaplotype[origLocus] = pGenotype;
      KLOG (LOGPARENTALPAIR, LOGDEBUG, "\t\t\t %d|%d \n",
	    pGenotype->allele[DAD], pGenotype->allele[MOM]);
      pProband->ppGenotypeList[origLocus] = pGenotype;
      pNextGenotype = pGenotype->pNext;
      /* temporarilly set the next pointer to NULL so to restrict
       * the genotype on the proband to current genotype only */
      pGenotype->pNext = NULL;
      pProband->pNumGenotype[origLocus] = 1;
      position = pGenotype->position;
      multiLocusIndex = multiLocusIndex2 + position;

      /* if this proband is a founder - no DAD also means no MOM */
      if (pProband->pParents[DAD] == NULL)
	{
	  if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
	    {
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  newWeightPolynomial = timesExp (2, (Polynomial *) weight, 1,
						  pGenotype->weightPolynomial,
						  1, 0);
		}
	      else
		{
		  newWeight = *(double *) weight *pGenotype->weight;
		}
#else
	      newWeight = *(double *) weight *pGenotype->weight;
#endif
	    }
	  else
	    {
	      /* Linkage Disequilibrium */
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  newWeightPolynomial = (Polynomial *) weight;
		}
	      else
		{
		  newWeight = *(double *) weight;
		}
#else
	      newWeight = *(double *) weight;
#endif
	      /* instead of using genotype weight at each locus
	       * we need to use haplotype frequencies */
	      /* see whether current locus is in LD with previous locus */
	      if (locus > 0)
		{
		  origLocus1 = locusList->pLocusIndex[locus - 1];
		  origLocus2 = locusList->pLocusIndex[locus];
		  pLocusTemp = originalLocusList.ppLocusList[origLocus1];
		  pLocus = originalLocusList.ppLocusList[origLocus2];
		  /* find the parameter values for these two loci */
		  pLDLoci = find_LD_loci (origLocus1, origLocus2);
		  KASSERT (pLDLoci != NULL,
			   "Can't find LD parameter between loci %d,%d.\n",
			   origLocus1, origLocus2);
		  /* now find the corresponding haplotype frequencies - 
		   * both paternal and maternal */
		  for (j = DAD; j <= MOM; j++)
		    {
		      alleleID1 =
			pProband->ppHaplotype[origLocus1]->allele[j];
		      alleleID2 =
			pProband->ppHaplotype[origLocus2]->allele[j];
		      pAlleleSet1 =
			pLocusTemp->ppAlleleSetList[alleleID1 - 1];
		      pAlleleSet2 = pLocus->ppAlleleSetList[alleleID2 - 1];
#ifndef NO_POLYNOMIAL
		      if (modelOptions.polynomial == TRUE)
			{
			  freqPolynomial = constantExp (0);
			}
		      else
			{
			  freq = 0;
			}
#else
		      freq = 0;
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
				  char vName[100];
				  sprintf (vName, "ppHaploFreq[%d][%d]",
					   allele1 - 1, allele2 - 1);
				  freqPolynomial =
				    plusExp (2, 1.0, freqPolynomial, 1.0,
					     variableExp (&pLDLoci->
							  ppHaploFreq[allele1
								      -
								      1]
							  [allele2 - 1],
							  NULL, 'D', vName),
					     1);
				}
			      else
				freq +=
				  pLDLoci->ppHaploFreq[allele1 - 1][allele2 -
								    1];
#else

			      freq +=
				pLDLoci->ppHaploFreq[allele1 - 1][allele2 -
								  1];
#endif
			    }
			}
#ifndef NO_POLYNOMIAL
		      if (modelOptions.polynomial == TRUE)
			{
			  newWeightPolynomial = timesExp (2,
							  newWeightPolynomial,
							  1,
							  freqPolynomial, 1,
							  1);
			}
		      else
			{
			  newWeight *= freq;
			}
#else
		      newWeight *= freq;
#endif
		      if (modelOptions.sexLinked && pProband->sex + 1 == MALE)
			break;
		    }		/* end of looping paternal and maternal haplotypes */
		}		/* if locus > 0 for LD case */

	    }			/* end of LD case */

	}			/* end of founder case */
      else
	{
	  /* for non-founders, just keep the previous weight */
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      newWeightPolynomial = (Polynomial *) weight;
	    }
	  else
	    {
	      newWeight = *(double *) weight;
	    }
#else
	  newWeight = *(double *) weight;
#endif
	}

#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	{
	  newPenetrancePolynomial = timesExp (2, (Polynomial *) penetrance, 1,
					      pGenotype->penetrancePolynomial,
					      1, 0);
	}
      else
	{
	  newPenetrance = *(double *) penetrance *pGenotype->penetrance;
	}
#else
      newPenetrance = *(double *) penetrance *pGenotype->penetrance;
#endif

      if (locus < locusList->numLocus - 1)
	{
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      loop_proband_genotype (pNucFam, pProband, peelingDirection,
				     locus + 1, multiLocusIndex,
				     newPenetrancePolynomial,
				     newWeightPolynomial);
	    }
	  else
	    {
	      loop_proband_genotype (pNucFam, pProband, peelingDirection,
				     locus + 1, multiLocusIndex,
				     &newPenetrance, &newWeight);
	    }
#else
	  loop_proband_genotype (pNucFam, pProband, peelingDirection,
				 locus + 1, multiLocusIndex,
				 &newPenetrance, &newWeight);
#endif
	}
      else
	{
	  /* we have got the entire multi-locus genotype for the proband */
	  compute_nuclear_family_likelihood (pNucFam, pProband,
					     peelingDirection);
	  /* store the likelihood in the corresponding flattened array */
	  if (pProband->pLikelihood[multiLocusIndex].touchedFlag == FALSE)
	    {
	      pProband->pLikelihood[multiLocusIndex].touchedFlag = TRUE;
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  pProband->pLikelihood[multiLocusIndex].
		    likelihoodPolynomial = newPenetrancePolynomial;
		}
	      else
		{
		  /* need to update the penetrance factors */
		  pProband->pLikelihood[multiLocusIndex].likelihood =
		    newPenetrance;
		}
#else
	      /* need to update the penetrance factors */
	      pProband->pLikelihood[multiLocusIndex].likelihood =
		newPenetrance;
#endif
	    }
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      pProband->pLikelihood[multiLocusIndex].likelihoodPolynomial =
		timesExp (2,
			  pProband->pLikelihood[multiLocusIndex].
			  likelihoodPolynomial, 1,
			  pNucFam->likelihoodPolynomial, 1, 1);
//        fprintf(stderr,"Likelihood for this entire multi-locus genotype %f %f\n",
//                evaluateValue(pNucFam->likelihoodPolynomial),
//                evaluateValue(pProband->pLikelihood[multiLocusIndex].likelihoodPolynomial));
	    }
	  else
	    {
	      pProband->pLikelihood[multiLocusIndex].likelihood *=
		pNucFam->likelihood;
	      /*  fprintf (stderr,
	         "Likelihood for this entire multi-locus genotype %f %f\n",
	         pNucFam->likelihood,
	         pProband->pLikelihood[multiLocusIndex].likelihood);
	       */
	    }
#else
	  pProband->pLikelihood[multiLocusIndex].likelihood *=
	    pNucFam->likelihood;
#endif
	  if (pProband->pParents[DAD] == NULL)
	    {
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  pProband->pLikelihood[multiLocusIndex].weightPolynomial =
		    newWeightPolynomial;
//          fprintf(stderr,"Weight=%e penetrance=%e\n",evaluateValue(newWeightPolynomial),
//                         evaluateValue(newPenetrancePolynomial));
		}
	      else
		{
		  pProband->pLikelihood[multiLocusIndex].weight = newWeight;
//          fprintf(stderr,"Weight=%e penetrance=%e\n",newWeight,
//                          newPenetrance);
		}
#else
	      pProband->pLikelihood[multiLocusIndex].weight = newWeight;
#endif
	    }
	}
      /*when we are done, need to restore the genotype linke list pointer */
      pGenotype->pNext = pNextGenotype;
      pGenotype = pNextGenotype;
    }

  return 0;
}

/* for now, we only use parental pair algorithm */
int
compute_nuclear_family_likelihood (NuclearFamily * pNucFam,
				   Person * pProband, int peelingDirection)
{
  int locus;
  //  HaplotypePair haplo;
  double weight[2] = { 1, 1 };
  double penetrance[2] = { 1, 1 };
#ifndef NO_POLYNOMIAL
  /* need to define some terms for the polynomail opertations */
  Polynomial *weightPolynomial[2];
  Polynomial *penetrancePolynomial[2];
  Polynomial *sumPolynomial;
#endif
  int multiLocusIndex[2] = { 0, 0 };
  double sum;

  int numHaplotypePair = 1;
  int numChild;
#if 0
  ParentalPair *pPair, *pNextPair;
  int i;
#endif

#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      weightPolynomial[0] = constantExp (1);
      weightPolynomial[1] = constantExp (1);
      penetrancePolynomial[0] = constantExp (1);
      penetrancePolynomial[1] = constantExp (1);
    }
#endif


  numChild = pNucFam->numChildren;
  /* first construct the parental pair for this nuclear family locus
   * by locus */
  numHaplotypePair = 1;
  for (locus = 0; locus < locusList->numLocus; locus++)
    {
      construct_parental_pair (pNucFam, pProband, locus);
      numHaplotypePair *= parentalPairSpace.pNumParentalPair[locus];
    }
  /* now we can construct haplotypes and get likelihood computed */

#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      sumPolynomial = constantExp (0);
    }
  else
    sum = 0;
#else
  sum = 0;
#endif


  KLOG (LOGPARENTALPAIR, LOGDEBUG, "Haplotype for nuclear family No. %d:\n",
	pNucFam->nuclearFamilyIndex);
  KLOG (LOGPARENTALPAIR, LOGDEBUG, "\t\t\t DAD(%s)\t\t\t MOM(%s)\n",
	pNucFam->pParents[DAD]->sID, pNucFam->pParents[MOM]->sID);
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      loop_parental_pair (pNucFam, pProband, 0, &parentalPairSpace,
			  multiLocusIndex, (void *) weightPolynomial,
			  (void *) penetrancePolynomial, &sumPolynomial);
//      fprintf(stderr,"Conditional likelihood for nuclear family %d is: %e\n",
//            pNucFam->nuclearFamilyIndex, evaluateValue(sumPolynomial));
    }
  else
    {
      loop_parental_pair (pNucFam, pProband, 0, &parentalPairSpace,
			  multiLocusIndex, (void *) weight,
			  (void *) penetrance, &sum);
      KLOG (LOGLIKELIHOOD, LOGDEBUG,
	    "Conditional likelihood for nuclear family %d is: %e\n",
	    pNucFam->nuclearFamilyIndex, sum);
    }
#else
  loop_parental_pair (pNucFam, pProband, 0, &parentalPairSpace,
		      multiLocusIndex, (void *) weight, (void *) penetrance,
		      &sum);
  KLOG (LOGLIKELIHOOD, LOGDEBUG,
	"Conditional likelihood for nuclear family %d is: %e\n",
	pNucFam->nuclearFamilyIndex, sum);
#endif

#if 0
  FREE ("haplo.ppParentalPairs", haplo.ppParentalPairs);
  FREE ("pNucFam->pHaplotypePairs", pNucFam->pHaplotypePairs);
  pNucFam->pHaplotypePairs = NULL;
  pNucFam->numHaplotypePairs = 0;

  /* free parental pair storage */
  for (locus = 0; locus < locusList->numLocus; locus++)
    {
      pNucFam->pNumParentalPair[locus] = 0;
      /* release memory of old parental pairs if exist already */
      pPair = pNucFam->ppParentalPair[locus];
      while (pPair != NULL)
	{
	  pNextPair = pPair->pNext;
	  for (i = 0; i < numChild; i++)
	    {
	      FREE ("pPair->pppChildGenoList[i]", pPair->pppChildGenoList[i]);
	    }
	  FREE ("pPair->pppChildGenoList", pPair->pppChildGenoList);
	  FREE ("pPair->pChildGenoLen", pPair->pChildGenoLen);
	  FREE ("pPair", pPair);
	  pPair = pNextPair;
	}
      pNucFam->ppParentalPair[locus] = NULL;
    }
#endif

#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      pNucFam->likelihoodPolynomial = sumPolynomial;
    }
  else
    pNucFam->likelihood = sum;
#else
  pNucFam->likelihood = sum;
#endif
  return 0;
}

/* using the parental pair algorithm to go through all possible 
 * parental pairs and calculate the likelihood of each parental pair 
 * nested looping is for the multi-locus */
int
loop_parental_pair (NuclearFamily * pNucFam,
		    Person * pProband,
		    int locus,
		    ParentalPairSpace * pHaplo,
		    int multiLocusIndex[2],
		    void *dWeight[2], void *dPenetrance[2], void *pParentSum)
{
  double sum;
  double prob;
  int chromosome[2] = { DAD, DAD };
  Locus *pLocus =
    originalLocusList.ppLocusList[locusList->pLocusIndex[locus]];
  Locus *pLocusTemp;
  ParentalPair *pPair;
  ParentalPair *pTempPair;
  int child;
  int k, l;
  int origLocus;
  int origLocus1, origLocus2;
  LDLoci *pLDLoci;
  int allele1, allele2;
  int alleleID1, alleleID2;
  AlleleSet *pAlleleSet1, *pAlleleSet2;
  int i, j;
  double freq = 0;
  int multiLocusIndex2[2];
  Person *pParent[2];
  int numGenotype[2];
  int xmissionIndex[2];
  double newWeight[2];		/* genotype weight */
  double newPenetrance[2];
  double childProduct = 1;
  double childPenetrance = 1;
  int numPair;
  int breakFlag = 0;

#ifndef NO_POLYNOMIAL
  Polynomial *freqPolynomial = NULL;
  Polynomial *sumPolynomial = NULL;
  Polynomial *probPolynomial = NULL;
  Polynomial *newWeightPolynomial[2] = { NULL, NULL };
  Polynomial *newPenetrancePolynomial[2] = { NULL, NULL };
  Polynomial *childProductPolynomial = NULL;
  Polynomial *childPenetrancePolynomial = NULL;
  if (modelOptions.polynomial == TRUE)
    {
      childPenetrancePolynomial = constantExp (1);
    }
#endif

  origLocus = locusList->pLocusIndex[locus];
  for (i = DAD; i <= MOM; i++)
    {
      pParent[i] = pNucFam->pParents[i];
      /* find the max number of possible genotypes for this parent */
      numGenotype[i] = pNucFam->pParents[i]->pSavedNumGenotype[origLocus];
      multiLocusIndex2[i] = multiLocusIndex[i] * numGenotype[i];
    }
  //  pPair = pNucFam->ppParentalPair[locus];
  numPair = -1;
  while ((numPair += 1) < pHaplo->pNumParentalPair[locus])
    {
      //      pHaplo->ppParentalPairs[locus] = pPair;
      pPair = &pHaplo->ppParentalPair[locus][numPair];
      pHaplo->pParentalPairInd[locus] = numPair;
      for (i = DAD; i <= MOM; i++)
	{
	  multiLocusIndex[i] =
	    multiLocusIndex2[i] + pPair->pGenotype[i]->position;
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {

	      newPenetrancePolynomial[i] =
		timesExp (2, (Polynomial *) dPenetrance[i], 1,
			  pPair->pGenotype[i]->penetrancePolynomial, 1, 0);

	    }
	  else
	    {
	      newPenetrance[i] = *((double *) dPenetrance + i) *
		pPair->pGenotype[i]->penetrance;
	    }
#else
	  newPenetrance[i] = *((double *) dPenetrance + i) *
	    pPair->pGenotype[i]->penetrance;
#endif


	  if (pParent[i] == pProband)
	    {
	      /* we are conditionaling on proband, so it's 1 */
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  newWeightPolynomial[i] = constantExp (1);
		}
	      else
		{
		  newWeight[i] = 1;
		}
#else
	      newWeight[i] = 1;
#endif
	    }
	  else
	    {
	      /* we don't have weight information of this individual yet */
	      if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
		{
#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      newWeightPolynomial[i] =
			timesExp (2, (Polynomial *) dWeight[i], 1,
				  pPair->pGenotype[i]->weightPolynomial, 1,
				  0);
		    }
		  else
		    {
		      newWeight[i] =
			*((double *) dWeight +
			  i) * pPair->pGenotype[i]->weight;
		    }
#else
		  newWeight[i] =
		    *((double *) dWeight + i) * pPair->pGenotype[i]->weight;
#endif
		}
	      else
		{		/* Linkage Disequilibrium */
#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      newWeightPolynomial[i] = (Polynomial *) dWeight[i];
		    }
		  else
		    {
		      newWeight[i] = *((double *) dWeight + i);
		    }
#else
		  newWeight[i] = *((double *) dWeight + i);
#endif
		  /* instead of using genotype weight at each locus
		   * we need to use haplotype frequencies */
		  /* see whether current locus is in LD with previous locus */
		  if (locus > 0)
		    {
		      origLocus1 = locusList->pLocusIndex[locus - 1];
		      origLocus2 = locusList->pLocusIndex[locus];
		      pLocusTemp = originalLocusList.ppLocusList[origLocus1];
		      /* find the parameter values for these two loci */
		      pLDLoci = find_LD_loci (origLocus1, origLocus2);
		      KASSERT (pLDLoci != NULL,
			       "Can't find LD parameter between loci %d,%d.\n",
			       origLocus1, origLocus2);
		      /* now find the corresponding haplotype frequency : 2 haplotypes 
		       * paternal haplotype & maternal haplotype */
		      pTempPair =
			&pHaplo->ppParentalPair[locus -
						1][pHaplo->
						   pParentalPairInd[locus -
								    1]];
		      for (j = DAD; j <= MOM; j++)
			{
			  alleleID2 = pPair->pGenotype[i]->allele[j];
			  alleleID1 = pTempPair->pGenotype[i]->allele[j];
			  pAlleleSet1 =
			    pLocusTemp->ppAlleleSetList[alleleID1 - 1];
			  pAlleleSet2 =
			    pLocus->ppAlleleSetList[alleleID2 - 1];
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    {
			      freqPolynomial = constantExp (0);
			    }
			  else
			    freq = 0;
#else
			  freq = 0;
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
				      char vName[100];
				      sprintf (vName, "ppHaploFreq[%d][%d]",
					       allele1 - 1, allele2 - 1);
				      freqPolynomial =
					plusExp (2, 1.0, freqPolynomial, 1.0,
						 variableExp (&pLDLoci->
							      ppHaploFreq
							      [allele1 -
							       1][allele2 -
								  1], NULL,
							      'D', vName), 1);
				    }
				  else
				    freq +=
				      pLDLoci->ppHaploFreq[allele1 -
							   1][allele2 - 1];
#else
				  freq +=
				    pLDLoci->ppHaploFreq[allele1 -
							 1][allele2 - 1];
#endif
				}
			    }
#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    {
			      newWeightPolynomial[i] =
				timesExp (2, newWeightPolynomial[i], 1,
					  freqPolynomial, 1, 0);
			    }
			  else
			    {
			      newWeight[i] *= freq;
			    }
#else
			  newWeight[i] *= freq;
#endif
			}	/* end of looping paternal and maternal haplotypes */
		    }		/* if locus > 0 for LD case */

		}		/*  LD case */
	    }			/* not a proband */
	}			/* looping dad and mom genotypes */


      if (locus < locusList->numLocus - 1)
	{
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {

	      loop_parental_pair (pNucFam, pProband, locus + 1, pHaplo,
				  multiLocusIndex,
				  (void *) newWeightPolynomial,
				  (void *) newPenetrancePolynomial,
				  pParentSum);
	    }
	  else
	    {
	      loop_parental_pair (pNucFam, pProband, locus + 1, pHaplo,
				  multiLocusIndex,
				  (void *) newWeight, (void *) newPenetrance,
				  pParentSum);
	    }
#else
	  loop_parental_pair (pNucFam, pProband, locus + 1, pHaplo,
			      multiLocusIndex,
			      (void *) newWeight, (void *) newPenetrance,
			      pParentSum);
#endif
	}
      else
	{
	  /* check to see whether we already have likelihood calculated
	   * for this individual. If this individual is  a connector to
	   * another nuclear family, then we should have already got the
	   * likelihood stored for this person */

#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      pHaplo->likelihoodPolynomial = constantExp (1);
	    }
	  else
	    {
	      pHaplo->likelihood = 1;
	    }
#else
	  pHaplo->likelihood = 1;
#endif

	  breakFlag = 0;
	  for (i = DAD; i <= MOM; i++)
	    {
	      /* we will update proband later */
	      if (pParent[i] != pProband
		  && pParent[i]->pLikelihood[multiLocusIndex[i]].
		  touchedFlag != TRUE)
		{
		  /* stored the weight */
		  pParent[i]->pLikelihood[multiLocusIndex[i]].touchedFlag =
		    TRUE;
#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      pParent[i]->pLikelihood[multiLocusIndex[i]].
			likelihoodPolynomial = newPenetrancePolynomial[i];
		    }
		  else
		    {
		      pParent[i]->pLikelihood[multiLocusIndex[i]].likelihood =
			newPenetrance[i];
		    }
#else
		  pParent[i]->pLikelihood[multiLocusIndex[i]].likelihood =
		    newPenetrance[i];
#endif
		}		/* getting the weight for this individual */

	      if (pParent[i] != pProband)
		{
#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      pHaplo->likelihoodPolynomial =
			timesExp (2,
				  pHaplo->likelihoodPolynomial,
				  1,
				  pParent[i]->pLikelihood[multiLocusIndex[i]].
				  likelihoodPolynomial, 1, 1);
		    }
		  else
		    {
		      pHaplo->likelihood *=
			pParent[i]->pLikelihood[multiLocusIndex[i]].
			likelihood;
		      if (newPenetrance[i] == 0)
			{
			  breakFlag = 1;
			}
		    }
#else
		  pHaplo->likelihood *=
		    pParent[i]->pLikelihood[multiLocusIndex[i]].likelihood;
		  if (newPenetrance[i] == 0)
		    {
		      breakFlag = 1;
		    }
#endif
		}
	      else
		{
		}
	      if (pParent[i]->pParents[DAD] == NULL && pParent[i] != pProband)
		{
#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      pHaplo->likelihoodPolynomial =
			timesExp (2, pHaplo->likelihoodPolynomial, 1,
				  newWeightPolynomial[i], 1, 1);
		    }
		  else
		    {
		      pHaplo->likelihood *= newWeight[i];
		    }
#else
		  pHaplo->likelihood *= newWeight[i];
#endif
		}
	    }			/* update for each parent */

	  if (breakFlag == 1)
	    continue;
	  /* we got a complete multi-locus haplotype for the parents 
	   * now we need to loop over children for multi-locus */
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      childProductPolynomial = constantExp (1);
	    }
	  else
	    childProduct = 1;
#else
	  childProduct = 1;
#endif

	  for (child = 0; child < pNucFam->numChildren; child++)
	    {
	      xmissionIndex[DAD] = 0;
	      xmissionIndex[MOM] = 0;
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  sumPolynomial = constantExp (0);
		  probPolynomial = constantExp (1);
		  childPenetrancePolynomial = constantExp (1);
		  loop_child_multi_locus_genotype (pNucFam->
						   ppChildrenList[child],
						   pProband, pHaplo, child, 0,
						   0,
						   childPenetrancePolynomial,
						   &sumPolynomial,
						   probPolynomial,
						   xmissionIndex, chromosome);
		  childProductPolynomial =
		    timesExp (2, childProductPolynomial, 1, sumPolynomial, 1,
			      1);

//          fprintf(stderr,"Child No. %d\n",child);
//          expPrinting(sumPolynomial);
//          fprintf(stderr,"\n  value=%f\n",evaluateValue(sumPolynomial));

		}
	      else
		{
		  sum = 0;
		  prob = 1;
		  childPenetrance = 1;
		  /* get transmission probability */
		  loop_child_multi_locus_genotype (pNucFam->
						   ppChildrenList[child],
						   pProband, pHaplo, child, 0,
						   0, &childPenetrance, &sum,
						   &prob,
						   xmissionIndex, chromosome);

		  childProduct *= sum;
//          fprintf(stderr,"Child No. %d\n",child);
//          fprintf(stderr,"likelihood=%f\n",sum);

		}
#else
	      sum = 0;
	      prob = 1;
	      childPenetrance = 1;
	      /* get transmission probability */
	      loop_child_multi_locus_genotype (pNucFam->ppChildrenList[child],
					       pProband, pHaplo, child, 0, 0,
					       &childPenetrance,
					       &sum, &prob,
					       xmissionIndex, chromosome);

	      childProduct *= sum;
#endif
	    }
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {

	      pHaplo->likelihoodPolynomial =
		timesExp (2, pHaplo->likelihoodPolynomial, 1,
			  childProductPolynomial, 1, 1);
	      *(Polynomial **) pParentSum =
		plusExp (2, 1.0, *(Polynomial **) pParentSum, 1.0,
			 pHaplo->likelihoodPolynomial, 1);
	    }
	  else
	    {
	      pHaplo->likelihood *= childProduct;
	      *(double *) pParentSum += pHaplo->likelihood;
	    }
#else
	  pHaplo->likelihood *= childProduct;
	  *(double *) pParentSum += pHaplo->likelihood;
#endif
#if 0
	  memcpy (&pNucFam->pHaplotypePairs[pNucFam->numHaplotypePairs],
		  pHaplo, sizeof (HaplotypePair));
	  pNucFam->numHaplotypePairs++;
#endif
	  //haploProb = get_haplo_probability(pHaplo);
	  for (k = 0; k < locusList->numLocus; k++)
	    {
	      pTempPair =
		&pHaplo->ppParentalPair[k][pHaplo->pParentalPairInd[k]];
	      KLOG (LOGPARENTALPAIR, LOGDEBUG,
		    "\t %4d-> %4d|%-4d \t %4d|%-4d\n",
		    locusList->pLocusIndex[k],
		    pTempPair->pGenotype[DAD]->allele[DAD],
		    pTempPair->pGenotype[DAD]->allele[MOM],
		    pTempPair->pGenotype[MOM]->allele[DAD],
		    pTempPair->pGenotype[MOM]->allele[MOM]);
	    }
	  KLOG (LOGPARENTALPAIR, LOGDEBUG, "\n");

#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
//      expPrinting(pHaplo->likelihoodPolynomial);
//      KLOG (LOGPARENTALPAIR, LOGDEBUG, "\t Likelihood = %e\n",
//              evaluateValue (pHaplo->likelihoodPolynomial));
	    }
	  else
	    {
	      KLOG (LOGPARENTALPAIR, LOGDEBUG, "\t Likelihood = %e\n",
		    pHaplo->likelihood);
	    }
#else
	  KLOG (LOGPARENTALPAIR, LOGDEBUG, "\t Likelihood = %e\n",
		pHaplo->likelihood);
#endif
	}			/* end of processing one parental pair */
      //      numPair++; /* move to next parental pair */
    }

  return 0;
}

int
loop_child_multi_locus_genotype (Person * pChild, Person * pProband,
				 ParentalPairSpace * pHaplo, int child,
				 int locus, int multiLocusIndex,
				 void *pPenetrance,
				 void *pSum, void *pProb,
				 int xmissionIndex[2], int chromosome[2])
{
  int i;
  Genotype *pGeno;
  ParentalPair *pParentalPair =
    &pHaplo->ppParentalPair[locus][pHaplo->pParentalPairInd[locus]];
  int parent;
  //  int alleleSetLen = originalLocusList.alleleSetLen;
  int origLocus = locusList->pLocusIndex[locus];
  Locus *pLocus =
    originalLocusList.ppLocusList[locusList->pLocusIndex[locus]];

  double newProb;

#ifndef NO_POLYNOMIAL
  Polynomial *newProbPolynomial = NULL;
#endif

  int newChromosome[2];
  int newChromosome2[2];
  int numGenotype;
  double newPenetrance;
#ifndef NO_POLYNOMIAL
  Polynomial *newPenetrancePolynomial = NULL;
#endif
  int multiLocusIndex2;
  int newXmissionIndex[2];

  numGenotype = pChild->pSavedNumGenotype[origLocus];
  multiLocusIndex2 = multiLocusIndex * numGenotype;
  xmissionIndex[0] *= 3;
  xmissionIndex[1] *= 3;
  for (i = 0; i < pParentalPair->pChildGenoLen[child]; i++)
    {
      pGeno = pParentalPair->pppChildGenoList[child][i];
      KLOG (LOGLIKELIHOOD, LOGDEBUG,
	    "\t child %s locus %4d -> %4d|%-4d \n",
	    pChild->sID, locusList->pLocusIndex[locus], pGeno->allele[DAD], pGeno->allele[MOM]);
      multiLocusIndex = multiLocusIndex2 + pGeno->position;
      newChromosome2[DAD] = chromosome[DAD];
      newChromosome2[MOM] = chromosome[MOM];
#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	{
	  if (pLocus->locusType == LOCUS_TYPE_TRAIT)
	    {
	      newPenetrancePolynomial =
		timesExp (2, (Polynomial *) pPenetrance, 1,
			  pGeno->penetrancePolynomial, 1, 0);

	    }
	  else
	    {
	      newPenetrancePolynomial = (Polynomial *) pPenetrance;
	    }
	  newProbPolynomial = (Polynomial *) pProb;

//          fprintf(stderr,"locus=%d newPenetrance=%f newProb=%f locus=%d locusList->numLocus=%d i=%d\n",
//                          locus,evaluateValue(newPenetrancePolynomial),evaluateValue(newProbPolynomial),
//                          locus,locusList->numLocus,i);

	}
      else
	{
	  if (pLocus->locusType == LOCUS_TYPE_TRAIT)
	    {
	      newPenetrance = *(double *) pPenetrance *pGeno->penetrance;
	    }
	  else
	    newPenetrance = *(double *) pPenetrance;
	  newProb = *(double *) pProb;

//          fprintf(stderr,"locus=%d newPenetrance=%f newProb=%f locus=%d locusList->numLocus=%d i=%d\n",
//                          locus,newPenetrance,newProb,
//                        locus,locusList->numLocus,i); 


	}
#else
      if (pLocus->locusType == LOCUS_TYPE_TRAIT)
	{
	  newPenetrance = *(double *) pPenetrance *pGeno->penetrance;
	}
      else
	newPenetrance = *(double *) pPenetrance;
      newProb = *(double *) pProb;
#endif
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
	    }			/* looping paternal and maternal chromsome */
	  if (locus < locusList->numLocus - 1)
	    {
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
//                  fprintf(stderr,"AAlocus=%d newPen=%f pSum=%f newProb=%f i=%d\n",
//                                  locus,evaluateValue(newPenetrancePolynomial),
//                                evaluateValue(*(Polynomial **)pSum),evaluateValue(newProbPolynomial),i);      
		  loop_child_multi_locus_genotype (pChild, pProband, pHaplo,
						   child, locus + 1,
						   multiLocusIndex,
						   newPenetrancePolynomial,
						   pSum, newProbPolynomial,
						   newXmissionIndex,
						   newChromosome2);
//                   fprintf(stderr,"BBlocus=%d newPen=%f pSum=%f newProb=%f\n",
//                                locus,evaluateValue(newPenetrancePolynomial),
//                                evaluateValue(*(Polynomial **)pSum),evaluateValue(newProbPolynomial));                   
		}
	      else
		{
//                   fprintf(stderr,"AAlocus=%d newPen=%f, pSum=%f newProb=%f i=%d\n",
//                                   locus, newPenetrance, *((double *)pSum),newProb,i);
		  loop_child_multi_locus_genotype (pChild, pProband, pHaplo,
						   child, locus + 1,
						   multiLocusIndex,
						   &newPenetrance, pSum,
						   &newProb,
						   newXmissionIndex,
						   newChromosome2);
//                   fprintf(stderr,"BBlocus=%d newPen=%f, pSum=%f newProb=%f\n",
//                                   locus,newPenetrance, *((double *)pSum),newProb);
		}
#else
	      loop_child_multi_locus_genotype (pChild, pProband, pHaplo,
					       child, locus + 1,
					       multiLocusIndex,
					       &newPenetrance, pSum, &newProb,
					       newXmissionIndex,
					       newChromosome2);
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
			      xmissionMatrix[newXmissionIndex[DAD]].
			      probPoly[DAD + 1], 1,
			      xmissionMatrix[newXmissionIndex[MOM]].
			      probPoly[MOM + 1], 1, 0);
//                 fprintf(stderr,"newProb=%f newPenetrance=%f\n",
//                               evaluateValue(newProbPolynomial),evaluateValue(newPenetrancePolynomial));
		}
	      else
		{
		  newProb =
		    xmissionMatrix[newXmissionIndex[DAD]].prob[DAD +
							       1] *
		    xmissionMatrix[newXmissionIndex[MOM]].prob[MOM + 1];
		  KLOG (LOGLIKELIHOOD, LOGDEBUG,
			"\t xmission prob: %f = %f * %f\n", newProb,
			xmissionMatrix[newXmissionIndex[DAD]].prob[DAD + 1],
			xmissionMatrix[newXmissionIndex[MOM]].prob[MOM + 1]);
//                 fprintf(stderr,"newProb=%f newPenetrance=%f\n",newProb,newPenetrance);
		}
#else
	      newProb = xmissionMatrix[newXmissionIndex[DAD]].prob[DAD + 1] *
		xmissionMatrix[newXmissionIndex[MOM]].prob[MOM + 1];
	      KLOG (LOGLIKELIHOOD, LOGDEBUG,
		    "\t xmission prob: %f = %f * %f\n",
		    newProb,
		    xmissionMatrix[newXmissionIndex[DAD]].prob[DAD + 1],
		    xmissionMatrix[newXmissionIndex[MOM]].prob[MOM + 1]);
#endif



	      /* we have completed one multilocus genotype for this child */
	      /* check whether we already have some information about this kid
	       * we should have if this kid is a connector to another nuclear
	       * family we have processed before */
	      if (pChild != pProband &&
		  pChild->pLikelihood[multiLocusIndex].touchedFlag != TRUE)
		{
		  pChild->pLikelihood[multiLocusIndex].touchedFlag = TRUE;
#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      pChild->pLikelihood[multiLocusIndex].
			likelihoodPolynomial = newPenetrancePolynomial;
		    }
		  else
		    {
		      pChild->pLikelihood[multiLocusIndex].likelihood =
			newPenetrance;
		    }
#else
		  pChild->pLikelihood[multiLocusIndex].likelihood =
		    newPenetrance;
#endif
		}
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
//                  fprintf(stderr,"pSum=%f, newProb=%f child.likelihood=%f pChild == pProband=%d\n",
//                                evaluateValue(*(Polynomial **) pSum),evaluateValue(newProbPolynomial),
//                                evaluateValue(pChild->pLikelihood[multiLocusIndex].likelihoodPolynomial),(pChild == pProband));
		  if (pChild == pProband)
		    {
		      *(Polynomial **) pSum =
			plusExp (2, 1.0, *(Polynomial **) pSum, 1.0,
				 newProbPolynomial, 1);

		    }
		  else
		    {
		      *(Polynomial **) pSum = plusExp (2, 1.0, *(Polynomial **) pSum, 1.0, timesExp (2, newProbPolynomial, 1, pChild->pLikelihood[multiLocusIndex].likelihoodPolynomial, 1, 0)	//end of timesExp
						       , 1);	//end of plusExp
		    }
//                  fprintf(stderr,"pSum=%f\n",evaluateValue(*(Polynomial **) pSum));
		}
	      else
		{

//                  fprintf(stderr,"pSum=%f, newProb=%f child.likelihood=%f pChild == pProband=%d\n",
//                                  *(double *) pSum,newProb,pChild->pLikelihood[multiLocusIndex].likelihood,pChild == pProband);
		  if (pChild == pProband)
		    *(double *) pSum += newProb;
		  else
		    *(double *) pSum +=
		      newProb *
		      pChild->pLikelihood[multiLocusIndex].likelihood;
//                   fprintf(stderr,"pSum=%f\n",*((double *)pSum));
		}
#else
	      if (pChild == pProband)
		*(double *) pSum += newProb;
	      else
		*(double *) pSum +=
		  newProb * pChild->pLikelihood[multiLocusIndex].likelihood;
#endif
	    }
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
			  int lastHetLoc, int prevPattern,
			  int loc)
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
	      if (prevPattern != 3 )
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
                              sprintf (vName1, "theta%d_%d", i,loc);
			      newProbPoly[i] =
				timesExp (2, newProbPoly[i], 1,
					     plusExp (2, 1.0, constantExp (1.0),
					 	        -1.0,
						        variableExp (&locusList->
								pPrevLocusDistance
								[i][loc], NULL,
								'D', vName1), 0),
					  1, 0);
			    }
			  else
			    {
			      newProb[i] *= (1 - locusList->pPrevLocusDistance[i][loc]);
			    }
#else
			  newProb[i] *= (1 - locusList->pPrevLocusDistance[i][loc]);
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
				sprintf (vName1, "theta%d_%d", i,loc);
				newProbPoly[i] = timesExp (2, newProbPoly[i], 1,
							      variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),1,0);
							      
                            }
                          else
			  {
                            newProb[i] *= locusList->pPrevLocusDistance[i][loc];
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
		      if(pattern==1)
			{
			  /* paternal inheritance for this locus
			     either no recombination from previous paternal strand 
			     or recombination from previous maternal strand */
#ifndef NO_POLYNOMIAL
                          if (modelOptions.polynomial == TRUE)
                            {
                                sprintf (vName1, "theta%d_%d", i,loc);
                                newProbPoly[i] = plusExp(2, 
                                  1.0, timesExp(2, (Polynomial *) prob[i], 
 						   1,
                                                   plusExp(2,1.0, constantExp(1.0),
							    -1.0, variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),0),
						   1,0),
				  1.0, timesExp(2, (Polynomial *) prob2[i],
						   1,
						   variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),
						   1,0),
				  0);
			    }
			  else
			    newProb[i] = *((double *)prob[i]) * 
			      (1 - locusList->pPrevLocusDistance[i][loc]) +
			      *((double *) prob2[i]) * locusList->pPrevLocusDistance[i][loc];
#else
                          newProb[i] = *((double *)prob[i]) *
                            (1 - locusList->pPrevLocusDistance[i][loc]) +
                            *((double *) prob2[i]) * locusList->pPrevLocusDistance[i][loc];
#endif 
			}
		      else
			{
			  /* has to be maternal */
#ifndef NO_POLYNOMIAL
                          if (modelOptions.polynomial == TRUE)
                            {
                                sprintf (vName1, "theta%d_%d", i,loc);
                                newProbPoly[i] = plusExp(2,
                                  1.0, timesExp(2, (Polynomial *) prob2[i], 
                                                   1,
                                                   plusExp(2,1.0, constantExp(1.0),
                                                            -1.0, variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),0),
                                                   1,0),
                                  1.0, timesExp(2, (Polynomial *) prob[i],
                                                   1,
                                                   variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),
                                                   1,0),
                                  0);

                            }
			  else			 
			    newProb[i] = *((double *)prob2[i]) *
			      (1 - locusList->pPrevLocusDistance[i][loc]) +
			      *((double *) prob[i]) * locusList->pPrevLocusDistance[i][loc];
			  
#else
			  newProb[i] = *((double *)prob2[i]) * 
			    (1 - locusList->pPrevLocusDistance[i][loc]) +
			    *((double *) prob[i]) * locusList->pPrevLocusDistance[i][loc];
#endif
			}
		    }
		      
		} /* end of prevPattern is homo and current pattern is het */
	    }/* end of prevHetLoc != -1 */
	  else
	    {
	      /* we don't have any het locus yet, this locus is the first het */
#ifndef NO_POLYNOMIAL
              if (modelOptions.polynomial == TRUE)
                {
                  for(i=0; i < 3; i++)
                    newProbPoly[i] = constantExp (0.5);
                }
              else
	      {
		for(i=0; i < 3; i++)
		  newProb[i] = 0.5;
              }
#else
	      for ( i=0; i < 3; i++)
		newProb[i] = 0.5; 
#endif
	    }
	  newLastHetLoc = loc; 


#ifndef NO_POLYNOMIAL
              if (modelOptions.polynomial == TRUE)
                {
		  for(i=0; i < 3; i++)
		    newHetProbPoly[i] = newProbPoly[i];
                }
	      else
                {
		  for(i=0; i < 3; i++)
		    newHetProbPtr[i] = newProbPtr[i];
                }
#else
	  for(i=0; i < 3; i++)
	    newHetProbPtr[i] = newProbPtr[i];
#endif


	} /* end of current pattern is not homo */
      else 
	{
 	  /* current pattern is homo */
	  if(lastHetLoc == -1)
	    /* nothing needs to be done for this locus */
	    ;
	  else
	    {
	      if(loc == totalLoci - 1)
		{
		  /* this is the last locus and it's homo, take the previous het locus */
		  for ( i=0; i < 3; i++)
		    {
#ifndef NO_POLYNOMIAL
	              if (modelOptions.polynomial == TRUE)
        	      {
			 newProbPoly[i]= (Polynomial *)hetProb[i];
		      }
 		      else
                         newProb[i] = *(double *)hetProb[i];
#else
		      newProb[i] = *(double *)hetProb[i];
#endif
		    }
		}
	      else
		{
		  if(prevPattern == 3)
		    {
		      /* previous locus pattern is homo */
		      for ( i=0; i < 3; i++)
			{

#ifndef NO_POLYNOMIAL
                          if (modelOptions.polynomial == TRUE)
                          {
                                sprintf (vName1, "theta%d_%d", i,loc);
                                newProbPoly[i] = plusExp(2,
                                  1.0, timesExp(2, (Polynomial *) prob[i],
                                                   1,
                                                   plusExp(2,1.0, constantExp(1.0),
                                                            -1.0, variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),0),
                                                   1,0),
                                  1.0, timesExp(2, (Polynomial *) prob2[i],
                                                   1,
                                                   variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),
                                                   1,0),
                                  0);
                                newProbPoly2[i] = plusExp(2,
                                  1.0, timesExp(2, (Polynomial *) prob2[i],
                                                   1,
                                                   plusExp(2,1.0, constantExp(1.0),
                                                            -1.0, variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),0),
                                                   1,0),
                                  1.0, timesExp(2, (Polynomial *) prob[i],
                                                   1,
                                                   variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),
                                                   1,0),
                                  0);
                            
                          }
                          else
			  {
                             newProb[i] = *(double *)prob[i] *
                               (1 - locusList->pPrevLocusDistance[i][loc]) +
                               *((double *) prob2[i]) * locusList->pPrevLocusDistance[i][loc];

                             newProb2[i] = *(double *)prob2[i] *
                               (1 - locusList->pPrevLocusDistance[i][loc]) +
                               *((double *) prob[i]) * locusList->pPrevLocusDistance[i][loc];
			  }
#else
			  newProb[i] = *(double *)prob[i] * 
			    (1 - locusList->pPrevLocusDistance[i][loc]) +
			    *((double *) prob2[i]) * locusList->pPrevLocusDistance[i][loc];
			  
			  newProb2[i] = *(double *)prob2[i] * 
			    (1 - locusList->pPrevLocusDistance[i][loc]) +
			    *((double *) prob[i]) * locusList->pPrevLocusDistance[i][loc];
#endif
			}
		    }
		  else 
		    {
		      for(i=0; i < 3; i++)
			{
			  if(prevPattern == 1)
			    {
#ifndef NO_POLYNOMIAL
                              if (modelOptions.polynomial == TRUE)
                              {
                                   sprintf (vName1, "theta%d_%d", i,loc);
                                   newProbPoly[i] = timesExp(2, (Polynomial *) prob[i],
                                                   1,
                                                   plusExp(2,1.0, constantExp(1.0),
                                                            -1.0, variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),0),
                                                   1,0);				   
                                   newProbPoly2[i] = timesExp(2, (Polynomial *) prob[i],
                                                   1,
                                                   variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),
                                                   1,0);
                              }
			      else
				{
	                              newProb[i] = *(double *)prob[i] *
	                                (1 - locusList->pPrevLocusDistance[i][loc]);
	                              newProb2[i] = *(double *)prob[i] *
	                                locusList->pPrevLocusDistance[i][loc];
				}
#else
			      newProb[i] = *(double *)prob[i] * 
				(1 - locusList->pPrevLocusDistance[i][loc]);
			      
			      newProb2[i] = *(double *)prob[i] * 
				locusList->pPrevLocusDistance[i][loc];
#endif
			    }
			  else
			    {
#ifndef NO_POLYNOMIAL
                              if (modelOptions.polynomial == TRUE)
                              {
                                   sprintf (vName1, "theta%d_%d", i,loc);
                                   newProbPoly2[i] = timesExp(2, (Polynomial *) prob[i],
                                                   1,
                                                   plusExp(2,1.0, constantExp(1.0),
                                                            -1.0, variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),0),
                                                   1,0);
                                   newProbPoly[i] = timesExp(2, (Polynomial *) prob[i],
                                                   1,
                                                   variableExp (&locusList->pPrevLocusDistance[i][loc],NULL,'D', vName1),
                                                   1,0);
			      }
			      else
			      {
                                  newProb2[i] = *(double *)prob[i] *
                                    (1 - locusList->pPrevLocusDistance[i][loc]);
                                  newProb[i] = *(double *)prob[i] *
                                    locusList->pPrevLocusDistance[i][loc];
			      }
#else
			      newProb2[i] = *(double *)prob[i] * 
				(1 - locusList->pPrevLocusDistance[i][loc]);
			      
			      newProb[i] = *(double *)prob[i] * 
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
		  pMatrix[newCellIndex].probPoly[i] = newProbPoly[i];
		}
	      else
		pMatrix[newCellIndex].prob[i] = newProb[i];
#else
	      pMatrix[newCellIndex].prob[i] = newProb[i];
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
					(void *) newProbPoly, (void *)newProbPoly2, 
					(void *) newHetProbPoly, 
					newCellIndex,
					newLastHetLoc,
					pattern, loc + 1);
	    }
	  else
	    populate_xmission_matrix (pMatrix, totalLoci,
				      (void *) newProbPtr, (void *)newProbPtr2, 
				      (void *) newHetProbPtr, 
				      newCellIndex,
				      newLastHetLoc,
				      pattern, loc + 1);
#else
	  populate_xmission_matrix (pMatrix, totalLoci,
				    newProbPtr, (void *)newProbPtr2, 
				    (void *) newHetProbPtr, 
				    newCellIndex,
				    newLastHetLoc,
				    pattern, loc + 1);
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
print_xmission_matrix(XMission *pMatrix, int totalLoci, int loc, int cellIndex, char *pID)
{
  int pattern; 
  int newCellIndex;
  int i;

  for(pattern=1; pattern <=3; pattern++)
    {
      newCellIndex = cellIndex*3 + pattern -1;
      if(pattern==1)
	{
	  pID[loc] = 'P';
	}
      else if (pattern==2)
	{
	  pID[loc] = 'M';
	}
      else
	{
	  pID[loc] = 'B'; 
	}
      if(loc != totalLoci -1)
	{
	  /* not complete multi-locus haplotype yet */
	  print_xmission_matrix(pMatrix, totalLoci, loc+1, newCellIndex, pID);
	}
      else
	{
	  /* print the xmission probability out */
	  for(i=0; i <= loc; i++)
	    {
	      fprintf(stderr, "%c", pID[i]);
	    }
	  fprintf(stderr, ":%f\n", pMatrix[newCellIndex].prob[0]);
	}
    }  
}
