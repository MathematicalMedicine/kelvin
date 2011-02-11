/**
#file iterationSupport.c

  The iterative half of the former kelvin main program.

  Iterates over a fixed grid of explicitly-specified trait space parameters 
  in order to average them out of the likelihood calculation. This is
  the original version of kelvin, and has been nearly completely replaced
  by a version that does dynamic integration to average-out the nuisance
  parameters (integrationSupport.c).

  Copyright &copy; 2010, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id: iterationSupport.c 15 2009-10-16 12:14:59Z whv001 $

*/
#include <memory.h>
#include <math.h>

#include "kelvin.h"
#include "kelvinGlobals.h"
#include "kelvinHandlers.h"     // For status request signal detection
#include "summary_result.h"
#include "kelvinWriteFiles.h"
#include "saveResults.h"
#include "trackProgress.h"
#include "ppl.h"

#include "iterationSupport.h"

extern struct swStopwatch *combinedComputeSW,  ///< Combined likelihood compute stopwatch
  *combinedBuildSW,      ///< Combined likelihood polynomial build stopwatch
  *overallSW;    ///< Overall stopwatch for the entire run.

/**

  Driver for iterative analysis.

  Original edifice of kelvin main program. Needs refactoring in the worst way, but
  our focus has moved on to the dynamic integration approach. This version is kept
  around for validation purposes as it can, for most types of analyses, operate on
  single set of trait space parameters just like all other linkage analysis programs.

  It can also do polyallelic LD and marker-to-marker analysis, capabilities not yet
  included in the dynamic integration approach.

  @author Yungui Huang - overall content.
  @author Bill Valentine-Cooper - progress tracking.

*/
void iterateMain ()
{
  ParamStruct paramSet;
  int numPositions;
  double log10AvgLR;

  char **markerNameList = NULL;

  int thresholdIdx = -1;
  double threshold = 0;
  double avgLR;
  double constraint;
  double log10_likelihood_null, log10_likelihood_alternative;
  int paramIdx = -1;
  int penIdx, gfreqInd, thetaInd;
  double pen_DD, pen_Dd, pen_dD, pen_dd;
  double mean_DD, mean_Dd, mean_dD, mean_dd;
  double SD_DD, SD_Dd, SD_dD, SD_dd;
  double theta[2];      /* theta */
  int breakFlag = FALSE;
  double gfreq = 0;     /* disease gene frequency */
  int ret;

  int status;

  double *marker1Pos, *marker2Pos;
  double *prevPos, *currPos;    /* for MP */
  double dist;
  double mkrFreq;
  double ppl;
  int i, j, k;
  int liabIdx;
  int mkrFreqIdx;
  int posIdx, pedIdx;
  TraitLocus *pTraitLocus = NULL;
  int prevTraitInd;
  int dprimeIdx;

  /* only for multipoint - we don't handle LD under multipoint yet */
  if (modelType->type == MP) {
    /* allocate space to save temporary results */
    CALCHOKE (markerNameList, sizeof (char *), (size_t) modelType->numMarkers, char **);
    if (modelType->trait == DT) {

      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        /* first dimension is gene freq */
        CALCHOKE (pPedigree->traitLikelihoodDT, sizeof (double *), (size_t) modelRange->ngfreq, double **);
        CALCHOKE (pPedigree->alternativeLikelihoodDT, sizeof (double *), (size_t) modelRange->ngfreq, double **);
        for (gfreqInd = 0; gfreqInd < modelRange->ngfreq; gfreqInd++) {
          /* second dimension is penetrance */
          CALCHOKE (pPedigree->traitLikelihoodDT[gfreqInd], sizeof (double), (size_t) modelRange->npenet, double *);
          CALCHOKE (pPedigree->alternativeLikelihoodDT[gfreqInd], sizeof (double), (size_t) modelRange->npenet, double *);
        }
      }
    } else {    /* QT */

      /* first dimension is pedigree */
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        /* second dimension is gene freq */
        CALCHOKE (pPedigree->traitLikelihoodQT, sizeof (double ***), (size_t) modelRange->ngfreq, double ****);
        CALCHOKE (pPedigree->alternativeLikelihoodQT, sizeof (double ***), (size_t) modelRange->ngfreq, double ****);
        for (gfreqInd = 0; gfreqInd < modelRange->ngfreq; gfreqInd++) {

          /* third dimension is mean */
          CALCHOKE (pPedigree->traitLikelihoodQT[gfreqInd], sizeof (double **), (size_t) modelRange->npenet, double ***);
          CALCHOKE (pPedigree->alternativeLikelihoodQT[gfreqInd], sizeof (double **), (size_t) modelRange->npenet, double ***);
          for (penIdx = 0; penIdx < modelRange->npenet; penIdx++) {
            /* fourth dimension is SD */
            CALCHOKE (pPedigree->traitLikelihoodQT[gfreqInd][penIdx], sizeof (double *), (size_t) modelRange->nparam, double **);
            CALCHOKE (pPedigree->alternativeLikelihoodQT[gfreqInd][penIdx], sizeof (double *), (size_t) modelRange->nparam, double **);
            for (paramIdx = 0; paramIdx < modelRange->nparam; paramIdx++) {
              /* 5th dimension is threshold */
              CALCHOKE (pPedigree->traitLikelihoodQT[gfreqInd][penIdx][paramIdx], sizeof (double), (size_t) modelRange->ntthresh, double *);
              CALCHOKE (pPedigree->alternativeLikelihoodQT[gfreqInd][penIdx][paramIdx], sizeof (double), (size_t) modelRange->ntthresh, double *);
            }   /* paramIdx */
          }     /* penIdx */
        }       /* gfreqInd */
      } /* pedIdx */

    }
  }

  if (fpIR != NULL) {
    memset (&dk_curModel, 0, sizeof (st_DKMaxModel));
    if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
      /* Assumes that dkelvin can only handle a single D' */
      CALCHOKE (dk_curModel.dprime, (size_t) 1, sizeof (double), double *);
    }
    CALCHOKE (dk_curModel.pen, (size_t) modelRange->nlclass, sizeof (st_DKMaxModelPenVector), void *);
    writeSurfaceFileHeader ();
  }

  /* find out the max we need to allocate */
  /* after genotype lists have been built, we want to pre-allocate parental pair work space
   * it used to be done dynamically, but it's too costly 
   * parental pair is based on each nuclear family */
  stat_parental_pair_workspace (&pedigreeSet);

  /* after genotype elimination and tracking the max work space needed 
   * for constructing parental pair */
  allocate_parental_pair_workspace (&parentalPairSpace, modelType->numMarkers + 1);

  /* conditional likelihood storage space for each individual */
  allocate_likelihood_space (&pedigreeSet, modelType->numMarkers + 1);

  if (modelOptions->markerAnalysis == FALSE || originalLocusList.ppLocusList[0]->locusType != LOCUS_TYPE_MARKER) {
    /* Assume the trait locus is the first one in the list */
    traitLocus = 0;
    pLocus = originalLocusList.ppLocusList[traitLocus];
    pTraitLocus = originalLocusList.ppLocusList[traitLocus]->pTraitLocus;
    pTrait = pTraitLocus->pTraits[traitLocus];
  }

  if (modelType->type == TP) {

    /* Two point. */
    if (originalLocusList.pLDLoci == NULL)
      CALCHOKE (originalLocusList.pLDLoci, (size_t) 1, sizeof (LDLoci), LDLoci *);
    pLDLoci = &originalLocusList.pLDLoci[0];
    originalLocusList.numLDLoci = 1;

    if (modelOptions->equilibrium == LINKAGE_EQUILIBRIUM) {
      /* fake some LD information to simplify looping */
      pLDLoci->numAllele1 = 2;
      MALCHOKE (pLDLoci->ppDPrime, sizeof (double *), double **);
      MALCHOKE (pLDLoci->ppDPrime[0], sizeof (double), double *);
      MALCHOKE (pLDLoci->ppDValue, sizeof (double *), double **);
      MALCHOKE (pLDLoci->ppDValue[0], sizeof (double), double *);
      MALCHOKE (pLDLoci->ppHaploFreq, sizeof (double *) * 2, double **);
      MALCHOKE (pLDLoci->ppHaploFreq[0], sizeof (double) * 2, double *);
      MALCHOKE (pLDLoci->ppHaploFreq[1], sizeof (double) * 2, double *);

      /* Initialize it */
      pLDLoci->ppDPrime[0][0] = 0;
    }

    analysisLocusList = &savedLocusList;
    savedLocusList.numLocus = 2;
    MALCHOKE (savedLocusList.pLocusIndex, sizeof (int) * savedLocusList.numLocus, int *);
    for (i = 0; i < 3; i++) {
      MALCHOKE (savedLocusList.pPrevLocusDistance[i], sizeof (double) * savedLocusList.numLocus, double *);
      MALCHOKE (savedLocusList.pNextLocusDistance[i], sizeof (double) * savedLocusList.numLocus, double *);
      savedLocusList.pPrevLocusDistance[i][0] = -1;
      savedLocusList.pNextLocusDistance[i][1] = -1;
    }

    if (modelOptions->polynomial == TRUE) {
      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);
      holdAllPolys ();
    }

    if (modelOptions->markerAnalysis == FALSE) {
      savedLocusList.traitLocusIndex = 0;
      savedLocusList.traitOrigLocus = 0;
    } else {
      savedLocusList.traitLocusIndex = -1;
      savedLocusList.traitOrigLocus = -1;
    }

    for (loc1 = 0; loc1 < originalLocusList.numLocus - 1; loc1++) {
      savedLocusList.pLocusIndex[0] = loc1;
      pLocus1 = originalLocusList.ppLocusList[loc1];
      if (modelOptions->markerAnalysis != FALSE && pLocus1->locusType != LOCUS_TYPE_MARKER)
        continue;
      if ((pLocus1->numAllele <= 1) || ((pLocus1->numAllele == 2) && ((pLocus1->pAlleleFrequency[0] <= ERROR_MARGIN) || (pLocus1->pAlleleFrequency[1] <= ERROR_MARGIN)))) {
        WARNING ("Biallelic marker %s has a minor allele frequency less than %g, skipping!", pLocus1->sName, ERROR_MARGIN);
        continue;
      }

      for (loc2 = loc1 + 1; loc2 < originalLocusList.numLocus; loc2++) {
        if (fpIR != NULL) {
          dk_curModel.posIdx = loc2;
        }
        pLocus2 = originalLocusList.ppLocusList[loc2];
        if (pLocus2->locusType != LOCUS_TYPE_MARKER)
          continue;
        if ((pLocus2->numAllele <= 1) || ((pLocus2->numAllele == 2) && ((pLocus2->pAlleleFrequency[0] <= ERROR_MARGIN) || (pLocus2->pAlleleFrequency[1] <= ERROR_MARGIN)))) {
          WARNING ("Biallelic marker %s has a minor allele frequency less than %g, skipping!", pLocus2->sName, ERROR_MARGIN);
          continue;
        }

        savedLocusList.pLocusIndex[1] = loc2;
        initialize_max_scale ();

        if (modelOptions->markerAnalysis == MM)
          SUBSTEP ((loc2 - 1) * 100 / (originalLocusList.numLocus - 1), "Starting w/loci %s(%d alleles) and %s(%d alleles", pLocus1->sName, pLocus1->numOriginalAllele, pLocus2->sName, pLocus2->numOriginalAllele);
        else
          SUBSTEP ((loc2 - 1) * 100 / (originalLocusList.numLocus - 1),
              "Starting w/loci %s(%d alleles) and %s(%d alleles) (%d of %d pairs)", pLocus1->sName, pLocus1->numOriginalAllele, pLocus2->sName, pLocus2->numOriginalAllele, loc2, originalLocusList.numLocus - 1);

        /* Find out number of alleles this marker locus has */
        if (modelOptions->equilibrium == LINKAGE_DISEQUILIBRIUM) {
          /* get the LD parameters */
          pLambdaCell = findLambdas (modelRange, pLocus1->numOriginalAllele, pLocus2->numOriginalAllele);
          reallocate_LD_loci (pLDLoci, pLocus1->numOriginalAllele, pLocus2->numOriginalAllele);

          // Create these variables ahead of likelihood polynomial build in hopes of preventing in-build creation.

          if (modelOptions->polynomial == TRUE) {
            char vName[128];
            int a0, a1;
            for (a0 = 0; a0 < pLocus1->numOriginalAllele; a0++) {
              for (a1 = 0; a1 < pLocus2->numOriginalAllele; a1++) {
                sprintf (vName, "ppHaploFreq_lA%d_rA%d", a0, a1);
                variableExp (&pLDLoci->ppHaploFreq[a0][a1], NULL, 'D', vName);
              }
            }
          }

          pLDLoci->locus1 = loc1;
          pLDLoci->locus2 = loc2;
          pLDLoci->numAllele1 = pLocus1->numOriginalAllele;
          pLDLoci->numAllele2 = pLocus2->numOriginalAllele;
          if (pLocus1->numOriginalAllele == 2 && pLocus2->numOriginalAllele == 2)
            R_square_flag = TRUE;
          else
            R_square_flag = FALSE;
        }
        // Build the polynomial now to avoid interfering with evaluation later
        if (modelOptions->polynomial == TRUE) {
          if (modelType->trait == DICHOTOMOUS) {
            sprintf (partialPolynomialFunctionName, "TD_LC%d_C%d_P%%s_%s_%s", modelRange->nlclass, pLocus2->pMapUnit->chromosome, pLocus1->sName, pLocus2->sName);
            swPushPhase ('k', "buildTD");
          } else {
            sprintf (partialPolynomialFunctionName, "TQ_LC%d_C%d_P%%s_%s_%s", modelRange->nlclass, pLocus2->pMapUnit->chromosome, pLocus1->sName, pLocus2->sName);
            swPushPhase ('k', "buildTQ");
          }
          swStart (combinedBuildSW);
	  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
	    Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
	    ret = build_likelihood_polynomial (pPedigree);
	  }
          swPopPhase ('k');
          swStop (combinedBuildSW);
        }
        swPushPhase ('k', "eval");      // Now we're full-tilt evaluation
        swStart (combinedComputeSW);

        loopMarkerFreqFlag = 0;
        if (modelRange->nafreq >= 2 && modelOptions->equilibrium == LINKAGE_DISEQUILIBRIUM && pLocus2->numOriginalAllele == 2) {
          loopMarkerFreqFlag = 1;
        } else if (modelRange->nafreq == 0) {
          /* add a fake one to facilitate loops and other handlings */
          addAlleleFreq (modelRange, pLocus2->pAlleleFrequency[0]);
        } else {
          modelRange->nafreq = 1;
          modelRange->afreq[0] = pLocus2->pAlleleFrequency[0];
        }

        /* allocate/initialize result storage */
        initialize_tp_result_storage ();

        /* we will force marker allele frequency loop to execute at least once */
        for (mkrFreqIdx = 0; mkrFreqIdx == 0 || mkrFreqIdx < modelRange->nafreq; mkrFreqIdx++) {
          mkrFreq = pLocus2->pAlleleFrequency[0];
          paramSet.mkrFreqIdx = mkrFreq;
          /* we should only loop over marker allele frequency under twopoint
           * and when markers are SNPs (only have two alleles) */
          if (loopMarkerFreqFlag) {
            mkrFreq = modelRange->afreq[mkrFreqIdx];
            /* update the locus */
            pLocus2->pAlleleFrequency[0] = mkrFreq;
            pLocus2->pAlleleFrequency[1] = 1 - mkrFreq;
            if (modelOptions->polynomial == TRUE);
            else
              update_locus (&pedigreeSet, loc2);
          }
          /* Loop over the penetrances, genefrequencies, thetas and call
           * the likelihood calculation, storing each value obtained to
           * disk. */

          for (gfreqInd = 0; (gfreqInd == 0 && modelOptions->markerAnalysis != FALSE) || gfreqInd < modelRange->ngfreq; gfreqInd++) {
            paramSet.gfreqIdx = gfreqInd;
            /* Here's a little bomb that should highlight if paramSet.gfreq is used
             * without being properly set.
             */
            paramSet.gfreq = -1;
            if (modelOptions->markerAnalysis == FALSE) {
              gfreq = modelRange->gfreq[gfreqInd];
              paramSet.gfreq = gfreq;

              if (fpIR != NULL)
                dk_curModel.dgf = gfreq;

              pLocus->pAlleleFrequency[0] = gfreq;
              pLocus->pAlleleFrequency[1] = 1 - gfreq;
              if (modelOptions->polynomial == TRUE);
              else
                update_locus (&pedigreeSet, loc1);
            }

            /* Clear Dprime combination impossible flag */
            memset (pLambdaCell->impossibleFlag, 0, sizeof (int) * pLambdaCell->ndprime);
            /* set up haplotype frequencies */
            dprime0Idx = -1;
	    if (modelOptions->equilibrium == LINKAGE_DISEQUILIBRIUM) {
	      for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
		if (isDPrime0 (pLambdaCell->lambda[dprimeIdx], pLambdaCell->m, pLambdaCell->n))
		  dprime0Idx = dprimeIdx;
		status = setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprimeIdx);
		if (status < 0)
		  pLambdaCell->impossibleFlag[dprimeIdx] = 1;
	      }
	      ASSERT ((modelOptions->equilibrium != LINKAGE_DISEQUILIBRIUM) || (dprime0Idx != -1), "The requisite zero D' was not found");
	    } else
	      dprime0Idx = 0;

            if (modelType->trait == DICHOTOMOUS) {

              for (penIdx = 0; (penIdx == 0 && modelOptions->markerAnalysis != FALSE) || penIdx < modelRange->npenet; penIdx++) {
                paramSet.penIdx = penIdx;
                if (modelOptions->markerAnalysis == FALSE && pLocus1->locusType == LOCUS_TYPE_TRAIT) {
                  for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++) {
                    pen_DD = modelRange->penet[liabIdx][0][penIdx];
                    pen_Dd = modelRange->penet[liabIdx][1][penIdx];
                    pen_dD = modelRange->penet[liabIdx][2][penIdx];
                    pen_dd = modelRange->penet[liabIdx][3][penIdx];
                    pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
                    pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
                    pTrait->penetrance[2][liabIdx][1][0] = pen_dD;
                    pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
                    pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
                    pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
                    pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_dD;
                    pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;

                    if (fpIR != NULL) {
                      dk_curModel.pen[liabIdx].DD = pen_DD;
                      dk_curModel.pen[liabIdx].Dd = pen_Dd;
                      dk_curModel.pen[liabIdx].dD = pen_dD;
                      dk_curModel.pen[liabIdx].dd = pen_dd;
                    }
                  }
                  if (modelOptions->polynomial == TRUE);
                  else
                    update_penetrance (&pedigreeSet, traitLocus);
                }
                /* get the likelihood at 0.5 first and LD=0 */
                if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
                  set_null_dprime (pLDLoci);
                  copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprime0Idx]);
                  copy_DValue (pLDLoci, pLambdaCell->DValue[dprime0Idx]);
                  ASSERT (pLambdaCell->impossibleFlag[dprime0Idx] == 0, "Haplotype frequency combination impossible at LE");
                }
                for (k = 0; k < 3; k++) {
                  analysisLocusList->pNextLocusDistance[k][0] = 0.5;
                  analysisLocusList->pPrevLocusDistance[k][1] = 0.5;
                }

                if (modelOptions->polynomial != TRUE)
                  status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);
                ret = compute_likelihood (&pedigreeSet);
                cL[0]++;        // TP DT NULL hypothesis

                if (swProgressRequestFlag) {
                  swProgressRequestFlag = FALSE;
                  DETAIL (0, "Likelihood calculation for NULL hypothesis %lu%% complete (~%lu min left)",
                      (cL[0] + cL[1]) * 100 / (eCL[0] + eCL[1]),
                      ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) * (eCL[0] + eCL[1]) / (cL[0] + cL[1]) - (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
                }
                if (ret == -2)
                  ERROR ("Negative likelihood for theta 0.5");

                if (ret == -1) {
                  DIAG (OVERALL, 1, {
                        fprintf (stderr, "dgf=%f\n", gfreq);
			for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++) {
                          pen_DD = modelRange->penet[liabIdx][0][penIdx];
			  pen_Dd = modelRange->penet[liabIdx][1][penIdx];
			  pen_dD = modelRange->penet[liabIdx][2][penIdx];
			  pen_dd = modelRange->penet[liabIdx][3][penIdx];
			  if (modelOptions->imprintingFlag)
			    fprintf (stderr, "Liab %d penentrance %f %f %f %f\n",
				     liabIdx + 1, pen_DD, pen_Dd, pen_dD, pen_dd);
                          else
			    fprintf (stderr, "Liab %d penentrance %f %f %f\n",
				     liabIdx + 1, pen_DD, pen_Dd, pen_dd);
			}
                      }
                  );
                  ERROR ("Zero likelihood for theta 0.5");
                }

                /* Save the results for NULL */
                for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                  /* save the likelihood at null */
                  Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                  pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
                }

                log10_likelihood_null = pedigreeSet.log10Likelihood;
                for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
                  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {

                    if (fpIR != NULL)
                      dk_curModel.dprime[0] = pLambdaCell->lambda[dprimeIdx][0][0];

                    copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);
                    if (pLambdaCell->impossibleFlag[dprimeIdx] != 0) {
                      // If we're going to bail at this point, add the progress count loop factor
                      cL[1] += modelRange->ntheta;
                      continue;
                    }
                    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
                    copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);
                    /* calculate R square if the marker is a SNP */
                    if (R_square_flag == TRUE)
                      R_square = calculate_R_square (pLocus1->pAlleleFrequency[0], pLocus2->pAlleleFrequency[0], pLDLoci->ppDValue[0][0]);
                    else
                      R_square = -1;
                    paramSet.R_square = R_square;
                  }
                  paramSet.dprimeIdx = dprimeIdx;
                  for (thetaInd = 0; thetaInd < modelRange->ntheta; thetaInd++) {
                    paramSet.thetaIdx = thetaInd;
                    if (modelOptions->mapFlag == SA) {
                      theta[0] = modelRange->theta[0][thetaInd];
                      theta[1] = modelRange->theta[1][thetaInd];
                      for (k = 0; k < 3; k++) {
                        analysisLocusList->pNextLocusDistance[k][0] = theta[0];
                        analysisLocusList->pPrevLocusDistance[k][1] = theta[0];
                      }
                    } else {
                      analysisLocusList->pNextLocusDistance[MAP_POS_MALE][0] = analysisLocusList->pPrevLocusDistance[MAP_POS_MALE][1] = modelRange->theta[0][thetaInd];
                      analysisLocusList->pNextLocusDistance[MAP_POS_FEMALE][0] = analysisLocusList->pPrevLocusDistance[MAP_POS_FEMALE][1] = modelRange->theta[1][thetaInd];
                    }
                    if (fpIR != NULL) {
                      dk_curModel.theta[0] = modelRange->theta[0][thetaInd];
                      dk_curModel.theta[1] = modelRange->theta[1][thetaInd];
                    }

                    if (modelOptions->polynomial == TRUE);
                    else
                      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);

                    // No new name for a polynomial here because we're reusing the existing one
                    ret = compute_likelihood (&pedigreeSet);
                    cL[1]++;    // TP DT alternative hypothesis (combined likelihood)
                    if (swProgressRequestFlag) {
                      swProgressRequestFlag = FALSE;
                      DETAIL (0, "Combined likelihood calculation %lu%% comp1ete (~%lu min left [%lu %lu])",
			      (cL[0] + cL[1]) * 100 / (eCL[0] + eCL[1]),
			      ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
			       (eCL[0] + eCL[1]) / (cL[0] + cL[1]) - (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60,
			      combinedComputeSW->swAccumWallTime, combinedBuildSW->swAccumWallTime);
                    }
                    record_tp_result (ret, &pedigreeSet, &paramSet, loc2);
                  }     /* end of theta loop */
                }       /* end of D prime loop */
                if (modelOptions->markerAnalysis != FALSE) {
                  /* marker to marker analysis, marker allele frequency is fixed */
                  gfreqInd = modelRange->ngfreq;
                  break;
                }
                if (modelOptions->markerAnalysis != FALSE) {
                  /* marker to marker analysis, penetrance stays at 1 */
                  break;
                }
              } /* end of penetrance loop */
            } /* end of DT */
            else {      /* start of QT or COMBINED */


              /* this should be MEAN + SD */
              for (paramIdx = 0; (paramIdx == 0 && modelType->distrib == QT_FUNCTION_CHI_SQUARE)
                  || (modelType->distrib != QT_FUNCTION_CHI_SQUARE && paramIdx < modelRange->nparam); paramIdx++) {
                paramSet.paramIdx = paramIdx;
                for (penIdx = 0; penIdx < modelRange->npenet; penIdx++) {
                  paramSet.penIdx = penIdx;
                  breakFlag = FALSE;
                  for (thresholdIdx = 0; thresholdIdx < modelRange->ntthresh; thresholdIdx++) {
                    paramSet.thresholdIdx = thresholdIdx;
                    if (modelOptions->markerAnalysis == FALSE) {
                      for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++) {
                        mean_DD = modelRange->penet[liabIdx][0][penIdx];
                        mean_Dd = modelRange->penet[liabIdx][1][penIdx];
                        mean_dD = modelRange->penet[liabIdx][2][penIdx];
                        mean_dd = modelRange->penet[liabIdx][3][penIdx];
                        SD_DD = modelRange->param[liabIdx][0][0][paramIdx];
                        SD_Dd = modelRange->param[liabIdx][1][0][paramIdx];
                        SD_dD = modelRange->param[liabIdx][2][0][paramIdx];
                        SD_dd = modelRange->param[liabIdx][3][0][paramIdx];
                        /* threshold for QT */
                        threshold = modelRange->tthresh[liabIdx][thresholdIdx];

                        if (fpIR != NULL) {
                          dk_curModel.pen[liabIdx].DD = mean_DD;
                          dk_curModel.pen[liabIdx].Dd = mean_Dd;
                          dk_curModel.pen[liabIdx].dD = mean_dD;
                          dk_curModel.pen[liabIdx].dd = mean_dd;
                          dk_curModel.pen[liabIdx].DDSD = SD_DD;
                          dk_curModel.pen[liabIdx].DdSD = SD_Dd;
                          dk_curModel.pen[liabIdx].dDSD = SD_dD;
                          dk_curModel.pen[liabIdx].ddSD = SD_dd;
                          dk_curModel.pen[liabIdx].threshold = threshold;
                        }
                        /* check against the hard coded constraint */
                        if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
                          constraint = (1 - gfreq) * (1 - gfreq) * mean_dd * SD_dd + 2 * gfreq * (1 - gfreq) * mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
                          if (constraint >= 3.0 || constraint <= -3.0) {
                            breakFlag = TRUE;
                            break;
                          }
                        }
                        pTrait->means[liabIdx][0][0] = mean_DD;
                        pTrait->means[liabIdx][0][1] = mean_Dd;
                        pTrait->means[liabIdx][1][0] = mean_dD;
                        pTrait->means[liabIdx][1][1] = mean_dd;
                        pTrait->stddev[liabIdx][0][0] = SD_DD;
                        pTrait->stddev[liabIdx][0][1] = SD_Dd;
                        pTrait->stddev[liabIdx][1][0] = SD_dD;
                        pTrait->stddev[liabIdx][1][1] = SD_dd;

                        /* threshold for QT */
                        pTrait->cutoffValue[liabIdx] = threshold;

                      } /* liability class Index */
                      if (breakFlag == TRUE)
                        continue;
                      if (modelOptions->polynomial == TRUE);
                      else
                        update_penetrance (&pedigreeSet, traitLocus);
                    }
                    /* marker to marker analysis */
                    /* get the likelihood at 0.5 first and LD=0 */
                    if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
                      set_null_dprime (pLDLoci);
                      copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprime0Idx]);
                      copy_DValue (pLDLoci, pLambdaCell->DValue[dprime0Idx]);
                      ASSERT (pLambdaCell->impossibleFlag[dprime0Idx] == 0, "Haplotype frequency combination impossible for LE\n");
                    }
                    for (k = 0; k < 3; k++) {
                      analysisLocusList->pNextLocusDistance[k][0] = 0.5;
                      analysisLocusList->pPrevLocusDistance[k][1] = 0.5;
                    }

                    if (modelOptions->polynomial != TRUE)
                      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);

                    ret = compute_likelihood (&pedigreeSet);
                    cL[2]++;    // TP QT NULL hypothesis likelihood
                    if (swProgressRequestFlag) {
                      swProgressRequestFlag = FALSE;
                      DETAIL (0, "Likelihood calculation for NULL hypothesis  %lu%% complete (~%lu min left)",
                          (cL[2] + cL[3]) * 100 / (eCL[2] + eCL[3]),
                          ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) * (eCL[2] + eCL[3]) / (cL[2] + cL[3]) - (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
                    }
                    if (ret == -2)
                      ERROR ("Negative likelihood for theta 0.5");

                    if (ret == -1) {
                      DIAG (OVERALL, 1, {
                            fprintf (stderr, "dgf=%f\n", gfreq);
                            for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++) {
                              pen_DD = modelRange->penet[liabIdx][0][penIdx];
			      pen_Dd = modelRange->penet[liabIdx][1][penIdx];
			      pen_dD = modelRange->penet[liabIdx][2][penIdx];
			      pen_dd = modelRange->penet[liabIdx][3][penIdx];
			      if (modelOptions->imprintingFlag)
				fprintf (stderr, "Liab %d penentrance %f %f %f %f\n",
					 liabIdx + 1, pen_DD, pen_Dd, pen_dD, pen_dd);
                              else
				fprintf (stderr, "Liab %d penentrance %f %f %f\n",
				     liabIdx + 1, pen_DD, pen_Dd, pen_dd);
			    }
                          }
                      );
                      ERROR ("Zero likelihood for theta of 0.5");
                    }
                    for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                      /* save the likelihood at null */
                      Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                      pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
                    }
                    log10_likelihood_null = pedigreeSet.log10Likelihood;
                    for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
                      paramSet.dprimeIdx = dprimeIdx;
                      if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {

                        if (fpIR != NULL)
                          dk_curModel.dprime[0] = pLambdaCell->lambda[dprimeIdx][0][0];

                        copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);
                        if (pLambdaCell->impossibleFlag[dprimeIdx] != 0) {
                          // If we're going to bail at this point, add the progress count loop factor
                          cL[3] += modelRange->ntheta;
                          continue;
                        }
                        copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
                        copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);
                      }
                      for (thetaInd = 0; thetaInd < modelRange->ntheta; thetaInd++) {
                        if (fpIR != NULL) {
                          dk_curModel.theta[0] = modelRange->theta[0][thetaInd];
                          dk_curModel.theta[1] = modelRange->theta[1][thetaInd];
                        }

                        paramSet.thetaIdx = thetaInd;
                        if (modelOptions->mapFlag == SA) {
                          theta[0] = modelRange->theta[0][thetaInd];
                          theta[1] = modelRange->theta[1][thetaInd];
                          for (k = 0; k < 3; k++) {
                            analysisLocusList->pNextLocusDistance[k][0] = theta[0];
                            analysisLocusList->pPrevLocusDistance[k][1] = theta[0];
                          }
                        } else {
                          analysisLocusList->pNextLocusDistance[MAP_POS_MALE][0] = analysisLocusList->pPrevLocusDistance[MAP_POS_MALE][1] = modelRange->theta[0][thetaInd];
                          analysisLocusList->pNextLocusDistance[MAP_POS_FEMALE][0] = analysisLocusList->pPrevLocusDistance[MAP_POS_FEMALE][1] = modelRange->theta[1][thetaInd];
                        }

                        if (modelOptions->polynomial != TRUE);
                        status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);

                        // No new name for a polynomial here because we're reusing the existing one
                        ret = compute_likelihood (&pedigreeSet);
                        cL[3]++;        // TP QT/CT alternative hypothesis
                        if (swProgressRequestFlag) {
                          swProgressRequestFlag = FALSE;
                          DETAIL (0, "Combined likelihood calculation %lu%% complete (~%lu min left)",
                              (cL[2] + cL[3]) * 100 / (eCL[2] + eCL[3]),
                              ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) * (eCL[2] + eCL[3]) / (cL[2] + cL[3]) - (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
                        }
                        record_tp_result (ret, &pedigreeSet, &paramSet, loc2);
                      } /* end of theta */
                    }   /* end of D prime */
                    if (modelOptions->markerAnalysis != FALSE)
                      break;
                  }     /* end of threshold loop */
                  if (modelOptions->markerAnalysis != FALSE)
                    break;
                }       /* end of penetrance loop */
                if (modelOptions->markerAnalysis != FALSE)
                  break;
              } /* end of parameter loop */
              if (modelOptions->markerAnalysis != FALSE)
                break;
            }   /* end of QT */
          }     /* end of gene freq */
          /* only loop marker allele frequencies when doing LD */
          if (modelOptions->equilibrium == LINKAGE_EQUILIBRIUM)
            break;
          /* we can only do SNPs when looping over marker allele frequency */
          if (pLocus2->numOriginalAllele > 2)
            break;
        }       /* end of marker allele frequency looping */

        get_average_LR (tp_result);     // Calculate the average BR per (D', theta) pair
        rescale_tp_result_dprime0 (dprime0Idx);
        rescale_tp_result (-1); // Rescale across (D', theta) pairs


        write2ptBRFile (loc1, loc2);
        write2ptMODFile (loc1, loc2, dprime0Idx);
        writePPLFileDetail (dprime0Idx);

        if (modelOptions->polynomial == TRUE)
          pedigreeSetPolynomialClearance (&pedigreeSet);

        if (modelOptions->markerAnalysis == ADJACENTMARKER)
          loc2 = originalLocusList.numLocus;

        free_tp_result_storage ();      // Free two point result storage

        swPopPhase ('k');
        swStop (combinedComputeSW);

      } /* end of looping second locus - loc2 */
      /* if we are doing trait marker, then we are done */
      /* Used to read: modelOptions->markerToMarker != TRUE which
       * is the same as markerAnalysis == FALSE as long as the old
       * markerToMarker and adjacentMarker flags were truly
       * orthogonal. Otherwise, it should be markerAnalysis !=
       * ADJACENTMARKER. */
      if (modelOptions->markerAnalysis == FALSE)
        loc1 = originalLocusList.numLocus;
    }   /* end of looping first locus - loc1 */

    SUBSTEP (0, "Finished analysis w/build time of %d, evaluate time of %d", combinedBuildSW->swAccumWallTime, combinedComputeSW->swAccumWallTime);

  } /* end of two point */
  else {        /* multipoint */

    int markerSetChanged;       /* Flag for multipoint analysis, did set of markers change? */
    int locusListChanged;       /* flag for multipoint analysis, did relative trait position or marker set change? */

    int prevFirstMarker;        /* first marker in the set for multipoint analysis */
    int prevLastMarker; /* last marker in the set for multipoint analysis */

    double relativePos;
    double traitPos;    /* trait position for multipoint analysis */

    /* marker set locus list for each position */
    markerLocusList.maxNumLocus = modelType->numMarkers;
    markerLocusList.numLocus = modelType->numMarkers;
    markerLocusList.traitOrigLocus = -1;
    markerLocusList.traitLocusIndex = -1;
    CALCHOKE (markerLocusList.pLocusIndex, (size_t) markerLocusList.maxNumLocus, sizeof (int), int *);
    for (k = 0; k < 3; k++) {
      CALCHOKE (markerLocusList.pPrevLocusDistance[k], (size_t) markerLocusList.maxNumLocus, sizeof (double), double *);
      CALCHOKE (markerLocusList.pNextLocusDistance[k], (size_t) markerLocusList.maxNumLocus, sizeof (double), double *);
    }

    /* assuming we always have trait in the analysis - this may not be true 
     * need to add code to process marker to marker analysis under multipoin
     */
    savedLocusList.numLocus = modelType->numMarkers + 1;
    savedLocusList.maxNumLocus = modelType->numMarkers + 1;
    CALCHOKE (savedLocusList.pLocusIndex, (size_t) savedLocusList.maxNumLocus, sizeof (int), int *);
    for (k = 0; k < 3; k++) {
      CALCHOKE (savedLocusList.pPrevLocusDistance[k], (size_t) savedLocusList.maxNumLocus, sizeof (double), double *);
      CALCHOKE (savedLocusList.pNextLocusDistance[k], (size_t) savedLocusList.maxNumLocus, sizeof (double), double *);
    }

    /* Allocate storage to calculate the trait likelihood independent of the trait position */
    traitLocusList.numLocus = 1;
    traitLocusList.maxNumLocus = 1;
    traitLocusList.traitLocusIndex = 0;
    traitLocusList.traitOrigLocus = traitLocus;
    CALCHOKE (traitLocusList.pLocusIndex, (size_t) traitLocusList.maxNumLocus, sizeof (int), int *);
    traitLocusList.pLocusIndex[0] = 0;
    for (k = 0; k < 3; k++) {
      CALCHOKE (traitLocusList.pPrevLocusDistance[k], (size_t) savedLocusList.maxNumLocus, sizeof (double), double *);
      CALCHOKE (traitLocusList.pNextLocusDistance[k], (size_t) savedLocusList.maxNumLocus, sizeof (double), double *);

      traitLocusList.pPrevLocusDistance[k][0] = -1;
      traitLocusList.pNextLocusDistance[k][0] = -1;
    }
    /* populate the trait xmission matrix */
    analysisLocusList = &traitLocusList;
    xmissionMatrix = traitMatrix;
    status = populate_xmission_matrix (traitMatrix, 1, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);
    if (modelOptions->polynomial == TRUE)
      holdAllPolys ();

    /* For trait likelihood */
    SUBSTEP (0, "Calculating trait likelihood");

    analysisLocusList = &traitLocusList;
    xmissionMatrix = traitMatrix;
    if (pTrait->type == DICHOTOMOUS) {
      /* load all saved trait likelihood */
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        if (modelOptions->saveResults == TRUE) {
          pPedigree->load_flag = restoreTrait (modelOptions->sexLinked, pPedigree->sPedigreeID, pPedigree->traitLikelihoodDT);
        } else
          pPedigree->load_flag = 0;
      }

      for (penIdx = 0; (penIdx == 0) || (modelOptions->dryRun == 0 && penIdx < modelRange->npenet); penIdx++) {
        for (liabIdx = 0; (liabIdx == 0) || (modelOptions->dryRun == 0 && liabIdx < modelRange->nlclass); liabIdx++) {
          pen_DD = modelRange->penet[liabIdx][0][penIdx];
          pen_Dd = modelRange->penet[liabIdx][1][penIdx];
          pen_dD = modelRange->penet[liabIdx][2][penIdx];
          pen_dd = modelRange->penet[liabIdx][3][penIdx];
          pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
          pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
          pTrait->penetrance[2][liabIdx][1][0] = pen_dD;
          pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
          pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
          pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
          pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_dD;
          pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;
        }

        if (modelOptions->polynomial == TRUE);
        else
          /* only need to update trait locus */
          update_penetrance (&pedigreeSet, traitLocus);

        /* Iterate over gene frequencies, but only one time thru if doing a dry-run. */
        for (gfreqInd = 0; (gfreqInd == 0) || (modelOptions->dryRun == 0 && gfreqInd < modelRange->ngfreq); gfreqInd++) {

          /* updated trait locus allele frequencies */
          gfreq = modelRange->gfreq[gfreqInd];
          pLocus->pAlleleFrequency[0] = gfreq;
          pLocus->pAlleleFrequency[1] = 1 - gfreq;

          if (modelOptions->polynomial == TRUE)
            sprintf (partialPolynomialFunctionName, "MDT_LC%d_C%d_P%%sSL%d", modelRange->nlclass, (originalLocusList.ppLocusList[1])->pMapUnit->chromosome, modelOptions->sexLinked);
          else
            update_locus (&pedigreeSet, traitLocus);

          /* Compute the likelihood for the trait */
          ret = compute_likelihood (&pedigreeSet);
          cL[4]++;      // MP DT trait likelihood
          if (swProgressRequestFlag) {
            swProgressRequestFlag = FALSE;
            DETAIL (0, "Trait likelihood calculations %lu%% complete", cL[4] * 100 / eCL[4]);
          }

          if (modelOptions->dryRun != 0)
            continue;

          if (ret == -2)
            ERROR ("Negative likelihood for trait");

          if (ret == -1) {
            DIAG (OVERALL, 1, {
                  fprintf (stderr, "dgf=%f\n", gfreq);
                  for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++) {
                    pen_DD = modelRange->penet[liabIdx][0][penIdx];
		    pen_Dd = modelRange->penet[liabIdx][1][penIdx];
		    pen_dD = modelRange->penet[liabIdx][2][penIdx];
		    pen_dd = modelRange->penet[liabIdx][3][penIdx];
		    if (modelOptions->imprintingFlag)
		      fprintf (stderr, "Liab %d penentrance %f %f %f %f\n",
			       liabIdx + 1, pen_DD, pen_Dd, pen_dD, pen_dd);
                    else
		      fprintf (stderr, "Liab %d penentrance %f %f %f\n",
			       liabIdx + 1, pen_DD, pen_Dd, pen_dd);
		  }
                }
            );
            ERROR ("Zero likelihood for trait");
          }
          /* save the results for NULL */
          for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
            /* save the likelihood at null */
            Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
            if (pPedigree->load_flag == 0) {    /*update only for the pedigrees which were add for this run */
              pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
              pPedigree->traitLikelihoodDT[gfreqInd][penIdx] = pPedigree->likelihood;
            }
          }

          log10_likelihood_null = pedigreeSet.log10Likelihood;
#if 0
          likelihoodDT[gfreqInd][penIdx] = log10_likelihood_null;
#endif
        }       /* gfreq */
      } /* pen */

      /* save all  trait likelihood which were created in this run */
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        /* save the likelihood at null */
        Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        if ((modelOptions->saveResults == TRUE) && (pPedigree->load_flag == 0)) {       /*save only for the pedigrees which were add for this run */
          pPedigree->load_flag = saveTrait (modelOptions->sexLinked, pPedigree->sPedigreeID, pPedigree->traitLikelihoodDT);
        } else {
          pPedigree->load_flag = 0;
        }
      }
    } else
      /* multipoint QT or COMBINED */
    {
      for (gfreqInd = 0; gfreqInd < modelRange->ngfreq; gfreqInd++) {
        gfreq = modelRange->gfreq[gfreqInd];
        pLocus->pAlleleFrequency[0] = gfreq;
        pLocus->pAlleleFrequency[1] = 1 - gfreq;

        update_locus (&pedigreeSet, traitLocus);
        /* this should be MEAN + SD */
        for (paramIdx = 0; paramIdx < modelRange->nparam; paramIdx++) {
          for (penIdx = 0; penIdx < modelRange->npenet; penIdx++) {
            breakFlag = FALSE;
            for (thresholdIdx = 0; thresholdIdx < modelRange->ntthresh; thresholdIdx++) {
              for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++) {
                mean_DD = modelRange->penet[liabIdx][0][penIdx];
                mean_Dd = modelRange->penet[liabIdx][1][penIdx];
                mean_dD = modelRange->penet[liabIdx][2][penIdx];
                mean_dd = modelRange->penet[liabIdx][3][penIdx];
                SD_DD = modelRange->param[liabIdx][0][0][paramIdx];
                SD_Dd = modelRange->param[liabIdx][1][0][paramIdx];
                SD_dD = modelRange->param[liabIdx][2][0][paramIdx];
                SD_dd = modelRange->param[liabIdx][3][0][paramIdx];
                threshold = modelRange->tthresh[liabIdx][thresholdIdx];

                /* check against the hard coded constraint */
                if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
                  constraint = (1 - gfreq) * (1 - gfreq) * mean_dd * SD_dd + 2 * gfreq * (1 - gfreq) * mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
                  DIAG (OVERALL, 1, {
                        fprintf (stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
				 constraint, gfreq, mean_DD, SD_DD, mean_Dd, SD_DD, mean_dd, SD_dd);
                      }
                  );
                  if (constraint >= 3.0 || constraint <= -3.0) {
                    breakFlag = TRUE;
                    break;
                  }
                }
                pTrait->means[liabIdx][0][0] = mean_DD;
                pTrait->means[liabIdx][0][1] = mean_Dd;
                pTrait->means[liabIdx][1][0] = mean_dD;
                pTrait->means[liabIdx][1][1] = mean_dd;
                pTrait->stddev[liabIdx][0][0] = SD_DD;
                pTrait->stddev[liabIdx][0][1] = SD_Dd;
                pTrait->stddev[liabIdx][1][0] = SD_dD;
                pTrait->stddev[liabIdx][1][1] = SD_dd;

                /* threshold for QT */
                pTrait->cutoffValue[liabIdx] = threshold;

              } /* liability class Index */
              if (breakFlag == TRUE)
                continue;
              if (modelOptions->polynomial == TRUE)
                sprintf (partialPolynomialFunctionName, "MQT_LC%d_C%d_P%%sSL%d", modelRange->nlclass, (originalLocusList.ppLocusList[1])->pMapUnit->chromosome, modelOptions->sexLinked);
              else
                update_penetrance (&pedigreeSet, traitLocus);
              ret = compute_likelihood (&pedigreeSet);
              cL[5]++;  // MP QT/CT trait likelihood
              if (swProgressRequestFlag) {
                swProgressRequestFlag = FALSE;
                DETAIL (0, "Trait likelihood calculations %lu%% complete", cL[5] * 100 / eCL[5]);
              }

              if (ret == -2)
                ERROR ("Negative likelihood for trait");

              if (ret == -1) {
                DIAG (OVERALL, 1, {
                      fprintf (stderr, "dgf=%f\n", gfreq);
                      for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++) {
                        pen_DD = modelRange->penet[liabIdx][0][penIdx];
			pen_Dd = modelRange->penet[liabIdx][1][penIdx];
			pen_dD = modelRange->penet[liabIdx][2][penIdx];
			pen_dd = modelRange->penet[liabIdx][3][penIdx];
			if (modelOptions->imprintingFlag)
			  fprintf (stderr, "Liab %d penentrance %f %f %f %f\n",
				   liabIdx + 1, pen_DD, pen_Dd, pen_dD, pen_dd);
                        else
			  fprintf (stderr, "Liab %d penentrance %f %f %f\n",
				   liabIdx + 1, pen_DD, pen_Dd, pen_dd);
		      }
                    }
                );
                ERROR ("Zero likelihood for trait");
              }

              for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                /* save the likelihood at null */
                Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
                pPedigree->traitLikelihoodQT[gfreqInd][penIdx]
                    [paramIdx][thresholdIdx] = pPedigree->likelihood;
              }

              log10_likelihood_null = pedigreeSet.log10Likelihood;
              if (isnan (log10_likelihood_null))
                ERROR ("Trait likelihood is NAN");
            }   /* thresholdIdx */
          }     /* penIdx */
        }       /* paramIdx */
      } /* gfreq */
    }   /* end of QT */

    SUBSTEP (0, "Trait likelihood calculations complete");

    if (modelOptions->polynomial == TRUE) {
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        /* save the likelihood at trait */
        Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        pPedigree->traitLikelihoodPolynomial = pPedigree->likelihoodPolynomial;
        pPedigree->traitLikelihoodPolyList = pPedigree->likelihoodPolyList;
      }
    }

    /* print out some statistics under dry run */
    if (modelOptions->dryRun != 0) {
      print_dryrun_stat (&pedigreeSet, -1);
    }

    /* get the trait locations we need to evaluate at */
    numPositions = modelRange->ntloc;
    CALCHOKE (mp_result, (size_t) numPositions, sizeof (SUMMARY_STAT), SUMMARY_STAT *);

    writeMPBRFileHeader ();
    writeMPMODFileHeader ();

    prevFirstMarker = -1;
    prevLastMarker = -1;
    prevTraitInd = -1;
    leftMarker = -1;

    /* Iterate over all positions in the analysis. */
    for (posIdx = 0; posIdx < numPositions; posIdx++) {
      if (fpIR != NULL) {
        dk_curModel.posIdx = posIdx;
      }

      mp_result[posIdx].max_penIdx = -1;
      /* positions listed are sex average positions */
      traitPos = modelRange->tloc[posIdx];
      /* Set the sex-averaged position first. The sex-specific positions will be updated 
       * once markers are selected since interpolation might be needed. */
      pTraitLocus->mapPosition[0] = traitPos;
      pTraitLocus->mapPosition[1] = traitPos;
      pTraitLocus->mapPosition[2] = traitPos;
      /* initialize the locusList */
      analysisLocusList = &savedLocusList;
      memset (analysisLocusList->pLocusIndex, 0, sizeof (int) * analysisLocusList->maxNumLocus);
      for (k = 0; k < 3; k++) {
        memset (&analysisLocusList->pPrevLocusDistance[k][0], 0, sizeof (double) * analysisLocusList->maxNumLocus);
        memset (&analysisLocusList->pNextLocusDistance[k][0], 0, sizeof (double) * analysisLocusList->maxNumLocus);
      }
      analysisLocusList->numLocus = 1;
      analysisLocusList->pLocusIndex[0] = traitLocus;
      for (k = 0; k < 3; k++) {
        analysisLocusList->pPrevLocusDistance[k][0] = -1;
        analysisLocusList->pNextLocusDistance[k][0] = -1;
      }
      /* select markers to be used for the multipoint analysis */
      add_markers_to_locuslist (analysisLocusList, modelType->numMarkers, &leftMarker, 0, originalLocusList.numLocus - 1, traitPos, 0);
      /* store the markers used */
      CALCHOKE (mp_result[posIdx].pMarkers, (size_t) modelType->numMarkers, sizeof (int), int *);
      k = 0;    /* marker index */
      for (i = 0; i < analysisLocusList->numLocus; i++) {
        j = analysisLocusList->pLocusIndex[i];
        if (originalLocusList.ppLocusList[j]->locusType == LOCUS_TYPE_MARKER) {
          mp_result[posIdx].pMarkers[k] = j;
          k++;
        } else {
          mp_result[posIdx].trait = i;
          traitIndex = i;
        }
      }
      analysisLocusList->traitLocusIndex = traitIndex;
      analysisLocusList->traitOrigLocus = traitLocus;

      SUBSTEP (((posIdx * 100) / numPositions), "Starting w/trait locus at %.9g (%d/%d positions)", traitPos, posIdx + 1, numPositions);

      markerSetChanged = FALSE;
      if (prevFirstMarker != mp_result[posIdx].pMarkers[0] || prevLastMarker != mp_result[posIdx].pMarkers[modelType->numMarkers - 1]) {
        /* marker set has changed */
        markerSetChanged = TRUE;
        markerLocusList.pLocusIndex[0] = mp_result[posIdx].pMarkers[0];
        prevPos = get_map_position (markerLocusList.pLocusIndex[0]);
        /* set up marker set locus list */
        for (k = 1; k < modelType->numMarkers; k++) {
          markerLocusList.pLocusIndex[k] = mp_result[posIdx].pMarkers[k];
          currPos = get_map_position (markerLocusList.pLocusIndex[k]);
          for (j = 0; j < 3; j++) {
            markerLocusList.pPrevLocusDistance[j][k] = markerLocusList.pNextLocusDistance[j][k - 1] = cm_to_recombination_fraction (currPos[j] - prevPos[j], map.mapFunction);
          }
          if (modelOptions->mapFlag == SA) {
            for (j = 1; j <= 2; j++) {
              markerLocusList.pPrevLocusDistance[j][k] = markerLocusList.pNextLocusDistance[j][k - 1] = markerLocusList.pPrevLocusDistance[0][k];
            }
          }
          prevPos = currPos;
        }       /* end of loop over the markers to set up locus list */

        {
          int length;
          char message[MAXLOGMSG + 1], *pMessage = message;
          pMessage += length = snprintf (message, MAXLOGMSG, "New markers: ");
          for (k = 0; k < modelType->numMarkers; k++)
            snprintf (pMessage, MAXLOGMSG - length, " %d(%.2f)", markerLocusList.pLocusIndex[k], *get_map_position (markerLocusList.pLocusIndex[k]));
          DETAIL (0, message);
        }

        analysisLocusList = &markerLocusList;
        xmissionMatrix = markerMatrix;
        if (modelOptions->polynomial == TRUE) {
          pedigreeSetPolynomialClearance (&pedigreeSet);
        }

        status = populate_xmission_matrix (markerMatrix, markerLocusList.numLocus, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);
        if (modelOptions->polynomial == TRUE)
          freePolys ();

        DIAG (XM, 1, {
            print_xmission_matrix (markerMatrix, markerLocusList.numLocus, 0, 0, tmpID);
	  }
        );

        /* Calculate likelihood for the marker set */
        DETAIL (0, "Determining marker set likelihood");

        for (k = 0; k < modelType->numMarkers; k++) {
          markerNameList[k] = (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[k]])->sName;
        }
        for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
          /* save the marker likelihood   */
          Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
          if (modelOptions->saveResults == TRUE) {
            pPedigree->load_flag = restoreMarker (pPedigree->sPedigreeID, (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome, modelType->numMarkers, markerNameList, &(pPedigree->markerLikelihood));
          } else {
            pPedigree->load_flag = 0;
          }
        }
        if (markerSetChanged) {
          swPushPhase ('k', "buildMM");
          /** Build a polynomial name including all involved marker ordinal numbers */
          char markerNo[8];
          if (modelOptions->polynomial == TRUE) {
            sprintf (partialPolynomialFunctionName, "MM_LC%d_C%d_P%%sM", modelRange->nlclass, (originalLocusList.ppLocusList[1])->pMapUnit->chromosome);
            for (k = 0; k < modelType->numMarkers; k++) {
              sprintf (markerNo, "_%d", markerLocusList.pLocusIndex[k]);
              strcat (partialPolynomialFunctionName, markerNo);
            }
          }
        } else
          swPushPhase ('k', "evalMM");
        ret = compute_likelihood (&pedigreeSet);
        cL[6]++;        // MP marker likelihood
        swPopPhase ('k');

        if (swProgressRequestFlag) {
          swProgressRequestFlag = FALSE;
          DETAIL (0, "Marker set likelihood evaluations %lu%% complete...\n", MAX (cL[6] * 100 / eCL[6], (posIdx + 1) * 100 / numPositions));
        }

        /* print out some statistics under dry run */
        if (modelOptions->dryRun != 0) {
          print_dryrun_stat (&pedigreeSet, -1);
        } else {
          /* save the results for marker likelihood */
          for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
            /* save the likelihood at null */
            Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];

            if (modelOptions->saveResults == TRUE) {
              if (pPedigree->load_flag == 0) {  /*save only for the pedigrees which were add for this run */
                pPedigree->markerLikelihood = pPedigree->likelihood;
                pPedigree->load_flag = saveMarker (pPedigree->sPedigreeID, (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome, modelType->numMarkers, markerNameList, &(pPedigree->markerLikelihood));
              }
            } else {
              pPedigree->markerLikelihood = pPedigree->likelihood;
            }
            pPedigree->load_flag = 0;
          }
          pedigreeSet.log10MarkerLikelihood = pedigreeSet.log10Likelihood;
        }
      }
      /* end of marker set change */
      prevFirstMarker = mp_result[posIdx].pMarkers[0];
      prevLastMarker = mp_result[posIdx].pMarkers[modelType->numMarkers - 1];
      if (markerSetChanged || prevTraitInd != mp_result[posIdx].trait)
        locusListChanged = TRUE;
      else
        locusListChanged = FALSE;
      prevTraitInd = mp_result[posIdx].trait;

      analysisLocusList = &savedLocusList;
      xmissionMatrix = altMatrix;
      /* interpolate trait postion for sex specific analysis */
      if (modelOptions->mapFlag == SS) {
        if (traitIndex == 0) {
          marker1Pos = get_map_position (analysisLocusList->pLocusIndex[1]);
          /* trait is the first one in the list */
          if (traitPos < ERROR_MARGIN && traitPos > -ERROR_MARGIN) {
            /* trait is at 0 */
            pTraitLocus->mapPosition[MAP_POS_MALE] = pTraitLocus->mapPosition[MAP_POS_FEMALE] = 0;
          } else {
            /* get the relative position on the sex average map */
            relativePos = traitPos / marker1Pos[0];
            pTraitLocus->mapPosition[MAP_POS_MALE] = relativePos * marker1Pos[MAP_POS_MALE];
            pTraitLocus->mapPosition[MAP_POS_FEMALE] = relativePos * marker1Pos[MAP_POS_FEMALE];
          }
          /* update the inter locus distance - sex averaged already done before */
          for (k = 1; k < 3; k++) {
            analysisLocusList->pNextLocusDistance[k][0] = analysisLocusList->pPrevLocusDistance[k][1] = cm_to_recombination_fraction (marker1Pos[k] - pTraitLocus->mapPosition[k], map.mapFunction);
          }
        } else if (traitIndex == modelType->numMarkers) {
          /* trait is the last one in the list */
          marker1Pos = get_map_position (analysisLocusList->pLocusIndex[modelType->numMarkers - 2]);
          marker2Pos = get_map_position (analysisLocusList->pLocusIndex[modelType->numMarkers - 1]);
          /* get the relative position on the sex average map */
          dist = marker2Pos[0] - marker1Pos[0];
          if (dist > ERROR_MARGIN) {
            relativePos = (traitPos - marker2Pos[0]) / dist;
            pTraitLocus->mapPosition[MAP_POS_MALE] = relativePos * (marker2Pos[MAP_POS_MALE] - marker1Pos[MAP_POS_MALE]) + marker2Pos[MAP_POS_MALE];
            pTraitLocus->mapPosition[MAP_POS_FEMALE] = relativePos * (marker2Pos[MAP_POS_FEMALE] - marker1Pos[MAP_POS_FEMALE]) + marker2Pos[MAP_POS_FEMALE];
          } else {
            pTraitLocus->mapPosition[MAP_POS_MALE] = traitPos - marker2Pos[0] + marker2Pos[MAP_POS_MALE];
            pTraitLocus->mapPosition[MAP_POS_FEMALE] = traitPos - marker2Pos[0] + marker2Pos[MAP_POS_FEMALE];
          }

          /* update the inter locus distance - sex averaged already done before */
          for (k = 1; k <= 2; k++) {
            analysisLocusList->pNextLocusDistance[k][traitIndex - 1] = analysisLocusList->pPrevLocusDistance[k][traitIndex] = cm_to_recombination_fraction (pTraitLocus->mapPosition[k] - marker2Pos[k], map.mapFunction);
          }

        } else {
          /* trait is in between two markers */
          marker1Pos = get_map_position (analysisLocusList->pLocusIndex[traitIndex - 1]);
          marker2Pos = get_map_position (analysisLocusList->pLocusIndex[traitIndex + 1]);
          /* get the relative position on the sex average map */
          dist = marker2Pos[0] - marker1Pos[0];
          if (dist > ERROR_MARGIN) {
            relativePos = (traitPos - marker1Pos[0]) / dist;
            pTraitLocus->mapPosition[MAP_POS_MALE] = relativePos * (marker2Pos[MAP_POS_MALE] - marker1Pos[MAP_POS_MALE]) + marker1Pos[MAP_POS_MALE];
            pTraitLocus->mapPosition[MAP_POS_FEMALE] = relativePos * (marker2Pos[MAP_POS_FEMALE] - marker1Pos[MAP_POS_FEMALE]) + marker1Pos[MAP_POS_FEMALE];
          } else {
            pTraitLocus->mapPosition[MAP_POS_MALE] = marker1Pos[MAP_POS_MALE];
            pTraitLocus->mapPosition[MAP_POS_FEMALE] = marker1Pos[MAP_POS_FEMALE];
          }
          /* update the inter locus distance - sex averaged already done before */
          for (k = 1; k < 3; k++) {
            analysisLocusList->pNextLocusDistance[k][traitIndex - 1] = analysisLocusList->pPrevLocusDistance[k][traitIndex] = cm_to_recombination_fraction (pTraitLocus->mapPosition[k] - marker1Pos[k], map.mapFunction);
            analysisLocusList->pNextLocusDistance[k][traitIndex] = analysisLocusList->pPrevLocusDistance[k][traitIndex + 1] = cm_to_recombination_fraction (marker2Pos[k] - pTraitLocus->mapPosition[k], map.mapFunction);
          }
        }
      }

      /* the locus list has been built, go on to the analysis 
       * multipoint DT */
      if (markerSetChanged || locusListChanged) {
        if (modelOptions->polynomial == TRUE) {
          pedigreeSetPolynomialClearance (&pedigreeSet);
          status = populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);
          DIAG (XM, 1, {
              print_xmission_matrix (altMatrix, savedLocusList.numLocus, 0, 0, tmpID);
	    }
          );
          if (modelOptions->polynomial == TRUE)
            freePolys ();
        }
      }

      if (modelOptions->polynomial != TRUE)
        status = populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);

      /* For Alternative */
      DETAIL (0, "Determining combined likelihood");

      if (pTrait->type == DICHOTOMOUS) {
        for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
          Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
          /* load stored alternative likelihood if they were already stored */
          if (modelOptions->saveResults == TRUE)
            pPedigree->load_flag = restoreAlternative (pPedigree->sPedigreeID, (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome, traitPos, pPedigree->alternativeLikelihoodDT);
          else
            pPedigree->load_flag = 0;
        }

        // Build the polynomial outside of the main loops
        {
          char markerNo[8];
          sprintf (partialPolynomialFunctionName, "MDA_LC%d_C%d_P%%sM", modelRange->nlclass, (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome);
          for (k = 0; k < modelType->numMarkers; k++) {
            if (traitPos <= *get_map_position (markerLocusList.pLocusIndex[k]) && (strstr (partialPolynomialFunctionName, "_T") == NULL))
              strcat (partialPolynomialFunctionName, "_T");
            sprintf (markerNo, "_%d", markerLocusList.pLocusIndex[k]);
            strcat (partialPolynomialFunctionName, markerNo);
          }
        }
        if (strstr (partialPolynomialFunctionName, "_T") == NULL)
          strcat (partialPolynomialFunctionName, "_T");
        swPushPhase ('k', "buildMDA");
        swStart (combinedBuildSW);
        ret = compute_likelihood (&pedigreeSet);
        swStop (combinedBuildSW);
        swPopPhase ('k');
        swPushPhase ('k', "evalMDA");
        swStart (combinedComputeSW);

        for (penIdx = 0; (penIdx == 0) || (modelOptions->dryRun == 0 && penIdx < modelRange->npenet); penIdx++) {
          paramSet.penIdx = penIdx;
          for (liabIdx = 0; (liabIdx == 0) || (modelOptions->dryRun == 0 && liabIdx < modelRange->nlclass); liabIdx++) {
            pen_DD = modelRange->penet[liabIdx][0][penIdx];
            pen_Dd = modelRange->penet[liabIdx][1][penIdx];
            pen_dD = modelRange->penet[liabIdx][2][penIdx];
            pen_dd = modelRange->penet[liabIdx][3][penIdx];
            pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
            pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
            pTrait->penetrance[2][liabIdx][1][0] = pen_dD;
            pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
            pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
            pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
            pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_dD;
            pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;

            if (fpIR != NULL) {
              dk_curModel.pen[liabIdx].DD = pen_DD;
              dk_curModel.pen[liabIdx].Dd = pen_Dd;
              dk_curModel.pen[liabIdx].dD = pen_dD;
              dk_curModel.pen[liabIdx].dd = pen_dd;
            }
          }

          if (modelOptions->polynomial != TRUE)
            update_penetrance (&pedigreeSet, traitLocus);       // Only need to update trait locus

          /* Iterate over gene frequencies -- just one loop for dry-runs. */
          for (gfreqInd = 0; (gfreqInd == 0) || (modelOptions->dryRun == 0 && gfreqInd < modelRange->ngfreq); gfreqInd++) {

            /* Updated trait locus allele frequencies */
            gfreq = modelRange->gfreq[gfreqInd];
            pLocus->pAlleleFrequency[0] = gfreq;
            pLocus->pAlleleFrequency[1] = 1 - gfreq;
            paramSet.gfreqIdx = gfreqInd;
            paramSet.gfreq = gfreq;

            if (fpIR != NULL)
              dk_curModel.dgf = gfreq;

            if (modelOptions->polynomial != TRUE)
              update_locus (&pedigreeSet, traitLocus);

            ret = compute_likelihood (&pedigreeSet);
            cL[7]++;    // MP DT alternative likelihood

            if (swProgressRequestFlag) {
              swProgressRequestFlag = FALSE;

              DETAIL (0, "Combined likelihood evaluations %lu%% complete (~%lu min left)", cL[7] * 100 / eCL[7],
		      ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) * eCL[7] / cL[7] -
		       (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);

            }
            /* Print out some statistics under dry run */
            if (modelOptions->dryRun != 0) {
              print_dryrun_stat (&pedigreeSet, traitPos);
            } else {

              if (ret == -2)
                ERROR ("Negative combined likelihood");

              log10_likelihood_alternative = pedigreeSet.log10Likelihood;

              /* add the result to the right placeholder */
              for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                if (pPedigree->load_flag == 0) {
                  pPedigree->alternativeLikelihoodDT[gfreqInd]
                      [penIdx] = pPedigree->likelihood;
                }
              }
              record_mp_result (ret, &pedigreeSet, &paramSet, posIdx);
            }   /* end of not dryRun */
          }     /* end of genFreq loop */
        }       /* end of penIdx loop */

        /* end of penetrance loop */
        /* save the alternative likelihood */
        for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
          Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
          if ((modelOptions->saveResults == TRUE) && (pPedigree->load_flag == 0)) {     /*save only for the pedigrees which were add for this run */
            pPedigree->load_flag = saveAlternative (pPedigree->sPedigreeID, (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome, traitPos, pPedigree->alternativeLikelihoodDT);
          }
          pPedigree->load_flag = 0;
        }

      } /* end of MP DT */
      else
        /* multipoint QT or COMBINED */
      {
	// Build the polynomial outside of the main loops
	{
	  char markerNo[8];
	  sprintf (partialPolynomialFunctionName, "MQA_LC%d_C%d_P%%sM", modelRange->nlclass, (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome);
	  for (k = 0; k < modelType->numMarkers; k++) {
	    if (traitPos <= *get_map_position (markerLocusList.pLocusIndex[k]) && (strstr (partialPolynomialFunctionName, "_T") == NULL))
	      strcat (partialPolynomialFunctionName, "_T");
	    sprintf (markerNo, "_%d", markerLocusList.pLocusIndex[k]);
	    strcat (partialPolynomialFunctionName, markerNo);
	  }
	}
	if (strstr (partialPolynomialFunctionName, "_T") == NULL)
	  strcat (partialPolynomialFunctionName, "_T");
	swPushPhase ('k', "buildMQA");
	swStart (combinedBuildSW);
	ret = compute_likelihood (&pedigreeSet);
	swStop (combinedBuildSW);
	swPopPhase ('k');
	swPushPhase ('k', "evalMQA");
	swStart (combinedComputeSW);

        for (gfreqInd = 0; gfreqInd < modelRange->ngfreq; gfreqInd++) {
          gfreq = modelRange->gfreq[gfreqInd];
          pLocus->pAlleleFrequency[0] = gfreq;
          pLocus->pAlleleFrequency[1] = 1 - gfreq;
          paramSet.gfreqIdx = gfreqInd;
          paramSet.gfreq = gfreq;

          if (fpIR != NULL)
            dk_curModel.dgf = gfreq;

          update_locus (&pedigreeSet, traitLocus);
          /* this should be MEAN + SD */
          for (paramIdx = 0; paramIdx < modelRange->nparam; paramIdx++) {
            paramSet.paramIdx = paramIdx;
            for (penIdx = 0; penIdx < modelRange->npenet; penIdx++) {
              paramSet.penIdx = penIdx;
              breakFlag = FALSE;
              for (thresholdIdx = 0; thresholdIdx < modelRange->ntthresh; thresholdIdx++) {
                paramSet.thresholdIdx = thresholdIdx;
                for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++) {
                  mean_DD = modelRange->penet[liabIdx][0][penIdx];
                  mean_Dd = modelRange->penet[liabIdx][1][penIdx];
                  mean_dD = modelRange->penet[liabIdx][2][penIdx];
                  mean_dd = modelRange->penet[liabIdx][3][penIdx];
                  SD_DD = modelRange->param[liabIdx][0][0][paramIdx];
                  SD_Dd = modelRange->param[liabIdx][1][0][paramIdx];
                  SD_dD = modelRange->param[liabIdx][2][0][paramIdx];
                  SD_dd = modelRange->param[liabIdx][3][0][paramIdx];
                  threshold = modelRange->tthresh[liabIdx][thresholdIdx];

                  if (fpIR != NULL) {
                    dk_curModel.pen[liabIdx].DD = mean_DD;
                    dk_curModel.pen[liabIdx].Dd = mean_Dd;
                    dk_curModel.pen[liabIdx].dD = mean_dD;
                    dk_curModel.pen[liabIdx].dd = mean_dd;
                    dk_curModel.pen[liabIdx].DDSD = SD_DD;
                    dk_curModel.pen[liabIdx].DdSD = SD_Dd;
                    dk_curModel.pen[liabIdx].dDSD = SD_dD;
                    dk_curModel.pen[liabIdx].ddSD = SD_dd;
                    dk_curModel.pen[liabIdx].threshold = threshold;
                  }

                  if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
                    /* check against the hard coded constraint */
                    constraint = (1 - gfreq) * (1 - gfreq) * mean_dd * SD_dd + 2 * gfreq * (1 - gfreq) * mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
                    /*      fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
                     * constraint, gfreq, mean_DD, SD_DD, 
                     * mean_Dd, SD_DD, 
                     * mean_dd, SD_dd);
                     */
                    if (constraint >= 3.0 || constraint <= -3.0) {
                      breakFlag = TRUE;
                      break;
                    }
                  }
                  pTrait->means[liabIdx][0][0] = mean_DD;
                  pTrait->means[liabIdx][0][1] = mean_Dd;
                  pTrait->means[liabIdx][1][0] = mean_dD;
                  pTrait->means[liabIdx][1][1] = mean_dd;
                  pTrait->stddev[liabIdx][0][0] = SD_DD;
                  pTrait->stddev[liabIdx][0][1] = SD_Dd;
                  pTrait->stddev[liabIdx][1][0] = SD_dD;
                  pTrait->stddev[liabIdx][1][1] = SD_dd;

                  /* threshold for QT */
                  pTrait->cutoffValue[liabIdx] = threshold;

                }       /* liability class Index */
                if (breakFlag == TRUE)
                  continue;
                if (modelOptions->polynomial == TRUE);
                else
                  update_penetrance (&pedigreeSet, traitLocus);
                /* ready for the alternative hypothesis */
                analysisLocusList = &savedLocusList;
                xmissionMatrix = altMatrix;
                if (modelOptions->polynomial == TRUE);
                else
                  status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);

		ret = compute_likelihood (&pedigreeSet);
		cL[8]++; // MP QT alternative likelihood
		if (swProgressRequestFlag) {
		  swProgressRequestFlag = FALSE;
		  DETAIL (((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) * eCL[8] / cL[8] -
			   (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60, "Combined likelihood evaluations %lu%% complete (~%lu min left)", cL[8] * 100 / eCL[8]);
		}

                log10_likelihood_alternative = pedigreeSet.log10Likelihood;
                /* add the result to the right placeholder */
                for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                  Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                  pPedigree->alternativeLikelihoodQT[gfreqInd]
                      [penIdx][paramIdx][thresholdIdx] = pPedigree->likelihood;
                }
                record_mp_result (ret, &pedigreeSet, &paramSet, posIdx);
              } /* end of threshold loop */
            }   /* end of penetrance loop */
          }     /* end of parameter loop */
        }       /* end of gene freq */

      } /* end of QT */

      swStop (combinedComputeSW);
      swPopPhase ('k');

      /* Print out average and log10(max) and maximizing parameters */
      avgLR = mp_result[posIdx].het_lr_total / (modelRange->nalpha * mp_result[posIdx].lr_count);
      log10AvgLR = log10 (avgLR) + mp_result[posIdx].scale;
      if (avgLR > 0.214) {
        if (log10AvgLR > 8)
          ppl = 1.00;
        else
          ppl = (avgLR * avgLR) / (-5.77 + 54 * avgLR + avgLR * avgLR);
      } else
        ppl = 0;

      writeMPBRFileDetail (posIdx, traitPos, ppl, avgLR);
      writeMPMODFileDetail (posIdx, traitPos);

    }   /* end of walking down the chromosome */
  }     /* end of multipoint */

  /* only for multipoint - deallocate memory  */
  if (modelType->type == MP) {
    /* allocate space to save temporary results */
    if (modelType->trait == DT) {
#if 0
      for (gfreqInd = 0; gfreqInd < modelRange->ngfreq; gfreqInd++) {
        free (likelihoodDT[gfreqInd]);
      }
      free (likelihoodDT);
#endif
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];

        for (gfreqInd = 0; gfreqInd < modelRange->ngfreq; gfreqInd++) {
          free (pPedigree->traitLikelihoodDT[gfreqInd]);
          free (pPedigree->alternativeLikelihoodDT[gfreqInd]);
        }
        free (pPedigree->traitLikelihoodDT);
        free (pPedigree->alternativeLikelihoodDT);
      }
      free (markerNameList);
    }
  }
  DIAG (OVERALL, 1, {
      dumpTrackingStats (cL, eCL);
    }
  );

  if (fpIR != NULL) {
    if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
      free (dk_curModel.dprime);
    free (dk_curModel.pen);
  }
}
