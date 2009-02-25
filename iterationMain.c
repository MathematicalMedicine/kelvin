  /* only for multipoint - we don't handle LD under multipoint yet */
  if (modelType.type == MP) {
    /* allocate space to save temporary results */
    markerNameList = (char **) calloc (sizeof (char *), modelType.numMarkers);
    if (modelType.trait == DT) {

      /* likelihoodDT is for homoLR */
      likelihoodDT = (double **) calloc (sizeof (double *), modelRange.ngfreq);
      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
        /* second dimension is penetrance */
        likelihoodDT[gfreqInd] = (double *) calloc (sizeof (double), modelRange.npenet);
      }
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        /* first dimension is gene freq */
        pPedigree->traitLikelihoodDT = (double **) calloc (sizeof (double *), modelRange.ngfreq);
        pPedigree->alternativeLikelihoodDT = (double **) calloc (sizeof (double *), modelRange.ngfreq);
        for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
          /* second dimension is penetrance */
          pPedigree->traitLikelihoodDT[gfreqInd] = (double *) calloc (sizeof (double), modelRange.npenet);
          pPedigree->alternativeLikelihoodDT[gfreqInd] = (double *) calloc (sizeof (double), modelRange.npenet);
        }
      }
    } else {    /* QT */

      /* first dimension is pedigree */
      likelihoodQT = (double *****) calloc (sizeof (double ****), pedigreeSet.numPedigree + 1);
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree + 1; pedIdx++) {
        /* second dimension is gene freq */
        likelihoodQT[pedIdx] = (double ****) calloc (sizeof (double ***), modelRange.ngfreq);
        for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {

          /* third dimension is mean */
          likelihoodQT[pedIdx][gfreqInd] = (double ***) calloc (sizeof (double **), modelRange.npenet);
          for (penIdx = 0; penIdx < modelRange.npenet; penIdx++) {
            /* fourth dimension is SD */
            likelihoodQT[pedIdx][gfreqInd][penIdx] = (double **) calloc (sizeof (double *), modelRange.nparam);
            for (paramIdx = 0; paramIdx < modelRange.nparam; paramIdx++) {
              /* 5th dimension is threshold */
              likelihoodQT[pedIdx][gfreqInd][penIdx][paramIdx] = (double *) calloc (sizeof (double), modelRange.ntthresh);
            }
          }
        }
      }
    }
  }

  /* find out the max we need to allocate */
  /* after genotype lists have been built, we want to pre-allocate parental pair work space
   * it used to be done dynamically, but it's too costly 
   * parental pair is based on each nuclear family */
  stat_parental_pair_workspace (&pedigreeSet);

  /* after genotype elimination and tracking the max work space needed 
   * for constructing parental pair */
  allocate_parental_pair_workspace (&parentalPairSpace, modelType.numMarkers + 1);

  /* conditional likelihood storage space for each individual */
  allocate_likelihood_space (&pedigreeSet, modelType.numMarkers + 1);

  /* Assume the trait locus is the first one in the list */
  traitLocus = 0;
  pLocus = originalLocusList.ppLocusList[traitLocus];
  pTraitLocus = originalLocusList.ppLocusList[traitLocus]->pTraitLocus;
  pTrait = pTraitLocus->pTraits[traitLocus];
  if (modelType.type == TP) {
    /* Two point. */
    if (originalLocusList.pLDLoci == NULL) {
      originalLocusList.pLDLoci = (LDLoci *) malloc (sizeof (LDLoci));
      memset (originalLocusList.pLDLoci, 0, sizeof (LDLoci));
    }
    pLDLoci = &originalLocusList.pLDLoci[0];
    originalLocusList.numLDLoci = 1;

    if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) {
      /* fake some LD information to simplify looping */
      pLDLoci->numAllele1 = 2;
      pLDLoci->ppDPrime = (double **) malloc (sizeof (double *));
      pLDLoci->ppDPrime[0] = (double *) malloc (sizeof (double));
      pLDLoci->ppDValue = (double **) malloc (sizeof (double *));
      pLDLoci->ppDValue[0] = (double *) malloc (sizeof (double));
      pLDLoci->ppHaploFreq = (double **) malloc (sizeof (double *) * 2);
      pLDLoci->ppHaploFreq[0] = (double *) malloc (sizeof (double) * 2);
      pLDLoci->ppHaploFreq[1] = (double *) malloc (sizeof (double) * 2);

      /* initialize it */
      pLDLoci->ppDPrime[0][0] = 0;
    }

    locusList = &savedLocusList;
    savedLocusList.numLocus = 2;
    savedLocusList.pLocusIndex = (int *) malloc (sizeof (int) * savedLocusList.numLocus);
    for (i = 0; i < 3; i++) {
      savedLocusList.pPrevLocusDistance[i] = (double *) malloc (sizeof (double) * savedLocusList.numLocus);
      savedLocusList.pNextLocusDistance[i] = (double *) malloc (sizeof (double) * savedLocusList.numLocus);
      savedLocusList.pPrevLocusDistance[i][0] = -1;
      savedLocusList.pNextLocusDistance[i][1] = -1;
    }

    if (modelOptions.polynomial == TRUE) {
      status =
        populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0,
                                  -1, -1, 0);
      holdAllPolys ();
    }

    total_count = modelRange.npenet * modelRange.ngfreq * modelRange.nalpha;

    if (modelOptions.markerAnalysis == FALSE) {
      savedLocusList.traitLocusIndex = 0;
      savedLocusList.traitOrigLocus = 0;
    } else {
      savedLocusList.traitLocusIndex = -1;
      savedLocusList.traitOrigLocus = -1;
    }

    for (loc1 = 0; loc1 < originalLocusList.numLocus - 1; loc1++) {
      savedLocusList.pLocusIndex[0] = loc1;
      pLocus1 = originalLocusList.ppLocusList[loc1];
      if (modelOptions.markerAnalysis != FALSE && pLocus1->locusType != LOCUS_TYPE_MARKER)
        continue;

      for (loc2 = loc1 + 1; loc2 < originalLocusList.numLocus; loc2++) {
        pLocus2 = originalLocusList.ppLocusList[loc2];
        if (pLocus2->locusType != LOCUS_TYPE_MARKER)
          continue;
        savedLocusList.pLocusIndex[1] = loc2;

#ifndef SIMPLEPROGRESS
	if (modelOptions.markerAnalysis == MM)
	  fprintf (stdout, "Starting w/loci %s(%d alleles) and %s(%d alleles\n", 
		   pLocus1->sName, pLocus1->numOriginalAllele, pLocus2->sName, pLocus2->numOriginalAllele);
	else
	  fprintf (stdout, "Starting w/loci %s(%d alleles) and %s(%d alleles) (%d of %d pairs)\n",
		   pLocus1->sName, pLocus1->numOriginalAllele, pLocus2->sName, pLocus2->numOriginalAllele,
		   loc2, originalLocusList.numLocus - 1);
#endif

        /* find out number of alleles this marker locus has */
        if (modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM) {
          /* get the LD parameters */
          pLambdaCell = findLambdas (&modelRange, pLocus1->numOriginalAllele, pLocus2->numOriginalAllele);
          reallocate_LD_loci (pLDLoci, pLocus1->numOriginalAllele, pLocus2->numOriginalAllele);

	  // Create these variables ahead of likelihood polynomial build in hopes of preventing in-build creation.

	  if (modelOptions.polynomial == TRUE) {
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

        loopMarkerFreqFlag = 0;
        if (modelRange.nafreq >= 2 && modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM && pLocus2->numOriginalAllele == 2) {
          loopMarkerFreqFlag = 1;
        } else if (modelRange.nafreq == 0) {
          /* add a fake one to facilitate loops and other handlings */
          addAlleleFreq (&modelRange, pLocus2->pAlleleFrequency[0]);
        } else {
          modelRange.nafreq = 1;
          modelRange.afreq[0] = pLocus2->pAlleleFrequency[0];
        }

        /* allocate/initialize result storage */
        initialize_tp_result_storage ();
	//	dumpTrackingStats(cL, eCL);

        /* we will force marker allele frequency loop to execute at least once */
        for (mkrFreqIdx = 0; mkrFreqIdx == 0 || mkrFreqIdx < modelRange.nafreq; mkrFreqIdx++) {
          mkrFreq = pLocus2->pAlleleFrequency[0];
          /* we should only loop over marker allele frequency under twopoint
           * and when markers are SNPs (only have two alleles) */
          if (loopMarkerFreqFlag) {
            mkrFreq = modelRange.afreq[mkrFreqIdx];
            /* update the locus */
            pLocus2->pAlleleFrequency[0] = mkrFreq;
            pLocus2->pAlleleFrequency[1] = 1 - mkrFreq;
            if (modelOptions.polynomial == TRUE);
            else
              update_locus (&pedigreeSet, loc2);
          }
          /* Loop over the penetrances, genefrequencies, thetas and call
           * the likelihood calculation, storing each value obtained to
           * disk. */
          for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
            gfreq = modelRange.gfreq[gfreqInd];
            // WHAT ON EARTH IS THIS ALL ABOUT? &&&
            if (1 && modelOptions.markerAnalysis == FALSE) {
              pLocus->pAlleleFrequency[0] = gfreq;
              pLocus->pAlleleFrequency[1] = 1 - gfreq;
              if (modelOptions.polynomial == TRUE);
              else
                update_locus (&pedigreeSet, loc1);
            }

            /* clear Dprime combination impossible flag */
            memset (pLambdaCell->impossibleFlag, 0, sizeof (int) * pLambdaCell->ndprime);
            /* set up haplotype frequencies */
            for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
              if (isDPrime0 (pLambdaCell->lambda[dprimeIdx], pLambdaCell->m, pLambdaCell->n))
                dprime0Idx = dprimeIdx;
              status = setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprimeIdx);
              if (status < 0)
                pLambdaCell->impossibleFlag[dprimeIdx] = 1;
            }

            if (modelType.trait == DICHOTOMOUS) {

              for (penIdx = 0; penIdx < modelRange.npenet; penIdx++) {
                if (modelOptions.markerAnalysis == FALSE && pLocus1->locusType == LOCUS_TYPE_TRAIT) {
                  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                    pen_DD = modelRange.penet[liabIdx][0][penIdx];
                    pen_Dd = modelRange.penet[liabIdx][1][penIdx];
                    pen_dD = modelRange.penet[liabIdx][2][penIdx];
                    pen_dd = modelRange.penet[liabIdx][3][penIdx];
                    pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
                    pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
                    pTrait->penetrance[2][liabIdx][1][0] = pen_dD;
                    pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
                    pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
                    pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
                    pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_dD;
                    pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;
                  }
                  if (modelOptions.polynomial == TRUE);
                  else
                    update_penetrance (&pedigreeSet, traitLocus);
                }
                /* get the likelihood at 0.5 first and LD=0 */
                if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
                  set_null_dprime (pLDLoci);
                  copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprime0Idx]);
                  copy_DValue (pLDLoci, pLambdaCell->DValue[dprime0Idx]);
                  KASSERT (pLambdaCell->impossibleFlag[dprime0Idx] == 0,
                           "Haplotype frequency combination impossible at LE. Exiting!\n");
                }
                for (k = 0; k < 3; k++) {
                  locusList->pNextLocusDistance[k][0] = 0.5;
                  locusList->pPrevLocusDistance[k][1] = 0.5;
                }

                if (modelOptions.polynomial == TRUE);
                else
                  status =
                    populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2,
                                              initialHetProbAddr, 0, -1, -1, 0);

                /* If we're not on the first iteration, it's not a polynomial build, so
                 * show progress at 1 minute intervals. Have a care to avoid division by zero. */
		sprintf (partialPolynomialFunctionName, "CL0_P%%s_%s_%s",
			 pLocus1->sName, pLocus2->sName);
                if (gfreqInd != 0 || penIdx != 0) {
		  pushStatus ('k', "evalCL0");
		  //                  swStart (combinedComputeSW);
                  compute_likelihood (&pedigreeSet);
                  cL[0]++;
		  //                  swStop (combinedComputeSW);
                  if (statusRequestSignal) {
                    statusRequestSignal = FALSE;
                    if (cL[0] > 1) {    // The first time thru we have no basis for estimation
                      fprintf (stdout, "%s %lu%% complete (~%lu min left)\r",
                               "Calculations", (cL[0] + cL[1]) * 100 / (eCL[0] + eCL[1]),
                               ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                                (eCL[0] + eCL[1]) / (cL[0] + cL[1]) -
                                (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
                      fflush (stdout);
                    }
                  }
		  popStatus ('k');
                } else { // This _is_ the first iteration
		  pushStatus ('k', "buildCL0");
		  swStart (combinedBuildSW);
		  compute_likelihood (&pedigreeSet);
		  cL[0]++;
		  swStop (combinedBuildSW);
		  fprintf (stdout, "%s %lu%% complete\r", "Calculations", (cL[0] + cL[1]) * 100 / (eCL[0] + eCL[1]));
		  fflush (stdout);
		  popStatus ('k');
		}
                if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
                  fprintf (stderr, "Theta 0.5 has likelihood 0\n");
                  fprintf (stderr, "dgf=%f\n", gfreq);
                  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                    pen_DD = modelRange.penet[liabIdx][0][penIdx];
                    pen_Dd = modelRange.penet[liabIdx][1][penIdx];
                    pen_dD = modelRange.penet[liabIdx][2][penIdx];
                    pen_dd = modelRange.penet[liabIdx][3][penIdx];
		    if (modelOptions.imprintingFlag)
		      fprintf (stderr, "Liab %d penentrance %f %f %f %f\n", liabIdx + 1, pen_DD, pen_Dd, pen_dD, pen_dd);
		    else
		      fprintf (stderr, "Liab %d penentrance %f %f %f\n", liabIdx + 1, pen_DD, pen_Dd, pen_dd);
                  }
                  exit (EXIT_FAILURE);
                }
                /* save the results for NULL */
                for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                  /* save the likelihood at null */
                  pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                  pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
                }

                log10_likelihood_null = pedigreeSet.log10Likelihood;
		//&&&		fprintf (stderr, "About to do %d iterations for ndprime at loc1/loc2: %d/%d\n",
		//			 pLambdaCell->ndprime, loc1, loc2);
                for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
                  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
                    copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);
                    if (pLambdaCell->impossibleFlag[dprimeIdx] != 0) {
		      // If we're going to bail at this point, add the progress count loop factor
		      cL[1] += modelRange.ntheta;
		      //&&&		      for (i=0; i<modelRange.ntheta; i++)
		      //			fprintf (stderr, "loc1,loc2,mkrFreqIdx,gfreqInd,penIdx,dprimeIdx,thetaInd:\t%d\t%d\t%d\t%d\t%d\t%d\tskipped\n",
		      //				 loc1,loc2,mkrFreqIdx,gfreqInd,penIdx,dprimeIdx);
                      continue;
		    }
                    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
                    copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);
                    /* calculate R square if the marker is a SNP */
                    if (R_square_flag == TRUE)
                      R_square =
                        calculate_R_square (pLocus1->
                                            pAlleleFrequency[0], pLocus2->pAlleleFrequency[0], pLDLoci->ppDValue[0][0]);
                    else
                      R_square = -1;
                  }

                  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++) {
                    if (modelOptions.mapFlag == SA) {
                      theta[0] = modelRange.theta[0][thetaInd];
                      theta[1] = modelRange.theta[1][thetaInd];
                      for (k = 0; k < 3; k++) {
                        locusList->pNextLocusDistance[k][0] = theta[0];
                        locusList->pPrevLocusDistance[k][1] = theta[0];
                      }
                    } else {
                      locusList->pNextLocusDistance[MAP_MALE][0] =
                        locusList->pPrevLocusDistance[MAP_MALE][1] = modelRange.theta[0][thetaInd];
                      locusList->pNextLocusDistance[MAP_FEMALE][0] =
                        locusList->pPrevLocusDistance[MAP_FEMALE][1] = modelRange.theta[1][thetaInd];
                    }

                    if (modelOptions.polynomial == TRUE);
                    else
                      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,
                                                         initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);

		    sprintf (partialPolynomialFunctionName, "CL1_P%%s_%s_%s",
			     pLocus1->sName, pLocus2->sName);
                    swStart (combinedComputeSW);
                    compute_likelihood (&pedigreeSet);
                    cL[1]++;
		    //&&&		    fprintf (stderr, "loc1,loc2,mkrFreqIdx,gfreqInd,penIdx,dprimeIdx,thetaInd:\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
		    //			     loc1,loc2,mkrFreqIdx,gfreqInd,penIdx,dprimeIdx,thetaInd);
                    swStop (combinedComputeSW);
                    if (statusRequestSignal) {
                      statusRequestSignal = FALSE;
                      if (cL[1] > 1) {  // The first time thru we have no basis for estimation
                        fprintf (stdout, "%s %lu%% complete (~%lu min left)\r",
                                 "Calculations", (cL[0] + cL[1]) * 100 / (eCL[0] + eCL[1]),
                                 ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                                  (eCL[0] + eCL[1]) / (cL[0] + cL[1]) -
                                  (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
                        fflush (stdout);
                      }
                    }

                    log10_likelihood_alternative = pedigreeSet.log10Likelihood;
                    if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
                      log10_likelihood_ratio = 0;
                    } else {
                      log10_likelihood_ratio = log10_likelihood_alternative - log10_likelihood_null;
                    }
                    /* check for overflow problem !!! */
                    if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
                      likelihood_ratio = DBL_MAX;
                      tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_total = DBL_MAX;
                    } else
                      /* check for underflow problem too !!! */
                    if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
                      likelihood_ratio = 0;
                    } else {
                      likelihood_ratio = pow (10.0, log10_likelihood_ratio);
                      tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_total += likelihood_ratio;
                    }
                    tp_result[dprimeIdx][thetaInd]
                      [mkrFreqIdx].lr_count++;

                    /* caculating the HET */
                    for (j = 0; j < modelRange.nalpha; j++) {
                      alphaV = modelRange.alpha[j];
                      alphaV2 = 1 - alphaV;
                      if (alphaV2 < 0)
                        alphaV2 = 0;
                      log10HetLR = 0;
                      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                        pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                        homoLR = pPedigree->likelihood / pedigreeSet.nullLikelihood[pedIdx];
                        tmp = log10 (alphaV * homoLR + (1 - alphaV));
                        log10HetLR += tmp * pPedigree->pCount[loc2]; // Use the pedigree weight from count file (CF)
                      }
                      if (log10HetLR >= DBL_MAX_10_EXP - 1) {
                        hetLR = DBL_MAX;
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total = DBL_MAX;
                      } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
                        hetLR = 0;
                      } else {
                        hetLR = pow (10, log10HetLR);
                        if (modelType.ccFlag)
                          /* scale it to prevent overflow */
                          tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total += hetLR / total_count;
                        else
                          tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total += hetLR;
                      }
                      if (tp_result[dprimeIdx][thetaInd]
                          [mkrFreqIdx].max_penIdx < 0 || hetLR > tp_result[dprimeIdx][thetaInd]
                          [mkrFreqIdx].max_lr) {
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_lr = hetLR;
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_alpha = alphaV;
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_gfreq = gfreq;
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_penIdx = penIdx;
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].R_square = R_square;
                        tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_mf = mkrFreq;
                      }
                    }   /* end of calculating HET LR */
                  }     /* end of theta loop */
                }       /* end of D prime loop */
                if (modelOptions.markerAnalysis != FALSE) {
                  /* marker to marker analysis, marker allele frequency is fixed */
                  gfreqInd = modelRange.ngfreq;
                  break;
                }
                if (modelOptions.markerAnalysis != FALSE) {
                  /* marker to marker analysis, penetrance stays at 1 */
                  break;
                }
              } /* end of penetrance loop */
            } /* end of TP */
            else
              /* should be QT or COMBINED - twopoint */
            {
              /* this should be MEAN + SD */
              for (paramIdx = 0; (paramIdx == 0 && modelType.distrib == QT_FUNCTION_CHI_SQUARE)
                   || (modelType.distrib != QT_FUNCTION_CHI_SQUARE && paramIdx < modelRange.nparam); paramIdx++) {
                for (penIdx = 0; penIdx < modelRange.npenet; penIdx++) {
                  breakFlag = FALSE;
                  for (thresholdIdx = 0; thresholdIdx < modelRange.ntthresh; thresholdIdx++) {
                    if (modelOptions.markerAnalysis == FALSE) {
                      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                        mean_DD = modelRange.penet[liabIdx][0][penIdx];
                        mean_Dd = modelRange.penet[liabIdx][1][penIdx];
                        mean_dD = modelRange.penet[liabIdx][2][penIdx];
                        mean_dd = modelRange.penet[liabIdx][3][penIdx];
                        SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
                        SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
                        SD_dD = modelRange.param[liabIdx][2][0][paramIdx];
                        SD_dd = modelRange.param[liabIdx][3][0][paramIdx];
                        /* threshold for QT */
                        threshold = modelRange.tthresh[liabIdx][thresholdIdx];
                        /* check against the hard coded constraint */
                        if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
                          constraint =
                            (1 - gfreq) * (1 - gfreq) * mean_dd * SD_dd + 2 * gfreq * (1 - gfreq) * mean_Dd * SD_Dd +
                            gfreq * gfreq * mean_DD * SD_DD;
                          if (constraint >= 3.0 || constraint <= -3.0) {
                            breakFlag = TRUE;
                            break;
                          }
                        }
                        pTrait->means[liabIdx][0][0] = mean_DD;
                        pTrait->means[liabIdx][0][1] = mean_Dd;
                        pTrait->means[liabIdx][1][0] = mean_Dd;
                        pTrait->means[liabIdx][1][1] = mean_dd;
                        pTrait->stddev[liabIdx][0][0] = SD_DD;
                        pTrait->stddev[liabIdx][0][1] = SD_Dd;
                        pTrait->stddev[liabIdx][1][0] = SD_Dd;
                        pTrait->stddev[liabIdx][1][1] = SD_dd;

                        /* threshold for QT */
                        pTrait->cutoffValue[liabIdx] = threshold;

                      } /* liability class Index */
                      if (breakFlag == TRUE)
                        continue;
                      if (modelOptions.polynomial == TRUE);
                      else
                        update_penetrance (&pedigreeSet, traitLocus);
                    }
                    /* marker to marker analysis */
                    /* get the likelihood at 0.5 first and LD=0 */
                    if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
                      set_null_dprime (pLDLoci);
                      copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprime0Idx]);
                      copy_DValue (pLDLoci, pLambdaCell->DValue[dprime0Idx]);

                      KASSERT (pLambdaCell->impossibleFlag[dprime0Idx] == 0,
                               "Haplotype frequency combination impossible at LE. Exiting!\n");
                    }
                    for (k = 0; k < 3; k++) {
                      locusList->pNextLocusDistance[k][0] = 0.5;
                      locusList->pPrevLocusDistance[k][1] = 0.5;
                    }
                    if (modelOptions.polynomial == TRUE);
                    else
                      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,
                                                         initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);

		    /* If we're not on the first iteration, it's not a polynomial build, so
		     * show progress at 1 minute intervals. Have a care to avoid division by zero. */

		    sprintf (partialPolynomialFunctionName, "CL2_P%%s_%s_%s",
			     pLocus1->sName, pLocus2->sName);
		    if (gfreqInd != 0 || penIdx != 0 || paramIdx != 0 || thresholdIdx != 0) {
		      pushStatus ('k', "evalCL2");
		      //		      swStart (combinedComputeSW);
		      compute_likelihood (&pedigreeSet);
		      cL[2]++;
		      //		      swStop (combinedComputeSW);
		      if (statusRequestSignal) {
			statusRequestSignal = FALSE;
			if (cL[2] > 1) {    // The first time thru we have no basis for estimation
			  fprintf (stdout, "%s %lu%% complete (~%lu min left)\r",
                               "Calculations", (cL[2] + cL[3]) * 100 / (eCL[2] + eCL[3]),
                               ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                                (eCL[2] + eCL[3]) / (cL[2] + cL[3]) -
                                (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
			  fflush (stdout);
			}
		      }
		      popStatus ('k');
		    } else { // This _is_ the first iteration
		      pushStatus ('k', "buildCL2");
		      swStart (combinedBuildSW);
		      compute_likelihood (&pedigreeSet);
		      cL[2]++;
		      swStop (combinedBuildSW);
		      fprintf (stdout, "%s %lu%% complete\r", "Calculations", (cL[2] + cL[3]) * 100 / (eCL[2] + eCL[3]));
		      fflush (stdout);
		      popStatus ('k');
		    }

                    if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
                      fprintf (stderr, "Theta 0.5 has likelihood 0\n");
                      fprintf (stderr, "dgf=%f\n", gfreq);
                      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                        pen_DD = modelRange.penet[liabIdx][0][penIdx];
                        pen_Dd = modelRange.penet[liabIdx][1][penIdx];
                        pen_dD = modelRange.penet[liabIdx][2][penIdx];
                        pen_dd = modelRange.penet[liabIdx][3][penIdx];
			if (modelOptions.imprintingFlag)
			  fprintf (stderr, "Liab %d penentrance %f %f %f %f\n", liabIdx + 1, pen_DD, pen_Dd, pen_dD, pen_dd);
			else
			  fprintf (stderr, "Liab %d penentrance %f %f %f\n", liabIdx + 1, pen_DD, pen_Dd, pen_dd);
                      }
                      exit (EXIT_FAILURE);
                    }
                    for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                      /* save the likelihood at null */
                      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                      pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
                    }
                    log10_likelihood_null = pedigreeSet.log10Likelihood;
                    for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
                      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
                        copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);
			if (pLambdaCell->impossibleFlag[dprimeIdx] != 0) {
			  // If we're going to bail at this point, add the progress count loop factor
			  cL[3] += modelRange.ntheta;
                          continue;
			}
                        copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
                        copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);
                      }
                      for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++) {

                        if (modelOptions.mapFlag == SA) {
                          theta[0] = modelRange.theta[0][thetaInd];
                          theta[1] = modelRange.theta[1][thetaInd];
                          for (k = 0; k < 3; k++) {
                            locusList->pNextLocusDistance[k][0] = theta[0];
                            locusList->pPrevLocusDistance[k][1] = theta[0];
                          }
                        } else {
                          locusList->pNextLocusDistance[MAP_MALE][0] =
                            locusList->pPrevLocusDistance[MAP_MALE][1] = modelRange.theta[0][thetaInd];
                          locusList->pNextLocusDistance[MAP_FEMALE][0] =
                            locusList->pPrevLocusDistance[MAP_FEMALE][1] = modelRange.theta[1][thetaInd];
                        }

                        if (modelOptions.polynomial == TRUE);
                        else
                          status =
                            populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2,
                                                      initialHetProbAddr, 0, -1, -1, 0);

                        strcpy (partialPolynomialFunctionName, "cL3_P%s");
			swStart (combinedComputeSW);
                        compute_likelihood (&pedigreeSet);
                        cL[3]++;
			swStop (combinedComputeSW);
			if (statusRequestSignal) {
			  statusRequestSignal = FALSE;
			  if (cL[3] > 1) {  // The first time thru we have no basis for estimation
			    fprintf (stdout, "%s %lu%% complete (~%lu min left)\r",
				     "Calculations", (cL[2] + cL[3]) * 100 / (eCL[2] + eCL[3]),
				     ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
				      (eCL[2] + eCL[3]) / (cL[2] + cL[3]) -
				      (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
			    fflush (stdout);
			  }
			}

                        log10_likelihood_alternative = pedigreeSet.log10Likelihood;
                        if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
                          log10_likelihood_ratio = 0;
                        } else {
                          log10_likelihood_ratio = log10_likelihood_alternative - log10_likelihood_null;
                        }
                        /* check for overflow problem !!! */
                        if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
                          likelihood_ratio = DBL_MAX;
                          tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_total = DBL_MAX;
                        } else
                          /* check for underflow problem too !!! */
                        if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
                          likelihood_ratio = 0;
                        } else {
                          likelihood_ratio = pow (10.0, log10_likelihood_ratio);
                          tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_total += likelihood_ratio;
                        }
                        tp_result[dprimeIdx][thetaInd]
                          [mkrFreqIdx].lr_count++;
                        /* caculating the HET */
                        for (j = 0; j < modelRange.nalpha; j++) {
                          alphaV = modelRange.alpha[j];
                          alphaV2 = 1 - alphaV;
                          if (alphaV2 < 0)
                            alphaV2 = 0;
                          log10HetLR = 0;
                          for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                            pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                            homoLR = pPedigree->likelihood / pedigreeSet.nullLikelihood[pedIdx];
                            log10HetLR += log10 (alphaV * homoLR + alphaV2);
                          }
                          if (log10HetLR >= DBL_MAX_10_EXP - 1) {
                            hetLR = DBL_MAX;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total = DBL_MAX;
                          } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
                            hetLR = 0;
                          } else {
                            adjustedHetLR = hetLR = pow (10, log10HetLR);
                            /* for threshold parameter, we need to make sure the weighting is even */
                            if (1 || modelType.distrib == QT_FUNCTION_CHI_SQUARE) {
                              if (modelRange.ntthresh == 1) {
                                adjustedHetLR *= 2 * (modelType.maxThreshold - modelType.minThreshold);
                              } else if (thresholdIdx == modelRange.ntthresh - 1) {
                                adjustedHetLR *=
                                  (2 * modelType.maxThreshold - threshold - modelRange.tthresh[0][thresholdIdx - 1]);
                              } else if (thresholdIdx == 0) {
                                adjustedHetLR *=
                                  (threshold + modelRange.tthresh[0][thresholdIdx + 1] - 2 * modelType.minThreshold);
                              } else {
                                adjustedHetLR *=
                                  (modelRange.tthresh[0][thresholdIdx + 1] - modelRange.tthresh[0][thresholdIdx - 1]);
                              }
                            }
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total += adjustedHetLR;
                          }
                          if (tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_penIdx < 0
                              || hetLR > tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_lr) {
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_lr = hetLR;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_alpha = alphaV;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_gfreq = gfreq;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_penIdx = penIdx;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_paramIdx = paramIdx;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_thresholdIdx = thresholdIdx;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].R_square = R_square;
                            tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_mf = mkrFreq;
                          }
                        }
                      } /* end of theta */
                    }   /* end of D prime */
                    if (modelOptions.markerAnalysis != FALSE)
                      break;
                  }     /* end of threshold loop */
                  if (modelOptions.markerAnalysis != FALSE)
                    break;
                }       /* end of penetrance loop */
                if (modelOptions.markerAnalysis != FALSE)
                  break;
              } /* end of parameter loop */
              if (modelOptions.markerAnalysis != FALSE)
                break;
            }   /* end of QT */
          }     /* end of gene freq */
          /* only loop marker allele frequencies when doing LD */
          if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
            break;
          /* we can only do SNPs when looping over marker allele frequency */
          if (pLocus2->numOriginalAllele > 2)
            break;
        }       /* end of marker allele frequency looping */

        /* calculate the average BR */
        get_average_LR (tp_result);

	write2ptBRFile ();
        write2ptMODFile ();
	writeMMFileDetail ();
	writePPLFileDetail ();

        /* need to clear polynomial */

        if (modelOptions.polynomial == TRUE) {
          pedigreeSetPolynomialClearance (&pedigreeSet);
        }


        if (modelOptions.markerAnalysis == ADJACENTMARKER)
          loc2 = originalLocusList.numLocus;

#ifndef SIMPLEPROGRESS
        fprintf (stdout, "\n");
#endif
	/* free two point result storage */
	free_tp_result_storage ();
      } /* end of looping second locus - loc2 */
      /* if we are doing trait marker, then we are done */
      /* Used to read: modelOptions.markerToMarker != TRUE which
       * is the same as markerAnalysis == FALSE as long as the old
       * markerToMarker and adjacentMarker flags were truly
       * orthogonal. Otherwise, it should be markerAnalysis !=
       * ADJACENTMARKER. */
      if (modelOptions.markerAnalysis == FALSE)
        loc1 = originalLocusList.numLocus;
    }   /* end of looping first locus - loc1 */
  } /* end of two point */
  else {        /* multipoint */

    /* marker set locus list for each position */
    markerLocusList.maxNumLocus = modelType.numMarkers;
    markerLocusList.numLocus = modelType.numMarkers;
    markerLocusList.traitOrigLocus = -1;
    markerLocusList.traitLocusIndex = -1;
    markerLocusList.pLocusIndex = (int *) calloc (markerLocusList.maxNumLocus, sizeof (int));
    for (k = 0; k < 3; k++) {
      markerLocusList.pPrevLocusDistance[k] = (double *) calloc (markerLocusList.maxNumLocus, sizeof (double));
      markerLocusList.pNextLocusDistance[k] = (double *) calloc (markerLocusList.maxNumLocus, sizeof (double));
    }

    /* assuming we always have trait in the analysis - this may not be true 
     * need to add code to process marker to marker analysis under multipoin
     */
    savedLocusList.numLocus = modelType.numMarkers + 1;
    savedLocusList.maxNumLocus = modelType.numMarkers + 1;
    savedLocusList.pLocusIndex = (int *) calloc (savedLocusList.maxNumLocus, sizeof (int));
    for (k = 0; k < 3; k++) {
      savedLocusList.pPrevLocusDistance[k] = (double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
      savedLocusList.pNextLocusDistance[k] = (double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
    }

    /* Allocate storage to calculate the trait likelihood independent of the trait position */
    traitLocusList.numLocus = 1;
    traitLocusList.maxNumLocus = 1;
    traitLocusList.traitLocusIndex = 0;
    traitLocusList.traitOrigLocus = traitLocus;
    traitLocusList.pLocusIndex = (int *) calloc (traitLocusList.maxNumLocus, sizeof (int));
    traitLocusList.pLocusIndex[0] = 0;
    for (k = 0; k < 3; k++) {
      traitLocusList.pPrevLocusDistance[k] = (double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
      traitLocusList.pNextLocusDistance[k] = (double *) calloc (savedLocusList.maxNumLocus, sizeof (double));

      traitLocusList.pPrevLocusDistance[k][0] = -1;
      traitLocusList.pNextLocusDistance[k][0] = -1;
    }
    /* populate the trait xmission matrix */
    locusList = &traitLocusList;
    xmissionMatrix = traitMatrix;
    status = populate_xmission_matrix (traitMatrix, 1, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);
    if (modelOptions.polynomial == TRUE)
      holdAllPolys ();

    /* For trait likelihood */
#ifndef SIMPLEPROGRESS
    fprintf (stdout, "Determining trait likelihood...\n");
#else
    fprintf (stdout, "Calculations 0%% complete\r");
    fflush (stdout);
#endif

    locusList = &traitLocusList;
    xmissionMatrix = traitMatrix;
    if (pTrait->type == DICHOTOMOUS) {
      /* load all saved trait likelihood */
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        if (modelOptions.saveResults == TRUE) {
          pPedigree->load_flag = restoreTrait (modelOptions.sexLinked, pPedigree->sPedigreeID, pPedigree->traitLikelihoodDT);
        } else
          pPedigree->load_flag = 0;
      }

      for (penIdx = 0; (penIdx == 0) || (modelOptions.dryRun == 0 && penIdx < modelRange.npenet); penIdx++) {
        for (liabIdx = 0; (liabIdx == 0) || (modelOptions.dryRun == 0 && liabIdx < modelRange.nlclass); liabIdx++) {
          pen_DD = modelRange.penet[liabIdx][0][penIdx];
          pen_Dd = modelRange.penet[liabIdx][1][penIdx];
          pen_dD = modelRange.penet[liabIdx][2][penIdx];
          pen_dd = modelRange.penet[liabIdx][3][penIdx];
          pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
          pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
          pTrait->penetrance[2][liabIdx][1][0] = pen_dD;
          pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
          pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
          pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
          pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_dD;
          pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;
        }

        if (modelOptions.polynomial == TRUE);
        else
          /* only need to update trait locus */
          update_penetrance (&pedigreeSet, traitLocus);

        /* Iterate over gene frequencies, but only one time thru if doing a dry-run. */
        for (gfreqInd = 0; (gfreqInd == 0) || (modelOptions.dryRun == 0 && gfreqInd < modelRange.ngfreq); gfreqInd++) {

          /* updated trait locus allele frequencies */
          gfreq = modelRange.gfreq[gfreqInd];
          pLocus->pAlleleFrequency[0] = gfreq;
          pLocus->pAlleleFrequency[1] = 1 - gfreq;

          if (modelOptions.polynomial == TRUE);
          else
            update_locus (&pedigreeSet, traitLocus);

          /* Compute the likelihood for the trait */
          sprintf (partialPolynomialFunctionName, "TL4_P%%sSL%d", modelOptions.sexLinked);
          compute_likelihood (&pedigreeSet);
          cL[4]++;
#ifndef SIMPLEPROGRESS
          if (cL[4] % MAX (1, eCL[4] / 5) == 1) {
            fprintf (stdout, "Trait likelihood evaluations %lu%% complete\r", cL[4] * 100 / eCL[4]);
            fflush (stdout);
          }
#endif
          if (modelOptions.dryRun != 0)
            continue;

          if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
            fprintf (stderr, "Trait has likelihood 0\n");
            fprintf (stderr, "dgf=%f\n", gfreq);
            for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
              pen_DD = modelRange.penet[liabIdx][0][penIdx];
              pen_Dd = modelRange.penet[liabIdx][1][penIdx];
              pen_dD = modelRange.penet[liabIdx][2][penIdx];
              pen_dd = modelRange.penet[liabIdx][3][penIdx];
	    if (modelOptions.imprintingFlag)
              fprintf (stderr, "Liab %d penentrance %f %f %f %f\n", liabIdx + 1, pen_DD, pen_Dd, pen_dD, pen_dd);
	    else
              fprintf (stderr, "Liab %d penentrance %f %f %f\n", liabIdx + 1, pen_DD, pen_Dd, pen_dd);
            }
            exit (EXIT_FAILURE);
          }
          /* save the results for NULL */
          for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
            /* save the likelihood at null */
            pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
            if (pPedigree->load_flag == 0) {    /*update only for the pedigrees which were add for this run */
              pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
              pPedigree->traitLikelihoodDT[gfreqInd][penIdx] = pPedigree->likelihood;
            }
          }

          log10_likelihood_null = pedigreeSet.log10Likelihood;
          likelihoodDT[gfreqInd][penIdx] = log10_likelihood_null;
        }       /* gfreq */
      } /* pen */

      /* save all  trait likelihood which were created in this run */
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        /* save the likelihood at null */
        pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        if ((modelOptions.saveResults == TRUE) && (pPedigree->load_flag == 0)) {        /*save only for the pedigrees which were add for this run */
          pPedigree->load_flag = saveTrait (modelOptions.sexLinked, pPedigree->sPedigreeID, pPedigree->traitLikelihoodDT);
        } else {
          pPedigree->load_flag = 0;
        }
      }
    } else
      /* multipoint QT or COMBINED */
    {
      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
        gfreq = modelRange.gfreq[gfreqInd];
        pLocus->pAlleleFrequency[0] = gfreq;
        pLocus->pAlleleFrequency[1] = 1 - gfreq;

        update_locus (&pedigreeSet, traitLocus);
        /* this should be MEAN + SD */
        for (paramIdx = 0; paramIdx < modelRange.nparam; paramIdx++) {
          for (penIdx = 0; penIdx < modelRange.npenet; penIdx++) {
            breakFlag = FALSE;
            for (thresholdIdx = 0; thresholdIdx < modelRange.ntthresh; thresholdIdx++) {
              for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                mean_DD = modelRange.penet[liabIdx][0][penIdx];
                mean_Dd = modelRange.penet[liabIdx][1][penIdx];
                mean_dD = modelRange.penet[liabIdx][2][penIdx];
                mean_dd = modelRange.penet[liabIdx][3][penIdx];
                SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
                SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
                SD_dD = modelRange.param[liabIdx][2][0][paramIdx];
                SD_dd = modelRange.param[liabIdx][3][0][paramIdx];
                threshold = modelRange.tthresh[liabIdx][thresholdIdx];

                /* check against the hard coded constraint */
                if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
                  constraint =
                    (1 - gfreq) * (1 -
                                   gfreq) * mean_dd *
                    SD_dd + 2 * gfreq * (1 - gfreq) * mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
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
                pTrait->means[liabIdx][1][0] = mean_Dd;
                pTrait->means[liabIdx][1][1] = mean_dd;
                pTrait->stddev[liabIdx][0][0] = SD_DD;
                pTrait->stddev[liabIdx][0][1] = SD_Dd;
                pTrait->stddev[liabIdx][1][0] = SD_Dd;
                pTrait->stddev[liabIdx][1][1] = SD_dd;

                /* threshold for QT */
                pTrait->cutoffValue[liabIdx] = threshold;

              } /* liability class Index */
              if (breakFlag == TRUE)
                continue;
              if (modelOptions.polynomial == TRUE);
              else
                update_penetrance (&pedigreeSet, traitLocus);
              sprintf (partialPolynomialFunctionName, "TL5_P%%sSL%d", modelOptions.sexLinked);
              compute_likelihood (&pedigreeSet);
              cL[5]++;
#ifndef SIMPLEPROGRESS
              if (cL[5] % MAX (1, eCL[5] / 5) == 1) {
                fprintf (stdout, "Trait likelihood evaluations %lu%% complete\r", cL[5] * 100 / eCL[5]);
                fflush (stdout);
              }
#endif
              if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
                fprintf (stderr, "Trait has likelihood 0\n");
                fprintf (stderr, "dgf=%f\n", gfreq);
                for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                  pen_DD = modelRange.penet[liabIdx][0][penIdx];
                  pen_Dd = modelRange.penet[liabIdx][1][penIdx];
                  pen_dD = modelRange.penet[liabIdx][2][penIdx];
                  pen_dd = modelRange.penet[liabIdx][3][penIdx];
		  if (modelOptions.imprintingFlag)
		    fprintf (stderr, "Liab %d penentrance %f %f %f %f\n", liabIdx + 1, pen_DD, pen_Dd, pen_dD, pen_dd);
		  else
		    fprintf (stderr, "Liab %d penentrance %f %f %f\n", liabIdx + 1, pen_DD, pen_Dd, pen_dd);
                }
                exit (EXIT_FAILURE);
              }

              for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                /* save the likelihood at null */
                pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
                likelihoodQT[pedIdx][gfreqInd][penIdx]
                  [paramIdx][thresholdIdx] = pPedigree->likelihood;
              }

              log10_likelihood_null = pedigreeSet.log10Likelihood;
              if (isnan (log10_likelihood_null))
                fprintf (stderr, "Trait likelihood is NAN.\n");
              likelihoodQT[pedigreeSet.numPedigree][gfreqInd][penIdx]
                [paramIdx][thresholdIdx] = log10_likelihood_null;
            }   /* thresholdIdx */
          }     /* penIdx */
        }       /* paramIdx */
      } /* gfreq */
    }   /* end of QT */

#ifndef SIMPLEPROGRESS
    fprintf (stdout, "Trait likelihood evaluations 100%% complete\n");
#endif

    if (modelOptions.polynomial == TRUE) {
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        /* save the likelihood at trait */
        pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        pPedigree->traitLikelihoodPolynomial = pPedigree->likelihoodPolynomial;
        pPedigree->traitLikelihoodPolyList = pPedigree->likelihoodPolyList;
      }
    }

    /* print out some statistics under dry run */
    if (modelOptions.dryRun != 0) {
      print_dryrun_stat (&pedigreeSet, -1);
    }

    /* get the trait locations we need to evaluate at */
    numPositions = modelRange.ntloc;
    mp_result = (SUMMARY_STAT *) calloc (numPositions, sizeof (SUMMARY_STAT));

    writeMPBRFileHeader ();
    writeMPMODFileHeader ();

    prevFirstMarker = -1;
    prevLastMarker = -1;
    prevTraitInd = -1;
    leftMarker = -1;

    /* Iterate over all positions in the analysis. */
    for (posIdx = 0; posIdx < numPositions; posIdx++) {
      /* positions listed are sex average positions */
      traitPos = modelRange.tloc[posIdx];
      /* Set the sex-averaged position first. The sex-specific positions will be updated 
       * once markers are selected since interpolation might be needed. */
      pTraitLocus->mapPosition[0] = traitPos;
      pTraitLocus->mapPosition[1] = traitPos;
      pTraitLocus->mapPosition[2] = traitPos;
      /* initialize the locusList */
      locusList = &savedLocusList;
      memset (locusList->pLocusIndex, 0, sizeof (int) * locusList->maxNumLocus);
      for (k = 0; k < 3; k++) {
        memset (&locusList->pPrevLocusDistance[k][0], 0, sizeof (double) * locusList->maxNumLocus);
        memset (&locusList->pNextLocusDistance[k][0], 0, sizeof (double) * locusList->maxNumLocus);
      }
      locusList->numLocus = 1;
      locusList->pLocusIndex[0] = traitLocus;
      for (k = 0; k < 3; k++) {
        locusList->pPrevLocusDistance[k][0] = -1;
        locusList->pNextLocusDistance[k][0] = -1;
      }
      /* select markers to be used for the multipoint analysis */
      add_markers_to_locuslist (locusList, modelType.numMarkers, &leftMarker, 0, originalLocusList.numLocus - 1, traitPos, 0);
      /* store the markers used */
      mp_result[posIdx].pMarkers = (int *) calloc (modelType.numMarkers, sizeof (int));
      k = 0;    /* marker index */
      for (i = 0; i < locusList->numLocus; i++) {
        j = locusList->pLocusIndex[i];
        if (originalLocusList.ppLocusList[j]->locusType == LOCUS_TYPE_MARKER) {
          mp_result[posIdx].pMarkers[k] = j;
          k++;
        } else {
          mp_result[posIdx].trait = i;
          traitIndex = i;
        }
      }
      locusList->traitLocusIndex = traitIndex;
      locusList->traitOrigLocus = traitLocus;

#ifndef SIMPLEPROGRESS
      /* Say where we're at with the trait locus and markers. */
      fprintf (stdout, "Starting w/trait locus at %.2f (%d/%d positions) with", traitPos, posIdx + 1, numPositions);
#endif

      markerSetChanged = FALSE;
      if (prevFirstMarker != mp_result[posIdx].pMarkers[0] ||
          prevLastMarker != mp_result[posIdx].pMarkers[modelType.numMarkers - 1]) {
        /* marker set has changed */
        markerSetChanged = TRUE;
        markerLocusList.pLocusIndex[0] = mp_result[posIdx].pMarkers[0];
        prevPos = get_map_position (markerLocusList.pLocusIndex[0]);
        /* set up marker set locus list */
        for (k = 1; k < modelType.numMarkers; k++) {
          markerLocusList.pLocusIndex[k] = mp_result[posIdx].pMarkers[k];
          currPos = get_map_position (markerLocusList.pLocusIndex[k]);
          for (j = 0; j < 3; j++) {
            markerLocusList.pPrevLocusDistance[j][k] =
              markerLocusList.pNextLocusDistance[j][k - 1] =
              cm_to_recombination_fraction (currPos[j] - prevPos[j], map.mapFunction);
          }
          if (modelOptions.mapFlag == SA) {
            for (j = 1; j <= 2; j++) {
              markerLocusList.pPrevLocusDistance[j][k] =
                markerLocusList.pNextLocusDistance[j][k - 1] = markerLocusList.pPrevLocusDistance[0][k];
            }
          }
          prevPos = currPos;
        }       /* end of loop over the markers to set up locus list */

#ifndef SIMPLEPROGRESS
        fprintf (stdout, " new markers");
        for (k = 0; k < modelType.numMarkers; k++)
          fprintf (stdout, " %d(%.2f)", markerLocusList.pLocusIndex[k], *get_map_position (markerLocusList.pLocusIndex[k]));
        fprintf (stdout, "\n");
#endif

	locusList = &markerLocusList;
        xmissionMatrix = markerMatrix;
        if (modelOptions.polynomial == TRUE) {
          pedigreeSetPolynomialClearance (&pedigreeSet);
        }

        /* save the polynomial flag */
        polynomialFlag = modelOptions.polynomial;
        status = populate_xmission_matrix (markerMatrix, markerLocusList.numLocus, initialProbAddr, initialProbAddr2,
                                           initialHetProbAddr, 0, -1, -1, 0);
        if (modelOptions.polynomial == TRUE)
          freePolys ();

        print_xmission_matrix (markerMatrix, markerLocusList.numLocus, 0, 0, tmpID);

#ifndef SIMPLEPROGRESS
        /* Calculate likelihood for the marker set */
        fprintf (stdout, "Determining marker set likelihood...\n");
#endif
        for (k = 0; k < modelType.numMarkers; k++) {
          markerNameList[k] = (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[k]])->sName;
        }
        for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
          /* save the marker likelihood   */
          pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
          if (modelOptions.saveResults == TRUE) {
            pPedigree->load_flag =
              restoreMarker (pPedigree->sPedigreeID,
                             (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
                             modelType.numMarkers, markerNameList, &(pPedigree->markerLikelihood));
          } else {
            pPedigree->load_flag = 0;
          }
        }
	if (markerSetChanged) {
	  pushStatus ('k', "buildML6");
	  char markerNo[8];
	  sprintf (partialPolynomialFunctionName, "ML6_P%%sC%dM",
		   (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome);
	  for (k = 0; k < modelType.numMarkers; k++) {
	    sprintf (markerNo, "_%d", markerLocusList.pLocusIndex[k]);
	    strcat (partialPolynomialFunctionName, markerNo);
	  }
	} else
	  pushStatus ('k', "evalML6");
        compute_likelihood (&pedigreeSet);
        cL[6]++;
	popStatus ('k');

#ifndef SIMPLEPROGRESS
        fprintf (stdout, "Marker set likelihood evaluations %lu%% complete...\n",
                 MAX (cL[6] * 100 / eCL[6], (posIdx + 1) * 100 / numPositions));
#endif

        modelOptions.polynomial = polynomialFlag;

        /* print out some statistics under dry run */
        if (modelOptions.dryRun != 0) {
          print_dryrun_stat (&pedigreeSet, -1);
        } else {
          /* save the results for marker likelihood */
          for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
            /* save the likelihood at null */
            pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];

            //      fprintf(stderr, "pedIdx=%d  markerpediLikehood %G\n", pedIdx, pPedigree->likelihood);
            if (modelOptions.saveResults == TRUE) {
              if (pPedigree->load_flag == 0) {  /*save only for the pedigrees which were add for this run */
                pPedigree->markerLikelihood = pPedigree->likelihood;
                pPedigree->load_flag =
                  saveMarker (pPedigree->sPedigreeID,
                              (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
                              modelType.numMarkers, markerNameList, &(pPedigree->markerLikelihood));
              }
            } else {
              pPedigree->markerLikelihood = pPedigree->likelihood;
            }
            pPedigree->load_flag = 0;
          }
          pedigreeSet.log10MarkerLikelihood = pedigreeSet.log10Likelihood;
        }
      } /* end of marker set change */
      else
#ifndef SIMPLEPROGRESS
        fprintf (stdout, " same markers\n");
#else
        ;
#endif
      prevFirstMarker = mp_result[posIdx].pMarkers[0];
      prevLastMarker = mp_result[posIdx].pMarkers[modelType.numMarkers - 1];
      if (markerSetChanged || prevTraitInd != mp_result[posIdx].trait)
        locusListChanged = TRUE;
      else
        locusListChanged = FALSE;
      prevTraitInd = mp_result[posIdx].trait;

      locusList = &savedLocusList;
      xmissionMatrix = altMatrix;
      /* interpolate trait postion for sex specific analysis */
      if (modelOptions.mapFlag == SS) {
        if (traitIndex == 0) {
          marker1Pos = get_map_position (locusList->pLocusIndex[1]);
          /* trait is the first one in the list */
          if (traitPos < ERROR_MARGIN && traitPos > -ERROR_MARGIN) {
            /* trait is at 0 */
            pTraitLocus->mapPosition[MAP_MALE] = pTraitLocus->mapPosition[MAP_FEMALE] = 0;
          } else {
            /* get the relative position on the sex average map */
            relativePos = traitPos / marker1Pos[0];
            pTraitLocus->mapPosition[MAP_MALE] = relativePos * marker1Pos[MAP_MALE];
            pTraitLocus->mapPosition[MAP_FEMALE] = relativePos * marker1Pos[MAP_FEMALE];
          }
          /* update the inter locus distance - sex averaged already done before */
          for (k = 1; k < 3; k++) {
            locusList->pNextLocusDistance[k][0] =
              locusList->pPrevLocusDistance[k][1] =
              cm_to_recombination_fraction (marker1Pos[k] - pTraitLocus->mapPosition[k], map.mapFunction);
          }
        } else if (traitIndex == modelType.numMarkers) {
          /* trait is the last one in the list */
          marker1Pos = get_map_position (locusList->pLocusIndex[modelType.numMarkers - 2]);
          marker2Pos = get_map_position (locusList->pLocusIndex[modelType.numMarkers - 1]);
          /* get the relative position on the sex average map */
          dist = marker2Pos[0] - marker1Pos[0];
          if (dist > ERROR_MARGIN) {
            relativePos = (traitPos - marker2Pos[0]) / dist;
            pTraitLocus->mapPosition[MAP_MALE] =
              relativePos * (marker2Pos[MAP_MALE] - marker1Pos[MAP_MALE]) + marker2Pos[MAP_MALE];
            pTraitLocus->mapPosition[MAP_FEMALE] =
              relativePos * (marker2Pos[MAP_FEMALE] - marker1Pos[MAP_FEMALE]) + marker2Pos[MAP_FEMALE];
          } else {
            pTraitLocus->mapPosition[MAP_MALE] = traitPos - marker2Pos[0] + marker2Pos[MAP_MALE];
            pTraitLocus->mapPosition[MAP_FEMALE] = traitPos - marker2Pos[0] + marker2Pos[MAP_FEMALE];
          }

          /* update the inter locus distance - sex averaged already done before */
          for (k = 1; k <= 2; k++) {
            locusList->pNextLocusDistance[k][traitIndex - 1] =
              locusList->pPrevLocusDistance[k][traitIndex] =
              cm_to_recombination_fraction (pTraitLocus->mapPosition[k] - marker2Pos[k], map.mapFunction);
          }

        } else {
          /* trait is in between two markers */
          marker1Pos = get_map_position (locusList->pLocusIndex[traitIndex - 1]);
          marker2Pos = get_map_position (locusList->pLocusIndex[traitIndex + 1]);
          /* get the relative position on the sex average map */
          dist = marker2Pos[0] - marker1Pos[0];
          if (dist > ERROR_MARGIN) {
            relativePos = (traitPos - marker1Pos[0]) / dist;
            pTraitLocus->mapPosition[MAP_MALE] =
              relativePos * (marker2Pos[MAP_MALE] - marker1Pos[MAP_MALE]) + marker1Pos[MAP_MALE];
            pTraitLocus->mapPosition[MAP_FEMALE] =
              relativePos * (marker2Pos[MAP_FEMALE] - marker1Pos[MAP_FEMALE]) + marker1Pos[MAP_FEMALE];
          } else {
            pTraitLocus->mapPosition[MAP_MALE] = marker1Pos[MAP_MALE];
            pTraitLocus->mapPosition[MAP_FEMALE] = marker1Pos[MAP_FEMALE];
          }
          /* update the inter locus distance - sex averaged already done before */
          for (k = 1; k < 3; k++) {
            locusList->pNextLocusDistance[k][traitIndex - 1] =
              locusList->pPrevLocusDistance[k][traitIndex] =
              cm_to_recombination_fraction (pTraitLocus->mapPosition[k] - marker1Pos[k], map.mapFunction);
            locusList->pNextLocusDistance[k][traitIndex] =
              locusList->pPrevLocusDistance[k][traitIndex + 1] =
              cm_to_recombination_fraction (marker2Pos[k] - pTraitLocus->mapPosition[k], map.mapFunction);
          }
        }
      }

      /* the locus list has been built, go on to the analysis 
       * multipoint DT */
      if (markerSetChanged || locusListChanged) {
        if (modelOptions.polynomial == TRUE) {
          pedigreeSetPolynomialClearance (&pedigreeSet);
          status =
            populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr,
                                      0, -1, -1, 0);
          print_xmission_matrix (altMatrix, savedLocusList.numLocus, 0, 0, tmpID);
          if (modelOptions.polynomial == TRUE)
            freePolys ();
        }
      }

      if (modelOptions.polynomial != TRUE)
        status =
          populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr, initialProbAddr2, initialHetProbAddr, 0, -1, -1, 0);
      /* For alternative */
#ifndef SIMPLEPROGRESS
      fprintf (stdout, "Determining combined likelihood...\n");
#endif

      if (pTrait->type == DICHOTOMOUS) {
        for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
          pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
          /* load stored alternative likelihood if they were already stored */
          if (modelOptions.saveResults == TRUE)
            pPedigree->load_flag =
              restoreAlternative (pPedigree->sPedigreeID,
                                  (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
                                  traitPos, pPedigree->alternativeLikelihoodDT);
          else
            pPedigree->load_flag = 0;
        }

        for (penIdx = 0; (penIdx == 0) || (modelOptions.dryRun == 0 && penIdx < modelRange.npenet); penIdx++) {
          for (liabIdx = 0; (liabIdx == 0) || (modelOptions.dryRun == 0 && liabIdx < modelRange.nlclass); liabIdx++) {
            pen_DD = modelRange.penet[liabIdx][0][penIdx];
            pen_Dd = modelRange.penet[liabIdx][1][penIdx];
            pen_dD = modelRange.penet[liabIdx][2][penIdx];
            pen_dd = modelRange.penet[liabIdx][3][penIdx];
            pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
            pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
            pTrait->penetrance[2][liabIdx][1][0] = pen_dD;
            pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
            pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
            pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
            pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_dD;
            pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;
          }

          if (modelOptions.polynomial != TRUE)
            update_penetrance (&pedigreeSet, traitLocus);       // Only need to update trait locus

          /* Iterate over gene frequencies -- just one loop for dry-runs. */
          for (gfreqInd = 0; (gfreqInd == 0) || (modelOptions.dryRun == 0 && gfreqInd < modelRange.ngfreq); gfreqInd++) {

            /* Updated trait locus allele frequencies */
            gfreq = modelRange.gfreq[gfreqInd];
            pLocus->pAlleleFrequency[0] = gfreq;
            pLocus->pAlleleFrequency[1] = 1 - gfreq;

            if (modelOptions.polynomial != TRUE)
              update_locus (&pedigreeSet, traitLocus);

            /* If we're not on the first iteration, it's not a polynomial build, so
             * show progress at 1 minute intervals. Have a care to avoid division by zero. */
	    char markerNo[8];
	    sprintf (partialPolynomialFunctionName, "CL7_P%%sC%dM",
		     (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome);
	    for (k = 0; k < modelType.numMarkers; k++) {
	      if (traitPos <= *get_map_position (markerLocusList.pLocusIndex[k]) &&
		  (strstr (partialPolynomialFunctionName, "_T") == NULL))
		strcat (partialPolynomialFunctionName, "_T");
	      sprintf (markerNo, "_%d", markerLocusList.pLocusIndex[k]);
	      strcat (partialPolynomialFunctionName, markerNo);
	    }
	    if (strstr (partialPolynomialFunctionName, "_T") == NULL)
	      strcat (partialPolynomialFunctionName, "_T");
            if (gfreqInd != 0 || penIdx != 0) {
	      pushStatus ('k', "evalCL7");
              swStart (combinedComputeSW);
              compute_likelihood (&pedigreeSet);
              cL[7]++;
              swStop (combinedComputeSW);
              if (statusRequestSignal) {
                statusRequestSignal = FALSE;
                if (cL[7] > 1) {        // The first time thru we have no basis for estimation
#ifndef SIMPLEPROGRESS
                  fprintf (stdout, "%s %lu%% complete (~%lu min left)\r",
                           "Combined likelihood evaluations", cL[7] * 100 / eCL[7],
                           ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                            eCL[7] / cL[7] - (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
#else
                  fprintf (stdout, "%s %lu%% complete (~%lu min left)\r",
                           "Calculations", (cL[6] + cL[7]) * 100 / (eCL[6] + eCL[7]),
                           ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                            (eCL[6] + eCL[7]) / (cL[6] + cL[7]) -
                            (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
#endif
                  fflush (stdout);
                }
              }
	      popStatus ('k');
            } else {     // This _is_ the first iteration
	      pushStatus ('k', "buildCL7");
              swStart (combinedBuildSW);
              compute_likelihood (&pedigreeSet);
              cL[7]++;
              swStop (combinedBuildSW);
#ifndef SIMPLEPROGRESS
              fprintf (stdout, "%s %lu%% complete\r", "Combined likelihood evaluations", cL[7] * 100 / eCL[7]);
#else
              fprintf (stdout, "%s %lu%% complete\r", "Calculations", (cL[6] + cL[7]) * 100 / (eCL[6] + eCL[7]));
#endif
              fflush (stdout);
	      popStatus ('k');
            }
            /* print out some statistics under dry run */
            if (modelOptions.dryRun != 0) {
              print_dryrun_stat (&pedigreeSet, traitPos);
            } else {

              log10_likelihood_alternative = pedigreeSet.log10Likelihood;
              if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99)
                log10_likelihood_ratio = 0;
              else
                log10_likelihood_ratio =
                  log10_likelihood_alternative - likelihoodDT[gfreqInd][penIdx] - pedigreeSet.log10MarkerLikelihood;
              /* check for overflow problem !!! */
              if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
                likelihood_ratio = DBL_MAX;
                mp_result[posIdx].lr_total += DBL_MAX;
              } else
                /* check for underflow problem too !!! */
              if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
                likelihood_ratio = 0;
              } else {
                likelihood_ratio = pow (10.0, log10_likelihood_ratio);
                mp_result[posIdx].lr_total += likelihood_ratio;
              }
              /* add the result to the right placeholder */
              mp_result[posIdx].lr_count++;
              for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                if (pPedigree->load_flag == 0) {
                  pPedigree->alternativeLikelihoodDT[gfreqInd]
                    [penIdx] = pPedigree->likelihood;
                }
              }
              /* caculating the Het */
              for (j = 0; j < modelRange.nalpha; j++) {
                alphaV = modelRange.alpha[j];
                alphaV2 = 1 - alphaV;
                if (alphaV2 < 0)
                  alphaV2 = 0;
                log10HetLR = 0;
                for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                  pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                  homoLR = pPedigree->alternativeLikelihoodDT[gfreqInd]
                    [penIdx] / (pPedigree->traitLikelihoodDT[gfreqInd][penIdx] * pPedigree->markerLikelihood);
                  /*              if (homoLR > 1.0e40 || homoLR < 1.0e-40) {
                   * fprintf(stderr, "homoLR %G, alt %G, trait %G, mrk %G\n",
                   * homoLR, pPedigree->alternativeLikelihoodDT[gfreqInd][penIdx],
                   * pPedigree->traitLikelihoodDT[gfreqInd][penIdx],
                   * pPedigree->markerLikelihood);
                   * } */
                  if (alphaV * homoLR + alphaV2 < 0)
                    fprintf (stderr, "HET LR less than 0. Check!!!\n");
                  log10HetLR += log10 (alphaV * homoLR + alphaV2);
                  // if (log10HetLR > 10 || log10HetLR < -40) {
                  /*if(gfreqInd ==0 && j==0){
                   * fprintf(stderr, "gf=%d pen=%d log10HetLR %G, homoLR %G, alt %G, trait %G, mrk %G\n",
                   * gfreqInd, penIdx,log10HetLR,
                   * homoLR, pPedigree->alternativeLikelihoodDT[gfreqInd][penIdx],
                   * pPedigree->traitLikelihoodDT[gfreqInd][penIdx],
                   * pPedigree->markerLikelihood);
                   * //  exit(0);
                   * } */
                }
                if (log10HetLR >= DBL_MAX_10_EXP - 1) {
                  hetLR = DBL_MAX;
                  mp_result[posIdx].het_lr_total = DBL_MAX;
                } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
                  hetLR = 0;
                } else {
                  hetLR = pow (10, log10HetLR);
                  mp_result[posIdx].het_lr_total += hetLR;
                }
                if (mp_result[posIdx].max_penIdx < 0 || hetLR > mp_result[posIdx].max_lr) {
                  mp_result[posIdx].max_lr = hetLR;
                  mp_result[posIdx].max_alpha = alphaV;
                  mp_result[posIdx].max_gfreq = gfreq;
                  mp_result[posIdx].max_penIdx = penIdx;
                }
              } /* end of calculating HET LR */
            }
          }     /* end of genFreq loop */
        }


        /* end of penetrance loop */
        /* save the alternative likelihood */
        for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
          pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
          if ((modelOptions.saveResults == TRUE) && (pPedigree->load_flag == 0)) {      /*save only for the pedigrees which were add for this run */
            pPedigree->load_flag =
              saveAlternative (pPedigree->sPedigreeID,
                               (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome, traitPos,
                               pPedigree->alternativeLikelihoodDT);
          }
          pPedigree->load_flag = 0;
        }

#ifndef SIMPLEPROGRESS
        fprintf (stdout, "%s %lu%% complete (~%lu min left)\n",
                 "Combined likelihood evaluations", cL[7] * 100 / eCL[7],
                 (combinedComputeSW->swAccumWallTime * eCL[7] / cL[7] - combinedComputeSW->swAccumWallTime) / 60);
#else
        fprintf (stdout, "%s %lu%% complete (~%lu min left)\r",
                 "Calculations", (cL[6] + cL[7]) * 100 / (eCL[6] + eCL[7]),
                 (combinedComputeSW->swAccumWallTime * (eCL[6] + eCL[7]) / (cL[6] + cL[7]) -
                  combinedComputeSW->swAccumWallTime) / 60);
#endif

      } /* end of TP */
      else
        /* multipoint QT or COMBINED */
      {
        for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
          gfreq = modelRange.gfreq[gfreqInd];
          pLocus->pAlleleFrequency[0] = gfreq;
          pLocus->pAlleleFrequency[1] = 1 - gfreq;

          update_locus (&pedigreeSet, traitLocus);
          /* this should be MEAN + SD */
          for (paramIdx = 0; paramIdx < modelRange.nparam; paramIdx++) {
            for (penIdx = 0; penIdx < modelRange.npenet; penIdx++) {
              breakFlag = FALSE;
              for (thresholdIdx = 0; thresholdIdx < modelRange.ntthresh; thresholdIdx++) {
                for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
                  mean_DD = modelRange.penet[liabIdx][0][penIdx];
                  mean_Dd = modelRange.penet[liabIdx][1][penIdx];
                  mean_dD = modelRange.penet[liabIdx][2][penIdx];
                  mean_dd = modelRange.penet[liabIdx][3][penIdx];
                  SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
                  SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
                  SD_dD = modelRange.param[liabIdx][2][0][paramIdx];
                  SD_dd = modelRange.param[liabIdx][3][0][paramIdx];
                  threshold = modelRange.tthresh[liabIdx][thresholdIdx];

                  if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
                    /* check against the hard coded constraint */
                    constraint =
                      (1 - gfreq) * (1 -
                                     gfreq) * mean_dd *
                      SD_dd + 2 * gfreq * (1 - gfreq) * mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
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
                  pTrait->means[liabIdx][1][0] = mean_Dd;
                  pTrait->means[liabIdx][1][1] = mean_dd;
                  pTrait->stddev[liabIdx][0][0] = SD_DD;
                  pTrait->stddev[liabIdx][0][1] = SD_Dd;
                  pTrait->stddev[liabIdx][1][0] = SD_Dd;
                  pTrait->stddev[liabIdx][1][1] = SD_dd;

                  /* threshold for QT */
                  pTrait->cutoffValue[liabIdx] = threshold;

                }       /* liability class Index */
                if (breakFlag == TRUE)
                  continue;
                if (modelOptions.polynomial == TRUE);
                else
                  update_penetrance (&pedigreeSet, traitLocus);
                /* ready for the alternative hypothesis */
                locusList = &savedLocusList;
                xmissionMatrix = altMatrix;
                if (modelOptions.polynomial == TRUE);
                else
                  status =
                    populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, initialProbAddr2,
                                              initialHetProbAddr, 0, -1, -1, 0);

                /* If we're not on the first iteration, it's not a polynomial build, so
                 * show progress at 1 minute intervals. Have a care to avoid division by zero. */
		char markerNo[8];
		sprintf (partialPolynomialFunctionName, "CL8_P%%sC%dM",
			 (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome);
		for (k = 0; k < modelType.numMarkers; k++) {
		  if (traitPos <= *get_map_position (markerLocusList.pLocusIndex[k]) &&
		      (strstr (partialPolynomialFunctionName, "_T") == NULL))
		    strcat (partialPolynomialFunctionName, "_T");
		  sprintf (markerNo, "_%d", markerLocusList.pLocusIndex[k]);
		  strcat (partialPolynomialFunctionName, markerNo);
		}
		if (strstr (partialPolynomialFunctionName, "_T") == NULL)
		  strcat (partialPolynomialFunctionName, "_T");
                if (gfreqInd != 0 || paramIdx != 0 || penIdx != 0) {
		  pushStatus ('k', "evalCL8");
                  swStart (combinedComputeSW);
                  compute_likelihood (&pedigreeSet);
                  cL[8]++;
                  swStop (combinedComputeSW);
                  if (statusRequestSignal) {
                    statusRequestSignal = FALSE;
                    if (cL[8] > 1) {    // The first time thru we have no basis for estimation
#ifndef SIMPLEPROGRESS
                      fprintf (stdout, "%s %lu%% complete (~%lu min left)\r",
                               "Combined likelihood evaluations", cL[8] * 100 / eCL[8],
                               ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                                eCL[8] / cL[8] - (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
#else
                      fprintf (stdout, "%s %lu%% complete (~%lu min left)\r",
                               "Calculations", (cL[6] + cL[8]) * 100 / (eCL[6] + eCL[8]),
                               ((combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime) *
                                (eCL[6] + eCL[8]) / (cL[6] + cL[8]) -
                                (combinedComputeSW->swAccumWallTime + combinedBuildSW->swAccumWallTime)) / 60);
#endif
                      fflush (stdout);
                    }
                  }
		  popStatus ('k');
                } else {  // This _is_ the first iteration
		  pushStatus ('k', "buildCL8");
                  swStart (combinedBuildSW);
                  compute_likelihood (&pedigreeSet);
                  cL[8]++;
                  swStop (combinedBuildSW);
#ifndef SIMPLEPROGRESS
                  fprintf (stdout, "%s %lu%% complete at %d\r", "Combined likelihood evaluations", cL[8] * 100 / eCL[8], nodeId);
#else
                  fprintf (stdout, "%s %lu%% complete\r", "Calculations", (cL[6] + cL[8]) * 100 / (eCL[6] + eCL[8]));
#endif
                  fflush (stdout);
		  popStatus ('k');
                }
                log10_likelihood_alternative = pedigreeSet.log10Likelihood;
                if (isnan (log10_likelihood_alternative))
                  fprintf (stderr, "ALT likelihood is NAN.\n");
                if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
                  log10_likelihood_ratio = 0;
                } else {
                  log10_likelihood_ratio = log10_likelihood_alternative - likelihoodQT[pedigreeSet.numPedigree][gfreqInd]
                    [penIdx][paramIdx][thresholdIdx] - pedigreeSet.log10MarkerLikelihood;
                }
                /* check for overflow problem !!! */
                if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
                  likelihood_ratio = DBL_MAX;
                  mp_result[posIdx].lr_total += DBL_MAX;
                } else
                  /* check for underflow problem too !!! */
                if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
                  likelihood_ratio = 0;
                } else {
                  likelihood_ratio = pow (10.0, log10_likelihood_ratio);
                  mp_result[posIdx].lr_total += likelihood_ratio;
                }
                /* add the result to the right placeholder */
                mp_result[posIdx].lr_count++;

                if (isnan (likelihood_ratio))
                  fprintf (stderr, "LR for the pedigree set is NAN.\n");
                /* caculating the HET */
                for (j = 0; j < modelRange.nalpha; j++) {
                  alphaV = modelRange.alpha[j];
                  alphaV2 = 1 - alphaV;
                  if (alphaV2 < 0)
                    alphaV2 = 0;
                  log10HetLR = 0;
                  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
                    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
                    homoLR =
                      pPedigree->likelihood / (likelihoodQT[pedIdx][gfreqInd][penIdx][paramIdx][thresholdIdx] *
                                               pPedigree->markerLikelihood);
                    log10HetLR += log10 (alphaV * homoLR + alphaV2);
                  }
                  if (log10HetLR >= DBL_MAX_10_EXP - 1) {
                    hetLR = DBL_MAX;
                    mp_result[posIdx].het_lr_total = DBL_MAX;
                  } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
                    hetLR = 0;
                  } else {
                    adjustedHetLR = hetLR = pow (10, log10HetLR);
                    /* for threshold parameter, we need to make sure the weighting is even */
                    if (1 || modelType.distrib == QT_FUNCTION_CHI_SQUARE) {
                      if (modelRange.ntthresh == 1) {
                        adjustedHetLR *= 2 * (modelType.maxThreshold - modelType.minThreshold);
                      } else if (thresholdIdx == modelRange.ntthresh - 1) {
                        adjustedHetLR *= (2 * modelType.maxThreshold - threshold - modelRange.tthresh[0][thresholdIdx - 1]);
                      } else if (thresholdIdx == 0) {
                        adjustedHetLR *= (threshold + modelRange.tthresh[0][thresholdIdx + 1] - 2 * modelType.minThreshold);
                      } else
                        adjustedHetLR *= modelRange.tthresh[0][thresholdIdx + 1] - modelRange.tthresh[0][thresholdIdx - 1];
                    }
                    mp_result[posIdx].het_lr_total += adjustedHetLR;
                  }
                  if (mp_result[posIdx].max_penIdx < 0 || hetLR > mp_result[posIdx].max_lr) {
                    mp_result[posIdx].max_lr = hetLR;
                    mp_result[posIdx].max_alpha = alphaV;
                    mp_result[posIdx].max_gfreq = gfreq;
                    mp_result[posIdx].max_penIdx = penIdx;
                    mp_result[posIdx].max_paramIdx = paramIdx;
                    mp_result[posIdx].max_thresholdIdx = thresholdIdx;
                  }
                }
              } /* end of threshold loop */
            }   /* end of penetrance loop */
          }     /* end of parameter loop */
        }       /* end of gene freq */

#ifndef SIMPLEPROGRESS
        fprintf (stdout, "%s %lu%% complete (~%lu min left)\n",
                 "Combined likelihood evaluations", cL[8] * 100 / eCL[8],
                 (combinedComputeSW->swAccumWallTime * eCL[8] / cL[8] - combinedComputeSW->swAccumWallTime) / 60);
#else
        fprintf (stdout, "%s %lu%% complete (~%lu min left)\r",
                 "Calculations", (cL[6] + cL[8]) * 100 / (eCL[6] + eCL[8]),
                 (combinedComputeSW->swAccumWallTime * (eCL[6] + eCL[8]) / (cL[6] + cL[8]) -
                  combinedComputeSW->swAccumWallTime) / 60);
#endif

      } /* end of QT */

      /* Print out average and log10(max) and maximizing parameters */
      if (modelType.trait == DT)
        avgLR = mp_result[posIdx].het_lr_total / (modelRange.nalpha * mp_result[posIdx].lr_count);
      else
        /* under QT CHISQ, threshold parameter has been evenly weighted */
        avgLR =
          mp_result[posIdx].het_lr_total / (modelRange.nalpha * (mp_result[posIdx].lr_count / modelRange.ntthresh) * 2 *
                                            (modelType.maxThreshold - modelType.minThreshold));

      if (avgLR > 0.214)
        ppl = (avgLR * avgLR) / (-5.77 + 54 * avgLR + avgLR * avgLR);
      else
        ppl = 0;

      writeMPBRFileDetail ();
      writeMPMODFileDetail ();

    }   /* end of walking down the chromosome */
  }     /* end of multipoint */
  fprintf (stdout, "\n");

  /* only for multipoint - deallocate memory  */
  if (modelType.type == MP) {
    /* allocate space to save temporary results */
    if (modelType.trait == DT) {

      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
        free (likelihoodDT[gfreqInd]);
      }
      free (likelihoodDT);
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];

        for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {
          free (pPedigree->traitLikelihoodDT[gfreqInd]);
          free (pPedigree->alternativeLikelihoodDT[gfreqInd]);
        }
        free (pPedigree->traitLikelihoodDT);
        free (pPedigree->alternativeLikelihoodDT);
      }
      free (markerNameList);
    }
  }
//  dumpTrackingStats(cL, eCL);
