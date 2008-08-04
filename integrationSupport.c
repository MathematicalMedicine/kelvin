
int
kelvin_dcuhre_integrate (double *integral, double *abserr, double vol_region)
{
  /* Local variables */
  //double a[15], b[15];
  int dim, return_val,i;
  //  dcuhre_state init_state;

  //extern /* Subroutine */ int ftest_();  

  localmax_value = 0.0;



  if (modelType.trait == DICHOTOMOUS) {

    dim = 1 + 3 * modelRange.nlclass;
    s = &init_state;
    initialize_state (s, xl, xu, dim);
    s->verbose = 1;
    s->nlclass = modelRange.nlclass;

    if (modelType.type == TP) {
      s->funsub = (U_fp) compute_hlod_2p_dt;
      s->mType = TP_DT;
    } else {
      s->funsub = (U_fp) compute_hlod_mp_dt;
      s->mType = MP_DT;
    }
  } else {			/*  QT or combined */

    dim = total_dim - 1;	// alpha
    if (modelType.type == TP) {
      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
	dim -= 2;		// theta and dprime
      } else {
	dim -= 1;		//theta
      }
    }

    s = &init_state;
    initialize_state (s, xl, xu, dim);
    s->verbose = 1;
    s->nlclass = modelRange.nlclass;

    s->maxcls = 20000;
   
    if (modelType.type == TP) {
      s->funsub = (U_fp) compute_hlod_2p_qt;
      s->mType = TP_DT;
    } else {
      s->funsub = (U_fp) compute_hlod_mp_qt;
      s->mType = MP_DT;
    }
  }

  for(i=0; i<s->nlclass; i++){
    s->vol_rate /= 6.0;
  }
  s->vol_rate *= vol_region;   /*This is the rate to convert to average function value*/

  fprintf (stderr, "Starting DCUHRE with dim=%d\n", dim);
  return_val = dcuhre_ (s);
  if (return_val > 0 && return_val < 20) {
    fprintf (stderr, "Ending program with error! ifail =%d \n", s->ifail);
  }

  fprintf (stderr,
	   "Final result =%15.10f  with error =%15.10f and neval = %d\n",
	   s->result, s->error, s->total_neval);
  fprintf (stderr, "End of DCUHRE with ifail =%d\n", s->ifail);


  *integral = s->result;
  *abserr = s->error;
  return return_val;

}


void
compute_hlod_mp_qt (double x[], double *f)
{

  int k, j;
  int pedIdx, liabIdx = 0, status;
  double constraint = 0.0;
  double mean_DD;
  double mean_Dd;
  double mean_dd;
  double SD_DD = 0.0;
  double SD_Dd = 0.0;
  double SD_dd = 0.0;
  double gfreq;
  double alphaV;

//  double  theta;
  double threshold = 0.0;
  double log10_likelihood_null, log10_likelihood_alternative,
    log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;
  double log10Likelihood;

  /* for null likelihood calculation */
  double product_likelihood = 1;	/* product of the likelihoods for all the pedigrees */
  double sum_log_likelihood = 0;	/* sum of the log10(likelihood) for all the pedigrees */

  Pedigree *pPedigree;

  //  checkpt ();

  int origLocus = locusList->pLocusIndex[0];

  if (locusList->numLocus > 1)
    origLocus = locusList->pLocusIndex[1];

  gfreq = x[0];
  pLocus->pAlleleFrequency[0] = gfreq;
  pLocus->pAlleleFrequency[1] = 1 - gfreq;

  update_locus (&pedigreeSet, traitLocus);

  //  checkpt ();
  if (print_point_flag)
    fprintf (fphlod, "dp=%f th=%f gf=%f ", fixed_dprime, fixed_theta, gfreq);


  j = 1;			// j=0 for gfrequency

  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    mean_DD = x[j];
    mean_Dd =
      (x[j + 1] - xl[j + 1]) * (x[j] - xl[j]) / (xu[j] - xl[j]) + xl[j + 1];
    mean_dd =
      (x[j + 2] - xl[j + 2]) * (x[j + 1] - xl[j + 1]) / (xu[j + 1] -
							 xl[j + 1]) * (x[j] -
								       xl[j])
      / (xu[j] - xl[j]) + xl[j + 2];
    if (print_point_flag)
      fprintf (fphlod, "pe %f %f %f ", mean_DD, mean_Dd, mean_dd);
    j += 3;
    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      SD_DD = x[j];		//      modelRange.param[liabIdx][0][0][paramIdx];
      SD_Dd = x[j + 1];		//        modelRange.param[liabIdx][1][0][paramIdx];
      SD_dd = x[j + 2];		//        modelRange.param[liabIdx][2][0][paramIdx];
      if (print_point_flag)
	fprintf (fphlod, "sd %f %f %f ", SD_DD, SD_Dd, SD_dd);
      j += 3;
    }
    /* threshold for QT */
    if (modelType.trait == CT) {
      threshold = x[j];		// modelRange.tthresh[liabIdx][thresholdIdx];
      if (print_point_flag)
	fprintf (fphlod, "thd=%f ", threshold);
      j++;
    }

    /* check against the hard coded constraint */
    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      constraint =
	(1.0 - gfreq) * (1.0 - gfreq) * mean_dd * SD_dd + 2 * gfreq * (1.0 -
								       gfreq)
	* mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
      /* fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
         constraint, gfreq, mean_DD, SD_DD, mean_Dd, SD_DD, mean_dd, SD_dd);
       */
      if (constraint >= 3.0 || constraint <= -3.0) {
	// fprintf(stderr,"Constraint is %f \n", constraint);
	num_out_constraint++;
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

  }				/* liability class Index */
  //  checkpt ();
  if (modelOptions.polynomial == TRUE);
  else
    /* only need to update trait locus */
    update_penetrance (&pedigreeSet, traitLocus);
  //  update_penetrance (&pedigreeSet, traitLocus);

  if (print_point_flag)
    fprintf (fphlod, "\n");

  /*This is a temporary checking. */
  if (j != s->ndim) {
    fprintf (stderr, "j=%d  while dim for BR is %d\n", j, s->ndim);
    exit (0);
  }

  /* for trait likelihood */
  locusList = &traitLocusList;
  xmissionMatrix = traitMatrix;

  //fprintf(stderr,"Null likelihood computation is starting for %d pedigrees \n",pedigreeSet.numPedigree);
  /* compute the null likelihood with   */
  pedigreeSet.likelihood = 1;
  pedigreeSet.log10Likelihood = 0;


  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {

    /* save the likelihood at null */
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];

    //fprintf(stderr,"pedIdx =%d ", pedIdx); 
    if (modelOptions.polynomial == TRUE) {
      KASSERT (pPedigree->traitLikelihoodPolynomial != NULL, "Error in  \n");
      /* evaluate likelihood */
      // fprintf (stderr, "evaluaing poly %d pedigree\n", pedIdx);

      evaluatePoly (pPedigree->traitLikelihoodPolynomial,
		    pPedigree->traitLikelihoodPolyList,
		    &pPedigree->likelihood);

      // fprintf (stderr, " is done %f with %d pedigrees\n",
      //               pPedigree->likelihood, pedigreeSet.numPedigree);
    } else {
      initialize_multi_locus_genotype (pPedigree);
      status = compute_pedigree_likelihood (pPedigree);
    }

    /*pPedigree->likelihood is now computed and now check it */
    if (pPedigree->likelihood == 0.0) {
      KLOG (LOGLIKELIHOOD, LOGWARNING,
	    "Pedigree %s has likelihood of 0 or too small.\n",
	    pPedigree->sPedigreeID);
      fprintf (stderr, "Pedigree %s has likelihood of 0 or too small.\n",
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
      sum_log_likelihood += log10Likelihood;
    }
    pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
    // fprintf (stderr,
    //     "null likelihood pedIdx=%d is done %20.18f with product =%20.16f\n",
    //     pedIdx, pPedigree->likelihood, product_likelihood);
  }

  pedigreeSet.likelihood = product_likelihood;
  pedigreeSet.log10Likelihood = sum_log_likelihood;
  log10_likelihood_null = pedigreeSet.log10Likelihood;
  KLOG (LOGLIKELIHOOD, LOGDEBUG, "Sum of log Likelihood is: %e\n",
	sum_log_likelihood);

  //  checkpt ();

  /* This is for alternative likelihood */
  locusList = &savedLocusList;
  xmissionMatrix = altMatrix;
  if (modelOptions.polynomial == TRUE);
  else
    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */

  KLOG (LOGLIKELIHOOD, LOGDEBUG, "Likelihood\n");
  compute_likelihood (&pedigreeSet);

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;
  if (isnan (log10_likelihood_alternative))
    fprintf (stderr, "ALT likelihood is NAN.\n");
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
  } else {
    log10_likelihood_ratio =
      log10_likelihood_alternative - log10_likelihood_null -
      pedigreeSet.log10MarkerLikelihood;
  }
  /* check for overflow problem !!! */
  if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
    likelihood_ratio = DBL_MAX;
  } else {
    /* check for underflow problem too !!! */
    if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
      likelihood_ratio = 0;
    } else {
      likelihood_ratio = pow (10.0, log10_likelihood_ratio);
      //mp_result[posIdx].lr_total += likelihood_ratio;
    }
  }
  if (isnan (likelihood_ratio))
    fprintf (stderr, "LR for the pedigree set is NAN.\n");

  //checkpt ();
  /* caculating the HET */
  for (j = 0; j < 6; j++) {
  // for (j = 0; j < 1; j++) {
    alphaV=alpha[j][0];
    alphaV2 = 1 - alphaV;
    if (alphaV2 < 0)
      alphaV2 = 0;

    log10HetLR = 0;
    for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
      homoLR =
	pPedigree->likelihood / (pedigreeSet.nullLikelihood[pedIdx] *
				 pPedigree->markerLikelihood);
      //fprintf(stderr,"j=%d pedIdx=%d  %20.18f %20.16f %20.16f %20.16f \n",j, pedIdx,pPedigree->likelihood,pedigreeSet.nullLikelihood[pedIdx] * pPedigree->markerLikelihood, homoLR ,log10HetLR);
      if (alphaV * homoLR + alphaV2 < 0)
	fprintf (stderr, "HET LR less than 0. Check!!!\n");
      log10HetLR += log10 (alphaV * homoLR + alphaV2);
    }
    if (log10HetLR >= DBL_MAX_10_EXP - 1) {
      hetLR = DBL_MAX;
    } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
      hetLR = 0;
    } else {
      hetLR = pow (10, log10HetLR);
    }
    if(j ==0){
      alpha_integral = hetLR ;//* alpha[j][1];
    }

    if (print_point_flag)
      fprintf (fphlod, "al=%f Hlod=%f\n", alphaV, hetLR);
    /*Update local maximum as necessary */
    if (hetLR > localmax_value) {
      localmax_value = hetLR;
      localmax_x[0] = gfreq;
      localmax_x[1] = alphaV;
      k = 2;
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	localmax_x[k] = x[k - 1];
	localmax_x[k + 1] =
	  (x[k] - xl[k]) * (x[k - 1] - xl[k - 1]) / (xu[k - 1] - xl[k - 1]) +
	  xl[k];
	localmax_x[k + 2] =
	  (x[k + 1] - xl[k + 1]) * (x[k] - xl[k]) / (xu[k] -
						     xl[k]) * (x[k - 1] -
							       xl[k -
								  1]) /
	  (xu[k - 1] - xl[k - 1]) + xl[k + 1];
	k += 3;
	if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	  localmax_x[k] = x[k - 1];
	  localmax_x[k + 1] = x[k];
	  localmax_x[k + 2] = x[k + 1];
	  k += 3;
	}
	/* threshold for QT */
	if (modelType.trait == CT) {
	  localmax_x[k] = x[k - 1];
	  k++;
	}
      }
    }
  }
  //checkpt ();
  avg_hetLR = alpha_integral;

  //  if (constraint >= 3.0 || constraint <= -3.0) {
  //fprintf (stderr, "Constraint is %f with avg hetLR =%f\n", constraint,
  //     avg_hetLR);
  //}


  /* Jacobian */
  k = 1;
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    avg_hetLR *=
      (x[k] - xl[k]) / (xu[k] - xl[k]) * (x[k] - xl[k]) / (xu[k] -
							   xl[k]) * (x[k +
								       1] -
								     xl[k +
									1]) /
      (xu[k + 1] - xl[k + 1]);
    k += 3;
    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      k += 3;
    }
    if (modelType.trait == CT) {
      k++;
    }
  }
  /*This is a temporary checking. */
  if (k != s->ndim) {
    fprintf (stderr, "k=%d  while dim for BR is %d\n", k, s->ndim);
    exit (0);
  }

  f[0] = avg_hetLR;

}






void
compute_hlod_mp_dt (double x[], double *f)
{

  int j;
  int pedIdx, liabIdx, status;

  double pen_DD, pen_Dd, pen_dd, gfreq, alphaV;
  double log10_likelihood_null, log10_likelihood_alternative,
    log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;
  double log10Likelihood;

  /* for null likelihood calculation */
  double product_likelihood = 1;	/* product of the likelihoods for all the pedigrees */
  double sum_log_likelihood = 0;	/* sum of the log10(likelihood) for all the pedigrees */

  // double *cw_center, *cw_hwidth;  /*pointers for centers and hwidths of the currrent working subregion*/

  Pedigree *pPedigree;
  int origLocus = locusList->pLocusIndex[0];

  if (locusList->numLocus > 1)
    origLocus = locusList->pLocusIndex[1];

  //fprintf(stderr,"in compute hlod x %G %G %G %G\n", x[0],x[1],x[2],x[3]);
  gfreq = x[0];

  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    pen_DD = x[3 * liabIdx + 1];
    pen_Dd = x[3 * liabIdx + 2] * x[3 * liabIdx + 1];
    pen_dd = x[3 * liabIdx + 3] * x[3 * liabIdx + 1] * x[3 * liabIdx + 2];
    pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
    pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
    pTrait->penetrance[2][liabIdx][1][0] = pen_Dd;
    pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
    pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
    pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
    pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_Dd;
    pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;

    //fprintf(stderr,"pene   DD=%G Dd=%G dd=%G\n", pen_DD, pen_Dd, pen_dd);

  }



  if (modelOptions.polynomial == TRUE);
  else
    /* only need to update trait locus */
    update_penetrance (&pedigreeSet, traitLocus);

  pLocus->pAlleleFrequency[0] = gfreq;
  pLocus->pAlleleFrequency[1] = 1 - gfreq;


  if (modelOptions.polynomial == TRUE);
  else
    update_locus (&pedigreeSet, traitLocus);

  /* for trait likelihood */
  locusList = &traitLocusList;
  xmissionMatrix = traitMatrix;

  // fprintf(stderr, "Null likelihood computation is starting for %d pedigrees \n",pedigreeSet.numPedigree);

  /* compute the null likelihood with   */
  pedigreeSet.likelihood = 1;
  pedigreeSet.log10Likelihood = 0;


  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {

    /* save the likelihood at null */
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];


    if (modelOptions.polynomial == TRUE) {
      KASSERT (pPedigree->traitLikelihoodPolynomial != NULL, "Error in  \n");
      /* evaluate likelihood */
      //fprintf(stderr, "evaluaing poly ");
      // printAllVariables();
      //expTermPrinting(stderr, pPedigree->traitLikelihoodPolynomial, 1);
      //fprintf(stderr, "\n");

      evaluatePoly (pPedigree->traitLikelihoodPolynomial,
		    pPedigree->traitLikelihoodPolyList,
		    &pPedigree->likelihood);
      //fprintf(stderr, " is done %f with %d pedigrees\n",pPedigree->likelihood, pedigreeSet.numPedigree);

      if (isnan (pPedigree->likelihood)) {
	fprintf (stderr, " likelihood is nan\n");
	exit (1);
      }

    } else {
      initialize_multi_locus_genotype (pPedigree);
      status = compute_pedigree_likelihood (pPedigree);
    }



    /*pPedigree->likelihood is now computed and now check it */
    if (pPedigree->likelihood == 0.0) {
      KLOG (LOGLIKELIHOOD, LOGWARNING,
	    "Pedigree %s has likelihood of 0 or too small.\n",
	    pPedigree->sPedigreeID);
      fprintf (stderr, "Pedigree %s has likelihood of 0 or too small.\n",
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
      sum_log_likelihood += log10Likelihood;
    }
    pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
    //fprintf(stderr,"null likelihood pedIdx=%d is done %20.18f with product =%20.16f\n",pedIdx,pPedigree->likelihood,product_likelihood );
  }


  pedigreeSet.likelihood = product_likelihood;
  pedigreeSet.log10Likelihood = sum_log_likelihood;
  log10_likelihood_null = pedigreeSet.log10Likelihood;
  KLOG (LOGLIKELIHOOD, LOGDEBUG, "Sum of log Likelihood is: %e\n",
	sum_log_likelihood);

  //fprintf(stderr,"Null likelihood = %20.15f\n", pedigreeSet.likelihood);

  /* This is for alternative likelihood */
  locusList = &savedLocusList;
  xmissionMatrix = altMatrix;
  compute_likelihood (&pedigreeSet);

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
  } else {
    log10_likelihood_ratio =
      log10_likelihood_alternative - log10_likelihood_null;
  }
  /* check for overflow problem !!! */
  if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
    likelihood_ratio = DBL_MAX;
  } else {
    /* check for underflow problem too !!! */
    if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
      likelihood_ratio = 0;
    } else {
      likelihood_ratio = pow (10.0, log10_likelihood_ratio);
      //mp_result[posIdx].lr_total += likelihood_ratio;
    }
  }
  /* add the result to the right placeholder */
  //mp_result[posIdx].lr_count++;
  //fprintf(stderr," %f %f %f %f %20.15f  ", gfreq,pen_DD, pen_Dd,pen_dd, log10_likelihood_alternative);
  /* caculating the HET */
  for (j = 0; j < 6; j++) {
    //for (j = 0; j < 1; j++) {
    alphaV = alpha[j][0];
    alphaV2 = 1 - alphaV;
    if (alphaV2 < 0)
      alphaV2 = 0;

    log10HetLR = 0;
    for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
      homoLR =
	pPedigree->likelihood / (pedigreeSet.nullLikelihood[pedIdx] *
				 pPedigree->markerLikelihood);
      //fprintf(stderr,"j=%d pedIdx=%d  %20.18f %20.16f %20.16f %20.16f \n",j, pedIdx,pPedigree->likelihood,pedigreeSet.nullLikelihood[pedIdx] * pPedigree->markerLikelihood, homoLR ,log10HetLR);
      if (alphaV * homoLR + alphaV2 < 0)
	fprintf (stderr, "HET LR less than 0. Check!!!\n");
      log10HetLR += log10 (alphaV * homoLR + alphaV2);
    }
    if (log10HetLR >= DBL_MAX_10_EXP - 1) {
      hetLR = DBL_MAX;
      //mp_result[posIdx].het_lr_total = DBL_MAX;
    } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
      hetLR = 0;
    } else {
      hetLR = pow (10, log10HetLR);
      //mp_result[posIdx].het_lr_total += hetLR;
    }

    //fprintf(stderr,"j=%d het LR=%15.10f ",j, hetLR);

    if( j==0){
      alpha_integral = hetLR; // * alpha[j][1];
    }

    if (hetLR > localmax_value) {
      localmax_value = hetLR;
      localmax_x[0] = gfreq;
      localmax_x[1] = alphaV;

      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	localmax_x[3 * liabIdx + 2] = x[3 * liabIdx + 1];
	localmax_x[3 * liabIdx + 3] = x[3 * liabIdx + 2] * x[3 * liabIdx + 1];
	localmax_x[3 * liabIdx + 4] =
	  x[3 * liabIdx + 3] * x[3 * liabIdx + 1] * x[3 * liabIdx + 2];
      }
    }

  }				/* end of calculating HET LR */




  avg_hetLR = alpha_integral;
  //fprintf(stderr,"avg hetLR =%15.10f \n", avg_hetLR);

  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    avg_hetLR *= x[3 * liabIdx + 1] * x[3 * liabIdx + 1] * x[3 * liabIdx + 2];
  }

  //  fprintf(stderr,"   hetLR =%f   ", avg_hetLR);


  *f = avg_hetLR;

  // return avg_hetLR;

}

void
compute_hlod_2p_qt (double x[], double *f)
{

  int k, j;
  int pedIdx, liabIdx = 0, status;
  double constraint = 0.0;
  double mean_DD;
  double mean_Dd;
  double mean_dd;
  double SD_DD = 0.0;
  double SD_Dd = 0.0;
  double SD_dd = 0.0;
  double gfreq;
  double alphaV;
  double theta;
  double threshold = 0.0;
  double log10_likelihood_null, log10_likelihood_alternative,
    log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;

  Pedigree *pPedigree;

  gfreq = x[0];
  theta = fixed_theta;		//x[5];  

  if (1 && modelOptions.markerAnalysis == FALSE) {
    pLocus->pAlleleFrequency[0] = gfreq;
    pLocus->pAlleleFrequency[1] = 1 - gfreq;

    if (modelOptions.polynomial == TRUE);
    else
      update_locus (&pedigreeSet, loc1);

  }

  /* this should be MEAN + SD */
  if (print_point_flag)
    fprintf (fphlod, "dp=%f th=%f gf=%f ", fixed_dprime, fixed_theta, gfreq);
  j = 1;
  if (modelOptions.markerAnalysis == FALSE) {
    for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
      mean_DD = x[j];
      mean_Dd =
	(x[j + 1] - xl[j + 1]) * (x[j] - xl[j]) / (xu[j] - xl[j]) + xl[j + 1];
      mean_dd =
	(x[j + 2] - xl[j + 2]) * (x[j + 1] - xl[j + 1]) / (xu[j + 1] -
							   xl[j +
							      1]) * (x[j] -
								     xl[j]) /
	(xu[j] - xl[j]) + xl[j + 2];
      if (print_point_flag)
	fprintf (fphlod, "pe %f %f %f ", mean_DD, mean_Dd, mean_dd);
      j += 3;
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	SD_DD = x[j];		//      modelRange.param[liabIdx][0][0][paramIdx];
	SD_Dd = x[j + 1];	//        modelRange.param[liabIdx][1][0][paramIdx];
	SD_dd = x[j + 2];	//        modelRange.param[liabIdx][2][0][paramIdx];
	if (print_point_flag)
	  fprintf (fphlod, "sd %f %f %f ", SD_DD, SD_Dd, SD_dd);
	j += 3;
      }
      /* threshold for QT */
      if (modelType.trait == CT) {
	threshold = x[j];	// modelRange.tthresh[liabIdx][thresholdIdx];
	if (print_point_flag)
	  fprintf (fphlod, "thd=%f ", threshold);
	j++;
      }

      /* check against the hard coded constraint */
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	constraint =
	  (1.0 - gfreq) * (1.0 - gfreq) * mean_dd * SD_dd + 2 * gfreq * (1 -
									 gfreq)
	  * mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
	/* fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
	   constraint, gfreq, mean_DD, SD_DD, mean_Dd, SD_DD, mean_dd, SD_dd);
	 */
	if (constraint >= 3.0 || constraint <= -3.0) {
	  // fprintf(stderr,"Constraint is %f \n", constraint);
	  num_out_constraint++;
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

    }				/* liability class Index */

    if (modelOptions.polynomial == TRUE);
    else
      update_penetrance (&pedigreeSet, traitLocus);


  }				/* marker to marker analysis */
  if (print_point_flag)
    fprintf (fphlod, "\n");

  /*This is a temporary checking. */
  if (j != s->ndim) {
    fprintf (stderr, "j=%d  while dim for BR is %d\n", j, s->ndim);
    exit (0);
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
    /* populate the matrix */
    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */

  KLOG (LOGLIKELIHOOD, LOGDEBUG, "NULL Likelihood\n");
  compute_likelihood (&pedigreeSet);



  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    fprintf (stderr, "Theta 0.5 has likelihood 0\n");
    fprintf (stderr, "dgf=%f\n", gfreq);
    for (j = 1; j < s->ndim; j++) {
      fprintf (stderr, " %f", x[j]);
    }
    fprintf (stderr, "\n");

    exit (-1);
  }

  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
    /* save the likelihood at null */
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
    pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
  }

  log10_likelihood_null = pedigreeSet.log10Likelihood;
  //for (dprimeIdx = 0;dprimeIdx < pLambdaCell->ndprime;dprimeIdx++){
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);
//    if (pLambdaCell->impossibleFlag[dprimeIdx] != 0)
//        continue;
    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
    copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);
  }

  if (modelOptions.mapFlag == SA) {
    for (k = 0; k < 3; k++) {
      locusList->pNextLocusDistance[k][0] = theta;
      locusList->pPrevLocusDistance[k][1] = theta;
    }
  } else {
    printf ("mapflag sould be SA\n");
    exit (-1);
  }

  if (modelOptions.polynomial == TRUE);
  else
    /* populate the matrix */
    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */

  KLOG (LOGLIKELIHOOD, LOGDEBUG, "ALT Likelihood\n");
  compute_likelihood (&pedigreeSet);

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
  } else {
    log10_likelihood_ratio =
      log10_likelihood_alternative - log10_likelihood_null;
  }

  /* check for overflow problem !!! */
  if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
    likelihood_ratio = DBL_MAX;
  } else {
    /* check for underflow problem too !!! */
    if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
      likelihood_ratio = 0;
    } else {
      likelihood_ratio = pow (10.0, log10_likelihood_ratio);
    }
  }

  /* caculating the HET */
  for (j = 1; j < 6; j++) {
    //for (j = 0; j < 1; j++) {
    alphaV = alpha[j][0];
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
    } else if (log10HetLR <= DBL_MIN_10_EXP + 1) {
      hetLR = 0;
    } else {
      hetLR = pow (10, log10HetLR);
    }
    
      alpha_integral += hetLR * alpha[j][1];
    

    if (print_point_flag)
      fprintf (fphlod, "al=%f Hlod=%f\n", alphaV, hetLR);
    /*Update local maximum as necessary */
    if (hetLR > localmax_value) {
      localmax_value = hetLR;
      localmax_x[0] = gfreq;
      localmax_x[1] = alphaV;
      k = 2;
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	localmax_x[k] = x[k - 1];
	localmax_x[k + 1] =
	  (x[k] - xl[k]) * (x[k - 1] - xl[k - 1]) / (xu[k - 1] - xl[k - 1]) +
	  xl[k];
	localmax_x[k + 2] =
	  (x[k + 1] - xl[k + 1]) * (x[k] - xl[k]) / (xu[k] -
						     xl[k]) * (x[k - 1] -
							       xl[k -
								  1]) /
	  (xu[k - 1] - xl[k - 1]) + xl[k + 1];
	k += 3;
	if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	  localmax_x[k] = x[k - 1];
	  localmax_x[k + 1] = x[k];
	  localmax_x[k + 2] = x[k + 1];
	  k += 3;
	}
	/* threshold for QT */
	if (modelType.trait == CT) {
	  localmax_x[k] = x[k - 1];
	  k++;
	}
      }
    }
  }

  avg_hetLR = alpha_integral;

  //  if (constraint >= 3.0 || constraint <= -3.0) {
  //fprintf (stderr, "Constraint is %f with avg hetLR =%f\n", constraint,
  //   avg_hetLR);
//}


  /* Jacobian */
  k = 1;
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    avg_hetLR *=
      (x[k] - xl[k]) / (xu[k] - xl[k]) * (x[k] - xl[k]) / (xu[k] -
							   xl[k]) * (x[k +
								       1] -
								     xl[k +
									1]) /
      (xu[k + 1] - xl[k + 1]);
    k += 3;
    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      k += 3;
    }
    if (modelType.trait == CT) {
      k++;
    }
  }
  /*This is a temporary checking. */
  if (k != s->ndim) {
    printf ("k=%d  while dim for BR is %d\n", k, s->ndim);
    exit (0);
  }

  *f = avg_hetLR;


}


void
compute_hlod_2p_dt (double x[], double *f)
{
//double compute_hlod(PedigreeSet *pedigreeSet,double x[], int loc1, int loc2, Locus *pLocus, Trait *pTrait, int traitLocus, int totalLoci, double * initialProbAddr[3], Locus *pLocus1){

/*  Limit of this function

   modelOptions.type := TP  (Two points)
   modelOptions.trait := DT (DICHOTOMOUS);
   modelOptions.equilibrium :=LINKAGE_EQUILIBRIUM
   modelOptions.polynomial := TRUE
   modelOptions->markerAnalysis := FALSE;
   modelOptions.mapFlag := SA 
   modelRange.nafreq :=1
   
*/

  int k, j;
  int pedIdx, liabIdx = 0, status;

  double pen_DD, pen_Dd, pen_dd, gfreq, alphaV, theta;
  double log10_likelihood_null, log10_likelihood_alternative,
    log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, tmp, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;

  Pedigree *pPedigree;

  gfreq = x[0];
  theta = fixed_theta;		//x[5];  
  //printf("Calculating hetLR with gf=%f DD=%f Dd=%f dd=%f theta=%f\n", gfreq, pen_DD,pen_Dd, pen_dd, fixed_theta);
  if (1 && modelOptions.markerAnalysis == FALSE) {
    pLocus->pAlleleFrequency[0] = gfreq;
    pLocus->pAlleleFrequency[1] = 1 - gfreq;

    if (modelOptions.polynomial == TRUE);
    else
      update_locus (&pedigreeSet, loc1);

  }

  if (modelOptions.markerAnalysis == FALSE
      && pLocus1->locusType == LOCUS_TYPE_TRAIT) {
    for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
      pen_DD = x[3 * liabIdx + 1];
      pen_Dd = x[3 * liabIdx + 2] * x[3 * liabIdx + 1];
      pen_dd = x[3 * liabIdx + 3] * x[3 * liabIdx + 1] * x[3 * liabIdx + 2];
      pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
      pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
      pTrait->penetrance[2][liabIdx][1][0] = pen_Dd;
      pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
      pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
      pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
      pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_Dd;
      pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;
    }
  }
  if (modelOptions.polynomial == TRUE);
  else
    update_penetrance (&pedigreeSet, traitLocus);

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
    /* populate the matrix */
    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */


  KLOG (LOGLIKELIHOOD, LOGDEBUG, "NULL Likelihood\n");
  compute_likelihood (&pedigreeSet);

  //printf("likelihood =%15.13f with theta 0.5 with %d pedigrees\n", pedigreeSet.likelihood, pedigreeSet.numPedigree);
  //scanf("%d ",&k);
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    fprintf (stderr, "Theta 0.5 has likelihood 0\n");
    fprintf (stderr, "dgf=%f\n", gfreq);
    exit (-1);
  }

  /* save the results for NULL */
  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
    /* save the likelihood at null */
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
    pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
  }

  log10_likelihood_null = pedigreeSet.log10Likelihood;

  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);
//      if (pLambdaCell->impossibleFlag[dprimeIdx] != 0)
//        continue;
    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
    copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);

    /* calculate R square if the marker is a SNP */
    if (R_square_flag == TRUE)
      R_square =
	calculate_R_square (pLocus1->pAlleleFrequency[0],
			    pLocus2->pAlleleFrequency[0],
			    pLDLoci->ppDValue[0][0]);
    else
      R_square = -1;
  }


  if (modelOptions.mapFlag == SA) {

    for (k = 0; k < 3; k++) {
      locusList->pNextLocusDistance[k][0] = theta;
      locusList->pPrevLocusDistance[k][1] = theta;

    }
  } else {
    printf ("mapflag sould be SA\n");
    exit (-1);
  }

  if (modelOptions.polynomial == TRUE);
  else
    /* populate the matrix */
    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */

  KLOG (LOGLIKELIHOOD, LOGDEBUG, "ALT Likelihood\n");
  compute_likelihood (&pedigreeSet);

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;

  //printf("likelihood =%15.13f with theta %f  %d pedigree\n", pedigreeSet.likelihood,fixed_theta, pedigreeSet.numPedigree);                                  
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
  } else {
    log10_likelihood_ratio =
      log10_likelihood_alternative - log10_likelihood_null;
  }


  /* check for overflow problem !!! */
  if (log10_likelihood_ratio >= DBL_MAX_10_EXP - 1) {
    likelihood_ratio = DBL_MAX;
  } else {
    /* check for underflow problem too !!! */
    if (log10_likelihood_ratio <= DBL_MIN_10_EXP + 1) {
      likelihood_ratio = 0;
    } else {
      likelihood_ratio = pow (10.0, log10_likelihood_ratio);
    }
  }
  /* caculating the HET */
  for (j = 0; j < 6; j++) {
    //for (j = 0; j < 1; j++) {
    alphaV = alpha[j][0];
    alphaV2 = 1 - alphaV;
    if (alphaV2 < 0)
      alphaV2 = 0;

    log10HetLR = 0;

    for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
      homoLR = pPedigree->likelihood / pedigreeSet.nullLikelihood[pedIdx];
      tmp = log10 (alphaV * homoLR + (1 - alphaV));
      log10HetLR += tmp * pPedigree->pCount[loc2];
    }

    if (log10HetLR >= __DBL_MAX_10_EXP__ - 1) {
      hetLR = __DBL_MAX__;
    } else if (log10HetLR <= __DBL_MIN_10_EXP__ + 1) {
      hetLR = 0;
    } else {
      hetLR = pow (10, log10HetLR);
    }

    if(j==0){
      alpha_integral += hetLR;// * alpha[j][1];
    }
    

    if (hetLR > localmax_value) {
      localmax_value = hetLR;
      localmax_x[0] = gfreq;
      localmax_x[1] = alphaV;

      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	localmax_x[3 * liabIdx + 2] = x[3 * liabIdx + 1];
	localmax_x[3 * liabIdx + 3] = x[3 * liabIdx + 2] * x[3 * liabIdx + 1];
	localmax_x[3 * liabIdx + 4] =
	  x[3 * liabIdx + 3] * x[3 * liabIdx + 1] * x[3 * liabIdx + 2];
      }
    }
    // fprintf(fphlod,"%f %f %f %f %f %f %f\n", log10(hetLR*x[1]*x[1]*x[2]), gfreq, pen_DD,pen_Dd, pen_dd, alphaV,fixed_theta);
  }				//end of calculating the HET         


  avg_hetLR = alpha_integral;
  //avg_hetLR= alpha_integral/modelRange.nalpha;

  //printf("avg hetLR =%15.10f with gf=%f DD=%f Dd=%f dd=%f theta=%f\n", avg_hetLR, gfreq, pen_DD,pen_Dd, pen_dd, fixed_theta);

  // avg_hetLR *=x[1]*x[1]*x[2];
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    avg_hetLR *= x[3 * liabIdx + 1] * x[3 * liabIdx + 1] * x[3 * liabIdx + 2];
  }


  f[0] = avg_hetLR;

  //  return 0.0;

}
