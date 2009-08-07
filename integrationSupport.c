#include "integrationSupport.h"
extern LDLoci *pLDLoci;
extern int R_square_flag;
extern double R_square;

int kelvin_dcuhre_integrate (double *integral, double *abserr, double vol_region, int *scale){

  /* INPUT
       vol_region : the volume of the rectangular region
           (( s->vol_rate is the rate to convert DCUHRE's output into the average))

     OUTPUT
       integral : average funtion value in the given region
       error    : estimate error in calculation of integral
       scale    : scale used in calculation of BR

     VARIABLE
       dim : number of variables in the middle layer of the 3-layer approach, which is dim of DCUHRE
       s   : sturcture to hold all global information in using DCUHRE
             and the only variable to change with different analysis
       localMOD : maximum in each BR  
                        maximization information is in localmax_x[]
       boost_rate : BR boosting BR^boost_rate

  */


  /* Local variables */
  int dim, return_val,i;
  double boost_rate=1.1; //1.1;

  /* 
  boost_rate=1.0;
  for(i=0;i<modelRange.nlclass;i++){
    boost_rate *= 1.3;
  }*/

  if(modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM) 
    boost_rate =1.0;
  //extern /* Subroutine */ int ftest_();  

  localMOD = DBL_MIN_10_EXP ;

  if (modelType.trait == DICHOTOMOUS) {

    dim = 1 + 3 * modelRange.nlclass;
    if(modelOptions.imprintingFlag)
      dim += modelRange.nlclass;	//dD

    s = &init_state;
    initialize_state (s, xl, xu, dim);

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
        if(modelOptions.mapFlag == SA)
	  dim -= 1;		//theta
        else
          dim -= 2;
      }
    }

    s = &init_state;
    initialize_state (s, xl, xu, dim);
   
    if (modelType.type == TP) {
      s->funsub = (U_fp) compute_hlod_2p_qt;
      s->mType = TP_DT;
    } else {
      s->funsub = (U_fp) compute_hlod_mp_qt;
      s->mType = MP_DT;
    }
  }
  if (modelOptions.maxIterations > -1) {
    s->maxcls = modelOptions.maxIterations;
  } else if (dim <10) {
    s->maxcls = 50000;
  } else {
    s->maxcls = 80* (int)pow(2.0,dim);
    //fprintf(stdout,"New maxcls is %d \n", s->maxcls);
  }

  s->verbose = 0;
  s->nlclass = modelRange.nlclass;

  for(i=0; i<s->nlclass; i++){
    if(modelOptions.imprintingFlag)
      s->vol_rate /= 16.0;
    else
      s->vol_rate /= 6.0; 
  }
  s->vol_rate *= vol_region;   /*This is the rate to convert to average function value*/

  if(s->verbose >0)
    fprintf (stderr,"Starting DCUHRE with dim=%d\n", dim);

  return_val = dcuhre_ (s);
  if(return_val >0){
    return return_val;
  }

  s->result /= s->vol_rate;  
  s->error /= s->vol_rate;

  if (return_val > 0 && return_val < 20) {
    fprintf (stderr, "Ending program with error! ifail =%d \n", s->ifail);
  }

  if(s->verbose >0){
    fprintf (stderr,
	   "Final result =%15.10f  with error =%15.10f and neval = %d\n",
	   s->result, s->error, s->total_neval);
    fprintf (stderr, "End of DCUHRE with ifail =%d\n", s->ifail);
  }

  /* BR boosting is done here */
  //fprintf(stderr, "Before boosting %e\n", s->result);
  if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM && modelType.trait == DT) {
  //  if (modelType.trait == DT) {  // boosting 
    //    for (i = 0; i < modelRange.nlclass; i++) 
      s->result = pow(10.0, (log10(s->result) * boost_rate));
    //fprintf(stderr, "After boosting %e\n", s->result);
  }

  *integral = s->result;
  *abserr = s->error;
  *scale = s->scale;
  return return_val;

}


void
compute_hlod_mp_qt (double x[], double *f, int *scale)
{

  int k, j;
  int pedIdx, liabIdx = 0, status,pen_size=3;
  double constraint = 0.0;
  double mean_DD=0.0, mean_Dd=0.0, mean_dD=0.0, mean_dd=0.0;
  double SD_DD=0.0,SD_Dd=0.0, SD_dD=0.0, SD_dd=0.0;
  double gfreq;
  double alphaV;
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

  int newscale, oldscale;  /* scaling related variables*/
  double newLog10HetLR;
  double oldsum;
  double oldsum_log10;
  double newsum_log10;


  if(modelOptions.imprintingFlag)
    pen_size=4;

  int origLocus = locusList->pLocusIndex[0];

  if (locusList->numLocus > 1)
    origLocus = locusList->pLocusIndex[1];

  gfreq = x[0];
  if(fpIR !=NULL)
    dk_curModel.dgf = gfreq;

  pLocus->pAlleleFrequency[0] = gfreq;
  pLocus->pAlleleFrequency[1] = 1 - gfreq;

  update_locus (&pedigreeSet, traitLocus);


  j = 1;			// j=0 for gfrequency
  if (modelType.trait == CT) {
    threshold= x[s->ndim -1];
  }
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    mean_DD = x[j];
    mean_Dd = (x[j+1]-xl[j+1])*(x[j]-xl[j])/(xu[j]-xl[j])+xl[j+1];

    if(modelOptions.imprintingFlag){
      mean_dD =(x[j+2]-xl[j+2])*(x[j]-xl[j])/(xu[j]-xl[j])+xl[j+2];
      mean_dd =(x[j+3]-xl[j+3])*(x[j+2]-xl[j+2])/(xu[j+2]- xl[j+2])*(x[j+1]-xl[j+1])/(xu[j+1]- xl[j+1])*(x[j]-xl[j])/(xu[j]-xl[j])+xl[j+3];
    }else{
      mean_dd =(x[j+2]-xl[j+2])*(x[j+1]-xl[j+1])/(xu[j+1]- xl[j+1])*(x[j]-xl[j])/(xu[j]-xl[j])+xl[j+2];
      mean_dD = mean_Dd;
    }
    j += pen_size;

    if(fpIR !=NULL){
      dk_curModel.pen[liabIdx].DD = mean_DD;
      dk_curModel.pen[liabIdx].Dd = mean_Dd;
      dk_curModel.pen[liabIdx].dD = mean_dD;
      dk_curModel.pen[liabIdx].dd = mean_dd;
    }


    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      /*SD_DD = x[j];		
      SD_Dd = x[j + 1];	
      if(modelOptions.imprintingFlag){
	SD_dD = x[j+2];
        SD_dd = x[j+3];
      }else{
	SD_dd = x[j+2];	
        SD_dD= SD_Dd;
      }
      j += pen_size;*/
      SD_DD= SD_Dd=SD_dD= SD_dd = x[j++];

      if(fpIR !=NULL){
        dk_curModel.pen[liabIdx].DDSD = SD_DD;
        dk_curModel.pen[liabIdx].DdSD = SD_Dd;
        dk_curModel.pen[liabIdx].dDSD = SD_dD;
        dk_curModel.pen[liabIdx].ddSD = SD_dd;
      }

    }
    if(fpIR !=NULL){
      if (modelType.trait == CT)
          dk_curModel.pen[liabIdx].threshold= threshold;
    }

    /* threshold for QT *
    if (modelType.trait == CT) {
      threshold = x[j];	// modelRange.tthresh[liabIdx][thresholdIdx];
      j++;
      }*/

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

  if (modelOptions.polynomial == TRUE);
  else
    /* only need to update trait locus */
    update_penetrance (&pedigreeSet, traitLocus);

  /*This is a temporary checking. *
  if (j != s->ndim) {
    fprintf (stderr, "j=%d  while dim for BR is %d\n", j, s->ndim);
    exit (EXIT_FAILURE);
    }*/

  /* for trait likelihood */
  locusList = &traitLocusList;
  xmissionMatrix = traitMatrix;

  /* compute the null likelihood with   */
  pedigreeSet.likelihood = 1;
  pedigreeSet.log10Likelihood = 0;

  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {

    /* save the likelihood at null */
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];

    if (modelOptions.polynomial == TRUE) {
      KASSERT (pPedigree->traitLikelihoodPolynomial != NULL, "Error in  \n");
      /* evaluate likelihood */
      evaluatePoly (pPedigree->traitLikelihoodPolynomial,
		    pPedigree->traitLikelihoodPolyList,
		    &pPedigree->likelihood);
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

      f[0]=1.0;
      return ;
      //       break;
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
  }

  pedigreeSet.likelihood = product_likelihood;
  pedigreeSet.log10Likelihood = sum_log_likelihood;
  log10_likelihood_null = pedigreeSet.log10Likelihood;
  KLOG (LOGLIKELIHOOD, LOGDEBUG, "Sum of log Likelihood is: %e\n",
	sum_log_likelihood);

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


  sprintf (partialPolynomialFunctionName, "MQA_C%d_P%%sM",
	   (originalLocusList.ppLocusList[1])->pMapUnit->chromosome);
  compute_likelihood (&pedigreeSet);
  cL[3]++;

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;
  if (isnan (log10_likelihood_alternative))
    fprintf (stderr, "ALT likelihood is NAN.\n");
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
    avg_hetLR=0.0;
  } else {
    log10_likelihood_ratio =
      log10_likelihood_alternative - log10_likelihood_null -
      pedigreeSet.log10MarkerLikelihood;
    // }
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

  /* caculating the HET */
  for (j = 0; j < 5; j++) {
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
      if (alphaV * homoLR + alphaV2 < 0)
	fprintf (stderr, "HET LR less than 0. Check!!!\n");
      log10HetLR += log10 (alphaV * homoLR + alphaV2);
    }

    if(fpIR !=NULL){
      dk_curModel.alpha = alphaV;
      fprintf(fpIR,"%6.3f", log10HetLR);

      fprintf(fpIR," %4.3f %4.3f",dk_curModel.alpha,dk_curModel.dgf);
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
        fprintf(fpIR," %4.3f %4.3f",dk_curModel.pen[liabIdx].DD,dk_curModel.pen[liabIdx].Dd);
        if(modelOptions.imprintingFlag){
          fprintf(fpIR," %4.3f %4.3f",dk_curModel.pen[liabIdx].dD,dk_curModel.pen[liabIdx].dd);
        }else{
          fprintf(fpIR," %4.3f",dk_curModel.pen[liabIdx].dd); 
        }
        if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
          fprintf(fpIR," %4.3f",dk_curModel.pen[liabIdx].DDSD); 
        }
      }
      if (modelType.trait == CT){
        fprintf(fpIR," %4.3f",dk_curModel.pen[0].threshold); 
      }
      fprintf(fpIR," %d\n",dk_curModel.posIdx);
    }

    oldsum=alpha_integral;
    oldscale= (*scale);
    newscale = 0;

    if (log10HetLR >= DBL_MAX_10_EXP - 1) {
      /* find the new scale, adjust the current sum */
      newscale=log10HetLR - (DBL_MAX_10_EXP - SCALE_RESERVE);
    } 
    if(newscale > oldscale) {
      /* need to use the newscale and adjust the old sum */
      if(oldsum>0) {
        oldsum_log10=log10(oldsum);
        newsum_log10=oldsum_log10+oldscale-newscale;
        if(newsum_log10<= DBL_MIN_10_EXP+1) {
          alpha_integral=0;
        }
        else {
          alpha_integral =pow(10, newsum_log10);
        }
      }
      *scale=newscale;
      oldscale=newscale;
    }else {
      /* use the old scale to adjust the new value */
      newscale = oldscale;
    }
    newLog10HetLR = log10HetLR - newscale; 

    if (newLog10HetLR <= DBL_MIN_10_EXP + 1) {
      hetLR = 0;
    } else {
      hetLR = pow (10, newLog10HetLR);
    }

    alpha_integral += hetLR * alpha[j][1];


    /*Update local maximum as necessary */
    if (log10HetLR > localMOD) {
      localMOD = log10HetLR;
      localmax_x[0] = gfreq;
      localmax_x[1] = alphaV;
      k = 2;
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	localmax_x[k] = x[k - 1];
	localmax_x[k + 1] =(x[k]-xl[k])*(x[k-1]-xl[k-1])/(xu[k-1]-xl[k-1])+xl[k];

        if(modelOptions.imprintingFlag){
	  localmax_x[k+2]=(x[k+1]-xl[k+1])*(x[k-1]-xl[k-1])/(xu[k-1]-xl[k-1])+xl[k+1];
          localmax_x[k+3]=(x[k+2]-xl[k+2])*(x[k+1]-xl[k+1])/(xu[k+1]-xl[k+1])*(x[k]-xl[k])/(xu[k]-xl[k])*(x[k-1]-xl[k-1])/(xu[k-1]-xl[k-1])+xl[k+2];
	}else{
	  localmax_x[k+2]=(x[k+1]-xl[k+1])*(x[k]-xl[k])/(xu[k]-xl[k])*(x[k-1]-xl[k-1])/(xu[k-1]-xl[k-1])+xl[k+1];
	}
	k += pen_size;
	if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	  /*localmax_x[k] = x[k - 1];
	  localmax_x[k + 1] = x[k];
	  localmax_x[k + 2] = x[k + 1];
          if(modelOptions.imprintingFlag)
	    localmax_x[k + 3] = x[k + 2];
	    k += pen_size;*/
          localmax_x[k] = x[k - 1];
          k++;
	}
	/* threshold for QT *
	if (modelType.trait == CT) {
	  localmax_x[k] = x[k - 1];
	  k++;
	  }*/
      }
      if (modelType.trait == CT) {
	localmax_x[k] = x[k - 1];
	k++;
      }
    }
  }

  avg_hetLR = alpha_integral;

  /* Jacobian */
  k = 1;
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    avg_hetLR *=(x[k]-xl[k])/(xu[k]-xl[k]);
    avg_hetLR *=(x[k]-xl[k])/(xu[k]-xl[k]);
    avg_hetLR *= (x[k+1]-xl[k+1])/(xu[k+1]-xl[k+1]);

    if(modelOptions.imprintingFlag){
      avg_hetLR *=(x[k]-xl[k])/(xu[k]-xl[k]);
      avg_hetLR *= (x[k+2]-xl[k+2])/(xu[k+2]-xl[k+2]);
    }

    k += pen_size;
    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      k ++; //+= pen_size;
    }
    /*if (modelType.trait == CT) {
      k++;
      }*/
  }
  /*This is a temporary checking. *
  if (k != s->ndim) {
    fprintf (stderr, "k=%d  while dim for BR is %d\n", k, s->ndim);
    exit (EXIT_FAILURE);
    }*/
  }
  f[0] = avg_hetLR;

}


void
compute_hlod_mp_dt (double x[], double *f, int *scale)
{

  int j;
  int pedIdx, liabIdx, status,pen_size=3, ret;

  double pen_DD, pen_Dd,pen_dD, pen_dd, gfreq, alphaV;
  double log10_likelihood_null, log10_likelihood_alternative,
    log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;
  double log10Likelihood;

  /* for null likelihood calculation */
  double product_likelihood = 1;	/* product of the likelihoods for all the pedigrees */
  double sum_log_likelihood = 0;	/* sum of the log10(likelihood) for all the pedigrees */

  Pedigree *pPedigree;

  int newscale, oldscale;  /* scaling related variables*/
  double newLog10HetLR;
  double oldsum;
  double oldsum_log10;
  double newsum_log10;



  int origLocus = locusList->pLocusIndex[0];

  if(modelOptions.imprintingFlag)
    pen_size=4;

  if (locusList->numLocus > 1)
    origLocus = locusList->pLocusIndex[1];

  //fprintf(stderr,"in compute hlod x %G %G %G %G\n", x[0],x[1],x[2],x[3]);
  gfreq = x[0];
  if(fpIR !=NULL)
    dk_curModel.dgf = gfreq;

  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    pen_DD = x[pen_size * liabIdx + 1];
    pen_Dd = x[pen_size * liabIdx + 2] * x[pen_size * liabIdx + 1];

    if(modelOptions.imprintingFlag){
      pen_dD = x[pen_size * liabIdx + 3] * x[pen_size * liabIdx + 1];
      pen_dd = x[pen_size * liabIdx + 4] * x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 2] *x[pen_size * liabIdx + 3];
    }else{
      pen_dd = x[pen_size * liabIdx + 3] * x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 2];
      pen_dD= pen_Dd;
    }
    if(fpIR !=NULL){
      dk_curModel.pen[liabIdx].DD = pen_DD;
      dk_curModel.pen[liabIdx].Dd = pen_Dd;
      dk_curModel.pen[liabIdx].dD = pen_dD;
      dk_curModel.pen[liabIdx].dd = pen_dd;
    }

    pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
    pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
    pTrait->penetrance[2][liabIdx][1][0] = pen_Dd;
    pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
    pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
    pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
    pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_Dd;
    pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;

  }

  if (modelOptions.polynomial == TRUE);
  else
    /* only need to update trait locus */
    update_penetrance (&pedigreeSet, traitLocus);

  pLocus->pAlleleFrequency[0] = gfreq;
  pLocus->pAlleleFrequency[1] = 1 - gfreq;


  if (modelOptions.polynomial != TRUE)
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
      ret=-1;
      product_likelihood = 0.0;
      sum_log_likelihood = -9999.99;
      break;
    } else if (pPedigree->likelihood < 0.0) {
      KASSERT (pPedigree->likelihood >= 0.0,
	       "Pedigree %s with NEGATIVE likelihood - This is CRAZY!!!.\n",
	       pPedigree->sPedigreeID);
      product_likelihood = 0.0;
      sum_log_likelihood = -9999.99;
      ret=-2;
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
  int k;
  char markerNo[8];
  sprintf (partialPolynomialFunctionName, "MDA_C%d_P%%sM",
	   (originalLocusList.ppLocusList[1])->pMapUnit->chromosome);
  for (k = 0; k < modelType.numMarkers; k++) {
    if (*get_map_position (traitLocus) <= *get_map_position (markerLocusList.pLocusIndex[k]) &&
	(strstr (partialPolynomialFunctionName, "_T") == NULL))
      strcat (partialPolynomialFunctionName, "_T");
    sprintf (markerNo, "_%d", markerLocusList.pLocusIndex[k]);
    strcat (partialPolynomialFunctionName, markerNo);
  }
  if (strstr (partialPolynomialFunctionName, "_T") == NULL)
    strcat (partialPolynomialFunctionName, "_T");
  cL[4]++;
  ret=compute_likelihood (&pedigreeSet);
  if(ret==-2){
    /* negative likelihood */
     fprintf(stderr, "Negative likelihood! Exiting!!\n");
     exit(EXIT_FAILURE);
  }


  log10_likelihood_alternative = pedigreeSet.log10Likelihood;
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
    avg_hetLR=0.0;
  } else {
    log10_likelihood_ratio =
      log10_likelihood_alternative - log10_likelihood_null;
    // }
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

  /* caculating the HET */
  for (j = 0; j < 5; j++) {
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

    if(fpIR !=NULL){
      dk_curModel.alpha = alphaV;
      fprintf(fpIR,"%6.3f", log10HetLR);

      fprintf(fpIR," %4.3f %4.3f",dk_curModel.alpha,dk_curModel.dgf);
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
        fprintf(fpIR," %4.3f %4.3f",dk_curModel.pen[liabIdx].DD,dk_curModel.pen[liabIdx].Dd);
        if(modelOptions.imprintingFlag){
          fprintf(fpIR," %4.3f %4.3f",dk_curModel.pen[liabIdx].dD,dk_curModel.pen[liabIdx].dd);
        }else{
          fprintf(fpIR," %4.3f",dk_curModel.pen[liabIdx].dd); 
        }
      }
      fprintf(fpIR," %d\n",dk_curModel.posIdx);
    }

    oldsum=alpha_integral;
    oldscale= (*scale);
    newscale = 0;

    if (log10HetLR >= DBL_MAX_10_EXP - 1) {
      /* find the new scale, adjust the current sum */
      newscale=log10HetLR - (DBL_MAX_10_EXP - SCALE_RESERVE);
    } 
    if(newscale > oldscale) {
      /* need to use the newscale and adjust the old sum */
      if(oldsum>0) {
        oldsum_log10=log10(oldsum);
        newsum_log10=oldsum_log10+oldscale-newscale;
        if(newsum_log10<= DBL_MIN_10_EXP+1) {
          alpha_integral=0;
        }
        else {
          alpha_integral =pow(10, newsum_log10);
        }
      }
      *scale=newscale;
      oldscale=newscale;
    }else {
      /* use the old scale to adjust the new value */
      newscale = oldscale;
    }
    newLog10HetLR = log10HetLR - newscale; 

    if (newLog10HetLR <= DBL_MIN_10_EXP + 1) {
      hetLR = 0;
    } else {
      hetLR = pow (10, newLog10HetLR);
    }

    alpha_integral += hetLR * alpha[j][1];

    if (log10HetLR > localMOD) {
      localMOD = log10HetLR;
      localmax_x[0] = gfreq;
      localmax_x[1] = alphaV;

      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	localmax_x[pen_size * liabIdx + 2] = x[pen_size * liabIdx + 1];
	localmax_x[pen_size * liabIdx + 3] = x[pen_size * liabIdx + 2] * x[pen_size * liabIdx + 1];
	localmax_x[pen_size * liabIdx + 4] = x[pen_size * liabIdx + 3] * x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 2];

        if(modelOptions.imprintingFlag){
   	  localmax_x[pen_size * liabIdx + 4] = x[pen_size * liabIdx + 3] * x[pen_size * liabIdx + 1];
	  localmax_x[pen_size * liabIdx + 5] =
	    x[pen_size * liabIdx + 4] * x[pen_size * liabIdx + 3] * x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 2];
	}
      }
    }

  }				/* end of calculating HET LR */

  avg_hetLR = alpha_integral;

  //Jacobian
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    if(modelOptions.imprintingFlag)
      avg_hetLR *= x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 1]* x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 2]* x[pen_size * liabIdx + 3];
    else
      avg_hetLR *= x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 2];
  }
  }
  *f = avg_hetLR;
}

void
compute_hlod_2p_qt (double x[], double *f, int *scale)
{

  int k, j, ret;
  int pedIdx, liabIdx = 0, status, pen_size=3;
  double constraint = 0.0;
  double mean_DD=0.0, mean_Dd=0.0, mean_dD=0.0, mean_dd=0.0;
  double SD_DD=0.0,SD_Dd=0.0, SD_dD=0.0, SD_dd=0.0;
  double gfreq;
  double alphaV;
  double thetaM, thetaF;
  double threshold = 0.0;
  double log10_likelihood_null, log10_likelihood_alternative,
    log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;

  Pedigree *pPedigree;

  int newscale, oldscale;  /* scaling related variables*/
  double newLog10HetLR;
  double oldsum;
  double oldsum_log10;
  double newsum_log10;

  if(modelOptions.imprintingFlag)
    pen_size=4;

  gfreq = x[0];
  if(fpIR !=NULL)
    dk_curModel.dgf = gfreq;

  if(modelOptions.mapFlag == SS){
    thetaM = fixed_thetaM;
    thetaF = fixed_thetaF;
  }else{
    thetaM = fixed_theta;
    thetaF = fixed_theta;
  }

  if (1 && modelOptions.markerAnalysis == FALSE) {
    pLocus->pAlleleFrequency[0] = gfreq;
    pLocus->pAlleleFrequency[1] = 1 - gfreq;

    if (modelOptions.polynomial == TRUE);
    else
      update_locus (&pedigreeSet, loc1);

  }

  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    status = setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprimeIdx);
    if(status<0)
      KASSERT (1,"Haplotype frequency combination impossible. Exiting!\n");
  }


  /* this should be MEAN + SD */
  j = 1;
  if (modelType.trait == CT) {
    threshold= x[s->ndim -1];
  }
  if (modelOptions.markerAnalysis == FALSE) {
    for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
      mean_DD = x[j];
      mean_Dd = (x[j+1]-xl[j+1])*(x[j]-xl[j])/(xu[j]-xl[j])+xl[j+1];

      if(modelOptions.imprintingFlag){
        mean_dD =(x[j+2]-xl[j+2])*(x[j]-xl[j])/(xu[j]-xl[j])+xl[j+2];
        mean_dd =(x[j+3]-xl[j+3])*(x[j+2]-xl[j+2])/(xu[j+2]- xl[j+2])*(x[j+1]-xl[j+1])/(xu[j+1]- xl[j+1])*(x[j]-xl[j])/(xu[j]-xl[j])+xl[j+3];
      }else{
        mean_dd =(x[j+2]-xl[j+2])*(x[j+1]-xl[j+1])/(xu[j+1]- xl[j+1])*(x[j]-xl[j])/(xu[j]-xl[j])+xl[j+2];
        mean_dD = mean_Dd;
      }
      j += pen_size;

      if(fpIR !=NULL){
        dk_curModel.pen[liabIdx].DD = mean_DD;
        dk_curModel.pen[liabIdx].Dd = mean_Dd;
        dk_curModel.pen[liabIdx].dD = mean_dD;
        dk_curModel.pen[liabIdx].dd = mean_dd;
      }

      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	/*SD_DD = x[j];		
	SD_Dd = x[j + 1];	
        if(modelOptions.imprintingFlag){
	  SD_dD = x[j+2];
          SD_dd = x[j+3];
	}else{
	  SD_dd = x[j+2];	
          SD_dD= SD_Dd;
	}
	j += pen_size;*/
        SD_DD= SD_Dd=SD_dD= SD_dd = x[j++];

        if(fpIR !=NULL){
          dk_curModel.pen[liabIdx].DDSD = SD_DD;
          dk_curModel.pen[liabIdx].DdSD = SD_Dd;
          dk_curModel.pen[liabIdx].dDSD = SD_dD;
          dk_curModel.pen[liabIdx].ddSD = SD_dd;
        }

      }
      if(fpIR !=NULL){
        if (modelType.trait == CT)
          dk_curModel.pen[liabIdx].threshold= threshold;
      }

      /* threshold for QT *
      if (modelType.trait == CT) {
	threshold = x[j];	// modelRange.tthresh[liabIdx][thresholdIdx];
	j++;
	}*/


      /* check against the hard coded constraint */
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	constraint =(1.0-gfreq)*(1.0-gfreq)*mean_dd*SD_dd+2*gfreq*(1-gfreq)*mean_Dd*SD_Dd+gfreq*gfreq*mean_DD*SD_DD;
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
      pTrait->means[liabIdx][1][0] = mean_dD;
      pTrait->means[liabIdx][1][1] = mean_dd;
      pTrait->stddev[liabIdx][0][0] = SD_DD;
      pTrait->stddev[liabIdx][0][1] = SD_Dd;
      pTrait->stddev[liabIdx][1][0] = SD_dD;
      pTrait->stddev[liabIdx][1][1] = SD_dd;

      /* threshold for QT */
      pTrait->cutoffValue[liabIdx] = threshold;

    }				/* liability class Index */

    if (modelOptions.polynomial == TRUE);
    else
      update_penetrance (&pedigreeSet, traitLocus);

  }				/* marker to marker analysis */

  /*This is a temporary checking. *
  if (j != s->ndim) {
    fprintf (stderr, "j=%d  while dim for BR is %d\n", j, s->ndim);
    exit (0);
    }*/

  /* get the likelihood at 0.5 first and LD=0 */
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {

    status = setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprime0Idx);
    if(status<0)
      KASSERT (1,"Haplotype frequency combination impossible. Exiting!\n");

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


  sprintf (partialPolynomialFunctionName, "TQ_C%d_P%%s_%s_%s",
	   pLocus2->pMapUnit->chromosome, 
	   pLocus1->sName, pLocus2->sName);
  cL[5]++;
  ret=compute_likelihood (&pedigreeSet);

  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    fprintf (stderr, "Theta 0.5 has likelihood 0 \n");
    fprintf (stderr, "dgf=%f\n", gfreq);
    for (j = 1; j < s->ndim; j++) {
      fprintf (stderr, " %f", x[j]);
    }
    fprintf(stderr, "mean %f %f %f %f SD %f %f %f %f\n",mean_DD, mean_Dd,mean_dD, mean_dd, SD_DD, SD_Dd, SD_dD,SD_dd );
    fprintf (stderr, "\n");

    //exit (EXIT_FAILURE);
    *f = 1.0;
    return ;
  }

  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
    /* save the likelihood at null */
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
    pedigreeSet.nullLikelihood[pedIdx] = pPedigree->likelihood;
  }

  log10_likelihood_null = pedigreeSet.log10Likelihood;

  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);
    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
    copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);
  }

  if (modelOptions.mapFlag == SA) {
    for (k = 0; k < 3; k++) {
      locusList->pNextLocusDistance[k][0] = thetaM;
      locusList->pPrevLocusDistance[k][1] = thetaF;
    }
  } else {
    locusList->pNextLocusDistance[MAP_MALE][0] =
    locusList->pPrevLocusDistance[MAP_MALE][1] = thetaM;
    locusList->pNextLocusDistance[MAP_FEMALE][0] =
    locusList->pPrevLocusDistance[MAP_FEMALE][1] = thetaF;
  }

  if (modelOptions.polynomial == TRUE);
  else
    /* populate the matrix */
    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */

  // No new name for a polynomial here because we're reusing the existing one
  cL[6]++;
  ret=compute_likelihood (&pedigreeSet);
  if(ret==-2){
    /* negative likelihood */
    fprintf(stderr, "Negative likelihood! Exiting!\n");
    exit(EXIT_FAILURE);
  }

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
    avg_hetLR=0.0;

    fprintf (stderr, "Theta %f has likelihood 0 \n", thetaM);
    fprintf (stderr, "dgf=%f\n", gfreq);
    for (j = 1; j < s->ndim; j++) {
      fprintf (stderr, " %f", x[j]);
    }
    fprintf(stderr, "mean %f %f %f %f SD %f %f %f %f\n",mean_DD, mean_Dd,mean_dD, mean_dd, SD_DD, SD_Dd, SD_dD,SD_dd );
    fprintf (stderr, "\n");

  } else {
    log10_likelihood_ratio =
      log10_likelihood_alternative - log10_likelihood_null;
    //  }

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
  for (j = 0; j < 5; j++) {
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

      if(isnan(homoLR)){
        printf("pedIdx =%d  homeLR=%e log10HLR=%e\n",pedIdx, homoLR, log10HetLR );
        exit(0);
      }
    }

    if(fpIR !=NULL){
      dk_curModel.alpha = alphaV;
      fprintf(fpIR,"%6.3f", log10HetLR);
      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
        fprintf(fpIR," %4.3f", dk_curModel.dprime[0]);
      }
      if(modelOptions.mapFlag == SA){
        fprintf(fpIR," %4.3f", dk_curModel.theta[0]);
      } else{
        fprintf(fpIR," %4.3f %4.3f", dk_curModel.theta[0], dk_curModel.theta[0] );
      }
      fprintf(fpIR," %4.3f %4.3f",dk_curModel.alpha,dk_curModel.dgf);
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
        fprintf(fpIR," %4.3f %4.3f",dk_curModel.pen[liabIdx].DD,dk_curModel.pen[liabIdx].Dd);
        if(modelOptions.imprintingFlag){
          fprintf(fpIR," %4.3f %4.3f",dk_curModel.pen[liabIdx].dD,dk_curModel.pen[liabIdx].dd);
        }else{
          fprintf(fpIR," %4.3f",dk_curModel.pen[liabIdx].dd); 
        }
        if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
          fprintf(fpIR," %4.3f",dk_curModel.pen[liabIdx].DDSD); 
	}
      }
      if (modelType.trait == CT){
        fprintf(fpIR," %4.3f",dk_curModel.pen[0].threshold); 
      }
      fprintf(fpIR," %d\n",dk_curModel.posIdx);
    }

    oldsum=alpha_integral;
    oldscale= (*scale);
    newscale = 0;

    if (log10HetLR >= DBL_MAX_10_EXP - 1) {
      /* find the new scale, adjust the current sum */
      newscale=log10HetLR - (DBL_MAX_10_EXP - SCALE_RESERVE);  
    } 
    if(newscale > oldscale) {
      /* need to use the newscale and adjust the old sum */
      if(oldsum>0) {
        oldsum_log10=log10(oldsum);
        newsum_log10=oldsum_log10+oldscale-newscale;
        if(newsum_log10<= DBL_MIN_10_EXP+1) {
          alpha_integral=0;
        }
        else {
          alpha_integral =pow(10, newsum_log10);
        }
      }
      *scale=newscale;
      oldscale=newscale;
    }else {
      /* use the old scale to adjust the new value */
      newscale = oldscale;
    }
    newLog10HetLR = log10HetLR - newscale; 

    if (newLog10HetLR <= DBL_MIN_10_EXP + 1) {
      hetLR = 0;
    } else {
      hetLR = pow (10, newLog10HetLR);
    }
    
    alpha_integral += hetLR * alpha[j][1];
    

    /*Update local maximum as necessary */
    if (log10HetLR > localMOD) {
      localMOD = log10HetLR;
      localmax_x[0] = gfreq;
      localmax_x[1] = alphaV;
      k = 2;
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {

	localmax_x[k] = x[k - 1];
	localmax_x[k + 1] =(x[k]-xl[k])*(x[k-1]-xl[k-1])/(xu[k-1]-xl[k-1])+xl[k];

        if(modelOptions.imprintingFlag){
	  localmax_x[k+2]=(x[k+1]-xl[k+1])*(x[k-1]-xl[k-1])/(xu[k-1]-xl[k-1])+xl[k+1];
          localmax_x[k+3]=(x[k+2]-xl[k+2])*(x[k+1]-xl[k+1])/(xu[k+1]-xl[k+1])*(x[k]-xl[k])/(xu[k]-xl[k])*(x[k-1]-xl[k-1])/(xu[k-1]-xl[k-1])+xl[k+2];
	}else{
	  localmax_x[k+2]=(x[k+1]-xl[k+1])*(x[k]-xl[k])/(xu[k]-xl[k])*(x[k-1]-xl[k-1])/(xu[k-1]-xl[k-1])+xl[k+1];
	}
	k += pen_size;
	if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	  /*localmax_x[k] = x[k - 1];
	  localmax_x[k + 1] = x[k];
	  localmax_x[k + 2] = x[k + 1];
          if(modelOptions.imprintingFlag)
	    localmax_x[k + 3] = x[k + 2];
	    k += pen_size;*/
          localmax_x[k] = x[k - 1];
          k++;
	}
	/* threshold for QT *
	if (modelType.trait == CT) {
	  localmax_x[k] = x[k - 1];
	  k++;
	  }*/
      }
      if (modelType.trait == CT) {
	localmax_x[k] = x[k - 1];
	k++;
      }
    }
  }

  avg_hetLR = alpha_integral;

  /* Jacobian */
  k = 1;
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    avg_hetLR *=(x[k]-xl[k])/(xu[k]-xl[k]);
    avg_hetLR *=(x[k]-xl[k])/(xu[k]-xl[k]);
    avg_hetLR *= (x[k+1]-xl[k+1])/(xu[k+1]-xl[k+1]);

    if(modelOptions.imprintingFlag){
      avg_hetLR *=(x[k]-xl[k])/(xu[k]-xl[k]);
      avg_hetLR *= (x[k+2]-xl[k+2])/(xu[k+2]-xl[k+2]);
    }

    k += pen_size;
    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      k ++; //= pen_size;
    }
    /*if (modelType.trait == CT) {
      k++;
      }*/
  }
  /*This is a temporary checking. *
  if (k != s->ndim) {
    printf ("k=%d  while dim for BR is %d\n", k, s->ndim);
    exit (EXIT_FAILURE);
    }*/

  }
  *f = avg_hetLR;
}

void compute_hlod_2p_dt (double x[], double *f, int *scale) {

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
  int pedIdx, liabIdx = 0, status, pen_size=3,ret;

  double pen_DD, pen_Dd, pen_dD,pen_dd, gfreq, alphaV, thetaM, thetaF;
  double log10_likelihood_null, log10_likelihood_alternative,
    log10_likelihood_ratio, likelihood_ratio;
  double log10HetLR, homoLR, alphaV2;
  double hetLR, avg_hetLR, alpha_integral=0.0,tmp;     
  Pedigree *pPedigree;

  int newscale, oldscale;  /* scaling related variables*/
  double newLog10HetLR;
  double oldsum;
  double oldsum_log10;
  double newsum_log10;

  
  if(modelOptions.imprintingFlag)
    pen_size=4;

  gfreq = x[0];
  // printf("Calculating hetLR with gf=%f  theta=%f\n", gfreq,fixed_theta);
  if(fpIR !=NULL)
    dk_curModel.dgf = gfreq;

  if(modelOptions.mapFlag == SS){
    thetaM = fixed_thetaM;
    thetaF = fixed_thetaF;
  }else{
    thetaM = fixed_theta;
    thetaF = fixed_theta;
  }



  if (1 && modelOptions.markerAnalysis == FALSE) {
    pLocus->pAlleleFrequency[0] = gfreq;
    pLocus->pAlleleFrequency[1] = 1 - gfreq;

    if (modelOptions.polynomial == TRUE);
    else
      update_locus (&pedigreeSet, loc1);

  }

  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    status = setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprimeIdx);
    if(status<0)
      KASSERT (1,"Haplotype frequency combination impossible. Exiting!\n");
  }

  if (modelOptions.markerAnalysis == FALSE
      && pLocus1->locusType == LOCUS_TYPE_TRAIT) {
    for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
      pen_DD = x[pen_size * liabIdx + 1];
      pen_Dd = x[pen_size * liabIdx + 2] * x[pen_size * liabIdx + 1];

      if(modelOptions.imprintingFlag){
        pen_dD = x[pen_size * liabIdx + 3] * x[pen_size * liabIdx + 1];
        pen_dd = x[pen_size * liabIdx + 4] * x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 2] *x[pen_size * liabIdx + 3];
      }else{
        pen_dd = x[pen_size * liabIdx + 3] * x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 2];
        pen_dD= pen_Dd;
      }
      if(fpIR !=NULL){
        dk_curModel.pen[liabIdx].DD = pen_DD;
        dk_curModel.pen[liabIdx].Dd = pen_Dd;
        dk_curModel.pen[liabIdx].dD = pen_dD;
        dk_curModel.pen[liabIdx].dd = pen_dd;
      }

      pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
      pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
      pTrait->penetrance[2][liabIdx][1][0] = pen_dD;
      pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
      pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
      pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
      pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_dD;
      pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;
    }
  }
  if (modelOptions.polynomial == TRUE);
  else
    update_penetrance (&pedigreeSet, traitLocus);

  /* get the likelihood at 0.5 first and LD=0 */
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    status = setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprime0Idx);
    if(status<0)
      KASSERT (1,"Haplotype frequency combination impossible. Exiting!\n");

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


  //  KLOG (LOGLIKELIHOOD, LOGDEBUG, "NULL Likelihood\n");
  sprintf (partialPolynomialFunctionName, "TD_C%d_P%%s_%s_%s",
	   pLocus2->pMapUnit->chromosome, 
	   pLocus1->sName, pLocus2->sName);
  cL[7]++;
  ret=compute_likelihood (&pedigreeSet);

  //printf("likelihood =%15.13f with theta 0.5 with %d pedigrees\n", pedigreeSet.likelihood, pedigreeSet.numPedigree);

  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    fprintf (stderr, "Theta 0.5 has likelihood 0\n");
    fprintf (stderr, "dgf=%f\n", gfreq);
    //exit (EXIT_FAILURE);
    f[0] = 1.0;
    return ;
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
      locusList->pNextLocusDistance[k][0] = thetaM;
      locusList->pPrevLocusDistance[k][1] = thetaF;
    }
  } else {
    locusList->pNextLocusDistance[MAP_MALE][0] =
      locusList->pPrevLocusDistance[MAP_MALE][1] =   thetaM;
    locusList->pNextLocusDistance[MAP_FEMALE][0] =
      locusList->pPrevLocusDistance[MAP_FEMALE][1] = thetaF;
  }

  if (modelOptions.polynomial == TRUE);
  else
    /* populate the matrix */
    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */

  //  KLOG (LOGLIKELIHOOD, LOGDEBUG, "ALT Likelihood\n");
  // No new name for a polynomial here because we're reusing the existing one
  cL[8]++;
  ret=compute_likelihood (&pedigreeSet);
  if (ret==-2){
    fprintf(stderr, "Negative likelihood! Exiting...\n");
    exit(EXIT_FAILURE);
  }
  log10_likelihood_alternative = pedigreeSet.log10Likelihood;


  //printf("likelihood =%15.13f with theta %f  %d pedigree\n", pedigreeSet.likelihood,fixed_theta, pedigreeSet.numPedigree);                                  
  if (pedigreeSet.likelihood == 0.0
      && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
    avg_hetLR=0.0;
  } else {
    log10_likelihood_ratio =
      log10_likelihood_alternative - log10_likelihood_null;

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
  for (j = 0; j < 5; j++) {
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
      // fprintf(stderr, "tmp=%15.10f cout=%d   log10= 15.10f \n",tmp, pPedigree->pCount[loc2], log10HetLR);
      // log10HetLR += tmp * sqrt(pPedigree->pCount[loc2] );  //kelvin log10 exponential for case-control
    }

    if(fpIR !=NULL){
      dk_curModel.alpha = alphaV;
      fprintf(fpIR,"%6.3f", log10HetLR);
      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
        fprintf(fpIR," %4.3f", dk_curModel.dprime[0]);
      }
      if(modelOptions.mapFlag == SA){
        fprintf(fpIR," %4.3f", dk_curModel.theta[0]);
      } else{
        fprintf(fpIR," %4.3f %4.3f", dk_curModel.theta[0], dk_curModel.theta[0] );
      }
      fprintf(fpIR," %4.3f %4.3f",dk_curModel.alpha,dk_curModel.dgf);
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
        fprintf(fpIR," %4.3f %4.3f",dk_curModel.pen[liabIdx].DD,dk_curModel.pen[liabIdx].Dd);
        if(modelOptions.imprintingFlag){
          fprintf(fpIR," %4.3f %4.3f",dk_curModel.pen[liabIdx].dD,dk_curModel.pen[liabIdx].dd);
        }else{
          fprintf(fpIR," %4.3f",dk_curModel.pen[liabIdx].dd); 
        }
      }
      fprintf(fpIR," %d\n",dk_curModel.posIdx);
    }


    oldsum=alpha_integral;
    oldscale= (*scale);
    newscale = 0;

    if (log10HetLR >= DBL_MAX_10_EXP - 1) {
      /* find the new scale, adjust the current sum */
      newscale=log10HetLR - (DBL_MAX_10_EXP - SCALE_RESERVE);
    }
    if(newscale > oldscale) {
      /* need to use the newscale and adjust the old sum */
      if(oldsum>0) {
	oldsum_log10=log10(oldsum);
	newsum_log10=oldsum_log10+oldscale-newscale;
	if(newsum_log10<= DBL_MIN_10_EXP+1) {
	  alpha_integral=0;
	}
	else {
	  alpha_integral =pow(10, newsum_log10);
	}
      }
      *scale=newscale;
      oldscale=newscale;
    }else {
      /* use the old scale to adjust the new value */
      newscale = oldscale;
    }
    newLog10HetLR = log10HetLR - newscale; 


    if (newLog10HetLR <= DBL_MIN_10_EXP + 1) {
      hetLR = 0;
    } else {
      hetLR = pow (10, newLog10HetLR);
    }

    alpha_integral += hetLR * alpha[j][1];

    if (log10HetLR > localMOD) {
      localMOD = log10HetLR;

      localmax_x[0] = gfreq;
      localmax_x[1] = alphaV;

      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	localmax_x[pen_size * liabIdx + 2] = x[pen_size * liabIdx + 1];
	localmax_x[pen_size * liabIdx + 3] = x[pen_size * liabIdx + 2] * x[pen_size * liabIdx + 1];
	localmax_x[pen_size * liabIdx + 4] =
	  x[pen_size * liabIdx + 3] * x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 2];

        if(modelOptions.imprintingFlag){
   	  localmax_x[pen_size * liabIdx + 4] = x[pen_size * liabIdx + 3] * x[pen_size * liabIdx + 1];
	  localmax_x[pen_size * liabIdx + 5] =
	    x[pen_size * liabIdx + 4] * x[pen_size * liabIdx + 3] * x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 2];
	}
      }
    }
    // fprintf(fphlod,"%f %f %f %f %f %f %f\n", log10(hetLR*x[1]*x[1]*x[2]), gfreq, pen_DD,pen_Dd, pen_dd, alphaV,fixed_theta);
  }				//end of calculating the HET         


  avg_hetLR = alpha_integral;

  //printf("avg hetLR =%15.10f with gf=%f DD=%f Dd=%f dd=%f theta=%f\n", avg_hetLR, gfreq, pen_DD,pen_Dd, pen_dd, fixed_theta);

  // Jacobian
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    if(modelOptions.imprintingFlag)
      avg_hetLR *= x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 1]* x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 2]* x[pen_size * liabIdx + 3];
    else
      avg_hetLR *= x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 1] * x[pen_size * liabIdx + 2];
  }

  }


  f[0] = avg_hetLR;


}
