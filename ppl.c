
/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#include "kelvin.h"

/* Calculate PPL - only for two point
 * Input:
 *   modelOptions.thetaCutoff r - (0, r) will have modelOptions.thetaWeight
 *                   while (r, 0.5) will have 1 - modelOptions.thetaWeight
 *
 */
double
calculate_PPL (SUMMARY_STAT ** result)
{
  int i, j;
  double w1, w2;
  double integral;
  double a;
  double theta1, theta2;
  double avgLR1, avgLR2;
  double thetam1, thetam2, thetaf1, thetaf2;
  double avgLR3, avgLR4;
  int numTheta = modelRange.ntheta;
  double PPL = 0;
  int num_thetam, num_thetaf;
  double rm, rf;
  double mavg, favg;		/* average thetas */
  double mdist, fdist;		/* distance to the cutoffs */
  double mdiff, fdiff;
  double I = 0;			/* integral for the entire region */
  double A = 0;			/* integral for current theta block */
  double B;			/* area of the theta block */
  int caseIdx = 0;

  integral = 0;
  if (modelOptions.mapFlag == SA) {
    w1 = modelOptions.thetaWeight / modelOptions.thetaCutoff[0];
    w2 =
      (1.0 - modelOptions.thetaWeight) / (0.5 - modelOptions.thetaCutoff[0]);
    for (i = 0; i < numTheta - 1; i++) {
      /* sex averaged theta */
      theta1 = modelRange.theta[0][i];
      theta2 = modelRange.theta[0][i + 1];
      avgLR1 = result[i][modelRange.nafreq].het_lr_avg;
      avgLR2 = result[i + 1][modelRange.nafreq].het_lr_avg;

      if (theta2 <= modelOptions.thetaCutoff[0]) {
	integral += 0.5 * w1 * (theta2 - theta1) * (avgLR1 + avgLR2);
      } else if (theta1 <= modelOptions.thetaCutoff[0]
		 && modelOptions.thetaCutoff[0] < theta2) {
	a =
	  avgLR1 + (avgLR2 - avgLR1) / (theta2 -
					theta1) *
	  (modelOptions.thetaCutoff[0] - theta1);
	integral +=
	  0.5 * w1 * (modelOptions.thetaCutoff[0] - theta1) * (avgLR1 + a);
	integral +=
	  0.5 * w2 * (theta2 - modelOptions.thetaCutoff[0]) * (avgLR2 + a);
      } else {
	/* modelOptions.thetaCutoff[0] < theta1 */
	integral += 0.5 * w2 * (theta2 - theta1) * (avgLR1 + avgLR2);
      }
    }

  } else {
    /* sex specific theta */
    /* theta cutoff */
    rm = modelOptions.thetaCutoff[0];
    rf = modelOptions.thetaCutoff[1];
    w1 = modelOptions.thetaWeight / (rm * rf);
    w2 = (1.0 - modelOptions.thetaWeight) / (0.5 * 0.5 - rm * rf);
    /* number of female/male thetas */
    num_thetam = modelRange.thetacnt[0];
    num_thetaf = modelRange.thetacnt[1];
    I = 0;
    for (i = 0; i < num_thetam - 1; i++) {
      thetam1 = modelRange.theta[0][i * num_thetam];
      thetam2 = modelRange.theta[0][(i + 1) * num_thetam];
      mavg = 0.5 * (thetam1 + thetam2);
      mdiff = thetam2 - thetam1;
      mdist = rm - thetam1;
      for (j = 0; j < num_thetaf - 1; j++) {
	thetaf1 = modelRange.theta[1][j];
	thetaf2 = modelRange.theta[1][j + 1];
	favg = 0.5 * (thetaf1 + thetaf2);
	fdiff = thetaf2 - thetaf1;
	fdist = rf - thetaf1;
	/* area of the theta block 
	 * it's divided into 4 sub blocks, each subblock is represented by 
	 * avgLR1(upper left), avgLR2 (upper right), avgLR3(lower left) or avgLR4 (lower right) 
	 * each block has an area of 0.25 of B
	 */
	B = mdiff * fdiff;

	avgLR1 = result[i * num_thetam + j][modelRange.nafreq].het_lr_avg;
	avgLR2 = result[i * num_thetam + j + 1][modelRange.nafreq].het_lr_avg;
	avgLR3 =
	  result[(i + 1) * num_thetam + j][modelRange.nafreq].het_lr_avg;
	avgLR4 =
	  result[(i + 1) * num_thetam + j + 1][modelRange.nafreq].het_lr_avg;

	/* divided into 10 cases */
	if (thetam2 < rm && thetaf2 < rf) {
	  /* case #1 */
	  /* entirely under skewed region */
	  caseIdx = 1;
	  A = w1 * 0.25 * B * (avgLR1 + avgLR2 + avgLR3 + avgLR4);
	} else if (rm <= thetam1 || rf <= thetaf1) {
	  /* case #2 */
	  caseIdx = 2;
	  /* entirely outside of skewed region */
	  A = w2 * 0.25 * B * (avgLR1 + avgLR2 + avgLR3 + avgLR4);

	} else if (thetam1 < rm && rm <= mavg && thetaf1 < rf && rf <= favg) {
	  /* case #3 */
	  caseIdx = 3;
	  A = avgLR1 * w1 * mdist * fdist +
	    avgLR1 * w2 * (0.25 * B - mdist * fdist) +
	    0.25 * w2 * B * (avgLR2 + avgLR3 + avgLR4);
	} else if (mavg < rm && rm <= thetam2 && thetaf1 < rf && rf < favg) {
	  /* case #4 */
	  caseIdx = 4;
	  A = avgLR1 * 0.5 * mdiff * (w1 * fdist + w2 * (favg - rf)) +
	    avgLR3 * w1 + (rm - mavg) * fdist +
	    avgLR3 * w2 * (0.25 * B - (rm - mavg) * fdist) +
	    w2 * 0.25 * B * (avgLR2 + avgLR4);
	} else if (thetam1 < rm && rm <= mavg && favg < rf && rf <= thetaf2) {
	  /* case #5 */
	  caseIdx = 5;
	  A = avgLR1 * 0.5 * fdiff * (w1 * mdist + w2 * (mavg - rm)) +
	    avgLR2 * w1 * (rf - favg) * mdist +
	    avgLR2 * w2 * (0.25 * B - (rf - favg) * mdist) +
	    w2 * 0.25 * B * (avgLR3 + avgLR4);
	} else if (mavg < rm && rm <= thetam2 && favg < rf && rf <= thetaf2) {
	  /* case #6 */
	  caseIdx = 6;
	  A = avgLR1 * w1 * 0.25 * B +
	    avgLR2 * 0.5 * mdiff * (w1 * (rf - favg) + w2 * (thetaf2 - rf)) +
	    avgLR3 * 0.5 * fdiff * (w1 * (rm - mavg) + w2 * (thetam2 - rm)) +
	    avgLR4 * w1 * (rm - mavg) * (rf - favg) +
	    avgLR4 * w2 * (0.25 * B - (rm - mavg) * (rf - favg));
	} else if (thetam2 < rm && thetaf1 < rf && rf <= favg) {
	  /* case #7 */
	  caseIdx = 7;
	  A =
	    (avgLR1 + avgLR3) * 0.5 * mdiff * (w1 * fdist +
					       w2 * (favg - rf)) + (avgLR2 +
								    avgLR4) *
	    w2 * 0.25 * B;

	} else if (thetam2 < rm && favg < rf && rf <= thetaf2) {
	  /* case #8 */
	  caseIdx = 8;
	  A = (avgLR1 + avgLR3) * w1 * 0.25 * B +
	    (avgLR2 + avgLR4) * 0.5 * mdiff * (w1 * (rf - favg) +
					       w2 * (thetaf2 - rf));
	} else if (thetaf2 < rf && thetam1 < rm && rm <= mavg) {
	  /* case #9 */
	  caseIdx = 9;
	  A =
	    (avgLR1 + avgLR2) * 0.5 * fdiff * (w1 * mdist +
					       w2 * (mavg - rm)) + (avgLR3 +
								    avgLR4) *
	    w2 * 0.25 * B;
	} else if (thetaf2 < rf && mavg < rm && rm <= thetam2) {
	  /* case #10 */
	  caseIdx = 10;
	  A = (avgLR1 + avgLR2) * w1 * 0.25 * B +
	    (avgLR3 + avgLR4) * 0.5 * fdiff * (w1 * (rm - mavg) +
					       w2 * (thetam2 - rm));
	}

	I += A;
      }				/* end of female theta loop */
    }				/* end of male theat loop */
    integral = I;
  }				/* end of sex specific processing */

  PPL =
    (modelOptions.prior * integral) / ((modelOptions.prior * integral) +
				       (1.0 - modelOptions.prior));
  return PPL;
}


/* Two-point LD analyses only. Loads the LDVals struct for later calculation
 * of the various LD-related PPL statistics.
 */
int
get_LDVals (SUMMARY_STAT ***result, LDVals *ldvals)
{
  int dprimeIdx, dprime0Idx, numdprimes=0;
  int thetaIdx;
  double theta1, theta2;
  double cutoff, weight;
  double lr1, lr2, cutlr;

  ldvals->ld_small_theta = ldvals->ld_big_theta = ldvals->ld_unlinked = 0;
  ldvals->le_small_theta = ldvals->le_big_theta = ldvals->le_unlinked = 0;
  cutoff = modelOptions.thetaCutoff[0];
  weight = modelOptions.thetaWeight;
  
  for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
    if (! (dprime0Idx =
	   isDPrime0 (pLambdaCell->lambda[dprimeIdx], pLambdaCell->m, pLambdaCell->n)))
      numdprimes++;
    
    for (thetaIdx = 1; thetaIdx < modelRange.ntheta; thetaIdx++) {
      theta1 = modelRange.theta[0][thetaIdx-1];
      theta2 = modelRange.theta[0][thetaIdx];
      lr1 = result[dprimeIdx][thetaIdx-1][modelRange.nafreq].het_lr_avg;
      lr2 = result[dprimeIdx][thetaIdx][modelRange.nafreq].het_lr_avg;
      /* isnan check on both LRs? */

      if (dprime0Idx) {
	/* D' == 0 */
	if (theta1 >= (cutoff - ERROR_MARGIN)) {
	  ldvals->le_big_theta += (lr1 + lr2) * (theta2 - theta1);
	} else if (theta2 <= (cutoff + ERROR_MARGIN)) {
	  ldvals->le_small_theta += (lr1 + lr2) * (theta2 - theta1);
	} else {
	  cutlr = lr1 + ((cutoff - theta1) / (theta2 - theta1)) * (lr2 - lr1);
	  ldvals->le_small_theta += (lr1 + cutlr) * (cutoff - theta1);
	  ldvals->le_big_theta += (cutlr + lr2) * (theta2 - cutoff);
	}
	if (theta2 >= (0.5 - ERROR_MARGIN))
	  ldvals->le_unlinked += lr2;
      } else {
	/* D' != 0 */
	if (theta1 >= (cutoff - ERROR_MARGIN)) {
	  ldvals->ld_big_theta += (lr1 + lr2) * (theta2 - theta1);
	} else if (theta2 <= (cutoff + ERROR_MARGIN)) {
	  ldvals->ld_small_theta += (lr1 + lr2) * (theta2 - theta1);
	} else {
	  cutlr = lr1 + ((cutoff - theta1) / (theta2 - theta1)) * (lr2 - lr1);
	  ldvals->ld_small_theta += (lr1 + cutlr) * (cutoff - theta1);
	  ldvals->ld_big_theta += (cutlr + lr2) * (theta2 - cutoff);
	}
	if (theta2 >= (0.5 - ERROR_MARGIN))
	  ldvals->ld_unlinked += lr2;
      }
    }
  }

  /* The 0.5 here is factored out of the calculation of area */
  ldvals->ld_small_theta *= 0.5 * (weight / cutoff);
  ldvals->le_small_theta *= 0.5 * (weight / cutoff);
  ldvals->ld_big_theta *= 0.5 * ((1 - weight) / (0.5 - cutoff));
  ldvals->le_big_theta *= 0.5 * ((1 - weight) / (0.5 - cutoff));

  /* Average LD values over the number of non-zero D's */
  ldvals->ld_small_theta /= numdprimes;
  ldvals->ld_big_theta /= numdprimes;
  ldvals->ld_unlinked /= numdprimes;

  return (0);
}


double calc_ldppl (LDVals *ldvals)
{
  double numerator;
  double denomRight;
  double prior = modelOptions.LDprior;
  double ldppl;

  numerator =
    ldvals->ld_small_theta * prior * 0.021 + 
    ldvals->ld_big_theta * prior *  0.0011+ 
    ldvals->le_small_theta * prior * 0.979 + 
    ldvals->le_big_theta * prior * 0.9989;
  denomRight =
    ldvals->le_unlinked * (1 - prior);
  ldppl = numerator / (numerator + denomRight);
  
  return (ldppl);
}


double calc_ppld_given_linkage (LDVals *ldvals)
{
  double numerator;
  double denomRight;
  double prior = modelOptions.LDprior;
  double ppld_given_l;

  numerator =
    ldvals->ld_small_theta * prior * 0.021 + 
    ldvals->ld_big_theta * prior * 0.0011;
  denomRight =
    ldvals->le_small_theta * prior * 0.979 + 
    ldvals->le_big_theta * prior * 0.9989;
  ppld_given_l = numerator / (numerator + denomRight);

  return (ppld_given_l);
}


double calc_ppld (LDVals *ldvals)
{
  double numerator;
  double denomRight;
  double prior = modelOptions.LDprior;
  double ppld; 

  numerator =
    ldvals->ld_small_theta * prior * 0.021 +
    ldvals->ld_big_theta * prior * 0.0011; 
  denomRight =
    ldvals->le_small_theta * prior * 0.979 + 
    ldvals->le_big_theta * prior * 0.9989 + 
    ldvals->le_unlinked * (1 - prior);
  ppld = numerator / (numerator + denomRight); 

  return (ppld);
}


double calc_ppld_and_linkage (LDVals *ldvals)
{
  double numerator;
  double denomRight;
  double prior = modelOptions.LDprior;
  double ppld_and_l;

  numerator =
    ldvals->ld_small_theta * prior * 0.021 + 
    ldvals->ld_big_theta * prior * 0.0011;
  denomRight =
    ldvals->le_small_theta * prior * 0.979 + 
    ldvals->le_big_theta * prior * 0.9989 + 
    ldvals->le_unlinked * (1 - prior);
  ppld_and_l = numerator / (numerator + denomRight);

  return (ppld_and_l);
}


/* per (Dprime, Theta) pair, calculate the average and max LR by
 * integrating out trait parameters first, then average over marker allele frequencies
 */
int
get_average_LR (SUMMARY_STAT *** result)
{
  int dprimeIdx;
  int thetaInd;
  int mkrFreqIdx;

  //  int maxFreqIdx;
  double total_lr;
  double max_lr_dprime_theta = -9999.99;
  double max_max_lr_dprime_theta = -9999.99;
  long long count;
  int max_mf = -1;
  int max_max_mf = -1;
  double max_lr_dprime = -9999.99;
  double max_max_lr_dprime = -9999.99;
  double max_lr = -9999.99;
  double max_max_lr = -9999.999;
  double dprime;
  double theta;
  double lr;
  int max_theta = -1;
  int max_max_theta = -1;
  int max_dprime = -1;
  int max_max_dprime = -1;
  double maxLR;

  fprintf(stderr, "YGH %e\n", tp_result[0][0][1].het_lr_avg);

  for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
    dprime = pLambdaCell->lambda[dprimeIdx][0][0];
    for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++) {
      /* sex averaged theta */
      theta = modelRange.theta[0][thetaInd];
      total_lr = 0;
      count = 0;
      for (mkrFreqIdx = 0; mkrFreqIdx < modelRange.nafreq; mkrFreqIdx++) {
	/* het_lr_total has trait parameters already integrated out - get the average */
	//              if(modelType.trait == DT || modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	if (isnan (tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total)) {
	  fprintf (stderr,
		   "het_lr_total is NAN at theta (%f,%f)(%d) and dprime %f(%d).\n",
		   modelRange.theta[0][thetaInd],
		   modelRange.theta[1][thetaInd], thetaInd, dprime,
		   dprimeIdx);
	} else
	  if (isnan (tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total -
 		     tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total))
	{
	  fprintf (stderr,
		   "het_lr_total is INF at theta (%f,%f)(%d) and dprime %f(%d).\n",
		   modelRange.theta[0][thetaInd],
		   modelRange.theta[1][thetaInd], thetaInd, dprime,
		   dprimeIdx);
	}
	if (modelType.trait == DT) {
	  if (modelType.ccFlag)
	    lr = tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total;
	  else
	    lr = tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total /
	      (tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_count *
	       modelRange.nalpha);
	} else
	  lr = tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total /
	    (tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_count /
	     modelRange.ntthresh * 2 * (modelType.maxThreshold -
					modelType.minThreshold) *
	     modelRange.nalpha);

	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_avg = lr;
	if(lr < 0) {
	  fprintf(stderr, "het_lr_avg negative %e! mkrFreqIdx: %d\n", lr, mkrFreqIdx);
	  exit(-1);
	}
	/* add up BR for each marker allele frequency */
	total_lr += lr;
	count++;
	/* keep track of the MOD for each marker allele frequency */
	maxLR = tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_lr;
	/* find max BR per (dprime, theta) pair */
	if (mkrFreqIdx == 0 || (lr > max_lr_dprime_theta)) {
	  max_lr_dprime_theta = lr;
	  max_mf = mkrFreqIdx;
	}
	/* find the max MOD per (dprime, theta) pair */
	if (mkrFreqIdx == 0 || (maxLR > max_max_lr_dprime_theta)) {
	  max_max_lr_dprime_theta = maxLR;
	  max_max_mf = mkrFreqIdx;
	}
	/* find max BR per dprime */
	if ((thetaInd == 0 && mkrFreqIdx == 0) || (lr > max_lr_dprime)) {
	  max_lr_dprime = lr;
	  max_theta = thetaInd;
	}
	/* find max MOD per dprime */
	if ((thetaInd == 0 && mkrFreqIdx == 0)
	    || (maxLR > max_max_lr_dprime)) {
	  max_max_lr_dprime = maxLR;
	  max_max_theta = thetaInd;
	}

	/* find overal max BR  */
	if ((thetaInd == 0 && mkrFreqIdx == 0 && dprimeIdx == 0) ||
	    (lr > max_lr)) {
	  max_lr = lr;
	  /* max theta would be in tp_result[max_dprime][#theta][#mf] */
	  max_dprime = dprimeIdx;
	}
	/* find overal max MOD  */
	if ((thetaInd == 0 && mkrFreqIdx == 0 && dprimeIdx == 0) ||
	    (maxLR > max_max_lr)) {
	  max_max_lr = maxLR;
	  /* max theta would be in tp_result[max_dprime][#theta][#mf] */
	  max_max_dprime = dprimeIdx;
	}

      }				/* end of looping marker allele frequencies */

      fprintf(stderr, "YGH %e\n", tp_result[0][0][1].het_lr_avg);

      /* recording the average and max */
      memcpy (&tp_result[dprimeIdx][thetaInd][modelRange.nafreq],
	      &tp_result[dprimeIdx][thetaInd][max_max_mf],
	      sizeof (SUMMARY_STAT));
      fprintf(stderr, "YGH %e\n", tp_result[0][0][1].het_lr_avg);

      tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_lr =
	max_max_lr_dprime_theta;
      tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_mf =
	modelRange.afreq[max_max_mf];
      /* this is the BR after integrating out marker allele frequencies */
      tp_result[dprimeIdx][thetaInd][modelRange.nafreq].het_lr_avg =
	total_lr / count;
      tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_br_lr =
	max_lr_dprime_theta;
      tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_br_mf =
	modelRange.afreq[max_mf];
    }				/* end of theta loop */
    /* recording the max per dprime */
    memcpy (&tp_result[dprimeIdx][modelRange.ntheta][modelRange.nafreq],
	    &tp_result[dprimeIdx][max_max_theta][max_max_mf],
	    sizeof (SUMMARY_STAT));
    tp_result[dprimeIdx][modelRange.ntheta][modelRange.nafreq].max_br_lr =
      max_lr_dprime;
    tp_result[dprimeIdx][modelRange.ntheta][modelRange.nafreq].max_br_mf =
      tp_result[dprimeIdx][max_theta][modelRange.nafreq].max_mf;
    tp_result[dprimeIdx][modelRange.ntheta][modelRange.nafreq].
      max_br_theta = modelRange.theta[0][thetaInd];
  }				/* end of d prime loop */

  /* record the max overal */
  memcpy (&tp_result[modelRange.ndprime][modelRange.ntheta]
	  [modelRange.nafreq],
	  &tp_result[max_max_dprime][max_max_theta][max_max_mf],
	  sizeof (SUMMARY_STAT));
  tp_result[dprimeIdx][modelRange.ntheta][modelRange.nafreq].max_br_lr =
    max_lr;
  tp_result[dprimeIdx][modelRange.ntheta][modelRange.nafreq].max_br_dprime =
    dprimeIdx;

  return 0;
}

/* allocate one extra D prime to store average LR for all D primes 
 * allocate one extra marker allele frequency to store average LR per (Dprime, Theta) pair */
int
initialize_tp_result_storage ()
{
  int i, j, k;
  int num;

  tp_result =
    (SUMMARY_STAT ***) calloc (pLambdaCell->ndprime + 1,
			       sizeof (SUMMARY_STAT **));
  for (i = 0; i < pLambdaCell->ndprime + 1; i++) {
    tp_result[i] =
      (SUMMARY_STAT **) calloc (modelRange.ntheta + 1,
				sizeof (SUMMARY_STAT *));
    for (j = 0; j < modelRange.ntheta + 1; j++) {
      num = modelRange.nafreq + 1;
      tp_result[i][j] =
	(SUMMARY_STAT *) calloc (num, sizeof (SUMMARY_STAT));
      for (k = 0; k < num; k++) {
	tp_result[i][j][k].max_penIdx = -1;
      }
    }
  }
  return 0;
}

int
free_tp_result_storage ()
{
  int i, j;

  for (i = 0; i < pLambdaCell->ndprime + 1; i++) {
    for (j = 0; j < modelRange.ntheta + 1; j++) {
      free (tp_result[i][j]);
    }
    free (tp_result[i]);
  }
  free (tp_result);
  tp_result = NULL;
  return 0;
}

double
calculate_R_square (double p1, double q1, double d)
{
  double p2 = 1.0 - p1;
  double q2 = 1.0 - q1;

  return (d * d / (p1 * p2 * q1 * q2));
}

void
free_likelihood_storage (PedigreeSet * pedSet)
{
  int pedIdx;
  int gfreqInd;
  int penIdx;
  int paramIdx;

  if (likelihoodDT == NULL && likelihoodQT == NULL)
    return;

  if (modelType.trait != DT) {	/*
				   for (pedIdx = 0; pedIdx < pedSet->numPedigree + 1; pedIdx++)
				   {
				   for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++)
				   {
				   * third dimension is penetrance *
				   free (likelihoodDT[pedIdx][gfreqInd]);
				   }
				   * second dimension is gene freq *
				   free (likelihoodDT[pedIdx]);
				   }
				   free (likelihoodDT);

				   }
				   else
       {                       *//* QT */
    for (pedIdx = 0; pedIdx < pedSet->numPedigree + 1; pedIdx++) {
      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++) {

	for (penIdx = 0; penIdx < modelRange.npenet; penIdx++) {
	  /* fourth dimension is SD */
	  for (paramIdx = 0; paramIdx < modelRange.nparam; paramIdx++) {
	    /* 5th dimension is threshold */
	    free (likelihoodQT[pedIdx][gfreqInd][penIdx][paramIdx]);

	  }
	  free (likelihoodQT[pedIdx][gfreqInd][penIdx]);
	}
	/* third dimension is mean */
	free (likelihoodQT[pedIdx][gfreqInd]);
      }
      /* second dimension is gene freq */
      free (likelihoodQT[pedIdx]);
    }
    free (likelihoodQT);
  }
}
