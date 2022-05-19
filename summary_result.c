/* Copyright (C) 2009, 2022 Mathematical Medicine LLC
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
#include <float.h> /* defines gcc DBL_MAX_EXP etc. */
#include <math.h> /* log10, pow */
#include <string.h> /* memcpy */

#include "utils/utils.h"
#include "summary_result.h" /* SUMMARY_STAT etc. */
#include "kelvinGlobals.h"
#include "pedlib/pedlib.h" /* for PedigreeSet, Pedigree - hetLR calculations */

int maxscale = 0; /* scale in terms of log10, scale of 0 means 1, scale of 1 means 10 */
/* three dimensional array for the two point summary results *
 * first dimension is the D prime, for LE, D prime=0 with just one element
 * in this dimension 
 * second dimension is theta values 
 * third dimension is marker allele frequency, for LE, only one element in this dimension */
SUMMARY_STAT ***tp_result;
/** One dimensional array, indexing by map position.  For multipoint,
 we don't know how to incorporate LD yet.  This map could be sex
 specific map or sex averaged map. For two point, we don't have to
 distinguish sex specific/avearge as we use theta relative to marker
 during analysis and after analysis (result) */
SUMMARY_STAT *mp_result;

/* storage for the NULL likelihood for the multipoint calculation under polynomial */
//double *****likelihoodQT = NULL;
//double **likelihoodDT = NULL;

/* allocate one extra D prime to store average LR for all D primes 
 * allocate one extra marker allele frequency to store average LR per (Dprime, Theta) pair */
int initialize_tp_result_storage ()
{
  int i, j, k;
  int num;

  CALCHOKE(tp_result, (size_t) pLambdaCell->ndprime + 1, sizeof (SUMMARY_STAT **), SUMMARY_STAT ***);
  for (i = 0; i < pLambdaCell->ndprime + 1; i++) {
    CALCHOKE(tp_result[i], (size_t) modelRange->ntheta + 1, sizeof (SUMMARY_STAT *), SUMMARY_STAT **);
    for (j = 0; j < modelRange->ntheta + 1; j++) {
      num = modelRange->nafreq + 1;
      CALCHOKE(tp_result[i][j], (size_t) num, sizeof (SUMMARY_STAT), SUMMARY_STAT *);
      for (k = 0; k < num; k++) {
	tp_result[i][j][k].max_penIdx = -1;
	tp_result[i][j][k].scale=0;
      }
    }
  }
  return 0;
}

int free_tp_result_storage ()
{
  int i, j;

  for (i = 0; i < pLambdaCell->ndprime + 1; i++) {
    for (j = 0; j < modelRange->ntheta + 1; j++) {
      free (tp_result[i][j]);
    }
    free (tp_result[i]);
  }
  free (tp_result);
  tp_result = NULL;
  return 0;
}

void initialize_max_scale() {
  maxscale = 0;
}

int record_tp_result(int callStatus, PedigreeSet *pedigreeSet, ParamStruct *param, int loc2) 
{
  int newscale, oldscale;
  double newLog10HetLR;
  double oldsum;
  double oldsum_log10;
  double newsum_log10;
  double homoLR, hetLR, log10HetLR;
  int dprimeIdx, thetaInd, mkrFreqIdx;
  int j,liabIdx;
  double alphaV, alphaV2;
  int pedIdx;
  Pedigree *pPedigree;
  double tmp;

  dprimeIdx = param->dprimeIdx;
  thetaInd = param->thetaIdx;
  mkrFreqIdx = param->mkrFreqIdx;

  if(callStatus == -2) {
    fprintf(stderr, "Negative likelihood! Exiting...\n");
    exit(EXIT_FAILURE);
  }

  if(callStatus == -1) {
    /* one of the pedigree has likelihood of 0 (under ALT)
     * this could be absolute 0 likelihood, or close to 0 with underflow problem 
     * either case for now, no need to calculate hetLR, as they all should be 0 */
  }
  else {
    /* caculating the HET */
    for (j = 0; (j==0 && modelOptions->markerAnalysis!=FALSE) || j < modelRange->nalpha; j++) {
      if(modelOptions->markerAnalysis != FALSE)
	alphaV = 1;
      else
	alphaV = modelRange->alpha[j];
      alphaV2 = 1 - alphaV;
      if (alphaV2 < 0)
	alphaV2 = 0;
      log10HetLR = 0;
      for (pedIdx = 0; pedIdx < pedigreeSet->numPedigree; pedIdx++) {
	pPedigree = pedigreeSet->ppPedigreeSet[pedIdx];
	homoLR = pPedigree->likelihood / pedigreeSet->nullLikelihood[pedIdx];
	tmp = log10 (alphaV * homoLR + (1 - alphaV));
	log10HetLR += tmp * pPedigree->pCount[loc2]; // Use the pedigree weight from count file (CF)
      }

      if(fpIR !=NULL){
        dk_curModel.alpha = alphaV;
        fprintf(fpIR,"%6.3f", log10HetLR);
        if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
          fprintf(fpIR," %4.3f", dk_curModel.dprime[0]);
        }
        if(modelOptions->mapFlag == SA){
          fprintf(fpIR," %4.3f", dk_curModel.theta[0]);
        } else{
          fprintf(fpIR," %4.3f %4.3f", dk_curModel.theta[0], dk_curModel.theta[1] );
        }
        fprintf(fpIR," %4.3f %4.3f",dk_curModel.alpha,dk_curModel.dgf);
        for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++) {
          fprintf(fpIR," %4.3f %4.3f",dk_curModel.pen[liabIdx].DD,dk_curModel.pen[liabIdx].Dd);
          if(modelOptions->imprintingFlag){
            fprintf(fpIR," %4.3f %4.3f",dk_curModel.pen[liabIdx].dD,dk_curModel.pen[liabIdx].dd);
          }else{
            fprintf(fpIR," %4.3f",dk_curModel.pen[liabIdx].dd); 
          }
          if (modelType->trait != DICHOTOMOUS && modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
            fprintf(fpIR," %4.3f",dk_curModel.pen[liabIdx].DDSD); 
	  }
        }
        if (modelType->trait == CT){
          fprintf(fpIR," %4.3f",dk_curModel.pen[0].threshold);  /* If each LC uses different threshold, this does not work*/
        }
        fprintf(fpIR," %d\n",dk_curModel.posIdx );
      }

      oldsum=tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total;
      oldscale=tp_result[dprimeIdx][thetaInd][mkrFreqIdx].scale;
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
	    tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total=0;
	  }
	  else {
	    tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total=pow(10, newsum_log10);
	  }
	}
	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].scale=newscale;
	oldscale=newscale;
      }
      else {
	/* use the old scale to adjust the new value */
	newscale = oldscale;
      }
      newLog10HetLR = log10HetLR - newscale; 
      if (newLog10HetLR <= DBL_MIN_10_EXP + 1) {
	hetLR = 0;
      } else {
	hetLR = pow(10, newLog10HetLR);
      }
      /* keep track of the maximum scale we have used so far 
       * we will need to rescale everything once we are done with 
       * populating tp_result */
      if(maxscale < newscale) {
	maxscale = newscale;
      }
      tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total += hetLR;
      if (tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_penIdx < 0 || 
	  log10HetLR > tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_log10_lr) {
	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_log10_lr = log10HetLR;
	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_alpha = alphaV;
	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_gfreq = param->gfreq;
	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_penIdx = param->penIdx;
	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].R_square = param->R_square;
	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_mf = param->mkrFreq;
	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_paramIdx = param->paramIdx;
	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_thresholdIdx = param->thresholdIdx;
      }
    } /* end of looping alpha */
  } /* likelihood of ALT is not 0 */
  tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_count++;
  return 0;
}

/* rescale tp results with the max scale per (D', theta) pair
 */
void rescale_tp_result(int maxScale)
{
  int dprimeIdx, thetaIdx;
  double hetLR, log10HetLR;
  double newHetLR, newLog10HetLR;
  int oldscale;
  double avgLR, log10AvgLR, newAvgLR, newLog10AvgLR;

  if(maxScale==-1)
    maxScale=maxscale;
  for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
    for (thetaIdx = 0; thetaIdx < modelRange->ntheta; thetaIdx++) {
      oldscale = tp_result[dprimeIdx][thetaIdx][modelRange->nafreq].scale;
      hetLR = tp_result[dprimeIdx][thetaIdx][modelRange->nafreq].het_lr_total;
      avgLR = tp_result[dprimeIdx][thetaIdx][modelRange->nafreq].het_lr_avg;
      if(oldscale != maxScale && hetLR > 0) {
	log10HetLR = log10(hetLR);
	newLog10HetLR = log10HetLR + oldscale - maxScale;
	if (newLog10HetLR < DBL_MIN_10_EXP+1) {
	  newHetLR = 0;
	}
	else {
	  newHetLR = pow(10, newLog10HetLR);
	}
	tp_result[dprimeIdx][thetaIdx][modelRange->nafreq].scale=maxScale;
	tp_result[dprimeIdx][thetaIdx][modelRange->nafreq].het_lr_total=newHetLR;
      } /* end of change needed */
      if(oldscale != maxScale && avgLR > 0) {
	log10AvgLR = log10(avgLR);
	newLog10AvgLR = log10AvgLR + oldscale - maxScale;
	if (newLog10AvgLR < DBL_MIN_10_EXP+1) {
	  newAvgLR = 0;
	}
	else {
	  newAvgLR = pow(10, newLog10AvgLR);
	}
	tp_result[dprimeIdx][thetaIdx][modelRange->nafreq].het_lr_avg=newAvgLR;
	tp_result[dprimeIdx][thetaIdx][modelRange->nafreq].scale=maxScale;
      } /* end of change needed */
    } /* end of theta loop */
  } /* end of D' loop */
}


void rescale_tp_result_dprime0(int dprime0Idx) 
{
  int maxScale, oldscale;
  int thetaIdx;
  int changeFlag = 0;
  double avgLR, log10AvgLR, newLog10AvgLR, newAvgLR; 

  /* find the max scale across thetas */
  maxScale = 0;
  for (thetaIdx = 0; thetaIdx < modelRange->ntheta; thetaIdx++) {
    oldscale = tp_result[dprime0Idx][thetaIdx][modelRange->nafreq].scale_orig;
    if(oldscale > maxScale) {
      maxScale = oldscale;
      if(changeFlag ==0 && thetaIdx > 0)
	changeFlag = 1;
    }
    tp_result[dprime0Idx][thetaIdx][modelRange->nafreq].scale_orig2 = oldscale;
    tp_result[dprime0Idx][thetaIdx][modelRange->nafreq].het_lr_avg_orig2 =
      tp_result[dprime0Idx][thetaIdx][modelRange->nafreq].het_lr_avg_orig;
  }

  
  if(changeFlag > 0) {
    for (thetaIdx = 0; thetaIdx < modelRange->ntheta; thetaIdx++) {
      if(oldscale != maxScale) {
	oldscale = tp_result[dprime0Idx][thetaIdx][modelRange->nafreq].scale_orig;
	avgLR = tp_result[dprime0Idx][thetaIdx][modelRange->nafreq].het_lr_avg_orig;
	log10AvgLR = log10(avgLR);
	newLog10AvgLR = log10AvgLR + oldscale - maxScale;
	if(newLog10AvgLR < DBL_MIN_10_EXP+1) {
	  newAvgLR = 0;
	}
	else {
	  newAvgLR = pow(10, newLog10AvgLR);
	}
	tp_result[dprime0Idx][thetaIdx][modelRange->nafreq].het_lr_avg_orig2 = newAvgLR ;
	tp_result[dprime0Idx][thetaIdx][modelRange->nafreq].scale_orig2 = maxScale;
      }
    }
  }
  
}

int record_mp_result(int callStatus, PedigreeSet *pedigreeSet, ParamStruct *param, int posIdx) 
{
  int newscale, oldscale;
  double newLog10HetLR;
  double oldsum;
  double oldsum_log10;
  double newsum_log10;
  double homoLR, hetLR, log10HetLR;
  int gfreqIdx, penIdx, paramIdx, thresholdIdx;
  int j,liabIdx;
  double alphaV, alphaV2;
  int pedIdx;
  Pedigree *pPedigree;
  double tmp;

  gfreqIdx = param->gfreqIdx;
  penIdx = param->penIdx;
  paramIdx = param->paramIdx;
  thresholdIdx = param->thresholdIdx;
  
  if(callStatus == -2) {
    fprintf(stderr, "Negative likelihood! Exiting...\n");
    exit(EXIT_FAILURE);
  }

  if(callStatus == -1) {
    /* one of the pedigree has likelihood of 0 (under ALT)
     * this could be absolute 0 likelihood, or close to 0 with underflow problem 
     * either case for now, no need to calculate hetLR, as they all should be 0 */
  }
  else {
    /* caculating the HET */
    for (j = 0; j < modelRange->nalpha; j++) {
      alphaV = modelRange->alpha[j];
      alphaV2 = 1 - alphaV;
      if (alphaV2 < 0)
	alphaV2 = 0;
      log10HetLR = 0;
      for (pedIdx = 0; pedIdx < pedigreeSet->numPedigree; pedIdx++) {
	pPedigree = pedigreeSet->ppPedigreeSet[pedIdx];
	if(modelType->trait == DT) {
	  homoLR = pPedigree->alternativeLikelihoodDT[gfreqIdx][penIdx] /
	    (pPedigree->traitLikelihoodDT[gfreqIdx][penIdx] *
	     pPedigree->markerLikelihood);
	}
	else {
	  homoLR = pPedigree->alternativeLikelihoodQT[gfreqIdx][penIdx][paramIdx][thresholdIdx] /
	    (pPedigree->traitLikelihoodQT[gfreqIdx][penIdx][paramIdx][thresholdIdx] *
	     pPedigree->markerLikelihood);
	}
	tmp = log10 (alphaV * homoLR + (1 - alphaV));
	log10HetLR += tmp * pPedigree->pCount[loc2]; // Use the pedigree weight from count file (CF)
      }

      if(fpIR !=NULL){
        dk_curModel.alpha = alphaV;
        fprintf(fpIR,"%6.3f", log10HetLR);

        fprintf(fpIR," %4.3f %4.3f",dk_curModel.alpha,dk_curModel.dgf);
        for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++) {
          fprintf(fpIR," %4.3f %4.3f",dk_curModel.pen[liabIdx].DD,dk_curModel.pen[liabIdx].Dd);
          if(modelOptions->imprintingFlag){
            fprintf(fpIR," %4.3f %4.3f",dk_curModel.pen[liabIdx].dD,dk_curModel.pen[liabIdx].dd);
          }else{
            fprintf(fpIR," %4.3f",dk_curModel.pen[liabIdx].dd); 
          }
          if (modelType->trait != DICHOTOMOUS && modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
            fprintf(fpIR," %4.3f",dk_curModel.pen[liabIdx].DDSD); 
	  }
        }
        if (modelType->trait == CT){
          fprintf(fpIR," %4.3f",dk_curModel.pen[0].threshold);  /* If each LC uses different threshold, this does not work*/
        }
        fprintf(fpIR," %d\n",dk_curModel.posIdx);
      }


      oldsum=mp_result[posIdx].het_lr_total;
      oldscale=mp_result[posIdx].scale;
      newscale=0;
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
	    mp_result[posIdx].het_lr_total=0;
	  }
	  else {
	    mp_result[posIdx].het_lr_total=pow(10, newsum_log10); 
	  }
	}
	mp_result[posIdx].scale=newscale;
	oldscale = newscale;
      }
      else {
	newscale = oldscale;
      }
      newLog10HetLR = log10HetLR - newscale; 
      if (newLog10HetLR <= DBL_MIN_10_EXP + 1) {
	hetLR = 0;
      } else {
	hetLR = pow (10, newLog10HetLR);
      }
      /* keep track of the maximum scale we have used so far 
       * we will need to rescale everything once we are done with 
       * populating tp_result */
      if(maxscale < newscale) {
	maxscale = newscale;
      }
      mp_result[posIdx].het_lr_total += hetLR;
      if (mp_result[posIdx].max_penIdx < 0 || 
	  log10HetLR > mp_result[posIdx].max_log10_lr) {
	mp_result[posIdx].max_log10_lr = log10HetLR;
	mp_result[posIdx].max_alpha = alphaV;
	mp_result[posIdx].max_gfreq = param->gfreq;
	mp_result[posIdx].max_penIdx = param->penIdx;
	mp_result[posIdx].max_paramIdx = param->paramIdx;
	mp_result[posIdx].max_thresholdIdx = param->thresholdIdx;
      }
    } /* end of looping alpha */
  } /* likelihood of ALT is not 0 */
  mp_result[posIdx].lr_count++;
  return 0;
}

/* per (Dprime, Theta) pair, calculate the average and max LR by
 * integrating out trait parameters first, then average over marker allele frequencies
 */
int get_average_LR (SUMMARY_STAT *** result)
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
  double lr;
  int max_theta = -1;
  int max_max_theta = -1;
  int max_max_dprime = -1;
  double maxLR;
  int maxScale, scale; 
  double newLog10LR, newHetLR, log10LR;

  for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
    dprime = pLambdaCell->lambda[dprimeIdx][0][0];
    for (thetaInd = 0; thetaInd < modelRange->ntheta; thetaInd++) {
      /* sex averaged theta */
      total_lr = 0;
      count = 0;
      if(modelRange->nafreq > 1) {
	/* find out the max scale for (D', theta) pair */
	maxScale = 0;
	for (mkrFreqIdx = 0; mkrFreqIdx < modelRange->nafreq; mkrFreqIdx++) {
	  scale = tp_result[dprimeIdx][thetaInd][mkrFreqIdx].scale;
	  if(scale > maxScale)
	    maxScale = scale; 
	}
	for (mkrFreqIdx = 0; mkrFreqIdx < modelRange->nafreq; mkrFreqIdx++) {
	  scale = tp_result[dprimeIdx][thetaInd][mkrFreqIdx].scale;
	  lr = tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total;
	  if(scale != maxScale && lr > 0) {
	    /* rescale */
	    log10LR = log10(lr);
	    newLog10LR = log10LR + scale - maxScale;
	    if(newLog10LR < DBL_MIN_10_EXP+1)
	      newHetLR = 0;
	    else
	      newHetLR = pow(10, newLog10LR); 
	    tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total = newHetLR;
	    tp_result[dprimeIdx][thetaInd][mkrFreqIdx].scale = maxScale;
	  }
	}
      }
      for (mkrFreqIdx = 0; mkrFreqIdx < modelRange->nafreq; mkrFreqIdx++) {
	/* het_lr_total has trait parameters already integrated out - get the average */
	//              if(modelType->trait == DT || modelType->distrib != QT_FUNCTION_CHI_SQUARE)
	if (isnan (tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total)) {
	  fprintf (stderr,
		   "het_lr_total is NAN at theta (%f,%f)(%d) and dprime %f(%d).\n",
		   modelRange->theta[0][thetaInd],
		   modelRange->theta[1][thetaInd], thetaInd, dprime,
		   dprimeIdx);
	} else
	  if (isnan (tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total -
 		     tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total))
	{
	  fprintf (stderr,
		   "het_lr_total is INF at theta (%f,%f)(%d) and dprime %f(%d).\n",
		   modelRange->theta[0][thetaInd],
		   modelRange->theta[1][thetaInd], thetaInd, dprime,
		   dprimeIdx);
	}
	if(modelOptions->markerAnalysis == FALSE)
	  lr = tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total /
	    (tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_count *
	     modelRange->nalpha);
	else
	  lr = tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_total /
	    tp_result[dprimeIdx][thetaInd][mkrFreqIdx].lr_count;

	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_avg = lr;
	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].het_lr_avg_orig = lr;
	tp_result[dprimeIdx][thetaInd][mkrFreqIdx].scale_orig = 
	  tp_result[dprimeIdx][thetaInd][mkrFreqIdx].scale; 
	/* add up BR for each marker allele frequency */
	total_lr += lr;
	count++;
	/* keep track of the MOD for each marker allele frequency */
	maxLR = tp_result[dprimeIdx][thetaInd][mkrFreqIdx].max_log10_lr;
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
	  /* max theta would be in tp_result[dprimeIdx][#theta][#mf] */
	}
	/* find overal max MOD  */
	if ((thetaInd == 0 && mkrFreqIdx == 0 && dprimeIdx == 0) ||
	    (maxLR > max_max_lr)) {
	  max_max_lr = maxLR;
	  /* max theta would be in tp_result[max_max_dprime][#theta][#mf] */
	  max_max_dprime = dprimeIdx;
	}

      }				/* end of looping marker allele frequencies */
      /* recording the average and max */
      memcpy (&tp_result[dprimeIdx][thetaInd][modelRange->nafreq],
	      &tp_result[dprimeIdx][thetaInd][max_max_mf],
	      sizeof (SUMMARY_STAT));
      tp_result[dprimeIdx][thetaInd][modelRange->nafreq].max_log10_lr =
	max_max_lr_dprime_theta;
      tp_result[dprimeIdx][thetaInd][modelRange->nafreq].max_mf =
	modelRange->afreq[max_max_mf];
      /* this is the BR after integrating out marker allele frequencies */
      tp_result[dprimeIdx][thetaInd][modelRange->nafreq].het_lr_avg =
	total_lr / count;
      tp_result[dprimeIdx][thetaInd][modelRange->nafreq].max_br_lr =
	max_lr_dprime_theta;
      tp_result[dprimeIdx][thetaInd][modelRange->nafreq].max_br_mf =
	modelRange->afreq[max_mf];
      tp_result[dprimeIdx][thetaInd][modelRange->nafreq].het_lr_avg_orig =
	total_lr/count; 
      tp_result[dprimeIdx][thetaInd][modelRange->nafreq].scale_orig =
	tp_result[dprimeIdx][thetaInd][modelRange->nafreq].scale;


    }				/* end of theta loop */
    /* recording the max per dprime */
    memcpy (&tp_result[dprimeIdx][modelRange->ntheta][modelRange->nafreq],
	    &tp_result[dprimeIdx][max_max_theta][max_max_mf],
	    sizeof (SUMMARY_STAT));
    tp_result[dprimeIdx][modelRange->ntheta][modelRange->nafreq].max_br_lr =
      max_lr_dprime;
    tp_result[dprimeIdx][modelRange->ntheta][modelRange->nafreq].max_br_mf =
      tp_result[dprimeIdx][max_theta][modelRange->nafreq].max_mf;
    tp_result[dprimeIdx][modelRange->ntheta][modelRange->nafreq].
      max_br_theta = modelRange->theta[0][thetaInd];
  }				/* end of d prime loop */

  /* record the max overal */
  memcpy (&tp_result[modelRange->ndprime][modelRange->ntheta]
	  [modelRange->nafreq],
	  &tp_result[max_max_dprime][max_max_theta][max_max_mf],
	  sizeof (SUMMARY_STAT));
  tp_result[dprimeIdx][modelRange->ntheta][modelRange->nafreq].max_br_lr =
    max_lr;
  tp_result[dprimeIdx][modelRange->ntheta][modelRange->nafreq].max_br_dprime =
    dprimeIdx;

  return 0;
}


