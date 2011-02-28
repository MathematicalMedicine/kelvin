/**
@file integrationSupport.c

  The new dynamic integration half of the former kelvin main program.

  Functions supporting maximized logarithm of odds (MOD) calculation
  via DCUHRE routines for approximation of a vector of definite
  integrals via recursive partitioning into subregions.

  A polynomial representation of the Elston-Stewart likelihood
  calculation for a pedigree structure is built in terms of trait
  space variables (position, frequency, penetrance) population
  admixture (alpha) and sometimes linkage disequilibrium (d').
  Plugging in values for these variables will allow calculation of the
  likelihood of the unvariant pedigree genotype given that set of
  values, and hence, the likelihood of that set of values.
  
  Our primary goal is to determine the peak likelihood of the set of
  pedigrees with the trait at various positions independent of other
  variables, one would think that we wouldn't particularly care about
  the likelihood at any of the other variables, since ultimately all
  we'd want is the peak likelihood for the pedigree at a range of
  trait positions independent of all other variabiles. Actually we do
  care about them since they can give us information about the nature
  of the disease.

  Copyright &copy; 2010, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

*/
#include <float.h>
#include <math.h>
#include <ctype.h>

#include "kelvin.h"
#include "kelvinGlobals.h"
#include "summary_result.h"
#include "dkelvinWriteFiles.h"
#include "kelvinWriteFiles.h"   // Just for writeSurfaceFileHeader
#include "trackProgress.h"
#include "ppl.h"

// This will go into a 'client include' for likelihood in the future.
#define compute_likelihood(pPedigreeList) compute_likelihood(__FILE__, __LINE__, pPedigreeList)
#define populate_xmission_matrix(pMatrix, totalLoci, prob, prob2, hetProb, cellIndex, lastHetLoc, prevPattern, loc) populate_xmission_matrix(__FILE__, __LINE__, pMatrix, totalLoci, prob, prob2, hetProb, cellIndex, lastHetLoc, prevPattern, loc)

#include "dcuhre.h"
#include "integrationGlobals.h"
#include "integrationSupport.h"
#include "integrationLocals.h"

// Application globals - globals that transcend the current file (module)

extern LDLoci *pLDLoci;

// Module-wide (file) globals

#ifdef STUDYDB
#include "database/StudyDB.h"
#include "database/databaseSupport.h"
extern struct StudyDB studyDB;

#define MAXLSTP 4096
double *oldTLoc; // Pointer to original list of trait loci, which we're going to ignore
double lociSetTransitionPositions[MAXLSTP];
double newTLoc[MAXLSTP];

int compare_doubles (const void *a, const void *b) {
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  return (*da > *db) - (*da < *db);
}

#endif

/**

  One-liner.

  More.

  @author Sang-Cheol Seok
  @par Reviewers:
     Bill Valentine-Cooper on 2009-10-30.<br>

  @par Global Inputs

  @par Global Outputs

  @return the new value of secondArg, unless firstArg wasn't positive.
  @retval -1 firstArg was less than zero
  @retval 0 firstArg was zero

*/
int kelvin_dcuhre_integrate (double *integralParam, double *abserrParam, double vol_region, int *scale)
{

  /* INPUT
   * vol_region : the volume of the rectangular region
   * (( s->vol_rate is the rate to convert DCUHRE's output into the average))
   * 
   * OUTPUT
   * integralParam : average funtion value in the given region
   * error    : estimate error in calculation of integral
   * scale    : scale used in calculation of BR
   * 
   * VARIABLE
   * dim : number of variables in the middle layer of the 3-layer approach, which is dim of DCUHRE
   * s   : sturcture to hold all global information in using DCUHRE
   * and the only variable to change with different analysis
   * localMOD : maximum in each BR  
   * maximization information is in localmax_x[]
   * boost_rate : BR boosting BR^boost_rate
   * 
   */


  /* Local variables */
  int dim, return_val, i;
  double boost_rate = 1.3;      //1.1;


  if (modelOptions->equilibrium == LINKAGE_DISEQUILIBRIUM)
    boost_rate = 1.0;
  //extern /* Subroutine */ int ftest_();  

  localMOD = DBL_MIN_10_EXP;

  if (modelType->trait == DICHOTOMOUS) {

    dim = 1 + 3 * modelRange->nlclass;
    if (modelOptions->imprintingFlag)
      dim += modelRange->nlclass;       //dD

    s = &init_state;
    initialize_state (s, xl, xu, dim);

    if (modelType->type == TP) {
      s->funsub = (U_fp) compute_hlod_2p_dt;
      s->mType = TP_DT;
    } else {
      s->funsub = (U_fp) compute_hlod_mp_dt;
      s->mType = MP_DT;
    }
  } else {      /*  QT or combined */

    dim = total_dim - 1;        // alpha
    if (modelType->type == TP) {
      if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
        dim -= 2;       // theta and dprime
      } else {
        if (modelOptions->mapFlag == SA)
          dim -= 1;     //theta
        else
          dim -= 2;
      }
    }

    s = &init_state;
    initialize_state (s, xl, xu, dim);

    if (modelType->type == TP) {
      s->funsub = (U_fp) compute_hlod_2p_qt;
      s->mType = TP_DT;
    } else {
      s->funsub = (U_fp) compute_hlod_mp_qt;
      s->mType = MP_DT;
    }
  }
  if (modelOptions->maxIterations > -1)
    s->maxcls = modelOptions->maxIterations;
  else if (dim < 10)
    s->maxcls = 10000000;       //50000;
  else
    s->maxcls = 100000 * (int) pow (2.0, dim);

  s->nlclass = modelRange->nlclass;
  s->aim_diff_suc = 3 * s->nlclass;

  for (i = 0; i < s->nlclass; i++) {
    if (modelOptions->imprintingFlag)
      s->vol_rate /= 16.0;
    else
      s->vol_rate /= 6.0;
  }
  s->vol_rate *= vol_region;    /*This is the rate to convert to average function value */

  DIAG (DCUHRE, 1, {
        fprintf (stderr, "Starting DCUHRE with dim=%d\n", dim);
      }
  );

  return_val = dcuhre_ (s);
  if (return_val > 0) {
    return return_val;
  }

  s->result /= s->vol_rate;
  s->error /= s->vol_rate;

  DIAG (DCUHRE, 1, {
      fprintf (stderr, "Final result =%15.10f  with error =%15.10f and neval = %d\n",
	       s->result, s->error, s->total_neval);
      fprintf (stderr, "End of DCUHRE with ifail =%d\n", s->ifail);
    }
    );

  /* BR boosting is done here */
  if (modelOptions->equilibrium == LINKAGE_EQUILIBRIUM && modelType->trait == DT) {
    //fprintf(stderr, "Before boosting %e\n", s->result);
    s->result = pow (10.0, (log10 (s->result) * boost_rate));
    //fprintf(stderr, "After boosting %e\n", s->result);
  }

  *integralParam = s->result;
  *abserrParam = s->error;
  *scale = s->scale;
  return return_val;

}


void compute_hlod_mp_qt (double x[], double *f, int *scale)
{

  int k, j;
  int pedIdx, liabIdxLocal = 0, statusLocal, pen_size = 3;
  double constraint = 0.0;
  double mean_DD = 0.0, mean_Dd = 0.0, mean_dD = 0.0, mean_dd = 0.0;
  double SD_DD = 0.0, SD_Dd = 0.0, SD_dD = 0.0, SD_dd = 0.0;
  double gfreq;
  double alphaV;
  double threshold = 0.0;
  double log10_likelihood_null, log10_likelihood_alternative, log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;
  double log10Likelihood;

  /* for null likelihood calculation */
  double product_likelihood = 1;        /* product of the likelihoods for all the pedigrees */
  double sum_log_likelihood = 0;        /* sum of the log10(likelihood) for all the pedigrees */

  Pedigree *pPedigreeLocal;

  int newscale, oldscale;       /* scaling related variables */
  double newLog10HetLR;
  double oldsum;
  double oldsum_log10;
  double newsum_log10;


  if (modelOptions->imprintingFlag)
    pen_size = 4;

  int origLocus = analysisLocusList->pLocusIndex[0];

  if (analysisLocusList->numLocus > 1)
    origLocus = analysisLocusList->pLocusIndex[1];

  gfreq = x[0];
  if (fpIR != NULL)
    dk_curModel.dgf = gfreq;

  pLocus->pAlleleFrequency[0] = gfreq;
  pLocus->pAlleleFrequency[1] = 1 - gfreq;

  update_locus (&pedigreeSet, traitLocus);


  j = 1;        // j=0 for gfrequency
  if (modelType->trait == CT) {
    threshold = x[s->ndim - 1];
  }
  for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
    mean_DD = x[j];
    mean_Dd = (x[j + 1] - xl[j + 1]) * (x[j] - xl[j]) / (xu[j] - xl[j]) + xl[j + 1];

    if (modelOptions->imprintingFlag) {
      mean_dD = (x[j + 2] - xl[j + 2]) * (x[j] - xl[j]) / (xu[j] - xl[j]) + xl[j + 2];
      mean_dd = (x[j + 3] - xl[j + 3]) * (x[j + 2] - xl[j + 2]) / (xu[j + 2] - xl[j + 2]) * (x[j + 1] - xl[j + 1]) / (xu[j + 1] - xl[j + 1]) * (x[j] - xl[j]) / (xu[j] - xl[j]) + xl[j + 3];
    } else {
      mean_dd = (x[j + 2] - xl[j + 2]) * (x[j + 1] - xl[j + 1]) / (xu[j + 1] - xl[j + 1]) * (x[j] - xl[j]) / (xu[j] - xl[j]) + xl[j + 2];
      mean_dD = mean_Dd;
    }
    j += pen_size;

    if (fpIR != NULL) {
      dk_curModel.pen[liabIdxLocal].DD = mean_DD;
      dk_curModel.pen[liabIdxLocal].Dd = mean_Dd;
      dk_curModel.pen[liabIdxLocal].dD = mean_dD;
      dk_curModel.pen[liabIdxLocal].dd = mean_dd;
    }


    if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
      /*SD_DD = x[j];           
       * SD_Dd = x[j + 1];      
       * if(modelOptions->imprintingFlag){
       * SD_dD = x[j+2];
       * SD_dd = x[j+3];
       * }else{
       * SD_dd = x[j+2];        
       * SD_dD= SD_Dd;
       * }
       * j += pen_size; */
      SD_DD = SD_Dd = SD_dD = SD_dd = x[j++];

      if (fpIR != NULL) {
        dk_curModel.pen[liabIdxLocal].DDSD = SD_DD;
        dk_curModel.pen[liabIdxLocal].DdSD = SD_Dd;
        dk_curModel.pen[liabIdxLocal].dDSD = SD_dD;
        dk_curModel.pen[liabIdxLocal].ddSD = SD_dd;
      }

    }
    if (fpIR != NULL) {
      if (modelType->trait == CT)
        dk_curModel.pen[liabIdxLocal].threshold = threshold;
    }

    /* threshold for QT *
     * if (modelType->trait == CT) {
     * threshold = x[j];        // modelRange->tthresh[liabIdxLocal][thresholdIdx];
     * j++;
     * } */

    /* check against the hard coded constraint */
    if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
      constraint = (1.0 - gfreq) * (1.0 - gfreq) * mean_dd * SD_dd + 2 * gfreq * (1.0 - gfreq)
          * mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
      /* fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
       * constraint, gfreq, mean_DD, SD_DD, mean_Dd, SD_DD, mean_dd, SD_dd);
       */
      if (constraint >= 3.0 || constraint <= -3.0) {
        // fprintf(stderr,"Constraint is %f \n", constraint);
        num_out_constraint++;
      }
    }

    pTrait->means[liabIdxLocal][0][0] = mean_DD;
    pTrait->means[liabIdxLocal][0][1] = mean_Dd;
    pTrait->means[liabIdxLocal][1][0] = mean_dD;
    pTrait->means[liabIdxLocal][1][1] = mean_dd;
    pTrait->stddev[liabIdxLocal][0][0] = SD_DD;
    pTrait->stddev[liabIdxLocal][0][1] = SD_Dd;
    pTrait->stddev[liabIdxLocal][1][0] = SD_dD;
    pTrait->stddev[liabIdxLocal][1][1] = SD_dd;

    /* threshold for QT */
    pTrait->cutoffValue[liabIdxLocal] = threshold;

  }     /* liability class Index */

  if (modelOptions->polynomial == TRUE);
  else
    /* only need to update trait locus */
    update_penetrance (&pedigreeSet, traitLocus);

  /* for trait likelihood */
  analysisLocusList = &traitLocusList;
  xmissionMatrix = traitMatrix;

  /* compute the null likelihood with   */
  pedigreeSet.likelihood = 1;
  pedigreeSet.log10Likelihood = 0;

  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {

    /* save the likelihood at null */
    pPedigreeLocal = pedigreeSet.ppPedigreeSet[pedIdx];

    if (modelOptions->polynomial == TRUE) {
      ASSERT (pPedigreeLocal->traitLikelihoodPolynomial != NULL, "Trait likelihood polynomial is NULL");
      /* evaluate likelihood */
      evaluatePoly (pPedigreeLocal->traitLikelihoodPolynomial, pPedigreeLocal->traitLikelihoodPolyList, &pPedigreeLocal->likelihood);
    } else {
      initialize_multi_locus_genotype (pPedigreeLocal);
      statusLocal = compute_pedigree_likelihood (pPedigreeLocal);
    }

    /*pPedigreeLocal->likelihood is now computed and now check it */
    if (pPedigreeLocal->likelihood == 0.0) {
      if (modelRange->atypicalQtTrait)
        WARNING ("Pedigree %s has likelihood of zero or nearly zero", pPedigreeLocal->sPedigreeID);
      else
        ERROR ("Pedigree %s has likelihood of zero or nearly zero", pPedigreeLocal->sPedigreeID);
      product_likelihood = 0.0;
      sum_log_likelihood = -9999.99;

      f[0] = 1.0;
      return;
      //       break;
    } else if (pPedigreeLocal->likelihood < 0.0) {
      ASSERT (pPedigreeLocal->likelihood >= 0.0, "Pedigree %s has negative likelihood", pPedigreeLocal->sPedigreeID);
      product_likelihood = 0.0;
      sum_log_likelihood = -9999.99;
      break;
    } else {
      if (pPedigreeLocal->pCount[origLocus] == 1) {
        product_likelihood *= pPedigreeLocal->likelihood;
        log10Likelihood = log10 (pPedigreeLocal->likelihood);
      } else {
        product_likelihood *= pow (pPedigreeLocal->likelihood, pPedigreeLocal->pCount[origLocus]);
        log10Likelihood = log10 (pPedigreeLocal->likelihood) * pPedigreeLocal->pCount[origLocus];
      }
      sum_log_likelihood += log10Likelihood;
    }
    pedigreeSet.nullLikelihood[pedIdx] = pPedigreeLocal->likelihood;
  }

  pedigreeSet.likelihood = product_likelihood;
  pedigreeSet.log10Likelihood = sum_log_likelihood;
  log10_likelihood_null = pedigreeSet.log10Likelihood;
  DIAG (OVERALL, 1, {
        fprintf (stderr, "Sum of log Likelihood is: %e\n", sum_log_likelihood);
      }
  );

  /* This is for alternative likelihood */
  analysisLocusList = &savedLocusList;
  xmissionMatrix = altMatrix;
  if (modelOptions->polynomial == TRUE);
  else
    statusLocal = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, /* probability */
        initialProbAddr2,       /* probability */
        initialHetProbAddr, 0,  /* cell index */
        -1, -1, /* last het locus & last het pattern (P-1 or M-2) */
        0);     /* current locus - start with 0 */


  char markerNo[8];
  sprintf (partialPolynomialFunctionName, "MQA_LC%d_C%d_P%%sM", modelRange->nlclass, (originalLocusList.ppLocusList[1])->pMapUnit->chromosome);
  for (k = 0; k < modelType->numMarkers; k++) {
    if (traitPos <= *get_map_position (markerLocusList.pLocusIndex[k]) && (strstr (partialPolynomialFunctionName, "_T") == NULL))
      strcat (partialPolynomialFunctionName, "_T");
    sprintf (markerNo, "_%d", markerLocusList.pLocusIndex[k]);
    strcat (partialPolynomialFunctionName, markerNo);
  }
  if (strstr (partialPolynomialFunctionName, "_T") == NULL)
    strcat (partialPolynomialFunctionName, "_T");
  compute_likelihood (&pedigreeSet);
  cL[3]++; // MP QT alternative likelihood

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;
  if (isnan (log10_likelihood_alternative))
    ERROR ("Alternative hypothesis likelihood is not a number");
  if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
    avg_hetLR = 0.0;
  } else {
    log10_likelihood_ratio = log10_likelihood_alternative - log10_likelihood_null - pedigreeSet.log10MarkerLikelihood;
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
      ERROR ("Likelihood ratio for the pedigree set is not a number");

    /* caculating the HET */
    for (j = 0; j < 5; j++) {
      alphaV = alpha[j][0];
      alphaV2 = 1 - alphaV;
      if (alphaV2 < 0)
        alphaV2 = 0;

      log10HetLR = 0;
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        pPedigreeLocal = pedigreeSet.ppPedigreeSet[pedIdx];
        homoLR = pPedigreeLocal->likelihood / (pedigreeSet.nullLikelihood[pedIdx] * pPedigreeLocal->markerLikelihood);
        //fprintf(stderr,"pedIdx=%d alternative=%e trait null=%e marker null=%e\n",pedIdx,pPedigreeLocal->likelihood ,pedigreeSet.nullLikelihood[pedIdx] , pPedigreeLocal->markerLikelihood);

        if (alphaV * homoLR + alphaV2 < 0)
          WARNING ("Heterogeneity likelihood ratio less than zero");
        log10HetLR += log10 (alphaV * homoLR + alphaV2);
      }

      if (fpIR != NULL) {
        dk_curModel.alpha = alphaV;
        fprintf (fpIR, "%6.3f", log10HetLR);

        fprintf (fpIR, " %9.8f %9.8f", dk_curModel.alpha, dk_curModel.dgf);
        for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
          fprintf (fpIR, " %9.8f %9.8f", dk_curModel.pen[liabIdxLocal].DD, dk_curModel.pen[liabIdxLocal].Dd);
          if (modelOptions->imprintingFlag) {
            fprintf (fpIR, " %9.8f %9.8f", dk_curModel.pen[liabIdxLocal].dD, dk_curModel.pen[liabIdxLocal].dd);
          } else {
            fprintf (fpIR, " %9.8f", dk_curModel.pen[liabIdxLocal].dd);
          }
          if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
            fprintf (fpIR, " %9.8f", dk_curModel.pen[liabIdxLocal].DDSD);
          }
        }
        if (modelType->trait == CT) {
          fprintf (fpIR, " %9.8f", dk_curModel.pen[0].threshold);
        }
        fprintf (fpIR, " %d\n", dk_curModel.posIdx);
      }

      oldsum = alpha_integral;
      oldscale = (*scale);
      newscale = 0;

      if (log10HetLR >= DBL_MAX_10_EXP - 1) {
        /* find the new scale, adjust the current sum */
        newscale = log10HetLR - (DBL_MAX_10_EXP - SCALE_RESERVE);
      }
      if (newscale > oldscale) {
        /* need to use the newscale and adjust the old sum */
        if (oldsum > 0) {
          oldsum_log10 = log10 (oldsum);
          newsum_log10 = oldsum_log10 + oldscale - newscale;
          if (newsum_log10 <= DBL_MIN_10_EXP + 1) {
            alpha_integral = 0;
          } else {
            alpha_integral = pow (10, newsum_log10);
          }
        }
        *scale = newscale;
        oldscale = newscale;
      } else {
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
        for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
          localmax_x[k] = x[k - 1];
          localmax_x[k + 1] = (x[k] - xl[k]) * (x[k - 1] - xl[k - 1]) / (xu[k - 1] - xl[k - 1]) + xl[k];

          if (modelOptions->imprintingFlag) {
            localmax_x[k + 2] = (x[k + 1] - xl[k + 1]) * (x[k - 1] - xl[k - 1]) / (xu[k - 1] - xl[k - 1]) + xl[k + 1];
            localmax_x[k + 3] = (x[k + 2] - xl[k + 2]) * (x[k + 1] - xl[k + 1]) / (xu[k + 1] - xl[k + 1]) * (x[k] - xl[k]) / (xu[k] - xl[k]) * (x[k - 1] - xl[k - 1]) / (xu[k - 1] - xl[k - 1]) + xl[k + 2];
          } else {
            localmax_x[k + 2] = (x[k + 1] - xl[k + 1]) * (x[k] - xl[k]) / (xu[k] - xl[k]) * (x[k - 1] - xl[k - 1]) / (xu[k - 1] - xl[k - 1]) + xl[k + 1];
          }
          k += pen_size;
          if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
            /*localmax_x[k] = x[k - 1];
             * localmax_x[k + 1] = x[k];
             * localmax_x[k + 2] = x[k + 1];
             * if(modelOptions->imprintingFlag)
             * localmax_x[k + 3] = x[k + 2];
             * k += pen_size; */
            localmax_x[k] = x[k - 1];
            k++;
          }
          /* threshold for QT *
           * if (modelType->trait == CT) {
           * localmax_x[k] = x[k - 1];
           * k++;
           * } */
        }
        if (modelType->trait == CT) {
          localmax_x[k] = x[k - 1];
          k++;
        }
      }
    }

    avg_hetLR = alpha_integral;

    /* Jacobian */
    k = 1;
    for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
      avg_hetLR *= (x[k] - xl[k]) / (xu[k] - xl[k]);
      avg_hetLR *= (x[k] - xl[k]) / (xu[k] - xl[k]);
      avg_hetLR *= (x[k + 1] - xl[k + 1]) / (xu[k + 1] - xl[k + 1]);

      if (modelOptions->imprintingFlag) {
        avg_hetLR *= (x[k] - xl[k]) / (xu[k] - xl[k]);
        avg_hetLR *= (x[k + 2] - xl[k + 2]) / (xu[k + 2] - xl[k + 2]);
      }

      k += pen_size;
      if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
        k++;    //+= pen_size;
      }
      /*if (modelType->trait == CT) {
       * k++;
       * } */
    }
  }
  f[0] = avg_hetLR;

}


void compute_hlod_mp_dt (double x[], double *f, int *scale)
{

  int j;
  int pedIdx, liabIdxLocal, statusLocal, pen_size = 3, ret;

  double pen_DD, pen_Dd, pen_dD, pen_dd, gfreq, alphaV;
  double log10_likelihood_null, log10_likelihood_alternative, log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;
  double log10Likelihood;

  /* for null likelihood calculation */
  double product_likelihood = 1;        /* product of the likelihoods for all the pedigrees */
  double sum_log_likelihood = 0;        /* sum of the log10(likelihood) for all the pedigrees */

  Pedigree *pPedigreeLocal;

  int newscale, oldscale;       /* scaling related variables */
  double newLog10HetLR;
  double oldsum;
  double oldsum_log10;
  double newsum_log10;



  int origLocus = analysisLocusList->pLocusIndex[0];

  if (modelOptions->imprintingFlag)
    pen_size = 4;

  if (analysisLocusList->numLocus > 1)
    origLocus = analysisLocusList->pLocusIndex[1];

  //fprintf(stderr,"in compute hlod x %G %G %G %G\n", x[0],x[1],x[2],x[3]);
  gfreq = x[0];
  if (fpIR != NULL)
    dk_curModel.dgf = gfreq;

  for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
    pen_DD = x[pen_size * liabIdxLocal + 1];
    pen_Dd = x[pen_size * liabIdxLocal + 2] * x[pen_size * liabIdxLocal + 1];

    if (modelOptions->imprintingFlag) {
      pen_dD = x[pen_size * liabIdxLocal + 3] * x[pen_size * liabIdxLocal + 1];
      pen_dd = x[pen_size * liabIdxLocal + 4] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 2] * x[pen_size * liabIdxLocal + 3];
    } else {
      pen_dd = x[pen_size * liabIdxLocal + 3] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 2];
      pen_dD = pen_Dd;
    }
    if (fpIR != NULL) {
      dk_curModel.pen[liabIdxLocal].DD = pen_DD;
      dk_curModel.pen[liabIdxLocal].Dd = pen_Dd;
      dk_curModel.pen[liabIdxLocal].dD = pen_dD;
      dk_curModel.pen[liabIdxLocal].dd = pen_dd;
    }

    pTrait->penetrance[2][liabIdxLocal][0][0] = pen_DD;
    pTrait->penetrance[2][liabIdxLocal][0][1] = pen_Dd;
    pTrait->penetrance[2][liabIdxLocal][1][0] = pen_dD;
    pTrait->penetrance[2][liabIdxLocal][1][1] = pen_dd;
    pTrait->penetrance[1][liabIdxLocal][0][0] = 1 - pen_DD;
    pTrait->penetrance[1][liabIdxLocal][0][1] = 1 - pen_Dd;
    pTrait->penetrance[1][liabIdxLocal][1][0] = 1 - pen_dD;
    pTrait->penetrance[1][liabIdxLocal][1][1] = 1 - pen_dd;

  }

  if (modelOptions->polynomial == TRUE);
  else
    /* only need to update trait locus */
    update_penetrance (&pedigreeSet, traitLocus);

  pLocus->pAlleleFrequency[0] = gfreq;
  pLocus->pAlleleFrequency[1] = 1 - gfreq;


  if (modelOptions->polynomial != TRUE)
    update_locus (&pedigreeSet, traitLocus);

  /* for trait likelihood */
  analysisLocusList = &traitLocusList;
  xmissionMatrix = traitMatrix;

  // fprintf(stderr, "Null likelihood computation is starting for %d pedigrees \n",pedigreeSet.numPedigree);

  /* compute the null likelihood with   */
  pedigreeSet.likelihood = 1;
  pedigreeSet.log10Likelihood = 0;

#ifdef STUDYDB
  if (toupper(*studyDB.role) == 'S') {
    compute_likelihood (&pedigreeSet);
  }
#endif

  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {

    /* save the likelihood at null */
    pPedigreeLocal = pedigreeSet.ppPedigreeSet[pedIdx];


    if (modelOptions->polynomial == TRUE) {
      ASSERT (pPedigreeLocal->traitLikelihoodPolynomial != NULL, "Pedigree trait likelihood is NULL");
      /* evaluate likelihood */
      evaluatePoly (pPedigreeLocal->traitLikelihoodPolynomial, pPedigreeLocal->traitLikelihoodPolyList, &pPedigreeLocal->likelihood);
      //fprintf(stderr, " is done %f with %d pedigrees\n",pPedigreeLocal->likelihood, pedigreeSet.numPedigree);

      if (isnan (pPedigreeLocal->likelihood))
        ERROR ("Null hypothesis likelihood is not a number");

    } else {
      initialize_multi_locus_genotype (pPedigreeLocal);
      statusLocal = compute_pedigree_likelihood (pPedigreeLocal);
    }

    /*pPedigreeLocal->likelihood is now computed and now check it */
    if (pPedigreeLocal->likelihood == 0.0) {
      WARNING ("Pedigree %s has likelihood of zero or nearly zero", pPedigreeLocal->sPedigreeID);
      ret = -1;
      product_likelihood = 0.0;
      sum_log_likelihood = -9999.99;
      break;
    } else if (pPedigreeLocal->likelihood < 0.0) {
      ASSERT (pPedigreeLocal->likelihood >= 0.0, "Pedigree %s has negative likelihood", pPedigreeLocal->sPedigreeID);
      product_likelihood = 0.0;
      sum_log_likelihood = -9999.99;
      ret = -2;
      break;
    } else {
      if (pPedigreeLocal->pCount[origLocus] == 1) {
        product_likelihood *= pPedigreeLocal->likelihood;
        log10Likelihood = log10 (pPedigreeLocal->likelihood);
      } else {
        product_likelihood *= pow (pPedigreeLocal->likelihood, pPedigreeLocal->pCount[origLocus]);
        log10Likelihood = log10 (pPedigreeLocal->likelihood) * pPedigreeLocal->pCount[origLocus];
      }
      sum_log_likelihood += log10Likelihood;
    }
    pedigreeSet.nullLikelihood[pedIdx] = pPedigreeLocal->likelihood;
    //fprintf(stderr,"null likelihood pedIdx=%d is done %20.18f with product =%20.16f\n",pedIdx,pPedigreeLocal->likelihood,product_likelihood );
  }


  pedigreeSet.likelihood = product_likelihood;
  pedigreeSet.log10Likelihood = sum_log_likelihood;
  log10_likelihood_null = pedigreeSet.log10Likelihood;
  DIAG (OVERALL, 1, {
        fprintf (stderr, "Sum of log Likelihood is: %e\n", sum_log_likelihood);
      }
  );

  /* This is for alternative likelihood */
  analysisLocusList = &savedLocusList;
  xmissionMatrix = altMatrix;
  int k;
  char markerNo[8];
  sprintf (partialPolynomialFunctionName, "MDA_LC%d_C%d_P%%sM", modelRange->nlclass, (originalLocusList.ppLocusList[1])->pMapUnit->chromosome);
  for (k = 0; k < modelType->numMarkers; k++) {
    if (*get_map_position (traitLocus) <= *get_map_position (markerLocusList.pLocusIndex[k]) && (strstr (partialPolynomialFunctionName, "_T") == NULL))
      strcat (partialPolynomialFunctionName, "_T");
    sprintf (markerNo, "_%d", markerLocusList.pLocusIndex[k]);
    strcat (partialPolynomialFunctionName, markerNo);
  }
  if (strstr (partialPolynomialFunctionName, "_T") == NULL)
    strcat (partialPolynomialFunctionName, "_T");
  ret = compute_likelihood (&pedigreeSet);
  cL[4]++; // MP DT alternative likelihood
  if (ret == -2)
    ERROR ("Negative likelihood for trait");

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;
  if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
    avg_hetLR = 0.0;
  } else {
    log10_likelihood_ratio = log10_likelihood_alternative - log10_likelihood_null;
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
    //fprintf(stderr," trait %e marker %e alter%e lr %e \n", log(pedigreeSet.nullLikelihood[0]) , log(pPedigreeLocal->markerLikelihood),log(pPedigreeLocal->likelihood), pPedigreeLocal->likelihood / (pedigreeSet.nullLikelihood[0] * pPedigreeLocal->markerLikelihood));
    /* caculating the HET */
    for (j = 0; j < 5; j++) {
      //for (j = 0; j < 1; j++) {
      alphaV = alpha[j][0];
      alphaV2 = 1 - alphaV;
      if (alphaV2 < 0)
        alphaV2 = 0;

      log10HetLR = 0;
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        pPedigreeLocal = pedigreeSet.ppPedigreeSet[pedIdx];
        homoLR = pPedigreeLocal->likelihood / (pedigreeSet.nullLikelihood[pedIdx] * pPedigreeLocal->markerLikelihood);
	//	fprintf(stderr,"j=%d pedIdx=%d  %e %e %e %e \n",j, pedIdx,pPedigreeLocal->likelihood,pedigreeSet.nullLikelihood[pedIdx] , pPedigreeLocal->markerLikelihood, homoLR);
        if (alphaV * homoLR + alphaV2 < 0)
          WARNING ("Heterogenous Likelihood Ratio is less than zero");
        log10HetLR += log10 (alphaV * homoLR + alphaV2);
      }

      if (fpIR != NULL) {
        dk_curModel.alpha = alphaV;
        fprintf (fpIR, "%6.3f", log10HetLR);

        fprintf (fpIR, " %9.8f %9.8f", dk_curModel.alpha, dk_curModel.dgf);
        for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
          fprintf (fpIR, " %9.8f %9.8f", dk_curModel.pen[liabIdxLocal].DD, dk_curModel.pen[liabIdxLocal].Dd);
          if (modelOptions->imprintingFlag) {
            fprintf (fpIR, " %9.8f %9.8f", dk_curModel.pen[liabIdxLocal].dD, dk_curModel.pen[liabIdxLocal].dd);
          } else {
            fprintf (fpIR, " %9.8f", dk_curModel.pen[liabIdxLocal].dd);
          }
        }
        fprintf (fpIR, " %d\n", dk_curModel.posIdx);
      }

      oldsum = alpha_integral;
      oldscale = (*scale);
      newscale = 0;

      if (log10HetLR >= DBL_MAX_10_EXP - 1) {
        /* find the new scale, adjust the current sum */
        newscale = log10HetLR - (DBL_MAX_10_EXP - SCALE_RESERVE);
      }
      if (newscale > oldscale) {
        /* need to use the newscale and adjust the old sum */
        if (oldsum > 0) {
          oldsum_log10 = log10 (oldsum);
          newsum_log10 = oldsum_log10 + oldscale - newscale;
          if (newsum_log10 <= DBL_MIN_10_EXP + 1) {
            alpha_integral = 0;
          } else {
            alpha_integral = pow (10, newsum_log10);
          }
        }
        *scale = newscale;
        oldscale = newscale;
      } else {
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

        for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
          localmax_x[pen_size * liabIdxLocal + 2] = x[pen_size * liabIdxLocal + 1];
          localmax_x[pen_size * liabIdxLocal + 3] = x[pen_size * liabIdxLocal + 2] * x[pen_size * liabIdxLocal + 1];
          localmax_x[pen_size * liabIdxLocal + 4] = x[pen_size * liabIdxLocal + 3] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 2];

          if (modelOptions->imprintingFlag) {
            localmax_x[pen_size * liabIdxLocal + 4] = x[pen_size * liabIdxLocal + 3] * x[pen_size * liabIdxLocal + 1];
            localmax_x[pen_size * liabIdxLocal + 5] = x[pen_size * liabIdxLocal + 4] * x[pen_size * liabIdxLocal + 3] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 2];
          }
        }
      }

    }   /* end of calculating HET LR */

    avg_hetLR = alpha_integral;

    //Jacobian
    for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
      if (modelOptions->imprintingFlag)
        avg_hetLR *= x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 2] * x[pen_size * liabIdxLocal + 3];
      else
        avg_hetLR *= x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 2];
    }
  }
  *f = avg_hetLR;
}

void compute_hlod_2p_qt (double x[], double *f, int *scale)
{

  int k, j, ret;
  int pedIdx, liabIdxLocal = 0, statusLocal, pen_size = 3;
  double constraint = 0.0;
  double mean_DD = 0.0, mean_Dd = 0.0, mean_dD = 0.0, mean_dd = 0.0;
  double SD_DD = 0.0, SD_Dd = 0.0, SD_dD = 0.0, SD_dd = 0.0;
  double gfreq;
  double alphaV;
  double thetaM, thetaF;
  double threshold = 0.0;
  double log10_likelihood_null, log10_likelihood_alternative, log10_likelihood_ratio, likelihood_ratio;
  double hetLR, log10HetLR, homoLR, alphaV2;
  double alpha_integral = 0.0, avg_hetLR;

  Pedigree *pPedigreeLocal;

  int newscale, oldscale;       /* scaling related variables */
  double newLog10HetLR;
  double oldsum;
  double oldsum_log10;
  double newsum_log10;

  if (modelOptions->imprintingFlag)
    pen_size = 4;

  gfreq = x[0];
  if (fpIR != NULL)
    dk_curModel.dgf = gfreq;

  if (modelOptions->mapFlag == SS) {
    thetaM = fixed_thetaM;
    thetaF = fixed_thetaF;
  } else {
    thetaM = fixed_theta;
    thetaF = fixed_theta;
  }

  if (1 && modelOptions->markerAnalysis == FALSE) {
    pLocus->pAlleleFrequency[0] = gfreq;
    pLocus->pAlleleFrequency[1] = 1 - gfreq;

    if (modelOptions->polynomial == TRUE);
    else
      update_locus (&pedigreeSet, loc1);

  }

  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
    statusLocal = setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprimeIdx);
    if (statusLocal < 0)
      ASSERT (1, "Haplotype frequency combination impossible");
  }


  /* this should be MEAN + SD */
  j = 1;
  if (modelType->trait == CT) {
    threshold = x[s->ndim - 1];
  }
  if (modelOptions->markerAnalysis == FALSE) {
    for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
      mean_DD = x[j];
      mean_Dd = (x[j + 1] - xl[j + 1]) * (x[j] - xl[j]) / (xu[j] - xl[j]) + xl[j + 1];

      if (modelOptions->imprintingFlag) {
        mean_dD = (x[j + 2] - xl[j + 2]) * (x[j] - xl[j]) / (xu[j] - xl[j]) + xl[j + 2];
        mean_dd = (x[j + 3] - xl[j + 3]) * (x[j + 2] - xl[j + 2]) / (xu[j + 2] - xl[j + 2]) * (x[j + 1] - xl[j + 1]) / (xu[j + 1] - xl[j + 1]) * (x[j] - xl[j]) / (xu[j] - xl[j]) + xl[j + 3];
      } else {
        mean_dd = (x[j + 2] - xl[j + 2]) * (x[j + 1] - xl[j + 1]) / (xu[j + 1] - xl[j + 1]) * (x[j] - xl[j]) / (xu[j] - xl[j]) + xl[j + 2];
        mean_dD = mean_Dd;
      }
      j += pen_size;

      if (fpIR != NULL) {
        dk_curModel.pen[liabIdxLocal].DD = mean_DD;
        dk_curModel.pen[liabIdxLocal].Dd = mean_Dd;
        dk_curModel.pen[liabIdxLocal].dD = mean_dD;
        dk_curModel.pen[liabIdxLocal].dd = mean_dd;
      }

      if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
        /*SD_DD = x[j];         
         * SD_Dd = x[j + 1];    
         * if(modelOptions->imprintingFlag){
         * SD_dD = x[j+2];
         * SD_dd = x[j+3];
         * }else{
         * SD_dd = x[j+2];      
         * SD_dD= SD_Dd;
         * }
         * j += pen_size; */
        SD_DD = SD_Dd = SD_dD = SD_dd = x[j++];

        if (fpIR != NULL) {
          dk_curModel.pen[liabIdxLocal].DDSD = SD_DD;
          dk_curModel.pen[liabIdxLocal].DdSD = SD_Dd;
          dk_curModel.pen[liabIdxLocal].dDSD = SD_dD;
          dk_curModel.pen[liabIdxLocal].ddSD = SD_dd;
        }

      }
      if (fpIR != NULL) {
        if (modelType->trait == CT)
          dk_curModel.pen[liabIdxLocal].threshold = threshold;
      }

      /* threshold for QT *
       * if (modelType->trait == CT) {
       * threshold = x[j];      // modelRange->tthresh[liabIdxLocal][thresholdIdx];
       * j++;
       * } */


      /* check against the hard coded constraint */
      if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
        constraint = (1.0 - gfreq) * (1.0 - gfreq) * mean_dd * SD_dd + 2 * gfreq * (1 - gfreq) * mean_Dd * SD_Dd + gfreq * gfreq * mean_DD * SD_DD;
        /* fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
         * constraint, gfreq, mean_DD, SD_DD, mean_Dd, SD_DD, mean_dd, SD_dd);
         */
        if (constraint >= 3.0 || constraint <= -3.0) {
          // fprintf(stderr,"Constraint is %f \n", constraint);
          num_out_constraint++;
        }
      }

      pTrait->means[liabIdxLocal][0][0] = mean_DD;
      pTrait->means[liabIdxLocal][0][1] = mean_Dd;
      pTrait->means[liabIdxLocal][1][0] = mean_dD;
      pTrait->means[liabIdxLocal][1][1] = mean_dd;
      pTrait->stddev[liabIdxLocal][0][0] = SD_DD;
      pTrait->stddev[liabIdxLocal][0][1] = SD_Dd;
      pTrait->stddev[liabIdxLocal][1][0] = SD_dD;
      pTrait->stddev[liabIdxLocal][1][1] = SD_dd;

      /* threshold for QT */
      pTrait->cutoffValue[liabIdxLocal] = threshold;

    }   /* liability class Index */

    if (modelOptions->polynomial == TRUE);
    else
      update_penetrance (&pedigreeSet, traitLocus);

  }



  /* marker to marker analysis */

  /* get the likelihood at 0.5 first and LD=0 */
  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {

    statusLocal = setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprime0Idx);
    if (statusLocal < 0)
      ASSERT (1, "Haplotype frequency combination impossible. Exiting!\n");

    set_null_dprime (pLDLoci);
    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprime0Idx]);
    copy_DValue (pLDLoci, pLambdaCell->DValue[dprime0Idx]);

    ASSERT (pLambdaCell->impossibleFlag[dprime0Idx] == 0, "Haplotype frequency combination impossible at LE. Exiting!\n");
  }

  for (k = 0; k < 3; k++) {
    analysisLocusList->pNextLocusDistance[k][0] = 0.5;
    analysisLocusList->pPrevLocusDistance[k][1] = 0.5;
  }

  if (modelOptions->polynomial == TRUE);
  else
    /* populate the matrix */
    statusLocal = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, /* probability */
        initialProbAddr2,       /* probability */
        initialHetProbAddr, 0,  /* cell index */
        -1, -1, /* last het locus & last het pattern (P-1 or M-2) */
        0);     /* current locus - start with 0 */


  sprintf (partialPolynomialFunctionName, "TQ_LC%d_C%d_P%%s_%s_%s", modelRange->nlclass, pLocus2->pMapUnit->chromosome, pLocus1->sName, pLocus2->sName);
  ret = compute_likelihood (&pedigreeSet);
  cL[5]++; // TP QT

  if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
    if (modelRange->atypicalQtTrait)
      WARNING ("Pedigree set has likelihood of zero or nearly zero");
    else
      ERROR ("Pedigree set has likelihood of zero or nearly zero");
    
    f[0] = 1.0;
    return;
  }

  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
    /* save the likelihood at null */
    pPedigreeLocal = pedigreeSet.ppPedigreeSet[pedIdx];
    pedigreeSet.nullLikelihood[pedIdx] = pPedigreeLocal->likelihood;
  }

  log10_likelihood_null = pedigreeSet.log10Likelihood;

  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
    copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);
    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
    copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);
  }

  if (modelOptions->mapFlag == SA) {
    for (k = 0; k < 3; k++) {
      analysisLocusList->pNextLocusDistance[k][0] = thetaM;
      analysisLocusList->pPrevLocusDistance[k][1] = thetaF;
    }
  } else {
    analysisLocusList->pNextLocusDistance[MAP_POS_MALE][0] = analysisLocusList->pPrevLocusDistance[MAP_POS_MALE][1] = thetaM;
    analysisLocusList->pNextLocusDistance[MAP_POS_FEMALE][0] = analysisLocusList->pPrevLocusDistance[MAP_POS_FEMALE][1] = thetaF;
  }

  if (modelOptions->polynomial == TRUE);
  else
    /* populate the matrix */
    statusLocal = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, /* probability */
        initialProbAddr2,       /* probability */
        initialHetProbAddr, 0,  /* cell index */
        -1, -1, /* last het locus & last het pattern (P-1 or M-2) */
        0);     /* current locus - start with 0 */

  // No new name for a polynomial here because we're reusing the existing one
  ret = compute_likelihood (&pedigreeSet);
  cL[6]++; // TP QT
  if (ret == -2)
    ERROR ("Negative likelihood for markers");

  log10_likelihood_alternative = pedigreeSet.log10Likelihood;
  if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
    avg_hetLR = 0.0;

    fprintf (stderr, "Theta %f has likelihood 0 \n", thetaM);
    fprintf (stderr, "dgf=%f\n", gfreq);
    for (j = 1; j < s->ndim; j++) {
      fprintf (stderr, " %f", x[j]);
    }
    fprintf (stderr, "mean %f %f %f %f SD %f %f %f %f\n", mean_DD, mean_Dd, mean_dD, mean_dd, SD_DD, SD_Dd, SD_dD, SD_dd);
    fprintf (stderr, "\n");

  } else {
    log10_likelihood_ratio = log10_likelihood_alternative - log10_likelihood_null;
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
        pPedigreeLocal = pedigreeSet.ppPedigreeSet[pedIdx];
        homoLR = pPedigreeLocal->likelihood / pedigreeSet.nullLikelihood[pedIdx];
        log10HetLR += log10 (alphaV * homoLR + alphaV2);

        if (isnan (homoLR)) {
          printf ("pedIdx =%d  homeLR=%e log10HLR=%e\n", pedIdx, homoLR, log10HetLR);
          exit (0);
        }
      }

      if (fpIR != NULL) {
        dk_curModel.alpha = alphaV;
        fprintf (fpIR, "%6.3f", log10HetLR);
        if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
          fprintf (fpIR, " %9.8f", dk_curModel.dprime[0]);
        }
        if (modelOptions->mapFlag == SA) {
          fprintf (fpIR, " %9.8f", dk_curModel.theta[0]);
        } else {
          fprintf (fpIR, " %9.8f %9.8f", dk_curModel.theta[0], dk_curModel.theta[1]);
        }
        fprintf (fpIR, " %9.8f %9.8f", dk_curModel.alpha, dk_curModel.dgf);
        for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
          fprintf (fpIR, " %9.8f %9.8f", dk_curModel.pen[liabIdxLocal].DD, dk_curModel.pen[liabIdxLocal].Dd);
          if (modelOptions->imprintingFlag) {
            fprintf (fpIR, " %9.8f %9.8f", dk_curModel.pen[liabIdxLocal].dD, dk_curModel.pen[liabIdxLocal].dd);
          } else {
            fprintf (fpIR, " %9.8f", dk_curModel.pen[liabIdxLocal].dd);
          }
          if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
            fprintf (fpIR, " %9.8f", dk_curModel.pen[liabIdxLocal].DDSD);
          }
        }
        if (modelType->trait == CT) {
          fprintf (fpIR, " %9.8f", dk_curModel.pen[0].threshold);
        }
        fprintf (fpIR, " %d\n", dk_curModel.posIdx);
      }

      oldsum = alpha_integral;
      oldscale = (*scale);
      newscale = 0;

      if (log10HetLR >= DBL_MAX_10_EXP - 1) {
        /* find the new scale, adjust the current sum */
        newscale = log10HetLR - (DBL_MAX_10_EXP - SCALE_RESERVE);
      }
      if (newscale > oldscale) {
        /* need to use the newscale and adjust the old sum */
        if (oldsum > 0) {
          oldsum_log10 = log10 (oldsum);
          newsum_log10 = oldsum_log10 + oldscale - newscale;
          if (newsum_log10 <= DBL_MIN_10_EXP + 1) {
            alpha_integral = 0;
          } else {
            alpha_integral = pow (10, newsum_log10);
          }
        }
        *scale = newscale;
        oldscale = newscale;
      } else {
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
        for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {

          localmax_x[k] = x[k - 1];
          localmax_x[k + 1] = (x[k] - xl[k]) * (x[k - 1] - xl[k - 1]) / (xu[k - 1] - xl[k - 1]) + xl[k];

          if (modelOptions->imprintingFlag) {
            localmax_x[k + 2] = (x[k + 1] - xl[k + 1]) * (x[k - 1] - xl[k - 1]) / (xu[k - 1] - xl[k - 1]) + xl[k + 1];
            localmax_x[k + 3] = (x[k + 2] - xl[k + 2]) * (x[k + 1] - xl[k + 1]) / (xu[k + 1] - xl[k + 1]) * (x[k] - xl[k]) / (xu[k] - xl[k]) * (x[k - 1] - xl[k - 1]) / (xu[k - 1] - xl[k - 1]) + xl[k + 2];
          } else {
            localmax_x[k + 2] = (x[k + 1] - xl[k + 1]) * (x[k] - xl[k]) / (xu[k] - xl[k]) * (x[k - 1] - xl[k - 1]) / (xu[k - 1] - xl[k - 1]) + xl[k + 1];
          }
          k += pen_size;
          if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
            /*localmax_x[k] = x[k - 1];
             * localmax_x[k + 1] = x[k];
             * localmax_x[k + 2] = x[k + 1];
             * if(modelOptions->imprintingFlag)
             * localmax_x[k + 3] = x[k + 2];
             * k += pen_size; */
            localmax_x[k] = x[k - 1];
            k++;
          }
          /* threshold for QT *
           * if (modelType->trait == CT) {
           * localmax_x[k] = x[k - 1];
           * k++;
           * } */
        }
        if (modelType->trait == CT) {
          localmax_x[k] = x[k - 1];
          k++;
        }
      }
    }

    avg_hetLR = alpha_integral;

    /* Jacobian */
    k = 1;
    for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
      avg_hetLR *= (x[k] - xl[k]) / (xu[k] - xl[k]);
      avg_hetLR *= (x[k] - xl[k]) / (xu[k] - xl[k]);
      avg_hetLR *= (x[k + 1] - xl[k + 1]) / (xu[k + 1] - xl[k + 1]);

      if (modelOptions->imprintingFlag) {
        avg_hetLR *= (x[k] - xl[k]) / (xu[k] - xl[k]);
        avg_hetLR *= (x[k + 2] - xl[k + 2]) / (xu[k + 2] - xl[k + 2]);
      }

      k += pen_size;
      if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
        k++;    //= pen_size;
      }
      /*if (modelType->trait == CT) {
       * k++;
       * } */
    }
    /*This is a temporary checking. *
     * if (k != s->ndim) {
     * printf ("k=%d  while dim for BR is %d\n", k, s->ndim);
     * exit (EXIT_FAILURE);
     * } */

  }
  *f = avg_hetLR;
}

void compute_hlod_2p_dt (double x[], double *f, int *scale)
{

//double compute_hlod(PedigreeSet *pedigreeSet,double x[], int loc1, int loc2, Locus *pLocus, Trait *pTrait, int traitLocus, int totalLoci, double * initialProbAddr[3], Locus *pLocus1){

/*  Limit of this function

   modelOptions->type := TP  (Two points)
   modelOptions->trait := DT (DICHOTOMOUS);
   modelOptions->equilibrium :=LINKAGE_EQUILIBRIUM
   modelOptions->polynomial := TRUE
   modelOptions->markerAnalysis := FALSE;
   modelOptions->mapFlag := SA 
   modelRange->nafreq :=1
   
*/

  int k, j;
  int pedIdx, liabIdxLocal = 0, statusLocal, pen_size = 3, ret;

  double pen_DD, pen_Dd, pen_dD, pen_dd, gfreq, alphaV, thetaM, thetaF;
  double log10_likelihood_null, log10_likelihood_alternative, log10_likelihood_ratio, likelihood_ratio;
  double log10HetLR, homoLR, alphaV2;
  double hetLR, avg_hetLR, alpha_integral = 0.0, tmp;
  Pedigree *pPedigreeLocal;

  int newscale, oldscale;       /* scaling related variables */
  double newLog10HetLR;
  double oldsum;
  double oldsum_log10;
  double newsum_log10;


  if (modelOptions->imprintingFlag)
    pen_size = 4;

  gfreq = x[0];
  // printf("Calculating hetLR with gf=%f  theta=%f\n", gfreq,fixed_theta);
  if (fpIR != NULL)
    dk_curModel.dgf = gfreq;

  if (modelOptions->mapFlag == SS) {
    thetaM = fixed_thetaM;
    thetaF = fixed_thetaF;
  } else {
    thetaM = fixed_theta;
    thetaF = fixed_theta;
  }



  if (1 && modelOptions->markerAnalysis == FALSE) {
    pLocus->pAlleleFrequency[0] = gfreq;
    pLocus->pAlleleFrequency[1] = 1 - gfreq;

    if (modelOptions->polynomial == TRUE);
    else
      update_locus (&pedigreeSet, loc1);

  }

  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
    statusLocal = setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprimeIdx);
    if (statusLocal < 0)
      ASSERT (1, "Haplotype frequency combination impossible. Exiting!\n");
  }

  if (modelOptions->markerAnalysis == FALSE && pLocus1->locusType == LOCUS_TYPE_TRAIT) {
    for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
      pen_DD = x[pen_size * liabIdxLocal + 1];
      pen_Dd = x[pen_size * liabIdxLocal + 2] * x[pen_size * liabIdxLocal + 1];

      if (modelOptions->imprintingFlag) {
        pen_dD = x[pen_size * liabIdxLocal + 3] * x[pen_size * liabIdxLocal + 1];
        pen_dd = x[pen_size * liabIdxLocal + 4] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 2] * x[pen_size * liabIdxLocal + 3];
      } else {
        pen_dd = x[pen_size * liabIdxLocal + 3] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 2];
        pen_dD = pen_Dd;
      }
      if (fpIR != NULL) {
        dk_curModel.pen[liabIdxLocal].DD = pen_DD;
        dk_curModel.pen[liabIdxLocal].Dd = pen_Dd;
        dk_curModel.pen[liabIdxLocal].dD = pen_dD;
        dk_curModel.pen[liabIdxLocal].dd = pen_dd;
      }

      pTrait->penetrance[2][liabIdxLocal][0][0] = pen_DD;
      pTrait->penetrance[2][liabIdxLocal][0][1] = pen_Dd;
      pTrait->penetrance[2][liabIdxLocal][1][0] = pen_dD;
      pTrait->penetrance[2][liabIdxLocal][1][1] = pen_dd;
      pTrait->penetrance[1][liabIdxLocal][0][0] = 1 - pen_DD;
      pTrait->penetrance[1][liabIdxLocal][0][1] = 1 - pen_Dd;
      pTrait->penetrance[1][liabIdxLocal][1][0] = 1 - pen_dD;
      pTrait->penetrance[1][liabIdxLocal][1][1] = 1 - pen_dd;
    }
  }
  if (modelOptions->polynomial == TRUE);
  else
    update_penetrance (&pedigreeSet, traitLocus);

  /* get the likelihood at 0.5 first and LD=0 */
  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
    statusLocal = setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprime0Idx);
    if (statusLocal < 0)
      ASSERT (1, "Haplotype frequency combination impossible. Exiting!\n");

    set_null_dprime (pLDLoci);
    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprime0Idx]);
    copy_DValue (pLDLoci, pLambdaCell->DValue[dprime0Idx]);
    ASSERT (pLambdaCell->impossibleFlag[dprime0Idx] == 0, "Haplotype frequency combination impossible at LE. Exiting!\n");
  }
  for (k = 0; k < 3; k++) {
    analysisLocusList->pNextLocusDistance[k][0] = 0.5;
    analysisLocusList->pPrevLocusDistance[k][1] = 0.5;
  }
  if (modelOptions->polynomial == TRUE);
  else
    /* populate the matrix */
    statusLocal = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, /* probability */
        initialProbAddr2,       /* probability */
        initialHetProbAddr, 0,  /* cell index */
        -1, -1, /* last het locus & last het pattern (P-1 or M-2) */
        0);     /* current locus - start with 0 */


  sprintf (partialPolynomialFunctionName, "TD_LC%d_C%d_P%%s_%s_%s", modelRange->nlclass, pLocus2->pMapUnit->chromosome, pLocus1->sName, pLocus2->sName);
  cL[7]++; // TP DT
  ret = compute_likelihood (&pedigreeSet);

  //printf("likelihood =%15.13f with theta 0.5 with %d pedigrees\n", pedigreeSet.likelihood, pedigreeSet.numPedigree);

  if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
    fprintf (stderr, "Theta 0.5 has likelihood 0\n");
    fprintf (stderr, "dgf=%f\n", gfreq);
    //exit (EXIT_FAILURE);
    f[0] = 1.0;
    return;
  }

  /* save the results for NULL */
  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
    /* save the likelihood at null */
    pPedigreeLocal = pedigreeSet.ppPedigreeSet[pedIdx];
    pedigreeSet.nullLikelihood[pedIdx] = pPedigreeLocal->likelihood;
  }

  log10_likelihood_null = pedigreeSet.log10Likelihood;


  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
    copy_dprime (pLDLoci, pLambdaCell->lambda[dprimeIdx]);

    copy_haploFreq (pLDLoci, pLambdaCell->haploFreq[dprimeIdx]);
    copy_DValue (pLDLoci, pLambdaCell->DValue[dprimeIdx]);

    /* calculate R square if the marker is a SNP */
    if (R_square_flag == TRUE)
      R_square = calculate_R_square (pLocus1->pAlleleFrequency[0], pLocus2->pAlleleFrequency[0], pLDLoci->ppDValue[0][0]);
    else
      R_square = -1;
  }


  if (modelOptions->mapFlag == SA) {
    for (k = 0; k < 3; k++) {
      analysisLocusList->pNextLocusDistance[k][0] = thetaM;
      analysisLocusList->pPrevLocusDistance[k][1] = thetaF;
    }
  } else {
    analysisLocusList->pNextLocusDistance[MAP_POS_MALE][0] = analysisLocusList->pPrevLocusDistance[MAP_POS_MALE][1] = thetaM;
    analysisLocusList->pNextLocusDistance[MAP_POS_FEMALE][0] = analysisLocusList->pPrevLocusDistance[MAP_POS_FEMALE][1] = thetaF;
  }

  if (modelOptions->polynomial == TRUE);
  else
    /* populate the matrix */
    statusLocal = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr, /* probability */
        initialProbAddr2,       /* probability */
        initialHetProbAddr, 0,  /* cell index */
        -1, -1, /* last het locus & last het pattern (P-1 or M-2) */
        0);     /* current locus - start with 0 */

  // No new name for a polynomial here because we're reusing the existing one
  cL[8]++; // TP DT
  ret = compute_likelihood (&pedigreeSet);
  if (ret == -2)
    ERROR ("Negative alternative likelihood");
  log10_likelihood_alternative = pedigreeSet.log10Likelihood;

  if (pedigreeSet.likelihood == 0.0 && pedigreeSet.log10Likelihood == -9999.99) {
    log10_likelihood_ratio = 0;
    avg_hetLR = 0.0;
  } else {
    log10_likelihood_ratio = log10_likelihood_alternative - log10_likelihood_null;

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
        pPedigreeLocal = pedigreeSet.ppPedigreeSet[pedIdx];
        homoLR = pPedigreeLocal->likelihood / pedigreeSet.nullLikelihood[pedIdx];
        tmp = log10 (alphaV * homoLR + (1 - alphaV));
        log10HetLR += tmp * pPedigreeLocal->pCount[loc2];
        // fprintf(stderr, "tmp=%15.10f cout=%d   log10= 15.10f \n",tmp, pPedigreeLocal->pCount[loc2], log10HetLR);
        // log10HetLR += tmp * sqrt(pPedigreeLocal->pCount[loc2] );  //kelvin log10 exponential for case-control
      }

      if (fpIR != NULL) {
        dk_curModel.alpha = alphaV;
        fprintf (fpIR, "%6.3f", log10HetLR);
        if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
          fprintf (fpIR, " %9.8f", dk_curModel.dprime[0]);
        }
        if (modelOptions->mapFlag == SA) {
          fprintf (fpIR, " %9.8f", dk_curModel.theta[0]);
        } else {
          fprintf (fpIR, " %9.8f %9.8f", dk_curModel.theta[0], dk_curModel.theta[1]);
        }
        fprintf (fpIR, " %9.8f %9.8f", dk_curModel.alpha, dk_curModel.dgf);
        for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
          fprintf (fpIR, " %9.8f %9.8f", dk_curModel.pen[liabIdxLocal].DD, dk_curModel.pen[liabIdxLocal].Dd);
          if (modelOptions->imprintingFlag) {
            fprintf (fpIR, " %9.8f %9.8f", dk_curModel.pen[liabIdxLocal].dD, dk_curModel.pen[liabIdxLocal].dd);
          } else {
            fprintf (fpIR, " %9.8f", dk_curModel.pen[liabIdxLocal].dd);
          }
        }
        fprintf (fpIR, " %d\n", dk_curModel.posIdx);
      }


      oldsum = alpha_integral;
      oldscale = (*scale);
      newscale = 0;

      if (log10HetLR >= DBL_MAX_10_EXP - 1) {
        /* find the new scale, adjust the current sum */
        newscale = log10HetLR - (DBL_MAX_10_EXP - SCALE_RESERVE);
      }
      if (newscale > oldscale) {
        /* need to use the newscale and adjust the old sum */
        if (oldsum > 0) {
          oldsum_log10 = log10 (oldsum);
          newsum_log10 = oldsum_log10 + oldscale - newscale;
          if (newsum_log10 <= DBL_MIN_10_EXP + 1) {
            alpha_integral = 0;
          } else {
            alpha_integral = pow (10, newsum_log10);
          }
        }
        *scale = newscale;
        oldscale = newscale;
      } else {
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

        for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
          localmax_x[pen_size * liabIdxLocal + 2] = x[pen_size * liabIdxLocal + 1];
          localmax_x[pen_size * liabIdxLocal + 3] = x[pen_size * liabIdxLocal + 2] * x[pen_size * liabIdxLocal + 1];
          localmax_x[pen_size * liabIdxLocal + 4] = x[pen_size * liabIdxLocal + 3] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 2];

          if (modelOptions->imprintingFlag) {
            localmax_x[pen_size * liabIdxLocal + 4] = x[pen_size * liabIdxLocal + 3] * x[pen_size * liabIdxLocal + 1];
            localmax_x[pen_size * liabIdxLocal + 5] = x[pen_size * liabIdxLocal + 4] * x[pen_size * liabIdxLocal + 3] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 2];
          }
        }
      }
      // fprintf(fphlod,"%f %f %f %f %f %f %f\n", log10(hetLR*x[1]*x[1]*x[2]), gfreq, pen_DD,pen_Dd, pen_dd, alphaV,fixed_theta);
    }   //end of calculating the HET         


    avg_hetLR = alpha_integral;

    //printf("avg hetLR =%15.10f with gf=%f DD=%f Dd=%f dd=%f theta=%f\n", avg_hetLR, gfreq, pen_DD,pen_Dd, pen_dd, fixed_theta);

    // Jacobian
    for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
      if (modelOptions->imprintingFlag)
        avg_hetLR *= x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 2] * x[pen_size * liabIdxLocal + 3];
      else
        avg_hetLR *= x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 1] * x[pen_size * liabIdxLocal + 2];
    }

  }

  f[0] = avg_hetLR;

}



/**

  Driver for dynamic integration analysis.

  Relatively recent cloning of the edifice of kelvin main program. Needs refactoring,
  and some progress has been made in that regard.

  @author Sang-Cheol Seok - overall content adapted from Yungui's iterative version..
  @author Bill Valentine-Cooper - progress tracking.

*/
void integrateMain ()
{

  int numPositions;
  int size_BR;
  int i, j, k;
  int liabIdxLocal, pedIdx, statusLocal;
  Pedigree *pPedigreeLocal;

  SUBSTEP (0, "Setting-up for integration-based analysis");

  /* total_dim is the number of all parameters in the 3-layer scheme
   * s->dim in dcuhre.c is the number of parameters in the middle layer alone */

  DETAIL (0, "Calculating dimensionality of outer and inner layers");
  total_dim = 2;        // alpha gf
  total_dim += 3 * modelRange->nlclass; //DD Dd dd
  if (modelOptions->imprintingFlag)
    total_dim += modelRange->nlclass;   //dD

  if (modelType->trait != DT) {
    if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
      total_dim += modelRange->nlclass;
    }
    if (modelType->trait == CT) {
      total_dim++;      //  One threshold for all LCs    //   = modelRange->nlclass;
    }
  }

  size_BR = total_dim;
  if (modelType->type == TP) {
    total_dim += 1;     // theta;
    if (modelOptions->mapFlag == SS)
      total_dim += 1;   // theta sex-specific case;

    if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
      total_dim += 1;   // dprime
  }

  DETAIL (0, "Outer dimension is %d, inner (BR) dimension is %d", total_dim, size_BR);
  DETAIL (0, "Allocating and initializing storage for analysis");

  MALCHOKE (xl, size_BR * sizeof (double), double *);
  MALCHOKE (xu, size_BR * sizeof (double), double *);
  for (i = 0; i < size_BR; i++) {
    xl[i] = 0;
    xu[i] = 1;
  }

  memset (&dk_globalmax, 0, sizeof (st_DKMaxModel));
  memset (&dk_dprime0max, 0, sizeof (st_DKMaxModel));
  memset (&dk_theta0max, 0, sizeof (st_DKMaxModel));
  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
    /* Assumes that dkelvin can only handle a single D' */
    CALCHOKE (dk_globalmax.dprime, (size_t) 1, sizeof (double), double *);
    CALCHOKE (dk_dprime0max.dprime, (size_t) 1, sizeof (double), double *);
    CALCHOKE (dk_theta0max.dprime, (size_t) 1, sizeof (double), double *);
  }
  CALCHOKE (dk_globalmax.pen, (size_t) modelRange->nlclass, sizeof (st_DKMaxModelPenVector), void *);
  CALCHOKE (dk_dprime0max.pen, (size_t) modelRange->nlclass, sizeof (st_DKMaxModelPenVector), void *);
  CALCHOKE (dk_theta0max.pen, (size_t) modelRange->nlclass, sizeof (st_DKMaxModelPenVector), void *);

  if (fpIR != NULL) {
    memset (&dk_curModel, 0, sizeof (st_DKMaxModel));
    if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
      /* Assumes that dkelvin can only handle a single D' */
      CALCHOKE (dk_curModel.dprime, (size_t) 1, sizeof (double), double *);
    }
    CALCHOKE (dk_curModel.pen, (size_t) modelRange->nlclass, sizeof (st_DKMaxModelPenVector), void *);
    writeSurfaceFileHeader ();
  }

  if (modelType->trait != DT) {
    /* Setting ranges for each variables. Default is [0,1] */
    k = 1;
    for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
      if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
        xl[k] = xl[k + 1] = xl[k + 2] = -3;
        xu[k] = xu[k + 1] = xu[k + 2] = 3;
        if (modelOptions->imprintingFlag) {
          xl[k + 3] = -3;
          xu[k + 3] = 3;
        }
      } else {
        xl[k] = modelRange->penetLimits[0][0];  //0.1;
        xl[k + 1] = modelRange->penetLimits[1][0];      //0.1;
        xl[k + 2] = modelRange->penetLimits[3][0];      //0.1;
        xu[k] = modelRange->penetLimits[0][1];  //30.0;
        xu[k + 1] = modelRange->penetLimits[1][1];      //30.0;
        xu[k + 2] = modelRange->penetLimits[3][1];      //30.0; //10.0; //23.0; //4.0;//30;

        if (modelOptions->imprintingFlag) {
          xl[k + 2] = modelRange->penetLimits[2][0];
          xu[k + 2] = modelRange->penetLimits[2][1];
          xl[k + 3] = modelRange->penetLimits[3][0];
          xu[k + 3] = modelRange->penetLimits[3][1];
        }
        //fprintf(stderr,"modelranage= %f %f %f %f %f %f %f %f\n",modelRange->penetLimits[0][0],modelRange->penetLimits[1][0],modelRange->penetLimits[2][0],modelRange->penetLimits[3][0],modelRange->penetLimits[0][1],modelRange->penetLimits[1][1],modelRange->penetLimits[2][1],modelRange->penetLimits[3][1]);
      }
      volume_region *= (xu[k] - xl[k]);
      volume_region *= (xu[k + 1] - xl[k + 1]);
      volume_region *= (xu[k + 2] - xl[k + 2]);
      if (modelOptions->imprintingFlag) {
        volume_region *= (xu[k + 3] - xl[k + 3]);
        k++;
      }
      k += 3;

      if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
        /*xl[k] = xl[k + 1] = xl[k + 2] = 0.7;
         * xu[k] = xu[k + 1] = xu[k + 2] = 1.0;//3.0;
         * if(modelOptions->imprintingFlag){
         * xl[k+3]= 0.7;
         * xu[k+3]= 1.0;
         * }
         * volume_region *= (xu[k] - xl[k]);
         * volume_region *= (xu[k + 1] - xl[k + 1]);
         * volume_region *= (xu[k + 2] - xl[k + 2]);
         * if(modelOptions->imprintingFlag){
         * volume_region *= (xu[k + 3] - xl[k + 3]);
         * k++;
         * }
         * k += 3;
         */
        xl[k] = 0.7;
        xu[k] = 1.0;
        volume_region *= (xu[k] - xl[k]);
        k++;
      }
      /*if (modelType->trait == CT) {
       * xl[k] = modelRange->tthresh[liabIdxLocal][0];//0.3;
       * xu[k] = modelRange->tthresh[liabIdxLocal][modelRange->ntthresh -1];// 23.0;
       * volume_region *= (xu[k] - xl[k]);
       * k++;
       * //   fprintf(stderr, " in CT\n ");
       * 
       * } */
    }   // retangular volume region is calculated and stored in volume_region
    if (modelType->trait == CT) {
      xl[k] = 0.0;      // modelRange->tthresh[liabIdxLocal][0];//0.3;
      xu[k] = 3.0;
      volume_region *= (xu[k] - xl[k]);
      k++;
    }
    // fprintf (stderr,"The number of dimension for calculation of BR is %d\n",k);
  }

  /*fpDK header */
  if (fpDK != NULL) {
    if (modelType->type == TP) {
      fprintf (fpDK, "num D1 Theta(M,F) numLR BR error scale MOD\n");
    } else {
      fprintf (fpDK, "traitPos ppl BR error numLR scale MOD\n");
    }
    fflush (fpDK);
  }

  /* only for multipoint - we don't handle LD under multipoint yet */
  /* DCUHRE do now use likelihoodDT or likelihoodQT to store null likelihoods */

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

  /* assume the trait locus is the first one in the list */
  traitLocus = 0;
  pLocus = originalLocusList.ppLocusList[traitLocus];
  TraitLocus *pTraitLocus = originalLocusList.ppLocusList[traitLocus]->pTraitLocus;
  pTrait = pTraitLocus->pTraits[traitLocus];
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

      /* initialize it */
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

    SUBSTEP (0, "Building transmission matrix");

    if (modelOptions->polynomial == TRUE) {
      /* populate the matrix */
      statusLocal = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,       /* probability */
          initialProbAddr2,     /* probability */
          initialHetProbAddr, 0,        /* cell index */
          -1,   /* last het locus */
          -1,   /* last  pattern (P-1 or M-2) */
          0);   /* current locus - start with 0 */

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

      if (modelOptions->markerAnalysis != FALSE && pLocus1->locusType != LOCUS_TYPE_MARKER) {
        DETAIL (0, "Skipping combination involving trait locus for marker analysis");
        continue;       // If we're doing a marker analysis, don't let the first locus be a trait
      }
      if ((pLocus1->numAllele <= 1) || ((pLocus1->numAllele == 2) && ((pLocus1->pAlleleFrequency[0] <= ERROR_MARGIN) || (pLocus1->pAlleleFrequency[1] <= ERROR_MARGIN)))) {
        WARNING ("Biallelic marker %s has a minor allele frequency less than %g, skipping!", pLocus1->sName, ERROR_MARGIN);
        continue;       // Skip MAF0 markers
      }

      for (loc2 = loc1 + 1; loc2 < originalLocusList.numLocus; loc2++) {
        if (fpIR != NULL) {
          dk_curModel.posIdx = loc2;
        }
        overallMOD = DBL_MIN_10_EXP + 1;        //0.0;  // global max 
        overallMin = DBL_MAX;
        dprime0_MOD = 0.0;      // max when D' == 0
        theta0_MOD = 0.0;       // max when Theta == 0
        /* Since dynamic sampling is unlikely to ever sample at Theta == 0,
         * we'll need to narrow down the Theta that's closest. Start by setting
         * the thetas for theta0max to a "large" value.
         */
        dk_theta0max.theta[0] = dk_theta0max.theta[1] = 0.5;
        /* Same thing for dprime0max as for theta0max, above */
        if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
          dk_dprime0max.dprime[0] = 1;

        pLocus2 = originalLocusList.ppLocusList[loc2];
        if (pLocus2->locusType != LOCUS_TYPE_MARKER)
          continue;
        if ((pLocus2->numAllele <= 1) || ((pLocus2->numAllele == 2) && ((pLocus2->pAlleleFrequency[0] <= ERROR_MARGIN) || (pLocus2->pAlleleFrequency[1] <= ERROR_MARGIN)))) {
          WARNING ("Biallelic marker %s has a minor allele frequency less than %g, skipping!", pLocus2->sName, ERROR_MARGIN);
          continue;
        }
        savedLocusList.pLocusIndex[1] = loc2;

        if (modelOptions->markerAnalysis == MM)
          SUBSTEP ((loc2 - 1) * 100 / (originalLocusList.numLocus - 1), "Starting w/loci %s(%d alleles) and %s(%d alleles)", pLocus1->sName, pLocus1->numOriginalAllele, pLocus2->sName, pLocus2->numOriginalAllele);
        else
          SUBSTEP ((loc2 - 1) * 100 / (originalLocusList.numLocus - 1),
              "Starting w/loci %s(%d alleles) and %s(%d alleles) (%d of %d pairs)", pLocus1->sName, pLocus1->numOriginalAllele, pLocus2->sName, pLocus2->numOriginalAllele, loc2, originalLocusList.numLocus - 1);

        /* Find out number of alleles this marker locus has *//* Check if this is okay with DCUHRE  ???????????? */
        if (modelOptions->equilibrium == LINKAGE_DISEQUILIBRIUM) {

          if (pLocus1->numOriginalAllele + pLocus2->numOriginalAllele > 4)
            ERROR ("Integration-based LD analysis not available for polyallelic loci");

          /* get the LD parameters */
          pLambdaCell = findLambdas (modelRange, pLocus1->numOriginalAllele, pLocus2->numOriginalAllele);
          reallocate_LD_loci (pLDLoci, pLocus1->numOriginalAllele, pLocus2->numOriginalAllele);

          /* Create these variables ahead of likelihood polynomial build to prevent
           * in-build creation. This allows polynomial compilation as the entire
           * likelihood is fully parameterized at the outset. */

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

        /* we will force marker allele frequency loop to execute at least once */
        for (mkrFreqIdx = 0; mkrFreqIdx == 0 || mkrFreqIdx < modelRange->nafreq; mkrFreqIdx++) {
          mkrFreq = pLocus2->pAlleleFrequency[0];
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
          //      fprintf (stderr, "mkrFreq =%d  model nafreq= %d \n",mkrFreqIdx, modelRange->nafreq);


          if (1 && modelOptions->markerAnalysis == FALSE) {

            if (modelOptions->polynomial == TRUE);
            else
              update_locus (&pedigreeSet, loc1);
          }

          /* clear Dprime combination impossible flag */
          memset (pLambdaCell->impossibleFlag, 0, sizeof (int) * pLambdaCell->ndprime);
          /* set up haplotype frequencies */
	  if(modelOptions->equilibrium == LINKAGE_DISEQUILIBRIUM) {
          for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
            if (isDPrime0 (pLambdaCell->lambda[dprimeIdx], pLambdaCell->m, pLambdaCell->n))
              dprime0Idx = dprimeIdx;
            /* statusLocal = setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprimeIdx);
             * if (statusLocal < 0) {
             * pLambdaCell->impossibleFlag[dprimeIdx] = 1;
             * }      moved to each HLR calculation function */
          }
	  }
	  else
	    dprime0Idx=0;

          /* for each D prime and theta, print out average and maximizing model information - MOD */
          if (modelOptions->markerAnalysis == FALSE)
            dk_write2ptBRHeader (loc1, loc2);

          /* analysis specific statistic initialization */
          if (modelOptions->mapFlag == SA) {
            num_BR = num_sample_Dp_theta;       // currently 271
          } else {      ///  This is for sec-specific analysis in four regions
            num_BR = num_sample_SS_theta;       // currenlty 260
          }
          max_scale = 0;
          CALCHOKE (BRscale, (size_t) num_BR, sizeof (int), int *);
          /*The main loop to Calculate BR(theta, dprime) or BR(thetaM, thetaF) */
          for (i = 0; i < num_BR; i++) {        /* num_BR = 271 for Sex-Average Analysis
                                                 * = 260 for Sex-Specific Analysis */

            if (modelOptions->mapFlag == SA) {
              fixed_dprime = dcuhre2[i][0];
              fixed_theta = dcuhre2[i][1];

              if (fpIR != NULL) {
                if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
                  dk_curModel.dprime[0] = fixed_dprime;
                dk_curModel.theta[0] = dk_curModel.theta[1] = fixed_theta;
              }
            } else {
              fixed_thetaM = thetaSS[i][0];
              fixed_thetaF = thetaSS[i][1];

              if (fpIR != NULL) {
                if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
                  dk_curModel.dprime[0] = 0.0;
                dk_curModel.theta[0] = fixed_thetaM;
                dk_curModel.theta[1] = fixed_thetaF;
              }
              //fprintf (stderr, "i=%d Dprime=%f theta=%f   loc1=%d  loc2=%d\n", i, fixed_thetaM, fixed_thetaF,loc1,loc2);
            }

            integral = 0.0;
            abserr = 0.0;

            if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {     // checking Dprime
              for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
                if (fabs (pLambdaCell->lambda[dprimeIdx][0][0] - fixed_dprime) < 0.0001) {
                  // fprintf(stderr,"dprimeIdx =%d with %15.13f which is matching wit fixed_dprime\n",dprimeIdx,pLambdaCell->lambda[dprimeIdx][0][0]);
                  break;
                }
              }
              if (dprimeIdx == pLambdaCell->ndprime) {
                // &&& This needs fixin'
                ERROR ( "dprimeIdx is %d for dprime=%f theta=%f\n", dprimeIdx,fixed_dprime, fixed_theta);
              }
            }
            num_out_constraint = 0;

            /* Call DCUHRE  Domain information is stored in global variables,  xl an xu */
            kelvin_dcuhre_integrate (&integral, &abserr, volume_region, &(BRscale[i]));
            ASSERT ((s->ifail == 0), "Dynamic integration failed with ifail of %d. Please increase the maxcls parameter in integrationSupport.c if ifail is 1. Others, check dchhre function in dcuhre.c", s->ifail);

            if (modelOptions->mapFlag == SA) {
              dcuhre2[i][3] = integral;
            } else {
              thetaSS[i][3] = integral;
            }
            if (BRscale[i] > max_scale) {
              max_scale = BRscale[i];
            }

            /* Dk specific results */
            if (fpDK != NULL) {
              if (modelOptions->mapFlag == SA) {
                fprintf (fpDK, "%d %6.4f %6.4f %6d %8.4f %8.4f %5d %8.4f\n", i, fixed_dprime, fixed_theta, s->total_neval, integral, abserr, BRscale[i], localMOD);
              } else {
                fprintf (fpDK, "%d %6.4f %6.4f %6d %8.4f %8.4f %5d %8.4f\n", i, fixed_thetaM, fixed_thetaF, s->total_neval, integral, abserr, BRscale[i], localMOD);
              }
              fflush (fpDK);
            }

            R_square = 0.0;     // tp_result[dprimeIdx][thetaInd][modelRange->nafreq].R_square;


            if (overallMOD < localMOD) {
              overallMOD = localMOD;
              if (modelOptions->mapFlag == SA) {
                if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
                  maxima_x[0] = dk_globalmax.dprime[0] = fixed_dprime;
                maxima_x[1] = dk_globalmax.theta[0] = dk_globalmax.theta[1] = fixed_theta;
              } else {
                if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
                  dk_globalmax.dprime[0] = 0.0;
                maxima_x[0] = dk_globalmax.theta[0] = fixed_thetaM;
                maxima_x[1] = dk_globalmax.theta[1] = fixed_thetaF;
              }
              dk_copyMaxModel (localmax_x, &dk_globalmax, size_BR);
              memcpy (&(maxima_x[2]), localmax_x, sizeof (double) * 18);
            }
	    if (overallMin > localMOD)
	      overallMin = localMOD;

            /* If LD, and (D' less then stored D', or BR larger than stored BR */
            if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
              double delta_dprime = dk_dprime0max.dprime[0] - fixed_dprime;
              if ((fabs (fixed_dprime) < fabs (dk_dprime0max.dprime[0]) - 1e-9) || ((fabs (delta_dprime) <= 1e-9) && (dprime0_MOD < localMOD))) {
                dprime0_MOD = localMOD;
                dk_dprime0max.dprime[0] = fixed_dprime;
                dk_dprime0max.theta[0] = dk_dprime0max.theta[1] = fixed_theta;
                dk_copyMaxModel (localmax_x, &dk_dprime0max, size_BR);
              }
            }

            if (modelOptions->mapFlag == SA) {
              double delta_theta = dk_theta0max.theta[0] - fixed_theta;
              /* If fixed_theta is closer to 0 than the current minimum theta,
               * or if fixed_theta is more or less the same and the new BR
               * is greater than the max BR */
              if ((fixed_theta < dk_theta0max.theta[0] - 1e-9) || ((fabs (delta_theta) <= 1e-9) && (theta0_MOD < localMOD))) {
                theta0_MOD = localMOD;
                if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
                  dk_theta0max.dprime[0] = fixed_dprime;
                dk_theta0max.theta[0] = dk_theta0max.theta[1] = fixed_theta;
                dk_copyMaxModel (localmax_x, &dk_theta0max, size_BR);
              }

            } else if (modelOptions->mapFlag == SS) {
              double delta_thetaM = dk_theta0max.theta[0] - fixed_thetaM, delta_thetaF = dk_theta0max.theta[1] - fixed_thetaF;
              if ((fixed_thetaM < dk_theta0max.theta[0] - 1e-9 && fabs (delta_thetaF) <= 1e-9) ||
                  (fixed_thetaF < dk_theta0max.theta[1] - 1e-9 && fabs (delta_thetaM) <= 1e-9) ||
                  (fixed_thetaM < dk_theta0max.theta[0] - 1e-9 && fixed_thetaF < dk_theta0max.theta[1] - 1e-9) || (fabs (delta_thetaM) <= 1e-9 && fabs (delta_thetaF) <= 1e-9 && theta0_MOD < localMOD)) {
                theta0_MOD = localMOD;
                if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
                  dk_theta0max.dprime[0] = 0.0;
                dk_theta0max.theta[0] = fixed_thetaM;
                dk_theta0max.theta[1] = fixed_thetaF;
                dk_copyMaxModel (localmax_x, &dk_theta0max, size_BR);
              }
            }
            if ((modelOptions->mapFlag == SA) && (modelOptions->equilibrium == LINKAGE_EQUILIBRIUM) && (i == 9)) {
              //fprintf (stderr,"End of LE case\n");
              i = num_BR;
            }
          }     /* end of for to calculate BR(theta, dprime) or BR(thetaM, thetaF) */

          /* analysis specific statistic initialization */
          if (modelOptions->mapFlag == SA) {
            le_small_theta = 0.0;
            le_big_theta = 0.0;
            ld_small_theta = 0.0;
            ld_big_theta = 0.0;
            ld_unlinked = 0.0;
            le_unlinked = 0.0;

          } else {      ///  This is for sec-specific analysis in four regions
            thetaSMSF = 0.0;    // 0< thetaM <0.05  0< thetaF <0.05  
            thetaBMSF = 0.0;    // 0.05< thetaM <0.5  0< thetaF <0.05  
            thetaSMBF = 0.0;    // 0< thetaM <0.05  0.05< thetaF <0.5  
            thetaBMBF = 0.0;    // 0.05< thetaM <0.05  0.05< thetaF <0.5 
          }
          for (i = 0; i < num_BR; i++) {        /* num_BR = 271 for Sex-Average Analysis
                                                 * = 260 for Sex-Specific Analysis */
            if (modelOptions->mapFlag == SA) {
              /* Use unifor scaling with max_scale */
              if ((BRscale[i] > max_scale) && (dcuhre2[i][3] > 0)) {
                newLog10BR = log10 (dcuhre2[i][3]) + BRscale[i] - max_scale;
                if (newLog10BR < DBL_MIN_10_EXP + 1) {
                  dcuhre2[i][3] = 0;
                } else {
                  dcuhre2[i][3] = pow (10, newLog10BR);
                }
              }

              if (i < 5) {
                le_small_theta += dcuhre2[i][3] * dcuhre2[i][2];
              } else if (i < 10) {
                le_big_theta += dcuhre2[i][3] * dcuhre2[i][2];
              } else if (i < 270) {
                if (dcuhre2[i][1] < modelOptions->thetaCutoff[0]) {
                  ld_small_theta += dcuhre2[i][3] * dcuhre2[i][2];
                } else {
                  ld_big_theta += dcuhre2[i][3] * dcuhre2[i][2];
                }
              } else {
                le_unlinked += dcuhre2[i][3] * dcuhre2[i][2];
              }
              if (modelOptions->markerAnalysis == FALSE)
                dk_write2ptBRData (dcuhre2[i][0], dcuhre2[i][1], dcuhre2[i][1], dcuhre2[i][3], max_scale);

              if ((modelOptions->equilibrium == LINKAGE_EQUILIBRIUM) && (i == 9)) {
                //fprintf (stderr,"End of LE case\n");
                i = 271;
              }
            } else {
              /* Use unifor scaling with max_scale */
              if ((BRscale[i] > max_scale) && (thetaSS[i][3] > 0)) {
                newLog10BR = log10 (thetaSS[i][3]) + BRscale[i] - max_scale;
                if (newLog10BR < DBL_MIN_10_EXP + 1) {
                  thetaSS[i][3] = 0;
                } else {
                  thetaSS[i][3] = pow (10, newLog10BR);
                }
              }
              if (i < 65) {
                thetaSMSF += thetaSS[i][3] * thetaSS[i][2];     // 0< thetaM <0.05  0< thetaF <0.05  
              } else if (i < 130) {
                thetaBMSF += thetaSS[i][3] * thetaSS[i][2];     // 0.05< thetaM <0.5  0< thetaF <0.05  
              } else if (i < 195) {
                thetaSMBF += thetaSS[i][3] * thetaSS[i][2];     // 0< thetaM <0.05  0.05< thetaF <0.5  
              } else {
                thetaBMBF += thetaSS[i][3] * thetaSS[i][2];     // 0.05< thetaM <0.05  0.05< thetaF <0.5 
              }
              if (modelOptions->markerAnalysis == FALSE)
                dk_write2ptBRData (0, thetaSS[i][0], thetaSS[i][1], thetaSS[i][3], max_scale);
            }
          }     /* end of for to calculate BR(theta, dprime) or BR(thetaM, thetaF) */
          free (BRscale);

          dk_write2ptMODHeader ();
	  if (overallMOD == 0 && overallMin == 0)
	    overallMOD = theta0_MOD = dprime0_MOD = -DBL_MAX;
          dk_write2ptMODData ("MOD(Overall)", overallMOD, &dk_globalmax);

          if (modelOptions->extraMODs) {
            dk_write2ptMODData ("MOD(Theta==0)", theta0_MOD, &dk_theta0max);
            if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
              dk_write2ptMODData ("MOD(D'==0)", dprime0_MOD, &dk_dprime0max);
          }

          /*Calculate ppl, ppld and ldppl */
          if (modelOptions->mapFlag == SS) {
            le_small_theta = thetaSMSF;
            le_big_theta = (0.09 * thetaBMSF + 0.09 * thetaSMBF + 0.81 * thetaBMBF) / 0.99;
          }
          ppl = modelOptions->thetaWeight * le_small_theta + (1 - modelOptions->thetaWeight) * le_big_theta;
          ppl = ppl / (ppl + (1 - modelOptions->prior) / modelOptions->prior);

	  fprintf (fpPPL, "%d", pLocus2->pMapUnit->chromosome);
	  if (modelOptions->markerAnalysis != FALSE)
	    fprintf (fpPPL, " %s %.4f %s %.4f",
		     pLocus1->sName, pLocus1->pMapUnit->mapPos[SEX_AVERAGED],
		     pLocus2->sName, pLocus2->pMapUnit->mapPos[SEX_AVERAGED]);
	  else
	    {
	      fprintf (fpPPL, " %s %s %.4f", pLocus1->sName, pLocus2->sName,
		       pLocus2->pMapUnit->mapPos[SEX_AVERAGED]);
	      if (modelOptions->physicalMap)
		fprintf (fpPPL, " %d", pLocus2->pMapUnit->basePairLocation);
	    }
	  fprintf (fpPPL, " %.*f", ppl >= .025 ? 2 : 3, KROUND (ppl, 3));
	  
          //printf("%f %f %f %f %f %f\n",le_small_theta,le_big_theta,ld_small_theta,ld_big_theta,le_unlinked ,modelOptions->thetaCutoff[0]);
          /* output LD-PPL now if needed */
          if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
            ldppl = 0.019 * (0.021 * ld_small_theta + 0.979 * le_small_theta) + 0.001 * (0.0011 * ld_big_theta + 0.9989 * le_big_theta);

            ldppl = ldppl / (ldppl + 0.98 * le_unlinked);

            ppld = 0.019 * 0.021 * ld_small_theta + 0.001 * 0.0011 * ld_big_theta;
            ppld = ppld / (ppld + 0.019 * 0.979 * le_small_theta + 0.001 * 0.9989 * le_big_theta + 0.98 * le_unlinked);
            //ppld= pow(0.02*(ld_small_theta *0.95 +ld_big_theta * 0.05), 0.3);
            //ppld *=  0.0004/0.020392;
            //ppld = ppld/(ppld + 1-0.0004/0.020392 * 0.02);

            ppldGl = 0.019 * 0.021 * ld_small_theta + 0.001 * 0.0011 * ld_big_theta;
            ppldGl = ppldGl / (ppldGl + 0.019 * 0.979 * le_small_theta + 0.001 * 0.9989 * le_big_theta);

            fprintf (fpPPL, " %.*f", ldppl >= .025 ? 2 : 4, KROUND (ldppl, 4));
            fprintf (fpPPL, " %.*f", ppldGl >= .025 ? 2 : 4, KROUND (ppldGl, 4));
            fprintf (fpPPL, " %.*f", ppld >= .025 ? 2 : 4, KROUND (ppld, 4));
          }
          fprintf (fpPPL, "\n");
          fflush (fpPPL);

          /* only loop marker allele frequencies when doing LD */
          if (modelOptions->equilibrium == LINKAGE_EQUILIBRIUM)
            break;
          /* we can only do SNPs when looping over marker allele frequency */
          if (pLocus2->numOriginalAllele > 2)
            break;
        }
        /* end of marker allele frequency looping */
        /* need to clear polynomial */
        if (modelOptions->polynomial) {
          pedigreeSetPolynomialClearance (&pedigreeSet);
        }


        if (modelOptions->markerAnalysis == ADJACENTMARKER)
          loc2 = originalLocusList.numLocus;
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
  }
  /* end of two point */
  else {
    /* multipoint */
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

    /* calculate the trait likelihood independent of the trait position */
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
    statusLocal = populate_xmission_matrix (traitMatrix, 1, initialProbAddr,    /* probability */
        initialProbAddr2,       /* probability */
        initialHetProbAddr, 0,  /* cell index */
        -1,     /* last he locus */
        -1,     /* last het pattern (P-1 or M-2) */
        0);     /* current locus - start with 0 */

    if (modelOptions->polynomial == TRUE) {
      holdAllPolys ();
    }

    /* for trait likelihood */
    analysisLocusList = &traitLocusList;
    xmissionMatrix = traitMatrix;
    if (pTrait->type == DICHOTOMOUS) {
      
      /*call compute_likelihood with dummy numbers to build polynomials */
      for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
        pTrait->penetrance[2][liabIdxLocal][0][0] = 0.7;
        pTrait->penetrance[2][liabIdxLocal][0][1] = 0.5;
        pTrait->penetrance[2][liabIdxLocal][1][0] = 0.5;
        pTrait->penetrance[2][liabIdxLocal][1][1] = 0.3;
        pTrait->penetrance[1][liabIdxLocal][0][0] = 1 - 0.7;
        pTrait->penetrance[1][liabIdxLocal][0][1] = 1 - 0.5;
        pTrait->penetrance[1][liabIdxLocal][1][0] = 1 - 0.5;
        pTrait->penetrance[1][liabIdxLocal][1][1] = 1 - 0.3;
      }

      if (modelOptions->polynomial == TRUE);
      else
        /* only need to update trait locus */
        update_penetrance (&pedigreeSet, traitLocus);

      pLocus->pAlleleFrequency[0] = 0.5;
      pLocus->pAlleleFrequency[1] = 1 - 0.5;

      if (modelOptions->polynomial == TRUE);
      else
        update_locus (&pedigreeSet, traitLocus);
      /* get the likelihood for the trait */
      sprintf (partialPolynomialFunctionName, "MDT_LC%d_C%d_P%%sSL%d", modelRange->nlclass, (originalLocusList.ppLocusList[1])->pMapUnit->chromosome, modelOptions->sexLinked);
      if (modelOptions->polynomial == TRUE) {
	cL[0]++; // MP DT trait likelihood
	compute_likelihood (&pedigreeSet);        /* This builds polynomials with dummy numbers */
      }
    } else {    // QT
      pLocus->pAlleleFrequency[0] = 0.5;
      pLocus->pAlleleFrequency[1] = 1 - 0.5;
      update_locus (&pedigreeSet, traitLocus);

      /*call compute_likelihood with dummy numbers to build polynomials */
      for (liabIdxLocal = 0; liabIdxLocal < modelRange->nlclass; liabIdxLocal++) {
        pTrait->means[liabIdxLocal][0][0] = 2.0;
        pTrait->means[liabIdxLocal][0][1] = 1.0;
        pTrait->means[liabIdxLocal][1][0] = 1.0;
        pTrait->means[liabIdxLocal][1][1] = 0.3;
        pTrait->stddev[liabIdxLocal][0][0] = 1.0;
        pTrait->stddev[liabIdxLocal][0][1] = 1.0;
        pTrait->stddev[liabIdxLocal][1][0] = 1.0;
        pTrait->stddev[liabIdxLocal][1][1] = 1.0;

        /* threshold for QT */
        pTrait->cutoffValue[liabIdxLocal] = 0.5;

      } /* liability class Index */
      if (modelOptions->polynomial == TRUE)
        sprintf (partialPolynomialFunctionName, "MQT_LC%d_C%d_P%%sSL%d", modelRange->nlclass, (originalLocusList.ppLocusList[1])->pMapUnit->chromosome, modelOptions->sexLinked);
      else
        update_penetrance (&pedigreeSet, traitLocus);
      compute_likelihood (&pedigreeSet);
      cL[1]++; // MP QT trait likelihood
    }

    /* Copy the polynomials built from above to traitLikelihoodPolynomials */
    if (modelOptions->polynomial == TRUE) {
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        /* save the likelihood at trait */
        pPedigreeLocal = pedigreeSet.ppPedigreeSet[pedIdx];
        pPedigreeLocal->traitLikelihoodPolynomial = pPedigreeLocal->likelihoodPolynomial;
        pPedigreeLocal->traitLikelihoodPolyList = pPedigreeLocal->likelihoodPolyList;
        pPedigreeLocal->likelihoodPolyList = NULL;
        pPedigreeLocal->likelihoodPolynomial = NULL;

        //fprintf(stderr,"Building traitPoly pedIdx =%d Null likelihood = %20.15f\n",pedIdx, pPedigreeLocal->likelihood);
        //fprintf(stderr,"pedIdx %d eType= %d\n", pedIdx, ((pPedigreeLocal->traitLikelihoodPolyList)->pList[0])->eType);
      }
    }


    /* get the trait locations we need to evaluate at */
    numPositions = modelRange->ntloc;

    /* Need to output the results */
    dk_writeMPBRHeader ();
    dk_writeMPMODHeader ();

    prevFirstMarker = -1;
    prevLastMarker = -1;
    prevTraitInd = -1;
    leftMarker = -1;

#ifdef STUDYDB

    /* Iterate over all positions in the analysis. */

    if (toupper(*studyDB.role) == 'S') {
      // We're a server! Completely suborn the trait loci vector in modelRange
      int i, j = 0;

      // Marker positions are themselves transition positions...
      for (i=1; i<originalLocusList.numLocus; i++)
	lociSetTransitionPositions[j++] = originalLocusList.ppLocusList[i]->pMapUnit->mapPos[0];

      // ...as are the midpoints of marker N and N+M, where M is the number of markers in the analysis
      for (i=0; i<(originalLocusList.numLocus - modelType->numMarkers); i++)
	lociSetTransitionPositions[j++] = lociSetTransitionPositions[i] +
	  ((lociSetTransitionPositions[i+modelType->numMarkers] - lociSetTransitionPositions[i])/2.0);

      // Sort the list of transition positions
      qsort (lociSetTransitionPositions, j, sizeof (double), compare_doubles);

      // Choose new trait positions between transition positions

      newTLoc[0] = lociSetTransitionPositions[0] - 1.0;
      for (i=1; i<j; i++) {
	newTLoc[i] = lociSetTransitionPositions[i-1] + ((lociSetTransitionPositions[i] - lociSetTransitionPositions[i-1]) / 2.0);
      }
      newTLoc[j] = lociSetTransitionPositions[j-1] + 1.0;

      oldTLoc = modelRange->tloc;
      modelRange->tloc = newTLoc;
      modelRange->ntloc = j+1;
      numPositions = j+1;
    }

    DIAG (ALTLSERVER, 1, {		 \
	for (i=0; i<numPositions; i++)					\
	  fprintf (stderr, "nTL[%d] is %.6g\n", i, newTLoc[i]);});

#endif

    CALCHOKE (mp_result, (size_t) numPositions, sizeof (SUMMARY_STAT), SUMMARY_STAT *);

    for (posIdx = 0; posIdx < numPositions; posIdx++) {

      if (fpIR != NULL)
        dk_curModel.posIdx = posIdx;

      /* positions listed are sex average positions */
      traitPos = modelRange->tloc[posIdx];

#ifdef STUDYDB
      
      int freeModels = 0;

      studyDB.driverPosIdx = posIdx;

      if (toupper(*studyDB.role) == 'S') {

	double lowPosition  = -99.99, highPosition = 9999.99;

	if (posIdx != 0)
	  lowPosition = lociSetTransitionPositions[posIdx - 1];
	if (posIdx != (modelRange->ntloc - 1))
	  highPosition = lociSetTransitionPositions[posIdx];

	// If we have models to work on, say how many, otherwise say we're skipping this position

	if ((freeModels = CountWork(lowPosition, highPosition)) == 0) {
	  SUBSTEP (posIdx * 100 / numPositions, "Skipping position %d (%.4gcM from %.4gcM to %.4gcM) of %d (no work)", posIdx + 1, traitPos, lowPosition, highPosition, numPositions);
	  continue;
	} else
	  SUBSTEP (posIdx * 100 / numPositions, "Starting with position %d (%.4gcM from %.4gcM to %.4gcM) of %d (%d available models)", posIdx + 1, traitPos, lowPosition, highPosition, numPositions, freeModels);
      }
#else
      SUBSTEP (posIdx * 100 / numPositions, "Starting with position %d (%.4gcM) of %d", posIdx + 1, traitPos, numPositions);
#endif

      /* set the sex average position first 
       * the sex specific positions will be updated once markers are selected
       * as interpolation might be needed
       */
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
      markerSetChanged = FALSE;
#ifdef STUDYDB
      if (TRUE) { // Marker set must change for every position because we don't know when it does for all study maps
#else
      if (prevFirstMarker != mp_result[posIdx].pMarkers[0]
          || prevLastMarker != mp_result[posIdx].pMarkers[modelType->numMarkers - 1]) {
#endif
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

        /* calculate likelihood for the marker set */
        analysisLocusList = &markerLocusList;
        xmissionMatrix = markerMatrix;
        if (modelOptions->polynomial == TRUE) {
          pedigreeSetPolynomialClearance (&pedigreeSet);
        }

        /* populate the matrix */
        statusLocal = populate_xmission_matrix (markerMatrix, markerLocusList.numLocus, initialProbAddr,        /* probability */
            initialProbAddr2,   /* probability */
            initialHetProbAddr, 0,      /* cell index */
            -1, /* last he locus */
            -1, /* last het pattern (P-1 or M-2) */
            0); /* current locus - start with 0 */

        if (modelOptions->polynomial == TRUE)
          freePolys ();

        DIAG (XM, 4, {
              print_xmission_matrix (markerMatrix, markerLocusList.numLocus, 0, 0, tmpID);
            }
        );

        char markerNo[8];
        sprintf (partialPolynomialFunctionName, "MM_LC%d_C%d_P%%sM", modelRange->nlclass, (originalLocusList.ppLocusList[1])->pMapUnit->chromosome);
        for (k = 0; k < modelType->numMarkers; k++) {
          sprintf (markerNo, "_%d", markerLocusList.pLocusIndex[k]);
          strcat (partialPolynomialFunctionName, markerNo);
        }
        cL[2]++; // MP marker likelihood
        compute_likelihood (&pedigreeSet);

        /* save the results for marker likelihood */
        for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
          /* save the likelihood at null */
          pPedigreeLocal = pedigreeSet.ppPedigreeSet[pedIdx];
          pPedigreeLocal->markerLikelihood = pPedigreeLocal->likelihood;
        }
        pedigreeSet.markerLikelihood = pedigreeSet.likelihood;
        pedigreeSet.log10MarkerLikelihood = pedigreeSet.log10Likelihood;
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
      /*end of SS model */

      /* the locus list has been built, go on to the analysis 
       * multipoint DT */
      if (markerSetChanged || locusListChanged) {
        if (modelOptions->polynomial == TRUE) {
          pedigreeSetPolynomialClearance (&pedigreeSet);
          /* populate the matrix */
          statusLocal = populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr,        /* probability */
              initialProbAddr2, /* probability */
              initialHetProbAddr, 0,    /* cell index */
              -1, -1,   /* last het locus & last het pattern (P-1 or M-2) */
              0);       /* current locus - start with 0 */
          DIAG (XM, 4, {
                print_xmission_matrix (altMatrix, savedLocusList.numLocus, 0, 0, tmpID);
              }
          );
          if (modelOptions->polynomial == TRUE)
            freePolys ();
        }
      }

      if (modelOptions->polynomial != TRUE);
      /* populate the matrix */
      statusLocal = populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr,    /* probability */
          initialProbAddr2,     /* probability */
          initialHetProbAddr, 0,        /* cell index */
          -1, -1,       /* last het locus & last het pattern (P-1 or M-2) */
          0);   /* current locus - start with 0 */

      /* multipoint DT */
      integral = 0.0;
      abserr = 0.0;
      num_out_constraint = 0;
      kelvin_dcuhre_integrate (&integral, &abserr, volume_region, &max_scale);
      ASSERT ((s->ifail == 0), "Dynamic integration failed with ifail of %d. Please increase the maxcls parameter in integrationSupport.c if ifail is 1. Others, check dchhre function in dcuhre.c", s->ifail);

      num_eval = s->total_neval;

      /* calculate imputed PPL and print the results */
      if (integral > 0.214) {
        if ((log10 (integral) + max_scale) > 8)
          ppl = 1.0;
        else
          ppl = (integral * integral) / (-5.77 + 54 * integral + integral * integral);
      } else
        ppl = 0;

      dk_writeMPBRData (posIdx, traitPos, ppl, integral, max_scale);
#ifdef STUDYDB
      if (studyDB.bogusLikelihoods > 0)
	fprintf (fpHet, "WARNING - Some positions have not been completely analyzed!\n");
#endif

      dk_copyMaxModel (localmax_x, &dk_globalmax, size_BR);
      dk_writeMPMODData (posIdx, traitPos, localMOD, &dk_globalmax);

      if (fpDK != NULL) {
        fprintf (fpDK, "%f  %6.4f %12.8f %12.8f %d %d %f\n", traitPos, ppl, integral, abserr, num_eval, max_scale, localMOD);
        fflush (fpDK);
      }

    }   /* end of walking down the chromosome */
  }     /* end of multipoint */

  DIAG (OVERALL, 1, {
        dumpTrackingStats (cL, eCL);
      }
  );

  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
    free (dk_globalmax.dprime);
    free (dk_dprime0max.dprime);
    free (dk_theta0max.dprime);
    if (fpIR != NULL)
      free (dk_curModel.dprime);
  }
  free (dk_globalmax.pen);
  free (dk_dprime0max.pen);
  free (dk_theta0max.pen);
  if (fpIR != NULL)
    free (dk_curModel.pen);

  free (xu);
  free (xl);
}
