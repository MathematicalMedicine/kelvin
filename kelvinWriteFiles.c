#include "kelvin.h"
//#include "kelvinLocals.h"
#include "kelvinGlobalsNew.h"
#include "summary_result.h"
#include "ppl.h"

void
writePPLFileHeader ()
{
  fprintf (fpPPL, "# Version %s\n", programVersion);
  if (modelOptions.markerAnalysis != FALSE)
    {
      fprintf (fpPPL, "Chr Marker1 Position1 Marker2 Position2 PPL");
      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
	{
	  fprintf (fpPPL, " PPL(LD) PPLD|L PPLD(L)");
	}
    }
  else
    {
      fprintf (fpPPL, "Chr Trait Marker Position PPL");
      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
	{
	  fprintf (fpPPL, " PPL(LD) PPLD|L PPLD(L)");
	}
    }
  fprintf (fpPPL, "\n");
  fflush (fpPPL);
}

void
writePPLFileDetail (int dprime0Idx)
{
  LDVals ldvals;
  double ldstat;
  float ppl;

  /* Chromosome, marker name, position, PPL */
  ppl = calculate_PPL (tp_result[dprime0Idx]);
  if (modelOptions.markerAnalysis != FALSE)
    {
      fprintf (fpPPL, "%d %s %.4f %s %.4f %.3f ",
	       pLocus2->pMapUnit->chromosome, pLocus1->sName,
	       pLocus1->pMapUnit->mapPos[SEX_AVERAGED], pLocus2->sName,
	       pLocus2->pMapUnit->mapPos[SEX_AVERAGED], ppl);
    }
  else
    {
      fprintf (fpPPL, "%d %s %s %.4f %.3f ",
	       pLocus2->pMapUnit->chromosome, pLocus1->sName,
	       pLocus2->sName, pLocus2->pMapUnit->mapPos[SEX_AVERAGED], ppl);
    }
  fflush (fpPPL);
  /* output LD-PPL now if needed */
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
    {
      /* load up ldvals first */
      get_LDVals (tp_result, &ldvals);
      ldstat = calc_ppl_allowing_ld (&ldvals, modelOptions.LDprior);
      fprintf (fpPPL, "%.*f ", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
      ldstat = calc_ppld_given_linkage (&ldvals, modelOptions.LDprior);
      fprintf (fpPPL, "%.*f ", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
      ldstat = calc_ppld_allowing_l (&ldvals, modelOptions.LDprior);
      fprintf (fpPPL, "%.*f ", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
    }
  fprintf (fpPPL, "\n");
  fflush (fpPPL);
}


void
write2ptBRFile (int loc1, int loc2)
{
  int i, j;
  int dprimeIdx, thetaInd;
  float theta[2];
  int exponent, scale;
  double base;
  double lr, log10LR;

  /*
  KASSERT (modelOptions.markerAnalysis == FALSE,
	   "Don't call write2ptBRFile for a marker-to-marker analysis\n");
  */

  /* Print the marker line: sequence number, trait name, marker name, and chromosome.
   * If the map is explicitly sex-specific, and the male and female position fields
   * contain non-negative values, print the average, female and male cM positions;
   * otherwise just print the average position. If the base pair field contains a
   * non-negative value, print it, too
   */
  fprintf (fpHet, "# Seq: %d Chr: %d Trait: %s Marker: %s", loc2,
	   pLocus2->pMapUnit->chromosome, pLocus1->sName, pLocus2->sName);
  if ((modelOptions.mapFlag == SEX_SPECIFIC)
      && (pLocus2->pMapUnit->mapPos[MAP_FEMALE] >= 0)
      && (pLocus2->pMapUnit->mapPos[MAP_MALE] >= 0))
    {
      fprintf (fpHet,
	       " AvgPosition: %.4f FemalePosition: %.4f MalePosition: %.4f",
	       pLocus2->pMapUnit->mapPos[MAP_SEX_AVERAGE],
	       pLocus2->pMapUnit->mapPos[MAP_FEMALE],
	       pLocus2->pMapUnit->mapPos[MAP_MALE]);
    }
  else
    {
      fprintf (fpHet, " Position: %.4f",
	       pLocus2->pMapUnit->mapPos[MAP_SEX_AVERAGE]);
    }
  if (pLocus2->pMapUnit->basePairLocation >= 0)
    fprintf (fpHet, " Phyiscal %d", pLocus2->pMapUnit->basePairLocation);
  fprintf (fpHet, "\n");

  /* For each D prime and theta, print out average and maximizing model information - MOD */

  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
    for (i = 0; i < pLocus1->numOriginalAllele - 1; i++)
      for (j = 0; j < pLocus2->numOriginalAllele - 1; j++)
	fprintf (fpHet, "D%1d%1d ", i + 1, j + 1);
  fprintf (fpHet, "Theta(M,F) BayesRatio\n");

  for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++)
    {
      for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
	{
	  if (tp_result[dprimeIdx][thetaInd][modelRange.nafreq].lr_count == 0)
	    continue;
	  theta[0] = modelRange.theta[0][thetaInd];
	  theta[1] = modelRange.theta[1][thetaInd];
	  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
	    {
	      for (i = 0; i < pLocus1->numOriginalAllele - 1; i++)
		for (j = 0; j < pLocus2->numOriginalAllele - 1; j++)
		  {
		    fprintf (fpHet, "%.2f ",
			     pLambdaCell->lambda[dprimeIdx][i][j]);
		  }
	    }
	  lr = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].het_lr_avg_orig;
	  if(lr > 0) {
	    log10LR=log10(lr);
	    exponent=floor(log10LR);
	    scale=tp_result[dprimeIdx][thetaInd][modelRange.nafreq].scale_orig;
	    base=pow(10, (log10LR-exponent));
	    fprintf (fpHet, "(%.4f,%.4f) %.6fe%+d\n", theta[0], theta[1],
		     base, exponent+scale);
	  }
	  else {
	    fprintf (fpHet, "(%.4f,%.4f) %.6fe%+d\n", theta[0], theta[1],
		     0.0, 0);
	  }
	}			/* theta loop */
    }				/* dprime loop */
}


void
writeMPBRFileHeader ()
{
  int k;

  /* Need to output the results */
  fprintf (fpHet, "Chr Position PPL BayesRatio");
  fprintf (fpHet, " MarkerList(0");
  for (k = 1; k < modelType.numMarkers; k++)
    fprintf (fpHet, ",%d", k);
  fprintf (fpHet, ")\n");
}


void
writeMPBRFileDetail (int posIdx, float traitPos, float ppl, double avgLR)
{
  int k;

  fprintf (fpHet, "%d %f %.*f %.6e",
	   (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->
	   pMapUnit->chromosome, traitPos, ppl >= .025 ? 2 : 3,
	   ppl >=
	   .025 ? rint (ppl * 100.) / 100. : rint (ppl * 1000.) / 1000.,
	   avgLR);

  /* print out markers used for this position */
  fprintf (fpHet, " (%d", mp_result[posIdx].pMarkers[0]);
  for (k = 1; k < modelType.numMarkers; k++)
    {
      fprintf (fpHet, ",%d", mp_result[posIdx].pMarkers[k]);
    }
  fprintf (fpHet, ")\n");
  fflush (fpHet);
}


void
writeMPMODFileHeader ()
{
  int liabIdx;

  if (fpMOD == NULL)
    return;

  fprintf (fpMOD, "Chr Position MOD Alpha DGF");

  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
    if (modelType.trait == DT)
      if (modelOptions.imprintingFlag)
	fprintf (fpMOD, " LC%dPV(DD,Dd,dD,dd)", liabIdx);
      else
	fprintf (fpMOD, " LC%dPV(DD,Dd,dd)", liabIdx);
    else
      {
	if (modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	  if (modelOptions.imprintingFlag)
	    fprintf (fpMOD,
		     " LC%dPV(DDMean,DdMean,dDMean,ddMean,DDSD,DdSD,dDSD,ddSD",
		     liabIdx);
	  else
	    fprintf (fpMOD, " LC%dPV(DDMean,DdMean,ddMean,DDSD,DdSD,ddSD",
		     liabIdx);
	else if (modelOptions.imprintingFlag)
	  fprintf (fpMOD, " LC%dPV(DDDF,DdDF,dDDF,ddDF", liabIdx);
	else
	  fprintf (fpMOD, " LC%dPV(DDDF,DdDF,ddDF", liabIdx);
	if (modelType.trait == CT)
	  fprintf (fpMOD, ",Thresh)");
	else
	  fprintf (fpMOD, ")");
      }

  fprintf (fpMOD, "\n");
  fflush (fpMOD);
  return;
}


void
writeMPMODFileDetail (int posIdx, float traitPos)
{
  int liabIdx;
  float pen_DD, pen_Dd, pen_dD, pen_dd;
  float SD_DD, SD_Dd, SD_dD, SD_dd;
  float threshold;
  double max;
  float gfreq, alphaV;
  int penIdx, paramIdx, thresholdIdx;

  if (fpMOD == NULL)
    return;

  max = mp_result[posIdx].max_log10_lr;
  gfreq = mp_result[posIdx].max_gfreq;
  alphaV = mp_result[posIdx].max_alpha;
  penIdx = mp_result[posIdx].max_penIdx;
  paramIdx = mp_result[posIdx].max_paramIdx;
  thresholdIdx = mp_result[posIdx].max_thresholdIdx;

  fprintf (fpMOD, "%d %f %.6f %f %f",
	   (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->
	   pMapUnit->chromosome, traitPos, max, alphaV, gfreq);
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
    {
      pen_DD = modelRange.penet[liabIdx][0][penIdx];
      pen_Dd = modelRange.penet[liabIdx][1][penIdx];
      pen_dD = modelRange.penet[liabIdx][2][penIdx];
      pen_dd = modelRange.penet[liabIdx][3][penIdx];
      if (modelOptions.imprintingFlag)
	fprintf (fpMOD, " (%.3f,%.3f,%.3f,%.3f", pen_DD, pen_Dd, pen_dD,
		 pen_dd);
      else
	fprintf (fpMOD, " (%.3f,%.3f,%.3f", pen_DD, pen_Dd, pen_dd);
      if (modelType.trait != DT
	  && modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	{
	  SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
	  SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
	  SD_dD = modelRange.param[liabIdx][2][0][paramIdx];
	  SD_dd = modelRange.param[liabIdx][3][0][paramIdx];
	  if (modelOptions.imprintingFlag)
	    fprintf (fpMOD, ",%.3f,%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dD,
		     SD_dd);
	  else
	    fprintf (fpMOD, ",%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dd);
	}
      if (modelType.trait == CT)
	{
	  threshold = modelRange.tthresh[liabIdx][thresholdIdx];
	  fprintf (fpMOD, ",%.3f)", threshold);
	}
      else
	fprintf (fpMOD, ")");
    }
  fprintf (fpMOD, "\n");
  return;
}


void
writeMaximizingModel (char *modelDescription, double myMOD, int myDPrimeIdx,
		      int myThetaIdx)
{
  float theta[2], gfreq, mkrFreq, alphaV, R_square;
  int penIdx, paramIdx, thresholdIdx, liabIdx;
  float pen_DD, pen_Dd, pen_dD, pen_dd;
  float SD_DD, SD_Dd, SD_dD, SD_dd;
  float threshold;
  int i, j; 

  if (modelOptions.extraMODs)
    fprintf (fpMOD, "%s ", modelDescription);
  fprintf (fpMOD, "%.4f", myMOD);

  theta[0] = modelRange.theta[0][myThetaIdx];
  theta[1] = modelRange.theta[1][myThetaIdx];
  gfreq = tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].max_gfreq;
  mkrFreq = tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].max_mf;
  alphaV = tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].max_alpha;
  penIdx = tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].max_penIdx;
  R_square = tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].R_square;
  paramIdx =
    tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].max_paramIdx;
  thresholdIdx =
    tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].max_thresholdIdx;
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
    for (i = 0; i < pLocus1->numOriginalAllele - 1; i++)
      for (j = 0; j < pLocus2->numOriginalAllele - 1; j++)
	fprintf (fpMOD, " %.2f", pLambdaCell->lambda[myDPrimeIdx][i][j]);
  fprintf (fpMOD, " (%.4f,%.4f)", theta[0], theta[1]);
  if (modelOptions.markerAnalysis != FALSE)
    fprintf (fpMOD, " %.3f", R_square);
  fprintf (fpMOD, " %.2f %.4f", alphaV, gfreq);

  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
    {
      pen_DD = modelRange.penet[liabIdx][0][penIdx];
      pen_Dd = modelRange.penet[liabIdx][1][penIdx];
      pen_dD = modelRange.penet[liabIdx][2][penIdx];
      pen_dd = modelRange.penet[liabIdx][3][penIdx];
      if (modelOptions.imprintingFlag)
	fprintf (fpMOD, " (%.3f,%.3f,%.3f,%.3f", pen_DD, pen_Dd, pen_dD,
		 pen_dd);
      else
	fprintf (fpMOD, " (%.3f %.3f %.3f", pen_DD, pen_Dd, pen_dd);
      if (modelType.trait != DT
	  && modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	{
	  SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
	  SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
	  SD_dD = modelRange.param[liabIdx][2][0][paramIdx];
	  SD_dd = modelRange.param[liabIdx][3][0][paramIdx];
	  if (modelOptions.imprintingFlag)
	    fprintf (fpMOD, ",%.3f,%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dD,
		     SD_dd);
	  else
	    fprintf (fpMOD, ",%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dd);
	}
      if (modelType.trait == CT)
	{
	  threshold = modelRange.tthresh[liabIdx][thresholdIdx];
	  fprintf (fpMOD, ",%.3f)", threshold);
	}
      else
	fprintf (fpMOD, ")");
    }
  fprintf (fpMOD, "\n");
  fflush (fpMOD);
}


void
write2ptMODFile (int loc1, int loc2)
{
  double log10_lr;
  int initialFlag;
  float max, max_at_theta0, max_at_dprime0;
  int maxDPrimeIdx, maxThetaIdx;
  int dprimeIdx, thetaInd;
  int maxDPrimeIdx_at_theta0;
  int maxTheta_at_dprime0;
  int i, j, liabIdx;
  float theta[2];
  int theta0Idx, dprime0Idx;
  
  if (fpMOD == NULL)
    return;

  if (modelOptions.markerAnalysis == FALSE)
    {
      fprintf (fpMOD, "# Seq: %d Chr: %d Trait: %s Marker: %s", loc2,
	       pLocus2->pMapUnit->chromosome, pLocus1->sName, pLocus2->sName);

      if ((modelOptions.mapFlag == SEX_SPECIFIC) &&
	  (pLocus2->pMapUnit->mapPos[MAP_FEMALE] >= 0) &&
	  (pLocus2->pMapUnit->mapPos[MAP_MALE] >= 0))
	{
	  fprintf (fpMOD,
		   " AvgPosition: %.4f FemalePosition: %.4f MalePosition: %.4f",
		   pLocus2->pMapUnit->mapPos[MAP_SEX_AVERAGE],
		   pLocus2->pMapUnit->mapPos[MAP_FEMALE],
		   pLocus2->pMapUnit->mapPos[MAP_MALE]);

	}
      else
	{
	  fprintf (fpMOD, " Position: %.4f",
		   pLocus2->pMapUnit->mapPos[MAP_SEX_AVERAGE]);
	}

    }
  else
    {
      fprintf (fpMOD,
	       "# Seq: %d Chr %d: Marker1: %s Position1: %.4f Marker2: %s Position2: %.4f",
	       loc2, pLocus2->pMapUnit->chromosome, pLocus1->sName,
	       pLocus1->pMapUnit->mapPos[MAP_SEX_AVERAGE], pLocus2->sName,
	       pLocus2->pMapUnit->mapPos[MAP_SEX_AVERAGE]);
    }

  if (pLocus2->pMapUnit->basePairLocation >= 0)
    fprintf (fpMOD, " Phyiscal %d", pLocus2->pMapUnit->basePairLocation);
  fprintf (fpMOD, "\n");

  if (modelOptions.extraMODs)
    fprintf (fpMOD, "Case ");
  fprintf (fpMOD, "MOD");

  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
    for (i = 0; i < pLocus1->numOriginalAllele - 1; i++)
      for (j = 0; j < pLocus2->numOriginalAllele - 1; j++)
	fprintf (fpMOD, " D%1d%1d", i + 1, j + 1);
  fprintf (fpMOD, " Theta(M,F)");
  if (modelOptions.markerAnalysis != FALSE)
    fprintf (fpMOD, " R2");
  fprintf (fpMOD, " Alpha DGF");

  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
    if (modelType.trait == DT)
      if (modelOptions.imprintingFlag)
	fprintf (fpMOD, " LC%dPV(DD,Dd,dD,dd)", liabIdx);
      else
	fprintf (fpMOD, " LC%dPV(DD,Dd,dd)", liabIdx);
    else
      {
	if (modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	  if (modelOptions.imprintingFlag)
	    fprintf (fpMOD,
		     " LC%dPV(DDMean,DdMean,dDMean,ddMean,DDSD,DdSD,dDSD,ddSD",
		     liabIdx);
	  else
	    fprintf (fpMOD, " LC%dPV(DDMean,DdMean,ddMean,DDSD,DdSD,ddSD",
		     liabIdx);
	else if (modelOptions.imprintingFlag)
	  fprintf (fpMOD, " LC%dPV(DDDF,DdDF,dDDF,ddDF", liabIdx);
	else
	  fprintf (fpMOD, " LC%dPV(DDDF,DdDF,ddDF", liabIdx);
	if (modelType.trait == CT)
	  fprintf (fpMOD, ",Thresh)");
	else
	  fprintf (fpMOD, ")");
      }
  fprintf (fpMOD, "\n");

  initialFlag = 1;
  max = -99999;
  max_at_theta0 = -99999;
  max_at_dprime0 = -99999;
  for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++)
    {
      //dprime = pLambdaCell->lambda[dprimeIdx][0][0];
      for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
	{
	  theta[0] = modelRange.theta[0][thetaInd];
	  theta[1] = modelRange.theta[1][thetaInd];
	  log10_lr = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_log10_lr;
	  if (initialFlag || log10_lr > max)
	    {
	      /* overall max */
	      max = log10_lr;
	      maxDPrimeIdx = dprimeIdx;
	      maxThetaIdx = thetaInd;
	    }
	  if (initialFlag
	      || (-ERROR_MARGIN <= theta[0] && theta[0] <= ERROR_MARGIN
		  && -ERROR_MARGIN <= theta[1] && theta[1] <= ERROR_MARGIN))
	    {
	      /* find the max for models with theta equal to 0 */
	      theta0Idx = thetaInd;
	      if (log10_lr > max_at_theta0)
		{
		  max_at_theta0 = log10_lr;
		  maxDPrimeIdx_at_theta0 = dprimeIdx;
		}
	    }
	  if (dprimeIdx == dprime0Idx)
	    {
	      if (initialFlag || maxTheta_at_dprime0 < 0
		  || log10_lr > max_at_dprime0)
		{
		  max_at_dprime0 = log10_lr;
		  maxTheta_at_dprime0 = thetaInd;
		}
	    }
	  initialFlag = 0;
	}
      initialFlag = 0;
    }

  /* Overall maximizing model - MOD */
  writeMaximizingModel ("MOD(Overall)", max, maxDPrimeIdx, maxThetaIdx);

  if (modelOptions.extraMODs)
    {
      /* Maximizing model at theta equal to 0 - MOD */
      writeMaximizingModel ("MOD(Theta==0)", max_at_theta0,
			    maxDPrimeIdx_at_theta0, theta0Idx);

      /* Maximizing model at d prime equal to 0 - MOD */
      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
	writeMaximizingModel ("MOD(D'==0)", max_at_dprime0, dprime0Idx,
			      maxTheta_at_dprime0);
    }

  return;
}
