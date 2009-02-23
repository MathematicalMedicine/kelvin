void writePPLFileHeader () {
  fprintf (fpPPL, "# Version %s\n", programVersion);
  if (modelOptions.markerAnalysis != FALSE) {
    fprintf (fpPPL, "Chr Marker1 Position1 Marker2 Position2 PPL");
    if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
      fprintf (fpPPL, " PPL(LD) PPLD|L PPLD(L) ");
    }
  } else {
    fprintf (fpPPL, "Chr Trait Marker Position PPL");
    if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
      fprintf (fpPPL, " PPL(LD) PPLD|L PPLD(L) ");
    }
  }
    fprintf (fpPPL, "\n");
    fflush (fpPPL);
}

void writePPLFileDetail () {
  LDVals ldvals;
  double ldstat;

  /* Chromosome, marker name, position, PPL */
  ppl = calculate_PPL (tp_result[dprime0Idx]);
  if (modelOptions.markerAnalysis != FALSE) {
    fprintf (fpPPL, "%d %s %.4f %s %.4f %.3f ",
	     pLocus2->pMapUnit->chromosome, pLocus1->sName, pLocus1->pMapUnit->mapPos[SEX_AVERAGED],
	     pLocus2->sName, pLocus2->pMapUnit->mapPos[SEX_AVERAGED], ppl);
  } else {
    fprintf (fpPPL, "%d %s %s %.4f %.3f ",
	     pLocus2->pMapUnit->chromosome, pLocus1->sName,
	     pLocus2->sName, pLocus2->pMapUnit->mapPos[SEX_AVERAGED], ppl);
  }
  fflush (fpPPL);
  /* output LD-PPL now if needed */
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
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


void write2ptBRFile() {

  /* For each D prime and theta, print out average and maximizing model information - MOD */

  fprintf (fpHet, "# %-d  %s %s \n", loc2, pLocus1->sName, pLocus2->sName);
  fprintf (fpHet, "Chr Position ");
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
    for (i = 0; i < pLocus1->numOriginalAllele - 1; i++)
      for (j = 0; j < pLocus2->numOriginalAllele - 1; j++)
	fprintf (fpHet, "D%1d%1d ", i + 1, j + 1);
  fprintf (fpHet, "Theta(M,F) BayesRatio MOD R2 Alpha DGF MF ");
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
    if (modelType.trait == DT)
      if (modelOptions.imprintingFlag)
	fprintf (fpHet, "LC%dPV(DD,Dd,dD,dd) ", liabIdx);
      else
	fprintf (fpHet, "LC%dPV(DD,Dd,dd) ", liabIdx);
    else
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	if (modelOptions.imprintingFlag)
	  fprintf (fpHet, "LC%dPV(DDMean,DdMean,dDMean,ddMean,DDSD,DdSD,dDSD,ddSD,Thresh) ", liabIdx);
	else
	  fprintf (fpHet, "LC%dPV(DDMean,DdMean,ddMean,DDSD,DdSD,ddSD,Thresh) ", liabIdx);
      else
	if (modelOptions.imprintingFlag)
	  fprintf (fpHet, "LC%dPV(DDDF,DdDF,dDDF,ddDF,Thresh) ", liabIdx);
	else
	  fprintf (fpHet, "LC%dPV(DDDF,DdDF,ddDF,Thresh) ", liabIdx);
  fprintf (fpHet, "\n");

  for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
    for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++) {
      if (tp_result[dprimeIdx][thetaInd]
	  [modelRange.nafreq].lr_count == 0)
	continue;
      theta[0] = modelRange.theta[0][thetaInd];
      theta[1] = modelRange.theta[1][thetaInd];
      max = log10 (tp_result[dprimeIdx][thetaInd]
		   [modelRange.nafreq].max_lr);
      gfreq = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_gfreq;
      alphaV = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_alpha;
      penIdx = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_penIdx;
      paramIdx = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_paramIdx;
      thresholdIdx = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_thresholdIdx;
      R_square = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].R_square;
      fprintf (fpHet, "%d %.4f ", pLocus2->pMapUnit->chromosome, pLocus2->pMapUnit->mapPos[SEX_AVERAGED]);
      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
	for (i = 0; i < pLocus1->numOriginalAllele - 1; i++)
	  for (j = 0; j < pLocus2->numOriginalAllele - 1; j++) {
	    fprintf (fpHet, "%.2f ", pLambdaCell->lambda[dprimeIdx][i][j]);
	  }
      }
      fprintf (fpHet, "(%.4f,%.4f) %.6e %.4f %.4f %.2f %.4f %.4f",
	       theta[0], theta[1],
	       tp_result[dprimeIdx][thetaInd][modelRange.nafreq].het_lr_avg, max,
	       tp_result[dprimeIdx][thetaInd][modelRange.nafreq].R_square, alphaV, gfreq,
	       tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_mf);
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	pen_DD = modelRange.penet[liabIdx][0][penIdx];
	pen_Dd = modelRange.penet[liabIdx][1][penIdx];
	pen_dD = modelRange.penet[liabIdx][2][penIdx];
	pen_dd = modelRange.penet[liabIdx][3][penIdx];
	if (modelOptions.imprintingFlag)
	  fprintf (fpHet, " (%.3f,%.3f,%.3f,%.3f", pen_DD, pen_Dd, pen_dD, pen_dd);
	else
	  fprintf (fpHet, " (%.3f,%.3f,%.3f", pen_DD, pen_Dd, pen_dd);
	if (modelType.trait != DT && modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	  SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
	  SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
	  SD_dD = modelRange.param[liabIdx][2][0][paramIdx];
	  SD_dd = modelRange.param[liabIdx][3][0][paramIdx];
	  if (modelOptions.imprintingFlag)
	    fprintf (fpHet, ",%.3f,%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dD, SD_dd);
	  else
	    fprintf (fpHet, ",%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dd);
	}
	if (modelType.trait != DT) {
	  threshold = modelRange.tthresh[liabIdx][thresholdIdx];
	  fprintf (fpHet, ",%.3f)", threshold);
	} else
	  fprintf (fpHet, ")");
      }
      fprintf (fpHet, "\n");
    }     /* theta loop */
  }       /* dprime loop */
}

void writeMPBRFileHeader () {

  /* Need to output the results */
  fprintf (fpHet, "Chr Position PPL BayesRatio MOD Alpha DGF ");
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
    if (modelType.trait == DT)
      if (modelOptions.imprintingFlag)
	fprintf (fpHet, "LC%dPV(DD,Dd,dD, dd) ", liabIdx);
      else
	fprintf (fpHet, "LC%dPV(DD,Dd,dd) ", liabIdx);
    else
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	if (modelOptions.imprintingFlag)
	  fprintf (fpHet, "LC%dPV(DDMean,DdMean,dDMean,ddMean,DDSD,DdSD,dDSD,ddSD,Thresh) ", liabIdx);
	else
	  fprintf (fpHet, "LC%dPV(DDMean,DdMean,ddMean,DDSD,DdSD,ddSD,Thresh) ", liabIdx);
      else
	if (modelOptions.imprintingFlag)
	  fprintf (fpHet, "LC%dPV(DDDF,DdDF,dDF,ddDF,Thresh) ", liabIdx);
	else
	  fprintf (fpHet, "LC%dPV(DDDF,DdDF,ddDF,Thresh) ", liabIdx);
  fprintf (fpHet, "MarkerList(0");
  for (k = 1; k < modelType.numMarkers; k++)
    fprintf (fpHet, ",%d", k);
  fprintf (fpHet, ")\n");
}

void writeMPBRFileDetail () {
  max = mp_result[posIdx].max_lr;
  gfreq = mp_result[posIdx].max_gfreq;
  alphaV = mp_result[posIdx].max_alpha;
  penIdx = mp_result[posIdx].max_penIdx;
  paramIdx = mp_result[posIdx].max_paramIdx;
  thresholdIdx = mp_result[posIdx].max_thresholdIdx;
  fprintf (fpHet, "%d %f %.*f %.6e %.6f %f %f",
	   (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
	   traitPos, ppl >= .025 ? 2 : 3, ppl >= .025 ? rint (ppl * 100.) / 100. : rint (ppl * 1000.) / 1000.,
	   avgLR, log10 (max), alphaV, gfreq);
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    pen_DD = modelRange.penet[liabIdx][0][penIdx];
    pen_Dd = modelRange.penet[liabIdx][1][penIdx];
    pen_dD = modelRange.penet[liabIdx][2][penIdx];
    pen_dd = modelRange.penet[liabIdx][3][penIdx];
    if (modelOptions.imprintingFlag)
      fprintf (fpHet, " (%.3f,%.3f,%.3f,%.3f", pen_DD, pen_Dd, pen_dD, pen_dd);
    else
      fprintf (fpHet, " (%.3f,%.3f,%.3f", pen_DD, pen_Dd, pen_dd);
    if (modelType.trait != DT && modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
      SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
      SD_dD = modelRange.param[liabIdx][2][0][paramIdx];
      SD_dd = modelRange.param[liabIdx][3][0][paramIdx];
      if (modelOptions.imprintingFlag)
	fprintf (fpHet, ",%.3f,%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dD, SD_dd);
      else
	fprintf (fpHet, ",%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dd);
    }
    if (modelType.trait != DT) {
      threshold = modelRange.tthresh[liabIdx][thresholdIdx];
      fprintf (fpHet, ",%.3f)", threshold);
    } else
      fprintf (fpHet, ")");
  }
  /* print out markers used for this position */
  fprintf (fpHet, " (%d", mp_result[posIdx].pMarkers[0]);
  for (k = 1; k < modelType.numMarkers; k++) {
    fprintf (fpHet, ",%d", mp_result[posIdx].pMarkers[k]);
  }
  fprintf (fpHet, ")\n");
  fflush (fpHet);
}

void writeMaximizingModel(char *modelDescription, double myMOD, int myDPrimeIdx, int myThetaIdx) {

  fprintf (fpTP, "# %s:\n", modelDescription);
  theta[0] = modelRange.theta[0][myThetaIdx];
  theta[1] = modelRange.theta[1][myThetaIdx];
  gfreq = tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].max_gfreq;
  mkrFreq = tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].max_mf;
  alphaV = tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].max_alpha;
  penIdx = tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].max_penIdx;
  R_square = tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].R_square;
  paramIdx = tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].max_paramIdx;
  thresholdIdx = tp_result[myDPrimeIdx][myThetaIdx][modelRange.nafreq].max_thresholdIdx;
  fprintf (fpTP,
	   "%d %s %.4f %.4f %.2f (%.4f,%.4f) %.3f %.2f %.4f %.4f",
	   pLocus2->pMapUnit->chromosome, pLocus2->sName,
	   pLocus2->pMapUnit->mapPos[SEX_AVERAGED], log10 (myMOD),
	   pLambdaCell->lambda[myDPrimeIdx][0][0], 
	   theta[0], theta[1], R_square, alphaV, gfreq, mkrFreq);
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    pen_DD = modelRange.penet[liabIdx][0][penIdx];
    pen_Dd = modelRange.penet[liabIdx][1][penIdx];
    pen_dD = modelRange.penet[liabIdx][2][penIdx];
    pen_dd = modelRange.penet[liabIdx][3][penIdx];
    if (modelOptions.imprintingFlag)
      fprintf (fpTP, " (%.3f,%.3f,%.3f,%.3f", pen_DD, pen_Dd, pen_dD, pen_dd);
    else
      fprintf (fpTP, " (%.3f %.3f %.3f", pen_DD, pen_Dd, pen_dd);
    if (modelType.trait != DT && modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
      SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
      SD_dD = modelRange.param[liabIdx][2][0][paramIdx];
      SD_dd = modelRange.param[liabIdx][3][0][paramIdx];
      if (modelOptions.imprintingFlag)
	fprintf (fpTP, ",%.3f,%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dD, SD_dd);
      else
	fprintf (fpTP, ",%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dd);
    }
    if (modelType.trait != DT) {
      threshold = modelRange.tthresh[liabIdx][thresholdIdx];
      fprintf (fpTP, ",%.3f)", threshold);
    } else
      fprintf (fpTP, ")");
  }
  fprintf (fpTP, "\n");
  fflush (fpTP);
}

void writeMMFileDetail() {

  fprintf (fpTP, "# %-d  %s %s\n", loc2, pLocus1->sName, pLocus2->sName);
  initialFlag = 1;
  max = -99999;
  max_at_theta0 = -99999;
  max_at_dprime0 = -99999;
  for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
    //dprime = pLambdaCell->lambda[dprimeIdx][0][0];
    for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++) {
      theta[0] = modelRange.theta[0][thetaInd];
      theta[1] = modelRange.theta[1][thetaInd];
      lr = tp_result[dprimeIdx][thetaInd][modelRange.nafreq].max_lr;
      if (initialFlag || lr > max) {
	/* overall max */
	max = lr;
	maxDPrimeIdx = dprimeIdx;
	maxThetaIdx = thetaInd;
      }
      if (initialFlag || (-ERROR_MARGIN <= theta[0] && theta[0] <= ERROR_MARGIN && -ERROR_MARGIN <= theta[1]
			  && theta[1] <= ERROR_MARGIN)) {
	/* find the max for models with theta equal to 0 */
	theta0Idx = thetaInd;
	if (lr > max_at_theta0) {
	  max_at_theta0 = lr;
	  maxDPrimeIdx_at_theta0 = dprimeIdx;
	}
      }
      if (dprimeIdx == dprime0Idx) {
	if (initialFlag || maxTheta_at_dprime0 < 0 || lr > max_at_dprime0) {
	  max_at_dprime0 = lr;
	  maxTheta_at_dprime0 = thetaInd;
	}
      }
      initialFlag = 0;
    }
    initialFlag = 0;
  }
  if (modelOptions.imprintingFlag)
    fprintf (fpTP, "Chr Marker Position MOD DPrime Theta(M,F) R2 Alpha DGF MF ");
  else
    fprintf (fpTP, "Chr Marker Position MOD DPrime Theta(M,F) R2 Alpha DGF MF ");
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
    if (modelType.trait == DT)
      if (modelOptions.imprintingFlag)
	fprintf (fpTP, "LC%dPV(DD,Dd,dD,dd) ", liabIdx);
      else
	fprintf (fpTP, "LC%dPV(DD,Dd,dd) ", liabIdx);
    else
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	if (modelOptions.imprintingFlag)
	  fprintf (fpTP, "LC%dPV(DDMean,DdMean,dDMean,ddMean,DDSD,DdSD,dDSD,ddSD,Thresh) ", liabIdx);
	else
	  fprintf (fpTP, "LC%dPV(DDMean,DdMean,ddMean,DDSD,DdSD,ddSD,Thresh) ", liabIdx);
      else
	if (modelOptions.imprintingFlag)
	  fprintf (fpTP, "LC%dPV(DDDF,DdDF,dDDF,ddDF,Thresh) ", liabIdx);
	else
	  fprintf (fpTP, "LC%dPV(DDDF,DdDF,ddDF,Thresh) ", liabIdx);
  fprintf (fpTP, "\n");
  
  /* Overall maximizing model - MOD */
  writeMaximizingModel ("Overall MOD maximizing model", max, maxDPrimeIdx, maxThetaIdx);
  
  /* Maximizing model at theta equal to 0 - MOD */
  writeMaximizingModel ("MOD maximizing model for theta=0", max_at_theta0, 
			    maxDPrimeIdx_at_theta0, theta0Idx);
  
  /* Maximizing model at d prime equal to 0 - MOD */
  writeMaximizingModel ("MOD maximizing model for dprime=0", max_at_dprime0, 
			    dprime0Idx, maxTheta_at_dprime0);
  
  /* find the overall maximizing theta and dprime - LR
   * with the other parameter integrated out */
  max = -9999.99;
  max_at_dprime0 = -9999.99;
  max_at_theta0 = -9999.99;
  for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
    for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++) {
      theta[0] = modelRange.theta[0][thetaInd];
      theta[1] = modelRange.theta[1][thetaInd];
      lr = tp_result[dprimeIdx][thetaInd][0].het_lr_avg;
      if (lr > max) {
	max = lr;
	maxThetaIdx = thetaInd;
	maxDPrimeIdx = dprimeIdx;
      }
      if (-ERROR_MARGIN <= theta[0] && theta[0] <= ERROR_MARGIN &&
	  -ERROR_MARGIN <= theta[1] && theta[1] <= ERROR_MARGIN) {
	if (lr > max_at_theta0) {
	  max_at_theta0 = lr;
	  maxDPrimeIdx_at_theta0 = dprimeIdx;
	}
      }
      if (dprime0Idx == dprimeIdx) {
	if (lr > max_at_dprime0) {
	  max_at_dprime0 = lr;
	  maxTheta_at_dprime0 = thetaInd;
	}
      }
    }
  }
  /* Overall maximizing model - LR */
  fprintf (fpTP, "# Overall LR maximizing model:\n");
  theta[0] = modelRange.theta[0][maxThetaIdx];
  theta[1] = modelRange.theta[1][maxThetaIdx];
  fprintf (fpTP,
	   "%d %s %.4f %.4f %.2f (%.4f,%.4f)\n",
	   pLocus2->pMapUnit->chromosome, pLocus2->sName,
	   pLocus2->pMapUnit->mapPos[SEX_AVERAGED], log10 (max),
	   pLambdaCell->lambda[maxDPrimeIdx][0][0], theta[0], theta[1]);
  
  /* Maximizing model at theta equal to 0 - LR */
  fprintf (fpTP, "# LR maximizing model for theta (0, 0):\n");
  fprintf (fpTP,
	   "%d %s %.4f %.4f %.2f (%.4f,%.4f)\n",
	   pLocus2->pMapUnit->chromosome, pLocus2->sName,
	   pLocus2->pMapUnit->mapPos[SEX_AVERAGED],
	   log10 (max_at_theta0), pLambdaCell->lambda[maxDPrimeIdx_at_theta0][0][0], 0.0, 0.0);
  
  /* Maximizing model at d prime equal to 0 - LR */
  fprintf (fpTP, "# LR maximizing model for dprime=0:\n");
  fprintf (fpTP,
	   "%d %s %.4f %.4f %.2f (%.4f,%.4f)\n",
	   pLocus2->pMapUnit->chromosome, pLocus2->sName,
	   pLocus2->pMapUnit->mapPos[SEX_AVERAGED], log10 (max_at_dprime0), 0.0, 0.0, 0.0);
}
