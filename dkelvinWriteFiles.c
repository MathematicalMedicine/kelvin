/* dkelvinWriteFiles.c */

void dk_write2ptBRHeader ()
{
  int i, j;

  fprintf (fpHet, "# Seq: %d Chr: %d Trait: %s Marker: %s", loc2,
	   pLocus2->pMapUnit->chromosome, pLocus1->sName, pLocus2->sName);
  
  if ((modelOptions.mapFlag == SEX_SPECIFIC) &&
      (pLocus2->pMapUnit->mapPos[MAP_FEMALE] >= 0) &&
      (pLocus2->pMapUnit->mapPos[MAP_MALE] >= 0)) {
    fprintf (fpHet, " AvgPosition: %.4f FemalePosition: %.4f MalePosition: %.4f",
	     pLocus2->pMapUnit->mapPos[MAP_SEX_AVERAGE],
	     pLocus2->pMapUnit->mapPos[MAP_FEMALE], pLocus2->pMapUnit->mapPos[MAP_MALE]);
    
  } else {
    fprintf (fpHet, " Position: %.4f", pLocus2->pMapUnit->mapPos[MAP_SEX_AVERAGE]);
  }
  
  if (pLocus2->pMapUnit->basePairLocation >= 0)
    fprintf (fpHet, " Phyiscal %d", pLocus2->pMapUnit->basePairLocation);
  fprintf (fpHet, "\n");
  
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
    for (i = 0; i < pLocus1->numOriginalAllele - 1; i++)
      for (j = 0; j < pLocus2->numOriginalAllele - 1; j++)
	fprintf (fpHet, "D%1d%1d ", i + 1, j + 1);
  fprintf (fpHet, "Theta(M,F) BayesRatio\n");    
  
  return;
}


void dk_write2ptBRData (double integral)
{
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    /* This bit is almost certainly wrong, given dkelvin's single D' limitation. */
    // int ii,jj;
    // for (ii = 0; ii < pLocus1->numOriginalAllele - 1; ii++)
    //   for (jj = 0; jj < pLocus2->numOriginalAllele - 1; jj++)
    //     fprintf (fpHet, "%.2f ", pLambdaCell->lambda[dprimeIdx][ii][jj]);
    fprintf (fpHet, "%.2f ", fixed_dprime);
  }
	      
  if (modelOptions.mapFlag == SA)
    fprintf (fpHet, "(%.4f,%.4f) %.6e\n", fixed_theta, fixed_theta, integral);
  else
    fprintf (fpHet, "(%.4f,%.4f) %.6e\n", fixed_thetaM, fixed_thetaF, integral);
  fflush (fpHet);
}


void dk_writeMPBRHeader ()
{
  int i;

  fprintf (fpHet, "Chr Position PPL BayesRatio");
  fprintf (fpHet, " MarkerList(0");
  for (i = 1; i < modelType.numMarkers; i++)
    fprintf (fpHet, ",%d", i);
  fprintf (fpHet, ")\n");
  fflush (fpHet);
}


void dk_writeMPBRData (double ppl, double br)
{
  int i;

  fprintf (fpHet, "%d %f %.*f %.6e",
	   (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
	   traitPos, ppl >= .025 ? 2 : 3, ppl >= .025 ? rint (ppl * 100.) / 100. : rint (ppl * 1000.) / 1000., br);
  /* print out markers used for this position */
  fprintf (fpHet, " (%d", mp_result[posIdx].pMarkers[0]);
  for (i = 1; i < modelType.numMarkers; i++) {
    fprintf (fpHet, ",%d", mp_result[posIdx].pMarkers[i]);
  }
  fprintf (fpHet, ")\n");
  fflush (fpHet);
}


void dk_writeMPMODHeader ()
{
  if (fpMOD == NULL)
    return;

  fprintf (fpMOD, "Chr Position MOD Alpha DGF");
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
    if (modelType.trait == DT)
      if (modelOptions.imprintingFlag)
	fprintf (fpMOD, " LC%dPV(DD,Dd,dD, dd)\n", liabIdx);
      else
	fprintf (fpMOD, " LC%dPV(DD,Dd,dd)\n", liabIdx);
    else {
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	if (modelOptions.imprintingFlag)
	  fprintf (fpMOD, " LC%dPV(DDMean,DdMean,dDMean,ddMean,DDSD,DdSD,dDSD,ddSD", liabIdx);
	else
	  fprintf (fpMOD, " LC%dPV(DDMean,DdMean,ddMean,DDSD,DdSD,ddSD", liabIdx);
      else
	if (modelOptions.imprintingFlag)
	  fprintf (fpMOD, " LC%dPV(DDDF,DdDF,dDF,ddDF", liabIdx);
	else
	  fprintf (fpMOD, " LC%dPV(DDDF,DdDF,ddDF", liabIdx);
      if (modelType.trait == CT)
	fprintf (fpMOD, ",Thesh)\n");
      else
	fprintf (fpMOD, ")\n");
    }
  fflush (fpMOD);
}


void dk_writeMPMODData (double value, st_DKMaxModel *model)
{
  if (fpMOD == NULL)
    return;

  fprintf (fpMOD, "%d %f %.6f %f %f",
	   (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->pMapUnit->chromosome,
	   traitPos, log10 (value), model->alpha, model->dgf);
  
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    if (! modelOptions.imprintingFlag)
      /* DD Dd dd or DDMean DdMean ddMean */
      fprintf (fpMOD, " (%.3f,%.3f,%.3f", model->pen[liabIdx].DD, model->pen[liabIdx].Dd,
	       model->pen[liabIdx].dd);
    else 
      /* DD Dd dD dd or DDMean DdMean dDMean ddMean */
      fprintf (fpMOD, " (%.3f,%.3f,%.3f,%.3f", model->pen[liabIdx].DD,
	       model->pen[liabIdx].Dd, model->pen[liabIdx].dD, model->pen[liabIdx].dd);
    
    if (modelType.trait != DICHOTOMOUS) {
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	if (! modelOptions.imprintingFlag)
	  /* DDSD DdSD ddSD */
	  fprintf (fpMOD, ",%.3f,%.3f,%.3f", model->pen[liabIdx].DDSD,
		   model->pen[liabIdx].DdSD, model->pen[liabIdx].ddSD);
	else 
	  /* DDSD DdSD dDSD ddSD */
	  fprintf (fpMOD, ",%.3f,%.3f,%.3f,%.3f", model->pen[liabIdx].DDSD,
		   model->pen[liabIdx].DdSD, model->pen[liabIdx].dDSD, model->pen[liabIdx].ddSD);
      }
      if (modelType.trait == CT)
	/* Theshold */
	fprintf (fpMOD, ",%.3f", model->pen[liabIdx].threshold);
    }
  }
  fprintf (fpMOD, ")\n");
  fflush (fpMOD);
}



void dk_write2ptMODHeader ()
{
  int i, j, liabIdx;

  if (fpMOD == NULL)
    return;

  if (modelOptions.markerAnalysis == FALSE) {
    fprintf (fpMOD, "# Seq: %d Chr: %d Trait: %s Marker: %s", loc2,
	     pLocus2->pMapUnit->chromosome, pLocus1->sName, pLocus2->sName);
    
    if ((modelOptions.mapFlag == SEX_SPECIFIC) &&
	(pLocus2->pMapUnit->mapPos[MAP_FEMALE] >= 0) &&
	(pLocus2->pMapUnit->mapPos[MAP_MALE] >= 0)) {
      fprintf (fpMOD, " AvgPosition: %.4f FemalePosition: %.4f MalePosition: %.4f",
	       pLocus2->pMapUnit->mapPos[MAP_SEX_AVERAGE],
	       pLocus2->pMapUnit->mapPos[MAP_FEMALE], pLocus2->pMapUnit->mapPos[MAP_MALE]);
      
    } else {
      fprintf (fpMOD, " Position: %.4f", pLocus2->pMapUnit->mapPos[MAP_SEX_AVERAGE]);
    }
    
  } else {
    fprintf (fpMOD, "# Seq: %d Chr %d: Marker1: %s Position1: %.4f Marker2: %s Position2: %.4f",
	     loc2, pLocus2->pMapUnit->chromosome,
	      pLocus1->sName, pLocus1->pMapUnit->mapPos[MAP_SEX_AVERAGE],
	      pLocus2->sName, pLocus2->pMapUnit->mapPos[MAP_SEX_AVERAGE]); 
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

  if (modelOptions.markerAnalysis != FALSE)
    fprintf (fpMOD, " Theta(M,F) R2 Alpha DGF MF");
  else
    fprintf (fpMOD, " Theta(M,F) Alpha DGF MF");
  
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
    if (modelType.trait == DT)
      if (modelOptions.imprintingFlag)
	fprintf (fpMOD, " LC%dPV(DD,Dd,dD,dd)\n", liabIdx);
      else
	fprintf (fpMOD, " LC%dPV(DD,Dd,dd)\n", liabIdx);
    else {
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	if (modelOptions.imprintingFlag)
	  fprintf (fpMOD, " LC%dPV(DDMean,DdMean,dDMean,ddMean,DDSD,DdSD,dDSD,ddSD", liabIdx);
	else
	  fprintf (fpMOD, " LC%dPV(DDMean,DdMean,ddMean,DDSD,DdSD,ddSD", liabIdx);
      else
	if (modelOptions.imprintingFlag)
	  fprintf (fpMOD, " LC%dPV(DDDF,DdDF,dDDF,ddDF", liabIdx);
	else
	  fprintf (fpMOD, " LC%dPV(DDDF,DdDF,ddDF", liabIdx);
      if (modelType.trait == CT)
	fprintf (fpMOD, ",Thresh)\n");
      else
	fprintf (fpMOD, ")\n");
    }
}


void dk_write2ptMODData (char *description, double value, st_DKMaxModel *model)
{
  if (fpMOD == NULL)
    return;

  /* MOD */
  if (modelOptions.extraMODs)
    fprintf (fpMOD, "%s ", description);
  fprintf (fpMOD, "%.4f", log10 (value));
  
  /* D' */
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {	
    /* This bit is almost certainly wrong, given dkelvin's single D' limitation. */
    // int ii,jj;
    // for (ii = 0; ii < pLocus1->numOriginalAllele - 1; ii++)
    //   for (jj = 0; jj < pLocus2->numOriginalAllele - 1; jj++)
    //     fprintf (fpMOD, " %.2f", pLambdaCell->lambda[dprimeIdx][ii][jj]);
    fprintf (fpMOD, " %.2f", model->dprime[0]);
  }

  if (modelOptions.markerAnalysis == FALSE)
    /* Theta Alpha DGF MF */
    fprintf (fpMOD, " (%.4f,%.4f) %.2f %.4f %.4f", model->theta[0], model->theta[1],
	     model->alpha, model->dgf, model->mf);
  else
    /* Theta R2 Alpha DGF MF */
    fprintf (fpMOD, " (%.4f,%.4f) %.3f %.2f %.4f %.4f", model->theta[0], model->theta[1],
	     model->r2, model->alpha, model->dgf, model->mf);
	    
  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
    if (! modelOptions.imprintingFlag)
      /* DD Dd dd or DDMean DdMean ddMean */
      fprintf (fpMOD, " (%.3f,%.3f,%.3f", model->pen[liabIdx].DD, model->pen[liabIdx].Dd,
	       model->pen[liabIdx].dd);
    else 
      /* DD Dd dD dd or DDMean DdMean dDMean ddMean */
      fprintf (fpMOD, " (%.3f,%.3f,%.3f,%.3f", model->pen[liabIdx].DD,
	       model->pen[liabIdx].Dd, model->pen[liabIdx].dD, model->pen[liabIdx].dd);
    
    if (modelType.trait != DICHOTOMOUS) {
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	if (! modelOptions.imprintingFlag)
	  /* DDSD DdSD ddSD */
	  fprintf (fpMOD, ",%.3f,%.3f,%.3f", model->pen[liabIdx].DDSD,
		   model->pen[liabIdx].DdSD, model->pen[liabIdx].ddSD);
	else 
	  /* DDSD DdSD dDSD ddSD */
	  fprintf (fpMOD, ",%.3f,%.3f,%.3f,%.3f", model->pen[liabIdx].DDSD,
		   model->pen[liabIdx].DdSD, model->pen[liabIdx].dDSD, model->pen[liabIdx].ddSD);
      }
      if (modelType.trait == CT)
	/* Theshold */
	fprintf (fpMOD, ",%.3f", model->pen[liabIdx].threshold);
    }
    fprintf (fpMOD, ")\n");
  }
  
  return;
}


void dk_copyMaxModel (double *arr, st_DKMaxModel *max)
{
  int j, idx;

  max->dgf = arr[0];
  max->alpha = arr[1];

  j=2;
  for (idx = 0; idx < modelRange.nlclass; idx++) {
    max->pen[idx].DD = arr[j++];
    max->pen[idx].Dd = arr[j++];
    if (! modelOptions.imprintingFlag) {
      max->pen[idx].dd = arr[j++];
    } else {
      max->pen[idx].dD = arr[j++];
      max->pen[idx].dd = arr[j++];
    }
    if (modelType.trait != DICHOTOMOUS) {
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	max->pen[idx].DDSD = arr[j++];
	max->pen[idx].DdSD = arr[j++];
	if (! modelOptions.imprintingFlag) {
	  max->pen[idx].ddSD = arr[j++];
	} else {
	  max->pen[idx].dDSD = arr[j++];
	  max->pen[idx].ddSD = arr[j++];
	}
      }
      if (modelType.trait == CT) 
	max->pen[idx].threshold = arr[j++];
    }
  }
}
