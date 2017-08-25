#include "kelvin.h"
#include "kelvinGlobals.h"
#include "summary_result.h"
#include "config/model_type.h"
#include <math.h>
#include <float.h>

/* dkelvinWriteFiles.c */
void dk_write2ptBRHeader (int loc1, int loc2)
{
  int i, j;

  if (fpHet == NULL)
    return;

  fprintf (fpHet, "# Seq: %d Chr: %d Trait: %s Marker: %s", loc2,
	   pLocus2->pMapUnit->chromosome, pLocus1->sName, pLocus2->sName);
  
  if ((modelOptions->mapFlag == SEX_SPECIFIC) &&
      (pLocus2->pMapUnit->mapPos[MAP_POS_FEMALE] >= 0) &&
      (pLocus2->pMapUnit->mapPos[MAP_POS_MALE] >= 0)) {
    fprintf (fpHet, " AvgPosition: %.4f FemalePosition: %.4f MalePosition: %.4f",
	     pLocus2->pMapUnit->mapPos[MAP_POS_SEX_AVERAGE],
	     pLocus2->pMapUnit->mapPos[MAP_POS_FEMALE], pLocus2->pMapUnit->mapPos[MAP_POS_MALE]);
    
  } else {
    fprintf (fpHet, " Position: %.4f", pLocus2->pMapUnit->mapPos[MAP_POS_SEX_AVERAGE]);
  }
  
  if (pLocus2->pMapUnit->basePairLocation >= 0)
    fprintf (fpHet, " Physical: %d", pLocus2->pMapUnit->basePairLocation);
  fprintf (fpHet, "\n");
  
  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
    for (i = 0; i < pLocus1->numOriginalAllele - 1; i++)
      for (j = 0; j < pLocus2->numOriginalAllele - 1; j++)
	fprintf (fpHet, "D%1d%1d ", i + 1, j + 1);
  fprintf (fpHet, "Theta(M,F) BayesRatio\n");    
  
  return;
}


void dk_write2ptBRData (double dprimevalue, double theta1, double theta2,double br, int max_scale)
{
  int exponent;
  double base, log10BR;

  log10BR=log10(br);
  exponent=floor(log10BR);
  base=pow(10, (log10BR-exponent));

  if (fpHet == NULL)
    return;

  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
    /* This bit is almost certainly wrong, given dkelvin's single D' limitation. */
    // int ii,jj;
    // for (ii = 0; ii < pLocus1->numOriginalAllele - 1; ii++)
    //   for (jj = 0; jj < pLocus2->numOriginalAllele - 1; jj++)
    //     fprintf (fpHet, "%.2f ", pLambdaCell->lambda[dprimeIdx][ii][jj]);
    fprintf (fpHet, "%.2f ", dprimevalue);
  }

  fprintf (fpHet, "(%.4f,%.4f) %.6fe%+.2d\n", theta1, theta2, base, exponent+max_scale);

  fflush (fpHet);
}


void dk_writeMPBRHeader ()
{
  int i;

  if (fpHet == NULL)
    return;

  fprintf (fpHet, "Chr Position");
  if (modelOptions->physicalMap)
    fprintf (fpHet, " Physical");
  fprintf (fpHet, " PPL BayesRatio MarkerList(0");
  for (i = 1; i < modelType->numMarkers; i++)
    fprintf (fpHet, ",%d", i);
  fprintf (fpHet, ")\n");
  fflush (fpHet);
}


void dk_writeMPBRData (int posIdx, float traitPos, double ppl, double br, int max_scale)
{
  int i;
  int exponent;
  double base, log10BR;

  log10BR=log10(br);
  exponent=floor(log10BR);
  base=pow(10, (log10BR-exponent));

  if (fpHet == NULL)
    return;

  fprintf (fpHet, "%d %f",
	   (originalLocusList.ppLocusList[1])->pMapUnit->chromosome, traitPos);
  if (modelOptions->physicalMap)
    fprintf (fpHet, " %d", interpolate_physical_location (traitPos));
  fprintf (fpHet, " %.*f %.6fe%+.2d", ppl >= .025 ? 2 : 3, KROUND (ppl, 3), base,
	   exponent+max_scale);
  /* print out markers used for this position */
  fprintf (fpHet, " (%d", mp_result[posIdx].pMarkers[0]);
  for (i = 1; i < modelType->numMarkers; i++) {
    fprintf (fpHet, ",%d", mp_result[posIdx].pMarkers[i]);
  }
  fprintf (fpHet, ")\n");
  fflush (fpHet);
}


void dk_writeMPMODHeader ()
{
  int liabIdx;

  if (fpMOD == NULL)
    return;

  fprintf (fpMOD, "Chr Position");
  if (modelOptions->physicalMap) 
    fprintf (fpMOD, " Physical");
  fprintf (fpMOD, " MOD Alpha DGF");

  for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++)
    if (modelType->trait == DT)
      if (modelOptions->imprintingFlag)
	fprintf (fpMOD, " LC%dPV(DD,Dd,dD,dd)", modelRange->lclassLabels[liabIdx]);
      else
	fprintf (fpMOD, " LC%dPV(DD,Dd,dd)", modelRange->lclassLabels[liabIdx]);
    else {
      if (modelType->distrib != QT_FUNCTION_CHI_SQUARE)
	if (modelOptions->imprintingFlag)
	  fprintf (fpMOD, " LC%dPV(DDMean,DdMean,dDMean,ddMean,DDSD,DdSD,dDSD,ddSD",
		   modelRange->lclassLabels[liabIdx]);
	else
	  fprintf (fpMOD, " LC%dPV(DDMean,DdMean,ddMean,DDSD,DdSD,ddSD",
		   modelRange->lclassLabels[liabIdx]);
      else
	if (modelOptions->imprintingFlag)
	  fprintf (fpMOD, " LC%dPV(DDDF,DdDF,dDF,ddDF", modelRange->lclassLabels[liabIdx]);
	else
	  fprintf (fpMOD, " LC%dPV(DDDF,DdDF,ddDF", modelRange->lclassLabels[liabIdx]);
      if (modelType->trait == CT)
	fprintf (fpMOD, ",Thesh)");
      else 
	fprintf (fpMOD, ")");
    }
  fprintf (fpMOD, "\n");
  fflush (fpMOD);
}


void dk_writeMPMODData (int posIdx, float traitPos, double value, st_DKMaxModel *model)
{
  int liabIdx;
  double Pen_DD, Pen_Dd, Pen_dD, Pen_dd;
  double SD_DD, SD_Dd, SD_dD, SD_dd;
  double threshold;

  if (fpMOD == NULL)
    return;

  fprintf (fpMOD, "%d %f", (originalLocusList.ppLocusList[mp_result[posIdx].pMarkers[0]])->
	   pMapUnit->chromosome, traitPos);
  if (modelOptions->physicalMap)
    fprintf (fpMOD, " %d", interpolate_physical_location (traitPos));
  fprintf (fpMOD, " %.6f %f %f", value, model->alpha, model->dgf);
  
  for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++) {
    if (! modelRange->lclassEmpty[liabIdx]) {
      Pen_DD = deNormalizeMean (modelType, model->pen[liabIdx].DD);
      Pen_Dd = deNormalizeMean (modelType, model->pen[liabIdx].Dd);
      Pen_dD = deNormalizeMean (modelType, model->pen[liabIdx].dD);
      Pen_dd = deNormalizeMean (modelType, model->pen[liabIdx].dd);
    } else 
      Pen_DD = Pen_Dd = Pen_dD = Pen_dd = -99.0;
    if (! modelOptions->imprintingFlag)
      /* DD Dd dd or DDMean DdMean ddMean */
      fprintf (fpMOD, " (%.3f,%.3f,%.3f", Pen_DD, Pen_Dd, Pen_dd);
    else 
      /* DD Dd dD dd or DDMean DdMean dDMean ddMean */
      fprintf (fpMOD, " (%.3f,%.3f,%.3f,%.3f", Pen_DD, Pen_Dd, Pen_dD, Pen_dd);
    
    if (modelType->trait != DICHOTOMOUS) {
      if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
	if (! modelRange->lclassEmpty[liabIdx]) {
	  SD_DD = deNormalizeStdev (modelType, model->pen[liabIdx].DDSD);
	  SD_Dd = deNormalizeStdev (modelType, model->pen[liabIdx].DdSD);
	  SD_dD = deNormalizeStdev (modelType, model->pen[liabIdx].dDSD);
	  SD_dd = deNormalizeStdev (modelType, model->pen[liabIdx].ddSD);
	} else 
	  SD_DD = SD_Dd = SD_dD = SD_dd = -99.0;
	if (! modelOptions->imprintingFlag)
	  /* DDSD DdSD ddSD */
	  fprintf (fpMOD, ",%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dd);
	else 
	  /* DDSD DdSD dDSD ddSD */
	  fprintf (fpMOD, ",%.3f,%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dD, SD_dd);
      }
      if (modelType->trait == CT) {
	threshold = (modelOptions->imprintingFlag) ? -99.0 :
	  deNormalizeMean (modelType, model->pen[liabIdx].threshold);
	/* Theshold */
	fprintf (fpMOD, ",%.3f", threshold);
      }
    }
    fprintf (fpMOD, ")");
  }

  fprintf (fpMOD, "\n");
  fflush (fpMOD);
}



void dk_write2ptMODHeader ()
{
  int i, j, liabIdx;

  if (fpMOD == NULL)
    return;

  if (modelOptions->markerAnalysis == FALSE) {
    fprintf (fpMOD, "# Seq: %d Chr: %d Trait: %s Marker: %s", loc2,
	     pLocus2->pMapUnit->chromosome, pLocus1->sName, pLocus2->sName);
    
    if ((modelOptions->mapFlag == SEX_SPECIFIC) &&
	(pLocus2->pMapUnit->mapPos[MAP_POS_FEMALE] >= 0) &&
	(pLocus2->pMapUnit->mapPos[MAP_POS_MALE] >= 0)) {
      fprintf (fpMOD, " AvgPosition: %.4f FemalePosition: %.4f MalePosition: %.4f",
	       pLocus2->pMapUnit->mapPos[MAP_POS_SEX_AVERAGE],
	       pLocus2->pMapUnit->mapPos[MAP_POS_FEMALE], pLocus2->pMapUnit->mapPos[MAP_POS_MALE]);
      
    } else {
      fprintf (fpMOD, " Position: %.4f", pLocus2->pMapUnit->mapPos[MAP_POS_SEX_AVERAGE]);
    }
    
  } else {
    fprintf (fpMOD, "# Seq: %d Chr %d: Marker1: %s Position1: %.4f Marker2: %s Position2: %.4f",
	     loc2, pLocus2->pMapUnit->chromosome,
	      pLocus1->sName, pLocus1->pMapUnit->mapPos[MAP_POS_SEX_AVERAGE],
	      pLocus2->sName, pLocus2->pMapUnit->mapPos[MAP_POS_SEX_AVERAGE]); 
  }
  
  if (pLocus2->pMapUnit->basePairLocation >= 0)
    fprintf (fpMOD, " Physical: %d", pLocus2->pMapUnit->basePairLocation);
  fprintf (fpMOD, "\n");

  if (modelOptions->extraMODs)
    fprintf (fpMOD, "Case ");
  fprintf (fpMOD, "MOD");
  
  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
    for (i = 0; i < pLocus1->numOriginalAllele - 1; i++)
      for (j = 0; j < pLocus2->numOriginalAllele - 1; j++)
	fprintf (fpMOD, " D%1d%1d", i + 1, j + 1);

  
  fprintf (fpMOD, " Theta(M,F)");
  if (modelOptions->markerAnalysis != FALSE) {
    if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
      fprintf (fpMOD, " R2");
  } else
    fprintf (fpMOD, " Alpha DGF");
  
  for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++)
    if (modelType->trait == DT)
      if (modelOptions->imprintingFlag)
	fprintf (fpMOD, " LC%dPV(DD,Dd,dD,dd)", modelRange->lclassLabels[liabIdx]);
      else
	fprintf (fpMOD, " LC%dPV(DD,Dd,dd)", modelRange->lclassLabels[liabIdx]);
    else {
      if (modelType->distrib != QT_FUNCTION_CHI_SQUARE)
	if (modelOptions->imprintingFlag)
	  fprintf (fpMOD, " LC%dPV(DDMean,DdMean,dDMean,ddMean,DDSD,DdSD,dDSD,ddSD",
		   modelRange->lclassLabels[liabIdx]);
	else
	  fprintf (fpMOD, " LC%dPV(DDMean,DdMean,ddMean,DDSD,DdSD,ddSD",
		   modelRange->lclassLabels[liabIdx]);
      else
	if (modelOptions->imprintingFlag)
	  fprintf (fpMOD, " LC%dPV(DDDF,DdDF,dDDF,ddDF", modelRange->lclassLabels[liabIdx]);
	else
	  fprintf (fpMOD, " LC%dPV(DDDF,DdDF,ddDF", modelRange->lclassLabels[liabIdx]);
      if (modelType->trait == CT)
	fprintf (fpMOD, ",Thresh)");
      else
	fprintf (fpMOD, ")");
    }
  fprintf (fpMOD, "\n");
}


void dk_write2ptMODData (char *description, double value, st_DKMaxModel *model)
{
  int liabIdx;
  double Pen_DD, Pen_Dd, Pen_dD, Pen_dd;
  double SD_DD, SD_Dd, SD_dD, SD_dd;
  double threshold = -99.0; // assigned to get rid of a spurious compiler warning

  if (fpMOD == NULL)
    return;

  /* MOD */
  if (modelOptions->extraMODs)
    fprintf (fpMOD, "%s ", description);
  if (value != -DBL_MAX)
    fprintf (fpMOD, "%.4f", value); ////log10 (value)); 6/4/2009
  else
    fprintf (fpMOD, "NoInf");
  
  /* D' */
  if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {	
    /* This bit is almost certainly wrong, given dkelvin's single D' limitation. */
    // int ii,jj;
    // for (ii = 0; ii < pLocus1->numOriginalAllele - 1; ii++)
    //   for (jj = 0; jj < pLocus2->numOriginalAllele - 1; jj++)
    //     fprintf (fpMOD, " %.2f", pLambdaCell->lambda[dprimeIdx][ii][jj]);
    fprintf (fpMOD, " %.2f", model->dprime[0]);
  }

  /* Theta */
  fprintf (fpMOD, " (%.4f,%.4f)", model->theta[0], model->theta[1]);
  if (modelOptions->markerAnalysis != FALSE) {
    if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM)
      /* R2 for LD marker-to-marker only */
      fprintf (fpMOD, " %.3f", model->r2);
  } else
    /* Alpha and DGF for trait-to-marker only */
    fprintf (fpMOD, " %.2f %.4f", model->alpha, model->dgf);
	    
  for (liabIdx = 0; liabIdx < modelRange->nlclass; liabIdx++) {
    if (! modelRange->lclassEmpty[liabIdx]) {
      Pen_DD = deNormalizeMean (modelType, model->pen[liabIdx].DD);
      Pen_Dd = deNormalizeMean (modelType, model->pen[liabIdx].Dd);
      Pen_dD = deNormalizeMean (modelType, model->pen[liabIdx].dD);
      Pen_dd = deNormalizeMean (modelType, model->pen[liabIdx].dd);
    } else 
      Pen_DD = Pen_Dd = Pen_dD = Pen_dd = -99.0;
    if (! modelOptions->imprintingFlag)
      /* DD Dd dd or DDMean DdMean ddMean */
      fprintf (fpMOD, " (%.3f,%.3f,%.3f", Pen_DD, Pen_Dd, Pen_dd);
    else 
      /* DD Dd dD dd or DDMean DdMean dDMean ddMean */
      fprintf (fpMOD, " (%.3f,%.3f,%.3f,%.3f", Pen_DD, Pen_Dd, Pen_dD, Pen_dd);
    
    if (modelType->trait != DICHOTOMOUS) {
      if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
	if (! modelRange->lclassEmpty[liabIdx]) {
	  SD_DD = deNormalizeStdev (modelType, model->pen[liabIdx].DDSD);
	  SD_Dd = deNormalizeStdev (modelType, model->pen[liabIdx].DdSD);
	  SD_dD = deNormalizeStdev (modelType, model->pen[liabIdx].dDSD);
	  SD_dd = deNormalizeStdev (modelType, model->pen[liabIdx].ddSD);
	} else 
	  SD_DD = SD_Dd = SD_dD = SD_dd = -99.0;
	if (! modelOptions->imprintingFlag)
	  /* DDSD DdSD ddSD */
	  fprintf (fpMOD, ",%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dd);
	else 
	  /* DDSD DdSD dDSD ddSD */
	  fprintf (fpMOD, ",%.3f,%.3f,%.3f,%.3f", SD_DD, SD_Dd, SD_dD, SD_dd);
      }
      if (modelType->trait == CT)
	threshold = (modelOptions->imprintingFlag) ? -99.0 :
	  deNormalizeMean (modelType, model->pen[liabIdx].threshold);
	/* Theshold */
	fprintf (fpMOD, ",%.3f", threshold);
    }
    fprintf (fpMOD, ")");
  }
  fprintf (fpMOD, "\n");
  return;
}


void dk_copyMaxModel (double *arr, st_DKMaxModel *max, int num)
{
  int j, idx;

  max->dgf = arr[0];
  max->alpha = arr[1];

  j=2;
  for (idx = 0; idx < modelRange->nlclass; idx++) {
    max->pen[idx].DD = arr[j++];
    max->pen[idx].Dd = arr[j++];
    if (! modelOptions->imprintingFlag) {
      max->pen[idx].dd = arr[j++];
    } else {
      max->pen[idx].dD = arr[j++];
      max->pen[idx].dd = arr[j++];
    }
    if (modelType->trait != DICHOTOMOUS) {
      if (modelType->distrib != QT_FUNCTION_CHI_SQUARE) {
	max->pen[idx].DDSD = arr[j];//arr[j++];
	max->pen[idx].DdSD = arr[j];//arr[j++];
	if (! modelOptions->imprintingFlag) {
	  max->pen[idx].ddSD = arr[j];//arr[j++];
	} else {
	  max->pen[idx].dDSD = arr[j];//arr[j++];
	  max->pen[idx].ddSD = arr[j];//arr[j++];
	}
      }
      if (modelType->trait == CT) 
	max->pen[idx].threshold = arr[num-1];  //arr[j++];
    }
  }
  /*if (modelType->trait == CT) 
    max->pen[idx].threshold = arr[j++];*/
}
