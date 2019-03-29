/* This file contains functions to read in marker map file and 
 * marker file 
 * 
 * Marker map file is in Marshfield format. For example, for chromosome 1
 * http://research.marshfieldclinic.org/genetics/Map_Markers/data/maps/map1.txt
 *   Marker      Dnumber      sex-ave(cM)      female(cM)        male(cM)
 * 1 AFM214yg7   D1S243              0.00             0.00             0.00
 *                           4.22             4.46             3.54
 * 2 AFM280we5   D1S468              4.22             4.46             3.54
 *                           4.63             2.94             6.67
 * For each marker except the last one, there are two lines:
 * MarkerNo  MarkerName DNumber SexAveragePos FemalePos MalePos
 * 		DistanceBetweenThisMarker and NextMarker (Avg Female Male)
 * Distance information will be ignored by this program
 * Occasionally the markerNo is lead by "*" or "x" or other special indicators
 * User should remove them before running this program. This is to force
 * users to think about the special meanings before running the program
 * 
 * Marker file follows the format below:
 * Total#ofMarkers
 * 3 <#ofAlleles> #MarkerName
 * <allele1> <allele2> ....
 * 3 <#ofAlleles> #MarkerName
 * <allele1> <allele2> ...
 * 
 * Using 3 is to be in sync of the linkage format marker spec section of the
 * datafile.dat
 * Marker Name has to match the name in map file
 * Markers should be listed based on their map order
 * Not all markers in the map file are required in the marker file
 *
 * So map file usually is a very stable file 
 * pedfile marker phenotype spec order needs to match the marker file, i.e.
 * based on map order
 *
 * Map file will be read in first, as it potentially has a bigger set of
 * markers
 * 
 * */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "gsl/gsl_sf_gamma.h"
#include "cdflib.h"

#include "pedlib.h"
#include "utils.h"		/* for logging */
#include "tools.h"
#include "allele_set_recoding.h"
#ifndef NO_POLYNOMIAL
#include "polynomial.h"
#endif
/* global variables */
Map map;
LocusList originalLocusList;
SubLocusList locusList;

/* internal functions */
MapUnit *add_map_unit (Map * map);
Locus *add_locus (LocusList * pLocusList, char *sName, int locusType,
		  int numAllele);

Genotype *add_genotype (Genotype ** ppList, int *pCount, int locusIndex,
			int allele1, int allele2);

Trait *add_trait (int trait, TraitLocus * pTraitLocus, int traitType);
void compute_penetrance (Person * pPerson, int locus,
			 int allele1, int allele2, void *pen);
void setup_LD_haplotype_freq (LDLoci * pLDLoci);

/* read in mapfile 
 * all marker information will be saved in the super marker list 
 * all markers should be on the same chromosome */
int
read_mapfile (char *sMapfileName)
{
  FILE *fp = NULL;
  int lineNo = 0;
  char line[MAX_LINE_LEN];
  char sTemp[MAX_LINE_LEN];
  int numRet;
  MapUnit *pMarker;
  int lineFlag = 1;

  fp = fopen (sMapfileName, "r");
  KASSERT (fp != NULL, "Can't open map file %s for read. Exiting!\n",
	   sMapfileName);
  /* skip the header line - assuming the leading word is "Marker" */
  while (fgets (line, MAX_LINE_LEN, fp))
    {
      lineNo++;
      numRet = sscanf (line, "%s", sTemp);
      if (numRet == 1 && !strcasecmp (sTemp, "Marker"))
	{
	  /* we can move on to the actual marker section */
	  break;
	}
    }

  KASSERT (feof (fp) == 0,
	   "No marker information at all in the map file %s (Total # of Lines=%d). Exiting!\n",
	   sMapfileName, lineNo);

  /* read a pair of lines for each marker, except the last marker */
  lineFlag = 1;
  while (fgets (line, MAX_LINE_LEN, fp))
    {
      lineNo++;
      if (is_line_blank (line))
	continue;
      switch (lineFlag)
	{
	case 1:
	  /* reading the marker line */
	  pMarker = add_map_unit (&map);
	  numRet = sscanf (line, "%d %s %s %lf %lf %lf",
			   &pMarker->no, pMarker->sName, pMarker->sDNumber,
			   &pMarker->mapPos[2],
			   &pMarker->mapPos[MOM], &pMarker->mapPos[DAD]);
	  KASSERT (numRet == 6,
		   "Can't get marker information from line %d in %s.",
		   lineNo, sMapfileName);
	  KASSERT (pMarker->no > 0,
		   "Marker number is not greater than 0 (%d) in file %s(%d).\n",
		   pMarker->no, sMapfileName, lineNo);
	  lineFlag = 2;
	  break;
	case 2:
	  /* ignore this line of inter-marker distance */
	  lineFlag = 1;
	  break;
	  /* lineFlag can not possibily be any other value, so no default is set */
	}

    }

  fclose (fp);
  return 0;
}

/* add a marker to the map */
MapUnit *
add_map_unit (Map * pMap)
{
  MapUnit *pMapUnit;

  if (pMap->maxUnit <= pMap->count)
    {
      /* need to allocate a chunk of memories now */
      pMap->maxUnit += DEF_LOCUS_MALLOC_INCREMENT;
      pMap->ppMapUnitList = (MapUnit **) REALLOC ("pMap->ppMapUnitList",
						  pMap->ppMapUnitList,
						  sizeof (MapUnit *) *
						  pMap->maxUnit);
    }

  pMapUnit = (MapUnit *) MALLOC ("pMapUnit", sizeof (MapUnit));
  pMap->ppMapUnitList[pMap->count] = pMapUnit;
  pMapUnit->mapIndex = pMap->count;

  /* increment the counter */
  pMap->count++;
  return pMapUnit;
}

/* locate marker map information by its name */
MapUnit *
find_map_unit (Map * pMap, char *sName)
{
  int i;

  for (i = 0; i < pMap->count; i++)
    {
      if (!strcasecmp (pMap->ppMapUnitList[i]->sName, sName))
	{
	  /* found it, return the pointer */
	  return pMap->ppMapUnitList[i];
	}
    }

  /* exhausted the list, failed to find it */
  return (NULL);
}

/* each marker Locus type has MapUnit pointer saved
 * each MapUnit has marker position information 
 * sexFlag: 3 DAD, 1-Mom, -1: average 
 * mapFunctionFlag: 0 - kosambi 1-Haldane
 * Return: distance in recombination fraction
 */
double
get_inter_marker_distance (MapUnit * pMarker1, MapUnit * pMarker2,
			   int sexFlag, int mapFunctionFlag)
{
  double pos1, pos2, distance;
  double theta;
  double temp;

  if (pMarker1 == NULL || pMarker2 == NULL || pMarker1 == pMarker2)
    return 0;
  pos1 = pMarker1->mapPos[sexFlag];
  pos2 = pMarker2->mapPos[sexFlag];
  if (pos2 > pos1)
    distance = pos2 - pos1;
  else
    distance = pos1 - pos2;
  /* now convert the map distance to recombination fraction */
  if (mapFunctionFlag == MAP_FUNCTION_KOSAMBI)
    {
      temp = exp (4 * distance);
      theta = 0.5 * (temp - 1) / (temp + 1);
    }
  else if (mapFunctionFlag == MAP_FUNCTION_HALDANE)
    {
      theta = 0.5 * (1 - exp (-2 * distance));
    }
  else
    /* something is wrong */
    theta = -1;

  return theta;
}

/* read datafile with specifications for each locus 
 * the format of the datafile is very similar to linkage format datafile.dat
 * */
int
read_datafile (char *sDatafileName)
{
  FILE *fp;
  int lineNo = 0;
  char line[MAX_LINE_LEN];
  int numRet;
  MapUnit *pMapUnit;
  int lineFlag = 1;
  int pos = 0;
  char *pLine;
  Locus *pLocus;
  Trait *pTrait = NULL;
  TraitLocus *pTraitLocus = NULL;
  int traitIdx;
  int locusType;		/* temporary place holder */
  int i, j;
  int numAllele;
  char sLocusName[MAX_LINE_LEN];
  int numGeno = 0;		/* possible genotypes */
  int allele1, allele2;
  double penetranceValue;
  int liabilityIndex;
  int numLocus;
  char sLocusType[MAX_LINE_LEN];
  char sTraitType[MAX_LINE_LEN];
  int numTrait;
  int traitType;
  double value;
  int trait1, trait2;
  LDLoci *pLDLoci;
  int numVariance;
  char sAnalysisType[MAX_LINE_LEN];
  char sPolyFlag[MAX_LINE_LEN];
  int locus1, locus2;
  int numAllele1, numAllele2;
  Locus *pLocus1, *pLocus2;
  int numLiab;


  fp = fopen (sDatafileName, "r");
  KASSERT (fp != NULL, "Can't open datafile %s for read. Exiting!\n",
	   sDatafileName);

  /* read analysis type  */
  while (fgets (line, MAX_LINE_LEN, fp))
    {
      lineNo++;
      /* skip the comment or blank lines */
      if (is_line_blank_or_comment (line))
	continue;
      numRet = sscanf (line, "%s %s", sAnalysisType, sPolyFlag);
      KASSERT (numRet >= 1,
	       "Can't get the analysis type (LD or LE) in file %s.\n",
	       sDatafileName);
      if (!strcasecmp (sAnalysisType, "LD"))
	modelOptions.equilibrium = LINKAGE_DISEQUILIBRIUM;
      else
	modelOptions.equilibrium = LINKAGE_EQUILIBRIUM;
      if (numRet > 1 && !strcasecmp (sPolyFlag, "NOPOLY"))
	modelOptions.polynomial = FALSE;
      else
	modelOptions.polynomial = TRUE;
      break;
    }
  /* read the total number of loci line */
  while (fgets (line, MAX_LINE_LEN, fp))
    {
      lineNo++;
      /* skip the comment or blank lines */
      if (is_line_blank_or_comment (line))
	continue;
      numRet = sscanf (line, "%d", &numLocus);
      KASSERT (numRet == 1,
	       "Can't get the total number of loci in file %s.\n",
	       sDatafileName);
      break;
    }
  KASSERT (feof (fp) == 0,
	   "Reach end of file %s(%d) before getting any locus information.\n",
	   sDatafileName, lineNo);

  /* lineFlag indicates what information we expect for the next line 
   * 1 - line of locus type, number of alleles, locus name */
  lineFlag = 1;
  while (get_nonblank_line (line, MAX_LINE_LEN, fp, &lineNo))
    {
      /* expect to read "<locusType> <#ofAlleles> #LocusName */
      numRet =
	sscanf (line, "%s %d # %s", sLocusType, &numAllele, sLocusName);
      if (numRet == 1 && !strcasecmp (sLocusType, "LD"))
	{
	  /* we are into LD parameter section 
	   * this should follow all trait and marker specification */
	  KASSERT (get_nonblank_line (line, MAX_LINE_LEN, fp, &lineNo) !=
		   NULL, "Can't get LD parameter values on line %s(%d).\n",
		   sDatafileName, lineNo);
	  numRet =
	    sscanf (line, "Loci=%d,%d DPrime=%n", &locus1, &locus2, &pos);
	  KASSERT (numRet == 2,
		   "Can't get two loci for the LD parameters.\n");
	  /* find the two loci first to determine how many D-primes 
	   * we expect */
	  pLocus1 = originalLocusList.ppLocusList[locus1 - 1];
	  pLocus1->LDLocus = locus2 - 1;
	  pLocus2 = originalLocusList.ppLocusList[locus2 - 1];
	  pLocus2->LDLocus = locus1 - 1;
	  numAllele1 = pLocus1->numAllele;
	  numAllele2 = pLocus2->numAllele;
	  if (originalLocusList.numLDLoci <= originalLocusList.maxNumLDLoci)
	    {
	      originalLocusList.maxNumLDLoci += DEF_LOCUS_MALLOC_INCREMENT;
	      originalLocusList.pLDLoci =
		(LDLoci *) REALLOC ("originalLocusList.pLDLoci",
				    originalLocusList.pLDLoci,
				    sizeof (LDLoci) *
				    originalLocusList.maxNumLDLoci);
	    }
	  pLDLoci = &originalLocusList.pLDLoci[originalLocusList.numLDLoci];
	  originalLocusList.numLDLoci += 1;
	  pLDLoci->locus1 = locus1 - 1;
	  pLDLoci->locus2 = locus2 - 1;
	  pLDLoci->ppDPrime = (double **) MALLOC ("pLDLoci->ppDPrime",
						  sizeof (double *) *
						  (numAllele1 - 1));
	  pLDLoci->ppDValue = (double **) MALLOC ("pLDLoci->ppDValue",
						  sizeof (double *) *
						  (numAllele1 - 1));
	  for (i = 0; i < numAllele1 - 1; i++)
	    {
	      pLDLoci->ppDPrime[i] =
		(double *) MALLOC ("pLDLoci->ppDPrime[i]",
				   sizeof (double) * (numAllele2 - 1));
	      pLDLoci->ppDValue[i] =
		(double *) MALLOC ("pLDLoci->ppDValue[i]",
				   sizeof (double) * (numAllele2 - 1));
	    }
	  pLDLoci->ppHaploFreq = (double **) MALLOC ("pLDLoci->ppHaploFreq",
						     sizeof (double *) *
						     numAllele1);
	  for (i = 0; i < numAllele1; i++)
	    {
	      pLDLoci->ppHaploFreq[i] =
		(double *) MALLOC ("pLDLoci->ppHaploFreq[i]",
				   sizeof (double) * numAllele2);
	    }


	  /* read D-prime one by one */
	  pLine = &line[pos];
	  for (i = 0; i < numAllele1 - 1; i++)
	    {
	      for (j = 0; j < numAllele2 - 1; j++)
		{
		  numRet =
		    sscanf (pLine, "%lf %n", &pLDLoci->ppDPrime[i][j], &pos);
		  KASSERT (numRet == 1,
			   "Can't read D Prime for allele pair (%d, %d) between locus (%d, %d)in file %s.\n",
			   i + 1, j + 1, locus1, locus2, sDatafileName,
			   lineNo);
		  pLine = &pLine[pos];
		}
	    }
	  /* set up D points */
	  setup_LD_haplotype_freq (pLDLoci);
	  continue;
	}

      KASSERT (numRet == 3,
	       "Can't get total number of alleles or locus name.\n");
      if (!strcasecmp (sLocusType, "marker"))
	{
	  locusType = LOCUS_TYPE_MARKER;
	  /* marker locus */
	  /* allocate space for the marker info first */
	  pLocus = add_locus (&originalLocusList, sLocusName,
			      locusType, numAllele);
	  /* locate this marker in the map */
	  pMapUnit = find_map_unit (&map, sLocusName);
	  KASSERT (pMapUnit != NULL,
		   "Can't find marker %s in map.\n", pLocus->sName);
	  pLocus->pMapUnit = pMapUnit;
	}
      else
	{
	  locusType = LOCUS_TYPE_TRAIT;
	  pLocus = add_locus (&originalLocusList, sLocusName,
			      locusType, numAllele);
	  pTraitLocus = pLocus->pTraitLocus;
	  numGeno = numAllele * (numAllele + 1) / 2;
	}

      /* now read allele frequencies */
      KASSERT (get_nonblank_line (line, MAX_LINE_LEN, fp, &lineNo) != NULL,
	       "Can't get allele frequencies for marker %s on line %s(%d).\n",
	       sLocusName, sDatafileName, lineNo);
      /* read each allele frequency one by one */
      pLine = line;
      for (i = 0; i < pLocus->numAllele; i++)
	{
	  numRet =
	    sscanf (pLine, "%lf %n", &pLocus->pAlleleFrequency[i], &pos);
	  KASSERT (numRet == 1,
		   "Can't read allele frequency for locus %s allele %d in file %s.\n",
		   pLocus->sName, i + 1, sDatafileName, lineNo);
	  pLine = &pLine[pos];
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      /* build polynomial for the allele frequencies */
	    }
#endif
	}

      construct_original_allele_set_list (originalLocusList.numLocus - 1);
      if (pLocus->locusType == LOCUS_TYPE_MARKER)
	{
	  //  construct_original_allele_set_list(originalLocusList.numLocus-1);
	  /* if we are dealing with marker, we are done with this locus */
	  continue;
	}

      /* need to find out number of trait variables we are dealing with 
       * Sometimes there are more than 1 trait variables associated 
       * with trait 
       * each trait varialble can be of affection status type or 
       * of quantitative trait type 
       * */
      KASSERT (get_nonblank_line (line, MAX_LINE_LEN, fp, &lineNo) != NULL,
	       "Can't get number of trait variable for locus %s on line %s(%d).\n",
	       sLocusName, sDatafileName, lineNo);
      /* expect to read number of trait variables */
      numRet = sscanf (line, "%d", &numTrait);
      KASSERT (numRet == 1,
	       "Can't get number of trait varaible in file %s(%d).\n",
	       sDatafileName, lineNo);
      /* although we could have only marker loci, no trait locus,
       * but if we do have a trait locus, then number of trait variables has 
       * to be greater than 0 */
      KASSERT (numTrait > 0,
	       "Number of trait variable in file %s(%d) is less than 1.\n",
	       sDatafileName, lineNo);
      pTraitLocus->numTrait = numTrait;

      /* allocate space for each trait variable */
      // pTraitLocus->pTraits = (Trait *)malloc(sizeof(Trait ) * numTrait);
      /* keep track of how many trait variables we have processed */
      traitIdx = 0;
      while (traitIdx < pTraitLocus->numTrait)
	{
	  /* read each trait variable */
	  KASSERT (get_nonblank_line (line, MAX_LINE_LEN, fp, &lineNo) !=
		   NULL,
		   "Can't get trait variable %d information for locus %d on line %s(%d).\n",
		   traitIdx + 1, sLocusName, sDatafileName, lineNo);
	  /* expecting <traitVariableType> [<NumberOfLiabilityClass>] */
	  numRet = sscanf (line, "%s %d %n", sTraitType, &numLiab, &pos);
	  KASSERT (numRet == 2,
		   "Can't get trait variable type or number of liablity class %s(%d).\n",
		   sDatafileName, lineNo);
	  KASSERT (numLiab > 0,
		   "Number of liability class in file %s(%d) is less than 1.\n",
		   sDatafileName, lineNo);
	  if (!strcasecmp (sTraitType, "Affect"))
	    {
	      /* dichotomouse trait - affection status */
	      traitType = DICHOTOMOUS;
	    }
	  else if (!strcasecmp (sTraitType, "QT"))
	    {
	      /* quantitative trait */
	      traitType = QUANTITATIVE;
	    }
	  else
	    {
	      traitType = COMBINED;
	    }
	  pTrait = add_trait (traitIdx, pTraitLocus, traitType);
	  pTrait->numLiabilityClass = numLiab;
	  pLine = &line[pos];
	  if (traitType == QUANTITATIVE || traitType == COMBINED)
	    {
	      pTrait->functionQT = 1;
	      pTrait->dfQT = 30;
	      pTrait->unknownTraitValue = -99.99;
	      numRet = sscanf (pLine, "%d %d %lf %lf %lf %lf",
			       &pTrait->functionQT,
			       &pTrait->dfQT,
			       &pTrait->unknownTraitValue,
			       &pTrait->lessCutoffFlag,
			       &pTrait->moreCutoffFlag, &pTrait->cutoffValue);
	      pTrait->functionQT = 1;
	      pTrait->dfQT = 30;
	      pTrait->unknownTraitValue = -99.99;
	    }
	  if (traitType == DICHOTOMOUS)
	    {
	      /* for affection status, continue to read penetrances 
	       * one penetrance line for each liability class - usually only one 
	       * number of penetrances on each line should be numGeno
	       * in most of cases, there are only two alleles, so 3 genotypes */
	      liabilityIndex = 0;
	      while (liabilityIndex < pTrait->numLiabilityClass)
		{
		  KASSERT (get_nonblank_line (line, MAX_LINE_LEN, fp, &lineNo)
			   != NULL,
			   "Can't get penetrance for liability %d on trait %s on line %s(%d).\n",
			   liabilityIndex + 1, sLocusName, sDatafileName,
			   lineNo);
		  allele1 = 1;
		  allele2 = 1;
		  i = 0;
		  pLine = line;
		  /* genotype order is: 11 12 13 ... 22 23... 33 ... */
		  while (i < numGeno)
		    {
		      numRet =
			sscanf (pLine, "%lf %n", &penetranceValue, &pos);
		      pTrait->
			penetrance[AFFECTION_STATUS_AFFECTED][liabilityIndex]
			[allele1 - 1][allele2 - 1] = penetranceValue;
		      /* do we really need to store the following two ??? */
		      pTrait->penetrance[AFFECTION_STATUS_UNAFFECTED]
			[liabilityIndex][allele1 - 1][allele2 - 1] =
			1 - penetranceValue;
		      pTrait->
			penetrance[AFFECTION_STATUS_UNKNOWN][liabilityIndex]
			[allele1 - 1][allele2 - 1] = 1;
		      /* the penetrance matrix is symmetric */
		      pTrait->
			penetrance[AFFECTION_STATUS_AFFECTED][liabilityIndex]
			[allele2 - 1][allele1 - 1] = penetranceValue;
		      pTrait->
			penetrance[AFFECTION_STATUS_UNAFFECTED]
			[liabilityIndex][allele2 - 1][allele1 - 1] =
			1 - penetranceValue;
		      pTrait->
			penetrance[AFFECTION_STATUS_UNKNOWN][liabilityIndex]
			[allele2 - 1][allele1 - 1] = 1;
		      KASSERT (numRet == 1,
			       "Can't get penetrance for genotype %d,%d on line %s(%d).\n",
			       allele1, allele2, sDatafileName, lineNo);
		      i++;
		      allele2++;
		      if (allele2 > numAllele)
			{
			  allele1++;
			  allele2 = allele1;
			}
		      pLine = &pLine[pos];
		    }
		  /* after we are done with this line, increment liability class counter */
		  liabilityIndex++;
		}		/* end of looping liability */
	    }			/* end of affection status */
	  if (traitType == QUANTITATIVE || traitType == COMBINED)
	    {
	      liabilityIndex = 0;
	      while (liabilityIndex < pTrait->numLiabilityClass)
		{
		  /* For quantitative trait, we expect to read means and variances 
		   * 1 line of means for each genotype 
		   * 1 line of variances for each genotype */
		  KASSERT (get_nonblank_line (line, MAX_LINE_LEN, fp, &lineNo)
			   != NULL,
			   "Can't get mean on trait %s on line %s(%d).\n",
			   sLocusName, sDatafileName, lineNo);
		  allele1 = 1;
		  allele2 = 1;
		  i = 0;
		  pLine = line;
		  /* genotype order is: 11 12 13 ... 22 23... 33 ... */
		  while (i < numGeno)
		    {
		      numRet = sscanf (pLine, "%lf %n", &value, &pos);
		      pTrait->means[liabilityIndex][allele1 - 1][allele2 -
								 1] = value;
		      pTrait->means[liabilityIndex][allele2 - 1][allele1 -
								 1] = value;
		      KASSERT (numRet == 1,
			       "Can't get mean for genotype %d,%d on line %s(%d).\n",
			       allele1, allele2, sDatafileName, lineNo);
		      i++;
		      allele2++;
		      if (allele2 > numAllele)
			{
			  allele1++;
			  allele2 = allele1;
			}
		      pLine = &pLine[pos];
		    }
		  /* variance */
		  KASSERT (get_nonblank_line (line, MAX_LINE_LEN, fp, &lineNo)
			   != NULL,
			   "Can't get variance on trait %s on line %s(%d).\n",
			   sLocusName, sDatafileName, lineNo);
		  allele1 = 1;
		  allele2 = 1;
		  i = 0;
		  pLine = line;
		  /* genotype order is: 11 12 13 ... 22 23... 33 ... */
		  while (i < numGeno)
		    {
		      numRet = sscanf (pLine, "%lf %n", &value, &pos);
		      pTrait->stddev[liabilityIndex][allele1 - 1][allele2 -
								  1] = value;
		      pTrait->stddev[liabilityIndex][allele2 - 1][allele1 -
								  1] = value;
		      KASSERT (numRet == 1,
			       "Can't get standard deviation for genotype %d,%d on line %s(%d).\n",
			       allele1, allele2, sDatafileName, lineNo);
		      i++;
		      allele2++;
		      if (allele2 > numAllele)
			{
			  allele1++;
			  allele2 = allele1;
			}
		      pLine = &pLine[pos];
		    }
		  liabilityIndex++;
		}		/* end of looping liability */
	    }			/* end of quantitative trait */
	  /* increment trait variable */
	  traitIdx++;
	}			/* end of trait variables */
      /* if we have more than 1 trait variables for the same locus, 
       * we should have covariance information */
      numVariance = pTraitLocus->numTrait * (pTraitLocus->numTrait - 1) / 2;
      i = 0;
      trait1 = 1;
      trait2 = 2;
      /* covariance order is 12 13 ... 23 ... 34 ... */
      while (i < numVariance)
	{
	  /* read a set of covariances one for each genotype between 
	   * trait1 and trait2 */
	  KASSERT (get_nonblank_line (line, MAX_LINE_LEN, fp, &lineNo) !=
		   NULL,
		   "Can't get covariance between (%d,%d) on locus %s on line %s(%d).\n",
		   sLocusName, trait1, trait2, sDatafileName, lineNo);
	  j = 0;
	  allele1 = 1;
	  allele2 = 1;
	  pLine = line;
	  /* genotype order is: 11 12 13 ... 22 23... 33 ... */
	  while (j < numGeno)
	    {
	      numRet = sscanf (pLine, "%lf %n", &value, &pos);
	      pTraitLocus->covariance[trait1][trait2][allele1 - 1][allele2 -
								   1] = value;
	      pTraitLocus->covariance[trait1][trait2][allele2 - 1][allele1 -
								   1] = value;
	      pTraitLocus->covariance[trait2][trait1][allele1 - 1][allele2 -
								   1] = value;
	      pTraitLocus->covariance[trait2][trait1][allele2 - 1][allele1 -
								   1] = value;
	      KASSERT (numRet == 1,
		       "Can't get covariance for genotype %d,%d on line %s(%d).\n",
		       allele1, allele2, sDatafileName, lineNo);
	      j++;
	      allele2++;
	      if (allele2 > numAllele)
		{
		  allele1++;
		  allele2 = allele1;
		}
	    }			/* end of genotype loop for covariance sets between two traits */
	}			/* end of covariance between traits */
    }				/* end of while */
  KASSERT (originalLocusList.numLocus == numLocus,
	   "Total Number of loci %d provided on the first line of datafile doesn't match actual number of loci %d specified in the datafile.\n",
	   numLocus, originalLocusList.numLocus);
  fclose (fp);
  return 0;
}

/* add a locus to a locus list */
Locus *
add_locus (LocusList * pLocusList, char *sName, int locusType, int numAllele)
{
  Locus *pLocus;

  if (pLocusList->maxNumLocus <= pLocusList->numLocus)
    {
      /* need to reallocate list */
      pLocusList->maxNumLocus += DEF_LOCUS_MALLOC_INCREMENT;
      pLocusList->ppLocusList = (Locus **) REALLOC ("pLocusList->ppLocusList",
						    pLocusList->ppLocusList,
						    sizeof (Locus *) *
						    pLocusList->maxNumLocus);
    }
  /* allocate space for the locus */
  pLocus = (Locus *) MALLOC ("pLocus", sizeof (Locus));
  memset (pLocus, 0, sizeof (Locus));

  /* add this locus to the list */
  pLocusList->ppLocusList[pLocusList->numLocus] = pLocus;
  pLocusList->numLocus++;

  /* set some known fields */
  strcpy (pLocus->sName, sName);
  pLocus->numAllele = numAllele;
  pLocus->numOriginalAllele = numAllele;
  /* if we are doing set recoding, we need to calculate the length of 
   * the allele set - how many integers we need to represent a allele in
   * bit format */
  //if(modelOptions.alleleSetRecodingFlag == TRUE)
  pLocus->alleleSetLen = numAllele / INT_BITS + 1;
  if (originalLocusList.alleleSetLen < pLocus->alleleSetLen)
    originalLocusList.alleleSetLen = pLocus->alleleSetLen;

  /* allocate space for frequency and count */
  pLocus->pAlleleFrequency = (double *) MALLOC ("pLocus->pAlleleFrequency",
						pLocus->numAllele *
						sizeof (double));
#ifndef NO_POLYNOMIAL
//  if(modelOptions.polynomial == TRUE) {
  pLocus->pAlleleFrequencyPolynomial =
    (Polynomial *) MALLOC ("pLocus->pAlleleFrequencyPolynomial",
			   pLocus->numAllele * sizeof (Polynomial));
//  }
#endif
  pLocus->pAlleleCount = (int *) MALLOC ("pLocus->pAlleleCount",
					 pLocus->numAllele * sizeof (int));

  /* if this is a trait locus, need to allocate more space */
  pLocus->locusType = locusType;
  if (locusType == LOCUS_TYPE_TRAIT)
    {
      pLocus->pTraitLocus = (TraitLocus *) MALLOC ("pLocus->pTraitLocus",
						   sizeof (TraitLocus));
      memset (pLocus->pTraitLocus, 0, sizeof (TraitLocus));
      pLocusList->numTraitLocus++;
    }

  return pLocus;
}

Trait *
add_trait (int trait, TraitLocus * pTraitLocus, int traitType)
{
  Trait *pTrait;

  /* allocate space */
  pTrait = (Trait *) MALLOC ("pTrait", sizeof (Trait));
  memset (pTrait, 0, sizeof (Trait));

  /* type can be either affection status or quantitative trait */
  pTrait->type = traitType;
  /* add to the list */
  pTraitLocus->pTraits[trait] = pTrait;

  return pTrait;
}

/* This function creates the base list of all possible phased genotype
 * pairs at the specified locus for the person
 * That is for a marker locus, if the phenotype of this marker is provided 
 * in pedfile without phase information (3,5), then (3,5) and (5,3) will 
 * be added to the genotype list 
 * for unknown phenotype, all combinations will be added based on number
 * of possible alleles. If for a locus, there are 3 alleles, then the 
 * following genotype list will be added for a person with unknown 
 * marker phenotypes:
 * (1,1) (1,2) (2,1) (2,2) (1,3) (3,1) (2,3) (3,2) (3,3)
 * even though some of the genotypes may not be compatible with 
 * inheritance (by looking at the parents or children)
 * The genotype elimination and set recoding (super alleles - by
 * combining nontransmitted alleles together) will be done later
 * */
int
create_baseline_marker_genotypes (int locus, Pedigree * pPedigree)
{
  int i;			/* person index */
  Person *pPerson;
  Genotype *pGenotype;
  Genotype *pGenotype2;
  int allele1, allele2;
  Locus *pLocus = originalLocusList.ppLocusList[locus];

  /* go through everyone in the pedigree */
  for (i = 0; i < pPedigree->numPerson; i++)
    {
      pPerson = pPedigree->ppPersonList[i];

      /* base on phenotype to make up the list */
      if (pPerson->pTypedFlag[locus] == 1)
	{
	  if (pPerson->pPhasedFlag[locus] == 1 ||
	      pPerson->pPhenotypeList[MOM][locus] ==
	      pPerson->pPhenotypeList[DAD][locus])
	    {
	      /* only 1 phased genotype for phased or homozygote */
	      pGenotype = add_genotype (&(pPerson->ppGenotypeList[locus]),
					&pPerson->pNumGenotype[locus],
					locus,
					pPerson->pPhenotypeList[DAD][locus],
					pPerson->pPhenotypeList[MOM][locus]);
	    }
	  else
	    {
	      /* add a pair of genotypes */
	      pGenotype = add_genotype (&pPerson->ppGenotypeList[locus],
					&pPerson->pNumGenotype[locus],
					locus,
					pPerson->pPhenotypeList[DAD][locus],
					pPerson->pPhenotypeList[MOM][locus]);
	      pGenotype2 = add_genotype (&pPerson->ppGenotypeList[locus],
					 &pPerson->pNumGenotype[locus],
					 locus,
					 pPerson->pPhenotypeList[MOM][locus],
					 pPerson->pPhenotypeList[DAD][locus]);
	      pGenotype->pDualGenotype = pGenotype2;
	      pGenotype2->pDualGenotype = pGenotype;

	    }
	}
      else
	{
	  /* marker phenotype is not known, add all possible combinations */
	  allele1 = 1;
	  allele2 = 1;
	  for (allele1 = 1; allele1 <= pLocus->numAllele; allele1++)
	    {
	      for (allele2 = allele1; allele2 <= pLocus->numAllele; allele2++)
		{
		  pGenotype = add_genotype (&pPerson->ppGenotypeList[locus],
					    &pPerson->pNumGenotype[locus],
					    locus, allele1, allele2);
		  if (allele2 != allele1)
		    {
		      pGenotype2 =
			add_genotype (&pPerson->ppGenotypeList[locus],
				      &pPerson->pNumGenotype[locus], locus,
				      allele2, allele1);
		      pGenotype->pDualGenotype = pGenotype2;
		      pGenotype2->pDualGenotype = pGenotype;
		    }
		}
	    }
	}
    }

  return 0;
}

/* similar to add marker baseline genotypes 
 * except it uses penetrance matrix to determine possible genotypes
 * in addition to disease status or quantitative trait values
 * */
int
create_baseline_trait_genotypes (int locus, Pedigree * pPedigree)
{
  int i;			/* person index */
  Person *pPerson;
  Genotype *pGenotype, *pGenotype2;
  int allele1, allele2;
  Locus *pLocus = originalLocusList.ppLocusList[locus];
  //TraitLocus *pTraitLocus = pLocus->pTraitLocus;
  double pen;
#ifndef NO_POLYNOMIAL
  Polynomial *penPolynomial;
#endif

  /* go through everyone in the pedigree */
  for (i = 0; i < pPedigree->numPerson; i++)
    {
      pPerson = pPedigree->ppPersonList[i];

      /* go through all possible genotype combinations */
      allele1 = 1;
      allele2 = 1;
      for (allele1 = 1; allele1 <= pLocus->numAllele; allele1++)
	{
	  for (allele2 = allele1; allele2 <= pLocus->numAllele; allele2++)
	    {
#ifndef NO_POLYNOMIAL		// including polynomial code
	      if (modelOptions.polynomial == TRUE)
		{
		  compute_penetrance (pPerson, locus, allele1, allele2,
				      &penPolynomial);
		}
	      else
		compute_penetrance (pPerson, locus, allele1, allele2, &pen);
#else
	      compute_penetrance (pPerson, locus, allele1, allele2, &pen);
#endif
	      /* if this genotype hasn't been rejected, then add it 
	       * under polynomial, penetrance is a variable (parameter), we
	       * can't check against one single value, so we include all
	       * trait genotypes */
#ifndef NO_POLYNOMIAL		// including polynomial code
	      if (modelOptions.polynomial == TRUE)
		{
		  pGenotype = add_genotype (&pPerson->ppGenotypeList[locus],
					    &pPerson->pNumGenotype[locus],
					    locus, allele1, allele2);
		  pGenotype->penetrancePolynomial = penPolynomial;
		  if (allele1 != allele2)
		    {
		      pGenotype2 =
			add_genotype (&pPerson->ppGenotypeList[locus],
				      &pPerson->pNumGenotype[locus], locus,
				      allele2, allele1);
		      pGenotype->pDualGenotype = pGenotype2;
		      pGenotype2->pDualGenotype = pGenotype;
		      pGenotype2->penetrancePolynomial = penPolynomial;
		    }
		}
	      else
		{
//              if ( pen > 0) {
		  pGenotype = add_genotype (&pPerson->ppGenotypeList[locus],
					    &pPerson->pNumGenotype[locus],
					    locus, allele1, allele2);
		  pGenotype->penetrance = pen;
		  if (allele1 != allele2)
		    {
		      pGenotype2 =
			add_genotype (&pPerson->ppGenotypeList[locus],
				      &pPerson->pNumGenotype[locus], locus,
				      allele2, allele1);
		      pGenotype->pDualGenotype = pGenotype2;
		      pGenotype2->pDualGenotype = pGenotype;
		      pGenotype2->penetrance = pen;
		    }
//              }
		}
#else
//      if ( pen > 0) {
	      pGenotype = add_genotype (&pPerson->ppGenotypeList[locus],
					&pPerson->pNumGenotype[locus],
					locus, allele1, allele2);
	      pGenotype->penetrance = pen;
	      if (allele1 != allele2)
		{
		  pGenotype2 = add_genotype (&pPerson->ppGenotypeList[locus],
					     &pPerson->pNumGenotype[locus],
					     locus, allele2, allele1);
		  pGenotype->pDualGenotype = pGenotype2;
		  pGenotype2->pDualGenotype = pGenotype;
		  pGenotype2->penetrance = pen;
		}
//      }
#endif
	    }
	}
    }

  return 0;
}

/* given a genotype, derive the penetrance either by looking up 
 * penetrance matrix under dichotomous trait case, or through some
 * function calls under quantitative trait and combined cases */

void
compute_penetrance (Person * pPerson, int locus, int allele1, int allele2,
		    void *pen)
{
  Locus *pLocus = originalLocusList.ppLocusList[locus];
  TraitLocus *pTraitLocus = pLocus->pTraitLocus;
  Trait *pTrait;
  int affectionStatus;
  int liabilityClass;
  int i, j;
//  double trait1, trait2;
//  double mean1, mean2;
//  double covariance;
  double trait, mean, stddev, temp, temp2;
  double sigma, df;
  int status;
  double bound;
  int tempInt;

  /* now go through all related traits and see whether 
   * this genotype is compatible with each trait value if
   * the trait value is known */
  /* for now, we assume we either have only 1 affection status trait
   * or 1 or more quantitative traits only 
   * i.e. there is no mixing between affection status trait and 
   * liability trait */
  i = 0;
  while (i < pTraitLocus->numTrait)
    {
      pTrait = pTraitLocus->pTraits[i];
      if (pTrait->numLiabilityClass > 1)
	liabilityClass = pPerson->ppLiabilityClass[locus][i];
      else
	liabilityClass = 1;
      if (pTrait->type == DICHOTOMOUS)
	{
	  if (pPerson->ppTraitKnown[locus][i] == TRUE)
	    {
	      affectionStatus = (int) pPerson->ppTraitValue[locus][i];

#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  if (affectionStatus == AFFECTION_STATUS_AFFECTED)
		    {
		      /* build the penetrance polynomial - single variable */
		      char vName[100];
		      sprintf (vName, "p[%d][%d][%d][%d]",
			       AFFECTION_STATUS_AFFECTED, liabilityClass - 1,
			       allele1 - 1, allele2 - 1);

		      *(Polynomial **) pen =
			variableExp (&pTrait->
				     penetrance[AFFECTION_STATUS_AFFECTED]
				     [liabilityClass - 1][allele1 -
							  1][allele2 - 1],
				     "");
		    }
		  else
		    {
		      /* unaffected: 1 - pen[affected] 
		       * 1 - pTrait->penetrance[AFFECTION_STATUS_AFFECTED][liabilityClass - 1]
		       *                                 [allele1-1][allele2-1] */
		      char vName[100];
		      sprintf (vName, "p[%d][%d][%d][%d]",
			       AFFECTION_STATUS_AFFECTED, liabilityClass - 1,
			       allele1 - 1, allele2 - 1);
		      *(Polynomial **) pen =
			plusExp (2, 1.0, constantExp (1), -1.0,
				 variableExp (&pTrait->
					      penetrance
					      [AFFECTION_STATUS_AFFECTED]
					      [liabilityClass - 1][allele1 -
								   1][allele2
								      - 1],
					      vName));
		    }
		}
	      else
		{
		  /* polynomial flag is not turned on */
		  *(double *) pen =
		    pTrait->penetrance[affectionStatus][liabilityClass - 1]
		    [allele1 - 1][allele2 - 1];
//          fprintf(stderr,"AffectionStatus=%d liabilityClass-1=%d allele1-1=%d allele2-1=%d pen=%f\n",
//                        affectionStatus,liabilityClass - 1,allele1-1,allele2-1,
//                        pTrait->penetrance[AFFECTION_STATUS_AFFECTED][liabilityClass - 1]
//                                          [allele1-1][allele2-1]);

		}
#else
	      /* no polynomial */
	      *(double *) pen =
		pTrait->penetrance[affectionStatus][liabilityClass - 1]
		[allele1 - 1][allele2 - 1];
#endif
	    }
	  else
	    {
	      affectionStatus = AFFECTION_STATUS_UNKNOWN;
	      /* when the affection status is unknown, the penetrance is set to 1 */

#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  *(Polynomial **) pen = constantExp (1);
		}
	      else
		*(double *) pen = 1;
#else
	      *(double *) pen = 1;
#endif
	    }
	}			/* end of affection status trait handling */
      else if (pTrait->type == QUANTITATIVE || pTrait->type == COMBINED)
	{
	  /* make sure we have all the trait values for this locus */
	  for (j = 0; j < pTraitLocus->numTrait; j++)
	    {
	      /* if the quantitative trait value is not known for any trait, then
	       * we can't calculate "exact" penetrance, just return 1 */
	      if (pPerson->ppTraitKnown[locus][j] == FALSE)
		{
		  *(double *) pen = 1;
		}
	    }

	  /*
	     pen = 0;
	     for (i = 0; i < pTraitLocus->numTrait; i++) {
	     trait1 = pPerson->ppTraitValue[locus][i];
	     mean1 = pTraitLocus->pTraits[i]->means[allele1-1][allele2-1];
	     for (j = 0; i < pTraitLocus->numTrait; i++) {
	     trait2 = pPerson->ppTraitValue[locus][j];
	     mean2 = pTraitLocus->pTraits[i]->means[allele1-1][allele2-1];

	     covariance = pTraitLocus->covariance[i][j][allele1-1][allele2-1];
	     pen += (trait1 - mean1) * (trait2 - mean2) * covariance;
	     }
	     } 
	   */

	  /* trait value is known */
	  trait = pPerson->ppTraitValue[locus][i];
	  mean = pTrait->means[liabilityClass - 1][allele1 - 1][allele2 - 1];
	  stddev =
	    pTrait->stddev[liabilityClass - 1][allele1 - 1][allele2 - 1];
	  /* working on the assumption of just 1 QT trait */
	  if (pTrait->functionQT == QT_FUNCTION_NORMAL)
	    {
	      if (pTrait->type == COMBINED &&
		  (trait >= pTrait->lessCutoffFlag - 0.000001 &&
		   trait <= pTrait->lessCutoffFlag + 0.000001))
		{
		  /* CDF of less than cutoff */
		  bound = -1;
		  tempInt = 1;
		  cdfnor (&tempInt, (double *) pen, &temp,
			  &pTrait->cutoffValue, &mean, &stddev, &status,
			  &bound);
		}
	      else if (pTrait->type == COMBINED &&
		       (trait >= pTrait->moreCutoffFlag - 0.000001 &&
			trait <= pTrait->moreCutoffFlag + 0.000001))
		{
		  /* CDF of more than cutoff */
		  bound = -1;
		  tempInt = 1;
		  cdfnor (&tempInt, &temp, (double *) pen,
			  &pTrait->cutoffValue, &mean, &stddev, &status,
			  &bound);
		}
	      else
		{
		  temp = (trait - mean) / stddev;
		  *(double *) pen =
		    exp (-0.5 * temp * temp) / (stddev * sqrt (2 * PI));
		}
	    }
	  else if (pTrait->functionQT == QT_FUNCTION_T)
	    {
	      if (pTrait->type == COMBINED &&
		  (trait >= pTrait->lessCutoffFlag - 0.000001 &&
		   trait <= pTrait->lessCutoffFlag + 0.000001))
		{
		  /* CDF of less than cutoff */
		  bound = -1;
		  tempInt = 1;
		  cdftnc (&tempInt, (double *) pen, &temp,
			  &pTrait->cutoffValue, &df, &mean, &status, &bound);
		}
	      else if (pTrait->type == COMBINED &&
		       (trait >= pTrait->moreCutoffFlag - 0.000001 &&
			trait <= pTrait->moreCutoffFlag + 0.000001))
		{
		  /* CDF of more than cutoff */
		  bound = -1;
		  tempInt = 1;
		  cdftnc (&tempInt, &temp, (double *) pen,
			  &pTrait->cutoffValue, &df, &mean, &status, &bound);
		}
	      else
		{
		  /* non-central T-distribution */
		  df = pTrait->dfQT;
		  /*
		     sigma = (trait - mean)/ stddev;
		     temp = gsl_sf_gamma((df+1)/2);
		     temp2 = gsl_sf_gamma(df/2);
		     *(double *)pen = temp / (sqrt(df*PI) * temp2 
		     * pow(1.0+sigma*sigma/(df-2), (df+1)/2)); 
		   */
		  sigma = 0;
		  temp = gsl_sf_gamma ((df + 1) / 2) / (sqrt (df * PI) *
							gsl_sf_gamma (df /
								      2));
		  temp2 = (trait - mean) / (stddev * sqrt ((df - 2) / df));
		  *(double *) pen =
		    temp / pow (1 + temp2 * temp2 / df, (df + 1) / 2);
		}
	    }			/* t-dsitribution */
	}			/* end of quantitative trait handling */
      i++;
    }				/* move to next trait for the same locus */
}

Genotype *
add_genotype (Genotype ** ppList, int *pCount, int locusIndex,
	      int allele1, int allele2)
{
  Genotype *pGenotype;
  int numInts;

  /* allocate space for the genotype */
  pGenotype = (Genotype *) MALLOC ("pGenotype", sizeof (Genotype));
  memset (pGenotype, 0, sizeof (Genotype));
  pGenotype->penetrance = 1;

  /* add this to the top of the genotype list */
  pGenotype->pNext = *ppList;
  *ppList = pGenotype;
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      pGenotype->penetrancePolynomial = constantExp (1);
    }
#endif

  /* increase the counter */
  (*pCount)++;

  /* set alleles */
  pGenotype->allele[DAD] = allele1;
  pGenotype->allele[MOM] = allele2;

  /* if we are required to do allele set recoding */
  //if(modelOptions.alleleSetRecodingFlag == TRUE) {
  /* allocate space for the allele bits - as potentially there are more
   * than 32 -1 possible alleles */
  numInts = originalLocusList.alleleSetLen;
  pGenotype->pAlleleBits[DAD] = (int *) MALLOC ("pGenotype->pAlleleBits[DAD]",
						sizeof (int) * numInts);
  pGenotype->pAlleleBits[MOM] = (int *) MALLOC ("pGenotype->pAlleleBits[MOM]",
						sizeof (int) * numInts);
  memset (pGenotype->pAlleleBits[DAD], 0, sizeof (int) * numInts);
  memset (pGenotype->pAlleleBits[MOM], 0, sizeof (int) * numInts);
  /* set the bits */
  set_allele_bit (allele1, pGenotype->pAlleleBits[DAD]);
  set_allele_bit (allele2, pGenotype->pAlleleBits[MOM]);
  //}

  return pGenotype;
}

int
remove_genotype (Genotype ** pHead, Genotype * pGenotype, int *pCount)
{
  Genotype *pPrev = NULL;
  Genotype *pCurr = NULL;

  if (*pHead == NULL)
    return -1;
  if (*pHead == pGenotype)
    {
      *pHead = pGenotype->pNext;
    }
  else
    {
      pCurr = (*pHead)->pNext;
      pPrev = (*pHead);
      while (pCurr != NULL && pCurr != pGenotype)
	{
	  pPrev = pCurr;
	  pCurr = pCurr->pNext;
	}
      /* failed to find the genotype */
      if (pCurr == NULL)
	return -1;
      pPrev->pNext = pGenotype->pNext;
    }

  /* free the space */
  free (pGenotype->pAlleleBits[DAD]);
  free (pGenotype->pAlleleBits[MOM]);
  free (pGenotype);
  /* decrement the counter */
  (*pCount)--;

  return 0;
}


/* This genotype list should be phased - for printing purpose, 
 * it doesn't matter */
void
print_person_locus_genotype_list (Person * pPerson, int locus)
{
  Genotype *pGenotype = pPerson->ppGenotypeList[locus];

  /* print out person lable */
  logMsg (LOGGENOELIM, LOGDEBUG, "    Person %s: ", pPerson->sID);
  while (pGenotype != NULL)
    {
      logMsg (LOGGENOELIM, LOGDEBUG, "(%d,%d) ", pGenotype->allele[DAD],
	      pGenotype->allele[MOM]);
      pGenotype = pGenotype->pNext;
    }
  logMsg (LOGGENOELIM, LOGDEBUG, "\n");

  return;
}

void
print_pedigree_locus_genotype_list (Pedigree * pPedigree, int locus)
{
  int i;

  logMsg (LOGGENOELIM, LOGDEBUG, "Pedigree %3s Locus %d: \n",
	  pPedigree->sPedigreeID, locus + 1);
  for (i = 0; i < pPedigree->numPerson; i++)
    {
      print_person_locus_genotype_list (pPedigree->ppPersonList[i], locus);
    }
  return;
}

/* This procedure set the genotype weight of each genotype for each
 * person in a pedigree 
 * for a homozygous genotype, the weight is p*p
 * for a phased heterozygous genotype, the weight is pq
 * This will facilitate likelihood calculation
 * 
 * */
int
set_genotype_weight (Pedigree * pPedigree)
{
  int locus;
  Person *pPerson;
  Genotype *pGenotype;
  int i;
  Locus *pLocus;
  double alleleFreq[2];
#ifndef NO_POLYNOMIAL
  Polynomial *alleleFreqPolynomial[2];
#endif

  for (i = 0; i < pPedigree->numPerson; i++)
    {
      pPerson = pPedigree->ppPersonList[i];
      memcpy (&pPerson->pSavedNumGenotype[0], &pPerson->pNumGenotype[0],
	      sizeof (int) * originalLocusList.numLocus);
      for (locus = 0; locus < originalLocusList.numLocus; locus++)
	{
	  pLocus = originalLocusList.ppLocusList[locus];
	  pGenotype = pPerson->ppGenotypeList[locus];
	  while (pGenotype)
	    {
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  if (pGenotype->allele[DAD] <= pLocus->numAllele)
		    {
		      char vName[100];
		      sprintf (vName, "f[%d][%d]", locus,
			       pGenotype->allele[DAD] - 1);
		      alleleFreqPolynomial[DAD] =
			variableExp (&pLocus->
				     pAlleleFrequency[pGenotype->allele[DAD] -
						      1], vName);
		    }
		  else
		    {
		      alleleFreqPolynomial[DAD] =
			pLocus->ppAlleleSetList[pGenotype->allele[DAD] -
						1]->sumFreqPolynomial;
		    }
		  if (pGenotype->allele[MOM] <= pLocus->numAllele)
		    {
		      char vName[100];
		      sprintf (vName, "f[%d][%d]", locus,
			       pGenotype->allele[MOM] - 1);
		      alleleFreqPolynomial[MOM] =
			variableExp (&pLocus->
				     pAlleleFrequency[pGenotype->allele[MOM] -
						      1], vName);
		    }
		  else
		    alleleFreqPolynomial[MOM] =
		      pLocus->ppAlleleSetList[pGenotype->allele[MOM] -
					      1]->sumFreqPolynomial;
		}
	      else
		{
		  if (pGenotype->allele[DAD] <= pLocus->numAllele)
		    {
		      alleleFreq[DAD] =
			pLocus->pAlleleFrequency[pGenotype->allele[DAD] - 1];
		    }
		  else
		    alleleFreq[DAD] =
		      pLocus->ppAlleleSetList[pGenotype->allele[DAD] -
					      1]->sumFreq;
		  if (pGenotype->allele[MOM] <= pLocus->numAllele)
		    {
		      alleleFreq[MOM] =
			pLocus->pAlleleFrequency[pGenotype->allele[MOM] - 1];
		    }
		  else
		    alleleFreq[MOM] =
		      pLocus->ppAlleleSetList[pGenotype->allele[MOM] -
					      1]->sumFreq;
		}
#else
	      if (pGenotype->allele[DAD] <= pLocus->numAllele)
		alleleFreq[DAD] =
		  pLocus->pAlleleFrequency[pGenotype->allele[DAD] - 1];
	      else
		alleleFreq[DAD] =
		  pLocus->ppAlleleSetList[pGenotype->allele[DAD] -
					  1]->sumFreq;
	      if (pGenotype->allele[MOM] <= pLocus->numAllele)
		alleleFreq[MOM] =
		  pLocus->pAlleleFrequency[pGenotype->allele[MOM] - 1];
	      else
		alleleFreq[MOM] =
		  pLocus->ppAlleleSetList[pGenotype->allele[MOM] -
					  1]->sumFreq;
#endif
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		{
		  /* build the polynomial sum */
		  pGenotype->weightPolynomial =
		    timesExp (2, alleleFreqPolynomial[DAD], 1,
			      alleleFreqPolynomial[MOM], 1);
		}
	      else
		{
		  pGenotype->weight = alleleFreq[DAD] * alleleFreq[MOM];
		}
#else
	      pGenotype->weight = alleleFreq[DAD] * alleleFreq[MOM];
#endif
	      pGenotype = pGenotype->pNext;
	    }
	}
    }

  return 0;
}

/* This procedure set the position of each genotype in this person's
 * genotype list for each person in a pedigree 
 * This will facilitate storing and retrieving multi locus 
 * conditional genotype likelihood
 * 
 * */
int
set_genotype_position (Pedigree * pPedigree)
{
  int locus;
  Person *pPerson;
  Genotype *pGenotype;
  int i;
  int position;
  Locus *pLocus;

  for (i = 0; i < pPedigree->numPerson; i++)
    {
      pPerson = pPedigree->ppPersonList[i];
      for (locus = 0; locus < originalLocusList.numLocus; locus++)
	{
	  pLocus = originalLocusList.ppLocusList[locus];
	  pGenotype = pPerson->ppGenotypeList[locus];
	  position = 0;
	  while (pGenotype)
	    {
	      pGenotype->position = position;
	      pGenotype = pGenotype->pNext;
	      position++;
	    }
	}
    }

  return 0;
}

int
allocate_multi_locus_genotype_storage (Pedigree * pPedigree)
{
  int locus;
  int origLocus;
  Person *pPerson;
  int i;
  int size;

  for (i = 0; i < pPedigree->numPerson; i++)
    {
      pPerson = pPedigree->ppPersonList[i];
      size = 1;
      for (locus = 0; locus < locusList.numLocus; locus++)
	{
	  origLocus = locusList.pLocusIndex[locus];
	  size *= pPerson->pNumGenotype[origLocus];
	}
      /* allocate space */
      pPerson->pLikelihood =
	(ConditionalLikelihood *) MALLOC ("pPerson->pLikelihood ",
					  sizeof (ConditionalLikelihood) *
					  size);
      memset (pPerson->pLikelihood, 0, sizeof (ConditionalLikelihood) * size);
      pPerson->numConditionals = size;
    }

  return 0;
}

int
free_multi_locus_genotype_storage (Pedigree * pPedigree)
{
  Person *pPerson;
  int i;

  for (i = 0; i < pPedigree->numPerson; i++)
    {
      pPerson = pPedigree->ppPersonList[i];
      if (pPerson->pLikelihood != NULL)
	{
	  /* free space */
	  FREE ("pPerson->pLikelihood", pPerson->pLikelihood);
	}
    }

  return 0;
}

int
initialize_multi_locus_genotype (Pedigree * pPedigree)
{
  Person *pPerson;
  int i, j;
  int size;

  for (i = 0; i < pPedigree->numPerson; i++)
    {
      pPerson = pPedigree->ppPersonList[i];
      size = pPerson->numConditionals;
      for (j = 0; j < size; j++)
	{
	  pPerson->pLikelihood[j].touchedFlag = FALSE;
#ifndef NO_POLYNOMIAL
	  if (modelOptions.polynomial == TRUE)
	    {
	      pPerson->pLikelihood[j].likelihoodPolynomial = constantExp (1);
	      pPerson->pLikelihood[j].weightPolynomial = constantExp (1);
	    }
	  else
	    {
	      pPerson->pLikelihood[j].likelihood = 1;
	      pPerson->pLikelihood[j].weight = 1;
	    }
#else
	  pPerson->pLikelihood[j].likelihood = 1;
	  pPerson->pLikelihood[j].weight = 1;
#endif
	}
      pPerson->numConditionals = size;
    }

  return 0;
}

void
setup_LD_haplotype_freq (LDLoci * pLDLoci)
{
  Locus *pLocus1, *pLocus2;
  int numAllele1, numAllele2;
  int i, j;
  double p1, p2, q1, q2;
  double DPrime, DValue, haploFreq;
  double maxD, minD, sum;

  pLocus1 = originalLocusList.ppLocusList[pLDLoci->locus1];
  pLocus2 = originalLocusList.ppLocusList[pLDLoci->locus2];
  numAllele1 = pLocus1->numOriginalAllele;
  numAllele2 = pLocus2->numOriginalAllele;
  for (i = 0; i < numAllele1 - 1; i++)
    {
      p1 = pLocus1->pAlleleFrequency[i];
      p2 = 1 - p1;
      sum = 0;
      for (j = 0; j < numAllele2 - 1; j++)
	{
	  q1 = pLocus2->pAlleleFrequency[j];
	  q2 = 1 - q1;

	  if (-p1 * q1 > -p2 * q2)
	    minD = -p1 * q1 + LD_E;
	  else
	    minD = -p2 * q2 + LD_E;

	  if (p1 * q2 < p2 * q1)
	    maxD = p1 * q2 - LD_E;
	  else
	    maxD = p2 * q1 - LD_E;

	  DPrime = pLDLoci->ppDPrime[i][j];
	  if (DPrime > 0)
	    DValue = DPrime * maxD;
	  else
	    DValue = DPrime * minD;
	  pLDLoci->ppDValue[i][j] = DValue;

	  haploFreq = p1 * q1 + DValue;

	  pLDLoci->ppHaploFreq[i][j] = haploFreq;

	  sum += haploFreq;
	}			/* end of looping the second marker allele frequencies */
      pLDLoci->ppHaploFreq[i][j] = p1 - sum;
    }				/* end of looping the first marker allele frequencies */

  /* fill out haplotype frequencies for column [**, n] */
  for (j = 0; j < numAllele2; j++)
    {
      q1 = pLocus2->pAlleleFrequency[j];
      sum = 0;
      for (i = 0; i < numAllele1 - 1; i++)
	{
	  sum += pLDLoci->ppHaploFreq[i][j];
	}

      pLDLoci->ppHaploFreq[i][j] = q1 - sum;


    }
}

/* find the place holder for the LD parameters between the given
 * two loci */
LDLoci *
find_LD_loci (int locus1, int locus2)
{
  int i;
  LDLoci *pLDLoci;

  for (i = 0; i < originalLocusList.numLDLoci; i++)
    {
      pLDLoci = &originalLocusList.pLDLoci[i];
      if ((pLDLoci->locus1 == locus1 && pLDLoci->locus2 == locus2) ||
	  (pLDLoci->locus1 == locus2 && pLDLoci->locus2 == locus1))
	return pLDLoci;
    }

  return NULL;
}

/* initialization procedure:
 * set recoding, genotype elimination */
int
initialize_loci (PedigreeSet * pPedigreeSet)
{
  int locus;
  int ped;
  Pedigree *pPedigree;
  Locus *pLocus;

  /* go through all loci in the original locus list */
  locus = 0;
  while (locus < originalLocusList.numLocus)
    {
      pLocus = originalLocusList.ppLocusList[locus];
      ped = 0;
      /* go through all peidgress in this set */
      while (ped < pPedigreeSet->numPedigree)
	{
	  pPedigree = pPedigreeSet->ppPedigreeSet[ped];
	  /* depending on the locus type, we call different function */

	  if (pLocus->locusType == LOCUS_TYPE_TRAIT)
	    create_baseline_trait_genotypes (locus, pPedigree);
	  else
	    create_baseline_marker_genotypes (locus, pPedigree);

	  logMsg (LOGGENOELIM, LOGDEBUG, "Baseline Genotype Lists:\n");
	  print_pedigree_locus_genotype_list (pPedigree, locus);


	  /* first step is do the set recoding: this should help speed up
	   * genotype elimination process */

	  //fprintf(stderr,"Before 1th allele_set_recoding\n");
	  allele_set_recoding (locus, pPedigree);
	  //fprintf(stderr,"After  allele_set_recoding\n");

	  /* do genotype elimination next */
	  pedigree_genotype_elimination (locus, pPedigree);
	  logMsg (LOGGENOELIM, LOGDEBUG,
		  "Genotype Lists after genotype elimination :\n");
	  print_pedigree_locus_genotype_list (pPedigree, locus);

	  //fprintf(stderr,"Before 2nd allele_set_recoding\n");
	  /* do the set recoding again */
	  allele_set_recoding (locus, pPedigree);
	  //fprintf(stderr,"After  allele_set_recoding\n");

	  logMsg (LOGGENOELIM, LOGDEBUG,
		  "Genotype Lists after set recoding :\n");
	  print_pedigree_locus_genotype_list (pPedigree, locus);

	  ped++;
	}
      locus++;
    }

  /* after genotype elimination, we can set genotype weight and position 
   * in the genotype list 
   * this will be done on the loci in the original locus list */
  /* If we ever need to loop over marker allele frequencies, 
   * the super allele frequencies, genotype weights should be set out of 
   * this routine and inside the looping */
  ped = 0;
  while (ped < pPedigreeSet->numPedigree)
    {
      pPedigree = pPedigreeSet->ppPedigreeSet[ped];
      set_genotype_weight (pPedigree);
      set_genotype_position (pPedigree);
      ped++;


/*
    int i,j,k;
    Person *pPerson;
    Genotype *g;
    for(i=0;i<pPedigree->numPerson;i++)
    {
       pPerson=pPedigree->ppPersonList[i];
       for(j=0;j<originalLocusList.numLocus;j++)
       {
           g=pPerson->ppGenotypeList[j]; 
           k=0;
           while(g)
           {

#ifndef NO_POLYNOMIAL  // including polynomial code
               if (modelOptions.polynomial == TRUE )
               {
               expPrinting(g->penetrancePolynomial);
               fprintf(stderr," Person %d locus %d genotype %d of %d (%d  %d)\n",
                       i+1,j+1,k+1,pPerson->pNumGenotype[j],
                       g->allele[0],
                       g->allele[1]);
               }
               else
               fprintf(stderr,"%f  Person %d locus %d genotype %d of %d (%d  %d)\n",
                       g->penetrance,
                       i+1,j+1,k+1,pPerson->pNumGenotype[j],
                       g->allele[0],
                       g->allele[1]);
#else
               fprintf(stderr,"%f  Person %d locus %d genotype %d of %d (%d  %d)\n",
                       g->penetrance,
                       i+1,j+1,k+1,pPerson->pNumGenotype[j],
                       g->allele[0],
                       g->allele[1]);
#endif

               g = g->pNext;
               k++;
           }
       }
    }
*/


    }
  return 0;
}


int
update_loci (PedigreeSet * pPedigreeSet)
{
  int ped;
  Pedigree *pPedigree;
  /* update genotype weight and position 
   * in the genotype list 
   * this will be done on the loci in the original locus list */
  ped = 0;
  while (ped < pPedigreeSet->numPedigree)
    {
      pPedigree = pPedigreeSet->ppPedigreeSet[ped];
      set_genotype_weight (pPedigree);
      set_genotype_position (pPedigree);
      ped++;
    }

  return 0;
}

int
update_penetrance (PedigreeSet * pPedigreeSet, int locus)
{
  Person *pPerson;
  Genotype *pGenotype;
  int i;
  double pen;
  int ped;
  Pedigree *pPedigree;
  /* update genotype weight and position 
   * in the genotype list 
   * this will be done on the loci in the original locus list */
  ped = 0;
  while (ped < pPedigreeSet->numPedigree)
    {
      pPedigree = pPedigreeSet->ppPedigreeSet[ped];
      for (i = 0; i < pPedigree->numPerson; i++)
	{
	  pPerson = pPedigree->ppPersonList[i];
	  pGenotype = pPerson->ppGenotypeList[locus];
	  while (pGenotype)
	    {
	      compute_penetrance (pPerson, locus,
				  pGenotype->allele[0], pGenotype->allele[1],
				  &pen);
	      pGenotype->penetrance = pen;
	      pGenotype = pGenotype->pNext;
	    }
	  pGenotype = pPerson->ppSavedGenotypeList[locus];
	  while (pGenotype)
	    {
	      compute_penetrance (pPerson, locus,
				  pGenotype->allele[0], pGenotype->allele[1],
				  &pen);
	      pGenotype->penetrance = pen;
	      pGenotype = pGenotype->pNext;
	    }
	  pGenotype = pPerson->ppShadowGenotypeList[locus];
	  while (pGenotype)
	    {
	      compute_penetrance (pPerson, locus,
				  pGenotype->allele[0], pGenotype->allele[1],
				  &pen);
	      pGenotype->penetrance = pen;
	      pGenotype = pGenotype->pNext;
	    }
	}
      ped++;
    }
  return 0;
}
