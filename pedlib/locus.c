
/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

/* This file contains functions to read in marker map file and 
 * marker file 
 * 
 * So map file usually is a very stable file 
 * pedfile marker phenotype spec order needs to match the marker file, i.e.
 * based on map order
 *
 * Map file will be read in first, as it potentially has a bigger set of
 * markers (never should be smaller set)
 * 
 * */

char *locusVersion = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"

#include "pedlib.h"
#include "utils.h"		/* for logging */
#include "tools.h"
#include "polynomial.h"

/* global variables */
Map map;
LocusList originalLocusList;
SubLocusList *locusList;

extern Polynomial *constant0Poly;
extern Polynomial *constant1Poly;

/* internal functions */
MapUnit *add_map_unit (Map * map);
Locus *add_locus (LocusList * pLocusList, char *sName, int locusType);
int add_allele (Locus * pLocus, char *sAlleleName, double freq);

Genotype *add_genotype (Genotype ** ppList, int *pCount, int locusIndex,
			int allele1, int allele2);

Trait *add_trait (int trait, TraitLocus * pTraitLocus, int traitType);
void compute_penetrance (Person * pPerson, int locus,
			 int allele1, int allele2, void *pen);
void free_allele_set (AlleleSet * pAlleleSet);
void free_person (Person * pPerson);
void free_nuclear_family (NuclearFamily * pNucFam);

/* hard removal of genotype flag */
int removeGenotypeFlag = TRUE;

void
set_removeGenotypeFlag (int flag)
{
  removeGenotypeFlag = flag;
}


/* read in mapfile 
 * all marker information will be saved in the super marker list 
 * all markers should be on the same chromosome */
int
read_mapfile (char *sMapfileName)
{
  FILE *fp = NULL;
  int lineNo = 0;
  char line[MAX_LINE_LEN];
  int numRet;
  MapUnit *pMarker;		/* current marker */
  MapUnit *pPrevMarker = NULL;	/* previous marker */
  char *p;
  int chr;
  char sName[MAX_LINE_LEN];
  double sexAvgPos, malePos, femalePos;
  int basePairLoc;

  fp = fopen (sMapfileName, "r");
  KASSERT (fp != NULL, "Can't open map file %s for read. Exiting!\n",
	   sMapfileName);
  /* default map function is KOSAMBI */
  map.mapFunction = MAP_FUNCTION_KOSAMBI;
  /* skip comment and blank lines 
   *  read in map function indication - required */
  while (fgets (line, MAX_LINE_LEN, fp)) {
    lineNo++;
    if (is_line_blank_or_comment (line))
      continue;
    /* assume we either have got our map function spec or we got marker data */
    p = strtok (line, " =\t");
    /* ignoring the case */
    if (p != NULL && !strcasecmp (p, "mapFunction")) {
      p = strtok (NULL, "\t =");
      if (p != NULL && (p[0] == 'h' || p[0] == 'H')) {
	map.mapFunction = MAP_FUNCTION_HALDANE;
      } else {
	map.mapFunction = MAP_FUNCTION_KOSAMBI;
      }
      /* we can move on to the actual marker section */
      break;
    } else if (p != NULL && !strncasecmp (p, "chr", 3)) {
      /* we get the header line */
      break;
    }
  }

  KASSERT (feof (fp) == 0,
	   "No marker information at all in the map file %s (Total # of Lines=%d). Exiting!\n",
	   sMapfileName, lineNo);

  /* read one line for each marker */
  pPrevMarker = NULL;
  while (fgets (line, MAX_LINE_LEN, fp)) {
    lineNo++;
    if (is_line_blank_or_comment (line))
      continue;
    malePos = -100;
    femalePos = -100;
    basePairLoc = -1;
    numRet =
      sscanf (line, "%d %s %lf %lf %lf %d", &chr, sName, &sexAvgPos,
	      &malePos, &femalePos, &basePairLoc);

    if (numRet == 0)
      /* assume this is a header line or comment line */
      continue;

    /* sex specific map and/or base pair location are not always known, so not required */
    KASSERT (numRet >= 3,
	     "Marker map file %s line %d seems not complete.\n",
	     sMapfileName, lineNo);
    KASSERT (chr > 0,
	     "Chromosome is not greater than 0 (%d) in file %s(%d).\n", chr,
	     sMapfileName, lineNo);

    pMarker = add_map_unit (&map);
    pMarker->chromosome = chr;
    strcpy (pMarker->sName, sName);
    pMarker->mapPos[MAP_SEX_AVERAGE] = sexAvgPos;
    pMarker->mapPos[MAP_FEMALE] = malePos;
    pMarker->mapPos[MAP_MALE] = femalePos;
    pMarker->basePairLocation = basePairLoc;
    /* make sure marker is specified based on map order */
    if (pPrevMarker != NULL
	&& (pPrevMarker->mapPos[MAP_SEX_AVERAGE] >
	    pMarker->mapPos[MAP_SEX_AVERAGE])) {
      /* this will kick the program out */
      KASSERT (1 == 0,
	       "Marker map given by %s is out of order between %s and %s.\n",
	       sMapfileName, pPrevMarker->sName, pMarker->sName);
    }

    /* we should move on to next marker */
    pPrevMarker = pMarker;
  }

  fclose (fp);
  return 0;
}

/* add a marker to the map */
MapUnit *
add_map_unit (Map * pMap)
{
  MapUnit *pMapUnit;

  if (pMap->maxUnit <= pMap->count) {
    /* need to allocate a chunk of memories now */
    pMap->maxUnit += DEF_LOCUS_MALLOC_INCREMENT;
    pMap->ppMapUnitList = (MapUnit **) REALLOC ("pMap->ppMapUnitList",
						pMap->ppMapUnitList,
						sizeof (MapUnit *) *
						pMap->maxUnit);
  }

  pMapUnit = (MapUnit *) MALLOC ("pMapUnit", sizeof (MapUnit));
  /* make sure the memory chunk is initialized */
  memset (pMapUnit, 0, sizeof (MapUnit));
  /* initialize the base pair location to negative - unknown */
  pMapUnit->basePairLocation = -1;
  pMap->ppMapUnitList[pMap->count] = pMapUnit;
  pMapUnit->mapIndex = pMap->count;

  /* increment the counter */
  pMap->count++;
  return pMapUnit;
}

void
free_map (Map * pMap)
{
  MapUnit *pMapUnit;
  int i;

  for (i = 0; i < pMap->count; i++) {
    pMapUnit = pMap->ppMapUnitList[i];
    free (pMapUnit);
  }
  free (pMap->ppMapUnitList);
  pMap->count = 0;
  pMap->maxUnit = 0;

  return;
}

/* locate marker map information by its name */
MapUnit *
find_map_unit (Map * pMap, char *sName)
{
  int i;

  for (i = 0; i < pMap->count; i++) {
    if (!strcasecmp (pMap->ppMapUnitList[i]->sName, sName)) {
      /* found it, return the pointer */
      return pMap->ppMapUnitList[i];
    }
  }

  /* exhausted the list, failed to find it */
  return (NULL);
}

/* each marker Locus type has MapUnit pointer saved
 * each MapUnit has marker position information 
 * mapFunctionFlag: 0 - kosambi 1-Haldane
 * distance - cM
 * if either of the loci is trait, then 
 * Return: distance in recombination fraction
 */
double
cm_to_recombination_fraction (double distance, int mapFunctionFlag)
{
  double temp;
  double theta;
  int traitFlag = 0;

  /* convert distance to Morgan */
  distance /= 100;

  if (distance >= 0.0 - ERROR_MARGIN && distance <= ERROR_MARGIN) {
    if (traitFlag)
      return 0;
    else
      return ERROR_MARGIN;
  }

  if (distance < 0)
    distance *= -1;

  /* now convert the map distance to recombination fraction */
  if (mapFunctionFlag == MAP_FUNCTION_KOSAMBI) {
    temp = exp (4 * distance);
    theta = 0.5 * (temp - 1) / (temp + 1);
  } else if (mapFunctionFlag == MAP_FUNCTION_HALDANE) {
    theta = 0.5 * (1 - exp (-2 * distance));
  } else
    /* something is wrong - we can only handle Kosambi or Haldane map functions */
    theta = -1;

  return theta;
}

/* read datafile with specifications for each locus */
int
read_datafile (char *sDatafileName)
{
  FILE *fp;
  int lineNo = 0;
  char line[MAX_LINE_LEN];
  int numRet;
  MapUnit *pMapUnit;
  Locus *pLocus;
  TraitLocus *pTraitLocus = NULL;
  int locusType;		/* temporary place holder */
  char sLocusName[MAX_LINE_LEN];
  char sLocusType[MAX_LINE_LEN];

  fp = fopen (sDatafileName, "r");
  KASSERT (fp != NULL, "Can't open datafile %s for read. Exiting!\n",
	   sDatafileName);

  while (fgets (line, MAX_LINE_LEN, fp)) {
    if (feof (fp))
      break;

    lineNo++;
    /* skip the comment or blank lines */
    if (is_line_blank_or_comment (line))
      continue;

    /* expect to read <LocusName> <LocusType> */
    numRet = sscanf (line, "%s %s", sLocusType, sLocusName);

    KASSERT (numRet == 2, "Can't get locus type or locus name.\n");
    if (!strcasecmp (sLocusType, "M")) {
      locusType = LOCUS_TYPE_MARKER;
      /* marker locus */
      /* allocate space for the marker info first */
      pLocus = add_locus (&originalLocusList, sLocusName, locusType);
      /* locate this marker in the map */
      pMapUnit = find_map_unit (&map, sLocusName);
      KASSERT (pMapUnit != NULL,
	       "Can't find marker %s in map.\n", pLocus->sName);
      pLocus->pMapUnit = pMapUnit;
    } else {
      locusType = LOCUS_TYPE_TRAIT;
      pLocus = add_locus (&originalLocusList, sLocusName, locusType);
      pTraitLocus = pLocus->pTraitLocus;
    }
  }				/* continue to next line */

  fclose (fp);
  return 0;
}				/* end of read_datafile() */

/* read marker frequency file with specifications for each marker */
int
read_markerfile (char *sMarkerfileName, int requiredMarkerCount)
{
  FILE *fp;
  int lineNo = 0;
  char line[MAX_LINE_LEN];
  int pos = 0;
  char *pLine;
  Locus *pLocus;
  int i;
  char sLocusName[MAX_LINE_LEN];
  int found;
  Locus *pTempLocus;
  int markerCount = 0;
  int allele = 1;
  double freq;
  char sAlleleName[MAX_LINE_LEN];


  fp = fopen (sMarkerfileName, "r");
  KASSERT (fp != NULL, "Can't open marker file %s for read. Exiting!\n",
	   sMarkerfileName);

  /* Sample marker file 
     M marker1
     F 0.1 0.2 0.3 0.4
     M marker2
     F 0.6
     F 0.4
     M marker3
     A 140 0.2
     A 141 0.2
     A 150 0.6
     M marker4
     F 0.1 0.2
     F 0.7
   */

  pLocus = NULL;
  while (fgets (line, MAX_LINE_LEN, fp)) {
    lineNo++;
    if (feof (fp))
      break;

    /* now read allele frequencies */
    if (sscanf (line, "M %s", sLocusName) == 1) {
      /* got to a marker line - find the locus in the locus list */
      markerCount++;
      found = FALSE;
      for (i = 0; i < originalLocusList.numLocus; i++) {
	pTempLocus = originalLocusList.ppLocusList[i];
	if (!strcasecmp (pTempLocus->sName, sLocusName)) {
	  /* found it */
	  found = TRUE;
	  pLocus = pTempLocus;
	  allele = 1;
	  break;
	}
      }
      KASSERT (found == TRUE, "Couldn't find marker %s in locus list.\n",
	       sLocusName);
    } else if (sscanf (line, "F %lf %n", &freq, &pos) == 1) {
      /* add the allele */
      sprintf (sAlleleName, "%d", allele);
      add_allele (pLocus, sAlleleName, freq);
      allele++;

      /* there could be more allele frequencies on the line */
      pLine = &line[pos];
      while (sscanf (pLine, "%lf %n", &freq, &pos) == 1) {
	sprintf (sAlleleName, "%d", allele);
	add_allele (pLocus, sAlleleName, freq);
	allele++;
	pLine = pLine + pos;
      }
    } else if (sscanf (line, "A %s %lf", sAlleleName, &freq) == 2) {
      add_allele (pLocus, sAlleleName, freq);
    }
  }				/* continue reading input */

  KASSERT (markerCount >= requiredMarkerCount, "Only %d of %d required markers found in file %s.\n",
	   markerCount, requiredMarkerCount, sMarkerfileName);

  fclose (fp);
  return 0;
}

/* add a locus to a locus list */
Locus *
add_locus (LocusList * pLocusList, char *sName, int locusType)
{
  Locus *pLocus;

  if (pLocusList->maxNumLocus <= pLocusList->numLocus) {
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

  /* if this is a trait locus, need to allocate more space */
  pLocus->locusType = locusType;
  if (locusType == LOCUS_TYPE_TRAIT) {
    pLocus->pTraitLocus = (TraitLocus *) MALLOC ("pLocus->pTraitLocus",
						 sizeof (TraitLocus));
    memset (pLocus->pTraitLocus, 0, sizeof (TraitLocus));
    pLocusList->numTraitLocus++;
  }

  return pLocus;
}

void
free_locus_list (LocusList * pLocusList)
{
  int i, j;
  Locus *pLocus;

  for (i = 0; i < pLocusList->numLocus; i++) {
    pLocus = pLocusList->ppLocusList[i];
    if (modelOptions.polynomial == TRUE) {
      free (pLocus->pAlleleFrequencyPolynomial);
    }
    if (pLocus->locusType == LOCUS_TYPE_TRAIT) {
      free (pLocus->pTraitLocus->pTraits[0]);
      free (pLocus->pTraitLocus);
    }
    free (pLocus->pAlleleFrequency);
    free (pLocus->pAlleleCount);
    for (j = 0; j < pLocus->numAllele; j++) {
      free (pLocus->ppAlleleNames[j]);
    }
    for (j = 0; j < pLocus->numAlleleSet; j++) {
      free_allele_set (pLocus->ppAlleleSetList[j]);
      free (pLocus->ppAlleleSetList[j]);
    }
    free (pLocus->ppAlleleNames);
    free (pLocus->ppAlleleSetList);
    free (pLocus);
  }
  free (pLocusList->pLDLoci);
  free (pLocusList->ppLocusList);
  pLocusList->numLocus = 0;
  pLocusList->maxNumLocus = 0;
}

void
free_allele_set (AlleleSet * pAlleleSet)
{
  free (pAlleleSet->pAlleleBits);
  free (pAlleleSet->pAlleles);
  /* sumFreqPolynomial? */
}

int
find_allele (int locus, char *sAlleleName)
{
  Locus *pLocus;
  int j;

  pLocus = originalLocusList.ppLocusList[locus];
  for (j = 0; j < pLocus->numAllele; j++) {
    if (!strcmp (pLocus->ppAlleleNames[j], sAlleleName))
      /* match found */
      return j + 1;
  }
  return -1;
}

/* add allele into a marker locus */
int
add_allele (Locus * pLocus, char *sAlleleName, double freq)
{
  int numAllele;

  if (pLocus == NULL)
    return -1;

  /* new number of alleles */
  numAllele = pLocus->numAllele + 1;

  /* reallocate alelel frequency space */
  pLocus->pAlleleFrequency = (double *) REALLOC ("pLocus->pAlleleFrequency",
						 pLocus->pAlleleFrequency,
						 numAllele * sizeof (double));

  /* if we are doing set recoding, we need to calculate the length of 
   * the allele set - how many integers we need to represent a allele in
   * bit format - be generous, as super alleles increases the number of allele sets we need */
  pLocus->alleleSetLen = numAllele / INT_BITS + 2;
  if (originalLocusList.alleleSetLen < pLocus->alleleSetLen)
    originalLocusList.alleleSetLen = pLocus->alleleSetLen;

  /* allocate space for frequency and count */
  if (modelOptions.polynomial == TRUE) {
    pLocus->pAlleleFrequencyPolynomial =
      (Polynomial *) REALLOC ("pLocus->pAlleleFrequencyPolynomial",
			      pLocus->pAlleleFrequencyPolynomial,
			      numAllele * sizeof (Polynomial));
  }
  /* actual count of the alleles in the pedigree */
  pLocus->pAlleleCount = (short *) REALLOC ("pLocus->pAlleleCount",
					    pLocus->pAlleleCount,
					    numAllele * sizeof (short));

  /* original allele names */
  pLocus->ppAlleleNames = (char **) REALLOC ("pLocus->ppAlleleNames",
					     pLocus->ppAlleleNames,
					     numAllele * sizeof (char *));
  /* add the name and frequency in */
  pLocus->ppAlleleNames[numAllele - 1] =
    (char *) MALLOC ("pLocus->ppAlleleNames[*]",
		     (strlen (sAlleleName) + 1) * sizeof (char));
  strcpy (pLocus->ppAlleleNames[numAllele - 1], sAlleleName);
  pLocus->pAlleleFrequency[numAllele - 1] = freq;
  pLocus->pAlleleCount[numAllele - 1] = 0;
  /* update allele counts */
  pLocus->numAllele = numAllele;
  pLocus->numOriginalAllele = numAllele;

  return 0;
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
  for (i = 0; i < pPedigree->numPerson; i++) {
    pPerson = pPedigree->ppPersonList[i];

    /* ignore the loop breaker duplicate */
    if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] == NULL)
      continue;

    /* base on phenotype to make up the list */
    if (pPerson->pTypedFlag[locus] == 1) {
      if (pPerson->pPhasedFlag[locus] == 1 ||
	  pPerson->pPhenotypeList[MOM][locus] ==
	  pPerson->pPhenotypeList[DAD][locus]) {
	/* only 1 phased genotype for phased or homozygote */
	pGenotype = add_genotype (&(pPerson->ppGenotypeList[locus]),
				  &pPerson->pNumGenotype[locus],
				  locus,
				  pPerson->pPhenotypeList[DAD][locus],
				  pPerson->pPhenotypeList[MOM][locus]);
      } else {
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
    } else {
      if (pPerson->pPhenotypeList[0][locus] == 0 &&
	  pPerson->pPhenotypeList[1][locus] == 0) {
	/* marker phenotype is not known, add all possible combinations */
	allele1 = 1;
	allele2 = 1;
	for (allele1 = 1; allele1 <= pLocus->numAllele; allele1++) {
	  for (allele2 = allele1; allele2 <= pLocus->numAllele; allele2++) {
	    /* for X chromosome, only add homozygous genotypes for MALE */
	    if (modelOptions.sexLinked && pPerson->sex + 1 == MALE
		&& allele1 != allele2)
	      continue;
	    pGenotype =
	      add_genotype (&pPerson->ppGenotypeList[locus],
			    &pPerson->pNumGenotype[locus], locus,
			    allele1, allele2);
	    if (allele2 != allele1) {
	      pGenotype2 =
		add_genotype (&pPerson->ppGenotypeList[locus],
			      &pPerson->pNumGenotype[locus],
			      locus, allele2, allele1);
	      pGenotype->pDualGenotype = pGenotype2;
	      pGenotype2->pDualGenotype = pGenotype;
	    }
	  }			/* end of allele2 loop */
	}			/* end of allele1 loop */
      } else {			/* must be half typed */

	if (pPerson->pPhenotypeList[0][locus] != 0)
	  allele1 = pPerson->pPhenotypeList[0][locus];
	else
	  allele1 = pPerson->pPhenotypeList[1][locus];
	for (allele2 = 1; allele2 <= pLocus->numAllele; allele2++) {
	  /* for X chromosome, only add homozygous genotypes for MALE */
	  if (modelOptions.sexLinked && pPerson->sex + 1 == MALE &&
	      allele1 != allele2)
	    continue;
	  pGenotype = add_genotype (&pPerson->ppGenotypeList[locus],
				    &pPerson->pNumGenotype[locus],
				    locus, allele1, allele2);
	  if (allele2 != allele1) {
	    pGenotype2 =
	      add_genotype (&pPerson->ppGenotypeList[locus],
			    &pPerson->pNumGenotype[locus], locus,
			    allele2, allele1);
	    pGenotype->pDualGenotype = pGenotype2;
	    pGenotype2->pDualGenotype = pGenotype;
	  }
	}			/* end of allele2 loop */
      }				/* half typed */
    }				/* not typed or not fully typed */
  }				/* loop of persons in a pedigree */

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
  double pen = 1.0;

  Polynomial *penPolynomial;

  /* go through everyone in the pedigree */
  for (i = 0; i < pPedigree->numPerson; i++) {
    pPerson = pPedigree->ppPersonList[i];

    /* ignore the loop breaker duplicate */
    if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] == NULL)
      continue;

    /* go through all possible genotype combinations */
    allele1 = 1;
    allele2 = 1;
    for (allele1 = 1; allele1 <= pLocus->numAllele; allele1++) {
      for (allele2 = allele1; allele2 <= pLocus->numAllele; allele2++) {
	/* for X chromosome, only add homozygous genotypes for MALE */
	if (modelOptions.sexLinked && pPerson->sex + 1 == MALE &&
	    allele1 != allele2)
	  continue;

	if (modelOptions.polynomial == TRUE) {
	  compute_penetrance (pPerson, locus, allele1, allele2,
			      &penPolynomial);
	}
	//              else
	//        compute_penetrance (pPerson, locus, allele1, allele2, &pen);
	/* if this genotype hasn't been rejected, then add it 
	 * under polynomial, penetrance is a variable (parameter), we
	 * can't check against one single value, so we include all
	 * trait genotypes */
	if (modelOptions.polynomial == TRUE) {
	  pGenotype = add_genotype (&pPerson->ppGenotypeList[locus],
				    &pPerson->pNumGenotype[locus],
				    locus, allele1, allele2);
	  pGenotype->penslot.penetrancePolynomial = penPolynomial;
	  if (allele1 != allele2) {
	    pGenotype2 =
	      add_genotype (&pPerson->ppGenotypeList[locus],
			    &pPerson->pNumGenotype[locus], locus,
			    allele2, allele1);
	    pGenotype->pDualGenotype = pGenotype2;
	    pGenotype2->pDualGenotype = pGenotype;
	    pGenotype2->penslot.penetrancePolynomial = penPolynomial;
	  }
	} else {
	  /* we can't add genotype base on pen here, as for the first input pen
	   * vector, some genotype might not be compatible, but as we move
	   * to next pen vector, it may be very well compatible */
	  //              if (pen > 0)
	  {
	    pGenotype =
	      add_genotype (&pPerson->ppGenotypeList[locus],
			    &pPerson->pNumGenotype[locus], locus,
			    allele1, allele2);
	    pGenotype->penslot.penetrance = pen;
	    if (allele1 != allele2) {
	      pGenotype2 =
		add_genotype (&pPerson->ppGenotypeList[locus],
			      &pPerson->pNumGenotype[locus],
			      locus, allele2, allele1);
	      pGenotype->pDualGenotype = pGenotype2;
	      pGenotype2->pDualGenotype = pGenotype;
	      pGenotype2->penslot.penetrance = pen;
	    }
	  }
	}
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
  double trait, mean, stddev, temp;
  double df;

  Polynomial *tempPoly = NULL, *tempPoly2 = NULL;

  tempPoly = NULL;
  tempPoly2 = NULL;

  if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] == NULL)
    pPerson = pPerson->pOriginalPerson;

  /* now go through all related traits and see whether 
   * this genotype is compatible with each trait value if
   * the trait value is known */
  /* for now, we assume we either have only 1 affection status trait
   * or 1 or more quantitative traits only 
   * i.e. there is no mixing between affection status trait and 
   * liability trait */
  i = 0;

  while (i < pTraitLocus->numTrait) {
    pTrait = pTraitLocus->pTraits[i];
    if (pTrait->numLiabilityClass > 1)
      liabilityClass = pPerson->ppLiabilityClass[locus][i];
    else
      liabilityClass = 1;

    if (pTrait->type == DICHOTOMOUS) {
      if (pPerson->ppTraitKnown[locus][i] == TRUE) {
	affectionStatus = (int) pPerson->ppTraitValue[locus][i];

	if (modelOptions.polynomial == TRUE) {
	  if (affectionStatus == AFFECTION_STATUS_AFFECTED) {
	    /* build the penetrance polynomial - single variable */
	    char vName[100];

	    sprintf (vName, "p[%d][%d][%d][%d]",
		     AFFECTION_STATUS_AFFECTED, liabilityClass - 1,
		     allele1 - 1, allele2 - 1);

	    *(Polynomial **) pen =
	      variableExp (&pTrait->penetrance[AFFECTION_STATUS_AFFECTED]
			   [liabilityClass - 1][allele1 -
						1][allele2 - 1],
			   NULL, 'D', "");
	  } else {
	    /* unaffected: 1 - pen[affected] 
	     * 1 - pTrait->penetrance[AFFECTION_STATUS_AFFECTED][liabilityClass - 1]
	     *                                 [allele1-1][allele2-1] */
	    char vName[100];

	    sprintf (vName, "p[%d][%d][%d][%d]",
		     AFFECTION_STATUS_AFFECTED, liabilityClass - 1,
		     allele1 - 1, allele2 - 1);
	    *(Polynomial **) pen =
	      plusExp (2, 1.0, constant1Poly, -1.0,
		       variableExp (&pTrait->
				    penetrance[AFFECTION_STATUS_AFFECTED]
				    [liabilityClass - 1][allele1 -
							 1][allele2
							    - 1],
				    NULL, 'D', vName), 0);
	  }
	} else {
	  /* polynomial flag is not turned on */
	  *(double *) pen =
	    pTrait->penetrance[affectionStatus][liabilityClass - 1]
	    [allele1 - 1][allele2 - 1];
//          fprintf(stderr,"AffectionStatus=%d liabilityClass-1=%d allele1-1=%d allele2-1=%d pen=%f\n",
//                        affectionStatus,liabilityClass - 1,allele1-1,allele2-1,
//                        pTrait->penetrance[AFFECTION_STATUS_AFFECTED][liabilityClass - 1]
//                                          [allele1-1][allele2-1]);

	}
      } else {
	affectionStatus = AFFECTION_STATUS_UNKNOWN;
	/* when the affection status is unknown, the penetrance is set to 1 */

	if (modelOptions.polynomial == TRUE) {
	  *(Polynomial **) pen = constant1Poly;
	} else
	  *(double *) pen = 1;
      }
    } /* end of affection status trait handling */
    else if (pTrait->type == QUANTITATIVE || pTrait->type == COMBINED) {

      /* make sure we have all the trait values for this locus */
      for (j = 0; j < pTraitLocus->numTrait; j++) {
	/* if the quantitative trait value is not known for any trait, then
	 * we can't calculate "exact" penetrance, just return 1 */
	if (pPerson->ppTraitKnown[locus][j] == FALSE) {
	  if (modelOptions.polynomial == TRUE) {
	    *(Polynomial **) pen = constant1Poly;
	  } else
	    *(double *) pen = 1;
	  return;
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
      stddev = pTrait->stddev[liabilityClass - 1][allele1 - 1][allele2 - 1];


      /* working on the assumption of just 1 QT trait */
      if (pTrait->functionQT == QT_FUNCTION_NORMAL) {
	if (pTrait->type == COMBINED &&
	    (trait >= pTrait->lessCutoffFlag - 0.000001 &&
	     trait <= pTrait->lessCutoffFlag + 0.000001)) {
	  /* the cutoff is number of SD, so it's already been converted to X in N(0,1) */

	  if (modelOptions.polynomial == TRUE) {
	    tempPoly = plusExp (2, 1.0,
				variableExp (&pTrait->
					     cutoffValue
					     [liabilityClass - 1],
					     NULL, 'D', "cutoff"),
				-1.0,
				variableExp (&pTrait->
					     means[liabilityClass -
						   1][allele1 -
						      1][allele2 -
							 1], NULL,
					     'D', "mean"), 0);
	    tempPoly = timesExp (2, tempPoly, 1,
				 variableExp (&pTrait->
					      stddev[liabilityClass
						     - 1][allele1 -
							  1][allele2
							     - 1],
					      NULL, 'D', "stddev"), -1, 1);
	    *(Polynomial **) pen =
	      functionCallExp (2, "gsl_cdf_ugaussian_P", tempPoly);

	  } else {
	    temp = pTrait->cutoffValue[liabilityClass - 1] - mean;
	    /* if standard deviation is 1, don't do the division calculation */
	    //                      if (0.999999 > stddev || stddev > 1.000001)
	    temp = temp / stddev;
	    *(double *) pen = gsl_cdf_ugaussian_P (temp);
	  }
	} else if (pTrait->type == COMBINED &&
		   (trait >= pTrait->moreCutoffFlag - 0.000001 &&
		    trait <= pTrait->moreCutoffFlag + 0.000001)) {
	  /* the cutoff is number of SD, so it's already been converted to X in N(0,1) */
	  if (modelOptions.polynomial == TRUE) {
	    tempPoly = plusExp (2, 1.0,
				variableExp (&pTrait->
					     cutoffValue
					     [liabilityClass - 1],
					     NULL, 'D', "cutoff"),
				-1.0,
				variableExp (&pTrait->
					     means[liabilityClass -
						   1][allele1 -
						      1][allele2 -
							 1], NULL,
					     'D', "mean"), 0);
	    tempPoly = timesExp (2, tempPoly, 1,
				 variableExp (&pTrait->
					      stddev[liabilityClass
						     - 1][allele1 -
							  1][allele2
							     - 1],
					      NULL, 'D', "stddev"), -1, 1);

	    *(Polynomial **) pen =
	      functionCallExp (2, "gsl_cdf_ugaussian_Q", tempPoly);

	  } else {
	    temp = pTrait->cutoffValue[liabilityClass - 1] - mean;
	    /* if standard deviation is 1, don't do the division calculation */
	    //  if (0.999999 > stddev || stddev > 1.000001)
	    temp /= stddev;
	    *(double *) pen = gsl_cdf_ugaussian_Q (temp);
	  }
	} else {		/* point PDF */
	  if (modelOptions.polynomial == TRUE) {
	    tempPoly =
	      timesExp (2,
			plusExp (2, 1.0,
				 variableExp (&pPerson->ppTraitValue[locus]
					      [i], NULL, 'D',
					      "trait"), -1.0,
				 variableExp (&pTrait->
					      means[liabilityClass -
						    1][allele1 -
						       1][allele2 -
							  1], NULL,
					      'D', "mean"), 0), 1,
			variableExp (&pTrait->
				     stddev[liabilityClass -
					    1][allele1 -
					       1][allele2 - 1],
				     NULL, 'D', "stddev"), -1, 0);

	    /* if lower bound exists and the point is right at the lower bound
	     * use CDF instead of pdf */
	    if (pTrait->minFlag && trait >= pTrait->min - 0.000001
		&& trait <= pTrait->min + 0.000001) {
	      *(Polynomial **) pen =
		functionCallExp (2, "gsl_cdf_ugaussian_P", tempPoly);
	    } else
	      if (pTrait->maxFlag && trait >= pTrait->max - 0.000001
		  && trait <= pTrait->max + 0.000001) {
	      *(Polynomial **) pen =
		functionCallExp (2, "gsl_cdf_ugaussian_Q", tempPoly);
	    } else {
	      *(Polynomial **) pen =
		functionCallExp (2, "gsl_ran_ugaussian_pdf", tempPoly);
	    }

	  } else {
	    temp = trait - mean;
	    /* if standard deviation is 1, don't do the division calculation */
	    //  if (0.999999 > stddev || stddev > 1.000001)
	    temp /= stddev;

	    /* if lower bound exists and the point is right at the lower bound
	     * use CDF instead of pdf */
	    if (pTrait->minFlag && trait >= pTrait->min - 0.000001
		&& trait <= pTrait->min + 0.000001) {
	      *(double *) pen = gsl_cdf_ugaussian_P (temp);
	    } else
	      if (pTrait->maxFlag && trait >= pTrait->max - 0.000001
		  && trait <= pTrait->max + 0.000001) {
	      *(double *) pen = gsl_cdf_ugaussian_Q (temp);
	    } else {
	      *(double *) pen = gsl_ran_ugaussian_pdf (temp);
	    }
	  }
	}
      } else if (pTrait->functionQT == QT_FUNCTION_T) {
	df = pTrait->dfQT;
	if (pTrait->type == COMBINED &&
	    (trait >= pTrait->lessCutoffFlag - 0.000001 &&
	     trait <= pTrait->lessCutoffFlag + 0.000001)) {

	  if (modelOptions.polynomial == TRUE) {
	    tempPoly = plusExp (2, 1.0,
				variableExp (&pTrait->
					     cutoffValue
					     [liabilityClass - 1],
					     NULL, 'D', "cutoff"),
				-1.0,
				variableExp (&pTrait->
					     means[liabilityClass -
						   1][allele1 -
						      1][allele2 -
							 1], NULL,
					     'D', "mean"), 0);
	    tempPoly = timesExp (2, tempPoly, 1,
				 timesExp (2,
					   variableExp (&pTrait->
							stddev
							[liabilityClass
							 - 1][allele1 - 1]
							[allele2 -
							 1], NULL,
							'D',
							"stddev"),
					   1, functionCallExp (2,
							       "sqrt",
							       timesExp
							       (2,
								plusExp
								(2,
								 1.0,
								 variableExp
								 (&pTrait->
								  dfQT,
								  NULL,
								  'D',
								  "df"),
								 -1.0,
								 constantExp
								 (2.0),
								 0),
								1,
								variableExp
								(&pTrait->
								 dfQT,
								 NULL,
								 'D',
								 "df"),
								-1,
								0)),
					   1, 0), -1, 1);

	    *(Polynomial **) pen =
	      functionCallExp (3, "gsl_cdf_tdist_P", tempPoly,
			       variableExp (&pTrait->dfQT, NULL, 'D', "df"));
	  } else {
	    temp = pTrait->cutoffValue[liabilityClass - 1] - mean;
	    temp /= stddev * sqrt ((df - 2) / df);
	    *(double *) pen = gsl_cdf_tdist_P (temp, df);
	  }
	} else if (pTrait->type == COMBINED &&
		   (trait >= pTrait->moreCutoffFlag - 0.000001 &&
		    trait <= pTrait->moreCutoffFlag + 0.000001)) {
	  if (modelOptions.polynomial == TRUE) {

	    tempPoly = plusExp (2, 1.0, variableExp (&pTrait->
						     cutoffValue
						     [liabilityClass
						      - 1], NULL,
						     'D', "cutoff"),
				-1.0,
				variableExp (&pTrait->
					     means[liabilityClass -
						   1][allele1 -
						      1][allele2 -
							 1], NULL,
					     'D', "mean"), 0);


	    tempPoly = timesExp (2, tempPoly, 1,
				 timesExp (2,
					   variableExp (&pTrait->
							stddev
							[liabilityClass
							 - 1][allele1 - 1]
							[allele2 -
							 1], NULL,
							'D',
							"stddev"),
					   1, functionCallExp (2,
							       "sqrt",
							       timesExp
							       (2,
								plusExp
								(2,
								 1.0,
								 variableExp
								 (&pTrait->
								  dfQT,
								  NULL,
								  'D',
								  "df"),
								 -1.0,
								 constantExp
								 (2.0),
								 0),
								1,
								variableExp
								(&pTrait->
								 dfQT,
								 NULL,
								 'D',
								 "df"),
								-1,
								0)),
					   1, 0), -1, 1);

	    *(Polynomial **) pen =
	      functionCallExp (3, "gsl_cdf_tdist_Q", tempPoly,
			       variableExp (&pTrait->dfQT, NULL, 'D', "df"));
	  } else {
	    temp = pTrait->cutoffValue[liabilityClass - 1] - mean;
	    temp /= stddev * sqrt ((df - 2) / df);
	    *(double *) pen = gsl_cdf_tdist_Q (temp, df);
	  }

	} else {		/* point pdf */
	  if (modelOptions.polynomial == TRUE) {
	    tempPoly = plusExp (2, 1.0,
				variableExp (&pPerson->ppTraitValue[locus]
					     [i], NULL, 'D',
					     "trait"),
				-1.0,
				variableExp (&pTrait->
					     means[liabilityClass -
						   1][allele1 -
						      1][allele2 -
							 1], NULL,
					     'D', "mean"), 0);

	    tempPoly = timesExp (2, tempPoly, 1,
				 timesExp (2,
					   variableExp (&pTrait->
							stddev
							[liabilityClass
							 - 1][allele1 - 1]
							[allele2 -
							 1], NULL,
							'D',
							"stddev"),
					   1, functionCallExp (2,
							       "sqrt",
							       timesExp
							       (2,
								plusExp
								(2,
								 1.0,
								 variableExp
								 (&pTrait->
								  dfQT,
								  NULL,
								  'D',
								  "df"),
								 -1.0,
								 constantExp
								 (2.0),
								 0),
								1,
								variableExp
								(&pTrait->
								 dfQT,
								 NULL,
								 'D',
								 "df"),
								-1,
								0)),
					   1, 0), -1, 1);

	    if (pTrait->minFlag && trait >= pTrait->min - 0.000001
		&& trait <= pTrait->min + 0.000001) {
	      *(Polynomial **) pen =
		functionCallExp (3, "gsl_cdf_tdist_P", tempPoly,
				 variableExp (&pTrait->dfQT, NULL,
					      'D', "df"));
	    } else
	      if (pTrait->maxFlag && trait >= pTrait->max - 0.000001
		  && trait <= pTrait->max + 0.000001) {
	      *(Polynomial **) pen =
		functionCallExp (3, "gsl_cdf_tdist_Q", tempPoly,
				 variableExp (&pTrait->dfQT, NULL,
					      'D', "df"));
	    } else {
	      *(Polynomial **) pen =
		functionCallExp (3, "gsl_ran_tdist_pdf", tempPoly,
				 variableExp (&pTrait->dfQT,
					      NULL, 'D', "df"));
	    }

	  } else {
	    temp = trait - mean;
	    temp /= stddev * sqrt ((df - 2) / df);

	    /* if lower bound exists and the point is right at the lower bound
	     * use CDF instead of pdf */
	    if (pTrait->minFlag && trait >= pTrait->min - 0.000001
		&& trait <= pTrait->min + 0.000001) {
	      *(double *) pen = gsl_cdf_tdist_P (temp, df);
	    } else
	      if (pTrait->maxFlag && trait >= pTrait->max - 0.000001
		  && trait <= pTrait->max + 0.000001) {
	      *(double *) pen = gsl_cdf_tdist_Q (temp, df);
	    } else {
	      *(double *) pen = gsl_ran_tdist_pdf (temp, df);
	    }
	  }
	}
      } /* t-dsitribution */
      else if (pTrait->functionQT == QT_FUNCTION_CHI_SQUARE) {
	if (pTrait->type == COMBINED &&
	    (trait >= pTrait->lessCutoffFlag - 0.000001 &&
	     trait <= pTrait->lessCutoffFlag + 0.000001)) {
	  if (modelOptions.polynomial == TRUE) {
	    *(Polynomial **) pen =
	      functionCallExp (3, "gsl_cdf_chisq_P",
			       variableExp (&pTrait->
					    cutoffValue
					    [liabilityClass - 1],
					    NULL, 'D', "cutoff"),
			       variableExp (&pTrait->
					    means[liabilityClass -
						  1][allele1 -
						     1][allele2 -
							1], NULL,
					    'D', "mean"));
	  } else {
	    /* mean is used to save the degree of freedom */
	    *(double *) pen =
	      gsl_cdf_chisq_P (pTrait->cutoffValue[liabilityClass - 1], mean);
	  }
	} else if (pTrait->type == COMBINED &&
		   (trait >= pTrait->moreCutoffFlag - 0.000001 &&
		    trait <= pTrait->moreCutoffFlag + 0.000001)) {
	  /* the cutoff is number of SD, so it's already been converted to X in N(0,1) */
	  if (modelOptions.polynomial == TRUE) {
	    *(Polynomial **) pen =
	      functionCallExp (3, "gsl_cdf_chisq_Q",
			       variableExp (&pTrait->
					    cutoffValue
					    [liabilityClass - 1],
					    NULL, 'D', "cutoff"),
			       variableExp (&pTrait->
					    means[liabilityClass -
						  1][allele1 -
						     1][allele2 -
							1], NULL,
					    'D', "mean"));
	  } else {
	    *(double *) pen =
	      gsl_cdf_chisq_Q (pTrait->cutoffValue[liabilityClass - 1], mean);
	  }
	} else {		/* point PDF */
	  if (modelOptions.polynomial == TRUE) {
	    if (pTrait->minFlag && trait >= pTrait->min - 0.000001
		&& trait <= pTrait->min + 0.000001) {
	      *(Polynomial **) pen =
		functionCallExp (3, "gsl_cdf_chisq_P",
				 variableExp (&pPerson->ppTraitValue[locus]
					      [i], NULL, 'D',
					      "trait"),
				 variableExp (&pTrait->
					      means[liabilityClass
						    - 1][allele1 - 1]
					      [allele2 - 1], NULL,
					      'D', "mean"));
	    } else
	      if (pTrait->maxFlag && trait >= pTrait->max - 0.000001
		  && trait <= pTrait->max + 0.000001) {
	      *(Polynomial **) pen =
		functionCallExp (3, "gsl_cdf_chisq_Q",
				 variableExp (&pPerson->ppTraitValue[locus]
					      [i], NULL, 'D',
					      "trait"),
				 variableExp (&pTrait->
					      means[liabilityClass
						    - 1][allele1 - 1]
					      [allele2 - 1], NULL,
					      'D', "mean"));
	    } else {
	      *(Polynomial **) pen =
		functionCallExp (3, "gsl_ran_chisq_pdf",
				 variableExp (&pPerson->ppTraitValue[locus]
					      [i], NULL, 'D',
					      "trait"),
				 variableExp (&pTrait->
					      means[liabilityClass
						    - 1][allele1 - 1]
					      [allele2 - 1], NULL,
					      'D', "mean"));

	    }
	  } else {
	    /* if lower bound exists and the point is right at the lower bound
	     * use CDF instead of pdf */
	    /*
	       if(pTrait->minFlag && trait >= pTrait->min - 0.000001 && 
	       trait <= pTrait->min + 0.000001)
	       {
	       *(double *) pen = gsl_cdf_chisq_P (trait, mean);
	       }
	       else 
	       if(pTrait->maxFlag && trait >= pTrait->max - 0.000001 && 
	       trait <= pTrait->max + 0.000001)
	       {
	       *(double *) pen = gsl_cdf_chisq_Q (trait, mean);
	       }
	       else
	     */
	    {
	      *(double *) pen = gsl_ran_chisq_pdf (trait, mean);
	    }
	  }
	}
      }
      /* end of chi square distribution */
    }				/* end of quantitative trait handling */
    i++;
  }				/* move to next trait for the same locus */

  return;
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
  pGenotype->penslot.penetrance = 1;

  /* add this to the top of the genotype list */
  pGenotype->pNext = *ppList;
  *ppList = pGenotype;
  if (modelOptions.polynomial == TRUE) {
    pGenotype->penslot.penetrancePolynomial = constant1Poly;
  }

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
  pGenotype->pAlleleBits[DAD] =
    (unsigned int *) MALLOC ("pGenotype->pAlleleBits[DAD]",
			     sizeof (unsigned int) * numInts);
  pGenotype->pAlleleBits[MOM] =
    (unsigned int *) MALLOC ("pGenotype->pAlleleBits[MOM]",
			     sizeof (unsigned int) * numInts);
  memset (pGenotype->pAlleleBits[DAD], 0, sizeof (unsigned int) * numInts);
  memset (pGenotype->pAlleleBits[MOM], 0, sizeof (unsigned int) * numInts);
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
  if (*pHead == pGenotype) {
    *pHead = pGenotype->pNext;
  } else {
    pCurr = (*pHead)->pNext;
    pPrev = (*pHead);
    while (pCurr != NULL && pCurr != pGenotype) {
      pPrev = pCurr;
      pCurr = pCurr->pNext;
    }
    /* failed to find the genotype */
    if (pCurr == NULL)
      return -1;
    pPrev->pNext = pGenotype->pNext;
  }

  /* free the space */
  if (removeGenotypeFlag == TRUE) {
    free (pGenotype->pAlleleBits[DAD]);
    free (pGenotype->pAlleleBits[MOM]);
    free (pGenotype);
  }

  /* decrement the counter */
  (*pCount)--;

  return 0;
}


/* This genotype list should be phased - for printing purpose, 
 * it doesn't matter */
void
print_person_locus_genotype_list (Person * pPerson, int locus)
{
#if DEBUG
  Genotype *pGenotype = pPerson->ppGenotypeList[locus];

  /* print out person lable */
  logMsg (LOGGENOELIM, LOGDEBUG,
	  "    Person %s locus %d num of geno %d: \n\t", pPerson->sID, locus,
	  pPerson->pNumGenotype[locus]);
  while (pGenotype != NULL) {
    logMsg (LOGGENOELIM, LOGDEBUG, "(%d,%d) ", pGenotype->allele[DAD],
	    pGenotype->allele[MOM]);
    pGenotype = pGenotype->pNext;
  }
  logMsg (LOGGENOELIM, LOGDEBUG, "\n");

#endif
  return;
}

void
print_pedigree_locus_genotype_list (Pedigree * pPedigree, int locus)
{
  int i;

  KLOG (LOGGENOELIM, LOGDEBUG, "Pedigree %3s Locus %d: \n",
	pPedigree->sPedigreeID, locus);
  for (i = 0; i < pPedigree->numPerson; i++) {
    print_person_locus_genotype_list (pPedigree->ppPersonList[i], locus);
  }
  return;
}

void
print_pedigree_locus_genotype_count (Pedigree * pPedigree, int locus)
{
  int i;
  Person *pPerson;

  fprintf (stderr, "Pedigree %3s:\n", pPedigree->sPedigreeID);
  for (i = 0; i < pPedigree->numPerson; i++) {
    pPerson = pPedigree->ppPersonList[i];
    fprintf (stderr, "    Person %s has %d genotypes.\n",
	     pPerson->sID, pPerson->pNumGenotype[locus]);
  }
  return;
}

/* This procedure set the genotype weight of each genotype for each
 * person in a pedigree 
 * for a homozygous genotype, the weight is p*p
 * for a phased heterozygous genotype, the weight is pq
 * This will facilitate likelihood calculation
 * locus - original locus index
 * 
 * */
int
set_genotype_weight (Pedigree * pPedigree, int locus)
{
  Person *pPerson;
  Genotype *pGenotype;
  int i;
  Locus *pLocus;
  double alleleFreq[2] = { 0, 0 };
  Polynomial *alleleFreqPolynomial[2];

  pLocus = originalLocusList.ppLocusList[locus];
  for (i = 0; i < pPedigree->numPerson; i++) {
    pPerson = pPedigree->ppPersonList[i];
    memcpy (&pPerson->pSavedNumGenotype[0], &pPerson->pNumGenotype[0],
	    sizeof (int) * originalLocusList.numLocus);
    pGenotype = pPerson->ppGenotypeList[locus];
    while (pGenotype) {
      if (modelOptions.polynomial == TRUE) {
	if (pGenotype->allele[DAD] <= pLocus->numAllele) {
	  char vName[100];

	  sprintf (vName, "f[%d][%d]", locus, pGenotype->allele[DAD] - 1);

	  if (pLocus->locusType == LOCUS_TYPE_MARKER)
	    alleleFreqPolynomial[DAD] =
	      constantExp (pLocus->
			   pAlleleFrequency[pGenotype->allele[DAD] - 1]);
	  else
	    alleleFreqPolynomial[DAD] =
	      variableExp (&pLocus->
			   pAlleleFrequency[pGenotype->allele[DAD] -
					    1], NULL, 'D', vName);
	} else {
	  alleleFreqPolynomial[DAD] =
	    pLocus->ppAlleleSetList[pGenotype->allele[DAD] -
				    1]->sumFreqPolynomial;
	}
	if (pGenotype->allele[MOM] <= pLocus->numAllele) {
	  char vName[100];

	  sprintf (vName, "f[%d][%d]", locus, pGenotype->allele[MOM] - 1);
	  if (pLocus->locusType == LOCUS_TYPE_MARKER)
	    alleleFreqPolynomial[MOM] =
	      constantExp (pLocus->
			   pAlleleFrequency[pGenotype->allele[MOM] - 1]);
	  else
	    alleleFreqPolynomial[MOM] =
	      variableExp (&pLocus->
			   pAlleleFrequency[pGenotype->allele[MOM] -
					    1], NULL, 'D', vName);
	} else
	  alleleFreqPolynomial[MOM] =
	    pLocus->ppAlleleSetList[pGenotype->allele[MOM] -
				    1]->sumFreqPolynomial;
      } else {
	if (pGenotype->allele[DAD] <= pLocus->numAllele) {
	  alleleFreq[DAD] =
	    pLocus->pAlleleFrequency[pGenotype->allele[DAD] - 1];
	} else
	  alleleFreq[DAD] =
	    pLocus->ppAlleleSetList[pGenotype->allele[DAD] - 1]->sumFreq;
	if (pGenotype->allele[MOM] <= pLocus->numAllele) {
	  alleleFreq[MOM] =
	    pLocus->pAlleleFrequency[pGenotype->allele[MOM] - 1];
	} else
	  alleleFreq[MOM] =
	    pLocus->ppAlleleSetList[pGenotype->allele[MOM] - 1]->sumFreq;
      }
      if (modelOptions.polynomial == TRUE) {
	if (modelOptions.sexLinked && pPerson->sex + 1 == MALE) {
	  pGenotype->wtslot.weightPolynomial = alleleFreqPolynomial[DAD];
	} else
	  /* build the polynomial sum */
	  pGenotype->wtslot.weightPolynomial =
	    timesExp (2, alleleFreqPolynomial[DAD], 1,
		      alleleFreqPolynomial[MOM], 1, 0);
      } else {
	if (modelOptions.sexLinked && pPerson->sex + 1 == MALE) {
	  pGenotype->wtslot.weight = alleleFreq[DAD];
	} else
	  pGenotype->wtslot.weight = alleleFreq[DAD] * alleleFreq[MOM];
      }
      pGenotype = pGenotype->pNext;
    }
  }				/* end of looping over persons in the given pedigree */

  return 0;
}

/* This procedure set the position of each genotype in this person's
 * genotype list for each person in a pedigree 
 * This will facilitate storing and retrieving multi locus 
 * conditional genotype likelihood
 * This function should be called after genotype elimination 
 * 
 * */
int
set_genotype_position (Pedigree * pPedigree, int locus)
{
  Person *pPerson;
  Genotype *pGenotype;
  int i;
  int position;
  Locus *pLocus;

  pLocus = originalLocusList.ppLocusList[locus];
  for (i = 0; i < pPedigree->numPerson; i++) {
    pPerson = pPedigree->ppPersonList[i];
    pGenotype = pPerson->ppGenotypeList[locus];
    position = 0;
    while (pGenotype) {
      pGenotype->position = position;
      pGenotype = pGenotype->pNext;
      position++;
    }
  }

  return 0;
}

/* Only allocated once */
int
allocate_multi_locus_genotype_storage (Pedigree * pPedigree, int numLocus)
{
  int locus;
  Person *pPerson;
  int i, j, k;
  int size;
  int *sortedList;
  int numGeno;

  /* sorted from max to min */
  sortedList = malloc (originalLocusList.numLocus * sizeof (int));
  memset (sortedList, 0, originalLocusList.numLocus * sizeof (int));

  for (i = 0; i < pPedigree->numPerson; i++) {
    pPerson = pPedigree->ppPersonList[i];
    pPerson->multiLocusAdjust = (int *) malloc (numLocus * sizeof (int));
    pPerson->numSavedGenotype2 = (int *) malloc (numLocus * sizeof (int));
    for (locus = 0; locus < originalLocusList.numLocus; locus++) {
      numGeno = pPerson->pSavedNumGenotype[locus];
      for (j = 0; j < locus; j++) {
	if (numGeno > sortedList[j])
	  break;
      }
      for (k = locus - 1; k >= j; k--) {
	sortedList[k + 1] = sortedList[k];
      }
      sortedList[j] = numGeno;
    }
    size = 1;
    for (locus = 0; locus < numLocus && locus < originalLocusList.numLocus;
	 locus++) {
      size *= sortedList[locus];
    }

    /* allocate space */
    pPerson->pLikelihood =
      (ConditionalLikelihood *) MALLOC ("pPerson->pLikelihood ",
					sizeof (ConditionalLikelihood) *
					size);
    //memset (pPerson->pLikelihood, 0, sizeof (ConditionalLikelihood) * size);
    pPerson->maxNumConditionals = size;
    pPerson->pTmpLikelihoodIndex = (int *) malloc (sizeof (int) * size);

    /* allocate loop breaker work space */
    if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] != NULL) {
      pPerson->loopBreakerStruct =
	(LoopBreaker *) calloc (1, sizeof (LoopBreaker));
      pPerson->loopBreakerStruct->maxNumGenotype = size;
      pPerson->loopBreakerStruct->genotype =
	(Genotype ***) malloc (sizeof (Genotype **) * size);
      for (j = 0; j < size; j++) {
	pPerson->loopBreakerStruct->genotype[j] =
	  (Genotype **) calloc (numLocus, sizeof (Genotype *));
      }
    }
  }

  free (sortedList);
  sortedList = NULL;
  return 0;
}

int
free_multi_locus_genotype_storage (Pedigree * pPedigree)
{
  Person *pPerson;
  int i;

  for (i = 0; i < pPedigree->numPerson; i++) {
    pPerson = pPedigree->ppPersonList[i];
    free (pPerson->multiLocusAdjust);
    free (pPerson->numSavedGenotype2);
    free (pPerson->pTmpLikelihoodIndex);
    if (pPerson->pLikelihood != NULL) {
      /* free space */
      FREE ("pPerson->pLikelihood", pPerson->pLikelihood);
      pPerson->pLikelihood = NULL;
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
  int numGeno = 1;
  int locus;
  int origLocus;
  ConditionalLikelihood *pConditional;

  for (i = 0; i < pPedigree->numPerson; i++) {
    pPerson = pPedigree->ppPersonList[i];
    pPerson->touchedFlag = 0;
    pPerson->numTmpLikelihood = 0;
    size = 1;
    for (locus = locusList->numLocus - 1; locus >= 0; locus--) {
      origLocus = locusList->pLocusIndex[locus];
      if (locus == locusList->numLocus - 1)
	pPerson->multiLocusAdjust[locus] = 1;
      else
	pPerson->multiLocusAdjust[locus] =
	  pPerson->multiLocusAdjust[locus + 1] * numGeno;
      numGeno = pPerson->pSavedNumGenotype[origLocus];
      pPerson->numSavedGenotype2[locus] = numGeno;
      size *= numGeno;
    }
    pPerson->numConditionals = size;

    /*      memset (pPerson->pLikelihood, 0,
       sizeof (ConditionalLikelihood) * pPerson->maxNumConditionals);
     */
    for (j = 0; j < size; j++) {
      pConditional = &pPerson->pLikelihood[j];
      pConditional->touchedFlag = 0;
      pConditional->tmpTouched = 0;
      if (modelOptions.polynomial == TRUE) {
	pConditional->lkslot.likelihoodPolynomial =
	  pConditional->tmpslot.tmpLikelihoodPolynomial = constant0Poly;
	pConditional->wtslot.weightPolynomial = constant1Poly;
      } else {
	pConditional->lkslot.likelihood =
	  pConditional->tmpslot.tmpLikelihood = 0;
	pConditional->wtslot.weight = 1;
      }
    }
  }

  return 0;
}

int
setup_LD_haplotype_freq (LDLoci * pLDLoci, LambdaCell * pCell, int dprimeIdx)
{
  Locus *pLocus1, *pLocus2;
  int numAllele1, numAllele2;
  int i, j;
  double p1, p2, q1, q2;
  double DPrime, DValue, haploFreq;
  double maxD, minD, sum;
  char buf1[MAX_LINE_LEN];
  char buf2[MAX_LINE_LEN];
  char *pBuf1;
  char *pBuf2;
  char *pTemp;

  pBuf1 = &buf1[0];
  pBuf2 = &buf2[0];
  buf1[0] = '\0';
  buf2[0] = '\0';
  pLocus1 = originalLocusList.ppLocusList[pLDLoci->locus1];
  pLocus2 = originalLocusList.ppLocusList[pLDLoci->locus2];
  numAllele1 = pLocus1->numOriginalAllele;
  numAllele2 = pLocus2->numOriginalAllele;
  for (i = 0; i < numAllele1 - 1; i++) {
    p1 = pLocus1->pAlleleFrequency[i];
    p2 = 1 - p1;
    sum = 0;
    for (j = 0; j < numAllele2 - 1; j++) {
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

      //      DPrime = pLDLoci->ppDPrime[i][j];
      DPrime = pCell->lambda[dprimeIdx][i][j];
      sprintf (pBuf1, "%s %4.2f", pBuf2, DPrime);
      /* switch pBuf1 and pBuf2 */
      pTemp = pBuf1;
      pBuf1 = pBuf2;
      pBuf2 = pTemp;

      if (DPrime > 0)
	DValue = DPrime * maxD;
      else
	DValue = -DPrime * minD;
      pCell->DValue[dprimeIdx][i][j] = DValue;

      haploFreq = p1 * q1 + DValue;

      pCell->haploFreq[dprimeIdx][i][j] = haploFreq;

      sum += haploFreq;
    }				/* end of looping the second marker allele frequencies */
    sprintf (pBuf1, "%5.3f %s", p1, pBuf2);
    if ((p1 - sum) >= 0 && (p1 - sum) < LD_E)
      pCell->haploFreq[dprimeIdx][i][j] = LD_E;
    else
      pCell->haploFreq[dprimeIdx][i][j] = p1 - sum;
    //      pLDLoci->ppHaploFreq[i][j] = p1 - sum;
    if ((p1 - sum) < 0) {
      KLOG (LOGINPUTFILE, LOGWARNING,
	    "Haplotype frequency is NEGATIVE - %s.\n", pBuf1);
      return -1;
    }
  }				/* end of looping the first marker allele frequencies */

  /* fill out haplotype frequencies for column [**, n] */
  for (j = 0; j < numAllele2; j++) {
    q1 = pLocus2->pAlleleFrequency[j];
    sum = 0;
    for (i = 0; i < numAllele1 - 1; i++) {
      //      sum += pLDLoci->ppHaploFreq[i][j];
      sum += pCell->haploFreq[dprimeIdx][i][j];
    }

    if ((q1 - sum) >= 0 && (q1 - sum) < LD_E)
      pCell->haploFreq[dprimeIdx][i][j] = LD_E;
    else
      pCell->haploFreq[dprimeIdx][i][j] = q1 - sum;
    if ((q1 - sum) < 0) {
      KLOG (LOGINPUTFILE, LOGWARNING,
	    "Haplotype frequency is NEGATIVE - %s.\n", pBuf1);
      return -1;
    }
  }

  return 0;
}

/* find the place holder for the LD parameters between the given
 * two loci */
LDLoci *
find_LD_loci (int locus1, int locus2)
{
  int i;
  LDLoci *pLDLoci;

  for (i = 0; i < originalLocusList.numLDLoci; i++) {
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
  int ret;

  set_removeGenotypeFlag (TRUE);
  if (modelOptions.polynomial == TRUE) {
    constant1Poly = constantExp (1);
    constant0Poly = constantExp (0);
  }

  /* go through all loci in the original locus list */
  locus = 0;
  while (locus < originalLocusList.numLocus) {
    pLocus = originalLocusList.ppLocusList[locus];
    ped = 0;
    /* go through all peidgress in this set */
    while (ped < pPedigreeSet->numPedigree) {
      pPedigree = pPedigreeSet->ppPedigreeSet[ped];
      /* depending on the locus type, we call different function */

      if (pLocus->locusType == LOCUS_TYPE_TRAIT)
	create_baseline_trait_genotypes (locus, pPedigree);
      else
	create_baseline_marker_genotypes (locus, pPedigree);

      KLOG (LOGGENOELIM, LOGDEBUG, "Baseline Genotype Lists:\n");
      print_pedigree_locus_genotype_list (pPedigree, locus);


      /* first step is do the set recoding: this should help speed up
       * genotype elimination process */

      //fprintf(stderr,"Before 1th allele_set_recoding\n");
      allele_set_recoding (locus, pPedigree);
      //fprintf(stderr,"After  allele_set_recoding\n");

      /* do genotype elimination next */
      ret = pedigree_genotype_elimination (locus, pPedigree);
      KASSERT (ret == 0, "Genotype incompatibility has been detected.\n");
      KLOG (LOGGENOELIM, LOGDEBUG,
	    "Genotype Lists after genotype elimination :\n");
      print_pedigree_locus_genotype_list (pPedigree, locus);

      //fprintf(stderr,"Before 2nd allele_set_recoding\n");
      /* do the set recoding again */
      allele_set_recoding (locus, pPedigree);
      //fprintf(stderr,"After  allele_set_recoding\n");

      /* set genotype weight now to save likelihood calculation time */
      set_genotype_weight (pPedigree, locus);
      set_genotype_position (pPedigree, locus);

      KLOG (LOGGENOELIM, LOGDEBUG, "Genotype Lists after set recoding :\n");
      print_pedigree_locus_genotype_list (pPedigree, locus);

      ped++;
    }
    locus++;
  }

  /* If we ever need to loop over marker allele frequencies, 
   * the super allele frequencies, genotype weights should be set out of 
   * this routine and inside the looping */
  ped = 0;
  while (ped < pPedigreeSet->numPedigree) {
    pPedigree = pPedigreeSet->ppPedigreeSet[ped];
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
               if (modelOptions.polynomial == TRUE )
               {
               expPrinting(g->penslot.penetrancePolynomial);
               fprintf(stderr," Person %d locus %d genotype %d of %d (%d  %d)\n",
                       i+1,j+1,k+1,pPerson->pNumGenotype[j],
                       g->allele[0],
                       g->allele[1]);
               }
               else
               fprintf(stderr,"%f  Person %d locus %d genotype %d of %d (%d  %d)\n",
                       g->penslot.penetrance,
                       i+1,j+1,k+1,pPerson->pNumGenotype[j],
                       g->allele[0],
                       g->allele[1]);
               g = g->pNext;
               k++;
           }
       }
    }
*/

  }				/* loop over pedigrees */

  /* populate the master genotype list */
  populate_saved_genotype_link (pPedigreeSet);

  set_removeGenotypeFlag (FALSE);

  return 0;
}


/* when gene frequency changes, founder's genotype weight needs to be updated 
 * This usually only happens to trait locus  
 * locus - original locus index
 */
int
update_locus (PedigreeSet * pPedigreeSet, int locus)
{
  int ped;
  Pedigree *pPedigree;

  /* update genotype weight and position 
   * in the genotype list 
   * this will be done on the loci in the original locus list */
  ped = 0;
  while (ped < pPedigreeSet->numPedigree) {
    pPedigree = pPedigreeSet->ppPedigreeSet[ped];
    set_genotype_weight (pPedigree, locus);
    //set_genotype_position (pPedigree, locus);
    ped++;
  }

  return 0;
}

/* When penetrance vector changes, each person's penetrance is affected, 
 * thus this function needs to be called */
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

  Polynomial *penPolynomial;

  ped = 0;
  while (ped < pPedigreeSet->numPedigree) {
    pPedigree = pPedigreeSet->ppPedigreeSet[ped];
    for (i = 0; i < pPedigree->numPerson; i++) {
      pPerson = pPedigree->ppPersonList[i];
      /* pass the loop breaker duplicates */
      if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] == NULL)
	continue;
      pGenotype = pPerson->ppSavedGenotypeList[locus];
      while (pGenotype) {
	if (modelOptions.polynomial == TRUE) {
	  compute_penetrance (pPerson, locus,
			      pGenotype->allele[0],
			      pGenotype->allele[1], &penPolynomial);

	  pGenotype->penslot.penetrancePolynomial = penPolynomial;
	} else {
	  compute_penetrance (pPerson, locus,
			      pGenotype->allele[0],
			      pGenotype->allele[1], &pen);

	  pGenotype->penslot.penetrance = pen;
	}


	pGenotype = pGenotype->pSavedNext;
      }
    }
    ped++;
  }
  return 0;
}


/* function to pick markers for multipoint analysis 
 * The strategy is we always want to pick two flanking markers first - the marker 
 * that's the closest among the markers on the left of the trait, and another marker
 * that's the closest among the markers on the right of the trait. 
 * Then we pick the rest of the markers strictly based on the distance. We want to 
 * pick the nearest ones. 
 * Usually we move trait from left to right in a fixed (an input parameter) interval
 */

/* pick the marker immediately to the left of the marker 
 * start - start marker index in the originalLocusList
 * trait - trait position in cM within the chromosome
 * search from left to right
 * mapFlag - male, female, or sex averaged 
 */
int
pick_left_flanking_marker (int start, int end, double traitPosition,
			   int mapFlag)
{
  int lastLeft = -1;
  Locus *pLocus;

  if (start < 0) {
    /* No prior experience. start from the beginning of the marker list */
    start = 0;
  }
  if (end < 0) {
    /* we are allowed to search to the end of the marker list */
    end = originalLocusList.numLocus - 1;
  }

  /* we are choosing marker on the left 
   * locus index starts from 0
   */
  while (start <= end) {
    pLocus = originalLocusList.ppLocusList[start];
    if (pLocus->locusType != LOCUS_TYPE_MARKER) {
      start++;
      continue;
    }
    /* use sex average cM position */
    if (pLocus->pMapUnit->mapPos[mapFlag] >= traitPosition) {
      break;
    }
    /* found one, remember the position 
     * if the marker is right on top of the trait position, we can use it as left 
     * flanking marker 
     */
    lastLeft = start;
    start++;
  }

  /* a return of -1 indicates no marker was found left of the trait */
  return lastLeft;
}

/* pick the marker immediately to the right of the marker 
 * start - start marker index in the originalLocusList
 * end - end marker index in the originalLocusList
 * Both start & end are inclusive
 * trait - trait position in cM within the chromosome
 * search from left to right - find the first one on the right of the trait
 * mapFlag - male, female, or sex averaged 
 */
int
pick_right_flanking_marker (int start, int end, double traitPosition,
			    int mapFlag)
{
  Locus *pLocus;

  if (start < 0) {
    /* No prior experience. start from the beginning of the marker list */
    start = 0;
  }
  if (end < 0) {
    /* we are allowed to search to the end of the marker list */
    end = originalLocusList.numLocus - 1;
  }

  while (start <= end) {
    pLocus = originalLocusList.ppLocusList[start];
    if (pLocus->locusType != LOCUS_TYPE_MARKER) {
      start++;
      continue;
    }
    /* use sex average cM position */
    if (pLocus->pMapUnit->mapPos[mapFlag] >= traitPosition) {
      /* found one, remember the position 
       * if the marker is right on top of the trait position, we can use it as right
       * flanking marker 
       */
      break;
    }
    /* continue the search */
    start++;
  }

  /* a return value bigger than end indicates no marker was found right of the trait */
  return start;
}

/* Pick the closet marker to the traitPosition
 * restricted to pick the firstLeft or the firstRight 
 * start and end provides the boundary 
 * a return of negative or value bigger than end indicates there is no marker available any more 
 */
int
pick_closest_marker (int *left, int *right, int start, int end,
		     double traitPosition, int mapFlag)
{
  Locus *pLocus;

  if (*right > end + 1)
    *right = end + 1;

  if (*left < -1)
    *left = -1;

  if (*left < start) {
    /* nothing on the left of the traitPosition, we probably are at the beginning of a
     * chromosome, so just pick the right marker */
    (*right)++;
    return *right - 1;
  }
  if (*right >= end + 1) {
    /* nothing on the right of the traitPosition, we probably are at the end of the 
     * chromosome, so just pick the left marker */
    (*left)--;
    return *left + 1;
  }
  pLocus = originalLocusList.ppLocusList[*left];
  while (pLocus->locusType == LOCUS_TYPE_TRAIT) {
    (*left)--;
    if (*left < 0) {
      (*right)++;
      return *right - 1;
    }
    pLocus = originalLocusList.ppLocusList[*left];
  }

  /* neither left or right is negative or out of range of the marker list */

  if (traitPosition -
      originalLocusList.ppLocusList[*left]->pMapUnit->mapPos[mapFlag] >
      originalLocusList.ppLocusList[*right]->pMapUnit->mapPos[mapFlag] -
      traitPosition) {
    /* the right distance is shorter, right marker is closer */
    (*right)++;
    return *right - 1;
  } else {
    /* the left distance is shorter/equal */
    (*left)--;
    return *left + 1;
  }

}

/* 
 * pick markers set up the locusList for multipoint analysis
 * numMarkers - number of markers to use for the multipoint analysis
 * leftMarker - since we are moving left to right, if you know the last leftMarker, that 
 *              will speed up the search, if you don't know, put -1
 * traitPosition - in cM
 * mapFlag - male, female, or sex averaged map
 * start, end - start and end marker locus index in the originalLocusList
 * pLocusList - should have contained the trait locus already
 */
int
add_markers_to_locuslist (SubLocusList * pLocusList,
			  int numMarkers, int *pLeftMarker, int start,
			  int end, double traitPosition, int mapFlag)
{
  /* take out the trait locus */
  int total = originalLocusList.numLocus - 1;
  int numMarkerSelected = 0;
  int rightMarker;
  int marker;
  int leftMarker;
  int ret;			/* return status */
  int i;

  /* number of markers to use can't exceed what we have */
  if (numMarkers > total) {
    /* out of our reach, got to say no */
    return -1;
  }
  if (numMarkers < 2) {
    /* we are not doing multipoint, why are we here? */
    return -2;
  }

  if (numMarkers == total) {
    /* no brainer - select all of them */
    for (i = 0; i < originalLocusList.numLocus; i++) {
      if (originalLocusList.ppLocusList[i]->locusType == LOCUS_TYPE_MARKER) {
	ret =
	  add_analysis_locus (pLocusList, i, DIRECTION_RIGHT,
			      map.mapFunction);
      }
    }
    return 0;

  } else {
    /* we got some work to do 
     * first - select flanking markers */
    *pLeftMarker =
      pick_left_flanking_marker (*pLeftMarker, end, traitPosition, mapFlag);
    if (*pLeftMarker >= 0) {
      numMarkerSelected++;
      ret =
	add_analysis_locus (pLocusList, *pLeftMarker, DIRECTION_RIGHT,
			    map.mapFunction);
    }
    /* left marker could be returned as -1, that means the search for right starts at 0 */
    rightMarker =
      pick_right_flanking_marker (*pLeftMarker + 1, end, traitPosition,
				  mapFlag);
    if (rightMarker <= end) {
      numMarkerSelected++;
      ret =
	add_analysis_locus (pLocusList, rightMarker, DIRECTION_RIGHT,
			    map.mapFunction);
    }

    /* move the markers bidirectionally and pick the rest of the markers */
    leftMarker = *pLeftMarker;
    leftMarker--;
    rightMarker++;
    while (numMarkerSelected < numMarkers) {
      /* left marker or right marker will be incremented if they were the marker chosen */
      marker = pick_closest_marker (&leftMarker, &rightMarker, start, end,
				    traitPosition, mapFlag);
      /* should not encounter this */
      if (marker < 0 || marker > end) {
	KLOG (LOGDEFAULT, LOGWARNING,
	      "We have run out of markers(%d out of %d) to choose for multipoint.\n",
	      numMarkerSelected, total);
	break;
      }
      numMarkerSelected++;
      ret =
	add_analysis_locus (pLocusList, marker, DIRECTION_LEFT,
			    map.mapFunction);
    }
  }

  return 0;
}

/* Add locus into the analysis locus list
 * this should only be used for multipoint analysis
 * pLocusList - where the list is stored at
 * locus - index to the originalLocusList starts from 0
 * directionFlag - left or right indicator as sometimes marker is right on top of trait 
 * to be robust, really should check whether we add duplicate locus in the list or not
 */
int
add_analysis_locus (SubLocusList * pLocusList, int locus, int directionFlag,
		    int mapFlag)
{
  int numLocus;
  int first;
  int last;
  int curr;
  double *firstPos;
  double *lastPos;
  double *thisPos;
  double *currPos;
  int i, j, k;
  int left;
  double *leftPos;

  /* if we can't return anything, what's the purpose of proceeding */
  if (pLocusList == NULL)
    return -1;

  numLocus = pLocusList->numLocus;
  if (numLocus + 1 > pLocusList->maxNumLocus) {
    return -1;
#if 0
    /* need to allocate/reallocate space */
    pLocusList->maxNumLocus += DEF_LOCUS_MALLOC_INCREMENT;
    pLocusList->pLocusIndex = (int *) REALLOC ("pLocusList->pLocusIndex",
					       pLocusList->pLocusIndex,
					       sizeof (int) *
					       locusList->maxNumLocus);
    /* 0 - sex averaged, 1 - male, 2 - female */
    for (i = 0; i < 3; i++) {
      pLocusList->pPrevLocusDistance[i] =
	(double *) REALLOC ("pLocusList->pPrevLocusDistance",
			    pLocusList->pPrevLocusDistance[i],
			    sizeof (double) * locusList->maxNumLocus);
      pLocusList->pNextLocusDistance[i] =
	(double *) REALLOC ("pLocusList->pNextLocusDistance",
			    pLocusList->pNextLocusDistance[i],
			    sizeof (double) * locusList->maxNumLocus);
    }
#endif
  }

  if (numLocus == 0) {
    pLocusList->pLocusIndex[0] = locus;
    /* 0 - sex averaged, 1 - male, 2 - female */
    for (i = 0; i < 3; i++) {
      pLocusList->pPrevLocusDistance[i][0] = -1;
      pLocusList->pNextLocusDistance[i][0] = -1;
    }
    pLocusList->numLocus++;
    /* done */
    return 0;
  }

  thisPos = get_map_position (locus);
  first = pLocusList->pLocusIndex[0];
  firstPos = get_map_position (first);
  if ((originalLocusList.ppLocusList[first]->locusType == LOCUS_TYPE_MARKER &&
       originalLocusList.ppLocusList[locus]->locusType == LOCUS_TYPE_MARKER &&
       locus < first) || thisPos[MAP_SEX_AVERAGE] < firstPos[MAP_SEX_AVERAGE]
      || (directionFlag == DIRECTION_LEFT
	  && (thisPos[MAP_SEX_AVERAGE] >
	      firstPos[MAP_SEX_AVERAGE] - ERROR_MARGIN)
	  && (thisPos[MAP_SEX_AVERAGE] <
	      firstPos[MAP_SEX_AVERAGE] + ERROR_MARGIN)))
  {
    /* add this locus to the beginning of the list, shift the rest */
    for (i = numLocus - 1; i >= 0; i--) {
      pLocusList->pLocusIndex[i + 1] = pLocusList->pLocusIndex[i];
      /* 0 - sex averaged, 1 - male, 2 - female */
      for (j = 0; j < 3; j++) {
	pLocusList->pPrevLocusDistance[j][i + 1] =
	  pLocusList->pPrevLocusDistance[j][i];
	pLocusList->pNextLocusDistance[j][i + 1] =
	  pLocusList->pNextLocusDistance[j][i];
      }
    }
    pLocusList->numLocus++;
    pLocusList->pLocusIndex[0] = locus;
    for (j = 0; j < 3; j++) {
      pLocusList->pPrevLocusDistance[j][0] = -1;

    }
    /* calculate the distance */
    pLocusList->pNextLocusDistance[0][0] =
      pLocusList->pPrevLocusDistance[0][1] =
      cm_to_recombination_fraction (firstPos[0] - thisPos[0],
				    map.mapFunction);
    if (modelOptions.mapFlag == SEX_AVERAGED) {
      pLocusList->pNextLocusDistance[1][0] =
	pLocusList->pPrevLocusDistance[1][1] =
	pLocusList->pNextLocusDistance[0][0];
      pLocusList->pNextLocusDistance[2][0] =
	pLocusList->pPrevLocusDistance[2][1] =
	pLocusList->pNextLocusDistance[0][0];
    } else {
      for (j = 1; j <= 2; j++) {
	pLocusList->pNextLocusDistance[j][0] =
	  pLocusList->pPrevLocusDistance[j][1] =
	  cm_to_recombination_fraction (firstPos[j] - thisPos[j],
					map.mapFunction);
      }
    }
    /* done */
    return 0;
  }

  last = pLocusList->pLocusIndex[numLocus - 1];
  lastPos = get_map_position (last);
  if ((originalLocusList.ppLocusList[last]->locusType == LOCUS_TYPE_MARKER &&
       originalLocusList.ppLocusList[locus]->locusType == LOCUS_TYPE_MARKER &&
       locus > last) || thisPos[MAP_SEX_AVERAGE] > lastPos[MAP_SEX_AVERAGE]
      || (directionFlag == DIRECTION_RIGHT
	  && (thisPos[MAP_SEX_AVERAGE] >
	      lastPos[MAP_SEX_AVERAGE] - ERROR_MARGIN)
	  && (thisPos[MAP_SEX_AVERAGE] <
	      lastPos[MAP_SEX_AVERAGE] + ERROR_MARGIN)))
  {
    /* add this locus to the end of the list */
    pLocusList->numLocus++;
    pLocusList->pLocusIndex[numLocus] = locus;
    for (j = 0; j < 3; j++) {
      pLocusList->pNextLocusDistance[j][numLocus] = -1;
    }

    /* calculate the distance */
    pLocusList->pPrevLocusDistance[0][numLocus] =
      pLocusList->pNextLocusDistance[0][numLocus - 1] =
      cm_to_recombination_fraction (thisPos[0] - lastPos[0], map.mapFunction);
    if (modelOptions.mapFlag == SEX_AVERAGED) {
      pLocusList->pPrevLocusDistance[1][numLocus] =
	pLocusList->pNextLocusDistance[1][numLocus - 1] =
	pLocusList->pPrevLocusDistance[0][numLocus];
      pLocusList->pPrevLocusDistance[2][numLocus] =
	pLocusList->pNextLocusDistance[2][numLocus - 1] =
	pLocusList->pPrevLocusDistance[0][numLocus];
    } else {
      for (j = 1; j <= 2; j++) {
	pLocusList->pPrevLocusDistance[j][numLocus] =
	  pLocusList->pNextLocusDistance[j][numLocus - 1] =
	  cm_to_recombination_fraction (thisPos[j] - lastPos[j],
					map.mapFunction);
      }
    }

    /* done */
    return 0;
  }

  /* We can't add it to the beginning or the end of the list
   * Find the right place first */
  left = first;
  leftPos = firstPos;
  i = 1;
  while (i < numLocus) {
    curr = pLocusList->pLocusIndex[i];
    currPos = get_map_position (curr);
    if ((originalLocusList.ppLocusList[curr]->locusType == LOCUS_TYPE_MARKER
	 && originalLocusList.ppLocusList[locus]->locusType ==
	 LOCUS_TYPE_MARKER && locus < curr)
	|| thisPos[MAP_SEX_AVERAGE] < currPos[MAP_SEX_AVERAGE]
	|| (directionFlag == DIRECTION_LEFT
	    && (thisPos[MAP_SEX_AVERAGE] >
		currPos[MAP_SEX_AVERAGE] - ERROR_MARGIN)
	    && (thisPos[MAP_SEX_AVERAGE] <
		currPos[MAP_SEX_AVERAGE] + ERROR_MARGIN))) {
      /* put it to the left the curr locus, shift the rest */
      for (j = numLocus - 1; j >= i; j--) {
	pLocusList->pLocusIndex[j + 1] = pLocusList->pLocusIndex[j];
	for (k = 0; k < 3; k++) {
	  pLocusList->pNextLocusDistance[k][j + 1] =
	    pLocusList->pNextLocusDistance[k][j];
	}
      }
      pLocusList->numLocus++;
      pLocusList->pLocusIndex[i] = locus;
      pLocusList->pPrevLocusDistance[0][i] =
	pLocusList->pNextLocusDistance[0][i - 1] =
	cm_to_recombination_fraction (thisPos[0] - leftPos[0],
				      map.mapFunction);
      /* calculate the distance */
      pLocusList->pNextLocusDistance[0][i] =
	pLocusList->pPrevLocusDistance[0][i + 1] =
	cm_to_recombination_fraction (currPos[0] - thisPos[0],
				      map.mapFunction);
      if (modelOptions.mapFlag == SEX_SPECIFIC) {
	for (k = 1; k <= 2; k++) {
	  pLocusList->pPrevLocusDistance[k][i] =
	    pLocusList->pNextLocusDistance[k][i - 1] =
	    cm_to_recombination_fraction (thisPos[k] - leftPos[k],
					  map.mapFunction);
	  /* calculate the distance */
	  pLocusList->pNextLocusDistance[k][i] =
	    pLocusList->pPrevLocusDistance[k][i + 1] =
	    cm_to_recombination_fraction (currPos[k] - thisPos[k],
					  map.mapFunction);

	}
      } else {
	pLocusList->pPrevLocusDistance[1][i] =
	  pLocusList->pNextLocusDistance[1][i - 1] =
	  pLocusList->pPrevLocusDistance[0][i];
	pLocusList->pPrevLocusDistance[2][i] =
	  pLocusList->pNextLocusDistance[2][i - 1] =
	  pLocusList->pPrevLocusDistance[0][i];
	pLocusList->pNextLocusDistance[1][i] =
	  pLocusList->pPrevLocusDistance[1][i + 1] =
	  pLocusList->pNextLocusDistance[0][i];
	pLocusList->pNextLocusDistance[2][i] =
	  pLocusList->pPrevLocusDistance[2][i + 1] =
	  pLocusList->pNextLocusDistance[0][i];
      }
      /* done */
      return 0;
    }
    i++;
    left = curr;
    leftPos = currPos;
  }
  return 0;
}

/* retrieve the map position for a marker locus or trait locus 
 */
double *
get_map_position (int locus)
{
  Locus *pLocus;

  if (locus < 0 && locus >= originalLocusList.numLocus) {
    return NULL;
  }
  pLocus = originalLocusList.ppLocusList[locus];
  if (pLocus->locusType == LOCUS_TYPE_MARKER) {
    return pLocus->pMapUnit->mapPos;
  } else {
    return pLocus->pTraitLocus->mapPosition;
  }

}

void
free_sub_locus_list (SubLocusList * pLocusList)
{
  int i;

  for (i = 0; i < 3; i++) {
    free (pLocusList->pPrevLocusDistance[i]);
    free (pLocusList->pNextLocusDistance[i]);
  }
  free (pLocusList->pLocusIndex);
  pLocusList->numLocus = 0;
  pLocusList->maxNumLocus = 0;

  /* polynomial free */
}

void
final_cleanup ()
{
  free_locus_list (&originalLocusList);
  free_map (&map);
}

void
free_pedigree_set (PedigreeSet * pPedigreeSet)
{
  int i, k;
  Pedigree *pPedigree;

  for (i = 0; i < pPedigreeSet->numPedigree; i++) {
    pPedigree = pPedigreeSet->ppPedigreeSet[i];
    //      free_multi_locus_genotype_storage (pPedigree);

    free (pPedigree->ppFounderList);

    /* free person list */
    for (k = 0; k < pPedigree->numPerson; k++) {
      free_person (pPedigree->ppPersonList[k]);
    }
    free (pPedigree->ppPersonList);

    /* free nuclear family list */
    for (k = 0; k < pPedigree->numNuclearFamily; k++) {
      free_nuclear_family (pPedigree->ppNuclearFamilyList[k]);
    }
    free (pPedigree->ppNuclearFamilyList);
    free (pPedigree->ppFounderNuclearFamilyList);
    free (pPedigree->loopBreakerList);
    free (pPedigree->pCount);
    free (pPedigree);
  }

  free (pPedigreeSet->nullLikelihood);
  free (pPedigreeSet->ppPedigreeSet);
  free (pPedigreeSet->pDonePerson);
}

void
free_person (Person * pPerson)
{
  int i;

  for (i = 0; i < originalLocusList.numTraitLocus; i++) {
    free (pPerson->ppOrigTraitValue[i]);
    free (pPerson->ppTraitValue[i]);
    free (pPerson->ppTraitKnown[i]);
    free (pPerson->ppLiabilityClass[i]);
  }
  free (pPerson->ppOrigTraitValue);
  free (pPerson->ppTraitValue);
  free (pPerson->ppTraitKnown);
  free (pPerson->ppLiabilityClass);

  free (pPerson->ppSpouseList);

  free (pPerson->pPhenotypeList[0]);
  free (pPerson->pPhenotypeList[1]);
  free (pPerson->pPhasedFlag);
  free (pPerson->pTypedFlag);

  if (pPerson->loopBreaker == 0 || pPerson->pParents[DAD] != NULL) {
    /* go through each locus for the genotypes */
    for (i = 0; i < originalLocusList.numLocus; i++) {
      /* delete each genotype in the list */
      while (pPerson->ppSavedGenotypeList[i] != NULL) {
	remove_genotype (&(pPerson->ppSavedGenotypeList[i]),
			 pPerson->ppSavedGenotypeList[i],
			 &pPerson->pSavedNumGenotype[i]);
      }
    }
    if (pPerson->loopBreaker == 1) {
      for (i = 0; i < pPerson->loopBreakerStruct->maxNumGenotype; i++)
	free (pPerson->loopBreakerStruct->genotype[i]);
      free (pPerson->loopBreakerStruct->genotype);
      free (pPerson->loopBreakerStruct);
    }
  }

  free (pPerson->ppGenotypeList);
  free (pPerson->pNumGenotype);
  free (pPerson->ppSavedGenotypeList);
  free (pPerson->pSavedNumGenotype);
  free (pPerson->ppShadowGenotypeList);
  free (pPerson->pShadowGenotypeListLen);
  free (pPerson->ppProbandGenotypeList);
  free (pPerson->pProbandNumGenotype);

  free (pPerson->ppHaplotype);

  free (pPerson->pTransmittedAlleles[DAD]);
  free (pPerson->pTransmittedAlleles[MOM]);
  free (pPerson->pNonTransmittedAlleles[DAD]);
  free (pPerson->pNonTransmittedAlleles[MOM]);

  free (pPerson->ppNuclearFamilyList);

  free (pPerson);
}

void
free_nuclear_family (NuclearFamily * pNucFam)
{
  NuclearFamilyConnector *pConnect, *pNext;

  pConnect = pNucFam->pUpConnectors;
  while (pConnect != NULL) {
    pNext = pConnect->pNextConnector;
    free (pConnect);
    pConnect = pNext;
  }

  pConnect = pNucFam->pDownConnectors;
  while (pConnect != NULL) {
    pNext = pConnect->pNextConnector;
    free (pConnect);
    pConnect = pNext;
  }

  free (pNucFam->hetFlag[DAD]);
  free (pNucFam->hetFlag[MOM]);
  free (pNucFam->tmpNumHet[DAD]);
  free (pNucFam->tmpNumHet[MOM]);
  free (pNucFam->relatedPPairStart);
  free (pNucFam->numRelatedPPair);
  free (pNucFam->totalRelatedPPair);
  free (pNucFam->ppChildrenList);
  free (pNucFam);
}


int
free_LD_loci (LDLoci * pLocus)
{
  int i;

  if (pLocus->ppDPrime == NULL)
    return 0;

  for (i = 0; i < pLocus->numAllele1 - 1; i++) {
    free (pLocus->ppDPrime[i]);
    free (pLocus->ppDValue[i]);
    free (pLocus->ppHaploFreq[i]);
  }
  free (pLocus->ppHaploFreq[pLocus->numAllele1 - 1]);

  free (pLocus->ppDPrime);
  free (pLocus->ppDValue);
  free (pLocus->ppHaploFreq);
  pLocus->ppDPrime = NULL;
  pLocus->ppDValue = NULL;
  pLocus->ppHaploFreq = NULL;

  return 0;
}

int
reallocate_LD_loci (LDLoci * pLocus, int m, int n)
{
  int i;

  /* if no size change, no need to realloc */
  if (pLocus->numAllele1 == m && pLocus->numAllele2 == n)
    return 0;

  /* size differs, free the old one first, then realloc */
  free_LD_loci (pLocus);
  pLocus->numAllele1 = m;
  pLocus->numAllele2 = n;

  pLocus->ppDPrime = (double **) calloc (m - 1, sizeof (double *));
  pLocus->ppDValue = (double **) calloc (m - 1, sizeof (double *));
  pLocus->ppHaploFreq = (double **) calloc (m, sizeof (double *));
  for (i = 0; i < m - 1; i++) {
    pLocus->ppDPrime[i] = (double *) calloc (n - 1, sizeof (double));
    pLocus->ppDValue[i] = (double *) calloc (n - 1, sizeof (double));
    pLocus->ppHaploFreq[i] = (double *) calloc (n, sizeof (double));
  }
  pLocus->ppHaploFreq[m - 1] = (double *) calloc (n, sizeof (double));

  return 0;
}

int
set_null_dprime (LDLoci * pLocus)
{
  int i;

  for (i = 0; i < pLocus->numAllele1 - 1; i++) {
    memset (pLocus->ppDPrime[i], 0,
	    sizeof (double) * (pLocus->numAllele2 - 1));
  }
  return 0;
}

int
copy_dprime (LDLoci * pLocus, double **pSrc)
{
  int i;

  for (i = 0; i < pLocus->numAllele1 - 1; i++) {
    memcpy (pLocus->ppDPrime[i], pSrc[i],
	    sizeof (double) * (pLocus->numAllele2 - 1));
  }
  return 0;
}

int
copy_haploFreq (LDLoci * pLocus, double **pSrc)
{
  int i;

  for (i = 0; i < pLocus->numAllele1; i++) {
    memcpy (pLocus->ppHaploFreq[i], pSrc[i],
	    sizeof (double) * (pLocus->numAllele2));
  }
  return 0;
}

int
copy_DValue (LDLoci * pLocus, double **pSrc)
{
  int i;

  for (i = 0; i < pLocus->numAllele1 - 1; i++) {
    memcpy (pLocus->ppDValue[i], pSrc[i],
	    sizeof (double) * (pLocus->numAllele2 - 1));
  }
  return 0;
}

int
isDPrime0 (double **ppDPrime, int m, int n)
{
  int i, j;

  for (i = 0; i < m - 1; i++)
    for (j = 0; j < n - 1; j++) {
      if (ppDPrime[i][j] < -ERROR_MARGIN || ppDPrime[i][j] > ERROR_MARGIN)
	return 0;
    }
  return 1;
}

/* find a locus in a locus list */
int
find_locus (LocusList * pLocusList, char *sName)
{
  Locus *pLocus;
  int i;

  for (i = 0; i < pLocusList->numLocus; i++) {
    pLocus = pLocusList->ppLocusList[i];
    if (!strcmp (pLocus->sName, sName)) {
      return i;
    }
  }
  return -1;
}

/* populate the master list of valid genotype list for given pedigree 
 * locus - index in original locus list */
void
populate_pedigree_saved_genotype_link (int locus, Pedigree * pPed)
{
  Person *pPerson;
  int i;
  Genotype *pGeno;

  for (i = 0; i < pPed->numPerson; i++) {
    pPerson = pPed->ppPersonList[i];
    if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] == NULL) {
      /* pass the duplicated loop breaker */
      continue;
    }
    pGeno = pPerson->ppGenotypeList[locus];
    pPerson->ppSavedGenotypeList[locus] = pGeno;
    pPerson->pSavedNumGenotype[locus] = pPerson->pSavedNumGenotype[locus];
    while (pGeno != NULL) {
      pGeno->pSavedNext = pGeno->pNext;
      pGeno = pGeno->pNext;
    }
  }
}

/* populate the master list of valid genotype list for all pedigrees */
void
populate_saved_genotype_link (PedigreeSet * pSet)
{
  int i, j;
  Pedigree *pPed;

  for (i = 0; i < originalLocusList.numLocus; i++) {
    for (j = 0; j < pSet->numPedigree; j++) {
      pPed = pSet->ppPedigreeSet[j];
      populate_pedigree_saved_genotype_link (i, pPed);
    }
  }
}

/* likelihood on pedigrees with loops are calculated with fixing loop breaker
 * genotype vector one at a time, followed by genotype elimination on the entire
 * pedigree. so after each calculation, the genotype list
 * needs to be restored for everyone from the saved master list */
void
restore_pedigree_genotype_link_from_saved (Pedigree * pPed)
{
  int i, j;
  int origLocus;
  Genotype *pGeno;
  Person *pPerson;

  for (j = 0; j < locusList->numLocus; j++) {
    origLocus = locusList->pLocusIndex[j];

    for (i = 0; i < pPed->numPerson; i++) {
      pPerson = pPed->ppPersonList[i];
      if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] == NULL) {
	/* pass the duplicated loop breaker */
	continue;
      }
      pGeno = pPerson->ppSavedGenotypeList[origLocus];
      pPerson->ppGenotypeList[origLocus] = pGeno;
      pPerson->pNumGenotype[origLocus] =
	pPerson->pSavedNumGenotype[origLocus];
      while (pGeno != NULL) {
	pGeno->pNext = pGeno->pSavedNext;
	pGeno = pGeno->pSavedNext;
      }
    }
  }
}
