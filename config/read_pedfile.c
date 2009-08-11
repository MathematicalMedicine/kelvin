/**
@file read_pedfile.c

 As for Vitesse, reads in post-makeped format pedigree files.

 Copyright 2008, Nationwide Children's Research Institute.  All rights
 reserved.  Permission is hereby given to use this software for
 non-profit educational purposes only.

 For the initial design, this library will read in post-makepd format
 file only. Standard post-makeped format has the following columns: 

 - PedID IndividualID FatherID MomID 1stChildID NextPaternalSiblingID
 - NextMaternalSiblingID Sex ProbandIndicator 
 - AffectionStatus(or QuantitativeTrait) 
 - Genotypes(a pair for each marker. One for a chromosome)
 - Pedigree ID (again)
 - Original Person ID (This can tell which individual is doubled. That is 
   two lines with different individual IDs with the same original Person IDs 
   at the end of the lines )

 Pre-makedped is simpler, only has:

 - PedID IndividualID FatherID MotherID Sex AffectionStatus(QuantitativeTrait)
 - Genotyping pairs for each marker

 The benefit of reading in post-makeped is that we don't have to deal with
 breaking loops as presumbly that some other program has been called 
 to break the loops or individuals have been doubled manually and
 proband is selected as well. Programs to set up post-makeped file and 
 break up loops are:

 <pre>
 makeped <pre-makeped pedfile> <pedfile.dat> n
   even though there are loops, still provide "n" when you run makeped

 unknown -l   (it needs pedfile.dat-post-makeped and datafile.dat to be
 in the same directory where you run unknown 
 </pre>

 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "../pedlib/pedlib.h"
#include "../pedlib/likelihood.h"
#include "../utils/utils.h"

// Functions listed below should not be called outside this library:
Pedigree *create_pedigree (PedigreeSet * pPedigreeSet, char *sPedLabel);
Person *create_person (Pedigree * pPed, char *sPersonLabel);
int add_founder (Pedigree * pPed, Person * pPerson);
NuclearFamily *create_nuclear_family (Pedigree * pPed);
NuclearFamilyConnector *create_nuclear_family_connector (Person * pPerson,
							 NuclearFamily *
							 pNucFam1,
							 int direction1,
							 NuclearFamily *
							 pNucFam2);
int add_person_nuclear_family (Person * pPerson, NuclearFamily * pNucFam);
int add_nuclear_family_parents (NuclearFamily * pNucFam,
				Person * pDad, Person * pMom);
int add_nuclear_family_child (NuclearFamily * pNucFam, Person * pChild);
int add_founder_nuclear_family (NuclearFamily * pNucFam);
int add_spouse (Person * pPerson, Person * pSpouse);
int read_person (char *sPedfileName, int lineNo, char *pLine,
		 Person * pPerson);
int setup_pedigree_ptrs (Pedigree * pPed);
int setup_nuclear_families (Pedigree * pPed);
int setup_loop_counts (Pedigree * pPed);
void add_loopbreaker (Pedigree * pPed, Person * pPerson);
Person *find_original_person (Pedigree * pPed, char *sPersonID);
Person *find_loop_breaker (Person * breaker);
void check_for_loop(Pedigree *);

/**

 Function to read the pedigree file - expected to be postmakeped format 
 and load data into PedigreeSet and Person structures.

*/
int
read_pedfile (char *sPedfileName, PedigreeSet * pPedigreeSet)
{
  FILE *fpPedfile = NULL;
  int lineNo = 0;
  char sPrevPedLabel[MAX_PED_LABEL_LEN];
  char sCurrPedLabel[MAX_PED_LABEL_LEN];
  char sCurrPersonLabel[MAX_PED_LABEL_LEN];
  int pos;
  int numRet;
  char *pLine = NULL;		/* current pointer to a string for strtok()  */
  Pedigree *pCurrPedigree = NULL;
  Person *pCurrPerson = NULL;
  int lastFlag = 0, i;

  /* Initialize pedigree set 
   * no longer doing initialization here, pPedigreeSet should have 
   * been created and initialized before this and pTrait field set */
  //memset(pPedigreeSet, 0, sizeof(PedigreeSet));

  /* open pedigree file */
  fpPedfile = fopen (sPedfileName, "r");
  KASSERT (fpPedfile != NULL,
	   "Can't open pedigree file %s for read. Exiting!", sPedfileName);

  sPrevPedLabel[0] = '\0';
  sPrevPedLabel[MAX_PED_LABEL_LEN - 1] = '\0';
  sCurrPedLabel[MAX_PED_LABEL_LEN - 1] = '\0';

  while (1) {
    /* lastFlag is to force some processing of the data we already read */
    if (lastFlag == 0
	&& fgetlongs (&flexBuffer, &flexBufferSize, fpPedfile) == NULL)
      lastFlag = 1;
    if (lastFlag == 0) {
      lineNo++;
      /* ignore blank or white space or comment lines */
      if (is_line_blank_or_comment (flexBuffer))
	continue;
      /* we got some data to process */
      numRet =
	sscanf (flexBuffer, "%s %s %n", sCurrPedLabel, sCurrPersonLabel,
		&pos);
      KASSERT (numRet == 2,
	       "Can't get pedigree Label from line no %d in pedfile %s.\n",
	       lineNo, sPedfileName);
      pLine = &flexBuffer[pos];
    }
    /* if this is not the first pedigree and it has a different pedigree
     * * Label than the previous one, it indicates a new pedigree starts now
     * * */
    if (strcmp (sPrevPedLabel, sCurrPedLabel) != 0) {
      // Make sure we haven't seen this pedigree before
      for (i = 0; i < pPedigreeSet->numPedigree; i++) {
	if (strcmp (sCurrPedLabel, pPedigreeSet->ppPedigreeSet[i]->sPedigreeID) == 0) {
	  KLOG (LOGPEDFILE, LOGFATAL, "Pedigree %s appears in separate places in the pedigree file, please correct and re-run.\n",
		sCurrPedLabel);
	}
      }
    }
    if (lastFlag == 1 ||strcmp (sPrevPedLabel, sCurrPedLabel) != 0) {
      /* a different ped Label indicates that we have got a new pedigree 
       * before we move onto the new pedigree, we need to set up 
       * * pointers for the current pedigree
       * * */
      /* If there is no proband found for this pedigree, default to the
       * first person and generate a warning message */
      if (pCurrPedigree && pCurrPedigree->pPeelingProband == NULL) {
	pCurrPedigree->pPeelingProband = *(pCurrPedigree->ppPersonList);
	KLOG (LOGPEDFILE, LOGWARNING,
	      "WARNING: No proband was given for this pedigree %s and proband is set to person %s.\n",
	      pCurrPedigree->sPedigreeID,
	      pCurrPedigree->pPeelingProband->sID);
	
      }
      if (pCurrPedigree) {
	pPedigreeSet->pDonePerson =
	  (int *) realloc (pPedigreeSet->pDonePerson,
			   sizeof (int) * pPedigreeSet->maxNumPerson);
	memset (pPedigreeSet->pDonePerson, 0,
		sizeof (int) * pPedigreeSet->maxNumPerson);
	
	/* set up pointers in the pedigree such as first child, 
	 * * * next paternal sibling, next maternal sibling */
	setup_pedigree_ptrs (pCurrPedigree);
	/* Really figure out if there is a loop */
	check_for_loop (pCurrPedigree);
	/* set up loop count, loop breaker count for this pedigree */
	setup_loop_counts (pCurrPedigree);
	/* set up nuclear families inside of this pedigree */
	setup_nuclear_families (pCurrPedigree);
	/* print out just for debug and verification purpose for now */
	print_nuclear_family (stdout, pCurrPedigree);
      }
      if (lastFlag == 1)
	/* we have processed the last pedigree in the file, done */
	return 0;
      /* create new pedigree */
      pCurrPedigree = create_pedigree (pPedigreeSet, sCurrPedLabel);
      strcpy (sPrevPedLabel, sCurrPedLabel);
    }

    // Make sure we haven't seen this person before
    for (i = 0; i < pCurrPedigree->numPerson; i++) {
      if (strcmp (sCurrPersonLabel, pCurrPedigree->ppPersonList[i]->sID) == 0) {
	KLOG (LOGPEDFILE, LOGFATAL, "Person %s appears more than once in pedigree %s, please correct and re-run.\n",
	      sCurrPersonLabel, sCurrPedLabel);
      }
    }

    /* each line should be a new person for this pedigree */
    pCurrPerson = create_person (pCurrPedigree, sCurrPersonLabel);

    /* read in this person's information from current line */
    read_person (sPedfileName, lineNo, pLine, pCurrPerson);

    /* mark this pedigree as having a loop if so */
    if (pCurrPerson->proband > 1) {
      pCurrPedigree->loopFlag = 1;
      pCurrPedigree->pPedigreeSet->loopFlag = 1;
    }

    if (pCurrPerson->proband == 1) {
      /* this person is selected as a peeling proband for this pedigree */
      pCurrPedigree->pPeelingProband = pCurrPerson;
    }
  }
  /* Finished reading pedigree file successfully */
  return 0;
}

int
read_person (char *sPedfileName, int lineNo, char *pLine, Person * pPerson)
{
  int numRet;			/* number of items read by sscanf() */
  int pos;			/* position of the unread in the buffer */
  int numMarker;		/* number of markers we have read genotypes */
  //  char tmpStr[MAX_LINE_LEN];

  //Pedigree *pPed = pPerson->pPedigree;
  Trait *pTrait;
  TraitLocus *pTraitLocus;
  int numTrait;
  int i, j;
  int ret;
  Locus *pLocus;
  /* input allele length is limited to 7 */
  char a1[8];
  char a2[8];

  //numRet = sscanf (pLine, "%s %s %n", pPerson->sDadID, pPerson->sMomID, &pos);
  numRet = sscanf (pLine, "%s %s %n",
		   pPerson->sParentID[DAD], pPerson->sParentID[MOM], &pos);
  KASSERT (numRet == 2,
	   "Line %d in pedfile %s doesn't have enough columns (parents). Is this a post-makeped file? \n",
	   lineNo, sPedfileName);
  /* ignore the next three columns: first child, next paternal sibling,
   * next maternal sibling, as we are going to build those pointers 
   * later once we have all the individuals in the pedigree */
  pLine = &pLine[pos];
  numRet = sscanf (pLine, "%s %s %s %n",
		   pPerson->sFirstChildID,
		   pPerson->sNextSibID[DAD], pPerson->sNextSibID[MOM], &pos);
  KASSERT (numRet == 3,
	   "Line %d in pedfile %s doesn't have enough columns (child & siblings). Is this a post-makeped file? \n",
	   lineNo, sPedfileName);
  pLine = &pLine[pos];
  numRet = sscanf (pLine, "%d %d %n", &pPerson->sex, &pPerson->proband, &pos);
  /* input use 1 as MALE and 2 as FEMALE, while internally we use 0 - MALE
   * and 1 - FEMALE */
  pPerson->sex -= 1;
  KASSERT (numRet == 2,
	   "Line %d in pedfile %s doesn't have enough columns (sex & proband). Is this a post-makeped file? \n",
	   lineNo, sPedfileName);
  if (pPerson->proband > 1)
    pPerson->loopBreaker = 1;

  pLine = &pLine[pos];

  /* Now we expect to read in disease locus related information 
   * if we ever got a disease locus 
   * i.e. disease locus information is expected to be listed 
   * before any marker locus */
  pPerson->ppOrigTraitValue =
    (double **) malloc (sizeof (double *) * originalLocusList.numTraitLocus);
  pPerson->ppTraitValue =
    (double **) malloc (sizeof (double *) * originalLocusList.numTraitLocus);
  pPerson->ppTraitKnown =
    (int **) malloc (sizeof (int *) * originalLocusList.numTraitLocus);
  pPerson->ppLiabilityClass =
    (int **) malloc (sizeof (int *) * originalLocusList.numTraitLocus);

  i = 0;
  while (i < originalLocusList.numTraitLocus) {
    pTraitLocus = originalLocusList.ppLocusList[i]->pTraitLocus;
    numTrait = pTraitLocus->numTrait;
    pPerson->ppOrigTraitValue[i] =
      (double *) malloc (sizeof (double) * numTrait);
    pPerson->ppTraitValue[i] = (double *) malloc (sizeof (double) * numTrait);
    pPerson->ppTraitKnown[i] = (int *) malloc (sizeof (int) * numTrait);
    memset (pPerson->ppTraitKnown[i], 0, sizeof (int) * numTrait); // Effectively sets to FALSE
    pPerson->ppLiabilityClass[i] = (int *) malloc (sizeof (int) * numTrait);
    j = 0;
    while (j < numTrait) {
      pTrait = pTraitLocus->pTraits[j];
      /* read trait value first - affection status or quantitative trait value */
      if (pTrait->type == DICHOTOMOUS) {
	numRet = sscanf (pLine, "%lf %n", &pPerson->ppTraitValue[i][j], &pos);
	KASSERT (numRet == 1,
		 "Failed to get affection status on line %d in file %s.\n",
		 lineNo, sPedfileName);
	if ((!isnan(pPerson->ppTraitValue[i][j])) &&
	    ((int) pPerson->ppTraitValue[i][j] !=
	     modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN]))
	  pPerson->ppTraitKnown[i][j] = TRUE;
      } else if (pTrait->type == QUANTITATIVE || pTrait->type == COMBINED) {
	numRet = sscanf (pLine, "%lf %n", &pPerson->ppOrigTraitValue[i][j], &pos);
	KASSERT (numRet == 1,
		 "Line %d in pedfile %s doesn't have enough columns (trait). Is this a post-makeped file? \n",
		 lineNo, sPedfileName);
	if ((!isnan(pPerson->ppOrigTraitValue[i][j])) &&
	    (pPerson->ppOrigTraitValue[i][j] != pTrait->unknownTraitValue)) {
	  pPerson->ppTraitKnown[i][j] = TRUE;
	  if ((pPerson->ppOrigTraitValue[i][j] != modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED]) &&
	      (pPerson->ppOrigTraitValue[i][j] != modelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED]))
	    /* Calculated the standardized quantitative trait value */
	    pPerson->ppTraitValue[i][j] = (pPerson->ppOrigTraitValue[i][j] - pTrait->sampleMean) / pTrait->sampleSD;
	  else
	    pPerson->ppTraitValue[i][j] = pPerson->ppOrigTraitValue[i][j];
	}
      } else {
	KASSERT (FALSE, "Unknown trait locus value type %d.\n", pTrait->type);
      }

      pLine = &pLine[pos];
      /* If the locus is defined as having liability class
       * we have to read the liability class data */
      if (pTrait->numLiabilityClass > 1) {
	numRet =
	  sscanf (pLine, "%d %n", &pPerson->ppLiabilityClass[i][j], &pos);
	KASSERT (numRet == 1,
		 "Line %d in pedfile %s doesn't have enough columns (LC). Is this a post-makeped file? \n",
		 lineNo, sPedfileName);
	pLine = &pLine[pos];
      }

      /* move on to next trait variable */
      j++;
    }
    /* move on to next trait locus */
    i++;
  }
  /* Next should be phenotypes for the markers 
   * */
  numMarker = i;
  while (numMarker < originalLocusList.numLocus) {
    pLocus = originalLocusList.ppLocusList[numMarker];

    /* read a pair of genotypes for the current marker */
    numRet = sscanf (pLine, "%s %s %n", a1, a2, &pos);
    if (numRet != 2) {
      /* phase of alleles might be known 
       * paternal | maternal     is assumed */
      numRet = sscanf (pLine, "%s | %s %n", a1, a2, &pos);
      pPerson->pPhasedFlag[numMarker] = 1;
    }
    KASSERT (numRet == 2,
	     "Line %d in pedfile %s doesn't have enough columns (Marker %d). Is this a post-makeped file? \n",
	     lineNo, sPedfileName, numMarker);
    /* check whether this locus is untyped */
    if (strcmp(a1, "0") == 0 || strcmp(a2, "0") == 0 || strcasecmp(a1, "X") == 0 || strcasecmp(a2, "X") == 0)
      pPerson->pTypedFlag[numMarker] = 0;
    else
      pPerson->pTypedFlag[numMarker] = 1;

    if (strcmp(a1, "0") != 0 && strcasecmp(a1, "X") != 0) {
      ret = find_allele (numMarker, a1);
      pPerson->pPhenotypeList[0][numMarker] = ret;
      KASSERT (ret >= 0,
	       "Line %d in pedfile %s contains a genotype with unknown allele %s at locus %s.\n",
	       lineNo, sPedfileName, a1, pLocus->sName);
    }
    else
      pPerson->pPhenotypeList[0][numMarker] = 0;

    if (strcmp(a2, "0") != 0 && strcasecmp(a2, "X") != 0) {
      ret = find_allele (numMarker, a2);
      pPerson->pPhenotypeList[1][numMarker] = ret;
      KASSERT (ret >= 0,
	       "Line %d in pedfile %s contains a genotype with unkown allele %s at locus %s.\n",
	       lineNo, sPedfileName, a2, pLocus->sName);
    }
    else
      pPerson->pPhenotypeList[1][numMarker] = 0;

    /* if this is X chromosome and this person is a male, then the genotype needs to be 
     * homozygous */
    if (modelOptions.sexLinked && pPerson->sex == 0
	&& pPerson->pPhenotypeList[0][numMarker] !=
	pPerson->pPhenotypeList[1][numMarker]) {
      KASSERT (0 == 1,
	       "Line %d in pedfile %s contains heterozygous genotype(%d, %d) for a male at locus %s while doing X chromosome analysis.\n",
	       lineNo, sPedfileName,
	       pPerson->pPhenotypeList[0][numMarker],
	       pPerson->pPhenotypeList[1][numMarker], pLocus->sName);

    }

    numMarker++;
    /* move the line buffer over to the next pair */
    pLine = &pLine[pos];
  }

  /* Next should be Ped: *** Per: *** 
   * We could find out whether an individual has been doubled 
   * from Per: field */
  numRet = sscanf (pLine, "Ped: %s Per: %s",
		   pPerson->pPedigree->sOriginalID, pPerson->sOriginalID);
  // If we don't find them, default to the 1st-column
  if (numRet <= 0) {
    strcpy (pPerson->pPedigree->sOriginalID, pPerson->pPedigree->sPedigreeID);
    strcpy (pPerson->sOriginalID, pPerson->sID);
  }
  KASSERT (numRet <= 0 || numRet == 2,
	   "Line %d in pedfile %s doesn't have enough columns (Old Ped, Per:). Is this a post-makeped file? \n",
	   lineNo, sPedfileName);

  return 0;
}

Pedigree *
create_pedigree (PedigreeSet * pPedigreeSet, char *sPedLabel)
{
  int num;
  Pedigree *ped;
  int oldNum;
  int i;

  oldNum = pPedigreeSet->numPedigree;
  if (pPedigreeSet->maxNumPedigree <= pPedigreeSet->numPedigree) {
    /* allocate more space */
    num = (pPedigreeSet->numPedigree + DEF_PED_MALLOC_INCREMENT);
    pPedigreeSet->ppPedigreeSet =
      (Pedigree **) realloc (pPedigreeSet->ppPedigreeSet,
			     num * sizeof (Pedigree *));
    pPedigreeSet->maxNumPedigree = num;
    pPedigreeSet->nullLikelihood = (double *)
      realloc (pPedigreeSet->nullLikelihood, num * sizeof (double));
  }
  /* increment number of pedigrees for the dataset */
  pPedigreeSet->numPedigree++;

  /* create new pedigree */
  ped = (Pedigree *) malloc (sizeof (Pedigree));
  memset (ped, 0, sizeof (Pedigree));

  pPedigreeSet->ppPedigreeSet[oldNum] = ped;
  ped->pPedigreeSet = pPedigreeSet;
  ped->pedigreeIndex = oldNum;
  /* only know this pedigree's ID so far  */
  strcpy (ped->sPedigreeID, sPedLabel);

  /* for casectrl */
  ped->pCount = (int *) malloc (originalLocusList.numLocus * sizeof (int));
  for (i = 0; i < originalLocusList.numLocus; i++) {
    ped->pCount[i] = 1;
  }

  return ped;
}

Person *
create_person (Pedigree * pPed, char *sID)
{
  int num;
  Person *pPerson;
  int oldNum;
  int size;

  oldNum = pPed->numPerson;
  if (pPed->maxNumPerson <= pPed->numPerson) {
    /* allocate more space */
    num = (pPed->numPerson + DEF_PED_MALLOC_INCREMENT);
    pPed->ppPersonList = (Person **) realloc (pPed->ppPersonList,
					      num * sizeof (Person *));
    pPed->maxNumPerson = num;
  }
  pPed->numPerson++;
  if (pPed->pPedigreeSet->maxNumPerson < pPed->numPerson)
    pPed->pPedigreeSet->maxNumPerson = pPed->numPerson;

  /* allocate space for this person */
  pPerson = (Person *) malloc (sizeof (Person));
  memset (pPerson, 0, sizeof (Person));

  pPed->ppPersonList[oldNum] = pPerson;
  pPerson->pPedigree = pPed;
  pPerson->personIndex = oldNum;

  /* allocate space for phenotype list for each locus */
  size = originalLocusList.numLocus * sizeof (int);
  pPerson->pPhenotypeList[0] = (int *) malloc (size);
  memset (pPerson->pPhenotypeList[0], 0, size);
  pPerson->pPhenotypeList[1] = (int *) malloc (size);
  memset (pPerson->pPhenotypeList[1], 0, size);
  pPerson->pPhasedFlag = (int *) malloc (size);
  memset (pPerson->pPhasedFlag, 0, size);
  pPerson->pTypedFlag = (int *) malloc (size);
  memset (pPerson->pTypedFlag, 0, size);

  /* allocate space for set recoding */
  size = originalLocusList.alleleSetLen;
  pPerson->pTransmittedAlleles[MOM] =
    (unsigned int *) malloc (sizeof (unsigned int) * size);
  pPerson->pTransmittedAlleles[DAD] =
    (unsigned int *) malloc (sizeof (unsigned int) * size);
  pPerson->pNonTransmittedAlleles[MOM] =
    (unsigned int *) malloc (sizeof (unsigned int) * size);
  pPerson->pNonTransmittedAlleles[DAD] =
    (unsigned int *) malloc (sizeof (unsigned int) * size);

  /* need to allocate some space for genotypes */
  pPerson->ppGenotypeList = (Genotype **) malloc (sizeof (Genotype *) *
						  originalLocusList.numLocus);
  memset (pPerson->ppGenotypeList, 0,
	  sizeof (Genotype *) * originalLocusList.numLocus);
  pPerson->pNumGenotype =
    (int *) malloc (sizeof (int) * originalLocusList.numLocus);
  memset (pPerson->pNumGenotype, 0,
	  sizeof (int) * originalLocusList.numLocus);
  pPerson->ppSavedGenotypeList =
    (Genotype **) malloc (sizeof (Genotype *) * originalLocusList.numLocus);
  memset (pPerson->ppSavedGenotypeList, 0,
	  sizeof (Genotype *) * originalLocusList.numLocus);
  pPerson->pSavedNumGenotype =
    (int *) malloc (sizeof (int) * originalLocusList.numLocus);
  memset (pPerson->pSavedNumGenotype, 0,
	  sizeof (int) * originalLocusList.numLocus);
  pPerson->ppProbandGenotypeList =
    (Genotype **) malloc (sizeof (Genotype *) * originalLocusList.numLocus);
  memset (pPerson->ppProbandGenotypeList, 0,
	  sizeof (Genotype *) * originalLocusList.numLocus);
  pPerson->pProbandNumGenotype =
    (int *) malloc (sizeof (int) * originalLocusList.numLocus);
  memset (pPerson->pProbandNumGenotype, 0,
	  sizeof (int) * originalLocusList.numLocus);
  pPerson->ppShadowGenotypeList =
    (Genotype **) malloc (sizeof (Genotype *) * originalLocusList.numLocus);
  memset (pPerson->ppShadowGenotypeList, 0,
	  sizeof (Genotype *) * originalLocusList.numLocus);
  pPerson->pShadowGenotypeListLen =
    (int *) malloc (sizeof (int) * originalLocusList.numLocus);
  memset (pPerson->pShadowGenotypeListLen, 0,
	  sizeof (int) * originalLocusList.numLocus);

  /* we have only read this person's ID so far */
  strcpy (pPerson->sID, sID);

  return pPerson;
}

/* This function searchs pedigrees that with the matching 
 * pedigree IDs - pedigree ID suppose to be unique */
Pedigree *
find_pedigree (PedigreeSet * pPedSet, char *sPedID)
{
  int i;

  for (i = 0; i < pPedSet->numPedigree; i++) {
    if (strcmp (pPedSet->ppPedigreeSet[i]->sPedigreeID, sPedID) == 0) {
      /* found the matching one */
      return pPedSet->ppPedigreeSet[i];
    }
  }

  /* failed to find matching one */
  return NULL;
}

/* This function searchs pedigrees that with the matching 
 * original pedigree IDs - pedigree ID suppose to be unique */
Pedigree *
find_original_pedigree (PedigreeSet * pPedSet, char *sPedID)
{
  int i;

  for (i = 0; i < pPedSet->numPedigree; i++) {
    if (strcmp (pPedSet->ppPedigreeSet[i]->sOriginalID, sPedID) == 0) {
      /* found the matching one */
      return pPedSet->ppPedigreeSet[i];
    }
  }

  /* failed to find matching one */
  return NULL;
}

/* This function searchs matching person in a pedigree 
 * Pedigree pointer is expected as input
 * it finds the person in the pedigree with the matching person ID */
Person *
find_person (Pedigree * pPed, char *sPersonID)
{
  int i;

  if (pPed == NULL)
    return NULL;

  for (i = 0; i < pPed->numPerson; i++) {
    if (strcmp (pPed->ppPersonList[i]->sID, sPersonID) == 0) {
      /* found it */
      return pPed->ppPersonList[i];
    }
  }

  /* failed to find matching one */
  return NULL;
}

/* This function searchs matching person in a pedigree 
 * Pedigree pointer is expected as input
 * it finds the person in the pedigree with the matching person ID */
Person *
find_original_person (Pedigree * pPed, char *sPersonID)
{
  int i;
  Person *pPerson;

  if (pPed == NULL)
    return NULL;

  for (i = 0; i < pPed->numPerson; i++) {
    pPerson = pPed->ppPersonList[i];
    if (strcmp (pPerson->sOriginalID, sPersonID) == 0
	&& pPerson->pParents[DAD] != NULL) {
      /* found it */
      return pPerson;
    }
  }

  /* failed to find matching one */
  return NULL;
}

/* This function searchs for another person that has the same proband ID as
 * the input person 
 */
Person *
find_loop_breaker (Person * breaker)
{
  int i;
  Person *pPerson;
  Pedigree *pPed = breaker->pPedigree;

  if (pPed == NULL)
    return NULL;

  for (i = 0; i < pPed->numPerson; i++) {
    pPerson = pPed->ppPersonList[i];
    if (pPerson != breaker && pPerson->proband == breaker->proband
	&& pPerson->pParents[DAD] != NULL) {
      /* found it */
      return pPerson;
    }
  }

  /* failed to find matching one */
  return NULL;
}

/* This function sets up pointers in the pedigree, such as
 * first child, next paternal or maternal siblings etc. 
 * We have decided to trust the post-makeped pedfile which 
 * provides the IDs of first child, next paternal or mater siblings
 * Later on, we may want to build in these information all by ourselves
 * when we want to take on the task of breaking loops as well
 *
 * This function also makes sure the same loop breaker points to
 * same sets of genotypes 
 * */
int
setup_pedigree_ptrs (Pedigree * pPed)
{
  int i;
  Person *pPerson;
  Person *pOrigPerson;

  /* nothing to do */
  if (pPed == NULL)
    return 0;

  /* go through each person to set up related pointers */
  for (i = 0; i < pPed->numPerson; i++) {
    pPerson = pPed->ppPersonList[i];
    /* set up parents links first */
    //if (strcmp (pPerson->sDadID, pPed->pPedigreeSet->sUnknownID) != 0) {
    if (strcmp (pPerson->sParentID[DAD], modelOptions.sUnknownPersonID) != 0) {
      /* dad is known */
      pPerson->pParents[DAD] = find_person (pPed, pPerson->sParentID[DAD]);
      /* mom is known */
      pPerson->pParents[MOM] = find_person (pPed, pPerson->sParentID[MOM]);
    } else {
      if (pPerson->loopBreaker == 0)
	/* this person is a true founder 
	 * doubled individual may appear to be a founder as no parents
	 * listed, but they can't possibly be founders - we never
	 * need to double true founders. they can't possibly 
	 * be involved in a loop situation */
	add_founder (pPed, pPerson);
    }

    /* set up first child link */
    if (strcmp (pPerson->sFirstChildID, modelOptions.sUnknownPersonID) != 0) {
      /* first child is known */
      pPerson->pFirstChild = find_person (pPed, pPerson->sFirstChildID);
    }
    /* set up next paternal sibling link */
    if (strcmp (pPerson->sNextSibID[DAD], modelOptions.sUnknownPersonID) != 0) {
      pPerson->pNextSib[DAD] = find_person (pPed, pPerson->sNextSibID[DAD]);
    } else
      pPerson->pNextSib[DAD] = NULL;

    /* set up next maternal sibling link */
    if (strcmp (pPerson->sNextSibID[MOM], modelOptions.sUnknownPersonID) != 0) {
      pPerson->pNextSib[MOM] = find_person (pPed, pPerson->sNextSibID[MOM]);
    } else
      pPerson->pNextSib[MOM] = NULL;

    /* If this person is a loop breaker, make sure it points to the same 
     * genotype set 
     * For a duplicated person, the parents are set to unknown
     * so we linked the genotype info to the original person who does
     * have genotype information 
     * */
    if (pPerson->loopBreaker >= 1 && pPerson->pParents[DAD] == NULL) {
      /* a duplicated person, find the original person */
      pOrigPerson = find_loop_breaker (pPerson);
      KASSERT (pOrigPerson != NULL,
	       "Can't find loop breaker's original information (original ID: %s current ID: %s).\n",
	       pPerson->sOriginalID, pPerson->sID);
      pPerson->pOriginalPerson = pOrigPerson;
      add_loopbreaker (pPed, pOrigPerson);

#if 0
      /* free space for the marker phenotype pair */
      free (pPerson->pPhenotypeList[0]);
      free (pPerson->pPhenotypeList[1]);
      free (pPerson->pPhasedFlag);
      free (pPerson->pTypedFlag);

      /* also need to free trait phenotype information !!!! - add this later */
      /* Points to the original person's marker genotype list */
      pPerson->pPhenotypeList[0] = pOrigPerson->pPhenotypeList[0];
      pPerson->pPhenotypeList[1] = pOrigPerson->pPhenotypeList[1];
      pPerson->pPhasedFlag = pOrigPerson->pPhasedFlag;
      pPerson->pTypedFlag = pOrigPerson->pTypedFlag;
#endif

      /* also need to point to the original person's trait information */
    }

  }

  return 0;
}

int
add_founder (Pedigree * pPed, Person * pPerson)
{

  pPerson->founderFlag = 1;
  if (pPed->numFounder <= pPed->maxNumFounder) {
    /* need to reallocate memory for the founder list */
    pPed->maxNumFounder += DEF_PED_MALLOC_INCREMENT;
    pPed->ppFounderList = (Person **) realloc (pPed->ppFounderList,
					       sizeof (Person *) *
					       pPed->maxNumFounder);
  }
  /* Add the person to the founder list */
  pPed->ppFounderList[pPed->numFounder] = pPerson;
  /* increment the founder list counter */
  pPed->numFounder++;

  /* update the pedigree set's record */
  if (pPed->pPedigreeSet->maxNumFounder < pPed->numFounder)
    pPed->pPedigreeSet->maxNumFounder = pPed->numFounder;

  return 0;

}

/* This function disects this pedigree into overlapping nuclear families
 * overlapping individuals are called connectors
 * Nuclear families are used in the peeling process 
 * */
int
setup_nuclear_families (Pedigree * pPed)
{
  int i, j, k;
  int *pDoneList;		/* keep track of whether we have done this child or not */
  NuclearFamily *pNucFam;
  Person *pPerson;

  /* allocate some working space 
   * This is to mark whether we are done with this child 
   * it really only applies if this person is a child */
  //  pDoneList = (int *) MALLOC ("pDoneList", sizeof (int) * pPed->numPerson);
  pDoneList = pPed->pPedigreeSet->pDonePerson;
  memset (pDoneList, 0, sizeof (int) * pPed->numPerson);

  pPed->numNuclearFamily = 0;
  /* search for nuclear families */
  for (i = 0; i < pPed->numPerson; i++) {
    /* we are only interested in nuclear families 
     * we will start with children as from children, the parents
     * can be found easily
     * while through one parent, if there are multiple marriages
     * it's a little messier to find all spouses */
    if (pPed->ppPersonList[i]->pParents[MOM] != NULL && pDoneList[i] == 0) {
      /* we have found a child that we haven't processed yet
       * Since a child can only be in one nuclear family as a child
       * create a new nuclear family */
      pNucFam = create_nuclear_family (pPed);
      /* we know the unique parents from this child, add them in
       * should also connect this nuclear family to other nuclear families
       * either parent has already been in */
      pPerson = pPed->ppPersonList[i];
      add_nuclear_family_parents (pNucFam, pPerson->pParents[DAD],
				  pPerson->pParents[MOM]);
      /* if both parents don't have parents, then this must be 
       * a founder pair */
      if (pPerson->pParents[DAD]->pParents[DAD] == NULL &&
	  pPerson->pParents[MOM]->pParents[DAD] == NULL) {
	/* create the founder pair and add to the founder pair list */
	add_founder_nuclear_family (pNucFam);
      }

      /* add this child to this nuclear family, should also connect 
       * this child to nuclear families in which this child is a parent */
      add_nuclear_family_child (pNucFam, pPerson);

      /* now we need to find all the children of this parent pairs 
       * we only need to scan through the rest of the persons, as 
       * if a child is indexed lower, we must have already processed
       * this nuclear family */
      for (k = i + 1; k < pPed->numPerson; k++) {
	if (pPed->ppPersonList[k]->pParents[MOM] == pNucFam->pParents[MOM]
	    && pPed->ppPersonList[k]->pParents[DAD] == pNucFam->pParents[DAD]) {
	  /* we found one with same pairs of parents, add it */
	  /* add this child to this nuclear family, should also connect 
	   * this child to nuclear families in which this child is a parent */
	  add_nuclear_family_child (pNucFam, pPed->ppPersonList[k]);
	  /* mark this child as done */
	  pDoneList[k] = 1;
	}
      }				/* end of loop searching other children with same parent pairs */
      /* mark this child as done 
       * it shouldn't matter, as we are going by the order, this child
       * will not get processed again anyway */
      pDoneList[i] = 1;

      /* set up peeling related data */
      if (pNucFam->pParents[MOM] == pPed->pPeelingProband ||
	  pNucFam->pParents[DAD] == pPed->pPeelingProband) {
	/* set up the first family we are going to do the peeling on */
	pPed->pPeelingNuclearFamily = pNucFam;
	pPed->peelingDirection = PEDIGREE_UP;
      }
      for (j = 0; j < pNucFam->numChildren; j++) {
	if (pNucFam->ppChildrenList[j] == pPed->pPeelingProband) {
	  pPed->pPeelingNuclearFamily = pNucFam;
	  pPed->peelingDirection = PEDIGREE_DOWN;
	}
      }
    }				/* end of processing any child */
  }				/* end of loop of looking at every person in this pedigree */

  return 0;
}

/* add mating nodes to person's structure */
int
add_spouse (Person * pPerson, Person * pSpouse)
{
  int oldNum = pPerson->numSpouse;

  if (oldNum >= pPerson->maxNumSpouse) {
    pPerson->maxNumSpouse += 1;
    /* need to reallocate (or allocate for the first time) */
    pPerson->ppSpouseList = (Person **) realloc (pPerson->ppSpouseList,
						 pPerson->maxNumSpouse *
						 sizeof (Person *));

  }
  /* add the spouse in */
  pPerson->ppSpouseList[oldNum] = pSpouse;
  /* increase the counter */
  pPerson->numSpouse++;

  return 0;
}

/* Add the nuclear family to Person's structure */
int
add_person_nuclear_family (Person * pPerson, NuclearFamily * pNucFam)
{
  int oldNum;

  oldNum = pPerson->numNuclearFamily;

  if (oldNum >= pPerson->maxNumNuclearFamily) {
    /* need to reallocate (or allocate for the first time) */
    pPerson->maxNumNuclearFamily += DEF_PED_MALLOC_INCREMENT;
    pPerson->ppNuclearFamilyList =
      (NuclearFamily **) realloc (pPerson->ppNuclearFamilyList,
				  pPerson->maxNumNuclearFamily *
				  sizeof (NuclearFamily *));

  }
  /* add the nuclear family in the list */
  pPerson->ppNuclearFamilyList[oldNum] = pNucFam;
  /* increase the counter */
  pPerson->numNuclearFamily++;

  return 0;
}

/* add the parents in for the nuclear family and also 
 * link them to existing nuclear families if they are connected to 
 * them */
int
add_nuclear_family_parents (NuclearFamily * pNucFam,
			    Person * pParent1, Person * pParent2)
{
  int i, j;
  Pedigree *pPed = pNucFam->pPedigree;
  NuclearFamily *pTempNucFam;

  /* first add the nuclear family to person's structure */
  add_person_nuclear_family (pParent1, pNucFam);
  add_person_nuclear_family (pParent2, pNucFam);

  /* add mom into dad's spouse list */
  add_spouse (pParent1, pParent2);
  /* add dad into mom's spouse list */
  add_spouse (pParent2, pParent1);

  /* add the parents in for the nuclear family */
  pNucFam->pParents[DAD] = pParent1;
  pNucFam->pParents[MOM] = pParent2;

  /* now find possible existing nuclear families that either parent is 
   * already in and link them */
  for (i = 0; i < pPed->numNuclearFamily; i++) {
    pTempNucFam = pPed->ppNuclearFamilyList[i];
    /* make sure we are not looking at the same nuclear family */
    if (pTempNucFam == pNucFam)
      /* it's likely the last one is the same, though it doesn't matter */
      continue;

    /* check whether dad is in the other nuclear family as a dad */
    if (pParent1 == pTempNucFam->pParents[DAD]) {
      /* multiple marriage with this dad 
       * in terms of peeling, it's considered as down direction 
       * add to the top of the linklist for the down families 
       * both families consider each other as down families */
      create_nuclear_family_connector (pParent1, pTempNucFam,
				       PEDIGREE_DOWN, pNucFam);
      create_nuclear_family_connector (pParent1, pNucFam, PEDIGREE_DOWN,
				       pTempNucFam);
    }
    /* check whether dad is in the other nuclear family as a child */
    for (j = 0; j < pTempNucFam->numChildren; j++) {
      if (pParent1 == pTempNucFam->ppChildrenList[j]) {
	/* the dad is also a child of the other nuclear family pTempNucFam 
	 * pTempNucFam is considered as a UP family for pNucFam 
	 * pNucFam is considered as a DOWN family for pTempNucFam */
	create_nuclear_family_connector (pParent1, pNucFam, PEDIGREE_UP,
					 pTempNucFam);
	create_nuclear_family_connector (pParent1, pTempNucFam,
					 PEDIGREE_DOWN, pNucFam);

	/* no need to search any more as a child can only in one nuclear family
	 * as a child */
	break;
      }
    }

    /* repeat the same process for mom */
    /* check whether mom is in the other nuclear family as a mom */
    if (pParent2 == pTempNucFam->pParents[MOM]) {
      /* multiple marriage with this mom 
       * in terms of peeling, it's considered as down direction 
       * add to the top of the linklist for the down families */
      create_nuclear_family_connector (pParent2, pTempNucFam,
				       PEDIGREE_DOWN, pNucFam);
      create_nuclear_family_connector (pParent2, pNucFam, PEDIGREE_DOWN,
				       pTempNucFam);
    }
    /* check whether mom is in the other nuclear family as a child */
    for (j = 0; j < pTempNucFam->numChildren; j++) {
      if (pParent2 == pTempNucFam->ppChildrenList[j]) {
	/* the mom is also a child of the other nuclear family pTempNucFam 
	 * pTempNucFam is considered as a UP family for pNucFam 
	 * pNucFam is considered as a DOWN family for pTempNucFam */
	create_nuclear_family_connector (pParent2, pNucFam,
					 PEDIGREE_UP, pTempNucFam);
	create_nuclear_family_connector (pParent2, pTempNucFam,
					 PEDIGREE_DOWN, pNucFam);

	/* no need to search any more as a child can only in one nuclear family
	 * as a child */
	break;
      }
    }
  }

  return 0;
}

int
add_nuclear_family_child (NuclearFamily * pNucFam, Person * pChild)
{
  //  Person **ppNew;
  int i;
  Pedigree *pPed = pNucFam->pPedigree;
  NuclearFamily *pTempNucFam;

  /* first add the nuclear family to person's structure */
  add_person_nuclear_family (pChild, pNucFam);

  if (pNucFam->numChildren >= pNucFam->maxNumChildren) {
    /* need to allocate new space */
    pNucFam->maxNumChildren += DEF_PED_MALLOC_INCREMENT;
#if 0
    ppNew = (Person **) malloc (sizeof (Person *) * pNucFam->maxNumChildren);
    /* copy old stuff over */
    for (i = 0; i < pNucFam->numChildren; i++)
      ppNew[i] = pNucFam->ppChildrenList[i];
    /* free the old space */
    free (pNucFam->ppChildrenList);
    /* take the new space */
    pNucFam->ppChildrenList = ppNew;
#endif
    pNucFam->ppChildrenList = (Person **) realloc (pNucFam->ppChildrenList,
						   sizeof (Person *) *
						   pNucFam->maxNumChildren);
  }

  /* add the child */
  pNucFam->ppChildrenList[pNucFam->numChildren] = pChild;
  /* increase the counter */
  pNucFam->numChildren++;

  /* Check for connected nuclear families 
   * It can't be a child in another family. It can only be a parent in
   * another nuclear family or more families */
  for (i = 0; i < pPed->numNuclearFamily; i++) {
    /* it can only be a mom or dad depends on the sex, but we really
     * don't care */
    pTempNucFam = pPed->ppNuclearFamilyList[i];
    /* make sure we are not looking at the same nuclear family */
    if (pTempNucFam == pNucFam)
      /* it's likely the last one is the same, though it doesn't matter */
      continue;
    if (pChild == pTempNucFam->pParents[DAD]
	|| pChild == pTempNucFam->pParents[MOM]) {
      /* found one connected one. It's possible there are more 
       * child in pNucFam, but parent in pTempNucFam
       * pTempNucFam is considered as a DOWN family to pNucFam
       * pNucFam is considered as a UP family to pTempNucFam*/
      create_nuclear_family_connector (pChild, pNucFam, PEDIGREE_DOWN,
				       pTempNucFam);
      create_nuclear_family_connector (pChild, pTempNucFam, PEDIGREE_UP,
				       pNucFam);
    }

  }

  return 0;
}

/* create a connector for one nuclear family
 * connectted person is pPerson
 * will be added to pUpFamilies or pDownFamilies depend on the "direction"
 * connect pNucFam2 to pNucFam1. pNucFam1 is the one we are setting
 * the connector */
NuclearFamilyConnector *
create_nuclear_family_connector (Person * pPerson,	/* connected person */
				 NuclearFamily * pNucFam1,	/* the family we are setting up */
				 int direction,	/* is Fam2 considered as UP or DOWN to Fam1 */
				 NuclearFamily * pNucFam2)
{
  NuclearFamilyConnector *pNucFamConnector;

  /* allocate space */
  pNucFamConnector =
    (NuclearFamilyConnector *) malloc (sizeof (NuclearFamilyConnector));

  if (direction == PEDIGREE_DOWN) {
    /* set up the down connector */
    pNucFamConnector->pConnectedPerson = pPerson;
    pNucFamConnector->pConnectedNuclearFamily = pNucFam2;
    pNucFamConnector->pNextConnector = pNucFam1->pDownConnectors;
    pNucFam1->pDownConnectors = pNucFamConnector;
  } else {
    /* set up the up connector  */
    pNucFamConnector->pConnectedPerson = pPerson;
    pNucFamConnector->pConnectedNuclearFamily = pNucFam2;
    pNucFamConnector->pNextConnector = pNucFam1->pUpConnectors;
    pNucFam1->pUpConnectors = pNucFamConnector;
  }
  return 0;
}

/* allocate space for a new nuclear family 
 * This function does nothing else than allocate space */
NuclearFamily *
create_nuclear_family (Pedigree * pPed)
{
  NuclearFamily *pNew;
  int oldNumNucFam;

  /* Get the current number of nuclear families */
  oldNumNucFam = pPed->numNuclearFamily;

  /* Check whether we need to allocate new block of space of list */
  if (oldNumNucFam >= pPed->maxNuclearFamily) {
    /* actually it should be never over it, otherwise we are in trouble
     * already 
     * Ok, we have used up all the space, need to reallocate a bigger
     * array and move old stuff over and free the old space */
    pPed->maxNuclearFamily += DEF_PED_MALLOC_INCREMENT;
    pPed->ppNuclearFamilyList =
      (NuclearFamily **) realloc (pPed->ppNuclearFamilyList,
				  sizeof (NuclearFamily *) *
				  pPed->maxNuclearFamily);
  }

  /* allocate space for the actual nuclear family */
  pNew = (NuclearFamily *) malloc (sizeof (NuclearFamily));
  memset (pNew, 0, sizeof (NuclearFamily));
#if 0
  pNew->ppParentalPair = (ParentalPair **) malloc (sizeof (ParentalPair *) *
						   originalLocusList.
						   numLocus);
  memset (pNew->ppParentalPair, 0,
	  sizeof (ParentalPair *) * originalLocusList.numLocus);
  pNew->pNumParentalPair =
    (int *) malloc (sizeof (int) * originalLocusList.numLocus);
#endif

  pPed->ppNuclearFamilyList[oldNumNucFam] = pNew;
  /* link this back to pedigree structure */
  pNew->pPedigree = pPed;
  /* set the index of this nuclear family */
  pNew->nuclearFamilyIndex = oldNumNucFam;
  /* increment number of nuclear families */
  pPed->numNuclearFamily++;

  return pNew;
}

/* This function prints out nuclear families within the given pedigree
 * this should help debug */
int
print_nuclear_family (FILE * fp, Pedigree * pPed)
{
  int i, j;
  NuclearFamily *pNucFam;
  NuclearFamilyConnector *pNucFamConnector;

  KLOG (LOGPEDFILE, LOGDEBUG, "Pedigree %s Nuclear Families (%d): \n",
	pPed->sPedigreeID, pPed->numNuclearFamily);
  for (i = 0; i < pPed->numNuclearFamily; i++) {
    pNucFam = pPed->ppNuclearFamilyList[i];
    KLOG (LOGPEDFILE, LOGDEBUG, "  %3d: Dad - %s Mom - %s \n",
	  i, pNucFam->pParents[DAD]->sID, pNucFam->pParents[MOM]->sID);
    KLOG (LOGPEDFILE, LOGDEBUG, "       Children - ");
    for (j = 0; j < pNucFam->numChildren; j++)
      KLOG (LOGPEDFILE, LOGDEBUG, "%s  ", pNucFam->ppChildrenList[j]->sID);
    KLOG (LOGPEDFILE, LOGDEBUG, "\n");

    /* now print out connector information */
    pNucFamConnector = pNucFam->pUpConnectors;
    KLOG (LOGPEDFILE, LOGDEBUG, "       Up Connectors:\n");
    while (pNucFamConnector != NULL) {
      KLOG (LOGPEDFILE, LOGDEBUG,
	    "\t Nuclear Family Index %3d, Connected Individual %s.\n",
	    pNucFamConnector->pConnectedNuclearFamily->
	    nuclearFamilyIndex, pNucFamConnector->pConnectedPerson->sID);
      pNucFamConnector = pNucFamConnector->pNextConnector;
    }
    pNucFamConnector = pNucFam->pDownConnectors;
    KLOG (LOGPEDFILE, LOGDEBUG, "       Down Connectors:\n");
    while (pNucFamConnector != NULL) {
      KLOG (LOGPEDFILE, LOGDEBUG,
	    "\t Nuclear Family Index %3d, Connected Individual %s.\n",
	    pNucFamConnector->pConnectedNuclearFamily->
	    nuclearFamilyIndex, pNucFamConnector->pConnectedPerson->sID);
      pNucFamConnector = pNucFamConnector->pNextConnector;
    }

  }

  /* print out top founder nuclear families */
  KLOG (LOGPEDFILE, LOGDEBUG, "    Founder Nuclear Families:\n");
  KLOG (LOGPEDFILE, LOGDEBUG, "      ");
  for (i = 0; i < pPed->numFounderNuclearFamily; i++)
    KLOG (LOGPEDFILE, LOGDEBUG, "%d ",
	  pPed->ppFounderNuclearFamilyList[i]->nuclearFamilyIndex);
  KLOG (LOGPEDFILE, LOGDEBUG, "\n");

  return 0;
}

/* founder nuclear families are those with both parents as founders 
 * this function adds the nuclear family to the list */
int
add_founder_nuclear_family (NuclearFamily * pNucFam)
{
  Pedigree *pPed = pNucFam->pPedigree;
  NuclearFamily **ppNew;
  int i;

  if (pPed->numFounderNuclearFamily >= pPed->maxNumFounderNuclearFamily) {
    /* need to reallocate space for the list */
    pPed->maxNumFounderNuclearFamily += DEF_PED_MALLOC_INCREMENT;
    ppNew = (NuclearFamily **) malloc (sizeof (NuclearFamily *) *
				       pPed->maxNumFounderNuclearFamily);
    /* copy over the old list */
    for (i = 0; i < pPed->numFounderNuclearFamily; i++)
      ppNew[i] = pPed->ppFounderNuclearFamilyList[i];
    /* free the old list */
    free (pPed->ppFounderNuclearFamilyList);
    /* assign the new list */
    pPed->ppFounderNuclearFamilyList = ppNew;
  }

  /* add the new founder nuclear family into the list */
  pPed->ppFounderNuclearFamilyList[pPed->numFounderNuclearFamily] = pNucFam;
  /* increment the counter */
  pPed->numFounderNuclearFamily++;

  return 0;
}

/* this function sets up number of loops and number of loop breakers 
 * Proband field of greater than 1 indicates this person is in a loop
 * multiple people with the same proband field value are the doubled, 
 * tripled ... individuals. 
 * if there are only two individuals with the same proband field value,
 * then this doubled individual has broken one loop 
 * if there are three individuals with the same proband field value, 
 * then this tripled individual has broken two loops */
int
setup_loop_counts (Pedigree * pPed)
{
  int i;
  int *pCount;
  Person *pPerson;
  int numTotal = 0;

  if (pPed->loopFlag == 0)
    /* no loop in this pedigree, don't even count anything */
    return 0;
  //  pPed->numLoopBreaker = 0;
  /* allocate some work space to count number of loops and number 
   * of loop breakers */
  pCount = pPed->pPedigreeSet->pDonePerson;
  /* initialize this workspace */
  memset (pCount, 0, sizeof (int) * pPed->numPerson);

  for (i = 0; i < pPed->numPerson; i++) {
    pPerson = pPed->ppPersonList[i];
    /* Assuming the proband column for loop breakers start at 2, then 3 etc. */
    if (pPerson->proband > 1 && pPerson->proband <= pPed->numPerson) {
      if (pCount[pPerson->proband] == 0)
	/* increase the count for this proband */
	pCount[pPerson->proband]++;
      numTotal++;
    }
  }
  /* get number of loops */
  pPed->numLoop = numTotal - pPed->numLoopBreaker;
  return 0;
}

/* NOT USED 
 * return indicator whether this pedigree set needs transmission matrices
 * There are four conditions, this function returns 1 indicating the need
 * for transmission matrices
 * 1. loop present in any pedigree
 * 2. any proband who is not founder (i.e. has parents)
 * 3. any proband who is founder, but spouse(s) is not a founder
 * 4. any pedigree has two sets of founders 
 * */
int
need_transmission_matrices (PedigreeSet * pPedSet)
{
  int i, j;
  Pedigree *pPed;
  Person *pProband;

  /* NuclearFamily *pNucFam; */

  /* condition 1 - at least one pedigree has loop */
  if (pPedSet->loopFlag == 1)
    return 1;

  /* go through each pedigree */
  for (i = 0; i < pPedSet->numPedigree; i++) {
    pPed = pPedSet->ppPedigreeSet[i];
    pProband = pPed->pPeelingProband;
    /* condition 2: proband who is not founder */
    if (pProband->pParents[DAD] != NULL) {
      return 1;
    }

    /* condition 3: proband is founder, but proband's spouse is not founder 
     * at this point, we know proband has to be a founder, otherwise
     * the function would have already returned */
    /*
       for(j=0; j < pProband->numNuclearFamily; j++)
       {
       pNucFam = pProband->ppNuclearFamilyList[j];
       if( (pNucFam->pDad == pProband && pNucFam->pMom->pDad != NULL) ||
       (pNucFam->pMom == pProband && pNucFam->pDad->pDad != NULL) )
       {
       return 1;
       }
       }
     */
    for (j = 0; j < pProband->numSpouse; j++) {
      if (pProband->ppSpouseList[j]->pParents[DAD] != NULL)
	return 1;
    }

    /* condition 4: two sets of founder */
    if (pPed->numFounderNuclearFamily > 1) {
      return 1;
    }
  }

  /* if the program reaches here, return false */
  return 0;

}

/** 

  Read count file (CF) - counts of different pedigree inheritance
  patterns.

  Builds a vector of pointers to existing pedigrees (in pPedigreeSet)
  for the header columns. For each row of counts by marker, updates
  each pedigree's pCount[i], where i is the locus index.

  Returns 0.

*/
int
read_ccfile (char *ccFileName, PedigreeSet * pPedigreeSet)
{
  FILE *fpFile = NULL;
  int lineNo = 0, pos, i, j, numRet, headerSeen = FALSE, numPed = 0, maxNumPed = 0;
  char *pLine;
  char sCurrPedLabel[MAX_PED_LABEL_LEN];
  char sCurrMarkerLabel[MAX_PED_LABEL_LEN];
  Pedigree **pPedSet = NULL; // Vector of pedigree pointers from header columns

  fpFile = fopen (ccFileName, "r");
  KASSERT (fpFile != NULL,
	   "Can't open case control count file %s for read. Exiting!",
	   ccFileName);
  while (fgetlongs (&flexBuffer, &flexBufferSize, fpFile) != NULL) {
    lineNo++;
    if (is_line_blank_or_comment (flexBuffer))
      continue;
    if (!headerSeen) {
      headerSeen = TRUE;
      /* The first header column should be the marker - ignore it */
      numRet = sscanf (flexBuffer, "%s %n", sCurrMarkerLabel, &pos);
      pLine = &flexBuffer[pos];
      while ((numRet = sscanf (pLine, "%s %n", sCurrPedLabel, &pos)) > 0) {
	pLine = &pLine[pos];
	Pedigree *pPed;
	pPed = find_original_pedigree (pPedigreeSet, sCurrPedLabel);
	KASSERT (pPed != NULL,
		 "Can't find pedigree %s in pedigree file, but we got pattern count.\n",
		 sCurrPedLabel);
	if (numPed + 1 > maxNumPed) {
	  maxNumPed += 6;
	  pPedSet = realloc (pPedSet, sizeof (Pedigree *) * maxNumPed);
	}
	pPedSet[numPed++] = pPed;
      }
      continue;
    }

    /* Read the count information for each SNP */
    numRet = sscanf (flexBuffer, "%s %n", sCurrMarkerLabel, &pos);
    KASSERT (numRet == 1,
	     "Can't get marker label from line no %d in count file %s.\n",
	     lineNo, ccFileName);
    pLine = &flexBuffer[pos];
    /* Find the locus */
    j = find_locus (&originalLocusList, sCurrMarkerLabel);
    KASSERT (j >= 0, "Can't get marker %s information.\n", sCurrMarkerLabel);
    i = 0;
    while (i < numPed) {
      numRet = sscanf (pLine, "%d %n", &(pPedSet[i]->pCount[j]), &pos);
      KASSERT (numRet == 1,
	       "Can't get count information for pedigree %s at locus %s.\n",
	       pPedSet[i]->sPedigreeID, sCurrMarkerLabel);
      pLine = &pLine[pos];
      i++;
    }
  }

  return 0;
}

/* add this person to the loop breaker list if not currently in */
void
add_loopbreaker (Pedigree * pPed, Person * pPerson)
{
  int num;
  int i;
  Person *pBreaker;

  num = pPed->numLoopBreaker;
  i = 0;
  while (i < num) {
    pBreaker = pPed->loopBreakerList[i];
    if (pBreaker == pPerson)
      return;
    i++;
  }
  pPed->loopBreakerList = (Person **) realloc (pPed->loopBreakerList,
					       (num + 1) * sizeof (Person *));
  pPed->loopBreakerList[num] = pPerson;
  pPed->numLoopBreaker++;
}


unsigned long long ancestryVectors[64];
unsigned long long aVDoneFlag = 1ULL << 63;
/**

  Return an ancestry vector for the specified individual.  Recursively
  compose it if it's not already defined. The ancestry vector is a bit
  vector with bits set for each ancestral personIndex. We are
  currently limited, therefore to 63 individuals per pedigree (we use
  the 64th bit to indicate if the vector is complete).  If, in the
  process of composing a vector, we find that the Mom's vector and
  Dad's vector have common bits set, then they've got a common
  ancestor, and that constitutes a loop.

*/
unsigned long long
getAncestryVector(int personIndex, Pedigree *pPed) {
  int momIndex = 0, dadIndex = 0, i;
  Person *pPerson, *pDad, *pMom;
  unsigned long long commonAncestors;
  
  /* If we already have the vector, return it. */
  if (ancestryVectors[personIndex] & aVDoneFlag)
    return ancestryVectors[personIndex];
  /* Since we don't have it, compose it from Mom's and Dad's. */
  ancestryVectors[personIndex] = 0;
  pPerson = pPed->ppPersonList[personIndex];
  if ((pMom = pPerson->pParents[MOM]) != NULL) {
    momIndex = pMom->personIndex;
    if (!(ancestryVectors[momIndex] & aVDoneFlag))
      getAncestryVector(momIndex, pPed);
    ancestryVectors[personIndex] = ancestryVectors[momIndex];
    ancestryVectors[personIndex] |= 1ULL << momIndex;
  }
  if ((pDad = pPerson->pParents[DAD]) != NULL) {
    dadIndex = pDad->personIndex;
    if (!(ancestryVectors[dadIndex] & aVDoneFlag))
      getAncestryVector(dadIndex, pPed);
    if ((commonAncestors = ((ancestryVectors[personIndex] & ancestryVectors[dadIndex]) & (~aVDoneFlag))) > 1) {
      pPed->currentLoopFlag = 1;
      fprintf (stderr, "Pedigree %s, individual %s has Mom %s and Dad %s with common ancestor(s) ",
	       pPed->sPedigreeID, pPerson->sID, pMom->sID, pDad->sID);
      for (i=0; i<64; i++) {
	if (commonAncestors & 1)
	  fprintf (stderr, "%s ", pPed->ppPersonList[i]->sID);
	commonAncestors = commonAncestors >> 1;
      }
      fprintf (stderr, "\n");
    }
    ancestryVectors[personIndex] |= ancestryVectors[dadIndex];
    ancestryVectors[personIndex] |= 1ULL << dadIndex;
  }
  /* Indicate that we have it now. */
  ancestryVectors[personIndex] |= aVDoneFlag;
  return ancestryVectors[personIndex];
}

/**

  Find only cosanguinity loops. This was never developed beyond a limitation
  of 63 individuals because it only detects cosanguinity loops, and we need
  to know about marriage loops as well. It is a kind of slick approach, though.

*/
void
check_for_common_ancestor (Pedigree *pPed) {
  int i;

  for (i=0; i<64; i++)
    ancestryVectors[i] = 0;

  /* For each parental pair, look for any common ancestors. Do it fast by
   just populating an ancestry bit matrix and watching to see if we've
   already set an ancestry bit earlier for this person. */

  if (pPed->numPerson > 63) {
    fprintf(stderr, "Current loop checking doesn't handle pedigrees of more than 63 individuals\n");
    return;
  }
  for (i = 0; i < pPed->numPerson; i++)
    getAncestryVector(i, pPed);
  return;
}


/**

  Find both marriage and cosanguinity loops.

  This approach relies upon the fact that all vertices in a loop must have
  at least two edges connecting them to other vertices.

  1. Convert the pedigree to a graph using individuals and "productive
     relationships" as vertices, and links between individuals and "productive
     relationships" as edges. For example if we have parents 1 & 2 with 
     children 3 & 4, then there are 5 vertices, one of which is the 
     "productive relationship" (we'll number it 5), and there are 4 edges,
     1-to-5 , 2-to-5, 3-to-5, and 4-to-5.
  2. Eliminate all singly-connected vertices. If none can be eliminated, then 
     you're done, and there is at least one loop.
  3. Recount vertex connections. If none are left, you're done and there are no 
     loops, otherwise go to 2.

*/
void
check_for_loop (Pedigree *pPed) {

  int tuple[256][3]; ///< Self, Mom, Dad
  Person *pPerson;
  int i, j, numPRs = 0, firstPR, secondPR, potentialPR, remainingPersons, removedSome,
    referenceCount;

  /// Copy individual and parent indexes to a structure we can destroy.
  for (i = 0; i < pPed->numPerson; i++) {
    pPerson = pPed->ppPersonList[i];
    tuple[i][0] = atoi(pPerson->sID);
    tuple[i][1] = atoi(pPerson->sParentID[MOM]);
    tuple[i][2] = atoi(pPerson->sParentID[DAD]);
    potentialPR = tuple[i][1] * 1000 + tuple[i][2];
    for (j = pPed->numPerson; j < (pPed->numPerson + numPRs); j++)
      if (tuple[j][0] == potentialPR) break;
    if (j == (pPed->numPerson + numPRs)) {
      tuple[j][0] = potentialPR; tuple[j][1] = tuple[i][1]; tuple[j][2] = tuple[i][2];
      //            printf("Adding %d %d %d at %d\n", tuple[j][0], tuple[j][1], tuple[j][2], j);
      numPRs++;
    }
  }
  //    printf("numPRs now %d\n", numPRs);

  remainingPersons = pPed->numPerson; removedSome = TRUE;
  while (remainingPersons > 0 && removedSome) {
    //        for (i = 0; i < (pPed->numPerson + numPRs); i++)
    //    printf("%d: %d %d %d\n", i, tuple[i][0], tuple[i][1], tuple[i][2]);
    removedSome = FALSE;
    // Try to remove individuals
    for (i = 0; i < pPed->numPerson; i++) { // Check everyone...
      if (tuple[i][0] == 0) continue; // ...still in the pedigree
      firstPR = 0; secondPR = 0; // Start with no links
      // Now look for the individual as a parent...
      for (j = 0; j < (pPed->numPerson + numPRs); j++) { // ...of anyone...
	if (tuple[j][0] == 0) continue; // ...still in the pedigree
	// Credit if subject is a parent...
	if ((tuple[i][0] == tuple[j][1]) || (tuple[i][0] == tuple[j][2])) {
	  potentialPR = tuple[j][1] * 1000 + tuple[j][2];
	  if (firstPR == 0)
	    firstPR = potentialPR;
	  else
	    if (potentialPR != firstPR) {
	      secondPR = potentialPR;
	      break;
	    }
	}
	// Credit if parental relationship still around
	if (j >= pPed->numPerson && tuple[i][1] == tuple[j][1] && tuple[i][2] == tuple[j][2]) {
	  potentialPR = tuple[j][1] * 1000 + tuple[j][2];
	  if (firstPR == 0)
	    firstPR = potentialPR;
	  else
	    if (potentialPR != firstPR) {
	      secondPR = potentialPR;
	      break;
	    }
	}
      }
      if (secondPR == 0) {
	//	printf("Removing %d, had only %d\n", tuple[i][0], firstPR);
	tuple[i][0] = 0; tuple[i][1] = 0; tuple[i][2] = 0;
	remainingPersons--;
	removedSome = TRUE;
      }//  else
      //	printf("%d has %d and %d\n", tuple[i][0], firstPR, secondPR);
    }
    // Try to remove pRs
    for (i = pPed->numPerson; i < (pPed->numPerson + numPRs); i++) { // Check all PRs..
      if (tuple[i][0] == 0) continue; // ...still in the pedigree
      // Now look for references to the PRs
      referenceCount = 0;
      for (j = 0; j < pPed->numPerson; j++) { // ...by anyone...
	if (tuple[j][0] == 0) continue; // ...still in the pedigree
	if (tuple[j][0] == tuple[i][1] || tuple[j][0] == tuple[i][2] ||
	    (tuple[j][1] == tuple[i][1] && tuple[j][2] == tuple[i][2])) {
	  referenceCount++;
	  //	  printf("PR %d got one from %d\n", i, tuple[j][0]);
	}
      }
      if (referenceCount < 2) {
	//	printf("Removing %d\n", tuple[i][0]);
	tuple[i][0] = 0; tuple[i][1] = 0; tuple[i][2] = 0;
	removedSome = TRUE;
      }
    }
  }
  if (remainingPersons > 0) {
    pPed->currentLoopFlag = 1;
    fprintf(stdout, "Pedigree %s, loop(s) found involving individuals ", pPed->sPedigreeID);
    for (i = 0; i< pPed->numPerson; i++)
      if (tuple[i][0] != 0)
	printf("%d ", tuple[i][0]);
    printf("\n");
  }
}


void getPedigreeTraitRange (PedigreeSet *pPedigreeSet, double *min , double *max)
{
  int va, vb;
  struct Person *person;

  *min = 999999999.00;
  *max = -999999999.00;
  for (va = 0; va < pPedigreeSet->numPedigree; va++)
    for (vb = 0 ; vb < pPedigreeSet->ppPedigreeSet[va]->numPerson; vb++) {
      person = pPedigreeSet->ppPedigreeSet[va]->ppPersonList[vb];
      if (person->ppOrigTraitValue[0][0] == modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN])
	continue;
      if (*min > person->ppOrigTraitValue[0][0])
	*min = person->ppOrigTraitValue[0][0];
      if (*max < person->ppOrigTraitValue[0][0])
	*max = person->ppOrigTraitValue[0][0];
    }
  return;
}