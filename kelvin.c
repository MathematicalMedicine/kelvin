/**********************************************************************
 * Kelvin - Linkage and Linkage Disequalibrium Analysis Program
 * Yungui Huang
 * Polynomial features - Hongling Wang
 * config.c and error logging modules - Alberto Maria Segre
 * Regex code - Nathan Burnette
 * RADSMM storage code - Martin Milder
 * 
 * Copyright 2007, Columbus Children's Research Institute.  All rights reserved.
 * Permission is hereby given to use and reproduce this software 
 * for non-profit educational purposes only.
 **********************************************************************/
#include "kelvin.h"
#include "likelihood.h"

/* Some default global values. */
char markerfile[KMAXFILENAMELEN + 1] = "markers.dat";
char mapfile[KMAXFILENAMELEN + 1] = "mapfile.dat";
char pedfile[KMAXFILENAMELEN + 1] = "pedfile.dat";
char datafile[KMAXFILENAMELEN + 1] = "datafile.dat";
char loopsfile[KMAXFILENAMELEN + 1] = "loops.dat";
char outfile[KMAXFILENAMELEN + 1] = "lods.out";
char avghetfile[KMAXFILENAMELEN + 1] = "avghet.out";
char avghomofile[KMAXFILENAMELEN + 1] = "avghomo.out";
char pplfile[KMAXFILENAMELEN + 1] = "ppl.out";
char ldPPLfile[KMAXFILENAMELEN + 1] = "ldppl.out";
FILE *fpHet = NULL;		/* average HET LR file */
//FILE *fpHomo = NULL;          /* average HOMO LR file */
FILE *fpPPL = NULL;		/* PPL output file */

/* Model datastructures. modelOptions is defined in the pedigree library. */
ModelType modelType;
ModelRange modelRange;
ModelOptions modelOptions;


/* number of D primes 
 * if there are more than 2 alleles in the marker/trait, number of D primes
 * and D prime ranges are assumed to be the same to reduce complexity 
 * for initial phase of this project */
int num_of_d_prime;
double *d_prime;
int num_of_theta;

/* three dimensional array for the two point summary results *
 * first dimension is the D prime, for LE, D prime=0 with just one element
 * in this dimension 
 * second dimension is theta values 
 * third dimension is marker allele frequency, for LE, only one element in this dimension */
SUMMARY_STAT ***tp_result;
/* two dimensional array per (dprime, theta) 
 * this will be used to calculate PPL
 */
/* storage for the NULL likelihood for the multipoint calculation under polynomial */
double ***likelihoodDT = NULL;
double *****likelihoodQT = NULL;

/* for multipoint, we use genetic map positions on a chromosome */
double *map_position;
int num_of_map_position;
/* one dimensional array, indexing by map position 
 * for multipoint, we don't know how to incorporate LD in yet 
 * This map could be sex specific map or sex averaged map 
 * For two point, we don't have to distinguish sex specific/avearge 
 * as we use theta relative to marker during analysis and after analysis
 * (result) */
SUMMARY_STAT *mp_result;
int numPositions;

XMission *nullMatrix;
XMission *altMatrix;

LambdaCell *pLambdaCell = NULL;
int prevNumDPrime = 0;
int loopMarkerFreqFlag = 0;



/**********************************************************************
 * Usage:
 *    kelvin [-s][-c] config.dat
 *
 * The config.dat file gives information about the specific linkage
 * analysis run. All information about, e.g., which markers to use,
 * what outputs to calculate, and so on, are stored in this
 * configuration file.
 *
 * -s : run serially
 * -c : restart from checkpoint
 **********************************************************************/
int
main (int argc, char *argv[])
{
  int i, j;
  int serial = FALSE;
  char configfile[KMAXFILENAMELEN] = "";
  char ckptfile[KMAXFILENAMELEN] = "";
  int breakFlag = FALSE;
  double alphaV, alphaV2;
  int loc1, loc2;

  PedigreeSet pedigreeSet;	/* Pedigrees. */
#if FALSE
  RADSMM_header_type header;	/* RADSMM */
#endif

  /* Start GAW */
  Pedigree *pPedigree;
  double pen_DD, pen_Dd, pen_dd;
  double mean_DD, mean_Dd, mean_dd;
  double SD_DD, SD_Dd, SD_dd;
  double gfreq;			/* disease gene frequency */
  double theta[2];		/* theta */
  int penIdx, liabIdx, gfreqInd, thetaInd;
  int paramIdx = -1;
  int dprimeIdx;
  //double likelihood_null, likelihood_alternative;
  double log10_likelihood_null, log10_likelihood_alternative;
  double likelihood_ratio;
  double log10_likelihood_ratio;
  Locus *pLocus;
  Locus *pLocus1, *pLocus2;
  Trait *pTrait;
  int pedIdx;
  double homoLR, hetLR;
  double adjustedHetLR = 0;
  double log10HetLR;
  double max;
  double constraint;
  LDLoci *pLDLoci = NULL;
  double traitPos;		/* trait position for multipoint analysis */
  TraitLocus *pTraitLocus;
  int traitLocus;
  int leftMarker = -1;
  int posIdx;
  int k;
  double avgLR;
  double ppl;
  double ldppl, ppld;
  int markerSetChanged;		/* flag for multipoint analysis */
  int locusListChanged;		/* flag for multipoint analysis */
  int prevFirstMarker;		/* first marker in the set for multipoint analysis */
  int prevLastMarker;		/* last marker in the set for multipoint analysis */
  int prevTraitInd;
  double *prevPos, *currPos;	/* for MP */
  int locus;
  int thresholdIdx = -1;
  double threshold = 0;
  int R_square_flag = FALSE;
  double R_square = 0;
  int dprime0Idx = 0;
  int mkrFreqIdx;
  double mkrFreq;
  int maxThetaIdx = 0;
  int maxDPrimeIdx = 0;
  int maxDPrimeIdx_at_theta0 = 0;
  double max_at_theta0;
  double lr;
  int theta0Idx = 0;
  double max_at_dprime0;
  int maxTheta_at_dprime0 = -1;
  double dprime;
  int totalLoci;
  int status;
  double *marker1Pos, *marker2Pos;
  double relativePos;
  int traitIndex = 0;
  double dist;
  int initialFlag = 0;

  SubLocusList savedLocusList;
  SubLocusList nullLocusList;

  clock_t time0, time1, time2;
  int numberOfCompute = 0;

#ifndef NO_POLYNOMIAL
  Polynomial *initialProbPoly[3];
#endif
  double initialProb[3];
  void *initialProbAddr[3];

  memset (&savedLocusList, 0, sizeof (savedLocusList));
  memset (&nullLocusList, 0, sizeof (nullLocusList));

  memset( &modelType, 0, sizeof(modelType));

#ifdef DEBUG
  //  mtrace();
#endif

  /* Initialize the logging system. */
  logInit ();
  /* logSet(LOGINPUTFILE, LOGDEBUG); */
  //logSet(LOGGENOELIM, LOGDEBUG);
  /* logSet(LOGPEELGRAPH, LOGFATAL); */
  //logSet (LOGLIKELIHOOD, LOGWARNING);
  //logSet(LOGLIKELIHOOD, LOGDEBUG); 
  //logSet(LOGPARENTALPAIR, LOGDEBUG); 
  //logSet(LOGSETRECODING, LOGDEBUG);

  /* Start by parsing command line arguments. Most essential: figure
   * out where the configuration file lives. */
  for (i = 1; i < argc; i++)
    {
      if (argv[i][0] == '-')
	switch (argv[i][1])
	  {
	  case '?':
	    /* Help */
	    fprintf (stdout, "Usage:\n");
	    fprintf (stdout, "  %s [-?][-s][-c <file>]\nwhere:\n", argv[0]);
	    fprintf (stdout, "      -? : this output;\n");
	    fprintf (stdout, "      -s : run serially;\n");
	    fprintf (stdout,
		     "      -c <file> : restart calculation from specified file.\n");
	    fprintf (stdout,
		     "Checkpoint data (for restarting) appears on stderr.\n");
	    exit (EXIT_FAILURE);
	    break;
	  case 's':
	    /* Run serially. */
	    serial = TRUE;
	    break;
	  case 'c':
	    /* Restart from checkpoint file. */
	    strncpy (ckptfile, argv[i + 1], KMAXFILENAMELEN);
	    break;
	  }
      else if (strlen (configfile) != 0)
	{
	  /* Unexpected argument; we already have a configuration file! Punt. */
	  KLOG (LOGDEFAULT, LOGFATAL,
		"Unexpected command line argument '%s'; aborting.\n",
		argv[i]);
	}
      else if (strlen (argv[i]) >= KMAXFILENAMELEN)
	{
	  /* Configuration file name too long! Punt. */
	  KLOG (LOGDEFAULT, LOGFATAL,
		"Configuration file name '%s' exceeds limit of %d; aborting.\n",
		argv[i], KMAXFILENAMELEN);
	}
      else
	{
	  /* Got a configuration file name. Copy it. */
	  strncpy (configfile, argv[i], KMAXFILENAMELEN);
	}
      i++;
    }

  /* Check to see if the configuration file name was specified. */
  KASSERT ((strlen (configfile) > 0),
	   "No configuration file specified; aborting.\n");

  /* set the default unknown person ID */
  modelOptions.sUnknownPersonID = malloc (sizeof (char) * 2);
  strcpy (modelOptions.sUnknownPersonID, "0");

  /* set default values for PPL calculations */
  /* LRs are weighted heavier for theta less than the cutoff */
  modelOptions.thetaCutoff[0] = 0.05;
  modelOptions.thetaCutoff[1] = 0.05;
  /* weight ofr theta less than the cutoff */
  modelOptions.thetaWeight = 0.95;
  /* prior probability of linkage */
  modelOptions.prior = 0.02;
  /* prior probability of LD given close linkage */
  modelOptions.LDprior = 0.02;

  /* Parse the configuration file. */
  KASSERT (readConfigFile (configfile, &modelType, &modelRange, &modelOptions)
	   != ERROR, "Error in configuration file; aborting.\n");

  /* For now, reject all models we can't deal with. So use KASSERT to
   * check that we're looking at, e.g., twopoint dichotomous models
   * with direct eval (not polynomial --> add this to model?) and
   * give appropriate error message otherwise. */
  KASSERT (modelRange.nalleles == 2, "Only biallelic traits supported.\n");

  /* the difference between QT and CT is whether we use threshold or not. Under CT -  yes to
   * threshold, under QT - no threshold */
  if (modelRange.ntthresh > 0 && modelType.trait != DT)
    {
      modelType.trait = CT;
    }
  if (modelType.trait == QT)
    {
      /* threshold value will not be used in any meaningful way, but we will use it for 
         the loop */
      modelRange.ntthresh = 1;
      if (modelRange.tthresh == NULL)
	{
	  modelRange.tthresh = (double **) malloc (sizeof (double *));
	  for (i = 0; i < modelRange.nlclass; i++)
	    {
	      modelRange.tthresh[i] = malloc (sizeof (double));
	    }
	}
    }

  /* open output files */
  fpHet = fopen (avghetfile, "w");
  KASSERT (fpHet != NULL, "Error in opening file %s for write.\n",
	   avghetfile);
  /*  fpHomo = fopen (avghomofile, "w");
     KASSERT (fpHomo != NULL, "Error in opening file %s for write.\n",
     avghomofile);
   */
  if (modelType.type == TP)
    {
      fpPPL = fopen (pplfile, "w");
      KASSERT (fpPPL != NULL, "Error in opening file %s for write.\n",
	       pplfile);
      fprintf (fpPPL, "%4s %15s %9s %6s ", "CHR", "MARKER", "cM", "PPL");
      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
	{
	  fprintf (fpPPL, "%6s %6s ", "LD-PPL", "PPLD");
	}
      fprintf (fpPPL, "\n");
      fflush (fpPPL);
    }

#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      polynomialInitialization ();
      fprintf (stderr,
	       "!!!!!!!!!!!The Computation is done in polynomial mode!!!!!!!!!!!!!!!\n");
    }
  else
    {
      fprintf (stderr, "Polynomial is off!\n");
    }
#else
  fprintf (stderr, "No polynomial available.\n");
#endif

  /* Read in the map file. */
  read_mapfile (mapfile);

  /* Initialize the locus list and read in the marker file. */
  memset (&originalLocusList, 0, sizeof (originalLocusList));
  /* read in what loci are in the pedigree file */
  read_datafile (datafile);

  /* The configuration has all the information about the disease trait if any */
  if (originalLocusList.numTraitLocus > 0)
    {
      /* we are not doing marker to marker analysis
       * Need to add the alleles into trait locus 
       * Assume the traitLoucs is 0 for now  - Need to fix this later */
      traitLocus = 0;
      pLocus = originalLocusList.ppLocusList[traitLocus];
      pTraitLocus = pLocus->pTraitLocus;
      add_allele (pLocus, "D", 0.5);
      add_allele (pLocus, "d", 0.5);
      /* fix number of trait variables at 1 for now */
      pTraitLocus->numTrait = 1;
      pTrait = add_trait (0, pTraitLocus, modelType.trait);
      pTrait->numLiabilityClass = modelRange.nlclass;
      if (modelType.trait == QT || modelType.trait == CT)
	{
	  modelType.min = (modelType.minOriginal - modelType.mean)/modelType.sd;
	  modelType.max = (modelType.maxOriginal - modelType.mean)/modelType.sd;
	  pTrait->minFlag = modelType.minFlag;
	  pTrait->maxFlag = modelType.maxFlag;
	  pTrait->min = modelType.min;
	  pTrait->max = modelType.max;
	  pTrait->functionQT = modelType.distrib;
	  if (modelType.distrib == QT_FUNCTION_T)
	    pTrait->dfQT = modelType.constants[0];
	  pTrait->sampleMean = modelType.mean;
	  pTrait->sampleSD = modelType.sd;
	  pTrait->unknownTraitValue =
	    modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN];
	  pTrait->lessCutoffFlag =
	    modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED];
	  pTrait->moreCutoffFlag =
	    modelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED];
	}
    }

  /* read in marker allele frequencies */
  read_markerfile (markerfile);

  /* build allele set information */
  for (locus = 0; locus < originalLocusList.numLocus; locus++)
    {
      construct_original_allele_set_list (locus);
    }

  /* Initialize the pedigree set datastructure and read in the pedigrees. */
  memset (&pedigreeSet, 0, sizeof (PedigreeSet));
  read_pedfile (pedfile, &pedigreeSet);

  fflush (stderr);
  fflush (stdout);

  /* allocate space for results */
  if (modelType.type == TP)
    {
      modelType.numMarkers = 1;
      totalLoci = 2;
      /* two point analysis */
      if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
	{
	  /* in order to simplify looping, even for LE, we add a fake LD parameter dprime=0, which
	   * is LE */
	  modelRange.ndprime = 1;
	  modelRange.dprime = (double *) calloc (1, sizeof (double));
	  modelRange.dprime[0] = 0;
	  pLambdaCell = findLambdas (&modelRange, 2, 2);
	}
    }
  else
    {
      /* we are doing multipoint analysis */
      totalLoci = modelType.numMarkers + originalLocusList.numTraitLocus;
      if (modelRange.tlmark == TRUE)
	{
	  /* add marker positions to the list of positions we want to conduct analysis */
	  for (i = 0; i < originalLocusList.numLocus; i++)
	    {
	      pLocus = originalLocusList.ppLocusList[i];
	      if (pLocus->locusType == LOCUS_TYPE_TRAIT)
		continue;
	      addTraitLocus (&modelRange,
			     pLocus->pMapUnit->mapPos[SEX_AVERAGED]);
	    }
	}
    }

  /* initialize some work space */
  initialize_parental_pair_workspace (&parentalPairSpace,
				      originalLocusList.numLocus);
  /* allocate transmission probability matrix */
  build_xmission_matrix (&nullMatrix, totalLoci);
  build_xmission_matrix (&altMatrix, totalLoci);
  xmissionMatrix = nullMatrix;

  /* initialize loci by doing genotype elimination, set recoding */
  initialize_loci (&pedigreeSet);

#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      for (k = 0; k < 3; k++)
	{
	  initialProbPoly[k] = constantExp (1.0);
	  initialProbAddr[k] = initialProbPoly[k];
	}
    }
  else
    {
      for (k = 0; k < 3; k++)
	{
	  initialProb[k] = 1.0;
	  initialProbAddr[k] = &initialProb[k];
	}
    }
#else
  for (k = 0; k < 3; k++)
    {
      initialProb[k] = 1.0;
      initialProbAddr[k] = &initialProb[k];
    }
#endif

  /* only for multipoint - we don't handle LD under multipoint yet */
  if (modelType.type == MP)
    {
      /* allocate space to save temporary results */
      if (modelType.trait == DT)
	{
	  /* first dimension is pedigree */
	  likelihoodDT =
	    (double ***) calloc (sizeof (double **),
				 pedigreeSet.numPedigree + 1);
	  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree + 1; pedIdx++)
	    {
	      /* second dimension is gene freq */
	      likelihoodDT[pedIdx] =
		(double **) calloc (sizeof (double *), modelRange.ngfreq);
	      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++)
		{
		  /* third dimension is penetrance */
		  likelihoodDT[pedIdx][gfreqInd] =
		    (double *) calloc (sizeof (double), modelRange.npenet);
		}
	    }
	}
      else
	{			/* QT */
	  /* first dimension is pedigree */
	  likelihoodQT =
	    (double *****) calloc (sizeof (double ****),
				   pedigreeSet.numPedigree + 1);
	  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree + 1; pedIdx++)
	    {
	      /* second dimension is gene freq */
	      likelihoodQT[pedIdx] =
		(double ****) calloc (sizeof (double ***), modelRange.ngfreq);
	      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++)
		{

		  /* third dimension is mean */
		  likelihoodQT[pedIdx][gfreqInd] =
		    (double ***) calloc (sizeof (double **),
					 modelRange.npenet);
		  for (penIdx = 0; penIdx < modelRange.npenet; penIdx++)
		    {
		      /* fourth dimension is SD */
		      likelihoodQT[pedIdx][gfreqInd][penIdx] =
			(double **) calloc (sizeof (double *),
					    modelRange.nparam);
		      for (paramIdx = 0; paramIdx < modelRange.nparam;
			   paramIdx++)
			{
			  /* 5th dimension is threshold */
			  likelihoodQT[pedIdx][gfreqInd][penIdx][paramIdx] =
			    (double *) calloc (sizeof (double),
					       modelRange.ntthresh);

			}

		    }
		}
	    }
	}
    }

  /* find out the max we need to allocate */
  /* after genotype lists have been built, we want to pre-allocate parental pair work space
   * it used to be done dynamically, but it's too costly 
   * parental pair is based on each nuclear family */
  stat_parental_pair_workspace (&pedigreeSet);

  /* after genotype elimination and tracking the max work space needed 
   * for constructing parental pair */
  allocate_parental_pair_workspace (&parentalPairSpace,
				    modelType.numMarkers + 1);

  /* conditional likelihood storage space for each individual */
  allocate_likelihood_space (&pedigreeSet, modelType.numMarkers + 1);

  /* The following segments are for RADSMM - store lods in binary format to a file */
#if FALSE
  /* Set up storage before you try to write LOD scores. */
  KASSERT ((RADSMM_setup_init (&header, 255) == 0),
	   "RADSMM initialization error.\n");

  /* Set up the analysis type. */
  KASSERT ((RADSMM_setup_type (&header,
			       ((modelType.type == TP) ? '2' : 'M'),
			       ((modelType.trait ==
				 DT) ? 'D' : ((modelType.trait ==
					       QT) ? 'Q' : 'C')),
			       ((modelOptions.equilibrium ==
				 LE) ? 'N' : 'Y')) == 0),
	   "RADSMM type initialization error.\n");

  /* Set up ordering of array dimensions (TODO: check this out). */
  KASSERT ((RADSMM_setup_ordering (&header, 'B') == 0),
	   "RADSMM ordering initialization error.\n");

  /* Set up the pedigrees. */
  KASSERT ((RADSMM_setup_pedigree
	    (&header, NULL, (long) pedigreeSet.numPedigree) == 0),
	   "RADSMM pedigree initialization error.\n");

  /* Set up the theta structures. */
  KASSERT ((RADSMM_setup_theta
	    (&header, modelRange.theta, (long) modelRange.ntheta,
	     ((modelRange.ngender == 1) ? 'D' : 'G'))),
	   "RADSMM theta initialization error.\n");

  /* Set up the liability classes. */
  KASSERT ((RADSMM_setup_LC (&header, modelRange.nlclass) == 0),
	   "RADSMM liability class initialization error.\n");

  /* Set up the penetrance arrays. TODO: will need work for
   * multiallelic diseases. Here, we just assume we're dealing with
   * DD, Dd, and dd. */
  for (i = 0; i < modelRange.nlclass; i++)
    KASSERT ((RADSMM_setup_penetrance (&header, i, modelRange.penet[i][0],
				       modelRange.penet[i][1],
				       modelRange.penet[i][2],
				       (long) modelRange.npenet) == 0),
	     "RADSMM penetrance initialization error.\n");

  /* Gene frequencies are next. */
  KASSERT ((RADSMM_setup_geneFreq (&header, modelRange.gfreq,
				   (long) modelRange.ngfreq) == 0),
	   "RADSMM gene frequency initialization error.\n");
#endif



  /**********/
//  exit(ERROR);
  /**********/


  time0 = clock ();
  time1 = clock ();


  if (modelType.trait == DT)
    fprintf (stderr, "Dichotomous Trait & ");
  else if (modelType.trait == QT)
    fprintf (stderr, "Quantitative Trait without threshoold & ");
  else
    fprintf (stderr, "Quantitative Trait with threshoold & ");
  fprintf (stderr, "%s\n",
	   (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) ? "LE" : "LD");

  fprintf (stderr, "Total number of markers in data: %d\n",
	   originalLocusList.numLocus - originalLocusList.numTraitLocus);
  fprintf (stderr, "Total number of trait locus in data: %d\n",
	   originalLocusList.numTraitLocus);
  fprintf (stderr, "Total number of families in pedigree file: %d\n",
	   pedigreeSet.numPedigree);
  fprintf (stderr, "Number of loci to use for analysis: %d\n",
	   modelType.numMarkers + originalLocusList.numTraitLocus);


  /* Initialize the connection to the infrastructure. */
  /* niceInit (); */

  /* assume the trait locus is the first one in the list */
  traitLocus = 0;
  pLocus = originalLocusList.ppLocusList[traitLocus];
  pTraitLocus = originalLocusList.ppLocusList[traitLocus]->pTraitLocus;
  pTrait = pTraitLocus->pTraits[traitLocus];
  if (modelType.type == TP)
    {
      /* Two point. */
      if (originalLocusList.pLDLoci == NULL)
	{
	  originalLocusList.pLDLoci = (LDLoci *) malloc (sizeof (LDLoci));
	  memset (originalLocusList.pLDLoci, 0, sizeof (LDLoci));
	}
      pLDLoci = &originalLocusList.pLDLoci[0];
      originalLocusList.numLDLoci = 1;

      if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
	{
	  /* fake some LD information to simplify looping */
	  pLDLoci->numAllele1 = 2;
	  pLDLoci->ppDPrime = (double **) malloc (sizeof (double *));
	  pLDLoci->ppDPrime[0] = (double *) malloc (sizeof (double));
	  pLDLoci->ppDValue = (double **) malloc (sizeof (double *));
	  pLDLoci->ppDValue[0] = (double *) malloc (sizeof (double));
	  pLDLoci->ppHaploFreq = (double **) malloc (sizeof (double *) * 2);
	  pLDLoci->ppHaploFreq[0] = (double *) malloc (sizeof (double) * 2);
	  pLDLoci->ppHaploFreq[1] = (double *) malloc (sizeof (double) * 2);

	  /* initialize it */
	  pLDLoci->ppDPrime[0][0] = 0;
	}

      locusList = &savedLocusList;
      savedLocusList.numLocus = 2;
      savedLocusList.pLocusIndex = (int *) malloc (sizeof (int) *
						   savedLocusList.numLocus);
      for (i = 0; i < 3; i++)
	{
	  savedLocusList.pPrevLocusDistance[i] =
	    (double *) malloc (sizeof (double) * savedLocusList.numLocus);
	  savedLocusList.pNextLocusDistance[i] =
	    (double *) malloc (sizeof (double) * savedLocusList.numLocus);

	  savedLocusList.pPrevLocusDistance[i][0] = -1;
	  savedLocusList.pNextLocusDistance[i][1] = -1;
	}

#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	{
	  /* populate the matrix */
	  status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
					     0,	/* cell index */
					     -1,	/* last het locus */
					     -1,	/* last het pattern (P-1 or M-2) */
					     0,	/* number of het loci */
					     0);	/* current locus - start with 0 */
	  makePolynomialStamp ();
	}
#endif


      for (loc1 = 0; loc1 < originalLocusList.numLocus - 1; loc1++)
	{
	  savedLocusList.pLocusIndex[0] = loc1;
	  pLocus1 = originalLocusList.ppLocusList[loc1];
	  if (modelOptions.markerAnalysis != FALSE
	      && pLocus1->locusType != LOCUS_TYPE_MARKER)
	    continue;

	  for (loc2 = loc1 + 1; loc2 < originalLocusList.numLocus; loc2++)
	    {
	      pLocus2 = originalLocusList.ppLocusList[loc2];
	      if (pLocus2->locusType != LOCUS_TYPE_MARKER)
		continue;
	      savedLocusList.pLocusIndex[1] = loc2;

	      /* find out number of alleles this marker locus has */
	      if (modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM)
		{
		  /* get the LD parameters */
		  pLambdaCell =
		    findLambdas (&modelRange, pLocus1->numOriginalAllele,
				 pLocus2->numOriginalAllele);
		  reallocate_LD_loci (pLDLoci, pLocus1->numOriginalAllele,
				      pLocus2->numOriginalAllele);
		  pLDLoci->locus1 = loc1;
		  pLDLoci->locus2 = loc2;
		  pLDLoci->numAllele1 = pLocus1->numOriginalAllele;
		  pLDLoci->numAllele2 = pLocus2->numOriginalAllele;
		  if (pLocus1->numOriginalAllele == 2
		      && pLocus2->numOriginalAllele == 2)
		    R_square_flag = TRUE;
		  else
		    R_square_flag = FALSE;
		}

	      loopMarkerFreqFlag = 0;
	      if (modelRange.nafreq >= 2
		  && modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM
		  && pLocus2->numOriginalAllele == 2)
		{
		  loopMarkerFreqFlag = 1;
		}
	      else if (modelRange.nafreq == 0)
		{
		  /* add a fake one to facilitate loops and other handlings */
		  addAlleleFreq (&modelRange, pLocus2->pAlleleFrequency[0]);
		}
	      else
		{
		  modelRange.nafreq = 1;
		  modelRange.afreq[0] = pLocus2->pAlleleFrequency[0];
		}

	      /* allocate result storage */
	      allocate_tp_result_storage ();

	      /* we will force marker allele frequency loop to execute at least once */
	      for (mkrFreqIdx = 0;
		   mkrFreqIdx == 0 || mkrFreqIdx < modelRange.nafreq;
		   mkrFreqIdx++)
		{
		  mkrFreq = pLocus2->pAlleleFrequency[0];
		  /* we should only loop over marker allele frequency under twopoint
		   * and when markers are SNPs (only have two alleles) */
		  if (loopMarkerFreqFlag)
		    {
		      mkrFreq = modelRange.afreq[mkrFreqIdx];
		      /* update the locus */
		      pLocus2->pAlleleFrequency[0] = mkrFreq;
		      pLocus2->pAlleleFrequency[1] = 1 - mkrFreq;
#ifndef NO_POLYNOMIAL
		      if (modelOptions.polynomial == TRUE)
			;
		      else
			update_locus (&pedigreeSet, loc2);
#else
		      update_locus (&pedigreeSet, loc2);
#endif
		    }
		  /* Loop over the penetrances, genefrequencies, thetas and call
		     the likelihood calculation, storing each value obtained to
		     disk. */
		  if (modelType.trait == DICHOTOMOUS)
		    {

		      for (penIdx = 0; penIdx < modelRange.npenet; penIdx++)
			{
			  if (modelOptions.markerAnalysis == FALSE
			      && pLocus1->locusType == LOCUS_TYPE_TRAIT)
			    {
			      for (liabIdx = 0; liabIdx < modelRange.nlclass;
				   liabIdx++)
				{
				  pen_DD =
				    modelRange.penet[liabIdx][0][penIdx];
				  pen_Dd =
				    modelRange.penet[liabIdx][1][penIdx];
				  pen_dd =
				    modelRange.penet[liabIdx][2][penIdx];
				  pTrait->penetrance[2][liabIdx][0][0] =
				    pen_DD;
				  pTrait->penetrance[2][liabIdx][0][1] =
				    pen_Dd;
				  pTrait->penetrance[2][liabIdx][1][0] =
				    pen_Dd;
				  pTrait->penetrance[2][liabIdx][1][1] =
				    pen_dd;
				  pTrait->penetrance[1][liabIdx][0][0] =
				    1 - pen_DD;
				  pTrait->penetrance[1][liabIdx][0][1] =
				    1 - pen_Dd;
				  pTrait->penetrance[1][liabIdx][1][0] =
				    1 - pen_Dd;
				  pTrait->penetrance[1][liabIdx][1][1] =
				    1 - pen_dd;
				}


#ifndef NO_POLYNOMIAL
			      if (modelOptions.polynomial == TRUE)
				;
			      else
				update_penetrance (&pedigreeSet, traitLocus);
#else
			      update_penetrance (&pedigreeSet, traitLocus);
#endif
			    }
			  for (gfreqInd = 0; gfreqInd < modelRange.ngfreq;
			       gfreqInd++)
			    {
			      gfreq = modelRange.gfreq[gfreqInd];
			      if (modelOptions.markerAnalysis == FALSE)
				{
				  pLocus->pAlleleFrequency[0] = gfreq;
				  pLocus->pAlleleFrequency[1] = 1 - gfreq;

#ifndef NO_POLYNOMIAL
				  if (modelOptions.polynomial == TRUE)
				    ;
				  else
				    update_locus (&pedigreeSet, loc1);
#else
				  update_locus (&pedigreeSet, loc1);
#endif
				}
			      /* get the likelihood at 0.5 first and LD=0 */
			      if (modelOptions.equilibrium !=
				  LINKAGE_EQUILIBRIUM)
				{
				  set_null_dprime (pLDLoci);
				  setup_LD_haplotype_freq (pLDLoci);
				}
			      for (k = 0; k < 3; k++)
				{
				  locusList->pNextLocusDistance[k][0] = 0.5;
				  locusList->pPrevLocusDistance[k][1] = 0.5;
				}

#ifndef NO_POLYNOMIAL
			      if (modelOptions.polynomial == TRUE)
				;
			      else
				/* populate the matrix */
				status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
								   0,	/* cell index */
								   -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
								   0,	/* number of het loci */
								   0);	/* current locus - start with 0 */

#else
			      /* populate the matrix */
			      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
								 0,	/* cell index */
								 -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
								 0,	/* number of het loci */
								 0);	/* current locus - start with 0 */
#endif

			      compute_likelihood (&pedigreeSet);


			      if (numberOfCompute == 0)
				{
				  time1 = clock ();
				}
			      numberOfCompute++;

			      if (pedigreeSet.likelihood == 0.0 &&
				  pedigreeSet.log10Likelihood == -9999.99)
				{
				  fprintf (stderr,
					   "Theta 0.5 has likelihood 0\n");
				  fprintf (stderr, "dgf=%f\n", gfreq);
				  for (liabIdx = 0;
				       liabIdx < modelRange.nlclass;
				       liabIdx++)
				    {
				      pen_DD =
					modelRange.penet[liabIdx][0][penIdx];
				      pen_Dd =
					modelRange.penet[liabIdx][1][penIdx];
				      pen_dd =
					modelRange.penet[liabIdx][2][penIdx];
				      fprintf (stderr,
					       "Liab %d penentrance %f %f %f\n",
					       liabIdx + 1, pen_DD, pen_Dd,
					       pen_dd);
				    }

				  exit (-1);
				}
			      /* save the results for NULL */
			      for (pedIdx = 0;
				   pedIdx < pedigreeSet.numPedigree; pedIdx++)
				{
				  /* save the likelihood at null */
				  pPedigree =
				    pedigreeSet.ppPedigreeSet[pedIdx];
				  pedigreeSet.nullLikelihood[pedIdx] =
				    pPedigree->likelihood;
				}

			      log10_likelihood_null =
				pedigreeSet.log10Likelihood;
			      for (dprimeIdx = 0;
				   dprimeIdx < pLambdaCell->ndprime;
				   dprimeIdx++)
				{
				  if (modelOptions.equilibrium !=
				      LINKAGE_EQUILIBRIUM)
				    {
				      copy_dprime (pLDLoci,
						   pLambdaCell->
						   lambda[dprimeIdx]);
				      setup_LD_haplotype_freq (pLDLoci);
				      /* calculate R square if the marker is a SNP */
				      if (R_square_flag == TRUE)
					R_square =
					  calculate_R_square (pLocus1->
							      pAlleleFrequency
							      [0],
							      pLocus2->
							      pAlleleFrequency
							      [0],
							      pLDLoci->
							      ppDValue[0][0]);
				      else
					R_square = -1;

				      if (-ERROR_MARGIN <=
					  pLambdaCell->lambda[dprimeIdx][0][0]
					  && pLambdaCell->
					  lambda[dprimeIdx][0][0] <=
					  ERROR_MARGIN)
					dprime0Idx = dprimeIdx;
				    }

				  for (thetaInd = 0;
				       thetaInd < modelRange.ntheta;
				       thetaInd++)
				    {
				      if (modelOptions.mapFlag == SA)
					{
					  theta[0] =
					    modelRange.theta[0][thetaInd];
					  theta[1] =
					    modelRange.theta[1][thetaInd];
					  for (k = 0; k < 3; k++)
					    {
					      locusList->
						pNextLocusDistance[k][0] =
						theta[0];
					      locusList->
						pPrevLocusDistance[k][1] =
						theta[0];
					    }
					}
				      else
					{
					  locusList->
					    pNextLocusDistance[MAP_MALE][0] =
					    locusList->
					    pPrevLocusDistance[MAP_MALE][1] =
					    modelRange.theta[0][thetaInd];
					  locusList->
					    pNextLocusDistance[MAP_FEMALE][0]
					    =
					    locusList->
					    pPrevLocusDistance[MAP_FEMALE][1]
					    = modelRange.theta[1][thetaInd];
					}

#ifndef NO_POLYNOMIAL
				      if (modelOptions.polynomial == TRUE)
					;
				      else
					/* populate the matrix */
					status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
									   0,	/* cell index */
									   -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
									   0,	/* number of het loci */
									   0);	/* current locus - start with 0 */
#else
				      /* populate the matrix */
				      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
									 0,	/* cell index */
									 -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
									 0,	/* number of het loci */
									 0);	/* current locus - start with 0 */
#endif

				      compute_likelihood (&pedigreeSet);

				      log10_likelihood_alternative =
					pedigreeSet.log10Likelihood;
				      if (pedigreeSet.likelihood == 0.0
					  && pedigreeSet.log10Likelihood ==
					  -9999.99)
					{
					  log10_likelihood_ratio = 0;
					}
				      else
					{
					  log10_likelihood_ratio =
					    log10_likelihood_alternative -
					    log10_likelihood_null;
					}
				      /* check for overflow problem !!! */
				      if (log10_likelihood_ratio >=
					  __DBL_MAX_10_EXP__ - 1)
					{
					  likelihood_ratio = __DBL_MAX__;
					  tp_result[dprimeIdx][thetaInd]
					    [mkrFreqIdx].lr_total =
					    __DBL_MAX__;
					}
				      else
					/* check for underflow problem too !!! */
				      if (log10_likelihood_ratio <=
					    __DBL_MIN_10_EXP__ + 1)
					{
					  likelihood_ratio = 0;
					}
				      else
					{
					  likelihood_ratio =
					    pow (10.0,
						 log10_likelihood_ratio);
					  tp_result[dprimeIdx][thetaInd]
					    [mkrFreqIdx].lr_total +=
					    likelihood_ratio;
					}
				      tp_result[dprimeIdx][thetaInd]
					[mkrFreqIdx].lr_count++;
				      //fprintf(stderr, "likelihood ratio: %e.\n", likelihood_ratio);

				      /* caculating the HET */
				      for (j = 0; j < modelRange.nalpha; j++)
					{
					  alphaV = modelRange.alpha[j];
					  alphaV2 = 1 - alphaV;
					  if (alphaV2 < 0)
					    alphaV2 = 0;
					  log10HetLR = 0;
					  for (pedIdx = 0;
					       pedIdx <
					       pedigreeSet.numPedigree;
					       pedIdx++)
					    {
					      pPedigree =
						pedigreeSet.
						ppPedigreeSet[pedIdx];
					      homoLR =
						pPedigree->likelihood /
						pedigreeSet.
						nullLikelihood[pedIdx];
					      log10HetLR +=
						log10 (alphaV * homoLR +
						       (1 - alphaV));
					    }
					  if (log10HetLR >=
					      __DBL_MAX_10_EXP__ - 1)
					    {
					      hetLR = __DBL_MAX__;
					      tp_result[dprimeIdx][thetaInd]
						[mkrFreqIdx].het_lr_total =
						__DBL_MAX__;
					    }
					  else
					    if (log10HetLR <=
						__DBL_MIN_10_EXP__ + 1)
					    {
					      hetLR = 0;
					    }
					  else
					    {
					      hetLR = pow (10, log10HetLR);
					      tp_result[dprimeIdx][thetaInd]
						[mkrFreqIdx].het_lr_total +=
						hetLR;
					    }
					  if (tp_result[dprimeIdx][thetaInd]
					      [mkrFreqIdx].max_penIdx < 0
					      || hetLR >
					      tp_result[dprimeIdx][thetaInd]
					      [mkrFreqIdx].max_lr)
					    {
					      tp_result[dprimeIdx][thetaInd]
						[mkrFreqIdx].max_lr = hetLR;
					      tp_result[dprimeIdx][thetaInd]
						[mkrFreqIdx].max_alpha =
						alphaV;
					      tp_result[dprimeIdx][thetaInd]
						[mkrFreqIdx].max_gfreq =
						gfreq;
					      tp_result[dprimeIdx][thetaInd]
						[mkrFreqIdx].max_penIdx =
						penIdx;
					      tp_result[dprimeIdx][thetaInd]
						[mkrFreqIdx].R_square =
						R_square;
					      tp_result[dprimeIdx][thetaInd]
						[mkrFreqIdx].max_mf = mkrFreq;
					    }
					}	/* end of calculating HET LR */
				    }	/* end of theta loop */
				}	/* end of D prime loop */
			      if (modelOptions.markerAnalysis != FALSE)
				{
				  /* marker to marker analysis, marker allele frequency is fixed */
				  gfreqInd = modelRange.ngfreq;
				  break;
				}
			    }	/* end of genFreq loop */
			  if (modelOptions.markerAnalysis != FALSE)
			    {
			      /* marker to marker analysis, penetrance stays at 1 */
			      break;
			    }
			}	/* end of penetrance loop */
		    }		/* end of TP */
		  else
		    /* should be QT or COMBINED - twopoint */
		    {
		      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq;
			   gfreqInd++)
			{
			  gfreq = modelRange.gfreq[gfreqInd];
			  if (modelOptions.markerAnalysis == FALSE)
			    {
			      pLocus->pAlleleFrequency[0] = gfreq;
			      pLocus->pAlleleFrequency[1] = 1 - gfreq;

			      update_locus (&pedigreeSet, traitLocus);
			    }
			  /* this should be MEAN + SD */
			  for (paramIdx = 0;
			       (paramIdx == 0
				&& modelType.distrib ==
				QT_FUNCTION_CHI_SQUARE)
			       || (modelType.distrib != QT_FUNCTION_CHI_SQUARE
				   && paramIdx < modelRange.nparam);
			       paramIdx++)
			    {
			      for (penIdx = 0; penIdx < modelRange.npenet;
				   penIdx++)
				{
				  breakFlag = FALSE;
				  for (thresholdIdx = 0;
				       thresholdIdx < modelRange.ntthresh;
				       thresholdIdx++)
				    {
				      if (modelOptions.markerAnalysis ==
					  FALSE)
					{
					  for (liabIdx = 0;
					       liabIdx < modelRange.nlclass;
					       liabIdx++)
					    {
					      mean_DD =
						modelRange.
						penet[liabIdx][0][penIdx];
					      mean_Dd =
						modelRange.
						penet[liabIdx][1][penIdx];
					      mean_dd =
						modelRange.
						penet[liabIdx][2][penIdx];
					      SD_DD =
						modelRange.
						param[liabIdx][0][0]
						[paramIdx];
					      SD_Dd =
						modelRange.
						param[liabIdx][1][0]
						[paramIdx];
					      SD_dd =
						modelRange.
						param[liabIdx][2][0]
						[paramIdx];
					      /* threshold for QT */
					      threshold =
						modelRange.
						tthresh[liabIdx]
						[thresholdIdx];


					      /* check against the hard coded constraint */
					      if (modelType.distrib !=
						  QT_FUNCTION_CHI_SQUARE)
						{
						  constraint =
						    pow (1.0 - gfreq,
							 2) * mean_dd *
						    SD_dd + 2 * gfreq * (1 -
									 gfreq)
						    * mean_Dd * SD_Dd +
						    pow (gfreq,
							 2) * mean_DD * SD_DD;
						  /*      fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
						     constraint, gfreq, mean_DD, SD_DD, 
						     mean_Dd, SD_DD, 
						     mean_dd, SD_dd);
						   */
						  if (constraint >= 3.0
						      || constraint <= -3.0)
						    {
						      breakFlag = TRUE;
						      break;
						    }
						}
					      pTrait->means[liabIdx][0][0] =
						mean_DD;
					      pTrait->means[liabIdx][0][1] =
						mean_Dd;
					      pTrait->means[liabIdx][1][0] =
						mean_Dd;
					      pTrait->means[liabIdx][1][1] =
						mean_dd;
					      pTrait->stddev[liabIdx][0][0] =
						SD_DD;
					      pTrait->stddev[liabIdx][0][1] =
						SD_Dd;
					      pTrait->stddev[liabIdx][1][0] =
						SD_Dd;
					      pTrait->stddev[liabIdx][1][1] =
						SD_dd;

					      /* threshold for QT */
					      pTrait->cutoffValue[liabIdx] =
						threshold;

					    }	/* liability class Index */
					  if (breakFlag == TRUE)
					    continue;
					  update_penetrance (&pedigreeSet,
							     traitLocus);
					}	/* marker to marker analysis */
				      /* get the likelihood at 0.5 first and LD=0 */
				      if (modelOptions.equilibrium !=
					  LINKAGE_EQUILIBRIUM)
					{
					  set_null_dprime (pLDLoci);
					  setup_LD_haplotype_freq (pLDLoci);
					}
				      for (k = 0; k < 3; k++)
					{
					  locusList->
					    pNextLocusDistance[k][0] = 0.5;
					  locusList->
					    pPrevLocusDistance[k][1] = 0.5;
					}

#ifndef NO_POLYNOMIAL
				      if (modelOptions.polynomial == TRUE)
					;
				      else
					/* populate the matrix */
					status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
									   0,	/* cell index */
									   -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
									   0,	/* number of het loci */
									   0);	/* current locus - start with 0 */
#else
				      /* populate the matrix */
				      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
									 0,	/* cell index */
									 -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
									 0,	/* number of het loci */
									 0);	/* current locus - start with 0 */
#endif

				      compute_likelihood (&pedigreeSet);


				      if (numberOfCompute == 0)
					{
					  time1 = clock ();
					}
				      numberOfCompute++;

				      if (pedigreeSet.likelihood == 0.0 &&
					  pedigreeSet.log10Likelihood ==
					  -9999.99)
					{
					  fprintf (stderr,
						   "Theta 0.5 has likelihood 0\n");
					  fprintf (stderr, "dgf=%f\n", gfreq);
					  for (liabIdx = 0;
					       liabIdx < modelRange.nlclass;
					       liabIdx++)
					    {
					      pen_DD =
						modelRange.
						penet[liabIdx][0][penIdx];
					      pen_Dd =
						modelRange.
						penet[liabIdx][1][penIdx];
					      pen_dd =
						modelRange.
						penet[liabIdx][2][penIdx];
					      fprintf (stderr,
						       "Liab %d penentrance %f %f %f\n",
						       liabIdx + 1, pen_DD,
						       pen_Dd, pen_dd);
					    }

					  exit (-1);
					}

				      for (pedIdx = 0;
					   pedIdx < pedigreeSet.numPedigree;
					   pedIdx++)
					{
					  /* save the likelihood at null */
					  pPedigree =
					    pedigreeSet.ppPedigreeSet[pedIdx];
					  pedigreeSet.nullLikelihood[pedIdx] =
					    pPedigree->likelihood;
					}

				      log10_likelihood_null =
					pedigreeSet.log10Likelihood;
				      for (dprimeIdx = 0;
					   dprimeIdx < pLambdaCell->ndprime;
					   dprimeIdx++)
					{
					  if (modelOptions.equilibrium !=
					      LINKAGE_EQUILIBRIUM)
					    {
					      copy_dprime (pLDLoci,
							   pLambdaCell->
							   lambda[dprimeIdx]);
					      setup_LD_haplotype_freq
						(pLDLoci);
					    }
					  for (thetaInd = 0;
					       thetaInd < modelRange.ntheta;
					       thetaInd++)
					    {

					      if (modelOptions.mapFlag == SA)
						{
						  theta[0] =
						    modelRange.
						    theta[0][thetaInd];
						  theta[1] =
						    modelRange.
						    theta[1][thetaInd];
						  for (k = 0; k < 3; k++)
						    {
						      locusList->
							pNextLocusDistance[k]
							[0] = theta[0];
						      locusList->
							pPrevLocusDistance[k]
							[1] = theta[0];
						    }
						}
					      else
						{
						  locusList->
						    pNextLocusDistance
						    [MAP_MALE][0] =
						    locusList->
						    pPrevLocusDistance
						    [MAP_MALE][1] =
						    modelRange.
						    theta[0][thetaInd];
						  locusList->
						    pNextLocusDistance
						    [MAP_FEMALE][0] =
						    locusList->
						    pPrevLocusDistance
						    [MAP_FEMALE][1] =
						    modelRange.
						    theta[1][thetaInd];
						}

#ifndef NO_POLYNOMIAL
					      if (modelOptions.polynomial ==
						  TRUE)
						;
					      else
						/* populate the matrix */
						status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
										   0,	/* cell index */
										   -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
										   0,	/* number of het loci */
										   0);	/* current locus - start with 0 */
#else
					      /* populate the matrix */
					      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
										 0,	/* cell index */
										 -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
										 0,	/* number of het loci */
										 0);	/* current locus - start with 0 */
#endif

					      compute_likelihood
						(&pedigreeSet);
					      log10_likelihood_alternative =
						pedigreeSet.log10Likelihood;
					      if (pedigreeSet.likelihood ==
						  0.0
						  && pedigreeSet.
						  log10Likelihood == -9999.99)
						{
						  log10_likelihood_ratio = 0;
						}
					      else
						{
						  log10_likelihood_ratio =
						    log10_likelihood_alternative
						    - log10_likelihood_null;
						}
					      /* check for overflow problem !!! */
					      if (log10_likelihood_ratio >=
						  __DBL_MAX_10_EXP__ - 1)
						{
						  likelihood_ratio =
						    __DBL_MAX__;
						  tp_result[dprimeIdx]
						    [thetaInd][mkrFreqIdx].
						    lr_total = __DBL_MAX__;
						}
					      else
						/* check for underflow problem too !!! */
					      if (log10_likelihood_ratio <=
						    __DBL_MIN_10_EXP__ + 1)
						{
						  likelihood_ratio = 0;
						}
					      else
						{
						  likelihood_ratio =
						    pow (10.0,
							 log10_likelihood_ratio);
						  tp_result[dprimeIdx]
						    [thetaInd][mkrFreqIdx].
						    lr_total +=
						    likelihood_ratio;
						}
					      tp_result[dprimeIdx][thetaInd]
						[mkrFreqIdx].lr_count++;
					      /* caculating the HET */
					      for (j = 0;
						   j < modelRange.nalpha; j++)
						{
						  alphaV =
						    modelRange.alpha[j];
						  alphaV2 = 1 - alphaV;
						  if (alphaV2 < 0)
						    alphaV2 = 0;
						  log10HetLR = 0;
						  for (pedIdx = 0;
						       pedIdx <
						       pedigreeSet.
						       numPedigree; pedIdx++)
						    {
						      pPedigree =
							pedigreeSet.
							ppPedigreeSet[pedIdx];
						      homoLR =
							pPedigree->
							likelihood /
							pedigreeSet.
							nullLikelihood
							[pedIdx];
						      log10HetLR +=
							log10 (alphaV *
							       homoLR +
							       alphaV2);
						    }
						  if (log10HetLR >=
						      __DBL_MAX_10_EXP__ - 1)
						    {
						      hetLR = __DBL_MAX__;
						      tp_result[dprimeIdx]
							[thetaInd]
							[mkrFreqIdx].
							het_lr_total =
							__DBL_MAX__;
						    }
						  else
						    if (log10HetLR <=
							__DBL_MIN_10_EXP__ +
							1)
						    {
						      hetLR = 0;
						    }
						  else
						    {
						      adjustedHetLR = hetLR =
							pow (10, log10HetLR);
						      /* for threshold parameter, we need to make sure the weighting is even */
						      if (1 || modelType.distrib ==
							  QT_FUNCTION_CHI_SQUARE)
							{
							  if (modelRange.
							      ntthresh == 1)
							    {
							      adjustedHetLR *=
								2 *
								(modelType.
								 max -
								 modelType.
								 min);
							    }
							  else
							    if (thresholdIdx
								==
								modelRange.
								ntthresh - 1)
							    {
							      adjustedHetLR *=
								(2 *
								 modelType.
								 max -
								 threshold -
								 modelRange.
								 tthresh[0]
								 [thresholdIdx
								  - 1]);
							    }
							  else
							    if (thresholdIdx
								== 0)
							    {
							      adjustedHetLR *=
								(threshold +
								 modelRange.
								 tthresh[0]
								 [thresholdIdx
								  + 1] -
								 2 *
								 modelType.
								 min);
							    }
							  else
							    {

							      adjustedHetLR *=
								(modelRange.
								 tthresh[0]
								 [thresholdIdx
								  + 1] -
								 modelRange.
								 tthresh[0]
								 [thresholdIdx
								  - 1]);
							    }
							}
						      tp_result[dprimeIdx]
							[thetaInd]
							[mkrFreqIdx].
							het_lr_total +=
							adjustedHetLR;
						    }
						  if (tp_result[dprimeIdx]
						      [thetaInd][mkrFreqIdx].
						      max_penIdx < 0
						      || hetLR >
						      tp_result[dprimeIdx]
						      [thetaInd][mkrFreqIdx].
						      max_lr)
						    {
						      tp_result[dprimeIdx]
							[thetaInd]
							[mkrFreqIdx].max_lr =
							hetLR;
						      tp_result[dprimeIdx]
							[thetaInd]
							[mkrFreqIdx].
							max_alpha = alphaV;
						      tp_result[dprimeIdx]
							[thetaInd]
							[mkrFreqIdx].
							max_gfreq = gfreq;
						      tp_result[dprimeIdx]
							[thetaInd]
							[mkrFreqIdx].
							max_penIdx = penIdx;
						      tp_result[dprimeIdx]
							[thetaInd]
							[mkrFreqIdx].
							max_paramIdx =
							paramIdx;
						      tp_result[dprimeIdx]
							[thetaInd]
							[mkrFreqIdx].
							max_thresholdIdx =
							thresholdIdx;
						      tp_result[dprimeIdx]
							[thetaInd]
							[mkrFreqIdx].
							R_square = R_square;
						      tp_result[dprimeIdx]
							[thetaInd]
							[mkrFreqIdx].max_mf =
							mkrFreq;
						    }
						}

					    }	/* end of theta */
					}	/* end of D prime */
				      if (modelOptions.markerAnalysis !=
					  FALSE)
					break;
				    }	/* end of threshold loop */
				  if (modelOptions.markerAnalysis != FALSE)
				    break;
				}	/* end of penetrance loop */
			      if (modelOptions.markerAnalysis != FALSE)
				break;
			    }	/* end of parameter loop */
			  if (modelOptions.markerAnalysis != FALSE)
			    break;
			}	/* end of gene freq */
		    }		/* end of QT */
		  /* only loop marker allele frequencies when doing LD */
		  if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
		    break;
		  /* we can only do SNPs when looping over marker allele frequency */
		  if (pLocus2->numOriginalAllele > 2)
		    break;
		}		/* end of marker allele frequency looping */

	      /* calculate the average BR */
	      get_average_LR (tp_result);

#if 0
	      fprintf (fpHomo, "# %-d  \"%s %s \" \n", loc2, pLocus2->sName,
		       pLocus1->sName);
	      for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime;
		   dprimeIdx++)
		{
		  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
		    {
		      /* male theta */
		      theta[0] = modelRange.theta[0][thetaInd];
		      /* female theta */
		      theta[1] = modelRange.theta[1][thetaInd];
		      if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
			{
			  fprintf (fpHomo, "\t (%f,%f)  %f(%d)\n",
				   theta[0], theta[1],
				   tp_result[dprimeIdx][thetaInd][0].
				   lr_total /
				   tp_result[dprimeIdx][thetaInd][0].lr_count,
				   tp_result[dprimeIdx][thetaInd][0].
				   lr_count);
			}
		      else
			{
			  fprintf (fpHomo, "\t %f (%f,%f)  %f(%d)\n",
				   pLambdaCell->lambda[dprimeIdx][0][0],
				   theta[0], theta[1],
				   tp_result[dprimeIdx][thetaInd][0].
				   lr_total /
				   tp_result[dprimeIdx][thetaInd][0].lr_count,
				   tp_result[dprimeIdx][thetaInd][0].
				   lr_count);
			}
		    }
		}
	      fprintf (fpHomo, "-	Total 1234(1234)\n");
	      fflush (fpHomo);
#endif
	      /* for each D prime and theta, print out average and maximizing model information - MOD */
	      fprintf (fpHet, "# %-d  %s %s \n", loc2, pLocus1->sName,
		       pLocus2->sName);
	      fprintf (fpHet,
		       "%5s %16s %6s %8s %8s %6s %4s %6s %6s %5s %5s %5s \n",
		       "DPrime",
		       "Theta(M, F)",
		       "COUNT", "AVG_LR", "MAX_HLOD", "R2", "ALPHA", "DGF",
		       "MF", "PEN_DD", "PEN_Dd", "PEN_dd");
	      for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime;
		   dprimeIdx++)
		{
		  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
		    {
		      theta[0] = modelRange.theta[0][thetaInd];
		      theta[1] = modelRange.theta[1][thetaInd];
		      max =
			log10 (tp_result[dprimeIdx][thetaInd]
			       [modelRange.nafreq].max_lr);
		      gfreq =
			tp_result[dprimeIdx][thetaInd][modelRange.nafreq].
			max_gfreq;
		      alphaV =
			tp_result[dprimeIdx][thetaInd][modelRange.nafreq].
			max_alpha;
		      penIdx =
			tp_result[dprimeIdx][thetaInd][modelRange.nafreq].
			max_penIdx;
		      paramIdx =
			tp_result[dprimeIdx][thetaInd][modelRange.nafreq].
			max_paramIdx;
		      thresholdIdx =
			tp_result[dprimeIdx][thetaInd][modelRange.nafreq].
			max_thresholdIdx;
		      R_square =
			tp_result[dprimeIdx][thetaInd][modelRange.nafreq].
			R_square;
		      fprintf (fpHet,
			       "%6.2f (%6.4f, %6.4f) %6d %8.4f %8.4f %6.4f %5.2f %6.4f %6.4f ",
			       pLambdaCell->lambda[dprimeIdx][0][0], theta[0],
			       theta[1],
			       modelRange.nalpha *
			       tp_result[dprimeIdx][thetaInd][modelRange.
							      nafreq].
			       lr_count,
			       tp_result[dprimeIdx][thetaInd][modelRange.
							      nafreq].
			       het_lr_avg, max,
			       tp_result[dprimeIdx][thetaInd][modelRange.
							      nafreq].
			       R_square, alphaV, gfreq,
			       tp_result[dprimeIdx][thetaInd][modelRange.
							      nafreq].max_mf);
		      for (liabIdx = 0; liabIdx < modelRange.nlclass;
			   liabIdx++)
			{
			  pen_DD = modelRange.penet[liabIdx][0][penIdx];
			  pen_Dd = modelRange.penet[liabIdx][1][penIdx];
			  pen_dd = modelRange.penet[liabIdx][2][penIdx];
			  fprintf (fpHet, " %5.3f %5.3f %5.3f ", pen_DD,
				   pen_Dd, pen_dd);
			  if (modelType.trait != DT
			      && modelType.distrib != QT_FUNCTION_CHI_SQUARE)
			    {
			      SD_DD =
				modelRange.param[liabIdx][0][0][paramIdx];
			      SD_Dd =
				modelRange.param[liabIdx][1][0][paramIdx];
			      SD_dd =
				modelRange.param[liabIdx][2][0][paramIdx];
			      fprintf (fpHet, " %5.3f %5.3f %5.3f ", SD_DD,
				       SD_Dd, SD_dd);
			    }
			  if (modelType.trait != DT)
			    {
			      threshold =
				modelRange.tthresh[liabIdx][thresholdIdx];
			      fprintf (fpHet, " %5.3f ", threshold);
			    }
			}
		      fprintf (fpHet, "\n");
		      fflush (fpHet);

		    }		/* theta loop */

		}		/* dprime loop */
	      fprintf (stderr, "# %-d  %s %s Max Het LR\n", loc2,
		       pLocus2->sName, pLocus1->sName);
	      initialFlag = 1;
	      max = -99999;
	      max_at_theta0 = -99999;
	      max_at_dprime0 = -99999;
	      for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime;
		   dprimeIdx++)
		{
		  dprime = pLambdaCell->lambda[dprimeIdx][0][0];
		  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
		    {
		      theta[0] = modelRange.theta[0][thetaInd];
		      theta[1] = modelRange.theta[1][thetaInd];
		      lr =
			tp_result[dprimeIdx][thetaInd][modelRange.nafreq].
			max_lr;
		      if (initialFlag || lr > max)
			{
			  /* overall max */
			  max = lr;
			  maxDPrimeIdx = dprimeIdx;
			  maxThetaIdx = thetaInd;
			}
		      if (initialFlag || (-ERROR_MARGIN <= theta[0]
					  && theta[0] <= ERROR_MARGIN
					  && -ERROR_MARGIN <= theta[1]
					  && theta[1] <= ERROR_MARGIN))
			{
			  /* find the max for models with theta equal to 0 */
			  theta0Idx = thetaInd;
			  if (lr > max_at_theta0)
			    {
			      max_at_theta0 = lr;
			      maxDPrimeIdx_at_theta0 = dprimeIdx;
			    }
			}
		      if (-ERROR_MARGIN <= dprime && dprime <= ERROR_MARGIN)
			{
			  /* find the max for models with D prime equal to 0 */
			  dprime0Idx = dprimeIdx;
			  if (initialFlag || maxTheta_at_dprime0 < 0
			      || lr > max_at_dprime0)
			    {
			      max_at_dprime0 = lr;
			      maxTheta_at_dprime0 = thetaInd;
			    }
			}

		      initialFlag = 0;
		    }
		  initialFlag = 0;
		}
	      fprintf (stderr,
		       "Chr     Marker   Position   MOD   DPrime Theta R2 ALPHA DGF MF PEN_DD PEN_Dd PEN_dd\n");
	      /* overal maximizing model - MOD */
	      fprintf (stderr, "# Overal MOD maximizing model:\n");
	      theta[0] = modelRange.theta[0][maxThetaIdx];
	      theta[1] = modelRange.theta[1][maxThetaIdx];
	      gfreq =
		tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].
		max_gfreq;
	      mkrFreq =
		tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].
		max_mf;
	      alphaV =
		tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].
		max_alpha;
	      penIdx =
		tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].
		max_penIdx;
	      R_square =
		tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].
		R_square;
	      paramIdx =
		tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].
		max_paramIdx;
	      thresholdIdx =
		tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].
		max_thresholdIdx;
	      fprintf (stderr,
		       "%4d %15s %8.4f %8.4f %5.2f (%6.4f, %6.4f) %5.3f %4.2f %6.4f %6.4f",
		       pLocus2->pMapUnit->chromosome, pLocus2->sName,
		       pLocus2->pMapUnit->mapPos[SEX_AVERAGED], log10 (max),
		       pLambdaCell->lambda[maxDPrimeIdx][0][0], theta[0],
		       theta[1], R_square, alphaV, gfreq, mkrFreq);
	      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
		{
		  pen_DD = modelRange.penet[liabIdx][0][penIdx];
		  pen_Dd = modelRange.penet[liabIdx][1][penIdx];
		  pen_dd = modelRange.penet[liabIdx][2][penIdx];
		  fprintf (stderr, " %5.3f %5.3f %5.3f ", pen_DD, pen_Dd,
			   pen_dd);
		  if (modelType.trait != DT
		      && modelType.distrib != QT_FUNCTION_CHI_SQUARE)
		    {
		      SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
		      SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
		      SD_dd = modelRange.param[liabIdx][2][0][paramIdx];
		      fprintf (stderr, " %5.3f %5.3f %5.3f ", SD_DD,
			       SD_Dd, SD_dd);
		    }
		  if (modelType.trait != DT)
		    {
		      threshold = modelRange.tthresh[liabIdx][thresholdIdx];
		      fprintf (stderr, " %5.3f ", threshold);
		    }
		}
	      fprintf (stderr, "\n");
	      fflush (stderr);

	      /* maximizing model at theta equal to 0 - MOD */
	      fprintf (stderr, "# MOD maximizing model for theta=0:\n");
	      gfreq =
		tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.
							     nafreq].
		max_gfreq;
	      mkrFreq =
		tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.
							     nafreq].max_mf;
	      alphaV =
		tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.
							     nafreq].
		max_alpha;
	      penIdx =
		tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.
							     nafreq].
		max_penIdx;
	      R_square =
		tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.
							     nafreq].R_square;
	      paramIdx =
		tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.
							     nafreq].
		max_paramIdx;
	      thresholdIdx =
		tp_result[maxDPrimeIdx_at_theta0][theta0Idx][modelRange.
							     nafreq].
		max_thresholdIdx;
	      fprintf (stderr,
		       "%4d %15s %8.4f %8.4f %5.2f %6.4f %5.3f %4.2f %6.4f %6.4f",
		       pLocus2->pMapUnit->chromosome, pLocus2->sName,
		       pLocus2->pMapUnit->mapPos[SEX_AVERAGED],
		       log10 (max_at_theta0),
		       pLambdaCell->lambda[maxDPrimeIdx_at_theta0][0][0], 0.0,
		       R_square, alphaV, gfreq, mkrFreq);
	      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
		{
		  pen_DD = modelRange.penet[liabIdx][0][penIdx];
		  pen_Dd = modelRange.penet[liabIdx][1][penIdx];
		  pen_dd = modelRange.penet[liabIdx][2][penIdx];
		  fprintf (stderr, " %5.3f %5.3f %5.3f ", pen_DD, pen_Dd,
			   pen_dd);
		  if (modelType.trait != DT
		      && modelType.distrib != QT_FUNCTION_CHI_SQUARE)
		    {
		      SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
		      SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
		      SD_dd = modelRange.param[liabIdx][2][0][paramIdx];
		      fprintf (stderr, " %5.3f %5.3f %5.3f ", SD_DD,
			       SD_Dd, SD_dd);
		    }
		  if (modelType.trait != DT)
		    {
		      threshold = modelRange.tthresh[liabIdx][thresholdIdx];
		      fprintf (stderr, " %5.3f ", threshold);
		    }
		}
	      fprintf (stderr, "\n");
	      fflush (stderr);

	      /* maximizing model at d prime equal to 0 - MOD */
	      fprintf (stderr, "# MOD maximizing model for dprime=0:\n");
	      gfreq =
		tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].
		max_gfreq;
	      mkrFreq =
		tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].
		max_mf;
	      alphaV =
		tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].
		max_alpha;
	      penIdx =
		tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].
		max_penIdx;
	      R_square =
		tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].
		R_square;
	      paramIdx =
		tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].
		max_paramIdx;
	      thresholdIdx =
		tp_result[dprime0Idx][maxTheta_at_dprime0][modelRange.nafreq].
		max_thresholdIdx;
	      fprintf (stderr,
		       "%4d %15s %8.4f %8.4f %5.2f %6.4f %5.3f %4.2f %6.4f %6.4f ",
		       pLocus2->pMapUnit->chromosome, pLocus2->sName,
		       pLocus2->pMapUnit->mapPos[SEX_AVERAGED],
		       log10 (max_at_dprime0), 0.0, 0.0, R_square, alphaV,
		       gfreq, mkrFreq);
	      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
		{
		  pen_DD = modelRange.penet[liabIdx][0][penIdx];
		  pen_Dd = modelRange.penet[liabIdx][1][penIdx];
		  pen_dd = modelRange.penet[liabIdx][2][penIdx];
		  fprintf (stderr, " %5.3f %5.3f %5.3f ", pen_DD, pen_Dd,
			   pen_dd);
		  if (modelType.trait != DT
		      && modelType.distrib != QT_FUNCTION_CHI_SQUARE)
		    {
		      SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
		      SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
		      SD_dd = modelRange.param[liabIdx][2][0][paramIdx];
		      fprintf (stderr, " %5.3f %5.3f %5.3f ", SD_DD,
			       SD_Dd, SD_dd);
		    }
		  if (modelType.trait != DT)
		    {
		      threshold = modelRange.tthresh[liabIdx][thresholdIdx];
		      fprintf (stderr, " %5.3f ", threshold);
		    }
		}
	      fprintf (stderr, "\n");
	      fflush (stderr);

	      /* find the overal maximizing theta and dprime - LR
	       * with the other parameter integrated out */
	      max = -9999.99;
	      max_at_dprime0 = -9999.99;
	      max_at_theta0 = -9999.99;
	      for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime;
		   dprimeIdx++)
		{
		  dprime = pLambdaCell->lambda[dprimeIdx][0][0];
		  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
		    {
		      theta[0] = modelRange.theta[0][thetaInd];
		      theta[1] = modelRange.theta[1][thetaInd];
		      lr = tp_result[dprimeIdx][thetaInd][0].het_lr_avg;
		      if (lr > max)
			{
			  max = lr;
			  maxThetaIdx = thetaInd;
			  maxDPrimeIdx = dprimeIdx;
			}
		      if (-ERROR_MARGIN <= theta[0]
			  && theta[0] <= ERROR_MARGIN
			  && -ERROR_MARGIN <= theta[1]
			  && theta[1] <= ERROR_MARGIN)
			{
			  if (lr > max_at_theta0)
			    {
			      max_at_theta0 = lr;
			      maxDPrimeIdx_at_theta0 = dprimeIdx;
			    }
			}
		      if (-ERROR_MARGIN <= dprime && dprime <= ERROR_MARGIN)
			{
			  /* find the max for models with D prime equal to 0 */
			  dprime0Idx = dprimeIdx;
			  if (lr > max_at_dprime0)
			    {
			      max_at_dprime0 = lr;
			      maxTheta_at_dprime0 = thetaInd;
			    }
			}

		    }
		}
	      /* overal maximizing model - LR */
	      fprintf (stderr, "# Overal LR maximizing model:\n");
	      theta[0] = modelRange.theta[0][maxThetaIdx];
	      theta[1] = modelRange.theta[1][maxThetaIdx];
	      //gfreq = tp_result[maxDPrimeIdx][maxThetaIdx][modelRange.nafreq].max_gfreq;
	      fprintf (stderr,
		       "%4d %15s %8.4f %8.4f %5.2f (%6.4f %6.4f)\n",
		       pLocus2->pMapUnit->chromosome, pLocus2->sName,
		       pLocus2->pMapUnit->mapPos[SEX_AVERAGED], log10 (max),
		       pLambdaCell->lambda[maxDPrimeIdx][0][0], theta[0],
		       theta[1]);
	      fflush (stderr);

	      /* maximizing model at theta equal to 0 - LR */
	      fprintf (stderr, "# LR maximizing model for theta (0, 0):\n");
	      fprintf (stderr,
		       "%4d %15s %8.4f %8.4f %5.2f %6.4f\n",
		       pLocus2->pMapUnit->chromosome, pLocus2->sName,
		       pLocus2->pMapUnit->mapPos[SEX_AVERAGED],
		       log10 (max_at_theta0),
		       pLambdaCell->lambda[maxDPrimeIdx_at_theta0][0][0],
		       0.0);
	      fflush (stderr);

	      /* maximizing model at d prime equal to 0 - LR */
	      fprintf (stderr, "# LR maximizing model for dprime=0:\n");
	      fprintf (stderr,
		       "%4d %15s %8.4f %8.4f %5.2f %6.4f\n",
		       pLocus2->pMapUnit->chromosome, pLocus2->sName,
		       pLocus2->pMapUnit->mapPos[SEX_AVERAGED],
		       log10 (max_at_dprime0), 0.0, 0.0);
	      fflush (stderr);


	      /* output PPL now */
	      /* chromosome, marker name, position, PPL */
	      ppl = calculate_PPL (tp_result[dprime0Idx]);
	      fprintf (fpPPL, "%4d %15s %9.4f %6.4f ",
		       pLocus2->pMapUnit->chromosome, pLocus2->sName,
		       pLocus2->pMapUnit->mapPos[SEX_AVERAGED], ppl);
	      fflush (fpPPL);
	      /* output LD-PPL now if needed */
	      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
		{
		  /* calculate the LD LR average first */
		  get_average_LD_LR (tp_result);
		  /* calculate the LD-PPL - posterior probability of linkage allowing for LD */
		  ldppl = calculate_PPL (tp_result[pLambdaCell->ndprime]);
		  /* now calculate the PPLD - posterior probability of LD given linkage */
		  ppld = calculate_PPLD (tp_result);
		  fprintf (fpPPL, "%6.4f %6.4f ", ldppl, ppld);

		}
	      fprintf (fpPPL, "\n");
	      fflush (fpPPL);

	      /* free two point result storage */
	      free_tp_result_storage (prevNumDPrime);
	      tp_result = NULL;
	      prevNumDPrime = pLambdaCell->ndprime;
	      /* need to clear polynomial */
#ifndef NO_POLYNOMIAL

	      if (modelOptions.polynomial == TRUE)
		{
		  pedigreeSetPolynomialClearance (&pedigreeSet);
		  partialPolynomialClearance ();
		}
#endif


	      if (modelOptions.markerAnalysis == ADJACENTMARKER)
		loc2 = originalLocusList.numLocus;
	    }			/* end of looping second locus - loc2 */
	  /* if we are doing trait marker, then we are done */
	  /* Used to read: modelOptions.markerToMarker != TRUE which
	     is the same as markerAnalysis == FALSE as long as the old
	     markerToMarker and adjacentMarker flags were truly
	     orthogonal. Otherwise, it should be markerAnalysis !=
	     ADJACENTMARKER. */
	  if (modelOptions.markerAnalysis == FALSE)
	    loc1 = originalLocusList.numLocus;
	}			/* end of looping first locus - loc1 */
    }				/* end of two point */
  else
    {
      /* multipoint */
      nullLocusList.maxNumLocus = modelType.numMarkers + 1;
      nullLocusList.pLocusIndex =
	(int *) calloc (nullLocusList.maxNumLocus, sizeof (int));
      for (k = 0; k < 3; k++)
	{
	  nullLocusList.pPrevLocusDistance[k] =
	    (double *) calloc (nullLocusList.maxNumLocus, sizeof (double));
	  nullLocusList.pNextLocusDistance[k] =
	    (double *) calloc (nullLocusList.maxNumLocus, sizeof (double));
	}
      nullLocusList.numLocus = modelType.numMarkers + 1;
      savedLocusList.maxNumLocus = modelType.numMarkers + 1;
      savedLocusList.pLocusIndex =
	(int *) calloc (savedLocusList.maxNumLocus, sizeof (int));
      for (k = 0; k < 3; k++)
	{
	  savedLocusList.pPrevLocusDistance[k] =
	    (double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
	  savedLocusList.pNextLocusDistance[k] =
	    (double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
	}


#ifndef NO_POLYNOMIAL
      if (modelOptions.polynomial == TRUE)
	{
	  locusList = &savedLocusList;
	  xmissionMatrix = altMatrix;
	  /* populate the matrix */
	  status = populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr,	/* probability */
					     0,	/* cell index */
					     -1,	/* last het locus */
					     -1,	/* last het pattern (P-1 or M-2) */
					     0,	/* number of het loci */
					     0);	/* current locus - start with 0 */
	  locusList = &nullLocusList;
	  xmissionMatrix = nullMatrix;
	  /* populate the matrix */
	  status = populate_xmission_matrix (nullMatrix, totalLoci, initialProbAddr,	/* probability */
					     0,	/* cell index */
					     -1,	/* last he locus */
					     -1,	/* last het pattern (P-1 or M-2) */
					     0,	/* number of het loci */
					     0);	/* current locus - start with 0 */
	  makePolynomialStamp ();
	}
#endif



      /* get the trait locations we need to evaluate at */
      numPositions = modelRange.ntloc;
      mp_result =
	(SUMMARY_STAT *) calloc (numPositions, sizeof (SUMMARY_STAT));
      /* Need to output the results */
      fprintf (fpHet,
	       "pos PPL avgLR(count) MOD Alpha dgf pen_vector markerList\n");
      prevFirstMarker = -1;
      prevLastMarker = -1;
      prevTraitInd = -1;
      leftMarker = -1;
      for (posIdx = 0; posIdx < numPositions; posIdx++)
	{
	  /* positions listed are sex average positions */
	  traitPos = modelRange.tloc[posIdx];
	  /* set the sex average position first 
	   * the sex specific positions will be updated once markers are selected
	   * as interpolation might be needed
	   */
	  pTraitLocus->mapPosition[0] = traitPos;
	  pTraitLocus->mapPosition[1] = traitPos;
	  pTraitLocus->mapPosition[2] = traitPos;
	  /* initialize the locusList */
	  locusList = &savedLocusList;
	  memset (locusList->pLocusIndex, 0,
		  sizeof (int) * locusList->maxNumLocus);
	  for (k = 0; k < 3; k++)
	    {
	      memset (&locusList->pPrevLocusDistance[k][0], 0,
		      sizeof (double) * locusList->maxNumLocus);
	      memset (&locusList->pNextLocusDistance[k][0], 0,
		      sizeof (double) * locusList->maxNumLocus);
	    }
	  locusList->numLocus = 1;
	  locusList->pLocusIndex[0] = traitLocus;
	  for (k = 0; k < 3; k++)
	    {
	      locusList->pPrevLocusDistance[k][0] = -1;
	      locusList->pNextLocusDistance[k][0] = -1;
	    }
	  /* select markers to be used for the multipoint analysis */
	  add_markers_to_locuslist (locusList, modelType.numMarkers,
				    &leftMarker, 0,
				    originalLocusList.numLocus - 1, traitPos,
				    0);

	  /* store the markers used */
	  mp_result[posIdx].pMarkers =
	    (int *) calloc (modelType.numMarkers, sizeof (int));
	  k = 0;		/* marker index */
	  for (i = 0; i < locusList->numLocus; i++)
	    {
	      j = locusList->pLocusIndex[i];
	      if (originalLocusList.ppLocusList[j]->locusType ==
		  LOCUS_TYPE_MARKER)
		{
		  mp_result[posIdx].pMarkers[k] = j;
		  k++;
		}
	      else
		{
		  mp_result[posIdx].trait = i;
		  traitIndex = i;
		}
	    }
	  markerSetChanged = FALSE;
	  if (prevFirstMarker != mp_result[posIdx].pMarkers[0] ||
	      prevLastMarker !=
	      mp_result[posIdx].pMarkers[modelType.numMarkers - 1])
	    {
	      /* marker set has changed */
	      markerSetChanged = TRUE;
	    }
	  prevFirstMarker = mp_result[posIdx].pMarkers[0];
	  prevLastMarker =
	    mp_result[posIdx].pMarkers[modelType.numMarkers - 1];
	  if (markerSetChanged || prevTraitInd != mp_result[posIdx].trait)
	    locusListChanged = TRUE;
	  else
	    locusListChanged = FALSE;
	  prevTraitInd = mp_result[posIdx].trait;

	  /* interpolate trait postion for sex specific analysis */
	  if (modelOptions.mapFlag == SS)
	    {
	      if (traitIndex == 0)
		{
		  marker1Pos = get_map_position (locusList->pLocusIndex[1]);
		  /* trait is the first one in the list */
		  if (traitPos < ERROR_MARGIN && traitPos > -ERROR_MARGIN)
		    {
		      /* trait is at 0 */
		      pTraitLocus->mapPosition[MAP_MALE] =
			pTraitLocus->mapPosition[MAP_FEMALE] = 0;
		    }
		  else
		    {
		      /* get the relative position on the sex average map */
		      relativePos = traitPos / marker1Pos[0];
		      pTraitLocus->mapPosition[MAP_MALE] =
			relativePos * marker1Pos[MAP_MALE];
		      pTraitLocus->mapPosition[MAP_FEMALE] =
			relativePos * marker1Pos[MAP_FEMALE];
		    }
		  /* update the inter locus distance - sex averaged already done before */
		  for (k = 1; k < 3; k++)
		    {
		      locusList->pNextLocusDistance[k][0] =
			locusList->pPrevLocusDistance[k][1] =
			cm_to_recombination_fraction (marker1Pos[k] -
						      pTraitLocus->
						      mapPosition[k],
						      map.mapFunction);
		    }
		}
	      else if (traitIndex == modelType.numMarkers)
		{
		  /* trait is the last one in the list */
		  marker1Pos =
		    get_map_position (locusList->
				      pLocusIndex[modelType.numMarkers - 2]);
		  marker2Pos =
		    get_map_position (locusList->
				      pLocusIndex[modelType.numMarkers - 1]);
		  /* get the relative position on the sex average map */
		  dist = marker2Pos[0] - marker1Pos[0];
		  if (dist > ERROR_MARGIN)
		    {
		      relativePos = (traitPos - marker2Pos[0]) / dist;
		      pTraitLocus->mapPosition[MAP_MALE] =
			relativePos * (marker2Pos[MAP_MALE] -
				       marker1Pos[MAP_MALE]) +
			marker2Pos[MAP_MALE];
		      pTraitLocus->mapPosition[MAP_FEMALE] =
			relativePos * (marker2Pos[MAP_FEMALE] -
				       marker1Pos[MAP_FEMALE]) +
			marker2Pos[MAP_FEMALE];
		    }
		  else
		    {
		      pTraitLocus->mapPosition[MAP_MALE] =
			traitPos - marker2Pos[0] + marker2Pos[MAP_MALE];
		      pTraitLocus->mapPosition[MAP_FEMALE] =
			traitPos - marker2Pos[0] + marker2Pos[MAP_FEMALE];
		    }

		  /* update the inter locus distance - sex averaged already done before */
		  for (k = 1; k <= 2; k++)
		    {
		      locusList->pNextLocusDistance[k][traitIndex - 1] =
			locusList->pPrevLocusDistance[k][traitIndex] =
			cm_to_recombination_fraction (pTraitLocus->
						      mapPosition[k] -
						      marker2Pos[k],
						      map.mapFunction);
		    }

		}
	      else
		{
		  /* trait is in between two markers */
		  marker1Pos =
		    get_map_position (locusList->pLocusIndex[traitIndex - 1]);
		  marker2Pos =
		    get_map_position (locusList->pLocusIndex[traitIndex + 1]);
		  /* get the relative position on the sex average map */
		  dist = marker2Pos[0] - marker1Pos[0];
		  if (dist > ERROR_MARGIN)
		    {
		      relativePos = (traitPos - marker1Pos[0]) / dist;
		      pTraitLocus->mapPosition[MAP_MALE] =
			relativePos * (marker2Pos[MAP_MALE] -
				       marker1Pos[MAP_MALE]) +
			marker1Pos[MAP_MALE];
		      pTraitLocus->mapPosition[MAP_FEMALE] =
			relativePos * (marker2Pos[MAP_FEMALE] -
				       marker1Pos[MAP_FEMALE]) +
			marker1Pos[MAP_FEMALE];

		    }
		  else
		    {
		      pTraitLocus->mapPosition[MAP_MALE] =
			marker1Pos[MAP_MALE];
		      pTraitLocus->mapPosition[MAP_FEMALE] =
			marker1Pos[MAP_FEMALE];
		    }
		  /* update the inter locus distance - sex averaged already done before */
		  for (k = 1; k < 3; k++)
		    {
		      locusList->pNextLocusDistance[k][traitIndex - 1] =
			locusList->pPrevLocusDistance[k][traitIndex] =
			cm_to_recombination_fraction (pTraitLocus->
						      mapPosition[k] -
						      marker1Pos[k],
						      map.mapFunction);
		      locusList->pNextLocusDistance[k][traitIndex] =
			locusList->pPrevLocusDistance[k][traitIndex + 1] =
			cm_to_recombination_fraction (marker2Pos[k] -
						      pTraitLocus->
						      mapPosition[k],
						      map.mapFunction);

		    }
		}

	    }

	  /* NULL hypothesis - the trait locus is not even on the map
	   * so we can set the trait at the left of the chromosome, with 0.5 recombination 
	   * to the first marker in the list 
	   * Build the locusList for the null hypothesis */
	  if (markerSetChanged == TRUE)
	    {
	      nullLocusList.pLocusIndex[0] = traitLocus;
	      for (k = 0; k < 3; k++)
		{
		  nullLocusList.pPrevLocusDistance[k][0] = -1;
		  nullLocusList.pNextLocusDistance[k][0] = 0.5;
		  nullLocusList.pPrevLocusDistance[k][1] = 0.5;
		}
	      nullLocusList.pLocusIndex[1] = mp_result[posIdx].pMarkers[0];
	      prevPos = get_map_position (nullLocusList.pLocusIndex[1]);
	      for (k = 1; k < modelType.numMarkers; k++)
		{
		  nullLocusList.pLocusIndex[k + 1] =
		    mp_result[posIdx].pMarkers[k];
		  currPos =
		    get_map_position (nullLocusList.pLocusIndex[k + 1]);
		  for (j = 0; j < 3; j++)
		    {
		      nullLocusList.pPrevLocusDistance[j][k + 1] =
			nullLocusList.pNextLocusDistance[j][k] =
			cm_to_recombination_fraction (currPos[j] - prevPos[j],
						      map.mapFunction);
		    }
		  if (modelOptions.mapFlag == SA)
		    {
		      for (j = 1; j <= 2; j++)
			{
			  nullLocusList.pPrevLocusDistance[j][k + 1] =
			    nullLocusList.pNextLocusDistance[j][k] =
			    nullLocusList.pPrevLocusDistance[0][k + 1];
			}
		    }
		  prevPos = currPos;
		}
	      for (j = 0; j < 3; j++)
		nullLocusList.pNextLocusDistance[j][modelType.numMarkers] =
		  -1;
	    }

	  /* the locus list has been built, go on to the analysis 
	   * multipoint DT */
	  if (pTrait->type == DICHOTOMOUS)
	    {
	      if (markerSetChanged == TRUE)
		{
		  locusList = &nullLocusList;
		  xmissionMatrix = nullMatrix;
#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      pedigreeSetPolynomialClearance (&pedigreeSet);
		      partialPolynomialClearance ();
		    }
		  else
		    {
		      /* populate the matrix */
		      status = populate_xmission_matrix (nullMatrix, totalLoci, initialProbAddr,	/* probability */
							 0,	/* cell index */
							 -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
							 0,	/* number of het loci */
							 0);	/* current locus - start with 0 */
		    }
#else
		  /* populate the matrix */
		  status = populate_xmission_matrix (nullMatrix, totalLoci, initialProbAddr,	/* probability */
						     0,	/* cell index */
						     -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
						     0,	/* number of het loci */
						     0);	/* current locus - start with 0 */
#endif

		  for (penIdx = 0; penIdx < modelRange.npenet; penIdx++)
		    {
		      for (liabIdx = 0; liabIdx < modelRange.nlclass;
			   liabIdx++)
			{
			  pen_DD = modelRange.penet[liabIdx][0][penIdx];
			  pen_Dd = modelRange.penet[liabIdx][1][penIdx];
			  pen_dd = modelRange.penet[liabIdx][2][penIdx];
			  pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
			  pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
			  pTrait->penetrance[2][liabIdx][1][0] = pen_Dd;
			  pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
			  pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
			  pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
			  pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_Dd;
			  pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;
			}


#ifndef NO_POLYNOMIAL
		      if (modelOptions.polynomial == TRUE)
			;
		      else
			/* only need to update trait locus */
			update_penetrance (&pedigreeSet, traitLocus);
#else
		      /* only need to update trait locus */
		      update_penetrance (&pedigreeSet, traitLocus);
#endif

		      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq;
			   gfreqInd++)
			{
			  /* updated trait locus allele frequencies */
			  gfreq = modelRange.gfreq[gfreqInd];
			  pLocus->pAlleleFrequency[0] = gfreq;
			  pLocus->pAlleleFrequency[1] = 1 - gfreq;


#ifndef NO_POLYNOMIAL
			  if (modelOptions.polynomial == TRUE)
			    ;
			  else
			    update_locus (&pedigreeSet, traitLocus);
#else
			  update_locus (&pedigreeSet, traitLocus);
#endif
			  /* get the likelihood at NULL hypothesis - use nullLocusList */
			  locusList = &nullLocusList;
			  xmissionMatrix = nullMatrix;
			  compute_likelihood (&pedigreeSet);
			  //printAllVariables(); 
			  //fprintf(stderr," Null Likelihood=%e log10Likelihood=%e \n",
			  //pedigreeSet.likelihood,pedigreeSet.log10Likelihood); 


			  if (pedigreeSet.likelihood == 0.0 &&
			      pedigreeSet.log10Likelihood == -9999.99)
			    {
			      fprintf (stderr, "NULL has likelihood 0\n");
			      fprintf (stderr, "dgf=%f\n", gfreq);
			      for (liabIdx = 0; liabIdx < modelRange.nlclass;
				   liabIdx++)
				{
				  pen_DD =
				    modelRange.penet[liabIdx][0][penIdx];
				  pen_Dd =
				    modelRange.penet[liabIdx][1][penIdx];
				  pen_dd =
				    modelRange.penet[liabIdx][2][penIdx];
				  fprintf (stderr,
					   "Liab %d penentrance %f %f %f\n",
					   liabIdx + 1, pen_DD, pen_Dd,
					   pen_dd);
				}

			      exit (-1);
			    }
			  /* save the results for NULL */
			  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree;
			       pedIdx++)
			    {
			      /* save the likelihood at null */
			      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
			      pedigreeSet.nullLikelihood[pedIdx] =
				pPedigree->likelihood;
			      likelihoodDT[pedIdx][gfreqInd][penIdx] =
				pPedigree->likelihood;
			    }

			  log10_likelihood_null = pedigreeSet.log10Likelihood;
			  likelihoodDT[pedigreeSet.
				       numPedigree][gfreqInd][penIdx] =
			    log10_likelihood_null;
			}	/* gfreq */
		    }		/* pen */

		}		/* end of marker set changed */

	      locusList = &savedLocusList;
	      xmissionMatrix = altMatrix;
	      if (markerSetChanged || locusListChanged)
		{
#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      pedigreeSetPolynomialClearance (&pedigreeSet);
		      partialPolynomialClearance ();
		    }
#endif
		}
#ifndef NO_POLYNOMIAL
	      if (modelOptions.polynomial == TRUE)
		;
	      else
		/* populate the matrix */
		status = populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr,	/* probability */
						   0,	/* cell index */
						   -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
						   0,	/* number of het loci */
						   0);	/* current locus - start with 0 */
#else
	      /* populate the matrix */
	      status = populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr,	/* probability */
						 0,	/* cell index */
						 -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
						 0,	/* number of het loci */
						 0);	/* current locus - start with 0 */

#endif


	      /* for alternative */
	      for (penIdx = 0; penIdx < modelRange.npenet; penIdx++)
		{
		  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
		    {
		      pen_DD = modelRange.penet[liabIdx][0][penIdx];
		      pen_Dd = modelRange.penet[liabIdx][1][penIdx];
		      pen_dd = modelRange.penet[liabIdx][2][penIdx];
		      pTrait->penetrance[2][liabIdx][0][0] = pen_DD;
		      pTrait->penetrance[2][liabIdx][0][1] = pen_Dd;
		      pTrait->penetrance[2][liabIdx][1][0] = pen_Dd;
		      pTrait->penetrance[2][liabIdx][1][1] = pen_dd;
		      pTrait->penetrance[1][liabIdx][0][0] = 1 - pen_DD;
		      pTrait->penetrance[1][liabIdx][0][1] = 1 - pen_Dd;
		      pTrait->penetrance[1][liabIdx][1][0] = 1 - pen_Dd;
		      pTrait->penetrance[1][liabIdx][1][1] = 1 - pen_dd;
		    }


#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    ;
		  else
		    /* only need to update trait locus */
		    update_penetrance (&pedigreeSet, traitLocus);
#else
		  update_penetrance (&pedigreeSet, traitLocus);
#endif
		  for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++)
		    {
		      /* updated trait locus allele frequencies */
		      gfreq = modelRange.gfreq[gfreqInd];
		      pLocus->pAlleleFrequency[0] = gfreq;
		      pLocus->pAlleleFrequency[1] = 1 - gfreq;


#ifndef NO_POLYNOMIAL
		      if (modelOptions.polynomial == TRUE)
			;
		      else
			update_locus (&pedigreeSet, traitLocus);
#else
		      update_locus (&pedigreeSet, traitLocus);
#endif

		      /* ready for the alternative hypothesis */
		      locusList = &savedLocusList;
		      xmissionMatrix = altMatrix;
		      compute_likelihood (&pedigreeSet);

		      //fprintf(stderr," Alternative Likelihood=%e log10Likelihood=%e\n",
		      //pedigreeSet.likelihood,pedigreeSet.log10Likelihood);


		      log10_likelihood_alternative =
			pedigreeSet.log10Likelihood;
		      if (pedigreeSet.likelihood == 0.0
			  && pedigreeSet.log10Likelihood == -9999.99)
			{
			  log10_likelihood_ratio = 0;
			}
		      else
			{
			  log10_likelihood_ratio =
			    log10_likelihood_alternative -
			    likelihoodDT[pedigreeSet.
					 numPedigree][gfreqInd][penIdx];

			}
		      /* check for overflow problem !!! */
		      if (log10_likelihood_ratio >= __DBL_MAX_10_EXP__ - 1)
			{
			  likelihood_ratio = __DBL_MAX__;
			  mp_result[posIdx].lr_total += __DBL_MAX__;
			}
		      else
			/* check for underflow problem too !!! */
		      if (log10_likelihood_ratio <= __DBL_MIN_10_EXP__ + 1)
			{
			  likelihood_ratio = 0;
			}
		      else
			{
			  likelihood_ratio =
			    pow (10.0, log10_likelihood_ratio);
			  mp_result[posIdx].lr_total += likelihood_ratio;
			}
		      /* add the result to the right placeholder */
		      mp_result[posIdx].lr_count++;

		      /* caculating the HET */
		      for (j = 0; j < modelRange.nalpha; j++)
			{
			  alphaV = modelRange.alpha[j];
			  alphaV2 = 1 - alphaV;
			  if (alphaV2 < 0)
			    alphaV2 = 0;
			  log10HetLR = 0;
			  for (pedIdx = 0;
			       pedIdx < pedigreeSet.numPedigree; pedIdx++)
			    {
			      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
			      homoLR = pPedigree->likelihood /
				likelihoodDT[pedIdx][gfreqInd][penIdx];
			      if (alphaV * homoLR + alphaV2 < 0)
				fprintf (stderr,
					 "HET LR less than 0. Check!!!\n");
			      log10HetLR += log10 (alphaV * homoLR + alphaV2);
			    }
			  if (log10HetLR >= __DBL_MAX_10_EXP__ - 1)
			    {
			      hetLR = __DBL_MAX__;
			      mp_result[posIdx].het_lr_total = __DBL_MAX__;
			    }
			  else if (log10HetLR <= __DBL_MIN_10_EXP__ + 1)
			    {
			      hetLR = 0;
			    }
			  else
			    {
			      hetLR = pow (10, log10HetLR);
			      mp_result[posIdx].het_lr_total += hetLR;
			    }
			  if (mp_result[posIdx].max_penIdx <
			      0 || hetLR > mp_result[posIdx].max_lr)
			    {
			      mp_result[posIdx].max_lr = hetLR;
			      mp_result[posIdx].max_alpha = alphaV;
			      mp_result[posIdx].max_gfreq = gfreq;
			      mp_result[posIdx].max_penIdx = penIdx;
			    }
			}	/* end of calculating HET LR */
		    }		/* end of genFreq loop */
		}		/* end of penetrance loop */
	    }			/* end of TP */
	  else
	    /* multipoint QT or COMBINED */
	    {
#ifndef NO_POLYNOMIAL
	      if (markerSetChanged)
		{
		  if (modelOptions.polynomial == TRUE)
		    {
		      pedigreeSetPolynomialClearance (&pedigreeSet);
		      partialPolynomialClearance ();
		    }
#endif
		  for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++)
		    {
		      gfreq = modelRange.gfreq[gfreqInd];
		      pLocus->pAlleleFrequency[0] = gfreq;
		      pLocus->pAlleleFrequency[1] = 1 - gfreq;

		      update_locus (&pedigreeSet, traitLocus);
		      /* this should be MEAN + SD */
		      for (paramIdx = 0; paramIdx < modelRange.nparam;
			   paramIdx++)
			{
			  for (penIdx = 0; penIdx < modelRange.npenet;
			       penIdx++)
			    {
			      breakFlag = FALSE;
			      for (thresholdIdx = 0;
				   thresholdIdx < modelRange.ntthresh;
				   thresholdIdx++)
				{
				  for (liabIdx = 0;
				       liabIdx < modelRange.nlclass;
				       liabIdx++)
				    {
				      mean_DD =
					modelRange.penet[liabIdx][0][penIdx];
				      mean_Dd =
					modelRange.penet[liabIdx][1][penIdx];
				      mean_dd =
					modelRange.penet[liabIdx][2][penIdx];
				      SD_DD =
					modelRange.
					param[liabIdx][0][0][paramIdx];
				      SD_Dd =
					modelRange.
					param[liabIdx][1][0][paramIdx];
				      SD_dd =
					modelRange.
					param[liabIdx][2][0][paramIdx];
				      threshold =
					modelRange.
					tthresh[liabIdx][thresholdIdx];

				      /* check against the hard coded constraint */
				      if (modelType.distrib !=
					  QT_FUNCTION_CHI_SQUARE)
					{
					  constraint =
					    pow (1.0 - gfreq,
						 2) * mean_dd * SD_dd +
					    2 * gfreq * (1 -
							 gfreq) * mean_Dd *
					    SD_Dd + pow (gfreq,
							 2) * mean_DD * SD_DD;
					  /*      fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
					     constraint, gfreq, mean_DD, SD_DD, 
					     mean_Dd, SD_DD, 
					     mean_dd, SD_dd);
					   */
					  if (constraint >= 3.0
					      || constraint <= -3.0)
					    {
					      breakFlag = TRUE;
					      break;
					    }
					}
				      pTrait->means[liabIdx][0][0] = mean_DD;
				      pTrait->means[liabIdx][0][1] = mean_Dd;
				      pTrait->means[liabIdx][1][0] = mean_Dd;
				      pTrait->means[liabIdx][1][1] = mean_dd;
				      pTrait->stddev[liabIdx][0][0] = SD_DD;
				      pTrait->stddev[liabIdx][0][1] = SD_Dd;
				      pTrait->stddev[liabIdx][1][0] = SD_Dd;
				      pTrait->stddev[liabIdx][1][1] = SD_dd;

				      /* threshold for QT */
				      pTrait->cutoffValue[liabIdx] =
					threshold;

				    }	/* liability class Index */
				  if (breakFlag == TRUE)
				    continue;
				  update_penetrance (&pedigreeSet,
						     traitLocus);

				  /* get the likelihood at NULL hypothesis - use nullLocusList */
				  locusList = &nullLocusList;
				  xmissionMatrix = nullMatrix;

#ifndef NO_POLYNOMIAL
				  if (modelOptions.polynomial == TRUE)
				    ;
				  else
				    /* populate the matrix */
				    status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
								       0,	/* cell index */
								       -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
								       0,	/* number of het loci */
								       0);	/* current locus - start with 0 */
#else
				  /* populate the matrix */
				  status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
								     0,	/* cell index */
								     -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
								     0,	/* number of het loci */
								     0);	/* current locus - start with 0 */
#endif
				  compute_likelihood (&pedigreeSet);
				  if (pedigreeSet.likelihood == 0.0 &&
				      pedigreeSet.log10Likelihood == -9999.99)
				    {
				      fprintf (stderr,
					       "NULL has likelihood 0\n");
				      fprintf (stderr, "dgf=%f\n", gfreq);
				      for (liabIdx = 0;
					   liabIdx < modelRange.nlclass;
					   liabIdx++)
					{
					  pen_DD =
					    modelRange.
					    penet[liabIdx][0][penIdx];
					  pen_Dd =
					    modelRange.
					    penet[liabIdx][1][penIdx];
					  pen_dd =
					    modelRange.
					    penet[liabIdx][2][penIdx];
					  fprintf (stderr,
						   "Liab %d penentrance %f %f %f\n",
						   liabIdx + 1, pen_DD,
						   pen_Dd, pen_dd);
					}

				      exit (-1);
				    }

				  for (pedIdx = 0;
				       pedIdx < pedigreeSet.numPedigree;
				       pedIdx++)
				    {
				      /* save the likelihood at null */
				      pPedigree =
					pedigreeSet.ppPedigreeSet[pedIdx];
				      pedigreeSet.nullLikelihood[pedIdx] =
					pPedigree->likelihood;
				      likelihoodQT[pedIdx][gfreqInd][penIdx]
					[paramIdx][thresholdIdx] =
					pPedigree->likelihood;
				    }

				  log10_likelihood_null =
				    pedigreeSet.log10Likelihood;
				  if (isnan (log10_likelihood_null))
				    fprintf (stderr,
					     "NULL likelihood is NAN.\n");

				  likelihoodQT[pedigreeSet.
					       numPedigree][gfreqInd][penIdx]
				    [paramIdx][thresholdIdx] =
				    log10_likelihood_null;
#ifndef NO_POLYNOMIAL
				}	/* thresholdIdx */
			    }	/* penIdx */
			}	/* paramIdx */
		    }		/* gfreq */
		}		/* markerSetChanged */

	      if (markerSetChanged || locusListChanged)
		{
		  if (modelOptions.polynomial == TRUE)
		    {
		      pedigreeSetPolynomialClearance (&pedigreeSet);
		      partialPolynomialClearance ();
		    }
		}

	      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++)
		{
		  gfreq = modelRange.gfreq[gfreqInd];
		  pLocus->pAlleleFrequency[0] = gfreq;
		  pLocus->pAlleleFrequency[1] = 1 - gfreq;

		  update_locus (&pedigreeSet, traitLocus);
		  /* this should be MEAN + SD */
		  for (paramIdx = 0; paramIdx < modelRange.nparam; paramIdx++)
		    {
		      for (penIdx = 0; penIdx < modelRange.npenet; penIdx++)
			{
			  breakFlag = FALSE;
			  for (thresholdIdx = 0;
			       thresholdIdx < modelRange.ntthresh;
			       thresholdIdx++)
			    {
			      for (liabIdx = 0; liabIdx < modelRange.nlclass;
				   liabIdx++)
				{
				  mean_DD =
				    modelRange.penet[liabIdx][0][penIdx];
				  mean_Dd =
				    modelRange.penet[liabIdx][1][penIdx];
				  mean_dd =
				    modelRange.penet[liabIdx][2][penIdx];
				  SD_DD =
				    modelRange.param[liabIdx][0][0][paramIdx];
				  SD_Dd =
				    modelRange.param[liabIdx][1][0][paramIdx];
				  SD_dd =
				    modelRange.param[liabIdx][2][0][paramIdx];
				  threshold =
				    modelRange.tthresh[liabIdx][thresholdIdx];

				  if (modelType.distrib !=
				      QT_FUNCTION_CHI_SQUARE)
				    {
				      /* check against the hard coded constraint */
				      constraint =
					pow (1.0 - gfreq,
					     2) * mean_dd * SD_dd +
					2 * gfreq * (1 -
						     gfreq) * mean_Dd *
					SD_Dd + pow (gfreq,
						     2) * mean_DD * SD_DD;
				      /*      fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
				         constraint, gfreq, mean_DD, SD_DD, 
				         mean_Dd, SD_DD, 
				         mean_dd, SD_dd);
				       */
				      if (constraint >= 3.0
					  || constraint <= -3.0)
					{
					  breakFlag = TRUE;
					  break;
					}
				    }
				  pTrait->means[liabIdx][0][0] = mean_DD;
				  pTrait->means[liabIdx][0][1] = mean_Dd;
				  pTrait->means[liabIdx][1][0] = mean_Dd;
				  pTrait->means[liabIdx][1][1] = mean_dd;
				  pTrait->stddev[liabIdx][0][0] = SD_DD;
				  pTrait->stddev[liabIdx][0][1] = SD_Dd;
				  pTrait->stddev[liabIdx][1][0] = SD_Dd;
				  pTrait->stddev[liabIdx][1][1] = SD_dd;

				  /* threshold for QT */
				  pTrait->cutoffValue[liabIdx] = threshold;

				}	/* liability class Index */
			      if (breakFlag == TRUE)
				continue;
			      update_penetrance (&pedigreeSet, traitLocus);

#endif

			      /* ready for the alternative hypothesis */
			      locusList = &savedLocusList;
			      xmissionMatrix = altMatrix;
#ifndef NO_POLYNOMIAL
			      if (modelOptions.polynomial == TRUE)
				;
			      else
				status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
								   0,	/* cell index */
								   -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
								   0,	/* number of het loci */
								   0);	/* current locus - start with 0 */
#else
			      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
								 0,	/* cell index */
								 -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
								 0,	/* number of het loci */
								 0);	/* current locus - start with 0 */
#endif
			      compute_likelihood (&pedigreeSet);
			      log10_likelihood_alternative =
				pedigreeSet.log10Likelihood;
			      if (isnan (log10_likelihood_alternative))
				fprintf (stderr, "ALT likelihood is NAN.\n");
			      if (pedigreeSet.likelihood == 0.0
				  && pedigreeSet.log10Likelihood == -9999.99)
				{
				  log10_likelihood_ratio = 0;
				}
			      else
				{
				  log10_likelihood_ratio =
				    log10_likelihood_alternative -
				    likelihoodQT[pedigreeSet.
						 numPedigree][gfreqInd]
				    [penIdx][paramIdx][thresholdIdx];
				}
			      /* check for overflow problem !!! */
			      if (log10_likelihood_ratio >=
				  __DBL_MAX_10_EXP__ - 1)
				{
				  likelihood_ratio = __DBL_MAX__;
				  mp_result[posIdx].lr_total += __DBL_MAX__;
				}
			      else
				/* check for underflow problem too !!! */
			      if (log10_likelihood_ratio <=
				    __DBL_MIN_10_EXP__ + 1)
				{
				  likelihood_ratio = 0;
				}
			      else
				{
				  likelihood_ratio =
				    pow (10.0, log10_likelihood_ratio);
				  mp_result[posIdx].lr_total +=
				    likelihood_ratio;
				}
			      /* add the result to the right placeholder */
			      mp_result[posIdx].lr_count++;

			      if (isnan (likelihood_ratio))
				fprintf (stderr,
					 "LR for the pedigree set is NAN.\n");
			      /* caculating the HET */
			      for (j = 0; j < modelRange.nalpha; j++)
				{
				  alphaV = modelRange.alpha[j];
				  alphaV2 = 1 - alphaV;
				  if (alphaV2 < 0)
				    alphaV2 = 0;
				  log10HetLR = 0;
				  for (pedIdx = 0;
				       pedIdx < pedigreeSet.numPedigree;
				       pedIdx++)
				    {
				      pPedigree =
					pedigreeSet.ppPedigreeSet[pedIdx];
				      homoLR =
					pPedigree->likelihood /
					likelihoodQT[pedIdx][gfreqInd][penIdx]
					[paramIdx][thresholdIdx];
				      log10HetLR +=
					log10 (alphaV * homoLR + alphaV2);
				    }
				  if (log10HetLR >= __DBL_MAX_10_EXP__ - 1)
				    {
				      hetLR = __DBL_MAX__;
				      mp_result[posIdx].het_lr_total =
					__DBL_MAX__;
				    }
				  else
				    if (log10HetLR <= __DBL_MIN_10_EXP__ + 1)
				    {
				      hetLR = 0;
				    }
				  else
				    {
				      adjustedHetLR = hetLR =
					pow (10, log10HetLR);
				      /* for threshold parameter, we need to make sure the weighting is even */
				      if (1 || modelType.distrib ==
					  QT_FUNCTION_CHI_SQUARE)
					{
					  if (modelRange.ntthresh == 1)
					    {
					      adjustedHetLR *=
						2 * (modelType.max -
						     modelType.min);
					    }
					  else if (thresholdIdx ==
						   modelRange.ntthresh - 1)
					    {
					      adjustedHetLR *=
						(2 * modelType.max -
						 threshold -
						 modelRange.
						 tthresh[0][thresholdIdx -
							    1]);
					    }
					  else if (thresholdIdx == 0)
					    {
					      adjustedHetLR *=
						(threshold +
						 modelRange.
						 tthresh[0][thresholdIdx +
							    1] -
						 2 * modelType.min);
					    }
					  else
					    {

					      adjustedHetLR *=
						modelRange.
						tthresh[0][thresholdIdx + 1] -
						modelRange.
						tthresh[0][thresholdIdx - 1];
					    }
					}
				      mp_result[posIdx].het_lr_total +=
					adjustedHetLR;
				    }
				  if (mp_result[posIdx].
				      max_penIdx < 0
				      || hetLR > mp_result[posIdx].max_lr)
				    {
				      mp_result[posIdx].max_lr = hetLR;
				      mp_result[posIdx].max_alpha = alphaV;
				      mp_result[posIdx].max_gfreq = gfreq;
				      mp_result[posIdx].max_penIdx = penIdx;
				      mp_result[posIdx].max_paramIdx =
					paramIdx;
				      mp_result[posIdx].max_thresholdIdx =
					thresholdIdx;
				    }
				}
			    }	/* end of threshold loop */

			}	/* end of penetrance loop */
		    }		/* end of parameter loop */
		}		/* end of gene freq */
	    }			/* end of QT */

	  /* print out average and log10(max) and maximizing parameters */
	  //	  if (modelType.trait == DT || modelType.distrib != QT_FUNCTION_CHI_SQUARE)
	  if (modelType.trait == DT)
	    avgLR =
	      mp_result[posIdx].het_lr_total / (modelRange.nalpha *
						mp_result[posIdx].lr_count);
	  else
	    /* under QT CHISQ, threshold parameter has been evenly weighted */
	    avgLR =
	      mp_result[posIdx].het_lr_total / (modelRange.nalpha *
						(mp_result[posIdx].lr_count /
						 modelRange.ntthresh) * 2 *
						(modelType.max -
						 modelType.min));

	  ppl = (avgLR * avgLR) / (-5.77 + 54 * avgLR + avgLR * avgLR);
	  max = mp_result[posIdx].max_lr;
	  gfreq = mp_result[posIdx].max_gfreq;
	  alphaV = mp_result[posIdx].max_alpha;
	  penIdx = mp_result[posIdx].max_penIdx;
	  paramIdx = mp_result[posIdx].max_paramIdx;
	  thresholdIdx = mp_result[posIdx].max_thresholdIdx;
	  fprintf (fpHet, "\t %f  %6.4f %12.8f(%d) %10.6f %f %f ",
		   traitPos, ppl, avgLR,
		   mp_result[posIdx].lr_count, log10 (max), alphaV, gfreq);
	  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
	    {
	      pen_DD = modelRange.penet[liabIdx][0][penIdx];
	      pen_Dd = modelRange.penet[liabIdx][1][penIdx];
	      pen_dd = modelRange.penet[liabIdx][2][penIdx];
	      fprintf (fpHet, " %f %f %f ", pen_DD, pen_Dd, pen_dd);
	      if (modelType.trait != DT
		  && modelType.distrib != QT_FUNCTION_CHI_SQUARE)
		{
		  SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
		  SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
		  SD_dd = modelRange.param[liabIdx][2][0][paramIdx];
		  fprintf (fpHet, " %5.3f %5.3f %5.3f ", SD_DD,
			   SD_Dd, SD_dd);
		}
	      if (modelType.trait != DT)
		{
		  threshold = modelRange.tthresh[liabIdx][thresholdIdx];
		  fprintf (fpHet, " %5.3f ", threshold);
		}
	    }
	  /* print out markers used for this position */
	  fprintf (fpHet, "(%d", mp_result[posIdx].pMarkers[0]);
	  for (k = 1; k < modelType.numMarkers; k++)
	    {
	      fprintf (fpHet, ",%d", mp_result[posIdx].pMarkers[k]);
	    }
	  fprintf (fpHet, ")\n");
	  fflush (fpHet);
	}			/* end of walking down the chromosome */
    }				/* end of multipoint */


  time2 = clock ();


#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
//   polyStatistics (NULL);
//   dismantle();
    }
#endif

  if (modelType.type == TP)
    {
      if (tp_result != NULL)
	{
	  /* two point */
	  for (i = 0; i < pLambdaCell->ndprime; i++)
	    {
	      free (tp_result[i]);
	    }
	  free (tp_result);
	}
    }
  else
    {
      /* multipoint */
      for (posIdx = 0; posIdx < numPositions; posIdx++)
	{
	  free (mp_result[posIdx].pMarkers);
	}
      free (mp_result);
    }

  /* free parental pair work space */
  free_parental_pair_workspace (&parentalPairSpace, modelType.numMarkers + 1);

  /* free transmission probability matrix */
  free (altMatrix);
  free (nullMatrix);

  if (modelType.type == TP)
    free_LD_loci (pLDLoci);

#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      pedigreeSetPolynomialClearance (&pedigreeSet);
      polynomialClearance ();
    }
#endif
  free_likelihood_storage (&pedigreeSet);
  free_likelihood_space (&pedigreeSet);
  free_pedigree_set (&pedigreeSet);
  free_sub_locus_list (&nullLocusList);
  free_sub_locus_list (&savedLocusList);
  free (modelOptions.sUnknownPersonID);
  final_cleanup ();


  fprintf (stderr, "Computation time:  %fs  %fs \n",
	   (double) (time1 - time0) / CLOCKS_PER_SEC,
	   (double) (time2 - time1) / CLOCKS_PER_SEC);


  /* close file pointers */
  if (modelType.type == TP)
    {
      fclose (fpPPL);
    }
  fclose (fpHet);
  //  fclose (fpHomo);
  return 0;
}
