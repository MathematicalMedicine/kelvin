/**********************************************************************
 * Multiprocessor Linkage Analysis
 * Alberto Maria Segre, Yungui Huang
 * RADSMM storage code Martin Milder
 * Regex code Nathan Burnette
 * 
 * Copyright 2006, The University of Iowa.  All rights reserved.
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
FILE *fpHomo = NULL;		/* average HOMO LR file */
FILE *fpPPL = NULL;		/* PPL output file */
FILE *fpLD = NULL;		/* LD PPL output file */

/* Model datastructures. modelOptions is defined in the pedigree library. */
ModelType modelType;
ModelRange modelRange;
ModelOptions modelOptions;

/* temporarily for GAW project */
/* summary statistics */
typedef struct SUMMARY_STAT
{
  /* for calculating average */
  double lr_total;		/* sum of log(LR) */
  int lr_count;			/* number of models for the (ld, theta) pair */
  double het_lr_avg;		/* average of HET log(LR) */
  double het_lr_total;		/* with heterogeneity - alpha, parameter */

  /* for max */
  double max_lr;		/* max het lr */
  double max_alpha;
  double max_gfreq;
  int max_penIdx;

  double R_square;		/* only applies in bi-allelic, i.e. SNPs */

  /* for MP marker list */
  int *pMarkers;
  double ppl;			/* imputed MP ppl based on average likelihood ratio */
  int trait;			/* trait locus index in the locus list for MP */
} SUMMARY_STAT;

/* number of D primes 
 * if there are more than 2 alleles in the marker/trait, number of D primes
 * and D prime ranges are assumed to be the same to reduce complexity 
 * for initial phase of this project */
int num_of_d_prime;
double *d_prime;
int num_of_theta;
double *theta;
/* two dimensional array for the two point summary results *
 * first dimension is the D prime, for LE, D prime=0 with just one element
 * in this dimension 
 * second dimension is theta values */
SUMMARY_STAT **tp_result;

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

LambdaCell *pLambdaCell = NULL;
int prevNumDPrime = 0;

void free_likelihood_storage ();
int free_tp_result_storage (int ndprime);
int allocate_tp_result_storage ();
double calculate_R_square (double p1, double q1, double d);
double calculate_PPL (SUMMARY_STAT * result,
		      double thetaCutoff, double weight, double prior);
int get_average_LD_LR (SUMMARY_STAT ** result, int numDPrime, int numTheta);

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
  double theta;			/* theta */
  int penIdx, liabIdx, gfreqInd, thetaInd;
  int paramIdx;
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
  int markerSetChanged;		/* flag for multipoint analysis */
  int locusListChanged;		/* flag for multipoint analysis */
  int prevFirstMarker;		/* first marker in the set for multipoint analysis */
  int prevLastMarker;		/* last marker in the set for multipoint analysis */
  int prevTraitInd;
  double prevPos, currPos;	/* for MP */
  int locus;
  int thresholdIdx;
  double threshold;
  int R_square_flag = FALSE;
  double R_square = 0;
  double thetaCutoff = 0.05;	/* for calculating PPL, should be configurable */
  double thetaWeight = 0.95;	/* for calculating PPL, should be configurable */
  double prior = 0.02;		/* priori probability of linkage */
  int dprime0Idx = 0;
  int mkrFreqIdx;
  double mkrFreq;
  int maxThetaIdx = 0;
  int maxDPrimeIdx = 0;
  int maxDPrimeIdx_at_theta0 = 0;
  double max_at_theta0;
  double lr;
  int theta0Idx = 0;

  SubLocusList savedLocusList;
  SubLocusList nullLocusList;

  clock_t time0, time1, time2;
  int numberOfCompute = 0;

  memset (&savedLocusList, 0, sizeof (savedLocusList));
  memset (&nullLocusList, 0, sizeof (nullLocusList));

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
  fpHomo = fopen (avghomofile, "w");
  KASSERT (fpHomo != NULL, "Error in opening file %s for write.\n",
	   avghomofile);
  if (modelType.type == TP)
    {
      fpPPL = fopen (pplfile, "w");
      KASSERT (fpPPL != NULL, "Error in opening file %s for write.\n",
	       pplfile);
      //setbuf(fpPPL, NULL);
      fprintf (fpPPL, "CHR MARKER   cM    PPL \n");
    }
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
    {
      fpLD = fopen (ldPPLfile, "w");
      KASSERT (fpPPL != NULL, "Error in opening file %s for write.\n",
	       ldPPLfile);
      fprintf (fpLD, "CHR MARKER   cM    PPL \n");
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
	  pTrait->functionQT = modelType.distrib;
	  if (modelType.distrib == T_DISTRIBUTION)
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

  /* initialize loci by doing genotype elimination, set recoding */
  initialize_loci (&pedigreeSet);

#ifndef NO_POLYNOMIAL
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

  if (modelOptions.polynomial == TRUE)
    {
      makePolynomialStamp ();
    }
#endif

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
      savedLocusList.pPrevLocusDistance = (double *) malloc (sizeof (double) *
							     savedLocusList.
							     numLocus);
      savedLocusList.pNextLocusDistance =
	(double *) malloc (sizeof (double) * savedLocusList.numLocus);

      savedLocusList.pPrevLocusDistance[0] = -1;
      savedLocusList.pNextLocusDistance[1] = -1;

      for (loc1 = 0; loc1 < originalLocusList.numLocus - 1; loc1++)
	{
	  savedLocusList.pLocusIndex[0] = loc1;
	  pLocus1 = originalLocusList.ppLocusList[loc1];
	  if (modelOptions.markerAnalysis == MARKERTOMARKER
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

	      /* allocate result storage */
	      allocate_tp_result_storage ();

	      /* we will force marker allele frequency loop to execute at least once */
	      for (mkrFreqIdx = 0;
		   mkrFreqIdx == 0 || mkrFreqIdx < modelRange.nafreq;
		   mkrFreqIdx++)
		{
		  /* we should only loop over marker allele frequency under twopoint
		   * and when markers are SNPs (only have two alleles) */
		  if (modelRange.nafreq >= 1
		      && modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM
		      && pLocus2->numOriginalAllele == 2)
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
		  if (pTrait->type == DICHOTOMOUS)
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
			      locusList->pNextLocusDistance[0] = 0.5;
			      locusList->pPrevLocusDistance[1] = 0.5;
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
				      theta = modelRange.theta[0][thetaInd];
				      locusList->pNextLocusDistance[0] =
					theta;
				      locusList->pPrevLocusDistance[1] =
					theta;

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
				      likelihood_ratio =
					pow (10.0, log10_likelihood_ratio);
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
					  hetLR = pow (10, log10HetLR);
					  tp_result[dprimeIdx][thetaInd].
					    het_lr_total += hetLR;
					  if (tp_result[dprimeIdx][thetaInd].
					      max_penIdx < 0
					      || hetLR >
					      tp_result[dprimeIdx][thetaInd].
					      max_lr)
					    {
					      tp_result[dprimeIdx][thetaInd].
						max_lr = hetLR;
					      tp_result[dprimeIdx][thetaInd].
						max_alpha = alphaV;
					      tp_result[dprimeIdx][thetaInd].
						max_gfreq = gfreq;
					      tp_result[dprimeIdx][thetaInd].
						max_penIdx = penIdx;
					      tp_result[dprimeIdx][thetaInd].
						R_square = R_square;
					    }
					}	/* end of calculating HET LR */
				      /* add the result to the right placeholder */
				      tp_result[dprimeIdx][thetaInd].
					lr_total += likelihood_ratio;
				      tp_result[dprimeIdx][thetaInd].
					lr_count++;
				      //fprintf(stderr, "likelihood ratio: %e.\n", likelihood_ratio);

				    }	/* end of theta loop */
				}	/* end of D prime loop */
			      if (modelOptions.markerAnalysis ==
				  MARKERTOMARKER)
				{
				  /* marker to marker analysis, marker allele frequency is fixed */
				  gfreqInd = modelRange.ngfreq;
				  break;
				}
			    }	/* end of genFreq loop */
			  if (modelOptions.markerAnalysis == MARKERTOMARKER)
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
					      constraint =
						pow (1.0 - gfreq,
						     2) * mean_dd * SD_dd +
						2 * gfreq * (1 -
							     gfreq) *
						mean_Dd * SD_Dd + pow (gfreq,
								       2) *
						mean_DD * SD_DD;
/*	  fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
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

				      locusList->pNextLocusDistance[0] = 0.5;
				      locusList->pPrevLocusDistance[1] = 0.5;
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
					      theta =
						modelRange.theta[0][thetaInd];
					      locusList->
						pNextLocusDistance[0] = theta;
					      locusList->
						pPrevLocusDistance[1] = theta;
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
					      likelihood_ratio =
						pow (10.0,
						     log10_likelihood_ratio);
					      /* caculating the HET */
					      for (j = 0;
						   j < modelRange.nalpha; j++)
						{
						  alphaV =
						    modelRange.alpha[j];
						  alphaV2 = 1 - alphaV;
						  if (alphaV2 < 0)
						    alphaV2 = 0;
						  hetLR = 0;
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
						      hetLR +=
							log10 (alphaV *
							       homoLR +
							       alphaV2);
						    }
						  hetLR = pow (10, hetLR);
						  tp_result[dprimeIdx]
						    [thetaInd].het_lr_total +=
						    hetLR;
						  if (tp_result[dprimeIdx]
						      [thetaInd].max_penIdx <
						      0
						      || hetLR >
						      tp_result[dprimeIdx]
						      [thetaInd].max_lr)
						    {
						      tp_result[dprimeIdx]
							[thetaInd].max_lr =
							hetLR;
						      tp_result[dprimeIdx]
							[thetaInd].max_alpha =
							alphaV;
						      tp_result[dprimeIdx]
							[thetaInd].max_gfreq =
							gfreq;
						      tp_result[dprimeIdx]
							[thetaInd].
							max_penIdx = penIdx;
						      tp_result[dprimeIdx]
							[thetaInd].R_square =
							R_square;
						    }
						}
					      /* add the result to the right placeholder */
					      tp_result[dprimeIdx][thetaInd].
						lr_total += likelihood_ratio;
					      tp_result[dprimeIdx][thetaInd].
						lr_count++;
					      //fprintf(stderr, "likelihood ratio: %e.\n", likelihood_ratio);

					    }	/* end of theta */
					}	/* end of D prime */
				      if (modelOptions.markerAnalysis ==
					  MARKERTOMARKER)
					break;
				    }	/* end of threshold loop */
				  if (modelOptions.markerAnalysis ==
				      MARKERTOMARKER)
				    break;
				}	/* end of penetrance loop */
			      if (modelOptions.markerAnalysis ==
				  MARKERTOMARKER)
				break;
			    }	/* end of parameter loop */
			  if (modelOptions.markerAnalysis == MARKERTOMARKER)
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
	      fprintf (fpHomo, "# %-d  \"%s %s \" \n", loc2, pLocus2->sName,
		       pLocus1->sName);
	      for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime;
		   dprimeIdx++)
		{
		  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
		    {
		      theta = modelRange.theta[0][thetaInd];
		      if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
			{
			  fprintf (fpHomo, "\t (%f,%f)  %f(%d)\n",
				   theta, theta,
				   tp_result[dprimeIdx][thetaInd].lr_total /
				   tp_result[dprimeIdx][thetaInd].lr_count,
				   tp_result[dprimeIdx][thetaInd].lr_count);
			}
		      else
			{
			  fprintf (fpHomo, "%f %d %f %f\n",
				   pLambdaCell->lambda[dprimeIdx][0][0],
				   tp_result[dprimeIdx][thetaInd].lr_count,
				   theta,
				   log10 (tp_result[dprimeIdx][thetaInd].
					  lr_total /
					  tp_result[dprimeIdx][thetaInd].
					  lr_count));
			}
		    }
		}
	      fprintf (fpHomo, "-	Total 1234(1234)\n");
	      fflush (fpHomo);

	      fprintf (fpHet, "# %-d  %s %s \n", loc2, pLocus1->sName,
		       pLocus2->sName);
	      fprintf (fpHet,
		       "DPrime Theta(M, F) COUNT AVG_LR MAX_HLOD R2 ALPHA DGF PEN_DD PEN_Dd PEN_dd\n");
	      for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime;
		   dprimeIdx++)
		{
		  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
		    {
		      theta = modelRange.theta[0][thetaInd];
		      max = log10 (tp_result[dprimeIdx][thetaInd].max_lr);
		      gfreq = tp_result[dprimeIdx][thetaInd].max_gfreq;
		      alphaV = tp_result[dprimeIdx][thetaInd].max_alpha;
		      penIdx = tp_result[dprimeIdx][thetaInd].max_penIdx;
		      tp_result[dprimeIdx][thetaInd].het_lr_avg =
			tp_result[dprimeIdx][thetaInd].
			het_lr_total / (modelRange.nalpha *
					tp_result[dprimeIdx]
					[thetaInd].lr_count);
		      fprintf (fpHet,
			       "%5.2f (%6.4f, %6.4f) %6d %8.4f %8.4f %6.4f %4.2f %6.4f ",
			       pLambdaCell->lambda[dprimeIdx][0][0], theta,
			       theta,
			       modelRange.nalpha *
			       tp_result[dprimeIdx][thetaInd].lr_count,
			       tp_result[dprimeIdx][thetaInd].het_lr_avg, max,
			       tp_result[dprimeIdx][thetaInd].R_square,
			       alphaV, gfreq);
		      for (liabIdx = 0; liabIdx < modelRange.nlclass;
			   liabIdx++)
			{
			  pen_DD = modelRange.penet[liabIdx][0][penIdx];
			  pen_Dd = modelRange.penet[liabIdx][1][penIdx];
			  pen_dd = modelRange.penet[liabIdx][2][penIdx];
			  fprintf (fpHet, " %5.3f %5.3f %5.3f ", pen_DD,
				   pen_Dd, pen_dd);
			}
		      fprintf (fpHet, "\n");
		      fflush (fpHet);

		    }

		}
	      fprintf (stderr, "# %-d  %s %s Max Het LR\n", loc2,
		       pLocus2->sName, pLocus1->sName);
	      max = -99999;
	      max_at_theta0 = -99999;
	      for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime;
		   dprimeIdx++)
		{
		  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
		    {
		      theta = modelRange.theta[0][thetaInd];
		      lr = tp_result[dprimeIdx][thetaInd].max_lr;
		      if (lr > max)
			{
			  /* overall max */
			  max = lr;
			  maxDPrimeIdx = dprimeIdx;
			  maxThetaIdx = thetaInd;
			}
		      if (-ERROR_MARGIN <= theta && theta <= ERROR_MARGIN)
			{
			  theta0Idx = thetaInd;
			  if (lr > max_at_theta0)
			    {
			      max_at_theta0 = lr;
			      maxDPrimeIdx_at_theta0 = dprimeIdx;
			    }
			}

		    }
		}
	      fprintf (stderr,
		       "Chr     Marker   Position   MOD   DPrime Theta R2 ALPHA DGF PEN_DD PEN_Dd PEN_dd\n");
	      theta = modelRange.theta[0][maxThetaIdx];
	      gfreq = tp_result[maxDPrimeIdx][maxThetaIdx].max_gfreq;
	      alphaV = tp_result[maxDPrimeIdx][maxThetaIdx].max_alpha;
	      penIdx = tp_result[maxDPrimeIdx][maxThetaIdx].max_penIdx;
	      R_square = tp_result[maxDPrimeIdx][maxThetaIdx].R_square;
	      fprintf (stderr,
		       "%4d %15s %8.4f %8.4f %5.2f %6.4f %5.3f %4.2f %6.4f",
		       pLocus2->pMapUnit->chromosome, pLocus2->sName,
		       pLocus2->pMapUnit->mapPos[SEX_AVERAGED], log10 (max),
		       pLambdaCell->lambda[maxDPrimeIdx][0][0], theta,
		       R_square, alphaV, gfreq);
	      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
		{
		  pen_DD = modelRange.penet[liabIdx][0][penIdx];
		  pen_Dd = modelRange.penet[liabIdx][1][penIdx];
		  pen_dd = modelRange.penet[liabIdx][2][penIdx];
		  fprintf (stderr, " %5.3f %5.3f %5.3f ", pen_DD, pen_Dd,
			   pen_dd);
		}
	      fprintf (stderr, "\n");
	      fflush (stderr);

	      gfreq = tp_result[maxDPrimeIdx_at_theta0][theta0Idx].max_gfreq;
	      alphaV = tp_result[maxDPrimeIdx_at_theta0][theta0Idx].max_alpha;
	      penIdx =
		tp_result[maxDPrimeIdx_at_theta0][theta0Idx].max_penIdx;
	      R_square =
		tp_result[maxDPrimeIdx_at_theta0][theta0Idx].R_square;
	      fprintf (stderr,
		       "%4d %15s %8.4f %8.4f %5.2f %6.4f %5.3f %4.2f %6.4f",
		       pLocus2->pMapUnit->chromosome, pLocus2->sName,
		       pLocus2->pMapUnit->mapPos[SEX_AVERAGED],
		       log10 (max_at_theta0),
		       pLambdaCell->lambda[maxDPrimeIdx_at_theta0][0][0], 0.0,
		       R_square, alphaV, gfreq);
	      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
		{
		  pen_DD = modelRange.penet[liabIdx][0][penIdx];
		  pen_Dd = modelRange.penet[liabIdx][1][penIdx];
		  pen_dd = modelRange.penet[liabIdx][2][penIdx];
		  fprintf (stderr, " %5.3f %5.3f %5.3f ", pen_DD, pen_Dd,
			   pen_dd);
		}
	      fprintf (stderr, "\n");
	      fflush (stderr);


#if 0
	      /* Print out maximizing model for each theta */

	      fprintf (stderr,
		       "Theta         LR      HLOD  alpha gfreq penDD.1 penDd.1 pendd.1 penDD.2 penDd.2 pendd.2 ... \n");
	      for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
		{
		  theta = modelRange.theta[0][thetaInd];
		  max = tp_result[0][thetaInd].max_lr;
		  gfreq = tp_result[0][thetaInd].max_gfreq;
		  alphaV = tp_result[0][thetaInd].max_alpha;
		  penIdx = tp_result[0][thetaInd].max_penIdx;
		  fprintf (stderr, "%f %10.6f %8.4f %f %f ",
			   theta, max, log10 (max), alphaV, gfreq);
		  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
		    {
		      pen_DD = modelRange.penet[liabIdx][0][penIdx];
		      pen_Dd = modelRange.penet[liabIdx][1][penIdx];
		      pen_dd = modelRange.penet[liabIdx][2][penIdx];
		      fprintf (stderr, " %f %f %f ", pen_DD, pen_Dd, pen_dd);
		    }
		  fprintf (stderr, "\n");
		}		/* end of looping theta for average output */
	      fflush (stderr);

#endif

	      /* output PPL now */
	      /* chromosome, marker name, position, PPL */
	      ppl =
		calculate_PPL (tp_result[dprime0Idx], thetaCutoff,
			       thetaWeight, prior);
	      fprintf (fpPPL, "%4d %s   %8.4f   %6.4f\n",
		       pLocus2->pMapUnit->chromosome, pLocus2->sName,
		       pLocus2->pMapUnit->mapPos[SEX_AVERAGED], ppl);
	      fflush (fpPPL);
	      /* output LD-PPL now if needed */
	      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
		{
		  /* calculate the LD LR average first */
		  get_average_LD_LR (tp_result, pLambdaCell->ndprime,
				     modelRange.ntheta);
		  ppl =
		    calculate_PPL (tp_result[pLambdaCell->ndprime],
				   thetaCutoff, thetaWeight, prior);
		  fprintf (fpLD, "%4d %s   %8.4f   %6.4f\n",
			   pLocus2->pMapUnit->chromosome, pLocus2->sName,
			   pLocus2->pMapUnit->mapPos[SEX_AVERAGED], ppl);
		  fflush (fpLD);

		}

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
      nullLocusList.pPrevLocusDistance =
	(double *) calloc (nullLocusList.maxNumLocus, sizeof (double));
      nullLocusList.pNextLocusDistance =
	(double *) calloc (nullLocusList.maxNumLocus, sizeof (double));
      nullLocusList.numLocus = modelType.numMarkers + 1;
      savedLocusList.maxNumLocus = modelType.numMarkers + 1;
      savedLocusList.pLocusIndex =
	(int *) calloc (savedLocusList.maxNumLocus, sizeof (int));
      savedLocusList.pPrevLocusDistance =
	(double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
      savedLocusList.pNextLocusDistance =
	(double *) calloc (savedLocusList.maxNumLocus, sizeof (double));


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
	  traitPos = modelRange.tloc[posIdx];
	  pTraitLocus->mapPosition = traitPos;
	  /* initialize the locusList */
	  locusList = &savedLocusList;
	  memset (locusList->pLocusIndex, 0,
		  sizeof (int) * locusList->maxNumLocus);
	  memset (locusList->pPrevLocusDistance, 0,
		  sizeof (double) * locusList->maxNumLocus);
	  memset (locusList->pNextLocusDistance, 0,
		  sizeof (double) * locusList->maxNumLocus);
	  locusList->numLocus = 1;
	  locusList->pLocusIndex[0] = traitLocus;
	  locusList->pPrevLocusDistance[0] = -1;
	  locusList->pNextLocusDistance[0] = -1;
	  /* select markers to be used for the multipoint analysis */
	  add_markers_to_locuslist (locusList, modelType.numMarkers,
				    &leftMarker, 0,
				    originalLocusList.numLocus - 1, traitPos,
				    modelOptions.mapFlag);

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
		mp_result[posIdx].trait = i;
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

	  /* NULL hypothesis - the trait locus is not even on the map
	   * so we can set the trait at the left of the chromosome, with 0.5 recombination 
	   * to the first marker in the list 
	   * Build the locusList for the null hypothesis */
	  if (markerSetChanged == TRUE)
	    {
	      nullLocusList.pLocusIndex[0] = traitLocus;
	      nullLocusList.pPrevLocusDistance[0] = -1;
	      nullLocusList.pNextLocusDistance[0] = 0.5;
	      nullLocusList.pLocusIndex[1] = mp_result[posIdx].pMarkers[0];
	      nullLocusList.pPrevLocusDistance[1] = 0.5;
	      prevPos =
		get_map_position (nullLocusList.pLocusIndex[1],
				  modelOptions.mapFlag);
	      for (k = 1; k < modelType.numMarkers; k++)
		{
		  nullLocusList.pLocusIndex[k + 1] =
		    mp_result[posIdx].pMarkers[k];
		  currPos =
		    get_map_position (nullLocusList.pLocusIndex[k + 1],
				      modelOptions.mapFlag);
		  nullLocusList.pPrevLocusDistance[k + 1] =
		    nullLocusList.pNextLocusDistance[k] =
		    cm_to_recombination_fraction (currPos - prevPos,
						  modelOptions.mapFlag);
		  prevPos = currPos;
		}
	      nullLocusList.pNextLocusDistance[modelType.numMarkers] = -1;
	    }

	  /* the locus list has been built, go on to the analysis 
	   * multipoint DT */

	  /* do NULL hypothesis first */
	  locusList = &nullLocusList;
	  //count_likelihood_space(&pedigreeSet);

	  if (pTrait->type == DICHOTOMOUS)
	    {
#ifndef NO_POLYNOMIAL
	      if (markerSetChanged == TRUE)
		{
		  if (modelOptions.polynomial == TRUE)
		    {
		      pedigreeSetPolynomialClearance (&pedigreeSet);
		      partialPolynomialClearance ();
		    }
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
			  //count_likelihood_space(&pedigreeSet);

			  /* NULL */
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
#ifndef NO_POLYNOMIAL
			      likelihoodDT[pedIdx][gfreqInd][penIdx] =
				pPedigree->likelihood;
#endif
			    }

			  log10_likelihood_null = pedigreeSet.log10Likelihood;
#ifndef NO_POLYNOMIAL
			  likelihoodDT[pedigreeSet.
				       numPedigree][gfreqInd][penIdx] =
			    log10_likelihood_null;
			}	/* gfreq */
		    }		/* pen */

		}		/* end of marker set changed */

	      if (markerSetChanged || locusListChanged)
		{
		  if (modelOptions.polynomial == TRUE)
		    {
		      pedigreeSetPolynomialClearance (&pedigreeSet);
		      partialPolynomialClearance ();
		    }
		}

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


		  if (modelOptions.polynomial == TRUE)
		    ;
		  else
		    /* only need to update trait locus */
		    update_penetrance (&pedigreeSet, traitLocus);

		  for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++)
		    {
		      /* updated trait locus allele frequencies */
		      gfreq = modelRange.gfreq[gfreqInd];
		      pLocus->pAlleleFrequency[0] = gfreq;
		      pLocus->pAlleleFrequency[1] = 1 - gfreq;


		      if (modelOptions.polynomial == TRUE)
			;
		      else
			update_locus (&pedigreeSet, traitLocus);
#endif

		      /* ready for the alternative hypothesis */
		      locusList = &savedLocusList;
		      //count_likelihood_space(&pedigreeSet);

		      compute_likelihood (&pedigreeSet);

//fprintf(stderr," Alternative Likelihood=%e log10Likelihood=%e\n",pedigreeSet.likelihood,pedigreeSet.log10Likelihood);


		      log10_likelihood_alternative =
			pedigreeSet.log10Likelihood;
		      if (pedigreeSet.likelihood == 0.0
			  && pedigreeSet.log10Likelihood == -9999.99)
			{
			  log10_likelihood_ratio = 0;
			}
		      else
			{
#ifndef NO_POLYNOMIAL
			  log10_likelihood_ratio =
			    log10_likelihood_alternative -
			    likelihoodDT[pedigreeSet.
					 numPedigree][gfreqInd][penIdx];

#else
			  log10_likelihood_ratio =
			    log10_likelihood_alternative -
			    log10_likelihood_null;
#endif
			}
		      likelihood_ratio = pow (10.0, log10_likelihood_ratio);
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
#ifndef NO_POLYNOMIAL
			      homoLR = pPedigree->likelihood /
				likelihoodDT[pedIdx][gfreqInd][penIdx];
#else
			      homoLR =
				pPedigree->likelihood /
				pedigreeSet.nullLikelihood[pedIdx];
#endif
			      log10HetLR += log10 (alphaV * homoLR + alphaV2);
			    }
			  hetLR = pow (10, log10HetLR);
			  mp_result[posIdx].het_lr_total += hetLR;
			  if (mp_result[posIdx].max_penIdx <
			      0 || hetLR > mp_result[posIdx].max_lr)
			    {
			      mp_result[posIdx].max_lr = hetLR;
			      mp_result[posIdx].max_alpha = alphaV;
			      mp_result[posIdx].max_gfreq = gfreq;
			      mp_result[posIdx].max_penIdx = penIdx;
			    }
			}	/* end of calculating HET LR */
		      /* add the result to the right placeholder */
		      mp_result[posIdx].lr_total += likelihood_ratio;
		      mp_result[posIdx].lr_count++;
		      //fprintf(stderr, "likelihood ratio: %e.\n", likelihood_ratio);

		    }		/* end of genFreq loop */
		}		/* end of penetrance loop */
	    }			/* end of TP */
	  else
	    /* multipoint should be QT or COMBINED */
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
#ifndef NO_POLYNOMIAL
				      likelihoodQT[pedIdx][gfreqInd][penIdx]
					[paramIdx][thresholdIdx] =
					pPedigree->likelihood;
#endif
				    }

				  log10_likelihood_null =
				    pedigreeSet.log10Likelihood;
				  if (isnan (log10_likelihood_null))
				    fprintf (stderr,
					     "NULL likelihood is NAN.\n");
#ifndef NO_POLYNOMIAL
				  likelihoodQT[pedigreeSet.
					       numPedigree][gfreqInd][penIdx]
				    [paramIdx][thresholdIdx] =
				    log10_likelihood_null;
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

				  /* check against the hard coded constraint */
				  constraint =
				    pow (1.0 - gfreq,
					 2) * mean_dd * SD_dd +
				    2 * gfreq * (1 -
						 gfreq) * mean_Dd * SD_Dd +
				    pow (gfreq, 2) * mean_DD * SD_DD;
				  /*      fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
				     constraint, gfreq, mean_DD, SD_DD, 
				     mean_Dd, SD_DD, 
				     mean_dd, SD_dd);
				   */
				  if (constraint >= 3.0 || constraint <= -3.0)
				    {
				      breakFlag = TRUE;
				      break;
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
#ifndef NO_POLYNOMIAL
				  log10_likelihood_ratio =
				    log10_likelihood_alternative -
				    likelihoodQT[pedigreeSet.
						 numPedigree][gfreqInd]
				    [penIdx][paramIdx][thresholdIdx];
#else
				  log10_likelihood_ratio =
				    log10_likelihood_alternative -
				    log10_likelihood_null;
#endif
				}
			      likelihood_ratio =
				pow (10.0, log10_likelihood_ratio);
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
				  hetLR = 0;
				  for (pedIdx = 0;
				       pedIdx < pedigreeSet.numPedigree;
				       pedIdx++)
				    {
				      pPedigree =
					pedigreeSet.ppPedigreeSet[pedIdx];
#ifndef NO_POLYNOMIAL
				      homoLR =
					pPedigree->likelihood /
					likelihoodQT[pedIdx][gfreqInd][penIdx]
					[paramIdx][thresholdIdx];
#else
				      homoLR =
					pPedigree->likelihood /
					pedigreeSet.nullLikelihood[pedIdx];
#endif
				      hetLR +=
					log10 (alphaV * homoLR + alphaV2);
				    }
				  hetLR = pow (10, hetLR);
				  mp_result[posIdx].het_lr_total += hetLR;
				  if (mp_result[posIdx].
				      max_penIdx < 0
				      || hetLR > mp_result[posIdx].max_lr)
				    {
				      mp_result[posIdx].max_lr = hetLR;
				      mp_result[posIdx].max_alpha = alphaV;
				      mp_result[posIdx].max_gfreq = gfreq;
				      mp_result[posIdx].max_penIdx = penIdx;
				    }
				}
			      /* add the result to the right placeholder */
			      mp_result[posIdx].lr_total += likelihood_ratio;
			      mp_result[posIdx].lr_count++;

			    }	/* end of threshold loop */

			}	/* end of penetrance loop */
		    }		/* end of parameter loop */
		}		/* end of gene freq */
	    }			/* end of QT */

	  /* print out average and log10(max) and maximizing parameters */
	  avgLR =
	    mp_result[posIdx].het_lr_total / (modelRange.nalpha *
					      mp_result[posIdx].lr_count);
	  ppl = (avgLR * avgLR) / (-5.77 + 54 * avgLR + avgLR * avgLR);
	  max = mp_result[posIdx].max_lr;
	  gfreq = mp_result[posIdx].max_gfreq;
	  alphaV = mp_result[posIdx].max_alpha;
	  penIdx = mp_result[posIdx].max_penIdx;
	  fprintf (fpHet, "\t %f  %6.4f %12.8f(%d) %10.6f %f %f ",
		   traitPos, ppl, avgLR,
		   mp_result[posIdx].lr_count, log10 (max), alphaV, gfreq);
	  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
	    {
	      pen_DD = modelRange.penet[liabIdx][0][penIdx];
	      pen_Dd = modelRange.penet[liabIdx][1][penIdx];
	      pen_dd = modelRange.penet[liabIdx][2][penIdx];
	      fprintf (fpHet, " %f %f %f ", pen_DD, pen_Dd, pen_dd);
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
      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
	{
	  fclose (fpLD);
	}
      fclose (fpPPL);
    }
  fclose (fpHet);
  fclose (fpHomo);
  return 0;
}

/* allocate one extra D prime to store average LR for all valid D primes */
int
allocate_tp_result_storage ()
{
  int i, j;

  tp_result =
    (SUMMARY_STAT **) calloc (pLambdaCell->ndprime + 1,
			      sizeof (SUMMARY_STAT *));
  for (i = 0; i < pLambdaCell->ndprime + 1; i++)
    {
      tp_result[i] =
	(SUMMARY_STAT *) calloc (modelRange.ntheta, sizeof (SUMMARY_STAT));
      for (j = 0; j < modelRange.ntheta; j++)
	{
	  tp_result[i][j].max_penIdx = -1;
	}
    }
  return 0;
}

int
free_tp_result_storage (int ndprime)
{
  int i;

  for (i = 0; i < ndprime; i++)
    {
      free (tp_result[i]);
    }
  free (tp_result);
  return 0;
}

double
calculate_R_square (double p1, double q1, double d)
{
  double p2 = 1.0 - p1;
  double q2 = 1.0 - q1;

  return (d * d / (p1 * p2 * q1 * q2));
}

void
free_likelihood_storage (PedigreeSet * pedSet)
{
  int pedIdx;
  int gfreqInd;
  int penIdx;
  int paramIdx;

  if (likelihoodDT == NULL && likelihoodQT == NULL)
    return;

  if (modelType.trait == DT)
    {
      for (pedIdx = 0; pedIdx < pedSet->numPedigree + 1; pedIdx++)
	{
	  for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++)
	    {
	      /* third dimension is penetrance */
	      free (likelihoodDT[pedIdx][gfreqInd]);
	    }
	  /* second dimension is gene freq */
	  free (likelihoodDT[pedIdx]);
	}
      free (likelihoodDT);
    }
  else
    {				/* QT */
      for (pedIdx = 0; pedIdx < pedSet->numPedigree + 1; pedIdx++)
	{
	  for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++)
	    {

	      for (penIdx = 0; penIdx < modelRange.npenet; penIdx++)
		{
		  /* fourth dimension is SD */
		  for (paramIdx = 0; paramIdx < modelRange.nparam; paramIdx++)
		    {
		      /* 5th dimension is threshold */
		      free (likelihoodQT[pedIdx][gfreqInd][penIdx][paramIdx]);

		    }
		  free (likelihoodQT[pedIdx][gfreqInd][penIdx]);
		}
	      /* third dimension is mean */
	      free (likelihoodQT[pedIdx][gfreqInd]);
	    }
	  /* second dimension is gene freq */
	  free (likelihoodQT[pedIdx]);
	}
      free (likelihoodQT);
    }
}


/* Calculate PPL 
 * Input:
 *   thetaCutoff r - (0, r) will have the given weight
 *                   while (r, 0.5) will have 1-weight
 *
 */
double
calculate_PPL (SUMMARY_STAT * result,
	       double thetaCutoff, double weight, double prior)
{
  int i;
  double w1, w2;
  double integral;
  double a;
  double theta1, theta2;
  double avgLR1, avgLR2;
  int numTheta = modelRange.ntheta;
  double PPL;

  integral = 0;
  w1 = weight / thetaCutoff;
  w2 = (1.0 - weight) / (0.5 - thetaCutoff);
  for (i = 0; i < numTheta - 1; i++)
    {
      /* sex averaged theta */
      theta1 = modelRange.theta[0][i];
      theta2 = modelRange.theta[0][i + 1];
      avgLR1 = result[i].het_lr_avg;
      avgLR2 = result[i + 1].het_lr_avg;

      if (theta2 <= thetaCutoff)
	{
	  integral += 0.5 * w1 * (theta2 - theta1) * (avgLR1 + avgLR2);
	}
      else if (theta1 <= thetaCutoff && thetaCutoff < theta2)
	{
	  a =
	    avgLR1 + (avgLR2 - avgLR1) / (theta2 - theta1) * (thetaCutoff -
							      theta1);
	  integral += 0.5 * w1 * (thetaCutoff - theta1) * (avgLR1 + a);
	  integral += 0.5 * w2 * (theta2 - thetaCutoff) * (avgLR2 + a);
	}
      else
	{
	  /* thetaCutoff < theta1 */
	  integral += 0.5 * w2 * (theta2 - theta1) * (avgLR1 + avgLR2);
	}
    }

  PPL = (prior * integral) / ((prior * integral) + (1.0 - prior));

  return PPL;
}

/* This only applies to two point 
 * Calculate the mean of LR at each theta point for all given D primes 
 */
int
get_average_LD_LR (SUMMARY_STAT ** result, int numDPrime, int numTheta)
{
  int i, j;
  double sumLR;
  double avgLR;

  for (i = 0; i < numTheta; i++)
    {
      sumLR = 0;
      for (j = 0; j < numDPrime; j++)
	{
	  sumLR += result[j][i].het_lr_avg;
	}
      avgLR = sumLR / numDPrime;
      result[numDPrime][i].het_lr_avg = avgLR;
    }
  return 0;
}
