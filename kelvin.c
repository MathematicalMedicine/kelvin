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
#include <mcheck.h>

/* Some default global values. */
char markerfile[KMAXFILENAMELEN + 1] = "markers.dat";
char mapfile[KMAXFILENAMELEN + 1] = "mapfile.dat";
char pedfile[KMAXFILENAMELEN + 1] = "pedfile.dat";
char datafile[KMAXFILENAMELEN + 1] = "datafile.dat";
char loopsfile[KMAXFILENAMELEN + 1] = "loops.dat";
char outfile[KMAXFILENAMELEN + 1] = "lods.out";
/* Model datastructures. modelOptions is defined in the pedigree library. */
ModelType modelType;
ModelRange modelRange;
ModelOptions modelOptions;
char sUnknownPedID[] = "0";

/* temporarily for GAW project */
/* summary statistics */
typedef struct SUMMARY_STAT
{
  /* for calculating average */
  double lr_total;		/* sum of log(LR) */
  int lr_count;			/* number of models for the (ld, theta) pair */
  double het_lr_total;		/* with heterogeneity - alpha, parameter */

  /* for max */
  double max_lr;		/* max het lr */
  double max_alpha;
  double max_gfreq;
  int max_penIdx;

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
  Trait *pTrait;
  int pedIdx;
  double homoLR, hetLR;
  double log10HetLR;
  double max;
  double constraint;
  double step_size;
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

  SubLocusList savedLocusList;
  SubLocusList nullLocusList;

  clock_t time0, time1, time2;
  int numberOfCompute = 0;

  memset(&savedLocusList, 0, sizeof(savedLocusList));
  memset(&nullLocusList, 0, sizeof(nullLocusList));

#ifdef DEBUG
  //  mtrace();
#endif

  /* Initialize the logging system. */
  logInit ();
  /* logSet(LOGINPUTFILE, LOGDEBUG); */
  //logSet(LOGGENOELIM, LOGDEBUG);
  /* logSet(LOGPEELGRAPH, LOGFATAL); */
  //logSet (LOGLIKELIHOOD, LOGWARNING);
  /* logSet(LOGLIKELIHOOD, LOGDEBUG); */
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
  modelOptions.sUnknownPersonID = sUnknownPedID;

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

#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      polynomialInitialization ();
      fprintf (stderr,
	       "!!!!!!!!!!!The Computation is done in polynomial mode!!!!!!!!!!!!!!!\n");
    }
#endif
  fprintf (stderr, "%s\n",
	   (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) ? "LE" : "LD");

  /* Read in the map file. */
  read_mapfile (mapfile);

  /* Initialize the locus list and read in the marker file. */
  memset (&originalLocusList, 0, sizeof (originalLocusList));
  /* read in what loci are in the pedigree file */
  read_datafile (datafile);

  /* The configuration has all the information about the disease trait if any */
  if (originalLocusList.numTraitLocus > 0)
    {
      /* Need to add the alleles into trait locus 
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
	  modelRange.nlambda = 1;
	  modelRange.nlamdim = 1;
	  modelRange.lambda = (double **) calloc (1, sizeof (double *));
	  modelRange.lambda[0] = (double *) calloc (1, sizeof (double));
	  modelRange.lambda[0][0] = 0;
	}
      else
	{
	  /* hard coded for GAW project -- They should be specified in kelvin.config 
	   * currently it is specified in datafile.dat
	   * for twopoint analysis only */
	  modelRange.nlambda = 21;
	  modelRange.nlamdim = 1;
	  /* dimension */
	  modelRange.lambda = (double **) calloc (1, sizeof (double *));
	  modelRange.lambda[0] =
	    (double *) calloc (modelRange.nlambda, sizeof (double));
	  step_size = 2.0 / (modelRange.nlambda - 1);
	  /* set up D prime values */
	  for (i = 0; i < modelRange.nlambda; i++)
	    {
	      modelRange.lambda[0][i] = -1.0 + step_size * i;
	    }
	}
      tp_result =
	(SUMMARY_STAT **) calloc (modelRange.nlambda,
				  sizeof (SUMMARY_STAT *));
      for (i = 0; i < modelRange.nlambda; i++)
	{
	  tp_result[i] =
	    (SUMMARY_STAT *) calloc (modelRange.ntheta,
				     sizeof (SUMMARY_STAT));
	  for (j = 0; j < modelRange.ntheta; j++)
	    {
	      tp_result[i][j].max_penIdx = -1;
	    }
	}
    }
  else
    {
      /* assuming we are doing multipoint analysis */
    }

  /* initialize some work space */
  initialize_parental_pair_workspace (&parentalPairSpace,
				      originalLocusList.numLocus);

  /* initialize loci by doing genotype elimination, set recoding */
  initialize_loci (&pedigreeSet);

#ifndef NO_POLYNOMIAL
  /* only for multipoint */
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
#endif

#ifndef NO_POLYNOMIAL
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

  /* after genotype elimination and tracking the max work space needed for constructing parental pair */
  allocate_parental_pair_workspace (&parentalPairSpace,
				    modelType.numMarkers + 1);


  allocate_likelihood_space(&pedigreeSet, modelType.numMarkers+1);

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

  /* Set up grid and invoke direct eval used as in Yungui's test
   * file. If we're doing 2 point, this is pretty straightforward. If
   * we're doing multipoint, this is the nflankloop behavior part. */
  /* assume the trait locus is the first one in the list */
  traitLocus = 0;
  pLocus = originalLocusList.ppLocusList[traitLocus];
  pTraitLocus = originalLocusList.ppLocusList[traitLocus]->pTraitLocus;
  pTrait = pTraitLocus->pTraits[traitLocus];
  if (modelType.type == TP)
    {
      /* Two point. */
      savedLocusList.numLocus = 2;
      savedLocusList.pLocusIndex = (int *) malloc (sizeof (int) *
						   savedLocusList.numLocus);
      savedLocusList.pPrevLocusDistance = (double *) malloc (sizeof (double) *
							     savedLocusList.
							     numLocus);
      savedLocusList.pNextLocusDistance =
	(double *) malloc (sizeof (double) * savedLocusList.numLocus);

      savedLocusList.pLocusIndex[0] = 0;
      savedLocusList.pLocusIndex[1] = 1;
      savedLocusList.pPrevLocusDistance[0] = -1;

      locusList = &savedLocusList;

      if (originalLocusList.pLDLoci == NULL)
	{
	  originalLocusList.pLDLoci = (LDLoci *) malloc (sizeof (LDLoci));
	  pLDLoci = &originalLocusList.pLDLoci[0];
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
      /* Loop over the penetrances, genefrequencies, thetas and call
         the likelihood calculation, storing each value optained to
         disk. */
      if (pTrait->type == DICHOTOMOUS)
	{

	  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
	    pLDLoci = &originalLocusList.pLDLoci[0];

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
		update_penetrance (&pedigreeSet, traitLocus);
#else
	      update_penetrance (&pedigreeSet, traitLocus);
#endif

	      for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++)
		{
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
		  /* get the likelihood at 0.5 first and LD=0 */
		  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
		    {
		      pLDLoci->ppDPrime[0][0] = 0;
		      setup_LD_haplotype_freq (pLDLoci);
		    }
		  locusList->pNextLocusDistance[0] = 0.5;
		  locusList->pPrevLocusDistance[1] = 0.5;
		  compute_likelihood (&pedigreeSet);


		  if (numberOfCompute == 0)
		    {
		      time1 = clock ();
		      numberOfCompute++;
		    }

		  if (pedigreeSet.likelihood == 0.0 &&
		      pedigreeSet.log10Likelihood == -9999.99)
		    {
		      fprintf (stderr, "Theta 0.5 has likelihood 0\n");
		      fprintf (stderr, "dgf=%f\n", gfreq);
		      for (liabIdx = 0; liabIdx < modelRange.nlclass;
			   liabIdx++)
			{
			  pen_DD = modelRange.penet[liabIdx][0][penIdx];
			  pen_Dd = modelRange.penet[liabIdx][1][penIdx];
			  pen_dd = modelRange.penet[liabIdx][2][penIdx];
			  fprintf (stderr,
				   "Liab %d penentrance %f %f %f\n",
				   liabIdx + 1, pen_DD, pen_Dd, pen_dd);
			}

		      exit (-1);
		    }
		  /* save the results for NULL */
		  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++)
		    {
		      /* save the likelihood at null */
		      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
		      pedigreeSet.nullLikelihood[pedIdx] =
			pPedigree->likelihood;
		    }

		  log10_likelihood_null = pedigreeSet.log10Likelihood;
		  for (dprimeIdx = 0; dprimeIdx < modelRange.nlambda;
		       dprimeIdx++)
		    {
		      if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
			{
			  pLDLoci->ppDPrime[0][0] =
			    modelRange.lambda[0][dprimeIdx];
			  setup_LD_haplotype_freq (pLDLoci);
			}
		      for (thetaInd = 0; thetaInd < modelRange.ntheta;
			   thetaInd++)
			{
			  theta = modelRange.theta[0][thetaInd];
			  locusList->pNextLocusDistance[0] = theta;
			  locusList->pPrevLocusDistance[1] = theta;

			  compute_likelihood (&pedigreeSet);

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
				   pedIdx < pedigreeSet.numPedigree; pedIdx++)
				{
				  pPedigree =
				    pedigreeSet.ppPedigreeSet[pedIdx];
				  homoLR =
				    pPedigree->likelihood /
				    pedigreeSet.nullLikelihood[pedIdx];
				  log10HetLR +=
				    log10 (alphaV * homoLR + (1 - alphaV));
				}
			      hetLR = pow (10, log10HetLR);
			      tp_result[dprimeIdx][thetaInd].het_lr_total +=
				hetLR;
			      if (tp_result[dprimeIdx][thetaInd].max_penIdx <
				  0
				  || hetLR >
				  tp_result[dprimeIdx][thetaInd].max_lr)
				{
				  tp_result[dprimeIdx][thetaInd].max_lr =
				    hetLR;
				  tp_result[dprimeIdx][thetaInd].max_alpha =
				    alphaV;
				  tp_result[dprimeIdx][thetaInd].max_gfreq =
				    gfreq;
				  tp_result[dprimeIdx][thetaInd].max_penIdx =
				    penIdx;
				}
			    }	/* end of calculating HET LR */
			  /* add the result to the right placeholder */
			  tp_result[dprimeIdx][thetaInd].lr_total +=
			    likelihood_ratio;
			  tp_result[dprimeIdx][thetaInd].lr_count++;
			  //fprintf(stderr, "likelihood ratio: %e.\n", likelihood_ratio);

			}	/* end of theta loop */
		    }		/* end of D prime loop */
		}		/* end of genFreq loop */
	    }			/* end of penetrance loop */
	}			/* end of TP */
      else
	/* should be QT or COMBINED - twopoint */
	{
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
			   thresholdIdx < modelRange.ntthresh; thresholdIdx++)
			{
			  for (liabIdx = 0; liabIdx < modelRange.nlclass;
			       liabIdx++)
			    {
			      mean_DD = modelRange.penet[liabIdx][0][penIdx];
			      mean_Dd = modelRange.penet[liabIdx][1][penIdx];
			      mean_dd = modelRange.penet[liabIdx][2][penIdx];
			      SD_DD =
				modelRange.param[liabIdx][0][0][paramIdx];
			      SD_Dd =
				modelRange.param[liabIdx][1][0][paramIdx];
			      SD_dd =
				modelRange.param[liabIdx][2][0][paramIdx];
			      /* threshold for QT */
			      threshold =
				modelRange.tthresh[liabIdx][thresholdIdx];


			      /* check against the hard coded constraint */
			      constraint =
				pow (1.0 - gfreq,
				     2) * mean_dd * SD_dd + 2 * gfreq * (1 -
									 gfreq)
				* mean_Dd * SD_Dd + pow (gfreq,
							 2) * mean_DD * SD_DD;
/*	  fprintf(stderr, "constraint: %f gfreq:%f DD (%f,%f) Dd(%f,%f) dd(%f,%f)\n",
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

			  /* get the likelihood at 0.5 first and LD=0 */
			  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
			    {
			      pLDLoci->ppDPrime[0][0] = 0;
			      setup_LD_haplotype_freq (pLDLoci);
			    }

			  locusList->pNextLocusDistance[0] = 0.5;
			  locusList->pPrevLocusDistance[1] = 0.5;
			  compute_likelihood (&pedigreeSet);


			  if (numberOfCompute == 0)
			    {
			      time1 = clock ();
			      numberOfCompute++;
			    }

			  if (pedigreeSet.likelihood == 0.0 &&
			      pedigreeSet.log10Likelihood == -9999.99)
			    {
			      fprintf (stderr,
				       "Theta 0.5 has likelihood 0\n");
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

			  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree;
			       pedIdx++)
			    {
			      /* save the likelihood at null */
			      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
			      pedigreeSet.nullLikelihood[pedIdx] =
				pPedigree->likelihood;
			    }

			  log10_likelihood_null = pedigreeSet.log10Likelihood;
			  for (dprimeIdx = 0; dprimeIdx < modelRange.nlambda;
			       dprimeIdx++)
			    {
			      if (modelOptions.equilibrium !=
				  LINKAGE_EQUILIBRIUM)
				{
				  pLDLoci->ppDPrime[0][0] =
				    modelRange.lambda[0][dprimeIdx];
				  setup_LD_haplotype_freq (pLDLoci);
				}
			      for (thetaInd = 0; thetaInd < modelRange.ntheta;
				   thetaInd++)
				{
				  theta = modelRange.theta[0][thetaInd];
				  locusList->pNextLocusDistance[0] = theta;
				  locusList->pPrevLocusDistance[1] = theta;
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
				      hetLR = 0;
				      for (pedIdx = 0;
					   pedIdx < pedigreeSet.numPedigree;
					   pedIdx++)
					{
					  pPedigree =
					    pedigreeSet.ppPedigreeSet[pedIdx];
					  homoLR =
					    pPedigree->likelihood /
					    pedigreeSet.
					    nullLikelihood[pedIdx];
					  hetLR +=
					    log10 (alphaV * homoLR + alphaV2);
					}
				      hetLR = pow (10, hetLR);
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
					}
				    }
				  /* add the result to the right placeholder */
				  tp_result[dprimeIdx][thetaInd].lr_total +=
				    likelihood_ratio;
				  tp_result[dprimeIdx][thetaInd].lr_count++;
				  //fprintf(stderr, "likelihood ratio: %e.\n", likelihood_ratio);

				}	/* end of theta */
			    }	/* end of D prime */
			}	/* end of threshold loop */
		    }		/* end of penetrance loop */
		}		/* end of parameter loop */
	    }			/* end of gene freq */
	}			/* end of QT */

      /* for GAW project only */
      fprintf (stderr, "#0  \"Average Homo LR\" \n");
      for (dprimeIdx = 0; dprimeIdx < modelRange.nlambda; dprimeIdx++)
	{
	  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
	    {
	      theta = modelRange.theta[0][thetaInd];
	      if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
		{
		  fprintf (stderr, "\t (%f,%f)  %f(%d)\n",
			   theta, theta,
			   tp_result[dprimeIdx][thetaInd].lr_total /
			   tp_result[dprimeIdx][thetaInd].lr_count,
			   tp_result[dprimeIdx][thetaInd].lr_count);
		}
	      else
		{
		  fprintf (stderr, "%f %d %f %f\n",
			   modelRange.lambda[0][dprimeIdx],
			   tp_result[dprimeIdx][thetaInd].lr_count,
			   theta,
			   log10 (tp_result[dprimeIdx][thetaInd].lr_total /
				  tp_result[dprimeIdx][thetaInd].lr_count));
		}
	    }
	}
      fprintf (stderr, "-	Total 1234(1234)\n");

      fprintf (stderr, "#1 \"Average Het LR\" \n");
      for (dprimeIdx = 0; dprimeIdx < modelRange.nlambda; dprimeIdx++)
	{
	  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
	    {
	      theta = modelRange.theta[0][thetaInd];
	      if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
		{
		  fprintf (stderr, "\t (%f,%f)  %12.8f(%d)\n",
			   theta, theta,
			   tp_result[dprimeIdx][thetaInd].het_lr_total /
			   (modelRange.nalpha *
			    tp_result[dprimeIdx][thetaInd].lr_count),
			   modelRange.nalpha *
			   tp_result[dprimeIdx][thetaInd].lr_count);
		}
	      else
		{
		  fprintf (stderr, "%f %d %f %12.8f\n",
			   modelRange.lambda[0][dprimeIdx],
			   modelRange.nalpha *
			   tp_result[dprimeIdx][thetaInd].lr_count, theta,
			   log10 (tp_result[dprimeIdx][thetaInd].
				  het_lr_total / (modelRange.nalpha *
						  tp_result[dprimeIdx]
						  [thetaInd].lr_count)));
		}
	    }

	}
      fprintf (stderr, "-	Total 1234(1234)\n");

      fprintf (stderr, "#2 Max Het LR \n");
      max = -99999;
      for (dprimeIdx = 0; dprimeIdx < modelRange.nlambda; dprimeIdx++)
	{
	  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
	    {
	      theta = modelRange.theta[0][thetaInd];
	      if (tp_result[dprimeIdx][thetaInd].max_lr > max)
		{
		  max = tp_result[dprimeIdx][thetaInd].max_lr;
		}
	    }
	}
      fprintf (stderr, "MOD: %10.6f\n", log10 (max));

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
	}

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
      fprintf (stderr, "#1 \"Average Het LR\" \n");
      fprintf (stderr,
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
			      //fprintf(stderr, "likelihood ratio: %e.\n", likelihood_ratio);

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
	  fprintf (stderr, "\t %f  %6.4f %12.8f(%d) %10.6f %f %f ",
		   traitPos, ppl, avgLR,
		   mp_result[posIdx].lr_count, log10 (max), alphaV, gfreq);
	  for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
	    {
	      pen_DD = modelRange.penet[liabIdx][0][penIdx];
	      pen_Dd = modelRange.penet[liabIdx][1][penIdx];
	      pen_dd = modelRange.penet[liabIdx][2][penIdx];
	      fprintf (stderr, " %f %f %f ", pen_DD, pen_Dd, pen_dd);
	    }
	  /* print out markers used for this position */
	  fprintf (stderr, "(%d", mp_result[posIdx].pMarkers[0]);
	  for (k = 1; k < modelType.numMarkers; k++)
	    {
	      fprintf (stderr, ",%d", mp_result[posIdx].pMarkers[k]);
	    }
	  fprintf (stderr, ")\n");
	  fflush (stderr);
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

  if(modelType.type == TP)
    {
      /* two point */
      for (i = 0; i < modelRange.nlambda; i++)
	{
	  free(tp_result[i]);
	}
      free(tp_result);
    }
  else 
    {
      /* multipoint */
      for (posIdx = 0; posIdx < numPositions; posIdx++)
	{
	  free(mp_result[posIdx].pMarkers);
	}
      free(mp_result);
    }

  /* free parental pair work space */
  free_parental_pair_workspace (&parentalPairSpace, modelType.numMarkers + 1);

#ifndef NO_POLYNOMIAL
		  if (modelOptions.polynomial == TRUE)
		    {
		      pedigreeSetPolynomialClearance (&pedigreeSet);
	              polynomialClearance();
                    }
#endif
  free_likelihood_space(&pedigreeSet);
  free_pedigree_set(&pedigreeSet);
  free_sub_locus_list(&nullLocusList);
  free_sub_locus_list(&savedLocusList);
  final_cleanup();
  

  fprintf (stderr, "Computation time:  %fs  %fs \n",
	   (double) (time1 - time0) / CLOCKS_PER_SEC,
	   (double) (time2 - time1) / CLOCKS_PER_SEC);


  exit (TRUE);
  return 0;
}
