#include "kelvinWriteFiles.h"
pthread_t statusThread;
int exitDueToLoop = FALSE; /* exit due to unbroken loop */
int k;
char messageBuffer[MAXSWMSG];

  overallSW = swCreate ("overall");
  combinedComputeSW = swCreate ("combinedComputeSW");
  combinedBuildSW = swCreate ("combinedBuildSW");

  /* Setup all of our signal handlers BEFORE we start any threads. */
  setupHandlers ();

  /* Start a thread with a timer to do the memory checks. It can afford
     to hang, while the main process cannot. */
  if (pthread_create ( &statusThread, NULL, monitorStatus, NULL))
    perror ("Failed to create status monitoring thread, no progress status will be displayed or written");

  /* Annouce ourselves for performance tracking. */

  pushStatus ('k', "NonSpecific");
  sprintf (messageBuffer, "kelvin %s built %s %s", programVersion, __DATE__, __TIME__);
  swLogMsg (messageBuffer);
  swLogMsg (kelvinVersion);
  swLogMsg (likelihoodVersion);
  swLogMsg (locusVersion);
  swLogMsg (polynomialVersion);
  sprintf (messageBuffer, "Compiler %s\n", __VERSION__);
  swLogMsg (messageBuffer);

#ifdef _OPENMP
  if ((envVar = getenv ("OMP_NUM_THREADS")) != NULL)
    threadCount = atoi (envVar);
  sprintf (messageBuffer, "OpenMP-enabled w/%d threads.", threadCount);
  swLogMsg (messageBuffer);
#endif
#ifdef FAKEEVALUATE
  swLogMsg ("Polynomial evaluation is being SKIPPED FOR TESTING, results will be wrong!");
#endif
#ifdef DMUSE
  swLogMsg ("Experimental static memory handling enabled!");
#endif
#ifdef DMTRACK
  swLogMsg ("Dynamic memory usage dumping is turned on, so performance will be poor!");
#endif
#ifdef GPROF
  sprintf (messageBuffer, "GNU profiler (gprof) run, use \"kill -%d %d\" to finish early.", SIGTERM, (int) getpid ());
  swLogMsg (messageBuffer);
#endif
#ifdef GCOV
  sprintf (messageBuffer, "GNU coverage analyzer (gcov) run, use \"kill -%d %d\" to finish early.", SIGTERM, (int) getpid ());
  swLogMsg (messageBuffer);
#endif
#ifdef USE_GSL
swLogMsg ("Using GNU Scientific Library (GSL) statistical functions instead of internal ones!");
#endif
  
  swStart (overallSW);
#ifdef GCCOPT
  sprintf (messageBuffer, "GCC optimization level %d enabled", GCCOPT);
  swLogMsg (messageBuffer);
#else
  swLogMsg ("GCC optimization disabled (or GCCOPT not defined)");
#endif
  fprintf (stdout, "To check status (at some risk), type CTRL-\\ or type \"kill -%d %d\".\n", SIGQUIT, (int) getpid ());

  // THESE SHOULD BE SOMEWHERE ELSE

  memset (&savedLocusList, 0, sizeof (savedLocusList));
  memset (&markerLocusList, 0, sizeof (markerLocusList));
  memset (&traitLocusList, 0, sizeof (traitLocusList));

  memset (&modelType, 0, sizeof (modelType));

  // Initialize the logging system.
  logInit ();

  /* Tolerate a request for help as the first argument; ignore the rest of the command 
   * line if help is requested. Otherwise, argv[1] better be the name of the configuration
   * file; the rest of the command line will be treated as override directives.
   */
  if (argc > 1) {
    if (! (strcmp (argv[1], "--help") && strcmp (argv[1], "-?"))) {
      fprintf (stderr, "usage: %s <conffile> [--directive arg1 arg2... [--directive...]]\n",
	       argv[0]);
      exit (0);
    } else
      strcpy (configfile, argv[1]);
  } else {
    fprintf (stderr, "usage: %s <conffile> [--directive arg1 arg2... [--directive...]]\n",
	     argv[0]);
    exit (-1);
  }

  /* Set modelRange, modelOptions and modelType to default values */
  initializeDefaults ();

  /* Parse the configuration file. */
  readConfigFile (argv[1]);

  /* If there's anything on the command line after the configuration file name, 
   * it must be override directives.
   */
  if (argc > 2) {
    parseCommandLine (argc-2, &argv[2]);
  }

  /* Make sure the config as read from the configuration file, and possibly modified on 
   * the command line, is legal. Then clean up the bits of memory allocated during
   * config parsing.
   */
  validateConfig ();
  finishConfig ();

  /* For now, reject all models we can't deal with. */
  KASSERT (modelRange.nalleles == 2, "Only biallelic traits supported.\n");

  if (modelOptions.polynomial == TRUE) {
    swLogMsg ("Computation is done in polynomial mode");
#ifdef POLYUSE_DL
    swLogMsg ("Dynamic libraries for polynomial evaluation will be used if found");
#endif
    polynomialInitialization (modelOptions.polynomialScale);
  } else {
    swLogMsg ("Computation is done in non-polynomial (direct evaluation) mode");
  }
  if (modelOptions.integration == TRUE) {
    swLogMsg ("Integration is done numerically (dkelvin)");
  } else {
    swLogMsg ("Integration is done with iteration (original kelvin)");
  }

  /* Read in the map file. */
  read_mapfile (modelOptions.mapfile);

  /* Initialize the locus list and read in the marker file. */
  memset (&originalLocusList, 0, sizeof (originalLocusList));
  /* read in what loci are in the pedigree file */
  read_datafile (modelOptions.datafile);

  /* The configuration has all the information about the disease trait if any */
  if (originalLocusList.numTraitLocus > 0) {
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
    if (modelType.trait == QT || modelType.trait == CT) {
      modelType.min = (modelType.minOriginal - modelType.mean) / modelType.sd;
      modelType.max = (modelType.maxOriginal - modelType.mean) / modelType.sd;
      pTrait->minFlag = modelType.minFlag;
      pTrait->maxFlag = modelType.maxFlag;
      pTrait->min = modelType.min;
      pTrait->max = modelType.max;
      pTrait->functionQT = modelType.distrib;
      if (modelType.distrib == QT_FUNCTION_T)
        pTrait->dfQT = modelType.constants[0];
      pTrait->sampleMean = modelType.mean;
      pTrait->sampleSD = modelType.sd;
      pTrait->unknownTraitValue = modelOptions.affectionStatus[AFFECTION_STATUS_UNKNOWN];
      pTrait->lessCutoffFlag = modelOptions.affectionStatus[AFFECTION_STATUS_UNAFFECTED];
      pTrait->moreCutoffFlag = modelOptions.affectionStatus[AFFECTION_STATUS_AFFECTED];
    }
  }

  /* read in marker allele frequencies */
  read_markerfile (modelOptions.markerfile, modelType.numMarkers);

  /* build allele set information */
  for (locus = 0; locus < originalLocusList.numLocus; locus++) {
    construct_original_allele_set_list (locus);
  }

  /* Initialize the pedigree set datastructure and read in the pedigrees. */
  memset (&pedigreeSet, 0, sizeof (PedigreeSet));
  read_pedfile (modelOptions.pedfile, &pedigreeSet);
  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
    if (pPedigree->currentLoopFlag)
      exitDueToLoop = TRUE;
  }
  KASSERT (exitDueToLoop == FALSE, "Not all loops in pedigrees are broken.\n");

  /* read in case control file if provided */
  if (strlen (modelOptions.ccfile) > 0)
    read_ccfile (modelOptions.ccfile, &pedigreeSet);
  flexBufferSize = 0;
  free (flexBuffer);
  fflush (stderr);
  fflush (stdout);

  if (modelType.trait == QT) {
    /* threshold value will not be used in any meaningful way, but we will use it for 
       the loop */
    modelRange.ntthresh = 1;
    modelType.minOriginal = 0;
    modelType.maxOriginal = 1;
    if (modelRange.tthresh == NULL) {
      modelRange.tthresh = (double **) malloc (sizeof (double *));
      for (i = 0; i < modelRange.nlclass; i++) {
	modelRange.tthresh[i] = malloc (sizeof (double));
      }
    }
  }

  /* allocate space for results */
  if (modelType.type == TP) {
    modelType.numMarkers = 1;
    totalLoci = 2;
    /* two point analysis */
    if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) {
      /* in order to simplify looping, even for LE, we add a fake LD parameter dprime=0, which
       * is LE */
      modelRange.ndprime = 1;
      modelRange.dprime = (double *) calloc (1, sizeof (double));
      modelRange.dprime[0] = 0;
      pLambdaCell = findLambdas (&modelRange, 2, 2);
      dprime0Idx = 0;
    }
  } else {
    /* we are doing multipoint analysis */
    totalLoci = modelType.numMarkers + originalLocusList.numTraitLocus;
    if (modelRange.tlocRangeStart >= 0) {
      double endofmap, tloc;
      endofmap = map.ppMapUnitList[map.count-1]->mapPos[SEX_AVERAGED] +
	modelRange.tlocRangeIncr;
      i = 0;
      while ((tloc = modelRange.tlocRangeStart + (i * modelRange.tlocRangeIncr)) <= endofmap) {
	addTraitLocus (&modelRange, tloc);
	i++;
      }
    }
    if (modelRange.tlmark == TRUE) {
      /* add marker positions to the list of positions we want to conduct analysis */
      for (i = 0; i < originalLocusList.numLocus; i++) {
        pLocus = originalLocusList.ppLocusList[i];
        if (pLocus->locusType == LOCUS_TYPE_TRAIT)
          continue;
        addTraitLocus (&modelRange, pLocus->pMapUnit->mapPos[SEX_AVERAGED]);
      }
    }
  }

  /* Estimate number of calls to each (appropriate) instance of compute_likelihood for
   * use in progress reporting, and display model information at this point since markers have
   * already been added to locus list */
  swLogMsg (estimateIterations (eCL));

  /* allocate storage for keeping track of het locus in nuclear families */
  allocate_nucfam_het (&pedigreeSet, totalLoci);

  /* initialize some work space */
  initialize_parental_pair_workspace (&parentalPairSpace, originalLocusList.numLocus);
  /* allocate transmission probability matrix */
  build_xmission_matrix (&nullMatrix, totalLoci);
  build_xmission_matrix (&altMatrix, totalLoci);
  build_xmission_matrix (&traitMatrix, 1);
  build_xmission_matrix (&markerMatrix, totalLoci - 1);
  xmissionMatrix = nullMatrix;
  tmpID = (char *) calloc (totalLoci, sizeof (char));

  /* initialize loci by doing genotype elimination, set recoding */
  initialize_loci (&pedigreeSet);

  if (modelOptions.dryRun != 0) {
    for (loc1 = 0; loc1 < originalLocusList.numLocus; loc1++) {
      fprintf (stderr, "Locus %d:\n", loc1);
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
        pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
        print_pedigree_locus_genotype_count (pPedigree, loc1);
      }
    }
  }
  if (modelOptions.polynomial == TRUE) {
    for (k = 0; k < 3; k++) {
      initialProbPoly[k] = constant1Poly;
      initialProbPoly2[k] = constant1Poly;
      initialProbAddr[k] = initialProbPoly[k];
      initialProbAddr2[k] = initialProbPoly2[k];
      initialHetProbAddr[k] = NULL;
    }
  } else {
    for (k = 0; k < 3; k++) {
      initialProb[k] = 1.0;
      initialProb2[k] = 1.0;
      initialProbAddr[k] = &initialProb[k];
      initialProbAddr2[k] = &initialProb2[k];
      initialHetProbAddr[k] = NULL;
    }
  }

  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
    pPedigree->load_flag = 0;   /* Initially 0 and changes to 1 when marker or 
                                 * alternative likelihood values are retrieved */
  }

  /* Open output files that get written across loops. */

  if(modelOptions.conditionalRun == 1 || modelOptions.loopCondRun == 1) {
    fpCond = fopen (modelOptions.condFile, "w");
    KASSERT (fpCond != NULL, "Error in opening file %s for write.\n", modelOptions.condFile); 
    //  fprintf( fpCond, "# Version %s\n", programVersion);
   }

  if (modelOptions.markerAnalysis == FALSE) {
    fpHet = fopen (modelOptions.avghetfile, "w");
    KASSERT (fpHet != NULL, "Error in opening file %s for write.\n", modelOptions.avghetfile);
    fprintf (fpHet, "# Version %s\n", programVersion);
  }

  if (modelType.type == TP) {
    fpPPL = fopen (modelOptions.pplfile, "w");
    KASSERT (fpPPL != NULL, "Error in opening file %s for write.\n", modelOptions.pplfile);
    writePPLFileHeader ();

    /*
    if (strlen (modelOptions.maxmodelfile) > 0) {
      fpTP = fopen (modelOptions.maxmodelfile, "w");
      KASSERT (fpTP != NULL, "Error in opening file %s for write.\n", modelOptions.maxmodelfile);
    }
    */
  }

  if (strlen (modelOptions.modfile) > 0) {
    fpMOD = fopen (modelOptions.modfile, "w");
    KASSERT (fpMOD != NULL, "Error in opening file %s for write.\n", modelOptions.modfile);
    fprintf (fpMOD, "# Version %s\n", programVersion);
  }

  if (strlen (modelOptions.intermediatefile) > 0) {
    fpIR = fopen (modelOptions.intermediatefile, "w");
    KASSERT (fpIR != NULL, "Error in opening file %s for write.\n", modelOptions.intermediatefile);
  }

 // DKelvin intermediate results are written here.
  if ((modelOptions.integration) && (strlen (modelOptions.dkelvinoutfile) > 0)) {
    fpDK= fopen(modelOptions.dkelvinoutfile, "w");
    KASSERT (fpDK != NULL, "Error in opening file %s for write.\n", modelOptions.dkelvinoutfile);
  }

