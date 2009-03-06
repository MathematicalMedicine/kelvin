  overallSW = swCreate ("overall");
  combinedComputeSW = swCreate ("combinedComputeSW");
  combinedBuildSW = swCreate ("combinedBuildSW");

  /* Setup all of our signal handlers BEFORE we start any threads. */
  setupHandlers ();

  /* Start a thread with a timer to do the memory checks. It can afford
     to hang, while the main process cannot. */
  pthread_t statusThread;
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
#ifdef __OPTIMIZE__
  swLogMsg ("GCC optimization enabled");
#else
  swLogMsg ("GCC optimization disabled");
#endif
  fprintf (stdout, "To check status (at some risk), type CTRL-\\ or type \"kill -%d %d\".\n", SIGQUIT, (int) getpid ());

  // THESE SHOULD BE SOMEWHERE ELSE

  memset (&savedLocusList, 0, sizeof (savedLocusList));
  memset (&markerLocusList, 0, sizeof (markerLocusList));
  memset (&traitLocusList, 0, sizeof (traitLocusList));

  memset (&modelType, 0, sizeof (modelType));

  // Initialize the logging system.
  logInit ();

  /* Start by parsing command line arguments. Most essential: figure
   * out where the configuration file lives. */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-')
      switch (argv[i][1]) {
      case '?':
        /* Help */
        fprintf (stdout, "Usage:\n");
        fprintf (stdout, "  %s [-?] <configuration file>\nwhere:\n", argv[0]);
        fprintf (stdout, "      -? : this output;\n");
        fprintf (stdout, "      <configuration file> : file containing run parameters.\n");
        exit (EXIT_FAILURE);
        break;
    } else if (strlen (configfile) != 0) {
      /* Unexpected argument; we already have a configuration file! Punt. */
      KLOG (LOGDEFAULT, LOGFATAL, "Unexpected command line argument '%s'; aborting.\n", argv[i]);
    } else if (strlen (argv[i]) >= KMAXFILENAMELEN) {
      /* Configuration file name too long! Punt. */
      KLOG (LOGDEFAULT, LOGFATAL, "Config file name '%s' too long (>=%d); aborting.\n", argv[i], KMAXFILENAMELEN);
    } else {
      /* Got a configuration file name. Copy it. */
      strncpy (configfile, argv[i], KMAXFILENAMELEN);
      char currentWorkingDirectory[MAXSWMSG - 32];
      getcwd (currentWorkingDirectory, sizeof (currentWorkingDirectory));
      sprintf (messageBuffer, "In %s w/%s", currentWorkingDirectory, configfile);
      swLogMsg (messageBuffer);
    }
    i++;
  }
  /* Check to see if the configuration file name was specified. */
  KASSERT ((strlen (configfile) > 0), "No configuration file specified; aborting.\n");

  /* Parse the configuration file. */
  KASSERT (readConfigFile (configfile)
           != ERROR, "Error in configuration file; aborting.\n");

  /* For now, reject all models we can't deal with. */
  KASSERT (modelRange.nalleles == 2, "Only biallelic traits supported.\n");

  /* The difference between QT and CT is whether we use threshold or not. Under CT there must 
   * be thresholds, under QT there should not. */
  if (modelRange.ntthresh > 0 && modelType.trait != DT) {
    modelType.trait = CT;
    KASSERT (modelType.minThreshold > -999999998 &&
             modelType.maxThreshold < 999999998,
             "Under QT threshold model, MIN and MAX of the QT threshold values need to be "
	     "provided through keywords T_MIN and T_MAX.\n");
  }

  if (modelOptions.polynomial == TRUE) {
    swLogMsg ("Computation is done in polynomial mode");
#ifdef POLYUSE_DL
    swLogMsg ("Dynamic libraries for polynomial evaluation will be used if found");
#endif
    polynomialInitialization ();
  } else {
    swLogMsg ("Computation is done in non-polynomial (direct evaluation) mode");
  }
  if (modelOptions.integration == TRUE) {
    swLogMsg ("Integration is done numerically (dkelvin)");
  } else {
    swLogMsg ("Integration is done with iteration (original kelvin)");
  }

  /* Read in the map file. */
  read_mapfile (mapfile);

  /* Initialize the locus list and read in the marker file. */
  memset (&originalLocusList, 0, sizeof (originalLocusList));
  /* read in what loci are in the pedigree file */
  read_datafile (datafile);

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
  read_markerfile (markerfile, modelType.numMarkers);

  /* build allele set information */
  for (locus = 0; locus < originalLocusList.numLocus; locus++) {
    construct_original_allele_set_list (locus);
  }

  /* Initialize the pedigree set datastructure and read in the pedigrees. */
  memset (&pedigreeSet, 0, sizeof (PedigreeSet));
  read_pedfile (pedfile, &pedigreeSet);
  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
    pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
    if (pPedigree->currentLoopFlag)
      exitDueToLoop = TRUE;
  }
  KASSERT (exitDueToLoop == FALSE, "Not all loops in pedigrees are broken.\n");

  /* read in case control file if provided */
  if (strlen (ccfile) > 0)
    read_ccfile (ccfile, &pedigreeSet);
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
    fpCond = fopen (condFile, "w");
    KASSERT (fpCond != NULL, "Error in opening file %s for write.\n", condFile); 
    //  fprintf( fpCond, "# Version %s\n", programVersion);
   }

  if (modelOptions.markerAnalysis == FALSE) {
    fpHet = fopen (avghetfile, "w");
    KASSERT (fpHet != NULL, "Error in opening file %s for write.\n", avghetfile);
    fprintf (fpHet, "# Version %s\n", programVersion);
  }

  fpMOD = fopen (modfile, "w");
  KASSERT (fpMOD != NULL, "Error in opening file %s for write.\n", modfile);
  fprintf (fpMOD, "# Version %s\n", programVersion);

  if (modelType.type == TP) {
    fpPPL = fopen (pplfile, "w");
    KASSERT (fpPPL != NULL, "Error in opening file %s for write.\n", pplfile);
    writePPLFileHeader ();

    if (strlen (maxmodelfile) > 0) {
      fpTP = fopen (maxmodelfile, "w");
      KASSERT (fpTP != NULL, "Error in opening file %s for write.\n", maxmodelfile);
    }
  }

  if (strlen (intermediatefile) > 0) {
    fpIR = fopen (intermediatefile, "w");
    KASSERT (fpIR != NULL, "Error in opening file %s for write.\n", intermediatefile);
  }

 // DKelvin intermediate results are written here.
  if ((modelOptions.integration) && (strlen (dkelvinoutfile) > 0)) {
    fpDK= fopen(dkelvinoutfile, "w");
    KASSERT (fpDK != NULL, "Error in opening file %s for write.\n", dkelvinoutfile);
  }
