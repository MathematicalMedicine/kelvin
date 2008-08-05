
  overallSW = swCreate ("overall");
  combinedComputeSW = swCreate ("combinedComputeSW");
  combinedBuildSW = swCreate ("combinedBuildSW");

  /* Setup all of our signal handlers BEFORE we start any threads. */
  setupHandlers ();

  /* Start a thread with a timer to do the memory checks. It can afford
     to hang, while the main process cannot. */
  pthread_t statusThread;
  if (pthread_create ( &statusThread, NULL, monitorStatus, NULL)) {
    perror ("Failed to create status monitoring thread");
    exit (EXIT_FAILURE);
  }

  /* Annouce ourselves for performance tracking. */
  char currentWorkingDirectory[MAXSWMSG - 32];

  sprintf (messageBuffer, "kelvin %s built %s %s", programVersion, __DATE__, __TIME__);
  swLogMsg (messageBuffer);
  swLogMsg (kelvinVersion);
  swLogMsg (likelihoodVersion);
  swLogMsg (locusVersion);
  swLogMsg (polynomialVersion);
  sprintf (messageBuffer, "Compiler %s, GSL %s\n", __VERSION__, GSL_VERSION);
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
  sprintf (messageBuffer, "GNU profiler (gprof) run, use \"kill -%d %d\" to finish early.", SIGTERM, getpid ());
  swLogMsg (messageBuffer);
#endif
#ifdef GCOV
  sprintf (messageBuffer, "GNU coverage analyzer (gcov) run, use \"kill -%d %d\" to finish early.", SIGTERM, getpid ());
  swLogMsg (messageBuffer);
#endif
  fprintf (stdout, "To check status (at some risk), type CTRL-\\ or type \"kill -%d %d\".\n", SIGQUIT, getpid ());
  swStart (overallSW);

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
      getcwd (currentWorkingDirectory, sizeof (currentWorkingDirectory));
      sprintf (messageBuffer, "In %s w/%s", currentWorkingDirectory, configfile);
      swLogMsg (messageBuffer);
    }
    i++;
  }
  /* Check to see if the configuration file name was specified. */
  KASSERT ((strlen (configfile) > 0), "No configuration file specified; aborting.\n");

  // IMHO, ALL OF THE FOLLOWING UP TO OPENING FILES BELONGS IN CONFIG.C

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

  /* set default for QT */
  modelType.minOriginal = -999999999.00;
  modelType.maxOriginal = 999999999.00;
  modelType.minThreshold = -999999999.00;
  modelType.maxThreshold = 999999999.00;

  /* Parse the configuration file. */
  KASSERT (readConfigFile (configfile, &modelType, &modelRange, &modelOptions)
           != ERROR, "Error in configuration file; aborting.\n");

  /* For now, reject all models we can't deal with. */
  KASSERT (modelRange.nalleles == 2, "Only biallelic traits supported.\n");

  /* the difference between QT and CT is whether we use threshold or not. Under CT -  yes to
   * threshold, under QT - no threshold */
  if (modelRange.ntthresh > 0 && modelType.trait != DT) {
    modelType.trait = CT;
    KASSERT (modelType.minThreshold > -999999998 &&
             modelType.maxThreshold < 999999998,
             "Under QT threshold model, MIN and MAX of the QT threshold values need to be "
	     "provided through keywords T_MIN and T_MAX.\n");
  }

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
