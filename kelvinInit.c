#include <limits.h>     // For things like PATH_MAX.
#include <float.h>      // Limits for floating point

#include <pthread.h>    // For memory-use tracking
#include <ctype.h>      // For toupper

#ifdef _OPENMP
#include <omp.h>
#endif

#include "kelvin.h"
#include "kelvinGlobals.h"
#include "kelvinHandlers.h"
#include "config/config.h"
#include "trackProgress.h"
#include "kelvinWriteFiles.h"
#include "utils/pageManagement.h"

#ifdef STUDYDB
#include "database/StudyDB.h"
extern struct StudyDB studyDB;
#include "database/databaseSupport.h"
#endif

struct swStopwatch *combinedComputeSW,  ///< Combined likelihood compute stopwatch
  *combinedBuildSW,      ///< Combined likelihood polynomial build stopwatch
  *overallSW,    ///< Overall stopwatch for the entire run.
  *singleModelSW;


char configfile[PATH_MAX];      ///< Configuration file read to populate all of this

int dprime0Idx = 0;

extern char *likelihoodVersion, *locusVersion, *polynomialVersion;
extern Polynomial *constant1Poly;

void kelvinInit (int argc, char *argv[])
{
  pthread_t monitorMemoryThread;
  int i, k;

  swDiagInit ();

  overallSW = swCreate ("overall");
  combinedComputeSW = swCreate ("combinedComputeSW");
  combinedBuildSW = swCreate ("combinedBuildSW");
  singleModelSW = swCreate ("singleModelSW");

  /* Setup all of our signal handlers. */
  setupHandlers ();

  /* Start a thread with a timer to do the memory checks. It can afford                                                 
   * to hang, while the main process cannot. */
  if (pthread_create (&monitorMemoryThread, NULL, monitorMemory, NULL))
    perror ("Failed to create memory monitoring thread, memory info will be displayed or written");

  /* Annouce ourselves for performance tracking. */

  swPushPhase ('k', "NonSpecific");
  INFO ("kelvin %s edit %s built %s %s on %s", programVersion, svnVersion, __DATE__, __TIME__, getenv("HOSTNAME"));
  INFO ("References for this version of kelvin:\n\n\t(1) Vieland VJ. Thermometers: Something for statistical geneticists to think\n\tabout. Human Hered, 61:144-156, 2006.\n\t(2) Huang Y, Segre A, O'Connell J, Valentine-Cooper W, Seok SC, Vieland VJ.\n\tKelvin:  A 2nd generation software package for computation of the PPL\n\tframework. [Abstract Program number 2336]. Presented at the annual meeting\n\tof The American Society of Human Genetics, November 2008, Philadelphia,\n\tPennsylvania. Available at http://www.ashg.org/2008meeting/abstracts/fulltext/\n");
  //  INFO ("%s", kelvinVersion);
  //  INFO ("%s", likelihoodVersion);
  //  INFO ("%s", locusVersion);
  //  INFO ("%s", polynomialVersion);
  INFO ("Compiler %s", __VERSION__);

#ifdef FAKEEVALUATE
  INFO ("Polynomial evaluation is being SKIPPED FOR TESTING, results will be wrong!");
#endif
#ifdef DMUSE
  INFO ("Experimental static memory handling enabled!");
#endif
#ifdef DMTRACK
  INFO ("Dynamic memory usage dumping is turned on, so performance will be poor!");
#endif

#ifdef USE_GSL
  INFO ("Using GNU Scientific Library (GSL) statistical functions instead of internal ones");
  #ifdef VERIFY_GSL
    #ifdef _OPENMP
      #undef _OPENMP
      #warning "Cannot use OpenMP when using internal statistical functions.");
      INFO ("OpenMP is DISABLED when using internal statistical functions.");
    #endif
  #else
    #ifdef _OPENMP
      INFO ("OpenMP-enabled w/maximum of %d thread(s).", omp_get_max_threads ());
    #endif
  #endif
#else
  INFO ("Using internal statistical functions instead of GNU Scientific Library (GSL)");
  #ifdef _OPENMP
    #undef _OPENMP
    #warning "Cannot use OpenMP when using internal statistical functions.");
    INFO ("OpenMP is DISABLED when using internal statistical functions.");
  #endif
#endif

  swStart (overallSW);
#ifdef GCCOPT
  INFO ("GCC optimization level %d enabled", GCCOPT);
#else
  INFO ("GCC optimization disabled (or GCCOPT not defined)");
#endif

#ifdef PTMALLOC3
  INFO ("Using alternative allocator ptmalloc3");
#endif

  INFO ("To check status (at some risk), type CTRL-\\ or type \"kill -%d %d\"", SIGQUIT, (int) getpid ());

  /* Tolerate a request for help as the first argument; ignore the rest of the command 
   * line if help is requested. Otherwise, argv[1] better be the name of the configuration
   * file; the rest of the command line will be treated as override directives.
   */
  if (argc > 1) {
    if (!(strcmp (argv[1], "--help") && strcmp (argv[1], "-?"))) {
      fprintf (stderr, "usage: %s <conffile> [--directive arg1 arg2... [--directive...]]\n", argv[0]);
      exit (0);
    } else
      strcpy (configfile, argv[1]);
  } else
    ERROR ("usage: %s <conffile> [--directive arg1 arg2... [--directive...]]\n", argv[0]);

  INFO ("Using configuration file %s", configfile);

  /* Set modelRange, modelOptions and modelType to default values */
  initializeDefaults ();

  SUBSTEP(0, "Processing analysis configuration");

  /* Parse the configuration file. */
  DETAIL(0, "Read and process config file %s", argv[1]);
  readConfigFile (argv[1]);

  /* If there's anything on the command line after the configuration file name, 
   * it must be override directives.
   */
  if (argc > 2)
    parseCommandLine (argc - 2, &argv[2]);

  /* Make sure the config as read from the configuration file, and possibly modified on 
   * the command line, is legal. Then clean up the bits of memory allocated during
   * config parsing.
   */
  DETAIL(0,"Validating configuration");
  validateConfig ();

  /* This fills in defaults for fields that are optional, and also copies the contents
   * of the static Model{Range,Options,Type} structures to their global page-allocated
   * counterparts so we can protect them from monkeying.
   */
  /* TODO: Protect all the dynamically allocated memory referenced through these structures too */
  modelOptions = (ModelOptions *) allocatePages (sizeof (ModelOptions));
  modelRange = (ModelRange *) allocatePages (sizeof (ModelRange));
  modelType = (ModelType *) allocatePages (sizeof (ModelType));

  fillConfigDefaults (modelRange, modelOptions, modelType);

  if (modelOptions->polynomial == TRUE) {
    INFO ("Computation is done in polynomial mode");
#ifdef POLYUSE_DL
    INFO ("Dynamic libraries for polynomial evaluation will be used if found");
#endif
    polynomialInitialization (modelOptions->polynomialScale);
  } else
    INFO ("Computation is done in non-polynomial (direct evaluation) mode");
  if (modelOptions->integration == TRUE)
    INFO ("Integration is done numerically (dkelvin)");
  else
    INFO ("Integration is done with iteration (original kelvin)");

  if (swProgressDelaySeconds > 0)
    swStartProgressWakeUps(swProgressDelaySeconds); // Make the configuration value have an effect

  /* Read in the map file. */
  DETAIL(0,"Read and process map file %s", modelOptions->mapfile);
  read_mapfile (modelOptions->mapfile, modelOptions->mapFlag == SS);

  /* Initialize the locus list and read in the marker file. */
  memset (&savedLocusList, 0, sizeof (savedLocusList));
  memset (&markerLocusList, 0, sizeof (markerLocusList));
  memset (&traitLocusList, 0, sizeof (traitLocusList));
  memset (&originalLocusList, 0, sizeof (originalLocusList));
  /* read in what loci are in the pedigree file */
  DETAIL(0,"Read and process locus file %s", modelOptions->datafile);
  read_datafile (modelOptions->datafile);

  if (originalLocusList.numTraitLocus == 0 && !modelOptions->markerAnalysis)
    ERROR("No trait information in %s, can't run trait analysis\n", 
	  modelOptions->datafile);

  /* Read in marker allele frequencies */
  DETAIL(0,"Read and process marker file %s", modelOptions->markerfile);
  read_markerfile (modelOptions->markerfile, modelType->numMarkers);
#ifdef DISTRIBUTION
  if (modelRange->microsats && modelOptions->equilibrium == LINKAGE_DISEQUILIBRIUM)
    ERROR("LD analysis not supported with microsatellite datasets\n");
#endif

  /* build allele set information */
  DETAIL(0, "Constructing allele set");
  for (locus = 0; locus < originalLocusList.numLocus; locus++)
    construct_original_allele_set_list (locus);

  /* Initialize the pedigree set datastructure and read in the pedigrees. */
  memset (&pedigreeSet, 0, sizeof (PedigreeSet));
  DETAIL(0,"Read and process pedigree file %s", modelOptions->pedfile);
  read_pedfile (modelOptions->pedfile, &pedigreeSet);

  if (!modelOptions->markerAnalysis) {
    /* We are not doing marker to marker analysis; the configuration
     * has all the information about the disease trait if any.
     * Assume the traitLocus is 0 for now  - Need to fix this later */
    traitLocus = 0;
    pLocus = originalLocusList.ppLocusList[traitLocus];
    /* Set the global pTrait */
    pTrait = pLocus->pTraitLocus->pTraits[0];
  }

  if (modelRange->nlclass > 1) {
    int va, vb=0;
    
    if (! pTrait->numLiabilityClass)
      ERROR ("Liability class analysis specified, but no class information in dataset.");

    for (va = 1; va <= modelRange->nlclass; va++) {
      if (pedigreeSet.liabilityClassCnt[va] != 0) {
	modelRange->lclassLabels[vb] = va;
	pedigreeSet.liabilityClassCnt[va] = ++vb;
      }
    }
    if (vb != modelRange->nlclass) {
      va = modelRange->nlclass - vb;
      WARNING("%d liability class%s empty. Dropping empty classes to improve performance", va, (va == 1) ? " is" : "es are");
      
      for (va = 1; va <= modelRange->nlclass; va++) {
	if (pedigreeSet.liabilityClassCnt[va] == 0)
	  WARNING("Dropping empty class %d", va);
      }
      if (pedigreeSet.liabilityClassCnt[vb] != vb)
	renumberLiabilityClasses (&pedigreeSet);
      modelRange->nlclass = vb;
    }
    if (modelRange->nlclass == 0)
      ERROR ("Liability class analysis specified, but all classes are empty.");
  }

  /* sort, uniquify and expand the trait model dimensions, subject to constraints */
  DETAIL(0,"Post-processing model and configuration data");
  finishConfig (modelRange, modelType);
  if (modelOptions->integration != TRUE && modelOptions->markerAnalysis == FALSE) {
    if (modelType->trait == DT && modelRange->npenet == 0)
      ERROR ("No Penetrance values left after application of Constraints");
    if (modelType->trait == QT || modelType->trait == CT) {
      if (modelType->distrib == QT_FUNCTION_T && modelRange->npenet == 0)
	ERROR ("No Mean values left after application of Constraints");
      if (modelType->distrib == QT_FUNCTION_T && modelRange->nparam == 0)
	ERROR ("No StandardDeviation values left after application of Constraints");
      if (modelType->distrib == QT_FUNCTION_CHI_SQUARE && modelRange->npenet == 0)
	ERROR ("No DegreesOfFreedom values left after application of Constraints");
      if (modelType->trait == CT && modelRange->ntthresh == 0)
	ERROR ("No Threshold values left after application of Constraints");
    }
  }

  /* Calculate sample mean and sample standard deviation for QT/CT T distrib, if needed */
  if (modelType->trait != DT && modelType->distrib == QT_FUNCTION_T && (modelType->mean == -DBL_MAX || modelType->sd == -DBL_MAX)) {
    double mean, stdev;
    
    getPedigreeSampleStdev (&pedigreeSet, &mean, &stdev);
    modelType->mean = mean;
    modelType->sd = stdev;
    DIAG (INPUTFILE, 3, {printf ("Sample Mean is %.4f, Standard Deviation is %.4f\n", mean, stdev);});
  }

  /* read in case control file if provided */
  if (strlen (modelOptions->ccfile) > 0) {
    DETAIL(0,"Reading case control file %s", modelOptions->ccfile);
    read_ccfile (modelOptions->ccfile, &pedigreeSet);
  }

  SUBSTEP(0, "Initializing analysis data");
  if (modelType->trait == QT || modelType->trait == CT) {
    modelType->min = (modelType->minOriginal - modelType->mean) / modelType->sd;
    modelType->max = (modelType->maxOriginal - modelType->mean) / modelType->sd;
    pTrait->minFlag = modelType->minFlag;
    pTrait->maxFlag = modelType->maxFlag;
    pTrait->min = modelType->min;
    pTrait->max = modelType->max;
    pTrait->functionQT = modelType->distrib;
    if (modelType->distrib == QT_FUNCTION_T)
      pTrait->dfQT = modelType->constants[0];
    pTrait->sampleMean = modelType->mean;
    pTrait->sampleSD = modelType->sd;
    pTrait->unknownTraitValue = modelOptions->affectionStatus[AFFECTION_STATUS_UNKNOWN];
    pTrait->lessCutoffFlag = modelOptions->affectionStatus[AFFECTION_STATUS_UNAFFECTED];
    pTrait->moreCutoffFlag = modelOptions->affectionStatus[AFFECTION_STATUS_AFFECTED];
    if (modelType->distrib == QT_FUNCTION_T)
      adjustQuantitativeTraits (&pedigreeSet);
    if (checkQtTraitRanges (&pedigreeSet) == -1)
      ERROR("Can't run analysis with illegal trait data");
  }

  /* FIXME: shouldn't this bit come BEFORE the !markerAnalysis block, above? */
  if (modelType->trait == QT) {
    /* threshold value will not be used in any meaningful way, but we will use it for 
     * the loop */
    modelRange->ntthresh = 1;
    modelType->minOriginal = 0;
    modelType->maxOriginal = 1;
    if (modelRange->tthresh == NULL) {
      MALCHOKE (modelRange->tthresh, sizeof (double *), double **);
      for (i = 0; i < modelRange->nlclass; i++) {
        MALCHOKE (modelRange->tthresh[i], sizeof (double), void *);
      }
    }
  }

  /* allocate space for results */
  if (modelType->type == TP) {
    modelType->numMarkers = 1;
    totalLoci = 2;
    /* two point analysis */
    if (modelOptions->equilibrium == LINKAGE_EQUILIBRIUM) {
      /* in order to simplify looping, even for LE, we add a fake LD parameter dprime=0, which
       * is LE */
      modelRange->ndprime = 1;
      CALCHOKE (modelRange->dprime, (size_t) 1, sizeof (double), double *);
      modelRange->dprime[0] = 0;
      pLambdaCell = findLambdas (modelRange, 2, 2);
    }
  } else {
    /* we are doing multipoint analysis */
    totalLoci = modelType->numMarkers + originalLocusList.numTraitLocus;
    if (modelRange->tlocRangeStart >= 0) {
      double endofmap, tloc;
      endofmap = map.ppMapUnitList[map.count - 1]->mapPos[MAP_POS_SEX_AVERAGE] + modelRange->tlocRangeIncr;
      i = 0;
      while ((tloc = modelRange->tlocRangeStart + (i * modelRange->tlocRangeIncr)) <= endofmap) {
        addTraitLocus (modelRange, tloc);
        i++;
      }
    }
    if (modelRange->tlmark == TRUE) {
      /* add marker positions to the list of positions we want to conduct analysis */
      for (i = 0; i < originalLocusList.numLocus; i++) {
        pLocus = originalLocusList.ppLocusList[i];
        if (pLocus->locusType == LOCUS_TYPE_TRAIT)
          continue;
        addTraitLocus (modelRange, pLocus->pMapUnit->mapPos[MAP_POS_SEX_AVERAGE]);
      }
    }
  }
#ifdef STUDYDB

  DETAIL(0,"Opening study database");
  initializeDB ();
  if (toupper(*studyDB.role) == 'S')
    SignOn (originalLocusList.ppLocusList[1]->pMapUnit->chromosome, "ES", modelType->numMarkers, programVersion);
  else if (toupper(*studyDB.role) == '2') {
    SignOn (originalLocusList.ppLocusList[1]->pMapUnit->chromosome, "es", 1, programVersion);
    DIAG (ALTLSERVER, 0, { fprintf (stderr, "Explicitly setting all trait and marker-set models for all pedigrees for this server instance to 1.0\n");});
    SetDummyNullLikelihood (); // Set all trait and marker-set models for all pedigrees for this server instance to 1.0.
  }
#endif

  /* Enable handling of segmentation faults/bus errors due to configuration monkeying */
  setupSegvHandler ();
  allowReadOnly (modelOptions, sizeof (ModelOptions));
//  allowReadOnly (modelRange, sizeof (ModelRange));
  allowReadOnly (modelType, sizeof (ModelType));

  /* Estimate number of calls to each (appropriate) instance of compute_likelihood for
   * use in progress reporting, and display model information at this point since markers have
   * already been added to locus list */
  estimateIterations (eCL);

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
  CALCHOKE (tmpID, (size_t) totalLoci, sizeof (char), char *);

  /* initialize loci by doing genotype elimination, set recoding */
  initialize_loci (&pedigreeSet);

  if (modelOptions->dryRun != 0) {
    for (loc1 = 0; loc1 < originalLocusList.numLocus; loc1++) {
      fprintf (stderr, "Locus %d:\n", loc1);
      for (i = 0; i < pedigreeSet.numPedigree; i++) {
        Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[i];
        print_pedigree_locus_genotype_count (pPedigree, loc1);
      }
    }
  }
  if (modelOptions->polynomial == TRUE) {
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

  for (i = 0; i < pedigreeSet.numPedigree; i++) {
    Pedigree *pPedigree = pedigreeSet.ppPedigreeSet[i];
    pPedigree->load_flag = 0;   /* Initially 0 and changes to 1 when marker or 
                                 * alternative likelihood values are retrieved */
  }

  /* Open output files that get written across loops. */

  SUBSTEP(0, "Opening cross-loop output files");

//  if (modelOptions->conditionalRun == 1 || modelOptions->loopCondRun == 1) {
//    fpCond = fopen (modelOptions->condFile, "w");
//    ASSERT (fpCond != NULL, "Error in opening file %s for write.\n", modelOptions->condFile);
//  }

  if (modelOptions->markerAnalysis == FALSE || modelOptions->forceAvghetFile == TRUE) {
    fpHet = fopen (modelOptions->avghetfile, "w");
    ASSERT (fpHet != NULL, "Error in opening file %s for write.\n", modelOptions->avghetfile);
    fprintf (fpHet, "# Version %s edit %s\n", programVersion, svnVersion);
  }

  if (modelType->type == TP) {
    fpPPL = fopen (modelOptions->pplfile, "w");
    ASSERT (fpPPL != NULL, "Error in opening file %s for write.\n", modelOptions->pplfile);
    writePPLFileHeader ();
  }

  if (strlen (modelOptions->modfile) > 0) {
    fpMOD = fopen (modelOptions->modfile, "w");
    ASSERT (fpMOD != NULL, "Error in opening file %s for write.\n", modelOptions->modfile);
    fprintf (fpMOD, "# Version %s edit %s\n", programVersion, svnVersion);
  }

  if (strlen (modelOptions->intermediatefile) > 0) {
    fpIR = fopen (modelOptions->intermediatefile, "w");
    ASSERT (fpIR != NULL, "Error in opening file %s for write.\n", modelOptions->intermediatefile);
  }
  // DKelvin intermediate results are written here.
  if ((modelOptions->integration) && (strlen (modelOptions->dkelvinoutfile) > 0)) {
    fpDK = fopen (modelOptions->dkelvinoutfile, "w");
    ASSERT (fpDK != NULL, "Error in opening file %s for write.\n", modelOptions->dkelvinoutfile);
  }

  R_square_flag = 0;
  R_square = FALSE;
  leftMarker = -1;
  traitIndex = 0;

}
