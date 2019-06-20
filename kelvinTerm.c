#include "utils/utils.h"
#include "utils/sw.h"
#include "kelvinGlobals.h"
#include "summary_result.h"
#include "ppl.h"

#ifdef STUDYDB
#include "database/StudyDB.h"
extern struct StudyDB studyDB;
#include "database/databaseSupport.h"
#endif

int posIdx;
extern LDLoci *pLDLoci;

void kelvinTerm () {

  int i;

  set_removeGenotypeFlag (TRUE);

  if (modelType->type == TP) {
    if (tp_result != NULL) {
      /* two point */
      for (i = 0; i < pLambdaCell->ndprime; i++) {
        free (tp_result[i]);
      }
      free (tp_result);
    }
  } else {
    /* multipoint */
    for (posIdx = 0; posIdx < modelRange->ntloc; posIdx++) {
      free (mp_result[posIdx].pMarkers);
    }
    free (mp_result);
  }

  /* free parental pair work space */
  free_parental_pair_workspace (&parentalPairSpace, modelType->numMarkers + 1);

  /* free transmission probability matrix */
  free (altMatrix);
  free (nullMatrix);

  if (modelType->type == TP)
    free_LD_loci (pLDLoci);

  if (modelOptions->polynomial == TRUE) {
    pedigreeSetPolynomialClearance (&pedigreeSet);
  }
  free_likelihood_storage (&pedigreeSet);
  free_likelihood_space (&pedigreeSet);
  free_pedigree_set (&pedigreeSet);
  free_sub_locus_list (&traitLocusList);
  free_sub_locus_list (&markerLocusList);
  free_sub_locus_list (&savedLocusList);
  free (modelOptions->sUnknownPersonID);
  final_cleanup ();

#ifdef STUDYDB
  if (toupper(*studyDB.role) == 'S') {
    DETAIL(0,"Signing out of database");
    SignOff (0);
  }
#endif

#ifdef SOURCEDIGRAPH
  if (modelOptions->polynomial == TRUE)
    dumpSourceParenting ();
#endif

  /* Final dump and clean-up for performance. */
  swStop (overallSW);
#ifndef DISTRIBUTION
  swDump (overallSW);
#endif
#ifdef POLYSTATISTICS
  if (modelOptions->polynomial == TRUE)
    polyStatistics ("End of run");
#endif
#ifdef DMTRACK
  swLogPeaks ("End of run");
  swDumpHeldTotals ();
  swDumpSources ();
  //  swDumpCrossModuleChunks ();
#endif
  STEP(0, "Finished run");

  if (modelOptions->dryRun != 0)
    INFO("Single model total parental pairs: %ld, peak per-polynomial parental pairs is %ld", grandTotalPairs, peakTotalPairs);

/* Close file pointers */
if (modelType->type == TP)
  fclose (fpPPL);

if (fpMOD != NULL)
  fclose (fpMOD);
if (fpHet != NULL)
  fclose (fpHet);
if (fpDry != NULL)
  fclose (fpDry);
}
