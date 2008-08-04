/**
@file trackProgress.c

  kelvin support routines for estimating and tracking progress thru analyses.

  In-progress.

*/
#include "kelvin.h"
#include "polynomial.h"
#include "trackProgress.h"

/** 

Make sure this thread is started AFTER signal handling has been setup.

Loop sleeping for 30 seconds and then:

- signal with SIGUSR1 to do allow a synchronous dump of statically-collected
statistitics.
- optionally write performance statistics to a file for graphing with gnuplot.
- optionally display dynamic statistics.

*/
void *monitorStatus () {

  int currentVMK, maximumVMK;
  time_t startTime;

  startTime = time (NULL);

  if ((maximumVMK = swGetMaximumVMK ()) != 0) {
#ifdef MEMGRAPH
    FILE *graphFile;
    char graphFileName[64];
    sprintf (graphFileName, "kelvin_%d_memory.dat", getpid ());
    if ((graphFile = fopen (graphFileName, "w")) == NULL) {
      perror ("Cannot open memory graph file!");
      exit (EXIT_FAILURE);
    }
#endif
    int wakeCount = 0;
    while (1) {
      sleep (30);
      wakeCount++;
      if (wakeCount % 2)
	kill (getppid (), SIGUSR1);   // Send a status-updating signal to parent.
      currentVMK = swGetCurrentVMK (getpid ());
#ifdef MEMGRAPH
      fprintf (graphFile, "%lu, %d, %d\n", time (NULL) - startTime, currentVMK, nodeId);
      fflush (graphFile);
#endif
#ifdef MEMSTATUS
      fprintf (stdout, "%lus, %dKb (%.1f%% of %.1fGb) at %d\n", time (NULL) - startTime, currentVMK,
	       currentVMK / (maximumVMK / 100.0), maximumVMK / (1024.0 * 1024.0), nodeId);
#endif
    }
  }
}

/**

  Dump out statistics for estimating the complexity of the pedigrees
  involved in the analysis.

*/
void print_dryrun_stat (PedigreeSet *pSet, ///< Pointer to set of pedigrees in analysis
			double pos ///< Position being analyzed.
)
{
  int pedIdx;
  long subTotalPairGroups, subTotalSimilarPairs;
  long totalPairGroups, totalSimilarPairs;
  NuclearFamily *pNucFam;
  Pedigree *pPedigree;
  int i;

  totalPairGroups = 0;
  totalSimilarPairs = 0;
  for (pedIdx = 0; pedIdx < pSet->numPedigree; pedIdx++) {
    /* save the likelihood at null */
    pPedigree = pSet->ppPedigreeSet[pedIdx];
    fprintf (stderr, "Ped %s(%d) has %d loops, %d nuclear families.\n",
             pPedigree->sPedigreeID, pedIdx, pPedigree->numLoop, pPedigree->numNuclearFamily);
    subTotalPairGroups = 0;
    subTotalSimilarPairs = 0;
    for (i = 0; i < pPedigree->numNuclearFamily; i++) {
      pNucFam = pPedigree->ppNuclearFamilyList[i];
      fprintf (stderr,
               "    Nuc %d w/ proband %s(%s) has %ld unique pp groups, %ld similar pp, total %ld.\n",
               i, pNucFam->pProband->sID,
               pNucFam->childProbandFlag ? "child" : "parent",
               pNucFam->totalNumPairGroups,
               pNucFam->totalNumSimilarPairs, pNucFam->totalNumPairGroups + pNucFam->totalNumSimilarPairs);
      subTotalPairGroups += pNucFam->totalNumPairGroups;
      subTotalSimilarPairs += pNucFam->totalNumSimilarPairs;
    }
    fprintf (stderr,
             "    Ped has total %ld unique pp groups, %ld similar pp, total %ld.\n",
             subTotalPairGroups, subTotalSimilarPairs, subTotalPairGroups + subTotalSimilarPairs);
    totalPairGroups += subTotalPairGroups;
    totalSimilarPairs += subTotalSimilarPairs;
  }
  fprintf (stderr,
           "POS %f has %ld unique pp groups, %ld similar pp, total %ld.\n",
           pos, totalPairGroups, totalSimilarPairs, totalPairGroups + totalSimilarPairs);
}

/**

  Log position and pedigree complexity statistics as produced by dry-run.

  At each position, after the determination of family pair groupingss, summarize and
  log the complexity data.

*/
void logPedigreeSetStatistics (PedigreeSet *pSet, ///< Pointer to pedigree set to be described ala dry-run.
	       int posIdx ///< Position for complexity analysis.
	       )
{
  char messageBuffer[MAXSWMSG];
  int pedIdx, i;
  NuclearFamily *pNucFam;
  Pedigree *pPedigree;
  int l = 0, nf = 0, pg = 0, sg = 0;

  for (pedIdx = 0; pedIdx < pSet->numPedigree; pedIdx++) {
    pPedigree = pSet->ppPedigreeSet[pedIdx];
    l += pPedigree->numLoop;
    nf += pPedigree->numNuclearFamily;
    for (i = 0; i < pPedigree->numNuclearFamily; i++) {
      pNucFam = pPedigree->ppNuclearFamilyList[i];
      pg += pNucFam->totalNumPairGroups;
      sg += pNucFam->totalNumSimilarPairs;
    }
  }
  sprintf (messageBuffer, "For %d pedigrees: unique groups:%d, similar groups:%d, polynomial terms:%d",
	   pSet->numPedigree, pg, sg, nodeId);
  swLogMsg (messageBuffer);
}

char analysisType[MAXSWMSG]; ///< Textual summary of analysis type built by dumpTrackingStats.
/**

  Derive textual summary of analysis type and build expected compute_likelihood() call counts.

*/
void dumpTrackingStats(ModelType modelType, ModelOptions modelOptions, ModelRange modelRange, int cl[], int eCL[])
{
  int i;
  fprintf (stderr, "compute_pedigree_likelihood counts: ");
  for (i=0; i<9; i++)
    fprintf (stderr, "%d=%d(%d) ", i, cl[i], eCL[i]);
  fprintf (stderr, "\n");
  fprintf (stderr, "modelRange. ntloc %d, npenet %d, nlclass %d, ngfreq %d, nafreq %d, "
	   "nparam %d, ntthresh %d, nalpha %d, ntheta %d, ndprime %d, originalLocusList.numLocus %d, "
	   "modelType.numMarkers %d\n",
	   modelRange.ntloc, modelRange.npenet, modelRange.nlclass, modelRange.ngfreq, 
	   modelRange.nafreq, modelRange.nparam, modelRange.ntthresh, modelRange.nalpha,
	   modelRange.ntheta, modelRange.ndprime, originalLocusList.numLocus, modelType.numMarkers);
}

// Construct string describing the type of analysis and determine evaluations required.
char *estimateIterations (ModelType modelType, ModelOptions modelOptions, ModelRange modelRange, int eCL[])
{
  //  int cL[9];
  //  dumpTrackingStats(modelType,  modelOptions,  modelRange, cL, eCL);
  if (modelOptions.markerAnalysis != FALSE) {
    /*
      Marker pair (not # in analysis, but locus list)
      Marker allele frequencies and penetrances stay at 1
      Theta and D' are still involved
    */
    eCL[0] = 0;
    eCL[1] = (originalLocusList.numLocus-2) * modelRange.ndprime * modelRange.ntheta;
    sprintf (analysisType, "Marker-to-marker, Linkage ");
    strcat (analysisType, (modelOptions.equilibrium == 
			    LINKAGE_EQUILIBRIUM) ? "Equilibrium." : "Disequilibrium.");
  } else { // not AM/MM
    if (modelType.type == TP) {
      /* 

      TP DT NULL hypothesis is cL[0], looped for marker pair, marker allele frequency (not really), modelRange. ngfreq, npenet
      TP DT alternative hypothesis is cL[1], looped for all of cL[0] and ndprime, ntheta
      TP QT NULL hypothesis is cL[2], looped for marker pair, marker allele frequency (not really),  modelRange. ngfreq, nparam, npenet, ntthresh
      TP QT alternative hypothesis is cL[3], looped for all of cL[2] and ndprime, ntheta

      */      
      sprintf (analysisType, "%da*%dgf*%dp*%dlc space, Trait-to-marker, Two-Point, ",
	       modelRange.nalpha, modelRange.ngfreq, modelRange.npenet, modelRange.nlclass);
      if (modelType.trait == DT) {
	strcat (analysisType, "Dichotomous Trait, ");
	eCL[0] = modelRange.ngfreq * modelRange.npenet * (originalLocusList.numLocus-1);
	eCL[1] = eCL[0] * modelRange.ndprime * modelRange.ntheta;
      } else { // TP not DT
	eCL[2] = modelRange.ngfreq * modelRange.npenet * (originalLocusList.numLocus-1) * modelRange.nparam * modelRange.ntthresh;
	eCL[3] = eCL[2] * modelRange.ndprime * modelRange.ntheta;
	if (modelType.trait == QT) { //QT
	  strcat (analysisType, "Quantitative Trait, ");
	} else { // TP not DT or QT, so CT
	  strcat (analysisType, "Quantitative Trait w/Threshold, ");
	}
	if (modelType.distrib == QT_FUNCTION_T) {
	  strcat (analysisType, "Truncated Student's T-Distribution, Linkage ");
	} else { // not T-Dist
	  if (modelType.distrib == QT_FUNCTION_CHI_SQUARE) {
	    strcat (analysisType, "Chi-Square Distribution, Linkage ");
	  } else { // not T-Dist or Chi-Sq Dist, so Normal Dist
	    strcat (analysisType, "Normal Distribution, Linkage ");
	  }
	}
      }
      strcat (analysisType, (modelOptions.equilibrium == 
			      LINKAGE_EQUILIBRIUM) ? "Equilibrium." : "Disequilibrium.");
    } else { // not TP, so multipoint

      /* Pedigree likelihood calculation looping for MP is for trait, marker, then alternative hypothesis.

	 modelRange.ntloc is all evaluation locations on the chromosome since we're past adding for TM.

      MP DT trait is cl[4], looped for 3(SA/SS)?, modelRange. npenet, nlclass, ngfreq
      MP QT/CT trait is cl[5], looped for 3(SA/SS)?, modelRange. ngfreq, nparam (QT dist), npenet, ntthresh
      MP marker is cl[6], looped for <= modelRange.ntloc
      MP DT alt is cl[7], looped for modelRange.ntloc, locusList->numLocus, modelRange. npenet, ngfreq
      MP QT/CT alt is cl[8], looped for modelRange.ntloc,  locusList->numLocus, 
              modelRange. ngfreq, nparam, npenet, ntthresh
      ...but locusList->numLocus is calculated, use modelRange.ntloc
      
      polynomials for each pedigree incorporate alpha?, # MP markers used in analysis.

      */
      sprintf (analysisType, "%da*%dgf*%dp*%dlc space, Trait-to-marker, Sex-%s Multipoint (w/%d loci), ",
	       modelRange.nalpha, modelRange.ngfreq, modelRange.npenet, modelRange.nlclass,
	       modelOptions.mapFlag == SS ? "Specific" : "Averaged",
	       modelType.numMarkers + originalLocusList.numTraitLocus);
      eCL[6] = modelRange.ntloc;
      if (modelType.trait == DT) {
	strcat (analysisType, "Dichotomous Trait.");
	eCL[4] = modelRange.npenet * modelRange.nlclass * modelRange.ngfreq;
	eCL[7] = modelRange.ntloc * modelRange.npenet * modelRange.ngfreq;
      } else { // SA/SS multipoint, but not DT
	eCL[5] = modelRange.npenet * modelRange.ntthresh * modelRange.ngfreq * 
	  modelRange.nparam;
	eCL[8] = modelRange.ntloc * modelRange.npenet * modelRange.ngfreq * 
	  modelRange.ntthresh * modelRange.nparam;
	if (modelType.trait == QT) {
	  strcat (analysisType, "Quantitative Trait, ");
	} else { // SA/SS multipoint but not DT or QT, so CT
	  strcat (analysisType, "Quantitative Trait w/Threshold, ");
	}
	if (modelType.distrib == QT_FUNCTION_T) {
	  strcat (analysisType, "Truncated Student's T-Distribution.");
	} else { // not T-Dist
	  if (modelType.distrib == QT_FUNCTION_CHI_SQUARE) {
	    strcat (analysisType, "Chi-Square Distribution.");
	  } else { // not T-Dist or Chi-Sq dist, so Normal Dist
	    strcat (analysisType, "Normal Distribution.");
	  }
	}
      }
    }
  }
  return (analysisType);
}
