/**
@file trackProgress.c

  kelvin support routines for estimating and tracking progress thru analyses.

  Estimating kelvin progress-to-completion is a lot more complicated that 
  you might expect.  Here are the problems:

  - There are multiple looping paths through the code, i.e. each different 
  type of analysis uses a different set of calls to compute_likelihood, and 
  each needs to be instrumented differently.

  - Multipoint analysis invokes compute_likelihood in three places:
    - trait likelihood (once for the entire run)
    - marker set likelihood (once for each distinct set of loci)
    - combined alternative and null likelihood (once for each position).

  We need to show progrees in the marker set loop as well as the combined
  one because marker set evaluation is our first real indication of 
  complexity, and sometimes it takes a long time.

  - Polynomial construction takes an arbitrarily large amount of time. While
  evaluation is iterative and therefore predictable, construction is a big
  unknown, even between positions in the same analysis, so when polynomial
  evaluation is requested, a progress graph can look extremely lumpy. There
  is currently no solution for this problem.

  - Integration approach (now the default) iterations are unpredictable
  because it re-partitions the trait space as much as needed to reduce error
  to tolerable amounts. We might know the upper-limit of iterations, but
  that can be wildly different from the actual number of iterations, and we
  don't want to scare people off by using that upper limit. There is currently
  no solution for this problem.

  - Position counting or loci pair counting can be used for progress tracking
  only so long as there are a reasonably large number of position or pairs
  being evaluated in a single run. No-one likes to see progress go from 0 to
  100% in one step at the end of the run, which is what would happen if a single
  position or pair is being evaluated.

  - Too little output during a slow run leaves the user in the dark about progress,
  while too much output during a complicated run overwhelms the user with details.
  The best approach would be to show the latest progress (assuming it has changed)
  at some interval that varies based upon how long the analysis has been running,
  or whenever the user requests status.

  - Finally, progress tracking for users must be simplistic and understandable,
  while internal users want the nitty-gritty details on what is going on at
  each step.

  I suggest that we no longer have a simple progress option, but
  rather display progress in a manner that gives as much detail as
  could be desired without either starving or glutting the user (or
  the log file). We can instrument the code to record internally as
  much detail as we can manage on latest step completed as it
  happens, but only display that information when the user requests,
  or at the run duration-related intervals that they specify. This makes
  output driven by the user's interest and the duration of the run
  instead of the complexity of the run or the duration of the
  individual steps.

  For example, assume we're doing a fixed grid multipoint run with
  polynomials for a large number of trait positions and a few simple
  pedigrees, so steps are performed quickly. The maximum level of 
  progress detail we normally provide would look like this:

    Building polynomial for trait likelihood, <minutes elapsed/estimated>
    Evaluating trait likelihood, <minutes elapsed/estimated>
    Building polynomial for new 3-marker set at trait position 1 of 140 (1.0cM), <minutes elapsed/estimated>
    Evaluating new 3-marker set likelihood at trait position 1 of 140 (1.0cM), <minutes elapsed/estimated>
    Building polynomial for combined likelihood at trait position 1 of 140 (1.0cM), <minutes elapsed/estimated>
    Evaluating combined likelihood at trait position 1 of 140 (1.0cM), <minutes elapsed/estimated>
    Evaluating combined likelihood at trait position 2 of 140 (2.0cM), <minutes elapsed/estimated>
    Evaluating combined likelihood at trait position 3 of 140 (3.0cM), <minutes elapsed/estimated>
    Evaluating combined likelihood at trait position 4 of 140 (4.0cM), <minutes elapsed/estimated>
    Building polynomial for new 3-marker set at trait position 5 of 140 (5.0cM), <minutes elapsed/estimated>
    Evaluating new 3-marker set likelihood at trait position 5 of 140 (5.0cM), <minutes elapsed/estimated>
    Building polynomial for combined likelihood at trait position 5 of 140 (5.0cM), <minutes elapsed/estimated>
    Evaluating combined likelihood at trait position 5 of 140 (5.0cM), <minutes elapsed/estimated>

    ...and so on for another thousand lines.

    ...and since it's fixed grid,  sprinkled throughout this at 2 minute intervals:
    "Total evaluations N% complete, <minutes elapsed/estimated>"

  This is just overwhelming for anyone, and worse, it can hide warning and
  error messages due to its sheer volume.

  But if we allow the user's needs to drive the output, it's much simpler.
  Say the user specifies (or accepts as default) progress displays at 1 
  minute intervals for the first 5 minutes, 5 minute intervals for the 
  next 30 minutes, and 30 minute intervals from then on. Here's what might
  be displayed for the same run if it takes 12 minutes:

    Total 1m elapsed, estimated 9m remaining, currently 17s into:
    Building polynomial for combined likelihood at trait position 14 of 140 (14.0cM).
    Total 2m elapsed, estimated 8m remaining, currently 8s into:
    Evaluating combined likelihood at trait position 30 of 140 (30.0cM).
    Total 3m elapsed, estimated 8m remaining, currently 12s into:
    Building polynomial for combined likelihood at trait position 42 of 140 (42.0cM).
    Total 4m elapsed, estimated 7m remaining, currently 6s into:
    Evaluating combined likelihood at trait position 56 of 140 (56.0cM).
    Total 5m elapsed, estimated 6m remaining, currently 1s into:
    Evaluating new 3-markers set likelihood at trait position 71 of 140 (71.0cM).
    Total 10m elapsed, estimated 1m remaining, currently 21s into:
    Building polynomial for combined likelihood at trait position 130 of 140 (130.0cM).

  If the run took less than a minute, none of this would be displayed, which
  is nice, because no-one would care about it.
  
  This approach is even more beneficial when we consider dynamic grid.
  Every significant step that the dynamic grid analysis takes, such as
  splitting a region, could be made available for reporting, but it
  wouldn't be overwhelming like it currently is even though its volume
  remains largely unpredicatable. Regardless of how an analysis was being 
  done, you would always get the same frequency of updates, and warning and
  error messages would still be obvious.

  For diagnostic purposes, we could configure the status display frequency
  to be 0 minutes to ensure that every step is displayed as it occurs.

  If we want an even simpler display, we could have a configuration parameter
  to eliminate the display of the current step while keeping the estimation
  output intervals.


  So here's what we can provide:

  For everything (mentioning "total" as opposed to stepwise):

  <minutes elapsed/estimated> as shown in the following examples will be either:
  "X total minutes elapsed, cannot estimate remaining time." (when no steps complete)
  "X total minutes elapsed, estimated Y remaining overall."

  For Simple Progress:

    Fixed or Dynamic Grid:

      Two-Point:

        Progress thru the N marker pairs will be displayed:

	"Building polynomial for marker pair N of M (Name1 and Name2), <minutes elapsed/estimated>"
	"Evaluating marker pair N of M (Name1 and Name2), <minutes elapsed/estimated>"

      Multi-Point:

        Progress thru the N positions will be displayed:

	"Building polynomial for trait likelihood, <minutes elapsed/estimated>"
	"Evaluating trait likelihood, <minutes elapsed/estimated>"
	"Building polynomial for new N-marker set at trait position N of M (Position), <minutes elapsed/estimated>"
	"Evaluating new N-marker set likelihood at trait position N of M (Position), <minutes elapsed/estimated>"
	"Building polynomial for combined likelihood at trait position N of M (Position), <minutes elapsed/estimated>"
	"Evaluating combined likelihood at trait position N of M (Position), <minutes elapsed/estimated>"

  For Detailed Progress:

    Fixed Grid:

      Two-Point:

        Progress thru the N marker pairs and evaluations every 2 minutes will be displayed:

	"Building polynomial for marker pair N of M (Name1 and Name2), <minutes elapsed/estimated>"
	"Evaluating marker pair N of M (Name1 and Name2), <minutes elapsed/estimated>"

	...and at 2 minute intervals:

	"<minutes elapsed/estimated>" (when no evaluations complete)
	"Total evaluations N% complete, <minutes elapsed/estimated>"

      Multi-Point:

        Progress thru trait, new marker sets, and the N positions will be displayed:

	"Building polynomial for trait likelihood, <minutes elapsed/estimated>"
	"Evaluating trait likelihood, <minutes elapsed/estimated>"
	"Building polynomial for new N-marker set at trait position N of M (Position), <minutes elapsed/estimated>"
	"Evaluating new N-marker set likelihood at trait position N of M (Position), <minutes elapsed/estimated>"
	"Building polynomial for combined likelihood at trait position N of M (Position), <minutes elapsed/estimated>"
	"Evaluating combined likelihood at trait position N of M (Position), <minutes elapsed/estimated>"
	...and at 2 minute intervals:
	"<minutes elapsed/estimated>" (when no evaluations complete)
	"Total evaluations N% complete, <minutes elapsed/estimated>"





	

  Fixed Grid in general:

  - Calculate iteration counts per compute_likelihood step based upon 
  trait space parameters and number of positions or loci pairs.

  Detailed Fixed Grid 2pt: calculate iteration counts per compute_likelihood step 
  based upon 
  trait space parameters. Show progress based solely upon percentage completion 
  of these counts. Lump polynomial build in with evaluation.

  Fixed Grid Multipoint:
  
  - 
  - Maintain stopwatches for overall time, polynomial build time, and evaluation time.
  - 

  Copyright &copy; 2009, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

*/
#include "kelvin.h"
#include "kelvinGlobals.h"
extern PedigreeSet pedigreeSet;
#include "kelvinHandlers.h"
#include "utils/polynomial.h"
#include "trackProgress.h"

/** 

Make sure this thread is started AFTER signal handling has been setup.

Loop sleeping for MONSTATDELAYSEC seconds and then:

- Optionally write performance statistics to a file for graphing with gnuplot.
- Optionally display dynamic statistics.
- Used to do thrashing checks and "program apoptosis", but that's commented out.
Thrashing was checked after some multiple of the MONSTATDELAYSEC periods.

  @author Bill Valentine-Cooper - overall content.
  @par Reviewers: None.
  @par Global Inputs

  - MEMGRAPH macro conditional. If defined, graph-friendly memory 
  utilization statistics are written to a file every MONSTATDELAYSEC seconds.

  - MEMSTATUS macro conditional. If defined, memory utilization statistics
  are written to stdout every MONSTATDELAYSEC seconds.

  @par Global Outputs

  - kelvin_pid_memory.dat if MEMGRAPH defined. This file has the current
  elapsed seconds, memory used in Kbytes, and polynomial nodeId appended
  to it every MONSTATDELAYSEC seconds. Can be plotted w/gnuplot as simply aa:

  <lit>gnuplot> plot "kelvin_pid_memory.dat" using 1:2 with lines;

  - memory status to stdout if MEMSTATUS defined.

  @return void

*/
void *monitorStatus ()
{

  long currentVMK, maximumPMK;
  time_t startTime;

  startTime = time (NULL);

  if ((maximumPMK = swGetMaximumPMK ()) != 0) {
#ifdef MEMGRAPH
    FILE *graphFile;
    char graphFileName[64];     // I define the file name, so this is long enough
    sprintf (graphFileName, "kelvin_%d_memory.dat", (int) getpid ());
    if ((graphFile = fopen (graphFileName, "w")) == NULL) {
      perror ("Cannot open memory graph file!");
      exit (EXIT_FAILURE);
    }
#endif
    int wakeCount = 0;
    while (1) {
      sleep (MONSTATDELAYSEC);
      wakeCount++;
      if (!(wakeCount % 2)) {
        /*      thrashingCheck (); */
        statusRequestSignal = TRUE;
//      kill (getpid (), SIGQUIT);   // Send a status-updating signal
      }
      currentVMK = swGetCurrentVMK (getpid ());
#ifdef MEMGRAPH
      fprintf (graphFile, "%lu, %d, %ld\n", time (NULL) - startTime, currentVMK, nodeId);
      fflush (graphFile);
#endif
#ifdef MEMSTATUS
      fprintf (stdout, "%lus, %dKb (%.1f%% of %.1fGb) at %d\n", time (NULL) - startTime, currentVMK, currentVMK / (maximumPMK / 100.0), maximumPMK / (1024.0 * 1024.0), nodeId);
#endif
    }
  }
  return NULL;
}

/**

  Dump out statistics for estimating the complexity of the pedigrees
  involved in the analysis.

*/
void print_dryrun_stat (PedigreeSet * pSet,     ///< Pointer to set of pedigrees in analysis
    double pos  ///< Position being analyzed.
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
    fprintf (stderr, "Ped %s(%d) has %d loops, %d nuclear families.\n", pPedigree->sPedigreeID, pedIdx, pPedigree->numLoop, pPedigree->numNuclearFamily);
    subTotalPairGroups = 0;
    subTotalSimilarPairs = 0;
    for (i = 0; i < pPedigree->numNuclearFamily; i++) {
      pNucFam = pPedigree->ppNuclearFamilyList[i];
      fprintf (stderr,
          "    Nuc %d w/ proband %s(%s) has %ld unique pp groups, %ld similar pp, total %ld.\n",
          i, pNucFam->pProband->sID, pNucFam->childProbandFlag ? "child" : "parent", pNucFam->totalNumPairGroups, pNucFam->totalNumSimilarPairs, pNucFam->totalNumPairGroups + pNucFam->totalNumSimilarPairs);
      subTotalPairGroups += pNucFam->totalNumPairGroups;
      subTotalSimilarPairs += pNucFam->totalNumSimilarPairs;
    }
    fprintf (stderr, "    Ped has total %ld unique pp groups, %ld similar pp, total %ld.\n", subTotalPairGroups, subTotalSimilarPairs, subTotalPairGroups + subTotalSimilarPairs);
    totalPairGroups += subTotalPairGroups;
    totalSimilarPairs += subTotalSimilarPairs;
  }
  fprintf (stderr, "POS %f has %ld unique pp groups, %ld similar pp, total %ld.\n", pos, totalPairGroups, totalSimilarPairs, totalPairGroups + totalSimilarPairs);
}

/**

  Log position and pedigree complexity statistics as produced by dry-run.

  At each position, after the determination of family pair groupings, summarize and
  log the complexity data.

*/
void logPedigreeSetStatistics (PedigreeSet * pSet,      ///< Pointer to pedigree set to be described ala dry-run.
    int posIdx  ///< Position for complexity analysis.
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
  sprintf (messageBuffer, "For %d pedigrees: unique groups:%d, similar groups:%d, polynomial terms:%ld", pSet->numPedigree, pg, sg, nodeId);
  swLogMsg (messageBuffer);
}

char analysisType[MAXSWMSG];    ///< Textual summary of analysis type built by dumpTrackingStats.
/**

  Derive textual summary of analysis type and build expected compute_likelihood() call counts.

*/
void dumpTrackingStats (unsigned long cl[], unsigned long eCL[])
{
  int i;
  fprintf (stderr, "compute_pedigree_likelihood counts: ");
  for (i = 0; i < 9; i++)
    fprintf (stderr, "%d=%lu(%lu) ", i, cl[i], eCL[i]);
  fprintf (stderr, "\n");
  fprintf (stderr, "modelRange-> ntloc %d, npenet %d, nlclass %d, ngfreq %d, nafreq %d, "
      "nparam %d, ntthresh %d, nalpha %d, ntheta %d, ndprime %d, originalLocusList.numLocus %d, "
      "modelType->numMarkers %d\n",
      modelRange->ntloc, modelRange->npenet, modelRange->nlclass, modelRange->ngfreq, modelRange->nafreq, modelRange->nparam, modelRange->ntthresh, modelRange->nalpha, modelRange->ntheta, modelRange->ndprime, originalLocusList.numLocus, modelType->numMarkers);
}

/** 
    Construct string describing the type of analysis and determine evaluations required. Note that
    the model numbers we're going to display along with the type of the analysis does not directly
    relate to the number of iterations we're concerned about. Iterations are about likelihood
    evaluation, while model numbers describe two things -- iterations over the model space which will
    show up as enumerated values in the results, and the model space itself which is everything that
    is variable but unenumerated in the results. As an example, alpha is definitely always a part of the 
    model space, but doesn't show up anywhere in the iterations for likelihood because it is iterated 
    over after combined likelihood is calculated. Another is pedigrees, which is always in the 
    iterations for likelihood but not a part of the model space (they don't vary) and not a part of the 
    enumerated results.

*/
char *estimateIterations (unsigned long eCL[])
{
  //  unsigned long cL[9];
  //  dumpTrackingStats(cL, eCL);
  int totalLoopsForDPrime = 0, loc1, loc2;
  Locus *pLocus1, *pLocus2;

  if (modelOptions->markerAnalysis != FALSE) {
    /*
     * Marker pair (not # in analysis, but locus list)
     * Marker allele frequencies and penetrances stay at 1
     * Theta and D' are still involved
     */

    for (loc1 = 0; loc1 < originalLocusList.numLocus - 1; loc1++) {
      pLocus1 = originalLocusList.ppLocusList[loc1];
      if (pLocus1->locusType != LOCUS_TYPE_MARKER)
        continue;
      for (loc2 = loc1 + 1; loc2 < originalLocusList.numLocus; loc2++) {
        pLocus2 = originalLocusList.ppLocusList[loc2];
        //      fprintf (stderr, "Adding loc1: %d (%d alleles) and loc2: %d (%d alleles)\n",
        //               loc1, pLocus1->numOriginalAllele, loc2, pLocus2->numOriginalAllele);
        if (pLocus2->locusType != LOCUS_TYPE_MARKER)
          continue;
        totalLoopsForDPrime += pow (modelRange->ndprime, (pLocus1->numOriginalAllele - 1) * (pLocus2->numOriginalAllele - 1));
        if (modelOptions->markerAnalysis == ADJACENTMARKER)
          loc1 = loc2;
      }
    }
    eCL[0] = 0;
    eCL[1] = totalLoopsForDPrime;
    sprintf (analysisType, "%dD' cases of %dAL*%dGF*%dpv(%dLC)' space for %d pedigree(s)\n" "Marker-to-marker Two-Point ", totalLoopsForDPrime, modelRange->nalpha, modelRange->ngfreq, modelRange->npenet, modelRange->nlclass, pedigreeSet.numPedigree);
    if (modelOptions->equilibrium == LINKAGE_EQUILIBRIUM)
      strcat (analysisType, "Equilibrium.");
    else
      strcat (analysisType, "Disequilibrium.");
  } else {      // not AM/MM
    if (modelType->type == TP) {
      /* 
       * 
       * TP DT NULL hypothesis is cL[0], looped for marker pair, marker allele frequency (not really), modelRange-> ngfreq, npenet
       * TP DT alternative hypothesis is cL[1], looped for all of cL[0] and allele pairs, ndprime, ntheta
       * TP QT NULL hypothesis is cL[2], looped for marker pair, marker allele frequency (not really),  modelRange-> ngfreq, nparam, npenet, ntthresh
       * TP QT alternative hypothesis is cL[3], looped for all of cL[2] and allele pairs, ndprime, ntheta
       * 
       */

      if (modelOptions->equilibrium != LINKAGE_EQUILIBRIUM) {
        pLocus1 = originalLocusList.ppLocusList[0];
        for (loc2 = 1; loc2 < originalLocusList.numLocus; loc2++) {
          pLocus2 = originalLocusList.ppLocusList[loc2];
          totalLoopsForDPrime += pow (modelRange->ndprime, (pLocus1->numOriginalAllele - 1) * (pLocus2->numOriginalAllele - 1));
        }
        // Divide by the iterations
        totalLoopsForDPrime /= (originalLocusList.numLocus - 1);
      } else
        totalLoopsForDPrime = 1;

      if (modelOptions->equilibrium == LINKAGE_EQUILIBRIUM)
        sprintf (analysisType, "%dTh*%d pair(s) of %dAL*%dGF*%dpv(%dLC) space for %d pedigree(s)\n"
            "Trait-to-marker Two-Point, ", modelRange->ntheta, (originalLocusList.numLocus - 1), modelRange->nalpha, modelRange->ngfreq, modelRange->npenet, modelRange->nlclass, pedigreeSet.numPedigree);
      else
        sprintf (analysisType, "%dTh*%dD' cases of %dAL*%dGF*%dp1*%dpv(%dLC)' space for %d pedigree(s)\n"
            "Trait-to-marker Two-Point, ", modelRange->ntheta, totalLoopsForDPrime, modelRange->nalpha, modelRange->ngfreq, modelRange->nparam, modelRange->npenet, modelRange->nlclass, pedigreeSet.numPedigree);
      if (modelType->trait == DT) {
        strcat (analysisType, "Dichotomous Trait, ");
        eCL[0] = (originalLocusList.numLocus - 1) * modelRange->ngfreq * modelRange->npenet;
        eCL[1] = eCL[0] * totalLoopsForDPrime * modelRange->ntheta;
      } else {  // TP not DT
        eCL[2] = (originalLocusList.numLocus - 1) * modelRange->ngfreq * modelRange->npenet * modelRange->nparam * modelRange->ntthresh;
        eCL[3] = eCL[2] * totalLoopsForDPrime * modelRange->ntheta;
        if (modelType->trait == QT) {    //QT
          strcat (analysisType, "Quantitative Trait, ");
        } else {        // TP not DT or QT, so CT
          strcat (analysisType, "Quantitative Trait w/Threshold, ");
        }
        if (modelType->distrib == QT_FUNCTION_T) {
          strcat (analysisType, "Student's T-Distribution, Linkage ");
        } else {        // not T-Dist
          if (modelType->distrib == QT_FUNCTION_CHI_SQUARE) {
            strcat (analysisType, "Chi-Square Distribution, Linkage ");
          } else {      // not T-Dist or Chi-Sq Dist, so Normal Dist
            strcat (analysisType, "Normal Distribution, Linkage ");
          }
        }
      }
      if (modelOptions->equilibrium == LINKAGE_EQUILIBRIUM)
        strcat (analysisType, "Equilibrium.");
      else
        strcat (analysisType, "Disequilibrium.");
    } else {    // not TP, so multipoint

      /* Pedigree likelihood calculation looping for MP is for trait, marker, then alternative hypothesis.
       * 
       * modelRange->ntloc is all evaluation locations on the chromosome since we're past adding for TM.
       * 
       * MP DT trait is cl[4], looped for 3(SA/SS)?, modelRange-> npenet, nlclass, ngfreq
       * MP QT/CT trait is cl[5], looped for 3(SA/SS)?, modelRange-> ngfreq, nparam (QT dist), npenet, ntthresh
       * MP marker is cl[6], looped for <= modelRange->ntloc
       * MP DT alt is cl[7], looped for modelRange->ntloc, locusList->numLocus, modelRange-> npenet, ngfreq
       * MP QT/CT alt is cl[8], looped for modelRange->ntloc,  locusList->numLocus, 
       * modelRange-> ngfreq, nparam, npenet, ntthresh
       * ...but locusList->numLocus is calculated, use modelRange->ntloc
       * 
       * polynomials for each pedigree incorporate alpha?, # MP markers used in analysis.
       * 
       */
      sprintf (analysisType, "%dTL of %dAL*%dGF*%dpv(%dLC) space for %d pedigree(s)\n"
          "Trait-to-marker, Sex-%s Multipoint (w/%d loci), ",
          modelRange->ntloc, modelRange->nalpha, modelRange->ngfreq, modelRange->npenet, modelRange->nlclass, pedigreeSet.numPedigree, modelOptions->mapFlag == SS ? "Specific" : "Averaged", modelType->numMarkers + originalLocusList.numTraitLocus);
      eCL[6] = modelRange->ntloc;
      if (modelType->trait == DT) {
        strcat (analysisType, "Dichotomous Trait.");
        eCL[4] = modelRange->npenet * modelRange->nlclass * modelRange->ngfreq;
        eCL[7] = modelRange->ntloc * modelRange->npenet * modelRange->ngfreq;
      } else {  // SA/SS multipoint, but not DT
        eCL[5] = modelRange->npenet * modelRange->ntthresh * modelRange->ngfreq * modelRange->nparam;
        eCL[8] = modelRange->ntloc * modelRange->npenet * modelRange->ngfreq * modelRange->ntthresh * modelRange->nparam;
        if (modelType->trait == QT) {
          strcat (analysisType, "Quantitative Trait, ");
        } else {        // SA/SS multipoint but not DT or QT, so CT
          strcat (analysisType, "Quantitative Trait w/Threshold, ");
        }
        if (modelType->distrib == QT_FUNCTION_T) {
          strcat (analysisType, "Student's T-Distribution.");
        } else {        // not T-Dist
          if (modelType->distrib == QT_FUNCTION_CHI_SQUARE) {
            strcat (analysisType, "Chi-Square Distribution.");
          } else {      // not T-Dist or Chi-Sq dist, so Normal Dist
            strcat (analysisType, "Normal Distribution.");
          }
        }
      }
    }
  }
  return (analysisType);
}
