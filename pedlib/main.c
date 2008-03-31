/**********************************************************************
 * Copyright 2007, Columbus Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "pedlib.h"
#include "utils.h"
#include "polynomial.h"

ModelOptions modelOptions;
char sUnknownPedID[] = "0";

time_t startTime;
time_t endTime;

int
main (int argc, char *argv[])
{
  PedigreeSet pedigreeSet;
  double theta = 0.2;
  double likelihood_null;
  double likelihood_alternative;
  double log10_null, log10_alternative;
  int i;

  if (argc < 4)
    {
      fprintf (stderr, "%s <postped> <mapfile> <markerfile>.\n", argv[0]);
      exit (-1);
    }
  if (argc > 4)
    theta = atof (argv[4]);

  /* initialize logging */
  logInit ();

  /* set up logging type levels 
   * all pedfile log messages with at least DEBUG level will be output */
  //logSet(LOGPEDFILE, LOGDEBUG);
  //logSet(LOGSETRECODING, LOGDEBUG);
  //logSet(LOGGENOELIM, LOGDEBUG);
  logSet (LOGGENOELIM, LOGFATAL);
  //logSet(LOGPARENTALPAIR, LOGDEBUG);
  logSet (LOGPEELGRAPH, LOGFATAL);
  //logSet(LOGLIKELIHOOD, LOGDEBUG);
  logSet (LOGMEMORY, LOGDEBUG);

  modelOptions.sUnknownPersonID = sUnknownPedID;
  modelOptions.equilibrium = LINKAGE_EQUILIBRIUM;

  /* read map file */
  read_mapfile (argv[2]);
  fprintf (stderr, "Total number of markers in map: %d\n", map.count);
  /* read data file */
  memset (&originalLocusList, 0, sizeof (originalLocusList));
  read_datafile (argv[3]);
  fprintf (stderr, "Total number of markers in data: %d\n",
	   originalLocusList.numLocus - originalLocusList.numTraitLocus);
  fprintf (stderr, "Total number of trait locus in data: %d\n",
	   originalLocusList.numTraitLocus);
  fprintf (stderr, "%s\n",
	   (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) ? "LE" : "LD");
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomialFlag == TRUE)
    {
      polynomialInitialization ();
    }
#endif

  /* initialize the space to make sure pointers are initialized to NULL */
  memset (&pedigreeSet, 0, sizeof (PedigreeSet));

  /* read postmakeped pedigree file and load pedigree related data structure */
  read_pedfile (argv[1], &pedigreeSet);

  /* initialize loci by doing genotype elimination, set recoding */
  initialize_loci (&pedigreeSet);

  /* set a subset of the locus list */
//  locusList.numLocus = 3;
  locusList.numLocus = 2;
  locusList.pLocusIndex = (int *) malloc (sizeof (int) * locusList.numLocus);
  locusList.pPrevLocusDistance = (double *) malloc (sizeof (double) *
						    locusList.numLocus);
  locusList.pNextLocusDistance = (double *) malloc (sizeof (double) *
						    locusList.numLocus);
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomialFlag == TRUE)
    {
      locusList.pPrevLocusDistance = (double *) malloc (sizeof (Polynomial) *
							locusList.numLocus);
      locusList.pNextLocusDistance = (double *) malloc (sizeof (Polynomial) *
							locusList.numLocus);
    }

#endif
  locusList.pLocusIndex[0] = 0;
  locusList.pLocusIndex[1] = 1;
  locusList.pLocusIndex[2] = 2;
  locusList.pPrevLocusDistance[0] = -1;
  locusList.pNextLocusDistance[0] = 0.5;
  locusList.pPrevLocusDistance[1] = 0.5;
  //locusList.pNextLocusDistance[1]= 0.0421;
  //locusList.pPrevLocusDistance[2]= 0.0421;
  locusList.pNextLocusDistance[locusList.numLocus - 1] = -1;

  printf ("Computing likelihood.\n");
  /* now compute likelihood */

  compute_likelihood (&pedigreeSet);


#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomialFlag == TRUE)
    {
      likelihood_null = evaluateValue (pedigreeSet.likelihoodPolynomial);
      locusList.pNextLocusDistance[0] = theta;
      locusList.pPrevLocusDistance[1] = theta;
      compute_likelihood (&pedigreeSet);
      likelihood_alternative =
	evaluateValue (pedigreeSet.likelihoodPolynomial);
    }
  else
    {
      likelihood_null = pedigreeSet.likelihood;
      log10_null = pedigreeSet.log10Likelihood;
      for (i = 1; i < 100000; i++)
	{
	  locusList.pNextLocusDistance[0] = theta;
	  locusList.pPrevLocusDistance[1] = theta;
	  compute_likelihood (&pedigreeSet);
	  likelihood_alternative = pedigreeSet.likelihood;
	  log10_alternative = pedigreeSet.log10Likelihood;
	}
    }
#else
  likelihood_null = pedigreeSet.likelihood;
  log10_null = pedigreeSet.log10Likelihood;
  time (&startTime);
  for (i = 1; i <= 10000; i++)
    {
      fprintf (stderr, "Iteration %d\n", i);
      locusList.pNextLocusDistance[0] = theta;
      locusList.pPrevLocusDistance[1] = theta;
      compute_likelihood (&pedigreeSet);
      likelihood_alternative = pedigreeSet.likelihood;
      log10_alternative = pedigreeSet.log10Likelihood;
    }
#endif
  printf ("Likelihood at 0.5: %e (log10: %e)\n",
	  likelihood_null, log10 (likelihood_null));
  printf ("Likelihood at %f: %e (log10: %e)\n", theta,
	  likelihood_alternative, log10 (likelihood_alternative));
  printf ("LOD score at %f: %f\n", theta, log10_alternative - log10_null);
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomialFlag == TRUE)
    {
      fprintf (stderr, "Computations are done with polynomials\n");
    }
  else
    {
      fprintf (stderr, "Computations are done without polynomials\n");
    }
#else
  fprintf (stderr, "Computations are done without polynomials\n");
#endif

  time (&endTime);
  fprintf (stderr, "Lapsed time: %llu sec\n",
	   (unsigned long long int) difftime (endTime, startTime));
  return 0;
}
