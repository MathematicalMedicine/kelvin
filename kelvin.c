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

/* Some default global values. */
char markerfile[KMAXFILENAMELEN + 1] = "markers.dat";
char mapfile[KMAXFILENAMELEN + 1] = "mapfile.dat";
char pedfile[KMAXFILENAMELEN + 1] = "pedfile.dat";
char loopsfile[KMAXFILENAMELEN + 1] = "loops.dat";
char outfile[KMAXFILENAMELEN + 1] = "lods.out";
/* Model datastructures. modelOptions is defined in the pedigree library. */
ModelType modelType;
ModelRange modelRange;
ModelOptions modelOptions;
char sUnknownPedID[] = "0";

/* temporarily for GAW project */
double theta_lr_total[100];
int theta_lr_count[100];
double theta_alpha_lr_total[100];
double alpha[21];
typedef struct
{
  double LR;
  double alpha;
  double gfreq;
  int penIdx;
  double theta;
} MaxLR;
MaxLR maxLR[100];
MaxLR maxMaxLR;

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
  //double likelihood_null, likelihood_alternative;
  double log10_likelihood_null, log10_likelihood_alternative;
  double likelihood_ratio;
  double log10_likelihood_ratio;
  Locus *pLocus;
  Trait *pTrait;
  int pedIdx;
  double homoLR, hetLR;
  double max;
  double constraint;


  clock_t time0, time1;


  for (i = 1; i <= 20; i++)
    {
      alpha[i] = i * 0.05;
    }
  for (i = 0; i < 100; i++)
    maxLR[i].penIdx = -99;
  /* End GAW */

  /* Initialize the logging system. */
  logInit ();
  /* logSet(LOGGENOELIM, LOGDEBUG); */
  /* logSet(LOGPEELGRAPH, LOGFATAL); */
  logSet (LOGLIKELIHOOD, LOGWARNING);
  /* logSet(LOGLIKELIHOOD, LOGDEBUG); */
  /* logSet(LOGPARENTALPAIR, LOGDEBUG); */

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

  /* Parse the configuration file. */
  KASSERT (readConfigFile (configfile, &modelType, &modelRange, &modelOptions)
	   != ERROR, "Error in configuration file; aborting.\n");
  fprintf (stderr, "model: %d\n", modelType.trait);

  /* For now, reject all models we can't deal with. So use KASSERT to
   * check that we're looking at, e.g., twopoint dichotomous models
   * with direct eval (not polynomial --> add this to model?) and
   * give appropriate error message otherwise. */
  KASSERT (modelType.type == TP, "Only two-point analysis supported.\n");
  KASSERT (modelType.trait != CT, "Combined traits unsupported.\n");
  KASSERT (modelRange.nalleles == 2, "Only biallelic traits supported.\n");
  KASSERT (modelOptions.equilibrium == LE,
	   "Only linkage equilibrium supported.\n");

  fflush (stderr);
  fflush (stdout);


  /* Read in the map file. */
  read_mapfile (mapfile);

  /* Initialize the locus list and read in the marker file. */
  memset (&originalLocusList, 0, sizeof (originalLocusList));
  read_datafile (markerfile);
  fprintf (stderr, "Total number of markers in data: %d\n",
	   originalLocusList.numLocus - originalLocusList.numTraitLocus);
  fprintf (stderr, "Total number of trait locus in data: %d\n",
	   originalLocusList.numTraitLocus);
  fprintf (stderr, "%s\n",
	   (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) ? "LE" : "LD");

  modelOptions.sUnknownPersonID = sUnknownPedID;
  /* Initialize the pedigree set datastructure and read in the pedigrees. */
  memset (&pedigreeSet, 0, sizeof (PedigreeSet));
  read_pedfile (pedfile, &pedigreeSet);


  time0 = clock ();

  //Polynomial initialization must be after reading markerfile because markfile has the
  //flag of POLY or NOPOLY.  Also, it must be before initialize_loci
#ifndef NO_POLYNOMIAL
  if (modelOptions.polynomial == TRUE)
    {
      polynomialInitialization ();
      fprintf (stderr,
	       "!!!!!!!!!!!The Computation is done in polynomial mode!!!!!!!!!!!!!!!\n");
    }
#endif




  /* initialize loci by doing genotype elimination, set recoding */
  initialize_loci (&pedigreeSet);
  memset (theta_lr_count, 0, sizeof (theta_lr_count));
  memset (theta_lr_total, 0, sizeof (theta_lr_total));

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

  /* Initialize the connection to the infrastructure. */
  /* niceInit (); */

  /* Set up grid and invoke direct eval used as in Yungui's test
   * file. If we're doing 2 point, this is pretty straightforward. If
   * we're doing multipoint, this is the nflankloop behavior part. */
  if (modelType.type == TP)
    {
      /* Two point. */

      /* Loop over the penetrances, genefrequencies, thetas and call
         the likelihood calculation, storing each value optained to
         disk. */
      pLocus = originalLocusList.ppLocusList[0];
      //pTrait = pLocus->pTraitLocus->pTraits[0];
      pTrait = originalLocusList.ppLocusList[0]->pTraitLocus->pTraits[0];
      locusList.numLocus = 2;
      locusList.pLocusIndex = (int *) malloc (sizeof (int) *
					      locusList.numLocus);
      locusList.pPrevLocusDistance = (double *) malloc (sizeof (double) *
							locusList.numLocus);
      locusList.pNextLocusDistance = (double *) malloc (sizeof (double) *
							locusList.numLocus);

      locusList.pLocusIndex[0] = 0;
      locusList.pLocusIndex[1] = 1;
      locusList.pPrevLocusDistance[0] = -1;


      /* For dichotomous trait, 2 point:
         for every pedigree
         for every marker
         for every penetrance
         for every genefrequency triple
         for every theta
         compute and store the liklihood */

      /* Start GAW */
      if (pTrait->type == DICHOTOMOUS)
	{

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
		update_penetrance (&pedigreeSet, 0);
#else
	      update_penetrance (&pedigreeSet, 0);
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
		    update_loci (&pedigreeSet);
#else
		  update_loci (&pedigreeSet);
#endif
		  /* get the likelihood at 0.5 first */
		  locusList.pNextLocusDistance[0] = 0.5;
		  locusList.pPrevLocusDistance[1] = 0.5;
		  compute_likelihood (&pedigreeSet);
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
			  fprintf (stderr, "Liab %d penentrance %f %f %f\n",
				   liabIdx + 1, pen_DD, pen_Dd, pen_dd);
			}

		      exit (-1);
		    }
		  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++)
		    {
		      /* save the likelihood at null */
		      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
		      pedigreeSet.nullLikelihood[pedIdx] =
			pPedigree->likelihood;
		    }

		  log10_likelihood_null = pedigreeSet.log10Likelihood;
		  for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
		    {
		      theta = modelRange.theta[0][thetaInd];
		      locusList.pNextLocusDistance[0] = theta;
		      locusList.pPrevLocusDistance[1] = theta;
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
		      likelihood_ratio = pow (10.0, log10_likelihood_ratio);
		      /* caculating the HET */
		      for (j = 1; j <= 20; j++)
			{
			  hetLR = 1;
			  for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree;
			       pedIdx++)
			    {
			      pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
			      homoLR =
				pPedigree->likelihood /
				pedigreeSet.nullLikelihood[pedIdx];
			      hetLR *= alpha[j] * homoLR + (1 - alpha[j]);
			    }
			  pedigreeSet.hetLR[j] = hetLR;
			  pedigreeSet.log10HetLR[j] = log10 (hetLR);
			  theta_alpha_lr_total[thetaInd] += hetLR;
			  if (maxLR[thetaInd].penIdx < 0
			      || hetLR > maxLR[thetaInd].LR)
			    {
			      maxLR[thetaInd].LR = hetLR;
			      maxLR[thetaInd].alpha = alpha[j];
			      maxLR[thetaInd].gfreq = gfreq;
			      maxLR[thetaInd].penIdx = penIdx;
			    }
			}
		      /* add the result to the right placeholder */
		      theta_lr_total[thetaInd] += likelihood_ratio;
		      theta_lr_count[thetaInd]++;
		      //fprintf(stderr, "likelihood ratio: %e.\n", likelihood_ratio);

		    }

		}
	    }
	}
      else			/* should be QT or COMBINED */
	{
	  for (gfreqInd = 0; gfreqInd < modelRange.ngfreq; gfreqInd++)
	    {
	      gfreq = modelRange.gfreq[gfreqInd];
	      pLocus->pAlleleFrequency[0] = gfreq;
	      pLocus->pAlleleFrequency[1] = 1 - gfreq;

	      update_loci (&pedigreeSet);
	      /* this should be MEAN + SD */
	      for (paramIdx = 0; paramIdx < modelRange.nparam; paramIdx++)
		{
		  for (penIdx = 0; penIdx < modelRange.npenet; penIdx++)
		    {
		      breakFlag = FALSE;
		      for (liabIdx = 0; liabIdx < modelRange.nlclass;
			   liabIdx++)
			{
			  mean_DD = modelRange.penet[liabIdx][0][penIdx];
			  mean_Dd = modelRange.penet[liabIdx][1][penIdx];
			  mean_dd = modelRange.penet[liabIdx][2][penIdx];
			  SD_DD = modelRange.param[liabIdx][0][0][paramIdx];
			  SD_Dd = modelRange.param[liabIdx][1][0][paramIdx];
			  SD_dd = modelRange.param[liabIdx][2][0][paramIdx];

			  /* check against the hard coded constraint */
			  constraint =
			    pow (1.0 - gfreq,
				 2) * mean_dd * SD_dd + 2 * gfreq * (1 -
								     gfreq) *
			    mean_Dd * SD_Dd + pow (gfreq,
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
			}
		      if (breakFlag == TRUE)
			continue;
		      update_penetrance (&pedigreeSet, 0);

		      /* get the likelihood at 0.5 first */
		      locusList.pNextLocusDistance[0] = 0.5;
		      locusList.pPrevLocusDistance[1] = 0.5;
		      compute_likelihood (&pedigreeSet);
		      if (pedigreeSet.likelihood == 0.0 &&
			  pedigreeSet.log10Likelihood == -9999.99)
			{
			  fprintf (stderr, "Theta 0.5 has likelihood 0\n");
			  fprintf (stderr, "dgf=%f\n", gfreq);
			  for (liabIdx = 0; liabIdx < modelRange.nlclass;
			       liabIdx++)
			    {
			      pen_DD = modelRange.penet[0][liabIdx][penIdx];
			      pen_Dd = modelRange.penet[1][liabIdx][penIdx];
			      pen_dd = modelRange.penet[2][liabIdx][penIdx];
			      fprintf (stderr,
				       "Liab %d penentrance %f %f %f\n",
				       liabIdx + 1, pen_DD, pen_Dd, pen_dd);
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
		      for (thetaInd = 0; thetaInd < modelRange.ntheta;
			   thetaInd++)
			{
			  theta = modelRange.theta[0][thetaInd];
			  locusList.pNextLocusDistance[0] = theta;
			  locusList.pPrevLocusDistance[1] = theta;
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
			  for (j = 1; j <= 20; j++)
			    {
			      hetLR = 1;
			      for (pedIdx = 0;
				   pedIdx < pedigreeSet.numPedigree; pedIdx++)
				{
				  pPedigree =
				    pedigreeSet.ppPedigreeSet[pedIdx];
				  homoLR =
				    pPedigree->likelihood /
				    pedigreeSet.nullLikelihood[pedIdx];
				  hetLR *= alpha[j] * homoLR + (1 - alpha[j]);
				}
			      pedigreeSet.hetLR[j] = hetLR;
			      pedigreeSet.log10HetLR[j] = log10 (hetLR);
			      theta_alpha_lr_total[thetaInd] += hetLR;
			      if (maxLR[thetaInd].penIdx < 0
				  || hetLR > maxLR[thetaInd].LR)
				{
				  maxLR[thetaInd].LR = hetLR;
				  maxLR[thetaInd].alpha = alpha[j];
				  maxLR[thetaInd].gfreq = gfreq;
				  maxLR[thetaInd].penIdx = penIdx;
				  //fprintf(stderr, "maxLR update: %f (Theta %f)\n", hetLR, theta);
				}
			    }
			  /* add the result to the right placeholder */
			  theta_lr_total[thetaInd] += likelihood_ratio;
			  theta_lr_count[thetaInd]++;
			  //fprintf(stderr, "likelihood ratio: %e.\n", likelihood_ratio);

			}
		    }		/* end of mean loop */
		}
	    }
	}

      /* for GAW project only */
      fprintf (stderr, "#0  \"Average Homo LR\" \n");
      for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
	{
	  theta = modelRange.theta[0][thetaInd];
	  fprintf (stderr, "\t (%f,%f)  %f(%d)\n",
		   theta, theta,
		   theta_lr_total[thetaInd] / theta_lr_count[thetaInd],
		   theta_lr_count[thetaInd]);
	}
      fprintf (stderr, "-	Total 1234(1234)\n");

      fprintf (stderr, "#1 \"Average Het LR\" \n");
      for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
	{
	  theta = modelRange.theta[0][thetaInd];
	  fprintf (stderr, "\t (%f,%f)  %f(%d)\n",
		   theta, theta,
		   theta_alpha_lr_total[thetaInd] / (20.0 *
						     theta_lr_count
						     [thetaInd]),
		   20 * theta_lr_count[thetaInd]);
	}
      fprintf (stderr, "-	Total 1234(1234)\n");

      fprintf (stderr, "#2 Max Het LR \n");
      for (thetaInd = 0; thetaInd < modelRange.ntheta; thetaInd++)
	{
	  theta = modelRange.theta[0][thetaInd];
	  if (maxLR[thetaInd].LR > max)
	    {
	      max = maxLR[thetaInd].LR;
	      memcpy (&maxMaxLR, &maxLR[thetaInd], sizeof (maxMaxLR));
	      maxMaxLR.theta = theta;
	    }
	}
      fprintf (stderr, "MOD: %10.6f\n", log10 (max));
      fprintf (stderr, "Maximizing parameters:\n");
      fprintf (stderr,
	       "Theta    LR      HLOD  alpha gfreq penDD.1 penDd.1 pendd.1 penDD.2 penDd.2 pendd.2 ... \n");
      fprintf (stderr, "\t (%f,%f)  %10.6f  %8.4f %f %f ", maxMaxLR.theta,
	       maxMaxLR.theta, maxMaxLR.LR, log10 (maxMaxLR.LR),
	       maxMaxLR.alpha, maxMaxLR.gfreq);
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++)
	{
	  pen_DD = modelRange.penet[0][liabIdx][maxMaxLR.penIdx];
	  pen_Dd = modelRange.penet[1][liabIdx][maxMaxLR.penIdx];
	  pen_dd = modelRange.penet[2][liabIdx][maxMaxLR.penIdx];
	  fprintf (stderr, " %f %f %f ", pen_DD, pen_Dd, pen_dd);
	}
      fprintf (stderr, "\n");

    }

  /*
     fprintf(stderr, "Theta         LR      HLOD  alpha gfreq penDD.1 penDd.1 pendd.1 penDD.2 penDd.2 pendd.2 ... \n");
     for(thetaInd=0; thetaInd < modelRange.ntheta; thetaInd++) 
     {
     theta = modelRange.theta[0][thetaInd];
     fprintf(stderr, "\t (%f,%f)  %10.6f  %8.4f %f %f ", 
     theta, theta,
     maxLR[thetaInd].LR, log10(maxLR[thetaInd].LR),
     maxLR[thetaInd].alpha, maxLR[thetaInd].gfreq);
     for(liabIdx=0; liabIdx < modelRange.nlclass; liabIdx++) 
     {
     pen_DD = modelRange.penet[0][liabIdx][penIdx];
     pen_Dd = modelRange.penet[1][liabIdx][penIdx];
     pen_dd = modelRange.penet[2][liabIdx][penIdx];
     fprintf(stderr, " %f %f %f ", pen_DD, pen_Dd, pen_dd);
     }
     fprintf(stderr, "\n");
     }
   */

  /* End GAW */

  /* For quantitative trait, 2 point:
     for every pedigree
     for every marker
     for every mean/sd tuple?
     for every theta
     compute and store the liklihood */

#if FALSE
  else
    {
      /* Multi point. */

      /* For dichotomous trait, multi point:
         for every pedigree
         for every marker
         for every penetrance
         for every genefrequency triple
         for every location (marker set is computed)
         compute and store the liklihood */


      /* For quantitative trait, multi point:
         for every pedigree
         for every marker
         for every mean/sd tuple?
         for every location (marker set is computed)
         compute and store the liklihood */

    }
#endif

  time1 = clock ();
  fprintf (stderr, "Computation time:  %fs ",
	   (double) (time1 - time0) / CLOCKS_PER_SEC);


  exit (TRUE);
}
