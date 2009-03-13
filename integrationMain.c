  /* total_dim is the number of all parameters in the 3-layer scheme
          s->dim in dcuhre.c is the number of parameters in the middle layer alone*/
  total_dim = 2;		// alpha gf
  total_dim += 3 * modelRange.nlclass;	//DD Dd dd
  if(modelOptions.imprintingFlag)
    total_dim += modelRange.nlclass;	//dD


  if (modelType.type == TP) {
    total_dim += 1;		// theta;
    if(modelOptions.mapFlag == SS)
      total_dim += 1;		// theta sex-specific case;

    if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
      total_dim += 1;		// dprime
    }
  }

  if (modelType.trait != DT) {
    if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
      total_dim += 3 * modelRange.nlclass;	//SD_DD SD_Dd SD_dd

      if(modelOptions.imprintingFlag)
        total_dim += modelRange.nlclass;   // SD_dD
    }
    if (modelType.trait == CT) {
      total_dim += modelRange.nlclass;
    }
  }

  memset (&dk_globalmax, 0, sizeof (st_DKMaxModel));
  memset (&dk_dprime0max, 0, sizeof (st_DKMaxModel));
  memset (&dk_theta0max, 0, sizeof (st_DKMaxModel));
  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    /* Assumes that dkelvin can only handle a single D' */
    dk_globalmax.dprime = (double *) calloc (1, sizeof (double));
    dk_dprime0max.dprime = (double *) calloc (1, sizeof (double));
    dk_theta0max.dprime = (double *) calloc (1, sizeof (double));
    KASSERT ((dk_globalmax.dprime != NULL) && (dk_dprime0max.dprime != NULL) &&
	     (dk_theta0max.dprime != NULL), "malloc failed");
  }
  dk_globalmax.pen = calloc (modelRange.nlclass, sizeof (st_DKMaxModelPenVector));
  dk_dprime0max.pen = calloc (modelRange.nlclass, sizeof (st_DKMaxModelPenVector));
  dk_theta0max.pen = calloc (modelRange.nlclass, sizeof (st_DKMaxModelPenVector));
  KASSERT ((dk_globalmax.pen != NULL) && (dk_dprime0max.pen != NULL) &&
	   (dk_theta0max.pen != NULL), "malloc failed");

  fprintf (stderr, "total dim =%d\n", total_dim);

  if (modelType.trait != DT) {
    /* Setting ranges for each variables. Default is [0,1] */
    k = 1;
    for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	xl[k] = xl[k + 1] = xl[k + 2] = -3;
	xu[k] = xu[k + 1] = xu[k + 2] = 3;
        if(modelOptions.imprintingFlag){
          xl[k+3]= -3;
          xu[k+3]= 3;
	}
      } else {
	xl[k] =    modelRange.penetLimits[0][0];//0.1;
        xl[k + 1] =modelRange.penetLimits[1][0];//0.1;
        xl[k + 2] =modelRange.penetLimits[3][0];//0.1;
	xu[k] =    modelRange.penetLimits[0][1];//30.0;
        xu[k + 1] =modelRange.penetLimits[1][1];//30.0;
        xu[k + 2] =modelRange.penetLimits[3][1];//30.0; //10.0; //23.0; //4.0;//30;

        if(modelOptions.imprintingFlag){
          xl[k+2]=  modelRange.penetLimits[2][0];
          xu[k+2]= modelRange.penetLimits[2][1];
          xl[k+3]=  modelRange.penetLimits[3][0];
          xu[k+3]= modelRange.penetLimits[3][1];
	}

      }
      volume_region *= (xu[k] - xl[k]);
      volume_region *= (xu[k + 1] - xl[k + 1]);
      volume_region *= (xu[k + 2] - xl[k + 2]);
      if(modelOptions.imprintingFlag){
        volume_region *= (xu[k + 3] - xl[k + 3]);
        k++;
      }
      k += 3;

      if (modelType.distrib != QT_FUNCTION_CHI_SQUARE) {
	xl[k] = xl[k + 1] = xl[k + 2] = 0.3;
	xu[k] = xu[k + 1] = xu[k + 2] = 1.0;//3.0;
        if(modelOptions.imprintingFlag){
          xl[k+3]= 0.3;
          xu[k+3]= 1.0;
	}
	volume_region *= (xu[k] - xl[k]);
	volume_region *= (xu[k + 1] - xl[k + 1]);
	volume_region *= (xu[k + 2] - xl[k + 2]);
        if(modelOptions.imprintingFlag){
          volume_region *= (xu[k + 3] - xl[k + 3]);
          k++;
        }
	k += 3;
      }
      if (modelType.trait == CT) {
	xl[k] = modelRange.tthresh[liabIdx][0];//0.3;
	xu[k] = modelRange.tthresh[liabIdx][modelRange.ntthresh -1];// 23.0;
	volume_region *= (xu[k] - xl[k]);
	k++;
	//   fprintf(stderr, " in CT\n ");

      }
    }  // retangular volume region is calculated and stored in volume_region

    fprintf (stderr,"The number of dimension for calculation of BR is %d\n",k);
  }

  /* only for multipoint - we don't handle LD under multipoint yet */
  /* DCUHRE do now use likelihoodDT or likelihoodQT to store null likelihoods */

  /* find out the max we need to allocate */
  /* after genotype lists have been built, we want to pre-allocate parental pair work space
   * it used to be done dynamically, but it's too costly 
   * parental pair is based on each nuclear family */
  stat_parental_pair_workspace (&pedigreeSet);

  /* after genotype elimination and tracking the max work space needed 
   * for constructing parental pair */
  allocate_parental_pair_workspace (&parentalPairSpace,
				    modelType.numMarkers + 1);

  /* conditional likelihood storage space for each individual */
  allocate_likelihood_space (&pedigreeSet, modelType.numMarkers + 1);

  if (modelType.trait == DT)
    fprintf (stderr, "Dichotomous Trait & ");
  else if (modelType.trait == QT)
    fprintf (stderr, "Quantitative Trait without threshold & ");
  else
    fprintf (stderr, "Quantitative Trait with threshold & ");

  /* assume the trait locus is the first one in the list */
  traitLocus = 0;
  pLocus = originalLocusList.ppLocusList[traitLocus];
  pTraitLocus = originalLocusList.ppLocusList[traitLocus]->pTraitLocus;
  pTrait = pTraitLocus->pTraits[traitLocus];
  if (modelType.type == TP) {
    /* Two point. */
    if (originalLocusList.pLDLoci == NULL) {
      originalLocusList.pLDLoci = (LDLoci *) malloc (sizeof (LDLoci));
      memset (originalLocusList.pLDLoci, 0, sizeof (LDLoci));
    }
    pLDLoci = &originalLocusList.pLDLoci[0];
    originalLocusList.numLDLoci = 1;

    if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) {
      /* fake some LD information to simplify looping */
      pLDLoci->numAllele1 = 2;
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

    locusList = &savedLocusList;
    savedLocusList.numLocus = 2;
    savedLocusList.pLocusIndex =(int *) malloc (sizeof (int) * savedLocusList.numLocus);

    for (i = 0; i < 3; i++) {
      savedLocusList.pPrevLocusDistance[i] =(double *) malloc (sizeof (double) * savedLocusList.numLocus);
      savedLocusList.pNextLocusDistance[i] =(double *) malloc (sizeof (double) * savedLocusList.numLocus);
      savedLocusList.pPrevLocusDistance[i][0] = -1;
      savedLocusList.pNextLocusDistance[i][1] = -1;
    }

    if (modelOptions.polynomial == TRUE) {
      /* populate the matrix */
      status = populate_xmission_matrix (xmissionMatrix, totalLoci, initialProbAddr,	/* probability */
					 initialProbAddr2,	/* probability */
					 initialHetProbAddr, 0,	/* cell index */
					 -1,	/* last het locus */
					 -1,	/* last  pattern (P-1 or M-2) */
					 0);	/* current locus - start with 0 */

      holdAllPolys ();
    }
    //total_count = modelRange.npenet * modelRange.ngfreq * modelRange.nalpha;

    if (modelOptions.markerAnalysis == FALSE) {
      savedLocusList.traitLocusIndex = 0;
      savedLocusList.traitOrigLocus = 0;
    } else {
      savedLocusList.traitLocusIndex = -1;
      savedLocusList.traitOrigLocus = -1;
    }

    for (loc1 = 0; loc1 < originalLocusList.numLocus - 1; loc1++) {

      savedLocusList.pLocusIndex[0] = loc1;
      pLocus1 = originalLocusList.ppLocusList[loc1];
      if (modelOptions.markerAnalysis != FALSE && pLocus1->locusType != LOCUS_TYPE_MARKER)
	continue;

      for (loc2 = loc1 + 1; loc2 < originalLocusList.numLocus; loc2++) {

	maximum_function_value = 0.0;  // global max 
        maximum_dprime0_value = 0.0;   // max when D' == 0
        maximum_theta0_value = 0.0;    // max when Theta == 0

	pLocus2 = originalLocusList.ppLocusList[loc2];
	if (pLocus2->locusType != LOCUS_TYPE_MARKER)
	  continue;
	savedLocusList.pLocusIndex[1] = loc2;

#ifndef SIMPLEPROGRESS
	if (modelOptions.markerAnalysis == MM)
	  fprintf (stdout, "Starting w/loci %s(%d alleles) and %s(%d alleles\n", 
		   pLocus1->sName, pLocus1->numOriginalAllele, pLocus2->sName, pLocus2->numOriginalAllele);
	else
	  fprintf (stdout, "Starting w/loci %s(%d alleles) and %s(%d alleles) (%d of %d pairs)\n",
		   pLocus1->sName, pLocus1->numOriginalAllele, pLocus2->sName, pLocus2->numOriginalAllele,
		   loc2, originalLocusList.numLocus - 1);
#endif

	/* find out number of alleles this marker locus has *//* Check if this is okay with DCUHRE  ???????????? */
	if (modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM) {
	  /* get the LD parameters */
	  pLambdaCell =findLambdas (&modelRange, pLocus1->numOriginalAllele,pLocus2->numOriginalAllele);
	  reallocate_LD_loci (pLDLoci, pLocus1->numOriginalAllele, pLocus2->numOriginalAllele);


	  // Create these variables ahead of likelihood polynomial build in hopes of preventing in-build creation.

	  if (modelOptions.polynomial == TRUE) {
	    char vName[128];
	    int a0, a1;
	    for (a0 = 0; a0 < pLocus1->numOriginalAllele; a0++) {
	      for (a1 = 0; a1 < pLocus2->numOriginalAllele; a1++) {
		sprintf (vName, "ppHaploFreq_lA%d_rA%d", a0, a1);
		variableExp (&pLDLoci->ppHaploFreq[a0][a1], NULL, 'D', vName);
	      }
	    }
	  }


	  pLDLoci->locus1 = loc1;
	  pLDLoci->locus2 = loc2;
	  pLDLoci->numAllele1 = pLocus1->numOriginalAllele;
	  pLDLoci->numAllele2 = pLocus2->numOriginalAllele;
	  if (pLocus1->numOriginalAllele == 2 && pLocus2->numOriginalAllele == 2)
	    R_square_flag = TRUE;
	  else
	    R_square_flag = FALSE;
	}

	loopMarkerFreqFlag = 0;
	if (modelRange.nafreq >= 2 && modelOptions.equilibrium == LINKAGE_DISEQUILIBRIUM && pLocus2->numOriginalAllele == 2) {
	  loopMarkerFreqFlag = 1;
	} else if (modelRange.nafreq == 0) {
	  /* add a fake one to facilitate loops and other handlings */
	  addAlleleFreq (&modelRange, pLocus2->pAlleleFrequency[0]);
	} else {
	  modelRange.nafreq = 1;
	  modelRange.afreq[0] = pLocus2->pAlleleFrequency[0];
	}

	/* we will force marker allele frequency loop to execute at least once */
	for (mkrFreqIdx = 0;  mkrFreqIdx == 0 || mkrFreqIdx < modelRange.nafreq;  mkrFreqIdx++) {
	  mkrFreq = pLocus2->pAlleleFrequency[0];
	  /* we should only loop over marker allele frequency under twopoint
	   * and when markers are SNPs (only have two alleles) */
	  if (loopMarkerFreqFlag) {
	    mkrFreq = modelRange.afreq[mkrFreqIdx];
	    /* update the locus */
	    pLocus2->pAlleleFrequency[0] = mkrFreq;
	    pLocus2->pAlleleFrequency[1] = 1 - mkrFreq;
	    if (modelOptions.polynomial == TRUE);
	    else
	      update_locus (&pedigreeSet, loc2);
	  }
	  //	  fprintf (stderr, "mkrFreq =%d  model nafreq= %d \n",mkrFreqIdx, modelRange.nafreq);


	  if (1 && modelOptions.markerAnalysis == FALSE) {

	    if (modelOptions.polynomial == TRUE);
	    else
	      update_locus (&pedigreeSet, loc1);
	  }

	  /* clear Dprime combination impossible flag */
	  memset (pLambdaCell->impossibleFlag, 0, sizeof (int) * pLambdaCell->ndprime);
	  /* set up haplotype frequencies */
	  for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
	    if (isDPrime0(pLambdaCell->lambda[dprimeIdx], pLambdaCell->m, pLambdaCell->n))
	      dprime0Idx = dprimeIdx;
	    status =
	      setup_LD_haplotype_freq (pLDLoci, pLambdaCell, dprimeIdx);
	    if (status < 0) {
	      pLambdaCell->impossibleFlag[dprimeIdx] = 1;
	    }
	  }

	  /* for each D prime and theta, print out average and maximizing model information - MOD */
	  if (modelOptions.markerAnalysis == FALSE)
	    dk_write2ptBRHeader ();

          /* analysis specific statistic initialization*/
          if(modelOptions.mapFlag == SA){
	    le_small_theta = 0.0;
	    le_big_theta = 0.0;
	    ld_small_theta = 0.0;
	    ld_big_theta = 0.0;
	    ld_unlinked = 0.0;
            le_unlinked=0.0;

            num_BR= num_sample_Dp_theta;  // currently 141

	  }else{      ///  This is for sec-specific analysis in four regions
            thetaSMSF= 0.0;  // 0< thetaM <0.05  0< thetaF <0.05  
            thetaBMSF= 0.0;  // 0.05< thetaM <0.5  0< thetaF <0.05  
            thetaSMBF= 0.0;  // 0< thetaM <0.05  0.05< thetaF <0.5  
            thetaBMBF= 0.0;  // 0.05< thetaM <0.05  0.05< thetaF <0.5 
            
            num_BR = num_sample_SS_theta;  // currenlty 260
	  }

	  /*The main loop to Calculate BR(theta, dprime) or BR(thetaM, thetaF)*/
   	  for (i = 0; i < num_BR; i++) {    /* num_BR = 141 for Sex-Average Analysis
					       = 260 for Sex-Specific Analysis */
            if(modelOptions.mapFlag == SA){ 
	      fixed_dprime = dcuhre2[i][0];
	      fixed_theta = dcuhre2[i][1];
	    }else{
              fixed_thetaM = thetaSS[i][0];
              fixed_thetaF = thetaSS[i][1];
	    }
   	    integral = 0.0;
	    abserr = 0.0;
	      //fprintf (stderr, "i=%d Dprime=%f theta=%f \n", i, fixed_dprime, fixed_theta);
	   
	    if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {  // checking Dprime
	      for (dprimeIdx = 0; dprimeIdx < pLambdaCell->ndprime; dprimeIdx++) {
		if (fabs (pLambdaCell->lambda[dprimeIdx][0][0] - fixed_dprime) < 0.0001) {
		  // fprintf(stderr,"dprimeIdx =%d with %15.13f which is matching wit fixed_dprime\n",dprimeIdx,pLambdaCell->lambda[dprimeIdx][0][0]);
		  break;
		}
	      }
	      if (dprimeIdx == pLambdaCell->ndprime) {
		fprintf (stderr, "dprimeIdx is %d\n", dprimeIdx);
		exit (0);
	      }
	    }
	    num_out_constraint = 0;

            /* Call DCUHRE  Domain information is stored in global variables,  xl an xu*/
	    kelvin_dcuhre_integrate (&integral, &abserr, volume_region);
	    
	    dcuhre2[i][3] = integral;


            /* Dk specific results*/
	    if (fpDK != NULL) {
	      if(modelOptions.imprintingFlag){
		fprintf(fpDK,"%d %6.4f %6.4f %6d %8.4f %8.4f %8.4f\n",i, fixed_thetaM,fixed_thetaF, s->total_neval, integral, abserr,log10 (localmax_value));
	      }else{
		fprintf(fpDK,"%d %6.4f %6.4f %6d %8.4f %8.4f %8.4f\n",i, fixed_dprime,fixed_theta, s->total_neval, integral, abserr,log10 (localmax_value));
	      }
	      fflush (fpDK);
	    }
       
            R_square = 0.0;// tp_result[dprimeIdx][thetaInd][modelRange.nafreq].R_square;

	    if (modelOptions.markerAnalysis == FALSE)
	      dk_write2ptBRData (integral);

	    if (maximum_function_value < localmax_value) {
	      maximum_function_value = localmax_value;
	      if (modelOptions.mapFlag == SA) {
		maxima_x[0] = dk_globalmax.dprime[0] = fixed_dprime;
		maxima_x[1] = dk_globalmax.theta[0] = dk_globalmax.theta[1] = fixed_theta;
	      } else {
		dk_globalmax.dprime[0] = 0.0;
		maxima_x[0] = dk_globalmax.theta[0] = fixed_thetaM;
		maxima_x[1] = dk_globalmax.theta[1] = fixed_thetaF;
              }
	      dk_copyMaxModel (localmax_x, &dk_globalmax);
	      memcpy (&(maxima_x[2]), localmax_x, sizeof (double) * 18);
	    }

	    if ((modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) &&
		(maximum_dprime0_value < localmax_value) && (fabs (fixed_dprime) < 1e-9)) {
	      maximum_dprime0_value = localmax_value;
	      dk_dprime0max.dprime[0] = fixed_dprime;
	      dk_dprime0max.theta[0] = dk_dprime0max.theta[1] = fixed_theta;
	      dk_copyMaxModel (localmax_x, &dk_dprime0max);
	    }
	    
	    if (maximum_theta0_value < localmax_value) {
	      if ((modelOptions.mapFlag == SA) && (fabs (fixed_theta) < 1e-9)) {
		maximum_theta0_value = localmax_value;
 		dk_theta0max.dprime[0] = fixed_dprime;
		dk_theta0max.theta[0] = dk_theta0max.theta[1] = fixed_theta;
		dk_copyMaxModel (localmax_x, &dk_theta0max);
	      } else if ((modelOptions.mapFlag == SS) && (fabs (fixed_thetaM) < 1e-9) &&
			 (fabs (fixed_thetaF) < 1e-9)) {
		maximum_theta0_value = localmax_value;
		dk_theta0max.dprime[0] = 0.0;
		dk_theta0max.theta[0] = fixed_thetaM;
		dk_theta0max.theta[1] = fixed_thetaF;
		dk_copyMaxModel (localmax_x, &dk_theta0max);
	      }
	    }

	    //fprintf (stderr, "tp result %f %f is %13.10f   \n",  fixed_theta, fixed_dprime, integral);

            if(modelOptions.mapFlag == SA){
            
	      if (i < 5) {
	        le_small_theta += integral * dcuhre2[i][2];
	      } else if (i < 10) {
	        le_big_theta += integral * dcuhre2[i][2];
	      } else if (i < 140){
	        if (fixed_theta < modelOptions.thetaCutoff[0]) {
		  ld_small_theta += integral * dcuhre2[i][2];
	        } else {
		  ld_big_theta += integral * dcuhre2[i][2];
	        }
	      } else {
                le_unlinked += integral * dcuhre2[i][2];
	      }
	      if ((modelOptions.equilibrium == LINKAGE_EQUILIBRIUM) && (i == 9)) {
	        //fprintf (stderr,"End of LE case\n");
	        i = 141;
	      }
	    }else{
              if(i<65){
                thetaSMSF += integral *thetaSS[i][2];  // 0< thetaM <0.05  0< thetaF <0.05  
	      }else if(i<130){
                thetaBMSF += integral *thetaSS[i][2];  // 0.05< thetaM <0.5  0< thetaF <0.05  
	      }else if(i<195){
                thetaSMBF += integral *thetaSS[i][2];  // 0< thetaM <0.05  0.05< thetaF <0.5  
	      }else {
                thetaBMBF += integral *thetaSS[i][2];  // 0.05< thetaM <0.05  0.05< thetaF <0.5 
	      }
	    }

	  }			/* end of for to calculate BR(theta, dprime) or BR(thetaM, thetaF)*/
	  dk_write2ptMODFile (maximum_function_value, &dk_globalmax);
	  dk_writeMAXHeader ();
	  dk_writeMAXData ("MOD(Overall)", maximum_function_value, &dk_globalmax);
	  dk_writeMAXData ("MOD(Theta==0)", maximum_theta0_value, &dk_theta0max);
	  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM)
	    dk_writeMAXData ("MOD(D'==0)", maximum_dprime0_value, &dk_dprime0max);

	  /*Calculate ppl, ppld and ldppl */
          if(modelOptions.mapFlag == SA){
	    ppl =  modelOptions.thetaWeight * le_small_theta + (1 -  modelOptions.thetaWeight) * le_big_theta;
	    ppl = ppl / (ppl + (1 - modelOptions.prior) / modelOptions.prior);
	  }else{
            ppl = modelOptions.thetaWeight* thetaSMSF;
            ppl += (1-modelOptions.thetaWeight)* 0.09/0.99 * thetaBMSF;
            ppl += (1-modelOptions.thetaWeight)* 0.09/0.99 * thetaSMBF;
            ppl += (1-modelOptions.thetaWeight)* 0.81/0.99 * thetaBMBF;
	    ppl = ppl / (ppl + (1 - modelOptions.prior) / modelOptions.prior);
	  }
          if (modelOptions.markerAnalysis != FALSE) {
            fprintf (fpPPL, "%d %s %.4f %s %.4f %.*f ",
	     pLocus2->pMapUnit->chromosome, pLocus1->sName, pLocus1->pMapUnit->mapPos[SEX_AVERAGED],
	     pLocus2->sName, pLocus2->pMapUnit->mapPos[SEX_AVERAGED],
	     ppl >= .025 ? 2 : 3, KROUND (ppl));
          } else {
            fprintf (fpPPL, "%d %s %s %.4f %.*f ",
	     pLocus2->pMapUnit->chromosome, pLocus1->sName,
	     pLocus2->sName, pLocus2->pMapUnit->mapPos[SEX_AVERAGED],
	     ppl >= .025 ? 2 : 3, KROUND (ppl));
          }


          /* output LD-PPL now if needed */
	  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
	    ldppl =0.019*(0.021*ld_small_theta+ 0.979*le_small_theta)+ 0.001*(0.011*ld_big_theta+ 0.9989*le_big_theta);

	    ldppl =
                ldppl / (ldppl + 0.98*le_unlinked);

	    ppld = 0.019*0.021*ld_small_theta+0.001*0.0011*ld_big_theta;
            ppld = ppld/(ppld+ 0.019*0.979*le_small_theta +0.001*0.9989*le_big_theta+ 0.98*le_unlinked);
            //ppld= pow(0.02*(ld_small_theta *0.95 +ld_big_theta * 0.05), 0.3);
            //ppld *=  0.0004/0.020392;
            //ppld = ppld/(ppld + 1-0.0004/0.020392 * 0.02);

	    ppldGl=0.019*0.021*ld_small_theta+0.001*0.0011*ld_big_theta;
            ppldGl = ppldGl/(ppldGl + 0.019*0.979*le_small_theta +0.001*0.9989*le_big_theta);

            fprintf (fpPPL, "%.*f ", ldppl >= .025 ? 2 : 4, KROUND (ldppl));
            fprintf (fpPPL, "%.*f ", ppldGl >= .025 ? 2 : 4, KROUND (ppldGl));
            fprintf (fpPPL, "%.*f ", ppld >= .025 ? 2 : 4, KROUND (ppld));
	  }
          fprintf (fpPPL, "\n");
          fflush (fpPPL);

	  if (fpDK != NULL) {
	    fprintf (fpDK, "Global max %8.4f %6.4f %6.4f %6.4f %6.4f ",
		     log10 (maximum_function_value), maxima_x[0],
		     maxima_x[1], maxima_x[3], maxima_x[2]);
	    for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	      fprintf (fpDK, "%6.4f %6.4f %6.4f ",
		       maxima_x[liabIdx * 3 + 4],
		       maxima_x[liabIdx * 3 + 5], maxima_x[liabIdx * 3 + 6]);
	    }
	    fprintf (fpDK, "\n");
	    fflush (fpDK);
	  }

	  /* only loop marker allele frequencies when doing LD */
	  if (modelOptions.equilibrium == LINKAGE_EQUILIBRIUM)
	    break;
	  /* we can only do SNPs when looping over marker allele frequency */
	  if (pLocus2->numOriginalAllele > 2)
	    break;
	}
	/* end of marker allele frequency looping */
	/* need to clear polynomial */

	if (modelOptions.polynomial) {
	  pedigreeSetPolynomialClearance (&pedigreeSet);
	}


	if (modelOptions.markerAnalysis == ADJACENTMARKER)
	  loc2 = originalLocusList.numLocus;
      }				/* end of looping second locus - loc2 */
      /* if we are doing trait marker, then we are done */
      /* Used to read: modelOptions.markerToMarker != TRUE which
         is the same as markerAnalysis == FALSE as long as the old
         markerToMarker and adjacentMarker flags were truly
         orthogonal. Otherwise, it should be markerAnalysis !=
         ADJACENTMARKER. */
      if (modelOptions.markerAnalysis == FALSE)
	loc1 = originalLocusList.numLocus;
    }				/* end of looping first locus - loc1 */
  }
  /* end of two point */
  else {
    /* multipoint */
    /* marker set locus list for each position */
    markerLocusList.maxNumLocus = modelType.numMarkers;
    markerLocusList.numLocus = modelType.numMarkers;
    markerLocusList.traitOrigLocus = -1;
    markerLocusList.traitLocusIndex = -1;
    markerLocusList.pLocusIndex =
      (int *) calloc (markerLocusList.maxNumLocus, sizeof (int));
    for (k = 0; k < 3; k++) {
      markerLocusList.pPrevLocusDistance[k] =
	(double *) calloc (markerLocusList.maxNumLocus, sizeof (double));
      markerLocusList.pNextLocusDistance[k] =
	(double *) calloc (markerLocusList.maxNumLocus, sizeof (double));
    }

    /* assuming we always have trait in the analysis - this may not be true 
     * need to add code to process marker to marker analysis under multipoin
     */
    savedLocusList.numLocus = modelType.numMarkers + 1;
    savedLocusList.maxNumLocus = modelType.numMarkers + 1;
    savedLocusList.pLocusIndex =
      (int *) calloc (savedLocusList.maxNumLocus, sizeof (int));
    for (k = 0; k < 3; k++) {
      savedLocusList.pPrevLocusDistance[k] =
	(double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
      savedLocusList.pNextLocusDistance[k] =
	(double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
    }

    /* calculate the trait likelihood independent of the trait position */
    traitLocusList.numLocus = 1;
    traitLocusList.maxNumLocus = 1;
    traitLocusList.traitLocusIndex = 0;
    traitLocusList.traitOrigLocus = traitLocus;
    traitLocusList.pLocusIndex =
      (int *) calloc (traitLocusList.maxNumLocus, sizeof (int));
    traitLocusList.pLocusIndex[0] = 0;
    for (k = 0; k < 3; k++) {
      traitLocusList.pPrevLocusDistance[k] =
	(double *) calloc (savedLocusList.maxNumLocus, sizeof (double));
      traitLocusList.pNextLocusDistance[k] =
	(double *) calloc (savedLocusList.maxNumLocus, sizeof (double));

      traitLocusList.pPrevLocusDistance[k][0] = -1;
      traitLocusList.pNextLocusDistance[k][0] = -1;
    }
    /* populate the trait xmission matrix */
    locusList = &traitLocusList;
    xmissionMatrix = traitMatrix;
    status = populate_xmission_matrix (traitMatrix, 1, initialProbAddr,	/* probability */
				       initialProbAddr2,	/* probability */
				       initialHetProbAddr, 0,	/* cell index */
				       -1,	/* last he locus */
				       -1,	/* last het pattern (P-1 or M-2) */
				       0);	/* current locus - start with 0 */

    if (modelOptions.polynomial == TRUE) {
      holdAllPolys ();
    }

    /* for trait likelihood */
    locusList = &traitLocusList;
    xmissionMatrix = traitMatrix;
    if (pTrait->type == DICHOTOMOUS) {

      /*call compute_likelihood with dummy numbers to build polynomials */
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	pTrait->penetrance[2][liabIdx][0][0] = 0.7;
	pTrait->penetrance[2][liabIdx][0][1] = 0.5;
	pTrait->penetrance[2][liabIdx][1][0] = 0.5;
	pTrait->penetrance[2][liabIdx][1][1] = 0.3;
	pTrait->penetrance[1][liabIdx][0][0] = 1 - 0.7;
	pTrait->penetrance[1][liabIdx][0][1] = 1 - 0.5;
	pTrait->penetrance[1][liabIdx][1][0] = 1 - 0.5;
	pTrait->penetrance[1][liabIdx][1][1] = 1 - 0.3;
      }

      if (modelOptions.polynomial == TRUE);
      else
	/* only need to update trait locus */
	update_penetrance (&pedigreeSet, traitLocus);

      pLocus->pAlleleFrequency[0] = 0.5;
      pLocus->pAlleleFrequency[1] = 1 - 0.5;

      if (modelOptions.polynomial == TRUE);
      else
	update_locus (&pedigreeSet, traitLocus);
      /* get the likelihood for the trait */
      KLOG (LOGLIKELIHOOD, LOGDEBUG, "Trait Likelihood\n");
      compute_likelihood (&pedigreeSet);	/* This builds polynomials with dummy numbers */

    } else {			// QT
      pLocus->pAlleleFrequency[0] = 0.5;
      pLocus->pAlleleFrequency[1] = 1 - 0.5;
      update_locus (&pedigreeSet, traitLocus);

      /*call compute_likelihood with dummy numbers to build polynomials */
      for (liabIdx = 0; liabIdx < modelRange.nlclass; liabIdx++) {
	pTrait->means[liabIdx][0][0] = 2.0;
	pTrait->means[liabIdx][0][1] = 1.0;
	pTrait->means[liabIdx][1][0] = 1.0;
	pTrait->means[liabIdx][1][1] = 0.3;
	pTrait->stddev[liabIdx][0][0] = 1.0;
	pTrait->stddev[liabIdx][0][1] = 1.0;
	pTrait->stddev[liabIdx][1][0] = 1.0;
	pTrait->stddev[liabIdx][1][1] = 1.0;

	/* threshold for QT */
	pTrait->cutoffValue[liabIdx] = 0.5;

      }				/* liability class Index */
      if (modelOptions.polynomial == TRUE);
      else
	update_penetrance (&pedigreeSet, traitLocus);
 
      compute_likelihood (&pedigreeSet);

    }

    /* coply the polynomials built from above to traitLikelihoodPolynomials */
    if (modelOptions.polynomial == TRUE) {
      for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
	/* save the likelihood at trait */
	pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
	pPedigree->traitLikelihoodPolynomial =
	  pPedigree->likelihoodPolynomial;
	pPedigree->traitLikelihoodPolyList = pPedigree->likelihoodPolyList;
	pPedigree->likelihoodPolyList = NULL;
	pPedigree->likelihoodPolynomial = NULL;

	//fprintf(stderr,"Building traitPoly pedIdx =%d Null likelihood = %20.15f\n",pedIdx, pPedigree->likelihood);
	//fprintf(stderr,"pedIdx %d eType= %d\n", pedIdx, ((pPedigree->traitLikelihoodPolyList)->pList[0])->eType);
      }
    }


    /* get the trait locations we need to evaluate at */
    numPositions = modelRange.ntloc;
    mp_result = (SUMMARY_STAT *) calloc (numPositions, sizeof (SUMMARY_STAT));

    /* Need to output the results */
    dk_writeMPBRHeader ();
    dk_writeMPMODHeader ();
    
    prevFirstMarker = -1;
    prevLastMarker = -1;
    prevTraitInd = -1;
    leftMarker = -1;

    /* Iterate over all positions in the analysis. */
    for (posIdx = 0; posIdx < numPositions; posIdx++) {
      /* positions listed are sex average positions */
      traitPos = modelRange.tloc[posIdx];
      /* set the sex average position first 
       * the sex specific positions will be updated once markers are selected
       * as interpolation might be needed
       */
      pTraitLocus->mapPosition[0] = traitPos;
      pTraitLocus->mapPosition[1] = traitPos;
      pTraitLocus->mapPosition[2] = traitPos;
      /* initialize the locusList */
      locusList = &savedLocusList;
      memset (locusList->pLocusIndex, 0,
	      sizeof (int) * locusList->maxNumLocus);
      for (k = 0; k < 3; k++) {
	memset (&locusList->pPrevLocusDistance[k][0], 0,
		sizeof (double) * locusList->maxNumLocus);
	memset (&locusList->pNextLocusDistance[k][0], 0,
		sizeof (double) * locusList->maxNumLocus);
      }
      locusList->numLocus = 1;
      locusList->pLocusIndex[0] = traitLocus;
      for (k = 0; k < 3; k++) {
	locusList->pPrevLocusDistance[k][0] = -1;
	locusList->pNextLocusDistance[k][0] = -1;
      }
      /* select markers to be used for the multipoint analysis */
      add_markers_to_locuslist (locusList, modelType.numMarkers,
				&leftMarker, 0,
				originalLocusList.numLocus - 1, traitPos, 0);

      /* store the markers used */
      mp_result[posIdx].pMarkers =
	(int *) calloc (modelType.numMarkers, sizeof (int));
      k = 0;			/* marker index */
      for (i = 0; i < locusList->numLocus; i++) {
	j = locusList->pLocusIndex[i];
	if (originalLocusList.ppLocusList[j]->locusType == LOCUS_TYPE_MARKER) {
	  mp_result[posIdx].pMarkers[k] = j;
	  k++;
	} else {
	  mp_result[posIdx].trait = i;
	  traitIndex = i;
	}
      }
      locusList->traitLocusIndex = traitIndex;
      locusList->traitOrigLocus = traitLocus;
      markerSetChanged = FALSE;
      if (prevFirstMarker != mp_result[posIdx].pMarkers[0]
	  || prevLastMarker !=
	  mp_result[posIdx].pMarkers[modelType.numMarkers - 1]) {
	/* marker set has changed */
	markerSetChanged = TRUE;
	markerLocusList.pLocusIndex[0] = mp_result[posIdx].pMarkers[0];
	prevPos = get_map_position (markerLocusList.pLocusIndex[0]);
	/* set up marker set locus list */
	for (k = 1; k < modelType.numMarkers; k++) {
	  markerLocusList.pLocusIndex[k] = mp_result[posIdx].pMarkers[k];
	  currPos = get_map_position (markerLocusList.pLocusIndex[k]);
	  for (j = 0; j < 3; j++) {
	    markerLocusList.pPrevLocusDistance[j][k] =
	      markerLocusList.pNextLocusDistance[j][k - 1] =
	      cm_to_recombination_fraction (currPos[j] - prevPos[j],
					    map.mapFunction);
	  }
	  if (modelOptions.mapFlag == SA) {
	    for (j = 1; j <= 2; j++) {
	      markerLocusList.pPrevLocusDistance[j][k] =
		markerLocusList.pNextLocusDistance[j][k - 1] =
		markerLocusList.pPrevLocusDistance[0][k];
	    }
	  }
	  prevPos = currPos;
	}			/* end of loop over the markers to set up locus list */

	/* calculate likelihood for the marker set */
	locusList = &markerLocusList;
	xmissionMatrix = markerMatrix;
	if (modelOptions.polynomial == TRUE) {
	  pedigreeSetPolynomialClearance (&pedigreeSet);
	}

	/* save the polynomial flag */
	polynomialFlag = modelOptions.polynomial;

	/* populate the matrix */
	status = populate_xmission_matrix (markerMatrix, markerLocusList.numLocus, initialProbAddr,	/* probability */
					   initialProbAddr2,	/* probability */
					   initialHetProbAddr, 0,	/* cell index */
					   -1,	/* last he locus */
					   -1,	/* last het pattern (P-1 or M-2) */
					   0);	/* current locus - start with 0 */

	if (modelOptions.polynomial == TRUE)
	  freePolys ();

	print_xmission_matrix (markerMatrix, markerLocusList.numLocus, 0, 0,
			       tmpID);

	compute_likelihood (&pedigreeSet);
	modelOptions.polynomial = polynomialFlag;

	/* save the results for marker likelihood */
	for (pedIdx = 0; pedIdx < pedigreeSet.numPedigree; pedIdx++) {
	  /* save the likelihood at null */
	  pPedigree = pedigreeSet.ppPedigreeSet[pedIdx];
	  pPedigree->markerLikelihood = pPedigree->likelihood;
	}
	pedigreeSet.markerLikelihood = pedigreeSet.likelihood;
	pedigreeSet.log10MarkerLikelihood = pedigreeSet.log10Likelihood;
      }
      /* end of marker set change */

      prevFirstMarker = mp_result[posIdx].pMarkers[0];
      prevLastMarker = mp_result[posIdx].pMarkers[modelType.numMarkers - 1];
      if (markerSetChanged || prevTraitInd != mp_result[posIdx].trait)
	locusListChanged = TRUE;
      else
	locusListChanged = FALSE;
      prevTraitInd = mp_result[posIdx].trait;

      locusList = &savedLocusList;
      xmissionMatrix = altMatrix;
      /* interpolate trait postion for sex specific analysis */
      if (modelOptions.mapFlag == SS) {
	if (traitIndex == 0) {
	  marker1Pos = get_map_position (locusList->pLocusIndex[1]);
	  /* trait is the first one in the list */
	  if (traitPos < ERROR_MARGIN && traitPos > -ERROR_MARGIN) {
	    /* trait is at 0 */
	    pTraitLocus->mapPosition[MAP_MALE] =
	      pTraitLocus->mapPosition[MAP_FEMALE] = 0;
	  } else {
	    /* get the relative position on the sex average map */
	    relativePos = traitPos / marker1Pos[0];
	    pTraitLocus->mapPosition[MAP_MALE] =
	      relativePos * marker1Pos[MAP_MALE];
	    pTraitLocus->mapPosition[MAP_FEMALE] =
	      relativePos * marker1Pos[MAP_FEMALE];
	  }
	  /* update the inter locus distance - sex averaged already done before */
	  for (k = 1; k < 3; k++) {
	    locusList->pNextLocusDistance[k][0] =
	      locusList->pPrevLocusDistance[k][1] =
	      cm_to_recombination_fraction (marker1Pos[k] -
					    pTraitLocus->
					    mapPosition[k], map.mapFunction);
	  }
	} else if (traitIndex == modelType.numMarkers) {
	  /* trait is the last one in the list */
	  marker1Pos =
	    get_map_position (locusList->
			      pLocusIndex[modelType.numMarkers - 2]);
	  marker2Pos =
	    get_map_position (locusList->
			      pLocusIndex[modelType.numMarkers - 1]);
	  /* get the relative position on the sex average map */
	  dist = marker2Pos[0] - marker1Pos[0];
	  if (dist > ERROR_MARGIN) {
	    relativePos = (traitPos - marker2Pos[0]) / dist;
	    pTraitLocus->mapPosition[MAP_MALE] =
	      relativePos * (marker2Pos[MAP_MALE] -
			     marker1Pos[MAP_MALE]) + marker2Pos[MAP_MALE];
	    pTraitLocus->mapPosition[MAP_FEMALE] =
	      relativePos * (marker2Pos[MAP_FEMALE] -
			     marker1Pos[MAP_FEMALE]) + marker2Pos[MAP_FEMALE];
	  } else {
	    pTraitLocus->mapPosition[MAP_MALE] =
	      traitPos - marker2Pos[0] + marker2Pos[MAP_MALE];
	    pTraitLocus->mapPosition[MAP_FEMALE] =
	      traitPos - marker2Pos[0] + marker2Pos[MAP_FEMALE];
	  }

	  /* update the inter locus distance - sex averaged already done before */
	  for (k = 1; k <= 2; k++) {
	    locusList->pNextLocusDistance[k][traitIndex - 1] =
	      locusList->pPrevLocusDistance[k][traitIndex] =
	      cm_to_recombination_fraction (pTraitLocus->
					    mapPosition[k] -
					    marker2Pos[k], map.mapFunction);
	  }
	} else {
	  /* trait is in between two markers */
	  marker1Pos =
	    get_map_position (locusList->pLocusIndex[traitIndex - 1]);
	  marker2Pos =
	    get_map_position (locusList->pLocusIndex[traitIndex + 1]);
	  /* get the relative position on the sex average map */
	  dist = marker2Pos[0] - marker1Pos[0];
	  if (dist > ERROR_MARGIN) {
	    relativePos = (traitPos - marker1Pos[0]) / dist;
	    pTraitLocus->mapPosition[MAP_MALE] =
	      relativePos * (marker2Pos[MAP_MALE] -
			     marker1Pos[MAP_MALE]) + marker1Pos[MAP_MALE];
	    pTraitLocus->mapPosition[MAP_FEMALE] =
	      relativePos * (marker2Pos[MAP_FEMALE] -
			     marker1Pos[MAP_FEMALE]) + marker1Pos[MAP_FEMALE];
	  } else {
	    pTraitLocus->mapPosition[MAP_MALE] = marker1Pos[MAP_MALE];
	    pTraitLocus->mapPosition[MAP_FEMALE] = marker1Pos[MAP_FEMALE];
	  }
	  /* update the inter locus distance - sex averaged already done before */
	  for (k = 1; k < 3; k++) {
	    locusList->pNextLocusDistance[k][traitIndex - 1] =
	      locusList->pPrevLocusDistance[k][traitIndex] =
	      cm_to_recombination_fraction (pTraitLocus->
					    mapPosition[k] -
					    marker1Pos[k], map.mapFunction);
	    locusList->pNextLocusDistance[k][traitIndex] =
	      locusList->pPrevLocusDistance[k][traitIndex + 1] =
	      cm_to_recombination_fraction (marker2Pos[k] -
					    pTraitLocus->
					    mapPosition[k], map.mapFunction);

	  }
	}

      }
      /*end of SS model */

      /* the locus list has been built, go on to the analysis 
       * multipoint DT */
      if (markerSetChanged || locusListChanged) {
	if (modelOptions.polynomial == TRUE) {
	  pedigreeSetPolynomialClearance (&pedigreeSet);
	  /* populate the matrix */
	  status = populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr,	/* probability */
					     initialProbAddr2,	/* probability */
					     initialHetProbAddr, 0,	/* cell index */
					     -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
					     0);	/* current locus - start with 0 */
	  print_xmission_matrix (altMatrix, savedLocusList.numLocus, 0, 0,
				 tmpID);
	  if (modelOptions.polynomial == TRUE)
	    freePolys ();
	}
      }

      if (modelOptions.polynomial != TRUE);
	/* populate the matrix */
	status = populate_xmission_matrix (altMatrix, totalLoci, initialProbAddr,	/* probability */
					   initialProbAddr2,	/* probability */
					   initialHetProbAddr, 0,	/* cell index */
					   -1, -1,	/* last het locus & last het pattern (P-1 or M-2) */
					   0);	/* current locus - start with 0 */

      /* multipoint DT */
      integral = 0.0;
      abserr = 0.0;
      num_out_constraint = 0;

      num_eval = kelvin_dcuhre_integrate (&integral, &abserr, volume_region);

      /* calculate imputed PPL and print the results */
      if (integral > 0.214)
	ppl =(integral * integral) / (-5.77 + 54 * integral + integral * integral);
      else
	ppl = 0;

      dk_writeMPBRData (ppl, integral);
      dk_copyMaxModel (localmax_x, &dk_globalmax);
      dk_writeMPMODData (localmax_value, &dk_globalmax);

      if (fpDK != NULL) {
	fprintf (fpDK, "%f  %6.4f %12.8f %12.8f %d  %f\n", traitPos, ppl,
		 integral, abserr, num_eval,log10 (localmax_value));
	fflush(fpDK);
      }

    }				/* end of walking down the chromosome */
  }				/* end of multipoint */

  if (modelOptions.equilibrium != LINKAGE_EQUILIBRIUM) {
    free (dk_globalmax.dprime);
    free (dk_dprime0max.dprime);
    free (dk_theta0max.dprime);
  }
  free (dk_globalmax.pen);
  free (dk_dprime0max.pen);
  free (dk_theta0max.pen);
