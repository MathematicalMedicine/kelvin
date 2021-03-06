/* Copyright (C) 2008, 2010, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */
  char **markerNameList = NULL;

  int initialFlag = 0;
  double max_at_dprime0;
  int maxTheta_at_dprime0 = -1;
  int theta0Idx = 0;
  double lr;
  double max_at_theta0;
  int maxDPrimeIdx = 0;
  int maxDPrimeIdx_at_theta0 = 0;
  int thresholdIdx = -1;
  double threshold = 0;
  double avgLR;
  double constraint;
  double max;
  double log10HetLR;
  double log10_likelihood_null, log10_likelihood_alternative;
  double likelihood_ratio;
  double log10_likelihood_ratio;
  int paramIdx = -1;
  double tmp;
  int maxThetaIdx = 0;
  double adjustedHetLR = 0;
  double homoLR, hetLR;
  int penIdx, gfreqInd, thetaInd;
  double pen_DD, pen_Dd, pen_dD, pen_dd;
  double mean_DD, mean_Dd, mean_dD, mean_dd;
  double SD_DD, SD_Dd, SD_dD, SD_dd;
  double theta[2];      /* theta */
  int breakFlag = FALSE;
  double alphaV, alphaV2;
  double gfreq; /* disease gene frequency */
  int ret;

int markerSetChanged; /* Flag for multipoint analysis, did set of markers change? */

int prevFirstMarker;		/* first marker in the set for multipoint analysis */
int prevLastMarker;		/* last marker in the set for multipoint analysis */

/* Variables became global from local */
PedigreeSet pedigreeSet;	/* Pedigrees. */

/** Storage for the NULL likelihood for the multipoint calculation under polynomial. */
double markerSetLikelihood;

/** For multipoint, we use genetic map positions on a chromosome. */
double *map_position;
int num_of_map_position;

int status;

Pedigree *pPedigree;
TraitLocus *pTraitLocus;
double *marker1Pos, *marker2Pos;
double *prevPos, *currPos;    /* for MP */
double dist;
double ldppl, ppld, ppldGl;
double mkrFreq;
double ppl;
double relativePos;
double traitPos;      /* trait position for multipoint analysis */
int i, j, k;
int liabIdx;
int locus;
int mkrFreqIdx;
int pedIdx;
int posIdx;
int prevTraitInd;
int dprimeIdx;

