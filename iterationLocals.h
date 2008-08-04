  struct swStopwatch *combinedComputeSW,        ///< Combined likelihood compute stopwatch
   *combinedBuildSW;    ///< Combined likelihood polynomial build stopwatch
  char **markerNameList = NULL;

  int initialFlag = 0;
  double max_at_dprime0;
  int maxTheta_at_dprime0 = -1;
  int theta0Idx = 0;
  double lr;
  double max_at_theta0;
  int maxDPrimeIdx = 0;
  int maxDPrimeIdx_at_theta0 = 0;
  double R_square = 0;
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
  int cL[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 },    ///< Actual # calls to each instance of compute_likelihood
    eCL[9] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0};   ///< Est. final # calls to each instance of compute_likelihood
  int breakFlag = FALSE;
  double alphaV, alphaV2;

