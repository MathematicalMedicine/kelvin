  struct swStopwatch *combinedComputeSW,        ///< Combined likelihood compute stopwatch
    *combinedBuildSW,    ///< Combined likelihood polynomial build stopwatch
*overallSW; 

  char configfile[KMAXFILENAMELEN] = "";

  unsigned long cL[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 },    ///< Actual # calls to each instance of compute_likelihood
    eCL[9] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0};   ///< Est. final # calls to each instance of compute_likelihood

FILE *fpCond = NULL;    ///< Conditional LR for genetic counseling, global due to likelihood.c write!
  FILE *fpHet = NULL;     ///< Average HET LR file (Bayes Ratio file) pointer
  FILE *fpPPL = NULL;     ///< PPL output file pointer
  FILE *fpMOD = NULL;     // MOD and maximizing model information
/* no longer needed, since MOD and MAX information goes to the same file now */
  //FILE *fpTP = NULL;      ///< Ancillary Two-point output, used to go to stderr
  FILE *fpIR = NULL;      ///< Intermediate results, used to go to stderr, normally dkelvin-only
  FILE *fpDK = NULL;      // DCHURE detail file

LambdaCell *pLambdaCell;
