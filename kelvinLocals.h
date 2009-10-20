unsigned long cL[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 },    ///< Actual # calls to each instance of compute_likelihood
    eCL[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};     ///< Est. final # calls to each instance of compute_likelihood

FILE *fpCond = NULL;    ///< Conditional LR for genetic counseling, global due to likelihood.c write!
FILE *fpHet = NULL;     ///< Average HET LR file (Bayes Ratio file) pointer
FILE *fpPPL = NULL;     ///< PPL output file pointer
FILE *fpMOD = NULL;     // MOD and maximizing model information
FILE *fpIR = NULL;      ///< Intermediate results, used to go to stderr, normally dkelvin-only
FILE *fpDK = NULL;      // DCHURE detail file

LambdaCell *pLambdaCell;
