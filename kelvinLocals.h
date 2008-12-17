  struct swStopwatch *combinedComputeSW,        ///< Combined likelihood compute stopwatch
   *combinedBuildSW;    ///< Combined likelihood polynomial build stopwatch

  Pedigree *pPedigree;
  Polynomial *initialProbPoly2[3];
  Polynomial *initialProbPoly[3];
  TraitLocus *pTraitLocus;
  char configfile[KMAXFILENAMELEN] = "";
  double *marker1Pos, *marker2Pos;
  double *prevPos, *currPos;    /* for MP */
  double dist;
  double initialProb[3];
  double ldppl, ppld, ppldGl,ppldAl;
  double mkrFreq;
  double ppl;
  double relativePos;
  double traitPos;      /* trait position for multipoint analysis */
  int exitDueToLoop = FALSE;
  int i, j, k;
  int leftMarker = -1;
  int liabIdx;
  int locus;
  int mkrFreqIdx;
  int pedIdx;
  int posIdx;
  int prevTraitInd;
  int status;
  int traitIndex = 0;

#ifdef _OPENMP
  char *envVar;
  int threadCount = 0;
#endif

  unsigned long cL[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 },    ///< Actual # calls to each instance of compute_likelihood
    eCL[9] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0};   ///< Est. final # calls to each instance of compute_likelihood
