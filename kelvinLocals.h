  Pedigree *pPedigree;
  Polynomial *initialProbPoly2[3];
  Polynomial *initialProbPoly[3];
  SubLocusList markerLocusList;
  SubLocusList savedLocusList;
  SubLocusList traitLocusList;
  TraitLocus *pTraitLocus;
  char configfile[KMAXFILENAMELEN] = "";
  double *marker1Pos, *marker2Pos;
  double *prevPos, *currPos;    /* for MP */
  double dist;
  double initialProb2[3];
  double initialProb[3];
  double ldppl, ppld;
  double mkrFreq;
  double ppl;
  double relativePos;
  double traitPos;      /* trait position for multipoint analysis */
  int dprime0Idx = 0;
  int dprimeIdx;
  int exitDueToLoop = FALSE;
  int i, j, k;
  int leftMarker = -1;
  int liabIdx;
  int loc1, loc2;
  int locus;
  int locusListChanged; /* flag for multipoint analysis */
  int markerSetChanged; /* flag for multipoint analysis */
  int mkrFreqIdx;
  int pedIdx;
  int polynomialFlag;
  int posIdx;
  int prevFirstMarker;  /* first marker in the set for multipoint analysis */
  int prevLastMarker;   /* last marker in the set for multipoint analysis */
  int prevTraitInd;
  int status;
  int totalLoci;
  int traitIndex = 0;
  int traitLocus;
  void *initialProbAddr2[3];
  void *initialProbAddr[3];

#ifdef _OPENMP
  char *envVar;
  int threadCount = 0;
#endif

