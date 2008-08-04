  int i, j, k;
  char configfile[KMAXFILENAMELEN] = "";
  int loc1, loc2;

  PedigreeSet pedigreeSet;      /* Pedigrees. */

  /* Start GAW */
  Pedigree *pPedigree;
  int liabIdx;
  int dprimeIdx;

  //double likelihood_null, likelihood_alternative;
  Locus *pLocus;
  Locus *pLocus1, *pLocus2;
  Trait *pTrait;
  int pedIdx;

  //  int pedID;
  LDLoci *pLDLoci = NULL;
  double traitPos;      /* trait position for multipoint analysis */
  TraitLocus *pTraitLocus;
  int traitLocus;
  int leftMarker = -1;
  int posIdx;
  double ppl;
  double ldppl, ppld;
  int markerSetChanged; /* flag for multipoint analysis */
  int locusListChanged; /* flag for multipoint analysis */
  int prevFirstMarker;  /* first marker in the set for multipoint analysis */
  int prevLastMarker;   /* last marker in the set for multipoint analysis */
  int prevTraitInd;
  double *prevPos, *currPos;    /* for MP */
  int locus;
  int R_square_flag = FALSE;
  int dprime0Idx = 0;
  int mkrFreqIdx;
  double mkrFreq;
  int totalLoci;
  int status;
  double *marker1Pos, *marker2Pos;
  double relativePos;
  int traitIndex = 0;
  double dist;
  int polynomialFlag;

  SubLocusList savedLocusList;
  SubLocusList traitLocusList;
  SubLocusList markerLocusList;

  Polynomial *initialProbPoly[3];
  Polynomial *initialProbPoly2[3];
  double initialProb[3];
  void *initialProbAddr[3];
  double initialProb2[3];
  void *initialProbAddr2[3];
  void *initialHetProbAddr[3];
  char *tmpID;
  int exitDueToLoop = FALSE;

#ifdef _OPENMP
  char *envVar;
  int threadCount = 0;
#endif

