
int liabIdx;
int mkrFreqIdx;
double mkrFreq;
int status;
int R_square_flag = FALSE;
double R_square = 0;
double ppl, ldppl, ppld, ppldGl;
int prevFirstMarker;		/* first marker in the set for multipoint analysis */
int prevLastMarker;		/* last marker in the set for multipoint analysis */
int prevTraitInd;
int leftMarker = -1;
int posIdx;
double traitPos;      /* trait position for multipoint analysis */
int j; 
int markerSetChanged; /* Flag for multipoint analysis, did set of markers change? */
int locusListChanged; /* flag for multipoint analysis, did relative trait position or marker set change? */
double *prevPos, *currPos;    /* for MP */
double *marker1Pos, *marker2Pos;
double dist;
int traitIndex = 0;
double relativePos;

int *BRscale;    /*scale per BR    added 6/17/2009*/
int max_scale; /*overall max scale*/
double newLog10BR;

  /* Variables for DCUHRE   added 1/2008 */
  double integral = 0.0, abserr = 0;
  int num_eval = 0;

  //int method_option=2;   * used to choose the integratio method
  // 1: kelvin_dcuhre with grid from the user
  // 2: kelvin_dcuhre with rule 13
  // */
  double le_small_theta = 0.0, le_big_theta = 0.0;
  double ld_small_theta = 0.0, ld_big_theta = 0.0;
  double ld_unlinked = 0.0;
  double le_unlinked =0.0;
  double volume_region = 1.0;

  double thetaSMSF=0.0,thetaBMSF=0.0, thetaSMBF=0.0, thetaBMBF=0.0;

  int num_BR;
  int num_sample_Dp_theta =141;
  int num_sample_SS_theta =260;

  /*Dcuhre rule points and weights for D' and theta or only theta */
  double dcuhre2[141][4] = { {0.0, 0.00234550385, 0.118463442, 0.0},
  {0.0, 0.0115382672, 0.239314336, 0.0},
  {0.0, 0.025, 0.284444444, 0.0},
  {0.0, 0.0384617328, 0.239314336, 0.0},
  {0.0, 0.0476544961, 0.118463442, 0.0},
  {0.0, 0.0711095347, 0.118463442, 0.0},
  {0.0, 0.153844405, 0.239314336, 0.0},
  {0.0, 0.275, 0.284444444, 0.0},
  {0.0, 0.396155595, 0.239314336, 0.0},
  {0.0, 0.478890465, 0.118463442, 0.0},
  {0.0000000000000, 0.0250000000000, 0.0084492309003, 0.0},
  {-0.2517129343453, 0.0250000000000, 0.0237714740190, 0.0},
  {0.2517129343453, 0.0250000000000, 0.0237714740190, 0.0},
  {-0.7013933644534, 0.0250000000000, 0.0294001617014, 0.0},
  {0.7013933644534, 0.0250000000000, 0.0294001617014, 0.0},
  {0.0000000000000, 0.0187071766414, 0.0237714740190, 0.0},
  {0.0000000000000, 0.0312928233586, 0.0237714740190, 0.0},
  {0.0000000000000, 0.0074651658887, 0.0294001617014, 0.0},
  {0.0000000000000, 0.0425348341113, 0.0294001617014, 0.0},
  {0.9590960631620, 0.0250000000000, 0.0066444364658, 0.0},
  {-0.9590960631620, 0.0250000000000, 0.0066444364658, 0.0},
  {0.0000000000000, 0.0489774015790, 0.0066444364658, 0.0},
  {0.0000000000000, 0.0010225984210, 0.0066444364658, 0.0},
  {0.9956010478552, 0.0250000000000, 0.0042536044255, 0.0},
  {-0.9956010478552, 0.0250000000000, 0.0042536044255, 0.0},
  {0.0000000000000, 0.0498900261964, 0.0042536044255, 0.0},
  {0.0000000000000, 0.0001099738036, 0.0042536044255, 0.0},
  {0.5000000000000, 0.0250000000000, 0.0000000000000, 0.0},
  {-0.5000000000000, 0.0250000000000, 0.0000000000000, 0.0},
  {0.0000000000000, 0.0375000000000, 0.0000000000000, 0.0},
  {0.0000000000000, 0.0125000000000, 0.0000000000000, 0.0},
  {0.1594544658298, 0.0289863616457, 0.0040664827466, 0.0},
  {-0.1594544658298, 0.0289863616457, 0.0040664827466, 0.0},
  {0.1594544658298, 0.0210136383543, 0.0040664827466, 0.0},
  {-0.1594544658298, 0.0210136383543, 0.0040664827466, 0.0},
  {0.3808991135940, 0.0345224778399, 0.0336223164632, 0.0},
  {-0.3808991135940, 0.0345224778399, 0.0336223164632, 0.0},
  {0.3808991135940, 0.0154775221601, 0.0336223164632, 0.0},
  {-0.3808991135940, 0.0154775221601, 0.0336223164632, 0.0},
  {0.6582769255267, 0.0414569231382, 0.0332008041365, 0.0},
  {-0.6582769255267, 0.0414569231382, 0.0332008041365, 0.0},
  {0.6582769255267, 0.0085430768618, 0.0332008041365, 0.0},
  {-0.6582769255267, 0.0085430768618, 0.0332008041365, 0.0},
  {0.8761473165029, 0.0469036829126, 0.0140936869250, 0.0},
  {-0.8761473165029, 0.0469036829126, 0.0140936869250, 0.0},
  {0.8761473165029, 0.0030963170874, 0.0140936869250, 0.0},
  {-0.8761473165029, 0.0030963170874, 0.0140936869250, 0.0},
  {0.9982431840532, 0.0499560796013, 0.0009770697703, 0.0},
  {-0.9982431840532, 0.0499560796013, 0.0009770697703, 0.0},
  {0.9982431840532, 0.0000439203987, 0.0009770697703, 0.0},
  {-0.9982431840532, 0.0000439203987, 0.0009770697703, 0.0},
  {0.9790222658168, 0.0412307108141, 0.0075319969436, 0.0},
  {-0.9790222658168, 0.0412307108141, 0.0075319969436, 0.0},
  {0.9790222658168, 0.0087692891859, 0.0075319969436, 0.0},
  {-0.9790222658168, 0.0087692891859, 0.0075319969436, 0.0},
  {0.6492284325645, 0.0494755566454, 0.0075319969436, 0.0},
  {-0.6492284325645, 0.0494755566454, 0.0075319969436, 0.0},
  {0.6492284325645, 0.0005244433546, 0.0075319969436, 0.0},
  {-0.6492284325645, 0.0005244433546, 0.0075319969436, 0.0},
  {0.8727421201131, 0.0339565366147, 0.0257718308672, 0.0},
  {-0.8727421201131, 0.0339565366147, 0.0257718308672, 0.0},
  {0.8727421201131, 0.0160434633853, 0.0257718308672, 0.0},
  {-0.8727421201131, 0.0160434633853, 0.0257718308672, 0.0},
  {0.3582614645881, 0.0468185530028, 0.0257718308672, 0.0},
  {-0.3582614645881, 0.0468185530028, 0.0257718308672, 0.0},
  {0.3582614645881, 0.0031814469972, 0.0257718308672, 0.0},
  {-0.3582614645881, 0.0031814469972, 0.0257718308672, 0.0},
  {0.5666666666667, 0.0301944444444, 0.0156250000000, 0.0},
  {-0.5666666666667, 0.0301944444444, 0.0156250000000, 0.0},
  {0.5666666666667, 0.0198055555556, 0.0156250000000, 0.0},
  {-0.5666666666667, 0.0198055555556, 0.0156250000000, 0.0},
  {0.2077777777778, 0.0391666666667, 0.0156250000000, 0.0},
  {-0.2077777777778, 0.0391666666667, 0.0156250000000, 0.0},
  {0.2077777777778, 0.0108333333333, 0.0156250000000, 0.0},
  {-0.2077777777778, 0.0108333333333, 0.0156250000000, 0.0},
  {0.0000000000000, 0.2750000000000, 0.0084492309003, 0.0},
  {-0.2517129343453, 0.2750000000000, 0.0237714740190, 0.0},
  {0.2517129343453, 0.2750000000000, 0.0237714740190, 0.0},
  {-0.7013933644534, 0.2750000000000, 0.0294001617014, 0.0},
  {0.7013933644534, 0.2750000000000, 0.0294001617014, 0.0},
  {0.0000000000000, 0.2183645897723, 0.0237714740190, 0.0},
  {0.0000000000000, 0.3316354102277, 0.0237714740190, 0.0},
  {0.0000000000000, 0.1171864929980, 0.0294001617014, 0.0},
  {0.0000000000000, 0.4328135070020, 0.0294001617014, 0.0},
  {0.9590960631620, 0.2750000000000, 0.0066444364658, 0.0},
  {-0.9590960631620, 0.2750000000000, 0.0066444364658, 0.0},
  {0.0000000000000, 0.4907966142114, 0.0066444364658, 0.0},
  {0.0000000000000, 0.0592033857886, 0.0066444364658, 0.0},
  {0.9956010478552, 0.2750000000000, 0.0042536044255, 0.0},
  {-0.9956010478552, 0.2750000000000, 0.0042536044255, 0.0},
  {0.0000000000000, 0.4990102357674, 0.0042536044255, 0.0},
  {0.0000000000000, 0.0509897642326, 0.0042536044255, 0.0},
  {0.5000000000000, 0.2750000000000, 0.0000000000000, 0.0},
  {-0.5000000000000, 0.2750000000000, 0.0000000000000, 0.0},
  {0.0000000000000, 0.3875000000000, 0.0000000000000, 0.0},
  {0.0000000000000, 0.1625000000000, 0.0000000000000, 0.0},
  {0.1594544658298, 0.3108772548117, 0.0040664827466, 0.0},
  {-0.1594544658298, 0.3108772548117, 0.0040664827466, 0.0},
  {0.1594544658298, 0.2391227451883, 0.0040664827466, 0.0},
  {-0.1594544658298, 0.2391227451883, 0.0040664827466, 0.0},
  {0.3808991135940, 0.3607023005587, 0.0336223164632, 0.0},
  {-0.3808991135940, 0.3607023005587, 0.0336223164632, 0.0},
  {0.3808991135940, 0.1892976994413, 0.0336223164632, 0.0},
  {-0.3808991135940, 0.1892976994413, 0.0336223164632, 0.0},
  {0.6582769255267, 0.4231123082435, 0.0332008041365, 0.0},
  {-0.6582769255267, 0.4231123082435, 0.0332008041365, 0.0},
  {0.6582769255267, 0.1268876917565, 0.0332008041365, 0.0},
  {-0.6582769255267, 0.1268876917565, 0.0332008041365, 0.0},
  {0.8761473165029, 0.4721331462132, 0.0140936869250, 0.0},
  {-0.8761473165029, 0.4721331462132, 0.0140936869250, 0.0},
  {0.8761473165029, 0.0778668537868, 0.0140936869250, 0.0},
  {-0.8761473165029, 0.0778668537868, 0.0140936869250, 0.0},
  {0.9982431840532, 0.4996047164120, 0.0009770697703, 0.0},
  {-0.9982431840532, 0.4996047164120, 0.0009770697703, 0.0},
  {0.9982431840532, 0.0503952835880, 0.0009770697703, 0.0},
  {-0.9982431840532, 0.0503952835880, 0.0009770697703, 0.0},
  {0.9790222658168, 0.4210763973270, 0.0075319969436, 0.0},
  {-0.9790222658168, 0.4210763973270, 0.0075319969436, 0.0},
  {0.9790222658168, 0.1289236026730, 0.0075319969436, 0.0},
  {-0.9790222658168, 0.1289236026730, 0.0075319969436, 0.0},
  {0.6492284325645, 0.4952800098088, 0.0075319969436, 0.0},
  {-0.6492284325645, 0.4952800098088, 0.0075319969436, 0.0},
  {0.6492284325645, 0.0547199901912, 0.0075319969436, 0.0},
  {-0.6492284325645, 0.0547199901912, 0.0075319969436, 0.0},
  {0.8727421201131, 0.3556088295323, 0.0257718308672, 0.0},
  {-0.8727421201131, 0.3556088295323, 0.0257718308672, 0.0},
  {0.8727421201131, 0.1943911704677, 0.0257718308672, 0.0},
  {-0.8727421201131, 0.1943911704677, 0.0257718308672, 0.0},
  {0.3582614645881, 0.4713669770255, 0.0257718308672, 0.0},
  {-0.3582614645881, 0.4713669770255, 0.0257718308672, 0.0},
  {0.3582614645881, 0.0786330229745, 0.0257718308672, 0.0},
  {-0.3582614645881, 0.0786330229745, 0.0257718308672, 0.0},
  {0.5666666666667, 0.3217500000000, 0.0156250000000, 0.0},
  {-0.5666666666667, 0.3217500000000, 0.0156250000000, 0.0},
  {0.5666666666667, 0.2282500000000, 0.0156250000000, 0.0},
  {-0.5666666666667, 0.2282500000000, 0.0156250000000, 0.0},
  {0.2077777777778, 0.4025000000000, 0.0156250000000, 0.0},
  {-0.2077777777778, 0.4025000000000, 0.0156250000000, 0.0},
  {0.2077777777778, 0.1475000000000, 0.0156250000000, 0.0},			     
  {-0.2077777777778, 0.1475000000000, 0.0156250000000, 0.0},
			     //			     {-0.949107912342759,0.5, 0.06474248308, 0.0},
			     //{-0.741531185599394,0.5, 0.13985269574, 0.0},
			     //{-0.405845151377397,0.5, 0.19090502525, 0.0},
			     {0.0,0.5,               1.0 , 0.0}//{0.0,0.5,               0.2089795184 , 0.0},
			     //{0.405845151377397,0.5, 0.19090502525, 0.0},
			     //{0.741531185599394,0.5, 0.13985269574, 0.0},
			     //{0.949107912342759,0.5, 0.06474248308, 0.0},
  };


  /*  DCUHRE Sex Specific theta points in 4 regions */
double thetaSS[260][4]= {{0.02500000000000, 0.02500000000000, 0.00844923090030, 0.0}, // 0< thetaM <0.05   0< thetaF <0.05  
{0.01870717664137, 0.02500000000000, 0.02377147401900, 0.0},
{0.03129282335863, 0.02500000000000, 0.02377147401900, 0.0},
{0.00746516588867, 0.02500000000000, 0.02940016170140, 0.0},
{0.04253483411133, 0.02500000000000, 0.02940016170140, 0.0},
{0.02500000000000, 0.01870717664140, 0.02377147401900, 0.0},
{0.02500000000000, 0.03129282335860, 0.02377147401900, 0.0},
{0.02500000000000, 0.00746516588870, 0.02940016170140, 0.0},
{0.02500000000000, 0.04253483411130, 0.02940016170140, 0.0},
{0.04897740157905, 0.02500000000000, 0.00664443646580, 0.0},
{0.00102259842095, 0.02500000000000, 0.00664443646580, 0.0},
{0.02500000000000, 0.04897740157900, 0.00664443646580, 0.0},
{0.02500000000000, 0.00102259842100, 0.00664443646580, 0.0},
{0.04989002619638, 0.02500000000000, 0.00425360442550, 0.0},
{0.00010997380362, 0.02500000000000, 0.00425360442550, 0.0},
{0.02500000000000, 0.04989002619640, 0.00425360442550, 0.0},
{0.02500000000000, 0.00010997380360, 0.00425360442550, 0.0},
{0.03750000000000, 0.02500000000000, 0.00000000000000, 0.0},
{0.01250000000000, 0.02500000000000, 0.00000000000000, 0.0},
{0.02500000000000, 0.03750000000000, 0.00000000000000, 0.0},
{0.02500000000000, 0.01250000000000, 0.00000000000000, 0.0},
{0.02898636164575, 0.02898636164570, 0.00406648274660, 0.0},
{0.02101363835425, 0.02898636164570, 0.00406648274660, 0.0},
{0.02898636164575, 0.02101363835430, 0.00406648274660, 0.0},
{0.02101363835425, 0.02101363835430, 0.00406648274660, 0.0},
{0.03452247783985, 0.03452247783990, 0.03362231646320, 0.0},
{0.01547752216015, 0.03452247783990, 0.03362231646320, 0.0},
{0.03452247783985, 0.01547752216010, 0.03362231646320, 0.0},
{0.01547752216015, 0.01547752216010, 0.03362231646320, 0.0},
{0.04145692313817, 0.04145692313820, 0.03320080413650, 0.0},
{0.00854307686183, 0.04145692313820, 0.03320080413650, 0.0},
{0.04145692313817, 0.00854307686180, 0.03320080413650, 0.0},
{0.00854307686183, 0.00854307686180, 0.03320080413650, 0.0},
{0.04690368291257, 0.04690368291260, 0.01409368692500, 0.0},
{0.00309631708743, 0.04690368291260, 0.01409368692500, 0.0},
{0.04690368291257, 0.00309631708740, 0.01409368692500, 0.0},
{0.00309631708743, 0.00309631708740, 0.01409368692500, 0.0},
{0.04995607960133, 0.04995607960130, 0.00097706977030, 0.0},
{0.00004392039867, 0.04995607960130, 0.00097706977030, 0.0},
{0.04995607960133, 0.00004392039870, 0.00097706977030, 0.0},
{0.00004392039867, 0.00004392039870, 0.00097706977030, 0.0},
{0.04947555664542, 0.04123071081410, 0.00753199694360, 0.0},
{0.00052444335458, 0.04123071081410, 0.00753199694360, 0.0},
{0.04947555664542, 0.00876928918590, 0.00753199694360, 0.0},
{0.00052444335458, 0.00876928918590, 0.00753199694360, 0.0},
{0.04123071081411, 0.04947555664540, 0.00753199694360, 0.0},
{0.00876928918589, 0.04947555664540, 0.00753199694360, 0.0},
{0.04123071081411, 0.00052444335460, 0.00753199694360, 0.0},
{0.00876928918589, 0.00052444335460, 0.00753199694360, 0.0},
{0.04681855300283, 0.03395653661470, 0.02577183086720, 0.0},
{0.00318144699717, 0.03395653661470, 0.02577183086720, 0.0},
{0.04681855300283, 0.01604346338530, 0.02577183086720, 0.0},
{0.00318144699717, 0.01604346338530, 0.02577183086720, 0.0},
{0.03395653661470, 0.04681855300280, 0.02577183086720, 0.0},
{0.01604346338530, 0.04681855300280, 0.02577183086720, 0.0},
{0.03395653661470, 0.00318144699720, 0.02577183086720, 0.0},
{0.01604346338530, 0.00318144699720, 0.02577183086720, 0.0},
{0.03916666666667, 0.03019444444440, 0.01562500000000, 0.0},
{0.01083333333333, 0.03019444444440, 0.01562500000000, 0.0},
{0.03916666666667, 0.01980555555560, 0.01562500000000, 0.0},
{0.01083333333333, 0.01980555555560, 0.01562500000000, 0.0},
{0.03019444444445, 0.03916666666670, 0.01562500000000, 0.0},
{0.01980555555556, 0.03916666666670, 0.01562500000000, 0.0},
{0.03019444444445, 0.01083333333330, 0.01562500000000, 0.0},
{0.01980555555556, 0.01083333333330, 0.01562500000000, 0.0},
{0.27500000000000, 0.02500000000000, 0.00844923090030, 0.0},// 0.05< thetaM <0.5   0< thetaF <0.05  
{0.21836458977231, 0.02500000000000, 0.02377147401900, 0.0},
{0.33163541022769, 0.02500000000000, 0.02377147401900, 0.0},
{0.11718649299799, 0.02500000000000, 0.02940016170140, 0.0},
{0.43281350700201, 0.02500000000000, 0.02940016170140, 0.0},
{0.27500000000000, 0.01870717664140, 0.02377147401900, 0.0},
{0.27500000000000, 0.03129282335860, 0.02377147401900, 0.0},
{0.27500000000000, 0.00746516588870, 0.02940016170140, 0.0},
{0.27500000000000, 0.04253483411130, 0.02940016170140, 0.0},
{0.49079661421145, 0.02500000000000, 0.00664443646580, 0.0},
{0.05920338578855, 0.02500000000000, 0.00664443646580, 0.0},
{0.27500000000000, 0.04897740157900, 0.00664443646580, 0.0},
{0.27500000000000, 0.00102259842100, 0.00664443646580, 0.0},
{0.49901023576742, 0.02500000000000, 0.00425360442550, 0.0},
{0.05098976423258, 0.02500000000000, 0.00425360442550, 0.0},
{0.27500000000000, 0.04989002619640, 0.00425360442550, 0.0},
{0.27500000000000, 0.00010997380360, 0.00425360442550, 0.0},
{0.38750000000000, 0.02500000000000, 0.00000000000000, 0.0},
{0.16250000000000, 0.02500000000000, 0.00000000000000, 0.0},
{0.27500000000000, 0.03750000000000, 0.00000000000000, 0.0},
{0.27500000000000, 0.01250000000000, 0.00000000000000, 0.0},
{0.31087725481170, 0.02898636164570, 0.00406648274660, 0.0},
{0.23912274518830, 0.02898636164570, 0.00406648274660, 0.0},
{0.31087725481170, 0.02101363835430, 0.00406648274660, 0.0},
{0.23912274518830, 0.02101363835430, 0.00406648274660, 0.0},
{0.36070230055865, 0.03452247783990, 0.03362231646320, 0.0},
{0.18929769944135, 0.03452247783990, 0.03362231646320, 0.0},
{0.36070230055865, 0.01547752216010, 0.03362231646320, 0.0},
{0.18929769944135, 0.01547752216010, 0.03362231646320, 0.0},
{0.42311230824351, 0.04145692313820, 0.03320080413650, 0.0},
{0.12688769175649, 0.04145692313820, 0.03320080413650, 0.0},
{0.42311230824351, 0.00854307686180, 0.03320080413650, 0.0},
{0.12688769175649, 0.00854307686180, 0.03320080413650, 0.0},
{0.47213314621315, 0.04690368291260, 0.01409368692500, 0.0},
{0.07786685378685, 0.04690368291260, 0.01409368692500, 0.0},
{0.47213314621315, 0.00309631708740, 0.01409368692500, 0.0},
{0.07786685378685, 0.00309631708740, 0.01409368692500, 0.0},
{0.49960471641197, 0.04995607960130, 0.00097706977030, 0.0},
{0.05039528358803, 0.04995607960130, 0.00097706977030, 0.0},
{0.49960471641197, 0.00004392039870, 0.00097706977030, 0.0},
{0.05039528358803, 0.00004392039870, 0.00097706977030, 0.0},
{0.49528000980878, 0.04123071081410, 0.00753199694360, 0.0},
{0.05471999019122, 0.04123071081410, 0.00753199694360, 0.0},
{0.49528000980878, 0.00876928918590, 0.00753199694360, 0.0},
{0.05471999019122, 0.00876928918590, 0.00753199694360, 0.0},
{0.42107639732701, 0.04947555664540, 0.00753199694360, 0.0},
{0.12892360267299, 0.04947555664540, 0.00753199694360, 0.0},
{0.42107639732701, 0.00052444335460, 0.00753199694360, 0.0},
{0.12892360267299, 0.00052444335460, 0.00753199694360, 0.0},
{0.47136697702545, 0.03395653661470, 0.02577183086720, 0.0},
{0.07863302297455, 0.03395653661470, 0.02577183086720, 0.0},
{0.47136697702545, 0.01604346338530, 0.02577183086720, 0.0},
{0.07863302297455, 0.01604346338530, 0.02577183086720, 0.0},
{0.35560882953232, 0.04681855300280, 0.02577183086720, 0.0},
{0.19439117046768, 0.04681855300280, 0.02577183086720, 0.0},
{0.35560882953232, 0.00318144699720, 0.02577183086720, 0.0},
{0.19439117046768, 0.00318144699720, 0.02577183086720, 0.0},
{0.40250000000001, 0.03019444444440, 0.01562500000000, 0.0},
{0.14749999999999, 0.03019444444440, 0.01562500000000, 0.0},
{0.40250000000001, 0.01980555555560, 0.01562500000000, 0.0},
{0.14749999999999, 0.01980555555560, 0.01562500000000, 0.0},
{0.32175000000000, 0.03916666666670, 0.01562500000000, 0.0},
{0.22825000000000, 0.03916666666670, 0.01562500000000, 0.0},
{0.32175000000000, 0.01083333333330, 0.01562500000000, 0.0},
{0.22825000000000, 0.01083333333330, 0.01562500000000, 0.0},
{0.02500000000000, 0.27500000000000, 0.00844923090030, 0.0}, // 0< thetaM <0.05   0.05< thetaF <0.5
{0.01870717664137, 0.27500000000000, 0.02377147401900, 0.0},
{0.03129282335863, 0.27500000000000, 0.02377147401900, 0.0},
{0.00746516588867, 0.27500000000000, 0.02940016170140, 0.0},
{0.04253483411133, 0.27500000000000, 0.02940016170140, 0.0},
{0.02500000000000, 0.21836458977260, 0.02377147401900, 0.0},
{0.02500000000000, 0.33163541022740, 0.02377147401900, 0.0},
{0.02500000000000, 0.11718649299830, 0.02940016170140, 0.0},
{0.02500000000000, 0.43281350700170, 0.02940016170140, 0.0},
{0.04897740157905, 0.27500000000000, 0.00664443646580, 0.0},
{0.00102259842095, 0.27500000000000, 0.00664443646580, 0.0},
{0.02500000000000, 0.49079661421100, 0.00664443646580, 0.0},
{0.02500000000000, 0.05920338578900, 0.00664443646580, 0.0},
{0.04989002619638, 0.27500000000000, 0.00425360442550, 0.0},
{0.00010997380362, 0.27500000000000, 0.00425360442550, 0.0},
{0.02500000000000, 0.49901023576760, 0.00425360442550, 0.0},
{0.02500000000000, 0.05098976423240, 0.00425360442550, 0.0},
{0.03750000000000, 0.27500000000000, 0.00000000000000, 0.0},
{0.01250000000000, 0.27500000000000, 0.00000000000000, 0.0},
{0.02500000000000, 0.38750000000000, 0.00000000000000, 0.0},
{0.02500000000000, 0.16250000000000, 0.00000000000000, 0.0},
{0.02898636164575, 0.31087725481130, 0.00406648274660, 0.0},
{0.02101363835425, 0.31087725481130, 0.00406648274660, 0.0},
{0.02898636164575, 0.23912274518870, 0.00406648274660, 0.0},
{0.02101363835425, 0.23912274518870, 0.00406648274660, 0.0},
{0.03452247783985, 0.36070230055910, 0.03362231646320, 0.0},
{0.01547752216015, 0.36070230055910, 0.03362231646320, 0.0},
{0.03452247783985, 0.18929769944090, 0.03362231646320, 0.0},
{0.01547752216015, 0.18929769944090, 0.03362231646320, 0.0},
{0.04145692313817, 0.42311230824380, 0.03320080413650, 0.0},
{0.00854307686183, 0.42311230824380, 0.03320080413650, 0.0},
{0.04145692313817, 0.12688769175620, 0.03320080413650, 0.0},
{0.00854307686183, 0.12688769175620, 0.03320080413650, 0.0},
{0.04690368291257, 0.47213314621340, 0.01409368692500, 0.0},
{0.00309631708743, 0.47213314621340, 0.01409368692500, 0.0},
{0.04690368291257, 0.07786685378660, 0.01409368692500, 0.0},
{0.00309631708743, 0.07786685378660, 0.01409368692500, 0.0},
{0.04995607960133, 0.49960471641170, 0.00097706977030, 0.0},
{0.00004392039867, 0.49960471641170, 0.00097706977030, 0.0},
{0.04995607960133, 0.05039528358830, 0.00097706977030, 0.0},
{0.00004392039867, 0.05039528358830, 0.00097706977030, 0.0},
{0.04947555664542, 0.42107639732690, 0.00753199694360, 0.0},
{0.00052444335458, 0.42107639732690, 0.00753199694360, 0.0},
{0.04947555664542, 0.12892360267310, 0.00753199694360, 0.0},
{0.00052444335458, 0.12892360267310, 0.00753199694360, 0.0},
{0.04123071081411, 0.49528000980860, 0.00753199694360, 0.0},
{0.00876928918589, 0.49528000980860, 0.00753199694360, 0.0},
{0.04123071081411, 0.05471999019140, 0.00753199694360, 0.0},
{0.00876928918589, 0.05471999019140, 0.00753199694360, 0.0},
{0.04681855300283, 0.35560882953230, 0.02577183086720, 0.0},
{0.00318144699717, 0.35560882953230, 0.02577183086720, 0.0},
{0.04681855300283, 0.19439117046770, 0.02577183086720, 0.0},
{0.00318144699717, 0.19439117046770, 0.02577183086720, 0.0},
{0.03395653661470, 0.47136697702520, 0.02577183086720, 0.0},
{0.01604346338530, 0.47136697702520, 0.02577183086720, 0.0},
{0.03395653661470, 0.07863302297480, 0.02577183086720, 0.0},
{0.01604346338530, 0.07863302297480, 0.02577183086720, 0.0},
{0.03916666666667, 0.32174999999960, 0.01562500000000, 0.0},
{0.01083333333333, 0.32174999999960, 0.01562500000000, 0.0},
{0.03916666666667, 0.22825000000040, 0.01562500000000, 0.0},
{0.01083333333333, 0.22825000000040, 0.01562500000000, 0.0},
{0.03019444444445, 0.40250000000030, 0.01562500000000, 0.0},
{0.01980555555556, 0.40250000000030, 0.01562500000000, 0.0},
{0.03019444444445, 0.14749999999970, 0.01562500000000, 0.0},
{0.01980555555556, 0.14749999999970, 0.01562500000000, 0.0},
{0.27500000000000, 0.27500000000000, 0.00844923090030, 0.0},// 0.05< thetaM <0.5   0.05< thetaF <0.5 
{0.21836458977231, 0.27500000000000, 0.02377147401900, 0.0},
{0.33163541022769, 0.27500000000000, 0.02377147401900, 0.0},
{0.11718649299799, 0.27500000000000, 0.02940016170140, 0.0},
{0.43281350700201, 0.27500000000000, 0.02940016170140, 0.0},
{0.27500000000000, 0.21836458977260, 0.02377147401900, 0.0},
{0.27500000000000, 0.33163541022740, 0.02377147401900, 0.0},
{0.27500000000000, 0.11718649299830, 0.02940016170140, 0.0},
{0.27500000000000, 0.43281350700170, 0.02940016170140, 0.0},
{0.49079661421145, 0.27500000000000, 0.00664443646580, 0.0},
{0.05920338578855, 0.27500000000000, 0.00664443646580, 0.0},
{0.27500000000000, 0.49079661421100, 0.00664443646580, 0.0},
{0.27500000000000, 0.05920338578900, 0.00664443646580, 0.0},
{0.49901023576742, 0.27500000000000, 0.00425360442550, 0.0},
{0.05098976423258, 0.27500000000000, 0.00425360442550, 0.0},
{0.27500000000000, 0.49901023576760, 0.00425360442550, 0.0},
{0.27500000000000, 0.05098976423240, 0.00425360442550, 0.0},
{0.38750000000000, 0.27500000000000, 0.00000000000000, 0.0},
{0.16250000000000, 0.27500000000000, 0.00000000000000, 0.0},
{0.27500000000000, 0.38750000000000, 0.00000000000000, 0.0},
{0.27500000000000, 0.16250000000000, 0.00000000000000, 0.0},
{0.31087725481170, 0.31087725481130, 0.00406648274660, 0.0},
{0.23912274518830, 0.31087725481130, 0.00406648274660, 0.0},
{0.31087725481170, 0.23912274518870, 0.00406648274660, 0.0},
{0.23912274518830, 0.23912274518870, 0.00406648274660, 0.0},
{0.36070230055865, 0.36070230055910, 0.03362231646320, 0.0},
{0.18929769944135, 0.36070230055910, 0.03362231646320, 0.0},
{0.36070230055865, 0.18929769944090, 0.03362231646320, 0.0},
{0.18929769944135, 0.18929769944090, 0.03362231646320, 0.0},
{0.42311230824351, 0.42311230824380, 0.03320080413650, 0.0},
{0.12688769175649, 0.42311230824380, 0.03320080413650, 0.0},
{0.42311230824351, 0.12688769175620, 0.03320080413650, 0.0},
{0.12688769175649, 0.12688769175620, 0.03320080413650, 0.0},
{0.47213314621315, 0.47213314621340, 0.01409368692500, 0.0},
{0.07786685378685, 0.47213314621340, 0.01409368692500, 0.0},
{0.47213314621315, 0.07786685378660, 0.01409368692500, 0.0},
{0.07786685378685, 0.07786685378660, 0.01409368692500, 0.0},
{0.49960471641197, 0.49960471641170, 0.00097706977030, 0.0},
{0.05039528358803, 0.49960471641170, 0.00097706977030, 0.0},
{0.49960471641197, 0.05039528358830, 0.00097706977030, 0.0},
{0.05039528358803, 0.05039528358830, 0.00097706977030, 0.0},
{0.49528000980878, 0.42107639732690, 0.00753199694360, 0.0},
{0.05471999019122, 0.42107639732690, 0.00753199694360, 0.0},
{0.49528000980878, 0.12892360267310, 0.00753199694360, 0.0},
{0.05471999019122, 0.12892360267310, 0.00753199694360, 0.0},
{0.42107639732701, 0.49528000980860, 0.00753199694360, 0.0},
{0.12892360267299, 0.49528000980860, 0.00753199694360, 0.0},
{0.42107639732701, 0.05471999019140, 0.00753199694360, 0.0},
{0.12892360267299, 0.05471999019140, 0.00753199694360, 0.0},
{0.47136697702545, 0.35560882953230, 0.02577183086720, 0.0},
{0.07863302297455, 0.35560882953230, 0.02577183086720, 0.0},
{0.47136697702545, 0.19439117046770, 0.02577183086720, 0.0},
{0.07863302297455, 0.19439117046770, 0.02577183086720, 0.0},
{0.35560882953232, 0.47136697702520, 0.02577183086720, 0.0},
{0.19439117046768, 0.47136697702520, 0.02577183086720, 0.0},
{0.35560882953232, 0.07863302297480, 0.02577183086720, 0.0},
{0.19439117046768, 0.07863302297480, 0.02577183086720, 0.0},
{0.40250000000001, 0.32174999999960, 0.01562500000000, 0.0},
{0.14749999999999, 0.32174999999960, 0.01562500000000, 0.0},
{0.40250000000001, 0.22825000000040, 0.01562500000000, 0.0},
{0.14749999999999, 0.22825000000040, 0.01562500000000, 0.0},
{0.32175000000000, 0.40250000000030, 0.01562500000000, 0.0},
{0.22825000000000, 0.40250000000030, 0.01562500000000, 0.0},
{0.32175000000000, 0.14749999999970, 0.01562500000000, 0.0},
{0.22825000000000, 0.14749999999970, 0.01562500000000, 0.0}};

/*********  end of local variable declaration   **********/
