#ifndef __iterationMain_h__
#define __iterationMain_h__

typedef struct ParamStruct {
  int gfreqIdx;
  int mkrFreqIdx;
  int dprimeIdx; /* could be multiple D' */
  int thetaIdx; /* could be multiple theta */
  int penIdx;   /* penetrance vector */
  int paramIdx; /* QT stdev index */
  int thresholdIdx; /* QT threshold index */

  float gfreq;
  float R_square;
  float mkrFreq;
  
} ParamStruct;

void iterateMain();

#endif
