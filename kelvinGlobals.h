#ifndef __kelvinGlobals_h__
#define __kelvinGlobals_h__

/**
@file kelvinGlobals.h

  Global variables and structures used by all dependent components of
  kelvin.

  Copyright &copy; 2009, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

*/

#include "config/model_type.h"
#include "config/model_range.h"
#include "config/model_options.h"
#include "pedlib/pedlib.h"
#include "pedlib/likelihood.h" // For XMission matrix

extern char *programVersion;
extern char *kelvinVersion;

extern struct swStopwatch *overallSW;
//char *messageBuffer;

extern LambdaCell *pLambdaCell;

int loopMarkerFreqFlag;
int total_count;

char *flexBuffer;
int flexBufferSize;

PedigreeSet pedigreeSet;        /* Pedigrees. */
Pedigree *pPedigree; // Probably an evil global, i.e. recycled local
int pedIdx; // Probably an evil global, i.e. recycled local
TraitLocus *pTraitLocus; // Probably an evil global, i.e. recycled local
int locus; // Probably an evil global, i.e. recycled local
Polynomial *initialProbPoly2[3];
double initialProb[3];
Polynomial *initialProbPoly[3];
//void *initialProbAddr[3];

/** Transmission matrices provide the pre-computed probability of
 inheritance of a a given combination of marker and trait alleles. */
XMission *nullMatrix;
XMission *altMatrix;
XMission *traitMatrix;
XMission *markerMatrix;

int loc1, loc2;
Locus *pLocus1;
Locus *pLocus2;
Locus *pLocus;
Trait *pTrait;
int traitLocus;
int totalLoci;
void *initialProbAddr[3];
char *tmpID;

extern FILE *fpCond;    ///< Conditional LR for genetic counseling, global due to likelihood.c write!
extern FILE *fpHet;     ///< Average HET LR file (Bayes Ratio file) pointer
extern FILE *fpPPL;     ///< PPL output file pointer
extern FILE *fpMOD;     // MOD and maximizing model information
extern FILE *fpIR;      ///< Intermediate results, used to go to stderr, normally dkelvin-only
extern FILE *fpDK;      // DCHURE detail file

/* Moved here from integrationGlobals.h 6/18/2009 
   Use dk_curModel for printing surface points in fpIR*/
typedef struct
{
  double DD, Dd, dD, dd, DDSD, DdSD, dDSD, ddSD, threshold;
} st_DKMaxModelPenVector;

typedef struct
{
  int posIdx;   // which stores loc2 for 2pt and posIdx for mp
  double *dprime, theta[2], alpha, dgf, mf, r2;
  st_DKMaxModelPenVector *pen;
} st_DKMaxModel;

st_DKMaxModel dk_globalmax, dk_dprime0max, dk_theta0max, dk_curModel;

SubLocusList traitLocusList;
SubLocusList markerLocusList;
SubLocusList savedLocusList;

double initialProb2[3];
void *initialProbAddr2[3];
void *initialHetProbAddr[3];

extern int dprime0Idx;

int R_square_flag;
double R_square;
int leftMarker;
int traitIndex;

typedef struct ParamStruct
{
  int gfreqIdx;
  int mkrFreqIdx;
  int dprimeIdx;        /* could be multiple D' */
  int thetaIdx; /* could be multiple theta */
  int penIdx;   /* penetrance vector */
  int paramIdx; /* QT stdev index */
  int thresholdIdx;     /* QT threshold index */

  float gfreq;
  float R_square;
  float mkrFreq;

} ParamStruct;

#endif
