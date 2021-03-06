#ifndef __kelvinGlobals_h__
#define __kelvinGlobals_h__

/**
@file kelvinGlobals.h

  Global variables and structures used by all dependent components of
  kelvin.
  
  Copyright (C) 2008, 2010, 2022 Mathematical Medicine LLC
  This program is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.
  You should have received a copy of the GNU General Public License along
  with this program. If not, see <https://www.gnu.org/licenses/>.

*/

#include "config/model.h"
#include "pedlib/pedlib.h"
#include "pedlib/likelihood.h" // For XMission matrix

extern char *programVersion;
extern char *svnVersion;
extern char *kelvinVersion;

extern struct swStopwatch *overallSW, *singleModelSW;

extern LambdaCell *pLambdaCell;

extern unsigned long grandTotalPairs;
extern unsigned long peakTotalPairs;

int loopMarkerFreqFlag;
int total_count;

PedigreeSet pedigreeSet;        /* Pedigrees. */
int locus; // Probably an evil global, i.e. recycled local
Polynomial *initialProbPoly2[3];
double initialProb[3];
Polynomial *initialProbPoly[3];

/** Transmission matrices provide the pre-computed probability of
 inheritance of a a given combination of marker and trait alleles. */
XMission *nullMatrix;
XMission *altMatrix;
XMission *traitMatrix;
XMission *markerMatrix;

/// Indicies into originalLocusList of loci in 2pt analysis, or first 2 in multipoint analysis
int loc1, loc2; 
/// Pointers to originalLocusList entries for loci in 2pth analysis, or first 2 in multipoint analysis
Locus *pLocus1, *pLocus2;

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
extern FILE *fpDry;     ///< Dry-run statistics output for sizing estimation   

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

st_DKMaxModel dk_globalmax, dk_dprime0max, dk_dprimeP1max, dk_dprimeN1max, dk_theta0max, dk_curModel;
st_DKMaxModel dk_localmax;   // 7/26/2018

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
