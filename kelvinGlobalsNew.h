#ifndef __kelvinGlobalsNew_h__
#define __kelvinGlobalsNew_h__

#include "config/model_type.h"
#include "config/model_range.h"
#include "config/model_options.h"
#include "pedlib/pedlib.h"
#include "pedlib/likelihood.h"

extern char *programVersion;
extern char *kelvinVersion;

extern struct swStopwatch *overallSW;
//char *messageBuffer;

/* Some default global values. */


extern LambdaCell *pLambdaCell;

int loopMarkerFreqFlag ;
int total_count;

char *flexBuffer;
int flexBufferSize;

PedigreeSet pedigreeSet;	/* Pedigrees. */
Pedigree *pPedigree;
int pedIdx;
//int totalLoci;
int i;
TraitLocus *pTraitLocus;
int locus;
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
/* no longer needed, since MOD and MAX information goes to the same file now */
  //FILE *fpTP = NULL;      ///< Ancillary Two-point output, used to go to stderr
extern FILE *fpIR;      ///< Intermediate results, used to go to stderr, normally dkelvin-only
extern FILE *fpDK;      // DCHURE detail file


/**********************************************************************
 * Progress meters.
 **********************************************************************/
extern struct swStopwatch *overallSW;


#endif
