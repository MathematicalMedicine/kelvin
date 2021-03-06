/** @file polynomial.h

 Copyright (C) 2007, 2008, 2010, 2022 Mathematical Medicine LLC
 This program is free software: you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the Free
 Software Foundation, either version 3 of the License, or (at your option)
 any later version.
 You should have received a copy of the GNU General Public License along
 with this program. If not, see <https://www.gnu.org/licenses/>.

*/
#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

#include <stdio.h>

/* The following dynamically-maintained variables are for debugging. See
   their definition in polynomial.c for descriptions. */

extern int maxHashLength, constantHashHits, variableHashHits,
  functionHashHits;
extern int sumReleaseableCount, sumNotReleaseableCount,
  sumReturnConstantCount, sumReturnSameConstantCount, sumReturn1TermCount, 
  sumHashHits, sumNewCount, sumListNewCount, sumListReplacementCount, 
  sumFreedCount, sum1stTermsFreedCount;
extern int productReleaseableCount, productNotReleaseableCount,
  productReturn0Count, productReturnConstantCount, productReturnSameConstantCount,
  productReturn1stTermCount, productReturn1TermSumCount, productHashHits, 
  productHashHitIsSumCount, productReturnNormalCount, productNon1FactorIsSumCount, 
  productListNewCount, productListReplacementCount, productFreedCount, 
  product1stTermsFreedCount;
extern int constantPLExpansions, variablePLExpansions, sumPCollectExpansions,
  sumPTermMergeExpansions, sumPListExpansions, productPCollectExpansions,
  productPTermMergeExpansions, productPListExpansions;
extern long nodeId;
extern int constantCount, variableCount, sumCount, productCount,
  functionCallCount;
extern int polyListSortingCount, 
  keepPolyCount, freePolysCount, holdPolyCount, holdAllPolysCount, 
  unHoldPolyCount, freeKeptPolysCount, freePolysAttemptCount;
extern long evaluatePolyCount, evaluateValueCount;
extern int containerExpansions;
extern unsigned long totalSPLLengths, totalSPLCalls, lowSPLCount, highSPLCount;
extern unsigned long initialHashSize;
extern int pendingExplicitDiscards;


/* These are the types of polynomials */
enum expressionType
{
  T_CONSTANT = 0,		// a constant value for example, 0, 1, 1.5 ...
  T_VARIABLE = 1,		// a variable such as x, y, z, ...
  T_SUM = 2,			// a sum such as 2x+3y+5.6z+...
  T_PRODUCT = 3,		// a product such as x^2*y^10*z^100
  T_FUNCTIONCALL = 4,		// a function call such as log10(x)
  T_EXTERNAL = 5,               // an external reference, e.g. compiled DL
  T_FREED = 6,			// freed, but the structure kept for diagnostics
  T_OFFLINE = 7			// not in memory currently, must be restored
};

#ifndef __POLYNOMIAL_INTERNAL_H__

/* Maximum size of a polynomial function name */
#define MAX_PFN_LEN 128
/* Maximum number of DLs supporting a single polynomial (up to 32K modules!) */
#define MAX_POLY_DL 32

typedef struct polynomial
{
  unsigned char eType;		// polynomial type
  double value;			// the value of the polynomial
  unsigned int id;		// unique id - 4 bytes
} Polynomial;

/* The following are convenience constants. */
extern Polynomial *constant0Poly;
extern Polynomial *constant1Poly;

/* Optimized list for polynomial evaluation.  When we evaluate a polynomial multiple
   times, it is more efficient to do a single traversal of the the polynomial tree
   to build a list of unique terms in dependency order, and then drive evaluation
   from the list. */

typedef struct polyList
{
  int listSize;			// size of the preallocated pList
  int listNext;			// next free position
  struct polynomial **pList;	// list of polynomials for evaluation
} polynomialList;

#endif

/* Prototypes */

// Initialization before polynomials care created and evaluated.
void polynomialInitialization (int polynomialScale);

// Constructor for a constant polynomial
Polynomial *constantExp (char *fileName, int lineNo, double con);

#define constantExp(con) constantExp(__FILE__, __LINE__, con)

// Constructor for a variable polynomial
Polynomial *variableExp (double *vD, int *vI, char vType,
				char name[10]);

// Constructor of a sum polynomial
Polynomial *plusExp (char *fileName, int lineNo, int num, ...);

#define plusExp(num, ...) plusExp(__FILE__, __LINE__, num, __VA_ARGS__)

// Constructor of a product polynomial
Polynomial *timesExp (char *fileName, int lineNo, int num, ...);

#define timesExp(num, ...) timesExp(__FILE__, __LINE__, num, __VA_ARGS__)

// Constructor of a functionCall polynomial
Polynomial *functionCallExp (int num, ...);

// Evaluate a polynomial recursively by traversal
double evaluateValue (Polynomial *p);

// Create and initialize an evaluation list for a polynomial
struct polyList *buildPolyList ();

// Sort and trim the subpolynomials on the evaluation list
void polyListSorting (Polynomial *p, struct polyList *l);

//Evaluate a polynomial with the help of the evaluation list 
void evaluatePoly (Polynomial *pp, struct polyList *l, double *pd);

// Recursively print a polynomial with full expansion and no annotation.
void expPrinting (Polynomial *p);

// Recursively print a polynomial to specified depth and with annotation as needed.
void expTermPrinting (FILE *outputFile, Polynomial *p, int depth);

// Dump dynamically collected statistics (no performance impact)
void polyDynamicStatistics (char *);

// Dump dynamic statistics and traverse structures for more information.
void polyStatistics (char *);

// Clear all data structures related to the polynomials
void polynomialClearance ();

// Flag the specified polynomial so it survives freePolys() calls
void keepPoly (Polynomial *);

// Flag the specified polynomial so it survives freePolys() & freeKeptPolys() calls
void holdPoly (Polynomial *);

// Unflag the specified polynomial so it won't survive free calls
void unHoldPoly (Polynomial *);

// Deallocate all polynomials not kept or held
void freePolys ();

// Deallocate all polynomials not held
void freeKeptPolys ();

// Flag all polynomials so they survive freePolys() & freeKeptPolys() calls
void holdAllPolys ();

// Print out all the polynomials for diagnostic purposes
void printAllVariables ();
void printAllPolynomials ();
void releaseExternalPoly (Polynomial *);

// Print a summary of the polynomial tree components (traversal)
void printSummaryPoly (Polynomial *p);


Polynomial *restoreExternalPoly (char * name);
void codePoly (Polynomial * p, struct polyList * l, char * name);

#ifdef SOURCEDIGRAPH
// Dump data for a Dot digraph of polynomial construction-to-use source lines
void dumpSourceParenting ();
#endif
// Dump data for a Dot digraph of polynomial tree structure
void writePolyDigraph (Polynomial *);

#endif

// Explicitly discard a polynomial
void discardPoly (Polynomial * p);
