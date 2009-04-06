/** @file polynomial.h

 Copyright 2008, Nationwide Children's Research Institute.  
 All rights reserved.
 Permission is hereby given to use this software 
 for non-profit educational purposes only.

*/
#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

#include <stdio.h>
#include "utils.h"		/* Kelvin utilities. */

// This is primarily for Cygwin
#ifndef RTLD_LOCAL
#define RTLD_LOCAL 0
#endif
#ifndef RTLD_GLOBAL
#define RTLD_GLOBAL RTLD_LOCAL
#endif

#if defined (DMTRACK) || defined (DMUSE)
#warning "Dynamic memory usage dumping is turned on, so performance will be poor!"
#define malloc(X) swMalloc((X), __FILE__, __LINE__)
#define calloc(X,Y) swCalloc((X),(Y), __FILE__, __LINE__)
#define realloc(X,Y) swRealloc((X),(Y), __FILE__, __LINE__)
#define free(X) swFree((X), __FILE__, __LINE__)
#endif

/* The following dynamically-maintained variables are for debugging. See
   their definition in polynomial.c for descriptions. */

extern int maxHashLength, constantHashHits, variableHashHits,
  functionHashHits;
extern int sumReleaseableCount, sumNotReleaseableCount,
  sumReturnConstantCount, sumReturn1TermCount, sumHashHits, sumNewCount,
  sumListNewCount, sumListReplacementCount, sumFreedCount,
  sum1stTermsFreedCount;
extern int productReleaseableCount, productNotReleaseableCount,
  productReturn0Count, productReturnConstantCount, productReturn1stTermCount,
  productReturn1TermSumCount, productHashHits, productHashHitIsSumCount,
  productReturnNormalCount, productNon1FactorIsSumCount, productListNewCount,
  productListReplacementCount, productFreedCount, product1stTermsFreedCount;
extern int constantPLExpansions, variablePLExpansions, sumPCollectExpansions,
  sumPTermMergeExpansions, sumPListExpansions, productPCollectExpansions,
  productPTermMergeExpansions, productPListExpansions;
extern int nodeId, constantCount, variableCount, sumCount, productCount,
  functionCallCount;
extern int polyListSortingCount, evaluatePolyCount, evaluateValueCount, 
  keepPolyCount, freePolysCount, holdPolyCount, holdAllPolysCount, 
  unHoldPolyCount, freeKeptPolysCount, freePolysAttemptCount;
extern int containerExpansions;
extern unsigned long totalSPLLengths, totalSPLCalls, lowSPLCount, highSPLCount;
extern unsigned long initialHashSize;

/* If CPU utilization % is below this in a check interval, we're thrashing
   and should take some kind of action (currently exit) */
#define THRASH_CPU 5

/* Maximum size of a polynomial function name */
#define MAX_PFN_LEN 128
/* Maximum number of DLs supporting a single polynomial (up to 32K modules!) */
#define MAX_POLY_DL 32

/* These are the types of polynomials */
enum expressionType
{
  T_CONSTANT = 0,		// a constant value for example, 0, 1, 1.5 ...
  T_VARIABLE = 1,		// a variable such as x, y, z, ...
  T_SUM = 2,			// a sum such as 2x+3y+5.6z+...
  T_PRODUCT = 3,		// a product such as x^2y^10z^100
  T_FUNCTIONCALL = 4,		// a function call such as log10(x)
  T_EXTERNAL = 5,               // an external reference, e.g. compiled DL
  T_FREED = 6,			// freed, but the structure kept for diagnostics
  T_OFFLINE = 7			// not in memory currently, must be restored
};

/* Constants are represented as polynomials, but have no subcomponent - just a value.
   Why do we use constant polynomials instead of just storing the constants directly?

  When we try to understand why we need constant polynomials, we need to consider 
  that there are two worlds:  polynomial world and non-polynomial world.  In the
  non-polynomial world, the result of mathematical computations is a value.  In the
  polynomial world, the result of mathematical computations is another polynomial,
  whether the result is a constant, a sum of a group of terms, or a product of a group 
  of terms.  Therefore, constant polynomials are for purpose of consistency in polynomial 
  calculations and polynomial expressions. There also aren't enough of them to constitute
  a real waste of memory. */

/* This structure represents the elements of a variable. A variable has an address 
   in memory (integer or double) for the client to populate, and a name.
   In general, the number of variable polynomials is very small, say 20. */

struct variablePoly
{
  char vType;			// variable type
  union
  {
    double *vAddrD;		// address for a double variable
    int *vAddrI;		// address for an integer variable
  } vAddr;
  char vName[100];		// variable name
};

/* This structure represents the elements of a sum. A sum is composed of a list
   of subpolynomial terms and their factors. */

struct sumPoly
{
  int num;			// number of terms - 4 bytes
#ifdef MIN_USE_SSD_DPS
  int iMTLIndex;                // -1 if not in-memory, otherwise the index to the iMTL.
#endif
  struct polynomial **sum;	// polynomial terms - 8 bytes
  double *factor;		// factors for polynomial terms - 8 bytes
}; // 20(24) bytes

/* This structure represents the elements of a product. A product is composed of a 
   list of subpolynomial terms and their exponents. */

struct productPoly
{
  int num;			// number of terms - 4 bytes
  struct polynomial **product;	// polynomial terms - 8 bytes
  int *exponent;		// exponents for polynomial terms - 4 bytes
}; // 16 bytes

/* This structure represents the elements of a function call. Each function call is 
   composed of the function name, and a list of parameters. */

struct functionPoly
{
  int num;			// number of parameters
  struct polynomial **para;	// parameters
  char *name;			// function name
};

/* This structure represents a term that has been moved outside of the context of 
   current polynomial memory structures. It is not created by the caller, but internally
   when a referenced polynomial of some other type is flushed from memory. */

struct externalPoly
{
  int fileOK;
  int entryOK;
  char polynomialFunctionName[MAX_PFN_LEN+1];
  double (*polynomialFunctionRoutine)();
  void *polynomialFunctionHandle[MAX_POLY_DL];
};

/* This structure represents a general polynomial. It is composed of
   a unique id, an polynomial type (eType), a value, and, for types
   other than constant polynomials, a pointer to the type-specific 
   polynomial structure. An attempt to reduce memory consumption by
   eliminating the linking pointer between the base polynomial and
   type-specific storage only reduced storage consumption by ~5% and
   increased runtime significantly, and so has been abandoned. */

typedef struct polynomial
{
  unsigned int id;		// unique id - 4 bytes
  int index;			// index in a polynomial list - 4 bytes
  int key;			// key of the polynomial - 4 bytes
  unsigned short count;		// hold reference count - 2 bytes
  unsigned char valid;		// preservation flag(s) - 1 byte
  unsigned char eType;		// polynomial type:  - 1 byte
#ifdef SOURCEDIGRAPH
  unsigned char source;		// index of entry in polySources
#endif
  double value;			// the value of the polynomial - 8 bytes
  union
  {
    struct variablePoly *v;	// variable
    struct sumPoly *s;		// sum
    struct productPoly *p;	// product
    struct functionPoly *f;	// function
    struct externalPoly *e;     // external
  } e;				// unused by constants - 8 bytes
  unsigned char oldEType;
} Polynomial; // 4 + 4 + 4 + 2 + 1 + 1 (+ 1) + 8 + 8 = 32 (33->40) bytes

/* Bit masks for the polynomial valid flag. */
#define VALID_EVAL_FLAG 1	// Used by tree traversal routines to limit to unique terms
#define VALID_KEEP_FLAG 2	// Kept. Weaker than HOLD, only kept until a freeKeptPolys() call
#define VALID_REF_FLAG 4	// Multiply-referenced
#define VALID_TOP_FLAG 8        // Explicitly kept or held (top polynomial in call)

/* Track the full source code module name and line number of calls to create
   polynomials so that we only have to store an index (unsigned char) with 
   the polynomial itself to know it's origin. */
#ifdef SOURCEDIGRAPH
#define MAXPOLYSOURCES 256	// maximum number of calls to plusExp or timesExp in source
#endif

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

/* hashStruct is used for facilitating the identification of redundant polynomials,
   which is the key to efficient storage and evaluation. Identification of redundancy
   is a two-step process. First a type-specific hash is performed on the polynomial to
   find a candidate bucket, and then all current entries in that bucket are compared
   to the potential new entry to see if it is redundant. */

struct hashStruct
{
  unsigned short num;		// number of polynomials in this hash bucket
  unsigned short length;	// length of the preallocated index list
  int *key;			// key of polynomials in this bucket
  int *index;			// indexes of unique polynomials matching this key
};

/** @ingroup PolyCons
    @{
  Hash and other storage sizes for claustrophobic environments. Will be
  scaled-up as needed in polynomialInitialization. Original scale-up was 
  10x. Currently 1 is the default if no other value is provided. Use larger
  values to gain run-time efficiency in large analyses, but watch as you
  run out of memory -- that's always been the real problem.

  All the following constants are prime numbers.

*/
#define MAX_POLYNOMIAL_KEY   2147483629
#define MIN_CONSTANT_HASH_SIZE 9991
#define MIN_VARIABLE_HASH_SIZE 97
#define MIN_SUM_HASH_SIZE      199993
#define MIN_PRODUCT_HASH_SIZE  199993
#define MIN_FUNCTIONCALL_HASH_SIZE 9991
#define MIN_HASH_TABLE_INCREASE 2
#define MIN_CONSTANT_LIST_INITIAL 1000
#define MIN_CONSTANT_LIST_INCREASE 1000
#define MIN_VARIABLE_LIST_INITIAL 5
#define MIN_VARIABLE_LIST_INCREASE 5
#define MIN_EXTERNAL_LIST_INITIAL 5
#define MIN_EXTERNAL_LIST_INCREASE 5
#define MIN_SUM_LIST_INITIAL 100000
#define MIN_SUM_LIST_INCREASE 10000
#define MIN_PRODUCT_LIST_INITIAL 100000
#define MIN_PRODUCT_LIST_INCREASE 10000
#define MIN_FUNCTIONCALL_LIST_INITIAL 100
#define MIN_FUNCTIONCALL_LIST_INCREASE 10
/*@}*/

/* Prototypes */

// Initialization before polynomials care created and evaluated.
void polynomialInitialization (int polynomialScale);

// Constructor for a constant polynomial
struct polynomial *constantExp (double con);

// Constructor for a variable polynomial
struct polynomial *variableExp (double *vD, int *vI, char vType,
				char name[10]);

// Constructor of a sum polynomial
struct polynomial *plusExp (char *fileName, int lineNo, int num, ...);

#define plusExp(num, ...) plusExp(__FILE__, __LINE__, num, __VA_ARGS__)

// Constructor of a product polynomial
struct polynomial *timesExp (char *fileName, int lineNo, int num, ...);

#define timesExp(num, ...) timesExp(__FILE__, __LINE__, num, __VA_ARGS__)

// Constructor of a functionCall polynomial
struct polynomial *functionCallExp (int num, ...);

// Evaluate a polynomial recursively by traversal
double evaluateValue (struct polynomial *p);

// Create and initialize an evaluation list for a polynomial
struct polyList *buildPolyList ();

// Sort and trim the subpolynomials on the evaluation list
void polyListSorting (struct polynomial *p, struct polyList *l);

//Evaluate a polynomial with the help of the evaluation list 
void evaluatePoly (struct polynomial *pp, struct polyList *l, double *pd);

// Recursively print a polynomial with full expansion and no annotation.
void expPrinting (struct polynomial *p);

// Recursively print a polynomial to specified depth and with annotation as needed.
void expTermPrinting (FILE * outputFile, struct polynomial *p, int depth);

// Dump dynamically collected statistics (no performance impact)
void polyDynamicStatistics (char *);

// Dump dynamic statistics and traverse structures for more information.
void polyStatistics (char *);

// Clear all data structures related to the polynomials
void polynomialClearance ();

// Flag the specified polynomial so it survives freePolys() calls
void keepPoly (struct polynomial *);

// Flag the specified polynomial so it survives freePolys() & freeKeptPolys() calls
void holdPoly (struct polynomial *);

// Unflag the specified polynomial so it won't survive free calls
void unHoldPoly (struct polynomial *);

// Deallocate all polynomials not kept or held
void freePolys ();

// Deallocate all polynomials not held
void freeKeptPolys ();

// Flag all polynomials so they survive freePolys() & freeKeptPolys() calls
void holdAllPolys ();

// Print out all the polynomials for diagnostic purposes
void printAllVariables ();
void printAllPolynomials ();

// Print a summary of the polynomial tree components (traversal)
void printSummaryPoly (struct polynomial *p);

#ifdef SOURCEDIGRAPH
// Dump data for a Dot digraph of polynomial construction-to-use source lines
void dumpSourceParenting ();
#endif
// Dump data for a Dot digraph of polynomial tree structure
void writePolyDigraph (Polynomial *);

#endif

void codePoly (Polynomial * p, struct polyList * l, char * name);
Polynomial *restoreExternalPoly (char * name);
int loadPolyDL (Polynomial * p);
Polynomial * exportPoly (Polynomial * p);
Polynomial * importPoly (Polynomial * p);
void importTermList (Polynomial * p);
void exportTermList (Polynomial * p, int writeFlag);
void deportTermList (Polynomial * p);
void thrashingCheck ();
void releaseExternalPoly (Polynomial *);
