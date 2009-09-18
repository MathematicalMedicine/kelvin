/** @file polynomial.h

 Copyright 2008, Nationwide Children's Research Institute.  
 All rights reserved.
 Permission is hereby given to use this software 
 for non-profit educational purposes only.

*/
#ifndef __POLYNOMIAL_INTERNAL_H__
#define __POLYNOMIAL_INTERNAL_H__

#include <stdio.h>
#ifdef USE_GMP
#include <gmp.h>                /* GNU Multi-Precision library. */
#endif
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

/* Maximum size of a polynomial function name */
#define MAX_PFN_LEN 128
/* Maximum number of DLs supporting a single polynomial (up to 32K modules!) */
#define MAX_POLY_DL 32

/* If CPU utilization % is below this in a check interval, we're thrashing
   and should take some kind of action (currently exit) */
#define THRASH_CPU 5

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
#ifdef DEPENDENCYFLAGGING
  unsigned long dependentNodes;
#endif
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
  unsigned char eType;		// polynomial type:  - 1 byte
  double value;			// the value of the polynomial - 8 bytes
  unsigned int id;		// unique id - 4 bytes
  int index;			// index in a polynomial list - 4 bytes
  int key;			// key of the polynomial - 4 bytes
  unsigned short count;		// hold reference count - 2 bytes
  unsigned char valid;		// preservation flag(s) - 1 byte
#ifdef SOURCEDIGRAPH
  unsigned char source;		// index of entry in polySources
#endif
#ifdef USE_GMP
  mpf_t mpfValue;               // the Multi-Precision Float value of the polynomial - lotsa bytes
#endif
  union
  {
    struct variablePoly *v;	// variable
    struct sumPoly *s;		// sum
    struct productPoly *p;	// product
    struct functionPoly *f;	// function
    struct externalPoly *e;     // external
  } e;				// unused by constants - 8 bytes
#ifdef DEPENDENCYFLAGGING
  unsigned long dependencyFlag; // Flags up to 64 variable dependencies
#endif
  //  unsigned char oldEType;
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
#define MIN_CONSTANT_HASH_SIZE 199993
#define MIN_VARIABLE_HASH_SIZE 97
#define MIN_SUM_HASH_SIZE      199993
#define MIN_PRODUCT_HASH_SIZE  199993
#define MIN_FUNCTIONCALL_HASH_SIZE 9991
#define MIN_HASH_TABLE_INCREASE 2
#define MIN_CONSTANT_LIST_INITIAL 100000
#define MIN_CONSTANT_LIST_INCREASE 10000
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

int loadPolyDL (Polynomial * p);
Polynomial * exportPoly (Polynomial * p);
Polynomial * importPoly (Polynomial * p);
void importTermList (Polynomial * p);
void exportTermList (Polynomial * p, int writeFlag);
void deportTermList (Polynomial * p);
void thrashingCheck ();
#ifdef DEPENDENCYFLAGGING
void dependencyFlagging (Polynomial * p);
#endif

#endif
