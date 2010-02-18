/**
@file polynomial.c

 Polynomial - build and evaluate arbitrarily complex polynomials.

 Written by Hongling Wang. Polynomial reduction optimizations,
 reorganization, compilation, tuning and documentation by Bill 
 Valentine-Cooper.

 Copyright 2009, Nationwide Children's Research Institute.  All rights
 reserved.  Permission is hereby given to use this software for
 non-profit educational purposes only.

  Builds and evaluates arbitrarily complex polynomials. Collects
  subpolynomials of the same type (addition or multiplication) into
  a single tier. This process will cause the tree to be alternating
  layers of additions and multiplications. Identifies redundant
  polynomials and maintains only one copy with multiple references.

  Polynomial has been validated against GiNaC for computational
  accuracy (and performance, for our specific purposes).

  SIMPLIFIED USAGE:

  1. Call polynomialInitialization before any other routines.
  2. Create constants and variables with calls to constantExp and
  variableExp. Use the returned Polynomials in future calls.
  3. Call plusExp to build a new Polynomial by adding subpolynomials
  multiplied by double factors. Use returned Polynomial in future calls.
  4. Call timesExp to build a new Polynomial by multiplying
  subpolynomials raised to integer powers. Use returned Polynomial in
  future calls.
  5. Set the values of the referenced variables.
  6. Evaluate your Polynomial by calling evaluatePoly or evaluateValue.
  Use evaluatePoly if you're going to re-evaluate numerous times, as
  the cost of having to call polyListSorting pays for itself quickly.
  Use evaluateValue if you're only going to evaluate a few times.

  See the test driver pt.c for a very simple example of usage.

  LIMITATIONS:

  Maximum constant accuracy is 9 digits due to integer comparison.
  Total polynomials ever seen (kept or not) is INT_MAX due to nodeId.

  CONDITIONALS:

  There are several compilation conditionals in the polynomial code.

  - FREEDEBUG - define this to enable debugging diagnosis of situations
  where polynomials that have been freed by one mechanism or another
  end up being referenced later on. This leads to segmentation faults
  and bus errors, and is typically caused by mis-management of calls
  to freePolys(), freeKeptPolys(), and unHoldPoly(). FREEDEBUG causes
  the code that actually frees the polynomial structure to instead
  modify it's eType to T_FREED, and store it's original type in its
  value. Most polynomial handling functions recognize this and will
  dump structure information that can be used to track down the
  problem. The key bit of information is the nodeId, which is unique
  across all types of polynomials. This can be used in conjunction
  with the environment variable polynomialLostNodeId to track down
  usage of a particular nodeId (under gdb, it will breakout of
  execution when defined and when referenced).

  - POLYSTATISTICS - dump extensive polynomial statistics every raw 8Mp
  and at significant milestones.

   - THRASH_CHECK - enables minute-by-minute checks of CPU utilization
  and causes kelvin to abort if it is getting less than THRASH_CPU
  percentage CPU during that interval. The intent is to prevent
  runs that are thrashing from wedging the host OS.

  - SOURCEDIGRAPH - define this to enable the generation of DOT-format
  data files that illustrate term source parenting. The data is
  collected for all polynomial build calls, and dumped to a file when
  the function dumpSourceParenting is invoked. This was used to
  determine if a dynamic hybrid direct-evaluation /
  polynomial-evaluation approach could be implemented to trade memory
  for execution time.  Unfortunately there were too few good
  "breakpoints". It still serves to illustrate the complexity of the
  code and resultant polynomials.

  - EVALUATESW - define this to enable the tracking and display of
  polynomial evaluation statistics. It _does_ affect performance as
  the stopwatch has to be turned on and off with each evaluation, and
  there can be millions of evaluations.

  - DEPENDENCYFLAGGING - define this to enable pre-evaluatePoly dumps
  of the percentage of polynomial nodes that reference each of the
  variables (up to 64) used. This was the starting point for an idea
  that didn't pan-out for pedigree likelihood polynomials, but might
  for other applications. The remaining code is as yet unimplemented.
  The idea is to build a flag for each polynomial node that indicates
  which variables affect its value. This could be used to avoid re-
  evaluation of portions of the tree that are unaffected by changes in
  specific variables. Unfortunately the larger the likelihood polynomial,
  the more of the tree covered by the variables. It was anywhere from
  85 to 99%, so they payoff just wasn't there. Try turning this on if
  you have a different application and see if it is worth finishing.

  - USE_GMP - define this to enable validation of polynomial double
  arithmetic. 192-bit mantissa arithmetic will be performed in parallel
  with normal double arithmetic to validate results for evaluation of 
  non-compiled polynomials only (using evaluatePoly, NOT evaluateValue).
  Note that any operations performed outside of the polynomial as well as
  constant definitions and function evaluations do not use GMP, and so
  could still lead to degredation of precision. We implemented this when 
  we discovered that a lot of pedigree likelihood evaluation steps were 
  additions of numbers with exponents differing by more than 18 (mantissa 
  precision for a 64-bit double), which causes the smaller addend to be 
  effectively zero. We needed to be sure that we weren't dropping weak 
  linkage signals when this happened, so entire analyses were run this 
  way. All runs indicated that there was no effect -- consistent 
  precision was maintained out to 15 places, so we have turned off this 
  code until such time as precision questions arise again.

  - USE_SSD - use the sSDHandler.c routines to allow data structures to
  significantly exceed physical memory.

  - POLYUSE_DL - look for, load and use dynamic libraries built earlier
  for the named polynomials in this run.

  - POLYCODE_DL - source code generated by codePoly will include statements
  to allow building of dynamically-loaded libraries.

  - POLYCOMP_DL - source code generated by codePoly will actually be compiled
  thru invocation of the compileDL.sh script in KELVIN_ROOT.

  - POLYCHECK_DL - internal polynomials will be build whether compiled ones
  exist or not, and internal vs compiled results will be compared.


  ENVIRONMENT VARIABLES:

  - polynomialLostNodeId - used with FREEDEBUG conditional to track
  down inappropriately freed polynomials. See FREEDEBUG.

  - polynomialDebugLevel - used to control level of diagnostic
  output. Not referenced in recursive or otherwise intense code,
  so not a performance problem. Will probably be integrated into
  the older disused diagnostic routines used in other parts of
  kelvin.

  COMMENTS ON PHILOSOPHY:

  - The hash tables are maintained for one reason only: polynomial reuse. If
  they start to adversely affect performance without providing a significant
  reduction in the number of unique polynomials, then we need to clean them
  up or use them less. That's the case with constant polys for sure, almost
  enough so to lead us to not treat constants as polynomials.

*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <stdarg.h>
#include <time.h>
#include <signal.h>
#include <float.h>

#include <dlfcn.h>

#ifdef USE_GSL
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#endif

#ifndef POLYUSE_DL
  #if defined (POLYCODE_DL)| defined (POLYCOMP_DL) | defined (POLYCHECK_DL)
    #warning "Cannot use POLYCODE_DL, POLYCOMP_DL or POLYCHECK_DL without POLYUSE_DL"
  #endif
  #undef POLYCODE_DL
  #undef POLYCOMP_DL
  #undef POLYCHECK_DL
#endif

#ifndef POLYCODE_DL
  #if defined (POLYCOMP_DL) | defined (POLYCHECK_DL)
    #warning "Cannot use POLYCOMP_DL or POLYCHECK_DL without POLYCODE_DL"
  #endif
  #undef POLYCOMP_DL
  #undef POLYCHECK_DL
#endif

#ifndef POLYCOMP_DL
  #ifdef POLYCHECK_DL
    #warning "Cannot use POLYCHECK_DL without POLYCOMP_DL"
  #endif
  #undef POLYCHECK_DL
#endif

#ifdef USE_SSD
#include "sSDHandler.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "polynomial_internal.h"
#include "polynomial.h"

#undef constantExp
#undef plusExp
#undef timesExp

#include "dists.h"

#include "sw.h"

int polynomialScale;     ///< Scaling factor for hash and other storage set in config.c

#ifdef THRASH_CHECK
#define THRASH_CPU 10
extern struct swStopwatch *overallSW;   ///< Total run statistics stopwatch defined in kelvin.c
#endif

/** This is a global variable used for giving each polynomial an unique ID
   so that we can know if two polynomials are the same just from their IDs */

long nodeId; ///< @ingroup vettedGlobals

/** This global variable is for status monitoring, and indicates the raw number
    of terms processed from an external point of view. Note that simplification
    and redundancy handling make the actual internal number of terms much less. */

long termCount; ///< @ingroup vettedGlobals

/**

  We have a polynomials list to save the polynomials in each polynomial category,
  which is named constantList, variableList, sumList, productList, and functionCallList.
  The length of these lists are saved in constantCount, variableCount, productCount,
  functionCallCount respectively.  Each of the lists is dynamically applied since we are
  not sure how many polynomials we may have.  Therefore, we apply to have a short list
  initially and increase its length gradually with the increase of the number of
  polynomials. The current length of the lists is recorded in constantListLength,
  variableListLength, sumListLength, productListLength, and functionCallListLength
  respectively.

  Why not free all the polynomials before building the likelihood polynomials for the
  next trait position? Before we even start to build the likelihood polynomials of
  pedigrees at the first trait position, some polynomials have already been built.
  These polynomials are useful for all the trait positions.  Therefore, we can't free
  everything before we start to build likelihood polynomials of pedigrees at the next
  trait position.

*/
static struct polynomial **constantList;
int constantCount;      ///< The number of constants in constantList @ingroup vettedGlobals
static int constantListLength; ///< The total space available for constants in constantList

static struct polynomial **variableList;       ///< The list of all variables
int variableCount;      ///< The number of variables in variableList @ingroup vettedGlobals
static int variableListLength; ///< The total space available for variables in variableList

struct polynomial **externalList;       ///< The list of all external functions
int externalCount;      ///< The number of externals in externalList @ingroup vettedGlobals
int externalListLength; ///< The total space available for externals in externalList

static struct polynomial **sumList;    ///< The list of all sums
int sumCount;   ///< The number of sums in sumList @ingroup vettedGlobals
static int sumListLength;      ///< The total space available for sums in sumList

static struct polynomial **productList;
int productCount; ///< The number of products in productList @ingroup vettedGlobals
static int productListLength;

static struct polynomial **functionCallList;
int functionCallCount; ///< The number of functions in functionCallList @ingroup vettedGlobals
static int functionCallListLength;

static char *eTypes[8] = { "C", "V", "S", "P", "F", "E", "U", "O" }; ///< Useful prefixes for polynomial types

/** @defgroup PolyCons Polynomial Scaling Constants
    @{
  Constants (not really) for static size of hashes and initial size of lists.
  These used to be constants, hence their case. Now they're scaled, but I kept the
  case because they act like constants. */

static int CONSTANT_HASH_SIZE, VARIABLE_HASH_SIZE, SUM_HASH_SIZE, PRODUCT_HASH_SIZE,
  FUNCTIONCALL_HASH_SIZE, HASH_TABLE_INCREASE,
  CONSTANT_LIST_INITIAL, CONSTANT_LIST_INCREASE,
  VARIABLE_LIST_INITIAL, VARIABLE_LIST_INCREASE, 
  EXTERNAL_LIST_INITIAL, EXTERNAL_LIST_INCREASE, 
  SUM_LIST_INITIAL, SUM_LIST_INCREASE, 
  PRODUCT_LIST_INITIAL, PRODUCT_LIST_INCREASE, 
  FUNCTIONCALL_LIST_INITIAL, FUNCTIONCALL_LIST_INCREASE = 0;
/*@}*/

/** @defgroup PolyPerf Polynomial Performance Monitoring Variables
    @{
  Variables for tracking internal polynomial memory usage. They are
  also externs in polynomial.h so you can monitor them from wherever
  you like. */

int maxHashLength = 0;  ///< Length of longest hash table collision list
int constantHashHits = 0,       ///< constantExp return is a pre-existing polynomial
    variableHashHits = 0,       ///< variableExp return is a pre-existing polynomial (surprise!)
    functionHashHits = 0;       ///< functionCallExp return is a pre-existing polynomial
int sumReleaseableCount = 0,    ///< Indicates if flag was set to release 1st term in sum
    sumNotReleaseableCount,     ///< ...or not.
    sumReturnConstantCount,     ///< plusExp return is actually a constant
    sumReturnSameConstantCount,     ///< plusExp return is actually the same constant as the 1st term
    sumReturn1TermCount,        ///< plusExp return is single-term polynomial
    sumHashHits,        ///< plusExp return is pre-existing polynomial
    sumNewCount,        ///< plusExp return is a new polynomial
    sumListNewCount,    ///< New sumList entry created
    sumListReplacementCount,    ///< Existing sumList entry used
    sumFreedCount = 0,  ///< Count of sumPolys freed
    sum1stTermsFreedCount = 0;  ///< Count of 1st-terms successfully freed
int productReleaseableCount = 0,        ///< Indicates if flag was set to release 1st term in product
    productNotReleaseableCount = 0,     ///< ...or not.
    productReturn0Count = 0,    ///< timesExp return is actually zero
    productReturnConstantCount = 0,     ///< timesExp return is actually a non-zero constant
    productReturnSameConstantCount = 0,     ///< timesExp return is actually the same constant as the 1st term
    productReturn1stTermCount = 0,      ///< timesExp return is actually the first term itself
    productReturn1TermSumCount = 0,     ///< timesExp return is a sum built from the first term itself
    productHashHits = 0,        ///< timesExp return is a pre-existing polynomial
    productHashHitIsSumCount = 0,       ///< timesExp return is a sum built from a pre-existing polynomial
    productReturnNormalCount = 0,       ///< timesExp return is a new product polynomial (not a sum)
    productNon1FactorIsSumCount = 0,    ///< timesExp return is a new sum of a new product
    productListNewCount = 0,    ///< New productList entry created
    productListReplacementCount = 0,    ///< Existing productList entry used
    productFreedCount = 0,      ///< Count of productPolys freed
    product1stTermsFreedCount = 0;      ///< Count of 1st-terms successfully freed
int constantPListExpansions = 0,        ///< Count of constantList expansions
    variablePListExpansions = 0,        ///< Count of variableList expansions
    sumPListExpansions = 0,     ///< Count of sumList expansions
    productPListExpansions = 0, ///< Count of productList expansions
    functionCallPListExpansions = 0;    ///< Count of functionCallList expansions
/// Count of calls to various statistically interesting functions
int polyListSortingCount = 0, evaluatePolyCount = 0, evaluateValueCount = 0, keepPolyCount = 0, freePolysCount = 0, holdPolyCount = 0, holdAllPolysCount = 0, unHoldPolyCount = 0, freeKeptPolysCount = 0, freePolysAttemptCount = 0;
int containerExpansions = 0;    ///< Count of expansions of any term-collection container.
unsigned long totalSPLLengths = 0, totalSPLCalls = 0, lowSPLCount = 0, highSPLCount = 0;
unsigned long initialHashSize = 0;      ///< Total initial size of hash table and collision lists
int pendingExplicitDiscards = 0; ///< Number of first-term explicit discards awaiting removal
/*@}*/

/// Flag set whenever we calculate our total data storage estimate
static int memoryLow = FALSE;

/// The svn version of this module displayed by kelvin @ingroup vettedGlobals
char *polynomialVersion = "$Id$";

/** Both of the following are set by initialization to value of environment variable of same name.
  They control diagnostic action in a manner not permitted by other approaches since they can
  be changed without rebuilding. */

static int polynomialDebugLevel = 0;   ///< Corresponds roughly to diagnostic output volume
static int polynomialLostNodeId = -1;  ///< For tracking down mis-freed polynomials

static struct swStopwatch *evaluatePolySW,     ///< Conditional stopwatch for evaluatePoly stats.
 *evaluateValueSW;      ///< Conditional stopwatch for evaluateValue stats.

#ifdef THRASH_CHECK
static unsigned long lastPDSAccumWallTime = 0, ///< last polyDynamicStatistics walltime for thrashing check.
    lastPDSAccumUserTime = 0;   ///< ...same but CPU time.
#endif

// The following are used for building sum polynomials...
/// Collecting variable terms...
static double *factor_v1;
static struct polynomial **p_v1;
static int containerLength_v1;
static int counter_v1;
// Collecting product terms...
static double *factor_p1;
static struct polynomial **p_p1;
static int containerLength_p1;
static int counter_p1;
// Collecting function call terms...
static double *factor_f1;
static struct polynomial **p_f1;
static int containerLength_f1;
static int counter_f1;
// Collecting all terms...
static double *factorSum;
static struct polynomial **pSum;
static int lengthSum;

// Used for building product polynomials...
// Collecting variable terms
static int *exponent_v2;
static struct polynomial **p_v2;
static int containerLength_v2;
static int counter_v2;
// Collecting sum terms...
static int *exponent_s2;
static struct polynomial **p_s2;
static int containerLength_s2;
static int counter_s2;
// Collecting function call terms...
static int *exponent_f2;
static struct polynomial **p_f2;
static int containerLength_f2;
static int counter_f2;
// Collecting all terms...
static int *exponentProd;
static struct polynomial **pProd;
static int lengthProd;

// Hash tables
static struct hashStruct *constantHash;        //... for constant polynomials
static struct hashStruct *variableHash;        //... for variable polynomials
static struct hashStruct *sumHash;     //...for sum polynomials
static struct hashStruct *productHash; //...for the product polynomials
static struct hashStruct *functionCallHash;    //...for the functionCall polynomials

#ifdef SOURCEDIGRAPH
struct polySource
{
  int entryNo;
  int lineNo;
  char moduleName[255];
  unsigned char eType;
  int totalCalls;
  int originalChildren[MAXPOLYSOURCES];
};

static struct polySource  polySources[MAXPOLYSOURCES];
static int polySourceCount = 0;
static int originalChildren[MAXPOLYSOURCES];
#endif

#ifdef USE_SSD
static int importedTerms = 0, exportedDroppedTerms = 0, exportedWrittenTerms = 0, peakInMemoryTerms = 0;
static long sSDDebt = 0; // When <= 0, we're OK, if not, must pay-off allocation debt with freeSSD calls
#define MAX_SSD_BUFFER ((MAX_DPC_MASK+1)*2)

#define IMTL_COUNT 32
struct inMemoryTermList {
  unsigned int lastNodeId; // This will enable recycling without retrieval
  struct chunkTicket cT;
  double buffer[MAX_SSD_BUFFER];
};
static struct inMemoryTermList  iMTL[IMTL_COUNT];
#endif



/**

  This was used for determining polynomial size, but it could not
  consider reuse. Now it's kept only because it's a debugged way
  of doing post-processing on all returned polynomials should that
  be useful in the future.

*/
Polynomial *polyReturnWrapper (Polynomial * p)
{
  // For testing purposes, we're going to export every polynomial!
  //  return exportPoly (p);
  return (p);
}

/**

  Clear the evaluation flag on the entire tree so we can mark where
  we've been and not retrace our steps regardless of redundancy.

  This code has been adapted for parallel execution by breaking-up
  the iteration over constants, sums and products.

*/
void clearValidEvalFlag ()
{
  int i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < constantCount; i++)
    constantList[i]->valid &= ~VALID_EVAL_FLAG;
  for (i = 0; i < variableCount; i++)
    variableList[i]->valid &= ~VALID_EVAL_FLAG;
  for (i = 0; i < externalCount; i++)
    externalList[i]->valid &= ~VALID_EVAL_FLAG;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < constantCount; i++)
    constantList[i]->valid &= ~VALID_EVAL_FLAG;
  for (i = 0; i < variableCount; i++)
    variableList[i]->valid &= ~VALID_EVAL_FLAG;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < sumCount; i++)
    sumList[i]->valid &= ~VALID_EVAL_FLAG;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < productCount; i++)
    productList[i]->valid &= ~VALID_EVAL_FLAG;
  for (i = 0; i < functionCallCount; i++)
    functionCallList[i]->valid &= ~VALID_EVAL_FLAG;
  return;
}

/**

  Recursively evaluate a polynomial without the benefit of a
  polynomial list.  This version of polynomial evaluation doesn't use
  the polynomial sorting list, so if you're just doing a few
  evaluations, this can be faster.

*/
double doEvaluateValue (Polynomial * p)
{
  int i;
  double result;
  struct sumPoly *sP;
  struct productPoly *pP;
  struct functionPoly *fp;
  double value0, value1;

  if (p->valid & VALID_EVAL_FLAG)
    return p->value;

  switch (p->eType) {
  case T_OFFLINE:
    importPoly(p);
    // Notice that there's no break here...

    // If a sub polynomial is a contant, return the value
  case T_CONSTANT:
    p->valid |= VALID_EVAL_FLAG;
    return p->value;
    // If a sub polynomial is a variable, return the value
  case T_VARIABLE:
    p->valid |= VALID_EVAL_FLAG;
    if (p->e.v->vType == 'D') {
      p->value = *(p->e.v->vAddr.vAddrD);
      return *(p->e.v->vAddr.vAddrD);
    } else if (p->e.v->vType == 'I') {
      p->value = *(p->e.v->vAddr.vAddrI);
      return *(p->e.v->vAddr.vAddrI);
    } else
      FATAL ("Unknown variable type '%c'", p->e.v->vType);

  case T_EXTERNAL:
#ifdef POLYUSE_DL
    if (!p->e.e->entryOK)
      if (!loadPolyDL (p))
	FATAL ("Cannot (re)load DL-based polynomial %s for evaluateValue", 
		 p->e.e->polynomialFunctionName);
#endif
    p->value = p->e.e->polynomialFunctionRoutine (1, variableList);
    return p->value;

    /* If a sub polynomial is a sum, evaluate the values of all the terms.
     * Add the values of the terms weighted by their coefficients. */
  case T_SUM:
    result = 0;
    importTermList(p);
    sP = p->e.s;
    for (i = 0; i < sP->num; i++) {
      if (sP->factor[i] == 1)
        result += doEvaluateValue (sP->sum[i]);
      else
        result += doEvaluateValue (sP->sum[i]) * sP->factor[i];
    }
    exportTermList(p, FALSE);
    p->value = result;
    p->valid |= VALID_EVAL_FLAG;
    return result;

    /* If a sub polynomial is a product, evaluate the values of all the terms.
     * Multiply the values of the terms powered by their exponents */
  case T_PRODUCT:
    result = 1;
    pP = p->e.p;
    for (i = 0; i < pP->num; i++) {
      switch (pP->exponent[i]) {
      case 0:
        break;
      case 1:
        result *= doEvaluateValue (pP->product[i]);
        break;
      default:
        result *= pow (doEvaluateValue (pP->product[i]), pP->exponent[i]);
        break;
      }
    }
    p->value = result;
    p->valid |= VALID_EVAL_FLAG;
    return result;

    /* If a sub polynomial is a function call, evaluate the values of all the parameters
     * and then the function according to the name of the function.  Therefore, a function
     * must have a unique name and a unique set of parameters. */
  case T_FUNCTIONCALL:
    fp = p->e.f;
    // log10 is for LOD score computation
    if (strcmp (fp->name, "log10") == 0)
      result = log10 (doEvaluateValue (fp->para[0]));
    else if (strcmp (fp->name, "log") == 0)
      result = log (doEvaluateValue (fp->para[0]));
    else if (strcmp (fp->name, "tanh") == 0)
      result = tanh (doEvaluateValue (fp->para[0]));
    else if (strcmp (fp->name, "atanh") == 0)
      result = atanh (doEvaluateValue (fp->para[0]));
    /* gsl_ran_tdist_pdf,gsl_cdf_tdist_Q,gsl_cdf_tdist_P,
     * gsl_ran_ugaussian_pdf,gsl_cdf_ugaussian_Q,gsl_cdf_ugaussian_P, and
     * gsl_cdf_chisq_P, gsl_cdf_chisq_Q, and gsl_ran_chisq_pdf
     * are for quantitative trait gene's penetrance computation */
    else if (strcmp (fp->name, "gsl_ran_tdist_pdf") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      //      result = gsl_ran_tdist_pdf (value0, value1);
      result = t_pdf_30 (value0, value1);
    } else if (strcmp (fp->name, "gsl_cdf_tdist_Q") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      //      result = gsl_cdf_tdist_Q (value0, value1);
      result = ((double) 1.0) - t_cdf (value0, value1);
    } else if (strcmp (fp->name, "gsl_cdf_tdist_P") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      //      result = gsl_cdf_tdist_P (value0, value1);
      result = t_cdf (value0, value1);
    } else if (strcmp (fp->name, "gsl_ran_ugaussian_pdf") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      //      result = gsl_ran_ugaussian_pdf (value0);
      result = gaussian_pdf (value0, (double) 0.0, (double) 1.0);
    } else if (strcmp (fp->name, "gsl_cdf_ugaussian_Q") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      //      result = gsl_cdf_ugaussian_Q (value0);
      result = ((double) 1.0) - gaussian_cdf (value0, (double) 0.0, (double) 1.0);
    } else if (strcmp (fp->name, "gsl_cdf_ugaussian_P") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      //      result = gsl_cdf_ugaussian_P (value0);
      result = gaussian_cdf (value0, (double) 0.0, (double) 1.0);
    } else if (strcmp (fp->name, "gsl_cdf_chisq_P") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      //      result = gsl_cdf_chisq_P (value0, value1);
      result = chisq_cdf (value0, value1);
    } else if (strcmp (fp->name, "gsl_cdf_chisq_Q") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      //      result = gsl_cdf_chisq_Q (value0, value1);
      result = ((double) 1.0) - chisq_cdf (value0, value1);
    } else if (strcmp (fp->name, "gsl_ran_chisq_pdf") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      //      result = gsl_ran_chisq_pdf (value0, value1);
      result = chisq_pdf (value0, value1);
    }
    // pow, exp, sqrt are standard functions
    else if (strcmp (fp->name, "pow") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      result = pow (value0, value1);
    } else if (strcmp (fp->name, "exp") == 0) {
      result = exp (doEvaluateValue (fp->para[0]));
    } else if (strcmp (fp->name, "sqrt") == 0) {
      result = sqrt (doEvaluateValue (fp->para[0]));
    } else
      ERROR ("Unknown function name %s in polynomial\n", fp->name);
    p->value = result;
    p->valid |= VALID_EVAL_FLAG;
    return result;

  default:
    FATAL ("In evaluateValue, unknown expression type: [%d]", p->eType);
  }
}

/**

  Single point-of-entry for evaluateValue. Does housekeeping and calls
  doEvaluateValue.

  This code has been adapted for parallel execution by stepping down
  one level in the tree and starting a separate thread for each of the
  subpolys at that point. After they finish, the root of the tree is
  processed to ensure that everything is done.

*/
double evaluateValue (Polynomial * p)
{
  double returnValue;
  int i;

  if ((evaluateValueCount++ & 0x3FF) == 0)
    DETAIL (0, "Evaluating polynomials, currently at %1.2g iterations", (double) evaluateValueCount);

#ifdef EVALUATESW
  swStart (evaluateValueSW);
#endif

  /* Clear all of the VALID_EVAL_FLAGs */
  clearValidEvalFlag ();

  if (p->eType == T_OFFLINE) importPoly (p);
  if (p->eType == T_SUM || p->eType == T_PRODUCT) {
    /* Step down a level where there are enough terms to keep the work interesting. */

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < p->e.s->num; i++) {
      doEvaluateValue (p->e.s->sum[i]);
    }
  }
  returnValue = doEvaluateValue (p);
#ifdef EVALUATESW
  swStop (evaluateValueSW);
#endif

  return returnValue;
}

/**

  A polynomial is zero if it is a constant polynomial and its value is 0.

*/
inline int isZeroExp (Polynomial * p)
{
  if (p->eType == T_OFFLINE) importPoly (p);
  if (p->eType != T_CONSTANT)
    return 0;
  if (p->value == 0.0)
    return 1;
  else
    return 0;
};

/**

  Binary search in an array of integers. Return 1 if found, 0 if
  not. location is the position for insertion of the target into the
  array. The polynomials are searched based on their keys.

*/
inline int binarySearch (int *array, int length, int target, int *location)
{
  int binaryStart, binaryEnd, binaryMiddle;

  if (length == 0) {
    *location = 0;
    return 0;
  }
  // Start from the two ends of the array, and bisect each time
  binaryStart = 0;
  binaryEnd = length - 1;
  while (binaryStart <= binaryEnd) {
    binaryMiddle = floor ((binaryStart + binaryEnd) / 2.0);
    //the target is found in the array
    if (target == array[binaryMiddle]) {
      *location = binaryMiddle;
      return 1;
    }
    // Target is not found in the array
    else if (target > array[binaryMiddle])
      binaryStart = binaryMiddle + 1;
    else
      binaryEnd = binaryMiddle - 1;
  }
  *location = binaryStart;
  return 0;
};

/**

  Search for a polynomial in the hash table with a key.  Because
  multiple polynomials may have the same key, the search in a hash
  table may return a list of polynomials that have the same key as the
  key of the target polynomial.

*/
inline int searchHashTable (struct hashStruct *hash, int *start, int *end, int key)
{
  int found;
  int location;

  if (binarySearch (hash->key, hash->num, key, &location) == 1) {
    *end = location + 1;
    *start = location - 1;
    while ((*start) >= 0 && hash->key[*start] == hash->key[location])
      (*start)--;
    while ((*end) < hash->num && hash->key[*end] == hash->key[location])
      (*end)++;
    if ((*start) < 0 || key != hash->key[*start])
      (*start)++;
    if ((*end) >= hash->num || key != hash->key[*end])
      (*end)--;
    found = 1;
  } else {
    *start = location;
    *end = location;
    found = 0;
  }
  return found;
};

/**
  Store a polynomial in a hash table. Polynomials in the same hash table index are
  sorted by their keys in increasing order.

*/
inline void insertHashTable (struct hashStruct *hash, int location, int key, int ourIndex)
{
  hash->num++;
  // If the hash collision list is full, allocate more memory
  if (hash->num > hash->length) {
    hash->length += HASH_TABLE_INCREASE;
    if (hash->length > maxHashLength)
      maxHashLength = hash->length;
    hash->key = realloc (hash->key, sizeof (int) * hash->length);
    hash->index = realloc (hash->index, sizeof (int) * hash->length);
    initialHashSize += HASH_TABLE_INCREASE * sizeof (int) * 2;
    if (hash->key == NULL || hash->index == NULL)
      FATAL ("Memory allocation for hash table failed!");
  }
  // Prepare a place in hash table for the new polynomial
  if (location <= hash->num - 2) {
    memmove (&hash->key[location + 1], &hash->key[location], sizeof (int) * (hash->num - location - 1));
    memmove (&hash->index[location + 1], &hash->index[location], sizeof (int) * (hash->num - location - 1));
  }
  //Insert the new polynomiual in the hash table
  hash->key[location] = key;
  hash->index[location] = ourIndex;

};

/// Delete the polynomial designated by location from a hash table
inline void deleteHashTable (struct hashStruct *hash, int location)
{
  if (location >= hash->num)
    FATAL ("Deletion in hash table failed!");
  memmove (&hash->key[location], &hash->key[location + 1], sizeof (int) * (hash->num - location - 1));
  memmove (&hash->index[location], &hash->index[location + 1], sizeof (int) * (hash->num - location - 1));
  hash->num--;
};

/**

  Compute the key of a constant polynomial from the normalized fraction and the exponent
  generated from function frexp, which converts a floating-point number to fractional
  and integral components.

*/
inline int keyConstantPolynomial (double tempD1, int tempI1)
{
  int key;

  // The key is a number between 0 and MAX_POLYNOMIAL_KEY
  if (tempD1 > 0)
    key = floor ((tempD1 - 0.5) * 2 * MAX_POLYNOMIAL_KEY) + tempI1;
  else
    key = floor ((tempD1 + 0.5) * 2 * MAX_POLYNOMIAL_KEY) + tempI1;
  return key;
};

/**

  Construct a constant polynomial. If one already exists, no new one is created, instead
   the existing constant polynomial is returned.

*/
Polynomial *constantExp (char *fileName, int lineNo, double con)
{
  int i;
  Polynomial *p;        // Newly built constant polynomial
  int key;      // Key of the new constant polynomial
  int hIndex, cIndex;   // Index in the hash table and the constant polynomial list
  int first, last, location;    /* First and last polynomial in a hash bucket for comparison,
                                 * and location of the newly-built polynomial in the bucket */
  int tempI1 = -1, tempI2 = -1; // Normalized integer component
  double tempD1, tempD2;        // Normalized fractional component

#ifdef SOURCEDIGRAPH
  memset (originalChildren, 0, sizeof (originalChildren));
#endif

#ifdef POLYSTATISTICS
  if ((termCount++ & 0x7FFFFF) == 0) ///< New term from an external standpoint.
    polyStatistics ("At 8M raw term multiple");
#endif

  // Compute a key for the constant
  tempD1 = frexp (con, &tempI1);
  key = keyConstantPolynomial (tempD1, tempI1);

  // Compute the index of this constant polynomial in the hash table of constant polynomials
  hIndex = key % CONSTANT_HASH_SIZE;
  if (hIndex < 0)
    hIndex += CONSTANT_HASH_SIZE;

  /* If the hash table item is not empty, determine if the constant is in the
   * constant polynomial list.  If it is not, find a position in the hash table to
   * save the key and the index in the constant list */
  if (constantHash[hIndex].num > 0) {
    /* If the key of this constant is equal to the keys of some polynomials in the
     * constant list, see if the value of the new constant is equal to the value
     * of an existing constant */
    if (searchHashTable (&constantHash[hIndex], &first, &last, key)) {
      for (i = first; i <= last; i++) {
        cIndex = constantHash[hIndex].index[i];
        tempD2 = frexp (constantList[cIndex]->value, &tempI2);
        // See if the two constants are the same
        if (tempI2 == tempI1 && (int) (tempD2 * 100000000) == (int) (tempD1 * 100000000)) {
          // If the two constants are the same, return the constant in the constant list
	  // Only set one or the other of VALID_REF_FLAG and VALID_NOTDISC_FLAG
	  if (constantList[cIndex]->valid & VALID_NOTDISC_FLAG)
	    // Not a pending discard, so now it's referenced
	    constantList[cIndex]->valid |= VALID_REF_FLAG;
	  else {
	    // Pending discard is recycled, so not pending anymore (and not referenced!)
	    constantList[cIndex]->valid |= VALID_NOTDISC_FLAG;
	    pendingExplicitDiscards--;
	  }
          constantHashHits++;
          return polyReturnWrapper (constantList[cIndex]);
        }
      }
    }
    location = last;
  } else
    location = 0;

  // Next, insert it into the constant list

  // Generate a constant polynomial
  MALCHOKE(p, sizeof (Polynomial), Polynomial *);
  p->eType = T_CONSTANT;
  p->value = con;
#ifdef USE_GMP
  mpf_init_set_d (p->mpfValue, p->value);
#endif

  // Check if the constant polynomial list is full.
  if (constantCount >= constantListLength) {
    constantListLength += CONSTANT_LIST_INCREASE;
    constantPListExpansions++;
    REALCHOKE (constantList, constantListLength * sizeof (Polynomial *), void *);
  }
  // Save the constant in the constant polynomial list
  constantList[constantCount] = p;
  p->index = constantCount;
  // Give the polynomial a unique ID
  p->id = nodeId;
  constantCount++;
  //  if (constantCount >= 100000)
  //    printf ("Hit 100K\n");
  nodeId++;
  p->key = key;
  p->valid = VALID_NOTDISC_FLAG;
  p->count = 0;
#ifdef SOURCEDIGRAPH
  p->source = findOrAddSource (fileName, lineNo, T_CONSTANT);
#endif

  // Record the constant polynomial in the hash table of constant polynomials
  insertHashTable (&constantHash[hIndex], location, key, constantCount - 1);

  return polyReturnWrapper (p);
};

/**

  Generate a key for a variable polynomial.

  The key of a variable polynomial is computed by the address of the variable.
  Different variables have different addresses in memory which provide nicely
  random keys.

*/
inline int keyVariablePolynomial (double *vD, int *vI, char vType)
{
  int key;

  /* Currently we can deal with integer and double variables. We use the address 
   * of the variable as the key of the variable polynomial. */
  if (vType == 'D')
    key = (long int) vD;
  else if (vType == 'I')
    key = (long int) vI;
  else
    FATAL ("Unknown variable type '%c'", vType);

  return key;
};

/**

  Builds a variable polynomial. Essentially just a name and a place to find
  it's value before evaluation.

*/

Polynomial *variableExp (double *vD, int *vI, char vType, char name[10])
{
  int i;
  Polynomial *p;        // Newly built variable polynomial
  struct variablePoly *vPoly;   // Variable structure for newly-built polynomial
  int key;      // Key of the newly-built variable polynomial
  int hIndex, vIndex;   // Index in the hash table and polynomial list
  int first, last, location;

#ifdef POLYSTATISTICS
  if ((termCount++ & 0x7FFFFF) == 0) ///< New term from an external standpoint.
    polyStatistics ("At 8M raw term multiple");
#endif

  // Create a key for this variable polynomial
  key = keyVariablePolynomial (vD, vI, vType);

  // Compute the index of this variable poly in the hash table of variable polynomials
  hIndex = key % VARIABLE_HASH_SIZE;
  if (hIndex < 0)
    hIndex += VARIABLE_HASH_SIZE;
  /* If the hash table item is not empty, determine if the variable is already in the
   * variable polynomial list.  If it is not, determine a position in the hash table to
   * save the key and an index in the variable list. Save this variable polynomial in
   * the variable polynomial list. */
  if (variableHash[hIndex].num > 0) {
    /* If the key of this variable is equal to the keys of some polynomials in the
     * variable list, check if the variable is already there. */
    if (searchHashTable (&variableHash[hIndex], &first, &last, key)) {
      for (i = first; i <= last; i++) {
        vIndex = variableHash[hIndex].index[i];
        // Compare if the two variables are the same
        if ((vType == 'D' && variableList[vIndex]->e.v->vAddr.vAddrD == vD)
            || (vType == 'I' && variableList[vIndex]->e.v->vAddr.vAddrI == vI)) {
          // If the two variables are the same, return the address from the list
          variableList[vIndex]->valid |= VALID_REF_FLAG;
          variableHashHits++;
          return polyReturnWrapper (variableList[vIndex]);
        }
      }
    }
    location = last;
  } else
    location = 0;

  // This variable polynomial doesn't exist, so we'll create it.
  MALCHOKE(p, sizeof (Polynomial), Polynomial *);

  MALCHOKE(p, sizeof (Polynomial), Polynomial *);
  MALCHOKE(vPoly, sizeof (struct variablePoly), struct variablePoly *);
  p->eType = T_VARIABLE;
  // Make sure we always have a name, either provided or based upon arrival order.
  if (strlen (name) == 0)
    sprintf (vPoly->vName, "u%d", variableCount);
  else
    strcpy (vPoly->vName, name);
  if (vType == 'D')
    vPoly->vAddr.vAddrD = vD;
  else
    vPoly->vAddr.vAddrI = vI;
  vPoly->vType = vType;
  p->e.v = vPoly;

  // If the polynomial list is full, get  more memory
  if (variableCount >= variableListLength) {
    variableListLength += VARIABLE_LIST_INCREASE;
    variablePListExpansions++;
    REALCHOKE (variableList, variableListLength * sizeof (Polynomial *), void *);
  }

  p->index = variableCount;
  p->id = nodeId;
  p->key = key;
  p->valid = VALID_NOTDISC_FLAG;
  p->count = 0;
#ifdef USE_GMP
  mpf_init (p->mpfValue);
#endif

  // Insert the variable polynomial in the variable polynomial list
  variableList[variableCount] = p;
  variableCount++;
  nodeId++;

  // Record the variable polynomial in the hash table of the variable polynomials
  insertHashTable (&variableHash[hIndex], location, key, variableCount - 1);

  return polyReturnWrapper (p);
};

/**

   Compute a key for a sum polynomial. There tends to be a huge number of sum
   and product polynomials in likelihood calculations, so it is important that we have a
   good hash strategy so that the polynomials are well distributed in the hash table.

*/
inline int keySumPolynomial (Polynomial ** p, double *factor, int counter)
{
  int key = 0;
  int j;
  double tempD1;
  int tempI1;

  /* Here we use the ids and keys of all terms in a sum polynomial, their coefficients, their
   * order in the sum to compute a key for the sum polynomial. */
  for (j = 0; j < counter; j++) {
    tempD1 = frexp (factor[j], &tempI1);
    if (p[j]->eType == T_OFFLINE) importPoly (p[j]);
    key += p[j]->id + p[j]->key * (j + 1) + (int) (tempD1 * 12345 + tempI1) * ((int) p[j]->eType + 1);
  }
  return key;
};

/*

  Search for a polynomial in a polynomial list. If the target polynomial is not in the list,
  return 0, otherwise return 1 and it's index in the location parameter. This function has
  been heavily optimized for degenerate cases due to it's intensive use.

*/
inline int searchPolynomialList (Polynomial ** p, int length, Polynomial * target, int *location)
{
  int binaryStart, binaryEnd, binaryMiddle;

  totalSPLLengths += length;
  totalSPLCalls++;

  if (length == 0) {
    lowSPLCount++;
    *location = 0;
    return 0;
  }
  if (target->index > p[length - 1]->index) {
    highSPLCount++;
    *location = length;
    return 0;
  }

  binaryStart = 0;
  binaryEnd = length - 1;
  while (binaryStart <= binaryEnd) {
    binaryMiddle = floor ((binaryStart + binaryEnd) / 2.0);
    if (target->index == p[binaryMiddle]->index) {
      *location = binaryMiddle;
      return 1;
    } else if (target->index > p[binaryMiddle]->index)
      binaryStart = binaryMiddle + 1;
    else
      binaryEnd = binaryMiddle - 1;
  }
  *location = binaryStart;
  return 0;
};

/**

  Collect the terms and their coefficients of a sum polynomial.  A plus operation can
  have any number of operands which are themselves polynomials.  We collect these operand
  polynomials and organize them in a unique manner as a basis for comparisons between the
  resulting polynomial and other polynomials.

*/
inline void collectSumTerms (double **factor,   ///< Container, pointer to a list of double factors
    Polynomial *** p,   ///< Container, pointer to a list of pointers to polynomials
    int *counter,       ///< Number of terms currently in the container lists
    int *containerLength,       ///< Current size of the container lists
    double f1,  ///< Factor of term to add to the container
    Polynomial * p1     ///< Pointer to polynomial of term to add to container
    )
{
  int location;

  // Search for the position where the new item should be inserted
  if (searchPolynomialList (*p, *counter, p1, &location) == 1) {
    // This item is currently in the sum, just merge their coefficients
    (*factor)[location] += f1;
    return;
  }
  // If container is full, apply for more memory
  if (*counter >= *containerLength - 1) {
    (*containerLength) += 50;
    containerExpansions++;
    REALCHOKE (*factor, (*containerLength) * sizeof (double), double *);
    REALCHOKE (*p, (*containerLength) * sizeof (Polynomial *), Polynomial **);
  }
  // This is a new item in the sum, insert it at the start or end
  if (location >= *counter) {
    (*p)[*counter] = p1;
    (*factor)[*counter] = f1;
    (*counter)++;
  }
  // Make space so we can insert the item where it goes in the middle of the sum
  else {
    memmove (&((*p)[location + 1]), &((*p)[location]), sizeof (Polynomial *) * ((*counter) - location));
    memmove (&((*factor)[location + 1]), &((*factor)[location]), sizeof (double) * ((*counter) - location));

    // Insert the new item
    (*p)[location] = p1;
    (*factor)[location] = f1;
    (*counter)++;
  }
}

#ifdef SOURCEDIGRAPH
/**

  qsort-referenced function to put sources in name order for source digraph.

*/
int compareSourcesByName (const void *left, const void *right)
{
  int result;
  struct polySource *pSLeft, *pSRight;

  pSLeft = (struct polySource *) left;
  pSRight = (struct polySource *) right;
  if ((result = strcmp (pSLeft->moduleName, pSRight->moduleName)) == 0) {
    result = pSLeft->lineNo - pSRight->lineNo;
  }
  return result;
}

/**

  qsort-referenced function to put sources in entryNo order for source digraph.

*/
int compareSourcesByEntryNo (const void *left, const void *right)
{
  struct polySource *pSLeft, *pSRight;

  pSLeft = (struct polySource *) left;
  pSRight = (struct polySource *) right;
  return (pSLeft->entryNo - pSRight->entryNo);
}

/**

  Find and return a polynomial term's source or add a new one.

*/
short findOrAddSource (char *fileName, int lineNo, unsigned char eType)
{
  struct polySource target, *result;
  int i;

  strcpy (target.moduleName, fileName);
  target.lineNo = lineNo;
  result = bsearch (&target, polySources, polySourceCount, sizeof (struct polySource), compareSourcesByName);
  if (result == NULL) {
    if (polySourceCount >= MAXPOLYSOURCES) {
      WARNING ("Exceeded maximum polynomial source count, no more locations can be monitored");
      return 0;
    }
    result = &polySources[polySourceCount];
    strcpy (result->moduleName, fileName);
    result->lineNo = lineNo;
    result->entryNo = polySourceCount;
    qsort (polySources, ++polySourceCount, sizeof (struct polySource), compareSourcesByName);
  }
  result->eType = eType;
  result->totalCalls++;
  for (i = 0; i < MAXPOLYSOURCES; i++) {
    result->originalChildren[i] += originalChildren[i];
  }
  return (result->entryNo);
}
#endif


/**

  Dump the specified polynomial in terms of constants, variables and operations.
  This produces a crazy amount of output if you're not careful.

*/
void expPrinting (Polynomial * p)
{
  int i;

  switch (p->eType) {
  case T_OFFLINE:
    importPoly(p);
    // Notice that there's no break here...

  case T_CONSTANT:
    fprintf (stderr, "%G", p->value);
    break;
  case T_VARIABLE:
    fprintf (stderr, "%s", p->e.v->vName);
    break;
  case T_EXTERNAL:
    fprintf (stderr, "%s()", p->e.e->polynomialFunctionName);
    break;
  case T_SUM:
    if (p->e.s->num > 1)
      fprintf (stderr, "(");
    importTermList (p);
    if (p->e.s->factor[0] != 1)
      fprintf (stderr, "%G*", p->e.s->factor[0]);
    expPrinting (p->e.s->sum[0]);
    for (i = 1; i < p->e.s->num; i++) {
      fprintf (stderr, "+");
      if (p->e.s->factor[i] != 1)
        fprintf (stderr, "%G*", p->e.s->factor[i]);
      expPrinting (p->e.s->sum[i]);
    }
    exportTermList (p, FALSE);
    if (p->e.s->num > 1)
      fprintf (stderr, ")");
    break;
  case T_PRODUCT:
    if (p->e.p->num > 1)
      fprintf (stderr, "(");
    expPrinting (p->e.p->product[0]);
    if (p->e.p->exponent[0] != 1)
      fprintf (stderr, "^%d", p->e.p->exponent[0]);
    for (i = 1; i < p->e.s->num; i++) {
      fprintf (stderr, "*");
      expPrinting (p->e.p->product[i]);
      if (p->e.p->exponent[i] != 1)
        fprintf (stderr, "^%d", p->e.p->exponent[i]);
    }
    if (p->e.p->num > 1)
      fprintf (stderr, ")");
    break;
  case T_FUNCTIONCALL:
    fprintf (stderr, "%s(", p->e.f->name);
    for (i = 0; i < p->e.f->num - 1; i++) {
      expPrinting (p->e.f->para[i]);
      fprintf (stderr, ", ");
    }
    expPrinting (p->e.f->para[p->e.f->num - 1]);
    fprintf (stderr, ")");
    break;
  case T_FREED:
    fprintf (stderr, "[FREED eType=%d id=%d index=%d key=%d count=%d valid=%d]\n", (int) p->value, p->id, p->index, p->key, p->count, p->valid);
  default:
    FATAL ("In expPrinting, unknown expression type %d", p->eType);
  }
}

/**
  Print the provided polynomial with sum and product terms in
  parentheses and preceeded by their type ('s' or 'p') and index
  down to depth indicated. At that level, simply refer to
  polynomial ids. This can produce useful output without excessive
  volume when used carefully.

*/
void expTermPrinting (FILE * output, Polynomial * p, int depth)
{
  int i;

  switch (p->eType) {
  case T_OFFLINE:
    importPoly(p);
    // Notice that there's no break here...

  case T_CONSTANT:
    fprintf (output, "%G", p->value);
    break;
  case T_VARIABLE:
    fprintf (output, "%s", p->e.v->vName);
    break;
  case T_EXTERNAL:
    fprintf (output, "%s()", p->e.e->polynomialFunctionName);
    break;
  case T_SUM:
    if (depth <= 0) {
      fprintf (output, "s%d", p->index);
      return;
    }
    fprintf (output, "(s%d:", p->index);
    if (p->e.s->factor[0] != 1)
      fprintf (output, "%G*", p->e.s->factor[0]);
    expTermPrinting (output, p->e.s->sum[0], depth - 1);
    importTermList (p);
    for (i = 1; i < p->e.s->num; i++) {
      fprintf (output, "+");
      if (p->e.s->factor[i] != 1)
        fprintf (output, "%G*", p->e.s->factor[i]);
      expTermPrinting (output, p->e.s->sum[i], depth - 1);
    }
    exportTermList (p, FALSE);
    fprintf (output, ")");
    break;
  case T_PRODUCT:
    if (depth <= 0) {
      fprintf (output, "p%d", p->index);
      return;
    }
    fprintf (output, "(p%d:", p->index);
    expTermPrinting (output, p->e.p->product[0], depth - 1);
    if (p->e.p->exponent[0] != 1)
      fprintf (output, "^%d", p->e.p->exponent[0]);
    for (i = 1; i < p->e.s->num; i++) {
      fprintf (output, "*");
      expTermPrinting (output, p->e.p->product[i], depth - 1);
      if (p->e.p->exponent[i] != 1)
        fprintf (output, "^%d", p->e.p->exponent[i]);
    }
    fprintf (output, ")");
    break;
  case T_FUNCTIONCALL:
    if (depth <= 0) {
      fprintf (output, "f%d", p->index);
      return;
    }
    fprintf (output, "%s(", p->e.f->name);
    for (i = 0; i < p->e.f->num - 1; i++) {
      expTermPrinting (output, p->e.f->para[i], depth - 1);
      fprintf (output, ", ");
    }
    expTermPrinting (output, p->e.f->para[p->e.f->num - 1], depth - 1);
    fprintf (output, ")");
    break;
  case T_FREED:
    fprintf (stderr, "[FREED eType=%d id=%d index=%d key=%d count=%d valid=%d]\n", (int) p->value, p->id, p->index, p->key, p->count, p->valid);
  default:
    FATAL ("In expTermPrinting, unknown expression type %d", p->eType);
  }
}


/** Explicitly discard a poly. Clear the VALID_NOTDISC_FLAG of valid to 
    indicate that the user explicitly discarded the polynomial. When we
    accumulate enough to be worth the effort, go purge them out. */
inline void discardPoly (Polynomial *p)
{
  unsigned short allFlags = VALID_EVAL_FLAG|VALID_KEEP_FLAG|VALID_REF_FLAG|VALID_TOP_FLAG|VALID_NOTDISC_FLAG;

  if (p->eType == T_CONSTANT && (p->value == 1 || p->value == 0 || p->value == 0.5)) {
    //    fprintf (stderr, "Not discarding a 1, 0 or 0.5 constant\n");
    return;
  }

  if ((p->valid & ~VALID_NOTDISC_FLAG) != 0) {
    FATAL ("Can't discard already flagged poly!");
  }

  p->valid &= ~VALID_NOTDISC_FLAG; // Explicitly discarded
  if (p->id == polynomialLostNodeId) {
    fprintf (stderr, "discardPoly sees id %d with valid %d and count %d\n", p->id, p->valid, p->count);
    raise (SIGINT);
  }
  //  fprintf (stderr, "Flagged polynomial %u as explicitly discarded, valid is %u, value is %g\n", p->id, p->valid, p->value);

  /*
#ifdef SOURCEDIGRAPH
  if (++pendingExplicitDiscards >= 1000) {
    dumpSourceParenting ();
    exit(EXIT_FAILURE);
  }
#endif
  */

  if (++pendingExplicitDiscards >= 100000) {

    swLogProgress(3 /* DETAIL + 1 */, 0, "Freeing %d explicitly discarded polynomials", pendingExplicitDiscards);
    /*
    int i, j;
    if (constantCount > 0) {
      fprintf (stderr, "All %d constants:\n", constantCount);
      for (i = 0; i < CONSTANT_HASH_SIZE; i++) {
	if (constantHash[i].num <= 0)
	  continue;
	for (j = 0; j < constantHash[i].num; j++) {
	  fprintf (stderr, "(%d %d) id=%d index=%d key=%d count=%d valid=%d constant: ", i, j, constantList[constantHash[i].index[j]]->id, constantHash[i].index[j], constantHash[i].key[j], constantList[constantHash[i].index[j]]->count, constantList[constantHash[i].index[j]]->valid);
	  expTermPrinting (stderr, constantList[constantHash[i].index[j]], 1);
	  fprintf (stderr, "\n");
	}
      }
      fprintf (stderr, "\n");
    }
    */
    if (polynomialDebugLevel >= 5)
      fprintf (stderr, "Free w/mask of all flags (%d) preserved...", allFlags);
    doFreePolys (allFlags);
    if (polynomialDebugLevel >= 5)
      fprintf (stderr, "\n");
    pendingExplicitDiscards = 0;
    //    exit (EXIT_SUCCESS);
  }
  return;
}

/**

  This function generates a sum polynomial.  It accepts a group of <factor, poly> pairs.
  The sum is represented as factor_1*poly_1+factor_2*poly_2+...+factor_n*poly_n.
  The sub polynomials are checked to see if it can be combined with other polynomials for
  simplification.  Also, if the sum is a constant, the result will be a constant.
  If the polynomial being created exists, the existing polynomial is returned.  Otherwise,
  a new polynomial is constructed.

  The parameters are:
  -# the number of <factor, poly> pairs in the parameter list
  -# a group of <factor, poly> pairs
  -# a flag whose value is set to be 0 for plus operations like
     p=p1 + p2 + p3 ... where p1, p2, p3 may not be freed and
     1 for plus operations like p = p  + p2 + p3 ... where the
     the parameter polynomial p will be replaced by the new sum
     polynomial and therefore the parameter polynomial p can be
     freed.

*/
Polynomial *plusExp (char *fileName, int lineNo, int num, ...)
{
  int i, k, l;
  va_list args; ///< Variable list of parameters
  struct sumPoly *sP;   ///< Element of a sum polynomial
  Polynomial *rp;       ///< Pointer of the new sum polynomial
  int counterSum;       ///< Counter of the number of terms in the new sum polynomial
  int flag;     //< 1st-term replacement flag
  double f1, f0;        ///< Some factors
  Polynomial *p1 = 0, *p0 = 0;  ///< Some polynomial terms
  double con = 0;       ///< Accumulates constant values appearing in the function parameters
  int key = 0;  ///< Key for the new polynomial
  int tempI, tempI2;    ///< Normalized powers (for comparison of two double numbers)
  double tempD, tempD2; ///< Normalized fractions( for comparison of two double numbers)
  int hIndex, sIndex;   ///< Index of the new polynomial in the hash table, sub index in a hash bucket
  int first,    ///< First polynomial for comparison
    last,       ///< Last polynomial for comparison
    location;   ///< Location of the new polynomial in a hash bucket 
  int p0SubHIndex = 0,  ///< Sub-index in a hash bucket
      p0HIndex = 0,     ///< Index in the hash table
      p0Index = 0,      ///< Index in the polynomial list
      p0Id = 0, p0Key, p0Valid;
  enum expressionType p0EType;

  termCount++;
  if ((termCount & 0x1FFFFF) == 0) {
    if ((termCount & 0xFFFFFF) == 0) {
#ifdef POLYSTATISTICS
      polyStatistics ("At 8M raw term multiple");
#endif
      swLogProgress(3 /* DETAIL + 1 */, 0, "Building polynomial, currently at %1.2g terms, %2.3gGb",
		    (double) sumCount + productCount + constantCount + variableCount + functionCallCount,
		    polyMemory() / (1024 * 1024 * 1024));
    } else
      swLogProgress(4 /* DETAIL + 2 */, 0, "Building polynomial, currently at %1.2g terms",
		    (double) sumCount + productCount + constantCount + variableCount + functionCallCount);
  }

  // Initialize the variables that keep tracks the number of terms collected for the new sum polynomial
  counter_v1 = 0;
  counter_p1 = 0;
  counter_f1 = 0;

  // Get the number of items for this sum from the parameter list of the function
  va_start (args, num);

#ifdef SOURCEDIGRAPH
  memset (originalChildren, 0, sizeof (originalChildren));
#endif

  // Iterate through all the items in the parameter list
  for (i = 0; i < num; i++) {
    f1 = va_arg (args, double); // Get the coefficient
    p1 = va_arg (args, Polynomial *);   // ...and the polynomial
    if (p1->eType == T_OFFLINE) importPoly (p1);

    //    if (p1->count == 0 && !((p1->valid & VALID_TOP_FLAG) == 0))
    //      fprintf (stderr, "UnHeld, non-top SUM %d referenced!\n", p1->id);

    p1->valid &= ~VALID_TOP_FLAG;       // No longer unreferenced

#ifdef SOURCEDIGRAPH
    originalChildren[p1->source]++;     // Bump the count of children of this subpoly source for this parent
#endif

    if (polynomialDebugLevel >= 60) {
      fprintf (stderr, "In plusExp factor=%f item No. %d of %d type=%d\n", f1, i + 1, num, p1->eType);
      if (polynomialDebugLevel >= 70) {
        expTermPrinting (stderr, p1, 1);
        fprintf (stderr, "\n");
      }
    }
    /* Record the first operand of the plus operation.  Often, the first operand and the result are the
     * same variable, e.g. p=p+..., therefore, when we build a new polynomial for the result of the
     * plus operation, the first operand polynomial of the plus operation can be freed. */
    if (i == 0) {
      f0 = f1;
      p0 = p1;
    }
    // If a term is constant 0, it has no effect on a sum, therefore we do nothing
    if (f1 == 0.0)
      continue;
    switch (p1->eType) {
    case T_OFFLINE:
      importPoly(p1);
      // Notice that there's no break here...

    case T_CONSTANT:   // If the term is a non-zero constant, we collect it
      con += p1->value * f1;
      break;
    case T_VARIABLE:   // The term is a variable
      collectSumTerms (&factor_v1, &p_v1, &counter_v1, &containerLength_v1, f1, p1);
      break;
    case T_PRODUCT:    // The item is a product
      collectSumTerms (&factor_p1, &p_p1, &counter_p1, &containerLength_p1, f1, p1);
      break;
    case T_FUNCTIONCALL:       // The term is a function call
      collectSumTerms (&factor_f1, &p_f1, &counter_f1, &containerLength_f1, f1, p1);
      break;
    case T_SUM:        // The term is a sum
      // Fold the components of the sum subpolynomial up into this sum polynomial
      importTermList (p1);
      for (l = 0; l < p1->e.s->num; l++) {
        switch (p1->e.s->sum[l]->eType) {
	case T_OFFLINE:
	  importPoly(p1->e.s->sum[l]);
	// Notice that there's no break here...

        case T_CONSTANT:       // Item is a constant, collected it into con
          con += f1 * p1->e.s->sum[l]->value * p1->e.s->factor[l];
          break;
        case T_VARIABLE:
          collectSumTerms (&factor_v1, &p_v1, &counter_v1, &containerLength_v1, f1 * p1->e.s->factor[l], p1->e.s->sum[l]);
          break;
        case T_PRODUCT:
          collectSumTerms (&factor_p1, &p_p1, &counter_p1, &containerLength_p1, f1 * p1->e.s->factor[l], p1->e.s->sum[l]);
          break;
        case T_FUNCTIONCALL:
          collectSumTerms (&factor_f1, &p_f1, &counter_f1, &containerLength_f1, f1 * p1->e.s->factor[l], p1->e.s->sum[l]);
          break;
        default:
          FATAL ("In plusExp, unknown expression type %d", p1->e.s->sum[l]->eType);
        }
      }
      exportTermList (p1, FALSE);
      break;
    case T_FREED:      // The term being referenced has been freed!
      expTermPrinting (stderr, p1, 1);
      fprintf (stderr, "\n");
      FATAL ("In plusExp, preceeding polynomial term %d was freed while still referenced", i);
      break;
    default:
      FATAL ("In plusExp, unknown polynomial type %d", p1->eType);
    }
  }
  flag = va_arg (args, int);

  va_end (args);

  if (flag == 0)
    sumNotReleaseableCount++;
  else
    sumReleaseableCount++;

  // Merge the collected polynomial terms to form a sum polynomial
  counterSum = counter_v1 + counter_p1 + counter_f1;
  if (counterSum + 1 > lengthSum) {
    lengthSum = counterSum + 1;
    REALCHOKE (factorSum, lengthSum * sizeof (double), double *);
    REALCHOKE (pSum, lengthSum * sizeof (Polynomial *), Polynomial **);
  }
  /* While we do end-up with zero factors for variables here, they get swallowed-up
   * into the next tier where they become non-zero, so it all comes out in the wash. */
  if (counter_v1 > 0) {
    memcpy (&factorSum[0], &factor_v1[0], sizeof (double) * counter_v1);
    memcpy (&pSum[0], &p_v1[0], sizeof (Polynomial *) * counter_v1);
  }
  if (counter_p1 > 0) {
    memcpy (&factorSum[counter_v1], &factor_p1[0], sizeof (double) * counter_p1);
    memcpy (&pSum[counter_v1], &p_p1[0], sizeof (Polynomial *) * counter_p1);
  }
  if (counter_f1 > 0) {
    memcpy (&factorSum[counter_v1 + counter_p1], &factor_f1[0], sizeof (double) * counter_f1);
    memcpy (&pSum[counter_v1 + counter_p1], &p_f1[0], sizeof (Polynomial *) * counter_f1);
  }
  // Handle the possible simple outcomes of collecting terms...
  if (counterSum == 0) {
    // After we go through all the items in the sum, we get only a constant.
    if (p0->eType == T_CONSTANT && p0->value == con) {
      sumReturnSameConstantCount++;
      return polyReturnWrapper (p0); // No net change at all from first term. Just return it.
    }
    sumReturnConstantCount++;
    if (flag != 0 && (p0->valid & ~VALID_NOTDISC_FLAG) == 0)
      discardPoly (p0);
    rp = constantExp (__FILE__, __LINE__, con);
    if (polynomialDebugLevel >= 60)
      fprintf (stderr, "Returning a constant %f\n", rp->value);
    return polyReturnWrapper (rp);
  } else if (con == 0.0 && counterSum == 1 && factorSum[0] == 1.0) {
    // We have only one term, no need to create a new polynomial
    rp = pSum[0];
    sumReturn1TermCount++;
    if (polynomialDebugLevel >= 60) {
      fprintf (stderr, "Returning a single term sum\n");
      expTermPrinting (stderr, rp, 1);
      fprintf (stderr, "\n");
    }
    return polyReturnWrapper (rp);
  } else {
    // We've got more than one term, so we're going to go ahead and create it.
    if (con != 0.0) {
      p1 = constantExp (__FILE__, __LINE__, con);
      factorSum[counterSum] = 1.0;
      pSum[counterSum] = p1;
      counterSum++;
    }
  }

  // Compute the key for this polynomial
  key = keySumPolynomial (pSum, factorSum, counterSum);

  hIndex = key % SUM_HASH_SIZE;
  if (hIndex < 0)
    hIndex += SUM_HASH_SIZE;

  /* If the hash table item is not empty, look for the sum in the sum polynomial 
   * list. If it is not there, determine a position in the hash table to
   * save the key and the index, and add this polynomial to the list */
  if (sumHash[hIndex].num > 0) {
    // If the key is already there, compare to polynomials under that key
    if (searchHashTable (&sumHash[hIndex], &first, &last, key)) {
      for (i = first; i <= last; i++) {
        sIndex = sumHash[hIndex].index[i];
        if (counterSum == sumList[sIndex]->e.s->num) {
          // Compare the two sums term-by-term
	  importTermList (sumList[sIndex]);
          for (k = 0; k < counterSum; k++) {
            // If the terms are identical, compare their factors
            if (pSum[k] == sumList[sIndex]->e.s->sum[k]) {
              tempD = frexp (factorSum[k], &tempI);
              tempD2 = frexp (sumList[sIndex]->e.s->factor[k], &tempI2);
              if (tempI2 != tempI || (int) (tempD2 * 100000000) != (int) (tempD * 100000000))
                break;
            } else
              break;
          }
	  exportTermList (sumList[sIndex], FALSE);
          if (k >= counterSum) {
            sumList[sIndex]->valid |= VALID_REF_FLAG;
            sumHashHits++;
            if (polynomialDebugLevel >= 60) {
              fprintf (stderr, "Returning an existing sum...\n");
              if (polynomialDebugLevel >= 70) {
                expPrinting (sumList[sIndex]);
                fprintf (stderr, "\n");
              }
            }
	    //	    if (flag != 0 && (p0->valid & ~VALID_NOTDISC_FLAG) == 0)
	    //	      discardPoly (p0); // Causes a crash, but shouldn't! Maybe sumList == p0?
            return polyReturnWrapper (sumList[sIndex]);
          }
        }
      }
      // Identical key, but not identical polynomial
      location = last;
    } else
      location = last;
  } else
    location = 0;

  // If the first polynomial in the parameter list can be freed, do so.
  if (flag != 0 && p0->eType == T_SUM && p0->valid == 0) {
    p0Index = p0->index;
    p0Id = p0->id;
    p0Key = p0->key;
    p0HIndex = p0->key % SUM_HASH_SIZE;
    if (p0HIndex < 0)
      p0HIndex += SUM_HASH_SIZE;
    p0EType = p0->eType;
    p0Valid = p0->valid;

    if (sumHash[p0HIndex].num <= 0) {
      expTermPrinting (stderr, p0, 1);
      fprintf (stderr, "\n");
      FATAL ("Preceeding polynomial %d is not in the sum hash table (empty list)", p0Id);
    }
    // Locate the polynomial in the hash table
    if (searchHashTable (&sumHash[p0HIndex], &first, &last, p0Key)) {
      for (i = first; i <= last; i++) {
        if (p0Index == sumHash[p0HIndex].index[i])
          break;
      }
      if (i <= last)
        p0SubHIndex = i;
      else {
        expTermPrinting (stderr, p0, 1);
        fprintf (stderr, "\n");
        FATAL ("Polynomial %d is not in the sum hash list (not in list)", p0Id);
      }
    } else {
      expTermPrinting (stderr, p0, 1);
      fprintf (stderr, "\n");
      FATAL ("Polynomial %d is not in the sum hash list (bad searchHashTable return)", p0Id);
    }

    // Free the memory of the first polynomial in the parameter list
    sum1stTermsFreedCount++;
#ifndef USE_SSD
    free (p0->e.s->sum);
    free (p0->e.s->factor);
#else
    deportTermList (p0);
#endif
    free (p0->e.s);
#ifdef FREEDEBUG
    p0->value = -sumList[i]->eType;
    p0->eType = T_FREED;
#else
    free (p0);
#endif
  } else {
    // Flag not set or not sum or actually valid
    p0EType = T_FREED;
    p0Valid = 7;
  }

  // Since the sum was not found in the sum list, a new polynomial is built
  MALCHOKE(rp, sizeof (Polynomial), Polynomial *);
  rp->eType = T_SUM;
  MALCHOKE(sP, sizeof (struct sumPoly), struct sumPoly *);
  sP->num = counterSum;
  // When usage is up to around 32Gb, these two mallocs can take seconds and even minutes!
  MALCHOKE(sP->factor, counterSum * sizeof (double),double *);
  MALCHOKE(sP->sum, counterSum * sizeof (Polynomial *),Polynomial **);
  for (i = 0; i < sP->num; i++) {
    sP->factor[i] = factorSum[i];
    sP->sum[i] = pSum[i];
    sP->sum[i]->valid |= VALID_REF_FLAG;
  }
  // Assign values to some attributes of the sum polynomial
#ifdef USE_SSD
  sP->iMTLIndex = -1;
#endif
  rp->e.s = sP;
  rp->key = key;
  rp->valid = VALID_NOTDISC_FLAG;
  rp->count = 0;
#ifdef USE_GMP
  mpf_init (rp->mpfValue);
#endif
#ifdef SOURCEDIGRAPH
  rp->source = findOrAddSource (fileName, lineNo, T_SUM);
#endif

  // Insert the new built polynomial in sum list

  /* If the first polynomial is the parameter list is freed, its position
   * in the polynomial list is occupied by the newly created sum polynomial. */
  if (flag != 0 && p0EType == T_SUM && p0Valid == 0) {
    sumListReplacementCount++;
    sumList[p0Index] = rp;
    sumList[p0Index]->index = p0Index;
    sumList[p0Index]->id = p0Id;
  } else {
    if (sumCount >= sumListLength) {
      sumListLength += SUM_LIST_INCREASE;
      sumPListExpansions++;
      REALCHOKE (sumList, sumListLength * sizeof (Polynomial *), void *);
    }
    sumListNewCount++;
    sumList[sumCount] = rp;
    sumList[sumCount]->index = sumCount;
    sumList[sumCount]->id = nodeId;
    if (nodeId == polynomialLostNodeId) {
      fprintf (stderr, "nodeId %d has valid of %d and count of %d\n", polynomialLostNodeId, rp->valid, rp->count);
      raise (SIGINT);
    }
    sumList[sumCount]->valid |= VALID_TOP_FLAG; // Currently unreferenced

    sumCount++;
    nodeId++;
  }

  // Insert the newly built polynomial into the Hash table
  if (flag != 0 && p0EType == T_SUM && p0Valid == 0) {
    if (p0SubHIndex != location || p0HIndex != hIndex) {
      insertHashTable (&sumHash[hIndex], location, key, sumList[p0Index]->index);
      // Delete the replaced polynomial from the hash table
      if (p0HIndex != hIndex)
        deleteHashTable (&sumHash[p0HIndex], p0SubHIndex);
      else if (p0SubHIndex < location)
        deleteHashTable (&sumHash[p0HIndex], p0SubHIndex);
      else
        deleteHashTable (&sumHash[p0HIndex], p0SubHIndex + 1);
    } else
      /* If the newly built polynomial and the one to be replaced are at the same position
       * in the hash table, just replace the key of the old polynomial with the key of the
       * new one. */
      sumHash[hIndex].key[location] = key;
  } else
    insertHashTable (&sumHash[hIndex], location, key, sumCount - 1);

  sumNewCount++;
  if (polynomialDebugLevel >= 60) {
    fprintf (stderr, "Returning a new sum...\n");
    if (polynomialDebugLevel >= 70) {
      expPrinting (rp);
      fprintf (stderr, "\n");
    }
  }
#ifdef USE_SSD
  if (sP->num >= MIN_USE_SSD_DPS)
    importedTerms++; // Because the following export is unpaired...
  if ((importedTerms-(exportedWrittenTerms+exportedDroppedTerms)) > peakInMemoryTerms)
    peakInMemoryTerms = importedTerms-(exportedWrittenTerms+exportedDroppedTerms);
#endif
  exportTermList (rp, TRUE);

  return polyReturnWrapper (rp);
};

/**

  This function computes a key for a product polynomial.  The number of
  product polynomials is usually very big.  Good keys for product polynomials
  are required for good distribution of product polynomials in product hash table.

*/
inline int keyProductPolynomial (Polynomial ** p, int *exponent, int counter)
{
  int key = 0;
  int j;

  for (j = 0; j < counter; j++) {
    if (p[j]->eType == T_OFFLINE) importPoly (p[j]);
    key += p[j]->id + p[j]->key * (j + 1) + (exponent[j] * 12345 * (int) (p[j]->eType + 1));
  }
  return key;

}

/**

  This function collects the terms for constructing a product polynomial.  The terms are sorted
  so that these terms can have a unique order in the constructed product polynomial.

*/
inline void collectProductTerms (int **exponent,        ///< Container, pointer to a list of int exponents
    Polynomial *** p,   ///< Container, pointer to a list of pointers to polynomials
    int *counter,       ///< Number of terms currently in the container lists
    int *containerLength,       ///< Current size of the container lists
    int e1,     ///< Exponent of term to add to container
    Polynomial * p1     ///< Pointer to polynomial of term to add to container
    )
{
  int location;

  // Search for a position for this component in the container
  if (searchPolynomialList (*p, *counter, p1, &location) == 1) {
    (*exponent)[location] += e1;
    if ((*exponent)[location] == 0)
      WARNING ("Zero polynomial exponent terms could be replaced with constant 1!");
    return;
  }

  if ((*counter) >= (*containerLength) - 1) {
    (*containerLength) += 50;
    REALCHOKE ((*exponent), (*containerLength) * sizeof (int), int *);
    REALCHOKE ((*p), (*containerLength) * sizeof (Polynomial *), Polynomial **);
  }
  // This is a new item in the product, insert it at the start or end
  if (location >= (*counter)) {
    (*p)[(*counter)] = p1;
    (*exponent)[(*counter)] = e1;
    (*counter)++;
  }
  // Make space so we can insert the item where it goes in the middle of the product
  else {
    memmove (&((*p)[location + 1]), &((*p)[location]), sizeof (Polynomial *) * ((*counter) - location));
    memmove (&((*exponent)[location + 1]), &((*exponent)[location]), sizeof (int) * ((*counter) - location));

    // Insert the new item.
    (*p)[location] = p1;
    (*exponent)[location] = e1;
    (*counter)++;
  }
};


/**

  This function creates a product polynomial.  It accepts a group of <exponent, poly>
  pairs as the operands of a times operation of poly_1^exponet1*poly_2^exponent_2*...
  *poly_n^exponent_n. It checks the operands to see whether they can be combined for 
  simplification. It is also possible the result is a polynomial other than type of a
  product if the the result isn't expressed by a product polynomial.

  The parameters are:
  -# the number of <poly, exponent> pairs in the parameter list
  -# a group of <poly, exponent> pairs

  -# a flag whose value is set to be 0 for times operations like
     p=p1 * p2 * p3 ... where p1, p2, p3 may not be freed
     and 1 for times operations like p = p * p2 * p3... where 
     the the parameter polynomial p will be replaced by the new
     product polynomial and therefore the parameter polynomial p 
     can be freed.

*/
Polynomial *timesExp (char *fileName, int lineNo, int num, ...)
{
  int i, k, l;
  int counterProd;      ///< Number of terms in the new product polynomial
  va_list args; ///< Variable list of parameters
  struct productPoly *pP;       ///< Product structure for the new polynomial
  Polynomial *rp;       ///< New product polynomial
  Polynomial *p1 = 0, *p0 = 0;  ///< Some polynomial terms
  int e1, e0;   ///< Component of a term
  int isZero = 0;       ///< If the result is zero
  double factor = 1;    ///< Collect the constant terms in the parameter list
  int key = 0;  ///< Key of the new polynomial
  int pIndex, hIndex;   ///< Index in the polynomial list, index in the hash table
  int first,    ///< First term for comparison
    last,       ///< Last term for comparison
    location;   ///< Position of the new polynomial in a hash bucket
  int flag;     ///< If the new product polynomial can replace an existing one
  int p0SubHIndex = 0,  ///< Sub index in a hash bucket
      p0HIndex = 0,     ///< Index in the hash table
      p0Index = 0,      ///< Index in the polynomial list
      p0Id = 0, ///< Unique Id
      p0Key,    ///< Key
      p0Valid;  ///< Valid byte
  enum expressionType p0EType;

  termCount++;
  if ((termCount & 0x1FFFFF) == 0) {
    if ((termCount & 0xFFFFFF) == 0) {
#ifdef POLYSTATISTICS
      polyStatistics ("At 8M raw term multiple");
#endif
      swLogProgress(3 /* DETAIL + 1 */, 0, "Building polynomial, currently at %1.2g terms, %2.3gGb",
		    (double) sumCount + productCount + constantCount + variableCount + functionCallCount,
		    polyMemory() / (1024 * 1024 * 1024));
    } else
      swLogProgress(4 /* DETAIL + 2 */, 0, "Building polynomial, currently at %1.2g terms",
		    (double) sumCount + productCount + constantCount + variableCount + functionCallCount);
  }

  // Initialize the containers for the operands of a times operation
  counter_v2 = 0;
  counter_s2 = 0;
  counter_f2 = 0;

#ifdef SOURCEDIGRAPH
  memset (originalChildren, 0, sizeof (originalChildren));
#endif

  // Get the number of operands for this times operation
  va_start (args, num);

  // Go through operand and its exponent of the product
  for (i = 0; i < num; i++) {
    p1 = va_arg (args, Polynomial *);
    if (p1->eType == T_OFFLINE) importPoly (p1);
    e1 = va_arg (args, int);

    //    if (p1->count == 0 && !((p1->valid & VALID_TOP_FLAG) == 0))
    //      fprintf (stderr, "UnHeld, non-top PRODUCT %d referenced!\n", p1->id);

    p1->valid &= ~VALID_TOP_FLAG;       // No longer unreferenced

#ifdef SOURCEDIGRAPH
    originalChildren[p1->source]++;     // Bump the count of children of this subpoly source for this parent
#endif

    if (polynomialDebugLevel >= 60) {
      fprintf (stderr, "In timesExp exponent=%d item No. %d of %d type=%d\n", e1, i + 1, num, p1->eType);
      if (polynomialDebugLevel >= 70) {
        expTermPrinting (stderr, p1, 1);
        fprintf (stderr, "\n");
      }
    }
    /* Record the first operand of the times operation.  Often, the first operand and the result are the
     * same variable, which refers to the case of p=p*..., therefore, when we build a new polynomial
     * for the result of the times operation, the polynomial that represents the first operand of the
     * times operation can be freed. */
    if (i == 0) {
      e0 = e1;
      p0 = p1;
    }
    // If any of the operand is zero, then the product is zero
    if (isZeroExp (p1)) {
      isZero = 1;
      break;
    }
    // If this operand is a constant, this constant and its exponent are accumulated into factor
    else if (p1->eType == T_CONSTANT) {
      factor *= pow (p1->value, e1);
      continue;
    } else {
      if (p1->eType == T_SUM && p1->e.s->num == 1) {    // A sum with only one item?
	Polynomial *q = p1;
	importTermList (q);
	factor *= pow (p1->e.s->factor[0], e1);
	p1 = p1->e.s->sum[0];
	exportTermList (q, FALSE);
      }
      /* If this operand is a variable, a sum that has more than one items, or a function call, 
       * we directly collect it.  If this operand is a product, we then collect each term of it. */
      switch (p1->eType) {
      case T_OFFLINE:
	importPoly(p1);
	// Notice that there's no break here...

      case T_VARIABLE:
        collectProductTerms (&exponent_v2, &p_v2, &counter_v2, &containerLength_v2, e1, p1);
        break;
      case T_SUM:
	importTermList (p1);
	collectProductTerms (&exponent_s2, &p_s2, &counter_s2, &containerLength_s2, e1, p1);
	exportTermList (p1, FALSE);
        break;
      case T_FUNCTIONCALL:
        collectProductTerms (&exponent_f2, &p_f2, &counter_f2, &containerLength_f2, e1, p1);
        break;
        // For a product term, we open it and check each of its terms
      case T_PRODUCT:
        for (l = 0; l < p1->e.p->num; l++) {
          switch (p1->e.p->product[l]->eType) {
	  case T_OFFLINE:
	    importPoly(p1->e.p->product[l]);
	    // Notice that there's no break here...

          case T_VARIABLE:
            collectProductTerms (&exponent_v2, &p_v2, &counter_v2, &containerLength_v2, e1 * p1->e.p->exponent[l], p1->e.p->product[l]);
            break;
          case T_SUM:
	    importTermList (p1->e.p->product[l]);
            collectProductTerms (&exponent_s2, &p_s2, &counter_s2, &containerLength_s2, e1 * p1->e.p->exponent[l], p1->e.p->product[l]);
	    exportTermList (p1->e.p->product[l], FALSE);
            break;
          case T_FUNCTIONCALL:
            collectProductTerms (&exponent_f2, &p_f2, &counter_f2, &containerLength_f2, e1 * p1->e.p->exponent[l], p1->e.p->product[l]);
            break;
          default:
            FATAL ("In timesExp, unknown expression type %d", p1->e.p->product[l]->eType);
          }
        }
        break;
      case T_FREED:    // The term being referenced has been freed!
        expTermPrinting (stderr, p1, 1);
	fprintf (stderr, "\n");
        FATAL ("In plusExp, preceeding polynomial term %d was freed while still referenced", i);
        break;
      default:
        FATAL ("In timesExp, unknown polynomial id %u type %d", p1->id, p1->eType);
        break;
      }
    }
  }
  flag = va_arg (args, int);

  va_end (args);

  if (flag == 0)
    productNotReleaseableCount++;
  else
    productReleaseableCount++;

  if (isZero) {
    // The product is zero, a zero polynomial is returned
    productReturn0Count++;
    if (flag != 0 && (p0->valid & ~VALID_NOTDISC_FLAG) == 0)
      discardPoly (p0);
    rp = constantExp (__FILE__, __LINE__, 0.0);
    if (polynomialDebugLevel >= 60)
      fprintf (stderr, "Returning a constant zero\n");
    return polyReturnWrapper (rp);
  }
  /* Combine all the items generated by collecting the operands of the times operation 
   * into one polynomial */
  counterProd = counter_v2 + counter_s2 + counter_f2;
  if (counterProd > lengthProd) {
    lengthProd = counterProd;
    REALCHOKE (exponentProd, lengthProd * sizeof (int), int *);
    REALCHOKE (pProd, lengthProd * sizeof (Polynomial *), Polynomial **);
  }
  if (counter_v2 > 0) {
    memcpy (&exponentProd[0], &exponent_v2[0], sizeof (int) * counter_v2);
    memcpy (&pProd[0], &p_v2[0], sizeof (Polynomial *) * counter_v2);
  }
  if (counter_s2 > 0) {
    memcpy (&exponentProd[counter_v2], &exponent_s2[0], sizeof (int) * counter_s2);
    memcpy (&pProd[counter_v2], &p_s2[0], sizeof (Polynomial *) * counter_s2);
  }
  if (counter_f2 > 0) {
    memcpy (&exponentProd[counter_v2 + counter_s2], &exponent_f2[0], sizeof (int) * counter_f2);
    memcpy (&pProd[counter_v2 + counter_s2], &p_f2[0], sizeof (Polynomial *) * counter_f2);
  }
  if (counterProd == 0) {
    // The product has 0 items, the result is a constant polynomial
    if (p0->eType == T_CONSTANT && p0->value == factor) {
      productReturnSameConstantCount++;
      return polyReturnWrapper (p0); // No net change at all from first term. Just return it.
    }
    productReturnConstantCount++;
    if (flag != 0 && (p0->valid & ~VALID_NOTDISC_FLAG) == 0)
      discardPoly (p0);
    rp = constantExp (__FILE__, __LINE__, factor);
    if (polynomialDebugLevel >= 60)
      fprintf (stderr, "Returning a constant %f\n", rp->value);
    return polyReturnWrapper (rp);

  } else        // If the result polynomial has only one term, it is not a product polynomial
  if (counterProd == 1 && exponentProd[0] == 1) {
    if (factor == 1.0) {
      // If the factor is 1, then the result polynomial is equal to its first term
      rp = pProd[0];
      if (polynomialDebugLevel >= 60)
        fprintf (stderr, "Returning first term from caller\n");
      productReturn1stTermCount++;
      return polyReturnWrapper (rp);
    } else {
      // If the factor is not 1, then the result polynomial is a sum polynomial
      rp = plusExp (fileName, lineNo, 1, factor, pProd[0], flag);
      if (polynomialDebugLevel >= 60)
        fprintf (stderr, "Returning via plusExp\n");
      productReturn1TermSumCount++;
      return polyReturnWrapper (rp);
    }
  } else {
    // The result polynomial is a product polynomial, compute the key
    key = keyProductPolynomial (pProd, exponentProd, counterProd);
    // Compute hash table index according to the key
    hIndex = key % PRODUCT_HASH_SIZE;
    if (hIndex < 0)
      hIndex += PRODUCT_HASH_SIZE;
    /* If the hash table item is not empty, look for the product in the product polynomial
     * list.  If it is not there, determine a position in the hash table to save the 
     * key and the index in the hash table. */
    if (productHash[hIndex].num > 0) {
      // If the key is already there, compare to polynomials under that key
      if (searchHashTable (&productHash[hIndex], &first, &last, key)) {
        for (i = first; i <= last; i++) {
          pIndex = productHash[hIndex].index[i];
          // Compare the two products term-by-term
          if (counterProd == productList[pIndex]->e.p->num) {
            for (k = 0; k < counterProd; k++)
              if (pProd[k] != productList[pIndex]->e.p->product[k]
                  || exponentProd[k] != productList[pIndex]->e.p->exponent[k])
                break;
            if (k >= counterProd) {
              // If the polynomials are the same compare the exponents
              if (factor == 1.0) {
                productList[pIndex]->valid |= VALID_REF_FLAG;
                productHashHits++;
                if (polynomialDebugLevel >= 60) {
                  fprintf (stderr, "Returning an existing product...\n");
                  if (polynomialDebugLevel >= 70) {
                    expPrinting (productList[pIndex]);
                    fprintf (stderr, "\n");
                  }
                }
		if (flag != 0 && (p0->valid & ~VALID_NOTDISC_FLAG) == 0)
		  discardPoly (p0);
                return polyReturnWrapper (productList[pIndex]);
              } else {
                /* This product polynomial already exists so we don't need to construct a
                 * new one, however, the factor is not 1, therefore the result polynomial is
                 * a sum polynomial */
                productHashHitIsSumCount++;
                productList[pIndex]->valid |= VALID_REF_FLAG;
                if (polynomialDebugLevel >= 60)
                  fprintf (stderr, "Returning existing product via plusExp\n");
                return plusExp (fileName, lineNo, 1, factor, polyReturnWrapper(productList[pIndex]), flag);
              }
            }
          }
        }
        // Identical key, but not identical polynomial
        location = last;
      }
      location = last;
    } else
      location = 0;

    // If the first operand can be freed, do so.
    if (flag != 0 && p0->eType == T_PRODUCT && p0->valid == 0) {
      /* We save the identity information of the polynomial to be freed
       * so that it can be used by the newly created polynomial. */
      p0Index = p0->index;
      p0Id = p0->id;
      p0Key = p0->key;
      p0HIndex = p0Key % PRODUCT_HASH_SIZE;
      if (p0HIndex < 0)
        p0HIndex += PRODUCT_HASH_SIZE;
      p0EType = p0->eType;
      p0Valid = p0->valid;

      // Determine the sub index of the old polynomial in the product hash table
      if (productHash[p0HIndex].num <= 0) {
        expTermPrinting (stderr, p0, 1);
        fprintf (stderr, "\n");
        FATAL ("Preceeding polynomial is not in the product hash list (1)");
      }
      if (searchHashTable (&productHash[p0HIndex], &first, &last, p0Key)) {
        for (i = first; i <= last; i++) {
          if (p0Index == productHash[p0HIndex].index[i])
            break;
        }
        if (i <= last)
          p0SubHIndex = i;
        else {
          expTermPrinting (stderr, p0, 1);
          fprintf (stderr, "\n");
          FATAL ("Preceeding polynomial is not in the product hash list (2)");
        }
      } else {
        expTermPrinting (stderr, p0, 1);
        fprintf (stderr, "\n");
        FATAL ("Preceeding polynomial is not in the product hash list (3)");
      }

      // Free the first operand
      product1stTermsFreedCount++;
      free (p0->e.p->exponent);
      free (p0->e.p->product);
      free (p0->e.p);
#ifdef FREEDEBUG
      p0->value = -productList[i]->eType;
      p0->eType = T_FREED;
#else
      free (p0);
#endif
    } else {
      p0EType = T_FREED;
      p0Valid = 7;
    }

    // Construct a new product polynomial from the terms saved in the container
    MALCHOKE(rp, sizeof (Polynomial), Polynomial *);
    rp->eType = T_PRODUCT;
    MALCHOKE(pP, sizeof (struct productPoly), struct productPoly *);
    pP->num = counterProd;
    MALCHOKE(pP->product, counterProd * sizeof (Polynomial *),Polynomial **);
    MALCHOKE(pP->exponent, counterProd * sizeof (int),int *);
    // Copy the collected terms and exponents into the product polynomial structure
    for (i = 0; i < counterProd; i++) {
      pP->product[i] = pProd[i];
      pP->exponent[i] = exponentProd[i];
      pP->product[i]->valid |= VALID_REF_FLAG;
    }
    rp->e.p = pP;
    rp->index = productCount;
    rp->id = nodeId;
    rp->key = key;
    rp->valid = VALID_NOTDISC_FLAG;
    rp->count = 0;
#ifdef USE_GMP
  mpf_init (rp->mpfValue);
#endif
#ifdef SOURCEDIGRAPH
    rp->source = findOrAddSource (fileName, lineNo, T_PRODUCT);
#endif

    // After the new polynomial is built, it is saved in the product polynomial list
    if (productCount >= productListLength) {
      productListLength += PRODUCT_LIST_INCREASE;
      productPListExpansions++;
      REALCHOKE (productList, productListLength * sizeof (Polynomial *), void *);
    }
    /* We either use a new position in the product list for this polynomial, or
     * replace an existing polynomial with the newly created polynomial at the position
     * occupied by the existing polynomial. */
    if (flag != 0 && p0EType == T_PRODUCT && p0Valid == 0) {
      // Use the resource of the freed polynomial for the newly-constructed polynomial
      productListReplacementCount++;
      productList[p0Index] = rp;
      productList[p0Index]->index = p0Index;
      productList[p0Index]->id = p0Id;
    } else {
      // Save the newly-constructed polynomial in the polynomial list
      productListNewCount++;
      productList[productCount] = rp;
      productList[productCount]->index = productCount;
      productList[productCount]->id = nodeId;
      if (polynomialDebugLevel >= 40) {
        fprintf (stderr, "Polynomial %ld, (product %d) added: ", nodeId, productCount);
        expTermPrinting (stderr, rp, 1);
        fprintf (stderr, "\n");
      }
      if (nodeId == polynomialLostNodeId) {
        fprintf (stderr, "nodeId %d has valid of %d and count of %d\n", polynomialLostNodeId, rp->valid, rp->count);
	raise (SIGINT);
      }
      productList[productCount]->valid |= VALID_TOP_FLAG;       // Currently unreferenced

      productCount++;
      nodeId++;
    }

    // The new polynomial is also recorded in the hash table
    if (flag != 0 && p0EType == T_PRODUCT && p0Valid == 0) {
      if (p0SubHIndex != location || p0HIndex != hIndex) {
        insertHashTable (&productHash[hIndex], location, key, p0Index);
        if (p0HIndex != hIndex) {
          deleteHashTable (&productHash[p0HIndex], p0SubHIndex);
        } else {
          if (p0SubHIndex < location) {
            deleteHashTable (&productHash[p0HIndex], p0SubHIndex);
          } else {
            deleteHashTable (&productHash[p0HIndex], p0SubHIndex + 1);
          }
        }
      }
      /* If the indexes and sub indexes of the old and new polynomials in the hash
       * table are identical, we don't need to do any more, just replace the key
       * of the old polynomial with the key of the new polynomial. */
      else
        productHash[hIndex].key[location] = key;
    } else
      // Just insert the newly-created polynomial into the hash table
      insertHashTable (&productHash[hIndex], location, key, productCount - 1);

    // If the factor is 1, return a product polynomial
    if (factor == 1.0) {
      productReturnNormalCount++;
      if (polynomialDebugLevel >= 60) {
        fprintf (stderr, "Returning a new product\n");
        if (polynomialDebugLevel >= 70) {
          expPrinting (rp);
          fprintf (stderr, "\n");
        }
      }
      return polyReturnWrapper (rp);
    } else {
      // If the factor is not 1, return a sum polynomial
      productNon1FactorIsSumCount++;
      if (polynomialDebugLevel >= 60)
        fprintf (stderr, "Returning new product via plusExp\n");
      return (plusExp (fileName, lineNo, 1, factor, polyReturnWrapper(rp), 0 /* Already handled flag */));
    }
  }
};

/**

   Create a key for a function call polynomial using the name of the 
   called function and the parameters.

*/
inline int keyFunctionCallPolynomial (char *fName, Polynomial ** p, int num)
{
  int key = 0;
  int i;

  // Compute a key from function name...
  for (i = 0; i < strlen (fName); i++)
    key += (int) fName[i];
  //...and parameters
  for (i = 0; i < num - 1; i++) {
    key = key + p[i]->key;
  }
  return key;
};

/**

  This function create a function call polynomial.  It accepts the function name 
  and parameters of the function.  If either the function or the parameters are 
  different, a new function call polynomial is created and returned. If everything 
  is the same as an existing function call polynomial, the existing one is returned.

  The parameters are:
  -# the number of parameters to the called function plus 1
  -# the name of the called function
  -# a group of parameters to the called function

*/
Polynomial *functionCallExp (int num, ...)
{
  int i, k;
  va_list args; ///< Variable list of parameters
  Polynomial *rp;       ///< New product polynomial
  struct functionPoly *fP;  ///< Function call structure for the new polynomial
  Polynomial **p;       ///< Parameters of the newly-built function call polynomial
  int fIndex, hIndex;   ///< Index in the polynomial list, index in the hash table
  char *fName;  ///< Function name
  int key = 0;      ///< Key of the new polynomial
  int first,    ///< First term for comparison
    last,       ///< Last term for comparison
    location;   ///< Position of the new polynomial in a hash bucket

#ifdef POLYSTATISTICS
  if ((termCount++ & 0x7FFFFF) == 0) ///< New term from an external standpoint.
    polyStatistics ("At 8M raw term multiple");
#endif

  // Get the number of parameters for this function call
  va_start (args, num);
  fName = va_arg (args, char *);

  MALCHOKE(p, (num - 1) * sizeof (Polynomial *),Polynomial **);

  // num is equal to 1 plus the number of parameters for the called function
  for (i = 0; i < num - 1; i++) {
    p[i] = va_arg (args, Polynomial *);
    if (p[i]->eType == T_OFFLINE) importPoly (p[i]);
    if (polynomialDebugLevel >= 60) {
      fprintf (stderr, "In functionCallExp item No. %d of %d type=%d: ", i + 1, num, p[i]->eType);
      expPrinting (p[i]);
      fprintf (stderr, "\n");
    }
  }
  va_end (args);

  // Compute a key for this polynomial
  key = keyFunctionCallPolynomial (fName, p, num);

  // Compute the index in the hash table according to the key

  hIndex = key % FUNCTIONCALL_HASH_SIZE;
  if (hIndex < 0)
    hIndex += FUNCTIONCALL_HASH_SIZE;

  /* Compare the key of this polynomial with the keys of other polynomials 
     already in the function list. */
  if (functionCallHash[hIndex].num > 0) {
    /* If the key of this function call is equal to the keys of polynomials in the
       function call list, see if the value of the new function call is equal to
       the value of an existing function call. */
    if (searchHashTable (&functionCallHash[hIndex], &first, &last, key)) {
      // Compare the new function call polynomials with a set of possible matching polynomials
      for (i = first; i <= last; i++) {
        fIndex = functionCallHash[hIndex].index[i];
        /* If the names of the called functions are identical and the numbers of parameters are 
	   identical, we compare their parameters. */
        if (strcmp (fName, functionCallList[fIndex]->e.f->name) == 0 && num - 1 == functionCallList[fIndex]->e.f->num) {
          // Compare the two function calls item by item
          for (k = 0; k < num - 1; k++)
            if (p[k] != functionCallList[fIndex]->e.f->para[k])
              break;
          // The two function calls are the same, return the existing one
          if (k >= num - 1) {
            free (p);
            functionHashHits++;
            return polyReturnWrapper (functionCallList[fIndex]);
          }
        }
      }
      location = last;
    } else
      location = last;
  } else
    location = 0;

  // The function call ws not found in the list, insert it in the list.
  // Build a new polynomial
  MALCHOKE(rp, sizeof (Polynomial), Polynomial *);
  MALCHOKE(fP, sizeof (struct functionPoly), struct functionPoly *);
  rp->eType = T_FUNCTIONCALL;
  fP->num = num - 1;
  MALCHOKE(fP->para, (num - 1) * sizeof (Polynomial *), Polynomial **);
  MALCHOKE(fP->name, strlen (fName) + 1, char *);
  strcpy (fP->name, fName);
  for (i = 0; i < fP->num; i++) {
    fP->para[i] = p[i];
  }
  rp->e.f = fP;
  rp->index = functionCallCount;
  rp->id = nodeId;
  rp->key = key;
  rp->valid = VALID_NOTDISC_FLAG;
  rp->count = 0;
#ifdef USE_GMP
  mpf_init (rp->mpfValue);
#endif

  // Insert the newly-built polynomial in function call list
  if (functionCallCount >= functionCallListLength) {
    functionCallListLength += FUNCTIONCALL_LIST_INCREASE;
    functionCallPListExpansions++;
    REALCHOKE (functionCallList, functionCallListLength * sizeof (Polynomial *), void *);
  }
  functionCallList[functionCallCount] = rp;
  functionCallCount++;
  if (polynomialDebugLevel >= 40)
    fprintf (stderr, "Polynomial %ld, (function %d) added\n", nodeId, functionCallCount);
  nodeId++;

  // Insert the polynomial in the hash table of function call polynomials
  insertHashTable (&functionCallHash[hIndex], location, key, functionCallCount - 1);

  free (p);

  return polyReturnWrapper (rp);

};

/**

   This function initialize a structure for creating a polynomial list which is used
   for polynomial evaluation.

*/
struct polyList *buildPolyList ()
{
  struct polyList *l;

  MALCHOKE(l, sizeof (struct polyList), struct polyList *);
  l->listSize = 1000;
  l->listNext = 0;
  MALCHOKE(l->pList, sizeof (Polynomial *) * l->listSize, Polynomial **);

  return l;
};

/**

  This function appends a polynomial onto a polynomial list.

*/
void polyListAppend (struct polyList *l, Polynomial * p)
{
  // Valid is a flag indicating that this polynomial appears on a sorting list.
  p->valid |= VALID_EVAL_FLAG;

  if (l->listNext >= l->listSize) {
    REALCHOKE (l->pList, sizeof (Polynomial *) * (l->listSize + 1000), void *);
    l->listSize = l->listSize + 1000;
  }
  l->pList[l->listNext] = p;
  l->listNext++;
};

/**

  Sort the polynomial list to prepare for polynomial evaluation.  Polynomials are
  evaluated in dependency order.  This function determines the order for evaluation.
  Shared terms are evaluated only once.  Together with polyListAppend(), this function
  constructs an evaluation list for evaluation of the polynomial.  For linkage computation,
  it allows shared polynomials within each likelihood polynomial and across likelihood
  polynomials of a set of pedigrees fully reused

*/
void doPolyListSorting (Polynomial * p, struct polyList *l)
{
  int i;

  switch (p->eType) {
  case T_OFFLINE:
    importPoly(p);
    // Notice that there's no break here...

    // If the polynomial is a constant, put it in the evaluation list
  case T_CONSTANT:
    polyListAppend (l, p);
    break;

    // If the polynomial is a variable, put it in the evaluation list
  case T_VARIABLE:
    if (!(p->valid & VALID_EVAL_FLAG)) {
      polyListAppend (l, p);
    }
    break;

    // If the polynomial is an external, put it in the evaluation list after all variables
  case T_EXTERNAL:
    if (p->valid & VALID_EVAL_FLAG)
      break;

    for (i = 0; i < variableCount; i++)
      if (!(variableList[i]->valid & VALID_EVAL_FLAG))
        doPolyListSorting (variableList[i], l);

    polyListAppend (l, p);
    break;

    // If the polynomial is a sum, put all the terms of the sum in the evaluation list
    //except constants and then put the sum in the evaluation list
  case T_SUM:

    if (p->valid & VALID_EVAL_FLAG)
      break;

    importTermList (p);
    for (i = 0; i < p->e.s->num; i++) {
      if (p->e.s->sum[i]->eType == T_OFFLINE) importPoly (p->e.s->sum[i]);
      if (p->e.s->sum[i]->eType != T_CONSTANT && (!(p->e.s->sum[i]->valid & VALID_EVAL_FLAG))) {
        doPolyListSorting (p->e.s->sum[i], l);
      }
    }
    exportTermList (p, FALSE);
    polyListAppend (l, p);
    break;

    /* If the polynomial is a product, put all the terms of the product in the
       evaluation list except constants and then put the product in the 
       evaluation list. */
  case T_PRODUCT:

    if (p->valid & VALID_EVAL_FLAG)
      break;

    for (i = 0; i < p->e.p->num; i++) {
      if (p->e.p->product[i]->eType == T_OFFLINE) importPoly (p->e.p->product[i]);
      if (p->e.p->product[i]->eType != T_CONSTANT && (!(p->e.p->product[i]->valid & VALID_EVAL_FLAG))) {
        doPolyListSorting (p->e.p->product[i], l);
      }
    }
    polyListAppend (l, p);
    break;

    /* If the polynomial is a functionCall, put the parameters in the evaluation list 
       and then put the functionCall in the evaluation list. */
  case T_FUNCTIONCALL:

    if (p->valid & VALID_EVAL_FLAG)
      break;

    for (i = 0; i < p->e.f->num; i++) {
      if (p->e.f->para[i]->eType == T_OFFLINE) importPoly (p->e.f->para[i]);
      if (p->e.f->para[i]->eType != T_CONSTANT && (!(p->e.f->para[i]->valid & VALID_EVAL_FLAG))) {
        doPolyListSorting (p->e.f->para[i], l);
      }
    }
    polyListAppend (l, p);
    break;

  default:
    break;
  }
}

/**

  Wrapper for doPolyListSorting to avoid recursing initialization and
  such.

*/
void polyListSorting (Polynomial * p, struct polyList *l)
{
  polyListSortingCount++;
#ifdef DEPENDENCYFLAGGING
  dependencyFlagging (p);
#endif
  //  writePolyDigraph(p);
  /* Clear all of the VALID_EVAL_FLAGs */
  clearValidEvalFlag ();
  doPolyListSorting (p, l);
}

/**

  This function compute the value of a polynomial.  It evaluate
  all the polynomials in the evaluation list for this polynomial.
  After the polynomials are completed, this is the only function
  to be executed for evaluation.  This function is not recursive,
  and is typically not used unless only a single evaluation is to
  be performed.

*/
void evaluatePoly (Polynomial * pp, struct polyList *l, double *pReturnValue)
{
  Polynomial *p;
  struct sumPoly *sP;
  struct productPoly *pP;
  register int i, j;
  double v, term;
#ifdef USE_GMP
  mpf_t mpfV, mpfTerm;
  mpf_init (mpfV);
  mpf_init (mpfTerm);
#endif

  if ((evaluatePolyCount++ & 0x3FF) == 0)
    DETAIL (0, "Evaluating polynomials, currently at %1.2g iterations", (double) evaluatePolyCount);

#ifdef EVALUATESW
  swStart (evaluatePolySW);
#endif

  if (polynomialDebugLevel >= 10)
    fprintf (stderr, "Starting evaluatePoly...\n");
  if (l->listNext == 0) {
    if (polynomialDebugLevel >= 10)
      fprintf (stderr, "...finished evaluatePoly with a zero!\n");
    *pReturnValue = pp->value;
#ifdef EVALUATESW
    swStop (evaluatePolySW);
#endif
    fprintf (stderr, "Evaluation %u of polynomial gives %.16g\n", evaluatePolyCount, pp->value);
    return;
  }

  for (j = 0; j <= l->listNext - 1; j++) {
    p = l->pList[j];
    switch (p->eType) {
    case T_OFFLINE:
      importPoly(p);
      // Notice that there's no break here...

    case T_CONSTANT:
#ifdef USE_GMP
      mpf_set_d (p->mpfValue, p->value);
#endif
      break;
      //Read the value of the variable
    case T_VARIABLE:
      if (p->e.v->vType == 'D')
        p->value = *(p->e.v->vAddr.vAddrD);
      else if (p->e.v->vType == 'I')
	p->value = *(p->e.v->vAddr.vAddrI);
      else
	FATAL ("Unknown variable type '%c'", p->e.v->vType);
#ifdef USE_GMP
      mpf_set_d (p->mpfValue, p->value);
#endif
      break;

    case T_EXTERNAL:
#ifdef POLYUSE_DL
      if (!p->e.e->entryOK)
	if (!loadPolyDL (p))
	  FATAL ("Cannot (re)load DL-based polynomial %s for evaluatePoly",
		   p->e.e->polynomialFunctionName);
      p->value = p->e.e->polynomialFunctionRoutine (1, variableList);
#endif
      break;

      //Sum up all the items in a sum
    case T_SUM:
      v = 0;
#ifdef USE_GMP
      mpf_set_ui (mpfV, 0UL);
#endif
      importTermList (p);
      sP = p->e.s;
      for (i = 0; i < sP->num; i++) {
        if (sP->factor[i] == 1)
          term = sP->sum[i]->value;
        else
          term = sP->sum[i]->value * sP->factor[i];
#ifdef USE_GMP
        if (sP->factor[i] == 1)
          mpf_set (mpfTerm, sP->sum[i]->mpfValue);
        else {
	  mpf_set_d (mpfTerm, sP->factor[i]);
          mpf_mul (mpfTerm, sP->sum[i]->mpfValue, mpfTerm);
	}
#endif
	
	// I think that this is a big deal in a signal processing context.
	/*
	int vPower, termPower;

	frexp(v, &vPower);
	frexp(term, &termPower);
	if ((v != 0) && (term != 0) && (abs(vPower - termPower) > 52)) {
	  fprintf (stderr, "Total precision loss (abs(%d-%d) = %d of 52 bits) in addition %d of %g and %g from  %d terms:\n", 
		   abs(vPower), abs(termPower), abs(vPower - termPower), i, v, term, sP->num);
	  for (j = 0; j < sP->num; j++)
	    fprintf (stderr, "%c %g*%g ",  j > 0 ? '+' : ' ', sP->factor[j], sP->sum[j]->value);
	  fprintf (stderr, "\n");
	}
	*/
	v += term;
#ifdef USE_GMP
	mpf_add (mpfV, mpfV, mpfTerm);
#endif
      }
      exportTermList (p, FALSE);
      p->value = v;
#ifdef USE_GMP
      mpf_set (p->mpfValue, mpfV);
#endif
      break;

      //Multiply all the items
    case T_PRODUCT:
      pP = p->e.p;
      v = 1;
      for (i = 0; i < pP->num; i++) {
        switch (pP->exponent[i]) {
        case 1:
          v *= pP->product[i]->value;
          break;
        case 2:
          v *= pP->product[i]->value * pP->product[i]->value;
          break;
        case 3:
          v *= pP->product[i]->value * pP->product[i]->value * pP->product[i]->value;
          break;
        case 4:
          v *= pP->product[i]->value * pP->product[i]->value * pP->product[i]->value * pP->product[i]->value;
          break;
        default:
          v *= pow (pP->product[i]->value, pP->exponent[i]);
          break;
        }
      }
      p->value = v;
#ifdef USE_GMP
      mpf_set_ui (mpfV, 1UL);
      for (i = 0; i < pP->num; i++) {
	mpf_pow_ui (mpfTerm, pP->product[i]->mpfValue, pP->exponent[i]);
	mpf_mul (mpfV, mpfV, mpfTerm);
      }
      mpf_set (p->mpfValue, mpfV);
#endif
      break;

      //Function calls are evaluated by calling the referred functions.
      //The referred function must be included in the linked library.
      //Otherwise, the program will exit.
    case T_FUNCTIONCALL:
      if (strcmp (p->e.f->name, "log10") == 0) {
        p->value = log10 (p->e.f->para[0]->value);
      } else if (strcmp (p->e.f->name, "log") == 0) {
        p->value = log (p->e.f->para[0]->value);
      } else if (strcmp (p->e.f->name, "tanh") == 0) {
        p->value = tanh (p->e.f->para[0]->value);
      } else if (strcmp (p->e.f->name, "atanh") == 0) {
        p->value = atanh (p->e.f->para[0]->value);
	/*
      } else if (strcmp (p->e.f->name, "gsl_ran_tdist_pdf") == 0) {
        p->value = gsl_ran_tdist_pdf (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_tdist_Q") == 0) {
        p->value = gsl_cdf_tdist_Q (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_tdist_P") == 0) {
        p->value = gsl_cdf_tdist_P (p->e.f->para[0]->value, p->e.f->para[1]->value);
	*/
      } else if (strcmp (p->e.f->name, "gsl_ran_tdist_pdf") == 0) {
        p->value = t_pdf_30 (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_tdist_Q") == 0) {
        p->value = ((double) 1.0) - t_cdf (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_tdist_P") == 0) {
        p->value = t_cdf (p->e.f->para[0]->value, p->e.f->para[1]->value);
	/*
      } else if (strcmp (p->e.f->name, "gsl_ran_ugaussian_pdf") == 0) {
        p->value = gsl_ran_ugaussian_pdf (p->e.f->para[0]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_ugaussian_Q") == 0) {
        p->value = gsl_cdf_ugaussian_Q (p->e.f->para[0]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_ugaussian_P") == 0) {
        p->value = gsl_cdf_ugaussian_P (p->e.f->para[0]->value);
	*/
      } else if (strcmp (p->e.f->name, "gsl_ran_ugaussian_pdf") == 0) {
        p->value = gaussian_pdf (p->e.f->para[0]->value, (double) 0.0, (double) 1.0);
      } else if (strcmp (p->e.f->name, "gsl_cdf_ugaussian_Q") == 0) {
        p->value = ((double) 1.0) - gaussian_cdf (p->e.f->para[0]->value, (double) 0.0, (double) 1.0);
      } else if (strcmp (p->e.f->name, "gsl_cdf_ugaussian_P") == 0) {
        p->value = gaussian_cdf (p->e.f->para[0]->value, (double) 0.0, (double) 1.0);
	/*
      } else if (strcmp (p->e.f->name, "gsl_cdf_chisq_P") == 0) {
        p->value = gsl_cdf_chisq_P (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_chisq_Q") == 0) {
        p->value = gsl_cdf_chisq_Q (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_ran_chisq_pdf") == 0) {
        p->value = gsl_ran_chisq_pdf (p->e.f->para[0]->value, p->e.f->para[1]->value);
	*/
      } else if (strcmp (p->e.f->name, "gsl_cdf_chisq_P") == 0) {
        p->value = chisq_cdf (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_chisq_Q") == 0) {
        p->value = ((double) 1.0) - chisq_cdf (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_ran_chisq_pdf") == 0) {
        p->value = chisq_pdf (p->e.f->para[0]->value, p->e.f->para[1]->value);

      } else if (strcmp (p->e.f->name, "pow") == 0) {
        p->value = pow (p->e.f->para[0]->value, p->e.f->para[1]->value);

      } else if (strcmp (p->e.f->name, "exp") == 0) {
        p->value = exp (p->e.f->para[0]->value);
      } else if (strcmp (p->e.f->name, "sqrt") == 0) {
        p->value = sqrt (p->e.f->para[0]->value);
      } else
        FATAL ("Unknown function name %s in polynomial", p->e.f->name);
#ifdef USE_GMP
      mpf_set_d (p->mpfValue, p->value);
#endif
      break;

    default:
      FATAL ("In evaluatePoly, unknown expression type: [%d]", p->eType);
      break;
    }
    if (isnan (p->value))
      ERROR ("In evaluatePoly, evaluated value of type %d as not a number (NaN)", p->eType);
    //    expTermPrinting (stderr, p, 1);
    //    fprintf (stderr, "\n%s[%d] = %g\n", eTypes[p->eType], p->id, p->value);

#ifdef USE_GMP
    //    fprintf (stderr, "...%s[%d] is %.18g, or ", eTypes[p->eType], p->id, p->value);
    //    mpf_out_str (stderr, 10, 50, p->mpfValue);
    //    fprintf (stderr, " aka %.18g\n", mpf_get_d (p->mpfValue));
#endif

  }
  if (polynomialDebugLevel >= 10)
    fprintf (stderr, "...finished evaluatePoly with %G\n", pp->value);
  *pReturnValue = pp->value;
#ifdef EVALUATESW
  swStop (evaluatePolySW);
#endif

#ifdef USE_GMP
  mpf_clear (mpfV);
  mpf_clear (mpfTerm);
  fprintf (stderr, "Evaluated %s[%d] to %.18g, or ", eTypes[pp->eType], pp->id, pp->value);
  mpf_out_str (stderr, 10, 50, pp->mpfValue);
  fprintf (stderr, " aka %.18g\n", mpf_get_d (pp->mpfValue));
  if ((fabs(mpf_get_d (pp->mpfValue) - pp->value) * 100.0 / pp->value) > .000001) {
    fprintf (stderr, "Error greater than 8 places, aborting!\n");
    exit (EXIT_FAILURE);
  }
#endif

  return;
}

/**

  This is the first function to call before using polynomials.  It allocates some memory
  and initializes important variables.

*/
void polynomialInitialization (int newPolynomialScale)
{
  int i;
  char *envVar;
  char messageBuffer[MAXSWMSG];

  // Get our maximum so we can react to depletion.
  maximumPMK = swGetMaximumPMK ();

#ifdef FREEDEBUG
#warning "freePoly protection is turned on, so memory will purposefully leak!"
  swLogMsg (stdout, "freePoly protection is turned on, so memory will purposefully leak!");
#endif

  if ((envVar = getenv ("polynomialDebugLevel")) != NULL) {
    polynomialDebugLevel = atoi (envVar);
  }
  if (polynomialDebugLevel > 0)
    fprintf (stderr, "polynomialDebugLevel is at %d\n", polynomialDebugLevel);

  /** polynomialScale sets the initial size of polynomial management lists, and 
  sets the permanent size of the polynomial hash tables. Kelvin parses a value 
  for it from the 'PE' directive. The default is 10, which corresponds to the 
  original maximum reasonable sizes for all data structures, and consumes just 
  under a gigabyte of memory.  Obviously small-memory environments benefit from 
  setting it to a lower value. */

  if (newPolynomialScale <= 0)
    polynomialScale = 1;
  else if (newPolynomialScale > 10)
    polynomialScale = 10;
  else
    polynomialScale = newPolynomialScale;
  sprintf (messageBuffer, "polynomialScale is %d (1-10, 1 is default)", polynomialScale);
  swLogMsg (stdout, messageBuffer);

  evaluatePolySW = swCreate ("evaluatePoly");
  evaluateValueSW = swCreate ("evaluateValue");

#ifdef USE_SSD
  initSSD ();
  for (i=0; i<IMTL_COUNT; i++) {
    iMTL[i].cT.doublePairCount = 0;
    iMTL[i].lastNodeId = 0;
  }
#endif

#ifdef USE_GMP
  mpf_set_default_prec (192);
#endif

  /* Scale all initial and growth sizes up by the polynomial scale (default 10) */

  CONSTANT_HASH_SIZE = MIN_CONSTANT_HASH_SIZE * polynomialScale;
  VARIABLE_HASH_SIZE = MIN_VARIABLE_HASH_SIZE * polynomialScale;
  SUM_HASH_SIZE = MIN_SUM_HASH_SIZE * polynomialScale;
  PRODUCT_HASH_SIZE = MIN_PRODUCT_HASH_SIZE * polynomialScale;
  FUNCTIONCALL_HASH_SIZE = MIN_FUNCTIONCALL_HASH_SIZE * polynomialScale;
  HASH_TABLE_INCREASE = MIN_HASH_TABLE_INCREASE * polynomialScale;
  maxHashLength = HASH_TABLE_INCREASE;  /* Set to actual default. */
  CONSTANT_LIST_INITIAL = MIN_CONSTANT_LIST_INITIAL * polynomialScale;
  CONSTANT_LIST_INCREASE = MIN_CONSTANT_LIST_INCREASE * polynomialScale;
  VARIABLE_LIST_INITIAL = MIN_VARIABLE_LIST_INITIAL * polynomialScale;
  VARIABLE_LIST_INCREASE = MIN_VARIABLE_LIST_INCREASE * polynomialScale;
  EXTERNAL_LIST_INITIAL = MIN_EXTERNAL_LIST_INITIAL * polynomialScale;
  EXTERNAL_LIST_INCREASE = MIN_EXTERNAL_LIST_INCREASE * polynomialScale;
  SUM_LIST_INITIAL = MIN_SUM_LIST_INITIAL * polynomialScale;
  SUM_LIST_INCREASE = MIN_SUM_LIST_INCREASE * polynomialScale;
  PRODUCT_LIST_INITIAL = MIN_PRODUCT_LIST_INITIAL * polynomialScale;
  PRODUCT_LIST_INCREASE = MIN_PRODUCT_LIST_INCREASE * polynomialScale;
  FUNCTIONCALL_LIST_INITIAL = MIN_FUNCTIONCALL_LIST_INITIAL * polynomialScale;
  FUNCTIONCALL_LIST_INCREASE = MIN_FUNCTIONCALL_LIST_INCREASE * polynomialScale;

  if ((envVar = getenv ("polynomialLostNodeId")) != NULL) {
    polynomialLostNodeId = atoi (envVar);
  }
  if (polynomialLostNodeId > 0)
    fprintf (stderr, "polynomialLostNodeId is %d\n", polynomialLostNodeId);

  nodeId = 0;
  termCount = 0;

  /* Allocate memory for polynomial list of each type of polynomials, set the counter of each
     type of polynomials to be 0. */
  constantListLength = CONSTANT_LIST_INITIAL;
  constantCount = 0;
  MALCHOKE(constantList, constantListLength * sizeof (Polynomial *), Polynomial **);

  variableListLength = VARIABLE_LIST_INITIAL;
  variableCount = 0;
  MALCHOKE(variableList, variableListLength * sizeof (Polynomial *), Polynomial **);

  externalListLength = EXTERNAL_LIST_INITIAL;
  externalCount = 0;
  MALCHOKE(externalList, externalListLength * sizeof (Polynomial *), Polynomial **);

  sumListLength = SUM_LIST_INITIAL;
  sumCount = 0;
  MALCHOKE(sumList, sumListLength * sizeof (Polynomial *), Polynomial **);

  productListLength = PRODUCT_LIST_INITIAL;
  productCount = 0;
  MALCHOKE(productList, productListLength * sizeof (Polynomial *), Polynomial **);

  functionCallListLength = FUNCTIONCALL_LIST_INITIAL;
  functionCallCount = 0;
  MALCHOKE(functionCallList, functionCallListLength * sizeof (Polynomial *), Polynomial **);

  // Allocate memory for hash tables, each type of polynomials has its own hash table
  MALCHOKE(constantHash, CONSTANT_HASH_SIZE * sizeof (struct hashStruct), struct hashStruct *);
  MALCHOKE(variableHash, VARIABLE_HASH_SIZE * sizeof (struct hashStruct), struct hashStruct *);
  MALCHOKE(sumHash, SUM_HASH_SIZE * sizeof (struct hashStruct), struct hashStruct *);
  MALCHOKE(productHash, PRODUCT_HASH_SIZE * sizeof (struct hashStruct), struct hashStruct *);
  MALCHOKE(functionCallHash, FUNCTIONCALL_HASH_SIZE * sizeof (struct hashStruct), struct hashStruct *);

  // Initialize the constant hash table, pre-allocate memory for recording polynomials
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < CONSTANT_HASH_SIZE; i++) {
    constantHash[i].num = 0;
    constantHash[i].length = HASH_TABLE_INCREASE;
    MALCHOKE(constantHash[i].index, constantHash[i].length * sizeof (int),int *);
    MALCHOKE(constantHash[i].key, constantHash[i].length * sizeof (int),int *);
  }

  // Initialize the variable hash table, pre-allocate memory for recording polynomials
  for (i = 0; i < VARIABLE_HASH_SIZE; i++) {
    variableHash[i].num = 0;
    variableHash[i].length = HASH_TABLE_INCREASE;
    MALCHOKE(variableHash[i].index, variableHash[i].length * sizeof (int),int *);
    MALCHOKE(variableHash[i].key, variableHash[i].length * sizeof (int),int *);
  }

  // Initialize the sum hash table, pre-allocate memory for recording polynomials
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < SUM_HASH_SIZE; i++) {
    sumHash[i].num = 0;
    sumHash[i].length = HASH_TABLE_INCREASE;
    MALCHOKE(sumHash[i].index, sumHash[i].length * sizeof (int),int *);
    MALCHOKE(sumHash[i].key, sumHash[i].length * sizeof (int),int *);
  }

  // Initialize the product hash table, pre-allocate memory for recording polynomials
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < PRODUCT_HASH_SIZE; i++) {
    productHash[i].num = 0;
    productHash[i].length = HASH_TABLE_INCREASE;
    MALCHOKE(productHash[i].index, productHash[i].length * sizeof (int),int *);
    MALCHOKE(productHash[i].key, productHash[i].length * sizeof (int),int *);
  }

  // Initialize the function call hash table, pre-allocate memory for recording polynomials
  for (i = 0; i < FUNCTIONCALL_HASH_SIZE; i++) {
    functionCallHash[i].num = 0;
    functionCallHash[i].length = HASH_TABLE_INCREASE;
    MALCHOKE(functionCallHash[i].index, functionCallHash[i].length * sizeof (int),int *);
    MALCHOKE(functionCallHash[i].key, functionCallHash[i].length * sizeof (int),int *);
  }

  initialHashSize = (sizeof (struct hashStruct) + (2 * HASH_TABLE_INCREASE * sizeof (int))) * (CONSTANT_HASH_SIZE + VARIABLE_HASH_SIZE + SUM_HASH_SIZE + PRODUCT_HASH_SIZE + FUNCTIONCALL_HASH_SIZE);

  // Acquire memory for containers for the terms of a sum polynomial

  // For variable polynomials
  containerLength_v1 = 100;
  MALCHOKE(factor_v1, containerLength_v1 * sizeof (double),double *);
  MALCHOKE(p_v1, containerLength_v1 * sizeof (Polynomial *),Polynomial **);

  // For product polynomials
  containerLength_p1 = 100;
  MALCHOKE(factor_p1, containerLength_p1 * sizeof (double),double *);
  MALCHOKE(p_p1, containerLength_p1 * sizeof (Polynomial *),Polynomial **);

  // For function call polynomials
  containerLength_f1 = 100;
  MALCHOKE(factor_f1, containerLength_f1 * sizeof (double),double *);
  MALCHOKE(p_f1, containerLength_f1 * sizeof (Polynomial *),Polynomial **);

  // Containers for organizing a sum polynomial
  lengthSum = 300;
  MALCHOKE(factorSum, lengthSum * sizeof (double),double *);
  MALCHOKE(pSum, lengthSum * sizeof (Polynomial *),Polynomial **);

  // Acquire memory for containers for the terms of a product polynomial

  // For variable polynomials
  containerLength_v2 = 100;
  MALCHOKE(exponent_v2, containerLength_v2 * sizeof (int),int *);
  MALCHOKE(p_v2, containerLength_v2 * sizeof (Polynomial *),Polynomial **);

  // For sum polynmials
  containerLength_s2 = 100;
  MALCHOKE(exponent_s2, containerLength_s2 * sizeof (int),int *);
  MALCHOKE(p_s2, containerLength_s2 * sizeof (Polynomial *),Polynomial **);

  // For function call polynomials
  containerLength_f2 = 100;
  MALCHOKE(exponent_f2, containerLength_f2 * sizeof (double),int *);
  MALCHOKE(p_f2, containerLength_f2 * sizeof (Polynomial *),Polynomial **);

  // Containers for organizing a product polynomial
  lengthProd = 300;
  MALCHOKE(exponentProd, lengthProd * sizeof (int),int *);
  MALCHOKE(pProd, lengthProd * sizeof (Polynomial *),Polynomial **);

}

/**

  Release all the memory occupied by the polynomials, polynomial lists, hash
  tables, etc. This has not been used for a while and so may not be entirely 
  up-to-date. I acknowledge that it should be resurrected for a variety of
  reasons, but it's so far from the top of the queue...

*/
void polynomialClearance ()
{
  int j;

  /* Hah! Exit works better. */
  return;

  nodeId = 0;
  // Clear constant polynomials
  for (j = 0; j < constantCount; j++)
    free (constantList[j]);
  constantCount = 0;
  // Clear variable polynomials
  for (j = 0; j < variableCount; j++) {
    free (variableList[j]->e.v);
    free (variableList[j]);
  }
  variableCount = 0;
  // Clear sum polynomials
  for (j = 0; j < sumCount; j++) {
    free (sumList[j]->e.s->sum);
    free (sumList[j]->e.s->factor);
    free (sumList[j]->e.s);
    free (sumList[j]);
  }
  sumCount = 0;
  // Clear product polynomials
  for (j = 0; j < productCount; j++) {
    free (productList[j]->e.p->product);
    free (productList[j]->e.p->exponent);
    free (productList[j]->e.p);
    free (productList[j]);
  }
  productCount = 0;
  // Clear function call polynomials
  for (j = 0; j < functionCallCount; j++) {
    free (functionCallList[j]->e.f->name);
    free (functionCallList[j]->e.f->para);
    free (functionCallList[j]->e.f);
    free (functionCallList[j]);
  }
  functionCallCount = 0;
  // Release memory occupied by constant hash table
  for (j = 0; j < CONSTANT_HASH_SIZE; j++) {
    if (constantHash[j].length > 0) {
      free (constantHash[j].index);
      free (constantHash[j].key);
    }
  }

  // Release memory occupied by variable hash table
  for (j = 0; j < VARIABLE_HASH_SIZE; j++) {
    if (variableHash[j].length > 0) {
      free (variableHash[j].index);
      free (variableHash[j].key);
    }
  }

  // Release memory occupied by sum hash table
  for (j = 0; j < SUM_HASH_SIZE; j++) {
    if (sumHash[j].length > 0) {
      free (sumHash[j].index);
      free (sumHash[j].key);
    }
  }

  // Release memory occupied by product hash table
  for (j = 0; j < PRODUCT_HASH_SIZE; j++) {
    if (productHash[j].length > 0) {
      free (productHash[j].index);
      free (productHash[j].key);
    }
  }

  // Release memory occupied by function call hash table
  for (j = 0; j < FUNCTIONCALL_HASH_SIZE; j++) {
    if (functionCallHash[j].length > 0) {
      free (functionCallHash[j].index);
      free (functionCallHash[j].key);
    }
  }

  // Release memory occupied by polynomial lists
  free (constantList);
  free (variableList);
  free (sumList);
  free (productList);
  free (functionCallList);

  // Release memory occupied by hash table
  free (constantHash);
  free (variableHash);
  free (sumHash);
  free (productHash);
  free (functionCallHash);

  // Release memory used for construction of sum polynomials
  free (factor_v1);
  free (p_v1);
  free (factor_p1);
  free (p_p1);
  free (factor_f1);
  free (p_f1);
  free (factorSum);
  free (pSum);

  // Release memory used for construction of product polynomials
  free (exponent_v2);
  free (p_v2);
  free (exponent_s2);
  free (p_s2);
  free (exponent_f2);
  free (p_f2);
  free (exponentProd);
  free (pProd);

}

/// Maximum expected number of tree tiers to track in summary of polynomial.
#define MAXPOLYTIERS 16
int polyTiers[MAXPOLYTIERS][5]; ///< Array for holding polynomial type counts by tier
int peakPolyTiers; ///< Maximum number of polynomial tiers expected (handled)
char *polyTypes[] = { "constant", "variable", "sum", "product", "function" }; ///< Polynomial type strings
/**

  Display a summarization of the specified polynomial as a number 
  of tiers with counts of each polynomial type at each tier.

*/
void doPrintSummaryPoly (Polynomial * p, int currentTier)
{
  int i;

  if (p->eType == T_OFFLINE) importPoly (p);
  if (p->valid & VALID_EVAL_FLAG)
    return;
  p->valid |= VALID_EVAL_FLAG;

  if (currentTier > peakPolyTiers)
    peakPolyTiers = currentTier;
  if (currentTier >= MAXPOLYTIERS)
    return;
  polyTiers[currentTier][p->eType]++;
  switch (p->eType) {
  case T_OFFLINE:
    importPoly(p);
    // Notice that there's no break here...

  case T_CONSTANT:
  case T_VARIABLE:
  case T_EXTERNAL:
    break;
  case T_SUM:
    importTermList (p);
    for (i = 0; i < p->e.s->num; i++) {
      doPrintSummaryPoly (p->e.s->sum[i], currentTier + 1);
    }
    exportTermList (p, FALSE);
    break;
  case T_PRODUCT:
    for (i = 0; i < p->e.p->num; i++) {
      doPrintSummaryPoly (p->e.p->product[i], currentTier + 1);
    }
    break;
  case T_FUNCTIONCALL:
    for (i = 0; i < p->e.f->num; i++) {
      doPrintSummaryPoly (p->e.f->para[i], currentTier + 1);
    }
    break;
  case T_FREED:
    expTermPrinting (stderr, p, 1);
    fprintf (stderr, "\n");
    FATAL ("In doPrintSummaryPoly, preceeding polynomial was freed");
    break;
  default:
    FATAL ("In doPrintSummaryPoly, unknown expression type: [%d]", p->eType);
  }
}

/**

  Wrapper for doPrintSummaryPoly to avoid recursing initialization and
  termination.

*/
void printSummaryPoly (Polynomial * p)
{
  int i, j;

  fprintf (stderr, "Summary of Remains of %ld Polynomials (from tree, unique terms only):\n", nodeId);

  clearValidEvalFlag ();
  memset (polyTiers, 0, sizeof (polyTiers));
  peakPolyTiers = 0;
  doPrintSummaryPoly (p, 0);
  for (i = 0; i <= peakPolyTiers; i++) {
    fprintf (stderr, "Tier %d:", i);
    int firstPrint = 1;

    for (j = 0; j < 5; j++) {
      if (polyTiers[i][j] != 0) {
        if (!firstPrint)
          fprintf (stderr, ",");
        firstPrint = 0;
        fprintf (stderr, " %d %s", polyTiers[i][j], polyTypes[j]);
        if (polyTiers[i][j] > 1)
          fprintf (stderr, "s");
      }
    }
    fprintf (stderr, "\n");
  }
  fprintf (stderr, "---\n");
}

#ifdef SOURCEDIGRAPH
/**

  This function produces a DOT-format (graphviz) data files named
  pSP.dot that illustrates term source parenting for the current
  set of polynomials. The data is collected for all polynomial build 
  calls.

*/
void dumpSourceParenting ()
{
  int i, j;
  qsort (polySources, polySourceCount, sizeof (struct polySource), compareSourcesByEntryNo);
  FILE *diGraph;

  if ((diGraph = fopen ("pSP.dot", "w")) == NULL) {
    perror ("Cannot open polynomial source parenting digraph file\n");
    exit (EXIT_FAILURE);
  }
  fprintf (diGraph, "digraph G {\n");
  for (i = 0; i < MAXPOLYSOURCES; i++)
    if (polySources[i].lineNo != 0)
      fprintf (diGraph, "%d [label=\"%s(%d)\"];\n", i, polySources[i].moduleName, polySources[i].lineNo);
  for (i = 0; i < MAXPOLYSOURCES; i++) {
    if (polySources[i].lineNo != 0) {
      for (j = 0; j < MAXPOLYSOURCES; j++) {
        if (i != j) {
          if (polySources[i].originalChildren[j] != 0)
            fprintf (diGraph, "%d -> %d [label=\"%d\"];\n", i, j, polySources[i].originalChildren[j]);
        }
      }
    }
  }
  fprintf (diGraph, "}\n");
  fclose (diGraph);
}

/**

  Full dump of all polynomial term-generating source locations along with the
  number of calls. Nice tuning diagnostic.

*/
void dumpPolySources ()
{
  int i;
  for (i = 0; i < polySourceCount; i++) {
    fprintf (stderr, "From %s line %d, %d type %d polynomials\n", polySources[i].moduleName, polySources[i].lineNo, polySources[i].totalCalls, polySources[i].eType);
  }
}
#endif

/**

  Produce a digraph for the specified polynomial itself.

*/
void doWritePolyDigraph (Polynomial * p, FILE * diGraph)
{
  int i;

  if (p->valid & VALID_EVAL_FLAG)
    return;
  p->valid |= VALID_EVAL_FLAG;

  switch (p->eType) {
  case T_OFFLINE:
    importPoly(p);
    // Notice that there's no break here...

  case T_CONSTANT:
    fprintf (diGraph, "%d [label=\"%G\"];\n", p->id, p->value);
    break;
  case T_VARIABLE:
    fprintf (diGraph, "%d [label=\"%s\"];\n", p->id, p->e.v->vName);
    break;
  case T_EXTERNAL:
    fprintf (diGraph, "%d [label=\"%s()\"];\n", p->id, p->e.e->polynomialFunctionName);
    break;
  case T_SUM:
    fprintf (diGraph, "%d [label=\"+\"];\n", p->id);
    importTermList (p);
    for (i = 0; i < p->e.s->num; i++) {
      doWritePolyDigraph (p->e.s->sum[i], diGraph);
      fprintf (diGraph, "%d -> %d;\n", p->id, p->e.s->sum[i]->id);
    }
    exportTermList (p, FALSE);
    break;
  case T_PRODUCT:
    fprintf (diGraph, "%d [label=\"*\"];\n", p->id);
    for (i = 0; i < p->e.p->num; i++) {
      doWritePolyDigraph (p->e.p->product[i], diGraph);
      fprintf (diGraph, "%d -> %d;\n", p->id, p->e.p->product[i]->id);
    }
    break;
  case T_FUNCTIONCALL:
    fprintf (diGraph, "%d [label=\"fn\"];\n", p->id);
    for (i = 0; i < p->e.f->num; i++) {
      doWritePolyDigraph (p->e.f->para[i], diGraph);
      fprintf (diGraph, "%d -> %d;\n", p->id, p->e.f->para[i]->id);
    }
    break;
  case T_FREED:
    expTermPrinting (stderr, p, 1);
    exit (EXIT_FAILURE);
    FATAL ("In doWritePolyDigraph, preceeding polynomial was freed");
    break;
  default:
    FATAL ("In doWritePolyDigraph, unknown expression type: [%d]", p->eType);
  }
}

/**

  Wrapper for doWritePolyDigraph to avoid recursing initialization and
  termination. Producesf files named pD_n.dot, where n is the polynomial
  root nodeId.

*/
void writePolyDigraph (Polynomial * p)
{
  FILE *diGraph;
  char fileName[PATH_MAX];

  sprintf (fileName, "pD_%d.dot", p->id);
  if ((diGraph = fopen (fileName, "w")) == NULL) {
    perror ("Cannot open polynomial digraph file\n");
    exit (EXIT_FAILURE);
  }
  fprintf (diGraph, "digraph G {\n");

  clearValidEvalFlag ();

  doWritePolyDigraph (p, diGraph);

  fprintf (diGraph, "}\n");
  fclose (diGraph);
}

void thrashingCheck ()
{
#ifdef THRASH_CHECK
  unsigned long deltaAccumWallTime, deltaAccumUserTime;
  char messageBuffer[MAXSWMSG];

  /* Now we check to see if we're thrashing... If we didn't get at least 10% CPU since the last
   * time we did statistics (provided that was at least a second ago), we should externalize polynomials,
   * or exit or take some evasive action. */

  swStop (overallSW);
  if (lastPDSAccumWallTime != 0) {
    deltaAccumWallTime = overallSW->swAccumWallTime - lastPDSAccumWallTime;
    deltaAccumUserTime = overallSW->swAccumRUSelf.ru_utime.tv_sec + overallSW->swAccumRUChildren.ru_utime.tv_sec - lastPDSAccumUserTime;
    if ((deltaAccumUserTime != 0) && (100 * deltaAccumUserTime / (deltaAccumWallTime ? deltaAccumWallTime : 1) < THRASH_CPU))
      WARNING ("Thrashing detected (utilization under %d%%), consider exiting!", THRASH_CPU);
  }
  swStart (overallSW);
  lastPDSAccumWallTime = overallSW->swAccumWallTime;
  lastPDSAccumUserTime = overallSW->swAccumRUSelf.ru_utime.tv_sec + overallSW->swAccumRUChildren.ru_utime.tv_sec;
#endif
  return;
}

void polyDynamicStatistics (char *title)
{

  fprintf (stderr, "Dynamic polynomial statistics (%s):\n", title);

#ifdef EVALUATESW
  if (evaluatePolyCount)
    swDump (evaluatePolySW);
  if (evaluateValueCount)
    swDump (evaluateValueSW);
#endif

  fprintf (stderr,
      "Counts/Hits: c=%d/%d(%2.1f:1), v=%d/%d, s=%d/%d(%2.1f:1), p=%d/%d(%2.1f:1), f=%d/%d\n",
      constantCount, constantHashHits,
      constantHashHits / (float) (constantCount ? constantCount : 1),
      variableCount, variableHashHits, sumCount, sumHashHits, sumHashHits / (float) (sumCount ? sumCount : 1), productCount, productHashHits, productHashHits / (float) (productCount ? productCount : 1), functionCallCount, functionHashHits);

  fprintf (stderr,
      "Lists: expansions(@size): c=%d(@%d), v=%d(@%d), s=%d(@%d), p=%d(@%d), f=%d(@%d), tcc=%d(@%lu)\n",
      constantPListExpansions, CONSTANT_LIST_INCREASE,
      variablePListExpansions, VARIABLE_LIST_INCREASE,
      sumPListExpansions, SUM_LIST_INCREASE, productPListExpansions, PRODUCT_LIST_INCREASE, functionCallPListExpansions, FUNCTIONCALL_LIST_INCREASE, containerExpansions, (unsigned long) sizeof (Polynomial *) + sizeof (double));

  fprintf (stderr,
      "...sizes: c=%lu, v=%lu, s=%lu, p=%lu, f=%lu\n",
      (unsigned long) sizeof (void *) * constantListLength,
      (unsigned long) sizeof (void *) * variableListLength, (unsigned long) sizeof (void *) * sumListLength, (unsigned long) sizeof (void *) * productListLength, (unsigned long) sizeof (void *) * functionCallListLength);

  fprintf (stderr,
	   "TermCount: %ld NodeId: %ld Hash: max len=%d, init size=%lu, SPL: eff=%lu%%, avg len=%lu\n", termCount, nodeId, maxHashLength, initialHashSize, 100 * (lowSPLCount + highSPLCount) / (totalSPLCalls ? totalSPLCalls : 1), totalSPLLengths / (totalSPLCalls ? totalSPLCalls : 1));

  fprintf (stderr,
      "Calls: pLS=%d eP=%d eV=%d hAP=%d kP=%d hP=%d uHP=%d fP=%d(%d) fKP=%d\n",
      polyListSortingCount, evaluatePolyCount, evaluateValueCount, holdAllPolysCount, keepPolyCount, holdPolyCount, unHoldPolyCount, freePolysAttemptCount, freePolysCount, freeKeptPolysCount);

  if (sumReleaseableCount == 0 && sumNotReleaseableCount == 0 && productReleaseableCount == 0 && productNotReleaseableCount == 0)
    return;

  fprintf (stderr, "Sum polys: try release=%d or not=%d return const=%d (same=%d) return 1-term=%d hash hits=%d\n", sumReleaseableCount, sumNotReleaseableCount, sumReturnConstantCount, sumReturnSameConstantCount, sumReturn1TermCount, sumHashHits);
  fprintf (stderr, "...really new=%d new on sumList=%d replaced on sumList=%d freed=%d 1st-term freed=%d\n", sumNewCount, sumListNewCount, sumListReplacementCount, sumFreedCount, sum1stTermsFreedCount);
  fprintf (stderr, "Product polys: try release=%d or not=%d return 0=%d return const=%d (same=%d) return 1st term=%d\n", productReleaseableCount, productNotReleaseableCount, productReturn0Count, productReturnConstantCount, productReturnSameConstantCount, productReturn1stTermCount);
  fprintf (stderr, "...actually a sum poly=%d hash hits=%d hash hit but sum=%d return normal (factor is 1)=%d\n", productReturn1TermSumCount, productHashHits, productHashHitIsSumCount, productReturnNormalCount);
  fprintf (stderr, "...factor not 1 now sum=%d new on productList=%d replaced on productList=%d\n", productNon1FactorIsSumCount, productListNewCount, productListReplacementCount);
  fprintf (stderr, "...freed=%d 1st-term freed=%d\n", productFreedCount, product1stTermsFreedCount);

  return;
}

/* 
   This is weighted to the high-end because although there
   are more short term lists, they are handled proportionately
   less often. This shifting approach is significantly faster
   than inequality comparison or masking when compiled without
   optimization, and only slightly slower when optimized.
*/
int hi16Bit (unsigned short value) {
  int i;
  for (i=15; i>0; i--)
    if (value & 0x8000)
      break;
    else
      value <<= 1;
  return i;
}

/**
 This function prints out polynomial statistic information.  It is mainly used for
 performance evaluation and debugging.

 The calculation of memory utilization is very accurate when there is no
 multithreading. Multithreading will add 15-20%.

*/
void polyStatistics (char *title)
{
  long constantSize, variableSize, sumSize, productSize, functionCallSize;
  double grandTotal;
  long sumTerms = 0, productTerms = 0, maxSumTerms = 0, maxProductTerms = 0;
  int sumTermCounts[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, 
    productTermCounts[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int constantHashSize = 0, constantHashPeak = 0, variableHashSize = 0, variableHashPeak = 0, sumHashSize = 0, sumHashPeak = 0, productHashSize = 0, productHashPeak = 0, functionCallHashSize = 0, functionCallHashPeak = 0;
  int i;

  polyDynamicStatistics (title);

  //  malloc_stats();
  /* VM calls have been known to block, so we're not going to do this in-line
   * if (swGetMaximumVMK () != 0) {
   * if (swGetCurrentVMK (getpid ()) > (0.95 * swGetMaximumVMK ())) {
   * fprintf (stderr,
   * "VM usage too high to permit list traversals for statistics.\n");
   * return;
   * }
   * }
   */

  fprintf (stderr, "Calculated polynomial statistics (%s):\n", title);

  constantSize = constantCount * sizeof (Polynomial);
  variableSize = variableCount * (sizeof (Polynomial) + sizeof (struct variablePoly));
  sumSize = sumCount * (sizeof (Polynomial) + sizeof (struct sumPoly));
  for (i = 0; i < sumCount; i++) {
    sumTermCounts[hi16Bit(sumList[i]->e.s->num)]++;
    sumTerms += sumList[i]->e.s->num;
    if (sumList[i]->e.s->num > maxSumTerms)
      maxSumTerms = sumList[i]->e.s->num;
  }
  productSize = productCount * (sizeof (Polynomial) + sizeof (struct productPoly));
  for (i = 0; i < productCount; i++) {
    productTermCounts[hi16Bit(productList[i]->e.p->num)]++;
    productTerms += productList[i]->e.p->num;
    if (productList[i]->e.p->num > maxProductTerms)
      maxProductTerms = productList[i]->e.p->num;
  }
  functionCallSize = functionCallCount * sizeof (Polynomial);

  fprintf (stderr, "Subterm count(avg): s=%ld(%ld), p=%ld(%ld), ", sumTerms, sumTerms / (sumCount ? sumCount : 1), productTerms, productTerms / (productCount ? productCount : 1));
  fprintf (stderr, "sizes (w/terms): c=%ld, s=%ld, p=%ld, f=%ld\n", constantSize, sumSize + (sumTerms * (sizeof (Polynomial *) + sizeof (double))), productSize + (productTerms * (sizeof (Polynomial *) + sizeof (int))), functionCallSize);
  fprintf (stderr, "Sum subterm dist <2^n: ");
  for (i=0; i<16; i++)
    fprintf (stderr, "%d=%d%s", i, sumTermCounts[i], i<15? "," : "\n");
  fprintf (stderr, "Product subterm dist <2^n: ");
  for (i=0; i<16; i++)
    fprintf (stderr, "%d=%d%s", i, productTermCounts[i], i<15? "," : "\n");

  constantHashSize = CONSTANT_HASH_SIZE * sizeof (struct hashStruct);
  for (i = 0; i < CONSTANT_HASH_SIZE; i++) {
    constantHashSize += constantHash[i].length * 2 * sizeof (int);
    if (constantHash[i].num > constantHashPeak)
      constantHashPeak = constantHash[i].num;
  }
  variableHashSize = VARIABLE_HASH_SIZE * sizeof (struct hashStruct);
  for (i = 0; i < VARIABLE_HASH_SIZE; i++) {
    variableHashSize += variableHash[i].length * 2 * sizeof (int);
    if (variableHash[i].num > variableHashPeak)
      variableHashPeak = variableHash[i].num;
  }
  sumHashSize = SUM_HASH_SIZE * sizeof (struct hashStruct);
  for (i = 0; i < SUM_HASH_SIZE; i++) {
    sumHashSize += sumHash[i].length * 2 * sizeof (int);
    if (sumHash[i].num > sumHashPeak)
      sumHashPeak = sumHash[i].num;
  }
  productHashSize = PRODUCT_HASH_SIZE * sizeof (struct hashStruct);
  for (i = 0; i < PRODUCT_HASH_SIZE; i++) {
    productHashSize += productHash[i].length * 2 * sizeof (int);
    if (productHash[i].num > productHashPeak)
      productHashPeak = productHash[i].num;
  }
  functionCallHashSize = FUNCTIONCALL_HASH_SIZE * sizeof (struct hashStruct);
  for (i = 0; i < FUNCTIONCALL_HASH_SIZE; i++) {
    functionCallHashSize += functionCallHash[i].length * 2 * sizeof (int);
    if (functionCallHash[i].num > functionCallHashPeak)
      functionCallHashPeak = functionCallHash[i].num;
  }
  fprintf (stderr, "Hash: size(peak length): c=%d(%d), v=%d(%d), s=%d(%d), p=%d(%d), f=%d(%d)\n",
      constantHashSize, constantHashPeak, variableHashSize, variableHashPeak, sumHashSize, sumHashPeak, productHashSize, productHashPeak, functionCallHashSize, functionCallHashPeak);

  grandTotal = constantHashSize + variableHashSize + sumHashSize + productHashSize +
      functionCallHashSize + constantSize + variableSize +
      sumSize + (sumTerms * (sizeof (Polynomial *) + sizeof (double))) +
      productSize + (productTerms * (sizeof (Polynomial *) + sizeof (int))) + functionCallSize + ((constantListLength + variableListLength + sumListLength + productListLength + functionCallListLength) * sizeof (void *));

  fprintf (stderr, "---Total data storage estimate: %.0fKb---\n", grandTotal / 1024);

  // Set a flag so we can take action as needed
  if ((grandTotal / 1024) > (maximumPMK / 4)) {
    memoryLow = TRUE;
#ifdef USE_SSD
    fprintf (stderr, "Maximum of %ldKb of physical memory is %.0f%% used, utilizing SSD.\n", 
	     maximumPMK, (grandTotal / 1024) / maximumPMK * 100);
#endif
  } else {
    memoryLow = FALSE;
#ifdef USE_SSD
    fprintf (stderr, "Maximum of %ldKb physical memory is at %.0f%% used.\n", 
	     maximumPMK, (grandTotal / 1024) / maximumPMK * 100);
#endif
  }

#ifdef DMTRACK
  swDumpHeldTotals ();
#endif
  return;
};

/**
 This function returns total polynomial memory utilization.

 The calculation of memory utilization is very accurate when there is no
 multithreading. Multithreading will add 15-20%.

*/
double polyMemory ()
{
  long constantSize, variableSize, sumSize, productSize, functionCallSize;
  double grandTotal;
  long sumTerms = 0, productTerms = 0;
  int constantHashSize = 0, variableHashSize = 0, sumHashSize = 0, 
    productHashSize = 0, functionCallHashSize = 0;
  int i;

  constantSize = constantCount * sizeof (Polynomial);
  variableSize = variableCount * (sizeof (Polynomial) + sizeof (struct variablePoly));

  sumSize = sumCount * (sizeof (Polynomial) + sizeof (struct sumPoly));
  for (i = 0; i < sumCount; i++)
    sumTerms += sumList[i]->e.s->num;

  productSize = productCount * (sizeof (Polynomial) + sizeof (struct productPoly));
  for (i = 0; i < productCount; i++)
    productTerms += productList[i]->e.p->num;

  functionCallSize = functionCallCount * sizeof (Polynomial);

  constantHashSize = CONSTANT_HASH_SIZE * sizeof (struct hashStruct);
  for (i = 0; i < CONSTANT_HASH_SIZE; i++)
    constantHashSize += constantHash[i].length * 2 * sizeof (int);

  variableHashSize = VARIABLE_HASH_SIZE * sizeof (struct hashStruct);
  for (i = 0; i < VARIABLE_HASH_SIZE; i++)
    variableHashSize += variableHash[i].length * 2 * sizeof (int);

  sumHashSize = SUM_HASH_SIZE * sizeof (struct hashStruct);
  for (i = 0; i < SUM_HASH_SIZE; i++)
    sumHashSize += sumHash[i].length * 2 * sizeof (int);

  productHashSize = PRODUCT_HASH_SIZE * sizeof (struct hashStruct);
  for (i = 0; i < PRODUCT_HASH_SIZE; i++)
    productHashSize += productHash[i].length * 2 * sizeof (int);

  functionCallHashSize = FUNCTIONCALL_HASH_SIZE * sizeof (struct hashStruct);
  for (i = 0; i < FUNCTIONCALL_HASH_SIZE; i++)
    functionCallHashSize += functionCallHash[i].length * 2 * sizeof (int);

  grandTotal = constantHashSize + variableHashSize + sumHashSize + productHashSize +
      functionCallHashSize + constantSize + variableSize +
      sumSize + (sumTerms * (sizeof (Polynomial *) + sizeof (double))) +
      productSize + (productTerms * (sizeof (Polynomial *) + sizeof (int))) + functionCallSize + ((constantListLength + variableListLength + sumListLength + productListLength + functionCallListLength) * sizeof (void *));

  return (grandTotal);
}

void printAllVariables ()
{
  int i, j;

  fprintf (stderr, "All variables (from hash):\n");

  if (variableCount > 0) {
    fprintf (stderr, "All %d variables:\n", variableCount);
    for (i = 0; i < VARIABLE_HASH_SIZE; i++) {
      if (variableHash[i].num <= 0)
        continue;
      for (j = 0; j < variableHash[i].num; j++) {
        fprintf (stderr, "(%d %d) index=%d key=%d count=%d valid=%d value=%.18g variable: ", i, j, variableHash[i].index[j], variableHash[i].key[j], variableList[variableHash[i].index[j]]->count, variableList[variableHash[i].index[j]]->valid, variableList[variableHash[i].index[j]]->value);
        expTermPrinting (stderr, variableList[variableHash[i].index[j]], 1);
        fprintf (stderr, "\n");
      }
    }
    fprintf (stderr, "\n");
  }
}

void printAllPolynomials ()
{
  int i, j;

  fprintf (stderr, "All Polynomials (from hash):\n");

  if (constantCount > 0) {
    fprintf (stderr, "All %d constants:\n", constantCount);
    for (i = 0; i < CONSTANT_HASH_SIZE; i++) {
      if (constantHash[i].num <= 0)
        continue;
      for (j = 0; j < constantHash[i].num; j++) {
        fprintf (stderr, "(%d %d) index=%d key=%d count=%d valid=%d constant: ", i, j, constantHash[i].index[j], constantHash[i].key[j], constantList[constantHash[i].index[j]]->count, constantList[constantHash[i].index[j]]->valid);
        expTermPrinting (stderr, constantList[constantHash[i].index[j]], 1);
        fprintf (stderr, "\n");
      }
    }
    fprintf (stderr, "\n");
  }

  printAllVariables ();

  if (sumCount > 0) {
    fprintf (stderr, "All %d sums:\n", sumCount);
    for (i = 0; i < SUM_HASH_SIZE; i++) {
      if (sumHash[i].num <= 0)
        continue;
      for (j = 0; j < sumHash[i].num; j++) {
        fprintf (stderr, "(%d %d) index=%d key=%d count=%d valid=%d value=%.18g sum: ", i, j, sumHash[i].index[j], sumHash[i].key[j], sumList[sumHash[i].index[j]]->count, sumList[sumHash[i].index[j]]->valid, sumList[sumHash[i].index[j]]->value);
        expTermPrinting (stderr, sumList[sumHash[i].index[j]], 1);
        fprintf (stderr, "\n");
      }
    }
    fprintf (stderr, "\n");
  }
  if (productCount > 0) {
    fprintf (stderr, "All %d products:\n", productCount);
    for (i = 0; i < PRODUCT_HASH_SIZE; i++) {
      if (productHash[i].num <= 0)
        continue;
      for (j = 0; j < productHash[i].num; j++) {
        fprintf (stderr, "(%d %d) index=%d key=%d count=%d valid=%d value=%.18g product: ", i, j, productHash[i].index[j], productHash[i].key[j], productList[productHash[i].index[j]]->count, productList[productHash[i].index[j]]->valid, productList[productHash[i].index[j]]->value);
        expTermPrinting (stderr, productList[productHash[i].index[j]], 1);
        fprintf (stderr, "\n");
      }
    }
    fprintf (stderr, "\n");
  }

  if (functionCallCount > 0) {
    fprintf (stderr, "All %d function calls:\n", functionCallCount);
    for (i = 0; i < FUNCTIONCALL_HASH_SIZE; i++) {
      if (functionCallHash[i].num <= 0)
        continue;
      for (j = 0; j < functionCallHash[i].num; j++) {
        fprintf (stderr, "(%d %d) index=%d key=%d count=%d valid=%d value=%.18g functionCall: ", i, j, functionCallHash[i].index[j], functionCallHash[i].key[j], functionCallList[functionCallHash[i].index[j]]->count, functionCallList[functionCallHash[i].index[j]]->valid, functionCallList[functionCallHash[i].index[j]]->value);
        expTermPrinting (stderr, functionCallList[functionCallHash[i].index[j]], 1);
        fprintf (stderr, "\n");
      }
    }
    fprintf (stderr, "\n");
  }
  fprintf (stderr, "---\n");
}

/* Mark all polynomials as held. */
void holdAllPolys ()
{
  int i, j;

  if (polynomialDebugLevel >= 8)
    fprintf (stderr, "holdAllPolys called for %d constants, %d variables, %d sums, %d products and %d functions\n",
	     constantCount, variableCount, sumCount, productCount, functionCallCount);
  holdAllPolysCount++;

#ifdef _OPENMP
#pragma omp sections
#endif
  {

#ifdef _OPENMP
#pragma omp section
#endif
    {
      if (constantCount > 0) {
        //      fprintf (stderr, "%d constants\n", constantCount);
        for (i = 0; i < CONSTANT_HASH_SIZE; i++) {
          if (constantHash[i].num <= 0)
            continue;
          for (j = 0; j < constantHash[i].num; j++)
            constantList[constantHash[i].index[j]]->count++;
        }
      }
    }

#ifdef _OPENMP
#pragma omp section
#endif
    {
      if (variableCount > 0) {
        //      fprintf (stderr, "%d variables:\n", variableCount);
        for (i = 0; i < VARIABLE_HASH_SIZE; i++) {
          if (variableHash[i].num <= 0)
            continue;
          for (j = 0; j < variableHash[i].num; j++) {
            //      fprintf (stderr, "%s\n", variableList[variableHash[i].index[j]]->e.v->vName);
            variableList[variableHash[i].index[j]]->count++;
          }
        }
      }
    }

#ifdef _OPENMP
#pragma omp section
#endif
    {
      if (sumCount > 0) {
        //      fprintf (stderr, "%d sums\n", sumCount);
        for (i = 0; i < SUM_HASH_SIZE; i++) {
          if (sumHash[i].num <= 0)
            continue;
          for (j = 0; j < sumHash[i].num; j++) {
            if (sumList[sumHash[i].index[j]]->id == polynomialLostNodeId)
              fprintf (stderr, "holdAllPolys for sums sees id %d and is bumping count\n", polynomialLostNodeId);
            sumList[sumHash[i].index[j]]->count++;
          }
        }
      }
    }

#ifdef _OPENMP
#pragma omp section
#endif
    {
      if (productCount > 0) {
        //      fprintf (stderr, "%d products\n", productCount);
        for (i = 0; i < PRODUCT_HASH_SIZE; i++) {
          if (productHash[i].num <= 0)
            continue;
          for (j = 0; j < productHash[i].num; j++) {
            if (productList[productHash[i].index[j]]->id == polynomialLostNodeId)
              fprintf (stderr, "holdAllPolys for products sees id %d and is bumping count\n", polynomialLostNodeId);
            productList[productHash[i].index[j]]->count++;
          }
        }
      }
    }

#ifdef _OPENMP
#pragma omp section
#endif
    {
      if (functionCallCount > 0) {
        //      fprintf (stderr, "%d function calls\n", functionCallCount);
        for (i = 0; i < FUNCTIONCALL_HASH_SIZE; i++) {
          if (functionCallHash[i].num <= 0)
            continue;
          for (j = 0; j < functionCallHash[i].num; j++)
            functionCallList[functionCallHash[i].index[j]]->count++;
        }
      }
    }
  }
//  printAllPolynomials();
//  fprintf (stderr, "---\n");
}

/* Flag all components of the provided polynomial for preservation. */
void doKeepPoly (Polynomial * p)
{
  int i;

  /* If we hit a term with the flag already set, we can stop, because
   * all of it's subpolys will have it set as well. This flag is only
   * cleared en-masse so there will never be a freeing of a subpoly
   * without those it contributes to being freed as well. */
  if (p->valid & VALID_KEEP_FLAG)
    return;

  if (p->id == polynomialLostNodeId)
    fprintf (stderr, "flagValids sees id %d and is flagging with %d\n", polynomialLostNodeId, VALID_KEEP_FLAG);

  p->valid |= VALID_KEEP_FLAG;
  switch (p->eType) {
  case T_OFFLINE:
    importPoly(p);
    // Notice that there's no break here...

  case T_CONSTANT:
  case T_VARIABLE:
  case T_EXTERNAL:
    break;
  case T_SUM:
    importTermList (p);
    for (i = 0; i < p->e.s->num; i++) {
      doKeepPoly (p->e.s->sum[i]);
    }
    exportTermList (p, FALSE);
    break;
  case T_PRODUCT:
  case T_FUNCTIONCALL:
    for (i = 0; i < p->e.p->num; i++) {
      doKeepPoly (p->e.p->product[i]);
    }
    break;
  case T_FREED:
    FATAL ("[FREED eType=%d id=%d index=%d key=%d count=%d valid=%d]",
	     (int) p->value, p->id, p->index, p->key, p->count, p->valid);
  default:
    FATAL ("In doKeepPoly, unknown expression type %d", p->eType);
  }
  return;
}

void keepPoly (Polynomial * p)
{
  if (p->valid & VALID_KEEP_FLAG)
    return;
  keepPolyCount++;
  if (polynomialDebugLevel >= 10)
    fprintf (stderr, "Into keepPoly\n");
  doKeepPoly (p);
  if (polynomialDebugLevel >= 10)
    fprintf (stderr, "Out of keepPoly\n");
  return;
}

/* Bump the hold count for this and all component subpolynomials. */
void doHoldPoly (Polynomial * p)
{
  int i;

  if (p->valid & VALID_EVAL_FLAG)
    return;
  p->valid |= VALID_EVAL_FLAG;

  if (p->id == polynomialLostNodeId)
    fprintf (stderr, "holdPoly sees id %d and is bumping hold count from %d\n", polynomialLostNodeId, p->count);

  switch (p->eType) {
  case T_OFFLINE:
    importPoly(p);
    // Notice that there's no break here...

  case T_CONSTANT:
  case T_VARIABLE:
  case T_EXTERNAL:
    p->count++;
    break;
  case T_SUM:
    importTermList (p);
    p->count++;
    for (i = 0; i < p->e.s->num; i++) {
      doHoldPoly (p->e.s->sum[i]);
    }
    exportTermList (p, FALSE);
    break;
  case T_PRODUCT:
  case T_FUNCTIONCALL:
    p->count++;
    for (i = 0; i < p->e.p->num; i++) {
      doHoldPoly (p->e.p->product[i]);
    }
    break;
  default:
    FATAL ("In holdPoly, unknown expression type %d", p->eType);
  }
  return;
}

void holdPoly (Polynomial * p)
{
  holdPolyCount++;
  if (polynomialDebugLevel >= 10)
    fprintf (stderr, "Into holdPoly\n");
  clearValidEvalFlag ();
  doHoldPoly (p);
  if (polynomialDebugLevel >= 10)
    fprintf (stderr, "Out of holdPoly\n");
  return;
}

/* Decrement the hold count for this and all component subpolynomials. */
void doUnHoldPoly (Polynomial * p)
{
  int i;

  if (p->valid & VALID_EVAL_FLAG)
    return;
  p->valid |= VALID_EVAL_FLAG;

  if (p->id == polynomialLostNodeId)
    fprintf (stderr, "UnHoldPoly sees id %d and is decrementing hold count from %d\n",
	     polynomialLostNodeId, p->count);

  switch (p->eType) {
  case T_OFFLINE:
    importPoly(p);
    // Notice that there's no break here...

  case T_CONSTANT:
  case T_VARIABLE:
  case T_EXTERNAL:
    p->count--;
    break;
  case T_SUM:
    p->count--;
    importTermList (p);
    for (i = 0; i < p->e.s->num; i++) {
      doUnHoldPoly (p->e.s->sum[i]);
    }
    exportTermList (p, FALSE);
    break;
  case T_PRODUCT:
  case T_FUNCTIONCALL:
    p->count--;
    for (i = 0; i < p->e.p->num; i++) {
      doUnHoldPoly (p->e.p->product[i]);
    }
    break;
  default:
    FATAL ("In holdPoly, unknown expression type %d", p->eType);
  }
  return;
}

void unHoldPoly (Polynomial * p)
{
  unHoldPolyCount++;
  clearValidEvalFlag ();
  doUnHoldPoly (p);
  return;
}

void doFreePolys (unsigned short keepMask)
{
  int i, j, k;
  Polynomial **newConstantList,
#ifdef VARIABLE_C_MESSES_UP_COMPILED_DL_ARGUMENTS
    **newVariableList,
#endif
    **newSumList,
    **newProductList;

  if (polynomialDebugLevel >= 5)
    fprintf (stderr, "Starting doFreePolys with keepMask of %uh\n", keepMask);
  if (polynomialDebugLevel >= 70) {
    fprintf (stderr, "...with the following:\n");
    printAllPolynomials ();
  }
#ifdef _OPENMP
#pragma omp sections
#endif
  {
#ifdef _OPENMP
#pragma omp section
#endif
    {
      /* First create a new polynomial-specific list to transcribe entries
       * we're keeping, and then go thru the old polynomial-specific list flagging
       * entries we're not keeping and replacing pointers with new indexes
       * for the ones we are keeping. */

      MALCHOKE(newConstantList, sizeof (Polynomial *) * (constantListLength),Polynomial **);
      k = 0;
      for (i = 0; i < constantCount; i++) {
	if (constantList[i]->eType == T_OFFLINE) importPoly (constantList[i]);
        if (constantList[i]->id == polynomialLostNodeId)
          fprintf (stderr, "doFreePolys sees id %d with valid %d and count %d during pass with mask %d\n",
		   polynomialLostNodeId, constantList[i]->valid, constantList[i]->count, keepMask);
        if ((constantList[i]->count > 0)
            || (constantList[i]->valid & keepMask)) {
          newConstantList[k] = constantList[i];
          constantList[i] = newConstantList[k];
          //      fprintf(stderr, "Index %d is now %d\n", newConstantList[k]->index, k);
          newConstantList[k]->index = k;
          k++;
        } else {
#ifndef FREEDEBUG
          free (constantList[i]);
#else
          // These are for debugging mis-freed pointers
          constantList[i]->value = constantList[i]->eType;
          constantList[i]->eType = T_FREED;
#endif
          constantList[i] = NULL;
        }
      }
      if (polynomialDebugLevel >= 5)
        fprintf (stderr, "%d/%d constant(s) ", k, constantCount);
      constantCount = k;

      /* Go thru the hash collapsing entries and freeing the corresponding
       * polynomials and terms. Zero polynomial-specific list entries so they
       * can be collapsed next. */
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
      for (j = 0; j < CONSTANT_HASH_SIZE; j++) {
        if (constantHash[j].num > 0) {
          k = 0;
          for (i = 0; i < constantHash[j].num; i++) {
            if (constantList[constantHash[j].index[i]] != NULL) {
              /* It's a keeper, slide it down and bump the count */
              //        fprintf(stderr, "Hash keeper index %d is now %d\n", constantHash[j].index[i],
              //              constantList[constantHash[j].index[i]]->index);
              constantHash[j].index[k] = constantList[constantHash[j].index[i]]->index;
              constantHash[j].key[k] = constantHash[j].key[i];
              k++;
            }
          }
          constantHash[j].num = k;
        }
      }
      free (constantList);
      constantList = newConstantList;
    }

#ifdef VARIABLE_C_MESSES_UP_COMPILED_DL_ARGUMENTS
#ifdef _OPENMP
#pragma omp section
#endif
    {
      MALCHOKE(newVariableList, sizeof (Polynomial *) * (variableListLength),Polynomial **);
      k = 0;
      for (i = 0; i < variableCount; i++) {
	if (variableList[i]->eType == T_OFFLINE) importPoly (variableList[i]);
        if (variableList[i]->id == polynomialLostNodeId)
          fprintf (stderr, "doFreePolys sees id %d with valid %d and count %d during pass with mask %d\n",
		   polynomialLostNodeId, variableList[i]->valid, variableList[i]->count, keepMask);
        if ((variableList[i]->count > 0)
            || (variableList[i]->valid & keepMask)) {
          newVariableList[k] = variableList[i];
          variableList[i] = newVariableList[k];
          //      fprintf(stderr, "Index %d is now %d\n", newVariableList[k]->index, k);
          newVariableList[k]->index = k;
          k++;
        } else {
#ifndef FREEDEBUG
          free (variableList[i]);
#ifdef USE_GMP
	  mpf_clear (variableList[i]->mpfValue);
#endif
#else
          // These are for debugging mis-freed pointers
          variableList[i]->value = variableList[i]->eType;
          variableList[i]->eType = T_FREED;
#endif
          variableList[i] = NULL;
        }
      }
      if (polynomialDebugLevel >= 5)
        fprintf (stderr, "%d/%d variable(s) ", k, variableCount);
      variableCount = k;

      /* Go thru the hash collapsing entries and freeing the corresponding
       * polynomials and terms. Zero polynomial-specific list entries so they
       * can be collapsed next. */
      for (j = 0; j < VARIABLE_HASH_SIZE; j++) {
        if (variableHash[j].num > 0) {
          k = 0;
          for (i = 0; i < variableHash[j].num; i++) {
            if (variableList[variableHash[j].index[i]] != NULL) {
              /* It's a keeper, slide it down and bump the count */
              //        fprintf(stderr, "Hash keeper index %d is now %d\n", variableHash[j].index[i],
              //              variableList[variableHash[j].index[i]]->index);
              variableHash[j].index[k] = variableList[variableHash[j].index[i]]->index;
              variableHash[j].key[k] = variableHash[j].key[i];
              k++;
            }
          }
          variableHash[j].num = k;
        }
      }
      free (variableList);
      variableList = newVariableList;
    }
#endif

#ifdef _OPENMP
#pragma omp section
#endif
    {
      MALCHOKE(newSumList, sizeof (Polynomial *) * (sumListLength),Polynomial **);
      k = 0;
      for (i = 0; i < sumCount; i++) {
	if (sumList[i]->eType == T_OFFLINE) importPoly (sumList[i]);
        if (sumList[i]->id == polynomialLostNodeId)
          fprintf (stderr, "doFreePolys sees id %d with valid %d and count %d during pass with mask %d\n",
		   polynomialLostNodeId, sumList[i]->valid, sumList[i]->count, keepMask);
        if ((sumList[i]->count > 0) || (sumList[i]->valid & keepMask)) {
          newSumList[k] = sumList[i];
          sumList[i] = newSumList[k];
          //      fprintf(stderr, "Index %d is now %d\n", newSumList[k]->index, k);
          newSumList[k]->index = k;
          k++;
        } else {
          sumFreedCount++;
#ifndef USE_SSD
          free (sumList[i]->e.s->sum);
          free (sumList[i]->e.s->factor);
#else
	  deportTermList (sumList[i]);
#endif
          free (sumList[i]->e.s);
#ifndef FREEDEBUG
#ifdef USE_GMP
	  mpf_clear (sumList[i]->mpfValue);
#endif
          free (sumList[i]);
#else
          // These are for debugging mis-freed pointers
          sumList[i]->value = sumList[i]->eType;
          sumList[i]->eType = T_FREED;
#endif
          sumList[i] = NULL;
        }
      }
      if (polynomialDebugLevel >= 5)
        fprintf (stderr, "%d/%d sum(s) ", k, sumCount);
      sumCount = k;

      /* Go thru the hash collapsing entries and freeing the corresponding
       * polynomials and terms. Zero polynomial-specific list entries so they
       * can be collapsed next. */
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
      for (j = 0; j < SUM_HASH_SIZE; j++) {
        if (sumHash[j].num > 0) {
          k = 0;
          for (i = 0; i < sumHash[j].num; i++) {
            if (sumList[sumHash[j].index[i]] != NULL) {
              /* It's a keeper, slide it down and bump the count */
              //        fprintf(stderr, "Hash keeper index %d is now %d\n", sumHash[j].index[i],
              //              sumList[sumHash[j].index[i]]->index);
              sumHash[j].index[k] = sumList[sumHash[j].index[i]]->index;
              sumHash[j].key[k] = sumHash[j].key[i];
              k++;
            }
          }
          sumHash[j].num = k;
        }
      }
      free (sumList);
      sumList = newSumList;
      if ((sumCount + (2 * SUM_LIST_INITIAL)) < sumListLength) {
        //      fprintf (stderr, "Reducing sumListLength from %d to %d\n",
        //               sumListLength, sumCount + SUM_LIST_INITIAL);
        sumListLength = sumCount + SUM_LIST_INITIAL;
        sumList = (Polynomial **) realloc (sumList, sizeof (Polynomial *) * (sumListLength));
      }
    }

#ifdef _OPENMP
#pragma omp section
#endif
    {
      MALCHOKE(newProductList, sizeof (Polynomial *) * (productListLength),Polynomial **);
      k = 0;
      for (i = 0; i < productCount; i++) {
	if (productList[i]->eType == T_OFFLINE) importPoly (productList[i]);
        if (productList[i]->id == polynomialLostNodeId)
          fprintf (stderr, "doFreePolys sees id %d with valid %d and count %d during pass with mask %d\n",
		   polynomialLostNodeId, productList[i]->valid, productList[i]->count, keepMask);
        if ((productList[i]->count > 0) || (productList[i]->valid & keepMask)) {
          newProductList[k] = productList[i];
          productList[i] = newProductList[k];
          //      fprintf(stderr, "Index %d is now %d\n", newProductList[k]->index, k);
          newProductList[k]->index = k;
          k++;
        } else {
          productFreedCount++;
          free (productList[i]->e.p->product);
          free (productList[i]->e.p->exponent);
          free (productList[i]->e.p);
#ifndef FREEDEBUG
#ifdef USE_GMP
	  mpf_clear (productList[i]->mpfValue);
#endif
          free (productList[i]);
#else
          // These are for debugging mis-freed pointers
          productList[i]->value = productList[i]->eType;
          productList[i]->eType = T_FREED;
#endif
          productList[i] = NULL;
        }
      }
      if (polynomialDebugLevel >= 5)
        fprintf (stderr, "%d/%d product(s) ", k, productCount);
      productCount = k;

      /* Go thru the hash collapsing entries and freeing the corresponding
       * polynomials and terms. Zero polynomial-specific list entries so they
       * can be collapsed next. */
#ifdef _OPENMP
#pragma omp parallel for private(i,j,k)
#endif
      for (j = 0; j < PRODUCT_HASH_SIZE; j++) {
        if (productHash[j].num > 0) {
          k = 0;
          for (i = 0; i < productHash[j].num; i++) {
            if (productList[productHash[j].index[i]] != NULL) {
              /* It's a keeper, slide it down and bump the count */
              //        fprintf(stderr, "Hash keeper index %d is now %d\n", productHash[j].index[i],
              //              productList[productHash[j].index[i]]->index);
              productHash[j].index[k] = productList[productHash[j].index[i]]->index;
              productHash[j].key[k] = productHash[j].key[i];
              k++;
            }
          }
          productHash[j].num = k;
        }
      }
      free (productList);
      productList = newProductList;
      if ((productCount + (2 * PRODUCT_LIST_INITIAL)) < productListLength) {
        //      fprintf (stderr, "Reducing productListLength from %d to %d\n",
        //               productListLength, productCount + PRODUCT_LIST_INITIAL);
        productListLength = productCount + PRODUCT_LIST_INITIAL;
        productList = (Polynomial **) realloc (productList, sizeof (Polynomial *) * (productListLength));
      }
    }
  }
  /* Reset building statistics. */
  constantHashHits = variableHashHits = functionHashHits = 0;
  sumReleaseableCount = sumNotReleaseableCount = sumReturnConstantCount = sumReturnSameConstantCount = sumReturn1TermCount = sumHashHits = sumNewCount = sumListNewCount = sumListReplacementCount = 0;
  productReleaseableCount = productNotReleaseableCount = productReturn0Count =
      productReturnConstantCount = productReturn1stTermCount = productReturn1TermSumCount = productHashHits = productHashHitIsSumCount = productReturnNormalCount = productNon1FactorIsSumCount = productListNewCount = productListReplacementCount = 0;
  totalSPLLengths = totalSPLCalls = lowSPLCount = highSPLCount = 0;

  if (polynomialDebugLevel >= 5)
    fprintf (stderr, "\n");
  if (polynomialDebugLevel >= 70) {
    fprintf (stderr, "...to the following:\n");
    printAllPolynomials ();
  }
  return;
}

/* Free all polynomials that aren't held or kept. Be a bit circumspect
   about freeing when we don't really need to, so that caching is
   improved, and even more, so we don't constantly go thru the work
   of freeing without much to show for it. */
void freePolys ()
{
  /*
#ifdef USE_SSD
  fprintf (stderr, "SSD sum terms imported/+written-dropped exported: %d/+%d-%d, peak in-memory: %d\n",
	   importedTerms, exportedWrittenTerms, exportedDroppedTerms, peakInMemoryTerms);
  statSSD ();
#endif
  */
  if (sumNotReleaseableCount + productNotReleaseableCount < 1024 * 512) {
    freePolysAttemptCount++;
    return;
  }
  freePolysCount++;
  if (polynomialDebugLevel >= 5)
    fprintf (stderr, "Free w/mask of VALID_KEEP_FLAG (%d) preserved...", VALID_KEEP_FLAG);
  doFreePolys (VALID_KEEP_FLAG);
  if (polynomialDebugLevel >= 5)
    fprintf (stderr, "\n");
  return;
}

/* Free all polynomials that aren't held. */
void freeKeptPolys ()
{
  freeKeptPolysCount++;
  keepPolyCount = 0;
  if (polynomialDebugLevel >= 5)
    fprintf (stderr, "Free w/mask of %d preserved...", 0);
  doFreePolys (0);
  if (polynomialDebugLevel >= 5)
    fprintf (stderr, "\n");
  return;
}

/*
  We need to track the pieces of the chunkTicket whether the term list is in 
  or out of memory. The new sumPoly.iMTLIndex will tell us the status of the
  term list: -1 means dynamically allocated in-memory (inital state), -2 means 
  not in-memory (exported), anything else means statically allocated in-memory.

  When the term list is out of memory, we need to keep the chunkTicket pieces
  in order to get it back. Since we're not using the sumPoly sum and factor 
  list pointers, they'll work for that, so the sum list pointer will be the
  chunkOffset and the factor list pointer will be the doublePairCount.

  When the term list is in memory, we need to keep the chunkTicket somewhere 
  else because sum and factor are used. Fortunately, being in-memory and having
  caching means that we can keep the chunkTicket in the cache entry.

*/  

void importTermList (Polynomial * p)
{
#ifdef USE_SSD
  int i;
  struct sumPoly *sP;

  if (p->eType != T_SUM)
    return;

  sP = (struct sumPoly *) p->e.s;

  if (sP->num < MIN_USE_SSD_DPS)
    return;

  // You can only import what you've exported before
  if (sP->iMTLIndex != -2)
    return;

  //  printf ("Import p%d of %d sum %lu and factor %lu\n", 
  //	  p->id, sP->num, (unsigned long) sP->sum, (unsigned long) sP->factor);

  importedTerms++;
  if ((importedTerms-(exportedWrittenTerms+exportedDroppedTerms)) > peakInMemoryTerms)
    peakInMemoryTerms = importedTerms-(exportedWrittenTerms+exportedDroppedTerms);

  // Actually bring it back in...

  for (i=0; i<=IMTL_COUNT; i++)
    if (iMTL[i].lastNodeId == p->id) {
      // Found it still hanging about!
      iMTL[i].cT.chunkOffset = (unsigned long) sP->sum;
      iMTL[i].cT.doublePairCount = (unsigned long) sP->factor;
      sP->iMTLIndex = i;
      sP->sum = (Polynomial **) &iMTL[i].buffer[0];
      sP->factor = (double *) &iMTL[i].buffer[sP->num];
      //      printf ("Import hit cache for p%d iMTLIndex %d of %d sum %lu and factor %lu\n",
      //	      p->id, sP->iMTLIndex, sP->num, (unsigned long) sP->sum, (unsigned long) sP->factor);
      return;
    }

  // Not still hanging about, so let's go get it...

  for (i=0; i<=IMTL_COUNT; i++)
    if (iMTL[i].cT.doublePairCount == 0)
      break;

  if (i >= IMTL_COUNT)
    FATAL ("Exhausted in-memory term lists for sum term lists!");

  iMTL[i].cT.chunkOffset = (unsigned long) sP->sum;
  iMTL[i].cT.doublePairCount = (unsigned long) sP->factor;
  iMTL[i].lastNodeId = p->id;
  getSSD (&iMTL[i].cT, iMTL[i].buffer);
  sP->iMTLIndex = i;
  sP->sum = (Polynomial **) &iMTL[i].buffer[0];
  sP->factor = (double *) &iMTL[i].buffer[sP->num];

  //  printf ("Imported p%d iMTLIndex %d of %d sum %lu and factor %lu\n",
  //	  p->id, sP->iMTLIndex, sP->num, (unsigned long) sP->sum, (unsigned long) sP->factor);

#endif
}

/**

  Write a list of sum terms out to SSD storage for later import.

*/
void exportTermList (Polynomial * p, int writeFlag)
{
#ifdef USE_SSD
  struct sumPoly *sP;

  // Only handle sum polynomials
  if (p->eType != T_SUM)
    return;

  sP = (struct sumPoly *) p->e.s;

  // Only handle largish sum polynomials
  if (sP->num < MIN_USE_SSD_DPS)
    return;

  // Can't export something already exported.
  if (sP->iMTLIndex == -2)
    return;

  //  printf ("Export (%d) p%d w/iMTLIndex %d of %d sum %lu and factor %lu\n", 
  //	  writeFlag, p->id, sP->iMTLIndex, sP->num, (unsigned long) sP->sum, (unsigned long) sP->factor);

  if (writeFlag) {
    // Actually write it out...

    struct chunkTicket *cT;
    if (sP->iMTLIndex == -1) {
      // ...and it's a new term list in allocated (discontiguous) memory

      // Only do exports when we're below some available physical memory threshold.
      if (memoryLow != TRUE)
	return;

      if (sSDDebt > 0)
	return; // ...no can do.

      double buffer[MAX_SSD_BUFFER];

      if (sP->num > (MAX_SSD_BUFFER / 2)) {
	WARNING ("MAX_SSD_BUFFER of %d is too small for request for %d double pairs",
		 MAX_SSD_BUFFER, sP->num);
	return; // ...leave it in memory.
      }

      // Avoid this painful mistake by changing sSDHandler's get and putSSD to have two separate pointers
      memcpy (&buffer[0], sP->sum, 8 * sP->num);
      memcpy (&buffer[sP->num], sP->factor, 8 * sP->num);
      if ((cT = putSSD (buffer, sP->num)) == NULL) {
	// We're out of business, put up the "Closed" sign until we get enough back.
	sSDDebt = sP->num << 2;
	WARNING ("In SSD debt for %ludps, reverting to memory allocation", sSDDebt);
	return;
      }
      free (sP->sum);
      free (sP->factor);
      sP->sum = (Polynomial **) cT->chunkOffset;
      sP->factor = (double *) cT->doublePairCount;
    } else {
      // ...and it's an old term list in a static buffer, so just rewrite it
      flushSSD (iMTL[sP->iMTLIndex].buffer, iMTL[sP->iMTLIndex].cT.chunkOffset, 
		iMTL[sP->iMTLIndex].cT.doublePairCount);
      sP->sum = (Polynomial **) iMTL[sP->iMTLIndex].cT.chunkOffset;
      sP->factor = (double *) iMTL[sP->iMTLIndex].cT.doublePairCount;
      iMTL[sP->iMTLIndex].cT.doublePairCount = 0;
    }
    //    printf ("Exported sum %lu and factor %lu\n", (unsigned long) sP->sum, (unsigned long) sP->factor);
    exportedWrittenTerms++;
  } else {
    if (sP->iMTLIndex == -1)
      // ...and it's a term list in allocated (discontiguous) memory, so this is a no-op
      return;
    // ...and it's an old term list in a static buffer, so just rewrite it
    sP->sum = (Polynomial **) iMTL[sP->iMTLIndex].cT.chunkOffset;
    sP->factor = (double *) iMTL[sP->iMTLIndex].cT.doublePairCount;
    iMTL[sP->iMTLIndex].cT.doublePairCount = 0;
    exportedDroppedTerms++;
  }
  sP->iMTLIndex = -2;

#endif
}

/*

  The other calls (import and export) are optional, i.e. they can choose not
  to do anything based upon some criteria. deportTermList is not option -- it
  will be called *instead* of an explicit free for the term lists, so it should
  always either do the explicit free or ditch the SSD term list.

*/

void deportTermList (Polynomial * p)
{
#ifdef USE_SSD
  struct sumPoly *sP;
  int i;

  if (p->eType != T_SUM)
    return;

  sP = (struct sumPoly *) p->e.s;

  //  printf ("Deport p%d w/iMTLIndex %d of %d sum %lu and factor %lu\n", 
  //	  p->id, sP->iMTLIndex, sP->num, (unsigned long) sP->sum, (unsigned long) sP->factor);

  if (sP->iMTLIndex != -1) {
    for (i=0; i<=IMTL_COUNT; i++)
      if (iMTL[i].lastNodeId == p->id)
	iMTL[i].lastNodeId = 0;
    struct chunkTicket *cT;
    MALCHOKE(cT, sizeof(struct chunkTicket),struct chunkTicket *);
    cT->chunkOffset = (unsigned long) sP->sum;
    cT->doublePairCount = (unsigned long) sP->factor;
    // Pay-off some (eventually all?) of our debt
    if (sSDDebt > 0) {
      sSDDebt -= cT->doublePairCount;
      if (sSDDebt <= 0)
	WARNING ("Paid-off SSD debt!");
    }
    freeSSD (cT);
  } else {
    free (sP->sum);
    free (sP->factor);
  }

#endif
}

/* Serialize and compress a polynomial and all subpolys for transfer
   Addresses can be left as-is since we'll create a translation table when we
   reload, and compression will eliminate any temporary inefficiency of size. */
Polynomial *exportPoly (Polynomial * p)
{
  //  p->oldEType = p->eType;
  return (p);
}

/* Decompress and deserialize a polynomial and all subpolys. Use an address
   translation table to maintain internal consistency. */
Polynomial *importPoly (Polynomial * p)
{
  //  p->eType = p->oldEType;
  return (p);
}

void releaseExternalPoly (Polynomial *rp)
{
  if (rp->eType == T_OFFLINE) importPoly (rp);
  if (rp->eType != T_EXTERNAL)
    FATAL ("releaseExternalPoly called with polynomial eType of %d", rp->eType);
#ifdef POLYUSE_DL
  int i;
  for (i=0; i<32; i++)
    if (rp->e.e->polynomialFunctionHandle[i] != NULL)
      dlclose (rp->e.e->polynomialFunctionHandle[i]);
    else
      break;
#endif
  swLogProgress(3 /* DETAIL + 1 */, 0, "Released polynomial DL %s", rp->e.e->polynomialFunctionName);
  rp->e.e->fileOK = FALSE;
  rp->e.e->entryOK = FALSE;
  return;
}

Polynomial *restoreExternalPoly (char *functionName)
{
  Polynomial *rp;
  struct externalPoly *eP;
  int i;

  // See if it already exists
  for (i = 0; i<externalCount; i++)
    if (strncmp (functionName, externalList[i]->e.e->polynomialFunctionName, MAX_PFN_LEN) == 0)
      return externalList[i];

  // New polynomial to create!

  MALCHOKE(rp, sizeof (Polynomial),Polynomial *);
  MALCHOKE(eP, sizeof (struct externalPoly),struct externalPoly *);
  rp->e.e = eP;
  strcpy (eP->polynomialFunctionName, functionName);
#ifdef POLYUSE_DL
  if (!loadPolyDL(rp)) {
    free (eP);
    free (rp);
    return NULL;
  }
#endif
  rp->index = externalCount;
  rp->id = nodeId;
  rp->eType = T_EXTERNAL;
  rp->valid = VALID_NOTDISC_FLAG;
  rp->count = 0;
  externalList[externalCount] = rp;
  externalCount++;
  if (externalCount >= externalListLength) {
    externalListLength += EXTERNAL_LIST_INCREASE;
    REALCHOKE (externalList, externalListLength * sizeof (Polynomial *), void *);
  }
  nodeId++;

  // These polynomials aren't hashed or listed
  return (rp);
}

#ifdef POLYUSE_DL
int loadPolyDL (Polynomial * p)
{
  int i = 0;
  char polynomialFileName[MAX_PFN_LEN+1];

  sprintf (polynomialFileName, "./%s.so", p->e.e->polynomialFunctionName);
  if ((p->e.e->polynomialFunctionHandle[0] = 
       dlopen (polynomialFileName, RTLD_LAZY|RTLD_LOCAL)) != NULL) {

    // Loaded! Do any supporting 1K clump DLs.
    for (i=0; i<=32; i++) {
      sprintf (polynomialFileName, "./%s_%dK.so", p->e.e->polynomialFunctionName, i);
      if ((p->e.e->polynomialFunctionHandle[i+1] = dlopen (polynomialFileName, RTLD_LAZY|RTLD_LOCAL)) == NULL)
	break;
    }

    // Now get the first entry point symbol
    if ((p->e.e->polynomialFunctionRoutine = 
	 dlsym (p->e.e->polynomialFunctionHandle[0], p->e.e->polynomialFunctionName)) != NULL) {
      // Found it!
      p->e.e->fileOK = TRUE;
      p->e.e->entryOK = TRUE;
      swLogProgress(3, 0, "Using %d DL(s) for %s", i+1, p->e.e->polynomialFunctionName);
    } else {
      WARNING ("dlsym() error [%s] for polynomial %s", dlerror(), p->e.e->polynomialFunctionName);
      for (i=0; i<32; i++)
	if (p->e.e->polynomialFunctionHandle[i] != NULL)
	  dlclose (p->e.e->polynomialFunctionHandle[i]);
	else
	  break;
      return FALSE;
    }
  } else {
    swLogProgress(3 /* DETAIL + 1 */, 0,
		  "dlopen() error [%s] for polynomial %s",
		  dlerror(), p->e.e->polynomialFunctionName);
    return FALSE;
  }
  return TRUE;
}
#endif

/**

  Largest single source file to produce.

  8192*496 is carefully tested to compile in just under 16Gb in our environment.
  Use 8192*128 for a 4Gb compile size, but be aware that large polynomials (over
  1,000 source files) might not compile due to command line size limitations.

  Dagnabit, the only problem is that while compilation seemed to only take that
  much space, I finally got the error:

cc1: out of memory allocating 18446744058046806904 bytes after a total of 289832960 bytes

  8192*384 works, though. And crazy as it may seem, renaming the source files to something
  short during compilation and linkage works as well.

*/
#define MAXSRCSIZE (8192*384)
/**

  Generate C surce code that, when compiled and executed, will perform the identical
  operations as evaluatePoly. Shockingly, much of the work here involved finding ways
  around limitations in the number of characters in a line of code, the size of source
  files, the number of routines in a single DL, etc.

*/
void codePoly (Polynomial * p, struct polyList *l, char *name)
{
  char srcFileName[PATH_MAX], srcCalledFileName[PATH_MAX], includeFileName[PATH_MAX];
  FILE *srcFile, *srcCalledFile = NULL, *includeFile;
  int i, j, srcSize = 0, fileCount = 0, totalSourceSize = 0, totalInternalSize = 0;
  int sumsUsed = 0, productsUsed = 0, functionCallsUsed = 0;
  Polynomial *result;
  struct sumPoly *sP;
  struct productPoly *pP;

  swPushPhase ('k', "encode poly");
  result = p;
  fprintf (stdout, "\nGenerating DL for %s...\n", name);

  sprintf (srcFileName, "%s.c", name);
  sprintf (includeFileName, "%s.h", name);
  if ((srcFile = fopen (srcFileName, "w")) == NULL) {
    perror ("Cannot open polynomial source file\n");
    exit (EXIT_FAILURE);
  }

  totalSourceSize += fprintf (srcFile, "#include <math.h>\n#include <stdlib.h>\n\n");
#ifdef USE_GSL
  totalSourceSize += fprintf (srcFile, "#include <gsl/gsl_randist.h>\n\n#include <gsl/gsl_cdf.h>\n\n#include <gsl/gsl_sf_gamma.h>\n\n");
#endif
  totalSourceSize += fprintf (srcFile, "#include \"%s.h\"\n\n", name);

  /* We need to be *very* considerate of what is global and what is not when dealing with
     multiple generated dynamic libraries. Nothing should be global that is not
     uniquely named amid all other currently loaded DLs, otherwise we'll definitely not be
     thread-safe, and we might not be able to work at all.
  */

  // Here's the global part...

#ifdef POLYCODE_DL
  totalSourceSize += fprintf (srcFile, "#include <dlfcn.h>\n#include \"polynomial.h\"\n#include \"dists.h\"\n\n");
  totalSourceSize += fprintf (srcFile, "char *baseFunctionName = \"%s\";\nint baseFunctionArgs = %d, dLFunctionCount = DLFUNCTIONCOUNT;\n\n", 
			      name, variableCount);
  totalSourceSize += fprintf (srcFile, "double %s (int num, struct polynomial **variableList) {\n", name);
#endif

  // And now we're local... A static here keeps these off the stack, which would otherwise frequently trash it.

  totalSourceSize += fprintf (srcFile, "static double V[VARIABLESUSED], S[SUMSUSED], " "P[PRODUCTSUSED], F[FUNCTIONCALLSUSED];\n\n");

  for (i = 0; i < variableCount; i++)
    totalSourceSize += fprintf (srcFile, "\tV[%d] = variableList[%d]->value; // %s\n", i, i, variableList[i]->e.v->vName);

  // Start by writing source lines to the first-tier DL
  srcCalledFile = srcFile;

  for (j = 0; j <= l->listNext - 1; j++) {

    if (srcSize >= MAXSRCSIZE) {

      if (fileCount > 0) {
	srcSize += fprintf (srcCalledFile, "}\n");
	fclose (srcCalledFile);
	totalSourceSize += srcSize;
      }
      srcSize = 0;

      sprintf (srcCalledFileName, "%s_%04d.c", name, fileCount);
      if ((srcCalledFile = fopen (srcCalledFileName, "w")) == NULL) {
        perror ("Cannot open polynomial source called file\n");
        exit (EXIT_FAILURE);
      }
      srcSize += fprintf (srcCalledFile, "#include <math.h>\n#include <stdio.h>\n\n");
      srcSize += fprintf (srcCalledFile, "double %s_%04d (double V[], double S[], double P[], double F[]) {\n", name, fileCount);
      fileCount++;
    }

    totalInternalSize += sizeof (Polynomial) + sizeof (Polynomial *);
    p = l->pList[j];
    switch (p->eType) {
    case T_OFFLINE:
      importPoly(p);
      // Notice that there's no break here...

    case T_CONSTANT:
      break;

    case T_VARIABLE:
      totalInternalSize += sizeof (struct variablePoly);
      p->value = p->index;
      break;

    case T_SUM:
      srcSize += fprintf (srcCalledFile, "\t%s[%d] = ", eTypes[p->eType], sumsUsed);
      sP = p->e.s;
      totalInternalSize += sizeof (struct sumPoly) + sP->num * (sizeof (Polynomial *) + sizeof (double *));
      importTermList (p);
      for (i = 0; i < sP->num; i++) {
	if (sP->sum[i]->eType == T_OFFLINE) importPoly (sP->sum[i]);
        if (i != 0)
          srcSize += fprintf (srcCalledFile, "+");
        if (sP->sum[i]->eType == T_CONSTANT)
          srcSize += fprintf (srcCalledFile, "%.*g", DBL_DIG, sP->sum[i]->value /* still the constant */ );
        else {
          if (sP->factor[i] == 1)
            srcSize += fprintf (srcCalledFile, "%s[%lu]", eTypes[sP->sum[i]->eType], (unsigned long) sP->sum[i]->value /* actually <whatever>Used */ );
          else
            srcSize += fprintf (srcCalledFile, "%.*g*%s[%lu]", DBL_DIG, sP->factor[i], eTypes[sP->sum[i]->eType], (unsigned long) sP->sum[i]->value /* actually <whatever>Used */ );
        }
      }
      exportTermList (p, FALSE);
      p->value = sumsUsed++;
      break;

    case T_PRODUCT:
      srcSize += fprintf (srcCalledFile, "\t%s[%d] = ", eTypes[p->eType], productsUsed);
      pP = p->e.p;
      totalInternalSize += sizeof (struct productPoly) + pP->num * (sizeof (Polynomial *) + sizeof (int *));
      for (i = 0; i < pP->num; i++) {
	if (pP->product[i]->eType == T_OFFLINE) importPoly (pP->product[i]);
        if (i != 0)
          srcSize += fprintf (srcCalledFile, "*");
        if (pP->product[i]->eType == T_CONSTANT)
          srcSize += fprintf (srcCalledFile, "%.*g", DBL_DIG, pP->product[i]->value /* still the constant */ );
        else {
          switch (pP->exponent[i]) {
          case 1:
            srcSize += fprintf (srcCalledFile, "%s[%lu]", eTypes[pP->product[i]->eType], (unsigned long) pP->product[i]->value);
            break;
          case 2:
            srcSize += fprintf (srcCalledFile, "%s[%lu]*%s[%lu]", eTypes[pP->product[i]->eType], (unsigned long) pP->product[i]->value, eTypes[pP->product[i]->eType], (unsigned long) pP->product[i]->value);
            break;
          case 3:
            srcSize += fprintf (srcCalledFile, "%s[%lu]*%s[%lu]*%s[%lu]", eTypes[pP->product[i]->eType],
                (unsigned long) pP->product[i]->value, eTypes[pP->product[i]->eType], (unsigned long) pP->product[i]->value, eTypes[pP->product[i]->eType], (unsigned long) pP->product[i]->value);
            break;
          case 4:
            srcSize += fprintf (srcCalledFile, "%s[%lu]*%s[%lu]*%s[%lu]*%s[%lu]", eTypes[pP->product[i]->eType],
                (unsigned long) pP->product[i]->value, eTypes[pP->product[i]->eType],
                (unsigned long) pP->product[i]->value, eTypes[pP->product[i]->eType], (unsigned long) pP->product[i]->value, eTypes[pP->product[i]->eType], (unsigned long) pP->product[i]->value);
            break;
          default:
            srcSize += fprintf (srcCalledFile, "pow(%s[%lu],%d)", eTypes[pP->product[i]->eType], (unsigned long) pP->product[i]->value, pP->exponent[i]);
            break;
          }
        }
      }
      //      srcSize += fprintf (srcCalledFile, ";\n\tprintf (\"%%g\\n\", %s[%d])", eTypes[p->eType], productsUsed);
      p->value = productsUsed++;
      break;

    case T_FUNCTIONCALL:
      totalInternalSize += sizeof (struct functionPoly) + (p->e.f->num * sizeof (Polynomial *));

      /// @todo Not sure what this was about...
      //      if (strcmp (p->e.f->name, "gsl_ran_ugaussian_pdf") == 0)
      //	srcSize += fprintf (srcCalledFile, "\tF[%d] = ugaussian_pdf(", functionCallsUsed);
      //      else
	srcSize += fprintf (srcCalledFile, "\tF[%d] = %s(", functionCallsUsed, p->e.f->name);
      for (i = 0; i < p->e.f->num; i++) {
	if (p->e.f->para[i]->eType == T_OFFLINE) importPoly (p->e.f->para[i]);
        if (i != 0)
          srcSize += fprintf (srcCalledFile, ",");
        if (p->e.f->para[i]->eType == T_CONSTANT)
          srcSize += fprintf (srcCalledFile, "%.*g", DBL_DIG, p->e.f->para[i]->value /* still the constant */ );
        else
          srcSize += fprintf (srcCalledFile, "%s[%lu]", eTypes[p->e.f->para[i]->eType], (unsigned long) p->e.f->para[i]->value);
      }
      srcSize += fprintf (srcCalledFile, ")");
      p->value = functionCallsUsed++;
      break;

    default:
      FATAL ("In codePoly, unknown expression type: [%d]", p->eType);
      break;
    }
    srcSize += fprintf (srcCalledFile, ";\n");
  }

  if (fileCount > 0) {
    totalSourceSize += fprintf (srcFile, "#include \"polyDLLoop.c\"\n\n");
    srcSize += fprintf (srcCalledFile, "}\n");
    fclose (srcCalledFile);
  }
  totalSourceSize += srcSize;

  if (result->eType == T_CONSTANT)
    totalSourceSize += fprintf (srcFile, "\n\treturn %.*g;\n}\n", DBL_DIG, result->value);
  else
    totalSourceSize += fprintf (srcFile, "\n\treturn %s[%lu];\n}\n", eTypes[result->eType], (unsigned long) result->value);

#ifdef POLYCODE_DL
  totalSourceSize += fprintf (srcFile, "#ifdef MAIN\n// Compile with gcc -o %s %s.c -DMAIN -I $KELVIN_ROOT/utils\n"
			      "double (*baseFunction)() = %s;\n#include \"polyDLMain.c\"\n#endif\n", name, name, name);
#endif
  fclose (srcFile);

  if ((includeFile = fopen (includeFileName, "w")) == NULL) {
    perror ("Cannot open polynomial include file\n");
    exit (EXIT_FAILURE);
  }
  totalSourceSize += fprintf (includeFile, "#define VARIABLESUSED %d\n", variableCount);
  totalSourceSize += fprintf (includeFile, "#define SUMSUSED %d\n", sumsUsed);
  totalSourceSize += fprintf (includeFile, "#define PRODUCTSUSED %d\n", productsUsed);
  totalSourceSize += fprintf (includeFile, "#define FUNCTIONCALLSUSED %d\n", functionCallsUsed);
#ifdef POLYCODE_DL
  totalSourceSize += fprintf (includeFile, "#define DLFUNCTIONCOUNT %d\n", fileCount);
#endif
  fclose (includeFile);
  swLogProgress(3 /* DETAIL + 1 */, 0,
		"Polynomial final internal size is %d, as code is %d",
		totalInternalSize, totalSourceSize);

#ifdef POLYCOMP_DL
  char command[PATH_MAX];
  swPushPhase ('k', "compile poly");
  sprintf (command, "$KELVIN_ROOT/compileDL.sh %s", name);
  int status;
  if ((status = system (command)) != 0) {
    perror ("system()");
    exit (EXIT_FAILURE);
  }
  swPopPhase ('k');
  /*
  sprintf (command, "echo internally %d, as code %d, and externally `wc -c %s.so`.", 
	   totalInternalSize, totalSourceSize, name);
  if ((status = system (command)) != 0) {
    perror ("system()");
    exit (EXIT_FAILURE);
  }
  */
#endif

  swPopPhase ('k');
  return;
}

#ifdef DEPENDENCYFLAGGING
unsigned long totalNodes = 0;
/**

  Traverse the polynomial specified racking-up dependency information
  in terms of variables.

*/
void doDependencyFlagging (Polynomial * p)
{
  int i;

  if (p->valid & VALID_EVAL_FLAG)
    return;
  p->valid |= VALID_EVAL_FLAG;

  switch (p->eType) {
  case T_OFFLINE:
    importPoly(p);
    // Notice that there's no break here...

  case T_CONSTANT:
    p->dependencyFlag = 0;
    break;
  case T_VARIABLE:
    p->dependencyFlag = (1UL << p->index);
    break;
  case T_EXTERNAL:
    FATAL ("Don't know how to flag dependencies for EXTERNAL");
    break;
  case T_SUM:
    totalNodes++;
    importTermList (p);
    for (i = 0; i < p->e.s->num; i++) {
      doDependencyFlagging (p->e.s->sum[i]);
      p->dependencyFlag |= p->e.s->sum[i]->dependencyFlag;
    }
    exportTermList (p, FALSE);
  for (i = 0; i < variableCount; i++)
    if (p->dependencyFlag & (1UL << i))
      variableList[i]->e.v->dependentNodes++;
    break;
  case T_PRODUCT:
    totalNodes++;
    for (i = 0; i < p->e.p->num; i++) {
      doDependencyFlagging (p->e.p->product[i]);
      p->dependencyFlag |= p->e.p->product[i]->dependencyFlag;
    }
  for (i = 0; i < variableCount; i++)
    if (p->dependencyFlag & (1UL << i))
      variableList[i]->e.v->dependentNodes++;
    break;
  case T_FUNCTIONCALL:
    totalNodes++;
    for (i = 0; i < p->e.f->num; i++) {
      p->dependencyFlag |= p->e.f->para[i]->dependencyFlag;
      doDependencyFlagging (p->e.f->para[i]);
    }
  for (i = 0; i < variableCount; i++)
    if (p->dependencyFlag & (1UL << i))
      variableList[i]->e.v->dependentNodes++;
    break;
  case T_FREED:
    expTermPrinting (stderr, p, 1);
    fprintf (stderr, "\n");
    FATAL ("In doDependencyFlagging, preceeding polynomial was freed");
    break;
  default:
    FATAL ("In doDependencyFlagging, unknown expression type: [%d]", p->eType);
  }
}

/**

  Wrapper for doDependencyFlagging to avoid recursing initialization and
  termination.

*/
void dependencyFlagging (Polynomial * p)
{
  int i;

  swLogProgress(3 /* DETAIL + 1 */, 0, "Resetting all variable dependency flags");
  totalNodes = 0;
  for (i=0; i<variableCount; i++)
    variableList[i]->e.v->dependentNodes = 0;

  clearValidEvalFlag ();
  doDependencyFlagging (p);

  swLogProgress(3 /* DETAIL + 1 */, 0, "There are %lu nodes", totalNodes);
  for (i=0; i<variableCount; i++)
    if (totalNodes != 0)
      INFO ("Variable %s referenced by %4.1d%% nodes", variableList[i]->e.v->vName,
	       (int) (variableList[i]->e.v->dependentNodes * 100.0 / totalNodes));
}
#endif
