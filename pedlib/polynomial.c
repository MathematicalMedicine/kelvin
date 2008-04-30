
#undef FREEDEBUG

/**********************************************************************
 * polynomial - build and evaluate arbitrarily complex polynomials.
 * Hongling Wang
 * 
 * Polynomial reduction optimizations - Bill Valentine-Cooper
 * 
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

/*
  Builds and evaluates arbitrarily complex polynomials. Collects
  subpolynomials of the same type (addition or multiplication) into
  a single tier. This process will cause the tree to be alternating
  layers of additions and multiplications. Identifies redundant 
  polynomials and maintains only one copy with multiple references.

  Usage:

  1. Call polynomialInitialization before any other routines.
  2. Create constants and variables with calls to constantExp and
  variableExp. Use the returned Polynomials in future calls.
  3. Call plusExp to build a new Polynomial by adding subpolynomials
  multiplied by double factors. Use returned Polynomial in future calls.
  4. Call timesExp to build a new Polynomial by multiplying
  subpolynomials raised to integer powers. Use returned Polynomial in
  future calls.
  5. Set the values of the referenced variables.
  6. Evaluate your Polynomial by calling evaluatePoly.
  
  Limits:

  Maximum constant accuracy is 9 digits due to integer comparison.
  Total polynomials ever seen (kept or not) is INT_MAX due to nodeId.
  
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "polynomial.h"
//#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"

/* This is a global variable used for giving each polynomial an unique ID
   so that we can know if two polynomials are the same just from their IDs */
int nodeId;

/* We have a polynomials list to save the polynomials in each polynomial category,
   which is named constantList, variableList, sumList, productList, and functionCallList.
   The length of these lists are saved in constantCount, variableCount, productCount,
   functionCallCount respectively.  Each of the lists is dynamically applied since we are
   not sure how many polynomials we may have.  Therefore, we apply to have a short list
   initially and increase its length gradually with the increase of the number of 
   polynomials. The current length of the lists is recorded in constantListLength,
   variableListLength, sumListLength, productListLength, and functionCallListLength
   respectively. */

/* Why not free all the polynomials before building the likelihood polynomials for the 
   next trait position? Before we even start to build the likelihood polynomials of 
   pedigrees at the first trait position, some polynomials have already been built. 
   These polynomials are useful for all the trait positions.  Therefore, we can't free 
   everything before we start to build likelihood polynomials of pedigrees at the next 
   trait position. */

struct polynomial **constantList;
int constantCount;
int constantListLength;
int constantCountStamp;
int constantCountStamp2;

struct polynomial **variableList;
int variableCount;
int variableListLength;
int variableCountStamp;
int variableCountStamp2;

struct polynomial **sumList;
int sumCount;
int sumListLength;
int sumCountStamp;
int sumCountStamp2;

struct polynomial **productList;
int productCount;
int productListLength;
int productCountStamp;
int productCountStamp2;

struct polynomial **functionCallList;
int functionCallCount;
int functionCallListLength;
int functionCallCountStamp;
int functionCallCountStamp2;

/* These used to be constants, hence their case. Now they're scaled, but I kept the
   case because they act like constants. */

int CONSTANT_HASH_SIZE, VARIABLE_HASH_SIZE, SUM_HASH_SIZE, PRODUCT_HASH_SIZE,
  FUNCTIONCALL_HASH_SIZE, HASH_TABLE_INCREASE, CONSTANT_LIST_INITIAL,
  CONSTANT_LIST_INCREASE, VARIABLE_LIST_INITIAL, VARIABLE_LIST_INCREASE,
  SUM_LIST_INITIAL, SUM_LIST_INCREASE, PRODUCT_LIST_INITIAL,
  PRODUCT_LIST_INCREASE, FUNCTIONCALL_LIST_INITIAL,
  FUNCTIONCALL_LIST_INCREASE = 0;

/* Variables for tracking internal polynomial memory usage. These dynamic ones are
   also externs in polynomial.h. */
int maxHashLength = 0;		/* Length of longest hash table collision list */
int constantHashHits = 0,	/* constantExp return is a pre-existing polynomial */
  variableHashHits = 0,		/* variableExp return is a pre-existing polynomial (surprise!) */
  functionHashHits = 0;		/* functionCallExp return is a pre-existing polynomial */
int sumReleaseableCount = 0,	/* Indicates if flag was set to release 1st term in sum */
  sumNotReleaseableCount = 0,	/* ...or not. */
  sumReturnConstantCount = 0,	/* plusExp return is actually a constant */
  sumReturn1TermCount = 0,	/* plusExp return is single-term polynomial */
  sumHashHits = 0,		/* plusExp return is pre-existing polynomial */
  sumNewCount = 0,		/* plusExp return is a new polynomial */
  sumListNewCount = 0,		/* New sumList entry created */
  sumListReplacementCount = 0,	/* Existing sumList entry used */
  sumFreedCount = 0,		/* Count of sumPolys freed */
  sum1stTermsFreedCount = 0;	/* Count of 1st-terms successfully freed */
int productReleaseableCount = 0,	/* Indicates if flag was set to release 1st term in product */
  productNotReleaseableCount = 0,	/* ...or not. */
  productReturn0Count = 0,	/* timesExp return is actually zero */
  productReturnConstantCount = 0,	/* timesExp return is actually a non-zero constant */
  productReturn1stTermCount = 0,	/* timesExp return is actually the first term itself */
  productReturn1TermSumCount = 0,	/* timesExp return is a sum built from the first term itself */
  productHashHits = 0,		/* timesExp return is a pre-existing polynomial */
  productHashHitIsSumCount = 0,	/* timesExp return is a sum built from a pre-existing polynomial */
  productReturnNormalCount = 0,	/* timesExp return is a new product polynomial (not a sum) */
  productNon1FactorIsSumCount = 0,	/* timesExp return is a new sum of a new product */
  productListNewCount = 0,	/* New productList entry created */
  productListReplacementCount = 0,	/* Existing productList entry used */
  productFreedCount = 0,	/* Count of productPolys freed */
  product1stTermsFreedCount = 0;	/* Count of 1st-terms successfully freed */
int constantPListExpansions = 0,	/* Count of constantList expansions */
  variablePListExpansions = 0,	/* Count of variableList expansions */
  sumPListExpansions = 0,	/* Count of sumList expansions */
  productPListExpansions = 0,	/* Count of productList expansions */
  functionCallPListExpansions = 0,	/* Count of functionCallList expansions */
  polyListSortingCount = 0, evaluatePolyCount = 0, evaluateValueCount = 0, keepPolyCount =
  0, freePolysCount = 0, holdPolyCount = 0, holdAllPolysCount =
  0, unHoldPolyCount = 0, freeKeptPolysCount = 0, freePolysAttemptCount = 0;
int containerExpansions = 0;	/* Count of expansions of any term-collection container.  */
unsigned long totalSPLLengths = 0, totalSPLCalls = 0, lowSPLCount =
  0, highSPLCount = 0;
unsigned long totalHashMemory = 0;	/* Total memory used for hash table and collision lists */

char *polynomialVersion = "0.34.2";	/* Make this meaningful since kelvin displays it. */

/* Both of the following are set by initialization to value of environment variable of same name.
   They control diagnostic action in a manner not permitted by other approaches since they can
   be changed without rebuilding. */
int polynomialDebugLevel = 0;	/* Corresponds roughly to diagnostic output volume */
int polynomialLostNodeId = -1;	/* For tracking down mis-freed polynomials */
extern int polynomialScale;	/* Scaling factor for hash and other storage */
extern struct swStopwatch *overallSW;

struct swStopwatch *evaluatePolySW, *evaluateValueSW;

unsigned long lastPDSAccumWallTime = 0, lastPDSAccumUserTime = 0;	/* For thrashing check. */

/* Clear the evaluation flag on the entire tree so we can mark
   where we've been and not retrace our steps regardless of
   redundancy. */
void
clearValidEvalFlag ()
{
  int i;

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

///////////////////////////////////////////////////////////////////////////////////////////////
//Recursively evaluate a polynomial.  This version of polynomial evaluation doesn't use
//polynomial sorting list.  It just compute the values of each sub polynomial recursively.  
//A reused sub polynomial maybe repeatedly evaluated.  The efficiency is lower than the
//function of evaluatePoly.  However, we don't need to build a sorting list of sub polynomials
//before the polynomial can be evaluated.
///////////////////////////////////////////////////////////////////////////////////////////////
double
doEvaluateValue (Polynomial * p)
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
    //If a sub polynomial is a contant, return the value
    //of this constant
  case T_CONSTANT:
    p->valid |= VALID_EVAL_FLAG;
    return p->value;
    //If a sub polynomial is a variable, return the value
    //of this variable
  case T_VARIABLE:
    p->valid |= VALID_EVAL_FLAG;
    if (p->e.v->vType == 'D') {
      p->value = *(p->e.v->vAddr.vAddrD);
      return *(p->e.v->vAddr.vAddrD);
    } else if (p->e.v->vType == 'I') {
      p->value = *(p->e.v->vAddr.vAddrI);
      return *(p->e.v->vAddr.vAddrI);
    } else {
      fprintf (stderr, "Wrong variable type, exit!\n");
      exit (1);
    }

    //If a sub polynomial is a sum, evaluate the values of all the terms.  Add the values of
    //the terms weighted by their coefficients.
  case T_SUM:
    result = 0;
    sP = p->e.s;
    for (i = 0; i < sP->num; i++) {
      if (sP->factor[i] == 1)
	result += doEvaluateValue (sP->sum[i]);
      else
	result += doEvaluateValue (sP->sum[i]) * sP->factor[i];
    }
    p->value = result;
    p->valid |= VALID_EVAL_FLAG;
    return result;

    //If a sub polynomial is a product, evaluate the values of all the terms.  Multiply the
    //values of the terms powered by their exponents
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

    //If a sub polynomial is a function call, evaluate the values of all the parameters
    //and then the function according to the name of the function.  Therefore, a function
    //must have a unique name and a unique set of parameters.
  case T_FUNCTIONCALL:
    fp = p->e.f;
    //log10 is for LOD score computation
    if (strcmp (fp->name, "log10") == 0)
      result = log10 (doEvaluateValue (fp->para[0]));
    //gsl_ran_tdist_pdf,gsl_cdf_tdist_Q,gsl_cdf_tdist_P,
    //gsl_ran_ugaussian_pdf,gsl_cdf_ugaussian_Q,gsl_cdf_ugaussian_P,
    //gsl_cdf_chisq_P, gsl_cdf_chisq_Q, and gsl_ran_chisq_pdf
    //are for quantitative trait gene's penetrance computation
    else if (strcmp (fp->name, "gsl_ran_tdist_pdf") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      result = gsl_ran_tdist_pdf (value0, value1);
    } else if (strcmp (fp->name, "gsl_cdf_tdist_Q") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      result = gsl_cdf_tdist_Q (value0, value1);
    } else if (strcmp (fp->name, "gsl_cdf_tdist_P") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      result = gsl_cdf_tdist_P (value0, value1);
    } else if (strcmp (fp->name, "gsl_ran_ugaussian_pdf") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      result = gsl_ran_ugaussian_pdf (value0);
    } else if (strcmp (fp->name, "gsl_cdf_ugaussian_Q") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      result = gsl_cdf_ugaussian_Q (value0);
    } else if (strcmp (fp->name, "gsl_cdf_ugaussian_P") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      result = gsl_cdf_ugaussian_P (value0);
    } else if (strcmp (fp->name, "gsl_cdf_chisq_P") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      result = gsl_cdf_chisq_P (value0, value1);
    } else if (strcmp (fp->name, "gsl_cdf_chisq_Q") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      result = gsl_cdf_chisq_Q (value0, value1);
    } else if (strcmp (fp->name, "gsl_ran_chisq_pdf") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      result = gsl_ran_chisq_pdf (value0, value1);
    }
    //pow, exp, sqrt are standard functions
    else if (strcmp (fp->name, "pow") == 0) {
      value0 = doEvaluateValue (fp->para[0]);
      value1 = doEvaluateValue (fp->para[1]);
      result = pow (value0, value1);

    } else if (strcmp (fp->name, "exp") == 0) {
      result = exp (doEvaluateValue (fp->para[0]));
    } else if (strcmp (fp->name, "sqrt") == 0) {
      result = sqrt (doEvaluateValue (fp->para[0]));
    } else {
      fprintf (stderr, "unknown function name %s in polynomials\n", fp->name);
      exit (1);
    }
    p->value = result;
    p->valid |= VALID_EVAL_FLAG;
    return result;

    //If the polynomial type is unknown, something must be wrong
  default:
    fprintf (stderr, "Error, Unknown expression type!!!!, exit(1)");
    exit (1);
  }
}
double
evaluateValue (Polynomial * p)
{
  double returnValue;
  int i;

  evaluateValueCount++;
#ifdef EVALUATESW
  swStart (evaluateValueSW);
#endif

  /* Clear all of the VALID_EVAL_FLAGs */
  clearValidEvalFlag ();

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

//////////////////////////////////////////////////////////////////////////////////
//A polynomial is zero if it is a constant polynomial and its value is 0.
//This is used when we build a product polynomial.  If the value of a term is 0,
//then the product will be 0.  We don't need to check other terms and this
//product expression can be build as 0.
/////////////////////////////////////////////////////////////////////////////////
int
isZeroExp (Polynomial * p)
{
  if (p->eType != T_CONSTANT)
    return 0;
  if (p->value == 0.0)
    return 1;
  else
    return 0;
};

/////////////////////////////////////////////////////////////////////////////////
//Binary search in an array of integers
//return 1: found
//       0: unfound
//       location: if not found, this is the position for insertion of the target
//                 in the array
//The polynomial are searched based on their keys.  The key of a polynomial
//is an integer
/////////////////////////////////////////////////////////////////////////////////
inline int
binarySearch (int *array, int length, int target, int *location)
{
  int binaryStart, binaryEnd, binaryMiddle;

  if (length == 0) {
    *location = 0;
    return 0;
  }
  //start from the two ends of the array, cut half of
  //the array each time
  binaryStart = 0;
  binaryEnd = length - 1;
  while (binaryStart <= binaryEnd) {
    binaryMiddle = floor ((binaryStart + binaryEnd) / 2.0);
    //the target is found in the array
    if (target == array[binaryMiddle]) {
      *location = binaryMiddle;
      return 1;
    }
    //target is not found in the array
    else if (target > array[binaryMiddle])
      binaryStart = binaryMiddle + 1;
    else
      binaryEnd = binaryMiddle - 1;
  }				//end of while
  *location = binaryStart;
  return 0;
};

//////////////////////////////////////////////////////////////////////////////////////////
//Search for a polynomial in the hash table with a key.  Because more than one polynomials
//may have the same key, the search in a hash table may return a list of polynomials that
//have the same key as the key of the target polynomial 
//////////////////////////////////////////////////////////////////////////////////////////
inline int
searchHashTable (struct hashStruct *hash, int *start, int *end, int key)
{
  int found;
  int location;

  //If the newly generated key is equal to the key of a polynomial that stays in
  //the same indexed position at the hash table
  if (binarySearch (hash->key, hash->num, key, &location) == 1) {
    *end = location + 1;
    *start = location - 1;
    //Search for more keys that are equal to the newly generated key
    while ((*start) >= 0 && hash->key[*start] == hash->key[location])
      (*start)--;
    //Search for more keys that are equal to the newly generated key
    while ((*end) < hash->num && hash->key[*end] == hash->key[location])
      (*end)++;
    if ((*start) < 0 || key != hash->key[*start])
      (*start)++;
    if ((*end) >= hash->num || key != hash->key[*end])
      (*end)--;
    found = 1;
  }
  //If the newly generated key is not found in hash table, we are sure that the new 
  //polynomial is not in the polynomial lists.
  else {
    *start = location;
    *end = location;
    found = 0;
  }
  return found;
};

//////////////////////////////////////////////////////////////////////////////////////////
//record a polynomial in a hash table.  The record of a polynomial in a hash table
//include the key of the polynomial and its position in the corresponding polynomial list.
//The polynomials recorded in the same hash table index are sorted by their keys in 
//the increasing order
//////////////////////////////////////////////////////////////////////////////////////////
inline void
insertHashTable (struct hashStruct *hash, int location, int key, int index)
{
  hash->num++;
  // If the hash collision list is full, allocate more memory
  if (hash->num > hash->length) {
    hash->length += HASH_TABLE_INCREASE;
    if (hash->length > maxHashLength)
      maxHashLength = hash->length;
    hash->key = realloc (hash->key, sizeof (int) * hash->length);
    hash->index = realloc (hash->index, sizeof (int) * hash->length);
    totalHashMemory += HASH_TABLE_INCREASE * sizeof (int) * 2;
    if (hash->key == NULL || hash->index == NULL) {
      fprintf (stderr, "Memory allocation for hash table failed!\n");
      exit (1);
    }
  }
  // Prepare a place in hash table for the new polynomial
  if (location <= hash->num - 2) {
    memmove (&hash->key[location + 1], &hash->key[location],
	     sizeof (int) * (hash->num - location - 1));
    memmove (&hash->index[location + 1], &hash->index[location],
	     sizeof (int) * (hash->num - location - 1));
  }
  //Insert the new polynomiual in the hash table
  hash->key[location] = key;
  hash->index[location] = index;

};

//////////////////////////////////////////////////////////////////////////////////////////
//delete the record of a polynomial from a hash table designated by the parameter of
//location, where the key and index are recorded in key list and index list respectively
//////////////////////////////////////////////////////////////////////////////////////////
inline void
deleteHashTable (struct hashStruct *hash, int location)
{
  if (location >= hash->num) {
    fprintf (stderr, "Deletion in hash table failed\n");
    exit (1);
  }

  memmove (&hash->key[location], &hash->key[location + 1],
	   sizeof (int) * (hash->num - location - 1));
  memmove (&hash->index[location], &hash->index[location + 1],
	   sizeof (int) * (hash->num - location - 1));
  hash->num--;
};

//////////////////////////////////////////////////////////////////////////////////////////
//compute the key of a constant polynomial from the normalized fraction and the exponent
//generated from function  frexp, which converts a floating-point number to fractional 
//and integral components 
//////////////////////////////////////////////////////////////////////////////////////////
inline int
keyConstantPolynomial (double tempD1, int tempI1)
{
  int key;

  //The key is a number between 0 and MAX_POLYNOMIAL_KEY
  if (tempD1 > 0)
    key = floor ((tempD1 - 0.5) * 2 * MAX_POLYNOMIAL_KEY) + tempI1;
  else
    key = floor ((tempD1 + 0.5) * 2 * MAX_POLYNOMIAL_KEY) + tempI1;

  return key;

};

/////////////////////////////////////////////////////////////////////////////////////////////
//Construct a constant polynomial according the the constant value provided by the parameter.
//If a constant polynomial which has the same constant value exists, then no new constant 
//polynomial is constructed.  Instead, the existing constant polynomial is returned.  In
//construction of likelihood polynomials, the number of constant polynomials generated is 
//pretty small compared with the number of sum and product polynomials and therefore is
//not a big burden for memory supply
/////////////////////////////////////////////////////////////////////////////////////////////
Polynomial *
constantExp (double con)
{
  int i;
  Polynomial *p;		//newly built constant polynomial
  int key;			//key of the new constant polynomial
  int hIndex, cIndex;		//index in hash table, index in constant polynomial list
  int first, last, location;	//first and last polynomial in a hash bucket for comparison, 

  //location of the newly built constant polynomial in the bucket
  int tempI1 = -1, tempI2 = -1;	//normalized integer component
  double tempD1, tempD2;	//normalized fractional component

  //Compute a key for the constant  
  tempD1 = frexp (con, &tempI1);
  key = keyConstantPolynomial (tempD1, tempI1);

  //Compute the index of this constant polynomial in the hash table of constant polynomials
  hIndex = key % CONSTANT_HASH_SIZE;
  if (hIndex < 0)
    hIndex += CONSTANT_HASH_SIZE;

  //If the hash table item is not empty, determine if the constant has been in the 
  //constant polynomial list.  If it is not, determine a position in the hash table to
  //save the key and the index in the constant list of this constant polynomial and 
  //save this constant polynomial in the constant polynomial list
  if (constantHash[hIndex].num > 0) {
    //if the key of this constant is equal to the keys of some polynomials in the 
    //constant list, compare if the value of the new constant is equal to
    //the value of an existing constant
    if (searchHashTable (&constantHash[hIndex], &first, &last, key)) {
      for (i = first; i <= last; i++) {
	cIndex = constantHash[hIndex].index[i];
	tempD2 = frexp (constantList[cIndex]->value, &tempI2);
	//Compare if the two constants are the same
	if (tempI2 == tempI1
	    && (int) (tempD2 * 100000000) == (int) (tempD1 * 100000000)) {
	  //if the two constants are the same, return the constant in the constant list 
	  constantList[cIndex]->valid |= VALID_REF_FLAG;
	  constantHashHits++;
	  return constantList[cIndex];
	}
      }
      location = last;
    }
    //If the newly generated key is unique, we are sure that this constant appears
    //the first time.
    else {
      location = last;
    }
  }
  //If the hash table item is empty, we are sure that this constant appears
  //the first time
  else {
    location = 0;
  }

  //next, insert it into the constant list

  //Generate a constant polynomial
  p = (Polynomial *) malloc (sizeof (Polynomial));
  if (p == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  p->eType = T_CONSTANT;
  p->value = con;

  //check if the constant polynomial list is full.  Apply for more items if it is full
  if (constantCount >= constantListLength) {
    constantListLength += CONSTANT_LIST_INCREASE;
    constantPListExpansions++;
    constantList =
      realloc (constantList, constantListLength * sizeof (Polynomial *));
    if (constantList == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }
  }
  //save the constant in the constant polynomial list
  constantList[constantCount] = p;
  p->index = constantCount;
  //give the polynomial a unique ID
  p->id = nodeId;
  constantCount++;
  nodeId++;
  if ((nodeId & 0x1FFFFF) == 0)
    polyStatistics ("At 2M poly multiple");
  p->key = key;
  p->valid = 0;
  p->count = 0;

  //Record the constant polynomial in the hash table of constant polynomials
  insertHashTable (&constantHash[hIndex], location, key, constantCount - 1);

  return p;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//The key of a variable polynomial is computed by the address of the variable.  Different variables have
//different addresses in memory which tend to provide different keys for different variable polynomials.
/////////////////////////////////////////////////////////////////////////////////////////////////////////
inline int
keyVariablePolynomial (double *vD, int *vI, char vType)
{
  int key;

  //currently we can deal with integer and double variables.
  //We use the address of the variable as the key of the variable polynomial
  if (vType == 'D')
    key = (long int) vD;
  else if (vType == 'I')
    key = (long int) vI;
  else {
    fprintf (stderr, "UNKNOWN variable type !");
    exit (1);
  }

  return key;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//This function builds a variable polynomial.  The construction of the product and sum polynomials starts
// from calling this function for constructing variable polynomials or indexing variable polynomials.  
// Therefore, the computation in this function must be very efficient.  The efficiency of this function 
//will affect the efficiency of the construction of almost every term in every polynomial
/////////////////////////////////////////////////////////////////////////////////////////////////////////
Polynomial *
variableExp (double *vD, int *vI, char vType, char name[10])
{
  int i;
  Polynomial *p;		//newly built variable polynomial
  struct variablePoly *vPoly;	//variable structure for newly built variable polynomial
  int key;			//key of the newly buily variable polynomial
  int hIndex, vIndex;		//index in the hash table, index in the polynomial list
  int first, last, location;

  //compuete a key for this variable polynomial
  key = keyVariablePolynomial (vD, vI, vType);


  //Compute the index of this variable polynomial in the hash table of variable polynomials
  hIndex = key % VARIABLE_HASH_SIZE;
  if (hIndex < 0)
    hIndex += VARIABLE_HASH_SIZE;
  //If the hash table item is not empty, determine if the variable has been in the 
  //variable polynomial list.  If it is not, determine a position in the hash table to
  //save the key and the index in the variable list of this variable polynomial and 
  //save this variable polynomial in the variable polynomial list
  if (variableHash[hIndex].num > 0) {

    //if the key of this variable is equal to the keys of some polynomials in the 
    //variable list, compare if the variable has already been there
    if (searchHashTable (&variableHash[hIndex], &first, &last, key)) {
      for (i = first; i <= last; i++) {
	vIndex = variableHash[hIndex].index[i];
	//Compare if the two variables are the same
	if ((vType == 'D' && variableList[vIndex]->e.v->vAddr.vAddrD == vD)
	    || (vType == 'I'
		&& variableList[vIndex]->e.v->vAddr.vAddrI == vI)) {
	  //if the two variables are the same, return the variable in the variable list 
	  variableList[vIndex]->valid |= VALID_REF_FLAG;
	  variableHashHits++;
	  return variableList[vIndex];
	}
      }
      location = last;
    }
    //If the newly generated key is unique, we are sure that this variable appears
    //the first time.
    else {
      location = last;
    }
  }
  //If the hash table item is empty, we are sure that this variable appears
  //the first time
  else {
    location = 0;
  }

  //This variable polynomial doesn't exist.  We create a new variable polynomial
  p = (Polynomial *) malloc (sizeof (Polynomial));
  vPoly = (struct variablePoly *) malloc (sizeof (struct variablePoly));
  if (p == NULL || vPoly == NULL)
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
  p->eType = T_VARIABLE;
  /* Make sure we always have a name, either provided or based upon arrival order. */
  if (strlen (name) == 0)
    sprintf (vPoly->vName, "u%d", variableCount);
  else
    strcpy (vPoly->vName, name);
  //  fprintf(stderr, "New variable %s\n", vPoly->vName);
  if (vType == 'D')
    vPoly->vAddr.vAddrD = vD;
  else
    vPoly->vAddr.vAddrI = vI;
  vPoly->vType = vType;
  p->e.v = vPoly;

  //If the polynomial list is full, apply for more memory
  if (variableCount >= variableListLength) {
    variableListLength += VARIABLE_LIST_INCREASE;
    variablePListExpansions++;
    variableList =
      realloc (variableList, variableListLength * sizeof (Polynomial *));
    if (variableList == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }
  }

  p->index = variableCount;
  p->id = nodeId;
  p->key = key;
  p->valid = 0;
  p->count = 0;

  //Insert the variable polynomial in the variable polynomial list
  variableList[variableCount] = p;
  variableCount++;
  nodeId++;
  if ((nodeId & 0x1FFFFF) == 0)
    polyStatistics ("At 2M poly multiple");

  //Record the variable polynomial in the hash table of the variable polynomials
  insertHashTable (&variableHash[hIndex], location, key, variableCount - 1);

  return p;
};

///////////////////////////////////////////////////////////////////////////////////////////////
//This function compute a key for a sum polynomial.  The number of sum polynomials for complex
//pedigrees and large number of pedigrees tend to be very big.  Therefore it is important that
//we have a good hash strategy so that the polynomials are well distributed in the hash table.
///////////////////////////////////////////////////////////////////////////////////////////////
inline int
keySumPolynomial (Polynomial ** p, double *factor, int counter)
{

  int key = 0;
  int j;
  double tempD1;
  int tempI1;

  //Here we use the ids and keys of all terms in a sum polynomial, their coefficients, their
  //order in the sum to compute a key for the sum polynomial 
  for (j = 0; j < counter; j++) {
    tempD1 = frexp (factor[j], &tempI1);
    key +=
      p[j]->id + p[j]->key * (j + 1) + (int) (tempD1 * 12345 +
					      tempI1) * ((int) p[j]->eType +
							 1);
  }

  return key;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//search for a polynomial in a polynomial list.  If the target polynomial is not in the polynomial list,
//return value:  0 not found
//               1 found
//If the target polynomial is found, then its index in the polynomial list is returned in parameter location
////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline int
searchPolynomialList (Polynomial ** p, int length,
		      Polynomial * target, int *location)
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
  }				//end of while
  *location = binaryStart;
  return 0;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Collect the terms and their coefficients of a sum polynomial.  A plus operation can have any number of
//operands which are polynomials.  We collect these operand polynomials and organize them in a unique order
//as a basis for comparisons between the resulted polynomial and other polynomials 
////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline void
collectSumTerms (double **factor, Polynomial *** p, int *counter,
		 int *containerLength, double f1, Polynomial * p1)
{
  int location;

  //search for the position where the new item should be inserted
  if (searchPolynomialList (*p, *counter, p1, &location) == 1) {
    //this item is currently in the sum, just merge their coefficients
    (*factor)[location] += f1;
    return;
  }
  //If container is full, apply for more memory
  if (*counter >= *containerLength - 1) {
    (*containerLength) += 50;
    containerExpansions++;
    *factor =
      (double *) realloc (*factor, (*containerLength) * sizeof (double));
    *p =
      (Polynomial **) realloc (*p,
			       (*containerLength) * sizeof (Polynomial *));
    if (*factor == NULL || *p == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }
  }
  //this is a new item in the sum, insert it at the start or end of the sum
  //  if (*counter == 0 || location >= *counter) {
  if (location >= *counter) {
    (*p)[*counter] = p1;
    (*factor)[*counter] = f1;
    (*counter)++;
  }
  //insert the item in the middle of the sum
  else {
    //Move the items backward
    memmove (&((*p)[location + 1]), &((*p)[location]),
	     sizeof (Polynomial *) * ((*counter) - location));
    memmove (&((*factor)[location + 1]), &((*factor)[location]),
	     sizeof (double) * ((*counter) - location));

    //insert the new item
    (*p)[location] = p1;
    (*factor)[location] = f1;
    (*counter)++;
  }
};


/////////////////////////////////////////////////////////////////////////////////////////////
//This function generates a sum polynomial.  It accepts a group of <factor, poly> pairs. 
//The sum is represented as factor_1*poly_1+factor_2*poly_2+...+factor_n*poly_n.         
//The sub polynomials are checked to see if it can be combined with other polynomials for
//simplification.  Also, if the sum is a constant, the result will be a constant.        
//If the polynomial being created exists, the existing polynomial is returned.  Otherwise,
//a new polynomial is constructed.
//The parameters are  1) the number of <factor, poly> pairs in the parameter list
//                    2) a group of <factor, poly> pairs
//                    3) a flag whose value is set to be 0 for plus operations like
//                       p=p1 + p2 + p3 ... where p1, p2, p3 may not be freed and
//                       1 for plus operations like p = p  + p2 + p3 ... where the
//                       the parameter polynomial p will be replaced by the new sum 
//                       polynomial and therefore the parameter polynomial p can be
//                       freed   
/////////////////////////////////////////////////////////////////////////////////////////////
Polynomial *
plusExp (int num, ...)
{
  int i, k, l;
  va_list args;			//parameters
  struct sumPoly *sP;		//element of a sum polynomial
  Polynomial *rp;		//pointer of the new sum polynomial
  int counterSum;		//counter of the number of terms in the new sum polynomial
  int flag;			//replacement flag
  double f1, f0;		//factoors
  Polynomial *p1 = 0, *p0 = 0;	//polynomial terms
  double con = 0;		//accumulate constant values appearing in the function parameters
  int key = 0;			//key for the new polynomial
  int tempI, tempI2;		//normalized powers (for comparison of two double numbers)
  double tempD, tempD2;		//normalized fractions( for comparison of two double numbers)
  int hIndex, sIndex;		//index of the new polynomial in the hash table, sub index in a hash bucket
  int first, last, location;	//first polynomial for comparison, last polynomial for comparison, location of the new polynomial in a hash bucket
  int p0SubHIndex = 0, p0HIndex = 0, p0Index = 0, p0Id = 0, p0Key, p0Valid;

  //sub index in a hash bucket, index in the hash table, index in the polynomial list, unique id, key, number of times refered of the
  //polynomial to be replaced
  enum expressionType p0EType;	//polynomial type of the polynomial to be replaced

  //Initialize the variables that keep tracks the number of terms collected for the new sum polynomial
  counter_v1 = 0;
  counter_p1 = 0;
  counter_f1 = 0;

  //get the number of items for this sum from the parameter list of the function
  va_start (args, num);

  //iterate through all the items in the parameter list
  for (i = 0; i < num; i++) {
    //get the coefficient
    f1 = va_arg (args, double);

    //get the polynomial
    p1 = va_arg (args, Polynomial *);

    if (polynomialDebugLevel >= 60) {
      fprintf (stderr, "In plusExp factor=%f item No. %d of %d type=%d\n", f1,
	       i + 1, num, p1->eType);
      if (polynomialDebugLevel >= 70) {
	expTermPrinting (stderr, p1, 1);
	fprintf (stderr, "\n");
      }
    }
    //Record the first operand of the plus operation.  Often, the first operand and the result are the
    //same variable, which refers to the case of p = p + ..., therefore, when we build a new polynomial
    //for the result of the plus operation, the polynomial that represents the first operand of the
    //plus operation can be freed.

    if (i == 0) {
      f0 = f1;
      p0 = p1;
    }
    //If a term is constant 0, it has no effect on a sum, therefore we do nothing
    if (f1 == 0.0) {
      continue;
    }
    switch (p1->eType) {
      //If the term is a non-zero constant, we collect it
    case T_CONSTANT:
      con += p1->value * f1;
      break;
      //The term is either a variable
    case T_VARIABLE:
      collectSumTerms (&factor_v1, &p_v1, &counter_v1,
		       &containerLength_v1, f1, p1);
      break;
      //The item is a product
    case T_PRODUCT:
      collectSumTerms (&factor_p1, &p_p1, &counter_p1,
		       &containerLength_p1, f1, p1);
      break;
      //The term is a function call
    case T_FUNCTIONCALL:
      collectSumTerms (&factor_f1, &p_f1, &counter_f1,
		       &containerLength_f1, f1, p1);
      break;
      //The term is a sum
    case T_SUM:
      //This item is a sum, we go through the items of this sum item
      //each item in this sum is called sum item
      for (l = 0; l < p1->e.s->num; l++) {
	switch (p1->e.s->sum[l]->eType) {
	  //If this sum item is a constant, collected it into con
	case T_CONSTANT:
	  con += f1 * p1->e.s->sum[l]->value * p1->e.s->factor[l];
	  break;
	  //Variable
	case T_VARIABLE:
	  collectSumTerms (&factor_v1, &p_v1, &counter_v1,
			   &containerLength_v1,
			   f1 * p1->e.s->factor[l], p1->e.s->sum[l]);
	  break;
	  //product
	case T_PRODUCT:
	  collectSumTerms (&factor_p1, &p_p1, &counter_p1,
			   &containerLength_p1,
			   f1 * p1->e.s->factor[l], p1->e.s->sum[l]);
	  break;
	case T_FUNCTIONCALL:
	  collectSumTerms (&factor_f1, &p_f1, &counter_f1,
			   &containerLength_f1,
			   f1 * p1->e.s->factor[l], p1->e.s->sum[l]);
	  break;
	default:
	  fprintf (stderr, "In plusExp, unknown expression type %d\n",
		   p1->e.s->sum[l]->eType);
	  exit (1);
	}
      }
      break;
    case T_FREED:
      fprintf (stderr, "In plusExp, polynomial term %d was freed:\n", i);
      expTermPrinting (stderr, p1, 1);
      exit (1);
      break;

    default:
      fprintf (stderr, "In plusExp, unknown polynomial type %d, exiting!\n",
	       p1->eType);
      raise (SIGUSR1);
      exit (1);
    }
  }
  flag = va_arg (args, int);

  va_end (args);

  if (flag == 0)
    sumNotReleaseableCount++;
  else
    sumReleaseableCount++;

  //merge the collected polynomial terms to form a sum polynomial
  counterSum = counter_v1 + counter_p1 + counter_f1;
  //make sure memory is enough for terms merging
  if (counterSum + 1 > lengthSum) {
    lengthSum = counterSum + 1;
    factorSum = (double *) realloc (factorSum, lengthSum * sizeof (double));
    pSum = (Polynomial **) realloc (pSum, lengthSum * sizeof (Polynomial *));
    if (factorSum == NULL || pSum == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }
  }
  /* While we do end-up with zero factors for variables here, they get swallowed-up
     into the next tier where they become non-zero, so it all comes out in the wash. */
  if (counter_v1 > 0) {
    memcpy (&factorSum[0], &factor_v1[0], sizeof (double) * counter_v1);
    memcpy (&pSum[0], &p_v1[0], sizeof (Polynomial *) * counter_v1);
  }
  if (counter_p1 > 0) {
    memcpy (&factorSum[counter_v1], &factor_p1[0],
	    sizeof (double) * counter_p1);
    memcpy (&pSum[counter_v1], &p_p1[0], sizeof (Polynomial *) * counter_p1);
  }
  if (counter_f1 > 0) {
    memcpy (&factorSum[counter_v1 + counter_p1], &factor_f1[0],
	    sizeof (double) * counter_f1);
    memcpy (&pSum[counter_v1 + counter_p1], &p_f1[0],
	    sizeof (Polynomial *) * counter_f1);
  }
  //After we go through all the items in the sum,
  //we get only a constant
  if (counterSum == 0) {
    rp = constantExp (con);
    sumReturnConstantCount++;
    return rp;
    if (polynomialDebugLevel >= 60)
      fprintf (stderr, "Returning a constant %f\n", rp->value);
    return rp;
  }
  //We have only one term for this polynomial.  There is no
  //need to create a new polynomial
  else if (con == 0.0 && counterSum == 1 && factorSum[0] == 1.0) {
    rp = pSum[0];
    sumReturn1TermCount++;
    if (polynomialDebugLevel >= 60) {
      fprintf (stderr, "Returning a single term sum\n");
      expTermPrinting (stderr, rp, 1);
      fprintf (stderr, "\n");
    }
    return rp;
  }
  //The result of plus operation is a more than 2 terms of sum polynomial
  else {
    if (con != 0.0) {
      p1 = constantExp (con);
      factorSum[counterSum] = 1.0;
      pSum[counterSum] = p1;
      counterSum++;
    }
  }

  //compute the key for this polynomial
  key = keySumPolynomial (pSum, factorSum, counterSum);

  hIndex = key % SUM_HASH_SIZE;
  if (hIndex < 0)
    hIndex += SUM_HASH_SIZE;

  //If the hash table item is not empty, determine if the sum has been in the 
  //sum polynomial list.  If it is not, determine a position in the hash table to
  //save the key and the index in the hash table and save this sum polynomial 
  //in the sum polynomial list
  if (sumHash[hIndex].num > 0) {
    //if the key of this sum is equal to the keys of some polynomials in the 
    //sum list, compare if the sum has already been there         
    if (searchHashTable (&sumHash[hIndex], &first, &last, key)) {
      for (i = first; i <= last; i++) {
	sIndex = sumHash[hIndex].index[i];
	if (counterSum == sumList[sIndex]->e.s->num) {
	  //compare the two sums item by item
	  for (k = 0; k < counterSum; k++) {
	    //If the cooresponding terms are identical,
	    //we compare their factors
	    if (pSum[k] == sumList[sIndex]->e.s->sum[k]) {
	      tempD = frexp (factorSum[k], &tempI);
	      tempD2 = frexp (sumList[sIndex]->e.s->factor[k], &tempI2);
	      if (tempI2 != tempI
		  || (int) (tempD2 * 100000000) != (int) (tempD * 100000000))
		break;
	    } else
	      break;
	  }
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
	    return sumList[sIndex];
	  }
	}			//end of if(counter==sumList[sIndex]->e.s->num)
      }				//end of for(i=binaryStart;i<=binaryEnd;i++)
      //identical key, but not identical polynomial
      location = last;
    }				//end of if(binaryStart>=0 && binarySt
    else {
      location = last;
    }
  }				//end of if(sumHash[hIndex].num>0)
  else {
    location = 0;
  }

  //If the first polynomial in the parameter list can be freed
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
      fprintf (stderr,
	       "This polynomial is not in the sum hash table (1), exiting!\n");
      expTermPrinting (stderr, p0, 1);
      fprintf (stderr, "\n");
      exit (1);
    }
    //locate the polynomial in the hash table
    if (searchHashTable (&sumHash[p0HIndex], &first, &last, p0Key)) {
      for (i = first; i <= last; i++) {
	if (p0Index == sumHash[p0HIndex].index[i])
	  break;
      }
      if (i <= last)
	p0SubHIndex = i;
      else {
	fprintf (stderr,
		 "This polynomial is not in the sum hash list (2), exiting!\n");
	expTermPrinting (stderr, p0, 1);
	fprintf (stderr, "\n");
	exit (1);
      }
    } else {
      fprintf (stderr,
	       "This polynomial is not in the sum hash list (3), exiting!\n");
      expTermPrinting (stderr, p0, 1);
      fprintf (stderr, "\n");
      exit (1);
    }
    //free the memory of the first polynomial in the parameter list
    sum1stTermsFreedCount++;
    free (p0->e.s->sum);
    free (p0->e.s->factor);
    free (p0->e.s);
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

  //If the sum is not found in the sum list, a new polynomial is built 
  rp = (Polynomial *) malloc (sizeof (Polynomial));
  if (rp == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  rp->eType = T_SUM;
  sP = (struct sumPoly *) malloc (sizeof (struct sumPoly));
  if (sP == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  sP->num = counterSum;
  sP->sum = (Polynomial **) malloc (counterSum * sizeof (Polynomial *));
  sP->factor = (double *) malloc (counterSum * sizeof (double));
  if (sP->sum == NULL || sP->factor == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  for (i = 0; i < sP->num; i++) {
    sP->factor[i] = factorSum[i];
    sP->sum[i] = pSum[i];
    sP->sum[i]->valid |= VALID_REF_FLAG;
  }
  //assign values to some attributes of the sum polynomial
  rp->e.s = sP;
  rp->key = key;
  rp->valid = 0;
  rp->count = 0;

  //Insert the new built polynomial in sum list

  //If the first polynomial is the parameter list is freed, its position
  //in the polynomial list is occupied by the newly created sum polynomial
  if (flag != 0 && p0EType == T_SUM && p0Valid == 0) {
    sumListReplacementCount++;
    sumList[p0Index] = rp;
    sumList[p0Index]->index = p0Index;
    sumList[p0Index]->id = p0Id;
  }
  //otherwise, insert the new polynomial is the sum polynomial list
  else {
    if (sumCount >= sumListLength) {
      sumListLength += SUM_LIST_INCREASE;
      sumPListExpansions++;
      sumList = realloc (sumList, sumListLength * sizeof (Polynomial *));
      if (sumList == NULL) {
	fprintf (stderr, "Memory allocation failure at %s line %d\n",
		 __FILE__, __LINE__);
	exit (1);
      }
    }
    sumListNewCount++;
    sumList[sumCount] = rp;
    sumList[sumCount]->index = sumCount;
    sumList[sumCount]->id = nodeId;
    if (nodeId == polynomialLostNodeId) {
      fprintf (stderr, "nodeId %d has valid of %d and count of %d\n",
	       polynomialLostNodeId, rp->valid, rp->count);
      expTermPrinting (stderr, rp, 16);
      fprintf (stderr,
	       "\nIf you're in gdb, continue from here to see more.\n");
      raise (SIGUSR1);
    }
    sumCount++;
    nodeId++;
    if ((nodeId & 0x1FFFFF) == 0)
      polyStatistics ("At 2M poly multiple");
  }

  //Insert the newly built polynomial into the Hash table
  if (flag != 0 && p0EType == T_SUM && p0Valid == 0) {
    if (p0SubHIndex != location || p0HIndex != hIndex) {
      //insert the newly built polynomial in the hash table
      insertHashTable (&sumHash[hIndex], location, key,
		       sumList[p0Index]->index);
      //delete the replaced polynomial from the hash table
      if (p0HIndex != hIndex) {
	deleteHashTable (&sumHash[p0HIndex], p0SubHIndex);

      } else			//p0SubHIndex!=binaryEnd
      {
	if (p0SubHIndex < location) {
	  deleteHashTable (&sumHash[p0HIndex], p0SubHIndex);
	} else {
	  deleteHashTable (&sumHash[p0HIndex], p0SubHIndex + 1);
	}
      }
    }
    //If the newly built polynomial and the one to be replaced are at the same position 
    //in the hash table, just replace the key of the old polynomial with the key of the
    //new polynomial
    else			//p0HIndex==hIndex && p0SubHIndex==binaryEnd
    {
      sumHash[hIndex].key[location] = key;
    }
  } else {
    insertHashTable (&sumHash[hIndex], location, key, sumCount - 1);
  }

  sumNewCount++;
  if (polynomialDebugLevel >= 60) {
    fprintf (stderr, "Returning a new sum...\n");
    if (polynomialDebugLevel >= 70) {
      expPrinting (rp);
      fprintf (stderr, "\n");
    }
  }
  return rp;
};

/////////////////////////////////////////////////////////////////////////////////
//This function compute a key for a product polynomial.  The number of
//product polynomial is usually very big.  Good keys for product polynomials
//are required for good distribution of product polynomials in product hash table
////////////////////////////////////////////////////////////////////////////////
inline int
keyProductPolynomial (Polynomial ** p, int *exponent, int counter)
{
  int key = 0;
  int j;

  for (j = 0; j < counter; j++)
    key +=
      p[j]->id + p[j]->key * (j + 1) +
      (exponent[j] * 12345 * (int) (p[j]->eType + 1));

  return key;

}

/////////////////////////////////////////////////////////////////////////////////////////////
//This function collects the terms for constructing a product polynomial.  The terms are sorted
//so that these terms can have a unique order in the constructed product polynomial.
////////////////////////////////////////////////////////////////////////////////////////////
inline void
collectProductTerms (int **exponent, Polynomial *** p, int *counter,
		     int *containerLength, int e1, Polynomial * p1)
{
  int location;

  //search for a position for this component in the container
  if (searchPolynomialList (*p, *counter, p1, &location) == 1) {
    (*exponent)[location] += e1;
    if ((*exponent)[location] == 0)
      fprintf (stderr,
	       "Zero exponent terms could be replaced with constant 1!\n");
    return;
  }

  if ((*counter) >= (*containerLength) - 1) {
    (*containerLength) += 50;
    (*exponent) =
      (int *) realloc ((*exponent), (*containerLength) * sizeof (int));
    (*p) =
      (Polynomial **) realloc ((*p),
			       (*containerLength) * sizeof (Polynomial *));
    if ((*exponent) == NULL || (*p) == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }
  }
  //If this item should be in the start or end
  //  if ((*counter) == 0 || location >= (*counter)) {
  if (location >= (*counter)) {
    (*p)[(*counter)] = p1;
    (*exponent)[(*counter)] = e1;
    (*counter)++;
  }
  //Otherwise, insert it in the container
  else {
    memmove (&((*p)[location + 1]), &((*p)[location]),
	     sizeof (Polynomial *) * ((*counter) - location));
    memmove (&((*exponent)[location + 1]), &((*exponent)[location]),
	     sizeof (int) * ((*counter) - location));
    //for(k=(*counter)-1;k>=j;k--)
    //{
    //  (*p)[k+1]        = (*p)[k];
    //  (*exponent)[k+1] = (*exponent)[k];
    //}
    //put this component in the container
    (*p)[location] = p1;
    (*exponent)[location] = e1;
    (*counter)++;
  }
};


////////////////////////////////////////////////////////////////////////////////////////
//This function create a product polynomial.  It accepts a group of <exponent, poly>  
//pairs as the operands of a times operation of                                       
//poly_1^exponet1*poly_2^exponent_2*...*poly_n^exponent_n.                            
//It checks the operands to see whether they can be combined for simplification.
//It is also possible the result is a polynomial other than type of a product if the
//the result can't be expressed by a product polynomial.
//The parameters are  1) the number of <poly, exponent> pairs in the parameter list
//                    2) a group of <poly, exponent> pairs
//                    3) a flag whose value is set to be 0 for times operations like
//                       p=p1 * p2 * p3 ... where p1, p2, p3 may not be freed and
//                       1 for times operations like p = p  * p2 * p3 ... where the
//                       the parameter polynomial p will be replaced by the new product
//                       polynomial and therefore the parameter polynomial p can be
//                       freed
////////////////////////////////////////////////////////////////////////////////////////

Polynomial *
timesExp (int num, ...)
{
  int i, k, l, counterProd;	//i,k and l are loop variables, counterProd is the number of terms in the new product polynomial
  va_list args;			//parameters
  struct productPoly *pP;	//product structure for the new polynomial
  Polynomial *rp;		//new product polynomial
  Polynomial *p1 = 0, *p0 = 0;	//pointer to a term
  int e1, e0;			//component of a term
  int isZero = 0;		//if the result is zero
  double factor = 1;		//collect the constant terms in the parameter list
  int key = 0;			//key of the new polynomial
  int pIndex, hIndex;		//index in the polynomial list, index in the hash table
  int first, last, location;	//first and last terms for comparison, position of the new polynomial in a hash bucket
  int flag;			//if the new product polynomial can replace an existing one
  int p0SubHIndex = 0, p0HIndex = 0, p0Index = 0, p0Id = 0, p0Key, p0Valid;

  //sub index in a hash bucket, index in the hash table, index in the polynomial list, unique id, key, number of times refered of the
  //polynomial to be replaced
  enum expressionType p0EType;

  //Initialize the containers for the operands of a times operation 
  counter_v2 = 0;
  counter_s2 = 0;
  counter_f2 = 0;

  //Get the number of operands for this times operation
  va_start (args, num);

  //go through operand and its exponent of the product
  for (i = 0; i < num; i++) {
    //get the polynomial
    p1 = va_arg (args, Polynomial *);

    //get the exponent
    e1 = va_arg (args, int);

    if (polynomialDebugLevel >= 60) {
      fprintf (stderr, "In timesExp exponent=%d item No. %d of %d type=%d\n",
	       e1, i + 1, num, p1->eType);
      if (polynomialDebugLevel >= 70) {
	expTermPrinting (stderr, p1, 1);
	fprintf (stderr, "\n");
      }
    }
    //Record the first operand of the times operation.  Often, the first operand and the result are the
    //same variable, which refers to the case of p = p * ..., therefore, when we build a new polynomial
    //for the result of the times operation, the polynomial that represents the first operand of the
    //times operation can be freed.      

    if (i == 0) {
      e0 = e1;
      p0 = p1;
    }
    //If any of the operand is zero, then the product is zero
    if (isZeroExp (p1)) {
      isZero = 1;
      break;
    }
    //If this operand is a constant, this constant and its exponent
    //is accumulated into factor
    else if (p1->eType == T_CONSTANT) {
      factor *= pow (p1->value, e1);
      continue;
    }
    //If the operand is variable/sum/product/function call
    else {
      //If this operand is a sum polynomial and this sum polynomial
      //has only one item
      if (p1->eType == T_SUM && p1->e.s->num == 1) {
	factor *= pow (p1->e.s->factor[0], e1);
	p1 = p1->e.s->sum[0];
      }
      //If this operand is a variable, a sum that has
      //more than one items, or a function call, we directly
      //collect it.  If this operand is a product, we then
      //collect each term of it.
      switch (p1->eType) {
      case T_VARIABLE:
	collectProductTerms (&exponent_v2, &p_v2, &counter_v2,
			     &containerLength_v2, e1, p1);
	break;

      case T_SUM:
	collectProductTerms (&exponent_s2, &p_s2, &counter_s2,
			     &containerLength_s2, e1, p1);
	break;

      case T_FUNCTIONCALL:
	collectProductTerms (&exponent_f2, &p_f2, &counter_f2,
			     &containerLength_f2, e1, p1);
	break;
	//For a product term, we open it and check each of its terms
      case T_PRODUCT:

	for (l = 0; l < p1->e.p->num; l++) {
	  switch (p1->e.p->product[l]->eType) {
	  case T_VARIABLE:
	    collectProductTerms (&exponent_v2, &p_v2, &counter_v2,
				 &containerLength_v2,
				 e1 * p1->e.p->exponent[l],
				 p1->e.p->product[l]);
	    break;

	  case T_SUM:
	    collectProductTerms (&exponent_s2, &p_s2, &counter_s2,
				 &containerLength_s2,
				 e1 * p1->e.p->exponent[l],
				 p1->e.p->product[l]);
	    break;

	  case T_FUNCTIONCALL:
	    collectProductTerms (&exponent_f2, &p_f2, &counter_f2,
				 &containerLength_f2,
				 e1 * p1->e.p->exponent[l],
				 p1->e.p->product[l]);
	    break;
	  default:
	    fprintf (stderr, "In timesExp, unknown expression type %d\n",
		     p1->e.p->product[l]->eType);
	    exit (1);
	  }			//end of switch
	}			//end of for
	break;

      case T_FREED:
	fprintf (stderr, "In plusExp, polynomial term %d was freed:\n", i);
	expTermPrinting (stderr, p1, 1);
	exit (1);
	break;

      default:
	fprintf (stderr,
		 "In timesExp, unknown polynomial type %d, exiting!\n",
		 p1->eType);
	raise (SIGUSR1);
	exit (1);
	break;
      }				//end of switch
    }				//end of else
  }				//end of for

  //Get the last parameter for this function.  It is a parameter to show if the
  //first operand of the times operation as a polynomial can be freed after the
  //result polynomial is constructed 
  flag = va_arg (args, int);

  va_end (args);

  //This is for performance checking use                                                        
  if (flag == 0)
    productNotReleaseableCount++;
  else
    productReleaseableCount++;

  //The product is zero, a zero polynomial is returned
  if (isZero) {
    rp = constantExp (0.0);
    productReturn0Count++;
    if (polynomialDebugLevel >= 60)
      fprintf (stderr, "Returning a constant zero\n");
    return rp;
  }
  //combine all the items generated by collecting the operands of the times operation into one polynomial
  counterProd = counter_v2 + counter_s2 + counter_f2;
  if (counterProd > lengthProd) {
    lengthProd = counterProd;
    exponentProd = (int *) realloc (exponentProd, lengthProd * sizeof (int));
    pProd =
      (Polynomial **) realloc (pProd, lengthProd * sizeof (Polynomial *));
    if (exponentProd == NULL || pProd == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }
  }
  if (counter_v2 > 0) {
    memcpy (&exponentProd[0], &exponent_v2[0], sizeof (int) * counter_v2);
    memcpy (&pProd[0], &p_v2[0], sizeof (Polynomial *) * counter_v2);
  }
  if (counter_s2 > 0) {
    memcpy (&exponentProd[counter_v2], &exponent_s2[0],
	    sizeof (int) * counter_s2);
    memcpy (&pProd[counter_v2], &p_s2[0], sizeof (Polynomial *) * counter_s2);
  }
  if (counter_f2 > 0) {
    memcpy (&exponentProd[counter_v2 + counter_s2], &exponent_f2[0],
	    sizeof (int) * counter_f2);
    memcpy (&pProd[counter_v2 + counter_s2], &p_f2[0],
	    sizeof (Polynomial *) * counter_f2);
  }
  //The product has 0 items, the result is a constant polynomial
  if (counterProd == 0) {
    rp = constantExp (factor);
    productReturnConstantCount++;
    if (polynomialDebugLevel >= 60)
      fprintf (stderr, "Returning a constant %f\n", rp->value);
    return rp;
  }
  //If the result polynomial has only one term, it is not a product polynomial
  else if (counterProd == 1 && exponentProd[0] == 1) {
    //If the factor is 1, then the result polynomial is equal to its first term
    if (factor == 1.0) {
      rp = pProd[0];
      if (polynomialDebugLevel >= 60)
	fprintf (stderr, "Returning first term from caller\n");
      productReturn1stTermCount++;
      return rp;
    }
    //If the factor is not 1, then the result polynomial is a sum polynomial
    else {
      rp = plusExp (1, factor, pProd[0], 0);
      if (polynomialDebugLevel >= 60)
	fprintf (stderr, "Returning via plusExp\n");
      productReturn1TermSumCount++;
      return rp;
    }
  }
  //The result polynomial is a product polynomial
  else {

    //compute the key for the product polynomial
    key = keyProductPolynomial (pProd, exponentProd, counterProd);


    //compute hash table index according to the key
    hIndex = key % PRODUCT_HASH_SIZE;
    if (hIndex < 0)
      hIndex += PRODUCT_HASH_SIZE;
    //If the hash table item is not empty, determine if the product has been in the 
    //product polynomial list.  If it is not, determine a position in the hash table to
    //save the key and the index in the hash table and save this product polynomial 
    //in the product polynomial list
    if (productHash[hIndex].num > 0) {
      //if the key of this product is equal to the keys of some polynomials in the 
      //product list, compare the new product polynomial with the polynomials having 
      //the identical key to determine if the product polynomial has already been constructed         
      if (searchHashTable (&productHash[hIndex], &first, &last, key)) {
	for (i = first; i <= last; i++) {
	  pIndex = productHash[hIndex].index[i];
	  //If the numer of items of the new polynomial is equal to the number of 
	  //items of another polynomial whose key is equal to that of the new polynomial,
	  //we need to compare them term by term to determine if they are identical
	  if (counterProd == productList[pIndex]->e.p->num) {
	    for (k = 0; k < counterProd; k++)
	      if (pProd[k] != productList[pIndex]->e.p->product[k]
		  || exponentProd[k] != productList[pIndex]->e.p->exponent[k])
		break;
	    //If the polynomials are the same
	    if (k >= counterProd) {
	      //this product polynomial has been existing and the factor is 1, return the
	      //existing polynomial
	      if (factor == 1.0) {
		//the attribute count is used as a sign to show that this polynomial is
		//refered in more than one places so that it can't be freed
		productList[pIndex]->valid |= VALID_REF_FLAG;
		productHashHits++;
		if (polynomialDebugLevel >= 60) {
		  fprintf (stderr, "Returning an existing product...\n");
		  if (polynomialDebugLevel >= 70) {
		    expPrinting (productList[pIndex]);
		    fprintf (stderr, "\n");
		  }
		}
		return productList[pIndex];
	      }
	      //this product polynomial has been existing so that we don't need to construct a 
	      //new one.  However, the factor is not 1, therefore the result polynomial is 
	      //a sum polynomial
	      else {
		productHashHitIsSumCount++;
		productList[pIndex]->valid |= VALID_REF_FLAG;
		if (polynomialDebugLevel >= 60)
		  fprintf (stderr,
			   "Returning existing product via plusExp\n");
		return plusExp (1, factor, productList[pIndex], 0);
	      }
	    }			//end of if
	  }			//end of if
	}			//end of for
	location = last;
      }				//end of if
      else {
	location = last;
      }
    }				//end of if
    else {
      location = 0;
    }

    //If the first operand of the times operation can be freed (flag==1)
    // and the first operand is a product that is refered only once
    if (flag != 0 && p0->eType == T_PRODUCT && p0->valid == 0) {
      //We save the identity information of the polynomial to be freed
      //so that the resource can be assigned to the newly created polynomial 
      p0Index = p0->index;
      p0Id = p0->id;
      p0Key = p0->key;
      //      p0HIndex = p0->key % PRODUCT_HASH_SIZE;
      p0HIndex = p0Key % PRODUCT_HASH_SIZE;
      if (p0HIndex < 0)
	p0HIndex += PRODUCT_HASH_SIZE;
      p0EType = p0->eType;
      p0Valid = p0->valid;

      //Determine the sub index of the old polynomial in the product hash table
      if (productHash[p0HIndex].num <= 0) {
	fprintf (stderr,
		 "This polynomial is not in the product hash list (1), exiting!\n");
	expTermPrinting (stderr, p0, 1);
	fprintf (stderr, "\n");
	exit (1);
      }
      if (searchHashTable (&productHash[p0HIndex], &first, &last, p0Key)) {
	for (i = first; i <= last; i++) {
	  if (p0Index == productHash[p0HIndex].index[i])
	    break;
	}
	if (i <= last)
	  p0SubHIndex = i;
	else {
	  fprintf (stderr,
		   "This polynomial is not in the product hash list (2), exiting!\n");
	  expTermPrinting (stderr, p0, 1);
	  fprintf (stderr, "\n");
	  exit (1);
	}
      } else {
	fprintf (stderr,
		 "This polynomial is not in the product hash list (3), exiting!\n");
	expTermPrinting (stderr, p0, 1);
	fprintf (stderr, "\n");
	exit (1);
      }

      //Free the first operand
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

    //Construct a new product polynomial from the terms
    //saved in the container
    rp = (Polynomial *) malloc (sizeof (Polynomial));
    if (rp == NULL)
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
    rp->eType = T_PRODUCT;
    pP = (struct productPoly *) malloc (sizeof (struct productPoly));
    if (pP == NULL)
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
    pP->num = counterProd;
    pP->product =
      (Polynomial **) malloc (counterProd * sizeof (Polynomial *));
    pP->exponent = (int *) malloc (counterProd * sizeof (int));
    if (pP->product == NULL || pP->exponent == NULL)
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
    //copy the collected polynomial terms and exponents into the
    //product polynomial structure
    for (i = 0; i < counterProd; i++) {
      pP->product[i] = pProd[i];
      pP->exponent[i] = exponentProd[i];
      pP->product[i]->valid |= VALID_REF_FLAG;
    }
    rp->e.p = pP;
    rp->index = productCount;
    rp->id = nodeId;
    rp->key = key;
    rp->valid = 0;
    rp->count = 0;

    //After the new polynomial is built, it is recorded in the product polynomial list
    if (productCount >= productListLength) {
      productListLength += PRODUCT_LIST_INCREASE;
      productPListExpansions++;
      productList =
	realloc (productList, productListLength * sizeof (Polynomial *));
      if (productList == NULL) {
	fprintf (stderr, "Memory allocation failure at %s line %d\n",
		 __FILE__, __LINE__);
	exit (1);
      }
    }
    //We either apply for a new position in the product list for this polynomial, or
    //replace an existing polynomial with the newly created polynomial at the position
    //occupied by the existing polynomial
    if (flag != 0 && p0EType == T_PRODUCT && p0Valid == 0) {
      //Assign the resource of the freed polynomial to the newly constructed polynomial
      productListReplacementCount++;
      productList[p0Index] = rp;
      productList[p0Index]->index = p0Index;
      productList[p0Index]->id = p0Id;
    } else {
      //save the newly constructed polynomial in the polynomial list
      productListNewCount++;
      productList[productCount] = rp;
      productList[productCount]->index = productCount;
      productList[productCount]->id = nodeId;
      if (polynomialDebugLevel >= 40) {
	fprintf (stderr,
		 "Polynomial %d, (product %d) added: ", nodeId, productCount);
	expTermPrinting (stderr, rp, 1);
	fprintf (stderr, "\n");
      }
      if (nodeId == polynomialLostNodeId) {
	fprintf (stderr, "nodeId %d has valid of %d and count of %d\n",
		 polynomialLostNodeId, rp->valid, rp->count);
	expTermPrinting (stderr, rp, 16);
	fprintf (stderr,
		 "\nIf you're in gdb, continue from here to see more.\n");
	raise (SIGUSR1);
      }
      productCount++;
      nodeId++;
      if ((nodeId & 0x1FFFFF) == 0)
	polyStatistics ("At 2M poly multiple");
    }

    //the new polynomial is also recorded in the hash table
    if (flag != 0 && p0EType == T_PRODUCT && p0Valid == 0) {

      //If the indexes or sub indexes of the new and old polynomials in the hash table
      //are different
      if (p0SubHIndex != location || p0HIndex != hIndex) {
	//Insert the new polynomial in the hash table of the product polynomial
	insertHashTable (&productHash[hIndex], location, key, p0Index);

	//Delete the old polynomial from the hash table
	if (p0HIndex != hIndex) {
	  deleteHashTable (&productHash[p0HIndex], p0SubHIndex);
	} else			//p0SubHIndex!=binaryEnd
	{

	  if (p0SubHIndex < location) {
	    deleteHashTable (&productHash[p0HIndex], p0SubHIndex);
	  } else {
	    deleteHashTable (&productHash[p0HIndex], p0SubHIndex + 1);
	  }
	}
      }
      //If the indexes and sub indexes of the old and new polynomials in the hash
      //table are identical, we don't need deletion and insertion any more, just 
      //replace the key of the old polynomial with the key of the new polynomial
      //in the hash table
      else			//p0HIndex==hIndex && p0SubHIndex==binaryEnd
      {
	productHash[hIndex].key[location] = key;
      }
    } else
      //There is no polynomial replacement problem, just insert the newly created polynomial
      //in the hash table
    {
      insertHashTable (&productHash[hIndex], location, key, productCount - 1);
    }

    //If the factor is 1, return a product polynomial
    if (factor == 1.0) {
      productReturnNormalCount++;
      if (polynomialDebugLevel >= 60) {
	fprintf (stderr, "Returning a new product\n");
	if (polynomialDebugLevel >= 70) {
	  expPrinting (rp);
	  fprintf (stderr, "\n");
	}
      }
      return rp;
    }
    //If the factor is not 1, return a sum polynomial
    else {
      productNon1FactorIsSumCount++;
      if (polynomialDebugLevel >= 60)
	fprintf (stderr, "Returning new product via plusExp\n");
      return plusExp (1, factor, rp, 0);
    }
  }
};

////////////////////////////////////////////////////////////////////////
//This function creates a key for a function call polynomial according
//to the name of the called function and the parameters               
////////////////////////////////////////////////////////////////////////
inline int
keyFunctionCallPolynomial (char *fName, Polynomial ** p, int num)
{
  int key = 0;
  int i;

  //compute a key from function name,
  for (i = 0; i < strlen (fName); i++)
    key += (int) fName[i];
  //and parameters
  for (i = 0; i < num - 1; i++) {
    key = key + p[i]->key;
  }

  return key;

};

////////////////////////////////////////////////////////////////////////////////////////
//This function create a function call.  It accepts the function name and parameters of
//the function.  Either a different function is called, or a same function is called but 
//with different parameters, a new function call polynomials is created.  If, on the other
//hand, a same function is called and with the same parameters, the existing function call
//polynomial is returned.
//The parameters include
//  (1) The number of parameters to the called function plus 1
//  (2) The name of the called function
//  (3) A group of parameters to the called function
////////////////////////////////////////////////////////////////////////////////////////
Polynomial *
functionCallExp (int num, ...)
{
  int i, k, hIndex, fIndex;	//hIndex: index in hash table, fIndex: index in polynomiall list
  char *fName;			//function name
  int key;			//key of the newly built function call polynomial
  Polynomial **p;		//parameters of the newly built function call polynomial
  int first, last, location;	//first and last polynomial in a hash bucket for comparison

  //location of the newly built variable polynomial in the hash bucket
  va_list args;			//arguments
  Polynomial *rp;		//newly built function call polynomial
  struct functionPoly *fP;	//function-call structure of the newly built function call polynomial

  //get the number of parameters for functionCallExp
  va_start (args, num);
  fName = va_arg (args, char *);

  p = (Polynomial **) malloc ((num - 1) * sizeof (Polynomial *));
  if (p == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  //num is equal to 1 plus the number of parameters for the called function
  for (i = 0; i < num - 1; i++) {
    for (i = 0; i < num - 1; i++) {
      if (polynomialDebugLevel >= 60) {
	fprintf (stderr,
		 "In functionCallExp item No. %d of %d type=%d: ",
		 i + 1, num, p[i]->eType);
	expPrinting (p[i]);
	fprintf (stderr, "\n");
      }
      p[i] = va_arg (args, Polynomial *);
    }
  }
  va_end (args);

  //compute a key for this polynomial
  key = keyFunctionCallPolynomial (fName, p, num);

  //compute the index in hash table according to the key

  hIndex = key % FUNCTIONCALL_HASH_SIZE;
  if (hIndex < 0)
    hIndex += FUNCTIONCALL_HASH_SIZE;

  //compare the key of this polynomial with the keys of polynomials in the
  //function list
  if (functionCallHash[hIndex].num > 0) {
    //if the key of this function call is equal to the keys of some polynomials in the 
    //function call list, compare if the value of the new function call is equal to
    //the value of an existing function call
    if (searchHashTable (&functionCallHash[hIndex], &first, &last, key)) {
      //compare the new function call polynomials with a set of possible matching polynomials
      for (i = first; i <= last; i++) {
	fIndex = functionCallHash[hIndex].index[i];
	//If the names of the called functions are identical and the numbers of parameters are identical,
	//we compare their parameters
	if (strcmp (fName, functionCallList[fIndex]->e.f->name) == 0
	    && num - 1 == functionCallList[fIndex]->e.f->paraNum) {
	  //compare the two function calls item by item
	  for (k = 0; k < num - 1; k++)
	    if (p[k] != functionCallList[fIndex]->e.f->para[k])
	      break;
	  //if the two function calls are the same, 
	  //return the function call in the function call list 
	  if (k >= num - 1) {
	    free (p);
	    functionHashHits++;
	    return functionCallList[fIndex];
	  }
	}			//end of if(num-1==functionCallList[fIndex]->e.f->paraNum)
      }				//end of for(i=binaryStart;i<=binaryEnd;i++)
      location = last;
    }				//end of if(binaryStart>=0 && binarySt
    else {
      location = last;
    }
  }				//end of if(functionCallHash[hIndex].num>0)
  else {
    location = 0;
  }


  //If the function call is not found in the list, insert it in the list
  //Build a new polynomial
  rp = (Polynomial *) malloc (sizeof (Polynomial));
  fP = (struct functionPoly *) malloc (sizeof (struct functionPoly));
  if (rp == NULL || fP == NULL)
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
  rp->eType = T_FUNCTIONCALL;
  fP->paraNum = num - 1;
  fP->para = (Polynomial **) malloc ((num - 1) * sizeof (Polynomial *));
  fP->name = (char *) malloc (strlen (fName) + 1);
  if (fP->para == NULL || fP->name == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  strcpy (fP->name, fName);
  for (i = 0; i < fP->paraNum; i++) {
    fP->para[i] = p[i];
  }
  rp->e.f = fP;
  rp->index = functionCallCount;
  rp->id = nodeId;
  rp->key = key;
  rp->valid = 0;
  rp->count = 0;

  //Insert the new built polynomial in function call list
  if (functionCallCount >= functionCallListLength) {
    functionCallListLength += FUNCTIONCALL_LIST_INCREASE;
    functionCallPListExpansions++;
    functionCallList =
      realloc (functionCallList,
	       functionCallListLength * sizeof (Polynomial *));
    if (functionCallList == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }
  }
  functionCallList[functionCallCount] = rp;
  functionCallCount++;
  if (polynomialDebugLevel >= 40)
    fprintf (stderr,
	     "Polynomial %d, (function %d) added\n",
	     nodeId, functionCallCount);
  nodeId++;
  if ((nodeId & 0x1FFFFF) == 0)
    polyStatistics ("At 2M poly multiple");

  //insert the polynomial in the hash table of the function call polynomials
  insertHashTable (&functionCallHash[hIndex], location, key,
		   functionCallCount - 1);

  free (p);

  return rp;

};

///////////////////////////////////////////////////////////////////////////////////////
//This function initialize a structure for creating a polynomial list which is used  
//for polynomial evaluation                                                          
///////////////////////////////////////////////////////////////////////////////////////
struct polyList *
buildPolyList ()
{
  struct polyList *l;

  /* Clear all of the VALID_EVAL_FLAGs */
  clearValidEvalFlag ();

  l = (struct polyList *) malloc (sizeof (struct polyList));
  if (l == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }

  l->listSize = 1000;
  l->listNext = 0;
  l->pList = (Polynomial **) malloc (sizeof (Polynomial *) * l->listSize);
  if (l->pList == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }

  return l;
};

//////////////////////////////////////////////////////////////////////////////////////
//This function appends a polynomial onto a polynomial list                         
//////////////////////////////////////////////////////////////////////////////////////
void
polyListAppend (struct polyList *l, Polynomial * p)
{
  //valid is a mark showing that this polynomial appears on a sorting list
  p->valid |= VALID_EVAL_FLAG;

  if (l->listNext >= l->listSize) {
    l->pList =
      realloc (l->pList, sizeof (Polynomial *) * (l->listSize + 1000));
    if (l->pList == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }
    l->listSize = l->listSize + 1000;
  }
  l->pList[l->listNext] = p;
  l->listNext++;
};

/////////////////////////////////////////////////////////////////////////////////////
// Sort the polynomial list to prepare for polynomial evaluation.  Polynomials are 
//evaluated in order.  This function is to determine a order for evaluation.  Also,
//shared terms are evaluated only once.  Together with polyListAppend(), this function
//construct a evaluation list for evaluation of the polynomial.  For linkage computation, 
//it allows shared polynomials within each likelihood polynomial and across likelihood 
//polynomials of a set of pedigrees fully reused
/////////////////////////////////////////////////////////////////////////////////////
void
doPolyListSorting (Polynomial * p, struct polyList *l)
{
  int i;

  switch (p->eType) {
    //If the polynomial is a constant, put it in the evaluation list
  case T_CONSTANT:
    polyListAppend (l, p);
    break;

    //If the polynomial is a variable, put it in the evaluation list
  case T_VARIABLE:
    if (!(p->valid & VALID_EVAL_FLAG)) {
      polyListAppend (l, p);
    }
    break;

    //If the polynomial is a sum, put all the terms of the sum in the evaluation list
    //except constants and then put the sum in the evaluation list
  case T_SUM:

    if (p->valid & VALID_EVAL_FLAG)
      break;

    for (i = 0; i < p->e.s->num; i++)
      if (p->e.s->sum[i]->eType != T_CONSTANT
	  && (!(p->e.s->sum[i]->valid & VALID_EVAL_FLAG))) {
	doPolyListSorting (p->e.s->sum[i], l);
      }
    polyListAppend (l, p);
    break;

    //If the polynomial is a product, put all the terms of the product in the 
    //evaluation list except constants and then put the product in the evaluation list
  case T_PRODUCT:

    if (p->valid & VALID_EVAL_FLAG)
      break;

    for (i = 0; i < p->e.p->num; i++)
      if (p->e.p->product[i]->eType != T_CONSTANT
	  && (!(p->e.p->product[i]->valid & VALID_EVAL_FLAG))) {
	doPolyListSorting (p->e.p->product[i], l);
      }
    polyListAppend (l, p);
    break;

    //If the polynomial is a functionCall, put the parameters in the evaluation list and then
    //put the functionCall in the evaluation list
  case T_FUNCTIONCALL:

    if (p->valid & VALID_EVAL_FLAG)
      break;

    for (i = 0; i < p->e.f->paraNum; i++)
      if (p->e.f->para[i]->eType != T_CONSTANT
	  && (!(p->e.f->para[i]->valid & VALID_EVAL_FLAG))) {
	doPolyListSorting (p->e.f->para[i], l);
      }
    polyListAppend (l, p);
    break;

  default:
    break;
  }
}

void polyListSorting (Polynomial * p, struct polyList *l) {
  polyListSortingCount++;
  doPolyListSorting(p, l);
}

///////////////////////////////////////////////////////////////////
//This function compute the value of a polynomial.  It evaluate  
//all the polynomials in the evaluation list for this polynomial.
//After the polynomials are completed, this is the only function 
//to be executed for evaluation.  This function is not recursive 
///////////////////////////////////////////////////////////////////
void
evaluatePoly (Polynomial * pp, struct polyList *l, double *pReturnValue)
{
  Polynomial *p;
  struct sumPoly *sP;
  struct productPoly *pP;
  register int i, j;
  double pV, re;
  int pE;

  evaluatePolyCount++;
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
    return;
  }

  for (j = 0; j <= l->listNext - 1; j++) {
    p = l->pList[j];
    switch (p->eType) {
    case T_CONSTANT:
      break;
      //Read the value of the variable
    case T_VARIABLE:
      if (p->e.v->vType == 'D')
	p->value = *(p->e.v->vAddr.vAddrD);
      else if (p->e.v->vType == 'I')
	p->value = *(p->e.v->vAddr.vAddrI);
      else {
	fprintf (stderr, "Wrong variable type, exit!\n");
	exit (1);
      }
      break;

      //Sum up all the items in a sum
    case T_SUM:
      p->value = 0;
      sP = p->e.s;
      for (i = 0; i < sP->num; i++) {
	if (sP->factor[i] == 1)
	  p->value += sP->sum[i]->value;
	else
	  p->value += sP->sum[i]->value * sP->factor[i];
      }
      break;

      //Multify all the items
    case T_PRODUCT:
      p->value = 1;
      pP = p->e.p;
      for (i = 0; i < pP->num; i++) {

	pV = pP->product[i]->value;
	pE = pP->exponent[i];

	switch (pE) {
	case 1:
	  re = pV;
	  break;
	case 2:
	  re = pV * pV;
	  break;
	case 3:
	  re = pV * pV * pV;
	  break;
	case 4:
	  re = pV * pV;
	  re = re * re;
	  break;
	case 5:
	  re = pV * pV;
	  re = re * re;
	  re = re * pV;
	  break;
	case 6:
	  re = pV * pV * pV;
	  re = re * re;
	  break;
	default:
	  re = pow (pV, pE);
	  break;
	}
	p->value *= re;
      }
      break;

      //Function calls are evaluated by calling the refered functions.
      //The refered function must be included in the linked library.  
      //Otherwise, the program will exit. 
    case T_FUNCTIONCALL:
      if (strcmp (p->e.f->name, "log10") == 0) {
	p->value = log10 (p->e.f->para[0]->value);
      } else if (strcmp (p->e.f->name, "gsl_ran_tdist_pdf") == 0) {
	p->value =
	  gsl_ran_tdist_pdf (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_tdist_Q") == 0) {
	p->value =
	  gsl_cdf_tdist_Q (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_tdist_P") == 0) {
	p->value =
	  gsl_cdf_tdist_P (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_ran_ugaussian_pdf") == 0) {
	p->value = gsl_ran_ugaussian_pdf (p->e.f->para[0]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_ugaussian_Q") == 0) {
	p->value = gsl_cdf_ugaussian_Q (p->e.f->para[0]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_ugaussian_P") == 0) {
	p->value = gsl_cdf_ugaussian_P (p->e.f->para[0]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_chisq_P") == 0) {
	p->value =
	  gsl_cdf_chisq_P (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_cdf_chisq_Q") == 0) {
	p->value =
	  gsl_cdf_chisq_Q (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "gsl_ran_chisq_pdf") == 0) {
	p->value =
	  gsl_ran_chisq_pdf (p->e.f->para[0]->value, p->e.f->para[1]->value);
      } else if (strcmp (p->e.f->name, "pow") == 0) {
	p->value = pow (p->e.f->para[0]->value, p->e.f->para[1]->value);

      } else if (strcmp (p->e.f->name, "exp") == 0) {
	p->value = exp (p->e.f->para[0]->value);
      } else if (strcmp (p->e.f->name, "sqrt") == 0) {
	p->value = sqrt (p->e.f->para[0]->value);
      } else {
	fprintf (stderr, "unknown function name %s in polynomials\n",
		 p->e.f->name);
	exit (1);
      }
      break;

    default:
      fprintf (stderr, "Error, unknown polynomial type!!!! exit(7)");
      break;
    }
    if (isnan (p->value)) {
      fprintf (stderr,
	       "In evaluatePoly, evaluated value of type %d as not a number\n",
	       p->eType);
      exit (1);
    }
  }
  if (polynomialDebugLevel >= 10)
    fprintf (stderr, "...finished evaluatePoly with %G\n", pp->value);
  *pReturnValue = pp->value;
#ifdef EVALUATESW
  swStop (evaluatePolySW);
#endif
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////
//print out an evaluation list.  This is used for performance evaluation and  debugging
/////////////////////////////////////////////////////////////////////////////////////////
void
printPolyList (struct polyList *l)
{
  Polynomial *p;
  register int i, j;

  for (j = l->listNext - 1; j >= 0; j--) {
    p = l->pList[j];
    fprintf (stderr, "No. %d  id %d", l->listNext - j, p->id);
    switch (p->eType) {
    case T_CONSTANT:
      fprintf (stderr, "  type:CONSTANT   value=%f\n", p->value);
      break;

    case T_VARIABLE:
      fprintf (stderr, "  type:VAR  name=%s  value=%f type=%c",
	       p->e.v->vName, p->value, p->e.v->vType);
      if (p->e.v->vType == 'D')
	fprintf (stderr, "address=%p\n", p->e.v->vAddr.vAddrD);
      else
	fprintf (stderr, "address=%p\n", p->e.v->vAddr.vAddrI);
      break;

    case T_SUM:
      fprintf (stderr, "  type:SUM ");
      for (i = 0; i < p->e.s->num; i++) {
	fprintf (stderr,
		 "       item %d type %d id %d factor %f item value%e  sum value %e \n",
		 i, p->e.s->sum[i]->eType, p->e.s->sum[i]->id,
		 p->e.s->factor[i], p->e.s->sum[i]->value, p->value);
      }
      break;

    case T_PRODUCT:
      fprintf (stderr, "  type:Product");
      for (i = 0; i < p->e.p->num; i++) {
	fprintf (stderr,
		 "       item %d type %d id %d item value %e exponent %d product value %e \n",
		 i, p->e.p->product[i]->eType, p->e.p->product[i]->id,
		 p->e.p->product[i]->value, p->e.p->exponent[i], p->value);
      }
      break;

    case T_FUNCTIONCALL:
      fprintf (stderr, "  type:CONSTANT CALL  name=%s  value=%e\n",
	       p->e.f->name, p->value);
      break;

    default:
      fprintf (stderr, "Error, unknown polynomial type!!!! exit(8)");
      break;
    }
    //    expPrinting(p);
    //    fprintf(stderr,"\n");

  }

}



/////////////////////////////////////////////////////////////////////////////////////////////
//This is the first function to call before we use polynomials.  It allocates some memory  
//and initiates important variables variables                                              
/////////////////////////////////////////////////////////////////////////////////////////////
void
polynomialInitialization ()
{
  int i;
  char *envVar;

#ifdef FREEDEBUG
#warning "freePoly protection is turned on, so memory will purposefully leak!"
  swLogMsg
    ("freePoly protection is turned on, so memory will purposefully leak!");
#endif

  if ((envVar = getenv ("polynomialDebugLevel")) != NULL) {
    polynomialDebugLevel = atoi (envVar);
  }
  if (polynomialDebugLevel > 0)
    fprintf (stderr, "polynomialDebugLevel is at %d\n", polynomialDebugLevel);

  if (polynomialScale <= 0)
    polynomialScale = 1;
  if (polynomialScale > 10)
    polynomialScale = 10;
  fprintf (stderr, "polynomialScale is %d (1-10, 1 is default)\n",
	   polynomialScale);

  evaluatePolySW = swCreate ("evaluatePoly");
  evaluateValueSW = swCreate ("evaluateValue");

  /* Scale all initial and growth sizes up by the polynomial scale (default 10) */

  CONSTANT_HASH_SIZE = MIN_CONSTANT_HASH_SIZE * polynomialScale;
  VARIABLE_HASH_SIZE = MIN_VARIABLE_HASH_SIZE * polynomialScale;
  SUM_HASH_SIZE = MIN_SUM_HASH_SIZE * polynomialScale;
  PRODUCT_HASH_SIZE = MIN_PRODUCT_HASH_SIZE * polynomialScale;
  FUNCTIONCALL_HASH_SIZE = MIN_FUNCTIONCALL_HASH_SIZE * polynomialScale;
  HASH_TABLE_INCREASE = MIN_HASH_TABLE_INCREASE * polynomialScale;
  maxHashLength = HASH_TABLE_INCREASE;	/* Set to actual default. */
  CONSTANT_LIST_INITIAL = MIN_CONSTANT_LIST_INITIAL * polynomialScale;
  CONSTANT_LIST_INCREASE = MIN_CONSTANT_LIST_INCREASE * polynomialScale;
  VARIABLE_LIST_INITIAL = MIN_VARIABLE_LIST_INITIAL * polynomialScale;
  VARIABLE_LIST_INCREASE = MIN_VARIABLE_LIST_INCREASE * polynomialScale;
  SUM_LIST_INITIAL = MIN_SUM_LIST_INITIAL * polynomialScale;
  SUM_LIST_INCREASE = MIN_SUM_LIST_INCREASE * polynomialScale;
  PRODUCT_LIST_INITIAL = MIN_PRODUCT_LIST_INITIAL * polynomialScale;
  PRODUCT_LIST_INCREASE = MIN_PRODUCT_LIST_INCREASE * polynomialScale;
  FUNCTIONCALL_LIST_INITIAL = MIN_FUNCTIONCALL_LIST_INITIAL * polynomialScale;
  FUNCTIONCALL_LIST_INCREASE =
    MIN_FUNCTIONCALL_LIST_INCREASE * polynomialScale;

  if ((envVar = getenv ("polynomialLostNodeId")) != NULL) {
    polynomialLostNodeId = atoi (envVar);
  }
  if (polynomialLostNodeId > 0)
    fprintf (stderr, "polynomialLostNodeId is %d\n", polynomialLostNodeId);

  nodeId = 0;

  // Allocate memory for polynomial list of each type of polynomials, set the counter of each
  //type of polynomials to be 0
  constantListLength = CONSTANT_LIST_INITIAL;
  constantCount = 0;
  constantList =
    (Polynomial **) malloc (constantListLength * sizeof (Polynomial *));
  if (constantList == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }

  variableListLength = 50;
  variableCount = 0;
  variableList =
    (Polynomial **) malloc (variableListLength * sizeof (Polynomial *));
  if (variableList == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }

  sumListLength = SUM_LIST_INITIAL;
  sumCount = 0;
  sumList = (Polynomial **) malloc (sumListLength * sizeof (Polynomial *));
  if (sumList == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }

  productListLength = PRODUCT_LIST_INITIAL;
  productCount = 0;
  productList =
    (Polynomial **) malloc (productListLength * sizeof (Polynomial *));
  if (productList == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }

  functionCallListLength = FUNCTIONCALL_LIST_INITIAL;
  functionCallCount = 0;
  functionCallList =
    (Polynomial **) malloc (functionCallListLength * sizeof (Polynomial *));
  if (functionCallList == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  // Allocate memory for hash tables, each type of polynomials has its own hash table
  constantHash =
    (struct hashStruct *) malloc (CONSTANT_HASH_SIZE *
				  sizeof (struct hashStruct));
  variableHash =
    (struct hashStruct *) malloc (VARIABLE_HASH_SIZE *
				  sizeof (struct hashStruct));
  sumHash =
    (struct hashStruct *) malloc (SUM_HASH_SIZE * sizeof (struct hashStruct));
  productHash =
    (struct hashStruct *) malloc (PRODUCT_HASH_SIZE *
				  sizeof (struct hashStruct));
  functionCallHash =
    (struct hashStruct *) malloc (FUNCTIONCALL_HASH_SIZE *
				  sizeof (struct hashStruct));
  if (constantHash == NULL || variableHash == NULL || sumHash == NULL
      || productHash == NULL || functionCallHash == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  //Initialize the constant hash table, pre-allocate memory for recording polynomials
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < CONSTANT_HASH_SIZE; i++) {
    constantHash[i].num = 0;
    constantHash[i].length = HASH_TABLE_INCREASE;
    constantHash[i].index =
      (int *) malloc (constantHash[i].length * sizeof (int));
    constantHash[i].key =
      (int *) malloc (constantHash[i].length * sizeof (int));
    if (constantHash[i].index == NULL || constantHash[i].key == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }
  }

  //Initialize the variable hash table, pre-allocate memory for recording polynomials
  for (i = 0; i < VARIABLE_HASH_SIZE; i++) {
    variableHash[i].num = 0;
    variableHash[i].length = HASH_TABLE_INCREASE;
    variableHash[i].index =
      (int *) malloc (variableHash[i].length * sizeof (int));
    variableHash[i].key =
      (int *) malloc (variableHash[i].length * sizeof (int));
    if (variableHash[i].index == NULL || variableHash[i].key == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }
  }

  //Initialize the sum hash table, pre-allocate memory for recording polynomials
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < SUM_HASH_SIZE; i++) {
    sumHash[i].num = 0;
    sumHash[i].length = HASH_TABLE_INCREASE;
    sumHash[i].index = (int *) malloc (sumHash[i].length * sizeof (int));
    sumHash[i].key = (int *) malloc (sumHash[i].length * sizeof (int));
    if (sumHash[i].index == NULL || sumHash[i].key == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }

  }

  //Initialize the product hash table, pre-allocate memory for recording polynomials
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < PRODUCT_HASH_SIZE; i++) {
    productHash[i].num = 0;
    productHash[i].length = HASH_TABLE_INCREASE;
    productHash[i].index =
      (int *) malloc (productHash[i].length * sizeof (int));
    productHash[i].key =
      (int *) malloc (productHash[i].length * sizeof (int));
    if (productHash[i].index == NULL || productHash[i].key == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }
  }

  //Initialize the function call hash table, pre-allocate memory for recording polynomials
  for (i = 0; i < FUNCTIONCALL_HASH_SIZE; i++) {
    functionCallHash[i].num = 0;
    functionCallHash[i].length = HASH_TABLE_INCREASE;
    functionCallHash[i].index =
      (int *) malloc (functionCallHash[i].length * sizeof (int));
    functionCallHash[i].key =
      (int *) malloc (functionCallHash[i].length * sizeof (int));
    if (functionCallHash[i].index == NULL || functionCallHash[i].key == NULL) {
      fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	       __LINE__);
      exit (1);
    }

  }

  totalHashMemory =
    (sizeof (struct hashStruct) +
     (2 * HASH_TABLE_INCREASE * sizeof (int))) * (CONSTANT_HASH_SIZE +
						  VARIABLE_HASH_SIZE +
						  SUM_HASH_SIZE +
						  PRODUCT_HASH_SIZE +
						  FUNCTIONCALL_HASH_SIZE);

  //Apply memory for containers to hold terms for a sum polynomial

  //For variable polynomials
  containerLength_v1 = 100;
  factor_v1 = (double *) malloc (containerLength_v1 * sizeof (double));
  p_v1 = (Polynomial **) malloc (containerLength_v1 * sizeof (Polynomial *));
  if (factor_v1 == NULL || p_v1 == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  //For product polynomials
  containerLength_p1 = 100;
  factor_p1 = (double *) malloc (containerLength_p1 * sizeof (double));
  p_p1 = (Polynomial **) malloc (containerLength_p1 * sizeof (Polynomial *));
  if (factor_p1 == NULL || p_p1 == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  //For function call polynomials
  containerLength_f1 = 100;
  factor_f1 = (double *) malloc (containerLength_f1 * sizeof (double));
  p_f1 = (Polynomial **) malloc (containerLength_f1 * sizeof (Polynomial *));
  if (factor_f1 == NULL || p_f1 == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  //Containers for organizing a sum polynomial
  lengthSum = 300;
  factorSum = (double *) malloc (lengthSum * sizeof (double));
  pSum = (Polynomial **) malloc (lengthSum * sizeof (Polynomial *));
  if (factorSum == NULL || pSum == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  //apply memory for container to hold terms for a product polynomial

  //For variable polynomials
  containerLength_v2 = 100;
  exponent_v2 = (int *) malloc (containerLength_v2 * sizeof (int));
  p_v2 = (Polynomial **) malloc (containerLength_v2 * sizeof (Polynomial *));
  if (exponent_v2 == NULL || p_v2 == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  //For sum polynmials
  containerLength_s2 = 100;
  exponent_s2 = (int *) malloc (containerLength_s2 * sizeof (int));
  p_s2 = (Polynomial **) malloc (containerLength_s2 * sizeof (Polynomial *));
  if (exponent_s2 == NULL || p_s2 == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  //For function call polynomials
  containerLength_f2 = 100;
  exponent_f2 = (int *) malloc (containerLength_f2 * sizeof (double));
  p_f2 = (Polynomial **) malloc (containerLength_f2 * sizeof (Polynomial *));
  if (exponent_f2 == NULL || p_f2 == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }
  //Containers for organizing a product polynomial
  lengthProd = 300;
  exponentProd = (int *) malloc (lengthProd * sizeof (int));
  pProd = (Polynomial **) malloc (lengthProd * sizeof (Polynomial *));
  if (exponentProd == NULL || pProd == NULL) {
    fprintf (stderr, "Memory allocation failure at %s line %d\n", __FILE__,
	     __LINE__);
    exit (1);
  }

}

///////////////////////////////////////////////////////////////////////////////
//Release all the memory occupied by the polynomials, polynomial lists, hash 
//tables
///////////////////////////////////////////////////////////////////////////////
void
polynomialClearance ()
{
  int j;

  /* Hah! Exit works better. */
  return;

  nodeId = 0;
  //clear constant polynomials
  for (j = 0; j < constantCount; j++)
    free (constantList[j]);
  constantCount = 0;
  //clear variable polynomials
  for (j = 0; j < variableCount; j++) {
    free (variableList[j]->e.v);
    free (variableList[j]);
  }
  variableCount = 0;
  //clear sum polynomials
  for (j = 0; j < sumCount; j++) {
    free (sumList[j]->e.s->sum);
    free (sumList[j]->e.s->factor);
    free (sumList[j]->e.s);
    free (sumList[j]);
  }
  sumCount = 0;
  //clear product polynomials
  for (j = 0; j < productCount; j++) {
    free (productList[j]->e.p->product);
    free (productList[j]->e.p->exponent);
    free (productList[j]->e.p);
    free (productList[j]);
  }
  productCount = 0;
  //clear function call polynomials
  for (j = 0; j < functionCallCount; j++) {
    free (functionCallList[j]->e.f->name);
    free (functionCallList[j]->e.f->para);
    free (functionCallList[j]->e.f);
    free (functionCallList[j]);
  }
  functionCallCount = 0;
  //release memory occupied by constant hash table
  for (j = 0; j < CONSTANT_HASH_SIZE; j++) {
    if (constantHash[j].length > 0) {
      free (constantHash[j].index);
      free (constantHash[j].key);
    }
  }

  //release memory occupied by variable hash table
  for (j = 0; j < VARIABLE_HASH_SIZE; j++) {
    if (variableHash[j].length > 0) {
      free (variableHash[j].index);
      free (variableHash[j].key);
    }
  }

  //release memory occupied by sum hash table
  for (j = 0; j < SUM_HASH_SIZE; j++) {
    if (sumHash[j].length > 0) {
      free (sumHash[j].index);
      free (sumHash[j].key);
    }
  }

  //release memory occupied by product hash table
  for (j = 0; j < PRODUCT_HASH_SIZE; j++) {
    if (productHash[j].length > 0) {
      free (productHash[j].index);
      free (productHash[j].key);
    }
  }

  //release memory occupied by function call hash table
  for (j = 0; j < FUNCTIONCALL_HASH_SIZE; j++) {
    if (functionCallHash[j].length > 0) {
      free (functionCallHash[j].index);
      free (functionCallHash[j].key);
    }
  }

  //release memory occupied by polynomial lists
  free (constantList);
  free (variableList);
  free (sumList);
  free (productList);
  free (functionCallList);

  //release memory occupied by hash table 
  free (constantHash);
  free (variableHash);
  free (sumHash);
  free (productHash);
  free (functionCallHash);

  //release memory used for construction of sum polynomials
  free (factor_v1);
  free (p_v1);
  free (factor_p1);
  free (p_p1);
  free (factor_f1);
  free (p_f1);
  free (factorSum);
  free (pSum);

  //release memory used for construction of product polynomials
  free (exponent_v2);
  free (p_v2);
  free (exponent_s2);
  free (p_s2);
  free (exponent_f2);
  free (p_f2);
  free (exponentProd);
  free (pProd);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
//This function prints out hash table contents.  This function is for performance evaluation 
///////////////////////////////////////////////////////////////////////////////////////////////////
void
printHashTables ()
{
  int i, j;
  int constantHashCount = 0, variableHashCount = 0, sumHashCount =
    0, productHashCount = 0, functionHashCount = 0;
  int constantHashSize = 0, variableHashSize = 0, sumHashSize =
    0, productHashSize = 0, functionHashSize = 0;

  fprintf (stderr, "Hash table usage for constants:\n");
  for (i = 0; i < CONSTANT_HASH_SIZE; i++) {
    if (constantHash[i].num > 0) {
      fprintf (stderr, "%d  ", constantHash[i].num);
      constantHashCount++;
      constantHashSize += constantHash[i].num * sizeof (int) * 2;
      for (j = 0; j < constantHash[i].num; j++) {
	fprintf (stderr, "(%d %d %d %d=%G) ", j, constantHash[i].index[j],
		 constantHash[i].key[j], constantHash[i].index[j],
		 constantList[constantHash[i].index[j]]->value);
      }
      fprintf (stderr, "\n");
    }
  }
  fprintf (stderr, "\n");

  exit (0);

  fprintf (stderr, "Hash table usage for variables:\n");
  for (i = 0; i < VARIABLE_HASH_SIZE; i++) {
    if (variableHash[i].num > 0) {
      fprintf (stderr, "%d  ", variableHash[i].num);
      variableHashCount++;
      variableHashSize += variableHash[i].num * sizeof (int) * 2;
      for (j = 0; j < variableHash[i].num; j++) {
	fprintf (stderr, "(%d %d %d) ", j, variableHash[i].index[j],
		 variableHash[i].key[j]);
      }
      fprintf (stderr, "\n");
    }
  }
  fprintf (stderr, "\n");

  fprintf (stderr, "Hash table usage for sums:\n");
  for (i = 0; i < SUM_HASH_SIZE; i++) {
    if (sumHash[i].num > 0) {
      fprintf (stderr, "%d  ", sumHash[i].num);
      sumHashCount++;
      sumHashSize += sumHash[i].num * sizeof (int) * 2;
      for (j = 0; j < sumHash[i].num; j++) {
	fprintf (stderr, "(%d %d %d) ", j, sumHash[i].index[j],
		 sumHash[i].key[j]);
      }
      fprintf (stderr, "\n");
    }
  }
  fprintf (stderr, "\n");

  fprintf (stderr, "Hash table usage for products:\n");
  for (i = 0; i < PRODUCT_HASH_SIZE; i++) {
    if (productHash[i].num > 0) {
      fprintf (stderr, "%d  ", productHash[i].num);
      productHashCount++;
      productHashSize += productHash[i].num * sizeof (int) * 2;
      for (j = 0; j < productHash[i].num; j++) {
	fprintf (stderr, "(%d %d %d) ", j, productHash[i].index[j],
		 productHash[i].key[j]);
      }
      fprintf (stderr, "\n");
    }
  }
  fprintf (stderr, "\n");

  fprintf (stderr, "Hash table usage for function calls:\n");
  for (i = 0; i < FUNCTIONCALL_HASH_SIZE; i++) {
    if (functionCallHash[i].num > 0) {
      fprintf (stderr, "%d  ", functionCallHash[i].num);
      functionHashCount++;
      functionHashSize += functionCallHash[i].num * sizeof (int) * 2;
      for (j = 0; j < functionCallHash[i].num; j++) {
	fprintf (stderr, "(%d %d %d) ", j, functionCallHash[i].index[j],
		 functionCallHash[i].key[j]);
      }
      fprintf (stderr, "\n");
    }
  }
  fprintf (stderr, "\n");

  fprintf (stderr,
	   "constantHashCount: %d   variableHashCount: %d   sumHashCount: %d    productHashCount: %d    functionHashCount: %d\n",
	   constantHashCount, variableHashCount, sumHashCount,
	   productHashCount, functionHashCount);
  fprintf (stderr,
	   "constantHashSize: %d   variableHashSize: %d   sumHashSize: %d    productHashSize: %d    functionHashSize: %d\n",
	   constantHashSize, variableHashSize, sumHashSize, productHashSize,
	   functionHashSize);

}

//////////////////////////////////////////////////////////////////////
//print out all variables and their values.  It is used for debugging
//////////////////////////////////////////////////////////////////////
void
printAllVariables ()
{
  int i;

  fprintf (stderr, "All the variables:\n");
  for (i = 0; i < variableCount; i++) {
    fprintf (stderr, "index=%d variable: ", i);
    expPrinting (variableList[i]);

    if (variableList[i]->e.v->vType == 'D') {
      fprintf (stderr, "  value: %f  %f\n",
	       *(variableList[i]->e.v->vAddr.vAddrD), variableList[i]->value);
    } else if (variableList[i]->e.v->vType == 'I') {
      fprintf (stderr, "  value: %d  %f\n",
	       *(variableList[i]->e.v->vAddr.vAddrI), variableList[i]->value);
    }
  }
  fprintf (stderr, "\n");
}

#define MAXPOLYTIERS 16
int polyTiers[MAXPOLYTIERS][5];
int peakPolyTiers;
char *polyTypes[] = { "constant", "variable", "sum", "product", "function" };

/* It might not be a sumPoly, but we can treat it as one. */
void
traversePoly (Polynomial * p, int currentTier)
{
  int i;

  if (p->valid & VALID_EVAL_FLAG)
    return;
  p->valid |= VALID_EVAL_FLAG;

  if (currentTier > peakPolyTiers)
    peakPolyTiers = currentTier;
  if (currentTier >= MAXPOLYTIERS)
    return;
  polyTiers[currentTier][p->eType]++;
  switch (p->eType) {
  case T_CONSTANT:
  case T_VARIABLE:
    break;
  case T_SUM:
    for (i = 0; i < p->e.s->num; i++) {
      traversePoly (p->e.s->sum[i], currentTier + 1);
    }
    break;
  case T_PRODUCT:
    for (i = 0; i < p->e.p->num; i++) {
      traversePoly (p->e.p->product[i], currentTier + 1);
    }
    break;
  case T_FUNCTIONCALL:
    for (i = 0; i < p->e.f->paraNum; i++) {
      traversePoly (p->e.f->para[i], currentTier + 1);
    }
    break;
  case T_FREED:
    fprintf (stderr,
	     "In traversePoly, evil caller is trying to use a polynomial that was freed:\n");
    expTermPrinting (stderr, p, 1);
    exit (1);
    break;
  default:
    fprintf (stderr,
	     "In traversePoly, unknown expression type: [%d], exiting!\n",
	     p->eType);
    exit (1);
  }
}

void
printSummaryPoly (Polynomial * p)
{
  int i, j;

  fprintf (stderr,
	   "Summary of Remains of %d Polynomials (from tree, unique terms only):\n",
	   nodeId);

  clearValidEvalFlag ();
  memset (polyTiers, 0, sizeof (polyTiers));
  peakPolyTiers = 0;
  traversePoly (p, 0);
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

void
expPrinting (Polynomial * p)
{
  int i;

  switch (p->eType) {
  case T_CONSTANT:
    fprintf (stderr, "%G", p->value);
    break;
  case T_VARIABLE:
    fprintf (stderr, "%s", p->e.v->vName);
    break;
  case T_SUM:
    if (p->e.s->num > 1)
      fprintf (stderr, "(");
    if (p->e.s->factor[0] != 1)
      fprintf (stderr, "%G*", p->e.s->factor[0]);
    expPrinting (p->e.s->sum[0]);
    for (i = 1; i < p->e.s->num; i++) {
      fprintf (stderr, "+");
      if (p->e.s->factor[i] != 1)
	fprintf (stderr, "%G*", p->e.s->factor[i]);
      expPrinting (p->e.s->sum[i]);
    }
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
    for (i = 0; i < p->e.f->paraNum - 1; i++) {
      expPrinting (p->e.f->para[i]);
      fprintf (stderr, ", ");
    }
    expPrinting (p->e.f->para[p->e.f->paraNum - 1]);
    fprintf (stderr, ")");
    break;
  case T_FREED:
    fprintf (stderr,
	     "[FREED eType=%d id=%d index=%d key=%d count=%d valid=%d]\n",
	     (int) p->value, p->id, p->index, p->key, p->count, p->valid);
  default:
    fprintf (stderr, "In expPrinting, unknown expression type %d, exiting\n",
	     p->eType);
    exit (1);
  }
}

/* Print the provided polynomial with sum and product terms in
   parentheses and preceeded by their type ('s' or 'p') and index
   down to depth indicated. At that level, simply refer to
   polynomial ids. */

void
expTermPrinting (FILE * output, Polynomial * p, int depth)
{
  int i;

  switch (p->eType) {
  case T_CONSTANT:
    fprintf (output, "%G", p->value);
    break;
  case T_VARIABLE:
    fprintf (output, "%s", p->e.v->vName);
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
    for (i = 1; i < p->e.s->num; i++) {
      fprintf (output, "+");
      if (p->e.s->factor[i] != 1)
	fprintf (output, "%G*", p->e.s->factor[i]);
      expTermPrinting (output, p->e.s->sum[i], depth - 1);
    }
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
    for (i = 0; i < p->e.f->paraNum - 1; i++) {
      expTermPrinting (output, p->e.f->para[i], depth - 1);
      fprintf (output, ", ");
    }
    expTermPrinting (output, p->e.f->para[p->e.f->paraNum - 1], depth - 1);
    fprintf (output, ")");
    break;
  case T_FREED:
    fprintf (stderr,
	     "[FREED eType=%d id=%d index=%d key=%d count=%d valid=%d]\n",
	     (int) p->value, p->id, p->index, p->key, p->count, p->valid);
  default:
    fprintf (stderr,
	     "In expTermPrinting, unknown expression type %d, exiting\n",
	     p->eType);
    raise (SIGUSR1);
    exit (1);
  }
}

void
polyDynamicStatistics (char *title)
{
  unsigned long deltaAccumWallTime, deltaAccumUserTime;

  fprintf (stderr, "Dynamic polynomial statistics (%s):\n", title);

  swDump (overallSW);
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
	   variableCount, variableHashHits, sumCount, sumHashHits,
	   sumHashHits / (float) (sumCount ? sumCount : 1), productCount,
	   productHashHits,
	   productHashHits / (float) (productCount ? productCount : 1),
	   functionCallCount, functionHashHits);

  fprintf (stderr,
	   "List expansions(@size): c=%d(@%d), v=%d(@%d), s=%d(@%d), p=%d(@%d), f=%d(@%d), tcc=%d(@%lu)\n",
	   constantPListExpansions, CONSTANT_LIST_INCREASE,
	   variablePListExpansions, VARIABLE_LIST_INCREASE,
	   sumPListExpansions, SUM_LIST_INCREASE, productPListExpansions,
	   PRODUCT_LIST_INCREASE, functionCallPListExpansions,
	   FUNCTIONCALL_LIST_INCREASE, containerExpansions,
	   (unsigned long) sizeof (Polynomial *) + sizeof (double));

  fprintf (stderr,
	   "List sizes: c=%lu, v=%lu, s=%lu, p=%lu, f=%lu\n",
	   (unsigned long) sizeof (void *) * (constantPListExpansions *
			      CONSTANT_LIST_INCREASE + CONSTANT_LIST_INITIAL),
	   (unsigned long) sizeof (void *) * (variablePListExpansions *
			      VARIABLE_LIST_INCREASE + VARIABLE_LIST_INITIAL),
	   (unsigned long) sizeof (void *) * (sumPListExpansions * SUM_LIST_INCREASE +
			      SUM_LIST_INITIAL),
	   (unsigned long) sizeof (void *) * (productPListExpansions * PRODUCT_LIST_INCREASE +
			      PRODUCT_LIST_INITIAL),
	   (unsigned long) sizeof (void *) * (functionCallPListExpansions *
			      FUNCTIONCALL_LIST_INCREASE +
			      FUNCTIONCALL_LIST_INITIAL));

  fprintf (stderr,
	   "NodeId: %d Hash: max len=%d, size=%lu, SPL: eff=%lu%%, avg len=%lu\n",
	   nodeId, maxHashLength, totalHashMemory,
	   100 * (lowSPLCount +
		  highSPLCount) / (totalSPLCalls ? totalSPLCalls : 1),
	   totalSPLLengths / (totalSPLCalls ? totalSPLCalls : 1));

  fprintf (stderr,
	   "Calls: pLS=%d eP=%d eV=%d hAP=%d kP=%d hP=%d uHP=%d fP=%d(%d) fKP=%d\n",
	   polyListSortingCount, evaluatePolyCount, evaluateValueCount, holdAllPolysCount,
	   keepPolyCount, holdPolyCount, unHoldPolyCount,
	   freePolysAttemptCount, freePolysCount, freeKeptPolysCount);

  if (sumReleaseableCount == 0 && sumNotReleaseableCount == 0 &&
      productReleaseableCount == 0 && productNotReleaseableCount == 0)
    return;

  fprintf (stderr,
	   "Sum polys: try release=%d or not=%d return constant=%d return 1-term=%d hash hits=%d\n",
	   sumReleaseableCount, sumNotReleaseableCount,
	   sumReturnConstantCount, sumReturn1TermCount, sumHashHits);
  fprintf (stderr,
	   "...really new=%d new on sumList=%d replaced on sumList=%d freed=%d 1st-term freed=%d\n",
	   sumNewCount, sumListNewCount, sumListReplacementCount,
	   sumFreedCount, sum1stTermsFreedCount);
  fprintf (stderr,
	   "Product polys: try release=%d or not=%d return 0=%d return constant=%d return 1st term=%d\n",
	   productReleaseableCount, productNotReleaseableCount,
	   productReturn0Count, productReturnConstantCount,
	   productReturn1stTermCount);
  fprintf (stderr,
	   "...actually a sum poly=%d hash hits=%d hash hit but sum=%d return normal (factor is 1)=%d\n",
	   productReturn1TermSumCount, productHashHits,
	   productHashHitIsSumCount, productReturnNormalCount);
  fprintf (stderr,
	   "...factor not 1 now sum=%d new on productList=%d replaced on productList=%d\n",
	   productNon1FactorIsSumCount, productListNewCount,
	   productListReplacementCount);
  fprintf (stderr, "...freed=%d 1st-term freed=%d\n", productFreedCount,
	   product1stTermsFreedCount);

  /* Now we check to see if we're thrashing... If we didn't get at least 25% CPU since the last
     time we did statistics (provided thats at least a second ago), we should bail. */

  if (lastPDSAccumWallTime != 0) {
    deltaAccumWallTime = overallSW->swAccumWallTime - lastPDSAccumWallTime;
    deltaAccumUserTime =
      overallSW->swAccumRU.ru_utime.tv_sec - lastPDSAccumUserTime;
    fprintf (stderr,
	     "Overall user CPU utilization was %lus for last period of %lus, or %lu%%\n",
	     deltaAccumUserTime, deltaAccumWallTime,
	     100 * deltaAccumUserTime / (deltaAccumWallTime ? deltaAccumWallTime : 1));
    if ((deltaAccumUserTime != 0) &&
	(100 * deltaAccumUserTime / (deltaAccumWallTime ? deltaAccumWallTime : 1) < 25)) {
      fprintf (stderr, "We're thrashing, time to abort!\n");
      exit (EXIT_FAILURE);
    }
  }
  lastPDSAccumWallTime = overallSW->swAccumWallTime;
  lastPDSAccumUserTime = overallSW->swAccumRU.ru_utime.tv_sec;

  return;
}

/*
 This function prints out polynomial statistic information.  It is mainly used for 
 performance evaluation and debugging. */
void
polyStatistics (char *title)
{
  long constantSize, variableSize, sumSize, productSize, functionCallSize;
  int sumTerms = 0, productTerms = 0, maxSumTerms = 0, maxProductTerms = 0;
  int constantHashCount = 0, constantHashSize = 0, constantHashPeak = 0,
    variableHashCount = 0, variableHashSize = 0, variableHashPeak = 0,
    sumHashCount = 0, sumHashSize = 0, sumHashPeak = 0,
    productHashCount = 0, productHashSize = 0, productHashPeak = 0,
    functionCallHashCount = 0, functionCallHashSize = 0,
    functionCallHashPeak = 0;
  int i;

  polyDynamicStatistics (title);

  /* VM calls have been known to block, so we're not going to do this in-line
  if (swGetMaximumVMK () != 0) {
    if (swGetCurrentVMK (getpid ()) > (0.95 * swGetMaximumVMK ())) {
      fprintf (stderr,
	       "VM usage too high to permit list traversals for statistics.\n");
      return;
    }
  }
  */

  fprintf (stderr, "Calculated polynomial statistics (%s):\n", title);

  constantSize = constantCount * sizeof (Polynomial);
  variableSize =
    variableCount * (sizeof (Polynomial) + sizeof (struct variablePoly));
  sumSize = sumCount * (sizeof (Polynomial) + sizeof (struct sumPoly));
  for (i = 0; i < sumCount; i++) {
    sumTerms += sumList[i]->e.s->num;
    if (sumList[i]->e.s->num > maxSumTerms)
      maxSumTerms = sumList[i]->e.s->num;
  }
  productSize =
    productCount * (sizeof (Polynomial) + sizeof (struct productPoly));
  for (i = 0; i < productCount; i++) {
    productTerms += productList[i]->e.p->num;
    if (productList[i]->e.p->num > maxProductTerms)
      maxProductTerms = productList[i]->e.p->num;
  }
  functionCallSize = functionCallCount * sizeof (Polynomial);

  fprintf (stderr, "Term count(avg): s=%d(%d), p=%d(%d), ",
	   sumTerms, sumTerms / (sumCount ? sumCount : 1), productTerms,
	   productTerms / (productCount ? productCount : 1));
  fprintf (stderr,
	   "sizes (w/terms): c=%ld, s=%ld, p=%ld, f=%ld\n",
	   constantSize,
	   sumSize + (sumTerms * (sizeof (Polynomial *) + sizeof (double))),
	   productSize +
	   (productTerms * (sizeof (Polynomial *) + sizeof (int))),
	   functionCallSize);

  if (constantCount > 0) {
    for (i = 0; i < CONSTANT_HASH_SIZE; i++) {
      constantHashCount++;
      constantHashSize += constantHash[i].num;
      if (constantHash[i].num > constantHashPeak)
	constantHashPeak = constantHash[i].num;
    }
  }
  if (variableCount > 0) {
    for (i = 0; i < VARIABLE_HASH_SIZE; i++) {
      variableHashCount++;
      variableHashSize += variableHash[i].num;
      if (variableHash[i].num > variableHashPeak)
	variableHashPeak = variableHash[i].num;
    }
  }
  if (sumCount > 0) {
    for (i = 0; i < SUM_HASH_SIZE; i++) {
      sumHashCount++;
      sumHashSize += sumHash[i].num;
      if (sumHash[i].num > sumHashPeak)
	sumHashPeak = sumHash[i].num;
    }
  }
  if (productCount > 0) {
    for (i = 0; i < PRODUCT_HASH_SIZE; i++) {
      productHashCount++;
      productHashSize += productHash[i].num;
      if (productHash[i].num > productHashPeak)
	productHashPeak = productHash[i].num;
    }
  }
  if (functionCallCount > 0) {
    for (i = 0; i < FUNCTIONCALL_HASH_SIZE; i++) {
      functionCallHashCount++;
      functionCallHashSize += functionCallHash[i].num;
      if (functionCallHash[i].num > functionCallHashPeak)
	functionCallHashPeak = functionCallHash[i].num;
    }
  }

  fprintf (stderr, "Hash length peak(avg): c=%d(%d), v=%d(%d), s=%d(%d), ",
	   constantHashPeak,
	   constantHashSize / (constantHashCount ? constantHashCount : 1),
	   variableHashPeak,
	   variableHashSize / (variableHashCount ? variableHashCount : 1),
	   sumHashPeak, sumHashSize / (sumHashCount ? sumHashCount : 1));
  fprintf (stderr, "p=%d(%d), f=%d(%d)\n", productHashPeak,
	   productHashSize / (productHashCount ? productHashCount : 1),
	   functionCallHashPeak,
	   functionCallHashSize /
	   (functionCallHashCount ? functionCallHashCount : 1));

  fprintf (stderr, "---\n");
  return;
};

void
printAllPolynomials ()
{
  int i, j;

  fprintf (stderr, "All Polynomials (from hash):\n");

  /*
  if (constantCount > 0) {
    fprintf (stderr, "All %d constants:\n", constantCount);
    for (i = 0; i < CONSTANT_HASH_SIZE; i++) {
      if (constantHash[i].num <= 0)
	continue;
      for (j = 0; j < constantHash[i].num; j++) {
	fprintf (stderr,
		 "(%d %d) index=%d key=%d count=%d valid=%d constant: ", i, j,
		 constantHash[i].index[j], constantHash[i].key[j],
		 constantList[constantHash[i].index[j]]->count,
		 constantList[constantHash[i].index[j]]->valid);
	expTermPrinting (stderr, constantList[constantHash[i].index[j]], 1);
	fprintf (stderr, "\n");
      }
    }
    fprintf (stderr, "\n");
  }
  */

  if (variableCount > 0) {
    fprintf (stderr, "All %d variables:\n", variableCount);
    for (i = 0; i < VARIABLE_HASH_SIZE; i++) {
      if (variableHash[i].num <= 0)
	continue;
      for (j = 0; j < variableHash[i].num; j++) {
	fprintf (stderr,
		 "(%d %d) index=%d key=%d count=%d valid=%d variable: ", i, j,
		 variableHash[i].index[j], variableHash[i].key[j],
		 variableList[variableHash[i].index[j]]->count,
		 variableList[variableHash[i].index[j]]->valid);
	expTermPrinting (stderr, variableList[variableHash[i].index[j]], 1);
	fprintf (stderr, "\n");
      }
    }
    fprintf (stderr, "\n");
  }
  if (sumCount > 0) {
    fprintf (stderr, "All %d sums:\n", sumCount);
    for (i = 0; i < SUM_HASH_SIZE; i++) {
      if (sumHash[i].num <= 0)
	continue;
      for (j = 0; j < sumHash[i].num; j++) {
	fprintf (stderr, "(%d %d) index=%d key=%d count=%d valid=%d sum: ",
		 i, j, sumHash[i].index[j], sumHash[i].key[j],
		 sumList[sumHash[i].index[j]]->count,
		 sumList[sumHash[i].index[j]]->valid);
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
	fprintf (stderr,
		 "(%d %d) index=%d key=%d count=%d valid=%d product: ", i, j,
		 productHash[i].index[j], productHash[i].key[j],
		 productList[productHash[i].index[j]]->count,
		 productList[productHash[i].index[j]]->valid);
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
	fprintf (stderr,
		 "(%d %d) index=%d key=%d count=%d valid=%d functionCall: ",
		 i, j, functionCallHash[i].index[j],
		 functionCallHash[i].key[j],
		 functionCallList[functionCallHash[i].index[j]]->count,
		 functionCallList[functionCallHash[i].index[j]]->valid);
	expTermPrinting (stderr,
			 functionCallList[functionCallHash[i].index[j]], 1);
	fprintf (stderr, "\n");
      }
    }
    fprintf (stderr, "\n");
  }
  fprintf (stderr, "---\n");
}

/* Mark all polynomials as held. */
void
holdAllPolys ()
{
  int i, j;

  holdAllPolysCount++;

#ifdef _OPENMP
#pragma omp sections
#endif
  {

#ifdef _OPENMP
#pragma omp section
#endif
    {
      fprintf (stderr, "Holding all current polynomials (via hashes):\n");
      if (constantCount > 0) {
	fprintf (stderr, "%d constants\n", constantCount);
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
	fprintf (stderr, "%d variables:\n", variableCount);
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
	fprintf (stderr, "%d sums\n", sumCount);
	for (i = 0; i < SUM_HASH_SIZE; i++) {
	  if (sumHash[i].num <= 0)
	    continue;
	  for (j = 0; j < sumHash[i].num; j++) {
	    if (sumList[sumHash[i].index[j]]->id == polynomialLostNodeId)
	      fprintf (stderr,
		       "holdAllPolys for sums sees id %d and is bumping count\n",
		       polynomialLostNodeId);
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
	fprintf (stderr, "%d products\n", productCount);
	for (i = 0; i < PRODUCT_HASH_SIZE; i++) {
	  if (productHash[i].num <= 0)
	    continue;
	  for (j = 0; j < productHash[i].num; j++) {
	    if (productList[productHash[i].index[j]]->id ==
		polynomialLostNodeId)
	      fprintf (stderr,
		       "holdAllPolys for products sees id %d and is bumping count\n",
		       polynomialLostNodeId);
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
	fprintf (stderr, "%d function calls\n", functionCallCount);
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
  fprintf (stderr, "---\n");
}

/* Flag all components of the provided polynomial for preservation. */
void
doKeepPoly (Polynomial * p)
{
  int i;

  /* If we hit a term with the flag already set, we can stop, because
     all of it's subpolys will have it set as well. This flag is only
     cleared en-masse so there will never be a freeing of a subpoly
     without those it contributes to being freed as well. */
  if (p->valid & VALID_KEEP_FLAG)
    return;

  if (p->id == polynomialLostNodeId)
    fprintf (stderr, "flagValids sees id %d and is flagging with %d\n",
	     polynomialLostNodeId, VALID_KEEP_FLAG);

  p->valid |= VALID_KEEP_FLAG;
  switch (p->eType) {
  case T_CONSTANT:
  case T_VARIABLE:
    break;
  case T_SUM:
  case T_PRODUCT:
  case T_FUNCTIONCALL:
    for (i = 0; i < p->e.s->num; i++) {
      doKeepPoly (p->e.s->sum[i]);
    }
    break;
  case T_FREED:
    fprintf (stderr,
	     "[FREED eType=%d id=%d index=%d key=%d count=%d valid=%d]\n",
	     (int) p->value, p->id, p->index, p->key, p->count, p->valid);
  default:
    fprintf (stderr, "In doKeepPoly, unknown expression type %d, exiting\n",
	     p->eType);
    raise (SIGUSR1);
    exit (1);
  }
  return;
}

void
keepPoly (Polynomial * p)
{
  keepPolyCount++;
  if (polynomialDebugLevel >= 10)
    fprintf (stderr, "Into keepPoly\n");
  doKeepPoly (p);
  if (polynomialDebugLevel >= 10)
    fprintf (stderr, "Out of keepPoly\n");
  return;
}

/* Bump the hold count for this and all component subpolynomials. */
void
doHoldPoly (Polynomial * p)
{
  int i;

  if (p->valid & VALID_EVAL_FLAG)
    return;
  p->valid |= VALID_EVAL_FLAG;

  if (p->id == polynomialLostNodeId)
    fprintf (stderr,
	     "holdPoly sees id %d and is bumping hold count from %d\n",
	     polynomialLostNodeId, p->count);

  switch (p->eType) {
  case T_CONSTANT:
  case T_VARIABLE:
    p->count++;
    break;
  case T_SUM:
  case T_PRODUCT:
  case T_FUNCTIONCALL:
    p->count++;
    for (i = 0; i < p->e.s->num; i++) {
      doHoldPoly (p->e.s->sum[i]);
    }
    break;
  default:
    fprintf (stderr, "In holdPoly, unknown expression type %d, exiting\n",
	     p->eType);
    exit (1);
  }
  return;
}

void
holdPoly (Polynomial * p)
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
void
doUnHoldPoly (Polynomial * p)
{
  int i;

  if (p->valid & VALID_EVAL_FLAG)
    return;
  p->valid |= VALID_EVAL_FLAG;

  if (p->id == polynomialLostNodeId)
    fprintf (stderr,
	     "UnHoldPoly sees id %d and is decrementing hold count from %d\n",
	     polynomialLostNodeId, p->count);

  switch (p->eType) {
  case T_CONSTANT:
  case T_VARIABLE:
    p->count--;
    break;
  case T_SUM:
  case T_PRODUCT:
  case T_FUNCTIONCALL:
    p->count--;
    for (i = 0; i < p->e.s->num; i++) {
      doUnHoldPoly (p->e.s->sum[i]);
    }
    break;
  default:
    fprintf (stderr, "In holdPoly, unknown expression type %d, exiting\n",
	     p->eType);
    exit (1);
  }
  return;
}

void
unHoldPoly (Polynomial * p)
{
  unHoldPolyCount++;
  clearValidEvalFlag ();
  doUnHoldPoly (p);
  return;
}

void
doFreePolys (unsigned short keepMask)
{
  int i, j, k;
  Polynomial **newConstantList, **newVariableList, **newSumList,
    **newProductList;

  if (polynomialDebugLevel >= 10)
    fprintf (stderr, "Starting doFreePolys\n");
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
         we're keeping, and then go thru the old polynomial-specific list flagging
         entries we're not keeping and replacing pointers with new indexes
         for the ones we are keeping. */

      newConstantList =
	(Polynomial **) malloc (sizeof (Polynomial *) * (constantListLength));
      if (newConstantList == NULL) {
	fprintf (stderr, "Memory allocation failure at %s line %d\n",
		 __FILE__, __LINE__);
	exit (1);
      }
      k = 0;
      for (i = 0; i < constantCount; i++) {
	if (constantList[i]->id == polynomialLostNodeId)
	  fprintf (stderr,
		   "doFreePolys sees id %d with valid %d and count %d during pass with mask %d\n",
		   polynomialLostNodeId, constantList[i]->valid,
		   constantList[i]->count, keepMask);
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
	fprintf (stderr, "Free w/mask of %d preserved %d of %d constant",
		 keepMask, k, constantCount);
      constantCount = k;

      /* Go thru the hash collapsing entries and freeing the corresponding
         polynomials and terms. Zero polynomial-specific list entries so they
         can be collapsed next. */
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
	      constantHash[j].index[k] =
		constantList[constantHash[j].index[i]]->index;
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

#ifdef _OPENMP
#pragma omp section
#endif
    {
      newVariableList =
	(Polynomial **) malloc (sizeof (Polynomial *) * (variableListLength));
      if (newVariableList == NULL) {
	fprintf (stderr, "Memory allocation failure at %s line %d\n",
		 __FILE__, __LINE__);
	exit (1);
      }
      k = 0;
      for (i = 0; i < variableCount; i++) {
	if (variableList[i]->id == polynomialLostNodeId)
	  fprintf (stderr,
		   "doFreePolys sees id %d with valid %d and count %d during pass with mask %d\n",
		   polynomialLostNodeId, variableList[i]->valid,
		   variableList[i]->count, keepMask);
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
#else
	  // These are for debugging mis-freed pointers
	  variableList[i]->value = variableList[i]->eType;
	  variableList[i]->eType = T_FREED;
#endif
	  variableList[i] = NULL;
	}
      }
      if (polynomialDebugLevel >= 5)
	fprintf (stderr, "Free w/mask of %d preserved %d of %d variable",
		 keepMask, k, variableCount);
      variableCount = k;

      /* Go thru the hash collapsing entries and freeing the corresponding
         polynomials and terms. Zero polynomial-specific list entries so they
         can be collapsed next. */
      for (j = 0; j < VARIABLE_HASH_SIZE; j++) {
	if (variableHash[j].num > 0) {
	  k = 0;
	  for (i = 0; i < variableHash[j].num; i++) {
	    if (variableList[variableHash[j].index[i]] != NULL) {
	      /* It's a keeper, slide it down and bump the count */
	      //        fprintf(stderr, "Hash keeper index %d is now %d\n", variableHash[j].index[i],
	      //              variableList[variableHash[j].index[i]]->index);
	      variableHash[j].index[k] =
		variableList[variableHash[j].index[i]]->index;
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

#ifdef _OPENMP
#pragma omp section
#endif
    {
      newSumList =
	(Polynomial **) malloc (sizeof (Polynomial *) * (sumListLength));
      if (newSumList == NULL) {
	fprintf (stderr, "Memory allocation failure at %s line %d\n",
		 __FILE__, __LINE__);
	exit (1);
      }
      k = 0;
      for (i = 0; i < sumCount; i++) {
	if (sumList[i]->id == polynomialLostNodeId)
	  fprintf (stderr,
		   "doFreePolys sees id %d with valid %d and count %d during pass with mask %d\n",
		   polynomialLostNodeId, sumList[i]->valid, sumList[i]->count,
		   keepMask);
	if ((sumList[i]->count > 0) || (sumList[i]->valid & keepMask)) {
	  newSumList[k] = sumList[i];
	  sumList[i] = newSumList[k];
	  //      fprintf(stderr, "Index %d is now %d\n", newSumList[k]->index, k);
	  newSumList[k]->index = k;
	  k++;
	} else {
	  sumFreedCount++;
	  free (sumList[i]->e.s->sum);
	  free (sumList[i]->e.s->factor);
	  free (sumList[i]->e.s);
#ifndef FREEDEBUG
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
	fprintf (stderr, "Free w/mask of %d preserved %d of %d sum", keepMask,
		 k, sumCount);
      sumCount = k;

      /* Go thru the hash collapsing entries and freeing the corresponding
         polynomials and terms. Zero polynomial-specific list entries so they
         can be collapsed next. */
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
	fprintf (stderr, "Reducing sumListLength from %d to %d\n",
		 sumListLength, sumCount + SUM_LIST_INITIAL);
	sumListLength = sumCount + SUM_LIST_INITIAL;
	sumList = (Polynomial **) realloc (sumList, sizeof (Polynomial *) *
					   (sumListLength));
      }
    }

#ifdef _OPENMP
#pragma omp section
#endif
    {
      newProductList =
	(Polynomial **) malloc (sizeof (Polynomial *) * (productListLength));
      if (newProductList == NULL) {
	fprintf (stderr, "Memory allocation failure at %s line %d\n",
		 __FILE__, __LINE__);
	exit (1);
      }
      k = 0;
      for (i = 0; i < productCount; i++) {
	if (productList[i]->id == polynomialLostNodeId)
	  fprintf (stderr,
		   "doFreePolys sees id %d with valid %d and count %d during pass with mask %d\n",
		   polynomialLostNodeId, productList[i]->valid,
		   productList[i]->count, keepMask);
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
	fprintf (stderr, " and %d of %d product polynomials\n", k,
		 productCount);
      productCount = k;

      /* Go thru the hash collapsing entries and freeing the corresponding
         polynomials and terms. Zero polynomial-specific list entries so they
         can be collapsed next. */
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
	      productHash[j].index[k] =
		productList[productHash[j].index[i]]->index;
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
	fprintf (stderr, "Reducing productListLength from %d to %d\n",
		 productListLength, productCount + PRODUCT_LIST_INITIAL);
	productListLength = productCount + PRODUCT_LIST_INITIAL;
	productList =
	  (Polynomial **) realloc (productList,
				   sizeof (Polynomial *) *
				   (productListLength));
      }
    }
  }
  /* Reset building statistics. */
  constantHashHits = variableHashHits = functionHashHits = 0;
  sumReleaseableCount = sumNotReleaseableCount = sumReturnConstantCount =
    sumReturn1TermCount = sumHashHits = sumNewCount = sumListNewCount =
    sumListReplacementCount = 0;
  productReleaseableCount = productNotReleaseableCount = productReturn0Count =
    productReturnConstantCount = productReturn1stTermCount =
    productReturn1TermSumCount = productHashHits = productHashHitIsSumCount =
    productReturnNormalCount = productNon1FactorIsSumCount =
    productListNewCount = productListReplacementCount = 0;
  totalSPLLengths = totalSPLCalls = lowSPLCount = highSPLCount = 0;

  if (polynomialDebugLevel >= 10)
    fprintf (stderr, "Reduced list of polynomials (and cleared hits)\n");
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
void
freePolys ()
{
  if (sumNotReleaseableCount + productNotReleaseableCount < 1024 * 512) {
    freePolysAttemptCount++;
    return;
  }
  freePolysCount++;
  doFreePolys (VALID_KEEP_FLAG);
  return;
}

/* Free all polynomials that aren't held. */
void
freeKeptPolys ()
{
  freeKeptPolysCount++;
  keepPolyCount = 0;
  doFreePolys (0);
  return;
}

/* Serialize and compress a polynomial and all subpolys for transfer
   Addresses can be left as-is since we'll create a translation table when we
   reload, and compression will eliminate any temporary inefficiency of size. */
void *
exportPoly (Polynomial * p)
{
  void *exportedPoly = NULL;

  return (exportedPoly);
}

/* Decompress and deserialize a polynomial and all subpolys. Use an address
   translation table to maintain internal consistency. */
Polynomial *
importPoly (void *exportedPoly)
{
  Polynomial *importedPoly = NULL;

  return (importedPoly);
}
