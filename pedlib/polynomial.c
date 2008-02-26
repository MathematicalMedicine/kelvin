#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include "polynomial.h"
//#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"

#include "../../diags/polynomial.c-head"
char *polynomialVersion = "0.0.29";
int polynomialDebugLevel = 0;

///////////////////////////////////////////////////////////////////////////////////////////////
//Recursively evaluate a polynomial.  This version of polynomail evaluation doesn't use
//polynomial sorting list.  It just compute the values of each sub polynomial recursively.  
//A reused sub polynomial maybe repeatedly evaluated.  The efficiency is lower than the
//function of evaluatePoly.  However, we don't need to build a sorting list of sub polynomials
//before the polynomial can be evaluated.
///////////////////////////////////////////////////////////////////////////////////////////////
double
evaluateValue (struct polynomial *p)
{
  int i;
  double result;
  struct sumPoly *sP;
  struct productPoly *pP;
  struct functionPoly *fp;
  double value0, value1;

  switch (p->eType) {
    //If a sub polynomial is a contant, return the value
    //of this constant
  case T_CONSTANT:
    return p->value;
    //If a sub polynomial is a variable, return the value
    //of this variable
  case T_VARIABLE:
    if (p->e.v->vType == 'D') {
      p->value = *(p->e.v->vAddrD);
      return *(p->e.v->vAddrD);
    } else if (p->e.v->vType == 'I') {
      p->value = *(p->e.v->vAddrI);
      return *(p->e.v->vAddrI);
    } else {
      fprintf (stderr, "Wrong variable type, exit!\n");
      exit (1);
    }

    //If a sub polynomial is a sum, evaluate the values of all the terms.  Add the values of
    //the terms weighted by their coefficients.
  case T_SUM:
    result = 0;
    sP = p->e.s;
    for (i = 0; i < sP->num; i++)
      result += evaluateValue (sP->sum[i]) * sP->factor[i];
    p->value = result;
    return result;

    //If a sub polynomial is a product, evaluate the values of all the terms.  Multiply the
    //values of the terms powered by their exponents
  case T_PRODUCT:
    result = 1;
    pP = p->e.p;
    for (i = 0; i < pP->num; i++) {
      result *= pow (evaluateValue (pP->product[i]), pP->exponent[i]);

    }
    p->value = result;
    return result;

    //If a sub polynomial is a function call, evaluate the values of all the parameters
    //and then the function according to the name of the function.  Therefore, a function
    //must have a unique name and a unique set of parameters.
  case T_FUNCTIONCALL:
    fp = p->e.f;
    if (strcmp (fp->name, "log10") == 0)
      result = log10 (evaluateValue (fp->para[0]));
    else if (strcmp (fp->name, "gsl_ran_tdist_pdf") == 0) {
      value0 = evaluateValue (fp->para[0]);
      value1 = evaluateValue (fp->para[1]);
      result = gsl_ran_tdist_pdf (value0, value1);
    } else if (strcmp (fp->name, "gsl_cdf_tdist_Q") == 0) {
      value0 = evaluateValue (fp->para[0]);
      value1 = evaluateValue (fp->para[1]);
      result = gsl_cdf_tdist_Q (value0, value1);
    } else if (strcmp (fp->name, "gsl_cdf_tdist_P") == 0) {
      value0 = evaluateValue (fp->para[0]);
      value1 = evaluateValue (fp->para[1]);
      result = gsl_cdf_tdist_P (value0, value1);
    } else if (strcmp (fp->name, "gsl_ran_ugaussian_pdf") == 0) {
      value0 = evaluateValue (fp->para[0]);
      result = gsl_ran_ugaussian_pdf (value0);
    } else if (strcmp (fp->name, "gsl_cdf_ugaussian_Q") == 0) {
      value0 = evaluateValue (fp->para[0]);
      result = gsl_cdf_ugaussian_Q (value0);
    } else if (strcmp (fp->name, "gsl_cdf_ugaussian_P") == 0) {
      value0 = evaluateValue (fp->para[0]);
      result = gsl_cdf_ugaussian_P (value0);
    } else if (strcmp (fp->name, "gsl_cdf_chisq_P") == 0) {
      value0 = evaluateValue (fp->para[0]);
      value1 = evaluateValue (fp->para[1]);
      result = gsl_cdf_chisq_P (value0, value1);
    } else if (strcmp (fp->name, "gsl_cdf_chisq_Q") == 0) {
      value0 = evaluateValue (fp->para[0]);
      value1 = evaluateValue (fp->para[1]);
      result = gsl_cdf_chisq_Q (value0, value1);
    } else if (strcmp (fp->name, "gsl_ran_chisq_pdf") == 0) {
      value0 = evaluateValue (fp->para[0]);
      value1 = evaluateValue (fp->para[1]);
      result = gsl_ran_chisq_pdf (value0, value1);
    } else if (strcmp (fp->name, "pow") == 0) {
      value0 = evaluateValue (fp->para[0]);
      value1 = evaluateValue (fp->para[1]);
      result = pow (value0, value1);

    } else if (strcmp (fp->name, "exp") == 0) {
      result = exp (evaluateValue (fp->para[0]));
    } else if (strcmp (fp->name, "sqrt") == 0) {
      result = sqrt (fp->para[0]->value);
    } else {
      fprintf (stderr, "unknown function name %s in polynomials\n", fp->name);
      exit (1);
    }
    p->value = result;
    return result;

    //If the polynomial type is unknown, something must be wrong
  default:
    fprintf (stderr, "Error, Unknown expression type!!!!, exit(1)");
    exit (1);
  }
};

/////////////////////////////////////////////////////////////////////////
//Count the numbers of polynomials belonging to different categories   
//in a polynomial sorting list.  This function gives us a quantitative 
//measure of computational complexity of a polynomial                  
/////////////////////////////////////////////////////////////////////////
void
countPoly (struct polyList *l, int *cCounter,
	   int *vCounter, int *sCounter, int *pCounter, int *fCounter)
{
  int j;
  struct polynomial *p;

  if (l->listNext == 0) {
    return;
  }
  for (j = l->listNext - 1; j >= 0; j--) {
    p = l->pList[j];
    switch (p->eType) {
    case T_CONSTANT:
      (*cCounter)++;
    case T_VARIABLE:
      (*vCounter)++;
      break;
    case T_SUM:
      (*sCounter)++;
      break;
    case T_PRODUCT:
      (*pCounter)++;
      break;
    case T_FUNCTIONCALL:
      (*fCounter)++;
      break;
    default:
      fprintf (stderr, "Error, unknown polynomial type!!! exit(2)");
      exit (1);
      break;
    }
  }

};

//////////////////////////////////////////////////////////////////////////////////
//A polynomial is zero if it is a constant polynomial and its value is 0.
//This is used when we build a product polynomial.  If the value of a term is 0,
//then the product will be 0.  We don't need to check other terms and this
//product expression can be build as 0.
/////////////////////////////////////////////////////////////////////////////////
int
isZeroExp (struct polynomial *p)
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
  binaryStart = 0;
  binaryEnd = length - 1;
  while (binaryStart <= binaryEnd) {
    binaryMiddle = floor ((binaryStart + binaryEnd) / 2.0);
    if (target == array[binaryMiddle]) {
      *location = binaryMiddle;
      return 1;
    } else if (target > array[binaryMiddle])
      binaryStart = binaryMiddle + 1;
    else
      binaryEnd = binaryMiddle - 1;
  }				//end of while
  *location = binaryStart;
  return 0;
};

//////////////////////////////////////////////////////////////////////////////////////////
//search for a polynomial in the hash table with a key.  Because more than one polynomials
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
  //If the hash table is full, allocate more memory
  hash->num++;
  if (hash->num > hash->length) {
    hash->length += HASH_TABLE_INCREASE;
    if (hash->length > maxHashListLength)
      maxHashListLength = hash->length;
    hash->key = realloc (hash->key, sizeof (int) * hash->length);
    hash->index = realloc (hash->index, sizeof (int) * hash->length);
    if (hash->key == NULL || hash->index == NULL) {
      fprintf (stderr, "Memory allocation for hash table failed!\n");
      exit (1);
    }
  }
  //Prepare a place in hash table for the new polynomial
  if (location <= hash->num - 2) {
    memmove (&hash->key[location + 1], &hash->key[location],
	     sizeof (int) * (hash->num - location - 1));
    memmove (&hash->index[location + 1], &hash->index[location],
	     sizeof (int) * (hash->num - location - 1));
  }
  //Insert the new polynomiual in the hash table
  hash->key[location] = key;
  hash->index[location] = index;


  //This piece of code is for performance evaluation
  if (hash->num > maxHashLength) {
    maxHash = hash;
    maxHashLength = hash->num;
  }

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

  //Make sure the key is positive and between 0 and MAX_POLYNOMIAL_KEY
  if (key >= 0)
    key = key % MAX_POLYNOMIAL_KEY;
  else
    key = (key + MAX_POLYNOMIAL_KEY) % MAX_POLYNOMIAL_KEY;

  if (key < 0)
    key = 0 - key;
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
struct polynomial *
constantExp (double con)
{
  int i;
  struct polynomial *p;
  int key;
  int hIndex, cIndex;
  int first, last, location;
  int tempI1, tempI2;
  double tempD1, tempD2;

  if (polynomialDebugLevel >= 6)
    fprintf (stderr, "In constantExp constant=%f\n", con);

  tempI1 = tempI2 = -1;


  //Compute a key for the constant  
  tempD1 = frexp (con, &tempI1);
  key = keyConstantPolynomial (tempD1, tempI1);

  //Compute the index of this constant polynomial in the hash table of constant polynomials
  hIndex = key % CONSTANT_HASH_SIZE;

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
	  constantList[cIndex]->count++;
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
  p = (struct polynomial *) malloc (sizeof (struct polynomial));
  if (p == NULL) {
    fprintf (stderr, "Memory allocation error: malloc returned NULL!");
    exit (1);
  }
  p->eType = T_CONSTANT;
  p->value = con;

  //check if the constant polynomial list is full.  Apply for more items if it is full
  if (constantCount >= constantListLength) {
    constantListLength += 1000;
    constantList =
      realloc (constantList,
	       constantListLength * sizeof (struct polynomial *));
  }
  //save the constant in the constant polynomial list
  constantList[constantCount] = p;
  p->index = constantCount;
  p->id = nodeId;
  if (polynomialDebugLevel >= 4)
    fprintf (stderr, "Polynomial %d, (constant %d) added\n", nodeId,
	     constantCount);
  constantCount++;
  nodeId++;
  p->key = key;
  //valid = 1 if the polynomial is used in evaluation
  //valid = 0 if the polynomial is not used in evaluation
  p->valid = 0;

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

  if (vType == 'I')
    key = (long int) vI;
  else if (vType == 'D')
    key = (long int) vD + 1;
  else {
    fprintf (stderr, "UNKNOWN variable type !");
    exit (1);
  }

  //Make sure the key is positive and between 0 and MAX_POLYNOMIAL_KEY
  if (key >= 0)
    key = key % MAX_POLYNOMIAL_KEY;
  else
    key = (key + MAX_POLYNOMIAL_KEY) % MAX_POLYNOMIAL_KEY;

  if (key < 0)
    key = 0 - key;
  return key;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//This function builds a variable polynomial.  The construction of the product and sum polynomials starts
// from calling this function for constructing variable polynomials or indexing variable polynomials.  
// Therefore, the computation in this function must be very efficient.  The efficiency of this function 
//will affect the efficiency of the construction of almost every term in every polynomial
/////////////////////////////////////////////////////////////////////////////////////////////////////////
struct polynomial *
variableExp (double *vD, int *vI, char vType, char name[10])
{
  int i;
  struct polynomial *p;
  struct variablePoly *vPoly;
  int key;
  int hIndex, vIndex;
  int first, last, location;

  //compuete a key for this variable polynomial
  key = keyVariablePolynomial (vD, vI, vType);

  if (polynomialDebugLevel >= 6)
    fprintf (stderr, "In variableExp name=%s\n", name);

  //Compute the index of this variable polynomial in the hash table of variable polynomials
  hIndex = key % VARIABLE_HASH_SIZE;

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
	if ((vType == 'D' && variableList[vIndex]->e.v->vAddrD == vD)
	    || (vType == 'I' && variableList[vIndex]->e.v->vAddrI == vI)) {
	  //if the two variables are the same, return the variable in the variable list 
	  variableList[vIndex]->count++;
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
  p = (struct polynomial *) malloc (sizeof (struct polynomial));
  vPoly = (struct variablePoly *) malloc (sizeof (struct variablePoly));
  if (p == NULL || vPoly == NULL)
    fprintf (stderr, "Memory allocation error: malloc returned NULL!");
  p->eType = T_VARIABLE;
  strcpy (vPoly->vName, name);
  if (vType == 'D')
    vPoly->vAddrD = vD;
  else
    vPoly->vAddrI = vI;
  vPoly->vType = vType;
  p->e.v = vPoly;

  //If the polynomial list is full, apply for more memory
  if (variableCount >= variableListLength) {
    variableListLength += 50;
    variableList =
      realloc (variableList,
	       variableListLength * sizeof (struct polynomial *));
    if (variableList == NULL) {
      fprintf (stderr, "Memory allocation error in variableExp(), exit!");
      exit (1);
    }
  }

  p->index = variableCount;
  p->id = nodeId;
  p->key = key;
  p->valid = 0;

  //Insert the variable polynomial in the variable polynomial list
  variableList[variableCount] = p;
  if (polynomialDebugLevel >= 4)
    fprintf (stderr, "Polynomial %d, (variable %d) added\n", nodeId,
	     variableCount);
  variableCount++;
  nodeId++;

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
keySumPolynomial (struct polynomial **p, double *factor, int counter)
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
    if (key >= 0)
      key = key % MAX_POLYNOMIAL_KEY;
    else
      key = key % MAX_POLYNOMIAL_KEY + MAX_POLYNOMIAL_KEY;
  }

  //Make sure that a key is always positive

  if (key < 0)
    key = 0 - key;

  return key;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//search for a polynomial in a polynomial list.  If the target polynomial is not in the polynomial list,
//return value:  0 not found
//               1 found
//If the target polynomial is found, then its index in the polynomial list is returned in parameter location
////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline int
searchPolynomialList (struct polynomial **p, int length,
		      struct polynomial *target, int *location)
{
  int binaryStart, binaryEnd, binaryMiddle;

  if (length == 0) {
    *location = 0;
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
//as a basis for comparisons between the resulted polynomail and other polynomials 
////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline void
collectSumTerms (double **factor, struct polynomial ***p, int *counter,
		 int *containerLength, double f1, struct polynomial *p1)
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
    *factor =
      (double *) realloc (*factor, (*containerLength) * sizeof (double));
    *p =
      (struct polynomial **) realloc (*p,
				      (*containerLength) *
				      sizeof (struct polynomial *));
    if (*factor == NULL || *p == NULL) {
      fprintf (stderr, "Momery allocation error!\n");
      exit (1);
    }
  }
  //this is a new item in the sum, insert it at the start or end of the sum
  if (*counter == 0 || location >= *counter) {
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
struct polynomial *
plusExp (int num, ...)
{
  int i, k, l;
  va_list args;
  struct sumPoly *sP;
  struct polynomial *rp;
  int counterSum;
  int flag;
  double f1, f0;
  struct polynomial *p1 = 0, *p0 = 0;
  double con = 0;
  int key = 0;
  int tempI, tempI2;
  double tempD, tempD2;
  int hIndex, sIndex;
  int first, last, location;
  int p0SubHIndex = 0, p0HIndex = 0, p0Index = 0, p0Id = 0, p0Key, p0Count;
  enum expressionType p0EType;

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
    p1 = va_arg (args, struct polynomial *);

    if (polynomialDebugLevel >= 6) {
      fprintf (stderr, "In plusExp factor=%f item No. %d of %d type=%d: ", f1,
	       i + 1, num, p1->eType);
      expPrinting (p1);
      fprintf (stderr, "\n");
    }
    //Record the first operand of the plus operation.  Often, the first operand and the result are the
    //same variable, which refers to the case of p = p + ..., therefore, when we build a new polynomial
    //for the result of the plus operation, the polynomial that represents the first operand of the
    //plus operation can be freed.

    if (i == 0) {
      f0 = f1;
      p0 = p1;

//        fprintf(stderr,"f0=%f p0=",f0);
//        expPrinting(p0);
//        fprintf(stderr,"\n");
    }
    //If a term is contant 0, it has no effect on a sum, therefore we do nothing
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
      collectSumTerms (&factor_v1, &p_v1, &counter_v1, &containerLength_v1,
		       f1, p1);
      break;
      //The item is a product
    case T_PRODUCT:
      collectSumTerms (&factor_p1, &p_p1, &counter_p1, &containerLength_p1,
		       f1, p1);
      break;
      //The term is a function call
    case T_FUNCTIONCALL:
      collectSumTerms (&factor_f1, &p_f1, &counter_f1, &containerLength_f1,
		       f1, p1);
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
			   &containerLength_v1, f1 * p1->e.s->factor[l],
			   p1->e.s->sum[l]);
	  break;
	  //product
	case T_PRODUCT:
	  collectSumTerms (&factor_p1, &p_p1, &counter_p1,
			   &containerLength_p1, f1 * p1->e.s->factor[l],
			   p1->e.s->sum[l]);
	  break;
	case T_FUNCTIONCALL:
	  collectSumTerms (&factor_f1, &p_f1, &counter_f1,
			   &containerLength_f1, f1 * p1->e.s->factor[l],
			   p1->e.s->sum[l]);
	  break;
	default:
	  fprintf (stderr, "UNKNOWN polynomial type, exit!!\n");
	  exit (1);
	}
      }
      break;
    default:
      fprintf (stderr, "UNKNOWN polynomial type, exit!!\n");
      exit (1);
    }

    //              for(k=0;k<counter;k++)
    //              {
    //                expPrinting(p[k]);
    //                fprintf(stderr," k=%d index=%d factor=%f\n",k,p[k]->index,factor[k]);
    //              }
  }
  flag = va_arg (args, int);

  va_end (args);

//flag=0;

  if (flag == 0)
    sum0++;
  else
    sum1++;

  counterSum = counter_v1 + counter_p1 + counter_f1;

  if (counterSum + 1 > lengthSum) {
    lengthSum = counterSum + 1;
    factorSum = (double *) realloc (factorSum, lengthSum * sizeof (double));
    pSum =
      (struct polynomial **) realloc (pSum,
				      lengthSum *
				      sizeof (struct polynomial *));
  }
//   j=0;

  if (counter_v1 > 0) {
    memcpy (&factorSum[0], &factor_v1[0], sizeof (double) * counter_v1);
    memcpy (&pSum[0], &p_v1[0], sizeof (struct polynomial *) * counter_v1);
  }
  if (counter_p1 > 0) {
    memcpy (&factorSum[counter_v1], &factor_p1[0],
	    sizeof (double) * counter_p1);
    memcpy (&pSum[counter_v1], &p_p1[0],
	    sizeof (struct polynomial *) * counter_p1);
  }
  if (counter_f1 > 0) {
    memcpy (&factorSum[counter_v1 + counter_p1], &factor_f1[0],
	    sizeof (double) * counter_f1);
    memcpy (&pSum[counter_v1 + counter_p1], &p_f1[0],
	    sizeof (struct polynomial *) * counter_f1);
  }
  //After we go through all the items in the sum,
  //we get only a constant
  if (counterSum == 0) {
    rp = constantExp (con);
    sum2++;
    return rp;
  } else if (con == 0.0 && counterSum == 1 && factorSum[0] == 1.0) {
    rp = pSum[0];
    sum3++;
    return rp;
  } else {
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
	    //                             sumList[sIndex]->count=2;
	    sumList[sIndex]->count++;
	    sum4++;
	    sumHashHits++;
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

  if (flag != 0 && p0->eType == T_SUM && p0->count == 1) {
    p0Index = p0->index;
    p0Id = p0->id;
    p0Key = p0->key;
    p0HIndex = p0->key % SUM_HASH_SIZE;
    p0EType = p0->eType;
    p0Count = p0->count;

    if (sumHash[p0HIndex].num <= 0) {
      fprintf (stderr,
	       "This polynomial is not in the sum hash list, exit(1) !\n");
      exit (1);
    }
    if (searchHashTable (&sumHash[p0HIndex], &first, &last, p0Key)) {
      for (i = first; i <= last; i++) {
	if (p0Index == sumHash[p0HIndex].index[i])
	  break;
      }
      if (i <= last)
	p0SubHIndex = i;
      else {
	fprintf (stderr,
		 "This polynomial is not in the sum hash list, exit(2)!\n");
	exit (1);
      }
    } else {
      fprintf (stderr,
	       "This polynomial is not in the sum hash list, exit(1) !\n");
      exit (1);
    }
    //free the memory of the first polynomial in the parameter list
    free (p0->e.s->sum);
    free (p0->e.s->factor);
    free (p0->e.s);
    free (p0);

  } else {
    p0EType = 999;
    p0Count = 999;
  }

  //If the sum is not found in the sum list, insert it in the sum list
  //Build a new polynomial
  rp = (struct polynomial *) malloc (sizeof (struct polynomial));
  if (rp == NULL) {
    fprintf (stderr, "Memory allocation error: malloc returned NULL!");
    exit (1);
  }
  rp->eType = T_SUM;
  sP = (struct sumPoly *) malloc (sizeof (struct sumPoly));
  if (sP == NULL) {
    fprintf (stderr, "Memory allocation error: malloc returned NULL!");
    exit (1);
  }
  sP->num = counterSum;
  sP->sum =
    (struct polynomial **) malloc (counterSum * sizeof (struct polynomial *));
  sP->factor = (double *) malloc (counterSum * sizeof (double));
  if (sP->sum == NULL || sP->factor == NULL) {
    fprintf (stderr, "Memory allocation error: malloc returned NULL!");
    exit (1);
  }
  for (i = 0; i < sP->num; i++) {
    sP->factor[i] = factorSum[i];
    sP->sum[i] = pSum[i];
  }
  rp->e.s = sP;
  rp->key = key;
  rp->valid = 0;
  rp->count = 1;


  //This piece of code is for performance test
  numSumTerms += counterSum;
  if (counterSum > maxSumLength)
    maxSumLength = counterSum;
  if (counterSum > sizeSumLength) {
    countSumLength =
      realloc (countSumLength, sizeof (int) * (counterSum + 30));
    for (i = sizeSumLength; i < counterSum + 30; i++) {
      countSumLength[i] = 0;
    }
    sizeSumLength = counterSum + 30;
  }
  countSumLength[counterSum - 1]++;



  //Insert the new built polynomial in sum list
  if (sumCount >= sumListLength) {
    sumListLength += 10000;
    sumList = realloc (sumList, sumListLength * sizeof (struct polynomial *));
  }
  if (sumList == NULL) {
    fprintf (stderr, "Memory allocation error in plusExp, exit!");
    exit (1);
  }
//      fprintf(stderr,"sumCount=%d sumListLength=%d key=%d hIndex=%d\n",sumCount, sumListLength, key, hIndex);


  if (flag != 0 && p0EType == T_SUM && p0Count == 1) {
    sum11++;
    sumList[p0Index] = rp;
    sumList[p0Index]->index = p0Index;
    sumList[p0Index]->id = p0Id;
  } else {
    sum00++;
    sumList[sumCount] = rp;
    sumList[sumCount]->index = sumCount;
    sumList[sumCount]->id = nodeId;
    if (polynomialDebugLevel >= 4)
      fprintf (stderr, "Polynomial %d, (sum %d) added\n", nodeId, sumCount);
    sumCount++;
    nodeId++;
  }


//         fprintf(stderr,"Before insert into sumHash\n");
//         for(i=0;i<sumHash[hIndex].num;i++)
//           fprintf(stderr,"hIndex=%d i=%d this key=%d key=%d index=%d\n",
//                           hIndex,i,key,sumHash[hIndex].key[i],
//                           sumHash[hIndex].index[i]);
//         fprintf(stderr,"\n");



  //Insert the newly built polynomial into the Hash table
  if (flag != 0 && p0EType == T_SUM && p0Count == 1) {
    if (p0HIndex != hIndex || p0SubHIndex != location) {
      insertHashTable (&sumHash[hIndex], location, key,
		       sumList[p0Index]->index);
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
    } else			//p0HIndex==hIndex && p0SubHIndex==binaryEnd
    {
      sumHash[hIndex].key[location] = key;
    }
  } else {
    insertHashTable (&sumHash[hIndex], location, key, sumCount - 1);
  }

  sum5++;
  return rp;
};

/////////////////////////////////////////////////////////////////////////////////
//This function compute a key for a product polynomial.  The number of
//product polynomial is usually very big.  Good keys for product polynomials
//are required for good distribution of product polynomials in product hash table
////////////////////////////////////////////////////////////////////////////////
inline int
keyProductPolynomial (struct polynomial **p, int *exponent, int counter)
{
  int key = 0;
  int j;

  for (j = 0; j < counter; j++) {
    key +=
      p[j]->id + p[j]->key * (j + 1) +
      (exponent[j] * 12345 * (int) (p[j]->eType)) % 10000;
    if (key >= 0)
      key = key % MAX_POLYNOMIAL_KEY;
    else
      key = key % MAX_POLYNOMIAL_KEY + MAX_POLYNOMIAL_KEY;

  }

  if (key < 0)
    key = 0 - key;
  return key;

}

/////////////////////////////////////////////////////////////////////////////////////////////
//This function collects the terms for constructing a product polynomial.  The terms are sorted
//so that these terms can have a unique order in the constructed product polynomial.
////////////////////////////////////////////////////////////////////////////////////////////
inline void
collectProductTerms (int **exponent, struct polynomial ***p, int *counter,
		     int *containerLength, int e1, struct polynomial *p1)
{
  int location;

  //search for a position for this component in the container
  if (searchPolynomialList (*p, *counter, p1, &location) == 1) {
    (*exponent)[location] += e1;
    return;
  }

  if ((*counter) >= (*containerLength) - 1) {
    (*containerLength) += 50;
    (*exponent) =
      (int *) realloc ((*exponent), (*containerLength) * sizeof (int));
    (*p) =
      (struct polynomial **) realloc ((*p),
				      (*containerLength) *
				      sizeof (struct polynomial *));
    if ((*exponent) == NULL || (*p) == NULL) {
      fprintf (stderr, "Momery allocation error in timesExp()!\n");
      exit (1);
    }
  }
  //If this item should be in the start or end
  if ((*counter) == 0 || location >= (*counter)) {
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

struct polynomial *
timesExp (int num, ...)
{
  int i, k, l, counterProd;
  va_list args;
  struct productPoly *pP;
  struct polynomial *rp;
  struct polynomial *p1 = 0, *p0 = 0;
  int e1, e0;
  int isZero = 0;
  double factor = 1;
  int key = 0;
  int pIndex, hIndex;
  int first, last, location;
  int flag;
  int p0SubHIndex = 0, p0HIndex = 0, p0Index = 0, p0Id = 0, p0Key, p0Count;
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
    p1 = va_arg (args, struct polynomial *);

    //get the exponent
    e1 = va_arg (args, int);


    if (polynomialDebugLevel >= 6) {
      fprintf (stderr, "In timesExp exponent=%d item No. %d of %d type=%d: ",
	       e1, i + 1, num, p1->eType);
      expPrinting (p1);
      fprintf (stderr, "\n");
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
	    //unknown product type
	  default:
	    fprintf (stderr, "unknown expression type\n");
	    exit (1);
	  }			//end of switch
	}			//end of for
	break;
	//unknown product type
      default:
	fprintf (stderr, "unknown expression type\n");
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
    product0++;
  else
    product1++;


  //The product is zero, a zero polynomial is returned
  if (isZero) {
    rp = constantExp (0.0);
    product2++;
    return rp;
  }
  //combine all the items generated by collecting the operands of the times operation into one polynomial
  counterProd = counter_v2 + counter_s2 + counter_f2;
  if (counterProd > lengthProd) {
    lengthProd = counterProd;
    exponentProd = (int *) realloc (exponentProd, lengthProd * sizeof (int));
    pProd =
      (struct polynomial **) realloc (pProd,
				      lengthProd *
				      sizeof (struct polynomial *));
  }
  if (counter_v2 > 0) {
    memcpy (&exponentProd[0], &exponent_v2[0], sizeof (int) * counter_v2);
    memcpy (&pProd[0], &p_v2[0], sizeof (struct polynomial *) * counter_v2);
  }
  if (counter_s2 > 0) {
    memcpy (&exponentProd[counter_v2], &exponent_s2[0],
	    sizeof (int) * counter_s2);
    memcpy (&pProd[counter_v2], &p_s2[0],
	    sizeof (struct polynomial *) * counter_s2);
  }
  if (counter_f2 > 0) {
    memcpy (&exponentProd[counter_v2 + counter_s2], &exponent_f2[0],
	    sizeof (int) * counter_f2);
    memcpy (&pProd[counter_v2 + counter_s2], &p_f2[0],
	    sizeof (struct polynomial *) * counter_f2);
  }
  //The product has 0 items, the result is a constant polynomial
  if (counterProd == 0) {
    rp = constantExp (factor);
    product3++;
    return rp;
  }
  //If the result polynomial has only one term, it is not a product polynomial
  else if (counterProd == 1 && exponentProd[0] == 1) {
    //If the factor is 1, then the result polynomial is equal to its first term
    if (factor == 1.0) {
      rp = pProd[0];
      product4++;
      return rp;
    }
    //If the factor is not 1, then the result polynomial is a sum polynomial
    else {
      rp = plusExp (1, factor, pProd[0], 0);
      product5++;
      return rp;
    }
  }
  //The result polynomial is a product polynomial
  else {

    //compute the key for the product polynomial
    key = keyProductPolynomial (pProd, exponentProd, counterProd);


    //compute hash table index according to the key
    hIndex = key % PRODUCT_HASH_SIZE;

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
		//                             productList[pIndex]->count=2;
		productList[pIndex]->count++;
		product6++;
		productHashHits++;
		return productList[pIndex];
	      }
	      //this product polynomail has been existing so that we don't need to construct a 
	      //new one.  However, the factor is not 1, therefore the result polynomial is 
	      //a sum polynomial
	      else {
		product7++;
		//                             productList[pIndex]->count=2;
		productList[pIndex]->count++;
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
    if (flag != 0 && p0->eType == T_PRODUCT && p0->count == 1) {
      //We save the identity information of the polynomial to be freed
      //so that the resource can be assigned to the newly created polynomial 
      p0Index = p0->index;
      p0Id = p0->id;
      p0Key = p0->key;
      p0HIndex = p0->key % PRODUCT_HASH_SIZE;
      p0EType = p0->eType;
      p0Count = p0->count;

      //Determine the sub index of the old polynomial in the product hash table
      if (productHash[p0HIndex].num <= 0) {
	fprintf (stderr,
		 "This polynomial is not in the product hash table, exit(1) !\n");
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
		   "This polynomial is not in the product hash list, exit(2)!\n");
	  exit (1);
	}
      } else {
	fprintf (stderr,
		 "This polynomial is not in the product hash list, exit(3) !\n");
	exit (1);
      }

      //Free the first operand
      free (p0->e.p->exponent);
      free (p0->e.p->product);
      free (p0->e.p);
      free (p0);

    } else {
      p0EType = 999;
      p0Count = 999;

    }


    //Construct a new product polynomial from the terms
    //saved in the container
    rp = (struct polynomial *) malloc (sizeof (struct polynomial));
    if (rp == NULL)
      fprintf (stderr, "Memory allocation error: malloc returned NULL!");
    rp->eType = T_PRODUCT;
    pP = (struct productPoly *) malloc (sizeof (struct productPoly));
    if (pP == NULL)
      fprintf (stderr, "Memory allocation error: malloc returned NULL!");
    pP->num = counterProd;
    pP->product =
      (struct polynomial **) malloc (counterProd *
				     sizeof (struct polynomial *));
    pP->exponent = (int *) malloc (counterProd * sizeof (int));
    if (pP->product == NULL || pP->exponent == NULL)
      fprintf (stderr, "Memory allocation error: malloc returned NULL!");
    //copy the collected polynomial terms and exponents into the
    //product polynomial structure
    for (i = 0; i < counterProd; i++) {
      pP->product[i] = pProd[i];
      pP->exponent[i] = exponentProd[i];
    }
    rp->e.p = pP;
    rp->index = productCount;
    rp->id = nodeId;
    rp->key = key;
    rp->valid = 0;
    rp->count = 1;



    //This piece of code is for performance evaluation
    numProductTerms += counterProd;
    if (counterProd > maxProductLength)
      maxProductLength = counterProd;
    if (counterProd > sizeProductLength) {
      countProductLength =
	realloc (countProductLength, sizeof (int) * (counterProd + 30));
      for (i = sizeProductLength; i < counterProd + 30; i++) {
	countProductLength[i] = 0;
      }
      sizeProductLength = counterProd + 30;
    }
    countProductLength[counterProd - 1]++;





    //After the new polynomial is built, it is recorded in the product polynomial list
    if (productCount >= productListLength) {
      productListLength += 10000;
      productList =
	realloc (productList,
		 productListLength * sizeof (struct polynomial *));
    }
    if (productList == NULL) {
      fprintf (stderr, "Memory allocation error in timesExp, exit!");
      exit (1);
    }
    //We either apply for a new position in the product list for this polynomial, or
    //replace an existing polynomial with the newly created polynomial at the position
    //occupied by the existing polynomial
    if (flag != 0 && p0EType == T_PRODUCT && p0Count == 1) {
      //Assign the resource of the freed polynomial to the newly constructed polynomial
      product11++;
      productList[p0Index] = rp;
      productList[p0Index]->index = p0Index;
      productList[p0Index]->id = p0Id;
    } else {
      //save the newly constructed polynomial in the polynomial list
      product00++;
      productList[productCount] = rp;
      productList[productCount]->index = productCount;
      productList[productCount]->id = nodeId;
      if (polynomialDebugLevel >= 4)
	fprintf (stderr, "Polynomial %d, (product %d) added\n", nodeId,
		 productCount);
      productCount++;
      nodeId++;
    }

    //the new polynomial is also recorded in the hash table
    if (flag != 0 && p0EType == T_PRODUCT && p0Count == 1) {

      //If the indexes or sub indexes of the new and old polynomials in the hash table
      //are different
      if (p0HIndex != hIndex || p0SubHIndex != location) {
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
      product8++;
      return rp;
    }
    //If the factor is not 1, return a sum polynomial
    else {
      product9++;
      return plusExp (1, factor, rp, 0);
    }
  }
};

////////////////////////////////////////////////////////////////////////
//This function creates a key for a function call polynomial according
//to the name of the called function and the parameters               
////////////////////////////////////////////////////////////////////////
inline int
keyFunctionCallPolynomial (char *fName, struct polynomial **p, int num)
{
  int key = 0;
  int i;

  for (i = 0; i < strlen (fName); i++)
    key += (int) fName[i];
  for (i = 0; i < num - 1; i++) {
    key = (key + p[i]->key) % FUNCTIONCALL_HASH_SIZE;
  }
  if (key < 0)
    key = 0 - key;

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
struct polynomial *
functionCallExp (int num, ...)
{
  int i, k, hIndex, fIndex;
  char *fName;
  int key;
  struct polynomial **p;
  int first, last, location;
  va_list args;
  struct polynomial *rp;
  struct functionPoly *fP;

  //get the number of parameters for functionCallExp
  va_start (args, num);
  fName = va_arg (args, char *);

  p =
    (struct polynomial **) malloc ((num - 1) * sizeof (struct polynomial *));

  //num is equal to 1 plus the number of parameters for the called function
  for (i = 0; i < num - 1; i++) {
    if (polynomialDebugLevel >= 6) {
      fprintf (stderr, "In functionCallExp item No. %d of %d type=%d: ",
	       i + 1, num, p[i]->eType);
      expPrinting (p[i]);
      fprintf (stderr, "\n");
    }
    p[i] = va_arg (args, struct polynomial *);
  }

  va_end (args);

  //compute a key for this polynomial
  key = keyFunctionCallPolynomial (fName, p, num);

  //compute the index in hash table according to the key

  hIndex = key % FUNCTIONCALL_HASH_SIZE;

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
  rp = (struct polynomial *) malloc (sizeof (struct polynomial));
  rp->eType = T_FUNCTIONCALL;
  fP = (struct functionPoly *) malloc (sizeof (struct functionPoly));
  if (rp == NULL || fP == NULL)
    fprintf (stderr, "Memory allocation error: malloc returned NULL!");
  fP->paraNum = num - 1;
  fP->para =
    (struct polynomial **) malloc ((num - 1) * sizeof (struct polynomial *));
  fP->name = (char *) malloc (strlen (fName) + 1);
  strcpy (fP->name, fName);


  if (fP->para == NULL) {
    fprintf (stderr, "Memory allocation error: malloc returned NULL!");
    exit (1);
  }
  for (i = 0; i < fP->paraNum; i++) {
    fP->para[i] = p[i];
  }
  rp->e.f = fP;
  rp->index = functionCallCount;
  rp->id = nodeId;
  rp->key = key;
  rp->valid = 0;

  //Insert the new built polynomial in function call list
  if (functionCallCount >= functionCallListLength) {
    functionCallListLength += 10000;
    functionCallList =
      realloc (functionCallList,
	       functionCallListLength * sizeof (struct polynomial *));
  }
  if (functionCallList == NULL) {
    fprintf (stderr, "Memory allocation error in functionCallExp, exit!");
    exit (1);
  }
  functionCallList[functionCallCount] = rp;
  if (polynomialDebugLevel >= 4)
    fprintf (stderr, "Polynomial %d, (function %d) added\n", nodeId,
	     functionCallCount);
  functionCallCount++;
  nodeId++;

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
  l = (struct polyList *) malloc (sizeof (struct polyList));
  l->listSize = 100;
  l->listNext = 0;
  l->pList =
    (struct polynomial **) malloc (sizeof (struct polynomial *) *
				   l->listSize);
  return l;
};

//////////////////////////////////////////////////////////////////////////////////////
//This function appends a polynomial onto a polynomial list                         
//////////////////////////////////////////////////////////////////////////////////////
void
polyListAppend (struct polyList *l, struct polynomial *p)
{
  //valid is a mark showing that this polynomial appears on a sorting list
  p->valid = 1;


  if (l->listNext >= l->listSize) {
    l->pList =
      realloc (l->pList, sizeof (struct polynomial *) * (l->listSize + 100));
    l->listSize = l->listSize + 100;
//        fprintf(stderr,"Length of polyList=%d size of polyList=%d\n",
//                l->listNext, sizeof(struct polynomial *)*l->listSize);
  }
  l->pList[l->listNext] = p;
  l->listNext++;
};

/////////////////////////////////////////////////////////////////////////////////////
// Sort the polynomial list to prepare for polynomial evaluation.  Polynomails are 
//evaluated in order.  This function is to determine a order for evaluation.  Also,
//shared terms are evaluated only once.  Together with polyListAppend(), this function
//construct a evaluation list for evaluation of the polynomial.  For linkage computation, 
//it allows shared polynomials within each likelihood polynomial and across likelihood 
//polynomials of a set of pedigrees fully reused
/////////////////////////////////////////////////////////////////////////////////////
void
polyListSorting (struct polynomial *p, struct polyList *l)
{
  int i;

  switch (p->eType) {
    //If the polynomial is a constant, put it in the evaluation list
  case T_CONSTANT:
    polyListAppend (l, p);
    break;

    //If the polynomial is a variable, put it in the evaluation list
  case T_VARIABLE:
    if (p->valid != 1) {
      polyListAppend (l, p);
    }
    break;

    //If the polynomial is a sum, put all the terms of the sum in the evaluation list
    //except constants and then put the sum in the evaluation list
  case T_SUM:

    if (p->valid == 1)
      break;

    for (i = 0; i < p->e.s->num; i++)
      if (p->e.s->sum[i]->eType != T_CONSTANT && p->e.s->sum[i]->valid != 1) {
	polyListSorting (p->e.s->sum[i], l);
      }
    polyListAppend (l, p);
    break;

    //If the polynomial is a product, put all the terms of the product in the 
    //evaluation list except constants and then put the product in the evaluation list
  case T_PRODUCT:

    if (p->valid == 1)
      break;

    for (i = 0; i < p->e.p->num; i++)
      if (p->e.p->product[i]->eType != T_CONSTANT
	  && p->e.p->product[i]->valid != 1) {
	polyListSorting (p->e.p->product[i], l);
      }
    polyListAppend (l, p);
    break;

    //If the polynomial is a functionCall, put the parameters in the evaluation list and then
    //put the functionCall in the evaluation list
  case T_FUNCTIONCALL:

    if (p->valid == 1)
      break;


    for (i = 0; i < p->e.f->paraNum; i++)
      if (p->e.f->para[i]->eType != T_CONSTANT && p->e.f->para[i]->valid != 1) {
	polyListSorting (p->e.f->para[i], l);
      }
    polyListAppend (l, p);
    break;

  default:
    break;
  }

};

/////////////////////////////////////////////////////////////////////////////////////////////////
//This function perform similar task as polyListSorting().  In linkage computation, it constructs
//an evaluation list of the likelihood polynomial of one pedigree which offers no reuse of shared
//polynomials.  Neither shared polynomials within the likelihood polynomial of one pedigree nor 
//those across the likelihood polynomials of a set of pedigrees are reused.   This function is for 
//performance evaluation             
/////////////////////////////////////////////////////////////////////////////////////////////////
void
polyListSorting2 (struct polynomial *p, struct polyList *l)
{
  int i;

  switch (p->eType) {
  case T_CONSTANT:
    polyListAppend (l, p);
    break;

  case T_VARIABLE:
    polyListAppend (l, p);
    break;

  case T_SUM:
    for (i = 0; i < p->e.s->num; i++)
      if (p->e.s->sum[i]->eType != T_CONSTANT) {
	polyListSorting2 (p->e.s->sum[i], l);
      }
    polyListAppend (l, p);
    break;

  case T_PRODUCT:
    for (i = 0; i < p->e.p->num; i++)
      if (p->e.p->product[i]->eType != T_CONSTANT) {
	polyListSorting2 (p->e.p->product[i], l);
      }
    polyListAppend (l, p);
    break;

  case T_FUNCTIONCALL:
    for (i = 0; i < p->e.f->paraNum; i++)
      if (p->e.f->para[i]->eType != T_CONSTANT) {
	polyListSorting (p->e.f->para[i], l);
      }
    polyListAppend (l, p);
    break;

  default:
    break;
  }

};

///////////////////////////////////////////////////////////////////
//This function compute the value of a polynomial.  It evaluate  
//all the polynomials in the evaluation list for this polynomial.
//After the polynomials are completed, this is the only function 
//to be executed for evaluation.  This function is not recursive 
///////////////////////////////////////////////////////////////////
double
evaluatePoly (struct polynomial *pp, struct polyList *l)
{
  struct polynomial *p;
  struct sumPoly *sP;
  struct productPoly *pP;
  register int i, j;
  double pV, re;
  int pE;

  if (l->listNext == 0) {
    return pp->value;
  }

  for (j = 0; j <= l->listNext - 1; j++) {
    p = l->pList[j];
//      expPrinting(p);
//      fprintf(stderr,"\n");
    switch (p->eType) {
//        case T_CONSTANT:
//             break;
      //Read the value of the variable
    case T_VARIABLE:
      if (p->e.v->vType == 'D')
	p->value = *(p->e.v->vAddrD);
      else if (p->e.v->vType == 'I')
	p->value = *(p->e.v->vAddrI);
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
	p->value += sP->sum[i]->value * sP->factor[i];
      }
      break;

      //Multify all the items
    case T_PRODUCT:
      p->value = 1;
      pP = p->e.p;
      for (i = 0; i < pP->num; i++)
//             p->value*=pow(p->e.p->product[i]->value,p->e.p->exponent[i]);
      {

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
	p->value = exp (evaluateValue (p->e.f->para[0]));
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
  }
  return pp->value;
//   return l->pList[l->listNext-1]->value;
}

/////////////////////////////////////////////////////////////////////////////////////////
//print out an evaluation list.  This is used for performance evaluation and  debugging
/////////////////////////////////////////////////////////////////////////////////////////
void
printPolyList (struct polyList *l)
{
  struct polynomial *p;
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
	fprintf (stderr, "address=%ld\n", (long) p->e.v->vAddrD);
      else
	fprintf (stderr, "address=%ld\n", (long) p->e.v->vAddrI);
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

  polynomialDebugLevel = atoi (getenv ("polynomialDebugLevel"));
  if (polynomialDebugLevel > 0)
    fprintf (stderr, "polynomialDebugLevel is at %d\n", polynomialDebugLevel);

  startTime = currentTime = clock ();

  //These variables are for performance evaluation
  maxHashLength = 0;
  sum0 = 0;
  sum1 = 0;
  sum2 = 0;
  sum3 = 0;
  sum4 = 0;
  sum5 = 0;
  sum00 = 0;
  sum11 = 0;
  product0 = 0;
  product1 = 0;
  product2 = 0;
  product3 = 0;
  product4 = 0;
  product5 = 0;
  product6 = 0;
  product7 = 0;
  product8 = 0;
  product9 = 0;
  product00 = 0;
  product11 = 0;
  numSumTerms = 0;
  numProductTerms = 0;
  maxSumLength = maxProductLength = 0;
  countSumLength = countProductLength = NULL;
  sizeSumLength = sizeProductLength = 0;






  //Each polynomial has a unique ID
  nodeId = 0;

  //allocate memory for polynomial list of each type of polynomials, set the counter of each
  //type of polynomials to be 0
  constantListLength = 1000;
  constantCount = 0;
  constantList =
    (struct polynomial **) malloc (constantListLength *
				   sizeof (struct polynomial *));

  variableListLength = 50;
  variableCount = 0;
  variableList =
    (struct polynomial **) malloc (variableListLength *
				   sizeof (struct polynomial *));

  sumListLength = 10000;
  sumCount = 0;
  sumList =
    (struct polynomial **) malloc (sumListLength *
				   sizeof (struct polynomial *));

  productListLength = 10000;
  productCount = 0;
  productList =
    (struct polynomial **) malloc (productListLength *
				   sizeof (struct polynomial *));

  functionCallListLength = 10000;
  functionCallCount = 0;
  functionCallList =
    (struct polynomial **) malloc (functionCallListLength *
				   sizeof (struct polynomial *));

  //allocate memory for hash tables, each type of polynomials has its own hash table
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


  //Initialize the constant hash table, pre-allocate memory for recording polynomials
  for (i = 0; i < CONSTANT_HASH_SIZE; i++) {
    constantHash[i].num = 0;
    constantHash[i].length = HASH_TABLE_INCREASE;
    constantHash[i].index =
      (int *) malloc (constantHash[i].length * sizeof (int));
    constantHash[i].key =
      (int *) malloc (constantHash[i].length * sizeof (int));
    if (constantHash[i].index == NULL || constantHash[i].key == NULL) {
      fprintf (stderr, "Memory allocation for hash table error!");
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
      fprintf (stderr, "Memory allocation for hash table error!");
      exit (1);
    }
  }

  //Initialize the sum hash table, pre-allocate memory for recording polynomials
  for (i = 0; i < SUM_HASH_SIZE; i++) {
    sumHash[i].num = 0;
    sumHash[i].length = HASH_TABLE_INCREASE;
    sumHash[i].index = (int *) malloc (sumHash[i].length * sizeof (int));
    sumHash[i].key = (int *) malloc (sumHash[i].length * sizeof (int));
    if (sumHash[i].index == NULL || sumHash[i].key == NULL) {
      fprintf (stderr, "Memory allocation for hash table error!");
      exit (1);
    }

  }

  //Initialize the product hash table, pre-allocate memory for recording polynomials
  for (i = 0; i < PRODUCT_HASH_SIZE; i++) {
    productHash[i].num = 0;
    productHash[i].length = HASH_TABLE_INCREASE;
    productHash[i].index =
      (int *) malloc (productHash[i].length * sizeof (int));
    productHash[i].key =
      (int *) malloc (productHash[i].length * sizeof (int));
    if (productHash[i].index == NULL || productHash[i].key == NULL) {
      fprintf (stderr, "Memory allocation for hash table error!");
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
      fprintf (stderr, "Memory allocation for hash table error!");
      exit (1);
    }

  }

  //Apply memory for containers to hold terms for a sum polynomial

  //For variable polynomials
  containerLength_v1 = 100;
  factor_v1 = (double *) malloc (containerLength_v1 * sizeof (double));
  p_v1 =
    (struct polynomial **) malloc (containerLength_v1 *
				   sizeof (struct polynomial *));
  if (factor_v1 == NULL || p_v1 == NULL) {
    fprintf (stderr, "Momery allocation error!\n");
    exit (1);
  }
  //For product polynomials
  containerLength_p1 = 100;
  factor_p1 = (double *) malloc (containerLength_p1 * sizeof (double));
  p_p1 =
    (struct polynomial **) malloc (containerLength_p1 *
				   sizeof (struct polynomial *));
  if (factor_p1 == NULL || p_p1 == NULL) {
    fprintf (stderr, "Momery allocation error!\n");
    exit (1);
  }
  //For function call polynomials
  containerLength_f1 = 100;
  factor_f1 = (double *) malloc (containerLength_f1 * sizeof (double));
  p_f1 =
    (struct polynomial **) malloc (containerLength_f1 *
				   sizeof (struct polynomial *));
  if (factor_f1 == NULL || p_f1 == NULL) {
    fprintf (stderr, "Momery allocation error!\n");
    exit (1);
  }
  //Containers for organizing a sum polynomial
  lengthSum = 300;
  factorSum = (double *) malloc (lengthSum * sizeof (double));
  pSum =
    (struct polynomial **) malloc (lengthSum * sizeof (struct polynomial *));


  //apply memory for container to hold terms for a product polynomial

  //For variable polynomials
  containerLength_v2 = 100;
  exponent_v2 = (int *) malloc (containerLength_v2 * sizeof (int));
  p_v2 =
    (struct polynomial **) malloc (containerLength_v2 *
				   sizeof (struct polynomial *));
  if (exponent_v2 == NULL || p_v2 == NULL) {
    fprintf (stderr, "Momery allocation error!\n");
    exit (1);
  }
  //For sum polynmials
  containerLength_s2 = 100;
  exponent_s2 = (int *) malloc (containerLength_s2 * sizeof (int));
  p_s2 =
    (struct polynomial **) malloc (containerLength_s2 *
				   sizeof (struct polynomial *));
  if (exponent_s2 == NULL || p_s2 == NULL) {
    fprintf (stderr, "Momery allocation error!\n");
    exit (1);
  }
  //For function call polynomials
  containerLength_f2 = 100;
  exponent_f2 = (int *) malloc (containerLength_f2 * sizeof (double));
  p_f2 =
    (struct polynomial **) malloc (containerLength_f2 *
				   sizeof (struct polynomial *));
  if (exponent_f2 == NULL || p_f2 == NULL) {
    fprintf (stderr, "Momery allocation error!\n");
    exit (1);
  }
  //Containers for organizing a product polynomial
  lengthProd = 300;
  exponentProd = (int *) malloc (lengthProd * sizeof (int));
  pProd =
    (struct polynomial **) malloc (lengthProd * sizeof (struct polynomial *));

}

/////////////////////////////////////////////////////////////////////////////////
//record the current status of the polynomials
/////////////////////////////////////////////////////////////////////////////////
void
makePolynomialStamp ()
{

//           nodeIdStamp            = nodeId;
  constantCountStamp = constantCount;
  variableCountStamp = variableCount;
  sumCountStamp = sumCount;
  productCountStamp = productCount;
  functionCallCountStamp = functionCallCount;

//           printAllPolynomials();  
}

///////////////////////////////////////////////////////////////////////////////
//record the current status of the polynomials
///////////////////////////////////////////////////////////////////////////////
void
makePolynomialStamp2 ()
{

//           nodeIdStamp2            = nodeId;
  constantCountStamp2 = constantCount;
  variableCountStamp2 = variableCount;
  sumCountStamp2 = sumCount;
  productCountStamp2 = productCount;
  functionCallCountStamp2 = functionCallCount;
//           printAllPolynomials();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//This function clears the polynomials that are built after a time stamp.  The time stamp recorded 
//the status of polynomials when the time stamp is made.  In genetic likelihood computations, this 
//function is used to clear the polynomials built for the likelihood polynomials of all pedigrees  
//at a specific trait position so that we may have enough memory for construction of likelihood    
//polynomials at next trait position                                                               
/////////////////////////////////////////////////////////////////////////////////////////////////////
void
partialPolynomialClearance ()
{
  int i, j, k;

//   nodeId=nodeIdStamp;

  //partially clear constant polynomials
  for (j = constantCountStamp; j < constantCount; j++)
    free (constantList[j]);
  constantCount = constantCountStamp;

  //partially clear variable polynomials
  for (j = variableCountStamp; j < variableCount; j++) {
    free (variableList[j]->e.v);
    free (variableList[j]);
  }

  variableCount = variableCountStamp;
  //partially clear sum polynomials
  for (j = sumCountStamp; j < sumCount; j++) {
    free (sumList[j]->e.s->sum);
    free (sumList[j]->e.s->factor);
    free (sumList[j]->e.s);
    free (sumList[j]);
  }
  sumCount = sumCountStamp;

  //partially clear product polynomials
  for (j = productCountStamp; j < productCount; j++) {
    free (productList[j]->e.p->product);
    free (productList[j]->e.p->exponent);
    free (productList[j]->e.p);
    free (productList[j]);
  }
  productCount = productCountStamp;

  //partially clear function call polynomials
  for (j = functionCallCountStamp; j < functionCallCount; j++) {
    free (functionCallList[j]->e.f->name);
    free (functionCallList[j]->e.f->para);
    free (functionCallList[j]->e.f);
    free (functionCallList[j]);
  }
  functionCallCount = functionCallCountStamp;

  //adjust the memory for constants
  constantListLength = constantCount + 10;
  constantList =
    realloc (constantList, constantListLength * sizeof (struct polynomial *));

  //adjust the memory for variables
  variableListLength = variableCount + 10;
  variableList =
    realloc (variableList, variableListLength * sizeof (struct polynomial *));

  //adjust the memory for sums
  sumListLength = sumCount + 10;
  sumList = realloc (sumList, sumListLength * sizeof (struct polynomial *));

  //adjust the memory for prodcuts
  productListLength = productCount + 10;
  productList =
    realloc (productList, productListLength * sizeof (struct polynomial *));

  //adjust the memory for functioncalls
  functionCallListLength = functionCallCount + 10;
  functionCallList =
    realloc (functionCallList,
	     functionCallListLength * sizeof (struct polynomial *));

  //rearrange the hash table of constants, delete the entries of polynomials in the hash table
  //that are built after the time stamp
  for (j = 0; j < CONSTANT_HASH_SIZE; j++) {
    if (constantHash[j].num > 0) {
      k = 0;
      for (i = 0; i < constantHash[j].num; i++) {
	if (constantHash[j].index[i] < constantCountStamp) {
	  constantHash[j].index[k] = constantHash[j].index[i];
	  constantHash[j].key[k] = constantHash[j].key[i];
	  k++;
	}
      }
      constantHash[j].num = k;
      constantHash[j].length = constantHash[j].num + HASH_TABLE_INCREASE;
      constantHash[j].key = realloc (constantHash[j].key,
				     sizeof (int) * constantHash[j].length);
      constantHash[j].index = realloc (constantHash[j].index,
				       sizeof (int) * constantHash[j].length);
      if (constantHash[j].key == NULL || constantHash[j].index == NULL) {
	fprintf (stderr, "memory allocation error, exit\n");
	exit (1);
      }

    }
  }


  //rearrange the hash table of variables
  for (j = 0; j < VARIABLE_HASH_SIZE; j++) {
    if (variableHash[j].num > 0) {
      k = 0;
      for (i = 0; i < variableHash[j].num; i++) {
	if (variableHash[j].index[i] < variableCountStamp) {
	  variableHash[j].index[k] = variableHash[j].index[i];
	  variableHash[j].key[k] = variableHash[j].key[i];
	  k++;
	}
      }
      variableHash[j].num = k;
      variableHash[j].length = variableHash[j].num + HASH_TABLE_INCREASE;
      variableHash[j].key = realloc (variableHash[j].key,
				     sizeof (int) * variableHash[j].length);
      variableHash[j].index = realloc (variableHash[j].index,
				       sizeof (int) * variableHash[j].length);
      if (variableHash[j].key == NULL || variableHash[j].index == NULL) {
	fprintf (stderr, "memory allocation error, exit\n");
	exit (1);
      }


    }
  }

  //rearrange hash table of sums
  for (j = 0; j < SUM_HASH_SIZE; j++) {
    if (sumHash[j].num > 0) {
      k = 0;
      for (i = 0; i < sumHash[j].num; i++) {
	if (sumHash[j].index[i] < sumCountStamp) {
	  sumHash[j].index[k] = sumHash[j].index[i];
	  sumHash[j].key[k] = sumHash[j].key[i];
	  k++;
	}
      }
      sumHash[j].num = k;
      sumHash[j].length = sumHash[j].num + HASH_TABLE_INCREASE;
      sumHash[j].key = realloc (sumHash[j].key,
				sizeof (int) * sumHash[j].length);
      sumHash[j].index = realloc (sumHash[j].index,
				  sizeof (int) * sumHash[j].length);
      if (sumHash[j].key == NULL || sumHash[j].index == NULL) {
	fprintf (stderr, "memory allocation error, exit\n");
	exit (1);
      }
    }
  }
  //rearrange hash table of products
  for (j = 0; j < PRODUCT_HASH_SIZE; j++) {
    if (productHash[j].num > 0) {
      k = 0;
      for (i = 0; i < productHash[j].num; i++) {
	if (productHash[j].index[i] < productCountStamp) {
	  productHash[j].index[k] = productHash[j].index[i];
	  productHash[j].key[k] = productHash[j].key[i];
	  k++;
	}
      }
      productHash[j].num = k;
      productHash[j].length = productHash[j].num + HASH_TABLE_INCREASE;
      productHash[j].key = realloc (productHash[j].key,
				    sizeof (int) * productHash[j].length);
      productHash[j].index = realloc (productHash[j].index,
				      sizeof (int) * productHash[j].length);
      if (productHash[j].key == NULL || productHash[j].index == NULL) {
	fprintf (stderr, "memory allocation error, exit\n");
	exit (1);
      }

    }
  }
  //rearrange hash table of function calls
  for (j = 0; j < FUNCTIONCALL_HASH_SIZE; j++) {
    if (functionCallHash[j].num > 0) {
      k = 0;
      for (i = 0; i < functionCallHash[j].num; i++) {
	if (functionCallHash[j].index[i] < functionCallCountStamp) {
	  functionCallHash[j].index[k] = functionCallHash[j].index[i];
	  functionCallHash[j].key[k] = functionCallHash[j].key[i];
	  k++;
	}
      }
      functionCallHash[j].num = k;
      functionCallHash[j].length =
	functionCallHash[j].num + HASH_TABLE_INCREASE;
      functionCallHash[j].key =
	realloc (functionCallHash[j].key,
		 sizeof (int) * functionCallHash[j].length);
      functionCallHash[j].index =
	realloc (functionCallHash[j].index,
		 sizeof (int) * functionCallHash[j].length);
      if (functionCallHash[j].key == NULL
	  || functionCallHash[j].index == NULL) {
	fprintf (stderr, "memory allocation error, exit\n");
	exit (1);
      }
    }
  }
  //clear the signs for evaluation referenceof the polynomials left after partial clearance 
  for (j = 0; j < constantCount; j++)
    constantList[j]->valid = 0;
  for (j = 0; j < variableCount; j++)
    variableList[j]->valid = 0;
  for (j = 0; j < sumCount; j++)
    sumList[j]->valid = 0;
  for (j = 0; j < productCount; j++)
    productList[j]->valid = 0;
  for (j = 0; j < functionCallCount; j++)
    functionCallList[j]->valid = 0;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This function is called for partial clearance of polynomials built after another time stamp.  In Kelvin,    
//it is used to clear the polynomialls built during the construction of the likelihood of one pedigree, while 
//the last function partialPolynomialClearance() used to clear the polynomials built during the construction  
//of the likelihood polynomials of all the pedigrees at the current trait position.                           
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
partialPolynomialClearance2 ()
{
  int i, j, k;
  int *sumIndex;
  int *productIndex;
  int *functionCallIndex;
  int sumRecount;
  int productRecount;
  int functionCallRecount;

  //
  sumIndex = (int *) malloc (sizeof (int) * (sumCount - sumCountStamp2 + 1));
  if (sumIndex == NULL) {
    fprintf (stderr, "sumIndex memory application failed!\n");
    exit (1);
  }

  productIndex =
    (int *) malloc (sizeof (int) * (productCount - productCountStamp2 + 1));
  if (productIndex == NULL) {
    fprintf (stderr, "productIndex memory application failed!\n");
    exit (1);
  }

  functionCallIndex =
    (int *) malloc (sizeof (int) *
		    (functionCallCount - functionCallCountStamp2 + 1));
  if (functionCallIndex == NULL) {
    fprintf (stderr, "functionCallIndex memory application failed!\n");
    exit (1);
  }
  //partially clear sums
  sumRecount = sumCountStamp2;
  for (j = sumCountStamp2; j < sumCount; j++) {
    if (sumList[j]->valid >= 1) {
      sumList[sumRecount] = sumList[j];
      sumList[sumRecount]->index = sumRecount;
      sumIndex[j - sumCountStamp2] = sumRecount;
      sumRecount++;
    } else {
      sumIndex[j - sumCountStamp2] = -1;
      free (sumList[j]->e.s->sum);
      free (sumList[j]->e.s->factor);
      free (sumList[j]->e.s);
      free (sumList[j]);
    }
  }
  sumCount = sumRecount;

//   fprintf(stderr,"sumList clearance completed!\n");

  //partially clear products 
  productRecount = productCountStamp2;
  for (j = productCountStamp2; j < productCount; j++) {
    if (productList[j]->valid >= 1) {
      productList[productRecount] = productList[j];
      productList[productRecount]->index = productRecount;
      productIndex[j - productCountStamp2] = productRecount;
      productRecount++;
    } else {
      productIndex[j - productCountStamp2] = -1;
      free (productList[j]->e.p->product);
      free (productList[j]->e.p->exponent);
      free (productList[j]->e.p);
      free (productList[j]);
    }
  }
  productCount = productRecount;

//   fprintf(stderr,"productList clearance completed!\n");

  //partially clear function-calls
  functionCallRecount = functionCallCountStamp2;
  for (j = functionCallCountStamp2; j < functionCallCount; j++) {
    if (functionCallList[j]->valid >= 1) {
      functionCallList[functionCallRecount] = functionCallList[j];
      functionCallList[functionCallRecount]->index = functionCallRecount;
      functionCallIndex[j - functionCallCountStamp2] = functionCallRecount;
      functionCallRecount++;
    } else {
      functionCallIndex[j - functionCallCountStamp2] = -1;
      free (functionCallList[j]->e.f->name);
      free (functionCallList[j]->e.f->para);
      free (functionCallList[j]->e.f);
      free (functionCallList[j]);
    }
  }
  functionCallCount = functionCallRecount;

//   fprintf(stderr,"function List clearance completed!\n");

  //rearrange hash table of sums
  for (j = 0; j < SUM_HASH_SIZE; j++) {
    if (sumHash[j].num > 0) {
      k = 0;
      for (i = 0; i < sumHash[j].num; i++) {
	if (sumHash[j].index[i] < sumCountStamp2) {
	  sumHash[j].index[k] = sumHash[j].index[i];
	  sumHash[j].key[k] = sumHash[j].key[i];
	  k++;
	} else if (sumIndex[sumHash[j].index[i] - sumCountStamp2] != -1) {
	  sumHash[j].index[k] =
	    sumIndex[sumHash[j].index[i] - sumCountStamp2];
	  sumHash[j].key[k] = sumHash[j].key[i];
	  k++;
	}
      }
      sumHash[j].num = k;
      sumHash[j].length = sumHash[j].num;
      sumHash[j].key = realloc (sumHash[j].key,
				sizeof (int) * sumHash[j].length);
      sumHash[j].index = realloc (sumHash[j].index,
				  sizeof (int) * sumHash[j].length);
    }
  }

  //  fprintf(stderr,"sumHash rearrangement completed!\n");

  //rearrange hash table of products
  for (j = 0; j < PRODUCT_HASH_SIZE; j++) {
    if (productHash[j].num > 0) {
      k = 0;
      for (i = 0; i < productHash[j].num; i++) {
	if (productHash[j].index[i] < productCountStamp2) {
	  productHash[j].index[k] = productHash[j].index[i];
	  productHash[j].key[k] = productHash[j].key[i];
	  k++;
	} else if (productIndex[productHash[j].index[i] - productCountStamp2]
		   != -1) {
	  productHash[j].index[k] =
	    productIndex[productHash[j].index[i] - productCountStamp2];
	  productHash[j].key[k] = productHash[j].key[i];
	  k++;
	}
      }
      productHash[j].num = k;
      productHash[j].length = productHash[j].num;
      productHash[j].key = realloc (productHash[j].key,
				    sizeof (int) * productHash[j].length);
      productHash[j].index = realloc (productHash[j].index,
				      sizeof (int) * productHash[j].length);
    }
  }

//   fprintf(stderr,"productHash rearrangement completed!\n");


  //rearrange hash table of function calls
  for (j = 0; j < FUNCTIONCALL_HASH_SIZE; j++) {
    if (functionCallHash[j].num > 0) {
      k = 0;
      for (i = 0; i < functionCallHash[j].num; i++) {
	if (functionCallHash[j].index[i] < functionCallCountStamp2) {
	  functionCallHash[j].index[k] = functionCallHash[j].index[i];
	  functionCallHash[j].key[k] = functionCallHash[j].key[i];
	  k++;
	} else
	  if (functionCallIndex
	      [functionCallHash[j].index[i] - functionCallCountStamp2] !=
	      -1) {
	  functionCallHash[j].index[k] =
	    functionCallIndex[functionCallHash[j].index[i] -
			      functionCallCountStamp2];
	  functionCallHash[j].key[k] = functionCallHash[j].key[i];
	  k++;
	}
      }
      functionCallHash[j].num = k;
      functionCallHash[j].length = functionCallHash[j].num;
      functionCallHash[j].key = realloc (functionCallHash[j].key,
					 sizeof (int) *
					 functionCallHash[j].length);
      functionCallHash[j].index =
	realloc (functionCallHash[j].index,
		 sizeof (int) * functionCallHash[j].length);
    }
  }

//   fprintf(stderr,"functionCallHash rearrangement completed!\n");

  free (sumIndex);
  free (productIndex);
  free (functionCallIndex);

};

///////////////////////////////////////////////////////////////////////////////
//Release all the memory occupied by the polynomials, polynomial lists, hash 
//tables
///////////////////////////////////////////////////////////////////////////////
void
polynomialClearance ()
{
  int j;

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
    if (constantHash[j].num > 0) {
      free (constantHash[j].index);
      free (constantHash[j].key);
    }
  }

  //release memory occupied by variable hash table
  for (j = 0; j < VARIABLE_HASH_SIZE; j++) {
    if (variableHash[j].num > 0) {
      free (variableHash[j].index);
      free (variableHash[j].key);
    }
  }

  //release memory occupied by sum hash table
  for (j = 0; j < SUM_HASH_SIZE; j++) {
    if (sumHash[j].num > 0) {
      free (sumHash[j].index);
      free (sumHash[j].key);
    }
  }

  //release memory occupied by product hash table
  for (j = 0; j < PRODUCT_HASH_SIZE; j++) {
    if (productHash[j].num > 0) {
      free (productHash[j].index);
      free (productHash[j].key);
    }
  }

  //release memory occupied by function call hash table
  for (j = 0; j < FUNCTIONCALL_HASH_SIZE; j++) {
    if (functionCallHash[j].num > 0) {
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


//////////////////////////////////////////////////////////////////////
//This function display the structure of all the polynomials.  It is
//used for debugging
//////////////////////////////////////////////////////////////////////
void
dismantle ()
{
  int i;

  for (i = 0; i < constantCount; i++)
    fprintf (stderr, "Constant index=%d  value=%10.8f key=%d\n",
	     constantList[i]->index, constantList[i]->value,
	     constantList[i]->key);
  for (i = 0; i < variableCount; i++)
    fprintf (stderr, "Variable index=%d  value=%10.8f key=%d\n",
	     variableList[i]->index, variableList[i]->value,
	     variableList[i]->key);
  for (i = 0; i < sumCount; i++) {
    fprintf (stderr, "Sum      index=%d  value=%10.8f key=%d\n",
	     sumList[i]->index, sumList[i]->value, sumList[i]->key);

/*
     for(j=0;j<sumList[i]->e.s->num;j++)
     {
        fprintf(stderr,"   %d  %ld  %10.8f  %10.8f\n",
              j,sumList[i]->e.s->sum[j]->index,sumList[i]->e.s->sum[j]->value, sumList[i]->e.s->factor[j]);
     }
*/
  }
  for (i = 0; i < productCount; i++) {
    fprintf (stderr, "Product  index=%d  value=%10.8f key=%d\n",
	     productList[i]->index, productList[i]->value,
	     productList[i]->key);

/*
     for(j=0;j<productList[i]->e.p->num;j++)
     {
        fprintf(stderr,"   %d  %ld  %10.8f  %d\n",
                j,productList[i]->e.p->product[j]->index,
                productList[i]->e.p->product[j]->value, 
                productList[i]->e.p->exponent[j]);
     }
*/
  }
  for (i = 0; i < functionCallCount; i++) {
    fprintf (stderr, "Function Call  id=%d  name: %s  value=%10.8f key=%d\n",
	     functionCallList[i]->index, functionCallList[i]->e.f->name,
	     functionCallList[i]->value, functionCallList[i]->key);

/*
     for(j=0;j<functionCallList[i]->e.f->paraNum;j++)
     {
        fprintf(stderr,"   %d  %ld  %10.8f \n",
                j,functionCallList[i]->e.f->para[j]->index,
                functionCallList[i]->e.f->para[j]->value);
     }
*/
  }
};

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
	fprintf (stderr, "(%d %d %d) ", j, constantHash[i].index[j],
		 constantHash[i].key[j]);
      }
      fprintf (stderr, "\n");
    }
  }
  fprintf (stderr, "\n");

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
      fprintf (stderr, "  value: %f  %f\n", *(variableList[i]->e.v->vAddrD),
	       variableList[i]->value);
    } else if (variableList[i]->e.v->vType == 'I') {
      fprintf (stderr, "  value: %d  %f\n", *(variableList[i]->e.v->vAddrI),
	       variableList[i]->value);
    }
  }
  fprintf (stderr, "\n");


}

//////////////////////////////////////////////////////////////////////
// print out a polynomial and its sorting list.  Used for debugging
//////////////////////////////////////////////////////////////////////
void
dismantlePolynomialAndSortingList (struct polynomial *p, struct polyList *l)
{
  int j;

  fprintf (stderr, "Polynomial Value:  %e\n", p->value);
  expPrinting (p);
  fprintf (stderr, "\n");

  for (j = 0; j <= l->listNext - 1; j++) {
    fprintf (stderr, "%4d  value: %e ", j, l->pList[j]->value);
    expPrinting (l->pList[j]);
    fprintf (stderr, "\n");
  }
};

////////////////////////////////////////////////////////////////////////////////
//append a polynomial to a list of polynomials. Used for performance evaluation.
//This function is paired with  polyListSorting3 for construction of evaluation
//list of a polynomial
////////////////////////////////////////////////////////////////////////////////
void
polyListAppend3 (struct polyList *l, struct polynomial *p, int signature)
{
  //valid is a mark showing that this polynomial appears on a sorting list
  p->valid = signature;

  if (l->listNext >= l->listSize) {
    l->pList =
      realloc (l->pList, sizeof (struct polynomial *) * (l->listSize + 100));
    l->listSize = l->listSize + 100;
  }
  l->pList[l->listNext] = p;
  l->listNext++;
};

///////////////////////////////////////////////////////////////////////////////////////
//Count the number of each type of polynomials in a list of a polynomials that determine
//the order of evaluation of a complex polynomial.  This is for performance evaluation
///////////////////////////////////////////////////////////////////////////////////////
void
polyListStatistics (struct polyList *l)
{
  int j;
  int numConstant = 0, numVariable = 0, numSum = 0, numProduct =
    0, numFunctionCall = 0;
  for (j = 0; j <= l->listNext - 1; j++) {
    switch (l->pList[j]->eType) {
    case T_CONSTANT:
      numConstant++;
      break;
    case T_VARIABLE:
      numVariable++;
      break;
    case T_SUM:
      numSum++;
      break;
    case T_PRODUCT:
      numProduct++;
      break;
    case T_FUNCTIONCALL:
      numFunctionCall++;
      break;
    default:
      fprintf (stderr, "Unknown expression type, exit!\n");
      exit (1);
    }
  }
  fprintf (stderr, "%d  %d  %d  %d  %d\n", numConstant, numVariable, numSum,
	   numProduct, numFunctionCall);
};

///////////////////////////////////////////////////////////////////////////////////////
//Building a evaluation order of a complex polynomial.  It is used for performance
//evaluation. For linkage computation, this function allows limited reuse of shared 
//polynomials, which means that shared polynomials within the likelihood polynomial of 
//each pedigree are reused while shared polynomials among likelihood polynomials of
//different pedigrees are not reused
///////////////////////////////////////////////////////////////////////////////////////
void
polyListSorting3 (struct polynomial *p, struct polyList *l, int signature)
{
  int i;

//   fprintf(stderr,"sorting3 valid=%d signature=%d \n",p->valid,signature);

  switch (p->eType) {
  case T_CONSTANT:
    polyListAppend3 (l, p, signature);
    break;

  case T_VARIABLE:
    if (p->valid != signature) {
      polyListAppend3 (l, p, signature);
    }
    break;

  case T_SUM:
    if (p->valid == signature)
      break;

    for (i = 0; i < p->e.s->num; i++)
      if (p->e.s->sum[i]->eType != T_CONSTANT
	  && p->e.s->sum[i]->valid != signature) {
	polyListSorting3 (p->e.s->sum[i], l, signature);
      }
    polyListAppend3 (l, p, signature);
    break;

  case T_PRODUCT:

    if (p->valid == signature)
      break;

    for (i = 0; i < p->e.p->num; i++)
      if (p->e.p->product[i]->eType != T_CONSTANT
	  && p->e.p->product[i]->valid != signature) {
	polyListSorting3 (p->e.p->product[i], l, signature);
      }
    polyListAppend3 (l, p, signature);
    break;

  case T_FUNCTIONCALL:

    if (p->valid == signature)
      break;

    for (i = 0; i < p->e.f->paraNum; i++)
      if (p->e.f->para[i]->eType != T_CONSTANT
	  && p->e.f->para[i]->valid != signature) {
	polyListSorting3 (p->e.f->para[i], l, signature);
      }
    polyListAppend3 (l, p, signature);
    break;

  default:
    break;
  }
//  for(i=0;i<l->listNext;i++)
//     if(l->pList[i]->eType==T_VARIABLE)
//       fprintf(stderr,"List element i=%d of %d l->pList[i].eType=%d  id=%d name=%s\n",
//                                  i,l->listNext,l->pList[i]->eType,l->pList[i]->id,l->pList[i]->vName);
//     else
//       fprintf(stderr,"List element i=%d of %d l->pList[i].eType=%d  id=%d \n",
//                                 i,l->listNext,l->pList[i]->eType,l->pList[i]->id);

};

#include "../../diags/polynomial.c-tail"
