
/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__
#include <stdarg.h>
#include <time.h>
#include <signal.h>

#include "sw.h"

//All the following constants are prime numbers
#define MAX_POLYNOMIAL_KEY   2147483629

/* Hash and other storage sizes for claustrophobic environments, must
   be scaled-up for larger models. Original was 10x, and that's the
   default if no other value is provided. */
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
#define MIN_SUM_LIST_INITIAL 100000
#define MIN_SUM_LIST_INCREASE 10000
#define MIN_PRODUCT_LIST_INITIAL 100000
#define MIN_PRODUCT_LIST_INCREASE 10000
#define MIN_FUNCTIONCALL_LIST_INITIAL 100
#define MIN_FUNCTIONCALL_LIST_INCREASE 10

#ifdef DMTRACK
#warning "Dynamic memory usage dumping is turned on, so performance will be poor!"
#define malloc(X) swMalloc((X), __FILE__, __LINE__)
#define calloc(X,Y) swCalloc((X),(Y), __FILE__, __LINE__)
#define realloc(X,Y) swRealloc((X),(Y), __FILE__, __LINE__)
#define free(X) swFree((X), __FILE__, __LINE__)
#endif

/* The following dynamically-maintained variables are for debugging. */
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

/* There are following categories of polynomials
   T_CONSTANT    : the polynomial represents a constant value for exampel, 0, 1, 1.5 ...
   T_VARIABLE    : the polynomial represents a variable such as x, y, z, ...
   T_SUM         : the polynomial represents a sum such as 2x+3y+5.6z+...
   T_PRODUCT     : the polynomial represents a product such as x^2y^10z^100
   T_FUNCTIONCALL: the polynomial represents a function call such as log10(x)
   T_FREED       : a polynomial that was freed but the structure kept for diagnostics.
*/
enum expressionType
{ T_CONSTANT = 0, T_VARIABLE = 1, T_SUM = 2, T_PRODUCT = 3, T_FUNCTIONCALL =
    4, T_FREED = 5
};

struct polynomial;

//This structure represents the elements of a variable,
//a variable has an address in memory, a name.
//In general, the number of variable polynomials is very small, say 20.
struct variablePoly
{
  char vType;			//variable type
  union
  {
    double *vAddrD;		//variable address for a double variable
    int *vAddrI;		//variable address for an integer variable
  } vAddr;
  char vName[100];		//variable name
};

//This structure represents the elements of a sum.
//A sum is composed of a number (saved in num) of items.
//Each item is a factor (saved in factor) times a polynomial (saved in sum).
struct sumPoly
{
  int num;			//number of terms
  struct polynomial **sum;	// polynomial terms
  double *factor;		//factors for polynomial terms

  //Answer to the question of why not using constant polynomials to express factors
  //and a similar question that there is inconsistency in expressing polynomials
  //because we have constant polynomials while we have double factors
  //
  //   
  //Here, factor is to express components of a sum polynomial.  Factors are constants.  
  //They can be expressed as constant polynomials.  However, constant polynomials 
  //need more memory than constants.  Therefore, here, we don't use constant polynomials.
  //Instead, we use constants firectly.


};

   //Following variables are for building sum polynomials
   //collect variable terms
double *factor_v1;
struct polynomial **p_v1;
int containerLength_v1;
int counter_v1;

   //collect product terms
double *factor_p1;
struct polynomial **p_p1;
int containerLength_p1;
int counter_p1;

   //collect function call terms
double *factor_f1;
struct polynomial **p_f1;
int containerLength_f1;
int counter_f1;

   //collect all terms
double *factorSum;
struct polynomial **pSum;
int lengthSum;


   //Following variables are for building product polynomials
   //collect variable terms
int *exponent_v2;
struct polynomial **p_v2;
int containerLength_v2;
int counter_v2;

   //collect sum terms
int *exponent_s2;
struct polynomial **p_s2;
int containerLength_s2;
int counter_s2;

   //collect function call terms
int *exponent_f2;
struct polynomial **p_f2;
int containerLength_f2;
int counter_f2;

   //collect all terms
int *exponentProd;
struct polynomial **pProd;
int lengthProd;

//This structure represents the elements of a product.
//A product is composed of a number (saved in num) of items.
//Each item includes an exponent (saved in exponent) and a 
//polynomial (saved in product).
struct productPoly
{
  int num;			//numer of terms
  struct polynomial **product;	//polynomial terms
  int *exponent;		//exponents for polynomial terms
};

//This structure represents the elements of a function call.
//Each function call is composed of the function name (saved in name)
// and a number (saved in paraNum) of parameters (saved in para)
struct functionPoly
{
  int paraNum;			//number of parameters
  struct polynomial **para;	//parameters
  char *name;			//function name
};

//This structure represents a general polynomial.
//A polynomial has a unique id, an polynomial type (eType),
//a value, and the pointer that refers the specific 
//polynomial structure depending on the polynomial type.

typedef struct polynomial
{
  unsigned int id;		//unique id
  int index;			//index in a polynomial list
  int key;			//key of the polynomial
  unsigned short count;		// Hold reference count
  unsigned char valid;		// Preservation flag(s)
  unsigned char eType;		//polynomial type: 
  double value;			//value saves the value of the polynomial
  //
  //Following is the answer to the question "why do we need constant polynomials instead of constants directly?"
  //
  //When we try to understand why we need constant polynomial, we may need to consider that there are two worlds:
  //polynomial world and non polynomial world.  In non polynomial world, the result of mathematical computations is
  //a value.  On the other hand, in polynomial world, the result of mathematical computations is a polynomial, no
  //matter the result is a constant, a sum of a group of terms, or a product of a group of terms.  Therefore, 
  //constant polynomials are for purpose of consistency in polynomial calculations and polynomial expressions.  
  //If the result of polynomial calculations is a constant, it has to be expressed as a constant polynomial because 
  //the calculations are in polynomial world. 
  union
  {
    struct variablePoly *v;	/*variable */
    struct sumPoly *s;		/*sum     */
    struct productPoly *p;	/*product */
    struct functionPoly *f;	/*function */
  } e;
} Polynomial;

#define VALID_EVAL_FLAG 1	/* Used by polyListSorting to manage inclusion */
#define VALID_KEEP_FLAG 2	/* Weaker than HOLD, only kept until a freeKeptPolys() call */
#define VALID_REF_FLAG 4	/* Weakest of all, but keeps 1st freeing on-track */

//List is for polynomail evaluation.  When we evaluate a polynomial,
//it is possible that we evaluate only some parts of it because the
//other parts have been evaluated when other polynomials are evaluated.
//List is a structure that records the sub polynomials that need evaluation
//and the evaluation order of these sub polynomials.
//listSize is how many sub polynomials this list can hold.
//listNext is the index for next sub polynomials to be insered in
//pList is the list of sub polynomials that should be evaluated
//for the evaluation of the polynomial
typedef struct polyList
{
  int listSize;			//size of the list preallocated
  int listNext;			//next free position
  struct polynomial **pList;	//list of polynomials for evaluation
} polynomialList;

//hashStruct is used for speeding up the creation of polynomials.
//Without hash struct, a polynomial has to compare with all the 
//currently existing polynomials in this category in order to 
//know if this polynomial has been created in store.
//num shows the total number of polynomials that fall into this hash item
//key saves keys of all the polynomials that fall into this hash item
//index saves indexes of all the polynomials that fall into this item
struct hashStruct
{
  unsigned short num;		//number of polynomials in a hash bucket
  unsigned short length;	//length of the hash bucket preallocated
  int *key;			//keys of polynomials
  int *index;			//indexes of polynomials
};

//Hash tables
//Each category of polynomials has a hash table
struct hashStruct *constantHash;	//Hash table for constant polynomials
struct hashStruct *variableHash;	//Hash table for variable polynomials
struct hashStruct *sumHash;	//Hash table for sum polynomials
struct hashStruct *productHash;	//Hash table for the product polynomials
struct hashStruct *functionCallHash;	//Hash table for the functionCall polynomials

//The function to evaluate a polynomial recursively
double evaluateValue (struct polynomial *p);

//Count the numbers of sub polynomials in all categories that resides in the evaluation 
//list of a polynomial
void countPoly (struct polyList *l, int *cCounter,
		int *vCounter, int *sCounter, int *pCounter, int *fCounter);
//Determine if a polynomial is 0
int isZeroExp (struct polynomial *p);

//Determine if a polynomial is 1
int isOneExp (struct polynomial *p);

//Determine if a polynomial is -1
int isMinusOneExp (struct polynomial *p);

//Constructor of a constant polynomial
struct polynomial *constantExp (double con);

//constructor of a variable polynomial
struct polynomial *variableExp (double *vD, int *vI, char vType,
				char name[10]);
//constructor of a sum polynomial
struct polynomial *plusExp (int num, ...);

//constructor of a product polynomial
struct polynomial *timesExp (int num, ...);

//constructor of a functionCall polynomial
struct polynomial *functionCallExp (int num, ...);

//Print a polynomial
void expPrinting (struct polynomial *p);

//Creation and Initialization of evaluation list for a polynomial
struct polyList *buildPolyList ();

//Append a polynomial to the evaluation list of a polynomial
void polyListAppend (struct polyList *l, struct polynomial *p);

//Append a polynomial to the evaluation list if the polynomial with this signature doesn't exist
void polyListAppend3 (struct polyList *l, struct polynomial *p,
		      int signature);
//Sort the sub polynomials on the evaluation list of a polynomial.
//Throw away the sub polynomials that doesn't need to be evaluated.
//Sharing is among different polynomials.
void polyListSorting (struct polynomial *p, struct polyList *l);

//Sort the sub polynomials on the evaluation list of a polynomial.
//Keep all the sub polynomials in this polynomial, no sharing
void polyListSorting2 (struct polynomial *p, struct polyList *l);

//Sort the sub polynomials on the evaluation list of a polynomial.
//Keep a separate evaluation list for this polynomial, sharing only inside this polynomial
void polyListSorting3 (struct polynomial *p, struct polyList *l,
		       int signature);
//Evaluate a polynomial with the help of the evaluation list 
void evaluatePoly (struct polynomial *pp, struct polyList *l, double *pd);

//Print a evaluation list
void printPolyList (struct polyList *l);

//record the current status of polynomials
void makePolynomialStamp ();
void makePolynomialStamp2 ();

//partially clear the polynomials
void partialPolynomialClearance ();
void partialPolynomialClearance2 ();

//Initialization before polynomials care created and evaluated.
void polynomialInitialization ();

//Clear all data structures related to the polynomials
void polynomialClearance ();

//Show polynomials' structures
void dismantle ();

//Statistics of polynomials
void polyStatistics (char *);
void polyDynamicStatistics (char *);

//print hash tables
void printHashTables ();

//print out all the variables
void printAllVariables ();

//print out all the polynomials
void printAllPolynomials ();
void dismantlePolynomialAndSortingList (struct polynomial *p,
					struct polyList *l);
//print out polylist
//void printPolyList(struct polyList *l)
void printSummaryPoly (struct polynomial *);

void keepPoly (struct polynomial *);
void holdPoly (struct polynomial *);
void unHoldPoly (struct polynomial *);
void freePolys ();
void freeKeptPolys ();
void holdAllPolys ();
void expTermPrinting (FILE *, struct polynomial *, int);
void printAllPolynomials ();

#endif
