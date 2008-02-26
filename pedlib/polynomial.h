#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

#define MAX_POLYNOMIAL_KEY   100000000
#define CONSTANT_HASH_SIZE 99991
#define VARIABLE_HASH_SIZE 97
#define SUM_HASH_SIZE      1999993
#define PRODUCT_HASH_SIZE  1999993
#define FUNCTIONCALL_HASH_SIZE 99991
#define HASH_TABLE_INCREASE 20
#include <stdarg.h>
#include <time.h>

#include "../../diags/polynomial.h-head"

clock_t startTime;
clock_t currentTime;
int maxHashLength;
int constantHashHits, variableHashHits, sumHashHits, productHashHits,
  functionHashHits;
struct hashStruct *maxHash;
int sum0, sum1, sum2, sum3, sum4, sum5, sum00, sum11;
int product0, product1, product2, product3, product4, product5, product6,
  product7, product8, product9, product00, product11;
int numSumTerms, numProductTerms;

int maxSumLength, maxProductLength;
int *countSumLength, *countProductLength;
int sizeSumLength, sizeProductLength;

//This is a global variable used for giving each polynomial an unique ID
//so that we can know if two polynomials are the same just from their IDs
int nodeId;
long nodeIdStamp;
long nodeIdStamp2;

//There are following categories of polynomails
//T_CONSTANT    : the polynomial represents a constant value for exampel, 0, 1, 1.5 ...
//T_VARIABLE    : the polynomial represents a variable such as x, y, z, ...
//T_SUM         : the polynomial represents a sum such as 2x+3y+5.6z+...
//T_PRODUCT     : the polynomial represents a product such as x^2y^10z^100
//T_FUNCTIONCALL: the polynomial represents a function call such as log10(x);
enum expressionType
{ T_CONSTANT = 0, T_VARIABLE = 1, T_SUM = 2, T_PRODUCT = 3, T_FUNCTIONCALL =
    4 };

struct polynomial;

//This structure represents the elements of a variable,
//a variavle has an address in memory, a name
struct variablePoly
{
  char vType;
  double *vAddrD;
  int *vAddrI;
  char vName[100];
};

//This structure represents the elements of a sum.
//A sum is composed of a number (saved in num) of items.
//Each item is a factor (saved in factor) times a polynomial (saved in sum).
struct sumPoly
{
  int num;
  struct polynomial **sum;
  double *factor;
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
  int num;
  struct polynomial **product;
  int *exponent;
};

//This structure represents the elements of a function call.
//Each function call is composed of the function name (saved in name)
// and a number (saved in paraNum) of parameters (saved in para)
struct functionPoly
{
  char *name;
  int paraNum;
  struct polynomial **para;
};

//This structure represents a general polynomial.
//A polynomial has a unique id, an polynomial type (eType),
//a value, and the pointer that refers the specific 
//polynomial structure depending on the polynomial type.

typedef struct polynomial
{
  int id;
  int index;
  int key;
  int count;
  enum expressionType eType;
  double value;
//  int                *count;
//  double              *values;
  union
  {
    struct variablePoly *v;	/*variable */
    struct sumPoly *s;		/*sum     */
    struct productPoly *p;	/*product */
    struct functionPoly *f;	/*function */
  } e;
  int valid;
} Polynomial;

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
  int listSize;
  int listNext;
  struct polynomial **pList;
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
  int num;
  int length;
  int *key;
  int *index;
};

//We have a polynomials list to save the polynomials in each polynomial category,
//which is named constantList, variableList, sumList, productList, and functionCallList.
//The length of these lists are saved in constantCount, variableCount, productCount,
//functionCallCount respectively.  Each of the lists is dynamically applied since we are
//not sure how many polynomials we may have.  Therefore, we apply to have a short list
//initially and increase its length gradually with the increase of the number of 
//polynomials. The current length of the lists is recorded in constantListLength,
//variableListLength, sumListLength, productListLength, and functionCallListLength
//respectively.
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
//count the number of polynomials in the evaluation list
void polyListStatistics (struct polyList *l);

//Evaluate a polynomial with the help of the evaluation list 
double evaluatePoly (struct polynomial *pp, struct polyList *l);

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

//statistics of polynomials
void polyStatistics ();

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

#endif
