#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

#define MAX_POLYNOMIAL_KEY   1000000
#define CONSTANT_HASH_SIZE 99
#define SUM_HASH_SIZE      10000
#define PRODUCT_HASH_SIZE  10000
#define FUNCTION_HASH_SIZE 10000
#include <stdarg.h>

//This is a global variable used for giving each polynomial an unique ID
//so that we can know if two polynomials are the same just from their IDs
int nodeId;

//There are following categories of polynomails
//T_CONSTANT    : the polynomial represents a constant value for exampel, 0, 1, 1.5 ...
//T_VARIABLE    : the polynomial represents a variable such as x, y, z, ...
//T_SUM         : the polynomial represents a sum such as 2x+3y+5.6z+...
//T_PRODUCT     : the polynomial represents a product such as x^2y^10z^100
//T_FUNCTIONCALL: the polynomial represents a function call such as log10(x);
enum expressionType {T_CONSTANT=0, T_VARIABLE=1, T_SUM=2, T_PRODUCT=3, T_FUNCTIONCALL=4};

struct polynomial;

//This structure represents the elements of a variable,
//a variavle has an address in memory, a name
struct variablePoly{
   double *vAddr;
   char   vName[10];
};

//This structure represents the elements of a sum.
//A sum is composed of a number (saved in num) of items.
//Each item is a factor (saved in factor) times a polynomial (saved in sum).
struct sumPoly{
   int num;
   struct polynomial **sum;
   double            *factor;
};

//This structure represents the elements of a product.
//A product is composed of a number (saved in num) of items.
//Each item includes an exponent (saved in exponent) and a 
//polynomial (saved in product).
struct productPoly{
   int num;
   struct polynomial **product;
   int                *exponent;
};

//This structure represents the elements of a function call.
//Each function call is composed of the function name (saved in name)
// and a number (saved in paraNum) of parameters (saved in para)
struct functionPoly{
   char *name;
   int              paraNum;
   struct polynomial **para;
};

//This structure represents a general polynomial.
//A polynomial has a unique id, an polynomial type (eType),
//a value, and the pointer that refers the specific 
//polynomial structure depending on the polynomial type.

typedef struct polynomial {
  long id;
  enum expressionType eType;
  double              value;
  union
  {
    struct variablePoly *v;    /*variable*/
    struct sumPoly      *s;    /*sum     */
    struct productPoly  *p;    /*product */
    struct functionPoly *f;    /*function*/ 
  } e;
  int                 valid;
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
typedef struct polyList {
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
struct hashStruct {
     int num;
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
int                constantCount;
int                constantListLength;

struct polynomial **variableList;
int               variableCount;
int               variableListLength;

struct polynomial **sumList;
int                sumCount;
int                sumListLength;

struct polynomial **productList;
int                productCount;
int                productListLength;

struct polynomial **functionCallList;
int                functionCallCount;
int                functionCallListLength;

//Hash tables
//Each category of polynomials has a hash table
struct hashStruct *constantHash;     //Hash table for constant polynomials
struct hashStruct *sumHash;          //Hash table for sum polynomials
struct hashStruct *productHash;      //Hash table for the product polynomials
struct hashStruct *functionHash;     //Hash table for the functionCall polynomials

//The function to evaluate a polynomial recursively
double evaluateValue(struct polynomial *p);

//Count the numbers of sub polynomials in all categories that resides in the evaluation 
//list of a polynomial
void countPoly(struct polyList *l, int *cCounter,
               int *vCounter,int *sCounter,int *pCounter, int *fCounter);
//Determine if a polynomial is 0
int isZeroExp(struct polynomial *p);
//Determine if a polynomial is 1
int isOneExp(struct polynomial *p);
//Determine if a polynomial is -1
int isMinusOneExp(struct polynomial *p);
//Print call the constants
void printConstantList();
//Constructor of a constant polynomial
struct polynomial *constantExp(double con);
//constructor of a variable polynomial
struct polynomial *variableExp(double *v, char name[10]);
//constructor of a sum polynomial
struct polynomial *plusExp(int num, ...);
//constructor of a product polynomial
struct polynomial *timesExp(int num,...);
//constructor of a functionCall polynomial
struct polynomial *functionCallExp(int num, ...);
//Print a polynomial
void expPrinting(struct polynomial *p);
//Creation and Initialization of evaluation list for a polynomial
struct polyList *buildPolyList();
//Append a sub polynomial to the evaluation list of a polynomial
void polyListAppend(struct polyList *l, struct polynomial *p);
//Sort the sub polynomials on the evaluation list of a polynomial.
//throw away the sub polynomials that doesn't need to be evaluated
void polyListSorting(struct polynomial *p, struct polyList *l);
//Sort the sub polynomials on the evaluation list of a polynomial.
//Keep all the sub polynomials in this polynomial
void polyListSorting2(struct polynomial *p, struct polyList *l);
//Evaluate a polynomial with the help of the evaluation list 
double evaluatePoly(struct polynomial *pp, struct polyList *l);
//Print a evaluation list
void printPolyList(struct polyList *l);
//Initialization before polynomials care created and evaluated.
void polynomialInitialization();
//Clear all data structures related to the polynomials
void clear();
//Show the statistics of polynomials
void dismantle();


#endif
