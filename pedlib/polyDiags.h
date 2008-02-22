#include "sw.h"
void dumpPStats();
#ifdef DMTRACK
#warning "Dynamic memory usage dumping is turned on, so performance will be poor!"
#define malloc(X) swMalloc((X), __FILE__, __LINE__)
#define calloc(X,Y) swCalloc((X),(Y), __FILE__, __LINE__)
#define realloc(X,Y) swRealloc((X),(Y), __FILE__, __LINE__)
#define free(X) swFree((X), __FILE__, __LINE__)
#endif

extern struct swStopwatch *overallSW;
extern volatile sig_atomic_t signalSeen;  /* Signalled dumps */

/* Variables for tracking internal polynomial memory usage. */
extern unsigned long maxHashListLength;
extern unsigned long constantPs, constantHashPs, variablePs, variableHashPs,
  sumPs, sumHashPs, productPs, productHashPs, functionPs, functionHashPs;
extern unsigned long peakConstantPs, peakVariablePs, peakSumPs, peakProductPs, peakFunctionPs;
extern unsigned long constantPLExpansions, variablePLExpansions, sumPCollectExpansions,
  sumPTermMergeExpansions, sumPListExpansions, productPCollectExpansions, 
  productPTermMergeExpansions, productPListExpansions;
extern unsigned long constantPsSize, variablePsSize, variablePsExpSize, sumPsSize, sumPColExpSize, 
  sumPTrmMrgExpSize, productPsSize, productPColExpSize, productPTrmMrgExpSize;
/* Display the internal polynomial memory usage statistics. */
void dumpPStats(char *);
void printAllPolynomials();
