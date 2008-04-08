/*

  Exercise the polynomial code by manually adding constants, variables, and
  two-operand sums and products. Print and evaluate the polynomials.

  If a file name is provided on the command line, it is read as keystrokes
  defining an initial set of polynomials and then the user is prompted for
  further input.

  Very primitive, but it gets the job done.

  Build with: gcc -o pt pt.c polynomial.c -L../utils/ -I../include/ -lutils -lm -lgsl -lgslcblas

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "polynomial.h"

int polynomialScale = 1;
struct swStopwatch *overallSW;

#define TRUE 1==1
#define FALSE !TRUE

struct handleInfo {
  char handle[32];
  double value;
  struct polynomial *pP;
} hIList[1024];
int hI = 0;

char iB[256];
int preExisting = FALSE;

int compareHandles (const void *left, const void *right) {
  int result;
  struct handleInfo *mHILeft, *mHIRight;
  mHILeft = (struct handleInfo *) left;
  mHIRight = (struct handleInfo *) right;
  return strcmp (mHILeft->handle, mHIRight->handle);
}

struct handleInfo *getHandle(FILE *inputFile, FILE *outputFile, char *type) {
  struct handleInfo target, *result;
  fprintf(outputFile, "%s > ", type);
  fgets(target.handle, sizeof(target.handle), inputFile);
  target.handle[strlen(target.handle)-1] = 0;
  result = bsearch (&target, hIList, hI, sizeof (struct handleInfo), compareHandles);
  if (result)
    preExisting = TRUE;
  else {
    strcpy(hIList[hI++].handle, target.handle);
    qsort (hIList, hI, sizeof (struct handleInfo), compareHandles);
    preExisting = FALSE;
    result = bsearch (&target, hIList, hI, sizeof (struct handleInfo), compareHandles);
  }
  return result;
}

void loopReading(FILE *inputFile, FILE *outputFile) {
  struct handleInfo *pHI, *pHIO1, *pHIO2;
  int i, freeFlag;
  double fO1, fO2;
  int eO1, eO2;
  char *promptString = "C/V/S/P/E/#/%%/?/q> ";

  fprintf(outputFile, promptString);
  while (fgets(iB, sizeof(iB), inputFile) != NULL) {
    switch (toupper(iB[0])) {
    case 'Q':			/* Quit! */
      exit(EXIT_SUCCESS);
    case '?':			/* Dump our work thus far */
      for (i=0;i<hI;i++) {
        fprintf(stderr,"%i: %s ", i, hIList[i].handle);
	expPrinting(hIList[i].pP);
	fprintf(stderr," or ");
	expTermPrinting(stderr, hIList[i].pP, 1);
	fprintf(stderr,"\n");
      }
      break;
    case '%':
      polyStatistics("Requested");
      break;
    case 'C':			/* Add a constant polynomial */
      pHI = getHandle(inputFile, outputFile,"constant poly name");
      if (!preExisting) {
        fprintf(outputFile,"Value> ");
        pHI->value = atof(fgets(iB, sizeof(iB), inputFile));
        pHI->pP = constantExp(pHI->value);
      } else
        fprintf(outputFile,"Can't overwrite existing poly\n");
      break;
    case 'V':			/* Add a variable polynomial */
      pHI = getHandle(inputFile, outputFile,"variable poly name (case-sensistive)");
      if (!preExisting)
        pHI->pP = variableExp(&pHI->value, 0, 'D', pHI->handle);
      else
        fprintf(outputFile,"Can't overwrite existing poly\n");
      break;
    case 'S':			/* Add a 2-operand sum polynomial */
      pHI = getHandle(inputFile, outputFile,"result poly name");
      fprintf(outputFile,"Factor of 1st operand> ");
      fO1 = atof(fgets(iB, sizeof(iB), inputFile));
      if (preExisting) {
        fprintf(outputFile,"Assuming 1st operand same as result\n");
        pHIO1 = pHI;
        freeFlag = TRUE;
      } else {
        pHIO1 = getHandle(inputFile, outputFile,"1st operand poly name");
        freeFlag = FALSE;
	if (!preExisting) {
	  fprintf(outputFile,"Can't use undefined poly name for 1st\n");
	  break;
	}
      }
      fprintf(outputFile,"Factor of 2nd operand> ");
      fO2 = atof(fgets(iB, sizeof(iB), inputFile));
      pHIO2 = getHandle(inputFile, outputFile,"2nd operand poly name");
      if (!preExisting) {
	printf("Can't use undefined poly name for 2nd\n");
	break;
      }
      pHI->pP = plusExp(2, fO1, pHIO1->pP, fO2, pHIO2->pP, freeFlag);
      break;
    case 'P':			/* Add a 2-operand product polynomial */
      pHI = getHandle(inputFile, outputFile,"result poly name");
      if (preExisting) {
        fprintf(outputFile,"Assuming 1st operand same as result\n");
        pHIO1 = pHI;
        freeFlag = TRUE;
      } else {
        pHIO1 = getHandle(inputFile, outputFile,"1st operand poly name");
        freeFlag = FALSE;
	if (!preExisting) {
	  printf("Can't use undefined poly name\n");
	  break;
	}
      }
      fprintf(outputFile,"Exponent of %s> ", pHIO1->handle);
      eO1 = atoi(fgets(iB, sizeof(iB), inputFile));
      pHIO2 = getHandle(inputFile, outputFile,"2nd operand poly name");
      if (!preExisting) {
	fprintf(outputFile,"Can't use undefined poly name\n");
	break;
      }
      fprintf(outputFile,"Exponent of %s> ", pHIO2->handle);
      eO2 = atoi(fgets(iB, sizeof(iB), inputFile));
      pHI->pP = timesExp(2, pHIO1->pP, eO1, pHIO2->pP, eO2, freeFlag);
      break;
    case 'F':
      break;
    case '#':			/* Comment line */
      fprintf(stderr, "%s", iB);
      break;
    case 'E':			/* Evaluate the polynomial */
      pHI = getHandle(inputFile, outputFile," poly to evaluate");
      if (!preExisting) {
	fprintf(outputFile, "Can't evaluate undefined poly\n");
	break;
      }
      for (i=0;i<hI;i++) {
	if (hIList[i].pP->eType == T_VARIABLE) {
	  fprintf(outputFile,"Value for %s> ", hIList[i].handle);
	  hIList[i].value = atof(fgets(iB, sizeof(iB), inputFile));
	}
      }
      fprintf(outputFile, "=%G\n", evaluateValue(pHI->pP));
    }
    fprintf(outputFile, promptString);
  }
  fprintf(outputFile,"\n");
  return;
}

int main(int argc, char *argv[]) {
  FILE *initFile, *nullFile;

  overallSW = swCreate("overall");
  swStart(overallSW);
  polynomialInitialization();

  if (argc > 1) {
    if ((initFile = fopen(argv[1],"r")) != NULL) {
      nullFile = fopen("/dev/null","w");
      loopReading(initFile, nullFile);
      close(initFile);
    } else {
      perror("Cannot open input file");
      return(1);
    }
  }
  loopReading(stdin, stderr);
  swStop(overallSW);
  swDump(overallSW);

  return(0);
}
