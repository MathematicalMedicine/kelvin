/*

  Exercise the polynomial code by manually adding constants, variables, and
  two-operand sums and products. Print and evaluate the polynomials.

  If a file name is provided on the command line, it is read as keystrokes
  defining an initial set of polynomials and then the user is prompted for
  further input.

  Very primitive, but it gets the job done.

  Build with: gcc -o pt pt.c -DPOLYUSE_DL -DPOLYCODE_DL polynomial.c -L../utils/ -I../include/ -lutils -lm -ldl

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "polynomial.h"
#include "sw.h"

int polynomialScale = 1;
struct swStopwatch *overallSW;

#define MAXHANDLE 32

struct handleInfo {
  char handle[MAXHANDLE];
  double *pValue;
  struct polynomial *pP;
} hIList[1024];
int hI = 0;

char iB[256];
int preExisting = FALSE;

FILE *responseFile = NULL;

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
  if (responseFile != NULL)
    fprintf(responseFile, "%s\n", target.handle);
  result = bsearch (&target, hIList, hI, sizeof (struct handleInfo), compareHandles);
  if (result)
    preExisting = TRUE;
  else {
    hIList[hI].pValue = (double *) malloc(sizeof(double));
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
  struct polyList *pL;
  char *promptString = "C/V/S/P/F/E/G/L/#/%%/?/Q> ";
  char polyName[16];

  fprintf(outputFile, promptString);
  while (fgets(iB, sizeof(iB), inputFile) != NULL) {
    if (responseFile != NULL)
      fprintf(responseFile, "%s", iB);
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
        *pHI->pValue = atof(fgets(iB, sizeof(iB), inputFile));
	if (responseFile != NULL)
	  fprintf(responseFile, "%s", iB);
        pHI->pP = constantExp(*pHI->pValue);
      } else
        fprintf(stderr,"Can't overwrite existing poly\n");
      break;
    case 'V':			/* Add a variable polynomial */
      pHI = getHandle(inputFile, outputFile,"variable poly name (case-sensistive)");
      if (!preExisting)
        pHI->pP = variableExp(pHI->pValue, 0, 'D', pHI->handle);
      else
        fprintf(stderr,"Can't overwrite existing poly\n");
      break;
    case 'S':			/* Add a 2-operand sum polynomial */
      pHI = getHandle(inputFile, outputFile,"result poly name");
      fprintf(outputFile,"Factor of 1st operand> ");
      fO1 = atof(fgets(iB, sizeof(iB), inputFile));
      if (responseFile != NULL)
	fprintf(responseFile, "%s", iB);
      if (preExisting) {
        fprintf(outputFile,"Assuming 1st operand same as result\n");
        pHIO1 = pHI;
        freeFlag = TRUE;
      } else {
        pHIO1 = getHandle(inputFile, outputFile,"1st operand poly name");
        freeFlag = FALSE;
	if (!preExisting) {
	  fprintf(stderr,"Can't use undefined poly name for 1st\n");
	  break;
	}
      }
      fprintf(outputFile,"Factor of 2nd operand> ");
      fO2 = atof(fgets(iB, sizeof(iB), inputFile));
      if (responseFile != NULL)
	fprintf(responseFile, "%s", iB);
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
      if (responseFile != NULL)
	fprintf(responseFile, "%s", iB);
      pHIO2 = getHandle(inputFile, outputFile,"2nd operand poly name");
      if (!preExisting) {
	fprintf(stderr,"Can't use undefined poly name\n");
	break;
      }
      fprintf(outputFile,"Exponent of %s> ", pHIO2->handle);
      eO2 = atoi(fgets(iB, sizeof(iB), inputFile));
      if (responseFile != NULL)
	fprintf(responseFile, "%s", iB);
      pHI->pP = timesExp(2, pHIO1->pP, eO1, pHIO2->pP, eO2, freeFlag);
      break;
    case 'F':			/* Add a single-operand function call polynomial */
      pHI = getHandle(inputFile, outputFile,"result poly name");
      if (preExisting) {
	printf("Can't overwrite existing poly\n");
	break;
      }
      pHIO1 = getHandle(inputFile, outputFile,"operand poly name");
      if (!preExisting) {
	printf("Can't use undefined poly name\n");
	break;
      }
      fprintf(outputFile,"Function name> ");
      fgets(iB, sizeof(iB), inputFile);
      iB[strlen(iB)-1] = 0;
      if (responseFile != NULL)
	fprintf(responseFile, "%s\n", iB);
      pHI->pP = functionCallExp(2, iB, pHIO1->pP);
      break;
    case '#':			/* Comment line */
      fprintf(stderr, "%s", iB);
      break;
    case 'E':			/* Evaluate the polynomial */
      pHIO1 = getHandle(inputFile, outputFile,"result variable poly name");
      if (!preExisting) {
	printf("Can't use undefined poly name\n");
	break;
      }
      if (pHIO1->pP->eType != T_VARIABLE) {
	printf("Can't assign to non-variable poly\n");
	break;
      }
      pHI = getHandle(inputFile, outputFile," poly to evaluate");
      if (!preExisting) {
	fprintf(stderr, "Can't evaluate undefined poly\n");
	break;
      }
      for (i=0;i<hI;i++) {
	if (hIList[i].pP->eType == T_VARIABLE) {
	  fprintf(outputFile,"Value for %s> ", hIList[i].handle);
	  if (hIList[i].pP == pHIO1->pP) {
	    fprintf(outputFile,"(return to keep %G) ", hIList[i].pP->value);
	    fgets(iB, sizeof(iB), inputFile);
	    if (strlen(iB) >= 2)
	      *hIList[i].pValue = atof(iB);
	    else
	      *hIList[i].pValue = hIList[i].pP->value;
	  } else
	    *hIList[i].pValue = atof(fgets(iB, sizeof(iB), inputFile));
	  if (responseFile != NULL)
	    fprintf(responseFile, "%s", iB);
	}
	//	fprintf(outputFile,"Value for %s is %G\n", hIList[i].handle, *hIList[i].pValue);
      }
      fprintf(stdout, "=%20.15g\n", (pHIO1->pP->value = evaluateValue(pHI->pP)));
      break;
    case 'G':			/* Graph the polynomial */
      pHI = getHandle(inputFile, outputFile," poly to graph");
      if (!preExisting) {
	fprintf(stderr, "Can't graph undefined poly\n");
	break;
      }
      writePolyDigraph(pHI->pP);
      fprintf(outputFile, "Dot-format digraph written to pD_%d.dot\n", pHI->pP->id);
      break;
    case 'L':			/* Generate polynomial sort list and C function */
      pHI = getHandle(inputFile, outputFile," poly to process");
      if (!preExisting) {
	fprintf(stderr, "Can't process undefined poly\n");
	break;
      }
      pL = buildPolyList();
      polyListSorting(pHI->pP, pL);
      sprintf (polyName, "P%d", nodeId - 1);
      codePoly(pHI->pP, pL, polyName);
      fprintf(outputFile, "Polynomial C function written to %s.c\n", polyName);
      break;
    default:
      fprintf(stderr, "One of:\n"
	      "C/V/S/P - create a Constant/Variable/Sum/Product poly\n"
	      "E - evaluate a poly\n"
	      "G - generate a digraph of a poly\n"
	      "L - generate poly list and separate C function\n"
	      "%% - display poly statistics\n"
	      "? - print all polys\n"
	      "Q or <EOF> - quit\n");
      break;
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
  if ((responseFile = fopen("pt.log","w")) == NULL) {
    perror("Cannot open response log file");
    return(1);
  }
  loopReading(stdin, stderr);
  fclose(responseFile);
  swStop(overallSW);
  swDump(overallSW);

  return(0);
}
