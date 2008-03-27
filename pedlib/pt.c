#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "polynomial.h"

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

struct handleInfo *getHandle(char *type) {
  struct handleInfo target, *result;
  printf("%s handle> ", type);
  gets(target.handle);
  result = bsearch (&target, hIList, hI, sizeof (struct handleInfo), compareHandles);
  if (result) {
    preExisting = TRUE;
    return result;
  } else {
    strcpy(hIList[hI].handle, target.handle);
    qsort (hIList, hI, sizeof (struct handleInfo), compareHandles);
    preExisting = FALSE;
    return &hIList[hI++];
  }
}

int main(int argc, char *argv[]) {
  struct handleInfo target, *result, *pHI, *pHIO1, *pHIO2;
  int i, freeFlag;
  double fO1, fO2;
  int eO1, eO2;

  polynomialInitialization();

  fprintf(stderr,"C/V/S/P/?> ");
  while (gets(iB) != NULL) {
    switch (toupper(iB[0])) {
    case '?':
      for (i=0;i<hI;i++) {
        fprintf(stderr,"%i: %s, ", i, hIList[i].handle);
        expPrinting(hIList[i].pP);
        fprintf(stderr," or ");
        expTermPrinting(stderr, hIList[i].pP, 1);
        fprintf(stderr,"\n");
      }
      break;
    case 'C':
      pHI = getHandle("constant poly name");
      if (!preExisting) {
        fprintf(stderr,"Value> ");
        pHI->value = atof(gets(iB));
        pHI->pP = constantExp(pHI->value);
      } else
        fprintf(stderr,"Can't overwrite existing poly\n");
      break;
    case 'V':
      pHI = getHandle("variable poly name");
      if (!preExisting)
        pHI->pP = variableExp(&pHI->value, 0, 'D', pHI->handle);
      else
        fprintf(stderr,"Can't overwrite existing poly\n");
      break;
    case 'S':
      pHI = getHandle("result poly name");
      fprintf(stderr,"Factor of 1st operand> ");
      fO1 = atof(gets(iB));
      if (preExisting) {
        fprintf(stderr,"Assuming 1st operand same as result\n");
        pHIO1 = pHI;
        freeFlag = TRUE;
      } else {
        pHIO1 = getHandle("1st operand poly name");
        freeFlag = FALSE;
      }
      fprintf(stderr,"Factor of 2nd operand> ");
      fO2 = atof(gets(iB));
      pHIO2 = getHandle("2nd operand poly name");
      if (!preExisting)
	printf("Can't use undefined poly name\n");
      else
        pHI->pP = plusExp(2, fO1, pHIO1->pP, fO2, pHIO2->pP, freeFlag);
      break;
    case 'P':
      pHI = getHandle("result poly name");
      if (preExisting) {
        fprintf(stderr,"Assuming 1st operand same as result\n");
        pHIO1 = pHI;
        freeFlag = TRUE;
      } else {
        pHIO1 = getHandle("1st operand poly name");
        freeFlag = FALSE;
      }
      fprintf(stderr,"Exponent of %s> ", pHIO1->handle);
      eO1 = atoi(gets(iB));
      pHIO2 = getHandle("2nd operand poly name");
      if (!preExisting)
	fprintf(stderr,"Can't use undefined poly name\n");
      else {
        fprintf(stderr,"Exponent of %s> ", pHIO2->handle);
        eO2 = atoi(gets(iB));
        pHI->pP = timesExp(2, pHIO1->pP, eO1, pHIO2->pP, eO2, freeFlag);
      }
      break;
    case 'F':
      break;
    }
    fprintf(stderr,"C/V/S/P/?> ");
  }
  return 0;
}
