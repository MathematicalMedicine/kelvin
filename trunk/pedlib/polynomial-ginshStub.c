#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "polynomial.h"

unsigned int nodeId;
char *polynomialVersion = "stub";
char *expressionTypeNames[6] = {"C", "V", "S", "P", "F", "-"};
Polynomial *variableList[8192];
int variableCount = 0;
FILE *outputFile;

void polynomialInitialization ()
{
  if ((outputFile = fopen("ginsh.in","w")) == NULL) {
    perror ("Cannot open polynomial output file\n");
    exit(EXIT_FAILURE);
  }
}

double evaluateValue (struct polynomial *p)
{
  double v = 0;
  struct polyList *l = NULL;

  evaluatePoly (p, l, &v);
  return v;
}

struct polyList *buildPolyList ()
{
  struct polyList *l = NULL;
  l = (struct polyList *) malloc(sizeof(struct polyList));
  l->pList = (int *) malloc(sizeof(int));
  return l;
}

void polyListSorting (struct polynomial *p, struct polyList *l)
{
}

void evaluatePoly (struct polynomial *p, struct polyList *l, double *value)
{
  int i;
  Polynomial *q;

  for (i=0; i<variableCount; i++) {
    q = variableList[i];
    if (q->vType == 'I')
      fprintf (outputFile, "V%lu = %d:\n", q->id, *q->iValue);
    else
      fprintf (outputFile, "V%lu = %G:\n", q->id, *q->dValue);
  }
  fprintf (outputFile, "%s%lu;\n", expressionTypeNames[p->eType], p->id);
  fflush (outputFile);
  *value = .005;
}

struct polynomial *constantExp (double con)
{
  struct polynomial *p = NULL;
  p = (Polynomial *) malloc(sizeof(Polynomial));
  p->id = nodeId++;
  p->eType = T_CONSTANT;
  fprintf (outputFile, "C%lu = %G:\n", p->id, con);
  return p;
}

struct polynomial *variableExp (double *vD, int *vI, char vType, char name[10])
{
  struct polynomial *p = NULL;
  p = (Polynomial *) malloc(sizeof(Polynomial));
  p->id = nodeId++;
  variableList[variableCount] = p;
  variableCount++;
  p->eType = T_VARIABLE;
  p->vType = vType;
  p->iValue = vI;
  p->dValue = vD;
  return p;
}
/*
struct polynomial *plusExp (int num, ...)
{
  Polynomial *p = NULL, *q;
  int i, flag;
  va_list args;
  double f;
  char term[128];
  char rightBuffer[256];

  rightBuffer[0] = '\0';
  va_start (args, num);
  for (i = 0; i < num; i++) {
    f = va_arg (args, double);
    q = va_arg (args, Polynomial *);
    if (i != 0)
      strcat (rightBuffer, "+");
    else
      p = q;
    sprintf (term, "%G*%s%lu", f, expressionTypeNames[q->eType], q->id);
    strcat (rightBuffer, term);
  }
  strcat (rightBuffer, ":\n");
  flag = va_arg (args, int);
  va_end (args);

  if (flag == 0) {
    p = (Polynomial *) malloc(sizeof(Polynomial));
    p->id = nodeId++;
    p->eType = T_SUM;
  }
  fprintf (outputFile, "%s%lu = %s", expressionTypeNames[p->eType], p->id, rightBuffer);

  //  if ((nodeId & 0x1FFFFF) == 0)
  if ((nodeId & 0x1F) == 0)
    fprintf (outputFile, "VC = %d;\n", nodeId);

  return p;
}
*/
struct polynomial *plusExp (int num, ...)
{
  Polynomial *p = NULL, *q;
  int i;
  va_list args;
  double f;

  p = (Polynomial *) malloc(sizeof(Polynomial));
  p->id = nodeId++;
  p->eType = T_SUM;
  fprintf (outputFile, "S%lu = ", p->id);
  va_start (args, num);
  for (i = 0; i < num; i++) {
    if (i != 0) fprintf (outputFile, "+");
    f = va_arg (args, double);
    q = va_arg (args, Polynomial *);
    fprintf (outputFile, "%G*%s%lu", f, expressionTypeNames[q->eType], q->id);
  }
  fprintf (outputFile, ":\n");
  va_end (args);

  if ((nodeId & 0x1F) == 0)
    fprintf (outputFile, "VC = %d;\n", nodeId);

  return p;
}
/*
struct polynomial *timesExp (int num, ...)
{
  Polynomial *p = NULL, *q;
  int i, flag;
  va_list args;
  int e;
  char term[128];
  char rightBuffer[256];

  rightBuffer[0] = '\0';
  va_start (args, num);
  for (i = 0; i < num; i++) {
    q = va_arg (args, Polynomial *);
    e = va_arg (args, int);
    if (i != 0)
      strcat (rightBuffer, "*");
    else
      p = q;
    sprintf (term, "%s%lu^%d", expressionTypeNames[q->eType], q->id, e);
    strcat (rightBuffer, term);
  }
  strcat (rightBuffer, ":\n");
  flag = va_arg (args, int);
  va_end (args);

  if (flag == 0) {
    p = (Polynomial *) malloc(sizeof(Polynomial));
    p->id = nodeId++;
    p->eType = T_PRODUCT;
  }
  fprintf (outputFile, "%s%lu = %s", expressionTypeNames[p->eType], p->id, rightBuffer);

  //  if ((nodeId & 0x1FFFFF) == 0)
  if ((nodeId & 0x1F) == 0)
    fprintf (outputFile, "VC = %d;\n", nodeId);

  return p;
}
*/

struct polynomial *timesExp (int num, ...)
{
  Polynomial *p = NULL, *q;
  int i;
  va_list args;
  int e;

  p = (Polynomial *) malloc(sizeof(Polynomial));
  p->id = nodeId++;
  p->eType = T_PRODUCT;
  fprintf (outputFile, "P%lu = ", p->id);
  va_start (args, num);
  for (i = 0; i < num; i++) {
    if (i != 0) fprintf (outputFile, "*");
    q = va_arg (args, Polynomial *);
    e = va_arg (args, int);
    fprintf (outputFile, "%s%lu^%d", expressionTypeNames[q->eType], q->id, e);
  }
  fprintf (outputFile, ":\n");
  va_end (args);

  if ((nodeId & 0x1F) == 0)
    fprintf (outputFile, "VC = %d;\n", nodeId);

  return p;
}


struct polynomial *functionCallExp (int num, ...)
{
  struct polynomial *p = NULL;
  p = (Polynomial *) malloc(sizeof(Polynomial));
  p->id = nodeId++;
  p->eType = T_FUNCTIONCALL;
  return p;
}


void expPrinting (struct polynomial *p)
{
}

void expTermPrinting (FILE *outputFile, struct polynomial *p, int depth)
{
}


void keepPoly (struct polynomial *p)
{
}

void holdPoly (struct polynomial *p)
{
}

void unHoldPoly (struct polynomial *p)
{
}

void freePolys ()
{
}

void freeKeptPolys ()
{
}

void holdAllPolys ()
{
}


void polyStatistics (char *title)
{
}

void polyDynamicStatistics (char *title)
{
}

