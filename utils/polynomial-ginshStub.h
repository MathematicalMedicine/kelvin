#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__
/* Copyright (C) 2008, 2010, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */

extern unsigned int nodeId;
extern char *polynomialVersion;

enum expressionType
{ T_CONSTANT = 0, T_VARIABLE = 1, T_SUM = 2, T_PRODUCT = 3, T_FUNCTIONCALL =
    4, T_FREED = 5
};

typedef struct polynomial
{
  unsigned long id;
  unsigned char eType;
  double factor;
  int exponent;
  double value;
  unsigned char vType;
  int *iValue;
  double *dValue;
  unsigned numTerms;
  struct polynomial **terms;
} Polynomial;

typedef struct polyList
{
  int *pList;
} polynomialList;

void polynomialInitialization ();

double evaluateValue (struct polynomial *p);
struct polyList *buildPolyList ();
void polyListSorting (struct polynomial *p, struct polyList *l);
void evaluatePoly (struct polynomial *p, struct polyList *l, double *value);

struct polynomial *constantExp (double con);
struct polynomial *variableExp (double *vD, int *vI, char vType, char name[10]);
struct polynomial *plusExp (int num, ...);
struct polynomial *timesExp (int num, ...);
struct polynomial *functionCallExp (int num, ...);

void expPrinting (struct polynomial *p);
void expTermPrinting (FILE *outputFile, struct polynomial *p, int depth);

void keepPoly (struct polynomial *p);
void holdPoly (struct polynomial *p);
void unHoldPoly (struct polynomial *p);
void freePolys ();
void freeKeptPolys ();
void holdAllPolys ();

void polyStatistics (char *title);
void polyDynamicStatistics (char *title);

#endif
