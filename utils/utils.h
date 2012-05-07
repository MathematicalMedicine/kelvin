/**
@file utils.h

  @see utils.c
  
  Copyright &copy; 2010, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

  @author Alberto Maria Segre
  @author Bill Valentine-Cooper

*/

#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdio.h>
#include <ctype.h>
#include "sw.h"

/**********************************************************************
 * Some convenient and commonly-used defines.
 **********************************************************************/

/**

  Macros for malloc and friends allows us to report module, line
  number, pointer name and cast for failures. Makes it very easy to
  diagnose problems without re-runs or excessive questioning.

  Casts aren't really necessary, but they can help clarify and constrain.

*/

#ifdef DMTRACK
//#warning "Dynamic memory usage dumping is turned on, so performance will be poor!"
#endif
#if defined (DMTRACK)
#define malloc(X) swMalloc((X), __FILE__, __LINE__)
#define calloc(X,Y) swCalloc((X),(Y), __FILE__, __LINE__)
#define realloc(X,Y) swRealloc((X),(Y), __FILE__, __LINE__)
#define free(X) swFree((X), __FILE__, __LINE__)
#endif

#define MALCHOKE(pChunk,chunkSize,chunkCast) \
do { \
  if ((pChunk = (chunkCast) malloc(chunkSize)) == 0) { \
    fprintf (stderr,"FATAL - ABORTING (%s:%d), malloc of %lu bytes for variable %s cast as (%s) failed!\n", \
	     (__FILE__),(__LINE__),(unsigned long)chunkSize,#pChunk,#chunkCast); \
    exit (EXIT_FAILURE); \
  } \
} while(0)

#define REALCHOKE(pChunk,chunkSize,chunkCast) \
do {									\
  if ((pChunk = (chunkCast) realloc(pChunk, chunkSize)) == 0) {		\
    fprintf (stderr,"FATAL - ABORTING (%s:%d), realloc to %lu bytes of variable %s cast as (%s) failed!\n", \
	     (__FILE__),(__LINE__),(unsigned long)chunkSize,#pChunk,#chunkCast); \
    exit (EXIT_FAILURE);						\
  } \
} while(0)

#define CALCHOKE(pChunk,chunkSize,chunkCount,chunkCast)			\
do { \
  if ((pChunk = (chunkCast) calloc(chunkSize, chunkCount)) == 0) {	\
    fprintf (stderr,"FATAL - ABORTING (%s:%d), calloc of %lu * %lu zeroed bytes of variable %s cast as (%s) failed!\n", \
	     (__FILE__),(__LINE__),(unsigned long)chunkSize,(unsigned long)chunkCount,#pChunk,#chunkCast); \
    exit (EXIT_FAILURE);						\
  } \
} while(0)

/**********************************************************************
 * NULL/FALSE/TRUE are probably already defined.
 **********************************************************************/
#ifndef NULL
#define NULL 0
#endif

#ifndef FALSE
#define FALSE 0
#define TRUE !FALSE
#endif
#define NOT(n) (((n) != FALSE)?0:1)

/* Unfortunately, ERROR is also TRUE. */
#ifndef ERROR
#define ERROR -1
#endif

/**********************************************************************
 * Max/Min/Abs macros. 
 **********************************************************************/
#ifndef MAX
#define MAX(m,n) ((m)>=(n)?(m):(n))
#endif
#ifndef MIN
#define MIN(m,n) ((m)<=(n)?(m):(n))
#endif
#ifndef ABS
#define ABS(m) ((m)>=0?(m):-(m))
#endif
#ifndef FABS
#define FABS(m) ((m)>=0.0?(m):-(m))
#endif

/**********************************************************************
 * Odd/Even macros. 
 **********************************************************************/
#define ODD(n) (((n)%2)==1)
#define EVEN(n) (((n)%2)==0)

/**********************************************************************
 * Pos/Neg macros. 
 **********************************************************************/
#define POS(n) ((n)>0)
#define NEG(n) ((n)<0)
#define SGN(n) ((n)<0?-1:1)
#define FPOS(n) ((n)>0.0)
#define FNEG(n) ((n)<0.0)
#define FSGN(n) ((n)<0.0?-1.0:1.0)

/**********************************************************************
 * Generate a random integer between 0 and n-1 (inclusive).
 **********************************************************************/
#define SEED(n) (srandom (n));
#define RANDINT(n) (random() % (n))
#define SEEDDBL(n) (srand48(n));
#define RANDDBL(n) (drand48() * (n))
#define RANDSGN() ((RANDINT(2)==0)?(-1):(1))

/* number of bits in an integer */
#define INT_BITS		(sizeof(int)*8)

#define KROUND(dbl,prec) dbl >= 0.025 ? rint (dbl * 100.0) / 100.0 : rint (dbl * pow (10.0,prec)) / pow (10.0,prec)

/* Routines for checking/manipulating lines of file input */

int is_line_blank_or_comment (char *line);
int is_line_blank (char *line);
int is_line_comment (char *line);
char *get_nonblank_line (char *pLine, int maxLen, FILE * fp, int *pLineNo);
char *strlower (char *str);
char *fgetlongs (char **buff, int *bufflen, FILE * fp);
int permuteLine (char *line, int maxlength);

#endif /* __UTILS_H__ */
