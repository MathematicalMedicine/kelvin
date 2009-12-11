/**
@file utils.h

  @see utils.c
  
  Copyright &copy; 2009, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

  @author Alberto Maria Segre
  @author Bill Valentine-Cooper

*/

#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdio.h>
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

#define MALCHOKE(pChunk,chunkSize,chunkCast)				\
  if ((pChunk = (chunkCast) malloc(chunkSize)) == 0) {			\
    fprintf (stderr,"malloc of %lu bytes for variable %s cast as (%s) failed at %s:%d!\n", \
	     (unsigned long)chunkSize,#pChunk,#chunkCast,(__FILE__),(__LINE__));	\
    exit (EXIT_FAILURE);						\
  }

#define REALCHOKE(pChunk,chunkSize,chunkCast)				\
  if ((pChunk = (chunkCast) realloc(pChunk, chunkSize)) == 0) {		\
    fprintf (stderr,"realloc to %lu bytes of variable %s cast as (%s) failed at %s:%d!\n", \
	     (unsigned long)chunkSize,#pChunk,#chunkCast,(__FILE__),(__LINE__));	\
    exit (EXIT_FAILURE);						\
  }

#define CALCHOKE(pChunk,chunkSize,chunkCount,chunkCast)			\
  if ((pChunk = (chunkCast) calloc(chunkSize, chunkCount)) == 0) {	\
    fprintf (stderr,"calloc of %lu * %lu zeroed bytes of variable %s cast as (%s) failed at %s:%d!\n", \
	     (unsigned long)chunkSize,(unsigned long)chunkCount,#pChunk,#chunkCast,(__FILE__),(__LINE__)); \
    exit (EXIT_FAILURE);						\
  }

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


/**********************************************************************
 * Kelvin log facility. 
 *
 * Levels are defined from 0 (most serious, and inevitably fatal) to
 * MAXLOGLEVELS (least serious). 
 *
 * Types of messages are defined as one of 32 orthogonal types.
 **********************************************************************/

/**********************************************************************
 * logFlag[] has to be global so the facility and severity can be
 * checked BEFORE invocation of logMsg so that parameters that are
 * function calls are not evaluated unnecessarily.
 **********************************************************************/
extern unsigned int logFlag[];
void logInit (void);
void logSet (unsigned int type, int level);
void logMsg (unsigned int type, int level, const char *format, ...);

/**********************************************************************
 * Levels of output verbosity: smaller is more severe. LOGFATAL errors
 * cause program termination, and are always produced on output.
 **********************************************************************/
#define MAXLOGLEVELS 5
#define LOGFATAL 0 /* report to user and abort */
#define LOGERROR 1 /* report to user but don't abort */
#define LOGWARNING 2 /* report to user */
#define LOGADVISE 3 /* allow filtered perusal by user */
#define LOGDEBUG 4 /* diagnostic runs only */

/**********************************************************************
 * Types of message output. We have 32 independent types of errors
 * allowed, each represented by a single bit in the 32 bit int
 * logFlag[].
 **********************************************************************/
#define LOGDEFAULT      1
#define LOGPEDFILE	(1 << 1)  /* processing pedigree file */
#define LOGSETRECODING	(1 << 2)  /* allele set recoding */
#define LOGGENOELIM	(1 << 3)  /* genotype elimination */
#define LOGPARENTALPAIR (1 << 4)  /* parental pair construction */
#define LOGPEELGRAPH    (1 << 5)  /* peel graph algorithm */
#define LOGLIKELIHOOD   (1 << 6)  /* likelihood caluclation */
#define LOGINPUTFILE    (1 << 7)  /* input configuration file */
#define LOGMEMORY	(1 << 8)  /* memory management */
#define LOGINTEGRATION	(1 << 9)  /* integration (dcuhre) */

/**********************************************************************
 * Macros for use in invoking the log function. These macro "wrappers"
 * provide for automatic inclusion of source filename and line number.
 **********************************************************************/

extern char *klog_prefix[]; ///< Prefixes allow message filtering
extern int klog_diagLevel[]; ///< Levels allow selective diagnostic detail per facility

/**********************************************************************
 * KLOG() reports an error of a given type at a given level. If the
 * level is LOGFATAL, will also abort.
 * 
 * Invocation:
 *   KLOG(type, level, formatString, args...)
 * No space allowed between KLOG and leading argument paren!
 **********************************************************************/
#define KLOGMSGLEN 2048
#define KLOG(FACILITY, LEVEL, ...) \
{ \
  char message[KLOGMSGLEN + 1], *pMessage = message; \
  pMessage += sprintf (message, "%s at %s:%d, ", klog_prefix[LEVEL], (__FILE__),(__LINE__)); \
  sprintf (pMessage, __VA_ARGS__); \
  if (LEVEL <= LOGWARNING) { \
    swLogMsg (stderr, message); \
    if (LEVEL == LOGFATAL) \
      exit (EXIT_FAILURE); \
  } else { \
    if (LEVEL == LOGADVISE) \
      swLogMsg (stdout, message); \
    else \
      if (LEVEL <= klog_diagLevel[FACILITY])	\
	  swLogMsg (stderr, message);		\
  } \
}

/**********************************************************************
 * KCHECK() checks a condition, and, if the condition fails, invokes
 * a TYPE error at level LEVEL. 
 * 
 * Invocation:
 *   KCHECK(condition, type, level, formatString, args...)
 * No space allowed between KCHECK and leading argument paren!
 **********************************************************************/
#define KCHECK(CONDITION, TYPE, LEVEL, ...)                           \
{                                                                     \
  if (!(CONDITION)) { \
    logMsg (TYPE, MAX(LOGERROR,LEVEL), "%s at %s:%d, ", klog_prefix[LEVEL], (__FILE__),(__LINE__)); \
    logMsg (TYPE, LEVEL, __VA_ARGS__);					\
  }									\
}

/**********************************************************************
 * KASSERT() checks a condition, and, if the condition fails, invokes
 * a LOGFATAL error and aborts. The TYPE of a KASSERT() is always
 * LOGDEFAULT, which is always turned on.
 * 
 * Invocation:
 *   KASSERT(condition, formatString, args...)
 * No space allowed between KASSERT and leading argument paren!
 **********************************************************************/
#define KASSERT(CONDITION, ...)                                       \
{                                                                     \
  if (!(CONDITION)) { \
    logMsg (LOGDEFAULT, LOGERROR, "%s at %s:%d, ", klog_prefix[LOGERROR], (__FILE__),(__LINE__)); \
    logMsg (LOGDEFAULT, LOGFATAL, __VA_ARGS__);				\
  }									\
}

#define KROUND(dbl) dbl >= 0.025 ? rint (dbl * 100.0) / 100.0 : rint (dbl * 10000.0) / 10000.0

/* Routines for checking/manipulating lines of file input */

int is_line_blank_or_comment (char *line);
int is_line_blank (char *line);
int is_line_comment (char *line);
char *get_nonblank_line (char *pLine, int maxLen, FILE * fp, int *pLineNo);
char *fgetlongs (char **buff, int *bufflen, FILE * fp);
int permuteLine (char *line, int maxlength);

#endif /* __UTILS_H__ */
