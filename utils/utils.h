/**********************************************************************
 * Kelvin utilities.
 * Alberto Maria Segre
 *
 * Copyright 2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software
 * for non-profit educational purposes only.
 **********************************************************************/

#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdio.h>

/**********************************************************************
 * Some convenient and commonly used defines.
 **********************************************************************/

#define MALCHOKE(pChunk,chunkSize,chunkCast)				\
  if ((pChunk = (chunkCast) malloc(chunkSize)) == 0) {			\
    fprintf (stderr, "malloc of %lu bytes failed at %s:%d!\n",		\
	     chunkSize, (__FILE__),(__LINE__));				\
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
void logInit ();
void logSet (unsigned int type, int level);
void logMsg (unsigned int type, int level, const char *format, ...);

/**********************************************************************
 * Levels of output verbosity: smaller is more severe. LOGFATAL errors
 * cause program termination, and are always produced on output.
 **********************************************************************/
#define MAXLOGLEVELS 5
#define LOGFATAL 0
#define LOGERROR 1
#define LOGWARNING 2
#define LOGADVISE 3
#define LOGDEBUG 4

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

/**********************************************************************
 * KLOG() invokes an error of a given type at a given level. If the
 * level is LOGFATAL, will also abort.
 * 
 * Invocation:
 *   KLOG(type, level, formatString, args...)
 * No space allowed between KLOG and leading argument paren!
 **********************************************************************/
#define KLOG(TYPE, LEVEL, ...)                                        \
{                                                                     \
  if ((LEVEL == 0) || (TYPE & logFlag[LEVEL - 1])) { \
    logMsg (TYPE, MAX(LOGERROR,LEVEL), "%s (%d): ", (__FILE__),(__LINE__)); \
    logMsg (TYPE, LEVEL, __VA_ARGS__);                                \
  }                                                                   \
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
  if (!(CONDITION))   				                      \
   {                                                                  \
      logMsg (TYPE, MAX(LOGERROR,LEVEL), "%s (%d): ", (__FILE__),(__LINE__));\
      logMsg (TYPE, LEVEL, __VA_ARGS__);                              \
    }							              \
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
  if (!(CONDITION))   				                      \
   {                                                                  \
      logMsg (LOGDEFAULT, LOGERROR, "%s (%d): ", (__FILE__),(__LINE__));\
      logMsg (LOGDEFAULT, LOGFATAL, __VA_ARGS__);                     \
    }							              \
}


/* Wrappers for malloc() and friends, with built in error checking
 * and logging.
 */
void *MALLOC (char *description, size_t size);
void FREE (char *description, void *ptr);
void *REALLOC (char *description, void *ptr, size_t size);


/* Routines for checking/manipulating lines of file input */

int is_line_blank_or_comment (char *line);
int is_line_blank (char *line);
int is_line_comment (char *line);
char *get_nonblank_line (char *pLine, int maxLen, FILE * fp, int *pLineNo);
char *fgetlongs (char **buff, int *bufflen, FILE * fp);
int permuteLine (char *line, int maxlength);


#endif /* __UTILS_H__ */
