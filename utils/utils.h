
/**********************************************************************
 * Kelvin utilities.
 * Alberto Maria Segre
 *
 * Copyright 2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software
 * for non-profit educational purposes only.
 **********************************************************************/

/**********************************************************************
 * Some convenient and commonly used defines.
 **********************************************************************/

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
 * cause program termination, and are always produced on output. Note
 * that MAXLOGLEVELS corresponds to the max integer value assigned (0
 * doesn't count).
 **********************************************************************/
#define MAXLOGLEVELS 4
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

/* processing pedigree file */
#define LOGPEDFILE	(1 << 1)

/* allele set recoding */
#define LOGSETRECODING	(1 << 2)

/* genotype elimination */
#define LOGGENOELIM	(1 << 3)

/* parental pair construction */
#define LOGPARENTALPAIR (1 << 4)

/* peel graph algorithm */
#define LOGPEELGRAPH    (1 << 5)

/* likelihood caluclation */
#define LOGLIKELIHOOD   (1 << 6)

/* input configuration file */
#define LOGINPUTFILE    (1 << 7)

/* memory management */
#define LOGMEMORY	(1 << 8)

/* integration (dcuhre) */
#define LOGINTEGRATION	(1 << 9)

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
{ \
  if ((LEVEL == 0) || (TYPE & logFlag[LEVEL - 1])) { \
    logMsg (TYPE, MAX(LOGERROR,LEVEL), "%s (%d): ", (__FILE__),(__LINE__)); \
    logMsg (TYPE, LEVEL, __VA_ARGS__);                                  \
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
