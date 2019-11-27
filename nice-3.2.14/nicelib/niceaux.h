/**********************************************************************
 * NICE Utilities
 * Alberto Maria Segre
 * General Definitions, Software Timers and Message Logging System.
 * 
 * Copyright (c) 1996-2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software 
 * for non-profit educational purposes only.
 **********************************************************************/
#ifndef NOALARM
#include <setjmp.h>		/* Needed for setjmp */
#endif

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
 * Angle conversion degrees/radians.
 **********************************************************************/
#ifndef PI
#define PI 3.14159265358979323846
#define TWOPI 6.28318530717958647692
#endif
#define RADS(degrees) ((degrees)*(PI/180.0))
#define DEGS(radians) ((radians)*(180.0/PI))

/**********************************************************************
 * Beats pow(), usually (unless expressions are deep).
 **********************************************************************/
#define SQUARE(x) ((x)*(x))

/**********************************************************************
 * Generate a random integer between 0 and n-1 (inclusive).
 **********************************************************************/
#define SEED(n) (srandom (n));
#define RANDINT(n) (random() % (n))
#define SEEDDBL(n) (srand48(n));
#define RANDDBL(n) (drand48() * (n))
#define RANDSGN() ((RANDINT(2)==0)?(-1):(1))

#ifndef NOALARM
/**********************************************************************
 * Software timers, used in the niceapi and, optionally, also at the
 * application level. This library supports multiplexing of multiple
 * signal handlers on only one signal.  
 *
 * Two types of software timers are supported by the alarm mechanism:
 * an interrupt, which invokes a handler function and then continues
 * execution (after optionally rescheudling the interrupt), and a
 * timeout, which aborts the current thread of execution and unwinds
 * back to where the timeout was set (nonlocal exit).
 **********************************************************************/
/**********************************************************************
 * The generic alarm handler almost follows signal() semantics; the
 * only difference between alarmhandler and sighandler_t is that the
 * signal number is unnecessary -- it would always be the same.
 **********************************************************************/
typedef void (*alarmhandler) ();

/**********************************************************************
 * Set an interrupt. Returns an alarm identifier (a positive integer),
 * which can be used later to cancel the pending interrupt, or ERROR.
 *
 * If recurs is TRUE, the interrupt, upon firing, is immediately
 * rescheduled for the original time interval; otherwise, the
 * interrupt is a one-shot deal.
 *
 * If returns is TRUE, the alarm handler returns normally; otherwise,
 * the alarm handler performs a nonlocal exit (and thus requires
 * some sort of special handling on the part of the calling program).
 **********************************************************************/
int niceSetInterrupt (unsigned int time, alarmhandler handler, int recurs);

/**********************************************************************
 * Set (or reset) a timeout. Returns an alarm identifier (a positive
 * integer), which can be used later to cancel the pending timeout,
 * FALSE or ERROR.
 *
 * A timeout is an alarm that causes a nonlocal exit when it expires:
 * execution is unwound to where the timeout was set, and the
 * niceSetTimeout () function exits with a zero (FALSE) value.
 *
 * Once a timeout has been set and cleared, it can be reset using
 * niceResetTimeout (), which would unwind to the previous call to
 * niceSetTimeout() on expiration.
 **********************************************************************/
int niceSetTimeout (unsigned int time);
int niceResetTimeout (unsigned int time);

/**********************************************************************
 * Remove specified interrupt or timeout. Returns the cancelled alarm
 * identifier (a positive integer) or ERROR.
 **********************************************************************/
int niceClearAlarm (int id);

/**********************************************************************
 * Commonly used timeout value in order to declare a protocol is
 * kaput. 120 seconds ought to be OK for either
 * application-application or daemon-daemon interactions, although
 * certainly some special cases call for the shorter (e.g., AUTHINT)
 * or longer (e.g., TXFRINT) values set elsewhere.
 **********************************************************************/
#define HUNGINT 120	/* Standard hung protocol threshold, in sec. */
#endif

/**********************************************************************
 * Message logging system. Based in part on a preliminary
 * implementation by James Hunsaker.
 *
 * Level constants. Smaller numbers represent more critical messages,
 * with the 0 level "always on" corresponding to a fatal error.
 *
 * Originally, there were 7 levels (DEBUG, INFO, NOTICE, WARNING, ERR,
 * CRIT, FATAL) of messages, but most were just not used at all, so we
 * reduced it to 3 levels, plus one level that is "always on."
 **********************************************************************/
#define MAXLOGLEVEL 4
#define NLOGFAT 0	/* SILENT: Fatal error; abort. */
#define NLOGERR 1	/* TERSE: Serious error (e.g., security problem). */
#define NLOGWRN 2	/* CHATTY: Unusual stuff (e.g., network outages). */
#define NLOGMSG 3	/* NOISY: Normal stuff. */

/* Destination constants. These are used to set bits (they are powers
 * of 2) so that you may specify more than one destination at a time
 * for log message. */
#define NLOGNULL   0
#define NLOGSTDOUT 1
#define NLOGSTDERR 2
#ifndef WIN32
#define NLOGSYSLOG 4   /* Ignored in Windows! */
#endif

/* Function prototypes. */
void niceLogSet (int level, int destination);
void niceLog (int level, const char *format, ...);
