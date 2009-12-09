/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#ifndef SW_H
#define SW_H

#include <sys/resource.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAXSWNAME 32
#define MAXSWMSG 220
#define MAXUDPMSG 230

void pushStatus (char program, char *currentStatus);
void popStatus (char program);

struct swStopwatch
{
  char swName[MAXSWNAME + 1];
  struct rusage swStartRUSelf;
  struct rusage swStartRUChildren;
  time_t swStartWallTime;
  struct rusage swAccumRUSelf;
  struct rusage swAccumRUChildren;
  time_t swAccumWallTime;
  int swStartedCount;
  int swRunning;
};

struct swStopwatch *swCreate (char *);
void swStart (struct swStopwatch *);
void swStop (struct swStopwatch *);
void swDumpM (struct swStopwatch *);
void swDump (struct swStopwatch *);
void swReset (struct swStopwatch *);

void swLogMsg (FILE *, char *);
int udpSend (char *, int, char *);

void *swMalloc (size_t, char *, int);
void *swCalloc (size_t, int, char *, int);
void *swRealloc (void *, size_t, char *, int);
void swFree (void *, char *, int);
void swDumpHeldTotals (void);
void swDumpBlockUse (void);
void swAddChunk (void *, size_t, int, char *, int);
size_t swDelChunk (void *, int, char *, int);
void swDumpSources (void);
void swDumpCrossModuleChunks (void);
void swLogPeaks (char *);
long swGetMaximumPMK (void);
long swGetCurrentVMK (pid_t);

/// Maximum amount of physical memory available in Kbytes
extern long maximumPMK;

#ifdef DMTRACK
extern double totalMalloc, totalFree, totalReallocOK, totalReallocMove,
  totalReallocFree, currentAlloc, peakAlloc;
extern int countMalloc, countFree, countReallocOK, countReallocMove,
  countReallocFree, maxListDepth, maxRecycles;
#endif

/**

  Logging Support

  Logging consists of 6 different message types:

  The first two types are always logged immediately to stderr, and are prefaced by a timestamp and
  a keyword that allows easy identification:

  FATAL - A fatal error prefixed by "FATAL ERROR - ABORTING (<module>:<line no>), ". Doesn't return.
  WARNING - A warning. Prefixed by "WARNING (<module>:<line no>), ". Returns success.

  The next three types are progress advisories that may or may not be displayed on stdout 
  depending upon their significance and the volume of progress advisories requested by the user.
  Only the latest advisory in each of these three categories is kept. Two different directives
  are used to determine which advisories are displayed. The first, LogProgressDelayMinutes <minutes>,
  determines the delay between progress advisories. If it is set to zero (0), they are displayed
  as they occur. The second, LogProgressLevel [step|substep|detail], determines the lowest level
  displayed.

  STEP - A major user-viewpoint step such as start, finish and marker set change.
  SUBSTEP - A lesser user/analyst-viewpoint step such as position change.
  DETAIL - A low-level analyst-viewpoint step such as trait/marker/alternative likelihood 
  phase, xmission matrix build or polynomial load.

  The last type is DIAG. These are
  used in conjunction with the facility-specific diagnostic level to determine if they are
  displayed to stdout.

*/

#define MAXLOGMSG 2048
    
#define FATAL(...) \
{ \
  int length; \
  char message[MAXLOGMSG + 1], *pMessage = message; \
  pMessage += length = snprintf (message, MAXLOGMSG, "FATAL ERROR - ABORTING (%s:%d), ", (__FILE__),(__LINE__)); \
  snprintf (pMessage, MAXLOGMSG - length,  __VA_ARGS__);			\
  swLogMsg (stderr, message);						\
  exit (EXIT_FAILURE); \
}

#define WARNING(...) \
{ \
  int length; \
  char message[MAXLOGMSG + 1], *pMessage = message; \
  pMessage += length = snprintf (message, MAXLOGMSG, "WARNING (%s:%d), ", (__FILE__),(__LINE__)); \
  snprintf (pMessage, MAXLOGMSG - length,  __VA_ARGS__);			\
  swLogMsg (stderr, message);						\
}

/* Use these to report progress (indents with tabs by level, doesn't show percent done if 0),
   or just call swLogProgress directly. */
#define STEP(PERCENTDONE, ...) swLogProgress(0, PERCENTDONE, __VA_ARGS__)
#define SUBSTEP(PERCENTDONE, ...) swLogProgress(1, PERCENTDONE, __VA_ARGS__)
#define DETAIL(PERCENTDONE, ...) swLogProgress(2, PERCENTDONE, __VA_ARGS__)

#endif
