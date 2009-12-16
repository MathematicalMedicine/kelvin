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
#include <signal.h>

#define MAXSWNAME 32
#define MAXSWMSG 220
#define MAXUDPMSG 230

void swPushPhase (char program, char *currentPhase);
void swPopPhase (char program);

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

  The first two types are always logged immediately to stderr, and are
  prefaced by a timestamp and a keyword that allows easy
  identification:

  FATAL - An integrity error prefixed by "FATAL - ABORTING (<module>:<line no>), ". Doesn't 
  return. These errors cannot be fixed by the user.
  ERROR - A user-induced error prefixed by "ERROR - EXITING, ". Doesn't return. Usually
  can be fixed by the user.
  WARNING - A warning requiring user attention. Prefixed by "WARNING, ". Returns success.

  The next simply prints the message to stdout prefaced by a timestamp:

  INFO - A normal informational message to stdout. Returns success.

  The next type is a progress advisory that always displays to stdout.

  STEP - A major user-viewpoint step such as start, finish and marker set change.

  The next two types are optional progress advisories that may be
  displayed on stdout depending upon their timing and the volume of
  progress advisories requested by the user.  Only the most recent
  advisory in each of these these categories is kept. Two global
  variables are used to determine which advisories are
  displayed. logProgressDelaySeconds determines the delay between
  optional progress advisories. If it is set to zero (0), they are
  displayed as they occur. logProgressLevel determines the lowest
  optional level displayed, with 1 representing a SUBSTEP and 2
  DETAIL. Note that you can have even more excruciating detail by
  invoking swLogProgress without a macro specifying levels beyond 2.

  SUBSTEP - A lesser user/analyst-viewpoint step such as position change.
  DETAIL - A low-level analyst-viewpoint step such as trait/marker/alternative likelihood 
  phase, xmission matrix build or polynomial load.

  The last type is DIAG. These are used in conjunction with the
  facility-specific diagnostic level to determine if they are
  displayed to stdout. They are not normally present in the
  distributed code, as they are all compilation conditionals. After
  all, they can densely populate sensitive loops to the point where
  performance is affected even though they're "turned off".

*/

#define MAXLOGMSG 2048
    
#define FATAL(...) \
do { \
  int length; \
  char message[MAXLOGMSG + 1], *pMessage = message; \
  pMessage += length = snprintf (message, MAXLOGMSG, "FATAL - ABORTING (%s:%d), ", (__FILE__),(__LINE__)); \
  snprintf (pMessage, MAXLOGMSG - length,  __VA_ARGS__); \
  swLogMsg (stderr, message); \
  exit (EXIT_FAILURE); \
} while(0)

// The rest of these don't really need to be macros.
#define ERROR(...) \
do { \
  int length; \
  char message[MAXLOGMSG + 1], *pMessage = message; \
  pMessage += length = snprintf (message, MAXLOGMSG, "ERROR - EXITING, "); \
  snprintf (pMessage, MAXLOGMSG - length,  __VA_ARGS__); \
  swLogMsg (stderr, message); \
  exit (EXIT_FAILURE); \
} while(0)

#define ASSERT(CONDITION, ...) if (!(CONDITION)) ERROR(__VA_ARGS__)

#define WARNING(...) \
do { \
  int length; \
  char message[MAXLOGMSG + 1], *pMessage = message; \
  pMessage += length = snprintf (message, MAXLOGMSG, "WARNING, "); \
  snprintf (pMessage, MAXLOGMSG - length,  __VA_ARGS__); \
  swLogMsg (stderr, message); \
} while(0)

#define INFO(...) \
do { \
  char message[MAXLOGMSG + 1]; \
  snprintf (message, MAXLOGMSG,  __VA_ARGS__); \
  swLogMsg (stdout, message); \
} while(0)


/* Use these to report progress (indents with tabs by level, doesn't show percent done if 0),
   or just call swLogProgress directly. */

extern volatile sig_atomic_t swProgressRequestFlag;
extern int swProgressDelaySeconds;
extern int swProgressLevel;

char *formatElapsedTime (unsigned int t, char *buffer);
void swLogProgress(int level, float percentDone, char *format, ...);
void swStartProgressWakeUps(int seconds);
#define STEP(PERCENTDONE, ...) { swLogProgress(0, PERCENTDONE, __VA_ARGS__); }
#define SUBSTEP(PERCENTDONE, ...) { swLogProgress(1, PERCENTDONE, __VA_ARGS__); }
#define DETAIL(PERCENTDONE, ...) { swLogProgress(2, PERCENTDONE, __VA_ARGS__); }

#define OVERALL 0
#define XM 0

/* The beauty of this lines in the fact that all diags go away completely if DISTRIBUTION is
   defined, and an entire chunk of code can be the diagnostic. */
#ifndef DISTRIBUTION
#define DIAG(ENV_LEVEL, DIAG_LEVEL, DIAG_CODE) if (ENV_LEVEL >= DIAG_LEVEL) { DIAG_CODE }
#else
#define DIAG(ENV_LEVEL, DIAG_LEVEL, DIAG_CODE) // There was a diag here
#endif

#endif
