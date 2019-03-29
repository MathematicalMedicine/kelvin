
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

#define MAXSWNAME 32
#define MAXSWMSG 220
#define MAXUDPMSG 230

struct swStopwatch
{
  char swName[MAXSWNAME + 1];
  struct rusage swStartRU;
  time_t swStartWallTime;
  struct rusage swAccumRU;
  time_t swAccumWallTime;
  int swStartedCount;
  int swRunning;
};

struct swStopwatch *swCreate (char *);
void swStart (struct swStopwatch *);
void swStop (struct swStopwatch *);
void swDump (struct swStopwatch *);
void swReset (struct swStopwatch *);

void swLogMsg (char *);
void udpSend (char *, int, char *);

void *swMalloc (size_t, char *, int);
void *swCalloc (size_t, size_t, char *, int);
void *swRealloc (void *, size_t, char *, int);
void swFree (void *, char *, int);
void swDumpBlockUse ();
void swAddChunk (void *, size_t, int, char *, int);
size_t swDelChunk (void *, int, char *, int);
void swDumpSources ();
void swDumpCrossModuleChunks ();
void swLogPeaks (char *);

#ifdef DMUSE
extern int used24s, used48s, used100s, missed24s, missed48s, missed100s;
#endif
#ifdef DMTRACK
extern double totalMalloc, totalFree, totalReallocOK, totalReallocMove,
  totalReallocFree, currentAlloc, peakAlloc;
extern int countMalloc, countFree, countReallocOK, countReallocMove,
  countReallocFree, maxListDepth, maxRecycles;
#endif

#endif
