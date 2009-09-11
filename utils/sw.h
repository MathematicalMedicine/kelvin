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

void swLogMsg (char *);
int udpSend (char *, int, char *);

void *swMalloc (size_t, char *, int);
void *swCalloc (size_t, int, char *, int);
void *swRealloc (void *, size_t, char *, int);
void swFree (void *, char *, int);
void swDumpHeldTotals ();
void swDumpBlockUse ();
void swAddChunk (void *, size_t, int, char *, int);
size_t swDelChunk (void *, int, char *, int);
void swDumpSources ();
void swDumpCrossModuleChunks ();
void swLogPeaks (char *);
long swGetMaximumPMK ();
long swGetCurrentVMK (pid_t);

/// Maximum amount of physical memory available in Kbytes
extern long maximumPMK;

#ifdef DMTRACK
extern double totalMalloc, totalFree, totalReallocOK, totalReallocMove,
  totalReallocFree, currentAlloc, peakAlloc;
extern int countMalloc, countFree, countReallocOK, countReallocMove,
  countReallocFree, maxListDepth, maxRecycles;
#endif

#endif
