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

/*
$Id$

by Bill Valentine-Cooper for the Center for Quantitative and Computational
Biology, Nationwide Children's Hospital Research Institute.

Provides stopwatches for tracking code performance. Different named stopwatches
can be started and stopped at the beginning and end of sections of your code so
as to track the total elapsed time, CPU time, context switches and paging
for that portion of your code. These values can be compared to the overall 
resource utilization of your code in order to identify sections that are 
performing abnormally.

Uses getrusage. Linux 2.6 populates only the following rusage structure members:
utime, stime, nvcsw, nivcsw, minflt, and majflt.

Usage:

struct swStopwatch *pSW = swCreate("test"); to create a new stopwatch
named "test".

swStart(pSW); to start the stopwatch running.

swStop(pSW); to stop (really just pause) the stopwatch.

swDump(pSW); display stopwatch information on stdout. Looks like:

Stopwatch overall(9) e:74s u:22.038650s s:32.414072s, vx:14, ivx:2636, sf:5396, hf:0

where "overall" is the name of the stopwatch specified in the call to swCreate,
the (9) means the stopwatch was started/stopped 9 times, and the other statistics
are:
e: elapsed time in seconds.
u: time spent executing user instructions.
s: time spent in operating system code on behalf of process.
vx: number of voluntary context switches (usually for some service)
ivx: number of involuntary context switches (time slice expiry or pre-emption)
sf: number of soft pagefaults (no actual IO)
hf: number of hard pagefaults (actual IO performed)

swReset(pSW); to reset the stopwatch to all zeroes.

You can reference and use the members of the structure pointer returned by
swCreate if you want to write your own output instead of relying upon swDump.

To do dynamic memory tracking, build with -DDMTRACK and put the following macros 
into the modules to be affected:

#ifdef DMTRACK
#warning "Dynamic memory usage dumping is turned on, so performance will be poor!"
#endif
#if defined (DMTRACK)
#define malloc(X) swMalloc((X), __FILE__, __LINE__)
#define calloc(X,Y) swCalloc((X),(Y), __FILE__, __LINE__)
#define realloc(X,Y) swRealloc((X),(Y), __FILE__, __LINE__)
#define free(X) swFree((X), __FILE__, __LINE__)
#endif

DO NOT USE THESE IN SW.C because you'll get recurive looping.

The test main() serves as an example of use. To build a standalone test 
version:

$ gcc -lm -o swTest sw.c -lpthread -lm -DMAIN
$ ./swTest

To build an object you can link into your program:

$ gcc -c sw.c

then #include sw.h in your source code, and link with sw.o.

*/

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <ctype.h>
#include <time.h>
#include <stdarg.h>
#include <pthread.h>
#ifdef TELLRITA
#include <netdb.h>
#include <netinet/in.h>
#endif
#include <sys/ipc.h>
#include <sys/shm.h>
#include <errno.h>
#include "sw.h"
#include "hashtab.h"
#include "lookupa.h"
#ifdef BACKTRACE
#include <execinfo.h> // For stack dump on demand
#endif

#ifndef MAX
#define MAX(m,n) ((m)>=(n)?(m):(n))
#endif
#ifndef MIN
#define MIN(m,n) ((m)>=(n)?(n):(m))
#endif

long currentVMK, maximumPMK = -1;


#ifdef linux
#include <linux/prctl.h>
int prctl (int, char *); // Strangely this is not included!
#endif

void swDumpStack (void) {
#ifdef BACKTRACE
  const size_t MAXSTACKDEPTH = 20;
  size_t stack_depth;
  void *stack_addrs[MAXSTACKDEPTH];
  char **stack_strings;

  stack_depth = backtrace(stack_addrs, MAXSTACKDEPTH);
  stack_strings = backtrace_symbols(stack_addrs, stack_depth);

  fprintf (stderr, "\nCall stack: ");
  size_t i;
  for (i = 2; i < stack_depth; i++) {
    char *wordBoundary;
    if ((wordBoundary = strchr( stack_strings[i], '+')) != '\0') {
      while ((isalnum(*wordBoundary) == 0) && (*wordBoundary != '_'))
	*(wordBoundary--) = '\0';
      while (isalnum(*wordBoundary) || (*wordBoundary == '_'))
	wordBoundary--;
      wordBoundary++;
      fprintf(stderr, "%s<-", wordBoundary);
      if (strcmp(wordBoundary, "main") == 0)
	break;
    }
  }
  fprintf (stderr, "\n");
  free(stack_strings); // malloc()ed by backtrace_symbols                                                             
  fflush(stderr);
#endif
}

/**

Show the current state of the process in some obvious manner -- currently
(on linux) set the process name.

*/
#define MAXPHASELEN 13
#define PHASE_STACK_DEPTH 32
char phaseStack[PHASE_STACK_DEPTH][MAXPHASELEN+1];
int phaseStackPosition = 0;
void swPushPhase (char program, char *currentPhase)
{
  char processName[16+1];
  if (phaseStackPosition == PHASE_STACK_DEPTH) {
    DIAG(0, 0, {
	WARNING("Phase stack overflow (not serious), phase not changed to [%s]", currentPhase);
	int i;
	for (i=0; i<PHASE_STACK_DEPTH; i++)
	  fprintf (stderr, "%d is %s\n", i, phaseStack[i]);
	fflush (stderr);
	});
    return;
  }
  strncpy (phaseStack[++phaseStackPosition], currentPhase, (size_t) MAXPHASELEN);
  sprintf (processName, "%c-%-.*s", program, MAXPHASELEN, currentPhase);
#ifdef PR_SET_NAME
  prctl (PR_SET_NAME, processName);
#endif
  return;
}
void swPopPhase (char program)
{
  char processName[16+1];
  if (phaseStackPosition == 0) {
    DIAG(0, 0, {
	WARNING("Phase stack underflow (not serious), phase not reverted");
      });
    return;
  }
  sprintf (processName, "%c-%-.*s", program, MAXPHASELEN, phaseStack[--phaseStackPosition]);
#ifdef PR_SET_NAME
  prctl (PR_SET_NAME, processName);
#endif
  return;
}

/* All we do to create a new stopwatch is allocate space for the structure, 
   zero the accumulated information and return a pointer to it. */
struct swStopwatch *
swCreate (char *swName)
{
  struct swStopwatch *newStopwatch;

  if ((newStopwatch =
       (struct swStopwatch *) malloc (sizeof (struct swStopwatch))) == 0) {
    fprintf (stderr,
	     "Unable to allocate memory for another stopwatch (%s)\n",
	     swName);
    fflush (stderr);
    exit (EXIT_FAILURE);
  } else {
    memset (newStopwatch, 0, sizeof (struct swStopwatch));
    strncpy (newStopwatch->swName, swName, (size_t) MAXSWNAME);
  }
  return newStopwatch;
}

/* Starting a stopwatch consists of 1) initializing the swStart structure 
   using a getrusage call for resource information and a time call for 
   elapsed time, and 2) turning on the swRunning flag. */
void
swStart (struct swStopwatch *theStopwatch)
{
  /* We need to tolerate multiple starts for recursive routines. */
  if (theStopwatch->swRunning)
    return;
  if (getrusage (RUSAGE_SELF, &theStopwatch->swStartRUSelf) != 0) {
    fprintf (stderr,
	     "Unable to get resource information for self using rusage call (%s)\n",
	     theStopwatch->swName);
    fflush (stderr);
    exit (EXIT_FAILURE);
  }
  if (getrusage (RUSAGE_CHILDREN, &theStopwatch->swStartRUChildren) != 0) {
    fprintf (stderr,
	     "Unable to get resource information for children using rusage call (%s)\n",
	     theStopwatch->swName);
    fflush (stderr);
    exit (EXIT_FAILURE);
  }
  theStopwatch->swStartWallTime = time (NULL);
  theStopwatch->swStartedCount += 1;
  theStopwatch->swRunning = 1;
  return;
}

/* Stopping a stopwatch consists of 1) getting current rusage and time info
   and using it to calculate deltas from swStart values, 2) bumping the 
   swAccumulated statistics with those deltas, and 3) turning off the
   swRunning flag. */
void
swStop (struct swStopwatch *theStopwatch)
{

  struct rusage stopStatsSelf, stopStatsChildren;
  time_t stopTime;
  unsigned long usec;

  if (!theStopwatch->swRunning)
    return;
  if (getrusage (RUSAGE_SELF, &stopStatsSelf) != 0) {
    fprintf (stderr,
	     "Unable to get resource information for self using rusage call (%s)\n",
	     theStopwatch->swName);
    fflush (stderr);
    exit (EXIT_FAILURE);
  }
  if (getrusage (RUSAGE_CHILDREN, &stopStatsChildren) != 0) {
    fprintf (stderr,
	     "Unable to get resource information for children using rusage call (%s)\n",
	     theStopwatch->swName);
    fflush (stderr);
    exit (EXIT_FAILURE);
  }
  stopTime = time (NULL);
  theStopwatch->swAccumWallTime +=
    difftime (stopTime, theStopwatch->swStartWallTime);

  usec =
    1000000 + theStopwatch->swAccumRUSelf.ru_utime.tv_usec +
    stopStatsSelf.ru_utime.tv_usec - theStopwatch->swStartRUSelf.ru_utime.tv_usec;
  theStopwatch->swAccumRUSelf.ru_utime.tv_sec +=
    stopStatsSelf.ru_utime.tv_sec - theStopwatch->swStartRUSelf.ru_utime.tv_sec +
    (usec / 1000000) - 1;
  theStopwatch->swAccumRUSelf.ru_utime.tv_usec = usec % 1000000;

  usec = 1000000 + theStopwatch->swAccumRUSelf.ru_stime.tv_usec +
    stopStatsSelf.ru_stime.tv_usec - theStopwatch->swStartRUSelf.ru_stime.tv_usec;
  theStopwatch->swAccumRUSelf.ru_stime.tv_sec += stopStatsSelf.ru_stime.tv_sec -
    theStopwatch->swStartRUSelf.ru_stime.tv_sec + (usec / 1000000) - 1;
  theStopwatch->swAccumRUSelf.ru_stime.tv_usec = usec % 1000000;

  theStopwatch->swAccumRUSelf.ru_nvcsw +=
    stopStatsSelf.ru_nvcsw - theStopwatch->swStartRUSelf.ru_nvcsw;
  theStopwatch->swAccumRUSelf.ru_nivcsw +=
    stopStatsSelf.ru_nivcsw - theStopwatch->swStartRUSelf.ru_nivcsw;
  theStopwatch->swAccumRUSelf.ru_minflt +=
    stopStatsSelf.ru_minflt - theStopwatch->swStartRUSelf.ru_minflt;
  theStopwatch->swAccumRUSelf.ru_majflt +=
    stopStatsSelf.ru_majflt - theStopwatch->swStartRUSelf.ru_majflt;

  usec =
    1000000 + theStopwatch->swAccumRUChildren.ru_utime.tv_usec +
    stopStatsChildren.ru_utime.tv_usec - theStopwatch->swStartRUChildren.ru_utime.tv_usec;
  theStopwatch->swAccumRUChildren.ru_utime.tv_sec +=
    stopStatsChildren.ru_utime.tv_sec - theStopwatch->swStartRUChildren.ru_utime.tv_sec +
    (usec / 1000000) - 1;
  theStopwatch->swAccumRUChildren.ru_utime.tv_usec = usec % 1000000;

  usec = 1000000 + theStopwatch->swAccumRUChildren.ru_stime.tv_usec +
    stopStatsChildren.ru_stime.tv_usec - theStopwatch->swStartRUChildren.ru_stime.tv_usec;
  theStopwatch->swAccumRUChildren.ru_stime.tv_sec += stopStatsChildren.ru_stime.tv_sec -
    theStopwatch->swStartRUChildren.ru_stime.tv_sec + (usec / 1000000) - 1;
  theStopwatch->swAccumRUChildren.ru_stime.tv_usec = usec % 1000000;

  theStopwatch->swAccumRUChildren.ru_nvcsw +=
    stopStatsChildren.ru_nvcsw - theStopwatch->swStartRUChildren.ru_nvcsw;
  theStopwatch->swAccumRUChildren.ru_nivcsw +=
    stopStatsChildren.ru_nivcsw - theStopwatch->swStartRUChildren.ru_nivcsw;
  theStopwatch->swAccumRUChildren.ru_minflt +=
    stopStatsChildren.ru_minflt - theStopwatch->swStartRUChildren.ru_minflt;
  theStopwatch->swAccumRUChildren.ru_majflt +=
    stopStatsChildren.ru_majflt - theStopwatch->swStartRUChildren.ru_majflt;

  theStopwatch->swRunning = 0;
  return;
}

void
swDumpOutput (struct swStopwatch *theStopwatch, char *appendText)
{
  char buffer[MAXUDPMSG];
  sprintf
    (buffer,
     "stopwatch %s(%d) e:%lus u:%lus s:%lus, vx:%lu, ivx:%lu, sf:%lu, hf:%lu%s",
     theStopwatch->swName, theStopwatch->swStartedCount,
     (unsigned long) theStopwatch->swAccumWallTime,
     (unsigned long) theStopwatch->swAccumRUSelf.ru_utime.tv_sec + theStopwatch->swAccumRUChildren.ru_utime.tv_sec,
     (unsigned long) theStopwatch->swAccumRUSelf.ru_stime.tv_sec + theStopwatch->swAccumRUChildren.ru_stime.tv_sec,
     (unsigned long) theStopwatch->swAccumRUSelf.ru_nvcsw + theStopwatch->swAccumRUChildren.ru_nvcsw,
     (unsigned long) theStopwatch->swAccumRUSelf.ru_nivcsw + theStopwatch->swAccumRUChildren.ru_nivcsw,
     (unsigned long) theStopwatch->swAccumRUSelf.ru_minflt + theStopwatch->swAccumRUChildren.ru_minflt,
     (unsigned long) theStopwatch->swAccumRUSelf.ru_majflt + theStopwatch->swAccumRUChildren.ru_majflt,
     appendText);
  INFO("%s", buffer);
  return;
}

void
swDumpM (struct swStopwatch *theStopwatch)
{
  char memoryBuffer[MAXUDPMSG];

  if (maximumPMK == -1) {
    maximumPMK = swGetMaximumPMK ();
  }
  if (maximumPMK != 0) {
    currentVMK = swGetCurrentVMK (getpid ());
    sprintf (memoryBuffer, ", vmG:%2.1f/%2.1f",
	     currentVMK / (1024.0 * 1024.0),
	     maximumPMK / (1024.0 * 1024.0));
  } else
    memoryBuffer[0] = '\0';

  if (theStopwatch->swRunning) {
    swStop (theStopwatch);
    swDumpOutput (theStopwatch, memoryBuffer);
    swStart (theStopwatch);
  } else {
    swDumpOutput (theStopwatch, memoryBuffer);
  }
  return;
}

void
swDump (struct swStopwatch *theStopwatch)
{
  if (theStopwatch->swRunning) {
    swStop (theStopwatch);
    swDumpOutput (theStopwatch, "");
    swStart (theStopwatch);
  } else {
    swDumpOutput (theStopwatch, "");
  }
  return;
}

/* Extremely lazy, please ignore. */
void
swReset (struct swStopwatch *theStopwatch)
{
  char swName[MAXSWNAME + 1];

  strncpy (swName, theStopwatch->swName, (size_t) MAXSWNAME);
  memset (theStopwatch, 0, sizeof (struct swStopwatch));
  strncpy (theStopwatch->swName, swName, (size_t) MAXSWNAME);
  return;
}

/* Dynamic memory debugging. */

struct swStopwatch *internalDMSW;

/* The proper way to do this would be another hash, but there are very
   few calls and even fewer modules, so a simple ordered list and a
   binary search should suffice. */
enum callTypes
  { cTMalloc, ctCalloc, cTReAlloc, cTReFree, cTFree };
#ifdef DMTRACK
static char *callTypeNames[] = { "malloc", "calloc", "realloc", "realloc", "free" };
#endif

#define MAXMEMCHUNKSOURCECOUNT 2048
struct memChunkSource
{
  unsigned char callType;	/* malloc, calloc, realloc, or free */
  char moduleName[255];		/* Calling module name */
  int entryNo;			/* The order created for linking to memChunks */
  int lineNo;			/* Calling module line number */
  int totalCalls;
  size_t totalBytes;
  size_t remainingBytes;
};
struct memChunkSource  memChunkSources[MAXMEMCHUNKSOURCECOUNT];
int memChunkSourceCount = 0;

struct memChunk
{
  void *chunkAddress;
  size_t chunkSize;
  int hashKey;
  int recycleCount;
  short allocSource;		/* Entry in memChunkSource list where this was allocated */
  short freeSource;		/* Entry in memChunkSource list where this was freed */
  struct memChunk *next;
};

#ifdef DMTRACK
#define DMHASHSIZE 27
htab *chunkHash;

/* Up to 1K in bytes */
#define MAXSMALLBLOCK 1024

/* Up to 1M in Kbytes */
#define MAXMEDIUMBLOCK 1024

/* Up to 2G in 1Mbytes */
#define MAXLARGEBLOCK 2*1024-1

/* 0 is for malloc, 1 is for realloc, 2 is for realloc failures. */
int smallBlocks[MAXSMALLBLOCK][3];
int mediumBlocks[MAXMEDIUMBLOCK][3];
int largeBlocks[MAXLARGEBLOCK][3];
double totalMalloc = 0, totalReallocOK = 0, totalReallocMove = 0,
  totalReallocFree = 0, totalFree = 0, currentAlloc = 0, peakAlloc = 0;
int countMalloc = 0, countReallocOK = 0, countReallocMove = 0,
  countReallocFree = 0, countFree = 0, maxListDepth = 0, maxRecycles = 0;
#endif
int firstMallocCall = 1;

void
swFirstDM (void)
{
  internalDMSW = swCreate ("internalDMSW");
  swStart (internalDMSW);
  firstMallocCall = 0;

#ifdef DMTRACK
  chunkHash = hcreate (DMHASHSIZE);
#endif

  return;
}

#ifdef DMTRACK

int
compareSourcesByName (const void *left, const void *right)
{
  int result;
  struct memChunkSource *mCSLeft, *mCSRight;

  mCSLeft = (struct memChunkSource *) left;
  mCSRight = (struct memChunkSource *) right;
  if ((result = strcmp (mCSLeft->moduleName, mCSRight->moduleName)) == 0) {
    result = mCSLeft->lineNo - mCSRight->lineNo;
  }
  return result;
}

int
compareSourcesByEntryNo (const void *left, const void *right)
{
  struct memChunkSource *mCSLeft, *mCSRight;

  mCSLeft = (struct memChunkSource *) left;
  mCSRight = (struct memChunkSource *) right;
  return (mCSLeft->entryNo - mCSRight->entryNo);
}

int
compareSourcesByTotalBytes (const void *left, const void *right)
{
  struct memChunkSource *mCSLeft, *mCSRight;

  mCSLeft = (struct memChunkSource *) left;
  mCSRight = (struct memChunkSource *) right;
  if (labs (mCSRight->totalBytes) - labs (mCSLeft->totalBytes) > 0)
    return 1;
  if (labs (mCSRight->totalBytes) - labs (mCSLeft->totalBytes) < 0)
    return -1;
  return 0;
}

int
compareSourcesByRemainingBytes (const void *left, const void *right)
{
  struct memChunkSource *mCSLeft, *mCSRight;

  mCSLeft = (struct memChunkSource *) left;
  mCSRight = (struct memChunkSource *) right;
  if (labs (mCSRight->remainingBytes) - labs (mCSLeft->remainingBytes) > 0)
    return 1;
  if (labs (mCSRight->remainingBytes) - labs (mCSLeft->remainingBytes) < 0)
    return -1;
  return 0;
}

/* Find it's source or add a new one. */
short
findOrAddSource (char *fileName, int lineNo, int callType, size_t chunkSize)
{
  struct memChunkSource target, *result;

  strcpy (target.moduleName, fileName);
  target.lineNo = lineNo;
  result =
    bsearch (&target, memChunkSources, memChunkSourceCount,
	     sizeof (struct memChunkSource), compareSourcesByName);
  if (result) {
    result->totalCalls++;
    result->totalBytes += chunkSize;
    if (callType != cTFree && callType != cTReFree)
      result->remainingBytes += chunkSize;
    return (result->entryNo);
  } else {
    if (memChunkSourceCount >= MAXMEMCHUNKSOURCECOUNT) {
      fprintf (stderr,
	       "Exceeded maximum chunkSourceCount, no more locations can be monitored\n");
      fflush (stderr);
      return 0;
    }
    if (callType != cTFree && callType != cTReFree)
      memChunkSources[memChunkSourceCount].remainingBytes = chunkSize;
    strcpy (memChunkSources[memChunkSourceCount].moduleName, fileName);
    memChunkSources[memChunkSourceCount].lineNo = lineNo;
    memChunkSources[memChunkSourceCount].callType = callType;
    memChunkSources[memChunkSourceCount].entryNo = memChunkSourceCount;
    memChunkSources[memChunkSourceCount].totalCalls = 1;
    memChunkSources[memChunkSourceCount].totalBytes = chunkSize;
    qsort (memChunkSources, ++memChunkSourceCount,
	   sizeof (struct memChunkSource), compareSourcesByName);
    return (memChunkSourceCount - 1);
  }
}

void
swAddChunk (void *chunkAddress, size_t chunkSize, int callType,
	    char *fileName, int lineNo)
{
  struct memChunk *newChunk, *oldChunk;
  int currentListDepth = 0;
  newChunk = (struct memChunk *) malloc (sizeof (struct memChunk));
  newChunk->chunkSize = chunkSize;
  newChunk->chunkAddress = chunkAddress;
  newChunk->recycleCount = 1;
  newChunk->allocSource =
    findOrAddSource (fileName, lineNo, callType, chunkSize);
  newChunk->next = NULL;
  newChunk->hashKey = lookup ((ub1 *) &chunkAddress, sizeof (chunkAddress), 0);
  if (hadd (chunkHash, (ub1 *) &newChunk->hashKey, sizeof (newChunk->hashKey),
       newChunk) == FALSE) {
    /* Collision - see if we have the address anywhere. */
    oldChunk = (struct memChunk *) hstuff (chunkHash);
    while (oldChunk->chunkAddress != chunkAddress) {
      if (oldChunk->next != NULL) {
	oldChunk = oldChunk->next;
	currentListDepth++;
	if (currentListDepth > maxListDepth) {
	  maxListDepth = currentListDepth;
	}
      } else {
	oldChunk->next = newChunk;
	return;
      }
    }
    if ((oldChunk->recycleCount % 2) != 0) {
      /*
         fprintf (stderr,
         "ERROR!! (%s: %d) - alloc of %u for size %u is in use for size %u, recycled %d times!\n",
         fileName, lineNo, chunkAddress, chunkSize, oldChunk->chunkSize,
         oldChunk->recycleCount);
       */
      oldChunk->recycleCount++;
    }
    oldChunk->recycleCount++;
    if (oldChunk->recycleCount > maxRecycles) {
      maxRecycles = oldChunk->recycleCount;
    }
    oldChunk->chunkSize = newChunk->chunkSize;
    oldChunk->allocSource = newChunk->allocSource;
    free (newChunk);
  }
  return;
}

size_t
swDelChunk (void *chunkAddress, int callType, char *fileName, int lineNo)
{
  struct memChunk *oldChunk;
  int hashKey;

  hashKey = lookup ((ub1 *) &chunkAddress, sizeof (chunkAddress), 0);
  if (hfind (chunkHash, (ub1 *) &hashKey, sizeof (hashKey)) == FALSE) {
    DIAG (0, 0, {
	WARNING("free() of address %lu that was never allocated (head)",
		(unsigned long) chunkAddress);
      });
    return 0;
  } else {
    oldChunk = (struct memChunk *) hstuff (chunkHash);
    while (oldChunk->chunkAddress != chunkAddress) {
      if (oldChunk->next != NULL) {
	oldChunk = oldChunk->next;
      } else {
	DIAG (0, 0, {
	    WARNING("free() of address %lu that was never allocated (not head)",
		    (unsigned long) chunkAddress);
	  });
	return 0;
      }
    }
    if ((oldChunk->recycleCount % 2) == 0) {
      DIAG (0, 0, {
	  ERROR("free() of address %lu that is no longer in use for size %u and recycled %d times",
		(unsigned long) chunkAddress, (unsigned int) oldChunk->chunkSize,
		oldChunk->recycleCount);
	});
      oldChunk->recycleCount++;
    } else {
      int i;
      for (i=0; i<memChunkSourceCount; i++)
        if (memChunkSources[i].entryNo == oldChunk->allocSource) {
	  DIAG (0, 3, {
	      if (callType == cTReAlloc || callType == cTReFree) {
		fprintf(stderr, "For %s: %d, subtracting %d from %d\n",
			fileName, lineNo, (int) oldChunk->chunkSize, (int) memChunkSources[i].remainingBytes);
		fflush (stderr);
	      }
	    });
	  memChunkSources[i].remainingBytes -= oldChunk->chunkSize;
	  break;
        }
    }
    oldChunk->freeSource =
      findOrAddSource (fileName, lineNo, callType, -oldChunk->chunkSize);
    oldChunk->recycleCount++;
    if (oldChunk->recycleCount > maxRecycles) {
      maxRecycles = oldChunk->recycleCount;
    }
    return (oldChunk->chunkSize);
  }
}

void
swDumpSources ()
{
  int i;

  fprintf (stderr, "Top 100 sources out of %d total:\n",
	   memChunkSourceCount);
  qsort (memChunkSources, memChunkSourceCount, sizeof (struct memChunkSource),
	 compareSourcesByTotalBytes);
  for (i=0; i<min(memChunkSourceCount, 100); i++) {
    fprintf (stderr, "At %s line %d, %s called %d times for %ld bytes\n",
	     memChunkSources[i].moduleName, memChunkSources[i].lineNo,
	     callTypeNames[memChunkSources[i].callType],
	     memChunkSources[i].totalCalls, memChunkSources[i].totalBytes);
  }
  fflush (stderr);
  qsort (memChunkSources, ++memChunkSourceCount,
         sizeof (struct memChunkSource), compareSourcesByName);
}

void
swDumpHeldTotals ()
{
  int i;
  fprintf (stderr, "Top 20 still-held (net) sources out of %d total:\n",
	   memChunkSourceCount);
  qsort (memChunkSources, memChunkSourceCount, sizeof (struct memChunkSource),
	 compareSourcesByRemainingBytes);
  for (i=0; i<min(memChunkSourceCount, 20); i++) {
    fprintf (stderr, "At %s line %d, %s called for net %ld bytes\n",
	     memChunkSources[i].moduleName, memChunkSources[i].lineNo,
	     callTypeNames[memChunkSources[i].callType],
	     memChunkSources[i].remainingBytes);
  }
  fflush (stderr);
  qsort (memChunkSources, ++memChunkSourceCount,
         sizeof (struct memChunkSource), compareSourcesByName);
}

void
swDumpCrossModuleChunks ()
{
  /* Need to watch that we only do freed blocks, otherwise we could pair a fresh alloc with an
     old free and send people off on a wild goose chase. */
  struct memChunk *chunk = NULL;
  struct memChunkSource target, *allocResult, *freeResult;

  fprintf (stderr, "Finding any cross-module alloc/free usage...\n");
  /* Change the sort order of memChunkSources to expedite cross-module checking */
  qsort (memChunkSources, memChunkSourceCount, sizeof (struct memChunkSource),
	 compareSourcesByEntryNo);
  /* Now traverse our chunkHash comparing alloc and free source modules */
  if (hfirst (chunkHash))
    do {
      chunk = (struct memChunk *) hstuff (chunkHash);
      if ((chunk->recycleCount % 2) == 0) {
	target.entryNo = chunk->allocSource;
	allocResult =
	  bsearch (&target, memChunkSources, memChunkSourceCount,
		   sizeof (struct memChunkSource), compareSourcesByEntryNo);
	if (allocResult) {
	  target.entryNo = chunk->freeSource;
	  freeResult =
	    bsearch (&target, memChunkSources, memChunkSourceCount,
		     sizeof (struct memChunkSource), compareSourcesByEntryNo);
	  if (freeResult) {
	    if (strcmp (allocResult->moduleName, freeResult->moduleName)) {
	      fprintf (stderr,
		       "Block of size %lu allocated in %s at line %d and freed in %s at line %d!\n",
		       chunk->chunkSize, allocResult->moduleName,
		       allocResult->lineNo, freeResult->moduleName,
		       freeResult->lineNo);
	    }
	  }
	}
      }
    } while (hnext (chunkHash));
  fflush (stderr);
  qsort (memChunkSources, ++memChunkSourceCount,
         sizeof (struct memChunkSource), compareSourcesByName);
}

void
swLogPeaks (char *reason)
{
  char messageBuffer[MAXSWMSG];
  if (firstMallocCall) {
    fprintf (stderr, "=> (%s): Before first memory allocation\n", reason);
    fflush (stderr);
    return;
  }
  swStop (internalDMSW);
  sprintf (messageBuffer,
	   "=> (%s): %lu seconds since 1st allocation, %g bytes in use, peak was %g\n",
	   reason, internalDMSW->swAccumWallTime, currentAlloc, peakAlloc);
  sprintf (messageBuffer,
	   "Count malloc:%d, free:%d, realloc OK:%d, realloc move:%d, realloc free:%d, max depth:%d, max recycles:%d",
	   countMalloc, countFree, countReallocOK, countReallocMove,
	   countReallocFree, maxListDepth, maxRecycles);
  INFO(messageBuffer);
  sprintf (messageBuffer,
	   "Size malloc:%g, free:%g, realloc OK:%g, realloc move:%g, realloc free:%g, current:%g, peak:%g",
	   totalMalloc, totalFree, totalReallocOK, totalReallocMove,
	   totalReallocFree, currentAlloc, peakAlloc);
  INFO(messageBuffer);
  swStart (internalDMSW);
  return;
}

void
swDumpChunks ()
{
  struct memChunk *newChunk;
  int i;

  if (hfirst (chunkHash))
    do {
      ub1 *hashKey;

      newChunk = (struct memChunk *) hstuff (chunkHash);
      hashKey = hkey (chunkHash);
      printf ("key ");
      for (i = 0; i < 4; i++) {
	printf ("%d ", hashKey[i]);
      }
      printf ("for address %lu and size %lu\n",
	      (unsigned long) newChunk->chunkAddress, newChunk->chunkSize);
    }
    while (hnext (chunkHash));
}

void
swDumpBlockUse ()
{
  char messageBuffer[MAXSWMSG];
  int i;

  sprintf (messageBuffer,
	   "Count malloc: %d, free: %d, realloc OK: %d, realloc move: %d, realloc free: %d, max depth: %d",
	   countMalloc, countFree, countReallocOK, countReallocMove,
	   countReallocFree, maxListDepth);
  INFO(messageBuffer);
  sprintf (messageBuffer,
	   "Size malloc: %g, free: %g, realloc OK: %g, realloc move: %g, realloc free: %g, current: %g, peak: %g",
	   totalMalloc, totalFree, totalReallocOK, totalReallocMove,
	   totalReallocFree, currentAlloc, peakAlloc);
  INFO(messageBuffer);

  fprintf (stderr, "Memory blocks allocated:\n    Size      Count\n");
  for (i = 0; i < MAXSMALLBLOCK; i++)
    if (smallBlocks[i][0] != 0)
      fprintf (stderr, "%8db  %8d\n", i, smallBlocks[i][0]);
  for (i = 0; i < MAXMEDIUMBLOCK; i++)
    if (mediumBlocks[i][0] != 0)
      fprintf (stderr, "%8dKb %8d\n", i, mediumBlocks[i][0]);
  for (i = 0; i < MAXLARGEBLOCK; i++)
    if (largeBlocks[i][0] != 0)
      fprintf (stderr, "%8dMb %8d\n", i, largeBlocks[i][0]);
  fprintf (stderr,
	   "\nMemory blocks reallocated in-place:\n    Size      Count\n");
  for (i = 0; i < MAXSMALLBLOCK; i++)
    if (smallBlocks[i][1] != 0)
      fprintf (stderr, "%8db  %8d\n", i, smallBlocks[i][1]);
  for (i = 0; i < MAXMEDIUMBLOCK; i++)
    if (mediumBlocks[i][1] != 0)
      fprintf (stderr, "%8dKb %8d\n", i, mediumBlocks[i][1]);
  for (i = 0; i < MAXLARGEBLOCK; i++)
    if (largeBlocks[i][1] != 0)
      fprintf (stderr, "%8dMb %8d\n", i, largeBlocks[i][1]);
  fprintf (stderr,
	   "\nMemory blocks moved when reallocated:\n    Size      Count\n");
  for (i = 0; i < MAXSMALLBLOCK; i++)
    if (smallBlocks[i][2] != 0)
      fprintf (stderr, "%8db  %8d\n", i, smallBlocks[i][2]);
  for (i = 0; i < MAXMEDIUMBLOCK; i++)
    if (mediumBlocks[i][2] != 0)
      fprintf (stderr, "%8dKb %8d\n", i, mediumBlocks[i][2]);
  for (i = 0; i < MAXLARGEBLOCK; i++)
    if (largeBlocks[i][2] != 0)
      fprintf (stderr, "%8dMb %8d\n", i, largeBlocks[i][2]);
  fflush (stderr);
  return;
}
#endif

void *
swCalloc (size_t size, int count, char *fileName, int lineNo)
{
  void *newBlock;
  size_t totalSize;

  totalSize = count * size;
  newBlock = swMalloc (totalSize, fileName, lineNo);
  memset (newBlock, 0, totalSize);
  return (newBlock);
}

void *
swMalloc (size_t size, char *fileName, int lineNo)
{
  void *newBlock;

#ifdef DMTRACK
  char messageBuffer[MAXSWMSG];
#endif

  if (firstMallocCall) {
    swFirstDM ();
  }

#ifdef DMTRACK
  countMalloc += 1;
  totalMalloc += size;
  currentAlloc += size;
  if (currentAlloc > peakAlloc)
    peakAlloc = currentAlloc;
  if (size < MAXSMALLBLOCK)
    smallBlocks[size][0]++;
  else {
    if (size < MAXMEDIUMBLOCK * MAXSMALLBLOCK)
      mediumBlocks[size / MAXSMALLBLOCK][0]++;
    else {
      if (size < MAXLARGEBLOCK * MAXMEDIUMBLOCK * MAXSMALLBLOCK)
	largeBlocks[size / (MAXMEDIUMBLOCK * MAXSMALLBLOCK)][0]++;
      else
	WARNING("Block of size %lu exceeds %dMb, not tracked", size, MAXLARGEBLOCK);
    }
  }
#endif
  newBlock = malloc (size);
#ifdef DMTRACK
  swAddChunk (newBlock, size, cTMalloc, fileName, lineNo);
#endif
  return (newBlock);
}

void *
swRealloc (void *pBlock, size_t newSize, char *fileName, int lineNo)
{
  void *newBlock;

#ifdef DMTRACK
  char messageBuffer[MAXSWMSG];
  int reallocFlag, oldSize;
#endif

  if (firstMallocCall) {
    swFirstDM ();
  }

#ifdef DMTRACK
  oldSize = swDelChunk (pBlock, cTReFree, fileName, lineNo);
  countReallocFree += 1;
  totalReallocFree += oldSize;
  currentAlloc -= oldSize;
#endif
  newBlock = realloc (pBlock, newSize);
#ifdef DMTRACK
  if (newBlock == pBlock) {
    reallocFlag = 1;
    countReallocOK += 1;
    totalReallocOK += newSize;
  } else {
    countReallocMove += 1;
    reallocFlag = 2;
    totalReallocMove += newSize;
  }
  currentAlloc += newSize;
  if (currentAlloc > peakAlloc)
    peakAlloc = currentAlloc;
  if (newSize < MAXSMALLBLOCK)
    smallBlocks[newSize][0]++;
  else {
    if (newSize < MAXMEDIUMBLOCK * MAXSMALLBLOCK)
      mediumBlocks[newSize / MAXSMALLBLOCK][reallocFlag]++;
    else {
      if (newSize < MAXLARGEBLOCK * MAXMEDIUMBLOCK * MAXSMALLBLOCK)
	largeBlocks[newSize /
		    (MAXMEDIUMBLOCK * MAXSMALLBLOCK)][reallocFlag]++;
      else {
	sprintf (messageBuffer,
		 "Block of size %lu exceeds %dMb, not tracked", newSize,
		 MAXLARGEBLOCK);
	INFO(messageBuffer);
      }
    }
  }
  swAddChunk (newBlock, newSize, cTReAlloc, fileName, lineNo);
#endif
  return (newBlock);
}

void
swFree (void *pBlock, char *fileName, int lineNo)
{
#ifdef DMTRACK
  size_t oldSize;
#endif

#ifdef DMTRACK
  oldSize = swDelChunk (pBlock, cTFree, fileName, lineNo);
  countFree += 1;
  totalFree += oldSize;
  currentAlloc -= oldSize;
#endif
  free (pBlock);
  return;
}

#ifdef TELLRITA
int
udpSend (char *hostName, int serverPort, char *message)
{
  int sockfd;
  struct sockaddr_in their_addr;	// connector's address information
  struct hostent *he;
  int numbytes;
  char *envVar;

  if ((he = gethostbyname (hostName)) == NULL) { /* Get the host info */
    if ((envVar = getenv ("swLogMsgHost")) != NULL) { /* Can't get out, try a local relay */
      if ((he = gethostbyname (envVar)) == NULL)
	return EXIT_FAILURE;
    } else {
      return EXIT_FAILURE;
    }
  }

  if ((sockfd = socket (AF_INET, SOCK_DGRAM, 0)) == -1) {
    perror ("socket");
    return EXIT_FAILURE;
  }

  their_addr.sin_family = AF_INET;	// host byte order
  their_addr.sin_port = htons (serverPort);	// short, network byte order
  their_addr.sin_addr = *((struct in_addr *) he->h_addr);
  memset (their_addr.sin_zero, '\0', sizeof their_addr.sin_zero);

  if ((numbytes = sendto (sockfd, message, strlen (message), 0,
			  (struct sockaddr *) &their_addr,
			  sizeof their_addr)) == -1) {
    perror ("sendto");
    return EXIT_FAILURE;
  }
  close (sockfd);

  return EXIT_SUCCESS;
}
#endif

#ifdef FULLLOG
FILE *fullLogFile = NULL;
char *fullLogFileName = "kelvin.full_log";
#endif

#define MAXDUPWARNING 4
char lastWarnings[MAXDUPWARNING][MAXLOGMSG + 1];
int lastLastWarning = 0;
void
swLogMsg (FILE *stream, char *message)
{
  time_t nowSec;
  struct tm *nowTm;
  int i;

  if (strncmp(message, "WARNING", 7) == 0) {
    for (i=0; i<MAXDUPWARNING; i++)
      if (strcmp(message, lastWarnings[i]) == 0)
	return;
    strcpy(lastWarnings[++lastLastWarning % MAXDUPWARNING], message);
  }

  nowSec = time (NULL);
  nowTm = localtime(&nowSec);
  fprintf (stream, "%02d/%02d/%02d %02d:%02d:%02d %s\n", nowTm->tm_year - 100, nowTm->tm_mon + 1, nowTm->tm_mday,
	   nowTm->tm_hour, nowTm->tm_min, nowTm->tm_sec, message);
  fflush (stream);
#ifdef TELLRITA
  char udpBuffer[MAXUDPMSG + 1];

  snprintf (udpBuffer, MAXUDPMSG, "PID: %d %s\n", (int) getpid (), message);
  // Don't need timestamp here because it'll show up with its own
  if (udpSend ("levi-montalcini.ccri.net", 4950, udpBuffer) == EXIT_SUCCESS)
    return;
#endif
#ifdef FULLLOG
  if (fullLogFile == 0)
    fullLogFile = fopen (fullLogFileName, "w+");
  fprintf (fullLogFile, "%02d/%02d/%02d %02d:%02d:%02d %s\n", nowTm->tm_year - 100, nowTm->tm_mon + 1, nowTm->tm_mday,
	   nowTm->tm_hour, nowTm->tm_min, nowTm->tm_sec, message);
  fflush (stream);
#endif
  return;
}

#define MAXPROGRESSLEVELS 8

struct progressMessage {
  time_t eventTime;
  int percentDone;
  char seen;
  char text[MAXLOGMSG + 1];
};

struct progressMessage progressLevels[MAXPROGRESSLEVELS + 1];
time_t logStartTime = 0;
int swProgressLevel = MAXPROGRESSLEVELS;
int logMaxUsedLevels = 0;

void
swLogTimedProgress() {
  int i, nothingShown = TRUE;

  for (i=0; i< MIN(logMaxUsedLevels + 1, swProgressLevel); i++)
    if (progressLevels[i].seen == FALSE) {
      swLogMsg (stderr, progressLevels[i].text);
      nothingShown = FALSE;
      progressLevels[i].seen = TRUE;
    }
  if (nothingShown)
    swLogMsg (stderr, "\t\t\t...running...");
}

char *formatElapsedTime (unsigned int t, char *buffer)
{
  unsigned int u /* Seconds less days */, v /* Seconds less days and hours */, w /* Seconds less minutes */;
  unsigned int d, h, m, s;
  // s = v - (w = (m = (v = (u - ((h = (u = (t - ((d = t/(24*60*60))*24*60*60))) / (60*60)) *60*60))) / 60) *60);
  d = t/(24*60*60);
  u = t - (d*24*60*60);
  h = u / (60*60);
  v = u - (h *60*60);
  m = v / 60;
  w = m *60;
  s = v - w;
  if (d != 0)
    sprintf (buffer, "%dd%dh%dm%ds", d, h, m, s);
  else if (h != 0)
    sprintf (buffer, "%dh%dm%ds", h, m, s);
    else if (m != 0)
      sprintf (buffer, "%dm%ds", m, s);
      else
	sprintf (buffer, "%ds", s);
  return buffer;
}

void 
swLogProgress(int level, float percentDone, char *format, ...) {
  int length;
  char tabs[MAXPROGRESSLEVELS];
  char *pMessage = progressLevels[level].text;
  char timeBuffer[32];
  va_list argp;

  /* This may look wrong, but we really do want to suppress any
     old lower-level messages if our new one is a higher level,
     because they're sub-steps or details for a previous step. */
  logMaxUsedLevels = level;

  if (logStartTime == 0)
    logStartTime = time (NULL);

  progressLevels[level].eventTime = time (NULL) - logStartTime;
  progressLevels[level].percentDone = percentDone;
  progressLevels[level].seen = FALSE;

  if (level >= MAXPROGRESSLEVELS)
    level = MAXPROGRESSLEVELS - 1;
  memset (tabs, '\t', MAXPROGRESSLEVELS);
  tabs[level] = '\0';
  if (progressLevels[level].percentDone > 0)
    pMessage += length = snprintf (progressLevels[level].text, MAXLOGMSG, "%s@%s (~%2d%%), ",
				   tabs, formatElapsedTime ((int) progressLevels[level].eventTime, timeBuffer),
				   progressLevels[level].percentDone);
  else
    pMessage += length = snprintf (progressLevels[level].text, MAXLOGMSG, "%s@%s, ",
				   tabs, formatElapsedTime ((int) progressLevels[level].eventTime, timeBuffer));
  
  va_start (argp, format);
  vsnprintf (pMessage, MAXLOGMSG - length, format, argp);
  va_end (argp);

  /* If we're not dumping everything the moment it happens, then
     skip any output */
  if (level > 0)
    if (swProgressDelaySeconds > 0 || level > swProgressLevel)
      return;

  swLogMsg (stderr, progressLevels[level].text);
  fflush (stderr);
  progressLevels[level].seen = TRUE;

}

/// Global flag to indicate that a progress advisory has been requested.
volatile sig_atomic_t swProgressRequestFlag = FALSE;

/// Global number of minutes to delay between progress updates.
int swProgressDelaySeconds;

/**
  Code for the progress thread.  We start by setting the request flag, 
  then we wait for the interval. If the flag is still set at the end of 
  the interval, we need to provide whatever progress information we can 
  ourselves. Otherwise, something in the program has given an update, and 
  we can go back to sleep for a while.
*/
void
*progressSignalHandler (void *threadId)
{
  int sleepSeconds;

  swProgressRequestFlag = TRUE;
  sleepSeconds = swProgressDelaySeconds;
  while (sleepSeconds > 0) {
    sleep (sleepSeconds - 1);
    swProgressRequestFlag = TRUE;
    sleep (1); // Give any main program loop a moment to provide status
    swProgressRequestFlag = FALSE;
    swLogTimedProgress ();
  }
  pthread_exit (NULL);
  return (void *) NULL; // Keep the cygwin GCC compiler happy
}

/// Timer thread to advise of progress
pthread_t progressWakeUpThread = 0;

/* Start a thread to wake-up and set a progress display flag at some interval. */
void
swStartProgressWakeUps (int seconds)
{
  swProgressDelaySeconds = seconds;

  INFO ("Further progress will be displayed at %d second intervals", seconds);
  if (progressWakeUpThread == 0)
    if (pthread_create (&progressWakeUpThread, NULL, progressSignalHandler, NULL))
      WARNING("Failed to create progress request thread, no interval-based progress advisories will be shown");
}

int procSMapWorks = TRUE;
long swGetCurrentVMK(pid_t pid) {
  char procSMap[PATH_MAX], buffer[128];
  FILE *fpProcSMap;
  long segmentSize = 0;

  currentVMK = 0;
  if (!procSMapWorks) {
    return currentVMK;
  }

  sprintf (procSMap, "/proc/%d/smaps", (int) pid);
  if ((fpProcSMap = fopen(procSMap, "r")) == NULL) {
    procSMapWorks = FALSE;
    return currentVMK;
  }  
  while (fgets(buffer, sizeof(buffer), fpProcSMap) != NULL)
    if (sscanf(buffer, "Size: %ld kB", &segmentSize) != 0)
      currentVMK += segmentSize;
  fclose(fpProcSMap);
  return currentVMK;
}

long swGetMaximumPMK(void) {
  char buffer[128];
  FILE *fpProcMemInfo;

  maximumPMK = 0;
  if ((fpProcMemInfo = fopen("/proc/meminfo", "r")) != NULL)
    while (fgets(buffer, sizeof(buffer), fpProcMemInfo) != NULL)
      if (sscanf(buffer, "MemTotal: %ld kB", &maximumPMK) != 0) {
	fclose(fpProcMemInfo);
	break;
      }
  return maximumPMK; 
}

volatile sig_atomic_t *envDiagLevel;
#if defined(DISTRIBUTION) || defined(__CYGWIN__)
#else
int swSharedDiagMemoryID;
#endif

/*

If you get a "bad system call" from cygwin using this, you need to do two things:

1. Make sure the environment variable CYGWIN is set to "server", and
2. Run (in another window) /usr/sbin/cygserver

*/

void swDiagInit(void) {
  int i;
  swProgressDelaySeconds = 120; ///< Default of two minutes delay between progress notifications
#if defined(DISTRIBUTION) || defined(__CYGWIN__) || !defined(USESHM)
  envDiagLevel = (int *) malloc(sizeof (int) * MAX_DIAG_FACILITY);
#else
  if ((swSharedDiagMemoryID = shmget ((key_t) getpid (), sizeof (int) * MAX_DIAG_FACILITY, 0666|IPC_CREAT)) == -1)
    ERROR ("Cannot create shared memory segment for diagnostics, %s", strerror(errno));

  INFO ("Use segment ID %d for diagnostic purposes", swSharedDiagMemoryID);
  if ((envDiagLevel = shmat (swSharedDiagMemoryID, NULL, 0)) == ((void *) -1))
    ERROR ("Cannot attach shared memory segment for diagnostics, %s", strerror(errno));
#endif
  // Not needed if we're setting them all to zero, since shmget does that.
  for (i=0; i<MAX_DIAG_FACILITY; i++)
    envDiagLevel[i] = 0;

  /* For the moment, to use, add a line like this right here: envDiagLevel[FACILITY] = N; where FACILITY is one of:

  OVERALL,
  LIKELIHOOD,
  READ_PEDFILE,
  ALLELE_SET_RECODING,
  GENOTYPE_ELIMINATION,
  PARENTAL_PAIR,
  CONFIG,
  INPUTFILE,
  XM,
  DCUHRE,
  POLYNOMIAL,
  ALTLSERVER,

  and N is the diagnostic level number, the higher the number the more diagnostics.

  e.g.:

  envDiagLevel[LIKELIHOOD] = 2;

  */

  return;
}

void swDiagTerm(void) {
#if defined(DISTRIBUTION) || defined(__CYGWIN__) || !defined(USESHM)
#else
  if (shmdt ((void *) envDiagLevel) != 0)
    ERROR ("Cannot detach shared memory segment, use ipcrm to clean-up");
  if (shmctl(swSharedDiagMemoryID, IPC_RMID, (struct shmid_ds *) NULL) != 0)
    ERROR ("Cannot mark shared memory segment for deletion, use ipcrm to clean-up");
#endif
  return;
}

/* Arrange the N elements of ARRAY in random order.  Only effective if N is much smaller than RAND_MAX;
   If this may not be the case, use a better random number generator. */
void swShuffle(int *array, size_t n)
{
  int stime;
  long ltime;

  ltime = time(NULL);
  stime = (unsigned) ltime/2;
  srand(stime);

  if (n > 1) {
    size_t i;
    for (i = 0; i < n - 1; i++) {
      size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
      int t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
}

#ifdef MAIN

#include <math.h>

struct swStopwatch *overallSW, *primeSW, *sqrtSW, *iterateSW, *fibSW, *restSW,
  *sleepSW, *napSW;

/* This is an example of recursive code that needs to be re-written in
   order to eliminate multiple exit points. */
int
fibonacci (int i)
{
  /* Code was...
     if (i <= 1)
     return 0;
     if (i == 2)
     return 1;
     return (fibonacci (i - 1) + fibonacci (i - 2));
     ...but now to ensure single entry/exit: */
  int j;

  /* See if we should dump our overall timer. We are considerate
     of both the executing code and the request for information. */
  if (swProgressRequestFlag) {
    DETAIL(0, "Down to %d in fibonacci recursion", i);
    swProgressRequestFlag = FALSE;
  }

  swStart (fibSW);
  if (i <= 1)
    j = 0;
  else if (i == 2)
    j = 1;
  else
    j = fibonacci (i - 1) + fibonacci (i - 2);
  swStop (fibSW);
  return j;
}

/* Test driver illustrates use of stopwatch and how to use user signals
   both in a timer and from another process to trigger a dump of statistics
   from within a cooperative chunk of code. It is possible to force a dump
   of statistics whether the code cooperates or not, but such an intrusion
   might well crash your process. */
int
main (int argc, char *argv[])
{
  int i, j, k, l, m, n;
  char messageBuffer[MAXSWMSG];

  swStartProgressWakeUps (10);

  STEP (0, "Initializing");

#ifdef DMTRACK
#warning "Dynamic memory usage dumping is turned on, so performance will be poor!"
  WARNING("Dynamic memory usage dumping is turned on, so performance will be poor!");
#endif

  /* Create the stopwatches we'll be using for this example... */

  SUBSTEP(0, "Creating stopwatches");
  overallSW = swCreate ("overall");
  primeSW = swCreate ("prime");
  iterateSW = swCreate ("iterations");
  sqrtSW = swCreate ("sqrt");
  fibSW = swCreate ("fibonacci");
  restSW = swCreate ("rest");
  sleepSW = swCreate ("sleep");
  napSW = swCreate ("nap");

  /* Start the overall stopwatch so we have the resource utilization info for
     the entire run of the program. */

  SUBSTEP(0, "Starting overall timer");
  swStart (overallSW);

  /* Now let's loop over our "hotspot" code (finding primes).
     primeSW will be the sum of sqrtSW + iterateSW + loop cost. */

  STEP(0, "Performing tests");
  SUBSTEP(0, "Testing progress advisory wake-ups");
  DETAIL(0, "Sleeping for 25 seconds");
  sleep (25);
  DETAIL(0, "Sleeping for 2 seconds");
  sleep (2);
  DETAIL(0, "Sleeping for 11 seconds");
  sleep (11);

  SUBSTEP(0, "Computing primes");
  INFO ("Computing primes between 1000000000 and 1000300000 responding to progress flag...");
  swStart (primeSW);
  k = 0;
  for (i = 1000000000; i < 1000300000; i++) {
    swStart (sqrtSW);
    k = sqrt (i) + 1;
    swStop (sqrtSW);

    /* See if we should dump our overall timer. We are considerate
       of both the executing code and the request for information by
       only modifying a volatile flag in the service routine, and then
       testing for that at our leisure in the executing code. */
    if (swProgressRequestFlag) {
      DETAIL(((i - 1000000000) * 100 / 300000), "Computing next prime after %d", i);
      swProgressRequestFlag = FALSE;
    }

    swStart (iterateSW);
    for (j = 2; j < k; j++)
      if (!(i % j))
	break;
    if (j == k) {
      k++;
      l = i;
    }
    swStop (iterateSW);
  }
  SUBSTEP(0, "Primes computed");
  INFO ("There are %d of them, and the last is %d", k, l);
  swStop (primeSW);

  SUBSTEP(0, "Computing fibonacci numbers");
  /* Now invoke a separate function like fibonacci. */
  INFO ("Computing fibonacci numbers up to the 32nd responding to progress flag...");
  for (i = 1; i <= 32; i++) {
    j = fibonacci (i);
    printf ("%d: %d\n", i, j);
  }

  /* Finally show some separate accumulation of elapsed time. restSW
     will be the sum of sleepSW and napSW. */

  SUBSTEP(0, "Testing timer accumulation and memory tracking");
  INFO ("Sleeping 28 seconds while marking memory in 200M and 1K chunks, ignoring progress flag...");
  swStart (restSW);
  char *p[7], *q[7];

  for (i = 0; i < 7; i++) {
    swStart (sleepSW);
    sleep (3);
    p[i] = (char *) swMalloc (256 * 1024 * 1024, __FILE__, __LINE__);
    memset (p[i], 'z', 256 * 1024 * 1024);
    swStop (sleepSW);
    swStart (napSW);
    sleep (1);
    q[i] = (char *) malloc (4 * 1024);
    memset (q[i], 'Z', 4 * 1024);
    swStop (napSW);
  }
  for (i = 0; i < 7; i++) {
    swFree (p[i], __FILE__, __LINE__);
    swFree (q[i], __FILE__, __LINE__);
  }
  swStop (restSW);

  p[0] = (char *) swMalloc (4 * 1024 * 1024, __FILE__, __LINE__);
  for (i = 3; i >= 0; i--) {
    INFO("realloc down %d", i);
    p[0] = swRealloc (p[0], i * 1024 * 1024, __FILE__, __LINE__);
  }

  STEP(0, "Displaying results");
  /* All done, so stop overall stopwatch and display results. */

  swStop (overallSW);

  INFO ("prime will be the sum of sqrt + iterate + unknown loop cost");
  swDumpM (sqrtSW);
  swDumpM (iterateSW);
  swDumpM (primeSW);

  swDumpM (fibSW);

  INFO ("rest will be almost exactly the sum of sleep and nap");
  swDumpM (sleepSW);
  swDumpM (napSW);
  swDumpM (restSW);

  INFO ("overall will be the sum of everything + unaccounted-for code");
  swDumpM (overallSW);

#ifdef DMTRACK
  swDumpSources ();
  swDumpCrossModuleChunks ();
  swDumpBlockUse ();
  swLogPeaks ("End of run");
#endif

  INFO("finished run");
  return 0;
}

#endif
