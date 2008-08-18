
/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

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
#if defined (DMTRACK) || defined (DMUSE)
#define malloc(X) swMalloc((X), __FILE__, __LINE__)
#define calloc(X,Y) swCalloc((X),(Y), __FILE__, __LINE__)
#define realloc(X,Y) swRealloc((X),(Y), __FILE__, __LINE__)
#define free(X) swFree((X), __FILE__, __LINE__)
#endif

DO NOT USE THESE IN SW.C because you'll get recurive looping.

The test main() serves as an example of use. To build a standalone test 
version:

$ gcc -lm -o swTest sw.c -DMAIN
$ ./swTest

To build an object you can link into your program:

$ gcc -c sw.c

then #include sw.h in your source code, and link with sw.o.

*/

#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#ifdef TELLRITA
#include <netdb.h>
#include <netinet/in.h>
#endif
#include "sw.h"
#include "hashtab.h"
#include "lookupa.h"

int currentVMK, maximumVMK = -1;

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
    exit (EXIT_FAILURE);
  } else {
    memset (newStopwatch, 0, sizeof (struct swStopwatch));
    strncpy (newStopwatch->swName, swName, MAXSWNAME);
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
    exit (EXIT_FAILURE);
  }
  if (getrusage (RUSAGE_CHILDREN, &theStopwatch->swStartRUChildren) != 0) {
    fprintf (stderr,
	     "Unable to get resource information for children using rusage call (%s)\n",
	     theStopwatch->swName);
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
    exit (EXIT_FAILURE);
  }
  if (getrusage (RUSAGE_CHILDREN, &stopStatsChildren) != 0) {
    fprintf (stderr,
	     "Unable to get resource information for children using rusage call (%s)\n",
	     theStopwatch->swName);
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
     theStopwatch->swAccumWallTime,
     theStopwatch->swAccumRUSelf.ru_utime.tv_sec + theStopwatch->swAccumRUChildren.ru_utime.tv_sec,
     theStopwatch->swAccumRUSelf.ru_stime.tv_sec + theStopwatch->swAccumRUChildren.ru_stime.tv_sec,
     theStopwatch->swAccumRUSelf.ru_nvcsw + theStopwatch->swAccumRUChildren.ru_nvcsw,
     theStopwatch->swAccumRUSelf.ru_nivcsw + theStopwatch->swAccumRUChildren.ru_nivcsw,
     theStopwatch->swAccumRUSelf.ru_minflt + theStopwatch->swAccumRUChildren.ru_minflt,
     theStopwatch->swAccumRUSelf.ru_majflt + theStopwatch->swAccumRUChildren.ru_majflt,
     appendText);
  swLogMsg (buffer);
  return;
}

void
swDumpM (struct swStopwatch *theStopwatch)
{
  char memoryBuffer[MAXUDPMSG];

  if (maximumVMK == -1) {
    maximumVMK = swGetMaximumVMK ();
    memoryBuffer[0] = 0;
  }
  if (maximumVMK != 0) {
    currentVMK = swGetCurrentVMK (getpid ());
    sprintf (memoryBuffer, ", vmG:%2.1f/%2.1f",
	     currentVMK / (1024.0 * 1024.0),
	     maximumVMK / (1024.0 * 1024.0));
  }

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

  strncpy (swName, theStopwatch->swName, MAXSWNAME);
  memset (theStopwatch, 0, sizeof (struct swStopwatch));
  strncpy (theStopwatch->swName, swName, MAXSWNAME);
  return;
}

/* Normally I wouldn't #ifdef-out anything but the essential calls, but
   these things (static allocation and memory tracking take up a lot of
   memory for their operation, and that's a waste if they're not in use. */

/* Static allocation alternative. */

/* Here's 84Mb of statically allocated 24-byte chunks (actually 28-bytes). */
#ifdef DMUSE
#define MAXFREE24S 7500000
struct free24
{
  unsigned char slot[24];
  struct free24 *next;
} free24s[MAXFREE24S];
struct free24 *nextFree24 = &free24s[0];
int used24s = 0, missed24s = 0;

/* Here's 156Mb of statically allocated 48-byte chunks (actually 52-bytes) */
#define MAXFREE48S 3000000
struct free48
{
  unsigned char slot[48];
  struct free48 *next;
} free48s[MAXFREE48S];
struct free48 *nextFree48 = &free48s[0];
int used48s = 0, missed48s = 0;

/* Here's 252Mb of statically allocated 80-byte chunks (actually 84-bytes) */
#define MAXFREE100S 9000000
struct free100
{
  unsigned char slot[100];
  struct free100 *next;
} free100s[MAXFREE100S];
struct free100 *nextFree100 = &free100s[0];
int used100s = 0, missed100s = 0;
#endif

/* Dynamic memory debugging. */

struct swStopwatch *internalDMSW;

/* The proper way to do this would be another hash, but there are very
   few calls and even fewer modules, so a simple ordered list and a
   binary search should suffice. */
enum callTypes
  { cTMalloc, ctCalloc, cTReAlloc, cTReFree, cTFree };
char *callTypeNames[] = { "malloc", "calloc", "realloc", "realloc", "free" };

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
} memChunkSources[MAXMEMCHUNKSOURCECOUNT];
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
swFirstDM ()
{
  internalDMSW = swCreate ("internalDMSW");
  swStart (internalDMSW);
  firstMallocCall = 0;

#ifdef DMTRACK
  chunkHash = hcreate (DMHASHSIZE);
#endif

  /* Experiment in memory self-management */
#ifdef DMUSE
  int i;

  for (i = 0; i < MAXFREE24S - 1; i++)
    free24s[i].next = &free24s[i + 1];
  free24s[MAXFREE24S - 1].next = NULL;
  for (i = 0; i < MAXFREE48S - 1; i++)
    free48s[i].next = &free48s[i + 1];
  free48s[MAXFREE48S - 1].next = NULL;
  for (i = 0; i < MAXFREE100S - 1; i++)
    free100s[i].next = &free100s[i + 1];
  free100s[MAXFREE100S - 1].next = NULL;
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
  newChunk->hashKey = lookup (&chunkAddress, sizeof (chunkAddress), 0);
  if (hadd
      (chunkHash, &newChunk->hashKey, sizeof (newChunk->hashKey),
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

  hashKey = lookup (&chunkAddress, sizeof (chunkAddress), 0);
  if (hfind (chunkHash, &hashKey, sizeof (hashKey)) == FALSE) {
    //    fprintf (stderr,
    //	     "WARNING!! (%s: %d) - free() of address %lu that was never allocated (head)!\n",
    //	     fileName, lineNo, (unsigned long) chunkAddress);
    return 0;
  } else {
    oldChunk = (struct memChunk *) hstuff (chunkHash);
    while (oldChunk->chunkAddress != chunkAddress) {
      if (oldChunk->next != NULL) {
	oldChunk = oldChunk->next;
      } else {
	//	fprintf (stderr,
	//		 "WARNING!! (%s: %d) = free() of address %lu that was never allocated (not head)!\n",
	//		 fileName, lineNo, (unsigned long) chunkAddress);
	return 0;
      }
    }
    if ((oldChunk->recycleCount % 2) == 0) {
      //      fprintf (stderr,
      //	       "ERROR!! (%s: %d) - free() of address %lu that is no longer in use for size %u and recycled %d times!\n",
      //	       fileName, lineNo, (unsigned long) chunkAddress, (unsigned int) oldChunk->chunkSize,
      //	       oldChunk->recycleCount);
      oldChunk->recycleCount++;
    } else {
      int i;
      for (i=0; i<memChunkSourceCount; i++)
        if (memChunkSources[i].entryNo == oldChunk->allocSource) {
//	  if (callType == cTReAlloc || callType == cTReFree)
	    //	    fprintf(stderr, "For %s: %d, subtracting %d from %d\n",
	    //		    fileName, lineNo, oldChunk->chunkSize, memChunkSources[i].remainingBytes);
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

  fprintf (stderr, "Top 20 sources out of %d total:\n",
	   memChunkSourceCount);
  qsort (memChunkSources, memChunkSourceCount, sizeof (struct memChunkSource),
	 compareSourcesByTotalBytes);
  for (i=0; i<min(memChunkSourceCount, 20); i++) {
    fprintf (stderr, "At %s line %d, %s called %d times for %ld bytes\n",
	     memChunkSources[i].moduleName, memChunkSources[i].lineNo,
	     callTypeNames[memChunkSources[i].callType],
	     memChunkSources[i].totalCalls, memChunkSources[i].totalBytes);
  }
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
  qsort (memChunkSources, ++memChunkSourceCount,
         sizeof (struct memChunkSource), compareSourcesByName);
}

void
swLogPeaks (char *reason)
{
  char messageBuffer[MAXSWMSG];
  if (firstMallocCall) {
    fprintf (stderr, "=> (%s): Before first memory allocation\n", reason);
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
  swLogMsg (messageBuffer);
  sprintf (messageBuffer,
	   "Size malloc:%g, free:%g, realloc OK:%g, realloc move:%g, realloc free:%g, current:%g, peak:%g",
	   totalMalloc, totalFree, totalReallocOK, totalReallocMove,
	   totalReallocFree, currentAlloc, peakAlloc);
  swLogMsg (messageBuffer);
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
  swLogMsg (messageBuffer);
  sprintf (messageBuffer,
	   "Size malloc: %g, free: %g, realloc OK: %g, realloc move: %g, realloc free: %g, current: %g, peak: %g",
	   totalMalloc, totalFree, totalReallocOK, totalReallocMove,
	   totalReallocFree, currentAlloc, peakAlloc);
  swLogMsg (messageBuffer);

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
  return;
}
#endif

void *
swCalloc (size_t count, size_t size, char *fileName, int lineNo)
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
#ifdef DMUSE
  struct free24 *your24;

  if (size <= 24) {
    used24s++;
    if (nextFree24 != NULL) {
      your24 = nextFree24;
      nextFree24 = your24->next;
      return (your24);
    } else {
      missed24s++;
      return (malloc (size));
    }
  }
  struct free48 *your48;

  if (size <= 48) {
    used48s++;
    if (nextFree48 != NULL) {
      your48 = nextFree48;
      nextFree48 = your48->next;
      return (your48);
    } else {
      missed48s++;
      return (malloc (size));
    }
  }
  struct free100 *your100;

  if (size <= 100) {
    used100s++;
    if (nextFree100 != NULL) {
      your100 = nextFree100;
      nextFree100 = your100->next;
      return (your100);
    } else {
      missed100s++;
      return malloc (size);
    }
  }
#endif

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
      else {
	sprintf (messageBuffer,
		 "Block of size %lu exceeds %dMb, not tracked", size,
		 MAXLARGEBLOCK);
	swLogMsg (messageBuffer);
      }
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
#ifdef DMUSE
  struct free24 *your24, *start24, *end24;

  start24 = &free24s[0];
  end24 = &free24s[MAXFREE24S - 1];
  if ((pBlock >= (void *) start24) && (pBlock <= (void *) end24)) {
    if (newSize <= 24)
      return (pBlock);
    your24 = (struct free24 *) pBlock;
    your24->next = nextFree24;
    nextFree24 = your24;
    return (swMalloc (newSize, fileName, lineNo));
  }
  struct free48 *your48, *start48, *end48;

  start48 = &free48s[0];
  end48 = &free48s[MAXFREE48S - 1];
  if ((pBlock >= (void *) start48) && (pBlock <= (void *) end48)) {
    if (newSize <= 48)
      return (pBlock);
    your48 = (struct free48 *) pBlock;
    your48->next = nextFree48;
    nextFree48 = your48;
    return (swMalloc (newSize, fileName, lineNo));
  }
  struct free100 *your100, *start100, *end100;

  start100 = &free100s[0];
  end100 = &free100s[MAXFREE100S - 1];
  if ((pBlock >= (void *) start100) && (pBlock <= (void *) end100)) {
    if (newSize <= 100)
      return (pBlock);

    your100 = (struct free100 *) pBlock;
    your100->next = nextFree100;
    nextFree100 = your100;
    return (swMalloc (newSize, fileName, lineNo));
  }
#endif

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
	swLogMsg (messageBuffer);
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

#ifdef DMUSE
  struct free24 *your24, *start24, *end24;

  start24 = &free24s[0];
  end24 = &free24s[MAXFREE24S - 1];
  if ((pBlock >= (void *) start24) && (pBlock <= (void *) end24)) {
    your24 = (struct free24 *) pBlock;
    your24->next = nextFree24;
    nextFree24 = your24;
    return;
  }
  struct free48 *your48, *start48, *end48;

  start48 = &free48s[0];
  end48 = &free48s[MAXFREE48S - 1];
  if ((pBlock >= (void *) start48) && (pBlock <= (void *) end48)) {
    your48 = (struct free48 *) pBlock;
    your48->next = nextFree48;
    nextFree48 = your48;
    return;
  }
  struct free100 *your100, *start100, *end100;

  start100 = &free100s[0];
  end100 = &free100s[MAXFREE100S - 1];
  if ((pBlock >= (void *) start100) && (pBlock <= (void *) end100)) {
    your100 = (struct free100 *) pBlock;
    your100->next = nextFree100;
    nextFree100 = your100;
    return;
  }
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

void
swLogMsg (char *message)
{
  char messageBuffer[MAXUDPMSG];

  sprintf (messageBuffer, "PID: %d, %s\n", (int) getpid (), message);
#ifdef TELLRITA
  if (udpSend ("levi-montalcini.ccri.net", 4950, messageBuffer) ==
      EXIT_FAILURE) messageBuffer[0] = 'p'; /* Yeah, it's embarassing */
#endif
  fprintf (stderr, messageBuffer);
  return;
}

int haveVMChildRunning = FALSE;
int swGetCurrentVMK(pid_t pid) {
  char commandString[128];
  FILE *gCFP;

  /* I've tried mallinfo and malloc_info (don't work), and mstats (DNE for Linux,
     On Crack for Mac). Help. */

  if ((maximumVMK == 0) || (haveVMChildRunning))
    return 0;

  haveVMChildRunning = TRUE;
  sprintf(commandString, "pmap %d 2>/dev/null | grep 'total' | cut -c 8-23", (int) pid);
  if ((gCFP = popen(commandString, "r")) != NULL) {
    fscanf(gCFP, "%d", &currentVMK);
    pclose(gCFP);
    haveVMChildRunning = FALSE;
  }
  return currentVMK;
}

int swGetMaximumVMK(void) {
  char commandString[] = "cat /proc/meminfo 2>/dev/null | grep 'MemTotal' | cut -c 11-22";
  /* For Mac  char commandString[] = "sysctl hw.memsize | cut -f 2 -d ' '";
     ...so now all we need is a reliable way to get current usage. */
  FILE *gCFP;
  
  if (maximumVMK != -1)
    return maximumVMK;

  if (haveVMChildRunning)
    return 0;

  haveVMChildRunning = TRUE;
  if ((gCFP = popen(commandString, "r")) != NULL) {
    fscanf(gCFP, "%d", &maximumVMK);
    pclose(gCFP);
    haveVMChildRunning = FALSE;
  }
  if (swGetCurrentVMK(getpid()) <= 0)
    maximumVMK = 0;
  return (maximumVMK);
}

#ifdef MAIN

#include <math.h>
#include <signal.h>

struct swStopwatch *overallSW, *primeSW, *sqrtSW, *iterateSW, *fibSW, *restSW,
  *sleepSW, *napSW;

/* Stuff for periodic timer signal. We use this to display statistics
   within cooperative code.  The same technique could be used to toggle
   debug output dynamically. */
volatile sig_atomic_t signalSeen = 0;
void
usr1SignalHandler (int signal)
{
  // signalSeen = 1;
  swLogPeaks ("Ding!");
}

void
quitSignalHandler (int signal)
{
  char messageBuffer[MAXSWMSG];

  swDump (overallSW);
#ifdef DMTRACK
  sprintf (messageBuffer,
	   "Count malloc: %d, free: %d, realloc OK: %d, realloc move: %d, realloc free: %d, max depth: %d, max recycles: %d",
	   countMalloc, countFree, countReallocOK, countReallocMove,
	   countReallocFree, maxListDepth, maxRecycles);
  swLogMsg (messageBuffer);
  sprintf (messageBuffer,
	   "Size malloc: %g, free: %g, realloc OK: %g, realloc move: %g, realloc free: %g, current: %g, peak: %g",
	   totalMalloc, totalFree, totalReallocOK, totalReallocMove,
	   totalReallocFree, currentAlloc, peakAlloc);
  swLogMsg (messageBuffer);
#endif
}

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
  if (signalSeen) {
    swDump (overallSW);
    signalSeen = 0;
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

  sprintf (messageBuffer, "At %s(%s %s)%d: %s", __FILE__, __DATE__, __TIME__,
	   __LINE__, "starting run");
  swLogMsg (messageBuffer);

#ifdef DMTRACK
#warning "Dynamic memory usage dumping is turned on, so performance will be poor!"
  swLogMsg
    ("Dynamic memory usage dumping is turned on, so performance will be poor!\n");
#endif

  /* Create the stopwatches we'll be using for this example... */

  overallSW = swCreate ("overall");
  primeSW = swCreate ("prime");
  iterateSW = swCreate ("iterations");
  sqrtSW = swCreate ("sqrt");
  fibSW = swCreate ("fibonacci");
  restSW = swCreate ("rest");
  sleepSW = swCreate ("sleep");
  napSW = swCreate ("nap");

  /* Fork a child that loops sleeping several seconds and then signalling 
     us with SIGUSR1. */

  pid_t childPID;

  childPID = fork ();
  if (childPID == 0) {
    while (1) {
      sleep (10);
      kill (getppid (), SIGUSR1);
    }				/* Does not return */
  }

  /* Setup signal handlers for SIGUSR1 and SIGQUIT (CTRL-\). */
  struct sigaction usr1Action, quitAction;
  sigset_t usr1BlockMask, quitBlockMask;

  sigfillset (&usr1BlockMask);
  usr1Action.sa_handler = usr1SignalHandler;
  usr1Action.sa_mask = usr1BlockMask;
  usr1Action.sa_flags = 0;
  sigaction (SIGUSR1, &usr1Action, NULL);
  sigfillset (&quitBlockMask);
  quitAction.sa_handler = quitSignalHandler;
  quitAction.sa_mask = quitBlockMask;
  quitAction.sa_flags = 0;
  sigaction (SIGQUIT, &quitAction, NULL);

  fprintf (stderr,
	   "To force a dump of stats, type CTRL-\\ (dangerous and terse but always works)\n");
  fprintf (stderr,
	   "or type \"kill -%d %d\" (safe and thorough, but requires program cooperation).\n",
	   SIGUSR1, getpid ());

  /* Start the overall stopwatch so we have the resource utilization info for
     the entire run of the program. */

  swStart (overallSW);

  /* Now let's loop over our "hotspot" code (finding primes).
     primeSW will be the sum of sqrtSW + iterateSW + loop cost. */

  printf
    ("Computing primes between 1000000000 and 1000200000 responding to SIGUSR1-set flag...\n");
  swStart (primeSW);
  k = 0;
  for (i = 1000000000; i < 1000200000; i++) {
    swStart (sqrtSW);
    k = sqrt (i) + 1;
    swStop (sqrtSW);

    /* See if we should dump our overall timer. We are considerate
       of both the executing code and the request for information by
       only modifying a volatile flag in the service routine, and then
       testing for that at our leisure in the executing code. */
    if (signalSeen) {
      swDump (overallSW);
      signalSeen = 0;
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
  printf ("There are %d of them, and the last is %d\n", k, l);
  swStop (primeSW);

  /* Now invoke a separate function like fibonacci. */
  printf
    ("Computing fibonacci numbers up to the 32nd responding to SIGUSR1-set flag...\n");
  for (i = 1; i <= 32; i++) {
    j = fibonacci (i);
    printf ("%d: %d\n", i, j);
  }

  /* Finally show some separate accumulation of elapsed time. restSW
     will be the sum of sleepSW and napSW. */

  printf
    ("Sleeping 21 seconds while marking memory in 2M and 1K chunks, ignoring SIGUSR1-set flag...\n");
  swStart (restSW);
  char *p[7], *q[7];

  for (i = 0; i < 7; i++) {
    swStart (sleepSW);
    sleep (2);
    p[i] = (char *) swMalloc (2 * 1024 * 1024, __FILE__, __LINE__);
    memset (p[i], 'z', 2 * 1024 * 1024);
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
    printf ("realloc down %d\n", i);
    p[0] = swRealloc (p[0], i * 1024 * 1024, __FILE__, __LINE__);
  }

  /* All done, so stop overall stopwatch and display results. */

  swStop (overallSW);

  printf ("prime will be the sum of sqrt + iterate + unknown loop cost\n");
  swDump (sqrtSW);
  swDump (iterateSW);
  swDump (primeSW);

  swDump (fibSW);

  printf ("rest will be almost exactly the sum of sleep and nap\n");
  swDump (sleepSW);
  swDump (napSW);
  swDump (restSW);

  printf ("overall will be the sum of everything + unaccounted-for code\n");
  swDump (overallSW);

#ifdef DMTRACK
  swDumpSources ();
  swDumpCrossModuleChunks ();
  swDumpBlockUse ();
  swLogPeaks ("End of run");
#endif

  swLogMsg ("finished run");
  return 0;
}

#endif
