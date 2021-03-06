/**
@file sSDHandler.c

 sSDHandler - manage storage on a Solid State Drive for polynomial term lists.

 We're out of memory, but have an SSD to help. Memory mapping a file
 on the SSD doesn't improve our situation because every page actually
 referenced has to be in physical memory. What we need is a "window"
 of fixed-memory size thru which we can access the full span of SSD
 space. That's normal file I/O, albeit with arbitrarily-sized chunks of
 data and an eye on efficient reuse. I searched and searched for a tool
 that would manage the space in the SSD efficiently and came up with
 nothing, so I wrote this.

 Philosophy:

 Allocate space in the SSD cache file in double pair chunks, i.e. some
 count of pairs of doubles.  Keep 16 lists of free chunks. List n
 contains only free chunks of size > 2^n double pairs.  Start with
 only one chunk in the 16th list that is the entire file. List heads
 are in memory, while subsequent list entries are at the starting
 offset for that entry (as pointed to by the previous entry's next
 pointer). Always push or pop to/from the head of the list. Don't keep
 back pointers, as it would double the number of writes, and would
 just make defragmentation easier, and defragmentation is not a
 priority. Find the leading high-order bit in the request size and
 left shift by one to determine the preferred list to use. The first
 chunk off of the preferred list will always work since they're all
 larger than the requested size. Use the entire chunk if it's a
 reasonably close fit, otherwise split it and return to the proper
 list. This works great for many initial allocations since we're just
 whittling-down the single large chunk in the 16th list, which is to
 say we're counting-down it's free space and starting address without
 ever writing any list information to the SSD. Returned entries and
 remainders smaller than their original list get pushed onto lesser
 lists. If the preferred list is empty, go up to the next larger list
 with entries. If all are empty, check the head of the list that is
 one smaller than the preferred one, as it may be large enough. If it
 isn't (but there are entries there), then pop all entries, sort them
 by descending size and push them back onto the list and try again. If
 it's still to small, do a garbage collection, and try again. If that
 fails, we're done.

 Garbage collection consists of popping all entries off of all lists
 and sorting them by their starting address. We then go thru this
 super list merging contiguous entries and pushing them back onto the
 appropriate lists. This eliminates "soap bubble" discontinuities, but
 not those where space is actually allocated between free entries.

 Usage:

 First call initSDD(), then call putSSD() to store a chunk of double pairs and get a ticket. Call
 getSSD() with the ticket to retrieve the chunk. Call freeSSD with the ticket to give-up the chunk.
 Call termSSD() when finished to display statistics. See the sample test driver at the end of
 the module.

 Performance:

 64K operations (putSSD, getSSD, freeSDD) of random sized chunks to to 64K in 47.5 seconds. 512K
 operations in 9 minutes.

 Caching aside (because I leave that to the OS), there is a fixed irreducible cost of 1 write 
 for a putSSD, 1 read for a getSSD. Anything beyond that is overhead. This algorithm's overhead is:

 - nothing for putSSD calls utilizing free lists where the remainder of the list head entry still
 belongs on the same list, e.g. 33 bytes taken from a 102 byte head entry on list 6 (2^6 = 64)
 leaving 69 bytes still in head entry of list 6.

 - nothing for putSSD calls using a single-entry list where the remainder of the list head entry
 is demoted but goes onto a completely empty list, e.g. 43 bytes taken from a 102 byte head entry
 on list 6 leaving 59 bytes which then moves to empty head entry for list 5 (2^5 = 32).

 - 1 write for putSSD calls using a single-entry list where the remainder of the list head entry
 is demoted and goes onto a non-empty list, e.g. 43 bytes taken from a 102 byte head entry
 on list 6 leaving 59 bytes which is then pushed onto list 5.

 - 1 read for putSSD calls using a multiple-entry list where the remainder of the list head entry
 is demoted but goes onto a completely empty list, e.g. 43 bytes taken from a 102 byte head entry
 on list 6 leaving 59 bytes which then moves to empty head entry for list 5 (2^5 = 32), and a
 new head entry is popped for list 6.

 - 1 read and 1 write for putSSD calls using a multiple-entry list where the remainder of the 
 list head entry is demoted and goes onto a non-empty list, e.g. 43 bytes taken from a 102 byte
 head entry on list 6 leaving 59 bytes which is then pushed onto list 5, and a new head entry
 is popped for list 6.

 - nothing for a freeSSD call where the returned entry goes onto a completely empty list.

 - 1 write for a freeSSD call where the returned entry gets pushed onto a non-empty list.

 Note:

 For diagnostic purposes, look in the cache file using something like:

 od --address-radix=d --read-bytes=<cO X 16> --format=f8 --output-duplicates --skip-bytes=<dPC X 16> /tmp/ssd/cache.dat 

 ...and compile test standalone version with: gcc -g -o sSDHandler -DMAIN sSDHandler.c

  @version $Id$

  @author Bill Valentine-Cooper.

  Copyright (C) 2008, 2010, 2022 Mathematical Medicine LLC
  This program is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.
  You should have received a copy of the GNU General Public License along
  with this program. If not, see <https://www.gnu.org/licenses/>.

*/

#include <limits.h>

#ifdef __APPLE__
#include <sys/param.h>
#include <sys/mount.h>
#else
#include <sys/vfs.h>
#include <sys/statfs.h>
#endif

#include <libgen.h>

#include "utils.h"
#include "sSDHandler.h"

// The default SSD path and size in Gb
char *defaultSSDFileName = "/ssd/cache.dat";
int defaultSSDFileSizeInGb = 28;
unsigned long maxSSDDPC;

#define DOUBLE_PAIR_SIZE (sizeof (double) * 2)

char messageBuffer[256];

FILE *sSDFD;

struct listEntry {
  unsigned long doublePairCount; /* Size in double pairs of this chunk, largest feasible is 
				    2^32-1*16 = ~64Gb, but might ought to be less than 2^16 (MAX_DPC_MASK) so
				    list 16 satisfies any request. 0 if end of list */
  unsigned long chunkOffset; // Offset in double pairs from beginning of file, max of 64Gb/16b = 4M.
  unsigned long nextFree; // Offset in double pairs of next listEntry in free list, 0 if DNE.
};
struct listEntry listHead[16];
int listDepth[16], getCallCount, putCallCount, freeCallCount;
unsigned long usedDPCs;
unsigned long handledDPCs;

/* 
   Returns index of leftmost set bit (right most being 0), max
   of 15.

   This is weighted to the high-end because although there
   are more short term lists, they are handled proportionately
   less often. This shifting approach is significantly faster
   than inequality comparison or masking when compiled without
   optimization, and only slightly slower when optimized.
*/
unsigned short high16Bit (unsigned long value) {
  unsigned short i;
  value >>= (MIN_USE_SSD_BITS - 1);
  if (value >= 0x8000)
    return 15;
  for (i=14; i>0; i--)
    if (value & 0x4000)
      break;
    else
      value <<= 1;
  return i;
}

/**

   Dump statistics on the use of the SSD to stderr.

*/
void statSSD () {
  int i;
  fprintf (stderr, "SSD list use: ");
  for (i=0; i<16; i++)
    fprintf (stderr, "%d:%d ", i, listDepth[i]);
  fprintf (stderr, "\nput/get/free calls: %d/%d/%d, total used: %lu, handled: %lu\n",
	   putCallCount, getCallCount, freeCallCount, usedDPCs, handledDPCs);
}

/**

   Finish using the SSD.

   Nothing fancy here. Just close the open file on the SSD.

  @par Global Inputs

  sSDFD - File descriptor for the open file on the SSD.

  @par Global Outputs

  @return void.

*/  
void termSSD () {
  fclose (sSDFD);
}

/**

   Prepare for using the SSD.

   Use environment variables to identify the name and size of the SSD file
   to open. Open the file and initialize free list pointers.

  @par Global Inputs

  sSDFileName - environment variable containing the full path and name of the 
  file to use on the solid state drive. If the environment variable is not
  defined, a default is used.

  sSDFileSizeInGb - environment variable containing a number which is the amount
  of space on the solid state drive to use. If the environment variable is not
  defined, the entire drive is used.

  @par Global Outputs

  The open file handle sSDFD.

  @return void.

*/
void initSSD() {
  int i;
  unsigned long sSDFileSizeInGb;
  char *envVar;
  char sSDFileName[PATH_MAX], *dirnameIsDamaged;
  struct statfs sSDStats;

  if ((envVar = getenv ("sSDFileName")) != NULL)
    strcpy (sSDFileName, envVar);
  else
    strcpy (sSDFileName, defaultSSDFileName);

  if ((envVar = getenv ("sSDFileSizeInGb")) != NULL)
    sSDFileSizeInGb = atol(envVar);
  else {
    dirnameIsDamaged = strdup (sSDFileName);
    if (
#ifdef __sun__
	statfs(dirname(dirnameIsDamaged), &sSDStats, 0, 0)
#else
	statfs(dirname(dirnameIsDamaged), &sSDStats)
#endif
	!= 0) {
      perror("statfs call failed on SSD path, using default");
      sSDFileSizeInGb = defaultSSDFileSizeInGb;
    } else
      sSDFileSizeInGb = (sSDStats.f_blocks * (sSDStats.f_bsize / 1024)) / (1024 * 1024) - 1;
  }

// Maximum number of double pairs on the 28Gb SSD is 28*1024*1024*1024/16, about 1792Mdps

  maxSSDDPC = sSDFileSizeInGb * 1024UL * 1024UL / 16UL * 1024UL;
  fprintf (stderr, "Using SSD located at %s of %luGb, or %ludps\n", 
	   sSDFileName, sSDFileSizeInGb, maxSSDDPC);

  if ((sSDFD = fopen(sSDFileName,"wb+")) == NULL) {
    perror ("Failed to open SSD cache file");
    exit (EXIT_FAILURE);
  }
  for (i=0; i<16; i++) {
    listHead[i].doublePairCount = listHead[i].chunkOffset = listDepth[i] = 0;
    listHead[i].nextFree = ULONG_MAX;
  }
  listHead[15].doublePairCount = maxSSDDPC;
  listDepth[15] = 1;
  getCallCount = putCallCount = freeCallCount = 0;
  return;
}

/**

   Simply replace the current list head with it's successor.

*/
void removeFreeListHead (unsigned short freeList) {

  // Handle failures first
  if (listHead[freeList].doublePairCount == 0) {
#ifdef DEBUG
    fprintf (stderr, "Attempt to remove head of empty L%d (Zero dpc, listDepth %d)!\n",
	     freeList, listDepth[freeList]);
#endif
    exit (EXIT_FAILURE);
  }
  if (listDepth[freeList] == 0) {
#ifdef DEBUG
    fprintf (stderr, "Attempt to remove head of empty L%d (Zero listDepth, %ludpc)!\n",
	     freeList, listHead[freeList].doublePairCount);
#endif
    exit (EXIT_FAILURE);
  }

  // Now handle last entries
  if (listHead[freeList].nextFree == ULONG_MAX) {
    listDepth[freeList] = 0;
    listHead[freeList].doublePairCount = 0;
#ifdef DEBUG
    fprintf (stderr, "Removed last head from L%d (according to nextFree)\n", freeList);
#endif
  } else {
    // Pull-in the next entry and shlorp-up the listEntry space!
    listDepth[freeList]--;
    if (fseek (sSDFD, listHead[freeList].nextFree * DOUBLE_PAIR_SIZE, SEEK_SET) != 0) {
      perror ("Failed to seek in SSD cache file");
      exit (EXIT_FAILURE);
    }
    if ((fread (&listHead[freeList], sizeof (struct listEntry), 1, sSDFD)) != 1) {
      perror ("Failed to read SSD cache file");
      exit (EXIT_FAILURE);
    }
    listHead[freeList].chunkOffset--;
    listHead[freeList].doublePairCount++;
#ifdef DEBUG
    fprintf (stderr, "Popped new head of L%d at co%lu of %ludps", freeList,
	     listHead[freeList].chunkOffset, listHead[freeList].doublePairCount);
    if (listHead[freeList].nextFree == ULONG_MAX)
      fprintf (stderr, ", next to end of list (one more, listDepth is %d!)\n", listDepth[freeList]);
    else
      fprintf (stderr, ", next at co%lu\n", listHead[freeList].nextFree);
#endif
  }
  return;
}

void insertFreeListHead (unsigned long chunkOffset, unsigned long doublePairCount) {
  unsigned short freeList;

  if (doublePairCount > 0)
    freeList = high16Bit (doublePairCount - 1UL);
  else
    freeList = 0;
#ifdef DEBUG
  fprintf (stderr, "Pushing co%lu of %ludps onto L%d\n",
	   chunkOffset, doublePairCount, freeList);
#endif
  if (listHead[freeList].doublePairCount == 0UL) {
    // First one, an easy insertion
#ifdef DEBUG
    fprintf (stderr, "First of list!\n");
#endif
    listHead[freeList].chunkOffset = chunkOffset;
    listHead[freeList].doublePairCount = doublePairCount;
    listHead[freeList].nextFree = ULONG_MAX;
    listDepth[freeList] = 1;
    return;
  }

  // Not the first one, but adjacent to it?
  if (chunkOffset+doublePairCount == listHead[freeList].chunkOffset) {
    // Add to beginning...
#ifdef DEBUG
    fprintf (stderr, "Immediately preceeds first of list!\n");
#endif
    listHead[freeList].chunkOffset = chunkOffset;
    listHead[freeList].doublePairCount += doublePairCount;
    return;
  }

  if (listHead[freeList].chunkOffset + listHead[freeList].doublePairCount == chunkOffset) {
    // Add to end...
#ifdef DEBUG
    fprintf (stderr, "Immediately follows first of list!\n");
#endif
    listHead[freeList].doublePairCount += doublePairCount;
    return;
  }

  // Not adjacent in any sense, just push it onto the free list...
  /* BTW - Combining contiguous free entries to defragment dynamically would require
     backpointers, which would double our read/write burden. */
  // Not the first one, need to flush the old first one to the SSD cache file
  
  listDepth[freeList]++;
  
  if (fseek (sSDFD, listHead[freeList].chunkOffset * DOUBLE_PAIR_SIZE, SEEK_SET) != 0) {
    sprintf (messageBuffer, "Failed to seek to %ludps in SSD cache file", listHead[freeList].chunkOffset);
    perror (messageBuffer);
    exit (EXIT_FAILURE);
  }
  listHead[freeList].chunkOffset++;
  listHead[freeList].doublePairCount--;
  if ((fwrite (&listHead[freeList], sizeof (struct listEntry), 1, sSDFD)) != 1) {
    perror ("Failed to write SSD cache file");
    exit (EXIT_FAILURE);
  }
  listHead[freeList].nextFree = listHead[freeList].chunkOffset - 1UL;
  listHead[freeList].chunkOffset = chunkOffset;
  listHead[freeList].doublePairCount = doublePairCount;

  return;
}

struct freeVectorEntry {
  unsigned long coStart;
  unsigned long coEnd;
};

static void swap_internal(char *a, char *b, size_t size)
{
  if (a != b)
    {
      char t;
      while (size--)
	{
	  t = *a;
	  *a++ = *b;
	  *b++ = t;
	}
    }
}

static void qsort_internal(char *begin, char *end, size_t size, int(*compar)(const void *, const void *))
{
  if (end > begin)
    {
      char *pivot = begin;
      char *l = begin + size, *r = end;

      while (l < r)
	{
	  if (compar(l, pivot) <= 0)
	    {
	      l += size;
	    }
	  else
	    {
	      r -= size;
	      swap_internal(l, r, size);
	    }
	}
      l -= size;
      swap_internal(begin, l, size);
      qsort_internal(begin, l, size, compar);
      qsort_internal(r, end, size, compar);
    }
}

void qqsort(void *base, size_t nmemb, size_t size, int(*compar)(const void *, const void *))
{
  qsort_internal((char *)base, (char *)base+nmemb*size, size, compar);
}

int compareFVE (const void *left, const void *right) {
  return (((struct freeVectorEntry *) left)->coStart - ((struct freeVectorEntry *) right)->coStart);
}

int compareSize (const void *left, const void *right) {
  return (
	  (
	   ((struct freeVectorEntry *) left)->coEnd -
	   ((struct freeVectorEntry *) left)->coStart
	  ) -
	  (
	   ((struct freeVectorEntry *) right)->coEnd -
	   ((struct freeVectorEntry *) right)->coStart
	  )
	 );
}

void garbageCollect () {
  /*
    This is not defragmentation, as we're not keeping track of chunkTickets we hand out, so we
    couldn't fix-up chunkOffsets. I think tracking usage will make you go blind.

    This is a complete rebuild of the freeLists. It is complicated, but can be done pretty quickly:

    1. Allocate a single vector of of offsets and endpoints of appropriate length (using 
    listDepths) for all freeLists.
    2. Follow each freeList thru the file collecting chunkOffsets and endpoints into the vector.
    3. Sort each vector by ascending chunkOffsets.
    4. Reset free listHead and listDepth.
    4. Collapse contiguous entries (next offset = previous endpoint)
    5. Inserting each fully-collapsed chunk into the appropriate freeList.

  */
  int i;
  unsigned long totalFreeBefore = 0, totalFreeAfter = 0, coStart = 0, coEnd = 0, freeVectorEntryCount;
  struct freeVectorEntry *freeVector;

  fprintf (stderr, "sSDHandler garbage collection started...\n");
  statSSD ();

  // Collect all freeList entries regardless of length

  freeVectorEntryCount = 0;
  for (i=0; i<16; i++)
    freeVectorEntryCount += listDepth[i];

  fprintf (stderr, "There are %lu total free list entries\n",  freeVectorEntryCount);
  if ((freeVector = (struct freeVectorEntry *) 
       malloc (freeVectorEntryCount * sizeof (struct freeVectorEntry))) == NULL) {
    fprintf (stderr, "Failed to allocate a %lu byte freeVector for garbage collection\n",
	     freeVectorEntryCount * sizeof (struct freeVectorEntry));
    exit (EXIT_FAILURE);
  }

  fprintf (stderr, "...extracting from list ");
  freeVectorEntryCount = 0;
  for (i=0; i<16; i++) {
    fprintf (stderr, "%d ", i);
    fflush (stderr);
    while (listHead[i].doublePairCount != 0) {
      freeVector[freeVectorEntryCount].coStart = listHead[i].chunkOffset;
      freeVector[freeVectorEntryCount].coEnd = listHead[i].chunkOffset + listHead[i].doublePairCount;
      freeVectorEntryCount++;
      totalFreeBefore += listHead[i].doublePairCount;
      removeFreeListHead (i);
    }
  }
  fprintf (stderr, "\nSorting %lu entries...", freeVectorEntryCount);
  fflush (stderr);

#ifdef DEBUG
  fprintf (stderr, "Pre-sort:\n");
  for (i=0; i<freeVectorEntryCount; i++)
    fprintf (stderr, "#%d: co%lu-co%lu\n", i, freeVector[i].coStart, freeVector[i].coEnd);
#endif

  qqsort (freeVector, freeVectorEntryCount, sizeof (struct freeVectorEntry), compareFVE);

#ifdef DEBUG
  fprintf (stderr, "Post-sort:\n");
  for (i=0; i<freeVectorEntryCount; i++) {
    fprintf (stderr, "#%d: co%lu-co%lu\n", i, freeVector[i].coStart, freeVector[i].coEnd);
    if ((i > 0) && (freeVector[i].coStart < freeVector[i-1].coEnd)) {
      fprintf (stderr, "Last two free chunks overlap!\n");
      exit (EXIT_FAILURE);
    }
  }
#endif

  fprintf (stderr, "\nRe-inserting...");

  // Reset listHead and listDepth

  for (i=0; i<16; i++) {
    listHead[i].doublePairCount = listHead[i].chunkOffset = listDepth[i] = 0;
    listHead[i].nextFree = ULONG_MAX;
  }

  // Connect and insert them.

  int scale = (freeVectorEntryCount > 100 ? 100 : 10);
  int j = MAX(1, (freeVectorEntryCount / scale));
  for (i=0; i<freeVectorEntryCount; i++) {
    if (i % j == 0) {
      fprintf (stderr, "\rRe-inserting...%lu%% done", i * 100 / freeVectorEntryCount);
      fflush (stderr);
    }
    coStart = freeVector[i].coStart;
    coEnd = freeVector[i].coEnd;
    while (i<(freeVectorEntryCount-1)) {
      if (coEnd == freeVector[i+1].coStart) {
	coEnd = freeVector[++i].coEnd;
      } else {
#ifdef DEBUG
	fprintf (stderr, "Unfree hole from co%lu to co%lu of %ludps\n",
		 coEnd, freeVector[i+1].coStart, (freeVector[i+1].coStart - coEnd));
#endif
	break;
      }      
    }
    totalFreeAfter += coEnd-coStart;
    insertFreeListHead (coStart, coEnd-coStart);
  }
  fprintf (stderr, "\rRe-inserting...100%% done\n");

  // Say how it went.
  fprintf (stderr, "Began with %ludps free and ended with %ludps\n",
	  totalFreeBefore, totalFreeAfter);
  statSSD ();
  free (freeVector);
}

/*
  Reorder by descending size one of the free lists so we can reliably snag largest
  entries when we're forced to drop below where we intend to operate.
*/

void reorderFreeList (int freeList) {
  unsigned long totalFreeBefore = 0, totalFreeAfter = 0, coStart = 0, coEnd = 0, freeVectorEntryCount;
  struct freeVectorEntry *freeVector;
  int i;

  fprintf (stderr, "sSDHandler reordering free list %d of %d entries...\n", freeList, listDepth[freeList]);

  freeVectorEntryCount = listDepth[freeList];

  if ((freeVector = (struct freeVectorEntry *) 
       malloc (freeVectorEntryCount * sizeof (struct freeVectorEntry))) == NULL) {
    fprintf (stderr, "Failed to allocate a %lu byte freeVector for garbage collection\n",
	     freeVectorEntryCount * sizeof (struct freeVectorEntry));
    exit (EXIT_FAILURE);
  }

  fprintf (stderr, "...extracting from list ");
  freeVectorEntryCount = 0;
  while (listHead[freeList].doublePairCount != 0) {
    freeVector[freeVectorEntryCount].coStart = listHead[freeList].chunkOffset;
    freeVector[freeVectorEntryCount].coEnd = listHead[freeList].chunkOffset + listHead[freeList].doublePairCount;
    freeVectorEntryCount++;
    totalFreeBefore += listHead[freeList].doublePairCount;
    removeFreeListHead (freeList);
  }
  fprintf (stderr, "\nSorting %lu entries...", freeVectorEntryCount);
  fflush (stderr);

  qqsort (freeVector, freeVectorEntryCount, sizeof (struct freeVectorEntry), compareSize);

  fprintf (stderr, "\nRe-inserting...");

  // Reset listHead and listDepth

  listHead[freeList].doublePairCount = listHead[freeList].chunkOffset = listDepth[freeList] = 0;
  listHead[freeList].nextFree = ULONG_MAX;

  // Connect and insert them.

  int scale = (freeVectorEntryCount > 100 ? 100 : 10);
  int j = MAX(1, (freeVectorEntryCount / scale));
  for (i=0; i<freeVectorEntryCount; i++) {
    if (i % j == 0) {
      fprintf (stderr, "\rRe-inserting...%lu%% done", i * 100 / freeVectorEntryCount);
      fflush (stderr);
    }
    coStart = freeVector[i].coStart;
    coEnd = freeVector[i].coEnd;
    while (i<(freeVectorEntryCount-1)) {
      if (coEnd == freeVector[i+1].coStart) {
	coEnd = freeVector[++i].coEnd;
      } else {
#ifdef DEBUG
	fprintf (stderr, "Unfree hole from co%lu to co%lu of %ludps\n",
		 coEnd, freeVector[i+1].coStart, (freeVector[i+1].coStart - coEnd));
#endif
	break;
      }      
    }
    totalFreeAfter += coEnd-coStart;
    insertFreeListHead (coStart, coEnd-coStart);
  }
  fprintf (stderr, "\rRe-inserting...100%% done\n");

  // Say how it went.
  fprintf (stderr, "Began with %ludps free and ended with %ludps\n",
	  totalFreeBefore, totalFreeAfter);
  statSSD ();
  free (freeVector);
}

/*
  Actually write the buffer.
*/

void flushSSD (double *buffer, unsigned long chunkOffset, unsigned long doublePairCount)
{
  putCallCount++;
  if ((putCallCount & 0x3FFFF) == 0x3FFFF)
    statSSD ();

  // Chunk in SSD file is assigned, put the buffer out there...
  if (fseek (sSDFD, chunkOffset * DOUBLE_PAIR_SIZE, SEEK_SET) != 0) {
    sprintf (messageBuffer, "Failed to seek to offset co%lu in SSD cache file",
	     chunkOffset);
    perror (messageBuffer);
    exit (EXIT_FAILURE);
  }
  if ((fwrite (buffer, DOUBLE_PAIR_SIZE, doublePairCount, sSDFD)) != doublePairCount) {
    sprintf (messageBuffer, "Failed to write %ludps at offset co%lu in SSD cache file",
	     doublePairCount, chunkOffset);
    perror (messageBuffer);
    exit (EXIT_FAILURE);
  }
}

/*
  Store a buffer of length doublePairCount and get a chunkTicket
  as your receipt.
*/

struct chunkTicket *doPutSSD (double *buffer, unsigned long myDPC, unsigned short freeList) {
  struct chunkTicket *newCT;
  unsigned short newFreeList;
  unsigned long leftOver;

  // Found a list, generate a chunkTicket...
  if ((newCT = (struct chunkTicket *) malloc (sizeof (struct chunkTicket))) == NULL) {
    perror ("Failed to allocate new chunkTicket");
    exit (EXIT_FAILURE);
  }
  newCT->chunkOffset = listHead[freeList].chunkOffset;
  newCT->doublePairCount = myDPC;

#ifdef DEBUG
  fprintf (stderr, " at offset co%lu OK\n", newCT->chunkOffset);
#endif

  // Handle leftovers...must be bigger than a freeList structure.
  if ((leftOver = listHead[freeList].doublePairCount - myDPC) > (MIN_USE_SSD_DPS + 1UL)) {
    // Put it on a free list...maybe the current one?
    if ((newFreeList = high16Bit (leftOver - 1UL)) == freeList) {
      // Keep the current list head, just smaller!
      listHead[freeList].chunkOffset += myDPC;
      listHead[freeList].doublePairCount -= myDPC;
    } else {
      // This space from this list head goes to the head of a lesser list
      insertFreeListHead (listHead[freeList].chunkOffset + myDPC, leftOver);
      removeFreeListHead (freeList);
    }
  } else {
    // No significant leftovers, give it all away.
    newCT->doublePairCount = listHead[freeList].doublePairCount;
    // This list head goes away completely due to no real leftovers
    removeFreeListHead (freeList);
  }
  flushSSD (buffer, newCT->chunkOffset, newCT->doublePairCount);
  usedDPCs += newCT->doublePairCount;
  handledDPCs += newCT->doublePairCount;
  return newCT;
}

/**

   Store a chunk of information on the SSD and get a ticket for
   it's retrieval.

   Accepts a pointer to a contiguous chunk of double-aligned memory
   up to maxSSDDPC*16 bytes long for storage on the SSD. Returns a
   pointer to a chunkTicket that can be passed to getSSD to retrieve
   that chunk of memory.

   @par Global Inputs

   none.

   @par Global Outputs

   none.

  @return A pointer to a chunkTicket structure redeemable for the
  chunk of data stored. Null if request cannot be fulfilled.

*/
struct chunkTicket *putSSD (
			    double *buffer, ///< Pointer to double-aligned buffer of data to be stored
			    unsigned long myDPC ///< Count of double pairs to be stored (X 16 = bytes)
			    ) {
  unsigned short freeList;
  unsigned short bestFreeList;

  if ((myDPC < MIN_USE_SSD_DPS) || (myDPC > maxSSDDPC)) {
    fprintf (stderr, "putSSD size of %ludps out-of-bounds\n", myDPC);
    return ((struct chunkTicket *) NULL);
  }

#ifdef DEBUG
  fprintf (stderr, "putSSD for %ludpc...", myDPC);
#endif

  // Use the head of the first available list that's large enough.
  bestFreeList  = MIN(15, high16Bit (myDPC));
  for (freeList = bestFreeList; freeList<16; freeList++) {
#ifdef DEBUG
    fprintf (stderr, "L%d has %ludpc ", freeList, listHead[freeList].doublePairCount);
#endif
    if (listHead[freeList].doublePairCount >= myDPC)
      return (doPutSSD (buffer, myDPC, freeList));
  }
  // Nothing left at all!!
  fprintf (stderr, "putSSD may be out of chunkage for a %ludps request (15 head has %ludps)!\n",
	   myDPC, listHead[15].doublePairCount);

  // Nothing, reorder the best free list and try again...
  /*
  reorderFreeList (bestFreeList);

  if (listHead[bestFreeList].doublePairCount >= myDPC) {
#ifdef DEBUG
    fprintf (stderr, "Best free list (%d) has %ludps\n", bestFreeList, listHead[bestFreeList].doublePairCount);
#endif
    return (doPutSSD (buffer, myDPC, bestFreeList));
  }
  */
  // That's about it...time to shuffle thru the trash.

  garbageCollect ();

  // Try again now...
  for (freeList = bestFreeList; freeList<16; freeList++) {
#ifdef DEBUG
    fprintf (stderr, "%d has %ludpc\n", freeList, listHead[freeList].doublePairCount);
#endif
    if (listHead[freeList].doublePairCount >= myDPC)
      return (doPutSSD (buffer, myDPC, freeList));
  }
  // Nothing left at all!!
  fprintf (stderr, "putSSD is completely out of chunkage for a %ludpc request!\n", myDPC);

  return ((struct chunkTicket *) NULL);
}

/**

  Present your chunkTicket to retrieve a copy of a chunk of data previously stored by putSSD.

  Accepts a pointer to a chunkTicket issued by putSSD and a pointer to a buffer.
  Copies the chunk of data associated with the chunkTicket into the provided buffer.
  The chunkTicket is still valid and can be used again to retrieve the same chunk of
  data as often as needed.

   @par Global Inputs

   none.

   @par Global Outputs

   none.

  @return void.

*/
void getSSD (
	     struct chunkTicket *myTicket, ///< A pointer to a chunk ticket issued by putSSD
	     double *buffer ///< The buffer into which the chunk data will be copied
	     ) {
  getCallCount++;
  if ((getCallCount & 0x3FFFF) == 0x3FFFF)
    statSSD ();

#ifdef DEBUG
  fprintf (stderr, "getSSD for %ludps at co%lu\n", myTicket->doublePairCount, myTicket->chunkOffset);
#endif
  // Simply seek and read...
  if (fseek (sSDFD, myTicket->chunkOffset * DOUBLE_PAIR_SIZE, SEEK_SET) != 0) {
    sprintf (messageBuffer, "Failed to seek to %ludps in SSD cache file", myTicket->chunkOffset);
    perror (messageBuffer);
    exit (EXIT_FAILURE);
  }
  if ((fread (buffer, DOUBLE_PAIR_SIZE, myTicket->doublePairCount, sSDFD)) != myTicket->doublePairCount) {
    perror ("Failed to read SSD cache file");
    exit (EXIT_FAILURE);
  }
  handledDPCs += myTicket->doublePairCount;
  return;
}

/**

  Release the resources associated with the chunkTicket, make it invalid.

  Accepts a pointer to a chunkTicket issued by putSSD. Identifies and releases all
  of the resources associated with that ticket for reuse, including the ticket itself.

   @par Global Inputs

   none.

   @par Global Outputs

   none.

  @return void.

*/
void freeSSD (
	      struct chunkTicket *myTicket ///< A pointer to a chunkTicket issued by putSSD.
) {
  freeCallCount++;
  if ((freeCallCount & 0x3FFFF) == 0x3FFFF)
    statSSD ();

#ifdef DEBUG
  fprintf (stderr, "freeSSD returning %ludps at co%lu\n", myTicket->doublePairCount, myTicket->chunkOffset);
#endif
  insertFreeListHead (myTicket->chunkOffset, myTicket->doublePairCount);
  usedDPCs -= myTicket->doublePairCount;
  free (myTicket);
  return;
}

#ifdef MAIN

// 536 is clean, 538 is not
#define TEST_ITERATIONS (1024 * 1024)

#include <time.h>

#define MAX_TICKET_MASK 0x7F

int main (int argc, char *argv[]) {
  
  struct chunkTicket *listOTickets[MAX_TICKET_MASK+2];
  int cTStatus[MAX_TICKET_MASK+2];
  int i, j, cT;
  unsigned long dPC;
  double buffer[(MAX_DPC_MASK + 1) * 4];

  initSSD ();

  fprintf (stderr, "Field of %d tickets of up to %ddps in size in file of %lddps\n",
	  MAX_TICKET_MASK, MAX_DPC_MASK, maxSSDDPC);

  //  srand(time(0)); // Let's go for repeatability

  for (i=0; i<=MAX_TICKET_MASK; i++)
    cTStatus[i] = 70;

  for (i=0; i<TEST_ITERATIONS; i++) {
    cT = rand() & MAX_TICKET_MASK;

#ifdef DEBUG
    fprintf (stderr, "%d, Ticket %d(%d):", i, cT, cTStatus[cT]);
#endif
    switch (cTStatus[cT]) {

    case 70: // Put it out there...
      dPC = rand() & MAX_DPC_MASK;
      dPC = MAX(MIN_USE_SSD_DPS, dPC);
      buffer[0] = dPC;
      for (j=1; j<dPC; j++) {
	buffer[j] = cT;
	buffer[dPC+j] = -cT;
      }

      listOTickets[cT] = putSSD (buffer, dPC);

      if (listOTickets[cT] != NULL)
	cTStatus[cT] = 71;
      else
	exit (EXIT_FAILURE);

      break;

    case 71: // Get it back and verify it...
      getSSD (listOTickets[cT], buffer);
      dPC = buffer[0];
      for (j=1; j<dPC; j++) {
	if ((buffer[j] != (double) cT) || (buffer[dPC+j] != (double) -cT)) {
	  fprintf (stderr, "At %dth position, wrote %g/%g at offset %ludps, got %g/%lu!\n",
		  j, (double) cT, (double) -cT, listOTickets[cT]->chunkOffset,
		  buffer[j], dPC+j);
	  fclose (sSDFD);
	  exit (EXIT_FAILURE);
	}
      }
      cTStatus[cT] = 72;
      break;

    case 72: // Release it...
      // Half the time we don't to this...
      if (rand() & 1) {
#ifdef DEBUG
	fprintf (stderr, "Skip free\n");
#endif
	break;
      }
      freeSSD (listOTickets[cT]);
      cTStatus[cT] = 70;
      listOTickets[cT] = 0;
      break;

    default:
      fprintf (stderr, "cTStatus has gone nuts with value of %d for index %d\n",
	       cTStatus[cT], cT);
      exit (EXIT_FAILURE);
    }
  }

  // Garbage collect to reduce to ticket count (hopefully)
  garbageCollect ();

  fprintf (stderr, "Freeing the rest...\n");
  for (i=0; i<=MAX_TICKET_MASK; i++) {
    if (cTStatus[i] != 70) {
#ifdef DEBUG
      fprintf (stderr, "%d, Ticket %d(%d):", i, cT, cTStatus[cT]);
#endif
      freeSSD (listOTickets[i]);
      cTStatus[i] = 70;
      listOTickets[i] = 0;
    }
  }

  // Now garbage collect, because it should be incredibly efficient!
  garbageCollect ();
  termSSD ();
}

#endif

#ifdef MAIN2

#define TEST_ITERATIONS (1024 * 1024)

#include <time.h>

int main (int argc, char *argv[]) {
  
  struct chunkTicket *ticket;
  double dPs[1024*1024];
  int i;

  initSSD ();

  for (i=0; i<=(8192*1024); i++)
    if ((ticket = putSSD (dPs, 1024UL)) == NULL)
      exit (EXIT_FAILURE);

  termSSD ();
}

#endif
