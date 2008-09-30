/**
@file sSDHandler.c

 sSDHandler - manage storage on a Solid State Drive for polynomial term lists.

 Written by Bill Valentine-Cooper.

 Copyright 2008, Nationwide Children's Research Institute.  All rights
 reserved.  Permission is hereby given to use this software for
 non-profit educational purposes only.

 We're out of memory, but have an SSD to help. Memory mapping a file on the SSD doesn't improve
 our situation because every page actually referenced has to be in physical memory. What we need
 is a "window" of fixed-memory size thru which we can access the full span of SSD space. That's
 normal file I/O. I searched and searched for a tool that would manage the space in the SSD 
 efficiently and came up with nothing, so I wrote this.

 Philosophy:

 Allocate space in the SSD cache file in double pair chunks, i.e. some count of pairs of doubles.
 Keep 16 lists of free chunks. List n contains only free chunks of size > 2^n double pairs.
 Start with only one chunk in the 16th list that is the entire file. List heads are in memory,
 while subsequent list entries are at the starting offset for that entry (as pointed to by the
 previous entry's next pointer). Always push or pop to/from the head of the list. Don't keep back
 pointers for now as it would double the number of writes, and would just make defragmentation
 easier, and defragmentation may not be necessary. Don't keep track of used memory -- leave that
 to the user.

 Usage:

 First call initSDD(), then call putSSD to store a chunk of double pairs and get a ticket. Call
 getSSD with the ticket to retrieve the chunk. Call freeSSD with the ticket to give-up the chunk.
 Call termSSD() when finished to display statistics.

 Performance:

 64K operations (putSSD, getSSD, freeSDD) of random sized chunks to to 64K in 47.5 seconds. 512K
 operations in 9 minutes.

 Caching aside (because I leave that to the OS), there is a fixed irreducible cost of 1 write 
 for a putSSD, 1 read for a getSSD. Anything beyond that is overhead. This algorithm's overhead is:

 -nothing for putSSD calls utilizing free lists where the remainder of the list head entry still
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

*/

#include "sSDHandler.h"

// The SSD path and an arbitrary filename
char *sSDFileName = "/tmp/ssd/cache.dat";
// Maximum number of double pairs on the 28Gb SSD is 28*1024*1024*1024/16, about 1792Mdps
#ifdef MAIN
#define MAX_SSD_DPC 512
#else
#define MAX_SSD_DPC (28 * 1024 * 1024 / 16 * 1024)
#endif

// Maximum size of a chunk of double pairs, 2^15=32K
#ifdef MAIN
#define MAX_DPC_MASK 0xF
#else
#define MAX_DPC_MASK 0x7FFF
#endif

#define DOUBLE_PAIR_SIZE (sizeof (double) * 2)

#define MIN(X,Y) (((X) > (Y)) ? (Y) : (X))
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))

char messageBuffer[256];

#ifdef MAIN
#define MAX_TICKET_MASK 0xF
struct chunkTicket *listOTickets[MAX_TICKET_MASK+1];
int cTStatus[MAX_TICKET_MASK+1];
#endif

FILE *sSDFD;

struct listEntry {
  unsigned long doublePairCount; /* Size in double pairs of this chunk, largest feasible is 
				    2^32-1*16 = ~64Gb, but really must be less than 2^16 (MAX_DPC_MASK) so
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
  if (value & 0xFFFF8000)
    return 15;
  for (i=14; i>0; i--)
    if (value & 0x4000)
      break;
    else
      value <<= 1;
  return i;
}

void statSSD () {
  int i;
  fprintf (stderr, "SSD list use: ");
  for (i=0; i<16; i++)
    fprintf (stderr, "%d:%d ", i, listDepth[i]);
  fprintf (stderr, "\nput/get/free calls: %d/%d/%d, total used: %lu, handled: %lu\n",
	   putCallCount, getCallCount, freeCallCount, usedDPCs, handledDPCs);
}

void termSSD () {
  fclose (sSDFD);
}

void initSSD() {
  int i;

  if ((sSDFD = fopen(sSDFileName,"wb+")) == NULL) {
    perror ("Failed to open SSD cache file");
    exit (EXIT_FAILURE);
  }
  for (i=0; i<16; i++)
    listHead[i].doublePairCount = listHead[i].chunkOffset = listHead[i].nextFree = listDepth[i] = 0;

  listHead[15].doublePairCount = MAX_SSD_DPC;
  listDepth[15] = 1;
  getCallCount = putCallCount = freeCallCount = 0;
  return;
}

void removeFreeListHead (unsigned short freeList) {
  listDepth[freeList]--;
  if (listHead[freeList].nextFree == 0) {
    // Out of list entries
    listHead[freeList].doublePairCount = 0;
  } else {
    // Pull-in the next entry and shlorp-up the listEntry space!
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
#ifdef DIAG
    fprintf (stderr, "Popped new head of list %d at %ludpc of %ludps, next at %ludpc\n", freeList,
	    listHead[freeList].chunkOffset,
	    listHead[freeList].doublePairCount,
	    listHead[freeList].nextFree);
#endif
  }
  return;
}

void insertFreeListHead (unsigned long chunkOffset, unsigned long doublePairCount) {
  unsigned short freeList;

  if (doublePairCount > 0)
    freeList = high16Bit (doublePairCount - 1);
  else
    freeList = 0;
  listDepth[freeList]++;
  if (listHead[freeList].doublePairCount == 0) {
    // First one, an easy insertion
    listHead[freeList].chunkOffset = chunkOffset;
    listHead[freeList].doublePairCount = doublePairCount;
    listHead[freeList].nextFree = 0;
  } else {
    /* BTW - Combining contiguous free entries to defragment dymaically would require
       backpointers, which would double our read/write burden. */
    // Not the first one, need to flush the old first one to the SSD cache file

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
    listHead[freeList].nextFree = listHead[freeList].chunkOffset - 1;
    listHead[freeList].chunkOffset = chunkOffset;
    listHead[freeList].doublePairCount = doublePairCount;

  }
  return;
}

struct freeVectorEntry {
  unsigned long startDPC;
  unsigned long endDPC;
};

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
  int i, j;
  unsigned long totalFreeBefore = 0, totalFreeAfter = 0, startDPC = 0, endDPC = 0, freeVectorEntryCount;
  struct freeVectorEntry *freeVector;

  fprintf (stderr, "sSDHandler garbage collection started...\n");
  statSSD ();

  // Collect all freeList entries regardless of length

  freeVectorEntryCount = 0;
  for (i=0; i<16; i++)
    freeVectorEntryCount += listDepth[i];

  if ((freeVector = (struct freeVectorEntry *) 
       malloc (freeVectorEntryCount * sizeof (struct freeVectorEntry))) == NULL) {
    fprintf (stderr, "Failed to allocate a %lu byte freeVector for garbage collection\n",
	     freeVectorEntryCount * sizeof (struct freeVectorEntry));
    exit (EXIT_FAILURE);
  }

  freeVectorEntryCount = 0;
  for (i=0; i<16; i++)
    while (listHead[i].doublePairCount != 0) {
      freeVector[freeVectorEntryCount].startDPC = listHead[i].chunkOffset;
      freeVector[freeVectorEntryCount].endDPC = listHead[i].chunkOffset + listHead[i].doublePairCount;
      freeVectorEntryCount++;
      totalFreeBefore += listHead[i].doublePairCount;
      removeFreeListHead (i);
    }

  /* Sort 'em. I'd use qsort, but it wants a vector of structure pointers, and I just allocated the
     whole bleedin' thing at once. */

  struct freeVectorEntry swap;
  for (i=0; i<freeVectorEntryCount; i++)
    for (j=i; j<freeVectorEntryCount-1; j++)
      if (freeVector[j].startDPC > freeVector[j+1].startDPC) {
	swap.startDPC = freeVector[j].startDPC;
	swap.endDPC = freeVector[j].endDPC;
	freeVector[j].startDPC = freeVector[j+1].startDPC;
	freeVector[j].endDPC = freeVector[j+1].endDPC;
	freeVector[j+1].startDPC = swap.startDPC;
	freeVector[j+1].endDPC = swap.endDPC;
      }

  // Reset listHead and listDepth

  for (i=0; i<16; i++)
    listHead[i].doublePairCount = listHead[i].chunkOffset = listHead[i].nextFree = listDepth[i] = 0;

  // Connect and insert them.

  for (i=0; i<freeVectorEntryCount; i++) {
    startDPC = freeVector[i].startDPC;
    endDPC = freeVector[i].endDPC;
    while ((i<(freeVectorEntryCount-1)) && (endDPC == freeVector[i+1].startDPC))
      endDPC = freeVector[++i].endDPC;
    totalFreeAfter += endDPC-startDPC;
    insertFreeListHead (startDPC, endDPC-startDPC);
  }

  // Say how it went.
  fprintf (stderr, "Began with %ludpc free and ended with %ludpc\n",
	  totalFreeBefore, totalFreeAfter);
  statSSD ();
  free (freeVector);
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

#ifdef DIAG
  fprintf (stderr, " at offset %ludpc OK\n", newCT->chunkOffset);
#endif

  // Handle leftovers...must be bigger than a freeList structure.
  if ((leftOver = listHead[freeList].doublePairCount - myDPC - 1) > MIN_USE_SSD) {
    // Something leftover, put it on a free list...maybe the current one?
    if ((newFreeList = high16Bit (leftOver)) == freeList) {
      // Keep the current list head, just smaller!
      listHead[freeList].chunkOffset += myDPC;
      listHead[freeList].doublePairCount -= myDPC;
    } else {
      // This space from this list head goes to the head of a lesser list
      insertFreeListHead (listHead[freeList].chunkOffset + myDPC, leftOver + 1);
      removeFreeListHead (freeList);
    }
  } else {
    // No significant leftovers, give it all away.
    newCT->doublePairCount = listHead[freeList].doublePairCount;
    // This list head goes away completely due to no real leftovers
    removeFreeListHead (freeList);
  }
  // Chunk in SSD file is assigned, put the buffer out there...
  if (fseek (sSDFD, newCT->chunkOffset * DOUBLE_PAIR_SIZE, SEEK_SET) != 0) {
    sprintf (messageBuffer, "Failed to seek to offset %ludpc in SSD cache file",
	     newCT->chunkOffset);
    perror (messageBuffer);
    exit (EXIT_FAILURE);
  }
  if ((fwrite (buffer, DOUBLE_PAIR_SIZE, newCT->doublePairCount, sSDFD)) != newCT->doublePairCount) {
    sprintf (messageBuffer, "Failed to write %ludps at offset %ludpc in SSD cache file",
	     newCT->doublePairCount, newCT->chunkOffset);
    perror (messageBuffer);
    exit (EXIT_FAILURE);
  }
  usedDPCs += newCT->doublePairCount;
  handledDPCs += newCT->doublePairCount;
  return newCT;
}

struct chunkTicket *putSSD (double *buffer, unsigned long myDPC) {
  unsigned short freeList;

  if ((myDPC < MIN_USE_SSD) || (myDPC > MAX_SSD_DPC)) {
    fprintf (stderr, "putSSD size of %ludpc out-of-bounds\n", myDPC);
    return ((struct chunkTicket *) NULL);
  }

  putCallCount++;
  if ((putCallCount & 0x3FFFF) == 0x3FFFF)
    statSSD ();

#ifdef DIAG
  fprintf (stderr, "putSSD for %ludpc...", myDPC);
#endif

  // Use the head of the first available list that's large enough.
  for (freeList = MIN(15, high16Bit (myDPC) +1); freeList<16; freeList++) {
#ifdef DIAG
    fprintf (stderr, "%d has %ludpc ", freeList, listHead[freeList].doublePairCount);
#endif
    if (listHead[freeList].doublePairCount != 0)
      return (doPutSSD (buffer, myDPC, freeList));
  }
  // Nothing left at all!!
  fprintf (stderr, "putSSD may be out of chunkage for a %ludpc request!\n", myDPC);
  garbageCollect ();

  // Try again now...
  for (freeList = MIN(15, high16Bit (myDPC) +1); freeList<16; freeList++) {
#ifdef DIAG
    fprintf (stderr, "%d has %ludpc ", freeList, listHead[freeList].doublePairCount);
#endif
    if (listHead[freeList].doublePairCount != 0)
      return (doPutSSD (buffer, myDPC, freeList));
  }
  // Nothing left at all!!
  fprintf (stderr, "putSSD is completely out of chunkage for a %ludpc request!\n", myDPC);

  return ((struct chunkTicket *) NULL);
}

/*
  Retrieve a buffer previously stored with putSSD by
  presenting your chunkTicket. chunkTicket is still valid.
*/
void getSSD (struct chunkTicket *myTicket, double *buffer) {
  getCallCount++;
#ifdef DIAG
  fprintf (stderr, "getSSD for %ludps at %ludpc\n", myTicket->doublePairCount, myTicket->chunkOffset);
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

/*
  Release the resources identified by the chunkTicket,
  which will no longer be valid.
*/
void freeSSD (struct chunkTicket *myTicket) {
  freeCallCount++;
#ifdef DIAG
  fprintf (stderr, "freeSSD returning %ludpc at %ludps\n", myTicket->doublePairCount, myTicket->chunkOffset);
#endif
  insertFreeListHead (myTicket->chunkOffset, myTicket->doublePairCount);
  usedDPCs -= myTicket->doublePairCount;
  free (myTicket);
  return;
}

#ifdef MAIN

#define TEST_ITERATIONS (1024 * 1024)

#include <time.h>

int main (int argc, char *argv[]) {
  
  int i, j, cT;
  unsigned long dPC;
  double buffer[(MAX_DPC_MASK * 2) + 1];

  initSSD ();

  fprintf (stderr, "Field of %d tickets of up to %ddps in size in file of %ddps\n",
	  MAX_TICKET_MASK, MAX_DPC_MASK, MAX_SSD_DPC);

  srand(time(0));

  for (i=0; i<=MAX_TICKET_MASK; i++)
    cTStatus[i] = 70;

  for (i=0; i<TEST_ITERATIONS; i++) {
    cT = rand() & MAX_TICKET_MASK;

    switch (cTStatus[cT]) {

    case 70: // Put it out there...
      dPC = rand() & MAX_DPC_MASK;
      dPC = MAX(MIN_USE_SSD, dPC);
#ifdef DIAG
      fprintf (stderr, "Ticket %d putSSD %ludpc\n", cT, dPC);
#endif
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
#ifdef DIAG
      fprintf (stderr, "Ticket %d getSSD %ludpc\n", cT, listOTickets[cT]->doublePairCount);
#endif
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
#ifdef DIAG
      fprintf (stderr, "Ticket %d freeSSD %ludpc\n", cT, listOTickets[cT]->doublePairCount);
#endif
      freeSSD (listOTickets[cT]);
      cTStatus[cT] = 70;
      break;

    default:
      fprintf (stderr, "cTStatus has gone nuts with value of %d for index %d\n",
	       cTStatus[cT], cT);
      exit (EXIT_FAILURE);
    }

  }
  termSSD ();
}

#endif
