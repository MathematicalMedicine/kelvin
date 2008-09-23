#include "sSDHandler.h"

FILE *sSDFD;
char *sSDFileName = "/tmp/ssd/cache.dat";
#define MAX_DPC_FWRITE 0x7FFF
#define MIN(X,Y) ((X) <= (Y) ? (X) : (Y))

struct chunkTicket {
  unsigned long doublePairCount; // Size in double pairs of this chunk, largest is 2^32-1*16 = ~64Gb
  unsigned long chunkOffset; // Offset in double pairs from beginning of file, max of 64Gb/16b = 4M.
};

/*
  The leftmost set bit of the doublePairCount picks the free list. The size of this (16 bytes) IS
  HARDCODED everywhere as 1 doublePair and 16 bytes.
*/
struct listEntry {
  unsigned long doublePairCount; // Size in double pairs of this chunk, largest is 2^32-1*16 = ~64Gb, 0 if end of list
  unsigned long chunkOffset; // Offset in double pairs from beginning of file, max of 64Gb/16b = 4M.
  unsigned long nextFree; // Offset in double pairs of next listEntry in free list, 0 if DNE.
};
struct listEntry listHead[16];
int listDepth[16];

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

int initSSD() {
  int i;

  if ((sSDFD = fopen(sSDFileName,"wb+")) == NULL) {
    perror ("Failed to open SSD cache file");
    exit (EXIT_FAILURE);
  }
  for (i=0; i<16; i++)
    listHead[i].doublePairCount = listHead[i].chunkOffset = listHead[i].nextFree = listDepth[i] = 0;

  listHead[15].doublePairCount = 0x7FFFFFFF; // For 32Gb / 16b, or 2Gb
  return;
}

void removeFreeListHead (unsigned short freeList) {
  if (listHead[freeList].nextFree == 0) {
    // Out of list entries
    listHead[freeList].doublePairCount = 0;
    printf ("Out of entries for list %d\n", freeList);
  } else {
    // Pull-in the next entry and shlorp-up the listEntry space!
    if (fseek (sSDFD, listHead[freeList].nextFree * 16, SEEK_SET) != 0) {
      perror ("Failed to seek in SSD cache file");
      exit (EXIT_FAILURE);
    }
    if ((fread (&listHead[freeList], sizeof (struct listEntry), 1, sSDFD)) != 1) {
      perror ("Failed to read SSD cache file");
      exit (EXIT_FAILURE);
    }
    listHead[freeList].chunkOffset--;
    listHead[freeList].doublePairCount++;
    listDepth[freeList]--;
    printf ("Popped new head of list %d at %ludpc of %ludps, next at %ludpc\n", freeList,
	    listHead[freeList].chunkOffset,
	    listHead[freeList].doublePairCount,
	    listHead[freeList].nextFree);
  }
  return;
}

void insertFreeListHead (unsigned long chunkOffset, unsigned long doublePairCount) {
  unsigned short freeList;

  freeList = high16Bit (doublePairCount - 1);
  printf ("insertFreeListHead in list %d for chunk at offset %lu of size %ludps\n",
	  freeList, chunkOffset, doublePairCount);
  if (listHead[freeList].doublePairCount == 0) {
    // First one, an easy insertion
    listHead[freeList].chunkOffset = chunkOffset;
    listHead[freeList].doublePairCount = doublePairCount;
    listHead[freeList].nextFree = 0;
  } else {
    // TBS - COMBINE CONTIGUOUS FREE ENTRIES TO DEFRAGMENT.
    // Not the first one, need to flush the old first one to the SSD cache file
    listDepth[freeList]++;
    printf ("Pushing old head of list %d at %ludpc of %ludps, next at %ludpc\n", freeList,
	    listHead[freeList].chunkOffset,
	    listHead[freeList].doublePairCount,
	    listHead[freeList].nextFree);
    if (fseek (sSDFD, listHead[freeList].chunkOffset * 16, SEEK_SET) != 0) {
      perror ("Failed to seek in SSD cache file");
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

/*
  Store a buffer of length doublePairCount and get a chunkTicket
  as your receipt.
*/
struct chunkTicket *putSSD (double *buffer, unsigned long myDPC) {
  struct chunkTicket *newCT;
  unsigned short freeList, newFreeList;
  unsigned long leftOver;

  printf ("putSSD for %ddps\n", myDPC);
  // Use the head of the first available list that's large enough.
  for (freeList = MIN(15, high16Bit (myDPC)+1); freeList<16; freeList++) {
    if (listHead[freeList].doublePairCount != 0) {
      printf ("...using freeList %d head of size %ludps\n", freeList, listHead[freeList].doublePairCount);
      // Found a list, generate a chunkTicket...
      if ((newCT = (struct chunkTicket *) malloc (sizeof (struct chunkTicket *))) == NULL) {
	perror ("Failed to allocate new chunkTicket");
	exit (EXIT_FAILURE);
      }
      newCT->chunkOffset = listHead[freeList].chunkOffset;
      newCT->doublePairCount = myDPC;
      printf ("...new chunk will be offset %ludpc\n", newCT->chunkOffset);
      // Handle leftovers...must be bigger than a freeList structure.
      if ((leftOver = listHead[freeList].doublePairCount - myDPC - 1) > 0) {
	printf ("...have %ludps left over\n", leftOver);
	printf ("...nCT->dPC is %ludps\n", newCT->doublePairCount);
	// Something leftover, put it on a free list...maybe the current one?
	if ((newFreeList = high16Bit (leftOver)) == freeList) {
	  // Keep the current list head, just smaller!
	  printf ("...keeping on current freeList %d\n", freeList);
	  printf ("...nCT->dPC is %ludps\n", newCT->doublePairCount);
	  listHead[freeList].chunkOffset += myDPC;
	  listHead[freeList].doublePairCount -= myDPC;
	} else {
	  // This space from this list head goes to the head of a lesser list
	  printf ("...moving to lesser freeList\n");
	  insertFreeListHead (listHead[freeList].chunkOffset + myDPC, leftOver + 1);
	  removeFreeListHead (freeList);
	}
      } else {
	// No significant leftovers, give it all away.
	printf ("...%ddps left over is not worth it, give it all away\n", leftOver);
	newCT->doublePairCount = listHead[freeList].doublePairCount;
	// This list head goes away completely due to no real leftovers
	removeFreeListHead (freeList);
      }
      // Chunk in SSD file is assigned, put the buffer out there...
      printf ("Seek to offset %ludps\n", newCT->chunkOffset * 16);
      if (fseek (sSDFD, newCT->chunkOffset * 16, SEEK_SET) != 0) {
	perror ("Failed to seek in SSD cache file");
	exit (EXIT_FAILURE);
      }
      printf ("Write %ludps of size 16\n", newCT->doublePairCount);
      if ((fwrite (buffer, 16, newCT->doublePairCount, sSDFD)) != newCT->doublePairCount) {
	perror ("Failed to write SSD cache file");
	exit (EXIT_FAILURE);
      }
      return newCT;
    }
  }
  // Nothing left at all!!
  fprintf (stderr, "putSSD is completely out of chunkage!\n");
  exit (EXIT_FAILURE);
}

/*
  Retrieve a buffer previously stored with putSSD by
  presenting your chunkTicket. chunkTicket is still valid.
*/
void getSSD (struct chunkTicket *myTicket, double *buffer) {
  printf ("getSSD for ticket with offset %ludpc and count %ddps\n", myTicket->chunkOffset, myTicket->doublePairCount);
  // Simply seek and read...
  printf ("Seek to offset %ludps\n", myTicket->chunkOffset);
  if (fseek (sSDFD, myTicket->chunkOffset * 16, SEEK_SET) != 0) {
    perror ("Failed to seek in SSD cache file");
    exit (EXIT_FAILURE);
  }
  printf ("Read %ludps of size 16\n", myTicket->doublePairCount);
  if ((fread (buffer, 16, myTicket->doublePairCount, sSDFD)) != myTicket->doublePairCount) {
    perror ("Failed to read SSD cache file");
    exit (EXIT_FAILURE);
  }
  return;
}

/*
  Release the resources identified by the chunkTicket,
  which will no longer be valid.
*/
void freeSSD (struct chunkTicket *myTicket) {
  printf ("freeSSD for ticket with offset %ludpc and count %ddps\n", myTicket->chunkOffset, myTicket->doublePairCount);
  insertFreeListHead (myTicket->chunkOffset, myTicket->doublePairCount);
  free (myTicket);
  return;
}


int main (int argc, char *argv[]) {
  struct chunkTicket *listOTickets[65536];
  int cTStatus[65536];

  double buffer[131072];
  int i, j, cT;
  unsigned long dPC;

  printf ("A listEntry takes %d bytes, and a nextFree takes %d\n",
	  sizeof (struct listEntry), sizeof (listHead[0].nextFree));

  initSSD ();

  srand(time(0));

  for (i=0; i<65536; i++) {
    cT = rand() & 0xFFFF;
    printf ("Working on ticket %d with status %d\n", cT, cTStatus[cT]);
    switch (cTStatus[cT]) {
    case 0: // Put it out there...
      dPC = rand() & 0x7FFF;
      for (j=0; j<dPC; j++) {
	buffer[j] = cT;
	buffer[dPC+j] = -cT;
      }
      listOTickets[cT] = putSSD (buffer, dPC);

      // Read it right back in to verify...
      getSSD (listOTickets[cT], buffer);
      printf ("Checking %ddps for cT %d\n", listOTickets[cT]->doublePairCount, cT);
      for (j=0; j<listOTickets[cT]->doublePairCount; j++) {
	if ((buffer[j] != (double) cT) || (buffer[listOTickets[cT]->doublePairCount+j] != (double) -cT)) {
	  printf ("At %dth position, wrote %g/%g at offset %ludps, got %g/%g!\n",
		  j, (double) cT, (double) -cT, listOTickets[cT]->chunkOffset,
		  buffer[j], buffer[listOTickets[cT]->doublePairCount+j]);
	  fclose (sSDFD);
	  exit (EXIT_FAILURE);
	}
      }
      printf ("OK!\n");

      cTStatus[cT] = 1;
      break;

    case 1: // Get it back and verify it...
      getSSD (listOTickets[cT], buffer);
      printf ("Checking %ddps for cT %d\n", listOTickets[cT]->doublePairCount, cT);
      for (j=0; j<listOTickets[cT]->doublePairCount; j++) {
	if ((buffer[j] != (double) cT) || (buffer[listOTickets[cT]->doublePairCount+j] != (double) -cT)) {
	  printf ("At %dth position, wrote %g/%g at offset %ludps, got %g/%g!\n",
		  j, (double) cT, (double) -cT, listOTickets[cT]->chunkOffset,
		  buffer[j], buffer[listOTickets[cT]->doublePairCount+j]);
	  fclose (sSDFD);
	  exit (EXIT_FAILURE);
	}
      }
      printf ("OK!\n");

      cTStatus[cT] = 2;
      break;

    case 2: // Release it...
      freeSSD (listOTickets[cT]);
      cTStatus[cT] = 0;
      break;
    }
  }

  for (i=0; i<16; i++)
    printf ("listDepth[%d] was %d\n", i, listDepth[i]);
  fclose (sSDFD);
  exit (EXIT_SUCCESS);
}
