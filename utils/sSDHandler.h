#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct chunkTicket {
  unsigned long doublePairCount; // Size in double pairs of this chunk, max of 2^15 = 32K.
  unsigned long chunkOffset; /* Offset in double pairs from beginning of file, max of 
				64Gb/16b = 4M, so largest usable SSD is 4M * 32K = 128G. */
};

void initSSD ();
void statSSD ();
void termSSD ();

// The smallest chunk to split, trying 257
#ifdef MAIN
#define MIN_USE_SSD 0
#else
#define MIN_USE_SSD 257
#endif

struct chunkTicket *putSSD (double *buffer, unsigned long myDPC);
void getSSD (struct chunkTicket *myTicket, double *buffer);
void freeSSD (struct chunkTicket *myTicket);

