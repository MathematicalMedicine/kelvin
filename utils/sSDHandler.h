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

// The smallest chunk to split, in double pairs trying 256
#ifdef MAIN
#define MIN_USE_SSD 2
#else
#define MIN_USE_SSD 256
#endif

// Maximum size of a chunk, in double pairs, 2^15=32K
#ifdef MAIN
#define MAX_DPC_MASK 0xF
#else
#define MAX_DPC_MASK 0x7FFF
#endif


struct chunkTicket *putSSD (double *buffer, unsigned long myDPC);
void getSSD (struct chunkTicket *myTicket, double *buffer);
void freeSSD (struct chunkTicket *myTicket);

