#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct chunkTicket {
  unsigned long doublePairCount; // Size in double pairs of this chunk
  unsigned long chunkOffset; /* Offset in double pairs from beginning of file, max of 
				64Gb/16b = 4M, so largest usable SSD is 4M * 32K = 128G. */
};

void initSSD ();
void statSSD ();
void termSSD ();

// The smallest chunk to split, in double pairs, trying 2^8 or 256 (2^MIN_USE_SSD_BITS)
#ifdef MAIN
#define MIN_USE_SSD_BITS 8
#else
#define MIN_USE_SSD_BITS 8
#endif
#define MIN_USE_SSD_DPS (1<<MIN_USE_SSD_BITS)

// Maximum size of a chunk, in double pairs, 2^15=32K, but we really can go beyond that
#ifdef MAIN
#define MAX_DPC_MASK 0xFFFF
#else
#define MAX_DPC_MASK 0x7FFFF
#endif


struct chunkTicket *putSSD (double *buffer, unsigned long myDPC);
void getSSD (struct chunkTicket *myTicket, double *buffer);
void freeSSD (struct chunkTicket *myTicket);
void flushSSD (double *buffer, unsigned long chunkOffset, unsigned long doublePairCount);

