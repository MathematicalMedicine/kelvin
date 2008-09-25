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

struct chunkTicket *putSSD (double *buffer, unsigned long myDPC);
void getSSD (struct chunkTicket *myTicket, double *buffer);
void freeSSD (struct chunkTicket *myTicket);

