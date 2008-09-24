#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// The SSD path and an arbitrary filename
char *sSDFileName = "/tmp/ssd/cache.dat";
// Maximum number of double pairs on the 28Gb SSD is 28*1024*1024*1024/16, about 1792Mdps
#define MAX_SSD_DPC (28 * 1024 * 1024 / 16 * 1024)

struct chunkTicket {
  unsigned long doublePairCount; // Size in double pairs of this chunk, max of 2^15 = 32K.
  unsigned long chunkOffset; /* Offset in double pairs from beginning of file, max of 
				64Gb/16b = 4M, so largest usable SSD is 4M * 32K = 128G. */
};

void initSSD ();
void termSSD ();

struct chunkTicket *putSSD (double *buffer, unsigned long myDPC);
void getSSD (struct chunkTicket *myTicket, double *buffer);
void freeSSD (struct chunkTicket *myTicket);

