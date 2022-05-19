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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct chunkTicket {
  unsigned long doublePairCount; // Size in double pairs of this chunk
  unsigned long chunkOffset; /* Offset in double pairs from beginning of file, max of 
				64Gb/16b = 4M, so largest usable SSD is 4M * 32K = 128G. */
};

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

void initSSD ();
struct chunkTicket *putSSD (double *buffer, unsigned long myDPC);
void getSSD (struct chunkTicket *myTicket, double *buffer);
void freeSSD (struct chunkTicket *myTicket);
void statSSD ();
void termSSD ();
void flushSSD (double *buffer, unsigned long chunkOffset, unsigned long doublePairCount);

