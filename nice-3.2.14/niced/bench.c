/**********************************************************************
 * NICE Benchmark Code
 * Alberto Maria Segre
 *
 * Copyright 2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software 
 * for non-profit educational purposes only.
 **********************************************************************/
#include "niced.h"

/**********************************************************************
 * Generate an array of BENCHMARKSZE unsigned short integers, randomly
 * permute their order, and then use insertion sort to sort
 * them. Repeat BENCHMARKREP times and use the elapsed time in
 * milliseconds as a crude relative measure of machine
 * load/performance and priority.
 **********************************************************************/
#define SWAP(i,j) {temp=array[i]; array[i]=array[j]; array[j]=temp;}
long int
bench ()
{
  unsigned int i, j, count;
  unsigned short int temp, array[BENCHMARKSZE];
  clock_t time;
  long int elapsed;

  /* Start generating values between 0 and USHRT_MAX. */
  for (i = 0; i < BENCHMARKSZE; i++)
    array[i] = (unsigned short int) (i % USHRT_MAX);

  /* Start the clock. */
  time = clock ();

  for (count = 0; count < BENCHMARKREP; count++)
    {
      /* Permute. */
      for (i = 0; i < BENCHMARKSZE; i++)
	for (j = i; j < BENCHMARKSZE; j++)
	  if (RANDINT (2))
	    SWAP (i, j);

      /* Insertion sort. */
      for (i = 1; i < BENCHMARKSZE; i++)
	{
	  j = i;
	  temp = array[j];
	  while (j > 0 && array[j - 1] > temp)
	    {
	      array[j] = array[j - 1];
	      j--;
	    }
	  array[j] = temp;
	}
    }

  /* Stop the clock. Note that resolution of timer is only about 10 usec. */
  elapsed = (((long int) ((clock () - time) * 1000)) / CLOCKS_PER_SEC);

  /* As somewhat of a hack, we're going to factor priority in here as
   * well. By using priority as a muiltiplier, we can percolate
   * machines of like load but running NICE at lower priority
   * downwards in the hierarchy. The problem is, there's no good way
   * to scale priority that I know of, so whatever we do is, ahem,
   * "heuristic." Our current heuristic makes "top" priority worth a
   * doubling in perceived speed, and a "bottom" priority worth a 50%
   * speed penalty. */
  return (elapsed + priority * (elapsed/40));
}
