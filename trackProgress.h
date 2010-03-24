#ifndef TRACKPROGRESS_H
#define TRACKPROGRESS_H

/**
@file trackProgress.h

  Functions related to tracking run progress and status display.

  Copyright &copy; 2010, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

*/
#include "pedlib/pedigree.h"

extern unsigned long cL[9], ///< Actual compute_likelihood call counts
  eCL[9]; ///< Estimated compute_likelihood call counts

/// Number of seconds to delay between updates of memory status
#define MONSTATDELAYSEC 30

void *monitorMemory ();
void print_dryrun_stat (PedigreeSet * pSet, double pos);
void logPedigreeSetStatistics (PedigreeSet * pSet, int posIdx);
void dumpTrackingStats(unsigned long cl[], unsigned long eCl[]);
void estimateIterations (unsigned long eCl[]);

#endif
