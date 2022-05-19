#ifndef TRACKPROGRESS_H
#define TRACKPROGRESS_H

/**
@file trackProgress.h

  Functions related to tracking run progress and status display.

  Copyright (C) 2008, 2010, 2022 Mathematical Medicine LLC
  This program is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.
  You should have received a copy of the GNU General Public License along
  with this program. If not, see <https://www.gnu.org/licenses/>.

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
