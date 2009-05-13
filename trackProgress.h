#include "pedlib/pedigree.h"

void *monitorStatus ();
void print_dryrun_stat (PedigreeSet * pSet, double pos);
void logPedigreeSetStatistics (PedigreeSet * pSet, int posIdx);
void dumpTrackingStats(unsigned long cl[], unsigned long eCl[]);
char *estimateIterations (unsigned long eCl[]);

extern unsigned long cL[9], eCL[9];
extern struct swStopwatch *combinedComputeSW, *combinedBuildSW, stopwatch;

