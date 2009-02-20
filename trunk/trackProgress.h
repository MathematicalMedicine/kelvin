#include "pedlib/pedigree.h"

void *monitorStatus ();
void print_dryrun_stat (PedigreeSet * pSet, double pos);
void logPedigreeSetStatistics (PedigreeSet * pSet, int posIdx);
void dumpTrackingStats(unsigned long cl[], unsigned long eCl[]);
char *estimateIterations (unsigned long eCl[]);
