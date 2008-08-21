#include "pedlib/pedigree.h"

void *monitorStatus ();
void print_dryrun_stat (PedigreeSet * pSet, double pos);
void logPedigreeSetStatistics (PedigreeSet * pSet, int posIdx);
void dumpTrackingStats(int cl[], int eCl[]);
char *estimateIterations (int eCl[]);
void pushStatus (char *);
void popStatus ();
