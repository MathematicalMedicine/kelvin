/**
@mainpage

  kelvin - Linkage and Linkage Disequilibrium Analysis Program.

  Implementation of the Elston-Stewart algorithm for linkage
  analysis. Currently supports two-point and multipoint analyses,
  dichotomous and quantitative traits, linkage equilibrium and
  disequilibrium, case control, and many other options.
  
  Copyright 2008, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  See http://hodgkin.ccri.net/software/kelvin/index.html for full
  documentation.

  AUTHORS

  - Yungui Huang - principle author
  - Hongling Wang - Polynomial features
  - Alberto Maria Segre - config and error logging modules
  - Nathan Burnette - Regex code
  - Sang-Cheol Seok - Dynamic grid
  - John Burian - configuration parser
  - Bill Valentine-Cooper - polynomial compilation, SSD, performance and tracking

@file kelvin.c

  Main driver and way too many loops.

  USAGE:
  <pre>
  kelvin <kelvin.conf>

  where <kelvin.conf> is a text file containing directives
  describing the locations of supporting files and specifying the
  nature of the analysis to be performed. See the provided documentation
  for details.
  </pre>
  LIMITATIONS

  Currently only handles biallelic traits.

  COMPILE-TIME CONDITIONALS

  There are numerous compilation conditions that affect kelvin. Many
  of them are diagnostic in nature.
  
  - DMTRACK - enable exhaustive memory management tracking. Keeps track
  of all allocations and frees by source code line, and can provide
  cross-references at any time by total utilization or retention. Can
  also report on cross-module operations, i.e. one module allocates and
  another module frees the same chunk of memory.

  - DMUSE
  - FAKEEVALUATE
  - GCOV
  - GPROF

  - MEMGRAPH - provides the same information as MEMSTATUS but in a
  format appropriate for graphing with a tool like gnuplot. The file
  is named kelvin_<pid>_memory.dat, where <pid> is the process id.
  These files will need to be cleaned-up periodically, so this 
  conditional is not on by default. This is only possible on systems
  with a supported pmap command.

  - MEMSTATUS - handy when there is concern about memory capacity,
  but can really clutter-up the display. Runs in a child process and
  reports on the main process' memory usage every 30 seconds. This
  is only possible on systems with a supported pmap command.

  - _OPENMP - this conditional is not defined by the user, but rather
  set by the compiler when the -fopenmp flag is specified to enable
  OpenMP multithreading. Aspects of polynomial build and evaluation
  are currently setup to handle multithreading.  As a caveat --
  polynomial evaluation performance is only improved with
  multi-threading if the analysis has more unique pedigrees than
  threads. Performance will actually degrade when there are only a few
  unique pedigrees because term optimization will cause their
  polynomials to share most of their terms, and then per-pedigree
  evaluation threads will get into extensive synchronization conflicts
  traversing those polynomials. You can monitor this by watching the
  System CPU time as opposed to User CPU time. System CPU time is
  almost exclusively wasted time due to thread synchronization
  conflicts. Note that it is not, however, necessary to disable
  multi-threading by rebuilding kelvin for these situations. Simply
  set the OMP_NUM_THREADS environment variable to 1, and no locking
  conflicts will occur.

  - POLYSTATISTICS - enable or disable display of extensive polynomial
  statistics at 2Mp creation intervals as well as at key points in
  pedigree processing. This does not affect performance or processing.

  - SIMPLEPROGRESS - suppresses more complicated progress reports in
  favor of a single percentage progress and estimated remaining time.
  This simple progress indication cannot be linear, but is easier to
  understand. This conditional is on by default.

  - SOURCEDIGRAPH

*/
#include <pthread.h>
#include "kelvin.h"
#include "kelvinHandlers.h"
#include "pedlib/likelihood.h"
#include "saveResults.h"
#include "ppl.h"
#ifdef USE_MPF
#include <gmp.h>                /* GNU Multi-Precision library. */
#endif

extern char *likelihoodVersion, *locusVersion, *polynomialVersion;
extern Polynomial *constant1Poly;

/**

  Global variables

*/
#include "dcuhre.h"
#include "integrationGlobals.h"
#include "kelvinGlobals.h"
#include "kelvinLocals.h"

/**

  Driver for all types of analyses.

  <pre>
  Usage:
     kelvin kelvin.conf
  </pre>
  The kelvin.conf file gives information about the specific linkage
  analysis run. All information about, e.g., which markers to use,
  what outputs to calculate, and so on, are stored in this
  configuration file.

*/
char *programVersion = "V0.38.0";       ///< Overall kelvin version set upon release.
char *kelvinVersion = "$Id$";        ///< svn's version for kelvin.c


int main (int argc, char *argv[])
{

  //  #include "iterationLocalsNew.h"
  //  #include "integrationLocals.h"

  //#include "kelvinWriteFiles.c"
  #include "dkelvinWriteFiles.c"

  #include "kelvinInit.c"

  if (modelOptions.integration) {
    #include "integrationMain.c"
  } else {
    int danglingWire=0;
    //#include "iterationMain.c"
    iterateMain();
  }

  #include "kelvinTerm.c"

  return EXIT_SUCCESS;
}

#include "integrationSupport.c"
