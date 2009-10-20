/**
@mainpage

  kelvin - Linkage and Linkage Disequilibrium Analysis Program.

  Implementation of the Elston-Stewart algorithm for linkage
  analysis. Currently supports two-point and multipoint analyses,
  dichotomous and quantitative traits, linkage equilibrium and
  disequilibrium, case control, and many other options.

  Copyright &copy; 2009, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @see http://hodgkin.ccri.net/software/kelvin/index.html for full
  documentation.

  @author Yungui Huang - principle author.
  @author Hongling Wang - Polynomial features.
  @author Alberto Maria Segre - config and error logging modules.
  @author Nathan Burnette - Regex code (now removed).
  @author Sang-Cheol Seok - Dynamic grid.
  @author John Burian - configuration parser, sequential update.
  @author Bill Valentine-Cooper - additional polynomial features,
  compilation, SSD, refactoring, performance and tracking.

@file kelvin.c

  Genetic linkage analysis with PPL that integrates over entire trait space.

  USAGE:
  <pre>
  kelvin <config file> [--directive arg1 arg2... [--directive...]]

  where <config file> is a text file containing directives
  describing the locations of supporting files and specifying the
  nature of the analysis to be performed. See the provided documentation
  for details. Any directives specified on the command line override
  those specified in the configuration file.
  </pre>
  LIMITATIONS

  Currently only handles biallelic traits.

  @version $Id$

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

/**

  Global variables

*/
#include "kelvinGlobals.h"
#include "kelvinLocals.h"
#include "kelvinInit.h"
#include "kelvinTerm.h"
#include "integrationSupport.h"
#include "iterationSupport.h"

/**

  Driver for all types of analyses.

*/
char *programVersion = "V0.38.0";       ///< Overall kelvin version set upon release.
char *kelvinVersion = "$Id$";        ///< svn's version for kelvin.c

int main (int argc, char *argv[])
{
  kelvinInit(argc, argv);

  if (modelOptions->integration) {
    integrateMain();
  } else {
    iterateMain();
  }

  kelvinTerm();

  return EXIT_SUCCESS;
}

