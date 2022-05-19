/* Copyright (C) 2008, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */
#include <ctype.h> // For isalnum()

#ifdef PTMALLOC3
#include "utils/malloc-2.8.3.h"
#endif

#include "kelvin.h"
#include "kelvinHandlers.h"
#include "utils/sw.h"
#include "config/config.h" // For the modelOptions->polynomial flag
#include "utils/polynomial.h" // To respond to the polynomial flag

extern struct swStopwatch *overallSW;

/**

  Handler for SIGQUIT.

  Typing ^\ or issuing a kill -s SIGQUIT gets a dump of statistics.

  P.S. - cygwin requires "stty quit ^C" first for this to work.

*/
void quitSignalHandler (int ourSignal)
{
#ifdef BACKTRACE
  swDumpStack ();
#endif
#ifdef POLYSTATISTICS
  if (modelOptions->polynomial == TRUE)
    polyDynamicStatistics ("Signal received");
  else
#endif
#ifndef DISTRIBUTION
    swDump (overallSW);
#endif
#ifdef DMTRACK
  swLogPeaks ("Signal");
#endif
#ifdef PTMALLOC3
  malloc_stats ();
#endif
  swLogTimedProgress ();
}

#if defined (GPROF) || (GCOV)
/*

  Handler for SIGTERM.

  If we're profiling or doing coverage analysis, we catch a SIGTERM to 
  allow early exit(). An orderly exit like this (with EXIT_SUCCESS
  status, ensures that the profileing or coverage information completely
  written and fit for analysis.

*/
void termSignalHandler (int ourSignal)
{
  INFO ("Terminating early for gprof or gcov!");
  exit (EXIT_SUCCESS);
}
#else
/// Handler for SIGTERM. is same as SIGINT
void termSignalHandler (int ourSignal)
{
  FATAL ("Terminating early via interrupt");
}
#endif

/// Handler for SIGINT
void intSignalHandler (int ourSignal)
{
  FATAL ("Terminating early via interrupt");
}

/// Handler for SIGUSR1
void usr1SignalHandler (int ourSignal)
{
  INFO ("This interrupt used to implement self-setting gdb breakpoints");
}

/**

  General-purpose exit handler

  Exit handler to clean-up after we hit any of our widely-distributed 
  exit points. Ensures that errant child processes are handled so
  we don't pester people with unnecessary messages.

*/
void exitKelvin ()
{
  INFO ("Cleaning-up for exit");
  swDiagTerm ();
}

void setupHandlers () {
  struct sigaction quitAction, intAction, usr1Action;
  sigset_t quitBlockMask, intBlockMask, usr1BlockMask;

  /* Add an exit handler to deal with wayward children. */

  if (atexit (exitKelvin)) {
    perror ("Could not register exit handler!");
    exit (EXIT_FAILURE);
  }

  /** Setup signal handlers. Be sure to do this before starting
   any threads, or interrupts might not be handled properly. */


  sigfillset (&usr1BlockMask);
  usr1Action.sa_handler = usr1SignalHandler;
  usr1Action.sa_mask = usr1BlockMask;
  usr1Action.sa_flags = 0;
  sigaction (SIGUSR1, &usr1Action, NULL);

#if defined (GPROF) || (GCOV)
  struct sigaction termAction;
  sigset_t termBlockMask;
#endif
  sigfillset (&quitBlockMask);
  quitAction.sa_handler = quitSignalHandler;
  quitAction.sa_mask = quitBlockMask;
  quitAction.sa_flags = 0;
  sigaction (SIGQUIT, &quitAction, NULL);

  sigfillset (&intBlockMask);
  intAction.sa_handler = intSignalHandler;
  intAction.sa_mask = intBlockMask;
  intAction.sa_flags = 0;
  sigaction (SIGINT, &intAction, NULL);

#if defined (GPROF) || (GCOV)
  sigfillset (&termBlockMask);
  termAction.sa_handler = termSignalHandler;
  termAction.sa_mask = termBlockMask;
  termAction.sa_flags = 0;
  sigaction (SIGTERM, &termAction, NULL);
#endif

}
