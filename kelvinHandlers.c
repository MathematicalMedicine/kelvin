#include "kelvin.h"
#include "kelvinHandlers.h"

volatile sig_atomic_t statusRequestSignal = FALSE;      ///< Status update requested via signal

/**

  Handler for SIGQUIT.

  Typing ^\ or issuing a kill -s SIGQUIT gets a dump of statistics.

  P.S. - cygwin requires "stty quit ^C" first for this to work.

*/
void quitSignalHandler (int signal)
{
  statusRequestSignal = TRUE;
#ifdef POLYSTATISTICS
  if (modelOptions.polynomial == TRUE)
    polyDynamicStatistics ("Signal received");
  else
#endif
    swDump (overallSW);
#ifdef DMTRACK
  swLogPeaks ("Signal");
#endif
}

#if defined (GPROF) || (GCOV)
/*

  Handler for SIGTERM.

  If we're profiling or doing coverage analysis, we catch a SIGTERM to 
  allow early exit(). An orderly exit like this (with EXIT_SUCCESS
  status, ensures that the profileing or coverage information completely
  written and fit for analysis.

*/
void termSignalHandler (int signal)
{
  fprintf (stderr, "Terminating early for gprof or gcov!\n");
  exit (EXIT_SUCCESS);
}
#endif

/// Handler for SIGINT
void intSignalHandler (int signal)
{
  fprintf (stderr, "Terminating early via interrupt!\n");
  exit (EXIT_FAILURE);
}

pid_t childPID = 0;     ///< For a child process producing timing (and memory?) statistics.
/**

  General-purpose exit handler

  Exit handler to clean-up after we hit any of our widely-distributed 
  exit points. Ensures that errant child processes are handled so
  we don't pester people with unnecessary messages.

*/
void exitKelvin ()
{
  swLogMsg ("Exiting");
  if (childPID != 0)
    kill (childPID, SIGKILL);   /* Sweep away any errant children */
}

void setupHandlers () {
  struct sigaction quitAction, intAction;
  sigset_t quitBlockMask, intBlockMask;

  /* Add an exit handler to deal with wayward children. */

  if (atexit (exitKelvin)) {
    perror ("Could not register exit handler!");
    exit (EXIT_FAILURE);
  }

  /** Setup signal handlers. Be sure to do this before starting
   any threads, or interrupts might not be handled properly. */

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
