#ifdef BACKTRACE
#include <execinfo.h> // For stack dump on demand
#endif
#include <ctype.h> // For isalnum()
#include "kelvin.h"
#include "kelvinHandlers.h"
#include "utils/sw.h"
#include "config/config.h" // For the modelOptions->polynomial flag
#include "utils/polynomial.h" // To respond to the polynomial flag

extern struct swStopwatch *overallSW;

volatile sig_atomic_t statusRequestSignal = FALSE;      ///< Status update requested via signal

void dumpStack (void) {
#ifdef BACKTRACE
  const size_t MAXSTACKDEPTH = 20;
  size_t stack_depth;
  void *stack_addrs[MAXSTACKDEPTH];
  char **stack_strings;

  stack_depth = backtrace(stack_addrs, MAXSTACKDEPTH);
  stack_strings = backtrace_symbols(stack_addrs, stack_depth);

  fprintf (stderr, "\nCall stack: ");
  size_t i;
  for (i = 2; i < stack_depth; i++) {
    char *wordBoundary;
    if ((wordBoundary = strchr( stack_strings[i], '+')) != '\0') {
      while ((isalnum(*wordBoundary) == 0) && (*wordBoundary != '_'))
	*(wordBoundary--) = '\0';
      while (isalnum(*wordBoundary) || (*wordBoundary == '_'))
	wordBoundary--;
      wordBoundary++;
      fprintf(stderr, "%s<-", wordBoundary);
      if (strcmp(wordBoundary, "main") == 0)
	break;
    }
  }
  fprintf (stderr, "\n");
  free(stack_strings); // malloc()ed by backtrace_symbols                                                             
#endif
  fflush(stderr);
}

/**

  Handler for SIGQUIT.

  Typing ^\ or issuing a kill -s SIGQUIT gets a dump of statistics.

  P.S. - cygwin requires "stty quit ^C" first for this to work.

*/
void quitSignalHandler (int ourSignal)
{
  dumpStack ();
  statusRequestSignal = TRUE;
#ifdef POLYSTATISTICS
  if (modelOptions->polynomial == TRUE)
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
void termSignalHandler (int ourSignal)
{
  fprintf (stderr, "Terminating early for gprof or gcov!\n");
  exit (EXIT_SUCCESS);
}
#endif

/// Handler for SIGINT
void intSignalHandler (int ourSignal)
{
  fprintf (stderr, "Terminating early via interrupt!\n");
  exit (EXIT_FAILURE);
}

/// Handler for SIGUSR1
void usr1SignalHandler (int ourSignal)
{
  fprintf (stderr, "This interrupt used to implement self-setting gdb breakpoints\n");
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
