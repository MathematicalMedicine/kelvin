/**********************************************************************
 * NICE System Daemon Query Controller
 * Alberto Maria Segre
 *
 * Copyright (c) 1997-2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software 
 * for non-profit educational purposes only.
 **********************************************************************/
#include "../config.h"
#include "niceaux.h"
#include "nicecom.h"
#include <stdio.h>		/* Needed for printf */
#include <stdlib.h>		/* EXIT_SUCCESS, EXIT_ERROR */
#include <string.h>		/* Needed for strncpy */
#include <unistd.h>		/* Needed for close */

#ifndef NOALARM
/* Used to timeout a query protocol. Note that this is a very short
 * alarm compared to the much longer HUNGING used internally in the
 * daemon and the API, simply because no user is going to want to sit
 * and wait forever. Also, timeouts here simply yield a "not
 * responding" message rather than doing something more fatal. */
#define HQRYINT 10		/* Hung threshold in sec. */
#endif

/* Port to query. */
unsigned short port = NICEPORT;	/* My port number. */

/* Formatting */
#define NICEQINDENT "   "	/* Indentation per generation of hierarchy. */
#define TARGET 0		/* Target host index. */

/* Internal prototypes. */
int niceQuery (int hid, int recursive, int *ncpus, int *navail);

/**********************************************************************
 * Query the NICE daemon running on specified machine. Can operate
 * recursively as well.
 **********************************************************************/
int
main (int argc, char *argv[])
{
  int i, host = FALSE;
  int recursive = FALSE;
  int ncpus = 0, navail = 0;
  int status;

  /* Next, check the command line arguments; there should be at most 1
   * of them, and anything else should cause an error. */
  for (i = 1; i < argc; i++)
    {
      if (strcmp (argv[i], "-?") == 0)
	{
	  printf ("Usage: %s [-n N] [-r] [-v] [host]\n", argv[0]);
	  printf ("where -? : this output\n");
	  printf ("      -n : specify port number\n");
	  printf ("      -r : query recursively\n");
	  printf ("      -v : version information.\n");
	  printf ("If host not specified, queries local daemon.\n");
	  exit (EXIT_FAILURE);
	}
      else if ((strcmp (argv[i], "-n") == 0) || (strcmp (argv[i], "--port") == 0))
	{
	  /* New port number. Use same limits as NICESAFEARGHD. */
	  i++;
	  sscanf (argv[i], "%5hd", &port);
	}
      else if (strcmp (argv[i], "-r") == 0)
	recursive = TRUE;	/* recursive specifies mode. */
      else if (strcmp (argv[i], "-v") == 0)
        printf ("NICE v%s\n", NICEVERSION);
      else if (!host)
	{
	  /* Copy specified host into host table. */
	  strncpy (nicehost, argv[i], NICEMAXSTRINGLEN);
	  host = TRUE;
	}
      else
	{
	  printf ("Usage: %s [-v] [-n N] [-r] [host]\n", argv[0]);
	  exit (EXIT_FAILURE);
	}
    }

  /* Check to see if a host was specified. If not, assume you're
   * interested in the local host. */
  if (!host)
    strncpy (nicehost, "localhost", NICEMAXSTRINGLEN);

  /* Get full host name from system. We're going to arbitrarily use
   * one of the host table indeces to store the target host's
   * parameters. We only need one such slot, so it doesn't much matter
   * which slot we use. */
  if ((qualifyHostname (nicehost, TARGET)) == ERROR)
    {
      fprintf (stderr, "Hostname error; aborting.\n");
      exit (EXIT_FAILURE);
    }

  /* Now go ahead and query the specified daemon. niceQuery can
   * return TRUE on success, FALSE if daemon doesn't respond, or ERROR
   * if recursive query failed once underway. */
  if ((status = niceQuery (TARGET, recursive, &ncpus, &navail)) == FALSE)
    {
      fprintf (stderr, "%s unresponsive.\n", nicehost);
      exit (EXIT_FAILURE);
    }
  else if (status == ERROR)
    {
      fprintf (stderr, "Error: query aborted.\n");
      exit (EXIT_FAILURE);
    }
  
  /* Print out stats. The number of available processors is calculated
   * as if you started a job on the given host, so we need to subtract
   * one from the number of available processors reported to account
   * for the root job. */
  if (recursive)
    fprintf (stderr, "%d/%d processors available.\n", navail, ncpus);

  /* Done. */
  exit (EXIT_SUCCESS);
}

/**********************************************************************
 * Print a status report after querying a daemon, optionally
 * recursively. niceQuery can return TRUE on success, FALSE if daemon
 * doesn't respond, or ERROR if recursive query failed once underway.
 **********************************************************************/
int
niceQuery (int hid, int recursive, int *ncpus, int *navail)
{
  int socket = ERROR, i, level;
  int cpus, avail, censor;
  char hname[NICEMAXSTRINGLEN + 1];	/* Host name */
  char hstat[NICEMAXSTRINGLEN + 1];	/* Host status */
  char message[NICEMAXSTRINGLEN + 1];
#ifndef NOALARM
  int queryAlarm;
#endif

#ifndef NOALARM
  /* Set a timer in case something blocks during the execution of the
   * protocol. A protocol is declared hung/blocked and attempts to
   * read or write are abandoned when it hangs for HQRYINT seconds. */
  if ((queryAlarm = niceSetTimeout (HQRYINT)) == FALSE)
    {
      /* The other end of the interaction has blocked. You shouldn't
       * need to clear the timer, which we know has already expired,
       * since that's how we got here. Note that here we return a
       * FALSE (rather than ERROR like in the API) to ensure the
       * proper "not responding" message is printed to the user. */
      close (socket);
      return (FALSE);
    }
#endif

  /* No daemon running on target machine, or target machine busy;
   * punt. We're not being very patient here; recall hailSocket may
   * return FALSE in the case of a temporary setback or ERROR when a
   * fatal condition is found. */
  if (hailSocket (hid, port, &socket) != TRUE)
    {
      if (socket != ERROR)
	close (socket);
      return (FALSE);
    }

  /* Send a status query. */
  if (sockprintf (socket, message, NICEQRYSTRDS, NQIDQUERY, NULLHOST) == ERROR)
    {
      /* Punt. */
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (queryAlarm);
#endif
      close (socket);
      return (ERROR);
    }

  /* Obtain information from designated host. */
  if (sockscanf (socket, message, NICEQRYSTRDDDSS, &level, &cpus, &avail,
		 hname, hstat) != 5)
    {
      /* Punt. */
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (queryAlarm);
#endif
      close (socket);
      return (ERROR);
    }

#ifndef NOALARM
  /* Clear the timer. Since its on such a short leash, and at this
   * point you're done waiting on your correspondent (unless the query
   * is recursive), clear the timer here. We'll reset the timer later
   * if indeed you are performing a recursive query. */
  niceClearAlarm (queryAlarm);
#endif

  /* If you're not recursive, you're done. Show output, tell the
   * server side you're done and then exit. */
  if (!recursive)
    {
      /* Print status: if the host name returned is "localhost" then
       * we should really "fix" the output to make it more
       * interpretable to the user. 
       *
       * Any recursively observed "localhost" designations will have
       * already been replaced by their host's IP address as they
       * percolate back up the hierachy.  Direct invocations of niceq
       * (either on current host or on some other named host) that
       * return "localhost" can instead simply return nicehost, since
       * that name was sufficient to get me to the appropriate niced
       * in the first place. */
      printf ("%s %s\n", 
	      ((strncmp (hname, "localhost", NICEMAXSTRINGLEN) == 0)?
	       nicehost:hname), hstat);

      /* No need to check status here, you're finished
       * anyway. Besides, you've already cleared the timer out, so
       * just blast away your ack and exit. */
      sockprintf (socket, message, NICEACKSTRD, FALSE);
      close (socket);
      return (TRUE);
    }

  /* Recursive query. Continue until you get the end-of-hierarchy mark
   * (level value will be ERROR; note level is initially 0, since
   * local host will be at level 0 by definition).
   *
   * We'll keep track of the number of cpus at each level and
   * accumulate the total in ncpus. To count the number of available
   * cpus in navail, we need to censor counts below any unavailable
   * node. We'll do this with the censor variable, which we'll reset
   * based on the level of the nodes as they fly by.
   *
   * One additional complication is the alarm. We can't know how many
   * machines are in the hierarchy a priori, so bounding everything
   * with an upper limit is probably a bad idea. Instead, set/clear
   * the alarm for each successive line of output, noting that each
   * line of output corresponds to a perhaps quite deep daemon/daemon
   * interaction (recall that the daemon side is on a longer timer as
   * well).
   * 
   * Crufty, but probably ok for queries. */
  censor = ERROR;
  while (level != ERROR)
    {
      /* Print current node status.
       *
       * If the host name returned is "localhost" then we should
       * really "fix" the output to make it more interpretable to the
       * user, but only if it is the top level hostname.
       *
       * Any recursively observed "localhost" designations will have
       * already been replaced by their host's IP address as they
       * percolate back up the hierachy.  Direct invocations of niceq
       * (either on current host or on some other named host) that
       * return "localhost" can instead simply return nicehost, since
       * that name was sufficient to get me to the appropriate niced
       * in the first place. */
      for (i = 0; i < level; i++) printf ("%s", NICEQINDENT); 
      printf ("%s %s\n", 
	      ((level == 0 && (strncmp (hname, "localhost", NICEMAXSTRINGLEN) == 0))?
	       nicehost:hname), hstat);

      /* Increment total number of processors. */
      *ncpus += cpus;

      /* Turn censoring on or off depending on availability and level
       * information. Censoring is turned on when you hit a node that
       * is either offline or fully populated with slaves, and turned
       * off once you pop out of that level. We'll use ERROR to
       * indicate no censoring, otherwise the level number at which
       * censoring takes effect. */
      if (censor != ERROR && level <= censor)
	censor = ERROR;
      if (censor == ERROR && !avail)
	censor = level;
      
      /* If you're not censored, increment number of available
       * processors. */
      if (censor == ERROR)
	*navail += cpus;

#ifndef NOALARM
      /* Reset the alarm. A protocol is declared hung/blocked and
       * attempts to read or write are abandoned when it hangs for
       * HQRYINT seconds. Returns from previous invocation of
       * niceSetTimeout () on expiration. */
      queryAlarm = niceResetTimeout (HQRYINT);
#endif
      /* Tell server side you want to go on. */
      if (sockprintf (socket, message, NICEACKSTRD, TRUE) == ERROR)
	{
	  /* Punt. */
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (queryAlarm);
#endif
	  close (socket);
	  return (ERROR);
	}
      
      /* Obtain information. */
      if (sockscanf (socket, message, NICEQRYSTRDDDSS, &level, &cpus, &avail,
		     hname, hstat) != 5)
	{
	  /* Punt. */
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (queryAlarm);
#endif
	  close (socket);
	  return (ERROR);
	}

#ifndef NOALARM
      /* Clear the timer before going on to the next line of
       * output. */
      niceClearAlarm (queryAlarm);
#endif
    }

  close (socket);
  return (TRUE);
}
