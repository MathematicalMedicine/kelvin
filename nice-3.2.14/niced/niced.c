/**********************************************************************
 * NICE System Daemon
 * Alberto Maria Segre
 *
 * Copyright (c) 1997-2005, The University of Iowa.  All rights
 * reserved.  Permission is hereby given to use and reproduce this
 * software for non-profit educational purposes only.
 **********************************************************************/
#include "niced.h"

/**********************************************************************
 * A bunch of global variables. These need to be global so they can be
 * accessed by the signal handler that traps kill or interrupt
 * signals; but the bottom line is that it is simpler to make these
 * global so we don't need to pass them around all over the place.
 **********************************************************************/
int mode = ACTIVE;		/* Daemon mode; one of ACTIVE or ORPHAN. */
int maxChildren = MAXCHILDREN;
int cardinality = 1;
int maxSlaves = MAXSLAVES;
int utcf = 0;			/* Universal time correction (minutes). */
int debug = FALSE;		/* Governs daemonization and output behavior. */
unsigned short int port = NICEPORT;	/* My port number. */
int priority = PRIORITYDFLT;	/* Daemon scheduling priority. */

/* The following variables are indexed by hid, like niceaddr[] from
 * nicecom.h. All values in nicestart[], niceend[] and nicelast[] are
 * expressed in local time, but are communicated in universal time. */
int niceload[MAXHOSTS];		/* Load on known daemons. */
int nicewall[MAXHOSTS];		/* Segregated subdomain. */
int nicedepth[MAXHOSTS];	/* Depth of highest opening in subtree. */
int nicestart[MAXHOSTS];	/* Scheduling information. */
int niceend[MAXHOSTS];
time_t nicelast[MAXHOSTS];	/* Last contact time. */
int nicebfact = MAXCHILDREN;	/* Hierarchy's global branching factor. */

/* Variables describing running slave applications. */
int slaves = 0;			/* Number of local slave processes. */
int spids[MAXSLAVES];		/* Slave PIDs. */

/**********************************************************************
 * This is the NICE system daemon. If a copy of the daemon is already
 * running on the current machine, the program exits and returns an
 * error.  Otherwise, it reports to the host specified in the
 * configuration file, and then monitors the NICE port for resource
 * requests, information about system status, or control signals.
 **********************************************************************/
int
main (int argc, char *argv[])
{
  int i, status;
  int lastload;			/* Last time load was measured. */
  int avail = TRUE;		/* Current availability. */
  char parent[NICEMAXSTRINGLEN + 1];	/* Parent at configuration time. */
  int socket;			/* Socket requiring service. */
  int qid;			/* Query identifier. */

  struct in_addr caddr;		/* Client IP address (from network). */
  char message[NICEMAXSTRINGLEN + 1];

  time_t now_t;			/* Current time: time_t format. */
  struct tm *now_tm;		/* Current time: tm format. */
  int now = 0;			/* Current time: mins after midnight format. */
#ifndef NOALARM
  int serverAlarm;
#endif
  struct sigaction sigact;	/* Signals. */

  /* Initialize host table. */
  for (i = 0; i < MAXHOSTS; i++)
    {
      /* 09/26/2006 */
      /* Moved to qualifyHostname() because of explicit malloc() due
       * to gcc v4.0 switch. */

      //inet_aton (NULLHOST, &niceaddr[i]);
      niceload[i] = nicedepth[i] = nicestart[i] = niceend[i] = 0;
      nicewall[i] = FALSE;
    }

  /* Initially, the parent specified in the configuration file will be
   * stored in both the ROOT and PARENT host table slots. When
   * registration is successfully completed, the PARENT host table
   * entry will be overwritten with the actual parent host
   * information. This is because the registration process may require
   * changing parent designation. */
  strncpy (parent, NULLHOST, NICEMAXSTRINGLEN);

  /* Next, obtain the nice configuration for the daemon. Configuration
   * information is contained on the command line as well as in the
   * configuration file, whose location may be specified on the
   * command line. Returns TRUE if configuration is successful.  
   *
   * Most configuration values are read into global variables, and so
   * need not be passed. This includes maxChildren, maxSlaves,
   * cardinality, parent, port, debug, nicewall, nicestart, and
   * niceend.
   *
   * The only local configuration variables that are set are the
   * parent name and the scheduling priority. */
  if (!niceConfig (argc, argv, parent, &priority))
    {
      printf ("Usage: %s [-c file] [-n port] [-p parent] [-b bfactor] [-d] [-?]\n",
	 argv[0]);
      exit (EXIT_FAILURE);
    }

  /* Trap any signals in order to make sure any necessary cleanup is
   * performed. Probably don't need to bother to trap all the other
   * possible signals; unlikely we need to bother cleaning up in the
   * event of a power or hardware failure! Recall SIGKILL and SIGSTOP
   * behavior can't be changed. TODO: probably should thoroughly check
   * all signals to see if this set is reasonable. */
  sigact.sa_handler = dieNeatly;
  sigemptyset (&sigact.sa_mask);
  sigact.sa_flags = 0;
  sigaction (SIGHUP, &sigact, NULL);
  sigaction (SIGINT, &sigact, NULL);
  sigaction (SIGQUIT, &sigact, NULL);
  sigaction (SIGTERM, &sigact, NULL);

  /* Get host name from system. On most Unix systems, this will be the
   * alphanumeric name of the machine. On DHCP configured nodes, this
   * will probably just be the dotted octet version of the IP
   * address. If the system is not networked, this might just be some
   * arbitrary string, or even the special string "localhost." */
  if ((gethostname (nicehost, NICEMAXSTRINGLEN)) == ERROR)
    niceLog (NLOGFAT, "Hostname error; aborting.\n");

  /* Make sure we have the fully-specified name (including domain) and
   * not just a short name as is returned by gethostname on some Unix
   * systems. qualifyHostname also stashes the host's IP address
   * (network byte order) in niceaddr[SELF].
   *
   * Note that if the machine does not really have a "real" hostname,
   * but just some arbitrary string, when you qualify the hostname it
   * might change into "localhost." (IS THIS TRUE?)
   *
   * Realize also that the address stashed in niceaddr[SELF] may not
   * be valid anywhere else on the network; it may indeed be
   * 127.0.0.1. See qualifyHostname() in nicecom.c for details. */
  if (qualifyHostname (nicehost, SELF) == ERROR)
    niceLog (NLOGFAT, "Hostname error; aborting.\n");

  /* Fix root host specification. We'll start by setting parent name
   * string and corresponding IP address in niceaddr[ROOT]. Later, the
   * registration process will establish niceaddr[PARENT]. */
  if (strncmp (parent, NULLHOST, NICEMAXSTRINGLEN) == 0)
    niceLog (NLOGFAT, "No parent specified; aborting.\n");
  else if (strncmp (parent, "localhost", NICEMAXSTRINGLEN) == 0)
    {
      /* Special parent name "localhost" was given. This is a synonym
       * for the current host, often used in a standalone
       * environment. We'll copy nicehost into the parent string and
       * its previously qualified address into the ROOT address
       * location. */
      strncpy (parent, nicehost, NICEMAXSTRINGLEN);
      memcpy (&niceaddr[ROOT], &niceaddr[SELF], ADDRBYTES);
    }
  else if (qualifyHostname (parent, ROOT) == ERROR)
    niceLog (NLOGFAT, "Unknown parent %s; aborting.\n", parent);

  /* Finally, open the NICE daemon port for monitoring. If you can't
   * grab the port, it means there is probably already a NICE daemon
   * running, so punt. */
  if ((nicesock = watchSocket (&port, port)) == ERROR)
    niceLog (NLOGFAT, "Port %hd unavailable; aborting.\n", port);

  /* OK, we're in this for the long haul. Unless you've said otherwise
   * on the command line, turn this process into a proper Unix
   * daemon. */
  if (debug == FALSE)
    daemonize ();

  /* Get a load estimate. Good to do this before you register, since
   * this way your parent has something to go on if you need to be
   * repositioned. The load estimate is based on a benchmark function
   * that is run at regular intervals, and weighted by process
   * priority. */
  niceload[SELF] = bench ();

  /* Unless you're the root node, report to your parent whose address
   * is already stored in the ROOT host table entry. You know you're
   * the root when you're your own parent, so check your address
   * against the specified parent address, which is already stored in
   * the ROOT host table location. */
  if (ADDRNEQ (niceaddr[SELF], niceaddr[ROOT]))
    {
      /* You're not the root, and you already have the desired
       * parent's FQDN, etc, stored in the host table ROOT location
       * (from qualifyHostname above). Copy it to the PARENT slot 
       * (if you are the root, then the parent slot must remain null). */
      memcpy (&niceaddr[PARENT], &niceaddr[ROOT], ADDRBYTES);

      /* Now try to register with the designated root node, now stored
       * in niceaddr[PARENT] location. An ERROR here means the root is
       * refusing you on grounds of failed authentication, in which
       * case you should punt. A FALSE here means that you either
       * failed to make contact (network problems?) or couldn't be
       * accomodated. Put yourself in orphan mode so that periodic
       * attempts are made to reregister with both the originally
       * designated parent (still in ROOT location) and whatever
       * alternate parent the registration process may have settled on
       * before aborting. */
      if ((status = niceRegister (PARENT)) == ERROR)
	exit (EXIT_FAILURE);
      else if (status == FALSE)
	mode = ORPHAN;
    }
  else
    {
      /* You're the root, so establish your maxChildren value
       * (established in niceConfig) as the hierarchy's target
       * branching factor. */
      nicebfact = maxChildren;

      /* TODO: Should we set niceaddr[SELF]=niceaddr[ROOT]=localaddr??
       *
       * Recall you recognize yourself as the root when
       * niceaddr[SELF]=niceaddr[ROOT]. Need to check the content of
       * niceaddr[] everywhere to ensure you never use the (known
       * unreliable) value in niceaddr[SELF]. (Is this note still
       * valid? has to do with use of localhost as the root
       * designation?) */
    }

  /* All systems are go at this point. Ensures load will get measured
   * again as soon as we enter the loop (now); note that we're on a
   * 24-hour clock here, so we have to wrap around appropriately. */
  lastload = nicestart[SELF] - (LOADINT + 1);
  lastload = (lastload < 0) ? lastload + MAXTIME : lastload;

  /* Set the seed on the random number generator to the current
   * process ID. */
  srandom (getpid ());

  /* Initialize the application downloading wish lists. */
  initWishes ();

  /* Be a good Unix citizen: lower the scheduling priority for the
   * daemon (i.e., "nice" the nice daemon!). */
  if (nice (priority) == ERROR)
    niceLog (NLOGFAT, "Priority change failed; aborting.\n");

  /* Tell the world you're up successfully. */
  niceLog (NLOGMSG, "Initialized%s\n",
	   (ADDRNULL (niceaddr[PARENT]) ? " (root)."
	    : ((mode == ORPHAN) ? " (orphan)." : ".")));

  /* Now start the monitoring loop. */
  while (TRUE)
    {
      /* Fetch current time as a time_t value. We could use
       * gettimeofday here (instead of time), followed by strftime
       * instead of ctime. While this would be simpler, it precludes
       * extending the configuration file to control availability
       * during different days of the week, for example. */
      time (&now_t);

      /* Convert time to an integer representing "minutes after
       * midnight". We can do this with the TIME macro, but first we
       * convert the time_t value into a broken-down time (tm)
       * structure, so we can easily extract minutes and hours. Also
       * reset the universal time correction factor using the timezone
       * variable, set as a side effect of the call to localtime. Note
       * that utcf may well change with daylight savings time, for
       * example. */
      now_tm = localtime (&now_t);
#ifndef __CYGWIN__
      utcf = (timezone / 60);
#else
      utcf = (_timezone / 60);
#endif
      now = TIME ((*now_tm).tm_hour, (*now_tm).tm_min);

      /* If you are in ORPHAN mode, then you need to attempt to
       * reconnect to the NICE hierarcy. */
      if (mode == ORPHAN)
	{
	  if ((status = niceRegister (PARENT)) == TRUE)
	    /* First, try to register with your last known parent. */
	    mode = ACTIVE;
	  else if (status == FALSE &&
		   ADDRNEQ (niceaddr[PARENT], niceaddr[ROOT]) &&
		   (status = niceRegister (ROOT)) == TRUE)
	    /* If that doesn't work, try the root, but only if its
	     * different than your last known parent. */
	    mode = ACTIVE;
	  else if (status == ERROR)
	    /* You got an authentication error on one of your
	     * attempts. Punt. */
	    exit (EXIT_FAILURE);
	}

      /* Check host availability for this iteration. First, we need to
       * check to see if there is a change in host availability. If
       * there is, then run the benchmark program to get an up to date
       * load figure, and report new load and availability to child
       * daemons. Otherwise, if the benchmark program hasn't been run
       * lately, run it and report the new load and current
       * availability to child daemons. */
      if ((avail == TRUE) &&
	  (((nicestart[SELF] > niceend[SELF]) &&
	    (now < nicestart[SELF]) && (now >= niceend[SELF])) ||
	   ((nicestart[SELF] < niceend[SELF]) &&
	    ((now < nicestart[SELF]) || (now >= niceend[SELF])))))
	{
	  /* Take host out of service; report change in status to
	   * child daemons. */
	  avail = FALSE;
	  niceload[SELF] = bench ();
	  lastload = now;
	  niceLog (NLOGMSG, "Going offline (%dx%d).\n",
		   cardinality, niceload[SELF]);

	  /* Force your slaves to sleep. */
	  pauseSlaves ();
	  /* Report status to children (if any). */
	  niceLoad ();
	}
      else if ((avail == FALSE) &&
	       (((nicestart[SELF] > niceend[SELF]) &&
		 ((now >= nicestart[SELF]) || (now < niceend[SELF]))) ||
		((nicestart[SELF] < niceend[SELF]) &&
		 (now >= nicestart[SELF]) && (now < niceend[SELF]))))
	{
	  /* Return host to service; report change in status to child
	   * daemons. */
	  avail = TRUE;
	  niceload[SELF] = bench ();
	  lastload = now;
	  niceLog (NLOGMSG, "Coming online (%dx%d).\n",
		   cardinality, niceload[SELF]);

	  /* Wake your slaves up. */
	  prodSlaves ();
	  /* Report status to children (if any). */
	  niceLoad ();
	}
      else if (((now >= lastload) && (lastload + LOADINT <= now)) ||
	       ((now < lastload) && (lastload + LOADINT - MAXTIME <= now)))
	{
	  /* No change in availability, but it is time to measure load
	   * and solicit information from your child daemons. This
	   * should happen regardless of whether or not you are
	   * online. */
	  niceload[SELF] = bench ();
	  lastload = now;

	  /* If you're online, your slaves should already be
	   * running. Yet it doesn't hurt to try waking them anyway,
	   * since someone may have forced them to sleep from the
	   * shell. Also, prodSlaves discovers any lapsed processes
	   * and update the daemon's internal state appropriately. */
	  if (avail == TRUE)
	    prodSlaves ();

	  /* Report status to children (if any). Do this even if we
	   * are offline; being offline only means you won't run any slave
	   * processes. */
	  niceLoad ();
	}

      /* Check the NICE port for incoming messages and act on them. We
       * do this always, even if the host isn't available, although
       * some messages may be handled differently when this is the
       * case. The SLEEPINT interval defaults to 5 seconds, and
       * regulates how often we do any bookkeeping chores when no
       * incoming messages are present. This value should be kept well
       * below the granularity of the config file's time
       * specifications. */
      while ((socket = checkSocket (nicesock, &caddr, SLEEPINT)) >= 0)
	{
	  /* Only executed if message received. Socket timeout would
	   * return ERROR from checkSocket, skipping this inner
	   * loop. */
#ifndef NOALARM
	  /* Set a timer in case something blocks during the execution
	   * of the protocol. A protocol is declared hung/blocked and
	   * attempts to read or write are abandoned when it hangs for
	   * HUNGINT seconds, typically a large number (as in the
	   * niceapi.c protocol timeout). Since our server is
	   * single-threaded, it is susceptible to denial-of-service
	   * attacks. So we first set a shorter alarm of AUTHINT
	   * seconds to authenticate, and then reset the alarm for
	   * HUNGINT seconds to let the protocol complete. */
	  if ((serverAlarm = niceSetTimeout (AUTHINT)) == FALSE)
	    {
	      /* The other end of the interaction has blocked. You
	       * shouldn't need to clear the timer, which we know has
	       * already expired, since that's how we got here. */
	      close (socket);
	      continue;
	    }
#endif
	  /* OK, obtain the QID for the incoming service message. We
	   * will need to decide how to dispatch according to
	   * QID. QIDs might come from another daemon, an application,
	   * or a control program. */
	  if (sockscanf (socket, message, NICEQRYSTRD, &qid) != 1)
	    {
#ifndef NOALARM
	      /* Clear the alarm. */
	      niceClearAlarm (serverAlarm);
#endif
	      /* Close the socket and go on to the next incoming message. */
	      close (socket);
	    }
	  else
	    {
	      /* Messages divide into roughly three groups: (1) those
	       * that come from known other damons, (2) those that
	       * come from applications on the local host, and (3)
	       * those that can come from anywhere. In the first group
	       * are most of the daemon-to-daemon protocols, except
	       * for NQIDREGISTER (since we can't know who that's
	       * coming from a priori), although it is possible to get
	       * a NQIDREGISTER from a previously known child that has
	       * just rebooted. In the second group are all the
	       * application protocols, and in the third group are
	       * things like niceq queries.
	       *
	       * First, try to identify the sender by index number if
	       * you can. This is cheap. */
	      i = 0;
	      while ((i <= nicehcnt) && ADDRNEQ (niceaddr[i], caddr))
		i++;

	      /* Next, handle messages that should only come from
	       * known correspondents but don't. */
	      if ((i > nicehcnt) && (qid == NQIDUNREGISTER))
		{
		  /* Recall NQIDUNREGISTER is asynchronous or
		   * nonblocking, so we needn't respond: just issue a
		   * warning message and close the socket. */
		  niceLog (NLOGWRN,
			   "Ignoring NQIDUNREGISTER from unknown host %s.\n",
			   inet_ntoa (caddr));
#ifndef NOALARM
		  /* Clear the alarm. */
		  niceClearAlarm (serverAlarm);
#endif
		  /* Close the socket and go on to the next incoming message. */
		  close (socket);
		}
	      else if ((i > nicehcnt)
		       && ((qid == NQIDPROMOTE) || (qid == NQIDORPHAN)
			   || (qid == NQIDLOAD) || (qid == NQIDSPAWN)
			   || (qid == NQIDTXFR)))
		{
		  /* Tell the supplicant to buzz off, then close the
		   * socket.  Since it is contacting you, and since
		   * these protocols are all parent initiated, its
		   * safe to assume this daemon probably thinks it is
		   * your parent. Either you've already dropped this
		   * parent and found another, or perhaps you rebooted
		   * and registered elsewhere. In any case, its safest
		   * to simply rebuff the connection and go about your
		   * business. Note that sending an ERROR here would
		   * allow the pseudoparent to continue thinking you
		   * were its child, albeit an unresponsive one: FALSE
		   * here is a rejection that forces the pseudoparent
		   * to drop you immediately. */
		  sockprintf (socket, message, NICEACKSTRD, FALSE);
#ifndef NOALARM
		  /* Clear the alarm. */
		  niceClearAlarm (serverAlarm);
#endif
		  /* Close the socket and go on to the next incoming message. */
		  close (socket);

		  /* Issue appropriate warning message. */
		  niceLog (NLOGWRN, "Ignoring %s from unknown host %s.\n",
			   ((qid == NQIDPROMOTE) ? "NQIDPROMOTE" :
			    (qid == NQIDORPHAN) ? "NQIDORPHAN" :
			    (qid == NQIDLOAD) ? "NQIDLOAD" :
			    (qid == NQIDSPAWN) ? "NQIDSPAWN" : "NQIDTXFR"),
			   inet_ntoa (caddr));
		}
	      else
		{
		  /* OK, all systems go. If this is a recognized host,
		   * then update their last contact time to the
		   * current (local) time. Then go ahead and dispatch
		   * on QID, being careful to close the appropriate
		   * socket within each handler. */
		  time (&nicelast[i]);

#ifndef NOALARM
		  /* You've properly authenticated, so clear the short
		   * AUTHINT alarm and reset it to the longer HUNGINT
		   * in order to allow enough time for the protocol
		   * handler to complete. */
		  niceClearAlarm (serverAlarm);
		  serverAlarm = niceResetTimeout (HUNGINT);
#endif
		  switch (qid)
		    {
		    case NQIDREGISTER:	/* daemon only */
		      serveNiceRegister (socket, caddr);
		      break;

		    case NQIDUNREGISTER:	/* daemon only, nonblocking */
		      serveNiceUnregister (socket, caddr);
		      break;
#ifndef NODYNAMIC
		    case NQIDPROMOTE:	/* daemon only */
		      serveNicePromote (socket, caddr);
		      break;
#endif
		    case NQIDORPHAN:	/* daemon only */
		      serveNiceOrphan (socket, caddr);
		      break;

		    case NQIDLOAD:	/* daemon only */
		      serveNiceLoad (socket, caddr);
		      break;

		    case NQIDSPAWN:	/* daemon only */
		      serveNiceSpawn (socket, caddr, avail);
		      break;

		    case NQIDTXFR:	/* daemon only */
#ifndef NOALARM
		      serveNiceTxfr (socket, caddr, avail, serverAlarm);
#else
		      serveNiceTxfr (socket, caddr, avail);
#endif
		      break;

		    case NQIDQUERY:	/* niceq */
		      serveNiceQuery (socket, avail);
		      break;

		    case NQIDPING:	/* niceapi */
		      serveNicePing (socket, avail);
		      break;

		    case NQIDSOLICIT:	/* niceapi */
		      serveNiceSolicit (socket, avail);
		      break;

		    case NQIDEXIT:	/* niceapi */
		      serveNiceExit (socket);
		      break;

		    default:	/* Unrecognized request. */
		      {
			/* Close the socket and print a warning message. */
			close (socket);
			niceLog (NLOGWRN,
				 "Received unknown QID %d from %s.\n", qid,
				 inet_ntoa (caddr));
		      }
		    }
#ifndef NOALARM
		  /* Finished. Clear the reset alarm before going on
		   * to the next incoming message. */
		  niceClearAlarm (serverAlarm);
#endif
		}
	    }
	}

      /* Check here for any unresponsive machines. This is where we
       * have to be careful if we wish to avoid a deadlock
       * condition. If the child gives up on its parent before the
       * parent gives up on the child, it is possible that the child
       * will attempt to reregister with that parent at the same time
       * the parent is trying to solicit a load update from the
       * child. The registration protocol is the only protocol
       * initiated by a child, and the liveness proof relies on the
       * notion that all protocols in this single-threaded server
       * model are initiated by the parent; so we must be absolutely
       * sure that there is no chance the parent could already know
       * about that child before the child tries to register. */
      time (&now_t);
      if (mode != ORPHAN && !ADDRNULL (niceaddr[PARENT]) &&
	  (difftime (now_t, nicelast[PARENT]) > PDRPINT))
	{
	  /* We have not heard from our parent recently. */
	  if (debug &&
	      (errhp =
	       gethostbyaddr ((char *) &niceaddr[PARENT], ADDRBYTES,
			      AF_INET)))
	    niceLog (NLOGWRN, "Orphaned by %s.\n", errhp->h_name);
	  else
	    niceLog (NLOGWRN, "Orphaned by %s.\n",
		     inet_ntoa (niceaddr[PARENT]));

	  /* Go into orphan mode if either my parent was the root or
	   * the root is also unreachable. Exit on authentication
	   * failure. */
	  if (ADDREQ (niceaddr[ROOT], niceaddr[PARENT]))
	    mode = ORPHAN;
	  else if ((status = niceRegister (ROOT)) == FALSE)
	    mode = ORPHAN;
	  else if (status == ERROR)
	    exit (EXIT_FAILURE);
	}

      /* Check children to see if any have lapsed. If they haven't
       * made contact recently, drop them. */
      i = 1;
      while (i <= nicehcnt)
	{
	  if (difftime (now_t, nicelast[i]) > CDRPINT)
	    dropChild (i);
	  else
	    i++;
	}

      /* Reap any stray zombie applications; this can happen due to
       * race conditions when applications exit. Reaped pids are
       * culled from spids[], which is then compacted. Normally,
       * serveNiceExit should reap the exiting process: reaping again
       * once each cycle is a low-cost failsafe measure. */
      i = 0;
      while (i < slaves)
	if (waitpid (spids[i], NULL, WNOHANG) != spids[i])
	  i++;
	else
	  spids[i] = spids[--slaves];

      /* Now cull any completed transfer processes from the wish
       * lists, and count how many active transfers are ongoing. Use
       * this to compute how many new transfers you can start. */
      i = MAXTXFRS - cullWishes ();

      /* Now start up i new transfer processes, but only if you are in
       * fact avail. Since you can only start source/put processes
       * (since all transfers are parent-initiated), we'll reserve
       * some available process for receiving an application should
       * your own parent suddenly elect to send you one. */
      if (i > 0 && avail == TRUE)
	{
	  i = i - countPendingGets ();
	  startPendingPuts (MAX (1, i));
	}
    }
}

/**********************************************************************
 * Call this procedure to turn yourself into a daemon. See p418 of
 * Steven's book "Advanced Programming in the UNIX Environment" for
 * more details.
 **********************************************************************/
void
daemonize ()
{
  int pid;

  /* Fork a process and have parent exit. Returns control to the shell
   * if this was called from the shell, and also makes sure the daemon
   * is not a process group leader. */
  if ((pid = fork ()) < 0)
    {
      /* Oops! Can't fork. Abort. */
      exit (EXIT_FAILURE);
    }
  else if (pid != 0)
    {
      /* Let the parent die. */
      exit (EXIT_SUCCESS);
    }
  else
    {
      /* Create a new session and become the session leader. */
      setsid ();

#if FALSE
      /* Make sure you are running in the NICE directory. */
      chdir (NBINDIR);
#endif

      /* Clear file mode creation mask. This allows daemon to open
       * files without any a priori constraints on file access. */
      umask (0);

      /* TODO: should probably also close stderr and stdout, but they
       * are currently useful for application debugging. */
      fclose (stdin);

      /* OK, we're done. */
      return;
    }
}

/**********************************************************************
 * dieNeatly is called when SIGINT or SIGKILL are received. It simply
 * invokes die with arguments that cause die to restructure the daemon
 * hierarchy.
 **********************************************************************/
void
dieNeatly (int signum)
{
  die (TRUE, TRUE);
}

/**********************************************************************
 * die is called by both dieNeatly and dieQuickly. It (optionally)
 * repairs the daemon hierarchy, then kills all local slave processes
 * and closes the nicesocket.
 **********************************************************************/
void
die (int restructure, int notifyParent)
{
  int i;

  /* Unregister this daemon unless it is the root daemon (has no
   * parent). */
  if (!ADDRNULL (niceaddr[PARENT]) && notifyParent && mode != ORPHAN)
    niceUnregister ();

  if (restructure && (nicehcnt > 0))
    {
      /* Tell your children they're orphans; trust them to do the
       * right thing and reattach themselves either with you (when you
       * come back up) or with their original root. Note: niceOrphan()
       * always succeeds, even if it doesn't manage to connect. This
       * shouldn't matter as I'm about to die anyway. */
      i = nicehcnt;
      while (i > 0)
	{
	  niceOrphan (i);
	  i--;
	}
    }

  /* Disable communication by closing the nicesock. This keeps your
   * slaves from trying to contact you back and tell you they are
   * exiting.  Much the same might be accomplished by sending SIGKILL
   * instead, but we don't want to do that since we do want the slave
   * to trap the signal and notify its own naggers it is exiting, even
   * if it can't get back to me. See niceExit for more details. */
  close (nicesock);

  /* Log output. */
  if (slaves > 0)
    niceLog (NLOGMSG, "Killing %d slaves.\n", slaves);

  /* Kill your running slave applications by sending SIGQUIT. */
  for (i = 0; i < slaves; i++)
    kill (spids[i], SIGQUIT);

  /* Kill your running file transfer applications by sending
   * SIGKILL. We can kill these pretty abruptly, and allow the other
   * end of the transfer to recognize that the process has
   * evaporated. Of course, we'll have to recognize the files on disk
   * are corrupt before we try to actually use them! */
  killActiveTransfers ();

  niceLog (NLOGMSG, "Exiting.\n");
  exit (EXIT_SUCCESS);
}

/**********************************************************************
 * Compute the overlap between two start/end time pairs.
 **********************************************************************/
int
computeOverlap (int start1, int end1, int start2, int end2)
{
  /* Gross, but effective. */
  if (start1 > end1)
    {
      /* Start late night, end early morning. */
      if (start2 > end2)
	return (MIN (end1, end2) + MAXTIME - MAX (start1, start2));
      else if (start2 < end2)
	return (MAX (end2 - start1, 0) + MAX (end1 - start2, 0));
      else
	return (MAXTIME - start1 + end1);
    }
  else if (start1 < end1)
    {
      /* Start early morning, end late night. */
      if (start2 > end2)
	return (MAX (end1 - start2, 0) + MAX (end2 - start1, 0));
      else if (start2 < end2)
	return (MIN (end1, end2) - MAX (start1, start2));
      else
	return (end1 - start1);
    }
  else
    {
      /* Always available; start1 = end1. */
      if (start2 > end2)
	return (end2 + (MAXTIME - start2));
      else if (start2 < end2)
	return (end2 - start2);
      else
	return (MAXTIME);
    }
}

/**********************************************************************
 * Find best child given a set of parameters and a definition of
 * "best." Defining best involves taking four factors into account:
 * service time overlaps, depth values, load factors, and whether or
 * not the child is a barrier node. While these are always considered
 * in the same priority order, you can define what makes a "good"
 * match for each by setting the Boolean variables highOvlp,
 * highDepth, highLoad, and noBarriers appropriately. If the defined
 * factors preclude finding anychild, returns ERROR.
 *
 * Used in serveNiceRegister() to decide where to forward a supplicant.
 * Used in serveNiceRegister() to decide who to demote.
 * Used in niceRegister() when changing branching factor.
 * Used in niceLoad() (via electChild()) to decide whom to promote.
 * Used in niceLoad() when changing branching factor.
 *********************************************************************/
int
bestMatch (int start, int end, int depth, int load, int noBarriers,
	   int highOvlp, int highDepth, int highLoad)
{
  int i, ovlp;
  int best = ERROR, bestovlp = 0, bestcnt = 0;
#if FALSE
  fprintf (stderr, "Computing [%d:%d,%d,%d]'s best (%s%s%s) match",
	   start, end, load, depth,
	   (highOvlp ? "t" : "f"), (highDepth ? "t" : "f"),
	   (highLoad ? "t" : "f"));
#endif
  /* Look at each child, and compare their overlap/load/latency with
   * the specified starttime/endtime. */
  i = 0;
  while (++i <= nicehcnt)
    {
      /* Skip barrier children if so requested. */
      if (noBarriers && nicewall[i])
	continue;

      ovlp = computeOverlap (start, end, nicestart[i], niceend[i]);
#if FALSE
      fprintf (stderr, " [%d,%d,%d]", ovlp, niceload[i], nicedepth[i]);
#endif
      /* Now check to see if the computed overlap is better than the
       * best overlap encountered so far, if any (if this is the first
       * you look at, then best will be ERROR and you will accept this
       * choice as the best so far). */
      if ((best == ERROR) ||
	  (highOvlp ? (ovlp > bestovlp) : (ovlp < bestovlp)) ||
	  (ovlp == bestovlp &&
	   ((highDepth ? (nicedepth[i] > nicedepth[best])
	     : (nicedepth[i] < nicedepth[best]))
	    || ((nicedepth[i] == nicedepth[best])
		&&
		((highLoad ? (niceload[i] > niceload[best])
		  : (niceload[i] < niceload[best])))))))
	{
	  best = i;
	  bestovlp = ovlp;
	  bestcnt = 1;
#if FALSE
	  fprintf (stderr, "*");
#endif
	}
      else if ((ovlp == bestovlp) &&
	       (nicedepth[i] == nicedepth[best]) &&
	       (niceload[i] == niceload[best]))
	{
#if FALSE
	  fprintf (stderr, "$");
#endif
	  bestcnt++;
	}
    }
  if (bestcnt > 1)
    {
      /* In the event that there was not a unique best child, break
       * ties at random. Start with the first equivalent child (best)
       * and scan up to nicehcnt, stopping when you get to the
       * randomly selected equivalent child. */
      bestcnt = RANDINT (bestcnt);
#if FALSE
      fprintf (stderr, " [%d-way tie]", bestcnt);
#endif
      i = best;
      while (bestcnt > 0 && ++i <= nicehcnt)
	{
	  if ((niceload[i] == niceload[best]) &&
	      (nicedepth[i] == nicedepth[best]) &&
	      (computeOverlap (start, end, nicestart[i], niceend[i]) ==
	       bestovlp))
	    {
	      best = i;
	      bestcnt--;
	    }
	}
    }
#if FALSE
  fprintf (stderr, " => %d\n", best);
#endif
  /* Return your best choice. */
  return (best);
}

#ifndef NODYNAMIC
/**********************************************************************
 * electChild() determines whether there is one of my children that
 * should be promoted. There are two general cases.
 *
 * If your parent has a higher vacancy, then the promotion is simply
 * an attempt to better fill out the tree. In this case, we would want
 * to promote the child with the greatest overlap, shallowest depth,
 * and lowest load. The depth constraint ensures that the hierarchy
 * doesn't change too fast; presumably a remaining deeper child would
 * then promote one of your grandchildren to fill in for your missing
 * child.
 *
 * In the second case, we are trying to see if there is a "better"
 * match between your child and your parent, with the goal of
 * replacing self in the registration process. Here, the current
 * heuristic finds the child of the current node that is the bestMatch
 * with the parent of the current node and only considers promoting
 * this child if it exceeds the appropriate promotion criteria.
 *
 * Never called by root or barrier node.
 *
 * Some day, we might also want to encourage random reconfiguration
 * with some small probability, to avoid local minima as best we can
 * (we might eventually do this by occasionally orphaning a child so
 * it reregisters with its originally specified ROOT).
 *
 * Used in niceLoad() to decide who to promote.
 **********************************************************************/
int
electChild ()
{
  int child, ovlp1, ovlp2;

  if (nicedepth[SELF] >= nicedepth[PARENT])
    {
      /* A sure promotion, because your parent has a higher level
       * vacancy. Find which of your children is most deserving. A
       * good match considers barrier nodes, has high overlap, is
       * shallow, and has low load. */
      child = bestMatch (nicestart[PARENT], niceend[PARENT],
			 nicedepth[PARENT], niceload[PARENT],
			 FALSE, TRUE, FALSE, FALSE);
      return (child);
    }
  else
    {
      /* A possible promotion, if your child can beat out you or one
       * of your siblings. Find which of your children best matches
       * your parent. A good match here considers barrier nodes, has
       * high overlap, is deep, and has low load. */
      child = bestMatch (nicestart[PARENT], niceend[PARENT],
			 nicedepth[PARENT], niceload[PARENT],
			 FALSE, TRUE, TRUE, FALSE);
      ovlp1 = computeOverlap (nicestart[PARENT], niceend[PARENT],
			      nicestart[child], niceend[child]);
      ovlp2 = computeOverlap (nicestart[SELF], niceend[SELF],
			      nicestart[child], niceend[child]);

      /* You should promote if selected child has a better overlap
       * with your parent, or if your child has equivalent overlap but
       * is less loaded or deeper than you are (recall your depth is the
       * minimum of your childrens' depth). */
      if ((ovlp1 > ovlp2) ||
	  (ovlp1 == ovlp2 &&
	   ((nicedepth[child] >= nicedepth[SELF]) ||
	    (niceload[child] < niceload[SELF]))))
	return (child);
    }

  /* No dice; ERROR indicates no candidate for promotion. */
  return (ERROR);
}
#endif

/**********************************************************************
 * Drop a child from your child list. Called from serveNiceUnregister
 * and niceLoad when promoting a child, among other places. Shifts the
 * last child into the ith position and decrements children. 
 *********************************************************************/
void
dropChild (int i)
{
  /* Log output. */
  if (debug &&
      (errhp = gethostbyaddr ((char *) &niceaddr[i], ADDRBYTES, AF_INET)))
    niceLog (NLOGMSG, "Dropping child %s\n", errhp->h_name);
  else
    niceLog (NLOGMSG, "Dropping child %s\n", inet_ntoa (niceaddr[i]));

  /* Update state by moving the last child to the ith position and
   * then decrementing nicehcnt. */
  if (i < nicehcnt)
    {
      niceload[i] = niceload[nicehcnt];
      nicewall[i] = nicewall[nicehcnt];
      nicedepth[i] = nicedepth[nicehcnt];
      nicelast[i] = nicelast[nicehcnt];
      nicestart[i] = nicestart[nicehcnt];
      niceend[i] = niceend[nicehcnt];
      memcpy (&niceaddr[i], &niceaddr[nicehcnt], ADDRBYTES);
    }

  /* Reset default values for last child record and decrement
   * nicehcnt. The only host table entry that needn't be reset is
   * nicelast[], since this is set explicitly at time of contact. */
  niceload[nicehcnt] = 0;
  nicewall[nicehcnt] = FALSE;
  nicedepth[nicehcnt] = 0;
  nicestart[nicehcnt] = 0;
  niceend[nicehcnt] = 0;
  memcpy (&niceaddr[nicehcnt], &nulladdr, ADDRBYTES);
  nicehcnt--;
}

/**********************************************************************
 **********************************************************************
 * Slave control functions.
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * Cause any running slaves to go to sleep.  Called when we exit the
 * legal time interval. Also updates the slave PID table, in the event
 * that some of these processes have terminated since last checked.
 **********************************************************************/
int
pauseSlaves ()
{
  int i = 0;

  while (i < slaves)
    if (kill (spids[i], SIGSTOP) == ERROR)
      {
	/* If the process no longer exists, an error is
	 * returned. Clear the PID table and reduce the slave count by
	 * one. */
	spids[i] = spids[slaves - 1];
	slaves--;
      }
    else			/* Go on to the next slave. */
      i++;

  return (i);
}

/**********************************************************************
 * Cause any sleeping slaves to wake up.  If slaves are awake, has
 * essentially no effect. Called when we reenter the legal time
 * interval, but is also called simply to update the slave PID table,
 * in the event that some of these processes have terminated (without
 * properly notifying us!) since last checked (that's why we start
 * with a waitpid()).
 **********************************************************************/
int
prodSlaves ()
{
  int i = 0;

  /* If the process no longer exists, waitpid will return the pid and
   * allow the zombie to expire. Otherwise, go ahead and send a
   * wake-up signal to the process. If kill returns an error, then the
   * process is not yours to awaken. In either case, reorganize the
   * PID table and reduce the slave count by one. */
  while (i < slaves)
    if ((waitpid (spids[i], NULL, WNOHANG) == spids[i]) ||
	(kill (spids[i], SIGCONT) == ERROR))
      spids[i] = spids[--slaves];	/* Move last slave and reduce total. */
    else
      i++;			/* Go on to the next slave. */

  return (i);
}

/*************************************************************************
 * Establish a socket connection to specified host's NICE port, and
 * update that host's last connection time entry so that we can keep
 * track of correspondants that have expired. Used when initiating any
 * daemon/daemon protocol. The function will set the socket identifier
 * as a side effect, and will return the status of the connection (the
 * first message returned by the server) which will always be TRUE (OK
 * to continue), FALSE (alternate continuation) or ERROR (message
 * rejected; socket closed, go away -- or, alternatively, failure to
 * connect because of timeout).
 *************************************************************************/
int
establishConnection (int hid, int qid, char *message, int *socket)
{
  int status;
#ifndef NOALARM
  int clientAlarm;
#endif

#ifndef NOALARM
  /* Set a timer in case something blocks while attempting to
   * connect. The initial connection should be relatively fast, like
   * when you authenticate, so a relatively short timeout is
   * appropriate. */
  if ((clientAlarm = niceSetTimeout (CONNINT)) == FALSE)
    {
      /* The other end of the interaction has timed out. You shouldn't
       * need to clear the timer, which we know has already expired,
       * since that's how we got here. What we should do is close the
       * socket and return ERROR, indicating that you were unable to
       * establish a connection. */
      close (*socket);
      return (ERROR);
    }
#endif

  /* Connect to the socket on the appropriate host, which corresponds
   * to niceaddr[i] in our internal hosts tables (or soon will, if
   * this is a new connection). Note that hailSocket() may return
   * ERROR (fatal or "hard" failure, e.g., connection refused) or
   * FALSE ("soft" failure, e.g., host busy) but that, in either case,
   * we'll treat this as if the connection was simply refused by the
   * target. In addition, hailSocket() might hang, in which case we
   * want to timeout CONNINT seconds later. */
  if (hailSocket (hid, port, socket) != TRUE)
    {
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      return (ERROR);
    }

  /* Got a connection.  Send the appropriate query identifier. */
  if (sockprintf (*socket, message, NICEQRYSTRD, qid) == ERROR)
    {
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      close (*socket);
      return (ERROR);
    }

  /* OK, all systems go. Mark this as a recent contact. */
  time (&nicelast[hid]);

  /* Now wait for an acknowledgement and return the appropriate
   * status. However, if the message you just sent was an unregister
   * message, there will be no response so just exit with appropriate
   * status. */
  if (qid == NQIDUNREGISTER)
    {
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      return (TRUE);
    }
  else if (sockscanf (*socket, message, NICEACKSTRD, &status) != 1)
    {
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      niceLog (NLOGWRN, "Connection failure [%d]\n", status);
      close (*socket);
      return (ERROR);
    }
  else
    {
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      return (status);
    }
}

/**********************************************************************
 * Obtain nice daemon configuration. Configuration information comes
 * from three sources: the command line, the configuration file, and
 * from compile-time defaults. Values on the command line should get
 * the highest priority. Note that only a few parameters can be set on
 * the command line, and all but the configuration file name
 * (obviously) can be read from a configuration file.
 **********************************************************************/
int
niceConfig (int argc, char *argv[], char *parent, int *priority)
{
  char cfgfile[NICEMAXSTRINGLEN + 1];	/* Config file name. */
  struct stat fstatus;
  FILE *cfgp;

  char line[NICEMAXSTRINGLEN + 1];
  int i, hours, minutes, verbose;

  /* The following variables need to be set and vetted for errors in
   * this function: parent[], nicestart[SELF], niceend[SELF],
   * maxChildren, maxSlaves, cardinality, priority, debug, and port
   * (some of these variables will already have some default value set
   * at compilation time). Our general strategy is as follows.
   *
   * First, check the command line for a configuration file name, or
   * use the default configuration file.
   *
   * Next, open and parse the configuration file, setting variables
   * appropriately as you go.
   *
   * Finally, process the remaining command line arguments, so as to
   * override any configuration file settings with the command line
   * settings. 
   *
   * Throughout, we want to avoid arbitrary length writes. So pay
   * special attention to sprintf (use snprintf) and sscanf to avoid
   * buffer overflows. */

  /* Set the default name of the configuration file (the configuration
   * directory is set in the top-level Makefile). */
  snprintf (cfgfile, NICEMAXSTRINGLEN, "%s/%s", CFGDIR, CFGFILE);

  /* Scan the command line arguments looking for the name of an
   * alternate configuration file to override the default. We'll also
   * take this opportunity to trap the -? argument, which, if found,
   * produces a brief usage message and exits without further ado. */
  for (i = 1; i < argc; i++)
    {
      if ((strcmp (argv[i], "-c") == 0)
	  || (strcmp (argv[i], "--config") == 0))
	{
	  /* Found one: set the configuration file name
	   * accordingly. NB: the last configuration file name
	   * actually encountered is used. */
	  strncpy (cfgfile, argv[i + 1], NICEMAXSTRINGLEN);
	  break;
	}
      else if ((strcmp (argv[i], "-?") == 0)
	       || (strcmp (argv[i], "--help") == 0))
	{
	  /* Print usage message then return ERROR, causing the daemon
	   * to exit immediately. */
	  printf
	    ("Usage: %s [-c file] [-n port] [-p parent] [-b bfactor] [-d] [-?]\n",
	     argv[0]);
	  printf ("where -? : this output\n");
	  printf ("      -c : specify configuration file\n");
	  printf ("      -n : specify port number\n");
	  printf ("      -p : specify parent host\n");
	  printf ("      -b : specity branching factor\n");
	  printf ("      -d : debug mode (no daemon).\n");
	  return (ERROR);
	}
    }
  /* Check to see if the configuration file is readable; if it isn't,
   * generate an error and exit (note that since we haven't configured
   * message logging yet, the best we can do is complain on stderr). */
  if (access (cfgfile, R_OK) == ERROR)
    {
      fprintf (stderr, "%s: file missing or unreadable\n", cfgfile);
      return (ERROR);
    }
  /* Make sure the file is a regular file; if it isn't, generate an
   * error and exit (note that since we haven't configured message
   * logging yet, the best we can do is complain on stderr). */
  stat (cfgfile, &fstatus);
  if (!S_ISREG (fstatus.st_mode))
    {
      fprintf (stderr, "%s: not regular file\n", cfgfile);
      return (ERROR);
    }
  /* Open the file for read; if this fails, generate an error and exit
   * (note that since we haven't configured message logging yet, the
   * best we can do is complain on stderr).  */
  if (!(cfgp = fopen (cfgfile, "r")))
    {
      fprintf (stderr, "%s: cannot open for read\n", cfgfile);
      return (ERROR);
    }

  /* Now read in the configuration file contents. While not EOF,
   * process each line one at a time. */
  while (fgets (line, NICEMAXSTRINGLEN, cfgp))
    {
      /* Flush comments. */
      sscanf (line, "#%*s");

      /* Name of parent as specified in configuration file. Use same
       * limits as NICESAFEARGS. */
      sscanf (line, "Parent %1024s", parent);

      /* Nice port number as specified in configuration file. Use same
       * limits as NICESAFEARGHD. */
      sscanf (line, "Port %5hd", &port);

      /* Time host becomes available, expressed in local time. */
      if (sscanf (line, "StartTime %2d:%2d", &hours, &minutes))
	nicestart[SELF] = TIME (hours, minutes);

      /* Time host becomes unavailable, expressed in local time. */
      if (sscanf (line, "EndTime %2d:%2d", &hours, &minutes))
	niceend[SELF] = TIME (hours, minutes);

      /* Initial max number of children for this daemon. Limited to
       * MAXCHILDREN = MAXHOSTS-3, since we need to reserve space for
       * SELF, PARENT, and ROOT in the host table. Also, the current
       * maxChildren value will eventually be obtained during the
       * registration process and maintained during niceLoad to allow
       * the root daemon to reconfigure the entire hierarchy. */
      sscanf (line, "MaxChildren %10d", &maxChildren);

      /* Max number of slave processes for this host. Limited to
       * MAXSLAVES, the size of the slave PID table. */
      sscanf (line, "MaxSlaves %10d", &maxSlaves);

      /* Are you a barrier machine? Used to segregate parts of the
       * hierarchy; especially useful when IP masquerading is being
       * used on a cluster. The default is FALSE or 0, so if no
       * barrier specification is found in the config file, the
       * machine is assumed not to be a barrier machine. */
      sscanf (line, "Barrier %1d", &nicewall[SELF]);

      /* Number of naggers to spawn per request. Also limited to
       * MAXSLAVES, the size of the slave PID table. If your
       * cardinality is <=0, then it means you're to act as a firewall
       * machine and fire up the proxy server instead. */
      sscanf (line, "Cardinality %10d", &cardinality);

      /* Daemon scheduling priority. Make sure it is within legal
       * range. */
      sscanf (line, "Priority %2d", priority);
      *priority = MAX (MIN (*priority, 19), -20);

      /* Verbosity of daemon output. */
      sscanf (line, "Verbose %1d", &verbose);
    }
  /* Close config file. */
  fclose (cfgp);

  /* Now process any additional command line arguments, so as to
   * override any parameters you've already set in the configuration
   * file. We're looking for port specifications (-n number), parent
   * specification (-p host), and the debug argument (-d) that
   * precludes daemonizing the niced process. We also need to trap any
   * bad arguments at this point. */
  for (i = 1; i < argc; i++)
    {
      if ((strcmp (argv[i], "-c") == 0)
	  || (strcmp (argv[i], "--config") == 0))
	{
	  /* Configuration file name has already been handled; skip
	   * file name and process the next argument. */
	  i++;
	}
      else if ((strcmp (argv[i], "-n") == 0)
	       || (strcmp (argv[i], "--port") == 0))
	{
	  /* New port number. Use same limits as NICESAFEARGHD. */
	  i++;
	  sscanf (argv[i], "%5hd", &port);
	}
      else if ((strcmp (argv[i], "-p") == 0)
	       || (strcmp (argv[i], "--parent") == 0))
	{
	  /* New parent specification. Use same limits as
	   * NICESAFEARGS. */
	  i++;
	  sscanf (argv[i], "%1024s", parent);
	}
      else if ((strcmp (argv[i], "-b") == 0)
	       || (strcmp (argv[i], "--bfactor") == 0))
	{
	  /* Initial max number of children for this daemon. Limited
	   * to MAXCHILDREN = MAXHOSTS-3, since we need to reserve
	   * space for SELF, PARENT, and ROOT in the host table. Also,
	   * the current maxChildren value will eventually be obtained
	   * during the registration process and maintained during
	   * niceLoad to allow the root daemon to reconfigure the
	   * entire hierarchy. Input size limit matches matches that
	   * for MaxChildren in config file parsing. */
	  i++;
	  sscanf (argv[i], "%10d", &maxChildren);
	}
      else if ((strcmp (argv[i], "-d") == 0)
	       || (strcmp (argv[i], "--debug") == 0))
	{
	  /* Don't daemonize; route log messages to stderr rather than
	   * syslog. */
	  debug = TRUE;
	}
      else if ((strcmp (argv[i], "-?") != 0)
	       && (strcmp (argv[i], "--help") != 0))
	{
	  /* Not a recognized option. Return an error; seems a shame
	   * to bomb out this late in the game, but it is harder to
	   * check for unrecognized options earlier. */
	  return (ERROR);
	}
    }

  /* Setup message logging. If you're not running as a system daemon,
   * report errors directly on STDERR (commonly used for debugging
   * purposes). Otherwise, send messages to system log. */
  niceLogSet (verbose, (debug ? NLOGSTDERR : NLOGSYSLOG));

  /* Finally, now that message logging is configured, we can check
   * and, if need be, adjust the final configuration while outputting
   * warning messages to their proper location.
   *
   * Max number of children for this daemon. Limited to MAXCHILDREN =
   * MAXHOSTS-3, since we need to reserve space for SELF, PARENT, and
   * ROOT in the host table. 
   *
   * We also need to ensure maxChildren is at least 1, or we shatter
   * the NICE hierarchy. */
  if (maxChildren > MAXCHILDREN)
    {
      niceLog (NLOGMSG, "Reducing MaxChildren to %d.\n", MAXCHILDREN);
      maxChildren = MAXCHILDREN;
    }
  else if (maxChildren < 1)
    maxChildren = 1;

  /* Max number of slave processes for this host. Limited to
   * MAXSLAVES, the size of the slave PID table. */
  if (maxSlaves > MAXSLAVES)
    {
      niceLog (NLOGMSG, "Reducing MaxSlaves to %d.\n", MAXSLAVES);
      maxSlaves = MAXSLAVES;
    }

  /* Number of naggers to spawn per request. Also limited to
   * MAXSLAVES, the size of the slave PID table. If your cardinality
   * is <=0, then it means you're to act as a firewall machine and
   * fire up the proxy server instead. In this case, ensure the daemon
   * only runs one child. */
  if (cardinality > MAXSLAVES)
    {
      niceLog (NLOGMSG, "Reducing cardinality to %d.\n", MAXSLAVES);
      cardinality = MAXSLAVES;
    }
  else if (cardinality <= 0)
    {
      niceLog (NLOGMSG, "Using NICE proxy server.\n");
      cardinality = 0;
      maxChildren = 1;
    }

  /* Finished! If you've gotten this far, return TRUE. */
  return (TRUE);
}

/**********************************************************************
 * NICE daemon protocols. 
 **********************************************************************/
/**********************************************************************
 * Register a NICE daemon with specified host. If the desired parent
 * is full, the request may be forwarded elsewhere. Returns TRUE if
 * registration was successful, FALSE if registration failed but
 * should be attempted once more, and ERROR if registration failed and
 * has no possibility of success (e.g., if some authentication problem
 * has occurred).
 *
 * Note that within daemond, this protocol is theoretically the only
 * possibility of deadlock, since it is initiated by a child (wheras
 * all other protocols are parent initiated). The alarms help to break
 * out of that possibility.
 **********************************************************************/
int
niceRegister (int hid)
{
  int socket;
  int status;
  char message[NICEMAXSTRINGLEN + 1];
  char ipaddr[NICEMAXSTRINGLEN + 1];	/* In dotted octet form. */
#ifndef NOALARM
  int clientAlarm;
#endif
#ifndef NODYNAMIC
  int child;
#endif

  if (debug &&
      (errhp = gethostbyaddr ((char *) &niceaddr[hid], ADDRBYTES, AF_INET)))
    niceLog (NLOGMSG, "Registering with %s: ", errhp->h_name);
  else
    niceLog (NLOGMSG, "Registering with %s: ", inet_ntoa (niceaddr[hid]));

  /* Open a connection with the prospective parent and send a
   * registration message. Here, establishConnection() will return
   * either TRUE if the connection can be made or FALSE if the host is
   * unreachable. */
  if (establishConnection (hid, NQIDREGISTER, message, &socket) == ERROR)
    {
      /* The parent may be unreachable, or there may be no daemon
       * running on the parent. In any case, generate an error and
       * return. Typically, the daemon will go into orphan mode and
       * continue to try to register. */
      niceLog (NLOGMSG, "failure.\n");
      return (FALSE);
    }

#ifndef NOALARM
  /* Set a timer in case something blocks during the execution of the
   * protocol.A protocol is declared hung/blocked and attempts to read
   * or write are abandoned when it hangs for HUNGINT seconds,
   * typically a large number. */
  if ((clientAlarm = niceSetTimeout (HUNGINT)) == FALSE)
    {
      /* The other end of the interaction has blocked. You shouldn't
       * need to clear the timer, which we know has already expired,
       * since that's how we got here. What we should do is close the
       * socket and return FALSE, indicating that something (nonfatal)
       * has gone awry, so that you can become an orphan and try
       * registering again later. */
      close (socket);
      return (FALSE);
    }
#endif

  /* Next, send your parameters. Some of these are fixed, while others
   * will subsequently need to be updated via the load reporting
   * protocol. All of them can be used by the parent to reassign you
   * if necessary. 
   *
   * TODO: Need to include some authentication parameters. */
  if (sockprintf (socket, message, NICEINFSTRSDDDDD,
		  NICEVERSION, niceload[SELF], nicewall[SELF],
		  nicedepth[SELF], UTIME (nicestart[SELF]),
		  UTIME (niceend[SELF])) == ERROR)
    {
      /* Lost link; return FALSE and become an orphan. */
      close (socket);

#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      /* Log output. */
      niceLog (NLOGMSG, "failure.\n");
      return (FALSE);
    }

  /* Prospective parent now returns his parameters. 
   *
   * TODO: Need to include some authentication parameters. */
  if (sockscanf (socket, message, NICEINFSTRSDDDDD,
		 &ipaddr, &niceload[hid], &nicewall[hid],
		 &nicestart[hid], &niceend[hid], &nicebfact) != 6)
    {
      /* Lost link; return FALSE and become an orphan. */
      close (socket);

#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      /* Log output. */
      niceLog (NLOGMSG, "failure.\n");
      return (FALSE);
    }

  /* Commit: close socket. */
  close (socket);

#ifndef NOALARM
  /* Clear the timer. */
  niceClearAlarm (clientAlarm);
#endif

  /* Now, depending on the value of ipaddr, we can figure out whether
   * the registration was successful or not. If ipaddr is the NULLHOST
   * address, then we've been rejected for authentication reasons and
   * should just exit altogether. If the IP address is LOCALHOST, the
   * registration is OK (note that we continue to use the IP address
   * established from the network). Otherwise, ipaddr is a forwarding
   * address for another attempt. */
  if (strncmp (ipaddr, NULLHOST, NICEMAXSTRINGLEN) == 0)
    {
      /* Log output. */
      niceLog (NLOGMSG, "authentication failed.\n");
      return (ERROR);
    }
  else if (strncmp (ipaddr, LOCALHOST, NICEMAXSTRINGLEN) == 0)
    {
      /* OK, found my parent. Mark the last contact time and make sure
       * your host table reflects the new parent information. We need
       * to do this because hid could have been ROOT rather than
       * parent when we invoked this function, such as, for example,
       * when registering with ROOT after being an orphan. */
      time (&nicelast[hid]);
      memcpy (&niceaddr[PARENT], &niceaddr[hid], ADDRBYTES);
      niceload[PARENT] = niceload[hid];
      nicewall[PARENT] = nicewall[hid];
      nicedepth[PARENT] = nicedepth[hid];
      nicestart[PARENT] = nicestart[hid];
      niceend[PARENT] = niceend[hid];
      nicelast[PARENT] = nicelast[hid];

      /* Next adjust your maxChildren value to reflect the hierarchy's
       * current branching factor. Limited to MAXCHILDREN =
       * MAXHOSTS-3, since we need to reserve space for SELF, PARENT,
       * and ROOT in the host table. */
      if (debug && maxChildren != MIN (nicebfact, MAXCHILDREN))
	niceLog (NLOGMSG, "Resetting branching factor to %d.\n",
		 MIN (nicebfact, MAXCHILDREN));
      maxChildren = MIN (nicebfact, MAXCHILDREN);
      if (nicehcnt > maxChildren)
	{
	  /* OK, we've suffered a branching factor restriction that
	   * will require we jettison a few extraneous children (note
	   * that nicehcnt must be <= maxChildren, and not strictly
	   * less than).
	   *
	   * If we're doing dynamic restructuring, we'll find the
	   * "best" children to demote by considering barrier nodes,
	   * low overlap, shallowness, and high load. Otherwise, we'll
	   * just drop a bunch until we're in compliance. */
	  while (nicehcnt > maxChildren)
	    {
#ifndef NODYNAMIC
	      child = bestMatch (nicestart[SELF], niceend[SELF],
				 nicedepth[SELF], niceload[SELF],
				 FALSE, FALSE, FALSE, TRUE);
	      niceOrphan ((child > 0
			   && child <= nicehcnt) ? child : nicehcnt);
#else
	      niceOrphan (nicehcnt);
#endif
	    }
	}

      /* Log output. */
      niceLog (NLOGMSG, "success.\n");
      return (TRUE);
    }
  else
    {
      /* If you get this far, it means that the prospective parent is
       * full and has directed you to register with ipaddr instead. */
      inet_aton (ipaddr, &niceaddr[PARENT]);

      /* Log output. */
      if (debug &&
	  (errhp = gethostbyaddr ((char *) &niceaddr[PARENT],
				  ADDRBYTES, AF_INET)))
	niceLog (NLOGMSG, "referred to %s.\n", errhp->h_name);
      else
	niceLog (NLOGMSG, "referred to %s.\n", ipaddr);

      /* Now try the new parent. If this fails (i.e., returns FALSE),
       * it is probably because the current parent selected a child
       * that was down or offline, so try ROOT again and hope for a
       * better outcome. */
      if ((status = niceRegister (PARENT)) == TRUE)
	return (TRUE);
      else if (status == FALSE && (status = niceRegister (ROOT) == TRUE))
	return (TRUE);
      else if (status == ERROR)
	return (ERROR);
    }
  return (FALSE);
}

/**********************************************************************
 * serveNiceRegister: Add a new child daemon. Updates internal hosts
 * tables based on local the given argument, which is what the
 * supplicant believes is its own address. TODO: the supplicant might
 * in fact be misconfigured, or might have more than one NIC, which
 * complicates verifying their stated IP address.
 **********************************************************************/
void
serveNiceRegister (int socket, struct in_addr caddr)
{
  int i, newload, newwall, newdepth, newstart, newend;
  int bestchild = 0;
  int status = FALSE;
#ifndef NODYNAMIC
  int ovlp1, ovlp2, worstchild = 0;
#endif
  char newvers[NICEMAXSTRINGLEN + 1];
  char message[NICEMAXSTRINGLEN + 1];

  /* OK, I'm listening. */
  if (sockprintf (socket, message, NICEACKSTRD, TRUE) == ERROR)
    {
      /* Link failed. Close socket and exit. */
      close (socket);
      return;
    }

  /* Must check to see if the new child might already be one of our
   * children. This can happen, for example, when the child dies,
   * reboots, and reregisters with us before we have a chance to
   * notice the child is unresponsive. The easy way to handle this is
   * to drop the child and then treat it exactly like a new
   * child. Important: deadlock could conceivably occur if we're not
   * careful! If we're executing this code as a result of the child
   * giving up on the parent (rather than e.g., rebooting) then you've
   * been really lucky that the parent hadn't initiated a load info
   * request in the interim, which would have led us to a deadlock
   * (the deadlock would eventually be broken with a timeout alarm, if
   * compiled with alarms). To avoid deadlocks, the parent should
   * always give up first in the event of a link failure. */
  i = 0;
  while (++i <= nicehcnt && !ADDREQ (niceaddr[i], caddr));
  if (i <= nicehcnt)
    dropChild (i);

  /* Obtain the supplicant's schedule information. */
  if (sockscanf (socket, message, NICEINFSTRSDDDDD, &newvers, &newload,
		 &newwall, &newdepth, &newstart, &newend) != 6)
    {
      /* Lost link; give it up. */
      close (socket);
      return;
    }

  /* Convert supplicant's times to local time. */
  newstart = LTIME (newstart);
  newend = LTIME (newend);

  /* Now we need to decide what to do with the supplicant. First,
   * perform whatever authentication is required to ensure the
   * supplicant is authorized to connect.  The default is that
   * authorization is denied.
   *
   * To be legitimate, we check the NICE version number to make sure
   * we have the same protocols. We also require that the child IP
   * address is reasonable (e.g., "localhost" is not acceptable). We
   * should probably do a reverse DNS lookup to ensure that the child
   * is who the network says they are, but this would be quite slow,
   * and if you can't trust the network anyway, the DNS lookup is
   * meaningless (i.e., we really need some stronger authorization
   * mechanism). */
  if ((strncmp (newvers, NICEVERSION, NICEMAXSTRINGLEN) == 0) &&
      (ADDRNEQ (caddr, localaddr)))
    status = TRUE;

  /* Once authentication passes, then there are three
   * scenarios. First, we might have room for the supplicant, in which
   * case we just accept it. Second, we may be out of room but the
   * supplicant may be better than an existing child we can demote
   * (DYNAMIC case only). Third, we may elect to forward the
   * supplicant elsewhere. */
  if (status && nicehcnt == maxChildren)
    {
#ifndef NODYNAMIC
      /* Determine whether the supplicant is more desirable than one
       * of your current children. First, find the child you would
       * most like to demote. A good candidate considers barrier
       * nodes, has low overlap, is shallow, and has high load. */
      worstchild = bestMatch (nicestart[SELF], niceend[SELF],
			      nicedepth[SELF], niceload[SELF],
			      FALSE, FALSE, FALSE, TRUE);
      ovlp1 = computeOverlap (nicestart[SELF], niceend[SELF],
			      nicestart[worstchild], niceend[worstchild]);
      ovlp2 = computeOverlap (nicestart[SELF], niceend[SELF],
			      newstart, newend);

      /* Figure out if you have a child to drop -- but wait to drop
       * him until later, when you commit to this new
       * registration. Note that here if ovlp is equal, load is more
       * important than depth, since fresh supplicants are likely to
       * be of low depth, at least initially. */
      if (!((ovlp1 < ovlp2) ||
	    (ovlp1 == ovlp2 &&
	     ((newload < niceload[worstchild]) ||
	      (newload == niceload[worstchild]
	       && newdepth > nicedepth[worstchild])))))
	{
	  /* Keep all of your children and forward supplicant
	   * instead. A good match excludes barrier nodes, has high
	   * overlap, is shallow, and has low load.
	   *
	   * Oops! Not so fast. If there's no place to send the
	   * supplicant because all of your children are barrier nodes
	   * (the supplicant isn't), prefer demoting an arguably
	   * better child and keeping the supplicant, essentially
	   * pushing the barrier nodes towards the leaves of the
	   * hierarchy. You really haven't got any choice. */
	  if (((bestchild = bestMatch (newstart, newend, newdepth, newload,
				       TRUE, TRUE, FALSE, FALSE)) == ERROR) &&
	      (newwall == FALSE))
	    bestchild = 0;
	  else
	    worstchild = 0;
	}
#else
      /* Decide where to forward the supplicant. A good match excludes
       * barrier nodes, has high overlap, is shallow, and has low
       * load. Result will be ERROR if there is no place to go. */
      bestchild =
	bestMatch (newstart, newend, newdepth, newload, TRUE, TRUE, FALSE,
		   FALSE);
#endif
    }
  /* Reply to the supplicant: must be able to indicate (1) acceptance,
   * (2) redirection, (3) rejection because I am full and have nowhere
   * to redirect, and (4) rejection for reasons of authentication.
   *
   * In practice:
   * status==FALSE, send NULLHOST => authentication error/exit.
   * bestchild==0, send LOCALHOST => accept
   * bestchild==-1, send ROOT => redirect to root
   * else, send bestchild IP address => redirect to child
   *
   * TODO: If there aren't enough "open slots" in the hierarchy for
   * all the nodes when barrier nodes are considered, it is possible
   * that you might end up continuously trying to register with the
   * root, in essence thrashing. Not clear how to avoid this.  */
  if (sockprintf (socket, message, NICEINFSTRSDDDDD,
		  ((status == FALSE) ? NULLHOST :
		   ((bestchild == 0) ? LOCALHOST :
		    ((bestchild == ERROR) ? inet_ntoa (niceaddr[ROOT]) :
		     inet_ntoa (niceaddr[bestchild])))),
		  niceload[SELF], nicewall[SELF],
		  nicestart[SELF], niceend[SELF], nicebfact) == ERROR)
    {
      /* Lost link; give it up. */
      close (socket);
      return;
    }

  /* Commit: close socket. */
  close (socket);

  /* If supplicant failed to authenticate, you're done. */
  if (status == FALSE)
    {
      /* Log output. */
      if (debug &&
	  (errhp = gethostbyaddr ((char *) &caddr, ADDRBYTES, AF_INET)))
	niceLog (NLOGWRN, "Authentication failure %s\n", errhp->h_name);
      else
	niceLog (NLOGWRN, "Authentication failure %s\n", inet_ntoa (caddr));
      return;
    }

#ifndef NODYNAMIC
  /* In the dynamic case, if you have elected to replace a child with
   * the supplicant, then you need to drop the worst child to make
   * room for the new one. Note: niceOrphan() always succeeds. */
  if (worstchild != 0)
    niceOrphan (worstchild);
#endif

  /* If you accepted the supplicant, you must copy their info into
   * your tables and mark the last contact time. If you redirected the
   * child, you needn't do anything else. */
  if (bestchild == 0)
    {
      nicehcnt++;
      memcpy (&niceaddr[nicehcnt], &caddr, ADDRBYTES);
      niceload[nicehcnt] = newload;
      nicewall[nicehcnt] = newwall;
      nicedepth[nicehcnt] = newdepth;
      nicestart[nicehcnt] = newstart;
      niceend[nicehcnt] = newend;

      /* Mark the last contact time. */
      time (&nicelast[nicehcnt]);

      /* Update your own depth estimate. */
      if (nicehcnt < maxChildren)
	nicedepth[SELF] = 0;
      else
	{
	  nicedepth[SELF] = nicedepth[nicehcnt] + 1;
	  for (i = 1; i < nicehcnt; i++)
	    nicedepth[SELF] = MIN (nicedepth[SELF], nicedepth[i] + 1);
	}

      /* Log output. */
      if (debug &&
	  (errhp = gethostbyaddr ((char *) &caddr, ADDRBYTES, AF_INET)))
	niceLog (NLOGMSG, "Accepting new child %s\n", errhp->h_name);
      else
	niceLog (NLOGMSG, "Accepting new child %s\n", inet_ntoa (caddr));
    }
  else
    {
      /* Log output. */
      if (debug &&
	  (errhp = gethostbyaddr ((char *) &niceaddr[bestchild],
				  ADDRBYTES, AF_INET)))
	niceLog (NLOGMSG, "Referring prospective child to %s\n",
		 errhp->h_name);
      else
	niceLog (NLOGMSG, "Referring prospective child to %s\n",
		 inet_ntoa (niceaddr[bestchild]));
    }
  return;
}

/**********************************************************************
 * Unregister with your parent. A NICE daemon would invoke this before
 * going down, gently. Aside from registering, this is the only
 * child-initiated protocol, and is therefore asynchronous
 * (non-blocking) by design in order to avoid deadlock conditions.
 **********************************************************************/
int
niceUnregister ()
{
  int socket, status;
  char message[NICEMAXSTRINGLEN + 1];

  /* Log output. */
  if (debug &&
      (errhp =
       gethostbyaddr ((char *) &niceaddr[PARENT], ADDRBYTES, AF_INET)))
    niceLog (NLOGMSG, "Unregistering with %s\n", errhp->h_name);
  else
    niceLog (NLOGMSG, "Unregistering with %s\n",
	     inet_ntoa (niceaddr[PARENT]));

  /* Send unregister message, but since this is an asynchronous
   * protocol you don't wait around for a reply (nor do you need to
   * set a timer if NOALARM is true). Here, establishConnection() will
   * either return TRUE (there's a socket to close) or ERROR (no
   * connection made, parent may be dead or unreachable, no socket to
   * close). */
  if ((status = establishConnection (PARENT, NQIDUNREGISTER, message,
				     &socket)) != ERROR)
    close (socket);

  /* Return status. */
  return (status);
}

/**********************************************************************
 * serveNiceUnregister: Remove a child daemon. Updates internal hosts
 * table. This is a non-blocking protocol; there is no acknowledgement.
 **********************************************************************/
void
serveNiceUnregister (int socket, struct in_addr caddr)
{
  int i = 1;			/* Only looking at children. */

  /* Close the socket immediately; this protocol is non-blocking. */
  close (socket);

  /* Find the child in question. */
  while ((i <= nicehcnt) && !(ADDREQ (niceaddr[i], caddr)))
    i++;

  if (i <= nicehcnt)
    {
      if (debug &&
	  (errhp = gethostbyaddr ((char *) &caddr, ADDRBYTES, AF_INET)))
	niceLog (NLOGMSG, "Child[%d] %s unregistered.\n", i, errhp->h_name);
      else
	niceLog (NLOGMSG, "Child[%d] %s unregistered.\n",
		 i, inet_ntoa (caddr));
      dropChild (i);
    }
  else
    {
      if (debug &&
	  (errhp = gethostbyaddr ((char *) &caddr, ADDRBYTES, AF_INET)))
	niceLog (NLOGWRN, "Child %s not found to unregister.\n",
		 errhp->h_name);
      else
	niceLog (NLOGWRN, "Child %s not found to unregister.\n",
		 inet_ntoa (caddr));
    }
  return;
}

/**********************************************************************
 * Tell a child NICE daemon to go register with a different parent. A
 * NICE daemon sends this message to a selected child during dynamic
 * restructuring. Return ERROR if the child can't be contacted, FALSE
 * if the child refuses the promotion, and TRUE if child accepts it.
 **********************************************************************/
#ifndef NODYNAMIC
int
nicePromote (int child, int parent)
{
  int socket, status;
  char message[NICEMAXSTRINGLEN + 1];
#ifndef NOALARM
  int clientAlarm;
#endif

  /* Open a connection with the child and send a promotion message. */
  status = establishConnection (child, NQIDPROMOTE, message, &socket);
  if (status == ERROR)
    {
      /* Child unresponsive; forget the whole idea. */
      if (debug &&
	  (errhp =
	   gethostbyaddr ((char *) &niceaddr[child], ADDRBYTES, AF_INET)))
	niceLog (NLOGWRN, "Can't promote %s.\n", errhp->h_name);
      else
	niceLog (NLOGWRN, "Can't promote %s.\n", inet_ntoa (niceaddr[child]));
      return (ERROR);
    }
  else if (parent != SELF && status == TRUE)
    {
      /* I had intended to promote this child, but I can't because it
       * believes I am also its root node (returns TRUE
       * status). Either I or my child are probably misconfigured
       * (shouldn't I be a barrier node? or should child point to
       * root?), but for now, forget the whole thing. */
      if (debug &&
	  (errhp =
	   gethostbyaddr ((char *) &niceaddr[child], ADDRBYTES, AF_INET)))
	niceLog (NLOGWRN, "Child %s refuses promotion.\n", errhp->h_name);
      else
	niceLog (NLOGWRN, "Child %s refuses promotion.\n",
		 inet_ntoa (niceaddr[child]));
      close (socket);
      return (FALSE);
    }

  /* All systems go. Log output.  Since inet_ntoa uses a static
   * string, we need to momentarily stash one of the addresses
   * somewhere in order to print out the error message. */
  if (debug &&
      (errhp = gethostbyaddr ((char *) &niceaddr[child], ADDRBYTES, AF_INET)))
    niceLog (NLOGMSG, "Promoting %s.\n", errhp->h_name);
  else
    niceLog (NLOGMSG, "Promoting %s.\n", inet_ntoa (niceaddr[child]));

#ifndef NOALARM
  /* Set a timer in case something blocks during the execution of the
   * remainder of the protocol. A protocol is declared hung/blocked
   * and attempts to read or write are abandoned when it hangs for
   * HUNGINT seconds, typically a large number. */
  if ((clientAlarm = niceSetTimeout (HUNGINT)) == FALSE)
    {
      /* The other end of the interaction has blocked. You shouldn't
       * need to clear the timer, which we know has already expired,
       * since that's how we got here. What we should do is close the
       * socket and return ERROR, indicating that something has gone
       * awry, and that the child will not be promoted.  */
      close (socket);
      return (ERROR);
    }
#endif

  /* Send child its new parent hostaddr. */
  if (sockprintf (socket, message, NICEINFSTRS,
		  inet_ntoa (niceaddr[parent])) == ERROR)
    {
      /* Lost link: give it up and trust the child. */
      close (socket);
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      return (ERROR);
    }

  /* OK; commit. Close socket, and wipe the child out of your internal
   * tables. Also update your own depth estimate to reflect the fact
   * you have an opening. */
  close (socket);
#ifndef NOALARM
  /* Clear the timer. */
  niceClearAlarm (clientAlarm);
#endif
  dropChild (child);
  nicedepth[SELF] = 0;

  /* Return success. */
  return (status);
}

/**********************************************************************
 * serveNicePromote: Handle promotion request. Updates internal hosts
 * table.
 **********************************************************************/
void
serveNicePromote (int socket, struct in_addr caddr)
{
  char message[NICEMAXSTRINGLEN + 1];
  char parent[NICEMAXSTRINGLEN + 1];

  /* OK, I'm listening. If my current parent is also my root node,
   * then I will reject a promotion message, since I do not wish to be
   * promoted outside my locally-specified hierarchy. */
  if (sockprintf (socket, message, NICEACKSTRD,
		  (ADDREQ (niceaddr[ROOT], niceaddr[PARENT]))) == ERROR)
    {
      /* Link failed. Close socket and exit. */
      close (socket);
      return;
    }
  else if (ADDREQ (niceaddr[ROOT], niceaddr[PARENT]))
    {
      /* Promotion rejected. Close socket and exit. */
      close (socket);
      return;
    }

  /* Get address of your new parent. */
  if (sockscanf (socket, message, NICEINFSTRS, parent) != 1)
    {
      /* Link failed. Close socket and exit. */
      close (socket);
      return;
    }

  /* All systems go; commit. Convert new parent address to addr
   * structure, overwriting the old parent address that we won't be
   * needing any longer. */
  inet_aton (parent, &niceaddr[PARENT]);

  /* Close the socket. */
  close (socket);

  if (debug &&
      (errhp =
       gethostbyaddr ((char *) &niceaddr[PARENT], ADDRBYTES, AF_INET)))
    niceLog (NLOGMSG, "Promoted: registering with new parent %s\n",
	     errhp->h_name);
  else
    niceLog (NLOGMSG, "Promoted: registering with new parent %s\n",
	     inet_ntoa (niceaddr[PARENT]));

  /* Now you're ready to register with the new parent, or, failing
   * that, with the root host. If both fail, become an ORPHAN and keep
   * trying to attach to one or the other until successful. 
   *
   * TODO: should probably exit if niceRegister returns ERROR, no? */
  if (niceRegister (PARENT) == FALSE &&
      (ADDREQ (niceaddr[PARENT], niceaddr[ROOT])
       || niceRegister (ROOT) == FALSE))
    mode = ORPHAN;

  /* Finished. */
  return;
}
#endif

/**********************************************************************
 * Make a child NICE daemon an orphan. Always succeeds!
 **********************************************************************/
void
niceOrphan (int child)
{
  int socket;
  char message[NICEMAXSTRINGLEN + 1];

  /* Open a connection with the child and send an orphan
   * message. Here, establishConnection() may return an ERROR (can't
   * reach the child), FALSE (child doesn't think you're its real
   * parent), or TRUE (acknowledge). While this isn't exactly an
   * asynchronous protocol (like niceUnregister), since we don't
   * really do anything after the connection is established, it isn't
   * necessary to set a timer. */
  if (establishConnection (child, NQIDORPHAN, message, &socket) == ERROR)
    {
      /* Child unresponsive; there's not much we can do about this. If
       * we're about to go down (one reason for invoking niceOrphan())
       * then it doesn't really much matter what we do here.  */
      if (debug &&
	  (errhp =
	   gethostbyaddr ((char *) &niceaddr[child], ADDRBYTES, AF_INET)))
	niceLog (NLOGWRN, "Abruptly orphaning %s.\n", errhp->h_name);
      else
	niceLog (NLOGWRN, "Abruptly orphaning %s.\n",
		 inet_ntoa (niceaddr[child]));
    }
  else
    {
      if (debug &&
	  (errhp =
	   gethostbyaddr ((char *) &niceaddr[child], ADDRBYTES, AF_INET)))
	niceLog (NLOGMSG, "Orphaning %s\n", errhp->h_name);
      else
	niceLog (NLOGMSG, "Orphaning %s\n", inet_ntoa (niceaddr[child]));
    }

  /* establishConnection() returns TRUE (listening) or FALSE (you
   * aren't my parent anyway).  Close socket. */
  close (socket);

  /* Always drop the child, regardless of whether or not you managed
   * to successfully tell the child it was going to do so. */
  dropChild (child);

  /* Finished. */
  return;
}

/**********************************************************************
 * serveNiceOrphan: services a request from your parent to turn
 * yourself into an orphan. This only happens when your parent is the
 * root daemon and is going down, hopefully temporarily. By turning
 * into an orphan, you speed up the process of restoring the hierarchy
 * when your parent comes back online.
 **********************************************************************/
void
serveNiceOrphan (int socket, struct in_addr caddr)
{
  char message[NICEMAXSTRINGLEN + 1];

  /* Make sure this is really your parent; then become an orphan.  No
   * need to check exit status, since, once we've gotten this far,
   * we're not going to act any differently if the link is lost. */
  if (ADDRNEQ (caddr, niceaddr[PARENT]))
    {
      sockprintf (socket, message, NICEACKSTRD, FALSE);
      close (socket);
      /* Log output. */
      niceLog (NLOGWRN, "Can't be abandoned by %s (!= %s).\n",
	       inet_ntoa (caddr), inet_ntoa (niceaddr[PARENT]));
      return;
    }

  /* Acknowledge. */
  sockprintf (socket, message, NICEACKSTRD, TRUE);

  /* Commit: close the socket. */
  close (socket);

  /* Log output. */
  if (debug &&
      (errhp =
       gethostbyaddr ((char *) &niceaddr[PARENT], ADDRBYTES, AF_INET)))
    niceLog (NLOGMSG, "Abandoned by %s.\n", errhp->h_name);
  else
    niceLog (NLOGMSG, "Abandoned by %s.\n", inet_ntoa (niceaddr[PARENT]));

  /* Set your mode to orphan. Once you've been orphaned, you'll
   * attempt to reconnect with your parent or the root node. Since
   * your parent just explicitly dropped you, it would no doubt be
   * better to prefer registering with the root; this is a design
   * decision. */
  mode = ORPHAN;
  memcpy (&niceaddr[PARENT], &niceaddr[ROOT], ADDRBYTES);
  return;
}

/**********************************************************************
 * Initiate exchange of load information with your children.
 **********************************************************************/
void
niceLoad ()
{
  int i, j, socket, status, wishes;
  char message[NICEMAXSTRINGLEN + 1];
  char name[NICEMAXSTRINGLEN + 1];
  char filename[NICEMAXSTRINGLEN + 1];
#ifndef NODYNAMIC
  int candidate;
#endif
#ifndef NOALARM
  int clientAlarm;
#endif

  /* Establish a connection with each child in turn and exchange load
   * and status information. */
  i = 1;
  while (i <= nicehcnt)		/* 0 < i <= nicehcnt */
    {
      /* Here, establishConnection() will return TRUE (child will send
       * load information), FALSE (you're not my real parent; drop
       * me!), or ERROR (can't reach child). */
      status = establishConnection (i, NQIDLOAD, message, &socket);
      if (status == FALSE)
	{
	  /* Your child rejects you: you are not its real parent. */
	  if (debug &&
	      (errhp = gethostbyaddr ((char *) &niceaddr[i],
				      ADDRBYTES, AF_INET)))
	    niceLog (NLOGWRN, "Spurned by pesudochild %s\n", errhp->h_name);
	  else
	    niceLog (NLOGWRN, "Spurned by pseudochild %s\n",
		     inet_ntoa (niceaddr[i]));

	  /* Close socket and drop child like a hot potato. Don't
	   * increment i, to ensure that you get the child that shifts
	   * into the dropped child's slot next time around. */
	  close (socket);
	  dropChild (i);
	}
      else if (status != ERROR)
	{
#ifndef NOALARM
	  /* Set a timer in case something blocks during the execution
	   * of the remainder of the protocol. A protocol is declared
	   * hung/blocked and attempts to read or write are abandoned
	   * when it hangs for HUNGINT seconds, typically a large
	   * number. */
	  if ((clientAlarm = niceSetTimeout (HUNGINT)) == FALSE)
	    {
	      /* The other end of the interaction has blocked. You
	       * shouldn't need to clear the timer, which we know has
	       * already expired, since that's how we got here. What
	       * we should do is close the socket and go on to the
	       * next child. */
	      close (socket);
	      continue;
	    }
#endif
	  /* Send your depth and load information, along with
	   * hierarchy's current branching factor. */
	  if (sockprintf (socket, message, NICEINFSTRDDD,
			  nicedepth[SELF], niceload[SELF],
			  nicebfact) == ERROR)
	    {
	      /* Lost link: use "continue" to skip this child. */
	      if (debug &&
		  (errhp = gethostbyaddr ((char *) &niceaddr[i],
					  ADDRBYTES, AF_INET)))
		niceLog (NLOGMSG, "No status report from %s\n",
			 errhp->h_name);
	      else
		niceLog (NLOGMSG, "No status report from %s\n",
			 inet_ntoa (niceaddr[i]));
	      close (socket);
#ifndef NOALARM
	      /* Clear the timer. */
	      niceClearAlarm (clientAlarm);
#endif
	      continue;
	    }

	  /* Obtain child's depth, load and availability information
	   * as long as the number of applications it would like to
	   * download. */
	  if (sockscanf (socket, message, NICEINFSTRDDD,
			 &nicedepth[i], &niceload[i], &wishes) != 3)
	    {
	      /* Lost link: use "continue" to skip this child. */
	      if (debug &&
		  (errhp = gethostbyaddr ((char *) &niceaddr[i],
					  ADDRBYTES, AF_INET)))
		niceLog (NLOGMSG, "No status report from %s\n",
			 errhp->h_name);
	      else
		niceLog (NLOGMSG, "No status report from %s\n",
			 inet_ntoa (niceaddr[i]));
	      close (socket);
#ifndef NOALARM
	      /* Clear the timer. */
	      niceClearAlarm (clientAlarm);
#endif
	      continue;
	    }

	  /* Fetch your child's currently inactive wishes, one by one,
	   * and place each wish on your putlist (if you have the
	   * requested application) or on your getlist (if you
	   * don't). */
	  j = 0;
	  while (j < wishes)
	    {
	      /* Acknowledge. On failure, you justexit the inner loop;
	       * since this is the last thing before you stop talking
	       * to this child, you needn't exit the outer loop
	       * explicitly, nor close the socket or clear the alarm
	       * (these things happen as you exit the outer loop). */
	      if (sockprintf (socket, message, NICEACKSTRD, TRUE) == ERROR)
		break;

	      /* Get child's Nth currently inactive get wishes. */
	      if (sockscanf (socket, message, NICEINFSTRS, &name) != 1)
		break;

	      /* If you have this application, make a put wish. If you
	       * don't, make a get wish (so that eventually you can
	       * pass it on). */
	      if (niceCheckApplication (name, filename) == TRUE)
		{
		  niceLog(NLOGMSG, "Adding PUT request for %s.\n", name); 
		  addPutWish (name, j);
		}
	      else
		{
		  niceLog(NLOGMSG, "Adding GET request for %s.\n", name); 
		  addGetWish (name, NULL, NULL);
		}
	      i++;
	    }

	  /* Finished with this child. Close the socket. */
	  close (socket);
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (clientAlarm);
#endif
	  /* Log output. */
	  if (debug &&
	      (errhp = gethostbyaddr ((char *) &niceaddr[i],
				      ADDRBYTES, AF_INET)))
	    niceLog (NLOGMSG, "Getting status (%d) from %s\n",
		     niceload[i], errhp->h_name);
	  else
	    niceLog (NLOGMSG, "Getting status (%d) from %s\n",
		     niceload[i], inet_ntoa (niceaddr[i]));
	  i++;
	}
      else
	{
	  if (debug &&
	      (errhp =
	       gethostbyaddr ((char *) &niceaddr[i], ADDRBYTES, AF_INET)))
	    niceLog (NLOGMSG, "No status report from %s\n", errhp->h_name);
	  else
	    niceLog (NLOGMSG, "No status report from %s\n",
		     inet_ntoa (niceaddr[i]));
	  i++;
	}
    }

  /* Your depth is recalculated each time you measure your childrens'
   * load. Initialize current depth calculation: if you have an open
   * slot, your depth is zero by definition. Otherwise, the initial
   * ERROR value means we're full up but don't yet know what each
   * child's current depth is. */
  nicedepth[SELF] = ((nicehcnt < maxChildren) ? 0 : ERROR);

  /* Update depth estimate from each child's current depth
   * information. After the loop is finished, your depth indicates
   * distance, in levels, to closest vacancy among your descendents,
   * with 0 meaning you have a vacancy. */
  for (i = 0; i < nicehcnt; i++)
    nicedepth[SELF] = ((nicedepth[SELF] == ERROR)
		       ? (nicedepth[i] + 1)
		       : MIN (nicedepth[SELF], nicedepth[i] + 1));

#ifndef NODYNAMIC
  /* If you are an internal non-barrier node and not currently
   * orphaned, take a moment to check the local situation and adjust
   * the hierarchy (by promoting one of your children) if need
   * be. What we're after is a more balanced or otherwise better
   * hierarchy. 
   *
   * Adjustments serve to ensure empty spaces higher up in the
   * hierarchy get filled. You know there's a space above you when
   * your parent's depth is less than (or equal to) your own.  Other
   * times, adjustment serves to help percolate better machines (in
   * terms of load or availability) up.  To avoid repeated
   * restructuring, we only adjust every RECONFODDS times you execute
   * niceLoad() (RECONFODDS = 13 and LOADINT = 9 means on average
   * about every 2 hours). 
   *
   * Since the normal situation would be for the parent's depth to be
   * one greater than mine, we increase the adjustment odds whenever
   * my parent has lower (or even equal) depth than I do. The bigger
   * the differences, the greater the odds.
   *
   * Note that electChild() identifies a promotable child if one can
   * be found, and that the criteria used to determine which child is
   * best for promotion depends on the circumstances of the promotion
   * itself. */
  if ((nicehcnt > 0) &&
      !(nicewall[SELF]) &&
      !(ADDRNULL (niceaddr[PARENT])) &&
      (mode != ORPHAN) &&
      (RANDINT (RECONFODDS) <
       MAX (1, (2 + nicedepth[SELF] - nicedepth[PARENT])))
      && (candidate = electChild () > 0))
    nicePromote (candidate, PARENT);
#endif

  /* Done. */
  return;
}

/**********************************************************************
 * serveNiceLoad: Service status request from your parent daemon.
 **********************************************************************/
void
serveNiceLoad (int socket, struct in_addr caddr)
{
  int i, status, wishes;
  char message[NICEMAXSTRINGLEN + 1];
#ifndef NODYNAMIC
  int child;
#endif

  /* Check if this is indeed your parent asking, then acknowledge. */
  if (ADDRNEQ (caddr, niceaddr[PARENT]))
    {
      if (debug &&
	  (errhp = gethostbyaddr ((char *) &caddr, ADDRBYTES, AF_INET)))
	niceLog (NLOGWRN, "Brushing off %s: my real parent is %s.\n",
		 errhp->h_name, inet_ntoa (niceaddr[PARENT]));
      else
	niceLog (NLOGWRN, "Brushing off %s: my real parent is %s.\n",
		 inet_ntoa (caddr), inet_ntoa (niceaddr[PARENT]));

      sockprintf (socket, message, NICEACKSTRD, FALSE);
      close (socket);
      return;
    }
  else if (sockprintf (socket, message, NICEACKSTRD, TRUE) == ERROR)
    {
      /* Lost link; give it up. */
      close (socket);
      return;
    }

  /* Get your parent's current depth and load information, along with
   * the hierarchy's current branching factor. The hierarchy's
   * branching factor is established by the root machine. */
  if (sockscanf (socket, message, NICEINFSTRDDD,
		 &nicedepth[PARENT], &niceload[PARENT], &nicebfact) != 3)
    {
      /* Lost link; give it up. */
      close (socket);
      return;
    }

  /* Now send your own depth and load information, along with the
   * length of your wishlist. */
  if (sockprintf (socket, message, NICEINFSTRDDD,
		  nicedepth[SELF], niceload[SELF],
		  (wishes = countPendingGets ())) == ERROR)
    {
      /* Lost link; give it up. */
      close (socket);
      return;
    }

  /* Send your parent your pending wishes, one by one. */
  i = 0;
  while (i < wishes)
    {
      /* Wait for acknowledgement. If it fails, just exit the inner
       * loop; since this is the last thing before you stop talking
       * with this parent, you needn't exit the outer loop explicitly,
       * nor close the socket or clear the alarm (these things happen
       * as you exit the outer loop). */
      if (sockscanf (socket, message, NICEACKSTRD, &status) != 1)
	break;

      /* Send your parent your ith pending get wish. This is woefully
       * inefficient, since each call to findNthPendingGet is O(N) in
       * the number of get wishes. But since these are limited, it's
       * probably not so bad. One pitfall is if we don't have an
       * appropriate number of pending gets (that is, somehow the get
       * wish list changes from the time we count to the time we
       * print). */
      if (sockprintf (socket, message, NICEINFSTRS,
		      findNthPendingGet (i)) == ERROR)
	break;
      i++;
    }

  /* Close the socket. */
  close (socket);

  /* Log output. */
  if (debug &&
      (errhp =
       gethostbyaddr ((char *) &niceaddr[PARENT], ADDRBYTES, AF_INET)))
    niceLog (NLOGMSG, "Serving status (%d) to %s\n", niceload[SELF],
	     errhp->h_name);
  else
    niceLog (NLOGMSG, "Serving status (%d) to %s\n",
	     niceload[SELF], inet_ntoa (niceaddr[PARENT]));

  /* Next adjust your maxChildren value to reflect the hierarchy's
   * current branching factor. Limited to MAXCHILDREN = MAXHOSTS-3,
   * since we need to reserve space for SELF, PARENT, and ROOT in the
   * host table. */
  if (debug && maxChildren != MIN (nicebfact, MAXCHILDREN))
    niceLog (NLOGMSG, "Resetting branching factor to %d.\n",
	     MIN (nicebfact, MAXCHILDREN));
  maxChildren = MIN (nicebfact, MAXCHILDREN);
  if (nicehcnt > maxChildren)
    {
      /* OK, we've suffered a branching factor restriction that will
       * require we jettison a few extraneous children (note that
       * nicehcnt must be <= maxChildren, and not strictly less than).
       *
       * If we're doing dynamic restructuring, we'll find the "best"
       * children to demote by considering barrier nodes, low overlap,
       * shallowness, and high load. Otherwise, we'll just drop a
       * bunch until we're in compliance. */
      while (nicehcnt > maxChildren)
	{
#ifndef NODYNAMIC
	  child = bestMatch (nicestart[SELF], niceend[SELF],
			     nicedepth[SELF], niceload[SELF],
			     FALSE, FALSE, FALSE, TRUE);
	  niceOrphan ((child > 0 && child <= nicehcnt) ? child : nicehcnt);
#else
	  niceOrphan (nicehcnt);
#endif
	}
    }

  /* Finished. */
  return;
}

/**********************************************************************
 * Request a new process be spawned on specified child machine.
 **********************************************************************/
void
niceSpawn (int child, char *executable, unsigned short int cport)
{
  int socket, status;
  char message[NICEMAXSTRINGLEN + 1];
  char appname[NICEMAXSTRINGLEN + 1];
  char filename[NICEMAXSTRINGLEN + 1];
#ifndef NOALARM
  int clientAlarm;
#endif
  if (debug &&
      (errhp = gethostbyaddr ((char *) &niceaddr[child], ADDRBYTES, AF_INET)))
    niceLog (NLOGMSG, "Spawning %s on %s\n", executable, errhp->h_name);
  else
    niceLog (NLOGMSG, "Spawning %s on %s\n",
	     executable, inet_ntoa (niceaddr[child]));

  /* Connect to the socket on the appropriate child. Here,
   * establishConnection () will return ERROR (no connection made),
   * FALSE (can't spawn any processes) or TRUE (go ahead). */
  status = establishConnection (child, NQIDSPAWN, message, &socket);
  if (status == ERROR)
    {
      /* No socket to close, give it up. */
      return;
    }
  else if (status == FALSE)
    {
      /* No help to be had. Close socket and exit. */
      close (socket);
      return;
    }

#ifndef NOALARM
  /* Set a timer in case something blocks during the execution of the
   * rest of the protocol. A protocol is declared hung/blocked and
   * attempts to read or write are abandoned when it hangs for HUNGINT
   * seconds, typically a large number. */
  if ((clientAlarm = niceSetTimeout (HUNGINT)) == FALSE)
    {
      /* The other end of the interaction has blocked. You shouldn't
       * need to clear the timer, which we know has already expired,
       * since that's how we got here. What we should do is close the
       * socket and return. */
      close (socket);
      return;
    }
#endif

  /* Send job information. Note we don't need to identify our host
   * address; only the port. The child will get the host address from
   * the network. Also, be sure to send only the executable name and
   * version, stripping off any architecture information that wouldn't
   * be valid on another machine anyway. */
  if (sockprintf (socket, message, NICEINFSTRSHU, executable, cport) == ERROR)
    {
      /* Lost link, give up. */
      close (socket);
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      return;
    }

  /* Wait for acknowledgement, then close socket, cleanup, and
   * exit. TRUE indicates the child is OK with the request. ERROR
   * indicates the request is in some way illegal (e.g., contains a
   * pathname), and FALSE indicates the child does not have the
   * required application.  */
  if (sockscanf (socket, message, NICEACKSTRD, &status) != 1)
    {
      /* Lost link, give up. */
      close (socket);
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      return;
    }

  /* If the child doesn't have the application but you do, make an
   * entry on the put list so that it can be sent on to the child at
   * the first opportunity. If you don't have the application, put it
   * on your own get list so that you can procure it to pass on to the
   * child. But before you can do all this, you need to get what the
   * child thinks is the full application name (including
   * architecture) so that you can obtain the right version. */
  if (status == FALSE)
    {
      /* Tell child to continue. */
      if (sockprintf (socket, message, NICEACKSTRD, TRUE) == ERROR)
	{
	  /* Lost link, give up. */
	  close (socket);
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (clientAlarm);
#endif
	  return;
	}

      /* Get application name. */
      if (sockscanf (socket, message, NICEINFSTRS, &appname) != 1)
	{
	  /* Lost link, give up. */
	  close (socket);
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (clientAlarm);
#endif
	  return;
	}

      /* Have you got the requested application? */
      if (niceCheckApplication (appname, filename))
	{
	  /* Found it! Add to your put list and respond TRUE. Your
	   * child will add the file to its own get list, marked for
	   * immediate download and automatic spawning on
	   * completion. */
	  niceLog (NLOGMSG, "Adding PUT request for %s.\n", appname);
	  addPutWish (appname, child);
	  sockprintf (socket, message, NICEACKSTRD, TRUE);
	}
      else
	{
	  /* No such luck. Add to your get list and respond
	   * FALSE. Your child will add the file to its own get list
	   * but without automatic spawning after download. */
	  niceLog (NLOGMSG, "Adding GET request for %s.\n", appname);
	  addGetWish (appname, NULL, NULL);
	  sockprintf (socket, message, NICEACKSTRD, FALSE);
	}
    }

  /* Close socket. */
  close (socket);
#ifndef NOALARM
  /* Clear the timer. */
  niceClearAlarm (clientAlarm);
#endif
  /* Finished. */
  return;
}

/**********************************************************************
 * serveNiceSpawn: services a request from parent daemon for a new
 * local application process, i.e., this is the procedure that spawns
 * a new job.
 **********************************************************************/
void
serveNiceSpawn (int socket, struct in_addr caddr, int avail)
{
  int i, status;
  char message[NICEMAXSTRINGLEN + 1];
  char executable[NICEMAXSTRINGLEN + 1];
  char filename[NICEMAXSTRINGLEN + 1];
  unsigned short int pport;
  char mport[NICEMAXSTRINGLEN + 1];

  /* First, we need to make sure that our status is OK to try to spawn
   * a new process locally. This includes checking to make sure its
   * really our parent asking. If it isn't, we punt. */
  if ((ADDRNEQ (caddr, niceaddr[PARENT])) ||
      (avail == FALSE) || (slaves == maxSlaves))
    {
      /* Reject request. Needn't check exit status, since we're just
       * going to return anyway. */
      sockprintf (socket, message, NICEACKSTRD, FALSE);
      niceLog (NLOGWRN, "Unable to spawn.\n");

      /* Punt. */
      close (socket);
      return;
    }

#if FALSE
  niceLog (NLOGMSG, "Attempting to spawn:\n");
  if (debug &&
      (errhp = gethostbyaddr ((char *) &niceaddr[0], ADDRBYTES, AF_INET)))
    fprintf (stderr, "+%s\n", errhp->h_name);
  else
    fprintf (stderr, "+%s\n", inet_ntoa (niceaddr[0]));
  if (debug &&
      (errhp = gethostbyaddr ((char *) &niceaddr[SELF], ADDRBYTES, AF_INET)))
    fprintf (stderr, " =%s\n", errhp->h_name);
  else
    fprintf (stderr, " =%s\n", inet_ntoa (niceaddr[SELF]));
  for (i = 1; i <= nicehcnt; i++)
    {
      if (debug &&
	  (errhp = gethostbyaddr ((char *) &niceaddr[i], ADDRBYTES, AF_INET)))
	fprintf (stderr, "  -%s\n", errhp->h_name);
      else
	fprintf (stderr, "  -%s\n", inet_ntoa (niceaddr[i]));
    }
#endif

  /* OK, still in the running. Acknowledge request. */
  if (sockprintf (socket, message, NICEACKSTRD, TRUE) == ERROR)
    {
      /* Lost link; give it up. */
      close (socket);
      return;
    }

  /* Get name of the executable to run. The exec name and port number
   * are important, but the hostname isn't, since I can get that from
   * caddr (besides, any self-reported network address might, in fact,
   * be wrong). 
   *
   * Here, the executable name will not include locally-dependent
   * components. So it will not have the full pathname or the
   * architecture information, but it should contain the executable
   * name and versioning information. */
  if (sockscanf (socket, message, NICEINFSTRSHU, executable, &pport) != 2)
    {
      /* Lost link; give it up. */
      close (socket);
      return;
    }

  /* Convert the port number into a string for execl. */
  snprintf (mport, NICEMAXSTRINGLEN, "%hu", pport);

  /* Scrutinize the application before launching, and ensure it is
   * locally available and properly signed. But before you do so, you
   * need to add the architecture designator, obtained locally, to the
   * exectuable string (which currently only contains the application
   * name and version information) after appending a separating
   * hyphen. */
  strncat (executable, "-", NICEMAXSTRINGLEN - strlen (executable));
  strncat (executable, NICEARCH, NICEMAXSTRINGLEN - strlen (executable));
  if ((status = niceCheckApplication (executable, filename)) == ERROR)
    {
      /* Illegal or malformed executable. Should never occur for
       * legitimate applications. */
      niceLog (NLOGWRN, "Illegal executable %s requested; fail.\n",
	       executable);

      /* No need to check exit status here, since we're just going
       * to return anyway. */
      sockprintf (socket, message, NICEACKSTRD, ERROR);
      close (socket);
      return;
    }
  else if (cardinality == 0)
    {
      /* Firewall machines with 0 cardinality are a problem. In the
       * future, might want to implement some sort of proxy forwarding
       * service to let applications within the firewalled network
       * interact with masters on the outside. */
      niceLog (NLOGWRN, "Can't fork %s; cardinality is zero.\n", executable);

      /* No need to check exit status here, since we're just going to
       * return anyway. We return TRUE, for now, because its not like
       * we're missing an application, we just can't launch it. */
      sockprintf (socket, message, NICEACKSTRD, TRUE);
      close (socket);
      return;
    }
  else if (status == FALSE)
    {
      /* No such executable. Add it to my wishlist and report failure. */
      if (sockprintf (socket, message, NICEACKSTRD, FALSE) == ERROR)
	{
	  /* Lost link; give it up. */
	  close (socket);
	  return;
	}

      /* Wait for ack. */
      if (sockscanf (socket, message, NICEACKSTRD, &status) != 1)
	{
	  /* Lost link; give it up. */
	  close (socket);
	  return;
	}

      /* Before sending the application name, we need to append the
       * local architecture, which is not in the executable string,
       * but is in the filename string, along with the path
       * structure. Just pass the application name. */
      if (sockprintf
	  (socket, message, NICEINFSTRS,
	   strrchr (filename, '/') + 1) == ERROR)
	{
	  /* Lost link; give it up. */
	  close (socket);
	  return;
	}

      /* Now depending on whether your parent has a copy available,
       * make an addition on the wishlist then exit. */
      if (sockscanf (socket, message, NICEACKSTRD, &status) != 1 ||
	  status == FALSE)
	{
	  /* No copy available; add to wishlist w/o launch
	   * option. When addr and port are both NULL, no attempt is
	   * made to launch the application after downloading. */
	  niceLog (NLOGMSG, "Adding GET request for %s.\n", (strrchr (filename, '/') + 1));
	  addGetWish ((strrchr (filename, '/') + 1), NULL, NULL);
	}
      else
	{
	  /* Parent has a copy; add to wishlist and mark for
	   * launch. When addr and port are set, application is
	   * immediately launched after download completes. */
	  niceLog (NLOGMSG, "Adding (live) GET request for %s.\n", (strrchr (filename, '/') + 1));
	  addGetWish ((strrchr (filename, '/') + 1), inet_ntoa (caddr),
		      mport);
	}
      /* Close socket and return. */
      close (socket);
      return;
    }

  /* All systems go. Acknowledge an attempt to spawn naggers. Whether
   * a nagger will actually be spawned depends on local conditions
   * (e.g., whether forking is successful). */
  if (sockprintf (socket, message, NICEACKSTRD, TRUE) == ERROR)
    {
      /* Lost link; give it up. */
      close (socket);
      return;
    }

  /* OK, committed. Close socket and fire up applications. */
  close (socket);

  /* Fire up. */
  i = spawnApplications (filename, 1, TRUE, caddr, mport);

  if (i == 0)
    niceLog (NLOGWRN, "Can't fork %s; fail.\n", filename);
  else
    niceLog (NLOGMSG, "Spawned %s %s %hu\n", executable, NICEFLAG, pport);
}

/**********************************************************************
 * Fire up a number of applications locally, depending on configured
 * cardinality of machine. This is also used in serveNiceSolicit()
 * when running naggers on a single machine (when cardinality is
 * greater than 1). 
 *
 * The filename name passed to spawnApplications() should have the
 * full pathname as well as versioning and architecture information
 * present explicitly in the string. spawnApplications() also assumes
 * that the executable is present and permitted (so check before
 * invoking).
 **********************************************************************/
int
spawnApplications (char *filename, int count, int fertile,
		   struct in_addr maddr, char *mport)
{
  int pid, i = 0;

  /* Spawn the appropriate number of slaves. */
  while (slaves < maxSlaves && i < count)
    {
      if ((pid = fork ()) < 0)
	/* Can't fork. Stop trying. */
	return (i);
      else if (pid == 0)
	{
#ifdef __CYGWIN__
	  /* Lower nagger's priority to match that of the daemon. For
	   * Unix, the nagger inherits its parent's priority, but
	   * that's not the case with Cygwin, so we need to set it
	   * explicitly. */
	  if (nice (priority) == ERROR)
	    niceLog (NLOGFAT, "Priority change failed; aborting.\n");
#endif
	  /* Are the copies spawned allowed to request new naggers
	   * themselves? */
	  if (fertile == TRUE)
	    execl (filename, filename, NICEFLAG, inet_ntoa (maddr), mport,
		   NULL);
	  else
	    execl (filename, filename, NICEFLAG, inet_ntoa (maddr), mport,
		   NICELEAF, NULL);

	  /* Should never execute. */
	  _exit (EXIT_FAILURE);
	}

      /* Update the slaves PID table. */
      spids[slaves] = pid;
      slaves++;
      /* Increment count of slaves spawned. */
      i++;
    }
  return (i);
}

/**********************************************************************
 * Fork a source process to transfer a copy of the specified
 * application to the designated child process. Invoked from
 * startPendingPuts (), which is in turn invoked at the end of main
 * loop. Note that startPendingPuts () will have already checked to
 * ensure it is OK to fork a new source process. Returns TRUE on
 * success, FALSE or ERROR on failure.
 **********************************************************************/
int
niceTxfr (PutWish *cell, char *filename, int length)
{
  int socket, status, spid;
  unsigned short int sport;
  char message[NICEMAXSTRINGLEN + 1];
#ifndef NOALARM
  int clientAlarm;
#endif

  if (debug &&
      (errhp = gethostbyaddr ((char *) &niceaddr[cell->child], ADDRBYTES, AF_INET)))
    niceLog (NLOGMSG, "Sending %s [%d] to %s\n", cell->name, length, errhp->h_name);
  else
    niceLog (NLOGMSG, "Sending %s [%d] to %s\n", cell->name, length, inet_ntoa (niceaddr[cell->child]));

  /* Connect to the socket on the appropriate child. Here,
   * establishConnection () will return FALSE/ERROR (no connection
   * made), or TRUE (go ahead). */
  status = establishConnection (cell->child, NQIDTXFR, message, &socket);
  if (status != TRUE)
    return (status);

#ifndef NOALARM
  /* Set a timer in case something blocks during the execution of the
   * rest of the protocol. A protocol is declared hung/blocked and
   * attempts to read or write are abandoned when it hangs for HUNGINT
   * seconds, typically a large number. */
  if ((clientAlarm = niceSetTimeout (HUNGINT)) == FALSE)
    {
      /* The other end of the interaction has blocked. You shouldn't
       * need to clear the timer, which we know has already expired,
       * since that's how we got here. What we should do is close the
       * socket and return. */
      close (socket);
      return (ERROR);
    }
#endif

  /* Next, send filename and length information to the child. */
  if (sockprintf (socket, message, NICEINFSTRSD, cell->name, length) == ERROR)
    {
      /* Punt. */
      close (socket);
#ifndef NOALARM
      niceClearAlarm (clientAlarm);
#endif
      return (ERROR);
    }

  /* Child will now open transfer stub socket and return the port
   * number to use for the transfer. If the port returned is 0 (or
   * message is not properly received), then that means the child
   * failed to get a socket; punt. */
  status = sockscanf (socket, message, NICEINFSTRHU, &sport);
  close (socket);
#ifndef NOALARM
  niceClearAlarm (clientAlarm);
#endif
  if ((status != 1) || (sport == 0))
    return (ERROR);

  /* All systems go. */
  niceLog (NLOGMSG, "Starting PUT for %s [%d] on %d.\n", cell->name, length, sport); 

#if TRUE
  /* All systems go; spawn a transfer stub process of your own,
   * retaining PID in putWish cell. */
  spid = fork ();
  if (spid < 0)
    return (ERROR);
  else if (spid == 0)
    {
      /* Source process. Clear the current niceTalk alarm in the
       * forked source process only. */
#ifndef NOALARM
      niceClearAlarm (clientAlarm);
#endif
      /* Push the file. When finished, just exit so you can be
       * culled. */
      niceApplicationSource (cell, sport, length);
      _exit (EXIT_SUCCESS);
    }

  /* Parent process after successful fork. Keep track of the source
   * pid and other transfer parameters so that it can be culled once
   * the file transfer is completed. */
  cell->pid = spid;
#else
  /* No forking for debugging. */
  niceApplicationSource (cell, sport, length);

  /* Mark putWish cell for cull. */
  cell->pid = ERROR;
#endif
  return (TRUE);
}

/**********************************************************************
 * Service a file transfer request. Forks a sink process that will
 * actually receive the file. Note that we pass the serverAlarm from
 * the main loop because it must be cleared in the forked subprocess!
 **********************************************************************/
#ifndef NOALARM
void
serveNiceTxfr (int socket, struct in_addr caddr, int avail, int alarm)
#else
void
serveNiceTxfr (int socket, struct in_addr caddr, int avail)
#endif
{
  int length, ssocket, spid;
  unsigned short int sport = port;
  GetWish *cell, *tmp;
  char message[NICEMAXSTRINGLEN + 1];
  char name[NICEMAXSTRINGLEN + 1];

  /* Check if this is indeed your parent asking, then acknowledge.
   * Also check if the host is available.
   *
   * TODO: Should check to make sure you don't exceed the number of
   * available transfer slots! */
  if (ADDRNEQ (caddr, niceaddr[PARENT]) || avail == FALSE)
    {
      /* Reject request. Needn't check exit status, since we're just
       * going to return anyway. */
      sockprintf (socket, message, NICEACKSTRD, FALSE);
      niceLog (NLOGWRN, "Unavailable for file transfer.\n");

      /* Punt. */
      close (socket);
      return;
    }

  /* OK, still in the running. Acknowledge request. */
  if (sockprintf (socket, message, NICEACKSTRD, TRUE) == ERROR)
    {
      /* Lost link; give it up. */
      close (socket);
      return;
    }

  /* Get name and length of application to download. This is also
   * where I would get, e.g., the signed MD5 hash of the file for
   * security. */
  if (sockscanf (socket, message, NICEINFSTRSD, name, &length) != 2)
    {
      /* Lost link; give it up. */
      close (socket);
      return;
    }

  /* Check for matching entry in your get list. If you can't find a
   * matching cell, you need to add one so that you can keep track of
   * the process parameters. This can happen when, e.g., a very old
   * get list item has fallen off your list (but not that of your
   * parent) and has finally become available. */
  if (!(tmp = (cell = findGetWish (name))))
    {
      niceLog (NLOGMSG, "Adding (unexpected) GET request for %s.\n", name);
      cell = addGetWish (name, NULL, NULL);
    }

  /* Next, we need to open a new socket for the sink process to listen
   * to. We need to do this before we fork off the sink process so
   * that we can report the port number to our parent. Just look for
   * an open port starting from the daemon's own port number. MAXPORT
   * is defined in nicecom.h as the largest unsigned short int. 
   * 
   * If you fail to get a socket, send a 0 port number to your parent
   * who will then know the transfer is off. Else, send the resulting
   * port number, then close the socket. */
  ssocket = watchSocket (&sport, MAXPORT);
  sockprintf (socket, message, NICEINFSTRHU, ((ssocket == ERROR)?0:sport));
  close (socket);
  if (ssocket == ERROR)
    return;

  /* All systems go. */
  niceLog (NLOGMSG, "Starting GET for %s [%d] on %d\n", name, length, sport);

#if TRUE
  /* Now you can fork your transfer stub. */
  spid = fork ();
  if (spid < 0)
    /* Fork failed. */
    return;
  else if (spid == 0)
    {
      /* Sink process. Clear the niced request-servicing alarm on the
       * forked sink process only. */
#ifndef NOALARM
      niceClearAlarm (alarm);
#endif
      /* Receive the file. When finished, just exit so you can be
       * culled. */
      niceApplicationSink (cell, ssocket, length);
      _exit (EXIT_SUCCESS);
    }

  /* Close the ssocket in the parent process; the forked stub will
   * close it later. */
  close (ssocket);

  /* OK, keep track of the sink process in the getWish cell. */
  cell->pid = spid;
#else
  /* No forking allowed for debugging. */
  niceApplicationSink (cell, ssocket, length);

  /* Mark getWish cell for cull. */
  cell->pid = ERROR;
#endif
  return;
}

/**********************************************************************
 * NICE control protocols.
 **********************************************************************/
/**********************************************************************
 * serveNiceQuery: Service request for a status report, including
 * daemon availability, load, and number of running slaves.  May
 * recursively query children and forward response. Invoked either by
 * niceq, the query program, or, recursively, by parent daemon.
 **********************************************************************/
void
serveNiceQuery (int psocket, int avail)
{
  int i, status = TRUE, level, csocket;
  int ccpus, cavail;
  char hname[NICEMAXSTRINGLEN + 1];	/* Host name */
  char hstat[NICEMAXSTRINGLEN + 1];	/* Host status */
  char message[NICEMAXSTRINGLEN + 1];

  /* Create a string describing your own availability, cardinality,
   * load and number of running slave processes. Use a plus sign to
   * indicate firewall nodes. Note that if we are orphaned, we'll
   * pretend we're root, at least temporarily. */
#if FALSE
  snprintf (hstat, NICEMAXSTRINGLEN, "[%s:%dx%d@%d.%d]%s",
	    ((avail == TRUE) ? "online" : "offline"),
	    cardinality, niceload[SELF], nicedepth[SELF], slaves,
	    ((ADDRNULL (niceaddr[PARENT])
	      || mode == ORPHAN) ? "*" : (nicewall[SELF] ? "+" : "")));
#else
  snprintf (hstat, NICEMAXSTRINGLEN, "[%s:%dx%d.%d]%s",
	    ((avail == TRUE) ? "online" : "offline"),
	    cardinality, niceload[SELF], slaves,
	    ((ADDRNULL (niceaddr[PARENT])
	      || mode == ORPHAN) ? "*" : (nicewall[SELF] ? "+" : "")));
#endif

  /* Send your name and status. By definition, your own level is
   * always 0. Note that if our name matches "localhost" it means we
   * don't really have a "real" network address, but probably just
   * some locally assigned string. While we'd really like to send our
   * IP address instead, its probably that niceaddr[SELF] is alo not a
   * meaningful address to the outside world (i.e., likely 127.0.0.1).
   * So what we will have to do is translate "localhost" to an IP
   * address on the receiving end of the nice query. */
  if (sockprintf (psocket, message, NICEQRYSTRDDDSS, 0,
		  cardinality, (avail && (slaves < maxSlaves)),
		  nicehost, hstat) == ERROR)
    {
      /* Lost link; give it up. */
      close (psocket);
      return;
    }

  /* Now wait for acknowledgement that will tell you if this is a
   * recursive query or not. */
  if ((sockscanf (psocket, message, NICEACKSTRD, &status) != 1) ||
      (status == FALSE))
    {
      /* Either you have a lost link, or you are not being queried
       * recursively. In either case, just close the socket and
       * exit. */
      close (psocket);
      return;
    }

  /* OK, this is a recursive query. Act as a filter passing up
   * information from your children. If a connection with your child
   * fails, report the child as unresponsive; if the connection with
   * your parent fails, give it up (no other choice, really). */
  i = 0;
  while (i < nicehcnt)
    {
      /* Preincrement i, since slaves are 1 indexed, but also because
       * we use continue extensively in this loop, so incrementing i
       * at the end would be problematic. 
       *
       * TODO: The problem is that the alarm set in main() closes
       * psocket, but not csocket. One solution might be to add a very
       * short alarm -- one that will expire BEFORE the pending alarm
       * from main() -- to close csocket and exit cleanly. But we
       * would have to make this alarm shorter than 1/nicehcnt the
       * length of the alarm in main(). */
      i++;
      if (hailSocket (i, port, &csocket) != TRUE)
	{
	  /* Child not answering (ERROR or FALSE from hailSocket), so
	   * inform your parent by crafting an artificial status
	   * string for the corresponding child, then go on to the
	   * next child. */
	  if (sockprintf (psocket, message, NICEQRYSTRDDDSS, 1, 0, FALSE,
			  inet_ntoa (niceaddr[i]), "[unresponsive]") == ERROR)
	    {
	      close (psocket);
	      return;
	    }
	  /* Wait for acknowledgement from parent. */
	  if (sockscanf (psocket, message, NICEACKSTRD, &status) != 1)
	    {
	      close (psocket);
	      return;
	    }
	  continue;
	}

      /* Got a socket to child, forward the status query. */
      if (sockprintf (csocket, message, NICEQRYSTRD, NQIDQUERY) == ERROR)
	{
	  /* Failed to send query; fake an answer for your parent and
	   * go on to next child. */
	  close (csocket);

	  /* Tell your parent child connection has evaporated. */
	  if (sockprintf (psocket, message, NICEQRYSTRDDDSS, 1, 0, FALSE,
			  inet_ntoa (niceaddr[i]), "[aborted]") == ERROR)
	    {
	      close (psocket);
	      return;
	    }
	  /* Wait for acknowledgement from parent. */
	  if (sockscanf (psocket, message, NICEACKSTRD, &status) != 1)
	    {
	      close (psocket);
	      return;
	    }

	  /* Try next child. */
	  continue;
	}

      /* Obtain information from the child. */
      if (sockscanf (csocket, message, NICEQRYSTRDDDSS, &level,
		     &ccpus, &cavail, hname, hstat) != 5)
	{
	  /* Failed to obtain information; fake an answer for your
	   * parent and go on to next child. */
	  close (csocket);

	  /* Tell your parent child connection has evaporated. */
	  if (sockprintf (psocket, message, NICEQRYSTRDDDSS, 1, 0, FALSE,
			  inet_ntoa (niceaddr[i]), "[aborted]") == ERROR)
	    {
	      close (psocket);
	      return;
	    }
	  /* Wait for acknowledgement from parent. */
	  if (sockscanf (psocket, message, NICEACKSTRD, &status) != 1)
	    {
	      close (psocket);
	      return;
	    }

	  /* Try next child. */
	  continue;
	}

      /* Keep reading from the current child and forwarding
       * information to your own parent until the child is done. */
      level = 0;
      while (level != ERROR)
	{
	  /* Send current node status to parent. If this fails, be
	   * sure to close the child socket that's still open before
	   * you give up. 
	   * 
	   * Careful: here's where you want to replace "localhost"
	   * received from your direct child with what you think is
	   * that child's IP address instead. TODO: Check how the 
	   * "abruptly orphaning" message discovers the hostname known
	   * by the outside world? */
	  if (sockprintf
	      (psocket, message, NICEQRYSTRDDDSS, (level + 1), ccpus, cavail,
	       ((strncmp (hname, "localhost", NICEMAXSTRINGLEN) == 0) ?
		inet_ntoa (niceaddr[i]) : hname), hstat) == ERROR)
	    {
	      close (csocket);
	      close (psocket);
	      return;
	    }
	  /* Wait for acknowledgement from parent. If this fails, be
	   * sure to close the child socket that's still open before
	   * you give up. */
	  if (sockscanf (psocket, message, NICEACKSTRD, &status) != 1)
	    {
	      close (csocket);
	      close (psocket);
	      return;
	    }

	  /* Tell child's server side you want to go on. */
	  if (sockprintf (csocket, message, NICEACKSTRD, TRUE) == ERROR)
	    {
	      /* Tell your parent child connection has evaporated. */
	      if (sockprintf (psocket, message, NICEQRYSTRDDDSS,
			      1, 0, FALSE, inet_ntoa (niceaddr[i]),
			      "[truncated]") == ERROR)
		{
		  close (psocket);
		  return;
		}
	      /* Wait for acknowledgement from parent. */
	      if (sockscanf (psocket, message, NICEACKSTRD, &status) != 1)
		{
		  close (psocket);
		  return;
		}
	      /* Try next child. */
	      break;
	    }
	  /* Obtain information from child. */
	  if (sockscanf
	      (csocket, message, NICEQRYSTRDDDSS, &level, &ccpus, &cavail,
	       hname, hstat) != 5)
	    {
	      /* Tell your parent child connection has evaporated. */
	      if (sockprintf (psocket, message, NICEQRYSTRDDDSS, 1, 0, FALSE,
			      inet_ntoa (niceaddr[i]),
			      "[truncated]") == ERROR)
		{
		  close (psocket);
		  return;
		}
	      /* Wait for acknowledgement from parent. */
	      if (sockscanf (psocket, message, NICEACKSTRD, &status) != 1)
		{
		  close (psocket);
		  return;
		}
	      /* Try next child. */
	      break;
	    }
	}
      /* Go on to next child. */
      close (csocket);
    }

  /* Send an end of protocol message. No need to check sockprintf
   * status because we're done anyway. */
  sockprintf (psocket, message, NICEQRYSTRDDDSS, ERROR, 0, FALSE,
	      NULLHOST, NULLHOST);

  /* Close the socket. */
  close (psocket);
  return;
}

/**********************************************************************
 **********************************************************************
 * NICE application protocols.
 **********************************************************************
 **********************************************************************/

/**********************************************************************
 * serveNicePing: Find out if a daemon is running and available.
 **********************************************************************/
void
serveNicePing (int socket, int avail)
{
  char message[NICEMAXSTRINGLEN + 1];

  /* Respond to the query with your current availability. Note there
   * is no need to check exit status, as this is all you're going to
   * do before you close the socket and exit anyway. */
  sockprintf (socket, message, NICEACKSTRD, avail);
  close (socket);
  return;
}

/**********************************************************************
 * serveNiceSolicit: Service request for new processes from a running
 * application via the API function niceInit(). Our policy is to spawn
 * (cardinality-1) naggers locally and then request naggers from each
 * child. This means to use all available machines, your application
 * has to be running on the same machine as the root daemon.
 **********************************************************************/
void
serveNiceSolicit (int socket, int avail)
{
  int i, status;
  char message[NICEMAXSTRINGLEN + 1];
  char executable[NICEMAXSTRINGLEN + 1];
  char filename[NICEMAXSTRINGLEN + 1];
  char sport[NICEMAXSTRINGLEN + 1];
  unsigned short int pport;
  char *archspec;

  if ((avail == FALSE) || (nicehcnt == 0 && cardinality == 1))
    {
      /* Reject the request. Either we're currently unavailable
       * because of time constraints, or we have no children on which
       * to spawn processes, and we are not allowed to spawn processes
       * locally (low cardinality). Note there is no need to check
       * sockprintf's exit status; even if the link were lost, we
       * would do the same thing anyway. */
      sockprintf (socket, message, NICEACKSTRD, FALSE);

      /* Close the socket. */
      close (socket);
      return;
    }

  /* Good to go; acknowledge, indicating the max number of naggers
   * that the soliciting process has a right to expect. We also
   * indicate whether or not we are a barrier node, so that the
   * application will know how to set the root of its own
   * descendents. */
  if (sockprintf (socket, message, NICEINFSTRDD, (nicehcnt + cardinality - 1),
		  nicewall[SELF]) == ERROR)
    {
      /* Link lost; give it up. */
      close (socket);
      return;
    }

  /* Get name of the executable to run. The name will include
   * versioning information, but no architecture or pathname
   * information. */
  if (sockscanf (socket, message, NICEINFSTRSHU, executable, &pport) != 2)
    {
      /* Link lost; give it up. */
      close (socket);
      return;
    }
  /* Close the socket. */
  close (socket);

  /* Scrutinize the application before launching, and ensure it is
   * locally available and properly signed. But before you do so, you
   * need to add the architecture designator, obtained locally, to the
   * exectuable string (which currently only contains the application
   * name and version information) after appending a separating
   * hyphen. */
  strncat (executable, "-", NICEMAXSTRINGLEN - strlen (executable));
  strncat (executable, NICEARCH, NICEMAXSTRINGLEN - strlen (executable));
  if ((status = niceCheckApplication (executable, filename)) == ERROR)
    {
      /* Someone tried to place a path structure. This should
       * never happen for a legitimate application. */
      niceLog (NLOGWRN, "Illegal executable %s requested; fail.\n",
	       executable);
      return;
    }
  else if (status == FALSE)
    {
      /* Can't find the executable. This should never happen because
       * the executable you are spawning is the same as the requesting
       * executable, which is already running on this machine! */
      niceLog (NLOGWRN, "Can't find executable %s; should never happen!\n",
	       executable);
      return;
    }

  /* Good to go. First, spawn any naggers you're allowed to spawn
   * locally. Here, we can use the executable name directly since we
   * know that executable is available on the local host. */
  if (cardinality > 1)
    {
      /* Fire up. Since you're the slave, you know that program will
       * be there, since we're running on the same machine. You just
       * need to construct strings with the name of the file and the
       * port number as a string. */
      snprintf (sport, NICEMAXSTRINGLEN, "%hu", pport);
      i = spawnApplications (filename, (cardinality - 1), FALSE,
			     niceaddr[SELF], sport);

      if (i == 0)
	niceLog (NLOGWRN, "Can't fork %s; fail.\n", filename);
      else
	niceLog (NLOGMSG, "[%d] %s %s %s %hu\n", slaves, executable,
		 NICEFLAG, inet_ntoa (niceaddr[SELF]), pport);
    }

  /* Careful: before you send the request to your children, strip the
   * architecture specification that you just added for the local
   * machine from the executable name before making your request by
   * overwriting the last hyphen with a NULL character, effectively
   * trimming the string. */
  archspec = strrchr (executable, '-');
  *archspec = 0;

#if FALSE
  if (nicehcnt > 0)
    {
      if (debug &&
	  (errhp = gethostbyaddr ((char *) &niceaddr[0], ADDRBYTES, AF_INET)))
	fprintf (stderr, "+%s\n", errhp->h_name);
      else
	fprintf (stderr, "+%s\n", inet_ntoa (niceaddr[0]));
      if (debug &&
	  (errhp =
	   gethostbyaddr ((char *) &niceaddr[SELF], ADDRBYTES, AF_INET)))
	fprintf (stderr, " =%s\n", errhp->h_name);
      else
	fprintf (stderr, " =%s\n", inet_ntoa (niceaddr[SELF]));
      for (i = 1; i <= nicehcnt; i++)
	{
	  if (debug &&
	      (errhp =
	       gethostbyaddr ((char *) &niceaddr[i], ADDRBYTES, AF_INET)))
	    fprintf (stderr, "  -%s\n", errhp->h_name);
	  else
	    fprintf (stderr, "  -%s\n", inet_ntoa (niceaddr[i]));
	}
    }
#endif

  /* Request new spawn on child daemons. */
  for (i = 1; i <= nicehcnt; i++)
    niceSpawn (i, executable, pport);

  /* Return. */
  return;
}

/**********************************************************************
 * serveNiceExit: Called when an application exits. Update your
 * internal tables. The primary raison d'etre here is to avoid leaving
 * dead applications around as zombies.
 **********************************************************************/
void
serveNiceExit (int socket)
{
  int i = 0;

  /* Close the socket. */
  close (socket);

  /* Our problem is that we don't know precisely which child has
   * exited, so we need to try all of them with an explicit WNOHANG
   * argument. Unfortunately, timing problems might conceivably leave
   * a zombie running when the application has already sent us an exit
   * signal but has yet to exit. Such zombies must be culled later. */
  while (i < slaves)
    if (waitpid (spids[i], NULL, WNOHANG) != spids[i])
      i++;
    else
      spids[i] = spids[--slaves];

  /* Return. */
  return;
}
