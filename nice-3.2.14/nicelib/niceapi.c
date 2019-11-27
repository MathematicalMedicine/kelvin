/**********************************************************************
 * NICE API Library
 * Alberto Maria Segre
 *
 * Copyright (c) 1997-2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software
 * for non-profit educational purposes only.
 **********************************************************************/
#include "niceaux.h"
#include "nicecom.h"
#include "niceapi.h"
#include "../config.h"
#include <unistd.h>		/* Needed for gethostname */
#include <stdio.h>		/* Needed for scanf, sscanf */
#include <stdlib.h>		/* EXIT_SUCCESS, EXIT_ERROR */
#include <string.h>		/* Needed for strncpy, strcmp */
#include <signal.h>		/* Needed for signal */
#include <limits.h>		/* USHRT_MAX */
#include <netinet/in.h>		/* Needed for in_addr */
#include <arpa/inet.h>		/* Needed for inet_aton */
#include <errno.h>		/* ETIMEDOUT, ENETUNREACH */

/**********************************************************************
 * A few global state variables. These are not extern!
 **********************************************************************/
/* The following variable is indexed by hid, like niceaddr[] from
 * nicecom.h. */
unsigned short int niceport[MAXHOSTS];	/* Parent, child (app) ports. */

/* Used for identification and authentication during protocols. */
int nicehid = -1;		/* My client id; not extern. */
int niceppid = 0;		/* My parent's pid; not extern. */
int nicerpid = 0;		/* My root's pid; not extern. */
int nicepid;			/* My own pid; not extern. */

#define MASTER		0x1	/* Status bits. Checked with &. */
#define SLAVE		0x2
#define NEITHER		0
int nicestatus = NEITHER;	/* Status; global, but not extern! */
int nicewall = FALSE;		/* Status; global, but not extern! */
#define LILENDIAN	0
#define BIGENDIAN	1
int niceendian;			/* Status; global, but not extern! */

/* Used to indicate we've been asked to exit. */
int die = FALSE;

/* Application query identifiers and message strings.  Define query
 * identifiers either here (if common to all applications) or in a
 * specific application. These query identifiers are used to index the
 * serverFn array directly, and therefore must be small ints, less
 * than NICEMAXQID (the largest query identifier supported). */
#define NQIDQUIT	-1
#define NQIDINIT	0

/* Used during handler definition and protocol execution. */
serverFn nicehandler[NICEMAXQID];	/* Vector of message handlers. */
int niceepoch[MAXHOSTS];	/* Vector of nag epoch indexes. */
#ifdef LIBDBG
int nicecalls[NICEMAXQID];	/* Vector of invocation counts. */
int nicewaits[NICEMAXQID]; 	/* Vector of deferred invocations. */
#endif

short int niceprops[NICEMAXQID];	/* Handler properties. */
#define EPOCHEQP	0x1		/* 00001 */
#define	EPOCHINC 	0x2		/* 00010 */
#define ROOTONLY	0x4		/* 00100 */
#define SAFETY		0x8		/* 01000 */
#define SYNCHPTCL	0x10		/* 10000 */

/* Used in dispatching signals: the server will elect one of the
 * following four options. */
#define NICEREFUSE -1
#define NICEACCEPT  0
#define NICEIGNORE  1
#define NICEDEFER   2

#ifndef NOALARM
#ifndef NODEFER
/* Used to ensure nicesocket is serviced at regular intervals. */
#define DEFERINT 30	/* Must be <75sec, typical TCP/IP kernel timeout. */
int deferAlarm;		/* Defer timer ID: retain to cancel. */
#endif
#endif

/* Used to queue up signals, or deferred asynchronous (but unsafe)
 * messages. Deferred synchronous unsafe messages are simply rebuffed, and
 * resubmitted by the client. But you can't do that with asynchronous
 * messages, since the client isn't waiting around to hear back from
 * you: hence these signals are cached internally on the server side, and 
 * executed later when the context is safe. */
typedef struct signalCell
{
  int qid;
  int qhid;
  int qpid;
  int qepoch;
  struct signalCell *next;
}
SignalCell;
/* Linked list of signals (i.e., deferred) asynchronous messages. */
SignalCell *signals = NULL;
SignalCell *freeSignals = NULL;

/* Function prototypes which needn't be known to the outside world. */
int niceEndian ();
void niceQuit (int signum);
#if FALSE
/* DEPRECATED: use SIGSTOP and SIGCONT instead. */
void niceSleep (int signum);
void niceAwake (int signum);
#endif
int niceSolicit (char *executable);
int nicePing (int hid);
void niceDefer ();
void niceServiceSocket (int defer);

/**********************************************************************
 * Set up an application to interface with NICE. First, parse any
 * input arguments, searching for NICE-specific information. Then make
 * sure the appropriate signal-handlers are set. Finally, execute the
 * initialization protocol (if you are a slave), and set up your own
 * socket to monitor. Returns ERROR or an upper bound on the number of
 * naggers the application can expect to hear from.
 **********************************************************************/
int
niceInit (int argc, char *argv[],
	  clientFn clientFn,	/* Sends init info. */
	  serverFn serverFn)	/* Gets init info. */
{
  int i, solicit = TRUE, status;
  char address[NICEMAXSTRINGLEN + 1];
  char executable[NICEMAXSTRINGLEN + 1];	/* My program name. */
  struct sigaction sigact; 			/* Signals. */

#ifdef LIBDBG
  fprintf (stderr, "NiceInit: ");
#endif
  /* First, determine name of our host machine and punt if you can't. */
  if ((gethostname (nicehost, NICEMAXSTRINGLEN)) == ERROR)
    {
#ifdef LIBDBG
      fprintf (stderr, "Hostname error; NICE unavailable.\n");
#endif
      return (ERROR);
    }
  /* Make sure we have the fully-specified name (including domain) and
   * not just a short name as is returned by gethostname on some Unix
   * systems, otherwise punt. While you're at it, get your address and
   * address length and cache them away for future use when making
   * socket connections. 
   *
   * Essentially, we are establishing a value for niceaddr[SELF] by
   * asking the machine what it thinks its IP address is. This is
   * different than the other values of niceaddr[], which are obtained
   * directly from the network layer, and reflect the fact that some
   * machine may have a different IP addresses which present
   * differently to other machines on the network because of routing
   * issues.
   *
   * Fortunately, we don't really ever use the niceaddr[SELF] value,
   * except for debugging information. Probably should just set it
   * explicitly to localaddr. */
  if ((qualifyHostname (nicehost, SELF)) == ERROR)
    {
#ifdef LIBDBG
      fprintf (stderr, "Hostname error; NICE unavailable.\n");
#endif
      return (ERROR);
    }
  /* Before you do anything else, find out if you have a NICE daemon
   * running and available on the local machine and punt if you
   * don't. */
  if (!nicePing (SELF))
    {
#ifdef LIBDBG
      fprintf (stderr, "Local NICE daemon unavailable.\n");
#endif
      return (ERROR);
    }

  /* Parse the NICE-specific command line arguments that follow the
   * NICEFLAG.  Look for niceaddr[PARENT] and niceport[PARENT]
   * (explicitly placed in that order by NICEAPI, along with an
   * optional NICELEAF designation for "extra" naggers launched on an
   * SMP machine). All things being equal, we would like to use the
   * same port as our parent, but if we are the root process we'll
   * instead start looking for a free port number one past the niced
   * port. 
   *
   * Input length limits derived from SHRT_MAX and INT_MAX. 
   *
   * TODO: Note that if someone tries to supply NICE-specific command
   * line arguments to the top-level application it will probably
   * eventually fail, since there won't be a copy of the application
   * running on the user-specified host and port. However, it would be
   * better to explicitly detect this situation and exit. */
  niceport[SELF] = NICEPORT + 1;
  for (i = 0; i < argc; i++)
    {
      if (strcmp (argv[i], NICEFLAG) == 0)
	{
	  nicestatus = nicestatus | SLAVE;
	  sscanf (argv[++i], NICESAFEARGS, address);
	  inet_aton (address, &niceaddr[PARENT]);
	  sscanf (argv[++i], NICESAFEARGHD, &niceport[PARENT]);
	  niceport[SELF] = niceport[PARENT];
	  if (i + 1 < argc && strcmp (argv[i+1], NICELEAF) == 0)
	    {
	      solicit = FALSE;
	      i++;
	    }
#if FALSE
	  fprintf (stderr, "Host address %s port %hu\n",
		   address, niceport[PARENT]);
#endif
	}
    }
  /* Get your process id number to use for authentication during
   * exchanges with your child processes. */
  nicepid = getpid ();

  /* Figure out what kind of endian we are (either LILENDIAN or
   * BIGENDIAN). */
  niceendian = niceEndian();

  /* If you are the root application process, then copy your host
   * information from nicepid, niceaddr[SELF] and niceport[SELF] into
   * nicerpid, niceaddr[ROOT] and niceport[ROOT]. These are the values
   * that will get handed down to all naggers deep in the hierarchy
   * (except for those hidden behind "barrier" nodes which will
   * instead record the barrier's address and port number). 
   *
   * This is the one place where we'd like to use niceaddr[SELF],
   * which gets copied to niceaddr[ROOT]. To use it safely, we had
   * better be darned sure this is a well-accessible IP
   * address. Unfortunately, this is not a safe bet; niceaddr[SELF]
   * might well be localaddr, which would be useless
   * elsewhere. Probably should simply set it to localaddr and stop
   * pretending!
   *
   * Note that we would like to check niceRoot() here, but we can't
   * because MASTER won't be set until the end of this routing. We'll
   * have to settle for recognizing the fact that we were not born a
   * slave to indicate we are intended to be the root. */
  if (!niceSlave ())
    {
      nicerpid = nicepid;
      niceaddr[ROOT] = niceaddr[SELF];
      niceport[ROOT] = niceport[SELF];
    }

  /* Register the initialization function's message handler. This is a
   * synchronous protocol, unaffected by epoch number. You have to do
   * this before you try to contact your master using niceTalk, or the
   * epoch numbers will be out of synch! Important: NQIDINIT's
   * serverFn is always assumed "unsafe" meaning that it will be
   * deferred if invoked from, e.g., a timer interrupt. */
  for (i = 0; i < NICEMAXQID; i++)
    nicehandler[i] = NULL;
  niceHandler (NQIDINIT, serverFn, TRUE, FALSE, FALSE, FALSE, FALSE);
  /* Initialize your parent's epoch count. */
  niceepoch[PARENT] = 0;

  /* Now set up a port for monitoring, and punt on failure.  Start at
   * same port used by your parent, if specified, and sweep up to
   * MAXPORT, the largest possible unsigned short int port
   * number. Somewhere we should find an unused port number; if we do,
   * then watchSocket will ensure niceport[SELF] is appropriately set
   * for the rest of the API to access. */
  if ((nicesock = (watchSocket (&niceport[SELF], MAXPORT))) == ERROR)
    {
      /* It bombed. You can't get an appropriate port to monitor, so
       * you may as well give up right away. */
#ifdef LIBDBG
      fprintf (stderr, "Can't get a port.\n");
#endif
      niceExit ();
    }
#ifdef LIBDBG
  fprintf (stderr, "[port %d] ", niceport[SELF]);
#endif

  /* If you are root, you need to update your ROOT port because it
   * might have changed while scanning for an open port to monitor. */
  if (!niceSlave())
    niceport[ROOT] = niceport[SELF];

  /* Next, we need make sure application processes trap and handle
   * most signals appropriately, closing the socket we just set up to
   * monitor. Here, we worry about the most common external signals.
   * Recall SIGKILL and SIGSTOP behavior can't be changed. TODO:
   * probably should thoroughly check all signals to see if this set
   * is reasonable. */
  sigact.sa_handler = niceQuit;
  sigemptyset (&sigact.sa_mask);
  sigact.sa_flags = 0;
  sigaction (SIGHUP, &sigact, NULL);
  sigaction (SIGINT, &sigact, NULL);
  sigaction (SIGQUIT, &sigact, NULL);
  sigaction (SIGTERM, &sigact, NULL);

#if FALSE
  /* The last signals we need to trap are the signals sent to an
   * application by the NICE daemon to tell it to go to sleep and to
   * wake up. To initialize this behavior properly, we can simply call
   * niceAwake which sets the signal handlers in the proper fashion
   * (i.e., implicitly assumes that a fresh slave is awake).
   *
   * DEPRECATED: use SIGSTOP and SIGCONT instead. */
  if (niceSlave ())
    niceAwake (0);
#endif

  /* Set up the program name in the program variable. Each program
   * should, ideally, initially have three hyphen-separated parts,
   * appname-version-arch, indicating the name of the application,
   * the application version (distinct from the NICE version), and the
   * machine architecture for which the executable is compiled.
   *
   * We can get this from argv[0] by flushing any directory
   * information. On the root processor, we can't guarantee that the
   * program is actually not invoked via, e.g., a soft link which
   * would have the effect of supplanting the real three-part name in
   * argv[0]. In this case, we need to use readlink () to resolve the
   * name before continuing; note the use of !niceSlave () since we
   * won't know we're the root quite yet.
   * 
   * All this effort is for security and versioning; we want to fire
   * up the matching version of the program but only from the approved
   * NICE BINDIR. Also, note that no application really has an
   * intrinsic notion of "self" compiled in -- it's all driven from
   * executable file name, which is up to the root processor admin to
   * maintain (easier, if Makefiles are appropriate) and up to the
   * Nice system itself as it propagates exectuables down the
   * hierarchy.  */
  if (!niceSlave () && ((i = readlink (argv[0], executable, NICEMAXSTRINGLEN)) > 0))
    {
      /* Resolve any soft links: copy real filename back into
       * argv[0]. Careful: readlink does not NULL terminate the
       * string. */
      executable[i] = (char) NULL;
      strncpy (argv[0], executable, NICEMAXSTRINGLEN);
    }
  if (strrchr (argv[0], '/'))
    /* Flush directory structure when copying from argv[0]. */
    strncpy (executable, (strrchr (argv[0], '/') + 1), NICEMAXSTRINGLEN);
  else
    /* No directory structure specified; copy directly from argv[0]. */
    strncpy (executable, argv[0], NICEMAXSTRINGLEN);

  /* Now explicitly blot out the architecture specification, which
   * follows the last hyphen (inclusive), leaving the version
   * information intact. Note that this is guaranteed to work on both
   * root and non-root (i.e., slave) machines, since we took the
   * trouble to resolve any soft links on the root machine. */
  executable[strrchr (executable, '-') - &(executable[0])] = 0;

  /* Report to your master if you are a nagger by exchanging the
   * authentication and optional problem-specific initialization
   * information. The latter corresponds to executing a protocol,
   * where the two sides of the protocol are specified by clientFn and
   * serverFn and the identifying QID is NQIDINIT. Need to do this
   * after you establish the nicesock so that niceport[SELF] is
   * properly reported to your parent. */
  if (niceSlave ())
    {
      /* A status of TRUE means the protocol has gone off as planned,
       * while FALSE means a nonfatal error has been encountered
       * (e.g., the parent is still there but too busy setting itself
       * up to respond yet). Alternatively if status is ERROR, then
       * the parent will never responed; perhaps this slave was
       * started by a master who has since expired.
       *
       * If you get FALSE, then you need to keep trying. The trick is
       * to always be in your master's queue so that when he is
       * eventually ready to service slaves you'll be waiting. */
#ifdef LIBDBG
      fprintf (stderr, "initializing...");
#endif
      while ((status = niceTalk (PARENT, NQIDINIT, clientFn)) == FALSE)
	{
#ifdef LIBDBG
	  fprintf (stderr, "retrying...");
#endif
	}

      /* It bombed.  Give it up. */
      if (status == ERROR)
	{
#ifdef LIBDBG
	  fprintf (stderr, "aborting.\n");
#endif
	  niceExit ();
	}
    }

  /* Unless you are a designated leaf process (happens when more than
   * one nagger is spawned on a machine) request help from your NICE
   * daemon, asking that volunteers register at your port. Return an
   * upper bound on how many volunteers one can expect. */
  if (solicit && (i = niceSolicit (executable)) > 0)
    nicestatus = nicestatus | MASTER;
  else
    i = 0;

#ifdef LIBDBG
  fprintf (stderr, "done.\n");
#endif
#ifndef NOALARM
#ifndef NODEFER
  /* One last thing to do before you return: establish a timer that
   * invokes niceDefer() in case you don't invoke niceCheck() often
   * enough. */
  deferAlarm = niceSetInterrupt (DEFERINT, niceDefer, FALSE);
#endif
#endif
  return (i);
}

/**********************************************************************
 * Tell a nagger to exit gracefully. This is the function used to trap
 * a quit signal (sent by a parent application, generated by the local
 * NICE daemon, or from the keyboard for the root process). It's up to
 * niceCheck to check this variable and invoke niceExit.
 **********************************************************************/
void
niceQuit (int signum)
{
  die = TRUE;
}

#if FALSE
/**********************************************************************
 * Take a nagger offline. This function is invoked in response to a
 * NICESLEEP signal, which is usually generated by the NICE daemon
 * running on the local machine.
 *
 * DEPRECATED: use SIGSTOP and SIGCONT instead.
 **********************************************************************/
void
niceSleep (int signum)
{
  /* Set signal handlers appropriately. */
  signal (NICESLEEP, SIG_IGN);
  signal (NICEAWAKE, niceAwake);

  /* Go to sleep until you get are NICEAWAKEned. */
  pause ();
}

/**********************************************************************
 * Wake a sleeping nagger up. This function is invoked in response to
 * a NICEAWAKE signal, which is usually generated by the NICE daemon
 * running on the local machine.
 *
 * DEPRECATED: use SIGSTOP and SIGCONT instead.
 **********************************************************************/
void
niceAwake (int signum)
{
  /* Set signal handlers appropriately. */
  signal (NICESLEEP, niceSleep);
  signal (NICEAWAKE, SIG_IGN);
}
#endif

/**********************************************************************
 * Close up connection with the nice environment. Involves some
 * various cleanup activities, followed by sending a message to your
 * local NICE daemon to tell him you're quitting.  No process ever
 * returns from this call.
 **********************************************************************/
void
niceExit ()
{
  int socket, i;
  char message[NICEMAXSTRINGLEN + 1];

#ifdef LIBDBG
  fprintf (stderr, "Exiting.\n");
  if (niceSlave())
    {
      fprintf (stderr, "niceTalk: ");
      for (i = 0; i < NICEMAXQID; i++)
	fprintf (stderr, "[%d:%d/%d] ", i, nicecalls[i] - nicewaits[i], nicecalls[i]);
      fprintf (stderr, "\n");
    }
#endif

  /* Tell your children you're quitting. */
  if (niceMaster ())
    for (i = 1; i <= nicehcnt; i++)
      niceTalk (i, NQIDQUIT, NULL);

  /* Close the socket we have been monitoring. Those children that
   * haven't gotten your message are probably waiting to talk to you,
   * in which case closing the socket should cause them to notice and
   * then die on their own. */
  if (niceParallel ())
    close (nicesock);

  /* Tell your daemon you are about to exit and then do it. If you
   * can't connect to your local NICE daemon, go ahead and die
   * directly. This will be the case when your daemon is killed and
   * sends you a SIGQUIT, invoking niceExit on your end.  While its
   * important to have niceExit recursively kill my naggers (because
   * when I expire they can't possibly do useful work), we shouldn't
   * be too concerned if we can't report back to the daemon; the
   * daemon will reap this zombie process eventually anyway. */
  if (niceSlave () && hailSocket (SELF, NICEPORT, &socket) == TRUE)
    {
      /* Prod the daemon to alert him that one of his slaves is
       * exiting, but don't wait around for an answer. */
      sockprintf (socket, message, NICEQRYSTRD, NQIDEXIT);

      /* Close the socket. */
      close (socket);
    }
#if FALSE
  else if (niceSlave ())
    fprintf (stderr, "Failed to alert my daemon of my demise...\n");
#endif

  /* Make sure you exit so your daemon (if you are nagger), who is now
   * waiting for you to die, can go on with its business. */
  exit (EXIT_SUCCESS);
}

/**********************************************************************
 * NICE status checking functions. These all check the value of the
 * nicestatus variable, and should be compiled inline.
 **********************************************************************/

/**********************************************************************
 * Return TRUE if you are a slave.
 **********************************************************************/
int
niceSlave ()
{
  return (nicestatus & SLAVE);
}

/**********************************************************************
 * Return TRUE if you are a master.
 **********************************************************************/
int
niceMaster ()
{
  return (nicestatus & MASTER);
}

/**********************************************************************
 * Return TRUE if you are a leaf process.
 **********************************************************************/
int
niceLeaf ()
{
  return (nicestatus == SLAVE);
}

/**********************************************************************
 * Return TRUE if you are the root process.
 **********************************************************************/
int
niceRoot ()
{
  return (nicestatus == MASTER);
}

/**********************************************************************
 * Return TRUE if running in parallel.
 **********************************************************************/
int
niceParallel ()
{
  return (nicestatus);
}

/**********************************************************************
 * Return TRUE if running standalone (serial). Note this will, by
 * design, also return TRUE if niceInit() was never invoked.
 **********************************************************************/
int
niceSerial ()
{
  return (!nicestatus);
}

/**********************************************************************
 * Endianess checking utility.
 **********************************************************************/
int
niceEndian ()
{
  short int word = 0x0001;
  char *byte = (char *) &word;
  return (byte[0] ? LILENDIAN : BIGENDIAN);
}

/**********************************************************************
 * NICE utilities.
 **********************************************************************/

/**********************************************************************
 * Called by application (client) to connect to designated host
 * (server) and initiate an exchange. Server is specified by the hid,
 * and the qid specifies the protocol to initiate. NQIDQUIT is handled
 * differently; all other qids are handled uniformly and according to
 * their properties.  Once the socket is successfully established and
 * the server acknowledges the query, the clientFn is invoked to
 * handle the rest of the protocol, with any additional arguments
 * passed to niceTalk forwarded on to the clientFn.
 *
 * Returns TRUE (positive integer) if successful, ERROR in the event
 * of an unrecoverable failure, or FALSE if client should reattempt.
 *
 * Important: every clientFn must also adhere to these return value
 * semantics, returning either TRUE for success or FALSE for deferral.
 * Also, when deferring, the client is advised to retry immediately;
 * thus, any use of timeout should be handled explicitly by the
 * application.
 **********************************************************************/
int
niceTalk (int hid, int qid, clientFn clientFn, ...)
{
  int socket, status, swap;
  char address[NICEMAXSTRINGLEN + 1];
  char message[NICEMAXSTRINGLEN + 1];
  va_list argp;
#ifndef NOALARM
  int clientAlarm;
#endif

  /* Before you even start, make sure that a client function is
   * provided only when appropriate. */
  if (((niceprops[qid] & SYNCHPTCL) && !clientFn) || 	/* Synch, no client! */
      (!(niceprops[qid] & SYNCHPTCL) && clientFn))	/* Asynch, w/client! */
    {
#ifdef LIBDBG
      fprintf (stderr, "Ill-formed client %d; aborting.\n", qid);
#endif
      return (ERROR);
    }

  /* Establish socket. If hailSocket fails, return the status of
   * hailSocket. This value can be ERROR (in the case that the failure
   * is fatal) or FALSE (in the case, say, that the target host is
   * busy, yielding a connection timeout). */
  if ((status = hailSocket (hid, niceport[hid], &socket)) != TRUE)
    {
#ifdef LIBDBG
      fprintf (stderr, "Hail socket failure (%d); aborting.\n", status);
#endif
      return (status);
    }

#ifndef NOALARM
  /* Set a timer in case something blocks during the execution of the
   * protocol. A protocol is declared hung/blocked and attempts to
   * read or write are abandoned when it hangs for HUNGINT seconds,
   * typically a large number. TODO: This should probably NOT be a
   * constant, but rather should be set based on some factor of the
   * observed RTT. */
  if ((clientAlarm = niceSetTimeout (HUNGINT)) == FALSE)
    {
      /* The other end of the interaction has blocked. You shouldn't
       * need to clear the timer, which we know has already expired,
       * since that's how we got here. */
#ifdef LIBDBG
      fprintf (stderr, "Blocked client [%d] protocol; aborting.\n", qid);
#endif
      close (socket);
      return (ERROR);
    }
#endif

  /* Send a query identifier along with endianess and authentication
   * information.
   *
   * If you are a child contacting your parent, send your nicehid
   * (obtained from your parent at initial contact) and niceppid (your
   * parent's pid) for comparison with the parent's pid.
   * 
   * If you are a parent contacting your child, send simply your own
   * pid, which your child can then compare to his copy.
   *
   * If you are an arbitrary node contacting the root, send ERROR as
   * the hid (since you won't be in the root's tables anyway) and
   * nicerpid (your root's pid, which, if you are behind a barrier
   * node, may not actually be the global root) which could be used
   * comparison with the server's pid, but would not offer much
   * assurance from an authentication point of view so we don't really
   * bother.
   *
   * Finally, note that you also need to send your epoch number for
   * synchronization purposes.  Note that, initially, niceppid=0 and
   * nicehid=-1, So if qid is NQIDINIT, these values will be
   * meaningless until your parent tells you (during the initial
   * exchange) what they are. */
  if (sockprintf (socket, message, NICEQRYSTRDDDDD, qid,
		  ((hid==PARENT)?nicehid:((hid==ROOT)?ERROR:PARENT)),
		  ((hid==PARENT)?niceppid:((hid==ROOT)?nicerpid:nicepid)),
		  niceepoch[hid], niceendian) == ERROR)
	{
#ifdef LIBDBG
	  fprintf (stderr, "niceTalk: write error [%s].\n", message);
#endif
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (clientAlarm);
#endif
	  close (socket);
	  return (ERROR);
	}
#ifdef LIBDBG
  /* Keep some statistics in debug mode. */
  if (qid >= 0 && qid < NICEMAXQID)
    nicecalls[qid]++;

  fprintf (stderr, "<-Sent QID %d to %d as %d",
	   qid, hid,
	   ((hid==PARENT)?nicehid:((hid==ROOT)?ERROR:PARENT)));
  if (qid != NQIDQUIT && niceprops[qid] & (EPOCHEQP | EPOCHINC))
    fprintf (stderr, " [epoch %d]", niceepoch[hid]);
#endif

  if ((qid != NQIDQUIT) && (niceprops[qid] & SYNCHPTCL))
    {
      /* This is a synchronous protocol. Wait for acknowledgement. At
       * this point, the server might elect to accept your credentials
       * (status NICEACCEPT), defer service (status NICEDEFER), reject
       * your credentials (status NICEREFUSE), or ignore your message
       * (status NICEIGNORE). There are lots of reasons to reject
       * credentials, but the action on this end should always be just
       * to exit with an ERROR status. In contrast, messages are
       * ignored when their epochs don't match. Not that this is also
       * where the server tells us if endianess-related byteswapping
       * should eventually occur along this connection. */
      if ((sockscanf (socket, message, NICEACKSTRDD, &status, &swap) != 2) ||
	  (status == NICEREFUSE))
	{
#ifdef LIBDBG
	  fprintf (stderr, ": credential error.\n");
#endif
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (clientAlarm);
#endif
	  close (socket);
	  return (ERROR);
	}
      else if (status == NICEDEFER)
	{
#ifdef LIBDBG
	  fprintf (stderr, ": service deferred.\n");
	  /* Keep statistics. */
	  if (qid >= 0 && qid < NICEMAXQID)
	    nicewaits[qid]++;
#endif
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (clientAlarm);
#endif
	  close (socket);
	  return (FALSE);
	}
      else if (status == NICEIGNORE)
	{
#ifdef LIBDBG
	  fprintf (stderr, ": ignored.\n");
#endif
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (clientAlarm);
#endif
	  close (socket);
	  /* Pretend everything worked smoothly! */
	  return (TRUE);
	}
#ifdef LIBDBG
      else if (status != NICEACCEPT)
	fprintf (stderr, " [something fishy!]");
#endif

      /* OK. You've been accepted by the server. If this is my first
       * interaction with my master, I need to tell him my port
       * information so he can make an appropriate entry in his host
       * table, returning the corresponding index value in his host
       * table, which I will use from now on to identify myself to my
       * parent. Note that he already has my IP address information
       * from the network.
       *
       * We might want to improve the authentication stuff a bit
       * later. As it is, we only get the network layer's IP address
       * on first contact, rather than checking it on each subsequent
       * contact. */
      if (qid == NQIDINIT)
	{
#ifdef LIBDBG
	  fprintf (stderr, " [new]");
#endif
	  if (sockprintf (socket, message, NICEINFSTRHU, niceport[SELF]) == ERROR)
	    {
#ifdef LIBDBG
	      fprintf (stderr, ": write error [%s].\n", message);
#endif
#ifndef NOALARM
	      /* Clear the timer. */
	      niceClearAlarm (clientAlarm);
#endif
	      close (socket);
	      return (ERROR);
	    }

	  /* Obtain the master's PID and the root's PID, which I can
	   * use for authentication in the future.  If the PID and
	   * index value returned are less than 0, it means the master
	   * failed to resolve the hostname I provided. Also, obtain
	   * from the master the address and port number of the root
	   * application host, which may be used in some protocols
	   * where a deeply nested nagger needs to report directly to
	   * the root processor.  Then procede with whatever protocol
	   * I originally intended. */
	  if (sockscanf (socket, message, NICEINFSTRSHUDDD, 
			 address, &niceport[ROOT], &nicehid, 
			 &niceppid, &nicerpid) != 5)
	    {
#ifdef LIBDBG
	      fprintf (stderr, ": authentication error [%s].\n", message);
#endif
#ifndef NOALARM
	      /* Clear the timer. */
	      niceClearAlarm (clientAlarm);
#endif
	      close (socket);
	      return (ERROR);
	    }

	  /* Obtain the root address. If your parent sent LOCALHOST,
	   * it means your parent is the root, so cache its IP
	   * address, from the network, in niceaddr[ROOT]. Anything
	   * other than LOCALHOST should be interpreted as the true
	   * root's address. The reason we have to go through this
	   * charade is because it is hard for a machine to truly know
	   * its own IP address, since that depends on proper
	   * configuration. */
	  if (strncmp (LOCALHOST, address, 9) == 0)
	    memcpy (&niceaddr[ROOT], &niceaddr[hid], sizeof (niceaddr[hid]));
	  else
	    inet_aton (address, &niceaddr[ROOT]);
#ifdef LIBDBG
	  fprintf (stderr, " [parent %s:%d][root %s:%d]",
		   inet_ntoa(niceaddr[PARENT]), niceport[PARENT],
		   inet_ntoa(niceaddr[ROOT]), niceport[ROOT]);
#endif
	}
    }

  /* Commit: at this point, authentication is OK, server has accepted
   * the message, and the client side function is correctly
   * specified. Increment the client epoch number if the protocol
   * requires. */
  status = TRUE;
  if (qid != NQIDQUIT && (niceprops[qid] & EPOCHINC))
    niceepoch[hid]++;

  /* If not an asynch protocol, invoke the clientFn (note that
   * NQIDQUIT is implicitly an async protocol). Important: client must
   * return TRUE (positive integer) for success or FALSE to indicate
   * another attempt should be undertaken immediately. */
  if ((qid != NQIDQUIT) && (niceprops[qid] & SYNCHPTCL))
    {
      /* Initialize the remaining arguments for the clientFn. */
      va_start (argp, clientFn);
#ifdef LIBDBG
      fprintf (stderr, ": dispatch.\n");
#endif
      status = (*clientFn) (socket, swap, message, argp);
      /* Close out the varargs. */
      va_end (argp);
    }
#ifdef LIBDBG
  else
    fprintf (stderr, ": done (no waiting).\n");
#endif

#ifndef NOALARM
  /* Clear the timer. */
  niceClearAlarm (clientAlarm);
#endif
  /* Close the socket. */
  close (socket);
  return (status);
}

/**********************************************************************
 * The next two functions are used to register a message handling
 * function (the first argument) so it can be invoked by niceCheck on
 * receipt of an appropriate message. The function niceHandler is
 * meant to be used by the application to associate a server function
 * with a specific QID.  The other arguments determine properties
 * associated with the protocol we are defining.
 *
 * The first property, SYNCHPTCL, if TRUE, means that the handler is
 * interactive, involving an exchange of messages between server and
 * client, with the server usually expected to go first. If FALSE, it
 * means that the client is not waiting around to hear from the
 * server; it is instead an asynchronous, ``send and forget'' type of
 * interaction.
 *
 * The second property, SAFETY, if TRUE, means that the handler can be
 * invoked safely anywhere within the server side's execution
 * flow. Typically, a server will invoke a handler by explicitly
 * calling niceCheck(), but handlers may also be invoked by other
 * means (e.g., by a timer expiring). If SAFETY is FALSE, the client
 * side will be repeatedly deferred until an explicit call to
 * niceCheck() is made. Important: SYNCHPTCL=FALSE implies that
 * SAFETY=TRUE, since asynchronous protocols cannot be deferred.
 *
 * The third property, ROOTONLY, if TRUE, means that only the root
 * server should handle this request from a nonchild. If FALSE, the
 * request is rejected. This can be useful for input/output
 * processing. Note: barrier nodes should also implement ROOTONLY
 * handlers, but these should only relay information to the "real"
 * root beyond the barrier.
 *
 * The last two properties, EPOCHINC and EPOCHEQP, govern the use of
 * the protocol ``epoch,'' where, for purposes of synchronization,
 * interactions between client and server are divided into ``epochs.''
 * A particular protocol might initiate a new epoch (EPOCHINC is true)
 * by incrementing the appropriate epoch count, and/or might be
 * rejected if the epoch number doesn't match (EPOCHEQP is true)
 * (happens when a message is delayed or received out of sequence).
 * The semantics of epochs are to check first, then increment if
 * necessary.
 **********************************************************************/
void
niceHandler (int qid, serverFn serverFn,
	     int synchronous, int safety, int rootOnly,
	     int epochInc, int epochEqp)
{
  if (qid < 0 ||
      qid >= NICEMAXQID ||
      nicehandler[qid] || (rootOnly && (epochInc || epochEqp)))
    {
      /* Illegal QID, QID already used, or inconsistent properties
       * (e.g., root processor can't keep track of epochs for
       * arbitrary nodes). */
      fprintf (stderr, "Bad QID %d.\n", qid);
      niceQuit (0);
      return;
    }

#ifdef LIBDBG
  nicecalls[qid] = nicewaits[qid] = 0;
#endif
  /* We still register a rootOnly handler for non root nodes, because
   * we can't know whether or not we're a barrier node or not. Recall
   * barrier nodes look like root nodes, but have roots of their own
   * and must therefore forward messages received from within their
   * private network to their external root. 
   *
   * Also, we used to require asynchronous handlers to be safe: now
   * they can be deferred, just like synchronous messages. */
  nicehandler[qid] = serverFn;
  /* Set up protocol property bits. */
  niceprops[qid] =
    (synchronous ? SYNCHPTCL : 0) |
    (safety ? SAFETY : 0) |
    (rootOnly ? ROOTONLY : 0) |
    (epochInc ? EPOCHINC : 0) | 
    (epochEqp ? EPOCHEQP : 0);
#ifdef LIBDBG
  fprintf (stderr, "Handler %d [%d%d%d:%d%d] => %hd\n",
	   qid, (synchronous) ? 1 : 0, (safety) ? 1 : 0, (rootOnly) ? 1 : 0,
	   (epochInc) ? 1 : 0, (epochEqp) ? 1 : 0, niceprops[qid]);
#endif
  return;
}

/**********************************************************************
 * Monitor nicesock, the application request line, and invoke the
 * appropriately registered handler function on receipt of a query.
 * There are two separate versions of niceServiceSocket(), the function
 * that actually does all the real work.
 **********************************************************************/
void
niceCheck ()
{
  SignalCell *cell, *last;

  if (niceParallel ())		/* Ignore unless MASTER or SLAVE */
    {
#ifndef NOALARM
#ifndef NODEFER
      /* OK, flush the defer timer that may be running. */
      niceClearAlarm (deferAlarm);
#endif
#endif
#ifdef LIBDBG
      fprintf (stderr, "niceCheck: servicing pending connections.\n");
#endif
      /* First, you need to service any previously deferred unsafe
       * asynchronous messages. Since we can't expect the clients to
       * resubmit them, they are instead queued internally and
       * processed (in the order they were received) only at "safe"
       * times, like now, indicated by explicit invocation of
       * niceCheck. */
      while ((cell = signals))
	{
	  /* Since signals are pushed onto the signals list, we need
	   * to peel these off from the end of the list. Inefficient,
	   * perhaps, but there should never be very many of these
	   * waiting to go. Start by queueing up the last cell. */
	  last = NULL;
	  while (cell->next)
	    {
	      last = cell;
	      cell = cell->next;
	    }
#ifdef LIBDBG
	  fprintf (stderr, "->Held QID %d from %d", cell->qid, cell->qhid);
#endif
	  /* Check the epochs if they are meaningful. Semantically,
	   * the epochs must be checked when the server is about to be
	   * invoked, and not when the message is received. If the
	   * epoch check is not met, ignore the message. */
	  if (!(niceprops[cell->qid] & EPOCHEQP) || cell->qepoch == niceepoch[cell->qhid])
	    {
	      /* Execute last cell: first, increment the server epoch
	       * number if the QID so specifies. */
	      if (niceprops[cell->qid] & EPOCHINC)
		niceepoch[cell->qhid]++;
	      
	      /* Invoke the handler. Note that we know this is an
	       * asynchronous protocol (or it would never have been
	       * deferred in this fashion) so we can send bogus
	       * socket, endian swap, and message arguments to the
	       * handler (there can never be any reply to the client,
	       * byte swapped or otherwise). Also, we know that any
	       * authentication issues (except epochs, which we just
	       * checked above) have been checked at deferral time, so
	       * they needn't be rechecked here. */
	      if (nicehandler[cell->qid])
		{
#ifdef LIBDBG
		  fprintf (stderr, ": handle.\n");
#endif
		  (*nicehandler[cell->qid]) (cell->qhid, ERROR, ERROR, NULL);
		}
	    }
#ifdef LIBDBG
	  else
	    fprintf (stderr, ": ignore.\n");
#endif
	  
	  /* Flush the (now serviced) signal cell. */
	  if (last)
	    last->next = NULL;
	  else
	    signals = NULL;
	  cell->next = freeSignals;
	  freeSignals = cell;
	}

      /* Service any new incoming messages. */
      niceServiceSocket (FALSE);
#ifndef NOALARM
#ifndef NODEFER
      /* Before you return, reestablish the timer that invokes
       * niceDefer() in case you don't invoke niceCheck() often
       * enough. */
      deferAlarm = niceSetInterrupt (DEFERINT, niceDefer, FALSE);
#endif
#endif
    }
}

void
niceDefer ()
{
  if (niceParallel ())		/* Ignore unless MASTER or SLAVE */
    {
      /* No need to flush the defer timer, since firing it is what got
       * you here in the first place. */
#ifdef LIBDBG
      fprintf (stderr, "niceDefer: deferring pending connections.\n");
#endif
      /* Service any incoming messages, deferring most of them. */
      niceServiceSocket (TRUE);
#ifndef NOALARM
#ifndef NODEFER
      /* Before you return, reestablish the timer that invokes
       * niceDefer() in case you don't invoke niceCheck() often
       * enough. */
      deferAlarm = niceSetInterrupt (DEFERINT, niceDefer, FALSE);
#endif
#endif
    }
}

void
niceServiceSocket (int defer)
{
  int socket, qid, qhid, qpid, qepoch, qendian;
  char message[NICEMAXSTRINGLEN + 1];
  struct in_addr saddr;
#ifndef NOALARM
  int serverAlarm;
#endif
  SignalCell *new;

  /* Here we'll use a zero timeout, since we don't want to block
   * and wait but only accept messages that are ready to go. Note
   * also that we keep servicing messages until all pending
   * messages are handled. This is so, if we should have to exit,
   * we can do so with minimal chance of getting into a deadlock
   * situation. 
   *
   * Note that checkSocket() also returns the source ipaddress as
   * reported by the network layer. */
  while ((socket = checkSocket (nicesock, &saddr, 0)) >= 0)
    {
#ifndef NOALARM
      /* Set a timer in case something blocks during the execution of
       * the protocol. A protocol is declared hung/blocked and
       * attempts to read or write are abandoned when it hangs for
       * HUNGINT seconds, typically a large number. TODO: This should
       * probably NOT be a constant, but rather should be set based on
       * some factor of the observed RTT. */
      if ((serverAlarm = niceSetTimeout (HUNGINT)) == FALSE)
	{
	  /* The other end of the interaction has blocked. You shouldn't
	   * need to clear the timer, which we know has already expired,
	   * since that's how we got here. */
#ifdef LIBDBG
	  fprintf (stderr, "Blocked server [%d] protocol; aborting.\n", qid);
#endif
	  close (socket);
	  continue;
	}
#endif
      /* Obtain the QID, authentication, and synchronization
       * information from the initial message. If the message is
       * garbled, flush it and continue servicing messages on this
       * or other sockets. */
      if (sockscanf (socket, message, NICEQRYSTRDDDDD,
		     &qid, &qhid, &qpid, &qepoch, &qendian) != 5)
	{
#ifdef LIBDBG
	  fprintf (stderr, "niceCheck: ill-formed query [%s].\n", message);
#endif
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (serverAlarm);
#endif
	  close (socket);
	  continue;
	}
#ifdef LIBDBG
      fprintf (stderr, "->Rcvd QID %d from %d", qid, qhid);
      if (niceprops[qid] & (EPOCHEQP | EPOCHINC))
	fprintf (stderr, " [epoch %d =? %d]",
		 qepoch, qhid >= 0 ? niceepoch[qhid] : -1);
#endif
      /* Service request. */
      if (qid == NQIDQUIT && 	/* Quit message */
	  qhid == PARENT &&	/* From my parent */
	  qpid == niceppid)	/* Legit PID */
	{
	  /* My master is telling me to exit. No acknowledgement is
	   * necessary, as NQIDQUIT is treated implicitly as
	   * asynchronous (SYNCHPTCL=false) and always safe
	   * (SAFETY=true). Also, EPOCHEQP is implicitly false, since
	   * a command to exit is respected whenever issued. The
	   * action here is to invoke niceQuit to set the global "die"
	   * flag. When the flag is set, other requests are rejected
	   * before exiting gracefully. */
	  niceQuit (0);
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (serverAlarm);
#endif
#ifdef LIBDBG
	  fprintf (stderr, ": quit.\n");
#endif
	  close (socket);
	  continue;
	}
      else if (die != TRUE && 		/* I'm not about to die */
	       qid == NQIDINIT && 	/* INIT request */
	       qhid < 0 && 		/* Unregistered child */
	       defer == FALSE &&	/* Not currently deferring NQIDINIT */
	       nicehcnt < MAXCHILDREN)	/* I have room */
	{
	  /* This is a child client connecting for the very first
	   * time, and I have room for it. Note that NQIDINIT is
	   * always considered "unsafe," which means it can be
	   * deferred if received at an inconvenient time.
	   *
	   * Give the new child the authentication information that
	   * will be used in future exchanges. First, tell the child
	   * we're listening. */
	  if (sockprintf (socket, message, NICEACKSTRDD, 
			  NICEACCEPT, (qendian!=niceendian)) == ERROR)
	    {
#ifdef LIBDBG
	      fprintf (stderr, ": write error [%s].\n", message);
#endif
#ifndef NOALARM
	      /* Clear the timer. */
	      niceClearAlarm (serverAlarm);
#endif
	      close (socket);
	      continue;
	    }

	  /* Next, get the child's port information. We already have
	   * the child's source IP address information from the
	   * network layer, as returned to us by checkSocket(). But we
	   * don't yet know what port the child is listening to.  We
	   * anticipate the increment of children on success in order
	   * to read the port number directly into the appropriate
	   * cell of niceport[]. */
	  if (sockscanf (socket, message, NICEINFSTRHU, 
			 &niceport[nicehcnt + 1]) != 1)
	    {
#ifdef LIBDBG
	      fprintf (stderr, ": bad slave.\n");
#endif
#ifndef NOALARM
	      /* Clear the timer. */
	      niceClearAlarm (serverAlarm);
#endif
	      close (socket);
	      continue;
	    }

	  /* Accept the new child: store the network-reported IP
	   * address into the niceaddr[] array (you've already cached
	   * the port address in the niceport array). */
	  nicehcnt++;
	  memcpy (&niceaddr[nicehcnt], &saddr, sizeof (saddr));
#ifdef LIBDBG
	  fprintf (stderr, " [new child #%d %s:%hu]",
		   nicehcnt, inet_ntoa (saddr), niceport[nicehcnt]);
#endif
	  /* Send the root application node's particulars. These are
	   * used by your descendent in order to contact the root
	   * directly with a ROOTONLY protocol, as well as for
	   * authentication (especially when in direct communication
	   * with root node).  Also, return index into hosts table (to
	   * be used by the child to identify itself in the future)
	   * and your pid (which will be used for authentication).
	   *
	   * One complicating factor: if you are a barrier node, your
	   * descendent should view you as their root, rather than the
	   * global root. That's because your descendants may not in
	   * fact have access to the "outside world" where the root
	   * resides. Instead, it's safer for them to send information
	   * destined to the root directly to you, the barrier node,
	   * for relaying. 
	   *
	   * Another complicating factor. If you're the root (or a
	   * barrier node), you can't be assured of knowing your own
	   * IP address in niceaddr[SELF]. That's because some
	   * machine's network configuration may return, for example,
	   * localhost to gethostbyname(), which would result in
	   * placing 127.0.0.1 in niceaddr[SELF]. So the convention
	   * is, if you're the root or a barrier, send LOCALHOST, and
	   * have the child get the IP address from the network layer.
	   *
	   * TODO: We should really notice if I perform IP
	   * masquerading. If I don't, my descendents should see me as
	   * their root. If I do, then my descendents can access the
	   * real root directly. */
	  if (sockprintf (socket, message, NICEINFSTRSHUDDD,
			  ((niceRoot() || (nicewall==TRUE))?LOCALHOST:inet_ntoa(niceaddr[ROOT])),
			  ((nicewall==TRUE)?niceport[SELF]:niceport[ROOT]),
			  nicehcnt, nicepid, 
			  ((nicewall==TRUE)?nicepid:nicerpid)) == ERROR)
	    {
#ifdef LIBDBG
	      fprintf (stderr, ": write error [%s].\n", message);
#endif
#ifndef NOALARM
	      /* Clear the timer. */
	      niceClearAlarm (serverAlarm);
#endif
	      close (socket);
	      continue;
	    }

	  /* Set the qhid variable so that the remainder of this
	   * handler knows the new qhid. */
	  qhid = nicehcnt;

	  /* Initialize the epoch counter for this child. */
	  niceepoch[qhid] = 0;

	  /* Now fall through and handle the rest of this
	   * protocol. Note that you haven't yet cancelled the
	   * alarm. */
	}
      else if (die != TRUE && 		/* I'm not about to die */
	       defer == TRUE && 	/* I'm deferring messages. */
	       !(niceprops[qid] & SAFETY)) /* Message is unsafe. */
	{
	  /* I'm deferring all of my incoming non-safe messages
	   * (including NQIDINIT, explicitly marked "unsafe" in
	   * niceInit). Send a rejection notice and return to
	   * processing other incoming messages. The receiving client
	   * should eventually try again. 
	   * 
	   * Note that unsafe asynchronous protocols are handled
	   * differently, internally, than synchronous protocols,
	   * because the client can't be told to try again later --
	   * its not listening to you after sending the message. Such
	   * signals are instead queued locally and then processed at
	   * a later "safe" time. 
	   *
	   * Message deferral occurs before we check for credential
	   * rejection or a deferred NQIDINIT message would be
	   * rejected based on its qhid. */
	  if (niceprops[qid] & SYNCHPTCL)
	    {
	      /* Unsafe synchronous signal. Tell client to try
	       * later; don't worry about endianess. */
	      sockprintf (socket, message, NICEACKSTRDD, NICEDEFER, ERROR);
	    }
	  else
	    {
	      /* Unsafe asynchronous signal. Get an empty cell to
	       * record the appropriate info. Note that since the
	       * signal is asynchronous, we needn't care about
	       * endianess. */
	      if (!freeSignals)
		new = (SignalCell *) malloc (sizeof(SignalCell));
	      else
		{
		  new = freeSignals;
		  freeSignals = new->next;
		}
	      /* Set up appropriate contents in the new cell. */
	      new->qid = qid;
	      new->qhid = qhid;
	      new->qpid = qpid;
	      new->qepoch = qepoch;
	      /* Push new cell onto signals list. */
	      new->next = signals;
	      signals = new;
	    }
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (serverAlarm);
#endif
#ifdef LIBDBG
	  fprintf (stderr, ": defer.\n");
#endif
	  close (socket);
	  continue;
	}
      else if ((die == TRUE) ||		/* I'm about to die */
	       (qid == NQIDINIT && 	/* I've no room for child */
		qhid < 0 &&
		nicehcnt == MAXCHILDREN) ||
	       (qid < 0) ||		/* Illegal QID */
	       (qid >= NICEMAXQID) ||
	       (!nicehandler[qid]) ||	/* No registered handler */
	       (qhid > 0 && qpid != nicepid) || 	/* Child calling. */
	       (qhid == 0 && qpid != niceppid) ||	/* Parent calling. */
	       (qhid < 0 && 		/* Root/barrier only. */
		!((niceRoot() || nicewall) && 
		  qpid == nicepid &&
		  (niceprops[qid] & ROOTONLY))))
	{
	  /* I'm either in the process of quitting, I have no room for
	   * the new child, or the incoming message is somehow
	   * unacceptable (an unregistered handler, out-of-range QID,
	   * or incorrect authentication). Send a credential rejection
	   * notice (only if this is a synchronous protocol; just
	   * ignore the message if it is asynchronous because the
	   * client isn't waiting around to hear from you anyway) and
	   * return to processing other incoming messages. The
	   * receiving client will exit if it is a synchronous
	   * protocol, and will continue in ignorance otherwise. Don't
	   * worry about endianess here. */
	  if (niceprops[qid] & SYNCHPTCL)
	    sockprintf (socket, message, NICEACKSTRDD, NICEREFUSE, ERROR);
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (serverAlarm);
#endif
#ifdef LIBDBG
	  fprintf (stderr, ": reject.\n");
#endif
	  close (socket);
	  continue;
	}
      else if ((niceprops[qid] & EPOCHEQP) && (qepoch != niceepoch[qhid]))
	{
	  /* Message received out of sequence; epoch number
	   * mismatch. Just ignore the message. There's no need to
	   * check sockprintf() status, as you're just going to close
	   * the socket and go on to the next message anyway. Don't
	   * worry about endianess here. */
	  sockprintf (socket, message, NICEACKSTRDD, NICEIGNORE, ERROR);
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (serverAlarm);
#endif
#ifdef LIBDBG
	  fprintf (stderr, ": ignore.\n");
#endif
	  close (socket);
	  continue;
	}
      else if (niceprops[qid] & SYNCHPTCL)
	{
	  /* A properly registered synchronous protocol. Tell your
	   * correspondent to go ahead. */
	  if (sockprintf (socket, message, NICEACKSTRDD, NICEACCEPT, (qendian!=niceendian)) == ERROR)
	    {
#ifdef LIBDBG
	      fprintf (stderr, ": write failure [%s].\n", message);
#endif
#ifndef NOALARM
	      /* Clear the timer. */
	      niceClearAlarm (serverAlarm);
#endif
	      close (socket);
	      continue;
	    }
	}

#ifdef LIBDBG
      fprintf (stderr, ": handle.\n");
#endif
      /* OK, ready to go on processing the query. First, increment
       * the server epoch number if the QID so specifies. */
      if (niceprops[qid] & EPOCHINC)
	niceepoch[qhid]++;

      /* Finally, go ahead and do whatever else it is the protocol
       * should do. Note that the handlers for synchronous protocols
       * require the socket and message buffer, while the asynchronous
       * handlers could as easily get bogus values for these
       * arguments. */
      if (nicehandler[qid])
	(*nicehandler[qid]) (qhid, socket, (qendian!=niceendian), message);

#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (serverAlarm);
#endif
      /* Close the socket. */
      close (socket);
    }
  /* Check the die status variable and exit gracefully if it is
   * set. */
  if (die == TRUE)
    {
      niceExit ();
      exit (EXIT_SUCCESS); /* Only executed by master. */
    }
}

/**********************************************************************
 * NICE daemon requests. The following functions are used to send the
 * appropriate NICE query to a daemon. These functions differ from
 * some of the other daemon request functions in that they are only
 * used by NICE applications and not by either control programs (e.g.,
 * niceq) or the NICE daemon itself. Note that interactions with the
 * NICE daemon are not as easy to code as the interactions between
 * applications; for example, we must manually manage the establishing
 * the socket and closing the socket. But these functions are not
 * intended to be modified by applications programmers, so they
 * needn't be so refined.
 **********************************************************************/

/**********************************************************************
 * Ping a NICE daemon (not just the host). Returns TRUE if the NICE
 * daemon is available, FALSE otherwise. Basically, try to open a
 * socket connection.
 **********************************************************************/
int
nicePing (int hid)
{
  int socket, status;
  char message[NICEMAXSTRINGLEN + 1];
#ifndef NOALARM
  int clientAlarm;
#endif

  /* No daemon running on hostname; generate an error and exit. Note
   * that we are not being terribly patient here, since hailSocket
   * could return FALSE (indicating a non fatal error) as easily as
   * ERROR. */
  if (hailSocket (hid, NICEPORT, &socket) != TRUE)
    return (FALSE);

#ifndef NOALARM
  /* Set a timer in case something blocks during the execution of the
   * protocol. A protocol is declared hung/blocked and attempts to
   * read or write are abandoned when it hangs for HUNGINT seconds,
   * typically a large number. TODO: This should probably NOT be a
   * constant, but rather should be set based on some factor of the
   * observed RTT. */
  if ((clientAlarm = niceSetTimeout (HUNGINT)) == FALSE)
    {
      /* The other end of the interaction has blocked. You shouldn't
       * need to clear the timer, which we know has already expired,
       * since that's how we got here. */
#ifdef LIBDBG
      fprintf (stderr, "Blocked client [%d] protocol; aborting.\n", NQIDPING);
#endif
      close (socket);
      return (ERROR);
    }
#endif

  /* Send a ping query. */
  if (sockprintf (socket, message, NICEQRYSTRD, NQIDPING) == ERROR)
    {
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      close (socket);
      return (FALSE);
    }

  if (sockscanf (socket, message, NICEACKSTRD, &status) != 1)
    status = FALSE;

  /* Close the socket and return status. */
#ifndef NOALARM
  /* Clear the timer. */
  niceClearAlarm (clientAlarm);
#endif
  close (socket);
  return (status);
}

/**********************************************************************
 * Ask your local NICE daemon to send some volunteers. Return a value
 * indicating the upper bound on the number of naggers you can expect,
 * or ERROR.  Note that even if it appears you can expect naggers, it
 * doesn't mean you will necessarily get any naggers, since the
 * executable you want to run may not be available on the child
 * machines.
 **********************************************************************/
int
niceSolicit (char *executable)
{
  int socket, count;
  char message[NICEMAXSTRINGLEN + 1];
#ifndef NOALARM
  int clientAlarm;
#endif

  /* Connect to NICE daemon socket on your own machine. Note that we
   * are not being terribly patient here, since hailSocket could
   * return FALSE (indicating a non fatal error) as easily as ERROR,
   * but since we are connecting to ourselves, network problems or a
   * connect timeout seem unlikely. */
  if (hailSocket (SELF, NICEPORT, &socket) != TRUE)
    return (ERROR);

#ifndef NOALARM
  /* Set a timer in case something blocks during the execution of the
   * protocol. A protocol is declared hung/blocked and attempts to
   * read or write are abandoned when it hangs for HUNGINT seconds,
   * typically a large number. TODO: This should probably NOT be a
   * constant, but rather should be set based on some factor of the
   * observed RTT. */
  if ((clientAlarm = niceSetTimeout (HUNGINT)) == FALSE)
    {
      /* The other end of the interaction has blocked. You shouldn't
       * need to clear the timer, which we know has already expired,
       * since that's how we got here. */
#ifdef LIBDBG
      fprintf (stderr, "Blocked client [%d] protocol; aborting.\n", NQIDSOLICIT);
#endif
      close (socket);
      return (ERROR);
    }
#endif

  /* Send a request for volunteers. */
  if (sockprintf (socket, message, NICEQRYSTRD, NQIDSOLICIT) == ERROR)
    {
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      close (socket);
      return (ERROR);
    }

  /* Wait for acknowledgement. */
  if (sockscanf (socket, message, NICEINFSTRDD, &count, &nicewall) != 2)
    {
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      close (socket);
      return (ERROR);
    }

  /* The acknowledgement gives you an upper bound on the number of
   * naggers you can expect. If it's zero, then you can just give up
   * now. */
  if (count <= 0)
    {
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      close (socket);
      return (count);
    }

  /* Send command line of job to run along with the port number you
   * expect to be contacted on. You needn't supply anything else (like
   * hname[SELF] or haddr[SELF]) since these will also be available on
   * the local niced that answers your solicitation. */
  if (sockprintf (socket, message, NICEINFSTRSHU,
		  executable, niceport[SELF]) == ERROR)
    {
#ifndef NOALARM
      /* Clear the timer. */
      niceClearAlarm (clientAlarm);
#endif
      close (socket);
      return (ERROR);
    }

#ifndef NOALARM
  /* Clear the timer. */
  niceClearAlarm (clientAlarm);
#endif
  close (socket);
  return (count);
}

/**********************************************************************
 * sockread() and sockwrite() provide streaming primitives not found
 * in nicecom. They're not in nicecom because they require accounting
 * for endianess, which is a function of the API, and not in niced at
 * all.
 **********************************************************************/
/**********************************************************************
 * sockwrite() returns the number of objects of specified size written
 * to socket or ERROR if, for example, the other end of the socket has
 * disappeared. The message buffer must have been declared to be
 * NICEMAXSTRINGLEN+1 bytes long to accomodate message terminator.
 *
 * Remember: Linux/i386 is little-endian but, e.g., HPUX/PA-RISC is
 * big-endian. The swap argument tells you whether you need to swap
 * bytes on this socket or not (byte swapping is always performed on
 * the write end).
 **********************************************************************/
int
sockwrite (void *data, size_t size, size_t count, int socket, int swap, char *message)
{
  int value, length;
  int i, j, k;
  unsigned char* mptr;
  unsigned char tmp;
#if FALSE
  /* Zero out the message string. */
  memset (message, (char) NULL, NICEMAXSTRINGLEN + 1);
#endif
  /* Construct the appropriate message First, figure out how many
   * objects you can pack in message buffer, then pack them in. */
  count = MIN (count, (NICEMAXSTRINGLEN / size));
  length = (size * count);
  memcpy (message, (char *) data, length);
  message[length] = NICEEOR;

  /* Swap bytes in place as determined by endianess and specified by
   * swap argument. Should work for whatever size we're dealing with.
   * We could have instead opted to use, e.g., ntohl(), but since
   * network byte order is big-endian and most Linux machines (those
   * that are i386 or derivative) are little-endian, this would entail
   * a lot of unnecessary swapping, resulting in a performance hit. So
   * instead we byte swap only when endian differences are noted. */
  if (swap)
    {
      mptr = (unsigned char *) message;
      for (i = 0; i < count; i++, mptr = mptr+size)
	{
	  j = 0;
	  k = size - 1;
	  while (j < k)
	    {
	      tmp = mptr[j];
	      mptr[j] = mptr[k];
	      mptr[k] = tmp;
	      j++;
	      k--;
	    }
	}
    }

  /* Write message line to socket. Keep trying unless you get an
   * unrecoverable error (interrupts are the only recoverable errors
   * considered here) or no bytes were actually written. Make sure you
   * explicitly send the null message terminator byte (hence 1 more
   * byte than length). */
  while ((value = send (socket, message, (length+1), 0)) < 0)
      {
#ifdef LIBDBG
	/* An error here usually means the other end of the socket has
	 * evaporated. */
	fprintf (stderr, "Socket %d (%dx%d): write error (%d) %s\n",
		 socket, count, size, errno, strerror (errno));
#endif
	return (ERROR);
      }

  /* Return value is number of objects of size size actually written,
   * not the number of bytes. */
  return (count);
}

/**********************************************************************
 * sockread() returns the number of objects of specified size bound
 * from the socket or ERROR if, for example, the other end of the
 * socket has disappeared.  The message buffer must have been declared
 * to be NICEMAXSTRINGLEN+1 bytes long.
 *
 * Remember: Linux/i386 is little-endian but, e.g., HPUX/PA-RISC is
 * big-endian. Since byte swapping is always performed on the write
 * end, swap is included here just for symmetry.
 **********************************************************************/
int
sockread (void *data, size_t size, size_t count, int socket, int swap, char *message)
{
  int value = 0;
  char *msgptr = &message[0];
#if FALSE
  /* Zero out the message buffer. */
  memset (message, (char) NULL, NICEMAXSTRINGLEN + 1);
#endif
  /* First, figure out how many objects you can pack in a message
   * buffer. */
  count = MIN (count, (NICEMAXSTRINGLEN / size));

  /* Next, read up to count elements from the socket. Note that you
   * may need to make multiple calls to recv in order to get the
   * entire message. You'll never get more than the message from the
   * stream, however, because we always alternate message and
   * response, so your correspondent will never stream more than one
   * message to you.
   * 
   * So: keep going until you get a string terminating byte, unless
   * you get an unrecoverable error (the only recoverable error here
   * is an interrupt). You might get an EOF (recv returns 0), which is
   * OK if the socket closes after the message, but not OK if the
   * message is not properly terminated. */
  while ((value = recv (socket, msgptr, (NICEMAXSTRINGLEN + 1 - (msgptr - message)), 0)) != 0)
    {
      if (value > 0)
	{
	  /* Advance msgptr to end of current message chunk. */
	  msgptr = msgptr + value;

	  /* Quit recv'ing if you got the message terminator in the
	   * appropriate position. */
	  if ((*(msgptr - 1) == NICEEOR) &&
	      ((msgptr - 1 - message) == count * size))
	    break;
	}
      else if (errno != EINTR)
	{
#ifdef LIBDBG
	  /* An error here usually means the other end of the socket has
	   * evaporated. */
	  fprintf (stderr, "Socket %d (%dx%d): read error (%d) %s.\n",
		   socket, count, size, errno, strerror (errno));
#endif
	  return (ERROR);
	}
    }

  /* OK, you've fallen out of the loop. Make sure that you've got a
   * non-zero-length message of proper length (modulo size) with a
   * terminating byte at the end; if not, generate an error. Careful:
   * you initialized message to all NULLs, so you need to be sure that
   * the transmitted terminator was indeed received in the byte before
   * msgptr. */
  if ((value = (msgptr - 1 - message)) == 0 || 
      *(msgptr - 1) != NICEEOR ||
      (value != count * size))
    {
#ifdef LIBDBG
      /* An error here might have different causes. */
      fprintf (stderr, "Socket %d (%dx%d = %d bytes) read error: %s.\n", 
	       socket, count, size, value,
	       (value==0)?"no message":((value != count*size)?"wrong size":"no EOR"));
#endif
      return (ERROR);
    }

  /* Store the message in data. */
  memcpy ((char*) data, message, value);

  /* Return number of objects of size size read, not the number of
   * bytes. */
  return (value / size);
}
