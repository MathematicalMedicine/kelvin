/**********************************************************************
 * NICE Test Program
 * Alberto Maria Segre
 *
 * Copyright (c) 1997-2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software 
 * for non-profit educational purposes only.
 **********************************************************************/
#include "niceaux.h"
#include "nicecom.h"
#include "niceapi.h"
#include <stdio.h>		/* Needed for printf */
#include <stdlib.h>		/* EXIT_SUCCESS, EXIT_ERROR */
#include <string.h>		/* Needed for strchr, strcmp */
#include <unistd.h>		/* Needed for getpid */

/* Number of iterations for this silly test program. 20 is good. */
#define ROUNDS 20

/* These functions define the protocol to establish the problem
 * specification on a new nagger. We'll use the INFOSTRING constant as
 * the string for socket transmission. */
#define INFOSTRING "%%TEST %d %s"
int initClient (int socket, int swap, char message[], va_list argp);
void initServer (int nid, int socket, int swap, char message[]);

/* This ID and functions define the protocol for "prodding" the
 * root. */
#define NQIDPROD 1
void prodded (int nid, int socket, int swap, char message[]);
int count = 1;		/* Total number of nodes. */

/* This ID and functions define the protocol for "tweaking" your
 * master, possibly recursively. */
#define NQIDTWEAK 2
int tweak (int socket, int swap, char message[], va_list argp);
void tweaked (int nid, int socket, int swap, char message[]);

/* This ID and function defines the protocol for "spanking" your child
 * on occasion when it is impertinent enough to tweak you. */
#define NQIDSPANK 3
void spanked (int nid, int socket, int swap, char message[]);

/* This ID and functions define the protocol for "pinching" the
 * root directly. */
#define NQIDPINCH 4
int pinch (int socket, int swap, char message[], va_list argp);
void pinched (int nid, int socket, int swap, char message[]);

/**********************************************************************
 * Silly sample NICE application. Spawns naggers on other
 * systems. Initially, naggers and masters exchange PIDs (as if this
 * were the problem specification). In addition, a nagger occasionally
 * chooses to "tweak" its master, who yelps on stdout and recursively
 * tweaks its own master, or "pinches" the root directly.
 *
 * Even if fairly silly, this application illustrates the steps
 * involved in parallelizing an application with NICE.
 *
 * First, you need to invoke the niceInit function when you start,
 * being sure to pass two basic "protocol defining" functions (here
 * initClient and initServer) which will be used to propagate the
 * initial problem description to a new nagger.
 *
 * Second, use niceTalk to initiate protocols with your master (0) or
 * any of your children (1...), passing an appropriate "protocol
 * defining" client function. Note that the appropriate handler
 * (server side of the "protocol defining" function) must have been
 * previously "registered" with the application interface, which will
 * automatically assign and return an appropriate message
 * identifier. For example, see the "tweaking" protocol given below.
 * 
 * Third, you must be sure to invoke niceCheck on occasion to handle
 * incoming messages.
 *
 * Finally, make sure to invoke niceExit when you exit.
 **********************************************************************/
int
main (int argc, char *argv[])
{
  int i, n, status;

  /* Connect to the NICE environment. */
  if ((n = niceInit (argc, argv, initClient, initServer)) == ERROR)
    {
      printf ("Initialization error.\n");
      return (ERROR);
    }

  /* Check status of process using niceSerial, niceRoot, etc. */
  if (niceSerial ())
    {
      printf ("Running serially.\n");
      return (FALSE);
    }
  else if (niceRoot ())
    printf ("Nagging (root node; %d slaves).\n", n);
  else if (niceMaster ())
    printf ("Nagging (internal node; %d slaves).\n", n);
  else if (niceSlave ())
    printf ("Nagging (leaf node; %d slaves).\n", n);
  else
    {
      printf ("No slaves available.\n");
      return (FALSE);
    }

  /* Register any additional protocol handlers and their appropriate
   * message identifier using niceHandler. The properties refer to the
   * protocols associated with, e.g., "tweaking" or "pinching". The
   * first implies this is a synchronous protocol (meaning client and
   * server will exchange information). The second means the protocol
   * can be safely executed anywhere in the server's thread of
   * execution. The third means protocol reports direct to root. For
   * the last two, we needn't worry about the epoch stuff for this
   * simple example. */
  niceHandler (NQIDPROD, prodded, FALSE, TRUE, TRUE, FALSE, FALSE);
  niceHandler (NQIDTWEAK, tweaked, TRUE, TRUE, FALSE, TRUE, FALSE);
  niceHandler (NQIDSPANK, spanked, FALSE, FALSE, FALSE, FALSE, TRUE);
  niceHandler (NQIDPINCH, pinched, TRUE, TRUE, TRUE, FALSE, FALSE);

  /* Seed the random number generator with your PID. */
  SEED (getpid ());

  /* Get the root process to count you. */
  if (!niceRoot ())
    {
      while ((status = niceTalk (ROOT, NQIDPROD, NULL)) == FALSE)
	fprintf (stderr, "Prodding deferred...\n");
      if (status == ERROR)
	{
	  fprintf (stderr, "Failed to prod root.\n");
	  niceExit ();
	}
    }

  /* Now go about your business, remembering to call niceCheck on 
   * occasion. Initiate protocols using niceTalk. */
  for (i = 0; i < ROUNDS; i++)
    {
      sleep (1);
      printf ("%d\n", i);
      /* Check and handle incoming messages. */
      niceCheck ();

      /* On occasion, tweak your master. */
      if (niceSlave () && RANDINT (10) == 0)
	{
	  while ((status = niceTalk (PARENT, NQIDTWEAK, tweak)) == FALSE)
	    fprintf (stderr, "Tweaking deferred...\n");
	  if (status == ERROR)
	    {
	      fprintf (stderr, "Tried and failed to tweak.\n");
	      niceExit ();
	    }
	}
      else if (niceSlave () && RANDINT (10) == 0)
	{
	  /* Note explicit use of args, since pinch is called recursively! */
	  while ((status = niceTalk (ROOT, NQIDPINCH, pinch, getpid(), nicehost)) == FALSE)
	    fprintf (stderr, "Pinching deferred...\n");
	  if (status == ERROR)
	    {
	      fprintf (stderr, "Tried and failed to pinch.\n");
	      niceExit ();
	    }
	}
    }

  /* Print out total number of nodes used. */
  if (niceRoot ())
    fprintf (stderr, "%d nodes total.\n", count);

  /* Terminate connection with the NICE environment; the exit is never
   * actually executed, but keeps the compiler happy. */
  niceExit ();
  exit (EXIT_SUCCESS);
}

/**********************************************************************
 * NICE Protocol defining functions.
 **********************************************************************/

/**********************************************************************
 * initClient defines the protocol for obtaining the problem
 * description from your master.
 *
 * The semantics of the return value are strictly regulated: an ERROR
 * will cause the client to exit, while a FALSE will cause the client
 * to repeat its attempt to obtain initialization information. TRUE
 * (or any other value) indicates successful completion of the
 * initialization protocol.
 **********************************************************************/
#define INITARRAYLEN 1024
#define INITARRAYTYPE double
int
initClient (int socket, int swap, char message[], va_list argp)
{
  int pid = getpid (), mpid;
  char mname[NICEMAXSTRINGLEN+1];
#ifndef NOSTREAM
  int i, j;
  INITARRAYTYPE initarray[INITARRAYLEN];
#endif
  
  /* Send your ID. */
  if (sockprintf (socket, message, INFOSTRING, pid, nicehost) == ERROR)
    {
      fprintf (stderr, "Bad message sent in initClient.\n");
      return (ERROR);
    }

  /* Obtain master's ID. */
  if (sockscanf (socket, message, INFOSTRING, &mpid, &mname) != 2)
    {
      fprintf (stderr, "Bad message rcvd in initClient:\n  =>%s.\n", message);
      return (ERROR);
    }

#ifndef NOSTREAM
  /* Initialize initarray. */
  for (i = 0; i < INITARRAYLEN; i++)
    initarray[i] = (INITARRAYTYPE) i;

  /* Send it. */
  i = 0;
  while (i < INITARRAYLEN)
    {
      if ((j = sockwrite (&initarray[i], sizeof (INITARRAYTYPE), INITARRAYLEN - i, 
			  socket, swap, message)) < 0)
	{
	  fprintf (stderr, "Error while sending initarray.\n");
	  return (ERROR);
	}
      i += j;

      /* Probably should ACK the message here... */
    }
#endif

  fprintf (stderr, "My master is %s [pid = %d]\n", mname, mpid);
  return (TRUE);
}

/**********************************************************************
 * initServer is the server side (master) of the initClient
 * protocol. It returns the problem description to the requesting
 * child process.
 **********************************************************************/
void
initServer (int nid, int socket, int swap, char message[])
{
  int pid = getpid (), cpid;
  char cname[NICEMAXSTRINGLEN+1];
#ifndef NOSTREAM
  int i, j;
  INITARRAYTYPE initarray[INITARRAYLEN];
#endif

  /* Obtain child ID. */
  if (sockscanf (socket, message, INFOSTRING, &cpid, &cname) != 2)
    {
      fprintf (stderr, "Bad message rcvd in initServer.\n");
      return;
    }

  /* Send your ID. */
  if (sockprintf (socket, message, INFOSTRING, pid, nicehost) == ERROR)
    {
      fprintf (stderr, "Bad message sent in initServer.\n");
      return;
    }

#ifndef NOSTREAM
  /* Obtain the test array. */
  i = 0;
  while (i < INITARRAYLEN)
    {
      if ((j = sockread (&initarray[i], sizeof(INITARRAYTYPE), INITARRAYLEN - i, 
			 socket, swap, message)) < 0)
	{
	  fprintf (stderr, "Error while receiving initarray.\n");
	  return;
	}
      i += j;

      /* Probably should ACK the message here... */
    }

  /* Check it. */
  for (i = 0; i < INITARRAYLEN; i++)
    if (initarray[i] != i)
      {
	fprintf (stderr, "Bad value received in initarray (%d).\n", i);
	return;
      }
  fprintf (stderr, "Initarray passed.\n");
#endif

  fprintf (stderr, "My #%d child is %s [pid = %d]\n", nid, cname, cpid);
  return;
}

/**********************************************************************
 * prodded is the asynchronous server-side function for the prod
 * protocol automatically invoked when a NQIDPROD message identifier
 * is received. It counts the descendents.
 **********************************************************************/
void
prodded (int nid, int socket, int swap, char message[])
{
  int status;

  if (niceRoot ())
    count++;
  else
    {
      while ((status = niceTalk (ROOT, NQIDPROD, NULL)) == FALSE)
	fprintf (stderr, "Prodding deferred...\n");
      if (status == ERROR)
	{
	  fprintf (stderr, "Tried and failed to relay prod.\n");
	  niceExit ();
	}
      else
	fprintf (stderr, "Relayed prodding to root.\n");
    }
  return;
}

/**********************************************************************
 * tweak defines the protocol for occasionally prodding your master.
 * The message identifier is NQIDTWEAK.
 **********************************************************************/
int
tweak (int socket, int swap, char message[], va_list argp)
{
  int pid = getpid ();

  fprintf (stderr, "Tweaking my parent...");

  /* Send your ID. */
  if (sockprintf (socket, message, INFOSTRING, pid, nicehost) == ERROR)
    {
      fprintf (stderr, "failure.\n");
      return (ERROR);
    }

  /* Success. */
  fprintf (stderr, "success!\n");
  return (TRUE);
}

/**********************************************************************
 * tweaked is the server-side function for the tweak protocol,
 * automatically invoked when a NQIDTWEAK message identifier is
 * received. It complains to stdout, possibly reprimands the tweaking
 * child, and then goes ahead and tweaks its own master.
 **********************************************************************/
void
tweaked (int nid, int socket, int swap, char message[])
{
  int cpid, status;
  char cname[NICEMAXSTRINGLEN+1];

  /* Obtain child ID. */
  if (sockscanf (socket, message, INFOSTRING, &cpid, &cname) != 2)
    {
      fprintf (stderr, "Bad message rcvd in tweaked.\n");
      return;
    }

  /* Success. */
  fprintf (stderr, "Getting tweaked by %s [pid = %d]\n", cname, cpid);

  /* Now with small probability spank the child that tweaked you. */
  if (RANDINT(10) <= 5)
    {
      fprintf (stderr, "Reprimanding child.\n");
      niceTalk (nid, NQIDSPANK, NULL);
    }
  /* Recursively tweak your own parent. */
 if (!niceRoot ())
    {
      while ((status = niceTalk (PARENT, NQIDTWEAK, tweak)) == FALSE)
	fprintf (stderr, "Tweaking deferred...\n");
      if (status == ERROR)
	{
	  fprintf (stderr, "Tried and failed to relay tweak.\n");
	  niceExit ();
	}
    }
}

/**********************************************************************
 * spanked is the server-side function for the spanking protocol,
 * which is the asynchronous message sent by a parent when it is tired
 * of being tweaked.
 **********************************************************************/
void
spanked (int nid, int socket, int swap, char message[])
{
  fprintf (stderr, "Ouch! Got spanked by my parent!\n");
}

/**********************************************************************
 * pinch defines the protocol for occasionally prodding the root.  The
 * message identifier is NQIDPINCH: note that we pass the pid and
 * hostname explicitly, unlike when tweaking.
 **********************************************************************/
int
pinch (int socket, int swap, char message[], va_list argp)
{
  int cpid;
  char *cname;
  fprintf (stderr, "Pinching the root...");

  /* Set up the arguments. */
  cpid = va_arg (argp, int);
  cname = va_arg (argp, char*); 

  /* Send your ID. */
  if (sockprintf (socket, message, INFOSTRING, cpid, cname) == ERROR)
    {
      fprintf (stderr, "failure.\n");
      return (ERROR);
    }

  /* Success. */
  fprintf (stderr, "success!\n");
  return (TRUE);
}

/**********************************************************************
 * pinched is the server-side function for the pinch protocol,
 * automatically invoked when a NQIDPINCH message identifier is
 * received. It complains to stdout.
 **********************************************************************/
void
pinched (int nid, int socket, int swap, char message[])
{
  int cpid, status;
  char cname[NICEMAXSTRINGLEN+1];

  /* Obtain child ID. */
  if (sockscanf (socket, message, INFOSTRING, &cpid, &cname) != 2)
    {
      fprintf (stderr, "Bad message rcvd in pinched.\n");
      return;
    }

  /* Now pinch your own root, if you are not the "real" root. This can
   * happen when you are, e.g., a barrier node. */
  if (!niceRoot ())
    {
      while ((status = niceTalk (ROOT, NQIDPINCH, pinch, cpid, cname)) == FALSE)
	fprintf (stderr, "Pinching deferred...\n");
      if (status == ERROR)
	{
	  fprintf (stderr, "Tried and failed to relay pinch.\n");
	  niceExit ();
	}
      else
	fprintf (stderr, "Relayed pinching from %s to root.\n", cname);
    }
  else
    {
      /* Success. */
      fprintf (stderr, "Getting pinched by %s [pid = %d]\n", cname, cpid);
    }
}
