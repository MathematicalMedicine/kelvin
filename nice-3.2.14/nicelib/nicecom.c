/**********************************************************************
 * NICE Communications Library
 * Alberto Maria Segre
 *
 * Copyright (c) 1997-2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software
 * for non-profit educational purposes only.
 **********************************************************************/
#include <stdio.h>		/* Needed for fprintf */
#include <stdarg.h>		/* Needed for va_list etc. */
#include <fcntl.h>		/* Needed for open */
#include <unistd.h>		/* Needed for unbuffered write */
#include <string.h>		/* Needed for strncpy, strcmp */
#include <netdb.h>		/* Needed for gethostbyname */
#include <limits.h>		/* USHRT_MAX */
#include <netinet/in.h>		/* Needed for in_addr */
#include <errno.h>		/* ETIMEDOUT, ENETUNREACH */
#ifdef _HPUX_SOURCE
#include <arpa/inet.h>		/* htons for HPUX native CC */
#endif
#include "niceaux.h"
#include "nicecom.h"

/**********************************************************************
 * Some global variables that must be accessible to a broad range of
 * functions.
 **********************************************************************/
char nicehost[NICEMAXSTRINGLEN + 1];	/* My host name. */
int nicesock;				/* My incoming messages. */
int nicehcnt = 0;			/* Number of hosts (< MAXHOSTS). */
struct in_addr *niceaddr;	/* Parent, child host addresses. */
struct in_addr nulladdr = {0};		/* Null address (for comparisons). */
struct in_addr localaddr = {16777343};	/* Localhost address (for comparisons). */

/**********************************************************************
 * Takes the name of a machine and ensures you have the correct full
 * hostname and not some shortened or alternate version. While you're
 * at it, cache the host's network address in the host table arrays at
 * hid.
 *
 * Note that the results are unpredictable when asking about your own
 * hostname (or address), since some machines have arbitrary
 * locally-assigned names that do not correspond to a network address
 * in any way. Furthermore, some machines may return, say, 127.0.0.1,
 * which would not be a valid address from anywhere else on the net.
 *
 * Returns TRUE if successful, ERROR if host not found or name can't
 * be resolved.
 **********************************************************************/
int
qualifyHostname (char *hostname, int hid)
{
  struct hostent *hp;
  int i;

  /* Do a lookup to get the hostentry from DNS. */
  if ((hp = gethostbyname (hostname)) == NULL)
    {
      /* Failure. */
      return (ERROR);
    }
  else if (hp->h_length != ADDRBYTES)
    {
      /* hp->h_length will tell you how many bytes are in the
       * address. This had better be consistent with the size of
       * in_addr structure, which is 32 bits for IPv4. */
#ifdef LIBDBG
      fprintf (stderr, "Host address length error (%d != %d).\n",
	       hp->h_length, ADDRBYTES);
#endif
      return (ERROR);
    }

  /* 09/26/2006 - modify for 64-bit */
   /* Before you go manipulating niceaddr[], make sure you've malloc'd
   * appropriate space. This wasn't needed before gcc v4.0, since we
   * could declare niceaddr[MAXHOSTS] at compile time, but now this
   * gives an "incomplete type" error and breaks. */
  if (!niceaddr)
    {
      niceaddr = (struct in_addr *) malloc (MAXHOSTS * sizeof (struct in_addr));
      for (i = 0; i < MAXHOSTS; i++)
	inet_aton (NULLHOST, &niceaddr[i]);
    }

  /* OK, got it. Now cache the FQDN in the hostname string and the
   * network address in the host table arrays if the address returned
   * is anything other than 127.0.0.1; if it is localaddr, then force
   * the name string to also read "localhost" regardless of whatever
   * name might be assigned locally. */
  memcpy (&niceaddr[hid], hp->h_addr, ADDRBYTES);
  if (ADDREQ(niceaddr[hid], localaddr))
    strncpy (hostname, "localhost", NICEMAXSTRINGLEN);
  else
    strncpy (hostname, hp->h_name, NICEMAXSTRINGLEN);

  /* Return success. */
  return (TRUE);
}

/**********************************************************************
 * Open a socket within a specified port address range and put it in
 * "listen" mode. Sets the port address actually used as a sideeffect.
 *
 * Returns socket descriptor (an integer >= 0) or ERROR.
 **********************************************************************/
int
watchSocket (unsigned short int *port, unsigned short int limit)
{
  int s, status;
  const int on = 1;
  struct sockaddr_in server;

  /* First make sure that the address structure is clear. */
  memset (&server, 0, sizeof (server));

  /* Create the socket. */
  if ((s = socket (AF_INET, SOCK_STREAM, 0)) < 0)
    {
#ifdef LIBDBG
      fprintf (stderr, "Socket creation error.\n");
#endif
      return (ERROR);
    }

  /* Set up the address structure and bind to the socket. */
  server.sin_family = AF_INET;

  /* Let the address be reused immediately should this process
   * die. This is useful especially for NICEPORT; otherwise, we have
   * to wait for the socket connection to time out. */
  setsockopt (s, SOL_SOCKET, SO_REUSEADDR, &on, sizeof (on));

  /* Now start sweeping the port range from *port all the way up to
   * limit looking for an address you can use. Note that if limit is
   * initially the same as *port, then that means we won't settle for
   * anything but the originally requested port number. Otherwise, we
   * start scanning at *port and end at limit. The port we actually
   * use is returned, indirectly, as the new value of *port. */
  status = ERROR;
  while ((*port <= limit) && (status == ERROR))
    {
      server.sin_port = htons (*port);	/* Port number */

      /* We need to cast the second argument to bind to avoid a
       * warning. We rely on the fact that sockaddr is used to type a
       * bunch of different types of address structures, including
       * TCP/IP's sockaddr_in. */
      status = bind (s, (struct sockaddr *) &server, sizeof (server));
      if (status == ERROR)
	(*port)++;
    }

  /* Finished the loop. Either you ran out of port addresses to try,
   * or you managed to bind the socket to one of them. The status
   * variable will tell you which is the case; the *port variable will
   * already be set appropriately. */
  if (status == ERROR)
    {
#ifdef LIBDBG
      fprintf (stderr, "Socket bind error.\n");
#endif
      close (s);
      return (ERROR);
    }

  /* Success. Put the socket in listen mode. SOMAXCONN limits the
   * number of connections, and is defined in <sys/socket.h>. */
  listen (s, SOMAXCONN);
  return (s);
}

/**********************************************************************
 * Open a socket on a specified host and set the appropriate socket
 * identifier. Returns TRUE if successful, FALSE if a "soft" error
 * occurred, or ERROR in the event of a "hard" or fatal failure.
 *
 * Note that the connect primitive may fail in several different
 * ways. The server might respond that there is no process waiting to
 * service a request, in which case the failure is a "hard"
 * failure. Alternatively, the server may timeout by failing to
 * respond within 75 seconds (this is apparently a kernel parameter),
 * either because the server is unreachable or because it is not
 * responding to the connect request; these are examples of "soft"
 * failures.
 **********************************************************************/
int
hailSocket (int hid, unsigned short int port, int *sock)
{
  struct sockaddr_in server;

  /* Create an Internet socket. */
  if ((*sock = socket (AF_INET, SOCK_STREAM, 0)) < 0)
    {
#ifdef LIBDBG
      fprintf (stderr, "Socket creation error (%d): %s.\n",
	       errno, strerror (errno));
#endif
      return (ERROR);
    }

  /* Make sure that the address structure is clear. */
  memset (&server, 0, sizeof (server));

  /* Set up the address structure for my destination host. */
  server.sin_family = AF_INET;
  server.sin_port = htons (port);
  memcpy (&server.sin_addr, &niceaddr[hid], ADDRBYTES);

  /* Try to connect. This can fail for a variety of reasons, some of
   * which (ETIMEDOUT or ENETUNREACH) has to do with the roughly 75
   * second timeout built into the connect primitive, while others are
   * more severe (like ECONNREFUSED). We rely on the fact that
   * sockaddr is used to type a bunch of different types of address
   * structures, including TCP/IP's sockaddr_in. */
  if (connect (*sock, (struct sockaddr *) &server, sizeof (server)) < 0)
    {
#ifdef LIBDBG
      fprintf (stderr, "Connect failure: %s.\n", strerror (errno));
#endif
      /* Close the failed socket connection. */
      close (*sock);
      if (errno == ETIMEDOUT || errno == ENETUNREACH)
	return (FALSE);
      else
	return (ERROR);
    }
  return (TRUE);
}

/**********************************************************************
 * Check socket you are watching. Return a socket that requires
 * service (a small nonegative integer), or ERROR if there is no
 * message waiting or if there is a socket error of some sort.
 *
 * Note by setting the select timeout interval to 0, we return
 * immediately if there is nothing waiting to be read.
 **********************************************************************/
int
checkSocket (int socket, struct in_addr *ipaddr, int timeout)
{
  fd_set ready;
  struct timeval to;
  struct sockaddr_in saddr;
  socklen_t saddrlen = sizeof(saddr);
  int result = 0;

  /* Zero out the bits corresponding to sockets in ready. */
  FD_ZERO (&ready);
  /* Add the socket you're interested in to ready set. */
  FD_SET (socket, &ready);
  to.tv_sec = timeout;
  to.tv_usec = 0;

  /* If select returns 0 then there are no pending messages; less than
   * 0 means there is a select error. */
  if ((result = select (socket + 1, &ready, 0, 0, &to)) <= 0)
    {
#ifdef LIBDBG
      /* Select error. */
      if (result < 0)
	fprintf (stderr, "Select error (%d) %s\n", errno, strerror (errno));
#endif
      return (ERROR);
    }

  /* Must be something on socket. Let's check to make sure socket is
   * in ready set (it should be) and then generate an appropriate
   * socket identifier. 
   *
   * We rely on the fact that sockaddr is used to type a bunch of
   * different types of address structures, including TCP/IP's
   * sockaddr_in. */
  if ((FD_ISSET (socket, &ready)) &&
      (result = accept (socket, (struct sockaddr *) &saddr, &saddrlen)) >= 0)
    {
      /* Since we don't want to rely on self-reported IP addresses, we're
       * going to get the source IP address from the network layer and set
       * it before returning. */
      memcpy (ipaddr, &(saddr.sin_addr), ADDRBYTES);
      return (result);
    }

#ifdef LIBDBG
  /* Accept error. */
  fprintf (stderr, "Accept error (%d) %s\n", errno, strerror (errno));
#endif
  return (ERROR);
}

/**********************************************************************
 * Socket interface. 
 *
 * As a debugging aid for protocol designers, any socket communication
 * that starts with the special string "%DBG" will be echoed if the
 * LIBDBG flag is set at compile time.
 **********************************************************************/
/**********************************************************************
 * sockprintf() returns the number of bytes written to socket or ERROR
 * if, for example, the other end of the socket has
 * disappeared. Avoids use of XDR by writing everything to a string
 * (where the string functions as the external representation
 * language). The message buffer must have been declared to be
 * NICEMAXSTRINGLEN+1 bytes long to accomodate string terminator.
 **********************************************************************/
int
sockprintf (int socket, char *message, char *format, ...)
{
  int length, value;
  va_list argp;
#if FALSE
  /* Zero out the message string. */
  memset (message, NICEEOR, NICEMAXSTRINGLEN + 1);
#endif
  /* Initialize argp and construct the appropriate message. Note that
   * vsnprintf will ensure that a terminating byte is written to
   * message, and as long as NICEEOR is the same as the string
   * terminator you're kosher. */
  va_start (argp, format);
  length = vsnprintf (message, NICEMAXSTRINGLEN, format, argp);
  va_end (argp);

  /* Write message line to socket. Keep trying until you send some
   * bytes or you get an unrecoverable error (interrupts are the only
   * recoverable errors considered here). Make sure you explicitly
   * send the terminating byte (hence 1 more than length bytes). */
  while ((value = send (socket, message, (length+1), 0)) < 0)
    if (errno != EINTR)
      {
#ifdef LIBDBG
	/* An error here usually means the other end of the socket has
	 * evaporated. */
	fprintf (stderr, "Socket (%d) write error: %d (%d: %s) %s\n",
		 socket, value, errno, strerror (errno), message);
#endif
	return (ERROR);
      }
#if FALSE
  /* Issue an error if you didn't get the whole message sent,
   * including the string terminator. Not actually needed, since send
   * will never send only part of a message. */
  if (value != length+1)
      {
#ifdef LIBDBG
	fprintf (stderr, "Socket (%d) send error (%d != %d) %s\n",
		 socket, (value-1), length, message);
#endif
	return (ERROR);
      }
#endif
#ifdef LIBDBG
  if (!strncmp (message, "%DBG", 4))	/* Remove? */
    fprintf (stderr, "\tsend (%d:%d) ``%s''\n", socket, length, message);
#elif FALSE
  fprintf (stderr, "\tsend (%d:%d) ``%s''\n", socket, length, message);
#endif

  /* Return number of bytes written, excluding the string
   * terminator. */
  return (length);
}

/**********************************************************************
 * sockscanf() returns the number of arguments bound from the socket
 * or ERROR if, for example, the other end of the socket has
 * disappeared.  Avoids use of XDR by writing everything to a string
 * (where the string functions as the external representation
 * language). The message buffer must have been declared to be
 * NICEMAXSTRINGLEN+1 bytes long to accomodate the (explicitly
 * transmitted) string terminator.
 **********************************************************************/
int
sockscanf (int socket, char *message, char *format, ...)
{
  int value = 0;
  char *msgptr = &message[0];
  va_list argp;
#if FALSE
  /* Zero out the message string. */
  memset (message, (char) NULL, NICEMAXSTRINGLEN + 1);
#endif
  /* Read the pending line from socket. Note that you may need to make
   * multiple calls to recv in order to get the entire message. You'll
   * never get more than the message from the stream, however, because
   * we always alternate message and response, so any legitimate
   * correspondent will never stream more than one message to you.
   * 
   * So: keep going until you get a string terminating byte, unless
   * you get an unrecoverable error (the only recoverable error here
   * is an interrupt). You might get an EOF (recv returns 0), which is
   * OK if the socket closes after the message, but not OK if the
   * message is not properly terminated. */
  while ((value = recv (socket, msgptr, (NICEMAXSTRINGLEN - (msgptr - message)), 0)) != 0)
    if (value > 0)
      {
	/* Advance msgptr to end of current message chunk. */
	msgptr = msgptr + value;

	/* Quit recv'ing if you got the string terminator. */
	if (*(msgptr - 1) == NICEEOR)
	  break;
      }
    else if (errno != EINTR)
      {
#ifdef LIBDBG
	/* Unrecoverable error. */
	fprintf (stderr, "Socket (%d) read error: %d (%d: %s).\n",
		 socket, (msgptr - 1 - message), errno, strerror (errno));
#endif
	return (ERROR);
      }

  /* OK, you've fallen out of the loop. Make sure that you've got a
   * non-zero-length message with a terminating byte at the end; if
   * not, generate an error. Careful: you initialized message to all
   * NULLs, so you need to be sure that the transmitted string
   * terminator was indeed received in the byte before msgptr. Won't
   * work if NICEEOR is different than the normal null byte string
   * terminator. */
  if ((value = (msgptr - 1 - message)) == 0 || *(msgptr - 1) != NICEEOR)
    {
#ifdef LIBDBG
      /* An error here usually means the other end of the socket has
       * evaporated. */
      fprintf (stderr, "Socket %d (%d bytes) read error: bad message.\n", 
	       socket, value);
#endif
      return (ERROR);
    }

#ifdef LIBDBG
  if (!strncmp (message, "%DBG", 4))	/* Remove? */
    fprintf (stderr, "\trecv (%d:%d) ``%s''\n", socket, value, message);
#elif FALSE
  fprintf (stderr, "\trecv (%d:%d) ``%s''\n", socket, value, message);
#endif

  /* Initialize argp and parse the message. */
  va_start (argp, format);
  value = vsscanf (message, format, argp);
  va_end (argp);

  /* Return number of items parsed. */
  return (value);
}
