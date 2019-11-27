/**********************************************************************
 * NICE Application Signing, Transfer and Installation
 * Alberto Maria Segre
 *
 * Based on a design and prototype impelementations by Mahesh Murthy,
 * Shashank Khandelwal, and Yan Liu.
 *
 * Copyright (c) 1997-2005, The University of Iowa.  All rights
 * reserved.  Permission is hereby given to use and reproduce this
 * software for non-profit educational purposes only.
 **********************************************************************/
#include "niced.h"

/**********************************************************************
 * Adds application downloading and signing functionality to the NICE
 * daemon.
 * 
 * Downloading is managed by two lists, getWishes and putWishes, which
 * contain information about the application transfer requests, or
 * "wishes." The niced acts on the contents of the lists once at the
 * end of each main loop, just after culling slave PIDs.  When an
 * application is available for download and there are not too many
 * active transfer processes, the parent initiates the download by
 * first forking a new file transfer process (simply to run
 * niceApplicationSource ()) and then sending a niceTxfr () message to
 * its intended child notifying it of the parameters of the forked
 * source process so that the child can fork a corresponding
 * niceApplicationSink () process.
 *
 * Download requests occur when a child does not have the application
 * it is asked to spawn (serveNiceSpawn). The child makes a wish entry
 * on its getWishes, and the parent makes an analogous wish entry on
 * its putWishes. Each wish is tagged with appropriate properties,
 * such as whether to launch the application once downloading is
 * complete, or, in the case of putWishes, whether the application is
 * available to download.
 *
 * During niceLoad exchanges, wishlist elements that are not currently
 * being downloaded are reported from child to parent to see if the
 * parent has since obtained the desired executable and should
 * therefore initiate a tranfer.
 *
 * Note that these lists are of limited size, so some applications
 * that never become available will eventually fall off the end of the
 * list. Also note that the number of active downloads is also
 * limited.
 **********************************************************************/

/**********************************************************************
 * Wishlists. 
 **********************************************************************/
GetWish *getWishes = NULL;
PutWish *putWishes = NULL;
GetWish *freeGetWishes = NULL;
PutWish *freePutWishes = NULL;

/**********************************************************************
 * Wishlist mechanism. Used by the niced to keep track of applications
 * it needs to download from its parent niced. Each entry is marked as
 * currently being downloaded and whether or not to launch on
 * completion.
 **********************************************************************/
void
initWishes ()
{
  int i;
  GetWish *g;
  PutWish *p;

  /* Generate and place MAXAPPLISTSIZE cells on each free list and
   * initialize their next pointers. Since we limit the maximum size
   * of the application wishlists, this can be done once and for all
   * at initialization time; there's never any need to realloc since
   * the cells are recycled as we go. */
  freeGetWishes = g = (GetWish *) malloc (sizeof (GetWish));
  for (i = 0; i < MAXAPPLISTSIZE - 1; i++, g = g->next)
    g->next = (GetWish *) malloc (sizeof (GetWish));
  freePutWishes = p = (PutWish *) malloc (sizeof (PutWish));
  for (i = 0; i < MAXAPPLISTSIZE - 1; i++, p = p->next)
    p->next = (PutWish *) malloc (sizeof (PutWish));
}

/**********************************************************************
 * Cull any completed transfer processes from the current wish lists,
 * and return how many remain active. Recall inactive cells have pid
 * of 0, while active cells have a real pid. When an active transfer
 * is completed, and its pid is culled, set the cell's pid to ERROR so
 * the cell can be recycled.
 **********************************************************************/
int
cullWishes ()
{
  int count = 0;
  GetWish *gtmp = getWishes, *gtmp2;
  PutWish *ptmp = putWishes, *ptmp2;
  char filename[NICEMAXSTRINGLEN + 1];
  struct in_addr addr;

  /* Scan the get list and try culling any active PIDs first. Anything
   * we successfully cull with waitpid we mark with an ERROR and then
   * recycle the cell in a later pass. */
  while (gtmp)
    {
      if (gtmp->pid > 0)
	{
	  if (waitpid (gtmp->pid, NULL, WNOHANG) == gtmp->pid)
	    {
	      /* Log output. */
	      niceLog (NLOGMSG, "Scheduled GET of %s complete.\n", gtmp->name);

	      /* Check that the application is good to go, and derive
	       * the apporpriate filename while you're at it. */
	      if (niceCheckApplication (gtmp->name, filename))
		{
		  /* Mark cell for recycling. */
		  gtmp->pid = ERROR;
	      
		  /* Spawn a copy of the application if port and addr are
		   * set (to indicate we wanted to spawn once the download
		   * was complete). */
		  if (gtmp->addr != "" && gtmp->port != "")
		    {
		      inet_aton (gtmp->addr, &addr);
		      if (spawnApplications (filename, 1, TRUE, addr, gtmp->port) > 0)
			niceLog (NLOGMSG, "Belatedly spawning %s (port %s).\n", 
				 gtmp->name, gtmp->port);
		      else
			niceLog (NLOGWRN, "Failed to spawn %s.\n", gtmp->name);
		    }
		}
	      else
		{
		  /* Something went wrong with the transfer. Don't
		   * cull the cell, but leave it so we can try again
		   * later. Remove any information regarding belated
		   * spawning since it will surely be too late for
		   * that (really? should check). */
		  memset (gtmp->addr, (char) NULL, NICEMAXSTRINGLEN + 1);
		  memset (gtmp->port, (char) NULL, NICEMAXSTRINGLEN + 1);
		  gtmp->pid = 0;
		}
	    }
	  else
	    count++;
	}
      gtmp = gtmp->next;
    }

  /* Scan again to recycle the list cells. We could do a better job of
   * this, but the wish lists are likely to be so short (at most one
   * or two elements, and usually empty) that it really doesn't
   * matter. */
  if ((gtmp = getWishes))
    {
      while (gtmp->next)
	if (gtmp->next->pid == ERROR)
	  {
	    gtmp2 = freeGetWishes;
	    freeGetWishes = gtmp->next;
	    gtmp->next = gtmp->next->next;
	    freeGetWishes->next = gtmp2;
	  }
	else
	  gtmp = gtmp->next;

      /* Check the head element too. */
      if (getWishes->pid == ERROR)
	{
	  gtmp = getWishes->next;
	  getWishes->next = freeGetWishes;
	  freeGetWishes = getWishes;
	  getWishes = gtmp;
	}
    }

  /* Do the same with the put list: cull completed transfers and count
   * still-active transfers. */
  while (ptmp)
    {
      if (ptmp->pid > 0)
	{
	  if (waitpid (ptmp->pid, NULL, WNOHANG) == ptmp->pid)
	    {
	      /* Log output. */
	      niceLog (NLOGMSG, "Scheduled PUT of %s to child #%d complete.\n", 
		       ptmp->name, ptmp->child);
	      ptmp->pid = ERROR;
	    }
	  else
	    count++;
	}
      ptmp = ptmp->next;
    }

  /* Scan again to recycle the list cells. We could do a better job of
   * this, but the wish lists are likely to be so short (at most one
   * or two elements, and usually empty) that it really doesn't
   * matter. */
  if ((ptmp = putWishes))
    {
      while (ptmp->next)
	if (ptmp->next->pid == ERROR)
	  {
	    ptmp2 = freePutWishes;
	    freePutWishes = ptmp->next;
	    ptmp->next = ptmp->next->next;
	    freePutWishes->next = ptmp2;
	  }
	else
	  ptmp = ptmp->next;

      /* Check the head element too. */
      if (putWishes->pid == ERROR)
	{
	  ptmp = putWishes->next;
	  putWishes->next = freePutWishes;
	  freePutWishes = putWishes;
	  putWishes = ptmp;
	}
    }

  /* Return number of still-active transfers found. */
  return (count);
}

/**********************************************************************
 * Add a wish list item to the get list. Returns pointer to cell if
 * successful, else NULL. If port and addr are provided, spawn the
 * application once it is downloaded.
 **********************************************************************/
GetWish*
addGetWish (char *name, char *addr, char *port)
{
  int i;
  GetWish *tmp;

  /* First, make sure all is in order: if some required arguments are
   * missing, fail and return NULL. */
  if (!name)
    return (NULL);

  /* Next, make sure that the cell is not already on the list. */
  if ((tmp = findGetWish (name)))
    {
      /* Found it! Reset addr and port on the premise that this is a
       * more recent request, but only if the requested transfer is
       * not already underway. */
      if ((tmp->pid == 0) && addr && port)
	{
	  strncpy (tmp->addr, addr, NICEMAXSTRINGLEN);
	  strncpy (tmp->port, port, NICEMAXSTRINGLEN);
	}
      return (tmp);
    }

  /* Push item onto the get wish list. If there are no more free get
   * cells, we need to recycle one of the older requests (but be
   * careful not to recycle active cells, i.e., those with ongoing
   * downloads). This ensures that requests for applications that
   * never become available for download from the parent are
   * eventually flushed.  */
  if (!freeGetWishes)
    {
      /* Find oldest (i.e., last) list cell that is not currently
       * active and recycle it. Since the tail elements of the list
       * might all be active, we need to keep track of how many
       * inactive cells we see during the scan so that we can cull the
       * last one.
       *
       * Note: getWishes is not NULL since free list is empty.
       * Note: i > 1 since MAXAPPLISTSIZE > MAXTXFRS+1. */
      tmp = getWishes;
      i = MAXAPPLISTSIZE - MAXTXFRS;	/* How many inactive cells to expect. */
      while (tmp->next->next)
	{
	  if (tmp->pid == 0)
	    i--;
	  if (i == 1 && tmp->next->pid == 0)
	    break;		/* Next cell is last inactive cell! */
	  tmp = tmp->next;
	}

      /* Recycle it. */
      freeGetWishes = tmp->next;
      tmp->next = freeGetWishes->next;
      freeGetWishes->next = NULL;
      memset (freeGetWishes->addr, (char) NULL, NICEMAXSTRINGLEN + 1);
      memset (freeGetWishes->port, (char) NULL, NICEMAXSTRINGLEN + 1);
      freeGetWishes->pid = 0;
    }

  /* Push a new (possibly just recovered) cell onto the get wish
   * list. */
  tmp = freeGetWishes;
  freeGetWishes = tmp->next;
  tmp->next = getWishes;
  getWishes = tmp;

  /* Populate the new cell with the appropriate information. */
  strncpy (tmp->name, name, NICEMAXSTRINGLEN);
  if (addr && port)
    {
      strncpy (tmp->addr, addr, NICEMAXSTRINGLEN);
      strncpy (tmp->port, port, NICEMAXSTRINGLEN);
    }
  tmp->pid = 0;

  /* Return succes. */
  return (tmp);
}

/**********************************************************************
 * Add a wish list item to the put list. Returns pointer to newly
 * added cell if successful, else NULL.
 **********************************************************************/
PutWish*
addPutWish (char *name, int child)
{
  int i;
  PutWish *tmp;

  /* First, make sure all is in order: if some required arguments are
   * NULL or out of bounds, fail and return NULL. */
  if (!name || child < 1 || child > nicehcnt + 1)
    return (NULL);

  /* Next, make sure that the cell is not already on the list. */
  if ((tmp = findPutWish (name, child)))
    {
      /* Found it! No need to add another copy. */
      return (tmp);
    }

  /* Push item onto the put wish list. If there are no more free put
   * cells, we need to recycle one of the older requests (but be
   * careful not to recycle active cells, i.e., those with ongoing
   * downloads). This ensures that requests for applications that
   * never become available for download from the parent are
   * eventually flushed.  */
  if (!freePutWishes)
    {
      /* Find oldest (i.e., last) list cell that is not currently
       * active and recycle it. Since the tail elements of the list
       * might all be active, we need to keep track of how many
       * inactive cells we see during the scan so that we can cull the
       * last one.
       *
       * Note: putWishes is not NULL since free list is empty.
       * Note: i > 1 since MAXAPPLISTSIZE > MAXTXFRS+1. */
      tmp = putWishes;
      i = MAXAPPLISTSIZE - MAXTXFRS;	/* How many inactive cells to expect. */
      while (tmp->next->next)
	{
	  if (tmp->pid == 0)
	    i--;
	  if (i == 1 && tmp->next->pid == 0)
	    break;		/* Next cell is last inactive cell! */
	  tmp = tmp->next;
	}

      /* Recycle it. */
      freePutWishes = tmp->next;
      tmp->next = freePutWishes->next;
      freePutWishes->next = NULL;
      freePutWishes->child = 0;
      freePutWishes->pid = 0;
    }

  /* Push a new (possibly just recovered) cell onto the put wish
   * list. */
  tmp = freePutWishes;
  freePutWishes = tmp->next;
  tmp->next = putWishes;
  putWishes = tmp;

  /* Populate the new cell with the appropriate information. */
  strncpy (tmp->name, name, NICEMAXSTRINGLEN);
  tmp->child = child;
  tmp->pid = 0;

  /* Return succes. */
  return (tmp);
}

/**********************************************************************
 * findGetWish () returns pointer to a get list cell indexed by
 * name.
 **********************************************************************/
GetWish*
findGetWish (char *name)
{
  GetWish *tmp = getWishes;
  while (tmp)
    if (strncmp (tmp->name, name, NICEMAXSTRINGLEN) == 0)
      return (tmp);
    else
      tmp = tmp->next;
  return (NULL);
}

/**********************************************************************
 * findPutWish () returns pointer to a put list cell indexed by
 * name and child.
 **********************************************************************/
PutWish*
findPutWish (char *name, int child)
{
  PutWish *tmp = putWishes;
  while (tmp)
    if ((strncmp (tmp->name, name, NICEMAXSTRINGLEN) == 0)
	&& tmp->child == child)
      return (tmp);
    else
      tmp = tmp->next;
  return (NULL);
}

/**********************************************************************
 * countPendingGets () returns how many gets are pending.
 **********************************************************************/
int
countPendingGets ()
{
  int count = 0;
  GetWish *tmp = getWishes;

  /* Step through the get list and count entries with 0 PID; these are
   * the currently inactive requests for downloads. */
  while (tmp)
    {
      if (tmp->pid == 0)
	count++;
      tmp = tmp->next;
    }
  return (count);
}

/**********************************************************************
 * findNthPendingGet () returns name of Nth pending get (0-indexed) in
 * wish list. Pretty silly way to do things, actually, but since the
 * get list is guaranteed to be short, it's probably OK.
 **********************************************************************/
char*
findNthPendingGet (int n)
{
  GetWish *tmp = getWishes;

  /* Step through the get list and return the Nth entry (0-indexed)
   * whose PID is 0; these are the currently inactive requests for
   * downloads. */
  while (tmp && n >= 0)
    {
      if (tmp->pid == 0 && n > 0)
	n--;
      else if (tmp->pid == 0)
	return (tmp->name);
      tmp = tmp->next;
    }
  return (NULL);
}

/**********************************************************************
 * startPendingPuts (int count) starts count file transfers (parent
 * initiated).
 **********************************************************************/
void
startPendingPuts (int count)
{
  PutWish *tmp = putWishes;
  struct stat info;
  char filename[NICEMAXSTRINGLEN + 1];

  /* Step through the put list and fire up transfers for the first
   * count cells with PID set to 0; these are currently inactive
   * entries. Once you start the transfer, update the cell and the
   * count of active transfers. */
  while (tmp && count > 0)
    {
      if (tmp->pid == 0)
	{
	  /* Got one to start. First, we need to determine the source
	   * filename and length. */
	  snprintf (filename, NICEMAXSTRINGLEN, "%s/%s", NBINDIR, tmp->name);

	  /* If successful in forking the put, decrement your license
	   * to start more puts by 1. Forked transfers, once complete,
	   * will be culled later, in the main loop, by PID which is
	   * stored in the putWish. */
	  if ((stat (filename, &info) != ERROR) &&
	      (niceTxfr (tmp, filename, (int) info.st_size) == TRUE))
	    count--;
	}
      tmp = tmp->next;
    }
}

/**********************************************************************
 * Kill any running file transfers on your way out the door. Invoked
 * at niced exit time.
 **********************************************************************/
void
killActiveTransfers ()
{
  GetWish *gtmp = getWishes;
  PutWish *ptmp = putWishes;
  while (gtmp)
    {
      if (gtmp->pid != 0)
	{
	  kill (gtmp->pid, SIGQUIT);
	  /* waitpid (gtmp->pid, NULL, WNOHANG); */
	}
      gtmp = gtmp->next;
    }
  while (ptmp)
    {
      if (ptmp->pid != 0)
	{
	  kill (ptmp->pid, SIGQUIT);
	  /* waitpid (ptmp->pid, NULL, WNOHANG); */
	}
      ptmp = ptmp->next;
    }
}

/**********************************************************************
 * Scrutinize application, e.g., prior to launching or to check if you
 * have a version for download. Instantiates local application name,
 * and checks the local copy to ensure it is properly signed. Takes
 * the application name (including version and architecture
 * information, but not including local pathname) and, if successful,
 * produces the fully-instantiated filename.
 *
 * Returns TRUE on success, FALSE if executable is not found or
 * improperly signed (TODO: should also mark such an exectuable for
 * replacement or deletion). ERROR is returned if the requested
 * executable is illegal or malformed (e.g., has a path component).
 **********************************************************************/
int
niceCheckApplication (char *application, char *filename)
{
  /* Derive the appropriate filename from the name of the
   * executable. Also trap and reject any funny path business. */
  if (strrchr (application, '/'))
    {
      /* Someone tried to place a path structure. This should never
       * happen for a legitimate application. */
      niceLog (NLOGERR, "Illegal executable %s requested; fail.\n", application);
      return (ERROR);
    }

  /* Complete filename with full pathname. */
  snprintf (filename, NICEMAXSTRINGLEN, "%s/%s", NBINDIR, application);

  /* Make sure the executable exists locally. */
  if (access (filename, X_OK) == ERROR)
    {
      /* No such exectuable. */
      niceLog (NLOGWRN, "Executable %s not found.\n", application);
      return (FALSE);
    }

  /* TODO: Here is where we should be checking, e.g., the local
   * application's signed certificate, and returning FALSE if it
   * doesn't pass muster. */

  /* Good to go. Local filename stored in filename. */
  return (TRUE);
}

/**********************************************************************
 * Register a new NICE application, checking that it is safe to run
 * and setting permissions appropriately.
 **********************************************************************/
int
niceRegisterApplication (char *filename)
{
  /* TODO: Here is where we should be checking, e.g., the local
   * application's signed certificate, and returning doing something
   * with the bad application if doesn't pass muster. 
   *
   * For now, just mark the file as readable and executable; i.e.,
   * assume its OK. */
  chmod (filename, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
  return (TRUE);
}

/**********************************************************************
 * Send application executable to my child.  This is the process
 * forked to fulfill a put wish; it gets the port number of the
 * child's analagous stub transfer process along with the putWish cell
 * containing the name of the file.
 **********************************************************************/
void
niceApplicationSource (PutWish *cell, unsigned short int port, int length)
{
  int i, j, socket;
#ifndef NOALARM
  int sourceAlarm;
#endif
  char filename[NICEMAXSTRINGLEN + 1];
  char message[NICEMAXSTRINGLEN + 1];
  FILE *file = NULL;

#ifndef NOALARM
  /* This is a lot like establishConnection () except we are
   * connecting to a different port.
   * 
   * Set a timer in case something blocks while attempting to
   * connect. The initial connection should be relatively fast, like
   * when you authenticate, so a relatively short timeout is
   * appropriate. */
  if ((sourceAlarm = niceSetTimeout (CONNINT)) == FALSE)
    {
      /* The other end of the interaction has timed out. You shouldn't
       * need to clear the timer, which we know has already expired,
       * since that's how we got here. What we should do is close the
       * socket and quit. */
      niceLog (NLOGERR, "Timeout [%ds] in niceApplicationSource().\n", CONNINT);
      return;
    }
#endif

  /* Contact your child at the specified port. Note that hailSocket()
   * may return ERROR (fatal or "hard" failure, e.g., connection
   * refused) or FALSE ("soft" failure, e.g., host busy) but that, in
   * either case, we'll treat this as if the connection was simply
   * refused by the target. */
  if (hailSocket (cell->child, port, &socket) != TRUE)
    {
      /* Punt. */
      niceLog (NLOGERR, "Can't contact sink in niceApplicationSource().\n");
      return;
    }

#ifndef NOALARM
  /* OK so far. Cancel the short alarm and reset it to a longer
   * interval for the file transfer. */
  niceClearAlarm (sourceAlarm);
#endif

  /* Open the actual application executable; first, you need to find
   * it in the appropriate directory. */
  snprintf (filename, NICEMAXSTRINGLEN, "%s/%s", NBINDIR, cell->name);
  if ((file = fopen (filename, "rb")) == NULL)
    {
      /* Punt. */
      niceLog (NLOGERR, "Can't open %s niceApplicationSource().\n", cell->name);
      return;
    }

#ifndef NOALARM
  /* Set a timer in case something blocks during the execution of the
   * protocol. A protocol is declared hung/blocked and attempts to
   * read or write are abandoned when it hangs for TXFRINT seconds,
   * typically a (very) large number. TODO: This should probably NOT
   * be a constant, but rather should be set based on some factor of
   * the observed RTT.
   *
   * Question: how do these work given that we're operating in a
   * forked process? Do we have our parent's alarms too? */
  if ((sourceAlarm = niceSetTimeout (TXFRINT)) == FALSE)
    {
      /* The other end of the interaction has blocked. You shouldn't
       * need to clear the timer, which we know has already expired,
       * since that's how we got here. */
      fclose (file);
      close (socket);
      niceLog (NLOGERR, "Timeout [%ds] in niceApplicationSource().\n", TXFRINT);
      return;
    }
#endif

  /* Start sending the application. This is a lot like
   * sockread/sockwrite, except we don't worry about endianess since
   * an architecture's own binaries will always have the proper
   * endianess. */
  i = 0;
  while (i < length)
    {
      /* Read next NICEMAXSTRINGLEN bytes from file. */
      if ((j = (int) fread (message, sizeof (char), NICEMAXSTRINGLEN, file)) < 0)
	{
	  /* Punt. */
	  fclose (file);
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (sourceAlarm);
#endif
	  close (socket);
	  niceLog (NLOGERR, "File read failure in niceApplicationSource().\n");
	  return;
	}

      /* Send next chunk of bytes over the socket. */
      while (write (socket, message, j) < 0)
	{
	  /* Punt. */
	  fclose (file);
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (sourceAlarm);
#endif
	  close (socket);
	  niceLog (NLOGERR, "Socket write failure in niceApplicationSource().\n");
	  return;
	}

      /* Increment count of bytes sent. */
      i += j;
      /* fprintf (stderr, "%d+ ", j); */
    }
  /* Close application. */
  fclose (file);
#ifndef NOALARM
  /* Clear the timer. */
  niceClearAlarm (sourceAlarm);
#endif
  close (socket);
  return;
}

/**********************************************************************
 * Receive the application executable from my parent. 
 **********************************************************************/
#ifdef NOALARM
#define TXFRINT	600	     /* Hung file transfer threshold, in sec. */
#endif
void
niceApplicationSink (GetWish *cell, int socket, int length)
{
  int i, j, tsocket;
#ifndef NOALARM
  int sinkAlarm;
#endif
  struct in_addr caddr;
  char filename[NICEMAXSTRINGLEN + 1];
  char message[NICEMAXSTRINGLEN + 1];
  FILE *file = NULL;

  /* Compute eventual name of application file, and open for
   * write. Recall the name I was given does not include the pathname,
   * since that is installation-dependent. */
  snprintf (filename, NICEMAXSTRINGLEN, "%s/%s", NBINDIR, cell->name);
  if ((file = fopen (filename, "wb")) == NULL)
    {
      /* Can't open file for write. Punt. */
      close (socket);
      niceLog (NLOGERR, "Can't open %s for write in niceApplicationSink().\n", cell->name);
      return;
    }

  /* Wait, a very long time, for an incoming message on the control
   * socket. Returns the transfer socket. If we time out here, just
   * give it up. */
  if ((tsocket = checkSocket (socket, &caddr, TXFRINT)) == ERROR)
    {
      /* Socket error. Give it up. */
      fclose (file);
      close (socket);
      niceLog (NLOGERR, "Request timeout [%ds] in niceApplicationSink().\n", TXFRINT);
      return;
    }
  /* We're done with the control socket, since we are only expecting a
   * single transfer request. */
  close (socket);

#ifndef NOALARM
  /* OK, we have a transfer socket, we've closed the control socket,
   * and are ready to transfer the executable.  Set a timer in case
   * something blocks during the execution of the protocol. A protocol
   * is declared hung/blocked and attempts to read or write are
   * abandoned when it hangs for TXFRINT seconds, typically a (very)
   * large number. TODO: This should probably NOT be a constant, but
   * rather should be set based on some factor of the observed RTT.
   *
   * Question: how do these work given that we're operating in a
   * forked process? Do we have our parent's alarms too? */
  if ((sinkAlarm = niceSetTimeout (TXFRINT)) == FALSE)
    {
      /* The other end of the interaction has blocked. You shouldn't
       * need to clear the timer, which we know has already expired,
       * since that's how we got here. */
      fclose (file);
      close (tsocket);
      niceLog (NLOGERR, "Timeout [%ds] in niceApplicationSink().\n", TXFRINT);
      return;
    }
#endif

  /* Start receiving the application. This is a lot like
   * sockread/sockwrite, except we don't worry about endianess since
   * an architectures own binaries will always have the proper
   * endianess, and since we are just buffering to a file control is a
   * bit simpler. 
   *
   * Also, by using read () we ensure that we never receive more than
   * NICEMAXSTRINGLEN bytes, although we may well receive fewer bytes
   * since TCP/IP may fragment things. */
  i = 0;
  while ((j = read (tsocket, message, NICEMAXSTRINGLEN)) != 0)
    {
      if (j > 0)
	{
	  /* Got some. Write it to disk, then go around again. */
	  if (fwrite (message, sizeof (char), j, file) != j)
	    {
	      /* Failed to write: punt. */
	      fclose (file);
#ifndef NOALARM
	      /* Clear the timer. */
	      niceClearAlarm (sinkAlarm);
#endif
	      close (tsocket);
	      niceLog (NLOGERR, "File write failure in niceApplicationSource().\n");
	      return;
	    }
	  i += j;
	  /* fprintf (stderr, "%d", j); */

	  /* Are we done? */
	  if (i >= length)
	    break;
	}
      else if (errno != EINTR && errno != EAGAIN)
	{
	  /* Bad reception: punt. */
	  fclose (file);
#ifndef NOALARM
	  /* Clear the timer. */
	  niceClearAlarm (sinkAlarm);
#endif
	  close (tsocket);
	  niceLog (NLOGERR, "Socket read failure in niceApplicationSink().\n");
	  return;
	}
      /* fprintf (stderr, "+ "); */
    }
  
  /* OK, done with file transfer. */
  fclose (file);
#ifndef NOALARM
  /* Clear the timer. */
  niceClearAlarm (sinkAlarm);
#endif
  close (tsocket);

  /* Before returning, register the application just received. This
   * involves setting the executable bits and so on, but should also
   * do much the same security check as niceCheckApplication (). */
  niceRegisterApplication (filename);
  return;
}

