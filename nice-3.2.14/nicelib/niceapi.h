/**********************************************************************
 * NICE System Library
 * Alberto Maria Segre
 * 
 * Copyright (c) 1997-2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software 
 * for non-profit educational purposes only.
 **********************************************************************/
#include <stdarg.h>		/* Needed for va_list etc. */

/* Maximum number of protocols in an application. */
#define NICEMAXQID 6		/* Max number of messages handled. */

/**********************************************************************
 * clientFn and serverFn are used to type all the protocol
 * functions. The "client" side is always the initiating end, and
 * always starts with a transmit (rather than a receive).
 **********************************************************************/
/* The clientFn is the initiating party, usually the nagger. Signature
 * is function(socket, swap, message, argp). */
typedef int (*clientFn) (int, int, char *, va_list);
/* The serverFn is the responding party, usually the master. Signature
 * is function(nid, socket, swap, message). */
typedef void (*serverFn) (int, int, int, char *);

/**********************************************************************
 * Function prototypes.
 **********************************************************************/
/* Used at the application level. */
int niceInit (int argc, char *argv[], clientFn clientFn, serverFn serverFn);
void niceHandler (int qid, serverFn serverFn, int synchronous, 
		  int safety, int rootOnly,
		  int epochInc, int epochEqp);
void niceCheck ();
int niceTalk (int serverid, int qid, clientFn clientFn, ...);

void niceExit ();
int niceSlave ();
int niceMaster ();
int niceLeaf ();
int niceRoot ();
int niceParallel ();
int niceSerial ();

/* Used in niceapi and at the application level, but not in niced or
 * niceq because of endianess issues. Instead, niced/niceq use
 * sockprintf/sockscanf exclusively. */
int sockwrite (void *data, size_t size, size_t count, int socket, int swap, char *message);
int sockread (void *data, size_t size, size_t count, int socket, int swap, char *message);
