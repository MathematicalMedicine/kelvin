/**********************************************************************
 * NICE Communications Library
 * Alberto Maria Segre
 * 
 * Copyright (c) 1997-2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software 
 * for non-profit educational purposes only.
 **********************************************************************/
#include <sys/types.h>		/* For size_t and more. */
#include <netinet/in.h>
#include <sys/socket.h>		/* Needed for socket, accept, bind, listen, etc */
#include <arpa/inet.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

#define MAXPORT USHRT_MAX	/* Largest allowable port number. */

#define NULLHOST "0.0.0.0"	/* Used to indicate a null host. */
#define LOCALHOST "127.0.0.1"	/* Used to indicate a null host. */
#define ADDRBYTES (sizeof (struct in_addr)) /* IP address bytes. */
#define ADDRCMP(ADDR1,ADDR2) (memcmp((ADDR1),(ADDR2),ADDRBYTES))
#define ADDREQ(ADDR1,ADDR2) (ADDRCMP(&(ADDR1),&(ADDR2)) == 0)
#define ADDRNEQ(ADDR1,ADDR2) (ADDRCMP(&(ADDR1),&(ADDR2)) != 0)
#define ADDRNULL(ADDR1) (ADDREQ(ADDR1,nulladdr))

#define NICEMAXSTRINGLEN 1024	/* Max chars in filenames, hostnames, etc. */
#define NICESAFEARGS    "%1024s" /* Used to input strings w/o overflow. */
#define NICESAFEARGHD   "%5hd"   /* Used to input strings w/o overflow. */
#define NICESAFEARGD    "%10d"   /* Used to input strings w/o overflow. */
#define NICESAFEARGLG   "%30lg"  /* Used to input strings w/o overflow. */
#define EOSTRING(STRING) (strchr((STRING), (char) NULL))	/* End of string pointer. */

/**********************************************************************
 * Some global variables that must be accessible to a broad range of
 * functions. These constitute our local network host table.
 **********************************************************************/
extern char nicehost[];		/* My host name. */
extern int nicesock;		/* My incoming messages. */
extern int nicehcnt;		/* Number of hosts (< MAXHOSTS). */
extern struct in_addr *niceaddr;	/* Parent, child host addresses. */
extern struct in_addr nulladdr;	/* Null address (for comparison). */
extern struct in_addr localaddr;	/* Localhost address (for comparisons). */

#define MAXHOSTS 20		/* Max hosts in host table. */
#define PARENT 0
#define SELF (MAXHOSTS-1)
#define ROOT (MAXHOSTS-2)	/* Used in niced only. */
/* While we don't use MAXHOSTS-2 in niceapi, directly, some
 * applications (e.g., mlip) may elect to use this location to store
 * root processor information (e.g., for sending results to be stored
 * to disk). Better stay safe and always use MAXCHILDREN as the upper
 * limit for either child niced's or child applications. */
#define MAXCHILDREN (MAXHOSTS-3)

/**********************************************************************
 * Used to communicate NICE application parameters.
 **********************************************************************/
#define NICEFLAG "-Nice"	/* Command-line argument string. */
#define NICELEAF "-NiceLeaf"	/* Command-line leaf designation. */

/**********************************************************************
 * NICE message strings. These are included here largely for
 * convenience, but also because many of them are used in
 * daemon/daemon, control/daemon, application/daemon, as well as
 * application/application contexts.
 **********************************************************************/
#define NICEQRYSTRD	"%%QRY %d"
#define NICEQRYSTRDS	"%%QRY %d %s"
#define NICEQRYSTRDD	"%%QRY %d %d"
#define NICEQRYSTRDSS	"%%QRY %d %s %s"
#define NICEQRYSTRDDDSS	"%%QRY %d %d %d %s %s"
#define NICEQRYSTRDDDD	"%%QRY %d %d %d %d"
#define NICEQRYSTRDDDDD	"%%QRY %d %d %d %d %d"
#define NICEQRYSTRSS	"%%QRY %s %s"

#define NICEACKSTRD	"%%ACK %d"
#define NICEACKSTRDD	"%%ACK %d %d"
#define NICEACKSTRDS	"%%ACK %d %s"

#define NICEINFSTRD	"%%INF %d"
#define NICEINFSTRDD	"%%INF %d %d"
#define NICEINFSTRDDD	"%%INF %d %d %d"
#define NICEINFSTRDDDD	"%%INF %d %d %d %d"
#define NICEINFSTRDDDDD	"%%INF %d %d %d %d %d"
#define NICEINFSTRDDDDDD "%%INF %d %d %d %d %d %d"

#define NICEINFSTRS	"%%INF %s"
#define NICEINFSTRSD	"%%INF %s %d"
#define NICEINFSTRSDD	"%%INF %s %d %d"
#define NICEINFSTRSDDDD "%%INF %s %d %d %d %d"
#define NICEINFSTRSDDDDD "%%INF %s %d %d %d %d %d"
#define NICEINFSTRSDDDDDD "%%INF %s %d %d %d %d %d %d"
#define NICEINFSTRSSD	"%%INF %s %s %d"
#define NICEINFSTRSSDD	"%%INF %s %s %d %d"
#define NICEINFSTRSS	"%%INF %s %s"

#define NICEINFSTRHU	"%%INF %hu"
#define NICEINFSTRSHU	"%%INF %s %hu"
#define NICEINFSTRSSSHU	"%%INF %s %s %s %hu"
#define NICEINFSTRSHUDD	"%%INF %s %hu %d %d"
#define NICEINFSTRSHUDDD "%%INF %s %hu %d %d %d"

/**********************************************************************
 * Daemon query identifiers. These are the identifiers used when
 * communicating with the nice daemon. They're here in nicecom.h
 * because they need to be the same in niced, niceq, as well as any
 * applications that talk to the daemon. Actual choice of numbers is
 * pretty much irrelevant, provided they don't conflict with the other
 * identifiers; we note here that the NICE API uses small ints, so
 * we'll pick numbers over 100 just to avoid a conflict.
 **********************************************************************/
/* Used only by niced. */
#define NQIDREGISTER	100
#define NQIDUNREGISTER	101
#define NQIDORPHAN	102
#define NQIDLOAD	103
#define NQIDPROMOTE	104
#define NQIDSPAWN	105
#define NQIDTXFR	106

/* Used only by niceq. */
#define NQIDQUERY	110

/* Used only by applications. */
#define NQIDPING	120
#define NQIDSOLICIT	121
#define NQIDEXIT	122

/* End of record marker; null character in a string. */
#define NICEEOR (char) '\0'

/**********************************************************************
 * Function prototypes.
 **********************************************************************/
/* Used in niceapi and niced/niceq only. */
int qualifyHostname (char *hostname, int hid);
int watchSocket (unsigned short int *port, unsigned short int limit);
int hailSocket (int hid, unsigned short int port, int *socket);
int checkSocket (int socket, struct in_addr *ipaddr, int timeout);

/* Used in niceapi, niced/niceq, and also at the application level. */
int sockprintf (int socket, char *message, char *format, ...);
int sockscanf (int socket, char *message, char *format, ...);
