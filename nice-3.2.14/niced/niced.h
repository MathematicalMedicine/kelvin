/**********************************************************************
 * NICE System Daemon
 * Alberto Maria Segre
 * 
 * Copyright (c) 1997-2005, The University of Iowa.  All rights
 * reserved.  Permission is hereby given to use and reproduce this
 * software for non-profit educational purposes only.
 **********************************************************************/
#include "niceaux.h"		/* Needed for alarms. */
#include "nicecom.h"		/* Needed for NICEMAXSTRINGLEN. */
#include "../config.h"
#include <stdio.h>		/* Needed for printf */
#include <stdlib.h>		/* EXIT_SUCCESS, EXIT_ERROR */
#include <stdarg.h>		/* vsscanf */
#include <signal.h>		/* Needed for signal */
#include <limits.h>		/* For USHRT_MAX */
#include <string.h>		/* Needed for strchr, strcmp */
#include <unistd.h>		/* Needed for gethostname */
#include <time.h>		/* For time, localtime */
#include <netdb.h>		/* Needed for gethostbyaddr */
#include <sys/resource.h>	/* For getpriority */
#include <sys/stat.h>		/* For umask */
#include <sys/wait.h>		/* For waitpid, WNOHANG */
#include <netinet/in.h>		/* Needed for in_addr */
#include <arpa/inet.h>		/* Needed for inet_aton */
#include <sys/socket.h>		/* For AF_INET */
#include <errno.h>		/* ETIMEDOUT, ENETUNREACH */

/* Time intervals. */
#define PRIORITYDFLT 0		/* Lower scheduling priority for daemon. */
#if FALSE		/* Alternate values for debugging. */
#define LOADINT 1		/* Benchmark interval = 1m */
#define CDRPINT 150		/* Child cull interval (>=2.5*loadint), 150s = 2.5m */
#define PDRPINT 270		/* Parent cull interval (>=4.5*loadint), 270s = 4.5m */
#define RECONFODDS 2		/* Reconfiguration odds (on avg once every 2m). */
#else			/* Normal values. */
#define LOADINT 9		/* Benchmark interval = 9m */
#define CDRPINT 1351		/* Child cull interval (>2.5*loadint), ~22.5m */
#define PDRPINT 2431		/* Parent cull interval (>4.5*loadint), ~40.5m */
#define RECONFODDS 13		/* Reconfiguration odds (on avg once every ~2h). */
#endif
#define ACTIVE 0		/* Daemon mode. */
#define ORPHAN 1

/* Socket timeouts. */
#define SLEEPINT 5		/* Timeout for socket monitoring (secs). */

/* Benchmark. */
#define BENCHMARKSZE 500	/* Size of benchmark problem. */
#define BENCHMARKREP 10		/* Number of benchmark cycles. */

/**********************************************************************
 * Time conversion (to "minutes after midnight" form). We'll store all
 * times in terms of local time, but will always communicate times to
 * other machines in terms of universal time. This way we can be
 * assured of correct behavior even as the machine's time changes wrt
 * to universal time (e.g., daylight savings time).
 **********************************************************************/
#define MAXTIME 1439		/* (24*60)-1 */
#define TIME(hours,minutes) ((hours)*60+(minutes))
#define UTIME(ltime) (((ltime) + utcf) % MAXTIME)
#define LTIME(utime) (((utime) - utcf) % MAXTIME)

/**********************************************************************
 * A bunch of global variables. These need to be global so they can be
 * accessed by the signal handler that traps kill or interrupt
 * signals; but the bottom line is that it is simpler to make these
 * global so we don't need to pass them around all over the place.
 **********************************************************************/
extern int mode;		/* Daemon mode; one of ACTIVE or ORPHAN. */
extern int maxChildren;
extern int cardinality;
extern int maxSlaves;
extern int utcf;		/* Universal time correction (minutes). */
extern int debug;		/* Governs daemonization and output behavior. */
extern unsigned short port;	/* My port number. */
extern int priority;		/* Daemon scheduling priority. */

/* The following variables are indexed by hid, like niceaddr[] from
 * nicecom.h. All values in nicestart[], niceend[] and nicelast[] are
 * expressed in local time, but are communicated in universal time. */
extern int niceload[MAXHOSTS];	/* Load on known daemons. */
extern int nicewall[MAXHOSTS];	/* Segregated subdomain. */
extern int nicedepth[MAXHOSTS];	/* Depth of highest opening in subtree. */
extern int nicestart[MAXHOSTS];	/* Scheduling information. */
extern int niceend[MAXHOSTS];
extern time_t nicelast[MAXHOSTS];	/* Last contact time. */
extern int nicebfact;  		/* Hierarchy's global branching factor. */

/* Variables describing running slave applications. */
extern int slaves;			/* Number of local slave processes. */
extern int spids[MAXSLAVES];		/* Slave PIDs. */

/* Temporary hostent used to construct error while debugging. */
struct hostent *errhp;

#ifndef NOALARM
/* Used to timeout either side of a blocked daemon-to-damon
 * protocol. We divide this into two parts: on the server side, the
 * time to authenticate, which should be short to guard against
 * denial-of-service attacks (which our single-threaded server cannot
 * really hope to resist) and the time to complete a protocol, which
 * should be longer (as in niceapi.c). On the client side, we have the
 * time to initially connect, which should be short to avoid any
 * possibility of deadlock, and the time to complete the
 * daemon-to-daemon protocol. We also want to have a (very long)
 * timeout interval for file transfers, so that they are simply not
 * left hanging. */
#define AUTHINT 5               /* Authentication timeout, in sec. */
#define CONNINT AUTHINT		/* Initial connection timeout, in sec. */
#define TXFRINT	600		/* Hung file transfer threshold, in sec. */
#endif

/**********************************************************************
 * NASTI stands for NICE Application Signing, Transfer and
 * Installation.  It adds application downloading and signing
 * functionality to NICE daemon, and is based on design and prototype
 * impelementations by Mahesh Murthy, Shashank Khandelwal, and Yan
 * Liu.
 **********************************************************************/

/* Note: MAXAPPLISTSIZE > MAXTXFRS + 1 required for transfers to
 * clear. */
#define MAXTXFRS 10		/* Maximum concurrent file transfers. */
#define MAXAPPLISTSIZE 20	/* Maximum size for file transfer lists. */

/* List cell structure. Each item retains the application name
 * (including version and architecture information), as well as the
 * PID of the transferring process (if active) and properties of the
 * transfer request (e.g., whether to launch on completion). Note that
 * some applications on these lists won't be binary compatible with
 * the current machine; rather, they're being requested so that they
 * can be propagated down to other machines below. */
typedef struct getWish
{
  char name[NICEMAXSTRINGLEN + 1];	/* application-version-architecture */
  char addr[NICEMAXSTRINGLEN + 1];	/* IP of requesting parent. */
  char port[NICEMAXSTRINGLEN + 1];	/* Port of requesting parent. */
  short pid;				/* PID of transfer, if active. */
  struct getWish *next;
} GetWish;

typedef struct putWish
{
  char name[NICEMAXSTRINGLEN + 1];	/* application-version-architecture */
  int child;				/* Index requesting child. */
  short pid;				/* PID of transfer, if active. */
  struct putWish *next;
} PutWish;

/**********************************************************************
 * Application wishlists.
 **********************************************************************/
extern GetWish *getWishes;
extern PutWish *putWishes;
extern GetWish *freeGetListCells;
extern PutWish *freePutListCells;

/**********************************************************************
 * Function prototypes.
 **********************************************************************/
void daemonize ();
void dieNeatly (int signum);
void die (int restructure, int notifyParent);
int pauseSlaves ();
int prodSlaves ();
void dropChild (int i);
int bestMatch (int start, int end, int depth, int load, int noBarriers,
	       int highOvlp, int highDepth, int highLoad);
#ifndef NODYNAMIC
int electChild ();
#endif
int establishConnection (int hid, int qid, char *message, int *socket);

int niceConfig  (int argc, char *argv[], char *parent, int *priority);
int niceRegister (int hid);
int niceUnregister ();
void niceOrphan (int child);
void niceLoad ();
#ifndef NODYNAMIC
int nicePromote (int child, int parent);
#endif
void niceSpawn (int child, char *executable, unsigned short int cport);
int niceTxfr (PutWish *cell, char *filename, int length);

void serveNiceRegister (int socket, struct in_addr caddr);
void serveNiceUnregister (int socket, struct in_addr caddr);
void serveNiceOrphan (int socket, struct in_addr caddr);
void serveNiceLoad (int socket, struct in_addr caddr);
#ifndef NODYNAMIC
void serveNicePromote (int socket, struct in_addr caddr);
#endif
void serveNiceSpawn (int socket, struct in_addr caddr, int avail);
int spawnApplications (char *filename, int count, int fertile, 
		       struct in_addr maddr, char *mport);
#ifndef NOALARM
void serveNiceTxfr (int socket, struct in_addr caddr, int avail, int alarm);
#else
void serveNiceTxfr (int socket, struct in_addr caddr, int avail);
#endif
void serveNiceQuery (int socket, int avail);
void serveNicePing (int socket, int avail);
void serveNiceSolicit (int socket, int avail);
void serveNiceExit (int socket);

long int bench ();

/**********************************************************************
 * NASTI function prototypes.
 **********************************************************************/
void initWishes ();
int cullWishes ();
GetWish* addGetWish (char *name, char *addr, char *port);
PutWish* addPutWish (char *name, int child);
GetWish* findGetWish (char *name);
PutWish* findPutWish (char *name, int child);
int getListEmpty ();
int countPendingGets ();
char* findNthPendingGet (int n);
void startPendingPuts (int count);
void killActiveTransfers ();

void niceApplicationSource (PutWish *cell, unsigned short int port, int length);
void niceApplicationSink (GetWish *cell, int socket, int length);

int niceCheckApplication (char *application, char *filename);
int niceRegisterApplication (char *name);
