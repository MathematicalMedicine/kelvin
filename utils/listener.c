/*
$Id$

by Bill Valentine-Cooper for the Center for Quantitative and Computational
Biology, Nationwide Children's Hospital Research Institute.

  1. Read the command line for <DB host name> <DB username> <DB password> <UDP socket>

  2. connect to database identified by command line parameters. DB should have
  a table RunLog in it as follows:

  Create table RunLog (fromNode varchar(32), changeDate timestamp, logEntry varchar(255));

  3. connect as a listener to UDP socket identified by command line parameter.

  4. As packets arrive on the socket, display them on stdout and log them to the database.

  5. Reconnect to DB or socket as needed.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
*/
#include <arpa/inet.h>
#include <netdb.h>
#include "sw.h"
#include "mysql.h"

char logEntry[MAXUDPMSG], fromNode[128];
unsigned long logEntryLength, fromNodeLength;
char dbHostName[64], dbName[64], dbUsername[64], dbPassword[64];
int listenerSocket;

MYSQL mysql;
MYSQL_STMT *insertRunLogStmt;

#define TRUE (0==0)
#define FALSE (!TRUE)

int doMySQLStuff() {

  /* Connect. */
  if (!mysql_real_connect(&mysql, dbHostName, dbUsername, dbPassword, NULL, 0, NULL, 0)) {
    fprintf(stderr, " Cannot connect to MySQL host [%s] (%s)\n", dbHostName, mysql_error(&mysql));
    return FALSE;
  }
  /* Change database. */
  if (mysql_select_db(&mysql, dbName)) {
    fprintf(stderr, "Cannot change MySQL db (%s)\n", mysql_error(&mysql));
    exit(EXIT_FAILURE);
  }
  /* Prepare a statement. */
  char *insertRunLog = "insert into RunLog (fromNode, logEntry) values (?, ?)";
  MYSQL_BIND boundValues[2];
  if (!(insertRunLogStmt = mysql_stmt_init(&mysql))) {
    fprintf(stderr, " Cannot initialize MySQL_STMT (%s)\n", mysql_error(&mysql));
    exit(EXIT_FAILURE);
  }    
  if (mysql_stmt_prepare(insertRunLogStmt, insertRunLog, strlen(insertRunLog))) {
    fprintf(stderr, " Cannot prepare MySQL_STMT (%s, %s)\n", insertRunLog, mysql_error(&mysql));
    exit(EXIT_FAILURE);
  }    
  /* Bind variables to the statement. */
  memset(boundValues, 0, sizeof(boundValues));
  boundValues[0].buffer_type = MYSQL_TYPE_STRING;
  boundValues[1].buffer_type = MYSQL_TYPE_STRING;
  boundValues[0].buffer = (char *) fromNode;
  boundValues[1].buffer = (char *) logEntry;
  boundValues[0].buffer_length = MAXUDPMSG;
  boundValues[1].buffer_length = MAXUDPMSG;
  boundValues[0].is_null = 0;
  boundValues[1].is_null = 0;
  boundValues[0].length = &fromNodeLength;
  boundValues[1].length = &logEntryLength;
  if (mysql_stmt_bind_param(insertRunLogStmt, boundValues)) {
    fprintf(stderr, " Cannot bind values to MySQL_STMT (%s %s)\n", insertRunLog, mysql_error(&mysql));
    exit(EXIT_FAILURE);
  }
  return TRUE;
}

int retryIntervals[] = {2, 5, 10, 30, 60, 300, 900, 3600};
#define MAXRETRYINTERVALS 8

void connectToDbOrSnooze() {
  int i, j, k, l, m, n;		/* Sure sign of a FORTRAN programmer */
  /* Get to the database or snooze a bit. */
  i = 0;
  fprintf(stderr, "(Re)connect to database...\n");
  while (doMySQLStuff() == FALSE) {
    fprintf(stderr, "Cannot connect and/or prepare, sleeping for %d seconds...\n",
	    retryIntervals[i]);
    sleep(retryIntervals[i]);
    if (++i >= MAXRETRYINTERVALS) i--;
  }
  fprintf(stderr, "OK\n");
}

int sockfd;
struct sockaddr_in my_addr;// my address information
struct sockaddr_in their_addr; // connector's address information
socklen_t addr_len;
int numbytes;

int doSocketListenerStuff() {
  if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) == -1) {
    perror("socket");
    exit(EXIT_FAILURE);
  }

  my_addr.sin_family = AF_INET; // host byte order
  my_addr.sin_port = htons(listenerSocket); // short, network byte order
  my_addr.sin_addr.s_addr = INADDR_ANY; // automatically fill with my IP
  memset(my_addr.sin_zero, '\0', sizeof my_addr.sin_zero);

  if (bind(sockfd, (struct sockaddr *)&my_addr, sizeof my_addr) == -1) {
    perror("bind");
    return FALSE;
  }
  addr_len = sizeof their_addr;
  return TRUE;
}

void connectToSocketOrSnooze() {
  int i, j, k, l, m, n;		/* Sure sign of a FORTRAN programmer */
  /* Get to the database or snooze a bit. */
  i = 0;
  fprintf(stderr, "(Re)connect to socket...\n");
  while (doSocketListenerStuff() == FALSE) {
    fprintf(stderr, "Cannot connect and/or prepare, sleeping for %d seconds...\n",
	    retryIntervals[i]);
    sleep(retryIntervals[i]);
    if (++i >= MAXRETRYINTERVALS) i--;
  }
  fprintf(stderr, "OK\n");
}

int main(int argc, char *argv[])
{
  struct hostent *he;
  /* Get command line parameters. */

  char *usage = "%s <DB host name> <DB name> <DB username> <DB password> <UDP socket, e.g. 4950>\n";
  if (argc<6) {
    fprintf(stderr, "%d is too few arguments!\n", argc+1);
    fprintf(stderr, usage, argv[0]);
    exit(EXIT_FAILURE);
  } else {
    strcpy(dbHostName, argv[1]);
    strcpy(dbName, argv[2]);
    strcpy(dbUsername, argv[3]);
    strcpy(dbPassword, argv[4]);
    if ((listenerSocket = atoi(argv[5])) == 0) {
      fprintf(stderr, "%d is not a valid UDP socket number\n", listenerSocket);
      fprintf(stderr, usage, argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  printf("Listening on socket %d and connecting to DB host %s, database %s as %s\n",
	 listenerSocket, dbHostName, dbName, dbUsername);

  /* Initialize structure for a MySQL connection. */
  if (!mysql_init(&mysql)) {
    perror("Cannot initialize MySQL\n");
    exit(EXIT_FAILURE);
  }
  connectToDbOrSnooze();
  connectToSocketOrSnooze();

  while (TRUE) {
    while ((numbytes = recvfrom(sockfd, logEntry, MAXUDPMSG-1 , 0,
			     (struct sockaddr *)&their_addr, &addr_len)) == -1) {
      connectToSocketOrSnooze();
    }

    logEntry[numbytes] = '\0';

    he = gethostbyaddr ((struct in_addr *)&their_addr.sin_addr, addr_len, AF_INET);
    if (he == NULL) /* Get the host info */
      strcpy(fromNode, inet_ntoa(their_addr.sin_addr));
    else
      strcpy(fromNode, he->h_name);
    printf("(%s): %s", fromNode, logEntry);

    logEntryLength = strlen(logEntry);
    fromNodeLength = strlen(fromNode);
    while (mysql_stmt_execute(insertRunLogStmt)) {
      connectToDbOrSnooze();
    }
  }

  /* We'd do these if we ever got here! */
  mysql_close(&mysql);
  close(sockfd);

  return 0;
}
