#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "../utils/sw.h"
#include "../utils/utils.h"

#include "StudyDB.h"

#ifndef MAIN
extern 
#endif
struct StudyDB studyDB;

int dBInitNotDone = TRUE, dBStmtsNotReady = TRUE;

void initializeDB () {

  /* Initialize structure for a MySQL connection. */
  if ((studyDB.connection = mysql_init(NULL)) == NULL)
    FATAL("Cannot initialize MySQL (%s)", strerror(errno));

  /* Connect. */
  if (!mysql_real_connect(studyDB.connection, studyDB.hostname, studyDB.username, studyDB.password, 
			  NULL, 0, NULL, CLIENT_MULTI_RESULTS /* Important discovery here */))
    ERROR("Cannot connect to MySQL on hostname [%s] as username [%s/%s] (%s)", studyDB.hostname, 
	  studyDB.username, studyDB.password, mysql_error(studyDB.connection));

  /* Change database. */
  if (mysql_select_db(studyDB.connection, studyDB.dBName))
    ERROR("Cannot change MySQL db (%s)", mysql_error(studyDB.connection));

  /* Verify our studyId. */
  sprintf (studyDB.strAdhocStatement, "Select Description from Studies where StudyId = %d", studyDB.inStudyId);
  if (mysql_query (studyDB.connection, studyDB.strAdhocStatement))
    ERROR("Cannot select study information (%s:%s)", studyDB.strAdhocStatement, mysql_error(studyDB.connection));
  if ((studyDB.resultSet = mysql_store_result (studyDB.connection)) == NULL)
    ERROR("Cannot retrieve study information (%s)", mysql_error(studyDB.connection));
  if (mysql_num_rows (studyDB.resultSet) == 0)
    ERROR("Study %d not found", studyDB.inStudyId);
  else {
    if ((studyDB.row = mysql_fetch_row (studyDB.resultSet)) == NULL)
      ERROR("Cannot fetch study information (%s)", mysql_error(studyDB.connection));
    INFO ("Storing/retrieving results under study %d (%s)", studyDB.inStudyId, studyDB.row[0]);
  }
  dBInitNotDone = FALSE;
}

// Don't use this willy-nilly since it does allocate space for the null and sign indicators.
#define BINDNUMERIC(WHERE, WHAT, TYPE) {	 \
  WHERE.buffer_type = TYPE; \
  WHERE.buffer = &(WHAT);		 \
  WHERE.buffer_length = sizeof (WHAT); \
  WHERE.length = 0; \
  WHERE.is_null = (my_bool *) calloc (1, sizeof (my_bool));	\
  WHERE.is_unsigned = 0; \
}

void prepareDBStatements () {

  if (dBInitNotDone)
    initializeDB ();

  // Prepare the GetDLOD call
  studyDB.stmtGetDLOD = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetDLOD, 0, sizeof(studyDB.bindGetDLOD));

  BINDNUMERIC (studyDB.bindGetDLOD[0], studyDB.inPedPosId, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetDLOD[1], studyDB.inAlpha, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[2], studyDB.inDGF, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[3], studyDB.inLC1BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[4], studyDB.inLC1BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[5], studyDB.inLC1LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[6], studyDB.inLC1LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[7], studyDB.inLC2BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[8], studyDB.inLC2BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[9], studyDB.inLC2LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[10], studyDB.inLC2LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[11], studyDB.inLC3BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[12], studyDB.inLC3BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[13], studyDB.inLC3LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[14], studyDB.inLC3LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[15], studyDB.inRegionId, MYSQL_TYPE_LONG);

  strncpy (studyDB.strGetDLOD, "call GetDLOD (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,@outRegionId,@outLOD)", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtGetDLOD, studyDB.strGetDLOD, strlen (studyDB.strGetDLOD)))
    ERROR("Cannot prepare GetDLOD call statement (%s)", mysql_error(studyDB.connection));
  if (mysql_stmt_bind_param (studyDB.stmtGetDLOD, studyDB.bindGetDLOD))
    ERROR("Cannot bind GetDLOD call statement (%s)", mysql_error(studyDB.connection));

  // Prepare the GetDLOD results call
  studyDB.stmtGetDLODResults = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetDLODResults, 0, sizeof(studyDB.bindGetDLODResults));

  BINDNUMERIC (studyDB.bindGetDLODResults[0], studyDB.outRegionId, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetDLODResults[1], studyDB.outLOD, MYSQL_TYPE_DOUBLE);

  strncpy (studyDB.strGetDLODResults, "Select @outRegionId, @outLOD", MAXSTMTLEN-1);
  if (mysql_stmt_prepare (studyDB.stmtGetDLODResults, studyDB.strGetDLODResults, strlen (studyDB.strGetDLODResults)))
    ERROR("Cannot prepare GetDLOD results select statement (%s)", mysql_error(studyDB.connection));
  if (mysql_stmt_bind_result (studyDB.stmtGetDLODResults, studyDB.bindGetDLODResults))
    ERROR("Cannot bind GetDLOD results select statement (%s)", mysql_error(studyDB.connection));
  dBStmtsNotReady = FALSE;
}

double GetDLOD (int inPedPosId, double inAlpha, double inDGF,
	      double inLC1BigPen, double inLC1BigLittlePen, double inLC1LittleBigPen, double inLC1LittlePen,
	      double inLC2BigPen, double inLC2BigLittlePen, double inLC2LittleBigPen, double inLC2LittlePen,
	      double inLC3BigPen, double inLC3BigLittlePen, double inLC3LittleBigPen, double inLC3LittlePen,
	      int inRegionId)
{
  int rows, columns;

  if (dBStmtsNotReady)
    prepareDBStatements ();

  studyDB.inPedPosId = inPedPosId;
  studyDB.inAlpha = inAlpha;
  studyDB.inDGF = inDGF;
  studyDB.inLC1BigPen = inLC1BigPen;
  studyDB.inLC1BigLittlePen = inLC1BigLittlePen;
  studyDB.inLC1LittleBigPen = inLC1LittleBigPen;
  studyDB.inLC1LittlePen = inLC1LittlePen;
  studyDB.inLC2BigPen = inLC2BigPen;
  studyDB.inLC2BigLittlePen = inLC2BigLittlePen;
  studyDB.inLC2LittleBigPen = inLC2LittleBigPen;
  studyDB.inLC2LittlePen = inLC2LittlePen;
  studyDB.inLC3BigPen = inLC3BigPen;
  studyDB.inLC3BigLittlePen = inLC3BigLittlePen;
  studyDB.inLC3LittleBigPen = inLC3LittleBigPen;
  studyDB.inLC3LittlePen = inLC3LittlePen;
  studyDB.inRegionId = inRegionId;

  if (mysql_stmt_execute (studyDB.stmtGetDLOD))
    ERROR("Cannot execute GetDLOD call statement (%s, %s)", 
	  mysql_stmt_error(studyDB.stmtGetDLOD), mysql_stmt_sqlstate(studyDB.stmtGetDLOD));

  if (mysql_stmt_execute (studyDB.stmtGetDLODResults))
    ERROR("Cannot execute GetDLOD results select statement (%s, %s)", 
	  mysql_stmt_error(studyDB.stmtGetDLODResults), mysql_stmt_sqlstate(studyDB.stmtGetDLODResults));
  if (mysql_stmt_store_result (studyDB.stmtGetDLODResults) != 0)
    ERROR("Cannot retrieve adhoc information (%s)", mysql_stmt_error(studyDB.stmtGetDLODResults));
  if (mysql_stmt_fetch (studyDB.stmtGetDLODResults) != 0)
    ERROR("Cannot fetch results (%s)", mysql_stmt_error(studyDB.stmtGetDLODResults));

  if (*studyDB.bindGetDLODResults[1].is_null) {
    INFO ("In RegionId %d, LOD is NULL", studyDB.outRegionId);
    return -1LL;
  } else {
    INFO ("In RegionId %d, LOD is %G", studyDB.outRegionId, studyDB.outLOD);
    return studyDB.outLOD;
  }
}

#ifdef MAIN

/*

gcc -g -o test databaseSupport.c ../utils/libklvnutls.a -DMAIN -lmysqlclient -I/usr/local/mysql/include -L/usr/local/mysql/lib/ 
./test 19 localhost LOD_dev whv001 foobar

*/

#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {
  
  studyDB.inStudyId = atoi (argv[1]);
  strcpy (studyDB.hostname, argv[2]);
  strcpy (studyDB.dBName, argv[3]);
  strcpy (studyDB.username, argv[4]);
  strcpy (studyDB.password, argv[5]);

  GetDLOD (10, .2, .3, .71, .42, .44, .13, .71, .42, .44, .13, .71, .42, .46, .13, 1);
  GetDLOD (10, .25, .3, .71, .42, .42, .13, .71, .42, .45, .13, .71, .42, .45, .13, 1);
  GetDLOD (17, .2, .3, .71, .42, .44, .13, .71, .42, .42, .13, .71, .42, .42, .13, 1);
  GetDLOD (9, .2, .3, .71, .42, .44, .13, .71, .42, .42, .13, .71, .42, .42, .13, 1);
  GetDLOD (6, .2, .3, .71, .42, .44, .13, -1, 0, 0, 0, -1, 0, 0, 0, 1);
  GetDLOD (7, .2, .3, .71, .42, .44, .13, .71, .42, .42, .13, -1, 0, 0, 0, 1);
  return EXIT_SUCCESS;
}
#endif
