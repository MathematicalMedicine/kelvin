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

void initializeDB () {

  /* Initialize structure for a MySQL connection. */
  if ((studyDB.connection = mysql_init(NULL)) == NULL)
    FATAL("Cannot initialize MySQL (%s)", strerror(errno));

  /* Connect. */
  if (!mysql_real_connect(studyDB.connection, studyDB.hostname, studyDB.username, studyDB.password, NULL, 0, NULL, 0))
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
}

#define BINDLONG(WHERE, WHAT) { \
  WHERE.buffer_type = MYSQL_TYPE_LONG; \
  WHERE.buffer = &(WHAT);	       \
  WHERE.buffer_length = sizeof (WHAT); \
  WHERE.length = 0; \
  WHERE.is_null = (my_bool *) 0; \
  WHERE.is_unsigned = 0; \
}
#define BINDDOUBLE(WHERE, WHAT) { \
  WHERE.buffer_type = MYSQL_TYPE_DOUBLE; \
  WHERE.buffer = &(WHAT);		 \
  WHERE.buffer_length = sizeof (WHAT); \
  WHERE.length = 0; \
  WHERE.is_null = (my_bool *) 0; \
  WHERE.is_unsigned = 0; \
}

void *prepareStatements () {

  // Prepare the GetDLOD call
  studyDB.stmtGetDLOD = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetDLOD, 0, sizeof(studyDB.bindGetDLOD));

  BINDLONG (studyDB.bindGetDLOD[0], studyDB.inPedPosId);
  BINDDOUBLE (studyDB.bindGetDLOD[1], studyDB.inAlpha);
  BINDDOUBLE (studyDB.bindGetDLOD[2], studyDB.inDGF);
  BINDDOUBLE (studyDB.bindGetDLOD[3], studyDB.inLC1BigPen);
  BINDDOUBLE (studyDB.bindGetDLOD[4], studyDB.inLC1BigLittlePen);
  BINDDOUBLE (studyDB.bindGetDLOD[5], studyDB.inLC1LittleBigPen);
  BINDDOUBLE (studyDB.bindGetDLOD[6], studyDB.inLC1LittlePen);
  BINDDOUBLE (studyDB.bindGetDLOD[7], studyDB.inLC2BigPen);
  BINDDOUBLE (studyDB.bindGetDLOD[8], studyDB.inLC2BigLittlePen);
  BINDDOUBLE (studyDB.bindGetDLOD[9], studyDB.inLC2LittleBigPen);
  BINDDOUBLE (studyDB.bindGetDLOD[10], studyDB.inLC2LittlePen);
  BINDDOUBLE (studyDB.bindGetDLOD[11], studyDB.inLC3BigPen);
  BINDDOUBLE (studyDB.bindGetDLOD[12], studyDB.inLC3BigLittlePen);
  BINDDOUBLE (studyDB.bindGetDLOD[13], studyDB.inLC3LittleBigPen);
  BINDDOUBLE (studyDB.bindGetDLOD[14], studyDB.inLC3LittlePen);
  BINDLONG (studyDB.bindGetDLOD[15], studyDB.inRegionId);

  strncpy (studyDB.strGetDLOD, "call GetDLOD (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtGetDLOD, studyDB.strGetDLOD, strlen (studyDB.strGetDLOD)))
    ERROR("Cannot prepare GetDLOD call statement (%s)", mysql_error(studyDB.connection));
  if (mysql_stmt_bind_param (studyDB.stmtGetDLOD, studyDB.bindGetDLOD))
    ERROR("Cannot bind GetDLODs call statement (%s)", mysql_error(studyDB.connection));
  fprintf (stderr, "BOUND!\n");
}

void GetDLOD (int inPedPosId, double inAlpha, double inDGF,
	      double inLC1BigPen, double inLC1BigLittlePen, double inLC1LittleBigPen, double inLC1LittlePen,
	      double inLC2BigPen, double inLC2BigLittlePen, double inLC2LittleBigPen, double inLC2LittlePen,
	      double inLC3BigPen, double inLC3BigLittlePen, double inLC3LittleBigPen, double inLC3LittlePen,
	      int inRegionId)
{
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
    ERROR("Cannot execute GetDLOD call statement (%s), mysql_stmt_error(studyDB.stmtGetDLOD)");

  if ((studyDB.resultSet = mysql_store_result (studyDB.connection)) == NULL)
    ERROR("Cannot retrieve LOD (%s)", mysql_error(studyDB.connection));
  if (mysql_num_rows (studyDB.resultSet) == 0)
    ERROR("LOD %d not found");
  else {
    if ((studyDB.row = mysql_fetch_row (studyDB.resultSet)) == NULL)
      ERROR("Cannot fetch study information (%s)", mysql_error(studyDB.connection));
    INFO ("Got LOD %G (%s)", studyDB.row[0]);
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

  initializeDB ();
  prepareStatements ();
  GetDLOD (7, .2, .3, .71, .42, .45, .13, .71, .42, .42, .13, .71, .42, .42, .13, 1);
}
#endif
