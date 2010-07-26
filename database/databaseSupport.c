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

// Don't use these willy-nilly since they do allocate space for length, is_null and error.
#define BINDNUMERIC(WHERE, WHAT, TYPE) {	 \
  WHERE.buffer_type = TYPE; \
  WHERE.buffer = &(WHAT);		 \
  WHERE.length = (unsigned long *) calloc (1, sizeof (unsigned long)); \
  WHERE.is_null = (my_bool *) calloc (1, sizeof (my_bool)); \
  WHERE.error = (my_bool *) calloc (1, sizeof (my_bool)); \
  WHERE.is_unsigned = 0; \
}
#define BINDSTRING(WHERE, WHAT, SIZE) {	 \
  WHERE.buffer_type = MYSQL_TYPE_STRING; \
  WHERE.buffer = &(WHAT);		 \
  WHERE.buffer_length = SIZE; \
  WHERE.length = (unsigned long *) calloc (1, sizeof (unsigned long)); \
  WHERE.is_null = (my_bool *) calloc (1, sizeof (my_bool)); \
  WHERE.error = (my_bool *) calloc (1, sizeof (my_bool)); \
  WHERE.is_unsigned = 0; \
}

void prepareDBStatements () {

  if (dBInitNotDone)
    initializeDB ();

  // Prepare the select for the pedigree/position
  studyDB.stmtGetPedPosId = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetPedPosId, 0, sizeof(studyDB.bindGetPedPosId));

  BINDNUMERIC (studyDB.bindGetPedPosId[0], studyDB.inStudyId, MYSQL_TYPE_LONG);
  BINDSTRING (studyDB.bindGetPedPosId[1], studyDB.inPedigreeSId, sizeof (studyDB.inPedigreeSId));
  BINDNUMERIC (studyDB.bindGetPedPosId[2], studyDB.inChromosomeNo, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetPedPosId[3], studyDB.inRefTraitPosCM, MYSQL_TYPE_DOUBLE);

  strncpy (studyDB.strGetPedPosId, "Select PedPosId from PedigreePositions where StudyId = ? AND PedigreeSId = ? AND ChromosomeNo = ? AND RefTraitPosCM = ?", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtGetPedPosId, studyDB.strGetPedPosId, strlen (studyDB.strGetPedPosId)))
    ERROR("Cannot prepare GetPedPosId call statement (%s)", mysql_error(studyDB.connection));
  if (mysql_stmt_bind_param (studyDB.stmtGetPedPosId, studyDB.bindGetPedPosId))
    ERROR("Cannot bind GetPedPosId call statement (%s)", mysql_error(studyDB.connection));

  BINDNUMERIC (studyDB.bindGetPedPosIdResults[0], studyDB.outPedPosId, MYSQL_TYPE_LONG);

  if (mysql_stmt_bind_result (studyDB.stmtGetPedPosId, studyDB.bindGetPedPosIdResults))
    ERROR("Cannot bind GetPedPosId results (%s)", mysql_error(studyDB.connection));

  // Prepare the GetDLOD call
  studyDB.stmtGetDLOD = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetDLOD, 0, sizeof(studyDB.bindGetDLOD));

  BINDNUMERIC (studyDB.bindGetDLOD[0], studyDB.inPedPosId, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetDLOD[1], studyDB.inDGF, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[2], studyDB.inLC1BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[3], studyDB.inLC1BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[4], studyDB.inLC1LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[5], studyDB.inLC1LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[6], studyDB.inLC2BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[7], studyDB.inLC2BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[8], studyDB.inLC2LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[9], studyDB.inLC2LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[10], studyDB.inLC3BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[11], studyDB.inLC3BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[12], studyDB.inLC3LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[13], studyDB.inLC3LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDLOD[14], studyDB.inRegionNo, MYSQL_TYPE_LONG);

  strncpy (studyDB.strGetDLOD, "call GetDLOD (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,@outRegionId,@outMarkerCount,@outLOD)", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtGetDLOD, studyDB.strGetDLOD, strlen (studyDB.strGetDLOD)))
    ERROR("Cannot prepare GetDLOD call statement (%s)", mysql_error(studyDB.connection));
  if (mysql_stmt_bind_param (studyDB.stmtGetDLOD, studyDB.bindGetDLOD))
    ERROR("Cannot bind GetDLOD call statement (%s)", mysql_error(studyDB.connection));

  // Prepare the GetDLOD results call
  studyDB.stmtGetDLODResults = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetDLODResults, 0, sizeof(studyDB.bindGetDLODResults));

  BINDNUMERIC (studyDB.bindGetDLODResults[0], studyDB.outRegionId, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetDLODResults[1], studyDB.outMarkerCount, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetDLODResults[2], studyDB.outLOD, MYSQL_TYPE_DOUBLE);

  strncpy (studyDB.strGetDLODResults, "Select @outRegionId, @outMarkerCount, @outLOD", MAXSTMTLEN-1);
  if (mysql_stmt_prepare (studyDB.stmtGetDLODResults, studyDB.strGetDLODResults, strlen (studyDB.strGetDLODResults)))
    ERROR("Cannot prepare GetDLOD results select statement (%s)", mysql_error(studyDB.connection));
  if (mysql_stmt_bind_result (studyDB.stmtGetDLODResults, studyDB.bindGetDLODResults))
    ERROR("Cannot bind GetDLOD results select statement (%s)", mysql_error(studyDB.connection));

  // Prepare the GetWork call
  studyDB.stmtGetWork = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetWork, 0, sizeof(studyDB.bindGetWork));

  BINDNUMERIC (studyDB.bindGetWork[0], studyDB.inServerId, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetWork[1], studyDB.inLowPosition, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWork[2], studyDB.inHighPosition, MYSQL_TYPE_DOUBLE);

  strncpy (studyDB.strGetWork, "call GetWork (?,?,?,@outPedPosId, @outPedigreeSId, @outPedTraitPosCM, @outLC1MPId, @outLC2MPId, @outLC3MPId)", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtGetWork, studyDB.strGetWork, strlen (studyDB.strGetWork)))
    ERROR("Cannot prepare GetWork call statement (%s)", mysql_error(studyDB.connection));
  if (mysql_stmt_bind_param (studyDB.stmtGetWork, studyDB.bindGetWork))
    ERROR("Cannot bind GetWork call statement (%s)", mysql_error(studyDB.connection));

  // Prepare the GetDTParts call
  studyDB.stmtGetDTParts = mysql_stmt_init (studyDB.connection);
  strncpy (studyDB.strGetDTParts, "call GetDTParts (@outLC1MPId, @outLC2MPId, @outLC3MPId)", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtGetDTParts, studyDB.strGetDTParts, strlen (studyDB.strGetDTParts)))
    ERROR("Cannot prepare GetDTParts call statement (%s)", mysql_error(studyDB.connection));

  // Prepare the GetWork results call
  studyDB.stmtGetWorkResults = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetWorkResults, 0, sizeof(studyDB.bindGetWorkResults));

  BINDNUMERIC (studyDB.bindGetWorkResults[0], studyDB.outPedPosId, MYSQL_TYPE_LONG);
  BINDSTRING (studyDB.bindGetWorkResults[1], studyDB.outPedigreeSId, sizeof (studyDB.outPedigreeSId));
  BINDNUMERIC (studyDB.bindGetWorkResults[2], studyDB.outPedTraitPosCM, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[3], studyDB.outDGF, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[4], studyDB.outLC1BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[5], studyDB.outLC1BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[6], studyDB.outLC1LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[7], studyDB.outLC1LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[8], studyDB.outLC2BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[9], studyDB.outLC2BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[10], studyDB.outLC2LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[11], studyDB.outLC2LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[12], studyDB.outLC3BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[13], studyDB.outLC3BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[14], studyDB.outLC3LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[15], studyDB.outLC3LittlePen, MYSQL_TYPE_DOUBLE);

  strncpy (studyDB.strGetWorkResults, "Select @outPedPosId, @outPedigreeSId, @outPedTraitPosCM, @outDGF,"
	   "@outLC1BigPen, @outLC1BigLittlePen, @outLC1LittleBigPen, @outLC1LittlePen, "
	   "@outLC2BigPen, @outLC2BigLittlePen, @outLC2LittleBigPen, @outLC2LittlePen, "
	   "@outLC3BigPen, @outLC3BigLittlePen, @outLC3LittleBigPen, @outLC3LittlePen", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtGetWorkResults, studyDB.strGetWorkResults, strlen (studyDB.strGetWorkResults)))
    ERROR("Cannot prepare GetWorkResults call statement (%s)", mysql_error(studyDB.connection));
  if (mysql_stmt_bind_param (studyDB.stmtGetWorkResults, studyDB.bindGetWorkResults))
    ERROR("Cannot bind GetWorkResults call statement (%s)", mysql_error(studyDB.connection));

  // Prepare the PutWork call
  studyDB.stmtPutWork = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindPutWork, 0, sizeof(studyDB.bindPutWork));

  BINDNUMERIC (studyDB.bindPutWork[0], studyDB.inMarkerCount, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindPutWork[1], studyDB.inLOD, MYSQL_TYPE_DOUBLE);

  strncpy (studyDB.strPutWork, "call PutWork (@outServerId, @outPedPosId, @outLC1MPId, @outLC2MPId, @outLC3MPId,?,?)", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtPutWork, studyDB.strPutWork, strlen (studyDB.strPutWork)))
    ERROR("Cannot prepare PutWork call statement (%s)", mysql_error(studyDB.connection));
  if (mysql_stmt_bind_param (studyDB.stmtPutWork, studyDB.bindPutWork))
    ERROR("Cannot bind PutWork call statement (%s)", mysql_error(studyDB.connection));

  dBStmtsNotReady = FALSE;
}

long GetPedPosId (char *inPedigreeSId, int inChromosomeNo, double inRefTraitPosCM)
{
  if (dBStmtsNotReady)
    prepareDBStatements ();

  strncpy (studyDB.inPedigreeSId, inPedigreeSId, 16);
  *studyDB.bindGetPedPosId[1].length = strlen(inPedigreeSId);
  studyDB.inChromosomeNo = inChromosomeNo;
  studyDB.inRefTraitPosCM = inRefTraitPosCM;

  if (mysql_stmt_execute (studyDB.stmtGetPedPosId))
    ERROR("Cannot execute PedPosId select statement w/%d, '%s', %d, %G (%s, %s)", 
	  studyDB.inStudyId, inPedigreeSId, inChromosomeNo, inRefTraitPosCM,
	  mysql_stmt_error(studyDB.stmtGetPedPosId), mysql_stmt_sqlstate(studyDB.stmtGetPedPosId));
  if (mysql_stmt_store_result (studyDB.stmtGetPedPosId) != 0)
    ERROR("Cannot retrieve PedPosId select results w/%d, '%s', %d, %G (%s)", 
	  studyDB.inStudyId, inPedigreeSId, inChromosomeNo, inRefTraitPosCM, 
	  mysql_stmt_error(studyDB.stmtGetPedPosId));
  if (mysql_stmt_fetch (studyDB.stmtGetPedPosId) != 0)
    ERROR("Cannot fetch PedPosId select results w/%d, '%s', %d, %G (%s %s)", 
	    studyDB.inStudyId, inPedigreeSId, inChromosomeNo, inRefTraitPosCM, 
	    mysql_stmt_error(studyDB.stmtGetPedPosId), mysql_stmt_sqlstate(studyDB.stmtGetPedPosId));
  return studyDB.outPedPosId;
}

double GetDLOD (int inPedPosId, double inDGF,
	      double inLC1BigPen, double inLC1BigLittlePen, double inLC1LittleBigPen, double inLC1LittlePen,
	      double inLC2BigPen, double inLC2BigLittlePen, double inLC2LittleBigPen, double inLC2LittlePen,
	      double inLC3BigPen, double inLC3BigLittlePen, double inLC3LittleBigPen, double inLC3LittlePen,
	      int inRegionNo)
{
  studyDB.inPedPosId = inPedPosId;
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
  studyDB.inRegionNo = inRegionNo;

  //  printf ("PPId:%u DGF:%G DD:%G Dd:%G dD:%G dd:%G\n", inPedPosId, inDGF, inLC1BigPen, inLC1BigLittlePen, inLC1LittleBigPen, inLC1LittlePen);

  if (mysql_stmt_execute (studyDB.stmtGetDLOD) != 0)
    ERROR("Cannot execute GetDLOD call statement w/%d (%s, %s)", inPedPosId,
	  mysql_stmt_error(studyDB.stmtGetDLOD), mysql_stmt_sqlstate(studyDB.stmtGetDLOD));

  if (mysql_stmt_execute (studyDB.stmtGetDLODResults) != 0)
    ERROR("Cannot execute GetDLOD results select statement (%s, %s)", 
	  mysql_stmt_error(studyDB.stmtGetDLODResults), mysql_stmt_sqlstate(studyDB.stmtGetDLODResults));
  if (mysql_stmt_store_result (studyDB.stmtGetDLODResults) != 0)
    ERROR("Cannot retrieve GetDLOD (%s)", mysql_stmt_error(studyDB.stmtGetDLODResults));
  if (mysql_stmt_fetch (studyDB.stmtGetDLODResults) != 0)
    ERROR("Cannot fetch results (%s)", mysql_stmt_error(studyDB.stmtGetDLODResults));

  if (*studyDB.bindGetDLODResults[2].is_null) {
    //    INFO ("In RegionId %d, LOD is NULL", studyDB.outRegionId);
    return -1LL;
  } else {
    //    INFO ("In RegionId %d, LOD is %G", studyDB.outRegionId, studyDB.outLOD);
    return studyDB.outLOD;
  }
}

int GetWork () {
  return FALSE;
}

#ifdef MAIN

/*

gcc -g -o test databaseSupport.c ../utils/libklvnutls.a -DMAIN -lmysqlclient -I/usr/local/mysql/include -L/usr/local/mysql/lib/ 
./test 19 localhost LOD_dev whv001 foobar

*/

#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {
  
  double myPedPosId;

  // Done for you by directive parsing:

  studyDB.inStudyId = atoi (argv[1]);
  strcpy (studyDB.hostname, argv[2]);
  strcpy (studyDB.dBName, argv[3]);
  strcpy (studyDB.username, argv[4]);
  strcpy (studyDB.password, argv[5]);

  // Annotated calls...

  myPedPosId = GetPedPosId (/* PedigreeSId-> */ "2", /* ChromosomeNo-> */ 44, /* RefTraitPosCM-> */ 5.1);
  GetDLOD (/* PedPosId-> */ myPedPosId, /* DGF-> */ .35,
	   /* LC1BigPen-> */ .71, /* LC1BigLittlePen-> */ .42, /* LC1LittleBigPen-> */ .44, /* LC1LittlePen-> */ .13, 
	   /* LC2BigPen-> */ .71, /* LC2BigLittlePen-> */ .42, /* LC2LittleBigPen-> */ .44, /* LC2LittlePen-> */ .13, 
	   /* LC3BigPen-> */ .71, /* LC3BigLittlePen-> */ .42, /* LC3LittleBigPen-> */ .46, /* LC3LittlePen-> */ .13, 
	   /* RegionNo-> */ 1);

  GetDLOD (myPedPosId, .35, .71, .42, .42, .13, .71, .42, .45, .13, .71, .42, .45, .13, 1);
  myPedPosId = GetPedPosId ("3", 44, 10.43210987);
  GetDLOD (myPedPosId, .3, .71, .42, .44, .13, .71, .42, .42, .13, .71, .42, .42, .13, 1);
  myPedPosId = GetPedPosId ("2", 44, 3.0);
  GetDLOD (myPedPosId, .3, .71, .42, .44, .13, .71, .42, .42, .13, .71, .42, .42, .13, 1);
  GetDLOD (6, .3, .71, .42, .44, .13, -1, 0, 0, 0, -1, 0, 0, 0, 1);
  GetDLOD (7, .3, .71, .42, .44, .13, .71, .42, .42, .13, -1, 0, 0, 0, 1);

  return EXIT_SUCCESS;
}
#endif
