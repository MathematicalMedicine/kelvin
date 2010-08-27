#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>

#include "../utils/sw.h"
#include "../utils/utils.h"

#include "StudyDB.h"

#ifndef MAIN
extern 
#endif
struct StudyDB studyDB;

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

int dBInitNotDone = TRUE, dBStmtsNotReady = TRUE;

void initializeDB () {

  /* Initialize structure for a MySQL connection. */
  if ((studyDB.connection = mysql_init(NULL)) == NULL)
    FATAL("Cannot initialize MySQL (%s)", strerror(errno));

  /* Connect. */
  if (!mysql_real_connect(studyDB.connection, studyDB.dBHostname, studyDB.username, studyDB.password, 
			  NULL, 0, NULL, CLIENT_MULTI_RESULTS /* Important discovery here */))
    ERROR("Cannot connect to MySQL on hostname [%s] as username [%s/%s] (%s)", studyDB.dBHostname, 
	  studyDB.username, studyDB.password, mysql_error(studyDB.connection));

  /* Change database. */
  if (mysql_select_db(studyDB.connection, studyDB.dBName))
    ERROR("Cannot change MySQL db (%s)", mysql_error(studyDB.connection));

  /* Verify our studyId. */
  sprintf (studyDB.strAdhocStatement, "Select Description from Studies where StudyId = %d", studyDB.studyId);
  if (mysql_query (studyDB.connection, studyDB.strAdhocStatement))
    ERROR("Cannot select study information (%s:%s)", studyDB.strAdhocStatement, mysql_error(studyDB.connection));
  if ((studyDB.resultSet = mysql_store_result (studyDB.connection)) == NULL)
    ERROR("Cannot retrieve study information (%s)", mysql_error(studyDB.connection));
  if (mysql_num_rows (studyDB.resultSet) == 0)
    ERROR("Study %d not found", studyDB.studyId);
  else {
    if ((studyDB.row = mysql_fetch_row (studyDB.resultSet)) == NULL)
      ERROR("Cannot fetch study information (%s)", mysql_error(studyDB.connection));
    DIAG (ALTLSERVER, 1, { fprintf (stderr, "Storing/retrieving results under study %d (%s)", studyDB.studyId, studyDB.row[0]);});
  }
  dBInitNotDone = FALSE;
}

void prepareDBStatements () {

  if (dBInitNotDone)
    initializeDB ();

  // Prepare the select for the pedigree/position
  studyDB.stmtGetPedPosId = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetPedPosId, 0, sizeof(studyDB.bindGetPedPosId));

  BINDNUMERIC (studyDB.bindGetPedPosId[0], studyDB.studyId, MYSQL_TYPE_LONG);
  BINDSTRING (studyDB.bindGetPedPosId[1], studyDB.pedigreeSId, sizeof (studyDB.pedigreeSId));
  BINDNUMERIC (studyDB.bindGetPedPosId[2], studyDB.chromosomeNo, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetPedPosId[3], studyDB.refTraitPosCM, MYSQL_TYPE_DOUBLE);

  strncpy (studyDB.strGetPedPosId, "Select PedPosId from PedigreePositions where StudyId = ? AND PedigreeSId = ? AND ChromosomeNo = ? AND RefTraitPosCM = ?", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtGetPedPosId, studyDB.strGetPedPosId, strlen (studyDB.strGetPedPosId)))
    ERROR("Cannot prepare GetPedPosId call statement (%s)", mysql_stmt_error(studyDB.stmtGetPedPosId));
  if (mysql_stmt_bind_param (studyDB.stmtGetPedPosId, studyDB.bindGetPedPosId))
    ERROR("Cannot bind GetPedPosId call statement (%s)", mysql_stmt_error(studyDB.stmtGetPedPosId));

  BINDNUMERIC (studyDB.bindGetPedPosIdResults[0], studyDB.pedPosId, MYSQL_TYPE_LONG);

  if (mysql_stmt_bind_result (studyDB.stmtGetPedPosId, studyDB.bindGetPedPosIdResults))
    ERROR("Cannot bind GetPedPosId results (%s)", mysql_stmt_error(studyDB.stmtGetPedPosId));

  // Prepare the GetDAltL call
  studyDB.stmtGetDAltL = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetDAltL, 0, sizeof(studyDB.bindGetDAltL));

  BINDNUMERIC (studyDB.bindGetDAltL[0], studyDB.pedPosId, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetDAltL[1], studyDB.dGF, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[2], studyDB.lC1BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[3], studyDB.lC1BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[4], studyDB.lC1LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[5], studyDB.lC1LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[6], studyDB.lC2BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[7], studyDB.lC2BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[8], studyDB.lC2LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[9], studyDB.lC2LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[10], studyDB.lC3BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[11], studyDB.lC3BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[12], studyDB.lC3LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[13], studyDB.lC3LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[14], studyDB.regionNo, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetDAltL[15], studyDB.parentRegionNo, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetDAltL[16], studyDB.parentRegionError, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetDAltL[17], studyDB.parentRegionSplitDir, MYSQL_TYPE_LONG);

  strncpy (studyDB.strGetDAltL, "call GetDAltL (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,@outRegionId,@outMarkerCount,@outAltL)", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtGetDAltL, studyDB.strGetDAltL, strlen (studyDB.strGetDAltL)))
    ERROR("Cannot prepare GetDAltL call statement (%s)", mysql_stmt_error(studyDB.stmtGetDAltL));
  if (mysql_stmt_bind_param (studyDB.stmtGetDAltL, studyDB.bindGetDAltL))
    ERROR("Cannot bind GetDAltL call statement (%s)", mysql_stmt_error(studyDB.stmtGetDAltL));

  // Prepare the GetDAltL results call
  studyDB.stmtGetDAltLResults = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetDAltLResults, 0, sizeof(studyDB.bindGetDAltLResults));

  BINDNUMERIC (studyDB.bindGetDAltLResults[0], studyDB.regionId, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetDAltLResults[1], studyDB.markerCount, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetDAltLResults[2], studyDB.lOD, MYSQL_TYPE_DOUBLE);

  strncpy (studyDB.strGetDAltLResults, "Select @outRegionId, @outMarkerCount, @outAltL", MAXSTMTLEN-1);
  if (mysql_stmt_prepare (studyDB.stmtGetDAltLResults, studyDB.strGetDAltLResults, strlen (studyDB.strGetDAltLResults)))
    ERROR("Cannot prepare GetDAltL results select statement (%s)", mysql_stmt_error(studyDB.stmtGetDAltLResults));
  if (mysql_stmt_bind_result (studyDB.stmtGetDAltLResults, studyDB.bindGetDAltLResults))
    ERROR("Cannot bind GetDAltL results select statement (%s)", mysql_stmt_error(studyDB.stmtGetDAltLResults));

  // Prepare the server sign-on
  studyDB.stmtSignOn = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindSignOn, 0, sizeof(studyDB.bindSignOn));

  BINDSTRING (studyDB.bindSignOn[0], studyDB.hostName, sizeof (studyDB.hostName));
  BINDNUMERIC (studyDB.bindSignOn[1], studyDB.processId, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindSignOn[2], studyDB.keepAliveFlag, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindSignOn[3], studyDB.studyId, MYSQL_TYPE_LONG);
  BINDSTRING (studyDB.bindSignOn[4], studyDB.pedigreeRegEx, sizeof (studyDB.pedigreeRegEx));
  BINDNUMERIC (studyDB.bindSignOn[5], studyDB.chromosomeNo, MYSQL_TYPE_LONG);
  BINDSTRING (studyDB.bindSignOn[6], studyDB.algorithm, sizeof (studyDB.algorithm));
  BINDNUMERIC (studyDB.bindSignOn[7], studyDB.markerCount, MYSQL_TYPE_LONG);
  BINDSTRING (studyDB.bindSignOn[8], studyDB.programVersion, sizeof (studyDB.programVersion));

  strncpy (studyDB.strSignOn, "Insert into Servers (HostName, ProcessId, KeepAliveFlag, StudyId, PedigreeRegEx, ChromosomeNo, Algorithm, MarkerCount, ProgramVersion) values (?,?,?,?,?,?,?,?,?)", MAXSTMTLEN-1);
  if (mysql_stmt_prepare (studyDB.stmtSignOn, studyDB.strSignOn, strlen (studyDB.strSignOn)))
    ERROR("Cannot prepare sign-on insert statement (%s)", mysql_stmt_error(studyDB.stmtSignOn));
  if (mysql_stmt_bind_param (studyDB.stmtSignOn, studyDB.bindSignOn))
    ERROR("Cannot bind sign-on insert statement (%s)", mysql_stmt_error(studyDB.stmtSignOn));

  // Prepare the CountWork call
  studyDB.stmtCountWork = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindCountWork, 0, sizeof(studyDB.bindCountWork));

  BINDNUMERIC (studyDB.bindCountWork[0], studyDB.serverId, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindCountWork[1], studyDB.lowPosition, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindCountWork[2], studyDB.highPosition, MYSQL_TYPE_DOUBLE);

  strncpy (studyDB.strCountWork, "call CountWork (?,?,?, @outWorkCount)", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtCountWork, studyDB.strCountWork, strlen (studyDB.strCountWork)))
    ERROR("Cannot prepare CountWork call statement (%s)", mysql_stmt_error(studyDB.stmtCountWork));
  if (mysql_stmt_bind_param (studyDB.stmtCountWork, studyDB.bindCountWork))
    ERROR("Cannot bind CountWork call statement (%s)", mysql_stmt_error(studyDB.stmtCountWork));

  // Prepare the CountWork results call
  studyDB.stmtCountWorkResults = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindCountWorkResults, 0, sizeof(studyDB.bindCountWorkResults));

  BINDNUMERIC (studyDB.bindCountWorkResults[0], studyDB.workCount, MYSQL_TYPE_LONG);

  strncpy (studyDB.strCountWorkResults, "Select @outWorkCount", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtCountWorkResults, studyDB.strCountWorkResults, strlen (studyDB.strCountWorkResults)))
    ERROR("Cannot prepare CountWorkResults call statement (%s)", mysql_stmt_error(studyDB.stmtCountWorkResults));
  if (mysql_stmt_bind_result (studyDB.stmtCountWorkResults, studyDB.bindCountWorkResults))
    ERROR("Cannot bind CountWorkResults call statement (%s)", mysql_stmt_error(studyDB.stmtCountWorkResults));

  // Prepare the GetWork call
  studyDB.stmtGetWork = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetWork, 0, sizeof(studyDB.bindGetWork));

  BINDNUMERIC (studyDB.bindGetWork[0], studyDB.serverId, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindGetWork[1], studyDB.lowPosition, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWork[2], studyDB.highPosition, MYSQL_TYPE_DOUBLE);

  strncpy (studyDB.strGetWork, "call GetWork (?,?,?,@outPedPosId, @outPedigreeSId, @outPedTraitPosCM, @outLC1MPId, @outLC2MPId, @outLC3MPId)", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtGetWork, studyDB.strGetWork, strlen (studyDB.strGetWork)))
    ERROR("Cannot prepare GetWork call statement (%s)", mysql_stmt_error(studyDB.stmtGetWork));
  if (mysql_stmt_bind_param (studyDB.stmtGetWork, studyDB.bindGetWork))
    ERROR("Cannot bind GetWork call statement (%s)", mysql_stmt_error(studyDB.stmtGetWork));

  // Prepare the GetDParts call
  studyDB.stmtGetDParts = mysql_stmt_init (studyDB.connection);
  strncpy (studyDB.strGetDParts, "call GetDParts (@outLC1MPId, @outLC2MPId, @outLC3MPId, @outDGF, "
	   "@outLC1BP, @outLC1BLP, @outLC1LBP, @outLC1LP,"
	   "@outLC2BP, @outLC2BLP, @outLC2LBP, @outLC2LP,"
	   "@outLC3BP, @outLC3BLP, @outLC3LBP, @outLC3LP)", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtGetDParts, studyDB.strGetDParts, strlen (studyDB.strGetDParts)))
    ERROR("Cannot prepare GetDParts call statement (%s)", mysql_stmt_error(studyDB.stmtGetDParts));

  // Prepare the GetWork results call
  studyDB.stmtGetWorkResults = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindGetWorkResults, 0, sizeof(studyDB.bindGetWorkResults));

  BINDNUMERIC (studyDB.bindGetWorkResults[0], studyDB.pedPosId, MYSQL_TYPE_LONG);
  BINDSTRING (studyDB.bindGetWorkResults[1], studyDB.pedigreeSId, sizeof (studyDB.pedigreeSId));
  BINDNUMERIC (studyDB.bindGetWorkResults[2], studyDB.pedTraitPosCM, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[3], studyDB.dGF, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[4], studyDB.lC1BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[5], studyDB.lC1BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[6], studyDB.lC1LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[7], studyDB.lC1LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[8], studyDB.lC2BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[9], studyDB.lC2BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[10], studyDB.lC2LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[11], studyDB.lC2LittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[12], studyDB.lC3BigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[13], studyDB.lC3BigLittlePen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[14], studyDB.lC3LittleBigPen, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindGetWorkResults[15], studyDB.lC3LittlePen, MYSQL_TYPE_DOUBLE);

  strncpy (studyDB.strGetWorkResults, "Select @outPedPosId, @outPedigreeSId, @outPedTraitPosCM, @outDGF,"
	   "@outLC1BP, @outLC1BLP, @outLC1LBP, @outLC1LP, "
	   "@outLC2BP, @outLC2BLP, @outLC2LBP, @outLC2LP, "
	   "@outLC3BP, @outLC3BLP, @outLC3LBP, @outLC3LP", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtGetWorkResults, studyDB.strGetWorkResults, strlen (studyDB.strGetWorkResults)))
    ERROR("Cannot prepare GetWorkResults call statement (%s)", mysql_stmt_error(studyDB.stmtGetWorkResults));
  if (mysql_stmt_bind_result (studyDB.stmtGetWorkResults, studyDB.bindGetWorkResults))
    ERROR("Cannot bind GetWorkResults call statement (%s)", mysql_stmt_error(studyDB.stmtGetWorkResults));

  // Prepare the PutWork call
  studyDB.stmtPutWork = mysql_stmt_init (studyDB.connection);
  memset (studyDB.bindPutWork, 0, sizeof(studyDB.bindPutWork));

  BINDNUMERIC (studyDB.bindPutWork[0], studyDB.serverId, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindPutWork[1], studyDB.markerCount, MYSQL_TYPE_LONG);
  BINDNUMERIC (studyDB.bindPutWork[2], studyDB.lOD, MYSQL_TYPE_DOUBLE);
  BINDNUMERIC (studyDB.bindPutWork[3], studyDB.runtimeCostSec, MYSQL_TYPE_LONG);

  strncpy (studyDB.strPutWork, "call PutWork (?, @outPedPosId, @outLC1MPId, @outLC2MPId, @outLC3MPId,?,?,?)", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtPutWork, studyDB.strPutWork, strlen (studyDB.strPutWork)))
    ERROR("Cannot prepare PutWork call statement (%s)", mysql_stmt_error(studyDB.stmtPutWork));
  if (mysql_stmt_bind_param (studyDB.stmtPutWork, studyDB.bindPutWork))
    ERROR("Cannot bind PutWork call statement (%s)", mysql_stmt_error(studyDB.stmtPutWork));

  dBStmtsNotReady = FALSE;
}

long GetPedPosId (char *pedigreeSId, int chromosomeNo, double refTraitPosCM)
{
  if (dBStmtsNotReady)
    prepareDBStatements ();

  strncpy (studyDB.pedigreeSId, pedigreeSId, 16);
  *studyDB.bindGetPedPosId[1].length = strlen(pedigreeSId);
  studyDB.chromosomeNo = chromosomeNo;
  studyDB.refTraitPosCM = refTraitPosCM;

  if (mysql_stmt_execute (studyDB.stmtGetPedPosId))
    ERROR("Cannot execute PedPosId select statement w/%d, '%s', %d, %G (%s, %s)", 
	  studyDB.studyId, pedigreeSId, chromosomeNo, refTraitPosCM,
	  mysql_stmt_error(studyDB.stmtGetPedPosId), mysql_stmt_sqlstate(studyDB.stmtGetPedPosId));
  if (mysql_stmt_store_result (studyDB.stmtGetPedPosId) != 0)
    ERROR("Cannot retrieve PedPosId select results w/%d, '%s', %d, %G (%s)", 
	  studyDB.studyId, pedigreeSId, chromosomeNo, refTraitPosCM, 
	  mysql_stmt_error(studyDB.stmtGetPedPosId));
  if (mysql_stmt_fetch (studyDB.stmtGetPedPosId) != 0)
    ERROR("Cannot fetch PedPosId select results w/%d, '%s', %d, %G (%s %s)", 
	    studyDB.studyId, pedigreeSId, chromosomeNo, refTraitPosCM, 
	    mysql_stmt_error(studyDB.stmtGetPedPosId), mysql_stmt_sqlstate(studyDB.stmtGetPedPosId));
  return studyDB.pedPosId;
}

double GetDAltL (int pedPosId, double dGF,
		double lC1BigPen, double lC1BigLittlePen, double lC1LittleBigPen, double lC1LittlePen,
		double lC2BigPen, double lC2BigLittlePen, double lC2LittleBigPen, double lC2LittlePen,
		double lC3BigPen, double lC3BigLittlePen, double lC3LittleBigPen, double lC3LittlePen,
		int regionNo, int parentRegionNo, double parentRegionError, int parentRegionSplitDir)
{
  studyDB.pedPosId = pedPosId;
  studyDB.dGF = dGF;
  studyDB.lC1BigPen = lC1BigPen;
  studyDB.lC1BigLittlePen = lC1BigLittlePen;
  studyDB.lC1LittleBigPen = lC1LittleBigPen;
  studyDB.lC1LittlePen = lC1LittlePen;
  studyDB.lC2BigPen = lC2BigPen;
  studyDB.lC2BigLittlePen = lC2BigLittlePen;
  studyDB.lC2LittleBigPen = lC2LittleBigPen;
  studyDB.lC2LittlePen = lC2LittlePen;
  studyDB.lC3BigPen = lC3BigPen;
  studyDB.lC3BigLittlePen = lC3BigLittlePen;
  studyDB.lC3LittleBigPen = lC3LittleBigPen;
  studyDB.lC3LittlePen = lC3LittlePen;
  studyDB.regionNo = regionNo;
  studyDB.parentRegionNo = parentRegionNo;
  studyDB.parentRegionError = parentRegionError;
  studyDB.parentRegionSplitDir = parentRegionSplitDir;

  if (mysql_stmt_execute (studyDB.stmtGetDAltL) != 0)
    ERROR("Cannot execute GetDAltL call statement w/%d (%s, %s)", pedPosId,
	  mysql_stmt_error(studyDB.stmtGetDAltL), mysql_stmt_sqlstate(studyDB.stmtGetDAltL));

  if (mysql_stmt_execute (studyDB.stmtGetDAltLResults) != 0)
    ERROR("Cannot execute GetDAltL results select statement (%s, %s)", 
	  mysql_stmt_error(studyDB.stmtGetDAltLResults), mysql_stmt_sqlstate(studyDB.stmtGetDAltLResults));
  if (mysql_stmt_store_result (studyDB.stmtGetDAltLResults) != 0)
    ERROR("Cannot retrieve GetDAltL (%s)", mysql_stmt_error(studyDB.stmtGetDAltLResults));
  if (mysql_stmt_fetch (studyDB.stmtGetDAltLResults) != 0)
    ERROR("Cannot fetch results (%s)", mysql_stmt_error(studyDB.stmtGetDAltLResults));

  if (*studyDB.bindGetDAltLResults[2].is_null) {
    DIAG (ALTLSERVER, 1, { fprintf (stderr, "In RegionId %d, AltL is NULL", studyDB.regionId);});
    return -1LL;
  } else {
    DIAG (ALTLSERVER, 1, { fprintf (stderr, "In RegionId %d, AltL is %G", studyDB.regionId, studyDB.lOD);});
    return studyDB.lOD;
  }
}

void SignOn (int chromosomeNo, char *algorithm, int markerCount, char *programVersion) {

  if (dBStmtsNotReady)
    prepareDBStatements ();

  // studyId and pedigreeRegEx are already set, but location info is needed and pedigreeRegEx needs size set
  gethostname (studyDB.hostName, 32);
  *studyDB.bindSignOn[0].length = strlen(studyDB.hostName);
  studyDB.processId = getpid();
  studyDB.keepAliveFlag = 1;
  *studyDB.bindSignOn[4].length = strlen(studyDB.pedigreeRegEx);
  studyDB.chromosomeNo = chromosomeNo;
  strncpy (studyDB.algorithm, algorithm, 2);
  *studyDB.bindSignOn[6].length = strlen(algorithm);
  studyDB.markerCount = markerCount;
  strncpy (studyDB.programVersion, programVersion, 32);
  *studyDB.bindSignOn[8].length = strlen(programVersion);

  if (mysql_stmt_execute (studyDB.stmtSignOn))
    ERROR("Cannot execute sign-on insert statement w/%d, '%s', %d, '%s', %d, '%s' (%s, %s)", 
	  studyDB.studyId, studyDB.pedigreeRegEx, chromosomeNo, algorithm, markerCount, programVersion,
	  mysql_stmt_error(studyDB.stmtSignOn), mysql_stmt_sqlstate(studyDB.stmtSignOn));

  /* Sigh...we really do need the serverId for logging, otherwise I'd let it be yet another
     temporary variable that flows between procedures under our connection context. */

  /* Verify our studyId. */
  sprintf (studyDB.strAdhocStatement, "Select LAST_INSERT_ID()");
  if (mysql_query (studyDB.connection, studyDB.strAdhocStatement))
    ERROR("Cannot select LAST_INSERT_ID() (serverId) (%s:%s)", studyDB.strAdhocStatement, mysql_error(studyDB.connection));
  if ((studyDB.resultSet = mysql_store_result (studyDB.connection)) == NULL)
    ERROR("Cannot retrieve LAST_SERVER_ID() (serverId) (%s)", mysql_error(studyDB.connection));
  if (mysql_num_rows (studyDB.resultSet) == 0)
    ERROR("serverId not found");
  else {
    if ((studyDB.row = mysql_fetch_row (studyDB.resultSet)) == NULL)
      ERROR("Cannot fetch LAST_SERVER_ID() (serverId) (%s)", mysql_error(studyDB.connection));
    studyDB.serverId = atoi(studyDB.row[0]);
    DIAG (ALTLSERVER, 1, { fprintf (stderr, "Signed on as serverId %d", studyDB.serverId);});
  }

}

int CountWork (double lowPosition, double highPosition)
{
  // serverId is already set
  studyDB.lowPosition = lowPosition;
  studyDB.highPosition = highPosition;

  // CountWork
  while (1) {
    if (mysql_stmt_execute (studyDB.stmtCountWork) != 0) {
      if (strcmp (mysql_stmt_sqlstate(studyDB.stmtCountWork), "40001") != 0) {
	ERROR("Cannot execute Count statement w/%G, %G, (%s, %s)", 
	      lowPosition, highPosition,
	      mysql_stmt_error(studyDB.stmtCountWork), mysql_stmt_sqlstate(studyDB.stmtCountWork));
      } else {
	swLogProgress(5, 0, "Retrying deadlock in 1 second");
	sleep(1);
	continue;
      }
    }
    break;
  }    
  if (mysql_stmt_execute (studyDB.stmtCountWorkResults))
    ERROR("Cannot execute CountWorkResults statement (%s, %s)", 
	  mysql_stmt_error(studyDB.stmtCountWorkResults), mysql_stmt_sqlstate(studyDB.stmtCountWorkResults));
  if (mysql_stmt_store_result (studyDB.stmtCountWorkResults) != 0)
    ERROR("Cannot retrieve CountWork results (%s)", mysql_stmt_error(studyDB.stmtCountWorkResults));
  if (mysql_stmt_fetch (studyDB.stmtCountWorkResults) != 0)
    ERROR("Cannot fetch results (%s)", mysql_stmt_error(studyDB.stmtCountWorkResults));

  return studyDB.workCount;
}

int GetDWork (double lowPosition, double highPosition, double *pedTraitPosCM, char *pedigreeSId, double *dGF,
	      double *lC1BigPen, double *lC1BigLittlePen, double *lC1LittleBigPen, double *lC1LittlePen,
	      double *lC2BigPen, double *lC2BigLittlePen, double *lC2LittleBigPen, double *lC2LittlePen,
	      double *lC3BigPen, double *lC3BigLittlePen, double *lC3LittleBigPen, double *lC3LittlePen)
{

  // serverId is already set
  studyDB.lowPosition = lowPosition;
  studyDB.highPosition = highPosition;

  // GetWork
  while (1) {
    if (mysql_stmt_execute (studyDB.stmtGetWork) != 0) {
      if ((strcmp (mysql_stmt_sqlstate(studyDB.stmtGetWork), "40001") != 0) &&
	  (strcmp (mysql_stmt_sqlstate(studyDB.stmtGetWork), "HY000") != 0)) {
	ERROR("Cannot execute Get statement w/%G, %G, (%s, %s)", 
	      lowPosition, highPosition,
	      mysql_stmt_error(studyDB.stmtGetWork), mysql_stmt_sqlstate(studyDB.stmtGetWork));
      } else {
	swLogProgress(5, 0, "Retrying deadlock in 1 second");
	sleep(1);
	continue;
      }
    }
    break;
  }    

  // Get DT parts - this uses the temporary variables in our session to get parameters from GetWork.
  if (mysql_stmt_execute (studyDB.stmtGetDParts))
    ERROR("Cannot execute GetDParts statement (%s, %s)", 
	  mysql_stmt_error(studyDB.stmtGetDParts), mysql_stmt_sqlstate(studyDB.stmtGetDParts));

  // GetWork results - a slew of parameters
  studyDB.pedPosId = 0;
  studyDB.pedTraitPosCM = studyDB.dGF = studyDB.lC1BigPen = studyDB.lC1BigLittlePen = 
    studyDB.lC1LittleBigPen = studyDB.lC1LittlePen = 0.0;
  if (mysql_stmt_execute (studyDB.stmtGetWorkResults))
    ERROR("Cannot execute GetWorkResults statement (%s, %s)", 
	  mysql_stmt_error(studyDB.stmtGetWorkResults), mysql_stmt_sqlstate(studyDB.stmtGetWorkResults));
  if (mysql_stmt_store_result (studyDB.stmtGetWorkResults) != 0)
    ERROR("Cannot retrieve GetWork results (%s)", mysql_stmt_error(studyDB.stmtGetWorkResults));
  if (mysql_stmt_fetch (studyDB.stmtGetWorkResults) != 0)
    ERROR("Cannot fetch results (%s)", mysql_stmt_error(studyDB.stmtGetWorkResults));
  if (*studyDB.bindGetWorkResults[0].is_null) {
    DIAG (ALTLSERVER, 1, { fprintf (stderr, "No more work! (bindGetWorkResults)");});
    return FALSE;
  } else {
    DIAG (ALTLSERVER, 1, { \
	fprintf (stderr, "Got work for PedPosId %d: pedigree %s, position %f, DGF %G, DD %G, Dd %G, dD %G, dd %G", \
		 studyDB.pedPosId, studyDB.pedigreeSId, studyDB.pedTraitPosCM, studyDB.dGF, \
		 studyDB.lC1BigPen, studyDB.lC1BigLittlePen, studyDB.lC1LittleBigPen, studyDB.lC1LittlePen);});
    strcpy (pedigreeSId, studyDB.pedigreeSId);
    *pedTraitPosCM = studyDB.pedTraitPosCM;
    *dGF = studyDB.dGF;
    *lC1BigPen = studyDB.lC1BigPen;
    *lC1BigLittlePen = studyDB.lC1BigLittlePen;
    *lC1LittleBigPen = studyDB.lC1LittleBigPen;
    *lC1LittlePen = studyDB.lC1LittlePen;
    *lC2BigPen = studyDB.lC2BigPen;
    *lC2BigLittlePen = studyDB.lC2BigLittlePen;
    *lC2LittleBigPen = studyDB.lC2LittleBigPen;
    *lC2LittlePen = studyDB.lC2LittlePen;
    *lC3BigPen = studyDB.lC3BigPen;
    *lC3BigLittlePen = studyDB.lC3BigLittlePen;
    *lC3LittleBigPen = studyDB.lC3LittleBigPen;
    *lC3LittlePen = studyDB.lC3LittlePen;
    return TRUE;
  }
}

void PutWork (int markerCount, double lOD, int runtimeCostSec)
{
  
  // serverId is already set
  studyDB.markerCount = markerCount;
  studyDB.lOD = lOD;
  studyDB.runtimeCostSec = runtimeCostSec;

  // PutWork
  while (1) {
    if (mysql_stmt_execute (studyDB.stmtPutWork) != 0) {
      if ((strcmp (mysql_stmt_sqlstate(studyDB.stmtPutWork), "40001") != 0) &&
	  (strcmp (mysql_stmt_sqlstate(studyDB.stmtPutWork), "HY000") != 0)) {
	ERROR("Cannot execute Put statement w/%d, %G, (%s, %s)", 
	      markerCount, lOD,
	      mysql_stmt_error(studyDB.stmtPutWork), mysql_stmt_sqlstate(studyDB.stmtPutWork));
      } else {
	swLogProgress(5, 0, "Retrying deadlock in 1 second");
	sleep(1);
	continue;
      }
    }
    break;
  }    
  DIAG (ALTLSERVER, 1, { fprintf (stderr, "Put work stored AltL of %.8g\n", lOD);});

}

#ifdef MAIN

/*

gcc -g -o test databaseSupport.c ../utils/libklvnutls.a -DMAIN -lmysqlclient -lpthread -I/usr/include/mysql -L/usr/lib64/mysql/
./test 1 mcclintock dev adminDev foobar '.*'


*/

#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {
  
  double pedPosId, pedTraitPosCM, dGF, lC1DD, lC1Dd, lC1dD, lC1dd, lC2DD, lC2Dd, lC2dD, lC2dd, lC3DD, lC3Dd, lC3dD, lC3dd, lOD;
  char pedigreeSId[33];

  // Done for you by directive parsing:

  studyDB.studyId = atoi (argv[1]);
  strcpy (studyDB.dBHostname, argv[2]);
  strcpy (studyDB.dBName, argv[3]);
  strcpy (studyDB.username, argv[4]);
  strcpy (studyDB.password, argv[5]);
  strcpy (studyDB.pedigreeRegEx, argv[6]);

  // Annotated calls...

  pedPosId = GetPedPosId (/* PedigreeSId-> */ "2", /* ChromosomeNo-> */ 40, /* RefTraitPosCM-> */ 5.1);
  GetDAltL (/* PedPosId-> */ pedPosId, /* DGF-> */ .35,
	   /* LC1DD-> */ .71, /* LC1Dd-> */ .42, /* LC1dD-> */ .44, /* LC1dd-> */ .13, 
	   /* LC2DD-> */ .71, /* LC2Dd-> */ .42, /* LC2dD-> */ .44, /* LC2dd-> */ .13, 
	   /* LC3DD-> */ .71, /* LC3Dd-> */ .42, /* LC3dD-> */ .46, /* LC3dd-> */ .13, 
	   /* RegionNo-> */ 1);

  GetDAltL (pedPosId, .35, .71, .42, .42, .13, .71, .42, .45, .13, .71, .42, .45, .13, 1);
  pedPosId = GetPedPosId ("3", 40, 10.43210987);
  GetDAltL (pedPosId, .3, .71, .42, .44, .13, .71, .42, .42, .13, .71, .42, .42, .13, 1);
  pedPosId = GetPedPosId ("2", 40, 3.0);
  GetDAltL (pedPosId, .3, .71, .42, .44, .13, .71, .42, .42, .13, .71, .42, .42, .13, 1);
  GetDAltL (6, .3, .71, .42, .44, .13, -1, 0, 0, 0, -1, 0, 0, 0, 1);
  GetDAltL (7, .3, .71, .42, .44, .13, .71, .42, .42, .13, -1, 0, 0, 0, 1);

  SignOn (40, "ES", 4, "Test driver");

  while (GetDWork (0, 10, &pedTraitPosCM, pedigreeSId, &dGF, &lC1DD, &lC1Dd, &lC1dD, &lC1dd,
		   &lC2DD, &lC2Dd, &lC2dD, &lC2dd, &lC3DD, &lC3Dd, &lC3dD, &lC3dd)) {
    lOD = ((double)(rand() % 9999)) / 1000.0;
    printf ("Trying AltL of %G\n", lOD);
    PutWork (5, lOD);
  }

  return EXIT_SUCCESS;
}
#endif
